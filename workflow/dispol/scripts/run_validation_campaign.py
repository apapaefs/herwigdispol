#!/usr/bin/env python3.10
"""
Run a full DIS validation campaign and collect the outputs.

This script orchestrates the standard Herwig validation runs for the
neutral-current DIS comparison:

  - setups:   GAMMA, Z, ALL
  - orders:   LO, POSNLO, NEGNLO
  - helicities:
      * GAMMA: 00, PP, PM
      * Z/ALL: 00, PP, PM, MP, MM

After the runs finish, it invokes the existing extraction scripts to produce
the validation summaries, merges shard YODA files back into one logical-run
YODA where possible, and stores a campaign manifest that includes both the raw
generated .out/.log/.yoda files and any merged YODA outputs.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import secrets
import shlex
import shutil
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, TypeVar

import yoda

from analyze_DIS_polarized import (
    SCALE_VARIATION_FACTORS,
    build_dis_polarized_objects,
    build_plot_scatter_objects,
    build_scale_envelope_plot_yoda,
    clone_objects_with_variation,
    filter_compatible_variation_plot_objects,
    harmonize_plot_yoda_bin_grids,
    sanitize_plot_yoda,
    write_yoda_gz as write_analysis_yoda_gz,
)
from extract_dis_out_results import parse_measurement
from poldis_top_to_yoda import convert_topdrawer_to_yoda
from rivet_scale_plot_postprocess import (
    rewrite_no_scale_ratio_plot_scripts,
    rewrite_scale_envelope_plot_scripts,
)

try:
    from prettytable import PrettyTable
except ImportError:
    PrettyTable = None


SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
REPO_ROOT = WORKFLOW_DIR.parent.parent
DEFAULT_BASE_DIR = WORKFLOW_DIR
GAMMA_HELICITIES = ("00", "PP", "PM")
FULL_HELICITIES = ("00", "PP", "PM", "MP", "MM")
SETUP_ORDER = ("GAMMA", "Z", "ALL")
PS_SETUP_ORDER = ("SPINVAL", "SPINCOMP", "SPINHAD")
SUPPORTED_SETUPS = SETUP_ORDER + ("CC",) + PS_SETUP_ORDER
ORDER_ORDER = ("LO", "POSNLO", "NEGNLO")
PROGRESS_MARKER_RE = re.compile(r"event>\s+(?P<current>init|\d+)(?:\s+(?P<total>\d+)|/(?P<total_alt>\d+))")

POWHEG_OPTION_LINE_PREFIXES = (
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:",
    "set /Herwig/Shower/ShowerHandler:HardEmission",
    "set /Herwig/Shower/ShowerHandler:LimitEmissions",
    "set /Herwig/Shower/ShowerHandler:UseConstituentMasses",
    "set /Herwig/Shower/KinematicsReconstructor:IFReconPOWHEG",
    "insert /Herwig/Analysis/Rivet:Analyses",
)
SCALE_VARIATION_ORDER = tuple(SCALE_VARIATION_FACTORS.keys())
SCALE_VARIATION_STEM_SUFFIXES = {
    "nominal": "",
    "ScaleFactorDown": "-SCALEDOWN",
    "ScaleFactorUp": "-SCALEUP",
}
PS_SETUP_FAMILIES = {
    "SPINVAL": ("RIVETPS-SPIN",),
    "SPINCOMP": ("RIVETPS-SPIN", "RIVETPS-NOSPIN", "RIVETPS-NOSPIN-UNPOL"),
    "SPINHAD": ("RIVETPS-SPIN", "RIVETPS-NOSPIN", "RIVETPS-NOSPIN-UNPOL"),
}
PS_FAMILY_LABELS = {
    "RIVETPS-SPIN": "Full",
    "RIVETPS-NOSPIN": "Born-Only",
    "RIVETPS-NOSPIN-UNPOL": "None",
}
PS_FAMILY_STYLES = {
    "RIVETPS-SPIN": ("red", "solid"),
    "RIVETPS-NOSPIN": ("blue", "dashed"),
    "RIVETPS-NOSPIN-UNPOL": ("green", "dotted"),
}
PS_ANALYSIS_NAME = "MC_DIS_PS"


def analyses_dir(base_dir: Path) -> Path:
    preferred = REPO_ROOT / "analyses" / "rivet" / "dis"
    return preferred


def cards_subdir(base_dir: Path) -> Path:
    candidate = base_dir / "cards"
    return candidate if candidate.exists() else base_dir


def card_relative_path(base_dir: Path, filename: str | Path) -> Path:
    path = Path(filename)
    if path.is_absolute():
        return path
    if path.parts and path.parts[0] == "cards":
        return path
    if cards_subdir(base_dir) != base_dir:
        return Path("cards") / path
    return path


def card_path(base_dir: Path, filename: str | Path) -> Path:
    path = card_relative_path(base_dir, filename)
    if path.is_absolute():
        return path
    return base_dir / path


def script_path(name: str) -> Path:
    return SCRIPT_DIR / name


def _analysis_component_matches(component: str, analysis_name: str) -> bool:
    if component == analysis_name:
        return True
    if not component.startswith(analysis_name) or len(component) == len(analysis_name):
        return False
    return component[len(analysis_name)] in ":;[("


def _yoda_path_matches_analysis_label(path: str, analysis_name: str, label: str) -> bool:
    parts = [part for part in path.split("/") if part]
    if not parts:
        return False
    if parts[0] == "RAW":
        if len(parts) < 3:
            return False
        return _analysis_component_matches(parts[1], analysis_name) and parts[2] == label
    if len(parts) < 2:
        return False
    return _analysis_component_matches(parts[0], analysis_name) and parts[1] == label


def print_stage(message: str) -> None:
    print(f"[stage] {message}", flush=True)


T = TypeVar("T")
U = TypeVar("U")


def run_parallel_ordered(
    items: Sequence[T],
    worker_fn: Callable[[T], U],
    max_workers: int,
) -> List[U]:
    if not items:
        return []
    worker_count = max(1, min(max_workers, len(items)))
    if worker_count == 1:
        return [worker_fn(item) for item in items]

    results: List[Optional[U]] = [None] * len(items)
    with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
        future_to_index = {
            executor.submit(worker_fn, item): index
            for index, item in enumerate(items)
        }
        for future in concurrent.futures.as_completed(future_to_index):
            index = future_to_index[future]
            results[index] = future.result()
    return [result for result in results if result is not None]


def extract_powheg_option_lines(in_path: Path) -> List[str]:
    if not in_path.exists():
        return []
    lines: List[str] = []
    for raw_line in in_path.read_text().splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if any(line.startswith(prefix) for prefix in POWHEG_OPTION_LINE_PREFIXES):
            lines.append(line)
    return lines


def write_powheg_options_summary(
    campaign_dir: Path,
    base_dir: Path,
    jobs: Sequence[JobSpec],
    analysis_variant: str,
) -> Optional[str]:
    unique_cards = sorted({spec.in_file for spec in jobs})
    output_path = campaign_dir / "powheg-options.txt"
    sections: List[str] = []
    for in_file in unique_cards:
        in_path = card_path(base_dir, in_file)
        option_lines = extract_powheg_option_lines(in_path)
        if not option_lines:
            continue
        sections.append(f"[{in_file}]")
        sections.extend(option_lines)
        sections.append("")

    variants = sorted({job.analysis_variant or analysis_variant for job in jobs if job.analysis_variant or analysis_variant})
    header = [
        f"# POWHEG-related card options for campaign '{campaign_dir.name}'",
        f"# analysis_variant: {analysis_variant or 'default'}",
        f"# analysis_variants: {', '.join(variants) if variants else 'default'}",
        f"# base_dir: {base_dir}",
        "",
    ]
    if not sections:
        header.append("# No matching POWHEG-related option lines were found in the selected .in files.")
        header.append("")
    output_path.write_text("\n".join(header + sections))
    return str(output_path)


@dataclass(frozen=True)
class JobSpec:
    setup: str
    order: str
    helicity: str
    stem: str
    run_file: str
    in_file: str
    events: int
    analysis_variant: str = ""
    scale_variation: str = "nominal"
    scale_factor: float = 1.0
    source_in_file: str = ""


@dataclass(frozen=True)
class ShardSpec:
    job: JobSpec
    shard_index: int
    shard_count: int
    tag: str
    seed: int
    events: int
    rerun_parent_tag: str = ""

    @property
    def stem(self) -> str:
        return self.job.stem


@dataclass
class JobResult:
    spec: ShardSpec
    command: List[str]
    returncode: Optional[int] = None
    started_at: Optional[float] = None
    ended_at: Optional[float] = None
    launcher_log: Optional[str] = None
    out_files: List[str] = field(default_factory=list)
    log_files: List[str] = field(default_factory=list)
    yoda_files: List[str] = field(default_factory=list)
    prepared: bool = False
    prepare_command: Optional[List[str]] = None
    prepare_returncode: Optional[int] = None

    @property
    def duration_s(self) -> Optional[float]:
        if self.started_at is None or self.ended_at is None:
            return None
        return self.ended_at - self.started_at


@dataclass
class YODAMergeResult:
    logical_run: str
    source_files: List[str]
    output_file: Optional[str]
    command: List[str]
    returncode: int
    copied: bool = False
    skipped: bool = False
    message: Optional[str] = None


@dataclass
class YODANLOResult:
    logical_run: str
    pos_file: Optional[str]
    neg_file: Optional[str]
    output_file: Optional[str]
    command: List[str]
    returncode: int
    skipped: bool = False
    message: Optional[str] = None


@dataclass
class RawPOWHEGYODAResult:
    logical_run: str
    channel: Optional[str]
    variant: Optional[str]
    source_log: Optional[str]
    source_out: Optional[str]
    output_file: Optional[str]
    summary_csv: Optional[str]
    command: List[str]
    returncode: int
    skipped: bool = False
    message: Optional[str] = None


def render_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    aligns: Optional[Sequence[str]] = None,
) -> List[str]:
    if not rows:
        return []

    string_rows = [[str(cell) for cell in row] for row in rows]
    alignments = list(aligns) if aligns is not None else ["l"] * len(headers)

    if PrettyTable is not None:
        table = PrettyTable()
        table.field_names = list(headers)
        for header, align in zip(headers, alignments):
            table.align[header] = "r" if align == "r" else "l"
        for row in string_rows:
            table.add_row(row)
        return table.get_string().splitlines()

    widths = []
    for idx, header in enumerate(headers):
        cell_width = max((len(row[idx]) for row in string_rows), default=0)
        widths.append(max(len(header), cell_width))

    def format_row(row: Sequence[str], is_header: bool = False) -> str:
        cells: List[str] = []
        for idx, cell in enumerate(row):
            align = alignments[idx] if idx < len(alignments) else "l"
            width = widths[idx]
            if align == "r" and not is_header:
                cells.append(cell.rjust(width))
            else:
                cells.append(cell.ljust(width))
        return "| " + " | ".join(cells) + " |"

    separator = "+-" + "-+-".join("-" * width for width in widths) + "-+"
    output = [separator, format_row(list(headers), is_header=True), separator]
    output.extend(format_row(row) for row in string_rows)
    output.append(separator)
    return output


def analysis_variant_from_args(args: argparse.Namespace) -> str:
    rivet = bool(getattr(args, "rivet", False))
    rivetfo = bool(getattr(args, "rivetfo", False))
    selected = int(rivet) + int(rivetfo)
    if selected > 1:
        raise ValueError("Use at most one of --rivet and --rivetfo.")
    if rivetfo:
        return "RIVETFO"
    if rivet:
        return "RIVET"
    return ""


def analysis_variant_args(analysis_variant: str) -> List[str]:
    if analysis_variant == "RIVET":
        return ["--rivet"]
    if analysis_variant == "RIVETFO":
        return ["--rivetfo"]
    return []


def resolved_analysis_variant_from_args(args: argparse.Namespace) -> str:
    resolved = getattr(args, "resolved_analysis_variant", None)
    if isinstance(resolved, str) and resolved:
        return resolved
    return analysis_variant_from_args(args)


def use_rivet_runtime(analysis_variant: str) -> bool:
    return analysis_variant in {"RIVET", "RIVETFO", "RIVETPS"} or analysis_variant.startswith("RIVETPS-")


def use_lo_runtime(analysis_variant: str) -> bool:
    return analysis_variant not in {"RIVETFO", "RIVETPS"} and not analysis_variant.startswith("RIVETPS-")


def use_raw_powheg_runtime(analysis_variant: str, enabled: bool) -> bool:
    return enabled and analysis_variant == "RIVETFO"


def is_cc_setup(setup: str) -> bool:
    return str(setup).upper() == "CC"


def is_ps_setup(setup: str) -> bool:
    return str(setup).upper() in PS_SETUP_ORDER


def contains_ps_setups(setups: Sequence[str]) -> bool:
    return any(is_ps_setup(setup) for setup in setups)


def contains_legacy_setups(setups: Sequence[str]) -> bool:
    return any(not is_ps_setup(setup) for setup in setups)


def validate_setup_family(setups: Sequence[str]) -> None:
    if contains_ps_setups(setups) and contains_legacy_setups(setups):
        raise ValueError("Do not mix SPINVAL/SPINCOMP/SPINHAD with GAMMA/Z/ALL/CC in the same invocation.")


def ps_families_for_setup(setup: str) -> Sequence[str]:
    try:
        return PS_SETUP_FAMILIES[setup]
    except KeyError as exc:
        raise ValueError(f"Unsupported PS setup {setup!r}") from exc


def ps_analysis_variants_by_setup(setups: Sequence[str]) -> Dict[str, Sequence[str]]:
    return {setup: ps_families_for_setup(setup) for setup in setups if is_ps_setup(setup)}


def supports_ps_family_comparison(setup: str) -> bool:
    return is_ps_setup(setup) and len(ps_families_for_setup(setup)) > 1


def multi_family_ps_setups(setups: Sequence[str]) -> List[str]:
    return [setup for setup in normalize_campaign_setups(setups) if supports_ps_family_comparison(setup)]


def has_ps_hadronization_toggle_inputs(setups: Sequence[str]) -> bool:
    selected = set(normalize_campaign_setups(setups))
    return {"SPINCOMP", "SPINHAD"}.issubset(selected)


def ps_family_label(analysis_variant: str) -> str:
    return PS_FAMILY_LABELS.get(analysis_variant, analysis_variant)


def herwig_legend_label(analysis_variant: str, setup: str = "") -> str:
    if analysis_variant in PS_FAMILY_LABELS:
        return ps_family_label(analysis_variant)
    if analysis_variant in {"RIVET", "RIVETFO"}:
        return "HERWIG POWHEG CC" if is_cc_setup(setup) else "HERWIG POWHEG NC"
    return "HERWIG NLO CC" if is_cc_setup(setup) else "HERWIG NLO NC"


def herwig_rivetplot_title(analysis_variant: str, setup: str = "") -> str:
    if analysis_variant in PS_FAMILY_LABELS:
        return ps_family_label(analysis_variant)
    if analysis_variant in {"RIVET", "RIVETFO"}:
        return "HERWIG POWHEG CC" if is_cc_setup(setup) else "HERWIG POWHEG NC"
    return "HERWIG NLO CC" if is_cc_setup(setup) else "HERWIG NLO NC"


def apply_setup_style_annotation(objects: Dict[str, object], setup: str) -> None:
    if not is_cc_setup(setup):
        return

    line_color = "blue"
    band_color = "blue"
    band_opacity = [0.35]
    line_style = "solid"

    for obj in objects.values():
        for key, value in (
            ("LineColor", line_color),
            ("MarkerColor", line_color),
            ("ErrorBandColor", band_color),
            ("ErrorBandOpacity", band_opacity),
            ("LineStyle", line_style),
        ):
            try:
                obj.setAnnotation(key, value)
            except Exception:
                pass


def apply_variant_style_annotation(objects: Dict[str, object], analysis_variant: str) -> None:
    style = PS_FAMILY_STYLES.get(analysis_variant)
    if style is None:
        return

    line_color, line_style = style
    apply_line_style_annotation(objects, line_color=line_color, line_style=line_style)


def apply_line_style_annotation(
    objects: Dict[str, object],
    line_color: str,
    line_style: str,
    error_band_color: Optional[str] = None,
) -> None:
    band_color = error_band_color or line_color
    for obj in objects.values():
        for key, value in (
            ("LineColor", line_color),
            ("MarkerColor", line_color),
            ("ErrorBandColor", band_color),
            ("LineStyle", line_style),
        ):
            try:
                obj.setAnnotation(key, value)
            except Exception:
                pass


def rewrite_plot_yoda_annotations(
    input_yoda: str,
    output_yoda: str,
    legend: str,
    line_color: str,
    line_style: str,
) -> str:
    objects = yoda.read(str(input_yoda))
    apply_legend_annotation(objects, legend)
    apply_line_style_annotation(objects, line_color=line_color, line_style=line_style)
    write_analysis_yoda_gz(objects, str(output_yoda))
    return str(output_yoda)


def raw_powheg_rivetplot_title() -> str:
    return "Raw_POWHEG"


RAW_POWHEG_CHANNELS = ("Compton", "BGF")


def raw_powheg_channel_rivetplot_title(channel: str) -> str:
    return f"Raw_POWHEG_{channel}"


def resolve_include_lo_requested(
    requested: Optional[bool],
    manifest: Optional[Dict[str, object]] = None,
    analysis_variant: str = "",
) -> bool:
    if requested is not None:
        return bool(requested)
    if manifest is not None:
        manifest_value = manifest.get("include_lo")
        if isinstance(manifest_value, bool):
            return manifest_value
    return use_lo_runtime(analysis_variant)


def resolve_extract_diagnostics_requested(
    requested: Optional[bool],
    manifest: Optional[Dict[str, object]] = None,
) -> bool:
    if requested is not None:
        return bool(requested)
    if manifest is not None:
        manifest_value = manifest.get("extract_diagnostics")
        if isinstance(manifest_value, bool):
            return manifest_value
    return False


def apply_legend_annotation(objects: Dict[str, object], legend: str) -> None:
    for obj in objects.values():
        try:
            obj.setAnnotation("Legend", legend)
        except Exception:
            pass


def normalize_campaign_setups(setups: Sequence[str]) -> List[str]:
    if not setups:
        return list(SETUP_ORDER)
    seen = set()
    ordered: List[str] = []
    ordering = {name: index for index, name in enumerate(SUPPORTED_SETUPS)}
    for setup in sorted(setups, key=lambda item: ordering.get(item, len(ordering))):
        if setup not in seen:
            ordered.append(setup)
            seen.add(setup)
    return ordered


def resolve_campaign_setups_requested(
    requested: Sequence[str],
    manifest: Optional[Dict[str, object]] = None,
) -> List[str]:
    if requested:
        return normalize_campaign_setups(requested)
    if manifest is not None:
        manifest_setups = manifest.get("campaign_setups")
        if isinstance(manifest_setups, list):
            setups = [str(item) for item in manifest_setups if isinstance(item, str)]
            if setups:
                return normalize_campaign_setups(setups)
    return normalize_campaign_setups(())


def resolve_analysis_name(args: argparse.Namespace, setups: Sequence[str]) -> str:
    analysis_name = str(getattr(args, "analysis_name", "MC_DIS_BREIT"))
    if contains_ps_setups(setups) and analysis_name == "MC_DIS_BREIT":
        analysis_name = PS_ANALYSIS_NAME
    args.analysis_name = analysis_name
    return analysis_name


def extend_results_extractor_scope(cmd: List[str], setups: Sequence[str]) -> None:
    cmd.append("--strict-tag")
    for setup in normalize_campaign_setups(setups):
        cmd.extend(["--setup", setup])


def resolve_scale_variations_requested(
    requested: Optional[bool],
    manifest: Optional[Dict[str, object]] = None,
    default: bool = False,
) -> bool:
    if requested is not None:
        return bool(requested)
    if manifest is not None:
        manifest_value = manifest.get("scale_variations")
        if isinstance(manifest_value, bool):
            return manifest_value
    return default


def selected_scale_variations(enabled: bool) -> Sequence[str]:
    if enabled:
        return SCALE_VARIATION_ORDER
    return ("nominal",)


def scale_variation_stem_suffix(scale_variation: str) -> str:
    try:
        return SCALE_VARIATION_STEM_SUFFIXES[scale_variation]
    except KeyError as exc:
        raise ValueError(f"Unsupported scale variation {scale_variation!r}") from exc


def build_logical_run_stem(
    setup: str,
    order: str,
    helicity: str,
    analysis_variant: str = "",
    scale_variation: str = "nominal",
) -> str:
    analysis_suffix = f"-{analysis_variant}" if analysis_variant else ""
    if order == "LO":
        base = f"DIS-POL-LO_{helicity}-{setup}"
    else:
        base = f"DIS-POL-POWHEG_{helicity}-{order}-{setup}"
    return f"{base}{analysis_suffix}{scale_variation_stem_suffix(scale_variation)}"


def generated_scale_variation_in_file(tag: str, stem: str) -> str:
    return str(Path("campaigns") / tag / "generated-inputs" / f"{stem}.in")


def build_jobs(
    lo_events: int,
    posnlo_events: int,
    negnlo_events: int,
    analysis_variant: str = "",
    analysis_variants_by_setup: Optional[Dict[str, Sequence[str]]] = None,
    setups: Sequence[str] = (),
    scale_variations: bool = False,
    campaign_tag: str = "",
    include_lo: Optional[bool] = None,
) -> List[JobSpec]:
    jobs: List[JobSpec] = []
    include_lo = resolve_include_lo_requested(include_lo, analysis_variant=analysis_variant)
    for setup in normalize_campaign_setups(setups):
        helicities = GAMMA_HELICITIES if setup == "GAMMA" else FULL_HELICITIES
        variants = (
            tuple(analysis_variants_by_setup.get(setup, (analysis_variant,)))
            if analysis_variants_by_setup is not None
            else (analysis_variant,)
        )
        if include_lo:
            for variant in variants:
                for hel in helicities:
                    nominal_stem = build_logical_run_stem(setup, "LO", hel, variant, "nominal")
                    for scale_variation in selected_scale_variations(scale_variations):
                        stem = build_logical_run_stem(setup, "LO", hel, variant, scale_variation)
                        jobs.append(
                            JobSpec(
                                setup=setup,
                                order="LO",
                                helicity=hel,
                                stem=stem,
                                run_file=f"{stem}.run",
                                in_file=(
                                    generated_scale_variation_in_file(campaign_tag, stem)
                                    if scale_variation != "nominal"
                                    else f"{stem}.in"
                                ),
                                events=lo_events,
                                analysis_variant=variant,
                                scale_variation=scale_variation,
                                scale_factor=SCALE_VARIATION_FACTORS[scale_variation],
                                source_in_file=f"{nominal_stem}.in",
                            )
                        )
        pieces = (("POSNLO", posnlo_events), ("NEGNLO", negnlo_events))
        for piece, events in pieces:
            for variant in variants:
                for hel in helicities:
                    nominal_stem = build_logical_run_stem(setup, piece, hel, variant, "nominal")
                    for scale_variation in selected_scale_variations(scale_variations):
                        stem = build_logical_run_stem(setup, piece, hel, variant, scale_variation)
                        jobs.append(
                            JobSpec(
                                setup=setup,
                                order=piece,
                                helicity=hel,
                                stem=stem,
                                run_file=f"{stem}.run",
                                in_file=(
                                    generated_scale_variation_in_file(campaign_tag, stem)
                                    if scale_variation != "nominal"
                                    else f"{stem}.in"
                                ),
                                events=events,
                                analysis_variant=variant,
                                scale_variation=scale_variation,
                                scale_factor=SCALE_VARIATION_FACTORS[scale_variation],
                                source_in_file=f"{nominal_stem}.in",
                            )
                        )
    return jobs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command")

    base_parent = argparse.ArgumentParser(add_help=False)
    base_parent.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR, help="Directory containing the DISPOL inputs and outputs.")

    tag_parent = argparse.ArgumentParser(add_help=False)
    tag_parent.add_argument("-t", "--tag", required=True, help="Herwig run tag to pass with -t.")

    analysis_name_parent = argparse.ArgumentParser(add_help=False)
    analysis_name_parent.add_argument("--analysis-name", default="MC_DIS_BREIT", help="Analysis path prefix.")

    setup_parent = argparse.ArgumentParser(add_help=False)
    setup_parent.add_argument(
        "--setup",
        action="append",
        choices=SUPPORTED_SETUPS,
        default=[],
        help="Setup(s) to operate on. May be given multiple times. In campaign/full, omitting this runs GAMMA, Z, and ALL; CC is available when requested explicitly. In analysis/postprocess/plotting, omitting this defaults to ALL.",
    )

    raw_powheg_parent = argparse.ArgumentParser(add_help=False)
    raw_powheg_parent.add_argument(
        "--raw-powheg",
        action="store_true",
        help="For RIVETFO, dump raw POWHEG momenta and build raw-POWHEG YODAs for side-by-side plotting.",
    )

    scale_variations_parent = argparse.ArgumentParser(add_help=False)
    scale_variations_parent.add_argument(
        "--scale-variations",
        action="store_true",
        default=None,
        help="Fan each DISPOL run out into nominal plus ScaleFactor up/down variations and use them for Rivet uncertainty envelopes.",
    )
    include_lo_parent = argparse.ArgumentParser(add_help=False)
    include_lo_parent.add_argument(
        "--include-lo",
        action="store_true",
        default=None,
        help="For --rivetfo workflows, also include LO runs instead of the default POWHEG-only set.",
    )

    diagnostics_parent = argparse.ArgumentParser(add_help=False)
    diagnostics_parent.add_argument(
        "--diagnostics",
        "--nlo-term-diagnostics",
        dest="diagnostics",
        action="store_true",
        default=None,
        help="Run extract_nlo_term_diagnostics.py after the main results extraction. Disabled by default.",
    )

    campaign_parent = argparse.ArgumentParser(add_help=False)
    campaign_parent.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="Maximum number of concurrent Herwig processes.")
    campaign_parent.add_argument("--shards", type=int, default=0, help="Number of shards per logical run (0 = choose automatically from --jobs).")
    campaign_parent.add_argument("--seed-base", type=int, default=100000, help="Base RNG seed used to generate distinct per-shard seeds.")
    campaign_parent.add_argument("--lo-events", type=int, default=1_000_000_000, help="Events per LO run. For --rivetfo, this is only used together with --include-lo.")
    campaign_parent.add_argument("--posnlo-events", type=int, default=1_000_000_000, help="Events per POSNLO run.")
    campaign_parent.add_argument("--negnlo-events", type=int, default=10_000_000, help="Events per NEGNLO run.")
    campaign_parent.add_argument("--no-prepare", action="store_true", help="Do not run 'Herwig read' for missing .run files.")
    campaign_parent.add_argument("--force-prepare", action="store_true", help="Always rerun 'Herwig read' from the .in cards, even if the .run files already exist.")
    campaign_parent.add_argument("--dry-run", action="store_true", help="Print the planned commands without executing them.")
    campaign_parent.add_argument("--keep-going", action="store_true", help="Continue scheduling jobs after a failure.")
    campaign_parent.add_argument(
        "--rerun-failed-random-seed",
        action="store_true",
        help="Reload the existing campaign manifest for this tag and rerun only unresolved failed shards with fresh random seeds.",
    )
    campaign_parent.add_argument(
        "--rerun-failed-shards",
        type=int,
        default=1,
        help="When rerunning failed shards with fresh random seeds, split each failed shard into this many replacement shards.",
    )
    campaign_parent.add_argument("--progress-interval", type=float, default=5.0, help="Seconds between progress refreshes.")
    campaign_parent.add_argument("--max-listed", type=int, default=12, help="Maximum number of logical runs to list in the live progress display.")
    campaign_parent.add_argument("--no-merge-yoda", action="store_true", help="Do not merge shard YODA files after the campaign.")
    campaign_parent.add_argument("--yoda-merge-tool", default="auto", help="Merge tool to use: auto, rivet-merge, or yodamerge.")
    campaign_parent.add_argument("--rivet", action="store_true", help="Use the -RIVET input/run cards and extract Rivet-mode outputs.")
    campaign_parent.add_argument("--rivetfo", action="store_true", help="Use the -RIVETFO input/run cards and extract Rivet-mode outputs.")
    postprocess_parent = argparse.ArgumentParser(add_help=False)
    postprocess_parent.add_argument("--dry-run", action="store_true", help="Print the planned commands without executing them.")
    postprocess_parent.add_argument("--no-merge-yoda", action="store_true", help="Do not rebuild merged YODA files.")
    postprocess_parent.add_argument("--yoda-merge-tool", default="auto", help="Merge tool to use: auto, rivet-merge, or yodamerge.")
    postprocess_parent.add_argument("--rivet", action="store_true", help="Use the -RIVET file naming convention.")
    postprocess_parent.add_argument("--rivetfo", action="store_true", help="Use the -RIVETFO file naming convention.")

    analyze_parent = argparse.ArgumentParser(add_help=False)
    postrun_jobs_parent = argparse.ArgumentParser(add_help=False)
    postrun_jobs_parent.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="Maximum number of concurrent post-run worker tasks.")
    analyze_io_parent = argparse.ArgumentParser(add_help=False)
    analyze_io_parent.add_argument("--rivet", action="store_true", help="Use the -RIVET file naming convention and raw-shard fallback.")
    analyze_io_parent.add_argument("--rivetfo", action="store_true", help="Use the -RIVETFO file naming convention and raw-shard fallback.")
    analyze_io_parent.add_argument("--yoda-merge-tool", default="auto", help="Merge tool to use for shard fallback: auto, rivet-merge, or yodamerge.")

    poldis_parent = argparse.ArgumentParser(add_help=False)
    poldis_parent.add_argument("--unpol", type=Path, help="POLDIS unpolarized .top file.")
    poldis_parent.add_argument("--pol", type=Path, help="POLDIS polarized .top file.")
    poldis_parent.add_argument("--order", default="NLO", choices=("NNLO", "NLO", "LO"), help="Which POLDIS dataset to export.")
    poldis_parent.add_argument("--out", type=Path, help="Output YODA path.")

    plot_parent = argparse.ArgumentParser(add_help=False)
    plot_parent.add_argument("--rivet-mkhtml-tool", default="rivet-mkhtml", help="Path or name of rivet-mkhtml.")
    plot_parent.add_argument("--ratiolabel", default="MC/POLDIS", help="Ratio label for rivet-mkhtml.")
    plot_parent.add_argument("--reflabel-prefix", default="POLDIS", help="Reference label prefix.")
    plot_parent.add_argument("--plot-dir", type=Path, help="Optional output plot directory (single-setup only).")

    subparsers.add_parser("campaign", parents=[base_parent, tag_parent, setup_parent, raw_powheg_parent, scale_variations_parent, include_lo_parent, diagnostics_parent, campaign_parent], help="Run the full validation campaign.")
    subparsers.add_parser("postprocess", parents=[base_parent, tag_parent, analysis_name_parent, setup_parent, raw_powheg_parent, scale_variations_parent, include_lo_parent, diagnostics_parent, postrun_jobs_parent, analyze_parent, postprocess_parent], help="Rebuild merged YODAs and summaries for an existing campaign tag without rerunning Herwig.")
    subparsers.add_parser("analyze-herwig", parents=[base_parent, tag_parent, analysis_name_parent, setup_parent, raw_powheg_parent, scale_variations_parent, include_lo_parent, postrun_jobs_parent, analyze_parent, analyze_io_parent], help="Build Herwig DIS polarized YODAs from campaign NLO outputs.")
    subparsers.add_parser("poldis-top", parents=[base_parent, analysis_name_parent, setup_parent, poldis_parent, analyze_io_parent], help="Convert POLDIS .top files into reference YODA.")
    subparsers.add_parser("rivetplot", parents=[base_parent, tag_parent, analysis_name_parent, setup_parent, raw_powheg_parent, scale_variations_parent, include_lo_parent, postrun_jobs_parent, analyze_parent, analyze_io_parent, poldis_parent, plot_parent], help="Run rivet-mkhtml comparing analyzed Herwig YODAs to the POLDIS reference.")
    full_parser = subparsers.add_parser("full", parents=[base_parent, tag_parent, analysis_name_parent, setup_parent, raw_powheg_parent, scale_variations_parent, include_lo_parent, diagnostics_parent, campaign_parent, analyze_parent, poldis_parent, plot_parent], help="Run campaign, analyze Herwig outputs, build POLDIS references, and make Rivet comparison plots.")
    full_parser.set_defaults(command="full")
    return parser


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    args = list(argv if argv is not None else sys.argv[1:])
    commands = {"campaign", "postprocess", "analyze-herwig", "poldis-top", "rivetplot", "full"}
    if not args:
        return build_parser().parse_args(["--help"])
    if args[0] in {"-h", "--help"}:
        return build_parser().parse_args(args)
    if args[0] not in commands:
        args = ["full", *args]
    return build_parser().parse_args(args)


def rivet_runtime_env(base_dir: Path) -> Dict[str, str]:
    env = os.environ.copy()
    analysis_dir = analyses_dir(base_dir).resolve()
    base_dir_resolved = base_dir.resolve()

    analysis_paths: List[str] = []
    if analysis_dir.exists():
        analysis_paths.append(str(analysis_dir))
    if base_dir_resolved.exists():
        analysis_paths.append(str(base_dir_resolved))
    existing_analysis = env.get("RIVET_ANALYSIS_PATH")
    if existing_analysis:
        analysis_paths.append(existing_analysis)
    if analysis_paths:
        env["RIVET_ANALYSIS_PATH"] = os.pathsep.join(analysis_paths)

    data_paths: List[str] = []
    if analysis_dir.exists():
        data_paths.append(str(analysis_dir))
    if base_dir_resolved.exists():
        data_paths.append(str(base_dir_resolved))
    existing_data = env.get("RIVET_DATA_PATH")
    if existing_data:
        data_paths.append(existing_data)
    if data_paths:
        env["RIVET_DATA_PATH"] = os.pathsep.join(data_paths)

    return env


def patch_raw_powheg_card(in_path: Path, analysis_variant: str, enable_raw_powheg: bool) -> bool:
    if analysis_variant != "RIVETFO" or not in_path.exists():
        return False

    text = in_path.read_text()
    anchor = "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission Yes"
    if anchor not in text:
        return False

    settings = [
        (
            "/Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics",
            "On" if enable_raw_powheg else "Off",
        ),
        (
            "/Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOTermDiagnostics",
            "No",
        ),
        (
            "/Herwig/MatrixElements/PowhegMEDISNCPol:DumpPOWHEGRawMomenta",
            "Yes" if enable_raw_powheg else "No",
        ),
    ]

    modified = False
    for setting, value in settings:
        replacement = f"set {setting} {value}"
        pattern = rf"^set {re.escape(setting)} \S+\s*$"
        if re.search(pattern, text, flags=re.MULTILINE):
            new_text = re.sub(pattern, replacement, text, flags=re.MULTILINE)
        else:
            new_text = text.replace(anchor, anchor + "\n" + replacement, 1)
        if new_text != text:
            modified = True
            text = new_text

    if not modified:
        return False
    in_path.write_text(text)
    return True


def set_or_insert_card_setting(text: str, setting: str, value: str, anchor: str) -> tuple[str, bool]:
    replacement = f"set {setting} {value}"
    pattern = rf"^set {re.escape(setting)} \S+\s*$"
    if re.search(pattern, text, flags=re.MULTILINE):
        new_text = re.sub(pattern, replacement, text, flags=re.MULTILINE)
    else:
        new_text = text.replace(anchor, anchor + "\n" + replacement, 1)
    return new_text, new_text != text


def patch_ps_spin_card(in_path: Path, spec: JobSpec) -> bool:
    if not spec.analysis_variant.startswith("RIVETPS-") or not in_path.exists():
        return False

    text = in_path.read_text()
    anchor = "set PowhegMEDISNCPol:Contribution PositiveNLO"
    if anchor not in text:
        anchor = "set PowhegMEDISNCPol:Contribution NegativeNLO"
    if anchor not in text:
        return False

    use_real_spin = "Yes" if spec.analysis_variant == "RIVETPS-SPIN" else "No"
    shower_spin = "No" if spec.analysis_variant == "RIVETPS-NOSPIN-UNPOL" else "Yes"
    diagnostics_enabled = spec.setup == "SPINVAL" and spec.analysis_variant == "RIVETPS-SPIN"

    modified = False
    settings = [
        ("/Herwig/MatrixElements/PowhegMEDISNCPol:UsePOWHEGRealSpinVertex", use_real_spin),
        ("/Herwig/Shower/ShowerHandler:SpinCorrelations", shower_spin),
        ("/Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics", "On" if diagnostics_enabled else "Off"),
        (
            "/Herwig/MatrixElements/PowhegMEDISNCPol:DiagnosePOWHEGRealSpinVertex",
            "Yes" if diagnostics_enabled else "No",
        ),
        (
            "/Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGRealSpinDiagMax",
            "200" if diagnostics_enabled else "0",
        ),
    ]
    for setting, value in settings:
        text, changed = set_or_insert_card_setting(text, setting, value, anchor)
        modified = modified or changed

    if not modified:
        return False
    in_path.write_text(text)
    return True


def format_scale_factor(scale_factor: float) -> str:
    return f"{scale_factor:.1f}"


def scale_variation_matrix_element(spec: JobSpec) -> str:
    return "MEDISNC" if spec.order == "LO" else "PowhegMEDISNCPol"


def resolve_scale_variation_snippet(base_dir: Path, snippet_name: str) -> Optional[Path]:
    direct_path = cards_subdir(base_dir) / "snippets" / snippet_name
    if direct_path.exists():
        return direct_path.resolve()

    preferred_candidates = (
        REPO_ROOT / "src" / "herwig" / "src" / "snippets" / snippet_name,
        REPO_ROOT / "src" / "herwig" / "share" / "Herwig" / "snippets" / snippet_name,
    )
    for candidate in preferred_candidates:
        if candidate.exists():
            return candidate.resolve()
    return None


def rewrite_scale_variation_card(source_text: str, base_dir: Path, spec: JobSpec) -> str:
    if spec.scale_variation == "nominal":
        return source_text

    matrix_element = scale_variation_matrix_element(spec)
    scale_factor_text = format_scale_factor(spec.scale_factor)
    lines = source_text.splitlines()
    insert_index: Optional[int] = None
    scale_factor_replaced = False

    scale_pattern = re.compile(
        rf"^(?P<indent>\s*set\s+(?:/Herwig/MatrixElements/)?{re.escape(matrix_element)}:ScaleFactor\s+)"
        rf"(?P<value>\S+)(?P<suffix>\s*(?:#.*)?)$"
    )
    matrix_setting_pattern = re.compile(
        rf"^set\s+(?:/Herwig/MatrixElements/)?{re.escape(matrix_element)}:"
    )
    matrix_insert_pattern = re.compile(
        rf"^insert\s+SubProcess:MatrixElements\[0\]\s+{re.escape(matrix_element)}\s*$"
    )
    saverun_pattern = re.compile(r"^(?P<indent>\s*saverun\s+)\S+(?P<suffix>\s+EventGenerator(?:\s*#.*)?)$")
    snippet_pattern = re.compile(r"^(?P<indent>\s*read\s+)snippets/(?P<snippet>\S+?)(?P<suffix>\s*(?:#.*)?)$")

    for index, line in enumerate(lines):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        match = scale_pattern.match(line)
        if match:
            lines[index] = f"{match.group('indent')}{scale_factor_text}{match.group('suffix')}"
            scale_factor_replaced = True
            insert_index = index + 1
            continue
        if matrix_setting_pattern.match(stripped):
            insert_index = index + 1
            continue
        if insert_index is None and matrix_insert_pattern.match(stripped):
            insert_index = index + 1
            continue
        match = saverun_pattern.match(line)
        if match:
            lines[index] = f"{match.group('indent')}{spec.stem}{match.group('suffix')}"
            continue
        match = snippet_pattern.match(line)
        if match:
            snippet_path = resolve_scale_variation_snippet(base_dir, match.group("snippet"))
            if snippet_path is not None:
                lines[index] = f"{match.group('indent')}{snippet_path.as_posix()}{match.group('suffix')}"

    if not scale_factor_replaced:
        if insert_index is None:
            raise RuntimeError(
                f"Could not find where to set {matrix_element}:ScaleFactor in {spec.source_in_file or spec.in_file}"
            )
        lines.insert(insert_index, f"set {matrix_element}:ScaleFactor {scale_factor_text}")

    text = "\n".join(lines)
    if source_text.endswith("\n"):
        text += "\n"
    return text


def materialize_scale_variation_card(base_dir: Path, spec: JobSpec, dry_run: bool = False) -> None:
    if spec.scale_variation == "nominal":
        return
    if not spec.source_in_file:
        raise RuntimeError(f"Missing source input card for scale variation job {spec.stem}")

    source_path = card_path(base_dir, spec.source_in_file)
    target_path = card_path(base_dir, spec.in_file)
    if not source_path.exists():
        raise FileNotFoundError(f"Missing source input file {source_path}")

    scale_factor_text = format_scale_factor(spec.scale_factor)
    rendered = rewrite_scale_variation_card(source_path.read_text(), base_dir, spec)
    if dry_run:
        print(f"[dry-run] generate-scale-card {source_path} -> {target_path} (ScaleFactor={scale_factor_text})")
        return

    target_path.parent.mkdir(parents=True, exist_ok=True)
    if target_path.exists() and target_path.read_text() == rendered:
        return
    target_path.write_text(rendered)


def ensure_run_file(
    base_dir: Path,
    spec: JobSpec,
    do_prepare: bool,
    dry_run: bool = False,
    analysis_variant: str = "",
    force_prepare: bool = False,
    enable_raw_powheg: bool = False,
) -> Optional[JobResult]:
    run_path = base_dir / spec.run_file
    in_path = card_path(base_dir, spec.in_file)
    materialize_scale_variation_card(base_dir, spec, dry_run=dry_run)
    card_patched = False
    effective_analysis_variant = spec.analysis_variant or analysis_variant
    if patch_raw_powheg_card(in_path, effective_analysis_variant, enable_raw_powheg):
        card_patched = True
    if patch_ps_spin_card(in_path, spec):
        card_patched = True
    if run_path.exists() and not force_prepare:
        if not card_patched:
            return None
        force_prepare = True
    if not do_prepare:
        raise FileNotFoundError(f"Missing run file {run_path}")
    if not in_path.exists() and not dry_run:
        raise FileNotFoundError(f"Missing input file {in_path}")
    prepare_cmd = ["Herwig", "read", str(card_relative_path(base_dir, spec.in_file))]
    dummy_spec = ShardSpec(job=spec, shard_index=1, shard_count=1, tag="", seed=0, events=spec.events)
    if dry_run:
        print(shlex.join(prepare_cmd))
        result = JobResult(spec=dummy_spec, command=[], prepared=True, prepare_command=prepare_cmd, prepare_returncode=0)
        return result
    env = rivet_runtime_env(base_dir) if use_rivet_runtime(effective_analysis_variant) else None
    proc = subprocess.run(prepare_cmd, cwd=base_dir, env=env, text=True)
    result = JobResult(spec=dummy_spec, command=[], prepared=True, prepare_command=prepare_cmd, prepare_returncode=proc.returncode)
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to prepare {spec.run_file} from {spec.in_file}")
    return result


def launcher_log_path(campaign_dir: Path, spec: ShardSpec) -> Path:
    safe = f"{spec.job.stem}-{spec.tag}.launcher.log"
    return campaign_dir / "launcher-logs" / safe


def collect_artifacts(base_dir: Path, stem: str, tag: str, seed: int) -> Dict[str, List[str]]:
    prefix = f"{stem}-S{seed}-{tag}"
    files = sorted(base_dir.glob(f"{prefix}*"))
    out_files: List[str] = []
    log_files: List[str] = []
    yoda_files: List[str] = []
    for path in files:
        name = path.name
        if name.endswith(".out"):
            out_files.append(str(path))
        elif name.endswith(".log"):
            log_files.append(str(path))
        elif ".yoda" in name:
            yoda_files.append(str(path))
    return {"out_files": out_files, "log_files": log_files, "yoda_files": yoda_files}


def choose_shards(logical_jobs: int, max_jobs: int, requested_shards: int) -> int:
    if requested_shards > 0:
        return requested_shards
    if logical_jobs <= 0:
        return 1
    return max(1, (max_jobs + logical_jobs - 1) // logical_jobs)


def split_events(total_events: int, shards: int) -> List[int]:
    if shards <= 1:
        return [total_events]
    base = total_events // shards
    remainder = total_events % shards
    parts = [base + (1 if index < remainder else 0) for index in range(shards)]
    return [part for part in parts if part > 0]


def build_shards(jobs: Sequence[JobSpec], base_tag: str, requested_shards: int, max_jobs: int, seed_base: int) -> List[ShardSpec]:
    shard_target = choose_shards(len(jobs), max_jobs, requested_shards)
    shards: List[ShardSpec] = []
    seed = seed_base
    for job in jobs:
        event_parts = split_events(job.events, shard_target)
        shard_count = len(event_parts)
        for index, events in enumerate(event_parts, start=1):
            tag = base_tag if shard_count == 1 else f"{base_tag}-s{index:03d}"
            shards.append(
                ShardSpec(
                    job=job,
                    shard_index=index,
                    shard_count=shard_count,
                    tag=tag,
                    seed=seed,
                    events=events,
                )
            )
            seed += 1
    return shards


def logical_run_label(spec: ShardSpec) -> str:
    return spec.job.stem


def logical_group_counts(shards: Sequence[ShardSpec]) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for shard in shards:
        counts[logical_run_label(shard)] = counts.get(logical_run_label(shard), 0) + 1
    return counts


def latest_progress_marker(result: JobResult) -> str:
    if not result.launcher_log:
        return "-"
    path = Path(result.launcher_log)
    if not path.exists():
        return "-"
    try:
        with path.open("rb") as handle:
            handle.seek(0, os.SEEK_END)
            size = handle.tell()
            handle.seek(max(0, size - 65536))
            chunk = handle.read().decode("utf-8", errors="replace")
    except OSError:
        return "-"

    markers = list(PROGRESS_MARKER_RE.finditer(chunk.replace("\r", "\n")))
    if not markers:
        return "-"
    last = markers[-1]
    current = last.group("current")
    total = last.group("total") or last.group("total_alt")
    if current == "init":
        return "init"
    if not total:
        return current
    try:
        cur_val = int(current)
        total_val = int(total)
    except ValueError:
        return f"{current}/{total}"
    if total_val <= 0:
        return f"{cur_val}/{total_val}"
    percent = 100.0 * cur_val / total_val
    return f"{cur_val}/{total_val} ({percent:.1f}%)"


def progress_lines(
    all_shards: Sequence[ShardSpec],
    pending: Sequence[ShardSpec],
    active_results: Sequence[JobResult],
    finished: Sequence[JobResult],
    max_listed: int,
    started_at: float,
) -> List[str]:
    total = len(all_shards)
    done = sum(1 for item in finished if item.returncode == 0)
    running = len(active_results)
    waiting = len(pending)
    failed = sum(1 for item in finished if item.returncode not in (0, None))
    elapsed = time.time() - started_at

    pending_counts = logical_group_counts(pending)
    running_specs = [item.spec for item in active_results]
    running_counts = logical_group_counts(running_specs)
    all_counts = logical_group_counts(all_shards)

    logical_rows = []
    logical_names = sorted(set(pending_counts) | set(running_counts))
    for name in logical_names[:max_listed]:
        logical_rows.append(
            [
                name,
                str(running_counts.get(name, 0)),
                str(pending_counts.get(name, 0)),
                str(all_counts.get(name, 0)),
            ]
        )

    running_rows = []
    for item in active_results[:max_listed]:
        shard = item.spec
        duration = 0.0 if item.started_at is None else time.time() - item.started_at
        running_rows.append(
            [
                shard.job.stem,
                shard.tag,
                latest_progress_marker(item),
                str(shard.events),
                str(shard.seed),
                f"{duration:.0f}s",
            ]
        )

    lines = [
        f"Tag: {all_shards[0].tag.rsplit('-s', 1)[0] if all_shards else 'n/a'}",
        f"Elapsed: {elapsed:.0f}s",
        f"Shards: completed {done}/{total} | running {running} | pending {waiting} | failed {failed}",
    ]
    if logical_rows:
        lines.append("")
        lines.append("Logical runs with work remaining:")
        lines.extend(render_table(["Run", "Running", "Pending", "Total"], logical_rows, aligns=("l", "r", "r", "r")))
    if running_rows:
        lines.append("")
        lines.append("Active shards:")
        lines.extend(render_table(["Run", "Tag", "Progress", "Events", "Seed", "Runtime"], running_rows, aligns=("l", "l", "r", "r", "r", "r")))
    return lines


def emit_progress(lines: Sequence[str], interactive: bool) -> None:
    text = "\n".join(lines)
    if interactive:
        sys.stdout.write("\x1b[2J\x1b[H")
        sys.stdout.write(text + "\n")
        sys.stdout.flush()
    else:
        sys.stdout.write(text + "\n\n")
        sys.stdout.flush()


def run_jobs(
    base_dir: Path,
    jobs: Sequence[ShardSpec],
    max_jobs: int,
    campaign_dir: Path,
    dry_run: bool,
    keep_going: bool,
    progress_interval: float,
    max_listed: int,
    analysis_variant: str,
) -> List[JobResult]:
    pending = list(jobs)
    active: List[tuple[subprocess.Popen[str], JobResult, object]] = []
    finished: List[JobResult] = []
    campaign_started = time.time()
    interactive = sys.stdout.isatty() and not dry_run
    last_progress = 0.0

    def start_job(spec: ShardSpec) -> None:
        cmd = ["Herwig", "run", spec.job.run_file, "-N", str(spec.events), "-t", spec.tag, "-s", str(spec.seed)]
        result = JobResult(spec=spec, command=cmd)
        log_path = launcher_log_path(campaign_dir, spec)
        result.launcher_log = str(log_path)
        if dry_run:
            print(shlex.join(cmd))
            result.returncode = 0
            result.started_at = result.ended_at = time.time()
            result.__dict__.update(collect_artifacts(base_dir, spec.job.stem, spec.tag, spec.seed))
            finished.append(result)
            return
        log_path.parent.mkdir(parents=True, exist_ok=True)
        handle = open(log_path, "w")
        result.started_at = time.time()
        effective_analysis_variant = spec.job.analysis_variant or analysis_variant
        env = rivet_runtime_env(base_dir) if use_rivet_runtime(effective_analysis_variant) else None
        proc = subprocess.Popen(cmd, cwd=base_dir, env=env, stdout=handle, stderr=subprocess.STDOUT, text=True)
        active.append((proc, result, handle))

    while pending or active:
        while pending and len(active) < max_jobs:
            start_job(pending.pop(0))
        if progress_interval >= 0.0:
            now = time.time()
            if now - last_progress >= progress_interval or (not pending and active):
                emit_progress(
                    progress_lines(jobs, pending, [result for _, result, _ in active], finished, max_listed, campaign_started),
                    interactive,
                )
                last_progress = now
        if not active:
            continue
        time.sleep(1.0)
        still_active: List[tuple[subprocess.Popen[str], JobResult, object]] = []
        for proc, result, handle in active:
            code = proc.poll()
            if code is None:
                still_active.append((proc, result, handle))
                continue
            handle.close()
            result.returncode = code
            result.ended_at = time.time()
            result.__dict__.update(collect_artifacts(base_dir, result.spec.job.stem, result.spec.tag, result.spec.seed))
            finished.append(result)
            if code != 0 and not keep_going:
                for proc2, result2, handle2 in still_active:
                    proc2.terminate()
                    handle2.close()
                    result2.returncode = proc2.wait()
                    result2.ended_at = time.time()
                    result2.__dict__.update(collect_artifacts(base_dir, result2.spec.job.stem, result2.spec.tag, result2.spec.seed))
                    finished.append(result2)
                emit_progress(
                    progress_lines(jobs, pending, [], finished, max_listed, campaign_started),
                    interactive,
                )
                return finished
        active = still_active
    emit_progress(progress_lines(jobs, pending, [], finished, max_listed, campaign_started), interactive)
    return finished


def run_extractor(cmd: List[str], output_path: Path, base_dir: Path, dry_run: bool) -> int:
    if dry_run:
        print(shlex.join(cmd))
        return 0
    proc = subprocess.run(cmd, cwd=base_dir, text=True, capture_output=True)
    output_path.write_text(proc.stdout)
    if proc.stderr:
        output_path.with_suffix(output_path.suffix + ".stderr").write_text(proc.stderr)
    return proc.returncode


def choose_yoda_merge_tool(preference: str, purpose: str = "generic") -> Optional[str]:
    if preference != "auto":
        return shutil.which(preference)
    if purpose == "shard":
        candidates = ("yodamerge", "rivet-merge")
    elif purpose == "nlo":
        candidates = ("yodamerge", "rivet-merge")
    else:
        candidates = ("yodamerge", "rivet-merge")
    for candidate in candidates:
        resolved = shutil.which(candidate)
        if resolved is not None:
            return resolved
    return None


def rivet_merge_env(base_dir: Path) -> Dict[str, str]:
    env = os.environ.copy()
    analysis_paths: List[str] = []
    analysis_dir = analyses_dir(base_dir).resolve()
    if analysis_dir.exists():
        analysis_paths.append(str(analysis_dir))
    base_dir_resolved = base_dir.resolve()
    if base_dir_resolved.exists():
        analysis_paths.append(str(base_dir_resolved))
    existing = env.get("RIVET_ANALYSIS_PATH")
    if existing:
        analysis_paths.append(existing)
    if analysis_paths:
        env["RIVET_ANALYSIS_PATH"] = os.pathsep.join(analysis_paths)
    return env


def merge_output_path(campaign_dir: Path, base_tag: str, logical_run: str, source_file: str) -> Path:
    suffixes = Path(source_file).suffixes
    ext = "".join(suffixes) if suffixes else ".yoda"
    return campaign_dir / "yoda" / f"{logical_run}-{base_tag}{ext}"


def merge_yoda_files(
    finished: Sequence[JobResult],
    base_dir: Path,
    campaign_dir: Path,
    base_tag: str,
    merge_tool_preference: str,
    max_workers: int,
    dry_run: bool,
) -> List[YODAMergeResult]:
    by_logical_run: Dict[str, List[str]] = {}
    for item in finished:
        if item.returncode != 0:
            continue
        if not item.yoda_files:
            continue
        by_logical_run.setdefault(item.spec.job.stem, []).extend(item.yoda_files)

    merge_dir = campaign_dir / "yoda"
    merge_dir.mkdir(parents=True, exist_ok=True)
    merge_tool = choose_yoda_merge_tool(merge_tool_preference, purpose="shard")
    if by_logical_run:
        worker_count = max(1, min(max_workers, len(by_logical_run)))
        mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
        print_stage(f"Merging shard YODA files for {len(by_logical_run)} logical runs ({mode})")

    tasks = sorted(by_logical_run.items())

    def merge_one(task: tuple[str, List[str]]) -> YODAMergeResult:
        logical_run, files = task
        unique_files = sorted({str(Path(path).resolve()) for path in files})
        if not unique_files:
            return YODAMergeResult(
                logical_run=logical_run,
                source_files=[],
                output_file=None,
                command=[],
                returncode=0,
                skipped=True,
                message="No YODA files found",
            )

        output_path = merge_output_path(campaign_dir, base_tag, logical_run, unique_files[0])
        if len(unique_files) == 1:
            command = ["cp", unique_files[0], str(output_path)]
            if dry_run:
                print(shlex.join(command))
                rc = 0
            else:
                shutil.copy2(unique_files[0], output_path)
                rc = 0
            return YODAMergeResult(
                logical_run=logical_run,
                source_files=unique_files,
                output_file=str(output_path),
                command=command,
                returncode=rc,
                copied=True,
                message="Single YODA file copied",
            )

        if merge_tool is None:
            return YODAMergeResult(
                logical_run=logical_run,
                source_files=unique_files,
                output_file=None,
                command=[],
                returncode=127,
                skipped=True,
                message="No YODA merge tool found (looked for rivet-merge, yodamerge)",
            )

        tool_name = Path(merge_tool).name
        command = [merge_tool]
        if tool_name == "rivet-merge":
            command.append("-e")
        command.extend([*unique_files, "-o", str(output_path)])
        if dry_run:
            print(shlex.join(command))
            rc = 0
        else:
            env = rivet_merge_env(base_dir) if tool_name == "rivet-merge" else None
            proc = subprocess.run(command, cwd=base_dir, env=env, text=True, capture_output=True)
            rc = proc.returncode
            stem = f"{logical_run}-{base_tag}"
            if proc.stdout:
                (merge_dir / f"{stem}.merge.stdout").write_text(proc.stdout)
            if proc.stderr:
                (merge_dir / f"{stem}.merge.stderr").write_text(proc.stderr)
        return YODAMergeResult(
            logical_run=logical_run,
            source_files=unique_files,
            output_file=str(output_path) if rc == 0 else None,
            command=command,
            returncode=rc,
            message=None if rc == 0 else f"Merge failed with rc={rc}",
        )

    return run_parallel_ordered(tasks, merge_one, 1 if dry_run else max_workers)


def build_nlo_yoda_files(
    yoda_merge_results: Sequence[YODAMergeResult],
    base_dir: Path,
    campaign_dir: Path,
    base_tag: str,
    merge_tool_preference: str,
    max_workers: int,
    dry_run: bool,
) -> List[YODANLOResult]:
    merged_by_run: Dict[str, str] = {}
    for item in yoda_merge_results:
        if item.returncode == 0 and item.output_file is not None:
            merged_by_run[item.logical_run] = item.output_file

    merge_tool = choose_yoda_merge_tool(merge_tool_preference, purpose="nlo")
    merge_dir = campaign_dir / "yoda"
    merge_dir.mkdir(parents=True, exist_ok=True)

    pos_runs = sorted(name for name in merged_by_run if "-POSNLO-" in name)
    if pos_runs:
        worker_count = max(1, min(max_workers, len(pos_runs)))
        mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
        print_stage(f"Building physical NLO YODAs from {len(pos_runs)} POSNLO/NEGNLO pairs ({mode})")

    def build_one(pos_run: str) -> YODANLOResult:
        neg_run = pos_run.replace("-POSNLO-", "-NEGNLO-", 1)
        nlo_run = pos_run.replace("-POSNLO-", "-NLO-", 1)
        pos_file = merged_by_run.get(pos_run)
        neg_file = merged_by_run.get(neg_run)
        if pos_file is None or neg_file is None:
            return YODANLOResult(
                logical_run=nlo_run,
                pos_file=pos_file,
                neg_file=neg_file,
                output_file=None,
                command=[],
                returncode=1,
                skipped=True,
                message=f"Missing {'POSNLO' if pos_file is None else 'NEGNLO'} merged YODA",
            )
        if merge_tool is None:
            return YODANLOResult(
                logical_run=nlo_run,
                pos_file=pos_file,
                neg_file=neg_file,
                output_file=None,
                command=[],
                returncode=127,
                skipped=True,
                message="No YODA merge tool found (looked for rivet-merge, yodamerge)",
            )

        source = pos_file
        ext = "".join(Path(source).suffixes) if Path(source).suffixes else ".yoda"
        output_path = merge_dir / f"{nlo_run}-{base_tag}{ext}"
        tool_name = Path(merge_tool).name
        if tool_name == "yodamerge":
            command = [merge_tool, f"{pos_file}:1", f"{neg_file}:-1", "-o", str(output_path)]
        else:
            command = [merge_tool, f"{pos_file}:x1", f"{neg_file}:x-1", "-o", str(output_path)]
        if dry_run:
            print(shlex.join(command))
            rc = 0
        else:
            env = rivet_merge_env(base_dir) if tool_name == "rivet-merge" else None
            proc = subprocess.run(command, cwd=base_dir, env=env, text=True, capture_output=True)
            rc = proc.returncode
            stem = f"{nlo_run}-{base_tag}"
            if proc.stdout:
                (merge_dir / f"{stem}.nlo.stdout").write_text(proc.stdout)
            if proc.stderr:
                (merge_dir / f"{stem}.nlo.stderr").write_text(proc.stderr)
        return YODANLOResult(
            logical_run=nlo_run,
            pos_file=pos_file,
            neg_file=neg_file,
            output_file=str(output_path) if rc == 0 else None,
            command=command,
            returncode=rc,
            message=None if rc == 0 else f"NLO merge failed with rc={rc}",
        )
    return run_parallel_ordered(pos_runs, build_one, 1 if dry_run else max_workers)


def skipped_nlo_yoda_results(yoda_merge_results: Sequence[YODAMergeResult], reason: str) -> List[YODANLOResult]:
    merged_by_run: Dict[str, str] = {}
    for item in yoda_merge_results:
        if item.returncode == 0 and item.output_file is not None:
            merged_by_run[item.logical_run] = item.output_file

    results: List[YODANLOResult] = []
    pos_runs = sorted(name for name in merged_by_run if "-POSNLO-" in name)
    for pos_run in pos_runs:
        neg_run = pos_run.replace("-POSNLO-", "-NEGNLO-", 1)
        nlo_run = pos_run.replace("-POSNLO-", "-NLO-", 1)
        results.append(
            YODANLOResult(
                logical_run=nlo_run,
                pos_file=merged_by_run.get(pos_run),
                neg_file=merged_by_run.get(neg_run),
                output_file=None,
                command=[],
                returncode=0,
                skipped=True,
                message=reason,
            )
        )
    return results


def write_manifest(
    campaign_dir: Path,
    args: argparse.Namespace,
    prepared: List[JobResult],
    finished: List[JobResult],
    extract_status: Dict[str, int],
    yoda_merge_results: Sequence[YODAMergeResult],
    yoda_nlo_results: Sequence[YODANLOResult],
    powheg_options_file: Optional[str] = None,
) -> None:
    resolved_analysis_variant = resolved_analysis_variant_from_args(args)
    include_lo = resolve_include_lo_requested(getattr(args, "include_lo", None), analysis_variant=resolved_analysis_variant)
    analysis_variants = sorted({item.spec.job.analysis_variant for item in finished if item.spec.job.analysis_variant})
    if not analysis_variants:
        analysis_variants = sorted({item.spec.job.analysis_variant for item in prepared if item.spec.job.analysis_variant})
    successful_finished = [item for item in finished if item.returncode == 0]
    failed_finished = unresolved_failed_results(finished)

    def manifest_job_result_entry(item: JobResult) -> Dict[str, object]:
        return {
            "spec": asdict(item.spec),
            "command": item.command,
            "returncode": item.returncode,
            "started_at": item.started_at,
            "ended_at": item.ended_at,
            "duration_s": item.duration_s,
            "launcher_log": item.launcher_log,
            "out_files": item.out_files,
            "log_files": item.log_files,
            "yoda_files": item.yoda_files,
        }

    manifest = {
        "tag": args.tag,
        "base_dir": str(args.base_dir),
        "jobs_limit": args.jobs,
        "shards_per_logical_run": args.shards,
        "seed_base": args.seed_base,
        "merge_yoda": not args.no_merge_yoda,
        "force_prepare": bool(getattr(args, "force_prepare", False)),
        "raw_powheg": bool(getattr(args, "raw_powheg", False)),
        "include_lo": include_lo,
        "extract_diagnostics": bool(getattr(args, "diagnostics", False)),
        "scale_variations": bool(getattr(args, "scale_variations", False)),
        "scale_variation_factors": dict(SCALE_VARIATION_FACTORS),
        "yoda_merge_tool": args.yoda_merge_tool,
        "rivet": bool(resolved_analysis_variant),
        "analysis_variant": resolved_analysis_variant,
        "analysis_variants": analysis_variants,
        "ps_mode": bool(analysis_variants and all(variant.startswith("RIVETPS-") for variant in analysis_variants)),
        "campaign_setups": normalize_campaign_setups(getattr(args, "setup", [])),
        "event_counts": {
            "LO": args.lo_events,
            "POSNLO": args.posnlo_events,
            "NEGNLO": args.negnlo_events,
        },
        "prepared": [
            {
                "spec": asdict(item.spec),
                "prepare_command": item.prepare_command,
                "prepare_returncode": item.prepare_returncode,
            }
            for item in prepared
        ],
        "finished": [manifest_job_result_entry(item) for item in successful_finished],
        "failed": [manifest_job_result_entry(item) for item in failed_finished],
        "finished_count": len(successful_finished),
        "failed_count": len(failed_finished),
        "extract_status": extract_status,
        "all_yoda_files": sorted({y for item in successful_finished for y in item.yoda_files}),
        "merged_yoda": [
            {
                "logical_run": item.logical_run,
                "source_files": item.source_files,
                "output_file": item.output_file,
                "command": item.command,
                "returncode": item.returncode,
                "copied": item.copied,
                "skipped": item.skipped,
                "message": item.message,
            }
            for item in yoda_merge_results
        ],
        "nlo_yoda": [
            {
                "logical_run": item.logical_run,
                "pos_file": item.pos_file,
                "neg_file": item.neg_file,
                "output_file": item.output_file,
                "command": item.command,
                "returncode": item.returncode,
                "skipped": item.skipped,
                "message": item.message,
            }
            for item in yoda_nlo_results
        ],
        "all_merged_yoda_files": sorted(
            {
                item.output_file
                for item in yoda_merge_results
                if item.output_file is not None
            }
        ),
        "all_nlo_yoda_files": sorted(
            {
                item.output_file
                for item in yoda_nlo_results
                if item.output_file is not None
            }
        ),
    }
    if powheg_options_file:
        manifest["powheg_options_file"] = powheg_options_file
    (campaign_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))


def load_manifest(campaign_dir: Path) -> Dict[str, object]:
    manifest_path = campaign_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    return json.loads(manifest_path.read_text())


def save_manifest(campaign_dir: Path, manifest: Dict[str, object]) -> None:
    (campaign_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))


def merge_nested_manifest_outputs(
    manifest: Dict[str, object],
    key: str,
    outputs: Dict[str, Dict[str, str]],
) -> None:
    existing = manifest.get(key, {})
    if not isinstance(existing, dict):
        existing = {}
    for family, family_outputs in outputs.items():
        family_existing = existing.get(family, {})
        if not isinstance(family_existing, dict):
            family_existing = {}
        family_existing.update(family_outputs)
        existing[family] = family_existing
    manifest[key] = existing


def update_manifest_outputs(
    campaign_dir: Path,
    analysis_outputs: Optional[Dict[str, str]] = None,
    analysis_plot_outputs: Optional[Dict[str, str]] = None,
    analysis_outputs_by_family: Optional[Dict[str, Dict[str, str]]] = None,
    analysis_plot_outputs_by_family: Optional[Dict[str, Dict[str, str]]] = None,
    raw_powheg_analysis_outputs: Optional[Dict[str, str]] = None,
    raw_powheg_analysis_plot_outputs: Optional[Dict[str, str]] = None,
    raw_powheg_channel_analysis_outputs: Optional[Dict[str, Dict[str, str]]] = None,
    raw_powheg_channel_analysis_plot_outputs: Optional[Dict[str, Dict[str, str]]] = None,
    reference_output: Optional[str] = None,
    plot_outputs: Optional[Dict[str, str]] = None,
    ps_hadronization_plot_outputs: Optional[Dict[str, str]] = None,
    raw_powheg_yoda_outputs: Optional[List[Dict[str, object]]] = None,
) -> None:
    manifest = load_manifest(campaign_dir)
    if analysis_outputs:
        existing = manifest.get("analysis_yoda", {})
        if not isinstance(existing, dict):
            existing = {}
        existing.update(analysis_outputs)
        manifest["analysis_yoda"] = existing
    if analysis_plot_outputs:
        existing = manifest.get("analysis_plot_yoda", {})
        if not isinstance(existing, dict):
            existing = {}
        existing.update(analysis_plot_outputs)
        manifest["analysis_plot_yoda"] = existing
    if analysis_outputs_by_family:
        merge_nested_manifest_outputs(manifest, "analysis_yoda_by_family", analysis_outputs_by_family)
    if analysis_plot_outputs_by_family:
        merge_nested_manifest_outputs(manifest, "analysis_plot_yoda_by_family", analysis_plot_outputs_by_family)
    if raw_powheg_analysis_outputs:
        existing = manifest.get("raw_powheg_analysis_yoda", {})
        if not isinstance(existing, dict):
            existing = {}
        existing.update(raw_powheg_analysis_outputs)
        manifest["raw_powheg_analysis_yoda"] = existing
    if raw_powheg_analysis_plot_outputs:
        existing = manifest.get("raw_powheg_analysis_plot_yoda", {})
        if not isinstance(existing, dict):
            existing = {}
        existing.update(raw_powheg_analysis_plot_outputs)
        manifest["raw_powheg_analysis_plot_yoda"] = existing
    if raw_powheg_channel_analysis_outputs:
        existing = manifest.get("raw_powheg_channel_analysis_yoda", {})
        if not isinstance(existing, dict):
            existing = {}
        for channel, outputs in raw_powheg_channel_analysis_outputs.items():
            channel_existing = existing.get(channel, {})
            if not isinstance(channel_existing, dict):
                channel_existing = {}
            channel_existing.update(outputs)
            existing[channel] = channel_existing
        manifest["raw_powheg_channel_analysis_yoda"] = existing
    if raw_powheg_channel_analysis_plot_outputs:
        existing = manifest.get("raw_powheg_channel_analysis_plot_yoda", {})
        if not isinstance(existing, dict):
            existing = {}
        for channel, outputs in raw_powheg_channel_analysis_plot_outputs.items():
            channel_existing = existing.get(channel, {})
            if not isinstance(channel_existing, dict):
                channel_existing = {}
            channel_existing.update(outputs)
            existing[channel] = channel_existing
        manifest["raw_powheg_channel_analysis_plot_yoda"] = existing
    if reference_output is not None:
        manifest["reference_yoda"] = reference_output
    if raw_powheg_yoda_outputs is not None:
        manifest["raw_powheg_yoda"] = raw_powheg_yoda_outputs
    if plot_outputs:
        existing_plots = manifest.get("plot_dirs", {})
        if not isinstance(existing_plots, dict):
            existing_plots = {}
        existing_plots.update(plot_outputs)
        manifest["plot_dirs"] = existing_plots
    if ps_hadronization_plot_outputs:
        existing_hadronization = manifest.get("ps_hadronization_plot_outputs", {})
        if not isinstance(existing_hadronization, dict):
            existing_hadronization = {}
        existing_hadronization.update(ps_hadronization_plot_outputs)
        manifest["ps_hadronization_plot_outputs"] = existing_hadronization
    if manifest:
        save_manifest(campaign_dir, manifest)


def resolve_campaign_dir(base_dir: Path, tag: str) -> Path:
    return base_dir / "campaigns" / tag


def raw_powheg_output_path(
    campaign_dir: Path,
    result: JobResult,
    channel: Optional[str] = None,
    variant: Optional[str] = None,
) -> Path:
    suffix = ".rawpowheg"
    if variant:
        suffix += f".{variant}"
    if channel:
        suffix += f".{channel}"
    suffix += ".yoda.gz"
    return campaign_dir / "raw-powheg" / f"{result.spec.job.stem}-S{result.spec.seed}-{result.spec.tag}{suffix}"


def raw_powheg_summary_csv_path(
    campaign_dir: Path,
    result: JobResult,
) -> Path:
    return campaign_dir / "raw-powheg" / (
        f"{result.spec.job.stem}-S{result.spec.seed}-{result.spec.tag}.rawpowheg.summary.csv"
    )


def run_raw_powheg_conversion_task(
    base_dir: Path,
    output_dir: Path,
    logical_run: str,
    channel: Optional[str],
    variant: Optional[str],
    log_file: Optional[str],
    out_file: Optional[str],
    output_path: Path,
    summary_csv: Optional[Path],
    command: List[str],
) -> RawPOWHEGYODAResult:
    proc = subprocess.run(command, cwd=base_dir, text=True, capture_output=True)
    if proc.stdout:
        (output_dir / f"{output_path.stem}.stdout").write_text(proc.stdout)
    if proc.stderr:
        (output_dir / f"{output_path.stem}.stderr").write_text(proc.stderr)
    return RawPOWHEGYODAResult(
        logical_run=logical_run,
        channel=channel,
        variant=variant,
        source_log=log_file,
        source_out=out_file,
        output_file=str(output_path),
        summary_csv=str(summary_csv) if summary_csv is not None else None,
        command=command,
        returncode=proc.returncode,
        skipped=proc.returncode != 0,
        message=None if proc.returncode == 0 else "raw POWHEG conversion failed",
    )


def build_raw_powheg_yoda_files(
    base_dir: Path,
    campaign_dir: Path,
    results: Sequence[JobResult],
    analysis_variant: str,
    max_workers: int = 1,
    dry_run: bool = False,
) -> List[RawPOWHEGYODAResult]:
    if analysis_variant != "RIVETFO":
        return []

    script = script_path("powheg_raw_momenta_to_yoda.py")
    output_dir = campaign_dir / "raw-powheg"
    output_dir.mkdir(parents=True, exist_ok=True)
    built: List[RawPOWHEGYODAResult] = []
    pending_tasks: List[tuple[str, Optional[str], Optional[str], Optional[str], Optional[str], Path, Optional[Path], List[str]]] = []

    for result in results:
        if result.spec.job.order not in {"POSNLO", "NEGNLO"}:
            continue
        raw_sources: List[str] = []
        seen_sources: set[str] = set()
        def add_source(path: Optional[str]) -> None:
            if not path:
                return
            resolved = str(Path(path))
            if resolved in seen_sources:
                return
            if not Path(resolved).exists():
                return
            seen_sources.add(resolved)
            raw_sources.append(resolved)

        for path in result.log_files:
            add_source(path)
        for path in result.out_files:
            add_source(path)
        if result.launcher_log:
            add_source(result.launcher_log)

        log_file = raw_sources[0] if raw_sources else None
        out_file = next((path for path in result.out_files if path.endswith(".out")), None)
        logical_run = result.spec.job.stem

        if not raw_sources or out_file is None:
            for channel, variant in [(None, None), *[(name, None) for name in RAW_POWHEG_CHANNELS]]:
                output_path = raw_powheg_output_path(campaign_dir, result, channel, variant)
                summary_csv = raw_powheg_summary_csv_path(campaign_dir, result) if channel is None else None
                built.append(
                    RawPOWHEGYODAResult(
                        logical_run=logical_run,
                        channel=channel,
                        variant=variant,
                        source_log=log_file,
                        source_out=out_file,
                        output_file=str(output_path),
                        summary_csv=str(summary_csv) if summary_csv is not None else None,
                        command=[],
                        returncode=1,
                        skipped=True,
                        message="Missing text artifacts or .out file for raw POWHEG conversion.",
                    )
                )
            continue

        measurement, _ = parse_measurement(Path(out_file).read_text(errors="replace"))
        if measurement is None:
            for channel, variant in [(None, None), *[(name, None) for name in RAW_POWHEG_CHANNELS]]:
                output_path = raw_powheg_output_path(campaign_dir, result, channel, variant)
                summary_csv = raw_powheg_summary_csv_path(campaign_dir, result) if channel is None else None
                built.append(
                    RawPOWHEGYODAResult(
                        logical_run=logical_run,
                        channel=channel,
                        variant=variant,
                        source_log=log_file,
                        source_out=out_file,
                        output_file=str(output_path),
                        summary_csv=str(summary_csv) if summary_csv is not None else None,
                        command=[],
                        returncode=1,
                        skipped=True,
                        message="Could not extract total cross section from .out file.",
                    )
                )
            continue
        conversion_targets = [(None, None), *[(name, None) for name in RAW_POWHEG_CHANNELS]]
        for channel, variant in conversion_targets:
            output_path = raw_powheg_output_path(campaign_dir, result, channel, variant)
            summary_csv = raw_powheg_summary_csv_path(campaign_dir, result) if channel is None else None
            legend = "Raw POWHEG" if channel is None else f"Raw POWHEG ({channel})"
            command = [
                "python3.10",
                str(script),
                "--output",
                str(output_path),
                "--analysis",
                "MC_DIS_BREIT",
                "--legend",
                legend,
                "--cross-section-pb",
                f"{measurement.value_pb:.16e}",
            ]
            if summary_csv is not None:
                command.extend(["--summary-csv", str(summary_csv)])
            if channel is not None:
                command.extend(["--channel", channel])
            for source in raw_sources:
                command.extend(["--log", source])
            if dry_run:
                print(shlex.join(command))
                built.append(
                    RawPOWHEGYODAResult(
                        logical_run=logical_run,
                        channel=channel,
                        variant=variant,
                        source_log=log_file,
                        source_out=out_file,
                        output_file=str(output_path),
                        summary_csv=str(summary_csv) if summary_csv is not None else None,
                        command=command,
                        returncode=0,
                    )
                )
                continue

            pending_tasks.append((logical_run, channel, variant, log_file, out_file, output_path, summary_csv, command))

    if pending_tasks and not dry_run:
        worker_count = max(1, min(max_workers, len(pending_tasks)))
        if worker_count == 1:
            for logical_run, channel, variant, log_file, out_file, output_path, summary_csv, command in pending_tasks:
                built.append(
                    run_raw_powheg_conversion_task(
                        base_dir,
                        output_dir,
                        logical_run,
                        channel,
                        variant,
                        log_file,
                        out_file,
                        output_path,
                        summary_csv,
                        command,
                    )
                )
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = [
                    executor.submit(
                        run_raw_powheg_conversion_task,
                        base_dir,
                        output_dir,
                        logical_run,
                        channel,
                        variant,
                        log_file,
                        out_file,
                        output_path,
                        summary_csv,
                        command,
                    )
                    for logical_run, channel, variant, log_file, out_file, output_path, summary_csv, command in pending_tasks
                ]
                for future in concurrent.futures.as_completed(futures):
                    built.append(future.result())

            order = {
                (
                    logical_run,
                    channel,
                    variant,
                    str(output_path),
                ): index
                for index, (logical_run, channel, variant, _log_file, _out_file, output_path, _summary_csv, _command) in enumerate(pending_tasks)
            }
            built.sort(key=lambda item: order.get((item.logical_run, item.channel, item.variant, str(item.output_file)), len(order)))

    return built


def collect_raw_powheg_inputs(
    campaign_dir: Path,
    tag: str,
    setup: str,
    order: str,
    channel: Optional[str] = None,
    variant: Optional[str] = None,
) -> Dict[str, List[Path]]:
    required_helicities = ("00", "PP", "PM") if setup == "GAMMA" else ("00", "PP", "PM", "MP", "MM")
    manifest = load_manifest(campaign_dir)
    found: Dict[str, List[Path]] = {hel: [] for hel in required_helicities}

    entries = manifest.get("raw_powheg_yoda", [])
    if isinstance(entries, list):
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            logical_run = str(entry.get("logical_run", ""))
            output_file = entry.get("output_file")
            returncode = entry.get("returncode")
            entry_channel = entry.get("channel")
            entry_variant = entry.get("variant")
            if not logical_run or not output_file or returncode != 0:
                continue
            if channel is None:
                if entry_channel not in (None, ""):
                    continue
            elif entry_channel != channel:
                continue
            if variant is None:
                if entry_variant not in (None, ""):
                    continue
            elif entry_variant != variant:
                continue
            for hel in required_helicities:
                prefix = f"DIS-POL-POWHEG_{hel}-{order}-{setup}"
                if logical_run.startswith(prefix):
                    found[hel].append(Path(str(output_file)))

    raw_dir = campaign_dir / "raw-powheg"
    if raw_dir.exists():
        for hel in required_helicities:
            if found[hel]:
                continue
            suffix = ".rawpowheg"
            if variant:
                suffix += f".{variant}"
            if channel:
                suffix += f".{channel}"
            pattern = f"DIS-POL-POWHEG_{hel}-{order}-{setup}*{tag}*{suffix}.yoda*"
            matches = sorted(
                path for path in raw_dir.glob(pattern)
                if path.name.endswith(".yoda") or path.name.endswith(".yoda.gz")
            )
            if matches:
                found[hel].extend(matches)

    missing = [hel for hel, paths in found.items() if not paths]
    if missing:
        channel_label = channel or "combined"
        if variant:
            channel_label = f"{channel_label}, variant={variant}"
        raise FileNotFoundError(
            f"Missing raw POWHEG YODA inputs for {order} setup {setup} ({channel_label}), helicities: {', '.join(missing)} in {campaign_dir / 'raw-powheg'}"
        )
    return {hel: sorted({path.resolve() for path in paths}) for hel, paths in found.items()}


def normalize_setups(setups: Sequence[str]) -> List[str]:
    if not setups:
        return ["ALL"]
    seen = set()
    ordered: List[str] = []
    ordering = {name: index for index, name in enumerate(("ALL", "GAMMA", "Z", "CC", *PS_SETUP_ORDER))}
    for setup in sorted(setups, key=lambda item: ordering.get(item, len(ordering))):
        if setup not in seen:
            ordered.append(setup)
            seen.add(setup)
    return ordered


def find_order_yoda_inputs(
    base_dir: Path,
    tag: str,
    setup: str,
    order: str,
    analysis_variant: Optional[str] = None,
    scale_variation: str = "nominal",
) -> Dict[str, Path]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    analysis_variant = analysis_variant if analysis_variant is not None else str(manifest.get("analysis_variant", ""))
    required_helicities = ("00", "PP", "PM") if setup == "GAMMA" else ("00", "PP", "PM", "MP", "MM")
    found: Dict[str, Path] = {}
    expected_runs = {
        hel: build_logical_run_stem(setup, order, hel, analysis_variant, scale_variation)
        for hel in required_helicities
    }

    if order == "NLO":
        manifest_entries = manifest.get("nlo_yoda", [])
    else:
        manifest_entries = manifest.get("merged_yoda", [])
    if isinstance(manifest_entries, list):
        for entry in manifest_entries:
            if not isinstance(entry, dict):
                continue
            logical_run = str(entry.get("logical_run", ""))
            output_file = entry.get("output_file")
            if not logical_run or not output_file:
                continue
            for hel in required_helicities:
                if logical_run == expected_runs[hel]:
                    found[hel] = Path(str(output_file))

    if len(found) == len(required_helicities):
        return found

    yoda_dir = campaign_dir / "yoda"
    if yoda_dir.exists():
        for hel in required_helicities:
            if hel in found:
                continue
            pattern = f"{expected_runs[hel]}-{tag}*.yoda*"
            matches = sorted(yoda_dir.glob(pattern))
            if matches:
                found[hel] = matches[0]

    missing = [hel for hel in required_helicities if hel not in found]
    if missing:
        raise FileNotFoundError(
            f"Missing {order} YODA inputs for setup {setup}, helicities: {', '.join(missing)} in {campaign_dir / 'yoda'}"
        )
    return found


def find_order_yoda_shards(
    base_dir: Path,
    tag: str,
    setup: str,
    order: str,
    analysis_variant: str,
    scale_variation: str = "nominal",
    scale_variations: bool = False,
) -> Dict[str, List[Path]]:
    required_helicities = ("00", "PP", "PM") if setup == "GAMMA" else ("00", "PP", "PM", "MP", "MM")
    found: Dict[str, List[Path]] = {hel: [] for hel in required_helicities}
    expected_runs = {
        hel: build_logical_run_stem(setup, order, hel, analysis_variant, scale_variation)
        for hel in required_helicities
    }
    results = collect_existing_yoda_results(
        base_dir,
        tag,
        analysis_variant,
        setups=[setup],
        scale_variations=scale_variations,
    )
    for item in results:
        stem = item.spec.job.stem
        for hel in required_helicities:
            if stem == expected_runs[hel]:
                for yoda_file in item.yoda_files:
                    found[hel].append(Path(yoda_file))
    missing = [hel for hel, paths in found.items() if not paths]
    if missing:
        raise FileNotFoundError(
            f"Missing raw shard YODA inputs for {order} setup {setup}, helicities: {', '.join(missing)} in {base_dir}"
        )
    return {hel: sorted({path.resolve() for path in paths}) for hel, paths in found.items()}


def yoda_has_analysis_objects(path: Path, analysis_name: str) -> bool:
    try:
        import yoda  # type: ignore
    except Exception:
        return False
    try:
        aos = yoda.read(str(path))
    except Exception:
        return False
    for label in ("Q2", "Pt", "Mjj"):
        if (
            f"/{analysis_name}/{label}" in aos
            or f"/RAW/{analysis_name}/{label}" in aos
            or any(_yoda_path_matches_analysis_label(obj_path, analysis_name, label) for obj_path in aos)
        ):
            return True
    return False


def merge_analysis_inputs(
    base_dir: Path,
    campaign_dir: Path,
    tag: str,
    setup: str,
    order: str,
    shard_inputs: Dict[str, List[Path]],
    merge_tool_preference: str,
    dry_run: bool,
) -> Dict[str, Path]:
    merge_tool = choose_yoda_merge_tool(merge_tool_preference, purpose="shard")
    if merge_tool is None:
        raise FileNotFoundError("No YODA merge tool found for analyze-herwig shard fallback")
    tool_name = Path(merge_tool).name
    tmp_dir = campaign_dir / "analysis" / "_inputs"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    merged: Dict[str, Path] = {}

    for hel, files in shard_inputs.items():
        output_path = tmp_dir / f"DIS-POL-POWHEG_{hel}-{order}-{setup}-{tag}.yoda"
        unique_files = [str(path) for path in sorted({path.resolve() for path in files})]
        if len(unique_files) == 1:
            if dry_run:
                print(f"[dry-run] cp {unique_files[0]} {output_path}")
            else:
                shutil.copy2(unique_files[0], output_path)
            merged[hel] = output_path
            continue

        command = [merge_tool]
        if tool_name == "rivet-merge":
            command.append("-e")
        command.extend([*unique_files, "-o", str(output_path)])
        if dry_run:
            print(shlex.join(command))
        else:
            env = rivet_merge_env(base_dir) if tool_name == "rivet-merge" else None
            proc = subprocess.run(command, cwd=base_dir, env=env, text=True, capture_output=True)
            stem = f"DIS-POL-POWHEG_{hel}-{order}-{setup}-{tag}"
            if proc.stdout:
                (tmp_dir / f"{stem}.stdout").write_text(proc.stdout)
            if proc.stderr:
                (tmp_dir / f"{stem}.stderr").write_text(proc.stderr)
            if proc.returncode != 0:
                raise RuntimeError(
                    f"Failed to merge shard YODAs for {hel} {order} {setup} with rc={proc.returncode}"
                )
        merged[hel] = output_path

    return merged


def analyze_herwig_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    analysis_variant: str = "",
    scale_variations: bool = False,
    merge_tool_preference: str = "auto",
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    analysis_variant = analysis_variant or str(manifest.get("analysis_variant", ""))
    analysis_dir = campaign_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    outputs: Dict[str, str] = {}
    plot_outputs: Dict[str, str] = {}
    selected = normalize_setups(setups)
    worker_count = max(1, min(max_workers, len(selected))) if selected else 1
    mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
    print_stage(f"Constructing analyzed Herwig YODAs for setups: {', '.join(selected)} ({mode})")

    def analyze_setup(setup: str) -> tuple[str, str, str]:
        legend_label = herwig_legend_label(analysis_variant, setup)
        variation_names = [name for name in selected_scale_variations(scale_variations) if name != "nominal"]

        if dry_run:
            output_path = analysis_dir / f"Herwig_{analysis_name}_{setup}_NLO_polarized.yoda.gz"
            plot_path = analysis_dir / f"Herwig_{analysis_name}_{setup}_NLO_polarized_plot.yoda.gz"
            print(f"[dry-run] analyze-herwig {setup} -> {output_path}")
            print(f"[dry-run] analyze-herwig-plot {setup} -> {plot_path}")
            return setup, str(output_path), str(plot_path)

        def build_objects_from_inputs(pos_inputs: Dict[str, object], neg_inputs: Dict[str, object], shard_mode: bool) -> Dict[str, object]:
            def one(mapping: Dict[str, object], helicity: str) -> Optional[object]:
                value = mapping.get(helicity)
                if value is None:
                    return None
                if shard_mode:
                    return [str(path) for path in value]
                return str(value)

            return build_dis_polarized_objects(
                setup=setup,
                zero_path=one(pos_inputs, "00"),
                pp_path=one(pos_inputs, "PP"),
                pm_path=one(pos_inputs, "PM"),
                mp_path=one(pos_inputs, "MP"),
                mm_path=one(pos_inputs, "MM"),
                zero_subtract_path=one(neg_inputs, "00"),
                pp_subtract_path=one(neg_inputs, "PP"),
                pm_subtract_path=one(neg_inputs, "PM"),
                mp_subtract_path=one(neg_inputs, "MP"),
                mm_subtract_path=one(neg_inputs, "MM"),
                analysis=analysis_name,
            )

        use_shard_fallback = False
        try:
            pos_inputs = find_order_yoda_inputs(base_dir, tag, setup, "POSNLO", scale_variation="nominal")
            neg_inputs = find_order_yoda_inputs(base_dir, tag, setup, "NEGNLO", scale_variation="nominal")
            variation_pos_inputs = {
                variation_name: find_order_yoda_inputs(base_dir, tag, setup, "POSNLO", scale_variation=variation_name)
                for variation_name in variation_names
            }
            variation_neg_inputs = {
                variation_name: find_order_yoda_inputs(base_dir, tag, setup, "NEGNLO", scale_variation=variation_name)
                for variation_name in variation_names
            }
            if not yoda_has_analysis_objects(pos_inputs["PP"], analysis_name):
                use_shard_fallback = True
            for variation_name in variation_names:
                if not yoda_has_analysis_objects(variation_pos_inputs[variation_name]["PP"], analysis_name):
                    use_shard_fallback = True
                    break
        except Exception:
            use_shard_fallback = True

        if use_shard_fallback:
            pos_inputs = find_order_yoda_shards(
                base_dir,
                tag,
                setup,
                "POSNLO",
                analysis_variant,
                scale_variation="nominal",
                scale_variations=scale_variations,
            )
            neg_inputs = find_order_yoda_shards(
                base_dir,
                tag,
                setup,
                "NEGNLO",
                analysis_variant,
                scale_variation="nominal",
                scale_variations=scale_variations,
            )
            variation_pos_inputs = {
                variation_name: find_order_yoda_shards(
                    base_dir,
                    tag,
                    setup,
                    "POSNLO",
                    analysis_variant,
                    scale_variation=variation_name,
                    scale_variations=scale_variations,
                )
                for variation_name in variation_names
            }
            variation_neg_inputs = {
                variation_name: find_order_yoda_shards(
                    base_dir,
                    tag,
                    setup,
                    "NEGNLO",
                    analysis_variant,
                    scale_variation=variation_name,
                    scale_variations=scale_variations,
                )
                for variation_name in variation_names
            }

        objects = build_objects_from_inputs(pos_inputs, neg_inputs, use_shard_fallback)
        apply_legend_annotation(objects, legend_label)
        apply_setup_style_annotation(objects, setup)
        output_path = analysis_dir / f"Herwig_{analysis_name}_{setup}_NLO_polarized.yoda.gz"
        write_analysis_yoda_gz(objects, str(output_path))
        plot_path = analysis_dir / f"Herwig_{analysis_name}_{setup}_NLO_polarized_plot.yoda.gz"
        plot_objects = build_plot_scatter_objects(objects)
        apply_legend_annotation(plot_objects, legend_label)
        apply_setup_style_annotation(plot_objects, setup)
        if variation_names:
            combined_plot_objects = dict(plot_objects)
            for variation_name in variation_names:
                variation_objects = build_objects_from_inputs(
                    variation_pos_inputs[variation_name],
                    variation_neg_inputs.get(variation_name, {}),
                    use_shard_fallback,
                )
                variation_plot_objects = build_plot_scatter_objects(variation_objects)
                apply_legend_annotation(variation_plot_objects, legend_label)
                apply_setup_style_annotation(variation_plot_objects, setup)
                variation_plot_objects, skipped_variation_paths = filter_compatible_variation_plot_objects(
                    plot_objects,
                    variation_plot_objects,
                )
                for skipped_path, reason in sorted(skipped_variation_paths.items()):
                    print(
                        f"[warn] Skipping {variation_name} plotted variation for {setup} "
                        f"at {skipped_path}: {reason}",
                        flush=True,
                    )
                combined_plot_objects.update(clone_objects_with_variation(variation_plot_objects, variation_name))
            plot_objects = combined_plot_objects
        write_analysis_yoda_gz(plot_objects, str(plot_path))
        return setup, str(output_path), str(plot_path)

    for setup, output_path, plot_path in run_parallel_ordered(
        selected,
        analyze_setup,
        1 if dry_run else max_workers,
    ):
        outputs[setup] = output_path
        plot_outputs[setup] = plot_path

    if outputs and not dry_run:
        update_manifest_outputs(campaign_dir, analysis_outputs=outputs, analysis_plot_outputs=plot_outputs)
    return outputs


def analyze_ps_herwig_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    scale_variations: bool = False,
    merge_tool_preference: str = "auto",
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, Dict[str, str]]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    analysis_dir = campaign_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    outputs: Dict[str, Dict[str, str]] = {}
    plot_outputs: Dict[str, Dict[str, str]] = {}
    selected = [setup for setup in normalize_campaign_setups(setups) if is_ps_setup(setup)]
    tasks = [(setup, family) for setup in selected for family in ps_families_for_setup(setup)]
    worker_count = max(1, min(max_workers, len(tasks))) if tasks else 1
    mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
    print_stage(
        "Constructing analyzed PS Herwig YODAs for "
        + ", ".join(f"{setup}/{family}" for setup, family in tasks)
        + f" ({mode})"
    )

    def analyze_task(task: tuple[str, str]) -> tuple[str, str, str, str]:
        setup, family = task
        legend_label = herwig_legend_label(family, setup)
        variation_names = [name for name in selected_scale_variations(scale_variations) if name != "nominal"]

        if dry_run:
            output_path = analysis_dir / f"Herwig_{analysis_name}_{family}_{setup}_NLO_polarized.yoda.gz"
            plot_path = analysis_dir / f"Herwig_{analysis_name}_{family}_{setup}_NLO_polarized_plot.yoda.gz"
            print(f"[dry-run] analyze-herwig-ps {setup}/{family} -> {output_path}")
            print(f"[dry-run] analyze-herwig-ps-plot {setup}/{family} -> {plot_path}")
            return setup, family, str(output_path), str(plot_path)

        def build_objects_from_inputs(
            pos_inputs: Dict[str, object],
            neg_inputs: Dict[str, object],
            shard_mode: bool,
        ) -> Dict[str, object]:
            def one(mapping: Dict[str, object], helicity: str) -> Optional[object]:
                value = mapping.get(helicity)
                if value is None:
                    return None
                if shard_mode:
                    return [str(path) for path in value]
                return str(value)

            return build_dis_polarized_objects(
                setup=setup,
                zero_path=one(pos_inputs, "00"),
                pp_path=one(pos_inputs, "PP"),
                pm_path=one(pos_inputs, "PM"),
                mp_path=one(pos_inputs, "MP"),
                mm_path=one(pos_inputs, "MM"),
                zero_subtract_path=one(neg_inputs, "00"),
                pp_subtract_path=one(neg_inputs, "PP"),
                pm_subtract_path=one(neg_inputs, "PM"),
                mp_subtract_path=one(neg_inputs, "MP"),
                mm_subtract_path=one(neg_inputs, "MM"),
                analysis=analysis_name,
            )

        use_shard_fallback = False
        try:
            pos_inputs = find_order_yoda_inputs(
                base_dir,
                tag,
                setup,
                "POSNLO",
                analysis_variant=family,
                scale_variation="nominal",
            )
            neg_inputs = find_order_yoda_inputs(
                base_dir,
                tag,
                setup,
                "NEGNLO",
                analysis_variant=family,
                scale_variation="nominal",
            )
            variation_pos_inputs = {
                variation_name: find_order_yoda_inputs(
                    base_dir,
                    tag,
                    setup,
                    "POSNLO",
                    analysis_variant=family,
                    scale_variation=variation_name,
                )
                for variation_name in variation_names
            }
            variation_neg_inputs = {
                variation_name: find_order_yoda_inputs(
                    base_dir,
                    tag,
                    setup,
                    "NEGNLO",
                    analysis_variant=family,
                    scale_variation=variation_name,
                )
                for variation_name in variation_names
            }
            if not yoda_has_analysis_objects(pos_inputs["PP"], analysis_name):
                use_shard_fallback = True
            for variation_name in variation_names:
                if not yoda_has_analysis_objects(variation_pos_inputs[variation_name]["PP"], analysis_name):
                    use_shard_fallback = True
                    break
        except Exception:
            use_shard_fallback = True

        if use_shard_fallback:
            pos_inputs = find_order_yoda_shards(
                base_dir,
                tag,
                setup,
                "POSNLO",
                family,
                scale_variation="nominal",
                scale_variations=scale_variations,
            )
            neg_inputs = find_order_yoda_shards(
                base_dir,
                tag,
                setup,
                "NEGNLO",
                family,
                scale_variation="nominal",
                scale_variations=scale_variations,
            )
            variation_pos_inputs = {
                variation_name: find_order_yoda_shards(
                    base_dir,
                    tag,
                    setup,
                    "POSNLO",
                    family,
                    scale_variation=variation_name,
                    scale_variations=scale_variations,
                )
                for variation_name in variation_names
            }
            variation_neg_inputs = {
                variation_name: find_order_yoda_shards(
                    base_dir,
                    tag,
                    setup,
                    "NEGNLO",
                    family,
                    scale_variation=variation_name,
                    scale_variations=scale_variations,
                )
                for variation_name in variation_names
            }

        objects = build_objects_from_inputs(pos_inputs, neg_inputs, use_shard_fallback)
        apply_legend_annotation(objects, legend_label)
        apply_variant_style_annotation(objects, family)
        output_path = analysis_dir / f"Herwig_{analysis_name}_{family}_{setup}_NLO_polarized.yoda.gz"
        write_analysis_yoda_gz(objects, str(output_path))

        plot_objects = build_plot_scatter_objects(objects)
        apply_legend_annotation(plot_objects, legend_label)
        apply_variant_style_annotation(plot_objects, family)
        if variation_names:
            combined_plot_objects = dict(plot_objects)
            for variation_name in variation_names:
                variation_objects = build_objects_from_inputs(
                    variation_pos_inputs[variation_name],
                    variation_neg_inputs.get(variation_name, {}),
                    use_shard_fallback,
                )
                variation_plot_objects = build_plot_scatter_objects(variation_objects)
                apply_legend_annotation(variation_plot_objects, legend_label)
                apply_variant_style_annotation(variation_plot_objects, family)
                variation_plot_objects, skipped_variation_paths = filter_compatible_variation_plot_objects(
                    plot_objects,
                    variation_plot_objects,
                )
                for skipped_path, reason in sorted(skipped_variation_paths.items()):
                    print(
                        f"[warn] Skipping {variation_name} plotted variation for {setup}/{family} "
                        f"at {skipped_path}: {reason}",
                        flush=True,
                    )
                combined_plot_objects.update(clone_objects_with_variation(variation_plot_objects, variation_name))
            plot_objects = combined_plot_objects
        plot_path = analysis_dir / f"Herwig_{analysis_name}_{family}_{setup}_NLO_polarized_plot.yoda.gz"
        write_analysis_yoda_gz(plot_objects, str(plot_path))
        return setup, family, str(output_path), str(plot_path)

    for setup, family, output_path, plot_path in run_parallel_ordered(
        tasks,
        analyze_task,
        1 if dry_run else max_workers,
    ):
        outputs.setdefault(family, {})[setup] = output_path
        plot_outputs.setdefault(family, {})[setup] = plot_path

    if outputs and not dry_run:
        update_manifest_outputs(
            campaign_dir,
            analysis_outputs_by_family=outputs,
            analysis_plot_outputs_by_family=plot_outputs,
        )
    return outputs


def analyze_raw_powheg_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    analysis_variant: str = "",
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    if analysis_variant != "RIVETFO":
        return {}

    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    if not manifest.get("raw_powheg"):
        return {}

    analysis_dir = campaign_dir / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    outputs: Dict[str, str] = {}
    plot_outputs: Dict[str, str] = {}
    channel_outputs: Dict[str, Dict[str, str]] = {channel: {} for channel in RAW_POWHEG_CHANNELS}
    channel_plot_outputs: Dict[str, Dict[str, str]] = {channel: {} for channel in RAW_POWHEG_CHANNELS}
    selected = normalize_setups(setups)
    tasks: List[tuple[str, Optional[str], Optional[str]]] = [(setup, None, None) for setup in selected]
    tasks.extend((setup, channel, None) for setup in selected for channel in RAW_POWHEG_CHANNELS)
    worker_count = max(1, min(max_workers, len(tasks))) if tasks else 1
    mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
    print_stage(f"Constructing analyzed raw POWHEG YODAs for {len(tasks)} setup/channel work units ({mode})")

    def analyze_raw_task(task: tuple[str, Optional[str], Optional[str]]) -> tuple[str, Optional[str], Optional[str], str, str]:
        setup, channel, variant = task
        if dry_run:
            if channel is None and variant is None:
                output_path = analysis_dir / f"RawPOWHEG_{analysis_name}_{setup}_NLO_polarized.yoda.gz"
                plot_path = analysis_dir / f"RawPOWHEG_{analysis_name}_{setup}_NLO_polarized_plot.yoda.gz"
                print(f"[dry-run] analyze-raw-powheg {setup} -> {output_path}")
                print(f"[dry-run] analyze-raw-powheg-plot {setup} -> {plot_path}")
                return setup, None, None, str(output_path), str(plot_path)
            channel_output_path = analysis_dir / f"RawPOWHEG_{channel}_{analysis_name}_{setup}_NLO_polarized.yoda.gz"
            channel_plot_path = analysis_dir / f"RawPOWHEG_{channel}_{analysis_name}_{setup}_NLO_polarized_plot.yoda.gz"
            print(f"[dry-run] analyze-raw-powheg-{channel.lower()} {setup} -> {channel_output_path}")
            print(f"[dry-run] analyze-raw-powheg-{channel.lower()}-plot {setup} -> {channel_plot_path}")
            return setup, channel, None, str(channel_output_path), str(channel_plot_path)

        if channel is None and variant is None:
            pos_inputs = collect_raw_powheg_inputs(campaign_dir, tag, setup, "POSNLO")
            neg_inputs = collect_raw_powheg_inputs(campaign_dir, tag, setup, "NEGNLO")
            objects = build_dis_polarized_objects(
                setup=setup,
                zero_path=[str(path) for path in pos_inputs["00"]],
                pp_path=[str(path) for path in pos_inputs["PP"]],
                pm_path=[str(path) for path in pos_inputs["PM"]],
                mp_path=[str(path) for path in pos_inputs.get("MP", [])] if "MP" in pos_inputs else None,
                mm_path=[str(path) for path in pos_inputs.get("MM", [])] if "MM" in pos_inputs else None,
                zero_subtract_path=[str(path) for path in neg_inputs["00"]],
                pp_subtract_path=[str(path) for path in neg_inputs["PP"]],
                pm_subtract_path=[str(path) for path in neg_inputs["PM"]],
                mp_subtract_path=[str(path) for path in neg_inputs.get("MP", [])] if "MP" in neg_inputs else None,
                mm_subtract_path=[str(path) for path in neg_inputs.get("MM", [])] if "MM" in neg_inputs else None,
                analysis=analysis_name,
            )
            apply_legend_annotation(objects, "Raw POWHEG")
            output_path = analysis_dir / f"RawPOWHEG_{analysis_name}_{setup}_NLO_polarized.yoda.gz"
            write_analysis_yoda_gz(objects, str(output_path))
            plot_objects = build_plot_scatter_objects(objects)
            apply_legend_annotation(plot_objects, "Raw POWHEG")
            plot_path = analysis_dir / f"RawPOWHEG_{analysis_name}_{setup}_NLO_polarized_plot.yoda.gz"
            write_analysis_yoda_gz(plot_objects, str(plot_path))
            return setup, None, None, str(output_path), str(plot_path)

        pos_channel_inputs = collect_raw_powheg_inputs(campaign_dir, tag, setup, "POSNLO", channel=channel)
        neg_channel_inputs = collect_raw_powheg_inputs(campaign_dir, tag, setup, "NEGNLO", channel=channel)
        channel_objects = build_dis_polarized_objects(
            setup=setup,
            zero_path=[str(path) for path in pos_channel_inputs["00"]],
            pp_path=[str(path) for path in pos_channel_inputs["PP"]],
            pm_path=[str(path) for path in pos_channel_inputs["PM"]],
            mp_path=[str(path) for path in pos_channel_inputs.get("MP", [])] if "MP" in pos_channel_inputs else None,
            mm_path=[str(path) for path in pos_channel_inputs.get("MM", [])] if "MM" in pos_channel_inputs else None,
            zero_subtract_path=[str(path) for path in neg_channel_inputs["00"]],
            pp_subtract_path=[str(path) for path in neg_channel_inputs["PP"]],
            pm_subtract_path=[str(path) for path in neg_channel_inputs["PM"]],
            mp_subtract_path=[str(path) for path in neg_channel_inputs.get("MP", [])] if "MP" in neg_channel_inputs else None,
            mm_subtract_path=[str(path) for path in neg_channel_inputs.get("MM", [])] if "MM" in neg_channel_inputs else None,
            analysis=analysis_name,
        )
        apply_legend_annotation(channel_objects, f"Raw POWHEG ({channel})")
        channel_output_path = analysis_dir / f"RawPOWHEG_{channel}_{analysis_name}_{setup}_NLO_polarized.yoda.gz"
        write_analysis_yoda_gz(channel_objects, str(channel_output_path))
        channel_plot_objects = build_plot_scatter_objects(channel_objects)
        apply_legend_annotation(channel_plot_objects, f"Raw POWHEG ({channel})")
        channel_plot_path = analysis_dir / f"RawPOWHEG_{channel}_{analysis_name}_{setup}_NLO_polarized_plot.yoda.gz"
        write_analysis_yoda_gz(channel_plot_objects, str(channel_plot_path))
        return setup, channel, None, str(channel_output_path), str(channel_plot_path)

    for setup, channel, variant, output_path, plot_path in run_parallel_ordered(
        tasks,
        analyze_raw_task,
        1 if dry_run else max_workers,
    ):
        if channel is None and variant is None:
            outputs[setup] = output_path
            plot_outputs[setup] = plot_path
        else:
            channel_outputs[channel][setup] = output_path
            channel_plot_outputs[channel][setup] = plot_path

    if (outputs or any(channel_outputs[channel] for channel in RAW_POWHEG_CHANNELS)) and not dry_run:
        update_manifest_outputs(
            campaign_dir,
            raw_powheg_analysis_outputs=outputs,
            raw_powheg_analysis_plot_outputs=plot_outputs,
            raw_powheg_channel_analysis_outputs=channel_outputs,
            raw_powheg_channel_analysis_plot_outputs=channel_plot_outputs,
        )
    return outputs


def use_cc_poldis_reference(setups: Sequence[str], analysis_variant: str) -> bool:
    selected = normalize_setups(setups)
    return analysis_variant == "RIVETFO" and selected == ["CC"]


def resolve_poldis_top_paths(
    base_dir: Path,
    unpol: Optional[Path],
    pol: Optional[Path],
    setups: Sequence[str],
    analysis_variant: str,
) -> tuple[Path, Path]:
    selected = normalize_setups(setups)
    has_cc = "CC" in selected
    has_nc = any(setup != "CC" for setup in selected)
    if analysis_variant == "RIVETFO" and has_cc and has_nc and (unpol is None or pol is None):
        raise ValueError(
            "Mixed NC+CC RivetFO comparisons require explicit --unpol and --pol POLDIS .top files."
        )

    if use_cc_poldis_reference(selected, analysis_variant):
        default_unpol = base_dir / "dijets_unpol_cc.top"
        default_pol = base_dir / "dijets_pol_cc.top"
    else:
        default_unpol = base_dir / "dijets_unpol.top"
        default_pol = base_dir / "dijets_pol.top"

    unpol_path = (unpol or default_unpol).resolve()
    pol_path = (pol or default_pol).resolve()
    missing = [path for path in (unpol_path, pol_path) if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing POLDIS .top file(s): " + ", ".join(str(path) for path in missing)
        )
    return unpol_path, pol_path


def build_poldis_reference(
    base_dir: Path,
    unpol: Optional[Path],
    pol: Optional[Path],
    analysis_name: str,
    order: str,
    out: Optional[Path],
    setups: Sequence[str] = (),
    analysis_variant: str = "",
    campaign_tag: Optional[str] = None,
    dry_run: bool = False,
) -> Path:
    unpol_path, pol_path = resolve_poldis_top_paths(base_dir, unpol, pol, setups, analysis_variant)

    if out is not None:
        output_path = out.resolve()
    elif campaign_tag:
        output_path = resolve_campaign_dir(base_dir, campaign_tag) / "refs" / f"POLDIS_{analysis_name}_ref.yoda.gz"
    else:
        output_path = (base_dir / f"POLDIS_{analysis_name}_ref.yoda.gz").resolve()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    print_stage(f"Building POLDIS reference YODA from {unpol_path.name} and {pol_path.name}")
    if dry_run:
        print(f"[dry-run] poldis-top {unpol_path} {pol_path} -> {output_path}")
        return output_path

    convert_topdrawer_to_yoda(
        unpol=str(unpol_path),
        pol=str(pol_path),
        out=str(output_path),
        analysis=analysis_name,
        order=order,
    )
    if campaign_tag:
        update_manifest_outputs(resolve_campaign_dir(base_dir, campaign_tag), reference_output=str(output_path))
    return output_path


def run_spin_diagnostic_extractor(
    base_dir: Path,
    tag: str,
    dry_run: bool = False,
) -> int:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    cmd = [
        "python3.10",
        str(script_path("extract_powheg_real_spin_diagnostics.py")),
        "--base-dir",
        str(base_dir),
        "-t",
        tag,
        "--strict-tag",
        "--json-out",
        str(campaign_dir / "spin_diagnostic.json"),
        "--csv-out",
        str(campaign_dir / "spin_diagnostic.csv"),
    ]
    print_stage("Extracting POWHEG real-spin diagnostics")
    return run_extractor(cmd, campaign_dir / "spin_diagnostic.txt", base_dir, dry_run)


def choose_rivet_mkhtml_tool(preference: str) -> Optional[str]:
    return shutil.which(preference)


def default_plot_dir(base_dir: Path, tag: str, setup: str, order: str) -> Path:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    return campaign_dir / f"plots_mc_vs_poldis_{setup.lower()}_{order.lower()}"


def default_ps_plot_dir(base_dir: Path, tag: str, setup: str) -> Path:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    return campaign_dir / f"plots_ps_spincomp_{setup.lower()}"


def default_ps_hadronization_plot_dir(base_dir: Path, tag: str) -> Path:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    return campaign_dir / "plots_ps_hadronization_full_toggle"


def run_ps_rivetplot_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    rivet_mkhtml_tool: str,
    plot_dir: Optional[Path] = None,
    scale_variations: Optional[bool] = None,
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    scale_variations_enabled = resolve_scale_variations_requested(scale_variations, manifest, default=False)
    analysis_outputs = manifest.get("analysis_yoda_by_family", {})
    if not isinstance(analysis_outputs, dict):
        analysis_outputs = {}
    analysis_plot_outputs = manifest.get("analysis_plot_yoda_by_family", {})
    if not isinstance(analysis_plot_outputs, dict):
        analysis_plot_outputs = {}

    tool = choose_rivet_mkhtml_tool(rivet_mkhtml_tool)
    if tool is None:
        raise FileNotFoundError(f"Could not find rivet-mkhtml tool '{rivet_mkhtml_tool}'")

    selected = multi_family_ps_setups(setups)
    if not selected:
        raise ValueError(
            "PS comparison plots are only defined for multi-family PS setups (currently SPINCOMP and SPINHAD)."
        )
    if plot_dir is not None and len(selected) != 1:
        raise ValueError("--plot-dir may only be used with a single --setup")

    outputs: Dict[str, str] = {}
    worker_count = max(1, min(max_workers, len(selected))) if selected else 1
    mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
    print_stage(f"Building PS spin-comparison plots for setups: {', '.join(selected)} ({mode})")

    def plot_setup(setup: str) -> tuple[str, str]:
        this_plot_dir = plot_dir.resolve() if plot_dir is not None else default_ps_plot_dir(base_dir, tag, setup)
        this_plot_dir.parent.mkdir(parents=True, exist_ok=True)
        sanitized_inputs_dir = campaign_dir / "analysis" / "_rivetplot_inputs_ps"
        sanitized_inputs_dir.mkdir(parents=True, exist_ok=True)

        def family_yoda(family: str) -> str:
            family_plot_map = analysis_plot_outputs.get(family, {})
            if not isinstance(family_plot_map, dict):
                family_plot_map = {}
            family_yoda_map = analysis_outputs.get(family, {})
            if not isinstance(family_yoda_map, dict):
                family_yoda_map = {}
            yoda_path = family_plot_map.get(setup) or family_yoda_map.get(setup)
            if not yoda_path:
                raise FileNotFoundError(
                    f"Missing analyzed PS Herwig YODA for {setup}/{family} in {campaign_dir / 'manifest.json'}"
                )
            return str(yoda_path)

        def sanitize_input_yoda(input_yoda: str, label: str, scale_envelope: bool = False) -> str:
            input_path = Path(input_yoda)
            output_name = input_path.name
            if output_name.endswith(".yoda.gz"):
                output_name = output_name[:-8] + f".{label}.sanitized.yoda.gz"
            elif output_name.endswith(".yoda"):
                output_name = output_name[:-5] + f".{label}.sanitized.yoda"
            else:
                output_name = output_name + f".{label}.sanitized.yoda"
            sanitized_path = sanitized_inputs_dir / output_name
            if not dry_run:
                if scale_envelope:
                    kept_count, dropped_count, dropped = build_scale_envelope_plot_yoda(input_path, sanitized_path)
                    print(
                        f"[stage] Built scale-envelope plot YODA for {label} in {setup}: kept {kept_count} objects, "
                        f"dropped {dropped_count} incompatible/empty objects -> {sanitized_path}",
                        flush=True,
                    )
                else:
                    kept_count, dropped_count, dropped = sanitize_plot_yoda(
                        input_path,
                        sanitized_path,
                        drop_band_annotations=not scale_envelope,
                    )
                    print(
                        f"[stage] Sanitized {label} plot YODA for {setup}: kept {kept_count} objects, "
                        f"dropped {dropped_count} incompatible/empty objects -> {sanitized_path}",
                        flush=True,
                    )
                if dropped_count:
                    for dropped_path, reason in sorted(dropped.items())[:20]:
                        print(f"[warn]   {dropped_path}: {reason}", flush=True)
                    if dropped_count > 20:
                        print(f"[warn]   ... and {dropped_count - 20} more", flush=True)
            return str(sanitized_path)

        sanitized_inputs: Dict[str, str] = {}
        for family in ps_families_for_setup(setup):
            sanitized_input = sanitize_input_yoda(
                family_yoda(family),
                f"{setup.lower()}_{family.lower()}",
                scale_envelope=scale_variations_enabled,
            )
            sanitized_inputs[family] = sanitized_input

        if not dry_run and not scale_variations_enabled and len(sanitized_inputs) > 1:
            harmonized_files, harmonized_objects = harmonize_plot_yoda_bin_grids(list(sanitized_inputs.values()))
            if harmonized_files or harmonized_objects:
                print_stage(
                    f"Harmonized PS plot bin grids for {setup}: "
                    f"updated {harmonized_objects} objects across {harmonized_files} files"
                )

        file_args: List[str] = []
        for family in ps_families_for_setup(setup):
            file_args.append(f"{sanitized_inputs[family]}:Title={ps_family_label(family)}")

        command: List[str] = [
            sys.executable,
            str(script_path("rivet_mkhtml_safe.py")),
            tool,
            *file_args,
            "--verbose",
            "-o",
            str(this_plot_dir),
        ]
        if dry_run:
            print(shlex.join(command))
        else:
            if this_plot_dir.exists():
                shutil.rmtree(this_plot_dir)
            this_plot_dir.mkdir(parents=True, exist_ok=True)
            env = os.environ.copy()
            if not scale_variations_enabled:
                env["DISPOL_FORCE_NO_ERROR_BANDS"] = "1"
            proc = subprocess.run(command, cwd=base_dir, text=True, capture_output=True, env=env)
            if proc.stdout:
                (this_plot_dir.parent / f"{this_plot_dir.name}.mkhtml.stdout").write_text(proc.stdout)
            if proc.stderr:
                (this_plot_dir.parent / f"{this_plot_dir.name}.mkhtml.stderr").write_text(proc.stderr)
            if proc.returncode != 0:
                raise RuntimeError(f"rivet-mkhtml failed for {setup} with rc={proc.returncode}")
            if scale_variations_enabled:
                patched_count, rerendered_count = rewrite_scale_envelope_plot_scripts(this_plot_dir)
                print_stage(
                    f"Postprocessed scale-envelope Rivet scripts for {setup}: "
                    f"patched {patched_count}, rerendered {rerendered_count}"
                )
            else:
                patched_count, rerendered_count = rewrite_no_scale_ratio_plot_scripts(this_plot_dir)
                print_stage(
                    f"Postprocessed no-scale Rivet ratio scripts for {setup}: "
                    f"patched {patched_count}, rerendered {rerendered_count}"
                )
        return setup, str(this_plot_dir)

    for setup, plot_output in run_parallel_ordered(
        selected,
        plot_setup,
        1 if dry_run else max_workers,
    ):
        outputs[setup] = plot_output

    if outputs and not dry_run:
        update_manifest_outputs(campaign_dir, plot_outputs=outputs)
    return outputs


def run_ps_hadronization_rivetplot_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    rivet_mkhtml_tool: str,
    scale_variations: Optional[bool] = None,
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    scale_variations_enabled = resolve_scale_variations_requested(scale_variations, manifest, default=False)
    analysis_outputs = manifest.get("analysis_yoda_by_family", {})
    if not isinstance(analysis_outputs, dict):
        analysis_outputs = {}
    analysis_plot_outputs = manifest.get("analysis_plot_yoda_by_family", {})
    if not isinstance(analysis_plot_outputs, dict):
        analysis_plot_outputs = {}

    tool = choose_rivet_mkhtml_tool(rivet_mkhtml_tool)
    if tool is None:
        raise FileNotFoundError(f"Could not find rivet-mkhtml tool '{rivet_mkhtml_tool}'")
    if not has_ps_hadronization_toggle_inputs(setups):
        raise ValueError(
            "PS hadronization toggle plots require both --setup SPINCOMP and --setup SPINHAD in the same campaign."
        )

    this_plot_dir = default_ps_hadronization_plot_dir(base_dir, tag)
    this_plot_dir.parent.mkdir(parents=True, exist_ok=True)
    sanitized_inputs_dir = campaign_dir / "analysis" / "_rivetplot_inputs_ps_hadronization"
    sanitized_inputs_dir.mkdir(parents=True, exist_ok=True)

    def family_yoda(setup: str, family: str) -> str:
        family_plot_map = analysis_plot_outputs.get(family, {})
        if not isinstance(family_plot_map, dict):
            family_plot_map = {}
        family_yoda_map = analysis_outputs.get(family, {})
        if not isinstance(family_yoda_map, dict):
            family_yoda_map = {}
        yoda_path = family_plot_map.get(setup) or family_yoda_map.get(setup)
        if not yoda_path:
            raise FileNotFoundError(
                f"Missing analyzed PS Herwig YODA for {setup}/{family} in {campaign_dir / 'manifest.json'}"
            )
        return str(yoda_path)

    def sanitize_input_yoda(input_yoda: str, label: str, scale_envelope: bool = False) -> str:
        input_path = Path(input_yoda)
        output_name = input_path.name
        if output_name.endswith(".yoda.gz"):
            output_name = output_name[:-8] + f".{label}.sanitized.yoda.gz"
        elif output_name.endswith(".yoda"):
            output_name = output_name[:-5] + f".{label}.sanitized.yoda"
        else:
            output_name = output_name + f".{label}.sanitized.yoda"
        sanitized_path = sanitized_inputs_dir / output_name
        if not dry_run:
            if scale_envelope:
                kept_count, dropped_count, dropped = build_scale_envelope_plot_yoda(input_path, sanitized_path)
                print(
                    f"[stage] Built hadronization plot YODA for {label}: kept {kept_count} objects, "
                    f"dropped {dropped_count} incompatible/empty objects -> {sanitized_path}",
                    flush=True,
                )
            else:
                kept_count, dropped_count, dropped = sanitize_plot_yoda(
                    input_path,
                    sanitized_path,
                    drop_band_annotations=not scale_envelope,
                )
                print(
                    f"[stage] Sanitized hadronization plot YODA for {label}: kept {kept_count} objects, "
                    f"dropped {dropped_count} incompatible/empty objects -> {sanitized_path}",
                    flush=True,
                )
            if dropped_count:
                for dropped_path, reason in sorted(dropped.items())[:20]:
                    print(f"[warn]   {dropped_path}: {reason}", flush=True)
                if dropped_count > 20:
                    print(f"[warn]   ... and {dropped_count - 20} more", flush=True)
        return str(sanitized_path)

    no_had_input = sanitize_input_yoda(
        family_yoda("SPINCOMP", "RIVETPS-SPIN"),
        "spincomp_rivetps_spin",
        scale_envelope=scale_variations_enabled,
    )
    had_input = sanitize_input_yoda(
        family_yoda("SPINHAD", "RIVETPS-SPIN"),
        "spinhad_rivetps_spin",
        scale_envelope=scale_variations_enabled,
    )

    if not dry_run and not scale_variations_enabled:
        harmonized_files, harmonized_objects = harmonize_plot_yoda_bin_grids([no_had_input, had_input])
        if harmonized_files or harmonized_objects:
            print_stage(
                "Harmonized PS hadronization-toggle plot bin grids: "
                f"updated {harmonized_objects} objects across {harmonized_files} files"
            )

    no_had_plot_input = no_had_input
    had_plot_input = had_input
    if not dry_run:
        no_had_plot_input = rewrite_plot_yoda_annotations(
            no_had_input,
            str(Path(no_had_input).with_name(Path(no_had_input).name.replace('.sanitized.', '.styled.'))),
            "Full (No Hadronization)",
            line_color="red",
            line_style="solid",
        )
        had_plot_input = rewrite_plot_yoda_annotations(
            had_input,
            str(Path(had_input).with_name(Path(had_input).name.replace('.sanitized.', '.styled.'))),
            "Full (Hadronization)",
            line_color="black",
            line_style="dashed",
        )

    command: List[str] = [
        sys.executable,
        str(script_path("rivet_mkhtml_safe.py")),
        tool,
        f"{no_had_plot_input}:Title=Full (No Hadronization)",
        f"{had_plot_input}:Title=Full (Hadronization)",
        "--verbose",
        "-o",
        str(this_plot_dir),
    ]
    if dry_run:
        print(shlex.join(command))
    else:
        if this_plot_dir.exists():
            shutil.rmtree(this_plot_dir)
        this_plot_dir.mkdir(parents=True, exist_ok=True)
        env = os.environ.copy()
        if not scale_variations_enabled:
            env["DISPOL_FORCE_NO_ERROR_BANDS"] = "1"
        proc = subprocess.run(command, cwd=base_dir, text=True, capture_output=True, env=env)
        if proc.stdout:
            (this_plot_dir.parent / f"{this_plot_dir.name}.mkhtml.stdout").write_text(proc.stdout)
        if proc.stderr:
            (this_plot_dir.parent / f"{this_plot_dir.name}.mkhtml.stderr").write_text(proc.stderr)
        if proc.returncode != 0:
            raise RuntimeError(f"rivet-mkhtml failed for SPINCOMP/SPINHAD hadronization toggle with rc={proc.returncode}")
        if scale_variations_enabled:
            patched_count, rerendered_count = rewrite_scale_envelope_plot_scripts(this_plot_dir)
            print_stage(
                "Postprocessed scale-envelope Rivet scripts for PS hadronization toggle: "
                f"patched {patched_count}, rerendered {rerendered_count}"
            )
        else:
            patched_count, rerendered_count = rewrite_no_scale_ratio_plot_scripts(this_plot_dir)
            print_stage(
                "Postprocessed no-scale Rivet ratio scripts for PS hadronization toggle: "
                f"patched {patched_count}, rerendered {rerendered_count}"
            )

    outputs = {"full_toggle": str(this_plot_dir)}
    if outputs and not dry_run:
        update_manifest_outputs(campaign_dir, ps_hadronization_plot_outputs=outputs)
    return outputs


def run_ps_plotting_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    rivet_mkhtml_tool: str,
    plot_dir: Optional[Path] = None,
    scale_variations: Optional[bool] = None,
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    outputs: Dict[str, str] = {}
    if multi_family_ps_setups(setups):
        outputs.update(
            run_ps_rivetplot_campaign(
                base_dir=base_dir,
                tag=tag,
                setups=setups,
                analysis_name=analysis_name,
                rivet_mkhtml_tool=rivet_mkhtml_tool,
                plot_dir=plot_dir,
                scale_variations=scale_variations,
                max_workers=max_workers,
                dry_run=dry_run,
            )
        )
    if has_ps_hadronization_toggle_inputs(setups):
        outputs.update(
            run_ps_hadronization_rivetplot_campaign(
                base_dir=base_dir,
                tag=tag,
                setups=setups,
                analysis_name=analysis_name,
                rivet_mkhtml_tool=rivet_mkhtml_tool,
                scale_variations=scale_variations,
                max_workers=max_workers,
                dry_run=dry_run,
            )
        )
    if not outputs:
        raise ValueError(
            "PS plotting requires a multi-family PS setup (SPINCOMP/SPINHAD) and/or the combined SPINCOMP+SPINHAD hadronization inputs."
        )
    return outputs



def run_rivetplot_campaign(
    base_dir: Path,
    tag: str,
    setups: Sequence[str],
    analysis_name: str,
    order: str,
    rivet_mkhtml_tool: str,
    ratiolabel: str,
    reflabel_prefix: str,
    plot_dir: Optional[Path] = None,
    scale_variations: Optional[bool] = None,
    max_workers: int = 1,
    dry_run: bool = False,
) -> Dict[str, str]:
    campaign_dir = resolve_campaign_dir(base_dir, tag)
    manifest = load_manifest(campaign_dir)
    scale_variations_enabled = resolve_scale_variations_requested(scale_variations, manifest, default=False)
    analysis_outputs = manifest.get("analysis_yoda", {})
    if not isinstance(analysis_outputs, dict):
        analysis_outputs = {}
    analysis_plot_outputs = manifest.get("analysis_plot_yoda", {})
    if not isinstance(analysis_plot_outputs, dict):
        analysis_plot_outputs = {}
    raw_powheg_analysis_outputs = manifest.get("raw_powheg_analysis_yoda", {})
    if not isinstance(raw_powheg_analysis_outputs, dict):
        raw_powheg_analysis_outputs = {}
    raw_powheg_analysis_plot_outputs = manifest.get("raw_powheg_analysis_plot_yoda", {})
    if not isinstance(raw_powheg_analysis_plot_outputs, dict):
        raw_powheg_analysis_plot_outputs = {}
    raw_powheg_channel_analysis_outputs = manifest.get("raw_powheg_channel_analysis_yoda", {})
    if not isinstance(raw_powheg_channel_analysis_outputs, dict):
        raw_powheg_channel_analysis_outputs = {}
    raw_powheg_channel_analysis_plot_outputs = manifest.get("raw_powheg_channel_analysis_plot_yoda", {})
    if not isinstance(raw_powheg_channel_analysis_plot_outputs, dict):
        raw_powheg_channel_analysis_plot_outputs = {}
    reference_output = manifest.get("reference_yoda")
    if not isinstance(reference_output, str) or not reference_output:
        raise FileNotFoundError(f"Missing POLDIS reference YODA in {campaign_dir / 'manifest.json'}")

    tool = choose_rivet_mkhtml_tool(rivet_mkhtml_tool)
    if tool is None:
        raise FileNotFoundError(f"Could not find rivet-mkhtml tool '{rivet_mkhtml_tool}'")

    selected = normalize_setups(setups)
    if plot_dir is not None and len(selected) != 1:
        raise ValueError("--plot-dir may only be used with a single --setup")

    outputs: Dict[str, str] = {}
    manifest_analysis_variant = str(manifest.get("analysis_variant", ""))
    raw_powheg_title = raw_powheg_rivetplot_title()
    raw_powheg_channel_titles = {
        channel: raw_powheg_channel_rivetplot_title(channel)
        for channel in RAW_POWHEG_CHANNELS
    }
    worker_count = max(1, min(max_workers, len(selected))) if selected else 1
    mode = f"parallel, max_workers={worker_count}" if not dry_run and worker_count > 1 else "serial"
    print_stage(f"Building Rivet comparison plots for setups: {', '.join(selected)} ({mode})")

    def plot_setup(setup: str) -> tuple[str, str]:
        herwig_title = herwig_rivetplot_title(manifest_analysis_variant, setup)
        herwig_yoda = analysis_plot_outputs.get(setup) or analysis_outputs.get(setup)
        if not herwig_yoda:
            raise FileNotFoundError(
                f"Missing analyzed Herwig YODA for setup {setup} in {campaign_dir / 'manifest.json'}"
            )
        this_plot_dir = (plot_dir.resolve() if plot_dir is not None else default_plot_dir(base_dir, tag, setup, order))
        this_plot_dir.parent.mkdir(parents=True, exist_ok=True)
        sanitized_inputs_dir = campaign_dir / "analysis" / "_rivetplot_inputs"
        sanitized_inputs_dir.mkdir(parents=True, exist_ok=True)

        def sanitize_input_yoda(input_yoda: str, label: str, scale_envelope: bool = False) -> str:
            input_path = Path(input_yoda)
            output_name = input_path.name
            if output_name.endswith(".yoda.gz"):
                output_name = output_name[:-8] + f".{label}.sanitized.yoda.gz"
            elif output_name.endswith(".yoda"):
                output_name = output_name[:-5] + f".{label}.sanitized.yoda"
            else:
                output_name = output_name + f".{label}.sanitized.yoda"
            sanitized_path = sanitized_inputs_dir / output_name
            if not dry_run:
                if scale_envelope:
                    kept_count, dropped_count, dropped = build_scale_envelope_plot_yoda(input_path, sanitized_path)
                    print(
                        f"[stage] Built scale-envelope plot YODA for {label} in {setup}: kept {kept_count} objects, "
                        f"dropped {dropped_count} incompatible/empty objects -> {sanitized_path}",
                        flush=True,
                    )
                else:
                    kept_count, dropped_count, dropped = sanitize_plot_yoda(input_path, sanitized_path)
                    print(
                        f"[stage] Sanitized {label} plot YODA for {setup}: kept {kept_count} objects, "
                        f"dropped {dropped_count} incompatible/empty objects -> {sanitized_path}",
                        flush=True,
                    )
                if dropped_count:
                    for dropped_path, reason in sorted(dropped.items())[:20]:
                        print(f"[warn]   {dropped_path}: {reason}", flush=True)
                    if dropped_count > 20:
                        print(f"[warn]   ... and {dropped_count - 20} more", flush=True)
            return str(sanitized_path)

        herwig_plot_input = sanitize_input_yoda(
            herwig_yoda,
            f"herwig_{setup.lower()}",
            scale_envelope=scale_variations_enabled,
        )
        herwig_file_arg = f"{herwig_plot_input}:Title={herwig_title}"
        command = [tool, herwig_file_arg]
        raw_powheg_yoda = raw_powheg_analysis_plot_outputs.get(setup) or raw_powheg_analysis_outputs.get(setup)
        if raw_powheg_yoda:
            raw_powheg_plot_input = sanitize_input_yoda(raw_powheg_yoda, f"raw_powheg_{setup.lower()}")
            command.append(f"{raw_powheg_plot_input}:Title={raw_powheg_title}")
        for channel in RAW_POWHEG_CHANNELS:
            channel_plot_map = raw_powheg_channel_analysis_plot_outputs.get(channel, {})
            if not isinstance(channel_plot_map, dict):
                channel_plot_map = {}
            channel_yoda_map = raw_powheg_channel_analysis_outputs.get(channel, {})
            if not isinstance(channel_yoda_map, dict):
                channel_yoda_map = {}
            channel_yoda = channel_plot_map.get(setup) or channel_yoda_map.get(setup)
            if channel_yoda:
                channel_plot_input = sanitize_input_yoda(
                    channel_yoda,
                    f"raw_powheg_{channel.lower()}_{setup.lower()}",
                )
                command.append(f"{channel_plot_input}:Title={raw_powheg_channel_titles[channel]}")
        reference_plot_input = sanitize_input_yoda(reference_output, f"reference_{setup.lower()}")
        command: List[str] = [
            sys.executable,
            str(script_path("rivet_mkhtml_safe.py")),
            tool,
            *command[1:],
        ]
        command.extend(
            [
                reference_plot_input,
                f"--reflabel={reflabel_prefix} {order}",
                f"--ratiolabel={ratiolabel}",
                "--verbose",
                "-o",
                str(this_plot_dir),
            ]
        )
        if dry_run:
            print(shlex.join(command))
        else:
            if this_plot_dir.exists():
                shutil.rmtree(this_plot_dir)
            this_plot_dir.mkdir(parents=True, exist_ok=True)
            proc = subprocess.run(command, cwd=base_dir, text=True, capture_output=True)
            if proc.stdout:
                (this_plot_dir.parent / f"{this_plot_dir.name}.mkhtml.stdout").write_text(proc.stdout)
            if proc.stderr:
                (this_plot_dir.parent / f"{this_plot_dir.name}.mkhtml.stderr").write_text(proc.stderr)
            if proc.returncode != 0:
                raise RuntimeError(f"rivet-mkhtml failed for {setup} with rc={proc.returncode}")
            if scale_variations_enabled:
                patched_count, rerendered_count = rewrite_scale_envelope_plot_scripts(this_plot_dir)
                print_stage(
                    f"Postprocessed scale-envelope Rivet scripts for {setup}: "
                    f"patched {patched_count}, rerendered {rerendered_count}"
                )
        return setup, str(this_plot_dir)

    for setup, plot_output in run_parallel_ordered(
        selected,
        plot_setup,
        1 if dry_run else max_workers,
    ):
        outputs[setup] = plot_output

    if outputs and not dry_run:
        update_manifest_outputs(campaign_dir, plot_outputs=outputs)
    return outputs


def collect_existing_yoda_results(
    base_dir: Path,
    base_tag: str,
    analysis_variant: str,
    setups: Sequence[str] = (),
    analysis_variants_by_setup: Optional[Dict[str, Sequence[str]]] = None,
    scale_variations: bool = False,
    include_lo: Optional[bool] = None,
) -> List[JobResult]:
    results: List[JobResult] = []
    tag_re = re.compile(rf"^(?P<stem>.+)-S(?P<seed>\d+)-(?P<tag>{re.escape(base_tag)}(?:-s[^.]+)?)")
    for job in build_jobs(
        0,
        0,
        0,
        analysis_variant,
        analysis_variants_by_setup=analysis_variants_by_setup,
        setups=setups,
        scale_variations=scale_variations,
        campaign_tag=base_tag,
        include_lo=include_lo,
    ):
        for path in sorted(base_dir.glob(f"{job.stem}-S*-{base_tag}*.yoda*")):
            match = tag_re.match(path.name)
            shard_tag = match.group("tag") if match else base_tag
            seed = int(match.group("seed")) if match else 0
            spec = ShardSpec(job=job, shard_index=1, shard_count=1, tag=shard_tag, seed=seed, events=0)
            results.append(
                JobResult(
                    spec=spec,
                    command=[],
                    returncode=0,
                    yoda_files=[str(path.resolve())],
                )
            )
    return results


def collect_existing_job_artifacts_from_manifest(campaign_dir: Path) -> List[JobResult]:
    manifest = load_manifest(campaign_dir)
    results: List[JobResult] = []
    seen_slots: set[tuple[str, str, int]] = set()
    for key in ("finished", "failed"):
        entries = manifest.get(key, [])
        if not isinstance(entries, list):
            continue
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            spec_data = entry.get("spec", {})
            if not isinstance(spec_data, dict):
                continue
            job_data = spec_data.get("job", {})
            if not isinstance(job_data, dict):
                continue
            job = JobSpec(
                setup=str(job_data.get("setup", "")),
                order=str(job_data.get("order", "")),
                helicity=str(job_data.get("helicity", "")),
                stem=str(job_data.get("stem", "")),
                run_file=str(job_data.get("run_file", "")),
                in_file=str(job_data.get("in_file", "")),
                events=int(job_data.get("events", 0)),
                analysis_variant=str(job_data.get("analysis_variant", "")),
                scale_variation=str(job_data.get("scale_variation", "nominal")),
                scale_factor=float(job_data.get("scale_factor", 1.0)),
                source_in_file=str(job_data.get("source_in_file", "")),
            )
            spec = ShardSpec(
                job=job,
                shard_index=int(spec_data.get("shard_index", 1)),
                shard_count=int(spec_data.get("shard_count", 1)),
                tag=str(spec_data.get("tag", "")),
                seed=int(spec_data.get("seed", 0)),
                events=int(spec_data.get("events", 0)),
                rerun_parent_tag=str(spec_data.get("rerun_parent_tag", "")),
            )
            slot = (spec.job.stem, spec.tag, spec.seed)
            if slot in seen_slots:
                continue
            seen_slots.add(slot)
            results.append(
                JobResult(
                    spec=spec,
                    command=entry.get("command", []) if isinstance(entry.get("command"), list) else [],
                    returncode=int(entry.get("returncode", 0)) if entry.get("returncode") is not None else None,
                    started_at=entry.get("started_at"),
                    ended_at=entry.get("ended_at"),
                    launcher_log=entry.get("launcher_log"),
                    out_files=[str(path) for path in entry.get("out_files", [])] if isinstance(entry.get("out_files"), list) else [],
                    log_files=[str(path) for path in entry.get("log_files", [])] if isinstance(entry.get("log_files"), list) else [],
                    yoda_files=[str(path) for path in entry.get("yoda_files", [])] if isinstance(entry.get("yoda_files"), list) else [],
                )
            )
    return results


def collect_existing_prepared_from_manifest(campaign_dir: Path) -> List[JobResult]:
    manifest = load_manifest(campaign_dir)
    entries = manifest.get("prepared", [])
    results: List[JobResult] = []
    if not isinstance(entries, list):
        return results
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        spec_data = entry.get("spec", {})
        if not isinstance(spec_data, dict):
            continue
        job_data = spec_data.get("job", {})
        if not isinstance(job_data, dict):
            continue
        job = JobSpec(
            setup=str(job_data.get("setup", "")),
            order=str(job_data.get("order", "")),
            helicity=str(job_data.get("helicity", "")),
            stem=str(job_data.get("stem", "")),
            run_file=str(job_data.get("run_file", "")),
            in_file=str(job_data.get("in_file", "")),
            events=int(job_data.get("events", 0)),
            analysis_variant=str(job_data.get("analysis_variant", "")),
            scale_variation=str(job_data.get("scale_variation", "nominal")),
            scale_factor=float(job_data.get("scale_factor", 1.0)),
            source_in_file=str(job_data.get("source_in_file", "")),
        )
        spec = ShardSpec(
            job=job,
            shard_index=int(spec_data.get("shard_index", 1)),
            shard_count=int(spec_data.get("shard_count", 1)),
            tag=str(spec_data.get("tag", "")),
            seed=int(spec_data.get("seed", 0)),
            events=int(spec_data.get("events", 0)),
            rerun_parent_tag=str(spec_data.get("rerun_parent_tag", "")),
        )
        results.append(
            JobResult(
                spec=spec,
                command=[],
                prepared=True,
                prepare_command=entry.get("prepare_command") if isinstance(entry.get("prepare_command"), list) else None,
                prepare_returncode=int(entry.get("prepare_returncode", 0)) if entry.get("prepare_returncode") is not None else None,
            )
        )
    return results


def shard_slot_key(spec: ShardSpec) -> tuple[str, str]:
    return (spec.job.stem, spec.tag)


def rerun_parent_slot_key(spec: ShardSpec) -> Optional[tuple[str, str]]:
    if not spec.rerun_parent_tag:
        return None
    return (spec.job.stem, spec.rerun_parent_tag)


def unresolved_failed_results(results: Sequence[JobResult]) -> List[JobResult]:
    latest_by_slot: Dict[tuple[str, str], JobResult] = {}
    successful_slots = set()
    replaced_slots = {
        parent_slot
        for item in results
        for parent_slot in [rerun_parent_slot_key(item.spec)]
        if parent_slot is not None
    }
    for item in results:
        slot = shard_slot_key(item.spec)
        latest_by_slot[slot] = item
        if item.returncode == 0:
            successful_slots.add(slot)
    unresolved = [
        item
        for slot, item in latest_by_slot.items()
        if slot not in replaced_slots and slot not in successful_slots and item.returncode not in (0, None)
    ]
    unresolved.sort(key=lambda item: (item.spec.job.stem, item.spec.tag, item.spec.seed))
    return unresolved


def build_failed_rerun_shards_with_random_seeds(
    existing_results: Sequence[JobResult],
    base_tag: str,
    rerun_shards: int = 1,
) -> List[ShardSpec]:
    failed = unresolved_failed_results(existing_results)
    if not failed:
        return []

    if rerun_shards < 1:
        raise ValueError("--rerun-failed-shards must be at least 1.")

    used_seeds = {item.spec.seed for item in existing_results if item.spec.seed > 0}
    used_tags_by_stem: Dict[str, set[str]] = {}
    for item in existing_results:
        used_tags_by_stem.setdefault(item.spec.job.stem, set()).add(item.spec.tag)

    rerun_specs: List[ShardSpec] = []
    for item in failed:
        stem_tags = used_tags_by_stem.setdefault(item.spec.job.stem, set())
        event_parts = split_events(item.spec.events, rerun_shards)
        shard_count = len(event_parts)
        for child_index, events in enumerate(event_parts, start=1):
            new_seed = 0
            while new_seed <= 0 or new_seed in used_seeds:
                new_seed = secrets.randbelow(2_000_000_000)
            used_seeds.add(new_seed)

            if shard_count == 1:
                new_tag = item.spec.tag
                rerun_parent_tag = ""
            else:
                if item.spec.tag == base_tag:
                    tag_prefix = f"{base_tag}-sr"
                else:
                    tag_prefix = f"{item.spec.tag}r"
                suffix_index = child_index
                new_tag = f"{tag_prefix}{suffix_index:03d}"
                while new_tag in stem_tags:
                    suffix_index += 1
                    new_tag = f"{tag_prefix}{suffix_index:03d}"
                rerun_parent_tag = item.spec.tag

            stem_tags.add(new_tag)
            rerun_specs.append(
                replace(
                    item.spec,
                    shard_index=child_index,
                    shard_count=shard_count,
                    tag=new_tag,
                    seed=new_seed,
                    events=events,
                    rerun_parent_tag=rerun_parent_tag,
                )
            )
    return rerun_specs


def update_manifest_postprocess(
    campaign_dir: Path,
    args: argparse.Namespace,
    extract_status: Dict[str, int],
    yoda_merge_results: Sequence[YODAMergeResult],
    yoda_nlo_results: Sequence[YODANLOResult],
) -> None:
    manifest = load_manifest(campaign_dir)
    resolved_analysis_variant = resolved_analysis_variant_from_args(args)
    include_lo = resolve_include_lo_requested(getattr(args, "include_lo", None), manifest, resolved_analysis_variant)
    analysis_variants = manifest.get("analysis_variants", [])
    if not isinstance(analysis_variants, list):
        analysis_variants = []
    manifest["tag"] = args.tag
    manifest["base_dir"] = str(args.base_dir.resolve())
    manifest["merge_yoda"] = not args.no_merge_yoda
    manifest["raw_powheg"] = bool(getattr(args, "raw_powheg", False))
    manifest["include_lo"] = include_lo
    manifest["extract_diagnostics"] = bool(getattr(args, "diagnostics", False))
    manifest["scale_variations"] = bool(getattr(args, "scale_variations", False))
    manifest["scale_variation_factors"] = dict(SCALE_VARIATION_FACTORS)
    manifest["yoda_merge_tool"] = args.yoda_merge_tool
    manifest["rivet"] = bool(resolved_analysis_variant)
    manifest["analysis_variant"] = resolved_analysis_variant
    manifest["analysis_variants"] = [str(item) for item in analysis_variants]
    manifest["ps_mode"] = bool(analysis_variants and all(str(item).startswith("RIVETPS-") for item in analysis_variants))
    manifest["extract_status"] = extract_status
    manifest["merged_yoda"] = [
        {
            "logical_run": item.logical_run,
            "source_files": item.source_files,
            "output_file": item.output_file,
            "command": item.command,
            "returncode": item.returncode,
            "copied": item.copied,
            "skipped": item.skipped,
            "message": item.message,
        }
        for item in yoda_merge_results
    ]
    manifest["nlo_yoda"] = [
        {
            "logical_run": item.logical_run,
            "pos_file": item.pos_file,
            "neg_file": item.neg_file,
            "output_file": item.output_file,
            "command": item.command,
            "returncode": item.returncode,
            "skipped": item.skipped,
            "message": item.message,
        }
        for item in yoda_nlo_results
    ]
    manifest["all_merged_yoda_files"] = sorted({item.output_file for item in yoda_merge_results if item.output_file is not None})
    manifest["all_nlo_yoda_files"] = sorted({item.output_file for item in yoda_nlo_results if item.output_file is not None})
    save_manifest(campaign_dir, manifest)


def run_campaign_command(args: argparse.Namespace) -> int:
    base_dir = args.base_dir.resolve()
    campaign_dir = resolve_campaign_dir(base_dir, args.tag)
    campaign_dir.mkdir(parents=True, exist_ok=True)
    analysis_variant = analysis_variant_from_args(args)
    enable_raw_powheg = bool(getattr(args, "raw_powheg", False))
    rerun_failed = bool(getattr(args, "rerun_failed_random_seed", False))
    rerun_failed_shards = int(getattr(args, "rerun_failed_shards", 1))
    manifest: Dict[str, object] = {}

    if rerun_failed_shards < 1:
        raise ValueError("--rerun-failed-shards must be at least 1.")
    if rerun_failed_shards != 1 and not rerun_failed:
        raise ValueError("--rerun-failed-shards is only meaningful together with --rerun-failed-random-seed.")

    existing_prepared: List[JobResult] = []
    existing_finished: List[JobResult] = []
    if rerun_failed:
        manifest = load_manifest(campaign_dir)
        if not manifest:
            raise FileNotFoundError(
                f"--rerun-failed-random-seed requested, but no existing manifest was found at {campaign_dir / 'manifest.json'}"
            )
        manifest_variant = str(manifest.get("analysis_variant", ""))
        if analysis_variant and manifest_variant and analysis_variant != manifest_variant:
            raise ValueError(
                f"Requested analysis variant {analysis_variant!r} does not match existing campaign variant {manifest_variant!r}"
            )
        analysis_variant = analysis_variant or manifest_variant
        enable_raw_powheg = bool(enable_raw_powheg or manifest.get("raw_powheg"))
    scale_variations = resolve_scale_variations_requested(getattr(args, "scale_variations", None), manifest, default=False)
    include_lo = resolve_include_lo_requested(getattr(args, "include_lo", None), manifest, analysis_variant)
    diagnostics_enabled = resolve_extract_diagnostics_requested(getattr(args, "diagnostics", None), manifest if rerun_failed else None)
    requested_setups = resolve_campaign_setups_requested(getattr(args, "setup", ()), manifest if rerun_failed else None)
    validate_setup_family(requested_setups)
    ps_mode = contains_ps_setups(requested_setups)
    analysis_variants_by_setup = ps_analysis_variants_by_setup(requested_setups) if ps_mode else None
    if ps_mode:
        if analysis_variant == "RIVETFO":
            raise ValueError("SPINVAL/SPINCOMP/SPINHAD are not supported with --rivetfo.")
        if enable_raw_powheg:
            raise ValueError("--raw-powheg is not supported for SPINVAL/SPINCOMP/SPINHAD.")
        if getattr(args, "include_lo", None):
            raise ValueError("SPINVAL/SPINCOMP/SPINHAD are NLO-only workflows; do not use --include-lo.")
        analysis_variant = "RIVETPS"
        include_lo = False
        diagnostics_enabled = False
    resolve_analysis_name(args, requested_setups)
    args.setup = requested_setups
    args.scale_variations = scale_variations
    args.include_lo = include_lo
    args.diagnostics = diagnostics_enabled
    args.resolved_analysis_variant = analysis_variant
    if rerun_failed:
        existing_prepared = collect_existing_prepared_from_manifest(campaign_dir)
        existing_finished = collect_existing_job_artifacts_from_manifest(campaign_dir)
        rerun_shards = build_failed_rerun_shards_with_random_seeds(
            existing_finished,
            args.tag,
            rerun_failed_shards,
        )
        expected_jobs = build_jobs(
            args.lo_events,
            args.posnlo_events,
            args.negnlo_events,
            analysis_variant,
            analysis_variants_by_setup=analysis_variants_by_setup,
            setups=requested_setups,
            scale_variations=scale_variations,
            campaign_tag=args.tag,
            include_lo=include_lo,
        )
        expected_shards = build_shards(
            expected_jobs,
            args.tag,
            args.shards,
            max(1, args.jobs),
            args.seed_base,
        )
        existing_slots = {shard_slot_key(item.spec) for item in existing_finished}
        missing_shards = [
            spec
            for spec in expected_shards
            if shard_slot_key(spec) not in existing_slots
        ]
        missing_jobs = sorted({spec.job.stem for spec in missing_shards})
        shards = [*missing_shards, *rerun_shards]
        jobs = sorted({item.job for item in shards}, key=lambda job: (job.setup, job.order, job.helicity, job.stem))
        unresolved_before = unresolved_failed_results(existing_finished)
        print_stage(
            f"Rerunning failed shards for campaign '{args.tag}' "
            f"({len(unresolved_before)} unresolved, {len(missing_jobs)} logical runs with missing shards, "
            f"{len(shards)} total replacement/new shards, rerun_shards={rerun_failed_shards}, "
            f"variant={analysis_variant or 'default'})"
        )
        if not shards:
            print_stage("No unresolved failed shards or missing logical runs remain; reusing existing campaign outputs")
    else:
        jobs = build_jobs(
            args.lo_events,
            args.posnlo_events,
            args.negnlo_events,
            analysis_variant,
            analysis_variants_by_setup=analysis_variants_by_setup,
            setups=requested_setups,
            scale_variations=scale_variations,
            campaign_tag=args.tag,
            include_lo=include_lo,
        )
        shards = build_shards(jobs, args.tag, args.shards, max(1, args.jobs), args.seed_base)
        print_stage(
            f"Preparing campaign '{args.tag}' ({len(jobs)} logical runs, {len(shards)} shards, variant={analysis_variant or 'default'})"
        )

    if enable_raw_powheg and analysis_variant != "RIVETFO":
        raise ValueError("--raw-powheg is only supported together with --rivetfo.")

    prepared: List[JobResult] = list(existing_prepared)
    do_prepare = args.force_prepare or not args.no_prepare
    for spec in jobs:
        prep = ensure_run_file(
            base_dir,
            spec,
            do_prepare,
            args.dry_run,
            analysis_variant,
            force_prepare=args.force_prepare,
            enable_raw_powheg=enable_raw_powheg,
        )
        if prep is not None:
            prepared.append(prep)

    powheg_options_file = None
    if not args.dry_run:
        print_stage("Recording POWHEG-related card options for this campaign")
        powheg_options_file = write_powheg_options_summary(campaign_dir, base_dir, jobs, analysis_variant)

    finished = run_jobs(
        base_dir,
        shards,
        max(1, args.jobs),
        campaign_dir,
        args.dry_run,
        args.keep_going,
        args.progress_interval,
        args.max_listed,
        analysis_variant,
    )

    finished = list(existing_finished) + finished
    failed = unresolved_failed_results(finished)
    extract_status: Dict[str, int] = {}
    yoda_merge_results: List[YODAMergeResult] = []
    yoda_nlo_results: List[YODANLOResult] = []
    raw_powheg_results: List[RawPOWHEGYODAResult] = []
    if not failed:
        if use_raw_powheg_runtime(analysis_variant, enable_raw_powheg):
            print_stage("Converting raw POWHEG log diagnostics to shard YODAs")
            raw_powheg_results = build_raw_powheg_yoda_files(
                base_dir,
                campaign_dir,
                finished,
                analysis_variant,
                max_workers=max(1, args.jobs),
                dry_run=args.dry_run,
            )
            raw_failures = [item for item in raw_powheg_results if item.returncode != 0 and not item.skipped]
            extract_status["raw_powheg"] = 0 if not raw_failures else 1
        else:
            extract_status["raw_powheg"] = -2
        if not args.no_merge_yoda:
            yoda_merge_results = merge_yoda_files(
                finished,
                base_dir,
                campaign_dir,
                args.tag,
                args.yoda_merge_tool,
                max(1, args.jobs),
                args.dry_run,
            )
            merge_failures = [item for item in yoda_merge_results if item.returncode != 0]
            if use_rivet_runtime(analysis_variant):
                nlo_skip_reason = "Skipped generic NLO YODA build for Rivet workflow; analyze-herwig subtracts POSNLO and NEGNLO directly."
                yoda_nlo_results = skipped_nlo_yoda_results(
                    yoda_merge_results,
                    nlo_skip_reason,
                )
                nlo_failures = []
            else:
                yoda_nlo_results = build_nlo_yoda_files(
                    yoda_merge_results,
                    base_dir,
                    campaign_dir,
                    args.tag,
                    args.yoda_merge_tool,
                    max(1, args.jobs),
                    args.dry_run,
                )
                nlo_failures = [item for item in yoda_nlo_results if item.returncode != 0 and not item.skipped]
            extract_status["yoda_merge"] = 0 if not merge_failures else 1
            extract_status["yoda_nlo"] = 0 if not nlo_failures else 1
        else:
            extract_status["yoda_merge"] = -2
            extract_status["yoda_nlo"] = -2
        if ps_mode:
            extract_status["results"] = -2
            extract_status["diagnostic"] = -2
            if "SPINVAL" in requested_setups:
                extract_status["spin_diagnostic"] = run_spin_diagnostic_extractor(base_dir, args.tag, args.dry_run)
            else:
                extract_status["spin_diagnostic"] = -2
        else:
            dis_cmd = [
                "python3.10",
                str(script_path("extract_dis_out_results.py")),
                "--base-dir",
                str(base_dir),
                "-t",
                args.tag,
                "--json-out",
                str(campaign_dir / "results.json"),
                "--csv-out",
                str(campaign_dir / "results.csv"),
            ]
            extend_results_extractor_scope(dis_cmd, requested_setups)
            dis_cmd.extend(analysis_variant_args(analysis_variant))
            print_stage("Extracting DIS run summaries")
            extract_status["results"] = run_extractor(dis_cmd, campaign_dir / "results.txt", base_dir, args.dry_run)
            if diagnostics_enabled:
                diag_cmd = [
                    "python3.10",
                    str(script_path("extract_nlo_term_diagnostics.py")),
                    "--base-dir",
                    str(base_dir),
                    "-t",
                    args.tag,
                    "--strict-tag",
                    "--json-out",
                    str(campaign_dir / "diagnostic.json"),
                    "--csv-out",
                    str(campaign_dir / "diagnostic.csv"),
                ]
                diag_cmd.extend(analysis_variant_args(analysis_variant))
                print_stage("Extracting NLO-term diagnostics")
                extract_status["diagnostic"] = run_extractor(diag_cmd, campaign_dir / "diagnostic.txt", base_dir, args.dry_run)
            else:
                extract_status["diagnostic"] = -2
            extract_status["spin_diagnostic"] = -2
    else:
        extract_status["raw_powheg"] = -1
        extract_status["yoda_merge"] = -1
        extract_status["yoda_nlo"] = -1
        extract_status["results"] = -1
        extract_status["diagnostic"] = -1
        extract_status["spin_diagnostic"] = -1

    write_manifest(
        campaign_dir,
        args,
        prepared,
        finished,
        extract_status,
        yoda_merge_results,
        yoda_nlo_results,
        powheg_options_file=powheg_options_file,
    )
    if raw_powheg_results and not args.dry_run:
        update_manifest_outputs(
            campaign_dir,
            raw_powheg_yoda_outputs=[
                {
                    "logical_run": item.logical_run,
                    "channel": item.channel,
                    "variant": item.variant,
                    "source_log": item.source_log,
                    "source_out": item.source_out,
                    "output_file": item.output_file,
                    "summary_csv": item.summary_csv,
                    "command": item.command,
                    "returncode": item.returncode,
                    "skipped": item.skipped,
                    "message": item.message,
                }
                for item in raw_powheg_results
            ],
        )

    print(f"Campaign tag: {args.tag}")
    print(f"Campaign directory: {campaign_dir}")
    print(f"Prepared run files: {len(prepared)}")
    print(f"Logical runs: {len(jobs)}")
    print(f"Shards launched: {len(shards)}")
    print(f"Completed shards: {sum(1 for item in finished if item.returncode == 0)} / {len(shards)}")
    print(f"Failed shards: {len(failed)}")
    if failed:
        print("Failed jobs:")
        for item in failed:
            print(f"  {item.spec.job.stem} [{item.spec.tag}]: rc={item.returncode} log={item.launcher_log}")
        return 1
    print("Extraction outputs:")
    if extract_status.get("results") not in (-2, -1):
        print(f"  {campaign_dir / 'results.txt'}")
        print(f"  {campaign_dir / 'results.csv'}")
        print(f"  {campaign_dir / 'results.json'}")
    if diagnostics_enabled:
        print(f"  {campaign_dir / 'diagnostic.txt'}")
        print(f"  {campaign_dir / 'diagnostic.csv'}")
        print(f"  {campaign_dir / 'diagnostic.json'}")
    if extract_status.get("spin_diagnostic") not in (-2, -1):
        print(f"  {campaign_dir / 'spin_diagnostic.txt'}")
        print(f"  {campaign_dir / 'spin_diagnostic.csv'}")
        print(f"  {campaign_dir / 'spin_diagnostic.json'}")
    if yoda_merge_results:
        print("Merged YODA outputs:")
        for item in yoda_merge_results:
            if item.output_file is not None:
                print(f"  {item.output_file}")
    if yoda_nlo_results:
        print("Physical NLO YODA outputs:")
        for item in yoda_nlo_results:
            if item.output_file is not None:
                print(f"  {item.output_file}")
    if raw_powheg_results:
        print("Raw POWHEG shard YODA outputs:")
        for item in raw_powheg_results:
            if item.output_file is not None and item.returncode == 0:
                print(f"  {item.output_file}")
    print(f"  {campaign_dir / 'manifest.json'}")
    return 0


def run_postprocess_command(args: argparse.Namespace) -> int:
    base_dir = args.base_dir.resolve()
    campaign_dir = resolve_campaign_dir(base_dir, args.tag)
    campaign_dir.mkdir(parents=True, exist_ok=True)
    manifest = load_manifest(campaign_dir)
    analysis_variant = analysis_variant_from_args(args) or str(manifest.get("analysis_variant", ""))
    requested_setups = resolve_campaign_setups_requested(args.setup, manifest)
    validate_setup_family(requested_setups)
    ps_mode = contains_ps_setups(requested_setups)
    if ps_mode:
        if analysis_variant == "RIVETFO":
            raise ValueError("SPINVAL/SPINCOMP/SPINHAD are not supported with --rivetfo.")
        analysis_variant = "RIVETPS"
    scale_variations = resolve_scale_variations_requested(getattr(args, "scale_variations", None), manifest, default=False)
    include_lo = resolve_include_lo_requested(getattr(args, "include_lo", None), manifest, analysis_variant)
    diagnostics_enabled = resolve_extract_diagnostics_requested(getattr(args, "diagnostics", None), manifest)
    if ps_mode:
        if include_lo:
            raise ValueError("SPINVAL/SPINCOMP/SPINHAD are NLO-only workflows; do not use --include-lo.")
        diagnostics_enabled = False
    resolve_analysis_name(args, requested_setups)
    args.setup = requested_setups
    args.scale_variations = scale_variations
    args.include_lo = include_lo
    args.diagnostics = diagnostics_enabled
    args.resolved_analysis_variant = analysis_variant
    if getattr(args, "raw_powheg", False) and analysis_variant != "RIVETFO":
        raise ValueError("--raw-powheg is only supported together with --rivetfo.")
    print_stage(f"Postprocessing existing campaign '{args.tag}'")

    finished = collect_existing_yoda_results(
        base_dir,
        args.tag,
        analysis_variant,
        setups=requested_setups,
        analysis_variants_by_setup=ps_analysis_variants_by_setup(requested_setups) if ps_mode else None,
        scale_variations=scale_variations,
        include_lo=include_lo,
    )
    if not finished:
        print(f"No shard YODA files found for tag {args.tag} in {base_dir}")
        return 1

    extract_status: Dict[str, int] = {}
    yoda_merge_results: List[YODAMergeResult] = []
    yoda_nlo_results: List[YODANLOResult] = []
    raw_powheg_results: List[RawPOWHEGYODAResult] = []

    if use_raw_powheg_runtime(analysis_variant, bool(getattr(args, "raw_powheg", False))):
        print_stage("Rebuilding raw POWHEG shard YODAs from existing logs")
        raw_sources = collect_existing_job_artifacts_from_manifest(campaign_dir)
        raw_powheg_results = build_raw_powheg_yoda_files(
            base_dir,
            campaign_dir,
            raw_sources,
            analysis_variant,
            max_workers=max(1, args.jobs),
            dry_run=args.dry_run,
        )
        raw_failures = [item for item in raw_powheg_results if item.returncode != 0 and not item.skipped]
        extract_status["raw_powheg"] = 0 if not raw_failures else 1
    else:
        extract_status["raw_powheg"] = -2

    if not args.no_merge_yoda:
        yoda_merge_results = merge_yoda_files(
            finished,
            base_dir,
            campaign_dir,
            args.tag,
            args.yoda_merge_tool,
            max(1, args.jobs),
            args.dry_run,
        )
        merge_failures = [item for item in yoda_merge_results if item.returncode != 0]
        if use_rivet_runtime(analysis_variant):
            nlo_skip_reason = "Skipped generic NLO YODA build for Rivet workflow; analyze-herwig subtracts POSNLO and NEGNLO directly."
            yoda_nlo_results = skipped_nlo_yoda_results(
                yoda_merge_results,
                nlo_skip_reason,
            )
            nlo_failures = []
        else:
            yoda_nlo_results = build_nlo_yoda_files(
                yoda_merge_results,
                base_dir,
                campaign_dir,
                args.tag,
                args.yoda_merge_tool,
                max(1, args.jobs),
                args.dry_run,
            )
            nlo_failures = [item for item in yoda_nlo_results if item.returncode != 0 and not item.skipped]
        extract_status["yoda_merge"] = 0 if not merge_failures else 1
        extract_status["yoda_nlo"] = 0 if not nlo_failures else 1
    else:
        extract_status["yoda_merge"] = -2
        extract_status["yoda_nlo"] = -2

    if ps_mode:
        extract_status["results"] = -2
        extract_status["diagnostic"] = -2
        if "SPINVAL" in requested_setups:
            extract_status["spin_diagnostic"] = run_spin_diagnostic_extractor(base_dir, args.tag, args.dry_run)
        else:
            extract_status["spin_diagnostic"] = -2
    else:
        dis_cmd = [
            "python3.10",
            str(script_path("extract_dis_out_results.py")),
            "--base-dir",
            str(base_dir),
            "-t",
            args.tag,
            "--json-out",
            str(campaign_dir / "results.json"),
            "--csv-out",
            str(campaign_dir / "results.csv"),
        ]
        extend_results_extractor_scope(dis_cmd, requested_setups)
        dis_cmd.extend(analysis_variant_args(analysis_variant))
        extract_status["results"] = run_extractor(dis_cmd, campaign_dir / "results.txt", base_dir, args.dry_run)
        if diagnostics_enabled:
            diag_cmd = [
                "python3.10",
                str(script_path("extract_nlo_term_diagnostics.py")),
                "--base-dir",
                str(base_dir),
                "-t",
                args.tag,
                "--strict-tag",
                "--json-out",
                str(campaign_dir / "diagnostic.json"),
                "--csv-out",
                str(campaign_dir / "diagnostic.csv"),
            ]
            diag_cmd.extend(analysis_variant_args(analysis_variant))
            print_stage("Extracting NLO-term diagnostics")
            extract_status["diagnostic"] = run_extractor(diag_cmd, campaign_dir / "diagnostic.txt", base_dir, args.dry_run)
        else:
            extract_status["diagnostic"] = -2
        extract_status["spin_diagnostic"] = -2

    if ps_mode:
        analyze_ps_herwig_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            scale_variations=scale_variations,
            merge_tool_preference=args.yoda_merge_tool,
            max_workers=max(1, args.jobs),
            dry_run=args.dry_run,
        )
    elif args.dry_run:
        analyze_herwig_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            analysis_variant=analysis_variant,
            scale_variations=scale_variations,
            merge_tool_preference=args.yoda_merge_tool,
            max_workers=max(1, args.jobs),
            dry_run=True,
        )
        if use_raw_powheg_runtime(analysis_variant, bool(getattr(args, "raw_powheg", False))):
            analyze_raw_powheg_campaign(
                base_dir=base_dir,
                tag=args.tag,
                setups=args.setup,
                analysis_name=args.analysis_name,
                analysis_variant=analysis_variant,
                max_workers=max(1, args.jobs),
                dry_run=True,
            )
    else:
        analyze_herwig_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            analysis_variant=analysis_variant,
            scale_variations=scale_variations,
            merge_tool_preference=args.yoda_merge_tool,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
        if use_raw_powheg_runtime(analysis_variant, bool(getattr(args, "raw_powheg", False))):
            analyze_raw_powheg_campaign(
                base_dir=base_dir,
                tag=args.tag,
                setups=args.setup,
                analysis_name=args.analysis_name,
                analysis_variant=analysis_variant,
                max_workers=max(1, args.jobs),
                dry_run=False,
            )

    update_manifest_postprocess(campaign_dir, args, extract_status, yoda_merge_results, yoda_nlo_results)
    if raw_powheg_results and not args.dry_run:
        update_manifest_outputs(
            campaign_dir,
            raw_powheg_yoda_outputs=[
                {
                    "logical_run": item.logical_run,
                    "channel": item.channel,
                    "variant": item.variant,
                    "source_log": item.source_log,
                    "source_out": item.source_out,
                    "output_file": item.output_file,
                    "summary_csv": item.summary_csv,
                    "command": item.command,
                    "returncode": item.returncode,
                    "skipped": item.skipped,
                    "message": item.message,
                }
                for item in raw_powheg_results
            ],
        )

    print(f"Postprocessed campaign tag: {args.tag}")
    print(f"Campaign directory: {campaign_dir}")
    print(f"Found shard YODA files: {len(finished)}")
    print("Outputs:")
    if extract_status.get("results") not in (-2, -1):
        print(f"  {campaign_dir / 'results.txt'}")
        print(f"  {campaign_dir / 'results.csv'}")
        print(f"  {campaign_dir / 'results.json'}")
    if diagnostics_enabled:
        print(f"  {campaign_dir / 'diagnostic.txt'}")
        print(f"  {campaign_dir / 'diagnostic.csv'}")
        print(f"  {campaign_dir / 'diagnostic.json'}")
    if extract_status.get("spin_diagnostic") not in (-2, -1):
        print(f"  {campaign_dir / 'spin_diagnostic.txt'}")
        print(f"  {campaign_dir / 'spin_diagnostic.csv'}")
        print(f"  {campaign_dir / 'spin_diagnostic.json'}")
    if yoda_merge_results:
        print("Merged YODA outputs:")
        for item in yoda_merge_results:
            if item.output_file is not None:
                print(f"  {item.output_file}")
    if yoda_nlo_results:
        print("Physical NLO YODA outputs:")
        for item in yoda_nlo_results:
            if item.output_file is not None:
                print(f"  {item.output_file}")
    if raw_powheg_results:
        print("Raw POWHEG shard YODA outputs:")
        for item in raw_powheg_results:
            if item.output_file is not None and item.returncode == 0:
                print(f"  {item.output_file}")
    return 0


def run_analyze_herwig_command(args: argparse.Namespace) -> int:
    base_dir = args.base_dir.resolve()
    campaign_dir = resolve_campaign_dir(base_dir, args.tag)
    manifest = load_manifest(campaign_dir)
    requested_setups = resolve_campaign_setups_requested(args.setup, manifest)
    validate_setup_family(requested_setups)
    ps_mode = contains_ps_setups(requested_setups)
    analysis_variant = analysis_variant_from_args(args) or str(manifest.get("analysis_variant", ""))
    if ps_mode:
        analysis_variant = "RIVETPS"
    scale_variations = resolve_scale_variations_requested(getattr(args, "scale_variations", None), manifest, default=False)
    resolve_analysis_name(args, requested_setups)
    args.setup = requested_setups
    args.scale_variations = scale_variations
    args.resolved_analysis_variant = analysis_variant
    if getattr(args, "raw_powheg", False) and analysis_variant != "RIVETFO":
        raise ValueError("--raw-powheg is only supported together with --rivetfo.")
    include_raw_powheg = use_raw_powheg_runtime(
        analysis_variant,
        bool(getattr(args, "raw_powheg", False) or manifest.get("raw_powheg")),
    )
    print_stage(f"Running analyze-herwig for campaign '{args.tag}'")
    if ps_mode:
        outputs_by_family = analyze_ps_herwig_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            scale_variations=scale_variations,
            merge_tool_preference=args.yoda_merge_tool,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
        print("Herwig analyzed YODA outputs:")
        for family in sorted(outputs_by_family):
            for setup in normalize_campaign_setups(args.setup):
                if setup in outputs_by_family[family]:
                    print(f"  {setup}/{family}: {outputs_by_family[family][setup]}")
        return 0

    outputs = analyze_herwig_campaign(
        base_dir=base_dir,
        tag=args.tag,
        setups=args.setup,
        analysis_name=args.analysis_name,
        analysis_variant=analysis_variant,
        scale_variations=scale_variations,
        merge_tool_preference=args.yoda_merge_tool,
        max_workers=max(1, args.jobs),
        dry_run=False,
    )
    if include_raw_powheg:
        analyze_raw_powheg_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            analysis_variant=analysis_variant,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
    print("Herwig analyzed YODA outputs:")
    for setup in normalize_setups(args.setup):
        print(f"  {setup}: {outputs[setup]}")
    return 0


def run_poldis_top_command(args: argparse.Namespace, campaign_tag: Optional[str] = None) -> int:
    base_dir = args.base_dir.resolve()
    label = campaign_tag or getattr(args, "tag", None) or "standalone"
    setups = normalize_campaign_setups(getattr(args, "setup", ()))
    validate_setup_family(setups)
    if contains_ps_setups(setups):
        raise ValueError("poldis-top is not used for SPINVAL/SPINCOMP/SPINHAD workflows.")
    resolve_analysis_name(args, setups)
    print_stage(f"Running poldis-top for campaign '{label}'")
    output_path = build_poldis_reference(
        base_dir=base_dir,
        unpol=args.unpol,
        pol=args.pol,
        analysis_name=args.analysis_name,
        order=args.order,
        out=args.out,
        setups=setups,
        analysis_variant=analysis_variant_from_args(args),
        campaign_tag=campaign_tag,
        dry_run=False,
    )
    print(f"POLDIS reference YODA: {output_path}")
    return 0


def run_rivetplot_command(args: argparse.Namespace) -> int:
    base_dir = args.base_dir.resolve()
    campaign_dir = resolve_campaign_dir(base_dir, args.tag)
    manifest = load_manifest(campaign_dir)
    requested_setups = resolve_campaign_setups_requested(args.setup, manifest)
    validate_setup_family(requested_setups)
    ps_mode = contains_ps_setups(requested_setups)
    analysis_variant = analysis_variant_from_args(args) or str(manifest.get("analysis_variant", ""))
    if ps_mode:
        analysis_variant = "RIVETPS"
    scale_variations = resolve_scale_variations_requested(getattr(args, "scale_variations", None), manifest, default=False)
    resolve_analysis_name(args, requested_setups)
    args.setup = requested_setups
    args.scale_variations = scale_variations
    args.resolved_analysis_variant = analysis_variant
    print_stage(f"Running rivetplot for campaign '{args.tag}'")
    if ps_mode:
        outputs = run_ps_plotting_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            rivet_mkhtml_tool=args.rivet_mkhtml_tool,
            plot_dir=args.plot_dir,
            scale_variations=scale_variations,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
    else:
        build_poldis_reference(
            base_dir=base_dir,
            unpol=args.unpol,
            pol=args.pol,
            analysis_name=args.analysis_name,
            order=args.order,
            out=args.out,
            setups=args.setup,
            analysis_variant=analysis_variant,
            campaign_tag=args.tag,
            dry_run=False,
        )
        outputs = run_rivetplot_campaign(
            base_dir=base_dir,
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            order=args.order,
            rivet_mkhtml_tool=args.rivet_mkhtml_tool,
            ratiolabel=args.ratiolabel,
            reflabel_prefix=args.reflabel_prefix,
            plot_dir=args.plot_dir,
            scale_variations=scale_variations,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
    print("Rivet comparison plot directories:")
    for setup in normalize_campaign_setups(args.setup):
        if setup in outputs:
            print(f"  {setup}: {outputs[setup]}")
    if "full_toggle" in outputs:
        print(f"  full_toggle: {outputs['full_toggle']}")
    return 0


def run_full_command(args: argparse.Namespace) -> int:
    print_stage(f"Running full workflow for campaign '{args.tag}'")
    rc = run_campaign_command(args)
    if rc != 0:
        return rc
    if args.dry_run:
        print("Skipping post-campaign analysis in dry-run mode.")
        return 0

    campaign_dir = resolve_campaign_dir(args.base_dir.resolve(), args.tag)
    manifest = load_manifest(campaign_dir)
    requested_setups = resolve_campaign_setups_requested(args.setup, manifest)
    validate_setup_family(requested_setups)
    ps_mode = contains_ps_setups(requested_setups)
    analysis_variant = analysis_variant_from_args(args) or str(manifest.get("analysis_variant", ""))
    if ps_mode:
        analysis_variant = "RIVETPS"
    scale_variations = resolve_scale_variations_requested(getattr(args, "scale_variations", None), manifest, default=False)
    resolve_analysis_name(args, requested_setups)
    args.setup = requested_setups
    args.scale_variations = scale_variations
    args.resolved_analysis_variant = analysis_variant
    include_raw_powheg = use_raw_powheg_runtime(
        analysis_variant,
        bool(getattr(args, "raw_powheg", False) or manifest.get("raw_powheg")),
    )

    if ps_mode:
        analyze_ps_herwig_campaign(
            base_dir=args.base_dir.resolve(),
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            scale_variations=scale_variations,
            merge_tool_preference=args.yoda_merge_tool,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
        if multi_family_ps_setups(requested_setups) or has_ps_hadronization_toggle_inputs(requested_setups):
            outputs = run_ps_plotting_campaign(
                base_dir=args.base_dir.resolve(),
                tag=args.tag,
                setups=args.setup,
                analysis_name=args.analysis_name,
                rivet_mkhtml_tool=args.rivet_mkhtml_tool,
                plot_dir=args.plot_dir,
                scale_variations=scale_variations,
                max_workers=max(1, args.jobs),
                dry_run=False,
            )
            print("Rivet comparison plot directories:")
            for setup in normalize_campaign_setups(args.setup):
                if setup in outputs:
                    print(f"  {setup}: {outputs[setup]}")
            if "full_toggle" in outputs:
                print(f"  full_toggle: {outputs['full_toggle']}")
        if "SPINVAL" in requested_setups:
            print("Spin diagnostic outputs:")
            print(f"  {campaign_dir / 'spin_diagnostic.txt'}")
            print(f"  {campaign_dir / 'spin_diagnostic.csv'}")
            print(f"  {campaign_dir / 'spin_diagnostic.json'}")
        return 0

    analyze_herwig_campaign(
        base_dir=args.base_dir.resolve(),
        tag=args.tag,
        setups=args.setup,
        analysis_name=args.analysis_name,
        analysis_variant=analysis_variant,
        scale_variations=scale_variations,
        merge_tool_preference=args.yoda_merge_tool,
        max_workers=max(1, args.jobs),
        dry_run=False,
    )
    if include_raw_powheg:
        analyze_raw_powheg_campaign(
            base_dir=args.base_dir.resolve(),
            tag=args.tag,
            setups=args.setup,
            analysis_name=args.analysis_name,
            analysis_variant=analysis_variant,
            max_workers=max(1, args.jobs),
            dry_run=False,
        )
    build_poldis_reference(
        base_dir=args.base_dir.resolve(),
        unpol=args.unpol,
        pol=args.pol,
        analysis_name=args.analysis_name,
        order=args.order,
        out=args.out,
        setups=args.setup,
        analysis_variant=analysis_variant,
        campaign_tag=args.tag,
        dry_run=False,
    )
    outputs = run_rivetplot_campaign(
        base_dir=args.base_dir.resolve(),
        tag=args.tag,
        setups=args.setup,
        analysis_name=args.analysis_name,
        order=args.order,
        rivet_mkhtml_tool=args.rivet_mkhtml_tool,
        ratiolabel=args.ratiolabel,
        reflabel_prefix=args.reflabel_prefix,
        plot_dir=args.plot_dir,
        scale_variations=scale_variations,
        max_workers=max(1, args.jobs),
        dry_run=False,
    )
    print("Rivet comparison plot directories:")
    for setup in normalize_setups(args.setup):
        print(f"  {setup}: {outputs[setup]}")
    return 0


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.command == "campaign":
        return run_campaign_command(args)
    if args.command == "postprocess":
        return run_postprocess_command(args)
    if args.command == "analyze-herwig":
        return run_analyze_herwig_command(args)
    if args.command == "poldis-top":
        return run_poldis_top_command(args)
    if args.command == "rivetplot":
        return run_rivetplot_command(args)
    if args.command == "full":
        return run_full_command(args)
    raise RuntimeError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
