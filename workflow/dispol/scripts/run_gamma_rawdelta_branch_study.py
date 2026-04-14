#!/usr/bin/env python3.10
"""
Run a focused GAMMA polarized-NLO representation study.

The helper:
  * generates temporary broad/interior GAMMA cards under campaigns/<tag>/<mode>/generated-inputs/
  * patches the PDF profile and polarized-NLO representation switches
  * prepares stale or missing .run files with `Herwig read`
  * launches Herwig shards
  * optionally ensures matching POLDIS references exist
  * writes one results.txt/json report with mode-to-mode and
    Herwig-minus-POLDIS comparisons
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import re
import shlex
import subprocess
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Mapping, Optional, Sequence

from extract_dis_out_results import Measurement, combine_measurements, parse_generated_event_count, parse_measurement


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_TAG = "plain53-gamma-rawdelta-branch"
DEFAULT_JOBS = 8
DEFAULT_POLDIS_JOBS = 4
DEFAULT_POLDIS_VARIANT_JOBS = 1
DEFAULT_SHARDS = 10
DEFAULT_EVENTS_PER_SHARD = 10_000_000
DEFAULT_SEED_BASE = 730_000
DEFAULT_POLDIS_EVENTS = 200_000_000

MODE = "gamma-rawdelta-branch-study"
GAMMA_HELICITIES = ("PP", "PM")
NLO_PIECES = ("POSNLO", "NEGNLO")
BRANCH_MODES = ("legacy", "eff", "raw")
DEFAULT_BRANCH_MODES = "eff,raw"
BRANCH_MODE_LABELS = {
    "legacy": "Legacy",
    "eff": "Uniform A",
    "raw": "Uniform B",
}
PDF_PROFILES = {
    "hybrid": {
        "unpolarized": "PDF4LHC15_nnlo_100_pdfas",
        "polarized_diff": "BDSSV24-NNLO",
    },
    "nnpdf_paired": {
        "unpolarized": "NNPDF40_nlo_pch_as_01180",
        "polarized_diff": "NNPDFpol20_nlo_as_01180",
    },
}


@dataclass(frozen=True)
class CutWindow:
    label: str
    q2_min: float
    q2_max: float
    y_min: float
    y_max: float


@dataclass(frozen=True)
class LogicalRun:
    key: str
    profile: str
    window: str
    branch_mode: str
    piece: str
    helicity: str
    source_card: str
    generated_stem: str
    generated_input_rel: Path
    run_file: str


@dataclass(frozen=True)
class ShardSpec:
    logical_run: LogicalRun
    shard_index: int
    shard_count: int
    tag: str
    seed: int
    events: int


@dataclass(frozen=True)
class ShardResult:
    spec: ShardSpec
    command: List[str]
    returncode: int
    duration_s: float
    launcher_log: Path


BROAD_WINDOW = CutWindow(label="broad", q2_min=49.0, q2_max=2500.0, y_min=0.2, y_max=0.6)
INTERIOR_WINDOW = CutWindow(label="interior", q2_min=100.0, q2_max=1000.0, y_min=0.3, y_max=0.5)
WINDOWS = {"broad": BROAD_WINDOW, "interior": INTERIOR_WINDOW}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used in Herwig shard tags.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DISPOL cards and outputs.",
    )
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Maximum concurrent Herwig processes.")
    parser.add_argument(
        "--poldis-jobs",
        type=int,
        default=DEFAULT_POLDIS_JOBS,
        help="Maximum concurrent POLDIS reference builds.",
    )
    parser.add_argument(
        "--poldis-variant-jobs",
        type=int,
        default=DEFAULT_POLDIS_VARIANT_JOBS,
        help="Maximum concurrent unpolarized/polarized jobs inside each POLDIS reference helper.",
    )
    parser.add_argument("--shards", type=int, default=DEFAULT_SHARDS, help="Shards per logical run.")
    parser.add_argument(
        "--events-per-shard",
        type=int,
        default=DEFAULT_EVENTS_PER_SHARD,
        help="Events per shard.",
    )
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE, help="Base RNG seed.")
    parser.add_argument("--poldis-events", type=int, default=DEFAULT_POLDIS_EVENTS, help="Events per POLDIS run.")
    parser.add_argument(
        "--poldis-totals-json",
        type=Path,
        help="Optional broad-window GAMMA override from an existing campaign-level poldis-totals.json.",
    )
    parser.add_argument(
        "--windows",
        default="broad,interior",
        help="Comma-separated window labels to run (broad,interior).",
    )
    parser.add_argument(
        "--pdf-profiles",
        default="hybrid,nnpdf_paired",
        help="Comma-separated PDF profiles to run (hybrid,nnpdf_paired).",
    )
    parser.add_argument(
        "--branch-modes",
        default=DEFAULT_BRANCH_MODES,
        help="Comma-separated polarized-NLO modes to run (legacy,eff,raw).",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without executing them.")
    parser.add_argument(
        "--skip-herwig",
        action="store_true",
        help="Reuse existing Herwig shard outputs and only run missing POLDIS references plus result collection.",
    )
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip Herwig and POLDIS execution and rebuild the report from existing outputs.",
    )
    parser.add_argument(
        "--skip-poldis",
        action="store_true",
        help="Do not auto-run missing POLDIS references; omit Herwig-vs-POLDIS comparisons when absent.",
    )
    return parser


def campaign_dir(base_dir: Path, tag: str) -> Path:
    return base_dir / "campaigns" / tag / MODE


def generated_input_dir(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "generated-inputs"


def launcher_log_dir(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "launcher-logs"


def monitor_dir(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "monitor"


def results_txt_path(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "results.txt"


def results_json_path(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "results.json"


def status_json_path(base_dir: Path, tag: str) -> Path:
    return monitor_dir(base_dir, tag) / "status.json"


def status_txt_path(base_dir: Path, tag: str) -> Path:
    return monitor_dir(base_dir, tag) / "status.txt"


def poldis_reference_json_path(base_dir: Path, tag: str, profile: str, window: CutWindow) -> Path:
    suffix = "poldis-gamma-interior" if window.label == "interior" else "poldis-gamma-broad"
    poldis_tag = f"{tag}-poldis-{profile}"
    return base_dir / "campaigns" / poldis_tag / suffix / "reference.json"


def build_gamma_source_card(piece: str, helicity: str) -> str:
    return f"DIS-POL-POWHEG_{helicity}-{piece}-GAMMA.in"


def stem_without_suffix(filename: str) -> str:
    return filename[:-3] if filename.endswith(".in") else filename


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def parse_csv_choices(text: str, allowed: Mapping[str, object], label: str) -> List[str]:
    values: List[str] = []
    seen: set[str] = set()
    for raw in text.split(","):
        token = raw.strip()
        if not token:
            continue
        if token not in allowed:
            raise ValueError(f"Unknown {label} {token!r}; allowed values are {', '.join(sorted(allowed))}")
        if token in seen:
            continue
        seen.add(token)
        values.append(token)
    if not values:
        raise ValueError(f"Need at least one {label}.")
    return values


def measurement_payload(measurement: Measurement) -> dict:
    return {"value_pb": measurement.value_pb, "error_pb": measurement.error_pb}


def measurement_from_payload(payload: Mapping[str, float]) -> Measurement:
    return Measurement(float(payload["value_pb"]), float(payload["error_pb"]))


def fmt_measurement(payload: Mapping[str, float]) -> str:
    return f"{payload['value_pb']:.10f} +- {payload['error_pb']:.10f} pb"


def sub_measurement(a: Measurement, b: Measurement) -> Measurement:
    return Measurement(a.value_pb - b.value_pb, math.sqrt(a.error_pb ** 2 + b.error_pb ** 2))


def half_diff(a: Measurement, b: Measurement) -> Measurement:
    delta = sub_measurement(a, b)
    return Measurement(delta.value_pb / 2.0, delta.error_pb / 2.0)


def load_poldis_totals_override(path: Path, profile: str, window: CutWindow) -> dict | None:
    if window.label != "broad":
        return None
    payload = json.loads(path.read_text())
    payload_window = payload.get("window")
    if isinstance(payload_window, str) and payload_window != "broad":
        raise RuntimeError(f"POLDIS totals override {path} is for window {payload_window!r}, expected 'broad'")
    setups = payload.get("setups")
    if not isinstance(setups, Mapping) or "GAMMA" not in setups:
        raise RuntimeError(f"POLDIS totals override {path} does not contain setups.GAMMA")
    gamma = setups["GAMMA"]
    assert isinstance(gamma, Mapping)
    override_profile = gamma.get("pdf_profile", payload.get("pdf_profile"))
    if override_profile is not None and override_profile != profile:
        raise RuntimeError(
            f"POLDIS totals override {path} uses pdf_profile={override_profile!r}, expected {profile!r}"
        )
    pdfs = gamma.get("pdfs", payload.get("pdfs"))
    if not isinstance(pdfs, Mapping):
        raise RuntimeError(f"POLDIS totals override {path} does not define PDF metadata")
    polarized = gamma.get("polarized")
    unpolarized = gamma.get("unpolarized")
    if not isinstance(polarized, Mapping) or not isinstance(unpolarized, Mapping):
        raise RuntimeError(f"POLDIS totals override {path} does not contain unpolarized/polarized totals")
    return {
        "tag": str(payload.get("tag", path.stem)),
        "mode": "poldis-totals-override",
        "window": asdict(window),
        "events": None,
        "pdf_profile": override_profile or profile,
        "pdfs": {
            "unpolarized": str(pdfs["unpolarized"]),
            "polarized_diff": str(pdfs["polarized_diff"]),
        },
        "source_totals_json": str(path),
        "source_reference_json": gamma.get("reference_json"),
        "source_reference_txt": gamma.get("reference_txt"),
        "unpolarized": unpolarized,
        "polarized": polarized,
    }


def build_runs(tag: str, profiles: Sequence[str], windows: Sequence[str], branch_modes: Sequence[str]) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    for profile in profiles:
        profile_token = "HYB" if profile == "hybrid" else "NNPDF"
        for window_name in windows:
            window_token = window_name.upper()
            for branch_mode in branch_modes:
                mode_token = branch_mode.upper()
                for piece in NLO_PIECES:
                    for helicity in GAMMA_HELICITIES:
                        source_card = build_gamma_source_card(piece, helicity)
                        source_stem = stem_without_suffix(source_card)
                        generated_stem = f"{source_stem}-RAWDELTA-{window_token}-{profile_token}-{mode_token}"
                        key = f"{profile}:{window_name}:{branch_mode}:{piece}_{helicity}"
                        runs.append(
                            LogicalRun(
                                key=key,
                                profile=profile,
                                window=window_name,
                                branch_mode=branch_mode,
                                piece=piece,
                                helicity=helicity,
                                source_card=source_card,
                                generated_stem=generated_stem,
                                generated_input_rel=Path("campaigns") / tag / MODE / "generated-inputs" / f"{generated_stem}.in",
                                run_file=f"{generated_stem}.run",
                            )
                        )
    return runs


def patch_card_text(source_text: str, generated_stem: str, window: CutWindow, profile: str, branch_mode: str) -> str:
    pdfs = PDF_PROFILES[profile]
    if branch_mode not in BRANCH_MODE_LABELS:
        raise ValueError(f"Unknown branch mode {branch_mode!r}")
    use_uniform = branch_mode != "legacy"
    raw_enabled = branch_mode == "raw"

    replacements = {
        "set /Herwig/Partons/DiffPDF:PDFName ": f"set /Herwig/Partons/DiffPDF:PDFName {pdfs['polarized_diff']}",
        "set /Herwig/Partons/LHAPDF:PDFName ": f"set /Herwig/Partons/LHAPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/HardLOPDF:PDFName ": f"set /Herwig/Partons/HardLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/ShowerNLOPDF:PDFName ": f"set /Herwig/Partons/ShowerNLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/HardNLOPDF:PDFName ": f"set /Herwig/Partons/HardNLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/ShowerLOPDF:PDFName ": f"set /Herwig/Partons/ShowerLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Cuts/NeutralCurrentCut:MinQ2 ": f"set /Herwig/Cuts/NeutralCurrentCut:MinQ2 {window.q2_min:.0f}.*GeV2",
        "set /Herwig/Cuts/NeutralCurrentCut:MaxQ2 ": f"set /Herwig/Cuts/NeutralCurrentCut:MaxQ2 {window.q2_max:.0f}.*GeV2",
        "set /Herwig/Cuts/NeutralCurrentCut:Miny ": f"set /Herwig/Cuts/NeutralCurrentCut:Miny {window.y_min}",
        "set /Herwig/Cuts/NeutralCurrentCut:Maxy ": f"set /Herwig/Cuts/NeutralCurrentCut:Maxy {window.y_max}",
    }
    branch_settings = {
        "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseUniformPolarizedNLORepresentation ": (
            "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseUniformPolarizedNLORepresentation Yes"
            if use_uniform
            else "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseUniformPolarizedNLORepresentation No"
        ),
        "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseRawFinitePolarizedNLODeltas ": (
            "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseRawFinitePolarizedNLODeltas Yes"
            if raw_enabled
            else "set /Herwig/MatrixElements/PowhegMEDISNCPol:UseRawFinitePolarizedNLODeltas No"
        ),
    }

    found = {prefix: False for prefix in replacements}
    branch_seen = {prefix: False for prefix in branch_settings}
    saverun_found = False
    use_q2_scale_found = False
    patched: List[str] = []

    for line in source_text.splitlines():
        stripped = line.strip()

        replaced = False
        for prefix, replacement in replacements.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                found[prefix] = True
                replaced = True
                break
        if replaced:
            continue

        for prefix, replacement in branch_settings.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                branch_seen[prefix] = True
                replaced = True
                break
        if replaced:
            continue

        if stripped.startswith("set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission "):
            patched.append(line)
            use_q2_scale_found = True
            for prefix, replacement in branch_settings.items():
                if not branch_seen[prefix]:
                    patched.append(replacement)
                    branch_seen[prefix] = True
            continue

        if stripped.startswith("saverun ") and stripped.endswith(" EventGenerator"):
            patched.append(f"saverun {generated_stem} EventGenerator")
            saverun_found = True
            continue

        patched.append(line)

    missing = [prefix.strip() for prefix, seen in found.items() if not seen]
    if missing:
        raise RuntimeError(f"Could not patch required card lines: {', '.join(missing)}")
    if not use_q2_scale_found:
        raise RuntimeError("Failed to find UseQ2ScaleInPOWHEGEmission insertion point.")
    if not all(branch_seen.values()):
        missing_switches = [prefix for prefix, seen in branch_seen.items() if not seen]
        raise RuntimeError(f"Failed to insert polarized-NLO switch lines: {missing_switches}")
    if not saverun_found:
        raise RuntimeError("Could not find saverun line to rewrite in generated card.")

    rendered = "\n".join(patched)
    if source_text.endswith("\n"):
        rendered += "\n"
    return rendered


def materialize_cards(base_dir: Path, runs: Sequence[LogicalRun]) -> None:
    for run in runs:
        source_path = base_dir / run.source_card
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source card {source_path}")
        target_path = base_dir / run.generated_input_rel
        rendered = patch_card_text(
            source_path.read_text(),
            generated_stem=run.generated_stem,
            window=WINDOWS[run.window],
            profile=run.profile,
            branch_mode=run.branch_mode,
        )
        target_path.parent.mkdir(parents=True, exist_ok=True)
        if not target_path.exists() or target_path.read_text() != rendered:
            target_path.write_text(rendered)


def ensure_run_files(base_dir: Path, runs: Sequence[LogicalRun], dry_run: bool) -> None:
    for run in runs:
        in_rel = run.generated_input_rel
        in_path = base_dir / in_rel
        run_path = base_dir / run.run_file
        cmd = ["Herwig", "read", str(in_rel)]
        if dry_run:
            print(shlex.join(cmd))
            continue
        if run_path.exists() and run_path.stat().st_mtime >= in_path.stat().st_mtime:
            continue
        proc = subprocess.run(cmd, cwd=base_dir, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Failed to prepare {run.run_file} from {in_path}")


def build_shards(
    runs: Sequence[LogicalRun],
    tag: str,
    shards: int,
    events_per_shard: int,
    seed_base: int,
) -> List[ShardSpec]:
    out: List[ShardSpec] = []
    seed = seed_base
    for run in runs:
        for shard_index in range(1, shards + 1):
            out.append(
                ShardSpec(
                    logical_run=run,
                    shard_index=shard_index,
                    shard_count=shards,
                    tag=f"{tag}-s{shard_index:03d}",
                    seed=seed,
                    events=events_per_shard,
                )
            )
            seed += 1
    return out


def launcher_log_path(base_dir: Path, tag: str, spec: ShardSpec) -> Path:
    return launcher_log_dir(base_dir, tag) / f"{spec.logical_run.generated_stem}-{spec.tag}.launcher.log"


def build_monitor_payload(
    base_dir: Path,
    tag: str,
    shards: Sequence[ShardSpec],
    results: Sequence[ShardResult],
    started_at: float,
    last_result: ShardResult | None = None,
    phase: str = "running",
) -> dict:
    total = len(shards)
    completed = len(results)
    failures = sum(1 for result in results if result.returncode != 0)
    pending = total - completed
    elapsed_s = time.time() - started_at
    mean_duration_s = sum(result.duration_s for result in results) / completed if completed else None
    eta_s = mean_duration_s * pending if mean_duration_s is not None else None

    per_run: Dict[str, dict] = {}
    for key in sorted({spec.logical_run.key for spec in shards}):
        run_specs = [spec for spec in shards if spec.logical_run.key == key]
        run_results = [result for result in results if result.spec.logical_run.key == key]
        per_run[key] = {
            "total": len(run_specs),
            "completed": len(run_results),
            "pending": len(run_specs) - len(run_results),
            "failures": sum(1 for result in run_results if result.returncode != 0),
        }

    payload = {
        "tag": tag,
        "mode": MODE,
        "phase": phase,
        "updated_unix_s": time.time(),
        "elapsed_s": elapsed_s,
        "total_shards": total,
        "completed_shards": completed,
        "pending_shards": pending,
        "failed_shards": failures,
        "mean_duration_s": mean_duration_s,
        "eta_s": eta_s,
        "per_run": per_run,
        "monitor_files": {
            "status_json": str(status_json_path(base_dir, tag)),
            "status_txt": str(status_txt_path(base_dir, tag)),
        },
    }
    if last_result is not None:
        payload["last_completed"] = {
            "run": last_result.spec.logical_run.key,
            "shard_index": last_result.spec.shard_index,
            "tag": last_result.spec.tag,
            "seed": last_result.spec.seed,
            "returncode": last_result.returncode,
            "duration_s": last_result.duration_s,
            "launcher_log": str(last_result.launcher_log),
        }
    return payload


def render_monitor_text(payload: Mapping[str, object]) -> str:
    lines = [
        f"Tag: {payload['tag']}",
        f"Mode: {payload['mode']}",
        f"Phase: {payload['phase']}",
        f"Completed: {payload['completed_shards']}/{payload['total_shards']}",
        f"Pending: {payload['pending_shards']}",
        f"Failures: {payload['failed_shards']}",
        f"Elapsed: {fmt_seconds(float(payload['elapsed_s']))}",
        f"ETA: {fmt_seconds(float(payload['eta_s'])) if payload['eta_s'] is not None else 'n/a'}",
        "",
        "Per run",
        "-------",
    ]
    per_run = payload["per_run"]
    assert isinstance(per_run, Mapping)
    for key in sorted(per_run):
        item = per_run[key]
        assert isinstance(item, Mapping)
        lines.append(
            f"{key}: completed={item['completed']}/{item['total']} "
            f"pending={item['pending']} failures={item['failures']}"
        )
    last_completed = payload.get("last_completed")
    if isinstance(last_completed, Mapping):
        lines.extend(
            [
                "",
                "Last completed",
                "--------------",
                f"{last_completed['run']}[{int(last_completed['shard_index']):03d}] "
                f"rc={last_completed['returncode']} duration={fmt_seconds(float(last_completed['duration_s']))}",
                f"log: {last_completed['launcher_log']}",
            ]
        )
    return "\n".join(lines).rstrip() + "\n"


def write_monitor_files(base_dir: Path, tag: str, payload: Mapping[str, object]) -> None:
    mon_dir = monitor_dir(base_dir, tag)
    mon_dir.mkdir(parents=True, exist_ok=True)
    status_json_path(base_dir, tag).write_text(json.dumps(payload, indent=2, sort_keys=True))
    status_txt_path(base_dir, tag).write_text(render_monitor_text(payload))


def run_one_shard(base_dir: Path, tag: str, spec: ShardSpec) -> ShardResult:
    cmd = [
        "Herwig",
        "run",
        spec.logical_run.run_file,
        "-N",
        str(spec.events),
        "-t",
        spec.tag,
        "-s",
        str(spec.seed),
    ]
    log_path = launcher_log_path(base_dir, tag, spec)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    started = time.time()
    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.run(cmd, cwd=base_dir, stdout=handle, stderr=subprocess.STDOUT, text=True)
    ended = time.time()
    return ShardResult(
        spec=spec,
        command=cmd,
        returncode=proc.returncode,
        duration_s=ended - started,
        launcher_log=log_path,
    )


def run_shards(base_dir: Path, tag: str, shards: Sequence[ShardSpec], jobs: int) -> List[ShardResult]:
    total = len(shards)
    started = time.time()
    results: List[ShardResult] = []
    write_monitor_files(base_dir, tag, build_monitor_payload(base_dir, tag, shards, results, started, phase="running"))

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        future_map = {executor.submit(run_one_shard, base_dir, tag, spec): spec for spec in shards}
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            phase = "failed" if any(item.returncode != 0 for item in results) else "running"
            if len(results) == total and phase != "failed":
                phase = "complete"
            write_monitor_files(
                base_dir,
                tag,
                build_monitor_payload(
                    base_dir,
                    tag,
                    shards,
                    results,
                    started,
                    last_result=result,
                    phase=phase,
                ),
            )
            if len(results) == 1 or len(results) % 10 == 0 or len(results) == total or result.returncode != 0:
                print(
                    f"[progress:{MODE}] {len(results)}/{total} shards complete, "
                    f"failures={sum(1 for item in results if item.returncode != 0)}, "
                    f"last={result.spec.logical_run.key}[{result.spec.shard_index:03d}] rc={result.returncode}",
                    flush=True,
                )

    bad = [item for item in results if item.returncode != 0]
    if bad:
        summary = "\n".join(
            f"  {item.spec.logical_run.key}[{item.spec.shard_index:03d}] rc={item.returncode} log={item.launcher_log}"
            for item in bad[:10]
        )
        raise RuntimeError(f"{len(bad)} shard(s) failed in {MODE}:\n{summary}")
    return results


def collect_explicit_sharded_measurements(
    base_dir: Path,
    stems: Mapping[str, str],
    tag: str,
    context: str,
) -> Dict[str, dict]:
    resolved: Dict[str, dict] = {}
    for key, stem in stems.items():
        pattern = f"{stem}-S*-{tag}-s*.out"
        paths = sorted(base_dir.rglob(pattern))
        if not paths:
            raise FileNotFoundError(
                f"Missing {context}: could not resolve any shards matching {pattern!r} under {base_dir}"
            )
        measurements: List[Measurement] = []
        generated_events: List[int | None] = []
        for path in paths:
            measurement, _unit = parse_measurement(path.read_text())
            if measurement is None:
                raise RuntimeError(f"Failed to parse cross section from {path}")
            measurements.append(measurement)
            generated_events.append(parse_generated_event_count(path.read_text()))
        resolved[key] = {
            "measurement": combine_measurements(measurements, generated_events),
            "paths": [str(path) for path in paths],
        }
    return resolved


def gamma_summary_from_flat_runs(flat_runs: Mapping[str, dict]) -> dict:
    pos_pp = flat_runs["POSNLO_PP"]["measurement"]
    pos_pm = flat_runs["POSNLO_PM"]["measurement"]
    neg_pp = flat_runs["NEGNLO_PP"]["measurement"]
    neg_pm = flat_runs["NEGNLO_PM"]["measurement"]
    pos_pol = half_diff(pos_pp, pos_pm)
    neg_pol = half_diff(neg_pp, neg_pm)
    nlo_pol = sub_measurement(pos_pol, neg_pol)
    return {
        "runs": {
            "POSNLO_PP": measurement_payload(pos_pp),
            "POSNLO_PM": measurement_payload(pos_pm),
            "NEGNLO_PP": measurement_payload(neg_pp),
            "NEGNLO_PM": measurement_payload(neg_pm),
        },
        "derived": {
            "POSNLO_pol": measurement_payload(pos_pol),
            "NEGNLO_pol": measurement_payload(neg_pol),
            "NLO_pol": measurement_payload(nlo_pol),
        },
    }


def run_poldis_reference(
    base_dir: Path,
    tag: str,
    profile: str,
    window: CutWindow,
    events: int,
    dry_run: bool,
    collect_only: bool,
    skip_poldis: bool,
    poldis_variant_jobs: int,
    poldis_totals_json: Path | None,
) -> dict | None:
    ref_path = poldis_reference_json_path(base_dir, tag, profile, window)
    if ref_path.exists():
        return json.loads(ref_path.read_text())
    if poldis_totals_json is not None:
        override = load_poldis_totals_override(poldis_totals_json.resolve(), profile, window)
        if override is not None:
            return override
    if skip_poldis:
        return None
    if collect_only:
        raise FileNotFoundError(f"Missing POLDIS reference {ref_path}")

    helper = base_dir / "run_poldis_gamma_window_reference.py"
    poldis_tag = f"{tag}-poldis-{profile}"
    cmd = [
        "python3.10",
        str(helper),
        "--base-dir",
        str(base_dir),
        "--tag",
        poldis_tag,
        "--window",
        window.label,
        "--pdf-profile",
        profile,
        "--events",
        str(events),
        "--jobs",
        str(poldis_variant_jobs),
    ]
    if dry_run:
        print(shlex.join(cmd))
        return None
    proc = subprocess.run(cmd, cwd=base_dir, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Failed to build POLDIS reference for {profile}:{window.label}")
    if not ref_path.exists():
        raise FileNotFoundError(f"Expected POLDIS reference {ref_path} after helper run")
    return json.loads(ref_path.read_text())


def collect_poldis_references(
    base_dir: Path,
    tag: str,
    profiles: Sequence[str],
    windows: Sequence[str],
    events: int,
    dry_run: bool,
    collect_only: bool,
    skip_poldis: bool,
    poldis_jobs: int,
    poldis_variant_jobs: int,
    poldis_totals_json: Path | None,
) -> Dict[str, Dict[str, dict | None]]:
    refs: Dict[str, Dict[str, dict | None]] = {profile: {} for profile in profiles}
    tasks = [(profile, window_name) for profile in profiles for window_name in windows]

    if not tasks:
        return refs

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, poldis_jobs)) as executor:
        future_map = {
            executor.submit(
                run_poldis_reference,
                base_dir=base_dir,
                tag=tag,
                profile=profile,
                window=WINDOWS[window_name],
                events=events,
                dry_run=dry_run,
                collect_only=collect_only,
                skip_poldis=skip_poldis,
                poldis_variant_jobs=poldis_variant_jobs,
                poldis_totals_json=poldis_totals_json,
            ): (profile, window_name)
            for profile, window_name in tasks
        }
        for future in concurrent.futures.as_completed(future_map):
            profile, window_name = future_map[future]
            refs[profile][window_name] = future.result()
            print(
                f"[progress:{MODE}] POLDIS {profile}:{window_name} ready "
                f"({sum(len(item) for item in refs.values())}/{len(tasks)})",
                flush=True,
            )

    return refs


def build_payload(
    base_dir: Path,
    tag: str,
    profiles: Sequence[str],
    windows: Sequence[str],
    branch_modes: Sequence[str],
    runs: Sequence[LogicalRun],
    shards: int,
    events_per_shard: int,
    poldis_events: int,
    poldis_jobs: int,
    poldis_variant_jobs: int,
    poldis_refs: Mapping[str, Mapping[str, dict | None]],
) -> dict:
    payload = {
        "tag": tag,
        "mode": MODE,
        "settings": {
            "profiles": list(profiles),
            "windows": list(windows),
            "branch_modes": list(branch_modes),
            "shards": shards,
            "events_per_shard": events_per_shard,
            "poldis_events": poldis_events,
            "poldis_jobs": poldis_jobs,
            "poldis_variant_jobs": poldis_variant_jobs,
        },
        "profiles": {},
    }

    run_index = {
        (run.profile, run.window, run.branch_mode, run.piece, run.helicity): run
        for run in runs
    }

    for profile in profiles:
        profile_payload = {
            "pdfs": dict(PDF_PROFILES[profile]),
            "windows": {},
        }
        for window_name in windows:
            window_payload: dict = {
                "window": asdict(WINDOWS[window_name]),
                "modes": {},
            }
            for branch_mode in branch_modes:
                stems = {
                    f"{piece}_{helicity}": run_index[(profile, window_name, branch_mode, piece, helicity)].generated_stem
                    for piece in NLO_PIECES
                    for helicity in GAMMA_HELICITIES
                }
                flat_runs = collect_explicit_sharded_measurements(
                    base_dir,
                    stems,
                    tag,
                    context=f"{profile}:{window_name}:{branch_mode}",
                )
                window_payload["modes"][branch_mode] = gamma_summary_from_flat_runs(flat_runs)

            ordered_pairs = (("eff", "legacy"), ("raw", "legacy"), ("raw", "eff"))
            mode_differences: dict = {}
            for upper_mode, lower_mode in ordered_pairs:
                if upper_mode not in branch_modes or lower_mode not in branch_modes:
                    continue
                upper = window_payload["modes"][upper_mode]["derived"]
                lower = window_payload["modes"][lower_mode]["derived"]
                mode_differences[f"{upper_mode}_minus_{lower_mode}"] = {
                    "label": f"{BRANCH_MODE_LABELS[upper_mode]} - {BRANCH_MODE_LABELS[lower_mode]}",
                    "POSNLO_pol": measurement_payload(
                        sub_measurement(
                            measurement_from_payload(upper["POSNLO_pol"]),
                            measurement_from_payload(lower["POSNLO_pol"]),
                        )
                    ),
                    "NEGNLO_pol": measurement_payload(
                        sub_measurement(
                            measurement_from_payload(upper["NEGNLO_pol"]),
                            measurement_from_payload(lower["NEGNLO_pol"]),
                        )
                    ),
                    "NLO_pol": measurement_payload(
                        sub_measurement(
                            measurement_from_payload(upper["NLO_pol"]),
                            measurement_from_payload(lower["NLO_pol"]),
                        )
                    ),
                }
            window_payload["mode_differences"] = mode_differences

            ref_payload = poldis_refs.get(profile, {}).get(window_name)
            if ref_payload is not None:
                poldis_nlo_pol = measurement_from_payload(ref_payload["polarized"]["NLO"])
                window_payload["poldis_reference"] = {
                    "path": str(
                        ref_payload.get(
                            "source_reference_json",
                            ref_payload.get(
                                "source_totals_json",
                                poldis_reference_json_path(base_dir, tag, profile, WINDOWS[window_name]),
                            ),
                        )
                    ),
                    "tag": str(ref_payload.get("tag", f"{tag}-poldis-{profile}")),
                    "pdf_profile": ref_payload.get("pdf_profile"),
                    "pdfs": ref_payload.get("pdfs"),
                    "polarized_NLO": ref_payload["polarized"]["NLO"],
                }
                window_payload["mode_minus_poldis_NLO_pol"] = {
                    branch_mode: measurement_payload(
                        sub_measurement(
                            measurement_from_payload(window_payload["modes"][branch_mode]["derived"]["NLO_pol"]),
                            poldis_nlo_pol,
                        )
                    )
                    for branch_mode in branch_modes
                }

            profile_payload["windows"][window_name] = window_payload
        payload["profiles"][profile] = profile_payload
    return payload


def render_results_text(payload: dict) -> str:
    lines = [
        "GAMMA polarized-NLO representation study",
        "========================================",
        f"Tag: {payload['tag']}",
        f"Profiles: {', '.join(payload['settings']['profiles'])}",
        f"Windows: {', '.join(payload['settings']['windows'])}",
        f"Branch modes: {', '.join(payload['settings']['branch_modes'])}",
        f"Shards per logical run: {payload['settings']['shards']}",
        f"Events per shard: {payload['settings']['events_per_shard']}",
        f"POLDIS events per run: {payload['settings']['poldis_events']}",
        f"Concurrent POLDIS runs: {payload['settings']['poldis_jobs']}",
        f"Concurrent variants per POLDIS run: {payload['settings']['poldis_variant_jobs']}",
    ]

    for profile in payload["settings"]["profiles"]:
        profile_payload = payload["profiles"][profile]
        lines.extend(
            [
                "",
                f"PDF profile: {profile}",
                "--------------------",
                f"Unpolarized PDF: {profile_payload['pdfs']['unpolarized']}",
                f"Polarized diff PDF: {profile_payload['pdfs']['polarized_diff']}",
            ]
        )
        for window_name in payload["settings"]["windows"]:
            window_payload = profile_payload["windows"][window_name]
            window = window_payload["window"]
            lines.extend(
                [
                    "",
                    f"{window_name.capitalize()} window: Q^2 in [{window['q2_min']:.0f}, {window['q2_max']:.0f}] GeV^2, y in [{window['y_min']:.1f}, {window['y_max']:.1f}]",
                ]
            )
            for branch_mode in payload["settings"]["branch_modes"]:
                derived = window_payload["modes"][branch_mode]["derived"]
                label = BRANCH_MODE_LABELS[branch_mode]
                lines.extend(
                    [
                        f"{label} POSNLO_pol: {fmt_measurement(derived['POSNLO_pol'])}",
                        f"{label} NEGNLO_pol: {fmt_measurement(derived['NEGNLO_pol'])}",
                        f"{label} NLO_pol:    {fmt_measurement(derived['NLO_pol'])}",
                    ]
                )
            for diff_key in ("eff_minus_legacy", "raw_minus_legacy", "raw_minus_eff"):
                if diff_key not in window_payload["mode_differences"]:
                    continue
                diff = window_payload["mode_differences"][diff_key]
                lines.extend(
                    [
                        f"{diff['label']} POSNLO_pol: {fmt_measurement(diff['POSNLO_pol'])}",
                        f"{diff['label']} NEGNLO_pol: {fmt_measurement(diff['NEGNLO_pol'])}",
                        f"{diff['label']} NLO_pol:    {fmt_measurement(diff['NLO_pol'])}",
                    ]
                )
            if "poldis_reference" in window_payload:
                lines.append(f"POLDIS NLO_pol:         {fmt_measurement(window_payload['poldis_reference']['polarized_NLO'])}")
                for branch_mode in payload["settings"]["branch_modes"]:
                    label = BRANCH_MODE_LABELS[branch_mode]
                    lines.append(
                        f"{label} - POLDIS:      {fmt_measurement(window_payload['mode_minus_poldis_NLO_pol'][branch_mode])}"
                    )
            else:
                lines.append("POLDIS reference: missing")
    return "\n".join(lines).rstrip() + "\n"


def write_results(base_dir: Path, tag: str, payload: dict) -> None:
    out_dir = campaign_dir(base_dir, tag)
    out_dir.mkdir(parents=True, exist_ok=True)
    text = render_results_text(payload)
    results_txt_path(base_dir, tag).write_text(text)
    results_json_path(base_dir, tag).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(text, end="")
    print(f"Wrote text results: {results_txt_path(base_dir, tag)}")
    print(f"Wrote JSON results: {results_json_path(base_dir, tag)}")


def run(args: argparse.Namespace) -> None:
    base_dir = args.base_dir.resolve()
    profiles = parse_csv_choices(args.pdf_profiles, PDF_PROFILES, "PDF profile")
    windows = parse_csv_choices(args.windows, WINDOWS, "window")
    branch_modes = parse_csv_choices(args.branch_modes, BRANCH_MODE_LABELS, "branch mode")
    runs = build_runs(args.tag, profiles, windows, branch_modes)
    shards = build_shards(runs, args.tag, args.shards, args.events_per_shard, args.seed_base)
    shard_results: List[ShardResult] = []
    started_at = time.time()

    if not args.collect_only and not args.skip_herwig:
        materialize_cards(base_dir, runs)
        ensure_run_files(base_dir, runs, dry_run=args.dry_run)
        if args.dry_run:
            for spec in shards:
                cmd = [
                    "Herwig",
                    "run",
                    spec.logical_run.run_file,
                    "-N",
                    str(spec.events),
                    "-t",
                    spec.tag,
                    "-s",
                    str(spec.seed),
                ]
                print(shlex.join(cmd))
        else:
            shard_results = run_shards(base_dir, args.tag, shards, args.jobs)

    if not args.dry_run:
        write_monitor_files(
            base_dir,
            args.tag,
            build_monitor_payload(
                base_dir,
                args.tag,
                shards,
                shard_results,
                started_at,
                phase="running-poldis",
            ),
        )

    poldis_refs = collect_poldis_references(
        base_dir=base_dir,
        tag=args.tag,
        profiles=profiles,
        windows=windows,
        events=args.poldis_events,
        dry_run=args.dry_run,
        collect_only=args.collect_only,
        skip_poldis=args.skip_poldis,
        poldis_jobs=args.poldis_jobs,
        poldis_variant_jobs=args.poldis_variant_jobs,
        poldis_totals_json=args.poldis_totals_json,
    )

    if args.dry_run:
        return

    payload = build_payload(
        base_dir=base_dir,
        tag=args.tag,
        profiles=profiles,
        windows=windows,
        branch_modes=branch_modes,
        runs=runs,
        shards=args.shards,
        events_per_shard=args.events_per_shard,
        poldis_events=args.poldis_events,
        poldis_jobs=args.poldis_jobs,
        poldis_variant_jobs=args.poldis_variant_jobs,
        poldis_refs=poldis_refs,
    )
    write_results(base_dir, args.tag, payload)
    write_monitor_files(
        base_dir,
        args.tag,
        build_monitor_payload(
            base_dir,
            args.tag,
            shards,
            shard_results,
            started_at,
            phase="complete",
        ),
    )


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.base_dir.resolve().exists():
        raise SystemExit(f"Base directory does not exist: {args.base_dir}")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
