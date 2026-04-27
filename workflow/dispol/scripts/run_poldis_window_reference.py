#!/usr/bin/env python3.10
"""
Build focused POLDIS total-cross-section references for a chosen DIS setup/window.

The helper:
  * templates user_dijet_rivetplots.f into a campaign-local work area
  * optionally splits the unpolarized/polarized runs into independent seeded shards
  * compiles and runs local poldis.x shard jobs in campaign-local work areas
  * parses the printed LO/NLO/NNLO inclusive totals
  * combines shard totals as event-weighted means with independent errors added
    in quadrature
  * writes reference.txt and reference.json
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import re
import shlex
import shutil
import subprocess
import threading
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Callable, Dict, List, Mapping, Optional, Sequence

from extract_dis_out_results import (
    Measurement,
    POLDIS_POL_REFS,
    POLDIS_UNPOL_REFS,
    combine_measurements,
)


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_ROOT_DIR = DEFAULT_BASE_DIR.parent.parent.parent.parent
DEFAULT_POLDIS_DIR = DEFAULT_ROOT_DIR / "POLDIS" / "POLDIS-public"
DEFAULT_TAG = "plain55-poldis-ref"
DEFAULT_EVENTS = 1_600_000_000
DEFAULT_PDF_PROFILE = "nnpdf_paired"
DEFAULT_JOBS = 1
DEFAULT_SHARDS = 96
DEFAULT_SEED_BASE = 14_217_136
SCALE_VARIATION_ORDER = ("nominal", "ScaleFactorDown", "ScaleFactorUp")
SCALE_VARIATION_MU2_FACTORS = {
    "nominal": "1D0",
    "ScaleFactorDown": "0.25D0",
    "ScaleFactorUp": "4D0",
}
REFERENCE_ERROR_MODES = ("stat", "scale", "stat+scale")
DEFAULT_REFERENCE_ERROR_MODE = "stat"
REFERENCE_ERROR_MODE_SUFFIXES = {
    "stat": "stat",
    "scale": "scale",
    "stat+scale": "stat_scale",
}

PDF_PROFILES: Dict[str, Dict[str, str]] = {
    "hybrid": {
        "unpolarized": "PDF4LHC15_nnlo_100_pdfas",
        "polarized_diff": "BDSSV24-NNLO",
    },
    "nnpdf_paired": {
        "unpolarized": "NNPDF40_nlo_pch_as_01180",
        "polarized_diff": "NNPDFpol20_nlo_as_01180",
    },
}

SETUP_TO_IBOSON = {
    "GAMMA": 0,
    "ALL": 1,
    "Z": 2,
    "CC": 3,
}

TOP_SUFFIX = {
    "GAMMA": "_GAM",
    "ALL": "",
    "Z": "_Z",
    "CC": "_W",
}


@dataclass(frozen=True)
class CutWindow:
    label: str
    q2_min: float
    q2_max: float
    y_min: float
    y_max: float


@dataclass(frozen=True)
class ShardSpec:
    label: str
    ipol: int
    scale_variation: str
    shard_index: int
    shard_count: int
    events: int
    seed: int


# Keep the production "broad" y/Q^2_max window, but match the rerun cards'
# raised Q^2 threshold so plain/RIVETFO POLDIS references stay aligned.
BROAD_WINDOW = CutWindow(label="broad", q2_min=100.0, q2_max=2500.0, y_min=0.2, y_max=0.6)
INTERIOR_WINDOW = CutWindow(label="interior", q2_min=100.0, q2_max=1000.0, y_min=0.3, y_max=0.5)
WINDOWS = {"broad": BROAD_WINDOW, "interior": INTERIOR_WINDOW}

TOTAL_RE = re.compile(
    r"^(?P<order>NNLO|NLO|LO)\s*=\s*(?P<value>[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*\+\-\s*"
    r"(?P<error>[+-]?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s*$"
)
PROGRESS_RE = re.compile(r"^\s*(?P<events>\d+),\s*ISEED=")

RUNTIME_FILES = ("poldis.f", "gbook.f", "jetalg.f")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used for the reference work area.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DISPOL helpers.",
    )
    parser.add_argument(
        "--poldis-dir",
        type=Path,
        default=DEFAULT_POLDIS_DIR,
        help="Directory containing the POLDIS source tree and compile helper.",
    )
    parser.add_argument(
        "--setup",
        choices=tuple(SETUP_TO_IBOSON),
        default="GAMMA",
        help="Electroweak setup to evaluate.",
    )
    parser.add_argument(
        "--window",
        choices=tuple(WINDOWS),
        default="broad",
        help="DIS window to evaluate.",
    )
    parser.add_argument(
        "--pdf-profile",
        choices=tuple(PDF_PROFILES),
        default=DEFAULT_PDF_PROFILE,
        help="Named PDF profile used for the unpolarized and polarized POLDIS runs.",
    )
    parser.add_argument(
        "--unpolarized-pdf",
        help="Optional override for the unpolarized POLDIS PDF set name.",
    )
    parser.add_argument(
        "--polarized-diff-pdf",
        help="Optional override for the polarized-difference POLDIS PDF set name.",
    )
    parser.add_argument(
        "--events",
        type=int,
        default=DEFAULT_EVENTS,
        help="Total number of POLDIS events per unpolarized/polarized reference after combining shards.",
    )
    parser.add_argument(
        "--shards",
        type=int,
        default=DEFAULT_SHARDS,
        help="Split each unpolarized/polarized reference run into this many independent seeded shards.",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=DEFAULT_SEED_BASE,
        help="Base random seed used to generate distinct per-shard seeds.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=DEFAULT_JOBS,
        help="Maximum concurrent POLDIS shard jobs inside this reference build.",
    )
    parser.add_argument(
        "--scale-variations",
        action="store_true",
        help="Also run common muR=muF scale variations matching Herwig ScaleFactor down/up by a factor of two.",
    )
    parser.add_argument(
        "--error-mode",
        choices=REFERENCE_ERROR_MODES,
        default=DEFAULT_REFERENCE_ERROR_MODE,
        help="Reference YODA error bars to write into reference.yoda.gz: statistical only, scale envelope only, or both in quadrature.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned compile/run commands only.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip compile/run and rebuild the parsed reference from existing logs.",
    )
    return parser


def work_dir(base_dir: Path, tag: str, setup: str, window: CutWindow) -> Path:
    return base_dir / "campaigns" / tag / f"poldis-{setup.lower()}-{window.label}"


def variant_dir(base_dir: Path, tag: str, setup: str, window: CutWindow, label: str) -> Path:
    return work_dir(base_dir, tag, setup, window) / label


def scale_variation_dir(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    label: str,
    scale_variation: str,
) -> Path:
    base = variant_dir(base_dir, tag, setup, window, label)
    if scale_variation == "nominal":
        return base
    return base / "scale-variations" / scale_variation


def variant_shards_dir(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    label: str,
    scale_variation: str,
) -> Path:
    return scale_variation_dir(base_dir, tag, setup, window, label, scale_variation) / "shards"


def variant_shard_dir(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    label: str,
    scale_variation: str,
    shard_index: int,
) -> Path:
    return variant_shards_dir(base_dir, tag, setup, window, label, scale_variation) / f"s{shard_index:03d}"


def variant_plan_path(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    label: str,
    scale_variation: str,
) -> Path:
    return scale_variation_dir(base_dir, tag, setup, window, label, scale_variation) / "shards.json"


def reference_txt_path(base_dir: Path, tag: str, setup: str, window: CutWindow) -> Path:
    return work_dir(base_dir, tag, setup, window) / "reference.txt"


def reference_json_path(base_dir: Path, tag: str, setup: str, window: CutWindow) -> Path:
    return work_dir(base_dir, tag, setup, window) / "reference.json"


def reference_yoda_path(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
) -> Path:
    return work_dir(base_dir, tag, setup, window) / "reference.yoda.gz"


def reference_yoda_variation_path(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    scale_variation: str = "nominal",
) -> Path:
    filename = (
        "reference.nominal.yoda.gz"
        if scale_variation == "nominal"
        else f"reference.{scale_variation}.yoda.gz"
    )
    return work_dir(base_dir, tag, setup, window) / filename


def reference_yoda_error_mode_path(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    error_mode: str,
) -> Path:
    suffix = REFERENCE_ERROR_MODE_SUFFIXES[error_mode]
    filename = f"reference.{suffix}.yoda.gz"
    return work_dir(base_dir, tag, setup, window) / filename


def monitor_dir(base_dir: Path, tag: str, setup: str, window: CutWindow) -> Path:
    return work_dir(base_dir, tag, setup, window) / "monitor"


def status_txt_path(base_dir: Path, tag: str, setup: str, window: CutWindow) -> Path:
    return monitor_dir(base_dir, tag, setup, window) / "status.txt"


def status_json_path(base_dir: Path, tag: str, setup: str, window: CutWindow) -> Path:
    return monitor_dir(base_dir, tag, setup, window) / "status.json"


def measurement_payload(measurement: Measurement) -> dict:
    return {"value_pb": measurement.value_pb, "error_pb": measurement.error_pb}


def measurement_from_payload(payload: Mapping[str, float]) -> Measurement:
    return Measurement(float(payload["value_pb"]), float(payload["error_pb"]))


def fmt_measurement(payload: Mapping[str, float]) -> str:
    return f"{payload['value_pb']:.10f} +- {payload['error_pb']:.10f} pb"


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def sub_measurement(a: Measurement, b: Measurement) -> Measurement:
    return Measurement(a.value_pb - b.value_pb, math.sqrt(a.error_pb ** 2 + b.error_pb ** 2))


def parse_poldis_totals(text: str, context: str) -> Dict[str, Measurement]:
    totals: Dict[str, Measurement] = {}
    for line in text.splitlines():
        match = TOTAL_RE.match(line.strip())
        if not match:
            continue
        totals[match.group("order")] = Measurement(
            value_pb=float(match.group("value")),
            error_pb=float(match.group("error")),
        )
    missing = [order for order in ("LO", "NLO", "NNLO") if order not in totals]
    if missing:
        raise RuntimeError(f"Failed to parse {context}: missing {', '.join(missing)} totals")
    return totals


def compile_command(poldis_dir: Path) -> list[str]:
    compile_script = poldis_dir / "compile_dijet_rivetplots"
    if not compile_script.exists():
        raise FileNotFoundError(f"Missing compile helper {compile_script}")
    command = compile_script.read_text().strip()
    if not command:
        raise RuntimeError(f"Compile helper {compile_script} is empty")
    return shlex.split(command)


def resolve_pdf_choice(profile: str, unpolarized_override: str | None, polarized_override: str | None) -> Dict[str, str]:
    if profile not in PDF_PROFILES:
        raise KeyError(f"Unknown PDF profile {profile!r}")
    resolved = dict(PDF_PROFILES[profile])
    if unpolarized_override:
        resolved["unpolarized"] = unpolarized_override
    if polarized_override:
        resolved["polarized_diff"] = polarized_override
    return resolved


def is_builtin_hybrid_pdf_choice(pdfs: Mapping[str, str]) -> bool:
    return (
        pdfs.get("unpolarized") == PDF_PROFILES["hybrid"]["unpolarized"]
        and pdfs.get("polarized_diff") == PDF_PROFILES["hybrid"]["polarized_diff"]
    )


def builtin_reference_pdf_profile(setup: str) -> str:
    return "hybrid" if str(setup).upper() == "CC" else "nnpdf_paired"


def uses_builtin_broad_reference_choice(setup: str, pdfs: Mapping[str, str]) -> bool:
    profile = builtin_reference_pdf_profile(setup)
    return (
        pdfs.get("unpolarized") == PDF_PROFILES[profile]["unpolarized"]
        and pdfs.get("polarized_diff") == PDF_PROFILES[profile]["polarized_diff"]
    )


def enabled_scale_variations(include_variations: bool) -> tuple[str, ...]:
    return SCALE_VARIATION_ORDER if include_variations else ("nominal",)


def scale_variation_task_label(scale_variation: str, label: str) -> str:
    if scale_variation == "nominal":
        return label
    return f"{scale_variation}-{label}"


def split_total_events(total_events: int, shard_count: int) -> List[int]:
    if shard_count < 1:
        raise ValueError("shard_count must be positive")
    if total_events < shard_count:
        raise ValueError(
            f"Cannot split {total_events} total events across {shard_count} shards; "
            "choose fewer shards or more events."
        )
    base, remainder = divmod(total_events, shard_count)
    return [base + (1 if index < remainder else 0) for index in range(shard_count)]


def variant_seed_offset(ipol: int) -> int:
    return 1_000_000 if ipol else 0


def build_variant_shards(
    *,
    label: str,
    ipol: int,
    scale_variation: str,
    total_events: int,
    shard_count: int,
    seed_base: int,
) -> List[ShardSpec]:
    event_counts = split_total_events(total_events, shard_count)
    specs: List[ShardSpec] = []
    for shard_index, shard_events in enumerate(event_counts, start=1):
        specs.append(
            ShardSpec(
                label=label,
                ipol=ipol,
                scale_variation=scale_variation,
                shard_index=shard_index,
                shard_count=shard_count,
                events=shard_events,
                seed=seed_base + variant_seed_offset(ipol) + shard_index - 1,
            )
        )
    return specs


def patch_user_source(
    source_text: str,
    *,
    ipol: int,
    setup: str,
    window: CutWindow,
    scale_variation: str,
    events: int,
    seed: int,
    unpolarized_pdf: str,
    polarized_diff_pdf: str,
) -> str:
    replacements = {
        "IPOL=": f"      IPOL={ipol}",
        "IBOSON=": f"      IBOSON={SETUP_TO_IBOSON[setup]}",
        "LEPCH=": "      LEPCH=-1D0",
        "IFRAME=": "      IFRAME=1",
        "ICH=": "      ICH=0",
        "INNLO=": "      INNLO=1",
        "NEV=": f"      NEV={events}",
        "YMIN=": f"      YMIN={window.y_min}",
        "YMAX=": f"      YMAX={window.y_max}",
        "Q2MIN=": f"      Q2MIN={window.q2_min:.0f}",
        "Q2MAX=": f"      Q2MAX={window.q2_max:.0f}",
        "R(1)=": f"      R(1)={seed}",
        "SCALEJET=Q2": f"      SCALEJET={SCALE_VARIATION_MU2_FACTORS[scale_variation]}*Q2",
    }
    seen = {key: False for key in replacements}
    seen_unpolarized_pdf = False
    seen_polarized_pdf = False
    patched = []
    pdf_branch: str | None = None
    for line in source_text.splitlines():
        stripped = line.strip()
        if stripped.startswith("IF (IPOL.EQ.0) THEN"):
            pdf_branch = "unpolarized"
            patched.append(line)
            continue
        if stripped.startswith("ELSEIF (IPOL.EQ.1) THEN"):
            pdf_branch = "polarized_diff"
            patched.append(line)
            continue
        if stripped.startswith("ENDIF"):
            pdf_branch = None
            patched.append(line)
            continue
        if stripped.startswith("call InitPDFsetByName("):
            if pdf_branch == "unpolarized":
                patched.append(f'         call InitPDFsetByName("{unpolarized_pdf}")')
                seen_unpolarized_pdf = True
                continue
            if pdf_branch == "polarized_diff":
                patched.append(f'         call InitPDFsetByName("{polarized_diff_pdf}")')
                seen_polarized_pdf = True
                continue
        replaced = False
        for prefix, replacement in replacements.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                seen[prefix] = True
                replaced = True
                break
        if not replaced:
            patched.append(line)
    missing = [prefix for prefix, ok in seen.items() if not ok]
    if missing:
        raise RuntimeError(f"Failed to patch POLDIS source lines: {missing}")
    if not seen_unpolarized_pdf or not seen_polarized_pdf:
        raise RuntimeError("Failed to patch POLDIS PDF set lines in user_dijet_rivetplots.f")
    rendered = "\n".join(patched)
    if source_text.endswith("\n"):
        rendered += "\n"
    return rendered


def expected_top_name(setup: str, ipol: int) -> str:
    prefix = "dijets_unpol" if ipol == 0 else "dijets_pol"
    suffix = TOP_SUFFIX[setup]
    return f"{prefix}{suffix}.top"


def _coerce_nonnegative_int(value: object, default: int = 0) -> int:
    try:
        return max(0, int(value))
    except (TypeError, ValueError):
        return default


def build_job_summary(
    *,
    phase: str,
    total_jobs: int = 1,
    job_counts: Optional[Mapping[str, object]] = None,
) -> Dict[str, int]:
    if job_counts is not None:
        running = _coerce_nonnegative_int(job_counts.get("running"))
        compiling = _coerce_nonnegative_int(job_counts.get("compiling"))
        pending = _coerce_nonnegative_int(job_counts.get("pending"))
        completed = _coerce_nonnegative_int(job_counts.get("completed"))
        total = _coerce_nonnegative_int(job_counts.get("total"))
        accounted = running + compiling + pending + completed
        if total <= 0:
            total = accounted
        elif accounted < total:
            pending += total - accounted
        elif accounted > total:
            total = accounted
        return {
            "running": running,
            "compiling": compiling,
            "pending": pending,
            "completed": completed,
            "active": running + compiling,
            "total": total,
        }

    total = max(0, int(total_jobs))
    running = 0
    compiling = 0
    pending = 0
    completed = 0
    phase_lower = str(phase).lower()
    if total > 0:
        if phase_lower == "complete" or phase_lower.startswith("completed-"):
            completed = total
        elif phase_lower.startswith("running-") or phase_lower == "collecting":
            running = min(1, total)
            pending = max(0, total - running)
        elif phase_lower.startswith("compiling-"):
            compiling = min(1, total)
            pending = max(0, total - compiling)
        else:
            pending = total
    return {
        "running": running,
        "compiling": compiling,
        "pending": pending,
        "completed": completed,
        "active": running + compiling,
        "total": total,
    }


def render_monitor_text(payload: Mapping[str, object]) -> str:
    lines = [
        f"Tag: {payload['tag']}",
        f"Setup: {payload['setup']}",
        f"Window: {payload['window']}",
        f"Phase: {payload['phase']}",
        f"Variant: {payload['variant']}",
        f"Elapsed: {fmt_seconds(float(payload['elapsed_s']))}",
    ]
    if payload.get("events_total") is not None:
        lines.extend(
            [
                f"Events done: {payload['events_done']}/{payload['events_total']}",
                f"Variant progress: {float(payload['variant_percent']):.2f}%",
                f"Overall progress: {float(payload['overall_percent']):.2f}%",
            ]
        )
    jobs = payload.get("jobs")
    if isinstance(jobs, dict) and jobs:
        lines.append(
            "Shard jobs: "
            f"running {jobs.get('running', 0)} | "
            f"compiling {jobs.get('compiling', 0)} | "
            f"pending {jobs.get('pending', 0)} | "
            f"completed {jobs.get('completed', 0)}"
        )
    active_variant = payload.get("active_variant")
    if active_variant:
        lines.append(f"Active variant: {active_variant}")
    if payload.get("last_progress_line"):
        lines.append(f"Last progress line: {payload['last_progress_line']}")
    if payload.get("log"):
        lines.append(f"Log: {payload['log']}")
    variants = payload.get("variants")
    if isinstance(variants, dict) and variants:
        lines.append("Parallel variants:")
        for label in sorted(variants):
            info = variants.get(label)
            if not isinstance(info, dict):
                continue
            phase = info.get("phase", "pending")
            percent = float(info.get("percent", 0.0))
            events_done = int(info.get("events_done", 0))
            events_total = info.get("events_total")
            if events_total is None:
                lines.append(f"  {label}: {phase}")
            else:
                lines.append(f"  {label}: {phase} {events_done}/{events_total} ({percent:.2f}%)")
            if info.get("last_progress_line"):
                lines.append(f"    last: {info['last_progress_line']}")
            if info.get("log"):
                lines.append(f"    log: {info['log']}")
    return "\n".join(lines).rstrip() + "\n"


def write_monitor_files(base_dir: Path, tag: str, setup: str, window: CutWindow, payload: Mapping[str, object]) -> None:
    mon_dir = monitor_dir(base_dir, tag, setup, window)
    mon_dir.mkdir(parents=True, exist_ok=True)
    status_json_path(base_dir, tag, setup, window).write_text(json.dumps(payload, indent=2, sort_keys=True))
    status_txt_path(base_dir, tag, setup, window).write_text(render_monitor_text(payload))


def build_monitor_payload(
    *,
    tag: str,
    setup: str,
    window: CutWindow,
    phase: str,
    variant: str,
    started_at: float,
    run_index: int,
    run_count: int,
    events_total: int | None = None,
    events_done: int = 0,
    last_progress_line: str | None = None,
    log_path: Path | None = None,
    total_jobs: int = 1,
    job_counts: Optional[Mapping[str, object]] = None,
) -> dict:
    variant_fraction = 0.0
    if events_total and events_total > 0:
        variant_fraction = min(max(events_done / events_total, 0.0), 1.0)
    overall_fraction = (run_index + variant_fraction) / max(run_count, 1)
    return {
        "tag": tag,
        "setup": setup,
        "window": window.label,
        "phase": phase,
        "variant": variant,
        "elapsed_s": time.time() - started_at,
        "events_total": events_total,
        "events_done": events_done,
        "variant_percent": 100.0 * variant_fraction,
        "overall_percent": 100.0 * overall_fraction,
        "last_progress_line": last_progress_line,
        "log": str(log_path) if log_path is not None else None,
        "jobs": build_job_summary(phase=phase, total_jobs=total_jobs, job_counts=job_counts),
    }


def build_parallel_monitor_payload(
    *,
    tag: str,
    setup: str,
    window: CutWindow,
    started_at: float,
    variant_states: Mapping[str, Mapping[str, object]],
    active_variant: str | None = None,
    job_counts: Optional[Mapping[str, object]] = None,
) -> dict:
    events_total = 0
    events_done = 0
    rendered_variants: Dict[str, object] = {}
    latest_state: Optional[Mapping[str, object]] = None
    latest_updated_at = float("-inf")

    for label in ("unpolarized", "polarized"):
        raw_state = variant_states.get(label, {})
        total_value = raw_state.get("events_total")
        total_int = int(total_value) if isinstance(total_value, int) else None
        done_value = raw_state.get("events_done", 0)
        done_int = int(done_value) if isinstance(done_value, int) else 0
        percent = 0.0
        if total_int and total_int > 0:
            percent = 100.0 * min(max(done_int / total_int, 0.0), 1.0)
            events_total += total_int
            events_done += min(max(done_int, 0), total_int)
        rendered_variants[label] = {
            "phase": raw_state.get("phase", "pending"),
            "events_total": total_int,
            "events_done": done_int,
            "percent": percent,
            "last_progress_line": raw_state.get("last_progress_line"),
            "log": str(raw_state["log"]) if raw_state.get("log") else None,
        }
        updated_at = float(raw_state.get("updated_at", 0.0))
        if updated_at >= latest_updated_at:
            latest_updated_at = updated_at
            latest_state = raw_state

    progress_fraction = (events_done / events_total) if events_total > 0 else 0.0
    latest_log = None
    latest_line = None
    if latest_state is not None:
        if latest_state.get("log"):
            latest_log = str(latest_state["log"])
        latest_line = latest_state.get("last_progress_line")

    return {
        "tag": tag,
        "setup": setup,
        "window": window.label,
        "phase": "running-parallel",
        "variant": active_variant or "parallel",
        "active_variant": active_variant,
        "elapsed_s": time.time() - started_at,
        "events_total": events_total if events_total > 0 else None,
        "events_done": events_done,
        "variant_percent": 100.0 * progress_fraction,
        "overall_percent": 100.0 * progress_fraction,
        "last_progress_line": latest_line,
        "log": latest_log,
        "variants": rendered_variants,
        "jobs": build_job_summary(phase="running-parallel", total_jobs=0, job_counts=job_counts),
    }


def materialize_variant_sources(
    poldis_dir: Path,
    target_dir: Path,
    *,
    ipol: int,
    setup: str,
    window: CutWindow,
    scale_variation: str,
    events: int,
    seed: int,
    unpolarized_pdf: str,
    polarized_diff_pdf: str,
) -> None:
    target_dir.mkdir(parents=True, exist_ok=True)
    for filename in RUNTIME_FILES:
        source_path = poldis_dir / filename
        if not source_path.exists():
            raise FileNotFoundError(f"Missing POLDIS runtime source {source_path}")
        shutil.copy2(source_path, target_dir / filename)

    template_path = poldis_dir / "user_dijet_rivetplots.f"
    if not template_path.exists():
        raise FileNotFoundError(f"Missing POLDIS user template {template_path}")
    rendered = patch_user_source(
        template_path.read_text(),
        ipol=ipol,
        setup=setup,
        window=window,
        scale_variation=scale_variation,
        events=events,
        seed=seed,
        unpolarized_pdf=unpolarized_pdf,
        polarized_diff_pdf=polarized_diff_pdf,
    )
    rendered_path = target_dir / "user_dijet_rivetplots.f"
    if not rendered_path.exists() or rendered_path.read_text() != rendered:
        rendered_path.write_text(rendered)


def write_variant_plan(
    *,
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    label: str,
    scale_variation: str,
    specs: Sequence[ShardSpec],
) -> None:
    plan_path = variant_plan_path(base_dir, tag, setup, window, label, scale_variation)
    plan_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "variant": label,
        "scale_variation": scale_variation,
        "shard_count": len(specs),
        "shards": [
            {
                "index": spec.shard_index,
                "events": spec.events,
                "seed": spec.seed,
                "relative_dir": str(Path("shards") / f"s{spec.shard_index:03d}"),
            }
            for spec in specs
        ],
    }
    plan_path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def load_variant_plan(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    label: str,
    scale_variation: str,
) -> Optional[dict]:
    plan_path = variant_plan_path(base_dir, tag, setup, window, label, scale_variation)
    if not plan_path.exists():
        return None
    return json.loads(plan_path.read_text())


def run_logged(cmd: Sequence[str], cwd: Path, log_path: Path, dry_run: bool) -> None:
    if dry_run:
        print(shlex.join(list(cmd)))
        return
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.run(cmd, cwd=cwd, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed in {cwd}: {shlex.join(list(cmd))}\nSee {log_path}")


def run_poldis_with_progress(
    *,
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    run_dir: Path,
    variant_label: str,
    events: int,
    run_index: int,
    run_count: int,
    started_at: float,
    dry_run: bool,
    write_monitor: bool = True,
    progress_callback: Optional[Callable[[Mapping[str, object]], None]] = None,
) -> None:
    cmd = ["./poldis.x"]
    log_path = run_dir / "run.log"
    if dry_run:
        print(shlex.join(cmd))
        return

    log_path.parent.mkdir(parents=True, exist_ok=True)
    payload = build_monitor_payload(
        tag=tag,
        setup=setup,
        window=window,
        phase=f"running-{variant_label}",
        variant=variant_label,
        started_at=started_at,
        run_index=run_index,
        run_count=run_count,
        events_total=events,
        events_done=0,
        log_path=log_path,
    )
    if write_monitor:
        write_monitor_files(base_dir, tag, setup, window, payload)
    if progress_callback is not None:
        progress_callback(payload)

    with log_path.open("w") as handle:
        handle.write(f"$ {' '.join(cmd)}\n")
        handle.flush()
        proc = subprocess.Popen(
            cmd,
            cwd=run_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        last_events = 0
        last_line = None
        for line in proc.stdout:
            handle.write(line)
            handle.flush()
            stripped = line.rstrip()
            match = PROGRESS_RE.match(stripped)
            if match:
                last_events = int(match.group("events"))
                last_line = stripped
                payload = build_monitor_payload(
                    tag=tag,
                    setup=setup,
                    window=window,
                    phase=f"running-{variant_label}",
                    variant=variant_label,
                    started_at=started_at,
                    run_index=run_index,
                    run_count=run_count,
                    events_total=events,
                    events_done=last_events,
                    last_progress_line=last_line,
                    log_path=log_path,
                )
                if write_monitor:
                    write_monitor_files(base_dir, tag, setup, window, payload)
                if progress_callback is not None:
                    progress_callback(payload)
        returncode = proc.wait()

    if returncode != 0:
        raise RuntimeError(f"Command failed in {run_dir}: {shlex.join(cmd)}\nSee {log_path}")

    payload = build_monitor_payload(
        tag=tag,
        setup=setup,
        window=window,
        phase=f"completed-{variant_label}",
        variant=variant_label,
        started_at=started_at,
        run_index=run_index + 1,
        run_count=run_count,
        events_total=events,
        events_done=events,
        last_progress_line=last_line,
        log_path=log_path,
    )
    if write_monitor:
        write_monitor_files(base_dir, tag, setup, window, payload)
    if progress_callback is not None:
        progress_callback(payload)


def ensure_compiled(poldis_dir: Path, run_dir: Path, dry_run: bool) -> None:
    cmd = compile_command(poldis_dir)
    binary = run_dir / "poldis.x"
    if not dry_run and binary.exists():
        source_mtime = max((run_dir / filename).stat().st_mtime for filename in (*RUNTIME_FILES, "user_dijet_rivetplots.f"))
        if binary.stat().st_mtime >= source_mtime:
            return
    run_logged(cmd, cwd=run_dir, log_path=run_dir / "compile.log", dry_run=dry_run)


def ensure_ran(run_dir: Path, dry_run: bool) -> None:
    run_logged(["./poldis.x"], cwd=run_dir, log_path=run_dir / "run.log", dry_run=dry_run)


def process_shard(
    *,
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    poldis_dir: Path,
    spec: ShardSpec,
    pdfs: Mapping[str, str],
    started_at: float,
    dry_run: bool,
    progress_callback: Optional[Callable[[str, int, Mapping[str, object]], None]] = None,
) -> str:
    run_dir = variant_shard_dir(base_dir, tag, setup, window, spec.label, spec.scale_variation, spec.shard_index)
    materialize_variant_sources(
        poldis_dir,
        run_dir,
        ipol=spec.ipol,
        setup=setup,
        window=window,
        scale_variation=spec.scale_variation,
        events=spec.events,
        seed=spec.seed,
        unpolarized_pdf=pdfs["unpolarized"],
        polarized_diff_pdf=pdfs["polarized_diff"],
    )
    variant_label = scale_variation_task_label(spec.scale_variation, spec.label)
    shard_label = f"{variant_label}-s{spec.shard_index:03d}"
    compile_payload = build_monitor_payload(
        tag=tag,
        setup=setup,
        window=window,
        phase=f"compiling-{shard_label}",
        variant=shard_label,
        started_at=started_at,
        run_index=0,
        run_count=1,
        log_path=run_dir / "compile.log",
    )
    if not dry_run and progress_callback is not None:
        progress_callback(variant_label, spec.shard_index, compile_payload)
    ensure_compiled(poldis_dir, run_dir, dry_run=dry_run)
    run_poldis_with_progress(
        base_dir=base_dir,
        tag=tag,
        setup=setup,
        window=window,
        run_dir=run_dir,
        variant_label=shard_label,
        events=spec.events,
        run_index=0,
        run_count=1,
        started_at=started_at,
        dry_run=dry_run,
        write_monitor=False,
        progress_callback=(
            None
            if progress_callback is None
            else (
                lambda payload, *, _label=variant_label, _index=spec.shard_index: progress_callback(_label, _index, payload)
            )
        ),
    )
    return shard_label


def load_variant_payload(
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    *,
    label: str,
    scale_variation: str,
    ipol: int,
) -> dict:
    run_dir = scale_variation_dir(base_dir, tag, setup, window, label, scale_variation)
    plan_path = run_dir / "shards.json"
    if plan_path.exists():
        plan = json.loads(plan_path.read_text())
        shard_payloads = []
        measurements_by_order: Dict[str, List[Measurement]] = {"LO": [], "NLO": [], "NNLO": []}
        shard_events: List[int] = []
        for shard_info in plan.get("shards", []):
            shard_dir = run_dir / str(shard_info["relative_dir"])
            log_path = shard_dir / "run.log"
            if not log_path.exists():
                raise FileNotFoundError(f"Missing POLDIS shard log {log_path}")
            totals = parse_poldis_totals(log_path.read_text(), context=str(log_path))
            top_path = shard_dir / expected_top_name(setup, ipol)
            shard_payload = {
                "index": int(shard_info["index"]),
                "events": int(shard_info["events"]),
                "seed": int(shard_info["seed"]),
                "log": str(log_path),
                "top": str(top_path) if top_path.exists() else None,
                "totals": {order: measurement_payload(totals[order]) for order in ("LO", "NLO", "NNLO")},
            }
            shard_payloads.append(shard_payload)
            shard_events.append(int(shard_info["events"]))
            for order in ("LO", "NLO", "NNLO"):
                measurements_by_order[order].append(totals[order])
        combined_totals = {
            order: measurement_payload(combine_measurements(measurements_by_order[order], shard_events))
            for order in ("LO", "NLO", "NNLO")
        }
        payload = {
            "log": str(run_dir / "shards"),
            "log_dir": str(run_dir / "shards"),
            "top": None,
            "shard_count": len(shard_payloads),
            "scale_variation": scale_variation,
            "totals": combined_totals,
            "shards": shard_payloads,
        }
        if len(shard_payloads) == 1:
            payload["log"] = shard_payloads[0]["log"]
            payload["top"] = shard_payloads[0]["top"]
        return payload

    log_path = run_dir / "run.log"
    if not log_path.exists():
        raise FileNotFoundError(f"Missing POLDIS run log {log_path}")
    totals = parse_poldis_totals(log_path.read_text(), context=str(log_path))
    top_path = run_dir / expected_top_name(setup, ipol)
    return {
        "log": str(log_path),
        "top": str(top_path) if top_path.exists() else None,
        "scale_variation": scale_variation,
        "totals": {order: measurement_payload(totals[order]) for order in ("LO", "NLO", "NNLO")},
    }


def build_broad_comparisons(setup: str, payload: dict) -> dict:
    comparisons: dict = {"unpolarized": {}, "polarized": {}}
    for order in ("LO", "NLO", "NNLO"):
        unpol_meas = measurement_from_payload(payload["unpolarized"][order])
        pol_meas = measurement_from_payload(payload["polarized"][order])
        ref_unpol = POLDIS_UNPOL_REFS[setup][order]
        ref_pol = POLDIS_POL_REFS[setup][order]
        comparisons["unpolarized"][order] = measurement_payload(sub_measurement(unpol_meas, ref_unpol))
        comparisons["polarized"][order] = measurement_payload(sub_measurement(pol_meas, ref_pol))
    return comparisons


def build_scale_variation_reference_objects(unpol_payload: Mapping[str, object], pol_payload: Mapping[str, object]) -> Dict[str, object]:
    from poldis_top_to_yoda import combine_ref_object_sets, convert_topdrawer_to_ref_objects

    unpol_shards = unpol_payload.get("shards")
    pol_shards = pol_payload.get("shards")
    if isinstance(unpol_shards, list) and isinstance(pol_shards, list):
        unpol_by_index = {
            int(shard["index"]): shard for shard in unpol_shards if isinstance(shard, dict) and shard.get("top")
        }
        pol_by_index = {
            int(shard["index"]): shard for shard in pol_shards if isinstance(shard, dict) and shard.get("top")
        }
        shard_indices = sorted(set(unpol_by_index) & set(pol_by_index))
        if not shard_indices:
            raise RuntimeError("No matching unpolarized/polarized shard Topdrawer pairs found.")

        ref_sets: List[Dict[str, object]] = []
        generated_events: List[int] = []
        for shard_index in shard_indices:
            unpol_shard = unpol_by_index[shard_index]
            pol_shard = pol_by_index[shard_index]
            unpol_events = int(unpol_shard.get("events", 0) or 0)
            pol_events = int(pol_shard.get("events", 0) or 0)
            if unpol_events <= 0 or pol_events <= 0:
                raise RuntimeError(f"Invalid shard event counts for shard {shard_index}")
            if unpol_events != pol_events:
                raise RuntimeError(
                    f"Mismatched unpolarized/polarized shard event counts for shard {shard_index}: "
                    f"{unpol_events} vs {pol_events}"
                )
            ref_sets.append(
                convert_topdrawer_to_ref_objects(
                    unpol=str(unpol_shard["top"]),
                    pol=str(pol_shard["top"]),
                )
            )
            generated_events.append(unpol_events)
        return combine_ref_object_sets(ref_sets, generated_events=generated_events)

    unpol_top = str(unpol_payload.get("top") or "")
    pol_top = str(pol_payload.get("top") or "")
    if not unpol_top or not pol_top:
        raise RuntimeError("Missing nominal unpolarized/polarized Topdrawer files for reference conversion.")
    return convert_topdrawer_to_ref_objects(unpol=unpol_top, pol=pol_top)


def build_reference_yodas(
    *,
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    variation_runs: Mapping[str, Mapping[str, Mapping[str, object]]],
    error_mode: str,
) -> tuple[str, Dict[str, str], Dict[str, str]]:
    from poldis_top_to_yoda import build_ref_object_error_mode, write_yoda_gz

    variation_ref_objects: Dict[str, Dict[str, object]] = {}
    variation_outputs: Dict[str, str] = {}
    for scale_variation, runs in variation_runs.items():
        ref_objects = build_scale_variation_reference_objects(runs["unpolarized"], runs["polarized"])
        variation_ref_objects[scale_variation] = ref_objects
        output_path = reference_yoda_variation_path(base_dir, tag, setup, window, scale_variation)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        write_yoda_gz(ref_objects, str(output_path))
        variation_outputs[scale_variation] = str(output_path)

    nominal_ref_objects = variation_ref_objects["nominal"]
    scale_down_ref_objects = variation_ref_objects.get("ScaleFactorDown")
    scale_up_ref_objects = variation_ref_objects.get("ScaleFactorUp")
    error_mode_outputs: Dict[str, str] = {}
    for candidate_mode in REFERENCE_ERROR_MODES:
        if candidate_mode != "stat" and (scale_down_ref_objects is None or scale_up_ref_objects is None):
            continue
        ref_objects = build_ref_object_error_mode(
            nominal_ref_objects,
            error_mode=candidate_mode,
            scale_down_ref_objects=scale_down_ref_objects,
            scale_up_ref_objects=scale_up_ref_objects,
        )
        output_path = reference_yoda_error_mode_path(base_dir, tag, setup, window, candidate_mode)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        write_yoda_gz(ref_objects, str(output_path))
        error_mode_outputs[candidate_mode] = str(output_path)

    if error_mode not in error_mode_outputs:
        raise RuntimeError(
            f"Requested reference error mode {error_mode!r} is unavailable for setup {setup} window {window.label}"
        )

    selected_path = reference_yoda_path(base_dir, tag, setup, window)
    selected_ref_objects = build_ref_object_error_mode(
        nominal_ref_objects,
        error_mode=error_mode,
        scale_down_ref_objects=scale_down_ref_objects,
        scale_up_ref_objects=scale_up_ref_objects,
    )
    selected_path.parent.mkdir(parents=True, exist_ok=True)
    write_yoda_gz(selected_ref_objects, str(selected_path))
    return str(selected_path), variation_outputs, error_mode_outputs


def build_reference_payload(
    *,
    base_dir: Path,
    tag: str,
    setup: str,
    window: CutWindow,
    events: int,
    shards: int,
    seed_base: int,
    pdf_profile: str,
    pdfs: Mapping[str, str],
    scale_variation_names: Sequence[str],
    error_mode: str,
) -> dict:
    variation_runs = {
        scale_variation: {
            "unpolarized": load_variant_payload(
                base_dir,
                tag,
                setup,
                window,
                label="unpolarized",
                scale_variation=scale_variation,
                ipol=0,
            ),
            "polarized": load_variant_payload(
                base_dir,
                tag,
                setup,
                window,
                label="polarized",
                scale_variation=scale_variation,
                ipol=1,
            ),
        }
        for scale_variation in scale_variation_names
    }
    unpol_payload = variation_runs["nominal"]["unpolarized"]
    pol_payload = variation_runs["nominal"]["polarized"]
    reference_yoda, reference_yoda_variations, reference_yoda_error_modes = build_reference_yodas(
        base_dir=base_dir,
        tag=tag,
        setup=setup,
        window=window,
        variation_runs=variation_runs,
        error_mode=error_mode,
    )
    payload = {
        "tag": tag,
        "mode": "poldis-window-reference",
        "setup": setup,
        "window": asdict(window),
        "events": events,
        "shards": shards,
        "seed_base": seed_base,
        "pdf_profile": pdf_profile,
        "scale_variations": list(scale_variation_names),
        "reference_yoda_error_mode": error_mode,
        "pdfs": {
            "unpolarized": pdfs["unpolarized"],
            "polarized_diff": pdfs["polarized_diff"],
        },
        "work_dir": str(work_dir(base_dir, tag, setup, window)),
        "reference_yoda": reference_yoda,
        "reference_yoda_variations": reference_yoda_variations,
        "reference_yoda_error_modes": reference_yoda_error_modes,
        "runs": {
            "unpolarized": unpol_payload,
            "polarized": pol_payload,
        },
        "scale_variation_runs": variation_runs,
        "unpolarized": unpol_payload["totals"],
        "polarized": pol_payload["totals"],
    }
    if (
        window.label == "broad"
        and uses_builtin_broad_reference_choice(setup, pdfs)
        and setup in POLDIS_POL_REFS
        and setup in POLDIS_UNPOL_REFS
    ):
        payload["comparisons_vs_builtin_hybrid_broad"] = build_broad_comparisons(setup, payload)
    return payload


def render_reference_text(payload: dict) -> str:
    scale_variations = payload.get("scale_variations", ["nominal"])
    lines = [
        "POLDIS window reference",
        "=======================",
        f"Tag: {payload['tag']}",
        f"Setup: {payload['setup']}",
        f"Window: {payload['window']['label']}",
        f"Total events per variant: {payload['events']}",
        f"Shards per variant: {payload.get('shards', 1)}",
        f"PDF profile: {payload['pdf_profile']}",
        f"Unpolarized PDF: {payload['pdfs']['unpolarized']}",
        f"Polarized diff PDF: {payload['pdfs']['polarized_diff']}",
        f"Scale variations: {', '.join(str(name) for name in scale_variations)}",
        f"Selected reference error mode: {payload.get('reference_yoda_error_mode', DEFAULT_REFERENCE_ERROR_MODE)}",
        f"Reference YODA: {payload['reference_yoda']}",
    ]
    if isinstance(payload.get("reference_yoda_error_modes"), dict):
        lines.append("Reference YODA error modes:")
        for name, path in payload["reference_yoda_error_modes"].items():
            lines.append(f"  {name}: {path}")
    if isinstance(payload.get("reference_yoda_variations"), dict) and len(payload["reference_yoda_variations"]) > 1:
        lines.append("Reference YODA variations:")
        for name, path in payload["reference_yoda_variations"].items():
            lines.append(f"  {name}: {path}")
    if int(payload.get("shards", 1)) > 1:
        lines.append("Shard combination: event-weighted mean with independent shard errors added in quadrature")
    lines.extend(
        [
            "",
            "Unpolarized totals",
            "------------------",
            f"LO:   {fmt_measurement(payload['unpolarized']['LO'])}",
            f"NLO:  {fmt_measurement(payload['unpolarized']['NLO'])}",
            f"NNLO: {fmt_measurement(payload['unpolarized']['NNLO'])}",
            f"log:  {payload['runs']['unpolarized']['log']}",
            f"top:  {payload['runs']['unpolarized']['top'] or 'multiple/missing'}",
            "",
            "Polarized totals",
            "----------------",
            f"LO:   {fmt_measurement(payload['polarized']['LO'])}",
            f"NLO:  {fmt_measurement(payload['polarized']['NLO'])}",
            f"NNLO: {fmt_measurement(payload['polarized']['NNLO'])}",
            f"log:  {payload['runs']['polarized']['log']}",
            f"top:  {payload['runs']['polarized']['top'] or 'multiple/missing'}",
        ]
    )
    comparisons = payload.get("comparisons_vs_builtin_hybrid_broad")
    if comparisons:
        lines.extend(["", "Difference vs built-in broad references", "-------------------------------------"])
        for order in ("LO", "NLO", "NNLO"):
            lines.append(
                f"Unpolarized {order} - built-in: {fmt_measurement(comparisons['unpolarized'][order])}"
            )
            lines.append(
                f"Polarized {order} - built-in:   {fmt_measurement(comparisons['polarized'][order])}"
            )
    return "\n".join(lines).rstrip() + "\n"


def write_reference(base_dir: Path, tag: str, setup: str, window: CutWindow, payload: dict) -> None:
    out_dir = work_dir(base_dir, tag, setup, window)
    out_dir.mkdir(parents=True, exist_ok=True)
    text = render_reference_text(payload)
    reference_txt_path(base_dir, tag, setup, window).write_text(text)
    reference_json_path(base_dir, tag, setup, window).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(text, end="")
    print(f"Wrote text reference: {reference_txt_path(base_dir, tag, setup, window)}")
    print(f"Wrote JSON reference: {reference_json_path(base_dir, tag, setup, window)}")
    print(
        f"Wrote selected reference YODA ({payload.get('reference_yoda_error_mode', DEFAULT_REFERENCE_ERROR_MODE)}): "
        f"{payload['reference_yoda']}"
    )
    for name, path in payload.get("reference_yoda_error_modes", {}).items():
        print(f"Wrote {name} error-mode reference YODA: {path}")
    for name, path in payload.get("reference_yoda_variations", {}).items():
        print(f"Wrote {name} reference YODA: {path}")


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    base_dir = args.base_dir.resolve()
    poldis_dir = args.poldis_dir.resolve()
    window = WINDOWS[args.window]
    pdfs = resolve_pdf_choice(args.pdf_profile, args.unpolarized_pdf, args.polarized_diff_pdf)
    started_at = time.time()
    scale_variation_names = enabled_scale_variations(args.scale_variations)
    if args.error_mode != "stat" and "ScaleFactorDown" not in scale_variation_names:
        raise ValueError("--error-mode scale/stat+scale requires --scale-variations so the POLDIS envelope exists.")
    variants = [
        {
            "scale_variation": scale_variation,
            "label": label,
            "ipol": ipol,
            "task_label": scale_variation_task_label(scale_variation, label),
        }
        for scale_variation in scale_variation_names
        for label, ipol in (("unpolarized", 0), ("polarized", 1))
    ]
    variant_specs = {
        variant["task_label"]: build_variant_shards(
            label=str(variant["label"]),
            ipol=int(variant["ipol"]),
            scale_variation=str(variant["scale_variation"]),
            total_events=args.events,
            shard_count=max(1, args.shards),
            seed_base=args.seed_base,
        )
        for variant in variants
    }

    if args.collect_only:
        write_monitor_files(
            base_dir,
            args.tag,
            args.setup,
            window,
            build_monitor_payload(
                tag=args.tag,
                setup=args.setup,
                window=window,
                phase="collecting",
                variant="all",
                started_at=started_at,
                run_index=0,
                run_count=len(variants),
                log_path=reference_json_path(base_dir, args.tag, args.setup, window),
            ),
        )
    else:
        if not args.dry_run:
            for variant in variants:
                specs = variant_specs[str(variant["task_label"])]
                write_variant_plan(
                    base_dir=base_dir,
                    tag=args.tag,
                    setup=args.setup,
                    window=window,
                    label=str(variant["label"]),
                    scale_variation=str(variant["scale_variation"]),
                    specs=specs,
                )

            shard_states: Dict[str, Dict[int, Dict[str, object]]] = {
                str(variant["task_label"]): {
                    spec.shard_index: {
                        "phase": "pending",
                        "events_total": spec.events,
                        "events_done": 0,
                        "last_progress_line": None,
                        "log": variant_shard_dir(
                            base_dir,
                            args.tag,
                            args.setup,
                            window,
                            str(variant["label"]),
                            str(variant["scale_variation"]),
                            spec.shard_index,
                        )
                        / "run.log",
                        "updated_at": started_at,
                    }
                    for spec in specs
                }
                for variant in variants
                for specs in [variant_specs[str(variant["task_label"])]]
            }
            monitor_lock = threading.Lock()

            def aggregate_variant_state(task_label: str) -> Dict[str, object]:
                states = shard_states[task_label]
                total_events = 0
                done_events = 0
                latest_state: Optional[Mapping[str, object]] = None
                latest_updated_at = float("-inf")
                has_running = False
                has_compiling = False
                all_completed = True
                for state in states.values():
                    events_total = int(state.get("events_total", 0) or 0)
                    events_done = int(state.get("events_done", 0) or 0)
                    total_events += events_total
                    done_events += min(max(events_done, 0), events_total)
                    phase = str(state.get("phase", "pending"))
                    if phase.startswith("running-"):
                        has_running = True
                    elif phase.startswith("compiling-"):
                        has_compiling = True
                    if not phase.startswith("completed-"):
                        all_completed = False
                    updated_at = float(state.get("updated_at", 0.0))
                    if updated_at >= latest_updated_at:
                        latest_updated_at = updated_at
                        latest_state = state
                if all_completed:
                    phase = f"completed-{task_label}"
                elif has_running:
                    phase = f"running-{task_label}"
                elif has_compiling:
                    phase = f"compiling-{task_label}"
                else:
                    phase = "pending"
                return {
                    "phase": phase,
                    "events_total": total_events,
                    "events_done": done_events,
                    "last_progress_line": latest_state.get("last_progress_line") if latest_state else None,
                    "log": latest_state.get("log") if latest_state else None,
                    "updated_at": latest_updated_at if latest_state is not None else started_at,
                }

            def flush_monitor(active_variant: Optional[str] = None) -> None:
                running_jobs = 0
                compiling_jobs = 0
                pending_jobs = 0
                completed_jobs = 0
                for states in shard_states.values():
                    for state in states.values():
                        phase = str(state.get("phase", "pending"))
                        if phase.startswith("running-"):
                            running_jobs += 1
                        elif phase.startswith("compiling-"):
                            compiling_jobs += 1
                        elif phase.startswith("completed-"):
                            completed_jobs += 1
                        else:
                            pending_jobs += 1
                write_monitor_files(
                    base_dir,
                    args.tag,
                    args.setup,
                    window,
                    build_parallel_monitor_payload(
                        tag=args.tag,
                        setup=args.setup,
                        window=window,
                        started_at=started_at,
                        variant_states={
                            str(variant["task_label"]): aggregate_variant_state(str(variant["task_label"]))
                            for variant in variants
                        },
                        active_variant=active_variant,
                        job_counts={
                            "running": running_jobs,
                            "compiling": compiling_jobs,
                            "pending": pending_jobs,
                            "completed": completed_jobs,
                            "total": running_jobs + compiling_jobs + pending_jobs + completed_jobs,
                        },
                    ),
                )

            def update_parallel_monitor(task_label: str, shard_index: int, payload: Mapping[str, object]) -> None:
                with monitor_lock:
                    shard_states[task_label][shard_index] = {
                        "phase": payload.get("phase", "pending"),
                        "events_total": payload.get("events_total"),
                        "events_done": payload.get("events_done", 0),
                        "last_progress_line": payload.get("last_progress_line"),
                        "log": payload.get("log"),
                        "updated_at": time.time(),
                    }
                    flush_monitor(active_variant=f"{task_label}-s{shard_index:03d}")

            flush_monitor()
        else:
            update_parallel_monitor = None

        shard_specs = [spec for specs in variant_specs.values() for spec in specs]
        worker_count = max(1, min(args.jobs, len(shard_specs)))
        if args.dry_run or worker_count <= 1:
            for spec in shard_specs:
                process_shard(
                    base_dir=base_dir,
                    tag=args.tag,
                    setup=args.setup,
                    window=window,
                    poldis_dir=poldis_dir,
                    spec=spec,
                    pdfs=pdfs,
                    started_at=started_at,
                    dry_run=args.dry_run,
                    progress_callback=update_parallel_monitor,
                )
        else:
            with concurrent.futures.ThreadPoolExecutor(max_workers=worker_count) as executor:
                futures = [
                    executor.submit(
                        process_shard,
                        base_dir=base_dir,
                        tag=args.tag,
                        setup=args.setup,
                        window=window,
                        poldis_dir=poldis_dir,
                        spec=spec,
                        pdfs=pdfs,
                        started_at=started_at,
                        dry_run=args.dry_run,
                        progress_callback=update_parallel_monitor,
                    )
                    for spec in shard_specs
                ]
                for future in concurrent.futures.as_completed(futures):
                    future.result()

    if args.dry_run:
        return 0

    payload = build_reference_payload(
        base_dir=base_dir,
        tag=args.tag,
        setup=args.setup,
        window=window,
        events=args.events,
        shards=max(1, args.shards),
        seed_base=args.seed_base,
        pdf_profile=args.pdf_profile,
        pdfs=pdfs,
        scale_variation_names=scale_variation_names,
        error_mode=args.error_mode,
    )
    write_reference(base_dir, args.tag, args.setup, window, payload)
    write_monitor_files(
        base_dir,
        args.tag,
        args.setup,
        window,
        build_monitor_payload(
            tag=args.tag,
            setup=args.setup,
            window=window,
            phase="complete",
            variant="all",
            started_at=started_at,
            run_index=len(variants),
            run_count=len(variants),
            log_path=reference_txt_path(base_dir, args.tag, args.setup, window),
        ),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
