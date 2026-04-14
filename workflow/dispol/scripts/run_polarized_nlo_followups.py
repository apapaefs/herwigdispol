#!/usr/bin/env python3.10
"""
Run focused polarized NLO follow-up studies for the DIS validation campaign.

Modes:
  - gamma-interior: interior-window polarized GAMMA NLO study
  - gamma-power-scan: broad-window polarized GAMMA SamplingPower scan
  - z-termdiag: Z NLO term-diagnostic runs and extraction
  - all: run the three studies in sequence

The script mirrors the operational style of run_z_lo_sigma0_check.py:
  * generate temporary .in files under campaigns/<tag>/<mode>/generated-inputs/
  * prepare missing or stale .run files with `Herwig read`
  * launch Herwig shards
  * maintain monitor/status.txt and monitor/status.json
  * collect mode-specific text/JSON reports
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import json
import math
import re
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence

from extract_dis_out_results import (
    Measurement,
    combine_measurements,
    load_run,
    parse_generated_event_count,
    parse_measurement,
)


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_TAG = "plain48-pnlo-followups"
DEFAULT_JOBS = 32
DEFAULT_SEED_BASE = 700_000

GAMMA_INTERIOR_MODE = "gamma-interior"
GAMMA_POWER_SCAN_MODE = "gamma-power-scan"
Z_TERMDIAG_MODE = "z-termdiag"
ALL_MODE = "all"
VALID_MODES = (GAMMA_INTERIOR_MODE, GAMMA_POWER_SCAN_MODE, Z_TERMDIAG_MODE, ALL_MODE)

GAMMA_CONTROL_TAG = "plain47"
GAMMA_INTERIOR_SUFFIX = "GNLOINT"
Z_TERMDIAG_SUFFIX = "ZTDIAG"

DEFAULT_GAMMA_INTERIOR_SHARDS = 100
DEFAULT_GAMMA_INTERIOR_EVENTS = 1_000_000
DEFAULT_GAMMA_POWER_SCAN_SHARDS = 20
DEFAULT_GAMMA_POWER_SCAN_EVENTS = 1_000_000
DEFAULT_Z_TERMDIAG_SHARDS = 1
DEFAULT_Z_TERMDIAG_EVENTS = 1_000_000
DEFAULT_POWER_SCAN_POWERS = "0.4,0.8"

SAMPLING_POWER_OBJECTS = ("MEDISNC", "MEDISNCPol", "PowhegMEDISNC", "PowhegMEDISNCPol")
FINITE_WIDTH_OBJECTS = SAMPLING_POWER_OBJECTS
GAMMA_HELICITIES = ("PP", "PM")
NLO_PIECES = ("POSNLO", "NEGNLO")


@dataclass(frozen=True)
class CutWindow:
    label: str
    q2_min: float
    q2_max: float
    y_min: float
    y_max: float


@dataclass(frozen=True)
class LogicalRun:
    mode: str
    key: str
    source_card: str
    requested_name: str
    generated_stem: str
    generated_input_rel: Path
    run_file: str
    piece: str
    helicity: str
    power: Optional[float] = None


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


BROAD_GAMMA_WINDOW = CutWindow(label="broad", q2_min=49.0, q2_max=2500.0, y_min=0.2, y_max=0.6)
INTERIOR_GAMMA_WINDOW = CutWindow(label="interior", q2_min=100.0, q2_max=1000.0, y_min=0.3, y_max=0.5)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=VALID_MODES, help="Study mode to run.")
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used in Herwig shard tags.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing DISPOL cards and outputs.",
    )
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Maximum concurrent Herwig processes.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without executing Herwig.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip Herwig read/run and only collect reports from existing outputs.",
    )
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE, help="Base RNG seed.")

    parser.add_argument(
        "--control-tag",
        default=GAMMA_CONTROL_TAG,
        help="Existing broad GAMMA control campaign tag reused for gamma-interior.",
    )
    parser.add_argument(
        "--poldis-reference",
        type=Path,
        help="Optional JSON file from run_poldis_gamma_window_reference.py for gamma-interior comparison.",
    )
    parser.add_argument(
        "--gamma-interior-shards",
        type=int,
        default=DEFAULT_GAMMA_INTERIOR_SHARDS,
        help="Shards per logical run for gamma-interior.",
    )
    parser.add_argument(
        "--gamma-interior-events-per-shard",
        type=int,
        default=DEFAULT_GAMMA_INTERIOR_EVENTS,
        help="Events per shard for gamma-interior.",
    )

    parser.add_argument(
        "--powers",
        default=DEFAULT_POWER_SCAN_POWERS,
        help="Comma-separated SamplingPower values to run for gamma-power-scan.",
    )
    parser.add_argument(
        "--include-existing",
        action="append",
        default=[],
        help="Overlay an existing campaign in the scan as power=tag, for example 0.2=plain45.",
    )
    parser.add_argument(
        "--gamma-power-scan-shards",
        type=int,
        default=DEFAULT_GAMMA_POWER_SCAN_SHARDS,
        help="Shards per logical run for gamma-power-scan.",
    )
    parser.add_argument(
        "--gamma-power-scan-events-per-shard",
        type=int,
        default=DEFAULT_GAMMA_POWER_SCAN_EVENTS,
        help="Events per shard for gamma-power-scan.",
    )

    parser.add_argument(
        "--z-termdiag-shards",
        type=int,
        default=DEFAULT_Z_TERMDIAG_SHARDS,
        help="Shards per logical run for z-termdiag.",
    )
    parser.add_argument(
        "--z-termdiag-events-per-shard",
        type=int,
        default=DEFAULT_Z_TERMDIAG_EVENTS,
        help="Events per shard for z-termdiag.",
    )
    parser.add_argument(
        "--nlo-audit-initial-samples",
        type=int,
        help="Optional override for NLOAuditInitialSamples in generated Z TERMDIAG cards.",
    )
    parser.add_argument(
        "--nlo-audit-sample-period",
        type=int,
        help="Optional override for NLOAuditSamplePeriod in generated Z TERMDIAG cards.",
    )
    parser.add_argument(
        "--nlo-term-diagnostic-period",
        type=int,
        help="Optional override for NLOTermDiagnosticPeriod in generated Z TERMDIAG cards.",
    )
    return parser


def build_gamma_requested_name(piece: str, helicity: str) -> str:
    return f"DIS-POL-POWHEG_{helicity}-{piece}-GAMMA.out"


def build_gamma_source_card(piece: str, helicity: str) -> str:
    return build_gamma_requested_name(piece, helicity).replace(".out", ".in")


def build_z_termdiag_requested_name(piece: str, helicity: str) -> str:
    return f"DIS-POL-POWHEG_{helicity}-{piece}-Z-TERMDIAG.out"


def build_z_termdiag_source_card(piece: str, helicity: str) -> str:
    return build_z_termdiag_requested_name(piece, helicity).replace(".out", ".in")


def stem_without_suffix(filename: str) -> str:
    return filename[:-3] if filename.endswith(".in") else filename


def campaign_mode_dir(base_dir: Path, tag: str, mode: str) -> Path:
    return base_dir / "campaigns" / tag / mode


def generated_input_dir(base_dir: Path, tag: str, mode: str) -> Path:
    return campaign_mode_dir(base_dir, tag, mode) / "generated-inputs"


def monitor_dir(base_dir: Path, tag: str, mode: str) -> Path:
    return campaign_mode_dir(base_dir, tag, mode) / "monitor"


def launcher_log_dir(base_dir: Path, tag: str, mode: str) -> Path:
    return campaign_mode_dir(base_dir, tag, mode) / "launcher-logs"


def status_json_path(base_dir: Path, tag: str, mode: str) -> Path:
    return monitor_dir(base_dir, tag, mode) / "status.json"


def status_txt_path(base_dir: Path, tag: str, mode: str) -> Path:
    return monitor_dir(base_dir, tag, mode) / "status.txt"


def results_txt_path(base_dir: Path, tag: str, mode: str) -> Path:
    return campaign_mode_dir(base_dir, tag, mode) / "results.txt"


def results_json_path(base_dir: Path, tag: str, mode: str) -> Path:
    return campaign_mode_dir(base_dir, tag, mode) / "results.json"


def launcher_log_path(base_dir: Path, tag: str, mode: str, spec: ShardSpec) -> Path:
    return launcher_log_dir(base_dir, tag, mode) / f"{spec.logical_run.generated_stem}-{spec.tag}.launcher.log"


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def build_monitor_payload(
    base_dir: Path,
    tag: str,
    mode: str,
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
    mean_duration_s = (
        sum(result.duration_s for result in results) / completed if completed else None
    )
    eta_s = mean_duration_s * pending if mean_duration_s is not None else None

    per_run: Dict[str, dict] = {}
    labels = [spec.logical_run.key for spec in shards]
    for label in sorted(set(labels)):
        run_specs = [spec for spec in shards if spec.logical_run.key == label]
        run_results = [result for result in results if result.spec.logical_run.key == label]
        run_completed = len(run_results)
        run_failures = sum(1 for result in run_results if result.returncode != 0)
        per_run[label] = {
            "total": len(run_specs),
            "completed": run_completed,
            "pending": len(run_specs) - run_completed,
            "failures": run_failures,
        }

    payload = {
        "tag": tag,
        "mode": mode,
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
            "status_json": str(status_json_path(base_dir, tag, mode)),
            "status_txt": str(status_txt_path(base_dir, tag, mode)),
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


def write_monitor_files(base_dir: Path, tag: str, mode: str, payload: Mapping[str, object]) -> None:
    mon_dir = monitor_dir(base_dir, tag, mode)
    mon_dir.mkdir(parents=True, exist_ok=True)
    status_json_path(base_dir, tag, mode).write_text(json.dumps(payload, indent=2, sort_keys=True))
    status_txt_path(base_dir, tag, mode).write_text(render_monitor_text(payload))


def run_one_shard(base_dir: Path, tag: str, mode: str, spec: ShardSpec) -> ShardResult:
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
    log_path = launcher_log_path(base_dir, tag, mode, spec)
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


def run_shards(base_dir: Path, tag: str, mode: str, shards: Sequence[ShardSpec], jobs: int) -> List[ShardResult]:
    total = len(shards)
    started = time.time()
    results: List[ShardResult] = []
    write_monitor_files(
        base_dir,
        tag,
        mode,
        build_monitor_payload(base_dir, tag, mode, shards, results, started, phase="running"),
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        future_map = {executor.submit(run_one_shard, base_dir, tag, mode, spec): spec for spec in shards}
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            phase = "failed" if any(item.returncode != 0 for item in results) else "running"
            if len(results) == total and phase != "failed":
                phase = "complete"
            write_monitor_files(
                base_dir,
                tag,
                mode,
                build_monitor_payload(
                    base_dir,
                    tag,
                    mode,
                    shards,
                    results,
                    started,
                    last_result=result,
                    phase=phase,
                ),
            )
            if len(results) == 1 or len(results) % 10 == 0 or len(results) == total or result.returncode != 0:
                print(
                    f"[progress:{mode}] {len(results)}/{total} shards complete, "
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
        raise RuntimeError(f"{len(bad)} shard(s) failed in {mode}:\n{summary}")
    return results


def build_shards(
    runs: Sequence[LogicalRun],
    tag: str,
    shards: int,
    events_per_shard: int,
    seed_base: int,
) -> List[ShardSpec]:
    specs: List[ShardSpec] = []
    seed = seed_base
    for run in runs:
        for shard_index in range(1, shards + 1):
            specs.append(
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
    return specs


def fmt_power_token(power: float) -> str:
    return f"{power:.6g}"


def power_suffix(power: float) -> str:
    code = int(round(power * 100.0))
    return f"SP{code:03d}"


def parse_power_list(text: str) -> List[float]:
    values: List[float] = []
    seen: set[float] = set()
    for raw in text.split(","):
        token = raw.strip()
        if not token:
            continue
        value = float(token)
        key = round(value, 12)
        if key in seen:
            continue
        seen.add(key)
        values.append(value)
    return values


def parse_existing_overlays(values: Sequence[str]) -> Dict[float, str]:
    overlays: Dict[float, str] = {}
    for raw in values:
        if "=" not in raw:
            raise ValueError(f"Expected power=tag, got {raw!r}")
        power_text, tag = raw.split("=", 1)
        power = round(float(power_text.strip()), 12)
        tag = tag.strip()
        if not tag:
            raise ValueError(f"Missing tag in include-existing value {raw!r}")
        if power in overlays:
            raise ValueError(f"Duplicate include-existing power {power_text!r}")
        overlays[power] = tag
    return overlays


def maybe_replace_line(stripped: str, prefix: str, replacement: str) -> str | None:
    if stripped.startswith(prefix):
        return replacement
    return None


def sampling_power_block(power: float) -> List[str]:
    return [f"set {obj}:SamplingPower {fmt_power_token(power)}" for obj in SAMPLING_POWER_OBJECTS]


def finite_width_block() -> List[str]:
    block: List[str] = []
    for obj in FINITE_WIDTH_OBJECTS:
        block.append(f"get {obj}:UseFiniteWidthSpacelikeZPropagator")
        block.append(f"set {obj}:UseFiniteWidthSpacelikeZPropagator Yes")
        block.append(f"get {obj}:UseFiniteWidthSpacelikeZPropagator")
    return block


def patch_card_text(
    source_text: str,
    generated_stem: str,
    window: CutWindow | None = None,
    sampling_power: float | None = None,
    enforce_finite_width: bool = False,
    z_diag_overrides: Mapping[str, int] | None = None,
) -> str:
    replacements = {}
    if window is not None:
        replacements = {
            "set /Herwig/Cuts/NeutralCurrentCut:MinQ2 ": f"set /Herwig/Cuts/NeutralCurrentCut:MinQ2 {window.q2_min:.0f}.*GeV2",
            "set /Herwig/Cuts/NeutralCurrentCut:MaxQ2 ": f"set /Herwig/Cuts/NeutralCurrentCut:MaxQ2 {window.q2_max:.0f}.*GeV2",
            "set /Herwig/Cuts/NeutralCurrentCut:Miny ": f"set /Herwig/Cuts/NeutralCurrentCut:Miny {window.y_min}",
            "set /Herwig/Cuts/NeutralCurrentCut:Maxy ": f"set /Herwig/Cuts/NeutralCurrentCut:Maxy {window.y_max}",
        }
    z_diag_overrides = dict(z_diag_overrides or {})
    source_has_sampling = "SamplingPower" in source_text
    source_has_finite_width = "UseFiniteWidthSpacelikeZPropagator" in source_text
    lines = source_text.splitlines()
    patched: List[str] = []
    found = {key: False for key in replacements}
    saverun_found = False
    inserted_sampling = False
    inserted_finite_width = False

    sampling_re = re.compile(r"^\s*#?\s*set\s+(\S+):SamplingPower\b")
    finite_width_seen = False

    for line in lines:
        stripped = line.strip()
        replaced = False
        for prefix, replacement in replacements.items():
            new_line = maybe_replace_line(stripped, prefix, replacement)
            if new_line is not None:
                patched.append(new_line)
                found[prefix] = True
                replaced = True
                break
        if replaced:
            continue

        if sampling_power is not None:
            match = sampling_re.match(line)
            if match and match.group(1) in SAMPLING_POWER_OBJECTS:
                patched.append(f"set {match.group(1)}:SamplingPower {fmt_power_token(sampling_power)}")
                replaced = True
        if replaced:
            continue

        if enforce_finite_width and "UseFiniteWidthSpacelikeZPropagator" in stripped:
            finite_width_seen = True
            if stripped.startswith("set "):
                lhs = stripped.split(None, 2)[1]
                patched.append(f"set {lhs.split()[0]} Yes")
            else:
                patched.append(line)
            continue

        if stripped.startswith("set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditInitialSamples ") and "NLOAuditInitialSamples" in z_diag_overrides:
            patched.append(
                f"set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditInitialSamples {z_diag_overrides['NLOAuditInitialSamples']}"
            )
            continue
        if stripped.startswith("set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditSamplePeriod ") and "NLOAuditSamplePeriod" in z_diag_overrides:
            patched.append(
                f"set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditSamplePeriod {z_diag_overrides['NLOAuditSamplePeriod']}"
            )
            continue
        if stripped.startswith("set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOTermDiagnosticPeriod ") and "NLOTermDiagnosticPeriod" in z_diag_overrides:
            patched.append(
                f"set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOTermDiagnosticPeriod {z_diag_overrides['NLOTermDiagnosticPeriod']}"
            )
            continue

        if stripped == "cd /Herwig/MatrixElements/":
            patched.append(line)
            if sampling_power is not None and not source_has_sampling and not inserted_sampling:
                patched.extend(sampling_power_block(sampling_power))
                inserted_sampling = True
            if enforce_finite_width and not source_has_finite_width and not inserted_finite_width:
                patched.extend(finite_width_block())
                inserted_finite_width = True
            continue

        if stripped.startswith("saverun ") and stripped.endswith(" EventGenerator"):
            patched.append(f"saverun {generated_stem} EventGenerator")
            saverun_found = True
            continue

        patched.append(line)

    missing = [prefix.strip() for prefix, ok in found.items() if not ok]
    if missing:
        raise RuntimeError(f"Could not patch required cut lines: {', '.join(missing)}")
    if not saverun_found:
        raise RuntimeError("Could not find saverun line to rewrite in generated card")
    if enforce_finite_width and source_has_finite_width and not finite_width_seen:
        raise RuntimeError("Expected to see finite-width propagator lines in source card but did not")

    rendered = "\n".join(patched)
    if source_text.endswith("\n"):
        rendered += "\n"
    return rendered


def materialize_cards(base_dir: Path, runs: Sequence[LogicalRun], patch_by_key: Mapping[str, dict]) -> None:
    for run in runs:
        source_path = base_dir / run.source_card
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source card {source_path}")
        target_path = base_dir / run.generated_input_rel
        rendered = patch_card_text(source_path.read_text(), run.generated_stem, **patch_by_key[run.key])
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
            measurement, _ = parse_measurement(path.read_text())
            if measurement is None:
                raise RuntimeError(f"Failed to parse cross section from {path}")
            measurements.append(measurement)
            generated_events.append(parse_generated_event_count(path.read_text()))
        resolved[key] = {
            "measurement": combine_measurements(measurements, generated_events),
            "paths": [str(path) for path in paths],
        }
    return resolved


def require_measurement(base_dir: Path, requested_name: str, preferred_tag: str, context: str):
    run = load_run(base_dir, requested_name, preferred_tag, strict_tag=True)
    if not run.exists or run.measurement is None:
        raise FileNotFoundError(
            f"Missing {context}: could not resolve {requested_name} with tag {preferred_tag!r} under {base_dir}"
        )
    return run


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


def gamma_requested_name_map() -> Dict[str, str]:
    out: Dict[str, str] = {}
    for piece in NLO_PIECES:
        for helicity in GAMMA_HELICITIES:
            out[f"{piece}_{helicity}"] = build_gamma_requested_name(piece, helicity)
    return out


def collect_gamma_requested_runs(base_dir: Path, tag: str, context: str) -> Dict[str, dict]:
    requested = gamma_requested_name_map()
    out: Dict[str, dict] = {}
    for key, requested_name in requested.items():
        run = require_measurement(base_dir, requested_name, tag, context)
        out[key] = {
            "measurement": run.measurement,
            "paths": [str(path) for path in run.paths],
        }
    return out


def render_gamma_summary_section(title: str, payload: dict) -> List[str]:
    runs = payload["runs"]
    derived = payload["derived"]
    return [
        title,
        "-" * len(title),
        f"POSNLO PP:    {fmt_measurement(runs['POSNLO_PP'])}",
        f"POSNLO PM:    {fmt_measurement(runs['POSNLO_PM'])}",
        f"POSNLO (PP-PM)/2: {fmt_measurement(derived['POSNLO_pol'])}",
        f"NEGNLO PP:    {fmt_measurement(runs['NEGNLO_PP'])}",
        f"NEGNLO PM:    {fmt_measurement(runs['NEGNLO_PM'])}",
        f"NEGNLO (PP-PM)/2: {fmt_measurement(derived['NEGNLO_pol'])}",
        f"NLO (PP-PM)/2:    {fmt_measurement(derived['NLO_pol'])}",
        "",
    ]


def load_poldis_reference(path: Path) -> dict:
    payload = json.loads(path.read_text())
    try:
        measurement_from_payload(payload["polarized"]["NLO"])
        measurement_from_payload(payload["unpolarized"]["NLO"])
    except Exception as exc:  # pragma: no cover - defensive validation
        raise RuntimeError(f"Invalid POLDIS reference JSON {path}: missing NLO totals") from exc
    return {
        "path": str(path.resolve()),
        "window": payload.get("window"),
        "unpolarized": {"NLO": payload["unpolarized"]["NLO"]},
        "polarized": {"NLO": payload["polarized"]["NLO"]},
    }


def build_gamma_interior_payload(
    base_dir: Path,
    tag: str,
    control_tag: str | None,
    interior_summary: dict,
    broad_summary: dict | None,
    poldis_reference: dict | None = None,
) -> dict:
    payload = {
        "tag": tag,
        "mode": GAMMA_INTERIOR_MODE,
        "control_tag": control_tag,
        "broad_window": asdict(BROAD_GAMMA_WINDOW),
        "interior_window": asdict(INTERIOR_GAMMA_WINDOW),
        "interior_run": interior_summary,
    }
    comparisons: dict = {}
    if broad_summary is not None:
        payload["broad_control"] = broad_summary
        comparisons["POSNLO_pol_minus_broad"] = measurement_payload(
            sub_measurement(
                measurement_from_payload(interior_summary["derived"]["POSNLO_pol"]),
                measurement_from_payload(broad_summary["derived"]["POSNLO_pol"]),
            )
        )
        comparisons["NEGNLO_pol_minus_broad"] = measurement_payload(
            sub_measurement(
                measurement_from_payload(interior_summary["derived"]["NEGNLO_pol"]),
                measurement_from_payload(broad_summary["derived"]["NEGNLO_pol"]),
            )
        )
        comparisons["NLO_pol_minus_broad"] = measurement_payload(
            sub_measurement(
                measurement_from_payload(interior_summary["derived"]["NLO_pol"]),
                measurement_from_payload(broad_summary["derived"]["NLO_pol"]),
            )
        )
    if poldis_reference is not None:
        payload["poldis_reference"] = poldis_reference
        comparisons["NLO_pol_minus_poldis"] = measurement_payload(
            sub_measurement(
                measurement_from_payload(interior_summary["derived"]["NLO_pol"]),
                measurement_from_payload(poldis_reference["polarized"]["NLO"]),
            )
        )
    if comparisons:
        payload["comparisons"] = comparisons
    return payload


def render_gamma_interior_text(payload: dict) -> str:
    lines = [
        "Polarized GAMMA NLO interior-window study",
        "=========================================",
        f"Tag: {payload['tag']}",
        f"Control tag: {payload.get('control_tag') or 'none'}",
        "",
        f"Broad window:    Q^2 in [{payload['broad_window']['q2_min']:.0f}, {payload['broad_window']['q2_max']:.0f}] GeV^2, "
        f"y in [{payload['broad_window']['y_min']:.1f}, {payload['broad_window']['y_max']:.1f}]",
        f"Interior window: Q^2 in [{payload['interior_window']['q2_min']:.0f}, {payload['interior_window']['q2_max']:.0f}] GeV^2, "
        f"y in [{payload['interior_window']['y_min']:.1f}, {payload['interior_window']['y_max']:.1f}]",
        "",
    ]
    if "broad_control" in payload:
        lines.extend(render_gamma_summary_section("Broad control", payload["broad_control"]))
    lines.extend(render_gamma_summary_section("Interior run", payload["interior_run"]))
    if "poldis_reference" in payload:
        reference = payload["poldis_reference"]
        lines.extend(
            [
                "POLDIS reference",
                "----------------",
                f"NLO unpolarized: {fmt_measurement(reference['unpolarized']['NLO'])}",
                f"NLO polarized:   {fmt_measurement(reference['polarized']['NLO'])}",
                f"Source: {reference['path']}",
                "",
            ]
        )
    if "comparisons" in payload:
        comparisons = payload["comparisons"]
        lines.extend(
            [
                "Comparisons",
                "-----------",
            ]
        )
        if "POSNLO_pol_minus_broad" in comparisons:
            lines.append(f"Interior POSNLO_pol - broad: {fmt_measurement(comparisons['POSNLO_pol_minus_broad'])}")
        if "NEGNLO_pol_minus_broad" in comparisons:
            lines.append(f"Interior NEGNLO_pol - broad: {fmt_measurement(comparisons['NEGNLO_pol_minus_broad'])}")
        if "NLO_pol_minus_broad" in comparisons:
            lines.append(f"Interior NLO_pol - broad:    {fmt_measurement(comparisons['NLO_pol_minus_broad'])}")
        if "NLO_pol_minus_poldis" in comparisons:
            lines.append(f"Interior NLO_pol - POLDIS:   {fmt_measurement(comparisons['NLO_pol_minus_poldis'])}")
    return "\n".join(lines).rstrip() + "\n"


def write_text_json_results(base_dir: Path, tag: str, mode: str, text: str, payload: dict) -> None:
    mode_dir = campaign_mode_dir(base_dir, tag, mode)
    mode_dir.mkdir(parents=True, exist_ok=True)
    results_txt_path(base_dir, tag, mode).write_text(text)
    results_json_path(base_dir, tag, mode).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(text, end="")
    print(f"Wrote text report: {results_txt_path(base_dir, tag, mode)}")
    print(f"Wrote JSON report: {results_json_path(base_dir, tag, mode)}")


def build_gamma_power_payload(entries: Sequence[dict], tag: str, include_existing: Mapping[float, str]) -> dict:
    return {
        "tag": tag,
        "mode": GAMMA_POWER_SCAN_MODE,
        "entries": list(entries),
        "include_existing": {fmt_power_token(power): include_existing[power] for power in sorted(include_existing)},
    }


def render_gamma_power_text(payload: dict) -> str:
    lines = [
        "Polarized GAMMA NLO SamplingPower scan",
        "======================================",
        f"Tag: {payload['tag']}",
        "",
        "Power  Source               POS_pol                  NEG_pol                  NLO_pol",
        "-----  -------------------  -----------------------  -----------------------  -----------------------",
    ]
    for entry in payload["entries"]:
        lines.append(
            f"{entry['power_label']:>5}  {entry['source']:19}  "
            f"{entry['derived']['POSNLO_pol']['value_pb']:>9.4f} +- {entry['derived']['POSNLO_pol']['error_pb']:<8.4f}  "
            f"{entry['derived']['NEGNLO_pol']['value_pb']:>9.4f} +- {entry['derived']['NEGNLO_pol']['error_pb']:<8.4f}  "
            f"{entry['derived']['NLO_pol']['value_pb']:>9.4f} +- {entry['derived']['NLO_pol']['error_pb']:<8.4f}"
        )
    return "\n".join(lines).rstrip() + "\n"


def write_gamma_power_csv(base_dir: Path, tag: str, entries: Sequence[dict]) -> Path:
    csv_path = campaign_mode_dir(base_dir, tag, GAMMA_POWER_SCAN_MODE) / "scan.csv"
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "power",
                "source",
                "source_tag",
                "pos_pp_pb",
                "pos_pp_err_pb",
                "pos_pm_pb",
                "pos_pm_err_pb",
                "neg_pp_pb",
                "neg_pp_err_pb",
                "neg_pm_pb",
                "neg_pm_err_pb",
                "pos_pol_pb",
                "pos_pol_err_pb",
                "neg_pol_pb",
                "neg_pol_err_pb",
                "nlo_pol_pb",
                "nlo_pol_err_pb",
            ]
        )
        for entry in entries:
            runs = entry["runs"]
            derived = entry["derived"]
            writer.writerow(
                [
                    entry["power"],
                    entry["source"],
                    entry["source_tag"],
                    runs["POSNLO_PP"]["value_pb"],
                    runs["POSNLO_PP"]["error_pb"],
                    runs["POSNLO_PM"]["value_pb"],
                    runs["POSNLO_PM"]["error_pb"],
                    runs["NEGNLO_PP"]["value_pb"],
                    runs["NEGNLO_PP"]["error_pb"],
                    runs["NEGNLO_PM"]["value_pb"],
                    runs["NEGNLO_PM"]["error_pb"],
                    derived["POSNLO_pol"]["value_pb"],
                    derived["POSNLO_pol"]["error_pb"],
                    derived["NEGNLO_pol"]["value_pb"],
                    derived["NEGNLO_pol"]["error_pb"],
                    derived["NLO_pol"]["value_pb"],
                    derived["NLO_pol"]["error_pb"],
                ]
            )
    return csv_path


def build_gamma_power_runs(tag: str, powers: Sequence[float]) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    for power in powers:
        suffix = power_suffix(power)
        for piece in NLO_PIECES:
            for helicity in GAMMA_HELICITIES:
                source_card = build_gamma_source_card(piece, helicity)
                source_stem = stem_without_suffix(source_card)
                generated_stem = f"{source_stem}-{suffix}"
                key = f"{fmt_power_token(power)}::{piece}_{helicity}"
                runs.append(
                    LogicalRun(
                        mode=GAMMA_POWER_SCAN_MODE,
                        key=key,
                        source_card=source_card,
                        requested_name=build_gamma_requested_name(piece, helicity),
                        generated_stem=generated_stem,
                        generated_input_rel=Path("campaigns") / tag / GAMMA_POWER_SCAN_MODE / "generated-inputs" / f"{generated_stem}.in",
                        run_file=f"{generated_stem}.run",
                        piece=piece,
                        helicity=helicity,
                        power=power,
                    )
                )
    return runs


def build_gamma_interior_runs(tag: str) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    for piece in NLO_PIECES:
        for helicity in GAMMA_HELICITIES:
            source_card = build_gamma_source_card(piece, helicity)
            source_stem = stem_without_suffix(source_card)
            generated_stem = f"{source_stem}-{GAMMA_INTERIOR_SUFFIX}"
            key = f"{piece}_{helicity}"
            runs.append(
                LogicalRun(
                    mode=GAMMA_INTERIOR_MODE,
                    key=key,
                    source_card=source_card,
                    requested_name=build_gamma_requested_name(piece, helicity),
                    generated_stem=generated_stem,
                    generated_input_rel=Path("campaigns") / tag / GAMMA_INTERIOR_MODE / "generated-inputs" / f"{generated_stem}.in",
                    run_file=f"{generated_stem}.run",
                    piece=piece,
                    helicity=helicity,
                )
            )
    return runs


def build_z_termdiag_runs(tag: str) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    for piece in NLO_PIECES:
        for helicity in GAMMA_HELICITIES:
            source_card = build_z_termdiag_source_card(piece, helicity)
            source_stem = stem_without_suffix(source_card)
            generated_stem = f"{source_stem}-{Z_TERMDIAG_SUFFIX}"
            key = f"{piece}_{helicity}"
            runs.append(
                LogicalRun(
                    mode=Z_TERMDIAG_MODE,
                    key=key,
                    source_card=source_card,
                    requested_name=build_z_termdiag_requested_name(piece, helicity),
                    generated_stem=generated_stem,
                    generated_input_rel=Path("campaigns") / tag / Z_TERMDIAG_MODE / "generated-inputs" / f"{generated_stem}.in",
                    run_file=f"{generated_stem}.run",
                    piece=piece,
                    helicity=helicity,
                )
            )
    return runs


def print_plan(base_dir: Path, tag: str, mode: str, runs: Sequence[LogicalRun], shards: Sequence[ShardSpec]) -> None:
    print(f"[plan:{mode}] base_dir={base_dir}")
    print(f"[plan:{mode}] logical_runs={len(runs)} shards={len(shards)}")
    print(f"[plan:{mode}] monitor_txt={status_txt_path(base_dir, tag, mode)}")
    print(f"[plan:{mode}] monitor_json={status_json_path(base_dir, tag, mode)}")
    for run in runs:
        print(
            f"[plan:{mode}] {run.key}: source={run.source_card} generated={run.generated_input_rel} run={run.run_file}"
        )
    for spec in shards[: min(12, len(shards))]:
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


def run_gamma_interior(args: argparse.Namespace, seed_base: int | None = None) -> None:
    base_dir = args.base_dir.resolve()
    tag = args.tag
    seed_base = args.seed_base if seed_base is None else seed_base
    runs = build_gamma_interior_runs(tag)
    patch_by_key = {run.key: {"window": INTERIOR_GAMMA_WINDOW} for run in runs}
    materialize_cards(base_dir, runs, patch_by_key)
    shards = build_shards(runs, tag, args.gamma_interior_shards, args.gamma_interior_events_per_shard, seed_base)

    if args.dry_run:
        print_plan(base_dir, tag, GAMMA_INTERIOR_MODE, runs, shards)
        ensure_run_files(base_dir, runs, dry_run=True)
        return

    if not args.collect_only:
        if args.control_tag:
            collect_gamma_requested_runs(base_dir, args.control_tag, "broad GAMMA control")
        ensure_run_files(base_dir, runs, dry_run=False)
        run_shards(base_dir, tag, GAMMA_INTERIOR_MODE, shards, args.jobs)

    generated = collect_explicit_sharded_measurements(
        base_dir,
        {run.key: run.generated_stem for run in runs},
        tag,
        context="gamma interior generated run",
    )
    interior_summary = gamma_summary_from_flat_runs(generated)
    broad_summary = None
    if args.control_tag:
        broad_summary = gamma_summary_from_flat_runs(
            collect_gamma_requested_runs(base_dir, args.control_tag, "broad GAMMA control")
        )
    poldis_reference = load_poldis_reference(args.poldis_reference.resolve()) if args.poldis_reference else None
    payload = build_gamma_interior_payload(
        base_dir,
        tag,
        args.control_tag,
        interior_summary,
        broad_summary,
        poldis_reference=poldis_reference,
    )
    write_text_json_results(base_dir, tag, GAMMA_INTERIOR_MODE, render_gamma_interior_text(payload), payload)


def run_gamma_power_scan(args: argparse.Namespace, seed_base: int | None = None) -> None:
    base_dir = args.base_dir.resolve()
    tag = args.tag
    seed_base = args.seed_base if seed_base is None else seed_base
    powers = parse_power_list(args.powers)
    existing = parse_existing_overlays(args.include_existing)
    power_keys = {round(power, 12) for power in powers}
    overlap = power_keys.intersection(existing.keys())
    if overlap:
        duplicate = ", ".join(fmt_power_token(value) for value in sorted(overlap))
        raise SystemExit(f"Duplicate powers between --powers and --include-existing: {duplicate}")

    runs = build_gamma_power_runs(tag, powers)
    patch_by_key = {run.key: {"sampling_power": run.power} for run in runs}
    materialize_cards(base_dir, runs, patch_by_key)
    shards = build_shards(
        runs,
        tag,
        args.gamma_power_scan_shards,
        args.gamma_power_scan_events_per_shard,
        seed_base,
    )

    if args.dry_run:
        print_plan(base_dir, tag, GAMMA_POWER_SCAN_MODE, runs, shards)
        ensure_run_files(base_dir, runs, dry_run=True)
        return

    if not args.collect_only and runs:
        ensure_run_files(base_dir, runs, dry_run=False)
        run_shards(base_dir, tag, GAMMA_POWER_SCAN_MODE, shards, args.jobs)

    entries: List[dict] = []

    if runs:
        generated_measurements = collect_explicit_sharded_measurements(
            base_dir,
            {run.key: run.generated_stem for run in runs},
            tag,
            context="gamma power-scan generated run",
        )
        grouped: Dict[float, Dict[str, dict]] = {}
        for run in runs:
            assert run.power is not None
            grouped.setdefault(run.power, {})[f"{run.piece}_{run.helicity}"] = generated_measurements[run.key]
        for power in sorted(grouped):
            summary = gamma_summary_from_flat_runs(grouped[power])
            entries.append(
                {
                    "power": power,
                    "power_label": fmt_power_token(power),
                    "source": "generated",
                    "source_tag": tag,
                    "runs": summary["runs"],
                    "derived": summary["derived"],
                }
            )

    for power in sorted(existing):
        summary = gamma_summary_from_flat_runs(
            collect_gamma_requested_runs(base_dir, existing[power], f"existing gamma scan power {fmt_power_token(power)}")
        )
        entries.append(
            {
                "power": power,
                "power_label": fmt_power_token(power),
                "source": f"existing:{existing[power]}",
                "source_tag": existing[power],
                "runs": summary["runs"],
                "derived": summary["derived"],
            }
        )

    entries.sort(key=lambda item: item["power"])
    payload = build_gamma_power_payload(entries, tag, existing)
    text = render_gamma_power_text(payload)
    write_text_json_results(base_dir, tag, GAMMA_POWER_SCAN_MODE, text, payload)
    csv_path = write_gamma_power_csv(base_dir, tag, entries)
    print(f"Wrote CSV report: {csv_path}")


def extractor_output_paths(base_dir: Path, tag: str) -> tuple[Path, Path, Path]:
    mode_dir = campaign_mode_dir(base_dir, tag, Z_TERMDIAG_MODE)
    return mode_dir / "extract.txt", mode_dir / "extract.json", mode_dir / "extract.csv"


def run_z_termdiag_extractor(base_dir: Path, tag: str, dry_run: bool) -> None:
    text_path, json_path, csv_path = extractor_output_paths(base_dir, tag)
    cmd = [
        sys.executable,
        str(base_dir / "extract_nlo_term_diagnostics.py"),
        "--base-dir",
        str(base_dir),
        "--tag",
        tag,
        "--strict-tag",
        "--setup",
        "Z",
        "--piece",
        "POSNLO",
        "--piece",
        "NEGNLO",
        "--helicity",
        "PP",
        "--helicity",
        "PM",
        "--name-filter",
        "TERMDIAG",
        "--z-spin-report",
        "--json-out",
        str(json_path),
        "--csv-out",
        str(csv_path),
    ]
    if dry_run:
        print(shlex.join(cmd))
        return
    mode_dir = campaign_mode_dir(base_dir, tag, Z_TERMDIAG_MODE)
    mode_dir.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(cmd, cwd=base_dir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"extract_nlo_term_diagnostics.py failed:\n{proc.stderr}")
    text_path.write_text(proc.stdout)
    print(proc.stdout, end="")
    print(f"Wrote extractor text report: {text_path}")
    print(f"Wrote extractor JSON report: {json_path}")
    print(f"Wrote extractor CSV report: {csv_path}")


def run_z_termdiag(args: argparse.Namespace, seed_base: int | None = None) -> None:
    base_dir = args.base_dir.resolve()
    tag = args.tag
    seed_base = args.seed_base if seed_base is None else seed_base
    overrides: Dict[str, int] = {}
    if args.nlo_audit_initial_samples is not None:
        overrides["NLOAuditInitialSamples"] = args.nlo_audit_initial_samples
    if args.nlo_audit_sample_period is not None:
        overrides["NLOAuditSamplePeriod"] = args.nlo_audit_sample_period
    if args.nlo_term_diagnostic_period is not None:
        overrides["NLOTermDiagnosticPeriod"] = args.nlo_term_diagnostic_period

    runs = build_z_termdiag_runs(tag)
    patch_by_key = {
        run.key: {
            "enforce_finite_width": True,
            "z_diag_overrides": overrides,
        }
        for run in runs
    }
    materialize_cards(base_dir, runs, patch_by_key)
    shards = build_shards(runs, tag, args.z_termdiag_shards, args.z_termdiag_events_per_shard, seed_base)

    if args.dry_run:
        print_plan(base_dir, tag, Z_TERMDIAG_MODE, runs, shards)
        ensure_run_files(base_dir, runs, dry_run=True)
        run_z_termdiag_extractor(base_dir, tag, dry_run=True)
        return

    if not args.collect_only:
        ensure_run_files(base_dir, runs, dry_run=False)
        run_shards(base_dir, tag, Z_TERMDIAG_MODE, shards, args.jobs)

    run_z_termdiag_extractor(base_dir, tag, dry_run=False)


def run_all(args: argparse.Namespace) -> None:
    base_seed = args.seed_base
    run_gamma_interior(args, seed_base=base_seed)

    powers = parse_power_list(args.powers)
    gamma_power_seed = base_seed + args.gamma_interior_shards * 4
    run_gamma_power_scan(args, seed_base=gamma_power_seed)

    generated_power_runs = len(powers) * 4
    z_seed = gamma_power_seed + args.gamma_power_scan_shards * generated_power_runs
    run_z_termdiag(args, seed_base=z_seed)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.base_dir.resolve().exists():
        raise SystemExit(f"Base directory does not exist: {args.base_dir}")

    if args.mode == GAMMA_INTERIOR_MODE:
        run_gamma_interior(args)
    elif args.mode == GAMMA_POWER_SCAN_MODE:
        run_gamma_power_scan(args)
    elif args.mode == Z_TERMDIAG_MODE:
        run_z_termdiag(args)
    elif args.mode == ALL_MODE:
        run_all(args)
    else:
        raise SystemExit(f"Unsupported mode: {args.mode}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
