#!/usr/bin/env python3.10
"""
Run a focused LO Z interior-window study and compare sigma0 vs sigma00.

The script:
  1. materializes interior-window copies of the standard LO Z cards
  2. prepares missing .run files with `Herwig read`
  3. launches 100 x 1e6-event shards for each LO Z helicity
  4. reuses broad-window plain45 shard outputs as the control sample
  5. collects Herwig and standalone reference results into results.txt/json
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import shlex
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

from extract_dis_out_results import (
    Measurement,
    combine_measurements,
    linear_combo,
    load_run,
    parse_generated_event_count,
    parse_measurement,
)


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
HELICITIES = ("00", "PP", "PM", "MP", "MM")
INTERIOR_SUFFIX = "ZINTCHK"
DEFAULT_TAG = "plain46-zlo-int"
DEFAULT_CONTROL_TAG = "plain45"
DEFAULT_JOBS = 32
DEFAULT_SHARDS = 100
DEFAULT_EVENTS_PER_SHARD = 1_000_000
DEFAULT_SEED_BASE = 600_000
FINITE_WIDTH_OBJECTS = ("MEDISNC", "MEDISNCPol", "PowhegMEDISNC", "PowhegMEDISNCPol")


@dataclass(frozen=True)
class CutWindow:
    label: str
    q2_min: float
    q2_max: float
    y_min: float
    y_max: float


@dataclass(frozen=True)
class LogicalRun:
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


@dataclass
class ShardResult:
    spec: ShardSpec
    command: List[str]
    returncode: int
    duration_s: float
    launcher_log: Path


BROAD_WINDOW = CutWindow(label="broad", q2_min=49.0, q2_max=2500.0, y_min=0.2, y_max=0.6)
INTERIOR_WINDOW = CutWindow(label="interior", q2_min=100.0, q2_max=1000.0, y_min=0.3, y_max=0.5)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used in Herwig shard tags.")
    parser.add_argument(
        "--control-tag",
        default=DEFAULT_CONTROL_TAG,
        help="Existing broad-window control shard tag to reuse for DIS-POL-LO_* -Z outputs.",
    )
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Maximum concurrent Herwig processes.")
    parser.add_argument("--shards", type=int, default=DEFAULT_SHARDS, help="Number of shards per helicity.")
    parser.add_argument(
        "--events-per-shard",
        type=int,
        default=DEFAULT_EVENTS_PER_SHARD,
        help="Events per shard for each helicity.",
    )
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE, help="Base RNG seed.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing the DISPOL cards and outputs.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned commands without executing Herwig.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip Herwig read/run and only collect results from existing shard outputs.",
    )
    return parser


def build_source_stem(helicity: str) -> str:
    return f"DIS-POL-LO_{helicity}-Z"


def build_generated_stem(helicity: str) -> str:
    return f"{build_source_stem(helicity)}-{INTERIOR_SUFFIX}"


def build_logical_runs(tag: str) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    for helicity in HELICITIES:
        stem = build_generated_stem(helicity)
        runs.append(
            LogicalRun(
                helicity=helicity,
                source_card=f"{build_source_stem(helicity)}.in",
                generated_stem=stem,
                generated_input_rel=Path("campaigns") / tag / "generated-inputs" / f"{stem}.in",
                run_file=f"{stem}.run",
            )
        )
    return runs


def patch_card_text(source_text: str, generated_stem: str, window: CutWindow) -> str:
    replacements = {
        "set /Herwig/Cuts/NeutralCurrentCut:MinQ2 ": f"set /Herwig/Cuts/NeutralCurrentCut:MinQ2 {window.q2_min:.0f}.*GeV2",
        "set /Herwig/Cuts/NeutralCurrentCut:MaxQ2 ": f"set /Herwig/Cuts/NeutralCurrentCut:MaxQ2 {window.q2_max:.0f}.*GeV2",
        "set /Herwig/Cuts/NeutralCurrentCut:Miny ": f"set /Herwig/Cuts/NeutralCurrentCut:Miny {window.y_min}",
        "set /Herwig/Cuts/NeutralCurrentCut:Maxy ": f"set /Herwig/Cuts/NeutralCurrentCut:Maxy {window.y_max}",
    }
    source_has_finite_width_block = "UseFiniteWidthSpacelikeZPropagator" in source_text
    finite_width_seen = False
    inserted_finite_width_block = False

    lines = source_text.splitlines()
    found = {key: False for key in replacements}
    saverun_found = False
    patched: List[str] = []

    for line in lines:
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
        if (
            stripped == "cd /Herwig/MatrixElements/"
            and not source_has_finite_width_block
            and not inserted_finite_width_block
        ):
            patched.append(line)
            for obj in FINITE_WIDTH_OBJECTS:
                patched.append(f"get {obj}:UseFiniteWidthSpacelikeZPropagator")
                patched.append(f"set {obj}:UseFiniteWidthSpacelikeZPropagator Yes")
                patched.append(f"get {obj}:UseFiniteWidthSpacelikeZPropagator")
            inserted_finite_width_block = True
            continue
        if "UseFiniteWidthSpacelikeZPropagator" in stripped:
            finite_width_seen = True
            if stripped.startswith("set "):
                lhs = stripped.split(None, 2)[1]
                patched.append(f"set {lhs.split()[0]} Yes")
            else:
                patched.append(line)
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
        raise RuntimeError("Could not find saverun line to rewrite in generated interior card")

    rendered = "\n".join(patched)
    if source_text.endswith("\n"):
        rendered += "\n"
    return rendered


def materialize_generated_cards(base_dir: Path, tag: str, window: CutWindow) -> List[LogicalRun]:
    runs = build_logical_runs(tag)
    for run in runs:
        source_path = base_dir / run.source_card
        target_path = base_dir / run.generated_input_rel
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source card {source_path}")
        rendered = patch_card_text(source_path.read_text(), run.generated_stem, window)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        if not target_path.exists() or target_path.read_text() != rendered:
            target_path.write_text(rendered)
    return runs


def ensure_run_files(base_dir: Path, runs: Sequence[LogicalRun], dry_run: bool) -> None:
    for run in runs:
        in_path = base_dir / run.generated_input_rel
        run_path = base_dir / run.run_file
        cmd = ["Herwig", "read", str(run.generated_input_rel)]
        if dry_run:
            print(shlex.join(cmd))
            continue
        if run_path.exists() and run_path.stat().st_mtime >= in_path.stat().st_mtime:
            continue
        proc = subprocess.run(cmd, cwd=base_dir, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f"Failed to prepare {run.run_file} from {in_path}")


def build_shards(runs: Sequence[LogicalRun], tag: str, shards: int, events_per_shard: int, seed_base: int) -> List[ShardSpec]:
    specs: List[ShardSpec] = []
    seed = seed_base
    for run in runs:
        for shard_index in range(1, shards + 1):
            shard_tag = f"{tag}-s{shard_index:03d}"
            specs.append(
                ShardSpec(
                    logical_run=run,
                    shard_index=shard_index,
                    shard_count=shards,
                    tag=shard_tag,
                    seed=seed,
                    events=events_per_shard,
                )
            )
            seed += 1
    return specs


def launcher_log_path(base_dir: Path, tag: str, spec: ShardSpec) -> Path:
    return base_dir / "campaigns" / tag / "launcher-logs" / f"{spec.logical_run.generated_stem}-{spec.tag}.launcher.log"


def monitor_dir(base_dir: Path, tag: str) -> Path:
    return base_dir / "campaigns" / tag / "monitor"


def status_json_path(base_dir: Path, tag: str) -> Path:
    return monitor_dir(base_dir, tag) / "status.json"


def status_txt_path(base_dir: Path, tag: str) -> Path:
    return monitor_dir(base_dir, tag) / "status.txt"


def build_monitor_payload(
    base_dir: Path,
    tag: str,
    control_tag: str,
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

    per_helicity: Dict[str, dict] = {}
    for helicity in HELICITIES:
        hel_specs = [spec for spec in shards if spec.logical_run.helicity == helicity]
        hel_results = [result for result in results if result.spec.logical_run.helicity == helicity]
        hel_completed = len(hel_results)
        hel_failures = sum(1 for result in hel_results if result.returncode != 0)
        hel_pending = len(hel_specs) - hel_completed
        hel_mean_duration = (
            sum(result.duration_s for result in hel_results) / hel_completed if hel_completed else None
        )
        per_helicity[helicity] = {
            "total": len(hel_specs),
            "completed": hel_completed,
            "pending": hel_pending,
            "failures": hel_failures,
            "mean_duration_s": hel_mean_duration,
        }

    payload = {
        "tag": tag,
        "control_tag": control_tag,
        "phase": phase,
        "updated_unix_s": time.time(),
        "elapsed_s": elapsed_s,
        "total_shards": total,
        "completed_shards": completed,
        "pending_shards": pending,
        "failed_shards": failures,
        "mean_duration_s": mean_duration_s,
        "eta_s": eta_s,
        "per_helicity": per_helicity,
        "monitor_files": {
            "status_json": str(status_json_path(base_dir, tag)),
            "status_txt": str(status_txt_path(base_dir, tag)),
        },
    }
    if last_result is not None:
        payload["last_completed"] = {
            "helicity": last_result.spec.logical_run.helicity,
            "generated_stem": last_result.spec.logical_run.generated_stem,
            "shard_index": last_result.spec.shard_index,
            "tag": last_result.spec.tag,
            "seed": last_result.spec.seed,
            "returncode": last_result.returncode,
            "duration_s": last_result.duration_s,
            "launcher_log": str(last_result.launcher_log),
        }
    return payload


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def render_monitor_text(payload: Mapping[str, object]) -> str:
    lines = [
        f"Tag: {payload['tag']}",
        f"Control tag: {payload['control_tag']}",
        f"Phase: {payload['phase']}",
        f"Completed: {payload['completed_shards']}/{payload['total_shards']}",
        f"Pending: {payload['pending_shards']}",
        f"Failures: {payload['failed_shards']}",
        f"Elapsed: {fmt_seconds(float(payload['elapsed_s']))}",
        f"ETA: {fmt_seconds(float(payload['eta_s'])) if payload['eta_s'] is not None else 'n/a'}",
        "",
        "Per helicity",
        "-----------",
    ]
    per_helicity = payload["per_helicity"]
    assert isinstance(per_helicity, Mapping)
    for helicity in HELICITIES:
        hel = per_helicity[helicity]
        assert isinstance(hel, Mapping)
        lines.append(
            f"{helicity}: completed={hel['completed']}/{hel['total']} pending={hel['pending']} "
            f"failures={hel['failures']} mean_duration={fmt_seconds(float(hel['mean_duration_s'])) if hel['mean_duration_s'] is not None else 'n/a'}"
        )
    last_completed = payload.get("last_completed")
    if isinstance(last_completed, Mapping):
        lines.extend(
            [
                "",
                "Last completed",
                "--------------",
                f"{last_completed['generated_stem']}[{int(last_completed['shard_index']):03d}] "
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
    return ShardResult(spec=spec, command=cmd, returncode=proc.returncode, duration_s=ended - started, launcher_log=log_path)


def run_shards(base_dir: Path, tag: str, control_tag: str, shards: Sequence[ShardSpec], jobs: int) -> List[ShardResult]:
    total = len(shards)
    completed = 0
    failures = 0
    results: List[ShardResult] = []
    started = time.time()

    write_monitor_files(
        base_dir,
        tag,
        build_monitor_payload(base_dir, tag, control_tag, shards, results, started, phase="running"),
    )

    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, jobs)) as executor:
        future_map = {executor.submit(run_one_shard, base_dir, tag, spec): spec for spec in shards}
        for future in concurrent.futures.as_completed(future_map):
            result = future.result()
            results.append(result)
            completed += 1
            if result.returncode != 0:
                failures += 1
            write_monitor_files(
                base_dir,
                tag,
                build_monitor_payload(
                    base_dir,
                    tag,
                    control_tag,
                    shards,
                    results,
                    started,
                    last_result=result,
                    phase="running" if failures == 0 and completed < total else "failed" if failures else "complete",
                ),
            )
            if completed == 1 or completed % 10 == 0 or completed == total or result.returncode != 0:
                elapsed = time.time() - started
                print(
                    f"[progress] {completed}/{total} shards complete, failures={failures}, "
                    f"elapsed={elapsed:.1f}s, last={result.spec.logical_run.generated_stem}[{result.spec.shard_index:03d}] rc={result.returncode}",
                    flush=True,
                )

    bad = [result for result in results if result.returncode != 0]
    if bad:
        write_monitor_files(
            base_dir,
            tag,
            build_monitor_payload(base_dir, tag, control_tag, shards, results, started, last_result=bad[-1], phase="failed"),
        )
        summary = "\n".join(
            f"  {item.spec.logical_run.generated_stem}[{item.spec.shard_index:03d}] rc={item.returncode} log={item.launcher_log}"
            for item in bad[:10]
        )
        raise RuntimeError(f"{len(bad)} shard(s) failed:\n{summary}")
    write_monitor_files(
        base_dir,
        tag,
        build_monitor_payload(base_dir, tag, control_tag, shards, results, started, last_result=results[-1] if results else None, phase="complete"),
    )
    return results


def require_measurement(base_dir: Path, requested_name: str, preferred_tag: str, context: str):
    run = load_run(base_dir, requested_name, preferred_tag, strict_tag=True)
    if not run.exists or run.measurement is None:
        raise FileNotFoundError(
            f"Missing {context}: could not resolve {requested_name} with tag {preferred_tag!r} under {base_dir}"
        )
    return run


def collect_helicity_measurements(base_dir: Path, names: Mapping[str, str], tag: str, context: str) -> Dict[str, dict]:
    resolved: Dict[str, dict] = {}
    for helicity, requested_name in names.items():
        run = require_measurement(base_dir, requested_name, tag, context)
        resolved[helicity] = {
            "measurement": run.measurement,
            "paths": [str(path) for path in run.paths],
            "requested_name": requested_name,
        }
    return resolved


def collect_explicit_sharded_measurements(
    base_dir: Path,
    stems: Mapping[str, str],
    tag: str,
    context: str,
) -> Dict[str, dict]:
    resolved: Dict[str, dict] = {}
    for helicity, stem in stems.items():
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

        resolved[helicity] = {
            "measurement": combine_measurements(measurements, generated_events),
            "paths": [str(path) for path in paths],
            "requested_name": f"{stem}.out",
        }
    return resolved


def build_requested_names_for_broad() -> Dict[str, str]:
    return {helicity: f"{build_source_stem(helicity)}.out" for helicity in HELICITIES}


def build_requested_names_for_interior() -> Dict[str, str]:
    return {helicity: f"{build_generated_stem(helicity)}.out" for helicity in HELICITIES}


def build_generated_stems_for_interior() -> Dict[str, str]:
    return {helicity: build_generated_stem(helicity) for helicity in HELICITIES}


def measurement_payload(measurement: Measurement) -> dict:
    return {"value_pb": measurement.value_pb, "error_pb": measurement.error_pb}


def compute_window_summary(helicity_runs: Mapping[str, dict]) -> dict:
    measurements = {helicity: helicity_runs[helicity]["measurement"] for helicity in HELICITIES}
    sigma0 = linear_combo([(0.25, measurements["PP"]), (0.25, measurements["PM"]), (0.25, measurements["MP"]), (0.25, measurements["MM"])])
    sigma_ll = linear_combo([(0.25, measurements["PP"]), (-0.25, measurements["PM"]), (-0.25, measurements["MP"]), (0.25, measurements["MM"])])
    sigma0_minus_00 = linear_combo([(1.0, sigma0), (-1.0, measurements["00"])])
    return {
        "helicities": {helicity: measurement_payload(measurements[helicity]) for helicity in HELICITIES},
        "derived": {
            "00": measurement_payload(measurements["00"]),
            "sigma0": measurement_payload(sigma0),
            "sigma_LL": measurement_payload(sigma_ll),
            "sigma0_minus_00": measurement_payload(sigma0_minus_00),
        },
    }


def load_calc_payload(stdout: str) -> dict:
    marker = "\nJSON\n----\n"
    if marker not in stdout:
        raise RuntimeError("calc_lo_dis_nc.py did not emit the expected JSON marker")
    return json.loads(stdout.split(marker, 1)[1])


def run_standalone_reference(base_dir: Path, window: CutWindow, herwig_lo_mode: bool) -> dict:
    cmd = [
        sys.executable,
        str(base_dir / "calc_lo_dis_nc.py"),
        "integrate",
        "--channel",
        "Z",
        "--flavor-scope",
        "total",
        "--q2-min",
        str(window.q2_min),
        "--q2-max",
        str(window.q2_max),
        "--y-min",
        str(window.y_min),
        "--y-max",
        str(window.y_max),
        "--json",
    ]
    if herwig_lo_mode:
        cmd.append("--herwig-lo-mode")
    proc = subprocess.run(cmd, cwd=base_dir, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Standalone calculator failed with rc={proc.returncode}:\n{proc.stderr}")
    payload = load_calc_payload(proc.stdout)
    totals = payload["integrate"]["totals"]
    if not totals:
        raise RuntimeError("Standalone calculator returned no integrated totals")
    row = totals[0]
    return {
        "mode": "herwig_lo" if herwig_lo_mode else "poldis",
        "sigma_unpol_pb": float(row["sigma_unpol_pb"]),
        "sigma_LL_pb": float(row["sigma_ll_pb"]),
    }


def delta_payload(value: Measurement, reference_value_pb: float) -> dict:
    delta = value.value_pb - reference_value_pb
    return {
        "value_pb": delta,
        "error_pb": value.error_pb,
        "reference_value_pb": reference_value_pb,
    }


def build_results_payload(
    base_dir: Path,
    tag: str,
    control_tag: str,
    broad_summary: dict,
    interior_summary: dict,
) -> dict:
    broad_poldis = run_standalone_reference(base_dir, BROAD_WINDOW, herwig_lo_mode=False)
    broad_herwig = run_standalone_reference(base_dir, BROAD_WINDOW, herwig_lo_mode=True)
    interior_poldis = run_standalone_reference(base_dir, INTERIOR_WINDOW, herwig_lo_mode=False)
    interior_herwig = run_standalone_reference(base_dir, INTERIOR_WINDOW, herwig_lo_mode=True)

    payload = {
        "tag": tag,
        "control_tag": control_tag,
        "broad_window": asdict(BROAD_WINDOW),
        "interior_window": asdict(INTERIOR_WINDOW),
        "broad_control": {
            "herwig": broad_summary,
            "standalone": {"poldis_like": broad_poldis, "herwig_like": broad_herwig},
        },
        "interior_run": {
            "herwig": interior_summary,
            "standalone": {"poldis_like": interior_poldis, "herwig_like": interior_herwig},
        },
    }

    for section_name, summary, poldis_ref, herwig_ref in (
        ("broad_control", broad_summary, broad_poldis, broad_herwig),
        ("interior_run", interior_summary, interior_poldis, interior_herwig),
    ):
        derived = summary["derived"]
        payload[section_name]["comparisons"] = {
            "00_vs_herwig_like_unpol": delta_payload(
                Measurement(**derived["00"]),
                herwig_ref["sigma_unpol_pb"],
            ),
            "sigma0_vs_herwig_like_unpol": delta_payload(
                Measurement(**derived["sigma0"]),
                herwig_ref["sigma_unpol_pb"],
            ),
            "sigma0_vs_poldis_like_unpol": delta_payload(
                Measurement(**derived["sigma0"]),
                poldis_ref["sigma_unpol_pb"],
            ),
            "sigma_LL_vs_poldis_like": delta_payload(
                Measurement(**derived["sigma_LL"]),
                poldis_ref["sigma_LL_pb"],
            ),
        }
    return payload


def fmt_measurement(value_pb: float, error_pb: float) -> str:
    return f"{value_pb:.10f} +- {error_pb:.10f} pb"


def fmt_scalar(value: float) -> str:
    return f"{value:.12f} pb"


def fmt_payload_measurement(payload: Mapping[str, float]) -> str:
    return fmt_measurement(payload["value_pb"], payload["error_pb"])


def render_window_section(title: str, payload: dict) -> List[str]:
    derived = payload["herwig"]["derived"]
    comparisons = payload["comparisons"]
    poldis_ref = payload["standalone"]["poldis_like"]
    herwig_ref = payload["standalone"]["herwig_like"]
    return [
        f"{title}",
        "-" * len(title),
        f"Herwig 00:        {fmt_payload_measurement(derived['00'])}",
        f"Herwig sigma0:    {fmt_payload_measurement(derived['sigma0'])}",
        f"Herwig sigma_LL:  {fmt_payload_measurement(derived['sigma_LL'])}",
        f"Herwig sigma0-00: {fmt_payload_measurement(derived['sigma0_minus_00'])}",
        "",
        f"POLDIS-like standalone unpolarized: {fmt_scalar(poldis_ref['sigma_unpol_pb'])}",
        f"POLDIS-like standalone sigma_LL:    {fmt_scalar(poldis_ref['sigma_LL_pb'])}",
        f"Herwig-like standalone unpolarized: {fmt_scalar(herwig_ref['sigma_unpol_pb'])}",
        f"Herwig-like standalone sigma_LL:    {fmt_scalar(herwig_ref['sigma_LL_pb'])}",
        "",
        "Comparisons",
        f"  00 - herwig_like_unpol:      {fmt_payload_measurement(comparisons['00_vs_herwig_like_unpol'])}",
        f"  sigma0 - herwig_like_unpol:  {fmt_payload_measurement(comparisons['sigma0_vs_herwig_like_unpol'])}",
        f"  sigma0 - poldis_like_unpol:  {fmt_payload_measurement(comparisons['sigma0_vs_poldis_like_unpol'])}",
        f"  sigma_LL - poldis_like_LL:   {fmt_payload_measurement(comparisons['sigma_LL_vs_poldis_like'])}",
        "",
    ]


def render_results_text(payload: dict) -> str:
    lines = [
        "Focused Z LO sigma0 vs 00 study",
        "================================",
        f"Tag: {payload['tag']}",
        f"Broad control tag: {payload['control_tag']}",
        "",
        f"Broad window:    Q^2 in [{payload['broad_window']['q2_min']:.0f}, {payload['broad_window']['q2_max']:.0f}] GeV^2, "
        f"y in [{payload['broad_window']['y_min']:.1f}, {payload['broad_window']['y_max']:.1f}]",
        f"Interior window: Q^2 in [{payload['interior_window']['q2_min']:.0f}, {payload['interior_window']['q2_max']:.0f}] GeV^2, "
        f"y in [{payload['interior_window']['y_min']:.1f}, {payload['interior_window']['y_max']:.1f}]",
        "",
    ]
    lines.extend(render_window_section("Broad control (plain45)", payload["broad_control"]))
    lines.extend(render_window_section("Interior run", payload["interior_run"]))
    broad_sigma0_minus_00 = payload["broad_control"]["herwig"]["derived"]["sigma0_minus_00"]
    interior_sigma0_minus_00 = payload["interior_run"]["herwig"]["derived"]["sigma0_minus_00"]
    lines.extend(
        [
            "Headline comparison",
            "-------------------",
            f"Broad sigma0 - 00:    {fmt_payload_measurement(broad_sigma0_minus_00)}",
            f"Interior sigma0 - 00: {fmt_payload_measurement(interior_sigma0_minus_00)}",
        ]
    )
    return "\n".join(lines).rstrip() + "\n"


def write_results(base_dir: Path, tag: str, payload: dict) -> None:
    campaign_dir = base_dir / "campaigns" / tag
    campaign_dir.mkdir(parents=True, exist_ok=True)
    text_path = campaign_dir / "results.txt"
    json_path = campaign_dir / "results.json"
    rendered = render_results_text(payload)
    text_path.write_text(rendered)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(rendered, end="")
    print(f"Wrote text report: {text_path}")
    print(f"Wrote JSON report: {json_path}")


def collect_results(base_dir: Path, tag: str, control_tag: str) -> dict:
    broad_runs = collect_helicity_measurements(
        base_dir,
        build_requested_names_for_broad(),
        control_tag,
        context="broad-window plain45 control",
    )
    interior_runs = collect_explicit_sharded_measurements(
        base_dir,
        build_generated_stems_for_interior(),
        tag,
        context="interior-window generated run",
    )
    broad_summary = compute_window_summary(broad_runs)
    interior_summary = compute_window_summary(interior_runs)
    return build_results_payload(base_dir, tag, control_tag, broad_summary, interior_summary)


def print_plan(base_dir: Path, tag: str, runs: Sequence[LogicalRun], shards: Sequence[ShardSpec]) -> None:
    print(f"[plan] base_dir={base_dir}")
    print(f"[plan] logical_runs={len(runs)} shards={len(shards)}")
    print(f"[plan] monitor_txt={status_txt_path(base_dir, tag)}")
    print(f"[plan] monitor_json={status_json_path(base_dir, tag)}")
    for run in runs:
        print(f"[plan] {run.generated_stem}: source={run.source_card} generated={run.generated_input_rel} run={run.run_file}")
    for spec in shards[: min(10, len(shards))]:
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


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    base_dir = args.base_dir.resolve()
    if not base_dir.exists():
        raise SystemExit(f"Base directory does not exist: {base_dir}")

    runs = materialize_generated_cards(base_dir, args.tag, INTERIOR_WINDOW)
    shards = build_shards(runs, args.tag, args.shards, args.events_per_shard, args.seed_base)

    if args.dry_run:
        print_plan(base_dir, args.tag, runs, shards)
        ensure_run_files(base_dir, runs, dry_run=True)
        return 0

    # Validate the broad control before launching new jobs.
    collect_helicity_measurements(
        base_dir,
        build_requested_names_for_broad(),
        args.control_tag,
        context="broad-window plain45 control",
    )

    if args.collect_only:
        payload = collect_results(base_dir, args.tag, args.control_tag)
        write_results(base_dir, args.tag, payload)
        return 0

    ensure_run_files(base_dir, runs, dry_run=False)
    run_shards(base_dir, args.tag, args.control_tag, shards, args.jobs)
    payload = collect_results(base_dir, args.tag, args.control_tag)
    write_results(base_dir, args.tag, payload)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
