#!/usr/bin/env python3.10
"""
Run multi-shard polarized GAMMA TERMDIAG comparisons across SamplingPower values.

The script:
  * generates temporary TERMDIAG cards under campaigns/<tag>/generated-inputs/
  * prepares stale or missing .run files with `Herwig read`
  * launches Herwig shards
  * maintains monitor/status.txt and monitor/status.json
  * extracts aggregated NLO term diagnostics per power
  * writes one comparison report summarizing which components move with SamplingPower
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
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence

from extract_dis_out_results import (
    Measurement,
    combine_measurements,
    parse_generated_event_count,
    parse_measurement,
)


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_TAG = "plain49-gamma-termdiag-power"
DEFAULT_JOBS = 8
DEFAULT_SHARDS = 10
DEFAULT_EVENTS_PER_SHARD = 1_000_000
DEFAULT_SEED_BASE = 710_000
DEFAULT_POWERS = "0.4,0.6"

MODE = "gamma-termdiag-power"
GAMMA_HELICITIES = ("PP", "PM")
NLO_PIECES = ("POSNLO", "NEGNLO")
SAMPLING_POWER_OBJECTS = ("MEDISNC", "MEDISNCPol", "PowhegMEDISNC", "PowhegMEDISNCPol")
PRIMARY_COMPONENTS = ("F_virt", "F_cq_even", "F_cq_odd", "F_cg_even", "F_cg_odd", "F_rq", "F_rg")
REAL_SPLIT_COMPONENTS = ("F_rq_even", "F_rq_odd", "F_rg_even", "F_rg_odd")
COMPONENT_ORDER = PRIMARY_COMPONENTS + REAL_SPLIT_COMPONENTS

DIAGNOSTIC_SETTINGS = {
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics On",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOAuditDiagnostics ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOAuditDiagnostics Yes",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOTermDiagnostics ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOTermDiagnostics Yes",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditInitialSamples ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditInitialSamples 50",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditSamplePeriod ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditSamplePeriod 0",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOTermDiagnosticPeriod ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOTermDiagnosticPeriod 10",
}


@dataclass(frozen=True)
class LogicalRun:
    key: str
    power: float
    piece: str
    helicity: str
    source_card: str
    requested_name: str
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


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tag", default=DEFAULT_TAG, help="Campaign tag used in Herwig shard tags.")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing DISPOL cards and outputs.",
    )
    parser.add_argument("--powers", default=DEFAULT_POWERS, help="Comma-separated SamplingPower values.")
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Maximum concurrent Herwig processes.")
    parser.add_argument("--shards", type=int, default=DEFAULT_SHARDS, help="Shards per logical run.")
    parser.add_argument(
        "--events-per-shard",
        type=int,
        default=DEFAULT_EVENTS_PER_SHARD,
        help="Events per shard.",
    )
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE, help="Base RNG seed.")
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without executing Herwig.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip Herwig read/run and rebuild reports from existing outputs.",
    )
    return parser


def campaign_dir(base_dir: Path, tag: str) -> Path:
    return base_dir / "campaigns" / tag


def generated_input_dir(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "generated-inputs"


def launcher_log_dir(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "launcher-logs"


def monitor_dir(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "monitor"


def status_json_path(base_dir: Path, tag: str) -> Path:
    return monitor_dir(base_dir, tag) / "status.json"


def status_txt_path(base_dir: Path, tag: str) -> Path:
    return monitor_dir(base_dir, tag) / "status.txt"


def results_txt_path(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "results.txt"


def results_json_path(base_dir: Path, tag: str) -> Path:
    return campaign_dir(base_dir, tag) / "results.json"


def extractor_text_path(base_dir: Path, tag: str, power: float) -> Path:
    return campaign_dir(base_dir, tag) / f"extract-{power_suffix(power).lower()}.txt"


def extractor_json_path(base_dir: Path, tag: str, power: float) -> Path:
    return campaign_dir(base_dir, tag) / f"extract-{power_suffix(power).lower()}.json"


def extractor_csv_path(base_dir: Path, tag: str, power: float) -> Path:
    return campaign_dir(base_dir, tag) / f"extract-{power_suffix(power).lower()}.csv"


def build_gamma_requested_name(piece: str, helicity: str) -> str:
    return f"DIS-POL-POWHEG_{helicity}-{piece}-GAMMA.out"


def build_gamma_source_card(piece: str, helicity: str) -> str:
    return build_gamma_requested_name(piece, helicity).replace(".out", ".in")


def stem_without_suffix(filename: str) -> str:
    return filename[:-3] if filename.endswith(".in") else filename


def fmt_power_token(power: float) -> str:
    return f"{power:.6g}"


def power_suffix(power: float) -> str:
    return f"SP{int(round(power * 100.0)):03d}"


def parse_power_list(text: str) -> List[float]:
    out: List[float] = []
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
        out.append(value)
    if not out:
        raise ValueError("Need at least one SamplingPower value.")
    return out


def sampling_power_block(power: float) -> List[str]:
    return [f"set {obj}:SamplingPower {fmt_power_token(power)}" for obj in SAMPLING_POWER_OBJECTS]


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def build_runs(tag: str, powers: Sequence[float]) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    for power in powers:
        suffix = power_suffix(power)
        for piece in NLO_PIECES:
            for helicity in GAMMA_HELICITIES:
                source_card = build_gamma_source_card(piece, helicity)
                source_stem = stem_without_suffix(source_card)
                generated_stem = f"{source_stem}-TERMDIAG-{suffix}"
                key = f"{fmt_power_token(power)}::{piece}_{helicity}"
                runs.append(
                    LogicalRun(
                        key=key,
                        power=power,
                        piece=piece,
                        helicity=helicity,
                        source_card=source_card,
                        requested_name=build_gamma_requested_name(piece, helicity),
                        generated_stem=generated_stem,
                        generated_input_rel=Path("campaigns") / tag / "generated-inputs" / f"{generated_stem}.in",
                        run_file=f"{generated_stem}.run",
                    )
                )
    return runs


def patch_card_text(source_text: str, generated_stem: str, power: float) -> str:
    sampling_re = re.compile(r"^\s*#?\s*set\s+(\S+):SamplingPower\b")
    lines = source_text.splitlines()
    patched: List[str] = []
    saverun_found = False
    use_q2_scale_found = False
    sampling_replacements = 0
    inserted_sampling_block = False
    diag_seen = {prefix: False for prefix in DIAGNOSTIC_SETTINGS}

    for line in lines:
        stripped = line.strip()

        match = sampling_re.match(line)
        if match and match.group(1) in SAMPLING_POWER_OBJECTS:
            patched.append(f"set {match.group(1)}:SamplingPower {fmt_power_token(power)}")
            sampling_replacements += 1
            continue

        diag_replaced = False
        for prefix, replacement in DIAGNOSTIC_SETTINGS.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                diag_seen[prefix] = True
                diag_replaced = True
                break
        if diag_replaced:
            continue

        if stripped == "cd /Herwig/MatrixElements/":
            patched.append(line)
            if "SamplingPower" not in source_text and not inserted_sampling_block:
                patched.extend(sampling_power_block(power))
                inserted_sampling_block = True
            continue

        if stripped.startswith("set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission "):
            patched.append(line)
            use_q2_scale_found = True
            for prefix, replacement in DIAGNOSTIC_SETTINGS.items():
                if not diag_seen[prefix]:
                    patched.append(replacement)
                    diag_seen[prefix] = True
            continue

        if stripped.startswith("saverun ") and stripped.endswith(" EventGenerator"):
            patched.append(f"saverun {generated_stem} EventGenerator")
            saverun_found = True
            continue

        patched.append(line)

    if "SamplingPower" in source_text and sampling_replacements == 0:
        raise RuntimeError("Expected SamplingPower lines in source card but did not rewrite any.")
    if "SamplingPower" not in source_text and not inserted_sampling_block:
        raise RuntimeError("Failed to insert SamplingPower block into generated card.")
    if not use_q2_scale_found:
        raise RuntimeError("Failed to find UseQ2ScaleInPOWHEGEmission insertion point.")
    if not all(diag_seen.values()):
        missing = [prefix for prefix, seen in diag_seen.items() if not seen]
        raise RuntimeError(f"Failed to insert diagnostic lines: {missing}")
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
        rendered = patch_card_text(source_path.read_text(), run.generated_stem, run.power)
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


def measurement_payload(measurement: Measurement) -> dict:
    return {"value_pb": measurement.value_pb, "error_pb": measurement.error_pb}


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


def run_extractor(base_dir: Path, tag: str, power: float, dry_run: bool) -> dict | None:
    script_path = base_dir / "extract_nlo_term_diagnostics.py"
    power_token = power_suffix(power)
    text_path = extractor_text_path(base_dir, tag, power)
    json_path = extractor_json_path(base_dir, tag, power)
    csv_path = extractor_csv_path(base_dir, tag, power)
    cmd = [
        "python3.10",
        str(script_path),
        "--base-dir",
        str(base_dir),
        "--tag",
        tag,
        "--setup",
        "GAMMA",
        "--piece",
        "POSNLO",
        "--piece",
        "NEGNLO",
        "--helicity",
        "PP",
        "--helicity",
        "PM",
        "--name-filter",
        power_token,
        "--json-out",
        str(json_path),
        "--csv-out",
        str(csv_path),
    ]
    if dry_run:
        print(shlex.join(cmd))
        return None

    proc = subprocess.run(cmd, cwd=base_dir, text=True, capture_output=True)
    if proc.returncode != 0:
        raise RuntimeError(
            f"extract_nlo_term_diagnostics.py failed for power {fmt_power_token(power)}\n"
            f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
        )
    text_path.parent.mkdir(parents=True, exist_ok=True)
    output_text = proc.stdout
    if proc.stderr:
        output_text += "\n[stderr]\n" + proc.stderr
    text_path.write_text(output_text)
    print(proc.stdout, end="")
    payload = json.loads(json_path.read_text())
    selected_runs = [
        entry
        for entry in payload.get("selected_runs", [])
        if entry.get("setup") == "GAMMA"
    ]
    if not selected_runs:
        raise RuntimeError(
            f"No GAMMA NLO_TERM diagnostics were found for power {fmt_power_token(power)} "
            f"(name-filter {power_token}). Check the run logs under {base_dir}."
        )
    return payload


def extract_component_block(payload: dict, top_key: str, piece: str, setup: str = "GAMMA") -> Dict[str, Optional[float]]:
    section = payload.get(top_key, {}).get(setup, {}).get(piece, {})
    out: Dict[str, Optional[float]] = {}
    for component, value_map in section.items():
        if not isinstance(value_map, Mapping):
            continue
        out[component] = value_map.get("pol_pb")
    return out


def build_component_payload(extractor_payload: dict) -> Dict[str, Dict[str, Optional[float]]]:
    out: Dict[str, Dict[str, Optional[float]]] = {}
    for piece in ("POSNLO", "NEGNLO", "NLO"):
        piece_payload: Dict[str, Optional[float]] = {}
        primary = extract_component_block(extractor_payload, "primary_estimates", piece)
        split = extract_component_block(extractor_payload, "real_split_estimates", piece)
        for component in PRIMARY_COMPONENTS:
            piece_payload[component] = primary.get(component)
        for component in REAL_SPLIT_COMPONENTS:
            piece_payload[component] = split.get(component)
        out[piece] = piece_payload
    return out


def build_closure_payload(extractor_payload: dict) -> Dict[str, dict]:
    out: Dict[str, dict] = {}
    for entry in extractor_payload.get("selected_runs", []):
        if entry.get("setup") != "GAMMA":
            continue
        key = f"{entry.get('piece')}_{entry.get('helicity')}"
        values = entry.get("values", {})
        out[key] = {
            "B_closure": values.get("B_closure"),
            "M_closure": values.get("M_closure"),
            "ME_self": values.get("ME_self"),
        }
    return out


def fmt_pb(value: Optional[float]) -> str:
    return "missing" if value is None else f"{value:+.6f} pb"


def fmt_unitless(value: Optional[float]) -> str:
    return "missing" if value is None else f"{value:+.6f}"


def build_power_payload(base_dir: Path, tag: str, power: float, runs: Sequence[LogicalRun], extractor_payload: dict) -> dict:
    stems = {f"{run.piece}_{run.helicity}": run.generated_stem for run in runs if round(run.power, 12) == round(power, 12)}
    flat_runs = collect_explicit_sharded_measurements(
        base_dir,
        stems,
        tag,
        context=f"GAMMA TERMDIAG power={fmt_power_token(power)}",
    )
    summary = gamma_summary_from_flat_runs(flat_runs)
    return {
        "power": power,
        "power_label": fmt_power_token(power),
        "source_stems": stems,
        "runs": summary["runs"],
        "derived": summary["derived"],
        "components": build_component_payload(extractor_payload),
        "closure": build_closure_payload(extractor_payload),
        "extractor_files": {
            "text": str(extractor_text_path(base_dir, tag, power)),
            "json": str(extractor_json_path(base_dir, tag, power)),
            "csv": str(extractor_csv_path(base_dir, tag, power)),
        },
    }


def build_comparison_payload(power_payloads: Mapping[str, dict], left_label: str, right_label: str) -> dict:
    comparison: Dict[str, dict] = {}
    for piece in ("POSNLO", "NEGNLO", "NLO"):
        rows: Dict[str, dict] = {}
        left_components = power_payloads[left_label]["components"][piece]
        right_components = power_payloads[right_label]["components"][piece]
        for component in COMPONENT_ORDER:
            left = left_components.get(component)
            right = right_components.get(component)
            rows[component] = {
                "left": left,
                "right": right,
                "delta": None if left is None or right is None else right - left,
            }
        comparison[piece] = rows
    return {
        "left_power": left_label,
        "right_power": right_label,
        "pieces": comparison,
    }


def render_component_table(component_map: Mapping[str, Optional[float]]) -> List[str]:
    width = max(len(name) for name in COMPONENT_ORDER)
    return [f"{name:<{width}}  {fmt_pb(component_map.get(name))}" for name in COMPONENT_ORDER]


def render_comparison_table(piece: str, rows: Mapping[str, dict], left_label: str, right_label: str) -> List[str]:
    lines = [
        f"{piece} component comparison",
        "-" * (len(piece) + 21),
        f"Component      {left_label:>12}  {right_label:>12}  delta({right_label}-{left_label})",
        "-------------  ------------  ------------  ----------------------",
    ]
    for component in COMPONENT_ORDER:
        row = rows[component]
        lines.append(
            f"{component:<13}  {fmt_pb(row['left']):>12}  {fmt_pb(row['right']):>12}  {fmt_pb(row['delta']):>22}"
        )
    return lines


def render_results_text(payload: dict) -> str:
    lines = [
        "Multi-shard GAMMA TERMDIAG SamplingPower comparison",
        "===================================================",
        f"Tag: {payload['tag']}",
        f"Powers: {', '.join(payload['power_labels'])}",
        f"Shards per logical run: {payload['shards']}",
        f"Events per shard: {payload['events_per_shard']}",
        "",
    ]
    for power_label in payload["power_labels"]:
        entry = payload["powers"][power_label]
        derived = entry["derived"]
        lines.extend(
            [
                f"SamplingPower {power_label}",
                "-" * (14 + len(power_label)),
                f"POSNLO_pol: {fmt_measurement(derived['POSNLO_pol'])}",
                f"NEGNLO_pol: {fmt_measurement(derived['NEGNLO_pol'])}",
                f"NLO_pol:    {fmt_measurement(derived['NLO_pol'])}",
                "",
                "Closure diagnostics",
                "-------------------",
            ]
        )
        for key in sorted(entry["closure"]):
            closure = entry["closure"][key]
            lines.append(
                f"{key}: B_closure={fmt_unitless(closure.get('B_closure'))} "
                f"M_closure={fmt_unitless(closure.get('M_closure'))} "
                f"ME_self={fmt_unitless(closure.get('ME_self'))}"
            )
        if not entry["closure"]:
            lines.append("No closure diagnostics found.")
        lines.append("")

        for piece in ("POSNLO", "NEGNLO", "NLO"):
            lines.append(f"{piece} polarized component estimates")
            lines.append("-" * (len(piece) + 31))
            lines.extend(render_component_table(entry["components"][piece]))
            lines.append("")

        files = entry["extractor_files"]
        lines.extend(
            [
                "Extractor artifacts",
                "-------------------",
                f"text: {files['text']}",
                f"json: {files['json']}",
                f"csv:  {files['csv']}",
                "",
            ]
        )

    comparison = payload.get("comparison")
    if isinstance(comparison, Mapping):
        left_label = comparison["left_power"]
        right_label = comparison["right_power"]
        lines.extend(
            [
                f"Cross-power comparison ({left_label} vs {right_label})",
                "-----------------------------------",
                "",
            ]
        )
        for piece in ("POSNLO", "NEGNLO", "NLO"):
            lines.extend(render_comparison_table(piece, comparison["pieces"][piece], left_label, right_label))
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_results(base_dir: Path, tag: str, payload: dict) -> None:
    out_dir = campaign_dir(base_dir, tag)
    out_dir.mkdir(parents=True, exist_ok=True)
    text = render_results_text(payload)
    results_txt_path(base_dir, tag).write_text(text)
    results_json_path(base_dir, tag).write_text(json.dumps(payload, indent=2, sort_keys=True))
    print(text, end="")
    print(f"Wrote text report: {results_txt_path(base_dir, tag)}")
    print(f"Wrote JSON report: {results_json_path(base_dir, tag)}")


def print_plan(base_dir: Path, tag: str, runs: Sequence[LogicalRun], shards: Sequence[ShardSpec], powers: Sequence[float]) -> None:
    print(f"[plan:{MODE}] base_dir={base_dir}")
    print(f"[plan:{MODE}] campaign_dir={campaign_dir(base_dir, tag)}")
    print(f"[plan:{MODE}] powers={', '.join(fmt_power_token(power) for power in powers)}")
    print(f"[plan:{MODE}] logical_runs={len(runs)} shards={len(shards)}")
    print(f"[plan:{MODE}] monitor_txt={status_txt_path(base_dir, tag)}")
    print(f"[plan:{MODE}] monitor_json={status_json_path(base_dir, tag)}")
    for run in runs:
        print(
            f"[plan:{MODE}] {run.key}: source={run.source_card} generated={run.generated_input_rel} run={run.run_file}"
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
    for power in powers:
        cmd = [
            "python3.10",
            str(base_dir / "extract_nlo_term_diagnostics.py"),
            "--base-dir",
            str(base_dir),
            "--tag",
            tag,
            "--setup",
            "GAMMA",
            "--piece",
            "POSNLO",
            "--piece",
            "NEGNLO",
            "--helicity",
            "PP",
            "--helicity",
            "PM",
            "--name-filter",
            power_suffix(power),
            "--json-out",
            str(extractor_json_path(base_dir, tag, power)),
            "--csv-out",
            str(extractor_csv_path(base_dir, tag, power)),
        ]
        print(shlex.join(cmd))


def run(args: argparse.Namespace) -> None:
    base_dir = args.base_dir.resolve()
    tag = args.tag
    powers = parse_power_list(args.powers)
    runs = build_runs(tag, powers)
    shards = build_shards(runs, tag, args.shards, args.events_per_shard, args.seed_base)

    materialize_cards(base_dir, runs)
    if not args.collect_only:
        ensure_run_files(base_dir, runs, args.dry_run)

    if args.dry_run:
        print_plan(base_dir, tag, runs, shards, powers)
        return

    if not args.collect_only:
        run_shards(base_dir, tag, shards, args.jobs)

    power_payloads: Dict[str, dict] = {}
    for power in powers:
        extractor_payload = run_extractor(base_dir, tag, power, dry_run=False)
        assert extractor_payload is not None
        power_payloads[fmt_power_token(power)] = build_power_payload(base_dir, tag, power, runs, extractor_payload)

    sorted_labels = [fmt_power_token(power) for power in sorted(powers)]
    payload = {
        "tag": tag,
        "mode": MODE,
        "powers": power_payloads,
        "power_labels": sorted_labels,
        "shards": args.shards,
        "events_per_shard": args.events_per_shard,
        "seed_base": args.seed_base,
    }
    if len(sorted_labels) >= 2:
        left_label = sorted_labels[0]
        right_label = sorted_labels[-1]
        payload["comparison"] = build_comparison_payload(power_payloads, left_label, right_label)

    write_results(base_dir, tag, payload)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
