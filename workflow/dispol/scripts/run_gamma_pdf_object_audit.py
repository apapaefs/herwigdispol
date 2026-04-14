#!/usr/bin/env python3.10
"""
Run a focused GAMMA NLO PDF/extractor object audit.

The helper:
  * generates broad-window PP/PM POSNLO/NEGNLO GAMMA cards with audit-only
    diagnostics enabled
  * prepares missing or stale .run files with `Herwig read`
  * launches the runs
  * parses NLO_AUDIT_OBJ and NLO_AUDIT_PDF records from the real Herwig logs
  * writes a compact results.txt/json summary
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
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


DEFAULT_BASE_DIR = Path(__file__).resolve().parent
DEFAULT_TAG = "plain50-gamma-pdf-audit"
DEFAULT_JOBS = 4
DEFAULT_SHARDS = 1
DEFAULT_EVENTS_PER_SHARD = 1_000_000
DEFAULT_SEED_BASE = 720_000

MODE = "pdf-audit"
GAMMA_HELICITIES = ("PP", "PM")
NLO_PIECES = ("POSNLO", "NEGNLO")
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

DIAGNOSTIC_SETTINGS = {
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:DISDiagnostics On",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOAuditDiagnostics ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOAuditDiagnostics Yes",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOTermDiagnostics ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:DumpNLOTermDiagnostics No",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditInitialSamples ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditInitialSamples 200",
    "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditSamplePeriod ": "set /Herwig/MatrixElements/PowhegMEDISNCPol:NLOAuditSamplePeriod 5000",
}


@dataclass(frozen=True)
class LogicalRun:
    key: str
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
        help="Directory containing the DISPOL cards and outputs.",
    )
    parser.add_argument("--jobs", type=int, default=DEFAULT_JOBS, help="Maximum concurrent Herwig processes.")
    parser.add_argument("--shards", type=int, default=DEFAULT_SHARDS, help="Shards per logical run.")
    parser.add_argument(
        "--events-per-shard",
        type=int,
        default=DEFAULT_EVENTS_PER_SHARD,
        help="Events per shard.",
    )
    parser.add_argument("--seed-base", type=int, default=DEFAULT_SEED_BASE, help="Base RNG seed.")
    parser.add_argument(
        "--pdf-profile",
        choices=sorted(PDF_PROFILES),
        default="hybrid",
        help="PDF profile to materialize into the generated audit cards.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print planned actions without executing Herwig.")
    parser.add_argument(
        "--collect-only",
        action="store_true",
        help="Skip Herwig read/run and rebuild the audit summary from existing logs.",
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


def build_gamma_requested_name(piece: str, helicity: str) -> str:
    return f"DIS-POL-POWHEG_{helicity}-{piece}-GAMMA.out"


def build_gamma_source_card(piece: str, helicity: str) -> str:
    return build_gamma_requested_name(piece, helicity).replace(".out", ".in")


def stem_without_suffix(filename: str) -> str:
    return filename[:-3] if filename.endswith(".in") else filename


def fmt_seconds(value: float | None) -> str:
    if value is None or not math.isfinite(value):
        return "n/a"
    total = int(round(value))
    hours, rem = divmod(total, 3600)
    minutes, seconds = divmod(rem, 60)
    return f"{hours:d}:{minutes:02d}:{seconds:02d}"


def build_runs(tag: str, pdf_profile: str) -> List[LogicalRun]:
    runs: List[LogicalRun] = []
    profile_token = "HYB" if pdf_profile == "hybrid" else "NNPDF"
    for piece in NLO_PIECES:
        for helicity in GAMMA_HELICITIES:
            source_card = build_gamma_source_card(piece, helicity)
            source_stem = stem_without_suffix(source_card)
            generated_stem = f"{source_stem}-PDFAUDIT-{profile_token}"
            key = f"{pdf_profile}:{piece}_{helicity}"
            runs.append(
                LogicalRun(
                    key=key,
                    piece=piece,
                    helicity=helicity,
                    source_card=source_card,
                    requested_name=build_gamma_requested_name(piece, helicity),
                    generated_stem=generated_stem,
                    generated_input_rel=Path("campaigns") / tag / MODE / "generated-inputs" / f"{generated_stem}.in",
                    run_file=f"{generated_stem}.run",
                )
            )
    return runs


def patch_card_text(source_text: str, generated_stem: str, pdf_profile: str) -> str:
    pdfs = PDF_PROFILES[pdf_profile]
    replacements = {
        "set /Herwig/Partons/DiffPDF:PDFName ": f"set /Herwig/Partons/DiffPDF:PDFName {pdfs['polarized_diff']}",
        "set /Herwig/Partons/LHAPDF:PDFName ": f"set /Herwig/Partons/LHAPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/HardLOPDF:PDFName ": f"set /Herwig/Partons/HardLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/ShowerNLOPDF:PDFName ": f"set /Herwig/Partons/ShowerNLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/HardNLOPDF:PDFName ": f"set /Herwig/Partons/HardNLOPDF:PDFName {pdfs['unpolarized']}",
        "set /Herwig/Partons/ShowerLOPDF:PDFName ": f"set /Herwig/Partons/ShowerLOPDF:PDFName {pdfs['unpolarized']}",
    }
    lines = source_text.splitlines()
    patched: List[str] = []
    saverun_found = False
    use_q2_scale_found = False
    diag_seen = {prefix: False for prefix in DIAGNOSTIC_SETTINGS}
    replacement_seen = {prefix: False for prefix in replacements}

    for line in lines:
        stripped = line.strip()

        replaced = False
        for prefix, replacement in replacements.items():
            if stripped.startswith(prefix):
                patched.append(replacement)
                replacement_seen[prefix] = True
                replaced = True
                break
        if replaced:
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

    missing_replacements = [prefix.strip() for prefix, seen in replacement_seen.items() if not seen]
    if missing_replacements:
        raise RuntimeError(f"Failed to rewrite PDF settings: {missing_replacements}")
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


def materialize_cards(base_dir: Path, runs: Sequence[LogicalRun], pdf_profile: str) -> None:
    for run in runs:
        source_path = base_dir / run.source_card
        if not source_path.exists():
            raise FileNotFoundError(f"Missing source card {source_path}")
        target_path = base_dir / run.generated_input_rel
        rendered = patch_card_text(source_path.read_text(), run.generated_stem, pdf_profile)
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
    write_monitor_files(
        base_dir,
        tag,
        build_monitor_payload(base_dir, tag, shards, results, started, phase="running"),
    )

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
            if len(results) == 1 or len(results) == total or result.returncode != 0:
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


def print_plan(base_dir: Path, tag: str, runs: Sequence[LogicalRun], shards: Sequence[ShardSpec], pdf_profile: str) -> None:
    print(f"[plan:{MODE}] base_dir={base_dir}")
    print(f"[plan:{MODE}] pdf_profile={pdf_profile}")
    print(f"[plan:{MODE}] logical_runs={len(runs)} shards={len(shards)}")
    print(f"[plan:{MODE}] monitor_txt={status_txt_path(base_dir, tag)}")
    print(f"[plan:{MODE}] monitor_json={status_json_path(base_dir, tag)}")
    for run in runs:
        print(
            f"[plan:{MODE}] {run.key}: source={run.source_card} generated={run.generated_input_rel} run={run.run_file}"
        )
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


def parse_bool_token(value: str) -> bool:
    token = value.strip()
    if token in {"1", "true", "True", "YES", "Yes"}:
        return True
    if token in {"0", "false", "False", "NO", "No"}:
        return False
    raise ValueError(f"Cannot parse boolean token {value!r}")


def parse_audit_line(prefix: str, line: str) -> Optional[dict]:
    if not line.startswith(prefix + " "):
        return None
    payload: Dict[str, str] = {}
    for token in line[len(prefix) + 1 :].strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        payload[key] = value
    required = {"run", "hel", "contrib", "n"}
    if not required.issubset(payload):
        return None
    payload["_kind"] = prefix
    return payload


def audit_key(payload: Mapping[str, str]) -> Tuple[str, str, str, int]:
    return (
        payload["run"],
        payload["hel"],
        payload["contrib"],
        int(payload["n"]),
    )


def find_runtime_logs(base_dir: Path, generated_stem: str, tag: str) -> List[Path]:
    pattern = f"{generated_stem}-S*-{tag}-s*.log"
    return sorted(
        path
        for path in base_dir.rglob(pattern)
        if path.name.endswith(".log") and not path.name.endswith("-EvtGen.log")
    )


def collect_runtime_logs(base_dir: Path, runs: Sequence[LogicalRun], tag: str) -> Dict[str, List[str]]:
    out: Dict[str, List[str]] = {}
    for run in runs:
        paths = find_runtime_logs(base_dir, run.generated_stem, tag)
        if not paths:
            raise FileNotFoundError(
                f"Missing runtime logs for {run.key}: no files matching {run.generated_stem}-S*-{tag}-s*.log"
            )
        out[run.key] = [str(path) for path in paths]
    return out


def load_audit_records(base_dir: Path, runs: Sequence[LogicalRun], tag: str) -> tuple[list[dict], Dict[str, List[str]]]:
    by_key: Dict[Tuple[str, str, str, int], dict] = {}
    runtime_logs = collect_runtime_logs(base_dir, runs, tag)
    for run in runs:
        for log_path_str in runtime_logs[run.key]:
            log_path = Path(log_path_str)
            for raw_line in log_path.read_text().splitlines():
                parsed = parse_audit_line("NLO_AUDIT_OBJ", raw_line)
                if parsed is None:
                    parsed = parse_audit_line("NLO_AUDIT_PDF", raw_line)
                if parsed is None:
                    continue
                key = audit_key(parsed)
                record = by_key.setdefault(
                    key,
                    {
                        "run": parsed["run"],
                        "hel": parsed["hel"],
                        "contrib": parsed["contrib"],
                        "n": int(parsed["n"]),
                        "logs": [],
                    },
                )
                if log_path_str not in record["logs"]:
                    record["logs"].append(log_path_str)
                if parsed["_kind"] == "NLO_AUDIT_OBJ":
                    record["obj"] = {
                        "lastExtractor": parsed.get("lastExtractor", "NULL"),
                        "eventExtractor": parsed.get("eventExtractor", "NULL"),
                        "sameExtractor": parse_bool_token(parsed.get("sameExtractor", "0")),
                        "lastIsPPE": parse_bool_token(parsed.get("lastIsPPE", "0")),
                        "eventIsPPE": parse_bool_token(parsed.get("eventIsPPE", "0")),
                        "beamPdf": parsed.get("beamPdf", "NULL"),
                        "extractorPdf": parsed.get("extractorPdf", "NULL"),
                        "diffPdf": parsed.get("diffPdf", "NULL"),
                    }
                else:
                    record["pdf"] = {
                        "xB": float(parsed["xB"]),
                        "xp": float(parsed["xp"]),
                        "Q2": float(parsed["Q2"]),
                        "mu2": float(parsed["mu2"]),
                        "jac": float(parsed["jac"]),
                        "Pz": float(parsed["Pz"]),
                        "hasDiffPdf": parse_bool_token(parsed.get("hasDiffPdf", "0")),
                    }
    return sorted(by_key.values(), key=lambda item: (item["run"], item["n"])), runtime_logs


def sample_payload(record: Mapping[str, object]) -> dict:
    obj = record.get("obj", {})
    pdf = record.get("pdf", {})
    assert isinstance(obj, Mapping)
    assert isinstance(pdf, Mapping)
    return {
        "run": record["run"],
        "hel": record["hel"],
        "contrib": record["contrib"],
        "n": record["n"],
        "lastExtractor": obj.get("lastExtractor"),
        "eventExtractor": obj.get("eventExtractor"),
        "beamPdf": obj.get("beamPdf"),
        "extractorPdf": obj.get("extractorPdf"),
        "diffPdf": obj.get("diffPdf"),
        "sameExtractor": obj.get("sameExtractor"),
        "hasDiffPdf": pdf.get("hasDiffPdf"),
        "xB": pdf.get("xB"),
        "xp": pdf.get("xp"),
        "Q2": pdf.get("Q2"),
        "mu2": pdf.get("mu2"),
    }


def summarize_records(records: Sequence[dict], runtime_logs: Mapping[str, List[str]]) -> dict:
    complete = [record for record in records if "obj" in record and "pdf" in record]
    obj_only = sum(1 for record in records if "obj" in record and "pdf" not in record)
    pdf_only = sum(1 for record in records if "pdf" in record and "obj" not in record)

    same_extractor_counts = Counter()
    extractor_pairs = Counter()
    pdf_triples = Counter()
    beam_pdf_mismatch_samples: List[dict] = []
    extractor_mismatch_samples: List[dict] = []
    diff_pdf_missing_samples: List[dict] = []

    for record in complete:
        obj = record["obj"]
        pdf = record["pdf"]
        assert isinstance(obj, Mapping)
        assert isinstance(pdf, Mapping)
        same_extractor = bool(obj["sameExtractor"])
        same_extractor_counts["true" if same_extractor else "false"] += 1
        extractor_pairs[(str(obj["lastExtractor"]), str(obj["eventExtractor"]))] += 1
        pdf_triples[(str(obj["beamPdf"]), str(obj["extractorPdf"]), str(obj["diffPdf"]))] += 1

        if str(obj["beamPdf"]) != str(obj["extractorPdf"]) and len(beam_pdf_mismatch_samples) < 5:
            beam_pdf_mismatch_samples.append(sample_payload(record))
        if not same_extractor and len(extractor_mismatch_samples) < 5:
            extractor_mismatch_samples.append(sample_payload(record))
        if not bool(pdf["hasDiffPdf"]) and len(diff_pdf_missing_samples) < 5:
            diff_pdf_missing_samples.append(sample_payload(record))

    return {
        "record_counts": {
            "total_joined_records": len(records),
            "complete_records": len(complete),
            "obj_only_records": obj_only,
            "pdf_only_records": pdf_only,
        },
        "same_extractor_counts": dict(same_extractor_counts),
        "distinct_extractor_pairs": [
            {"count": count, "lastExtractor": pair[0], "eventExtractor": pair[1]}
            for pair, count in extractor_pairs.most_common()
        ],
        "distinct_pdf_triples": [
            {"count": count, "beamPdf": triple[0], "extractorPdf": triple[1], "diffPdf": triple[2]}
            for triple, count in pdf_triples.most_common()
        ],
        "flags": {
            "beamPdf_ne_extractorPdf": bool(beam_pdf_mismatch_samples),
            "sameExtractor_false": bool(extractor_mismatch_samples),
            "hasDiffPdf_false": bool(diff_pdf_missing_samples),
        },
        "samples": {
            "beamPdf_ne_extractorPdf": beam_pdf_mismatch_samples,
            "sameExtractor_false": extractor_mismatch_samples,
            "hasDiffPdf_false": diff_pdf_missing_samples,
        },
        "runtime_logs": runtime_logs,
    }


def render_pair_lines(title: str, rows: Sequence[Mapping[str, object]], lhs: str, rhs: str, limit: int = 10) -> List[str]:
    lines = [title, "-" * len(title)]
    if not rows:
        lines.append("none")
        return lines
    for row in rows[:limit]:
        lines.append(
            f"{int(row['count']):>6}  {lhs}={row[lhs]}  {rhs}={row[rhs]}"
        )
    return lines


def render_triple_lines(title: str, rows: Sequence[Mapping[str, object]], limit: int = 10) -> List[str]:
    lines = [title, "-" * len(title)]
    if not rows:
        lines.append("none")
        return lines
    for row in rows[:limit]:
        lines.append(
            f"{int(row['count']):>6}  beamPdf={row['beamPdf']}  extractorPdf={row['extractorPdf']}  diffPdf={row['diffPdf']}"
        )
    return lines


def render_sample_lines(title: str, rows: Sequence[Mapping[str, object]]) -> List[str]:
    lines = [title, "-" * len(title)]
    if not rows:
        lines.append("none")
        return lines
    for row in rows:
        lines.append(
            f"run={row['run']} hel={row['hel']} contrib={row['contrib']} n={row['n']} "
            f"xB={row['xB']} xp={row['xp']} Q2={row['Q2']} mu2={row['mu2']}"
        )
        lines.append(
            f"  lastExtractor={row['lastExtractor']} eventExtractor={row['eventExtractor']}"
        )
        lines.append(
            f"  beamPdf={row['beamPdf']} extractorPdf={row['extractorPdf']} diffPdf={row['diffPdf']}"
        )
    return lines


def render_results_text(payload: dict) -> str:
    counts = payload["record_counts"]
    flags = payload["flags"]
    same_counts = payload["same_extractor_counts"]
    lines = [
        "GAMMA PDF object audit",
        "======================",
        f"Tag: {payload['tag']}",
        f"PDF profile: {payload['pdf_profile']}",
        f"Unpolarized PDF: {payload['pdfs']['unpolarized']}",
        f"Polarized diff PDF: {payload['pdfs']['polarized_diff']}",
        "",
        "Record counts",
        "-------------",
        f"joined records:   {counts['total_joined_records']}",
        f"complete records: {counts['complete_records']}",
        f"obj-only records: {counts['obj_only_records']}",
        f"pdf-only records: {counts['pdf_only_records']}",
        "",
        "Flags",
        "-----",
        f"beamPdf != extractorPdf ever: {'Yes' if flags['beamPdf_ne_extractorPdf'] else 'No'}",
        f"sameExtractor = false ever:   {'Yes' if flags['sameExtractor_false'] else 'No'}",
        f"hasDiffPdf = false ever:      {'Yes' if flags['hasDiffPdf_false'] else 'No'}",
        "",
        "sameExtractor counts",
        "--------------------",
        f"true:  {same_counts.get('true', 0)}",
        f"false: {same_counts.get('false', 0)}",
        "",
    ]
    lines.extend(render_pair_lines("Distinct extractor pairs", payload["distinct_extractor_pairs"], "lastExtractor", "eventExtractor"))
    lines.append("")
    lines.extend(render_triple_lines("Distinct beam/extractor/diff PDF triples", payload["distinct_pdf_triples"]))
    lines.append("")
    lines.extend(render_sample_lines("beamPdf != extractorPdf samples", payload["samples"]["beamPdf_ne_extractorPdf"]))
    lines.append("")
    lines.extend(render_sample_lines("sameExtractor = false samples", payload["samples"]["sameExtractor_false"]))
    lines.append("")
    lines.extend(render_sample_lines("hasDiffPdf = false samples", payload["samples"]["hasDiffPdf_false"]))
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


def run(args: argparse.Namespace) -> None:
    base_dir = args.base_dir.resolve()
    tag = args.tag
    runs = build_runs(tag, args.pdf_profile)
    materialize_cards(base_dir, runs, args.pdf_profile)
    shards = build_shards(runs, tag, args.shards, args.events_per_shard, args.seed_base)

    if args.dry_run:
        print_plan(base_dir, tag, runs, shards, args.pdf_profile)
        ensure_run_files(base_dir, runs, dry_run=True)
        return

    if not args.collect_only:
        ensure_run_files(base_dir, runs, dry_run=False)
        run_shards(base_dir, tag, shards, args.jobs)

    records, runtime_logs = load_audit_records(base_dir, runs, tag)
    summary = summarize_records(records, runtime_logs)
    payload = {
        "tag": tag,
        "mode": MODE,
        "pdf_profile": args.pdf_profile,
        "pdfs": dict(PDF_PROFILES[args.pdf_profile]),
        **summary,
    }
    write_results(base_dir, tag, payload)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.base_dir.resolve().exists():
        raise SystemExit(f"Base directory does not exist: {args.base_dir}")
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
