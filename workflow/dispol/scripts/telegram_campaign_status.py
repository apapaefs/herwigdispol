#!/usr/bin/env python3.10
"""Send periodic DISPOL campaign status updates to Telegram.

This helper is intentionally out-of-band and read-only with respect to the
campaign workflow. It reconstructs the expected shard layout using the same
logic as `recover_campaign_manifest.py`, inspects launcher logs and produced
artifacts under `--base-dir`, and sends a compact status summary via the
Telegram Bot API.

It does not modify the running campaign, does not invoke `Herwig`, and does
not write into the workflow tree unless the user explicitly enables an
external state file via `--only-on-change`.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import socket
import subprocess
import sys
import tempfile
import time
import urllib.parse
import urllib.request
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Sequence

DEFAULT_BASE_DIR = Path(__file__).resolve().parent

from recover_campaign_manifest import (
    SUPPORTED_SETUPS,
    build_artifact_index,
    build_jobs,
    build_launcher_log_index,
    build_shards,
    normalize_campaign_setups,
    resolve_analysis_configuration,
)


PROGRESS_MARKER_RE = re.compile(r"event>\s+(?P<current>init|\d+)(?:\s+(?P<total>\d+)|/(?P<total_alt>\d+))")
SCREEN_SESSION_RE = re.compile(r"^\s*(?P<name>\d+\.[^\s]+)\s+\((?P<state>[^)]+)\)\s*$")
SCREEN_TAG_RE = re.compile(r"^Tag:\s*(?P<tag>.+)$")
SCREEN_ELAPSED_RE = re.compile(r"^Elapsed:\s*(?P<elapsed>.+)$")
SCREEN_SHARDS_RE = re.compile(
    r"^Shards:\s*completed\s+(?P<completed>\d+)/(?P<total>\d+)\s+\|\s+running\s+(?P<running>\d+)\s+\|\s+pending\s+(?P<pending>\d+)\s+\|\s+failed\s+(?P<failed>\d+)\s*$"
)
PROGRESS_FRACTION_RE = re.compile(r"(?P<current>\d+)/(?P<total>\d+)")
SCREEN_RUNNING_EQUIVALENT_FRACTION = 0.5


def format_duration(total_seconds: int) -> str:
    total_seconds = max(0, int(total_seconds))
    days, remainder = divmod(total_seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, seconds = divmod(remainder, 60)
    parts = []
    if days:
        parts.append(f"{days}d")
    if hours or days:
        parts.append(f"{hours}h")
    if minutes or hours or days:
        parts.append(f"{minutes}m")
    parts.append(f"{seconds}s")
    return " ".join(parts)


def parse_elapsed_seconds(text: str) -> Optional[int]:
    stripped = text.strip()
    if not stripped:
        return None
    if stripped.endswith("s") and stripped[:-1].isdigit():
        return int(stripped[:-1])

    total = 0
    consumed = ""
    for match in re.finditer(r"(\d+)([dhms])", stripped):
        value = int(match.group(1))
        unit = match.group(2)
        consumed += match.group(0)
        if unit == "d":
            total += value * 86400
        elif unit == "h":
            total += value * 3600
        elif unit == "m":
            total += value * 60
        elif unit == "s":
            total += value
    if consumed:
        return total
    return None


def estimate_completion(
    elapsed_seconds: Optional[int],
    completed: int,
    total: int,
    active_equivalent: float = 0.0,
) -> Dict[str, object]:
    if elapsed_seconds is None or elapsed_seconds <= 0 or completed < 0 or total <= 0 or completed > total:
        return {}
    effective_completed = completed + max(0.0, active_equivalent)
    if effective_completed <= 0 or effective_completed > total:
        return {}
    remaining = max(0.0, total - effective_completed)
    estimated_total = elapsed_seconds * total / effective_completed
    remaining_seconds = max(0, int(round(estimated_total - elapsed_seconds)))
    return {
        "elapsed_seconds": int(elapsed_seconds),
        "remaining_seconds": remaining_seconds,
        "estimated_total_seconds": max(int(round(estimated_total)), int(elapsed_seconds)),
        "completion_fraction": effective_completed / total,
        "effective_completed": effective_completed,
    }


def latest_progress_marker(path: Optional[Path]) -> str:
    if path is None or not path.exists():
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


def load_manifest(campaign_dir: Path) -> Dict[str, object]:
    manifest_path = campaign_dir / "manifest.json"
    if not manifest_path.exists():
        return {}
    return json.loads(manifest_path.read_text())


def maybe_int(value: object, default: int = 0) -> int:
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def maybe_float(value: object, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def manifest_shard_slot(spec: object) -> Optional[tuple[str, str]]:
    if not isinstance(spec, dict):
        return None
    job = spec.get("job", {})
    if not isinstance(job, dict):
        return None
    stem = str(job.get("stem", ""))
    tag = str(spec.get("rerun_parent_tag") or spec.get("tag") or "")
    if not stem or not tag:
        return None
    return (stem, tag)


def load_target_shards_from_manifest(manifest: Dict[str, object]) -> List[Dict[str, object]]:
    entries = manifest.get("target_shards", [])
    if not isinstance(entries, list):
        return []
    shards: List[Dict[str, object]] = []
    seen: set[tuple[str, str, int]] = set()
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        job = entry.get("job", {})
        if not isinstance(job, dict) or "stem" not in job:
            continue
        slot = (
            str(job.get("stem", "")),
            str(entry.get("tag", "")),
            maybe_int(entry.get("seed")),
        )
        if not slot[0] or not slot[1] or slot in seen:
            continue
        seen.add(slot)
        shards.append(entry)
    return shards


def progress_fraction_from_marker(marker: str) -> Optional[float]:
    match = PROGRESS_FRACTION_RE.search(marker)
    if not match:
        return None
    current = maybe_int(match.group("current"), default=-1)
    total = maybe_int(match.group("total"), default=-1)
    if current < 0 or total <= 0:
        return None
    return max(0.0, min(1.0, current / total))


def extract_active_progress_samples(lines: Sequence[str]) -> List[float]:
    samples: List[float] = []
    in_active_section = False
    for raw_line in lines:
        line = raw_line.rstrip()
        stripped = line.strip()
        if stripped == "Active shards:":
            in_active_section = True
            continue
        if not in_active_section:
            continue
        if not stripped:
            break
        if stripped.startswith("+"):
            continue
        if not stripped.startswith("|"):
            if samples:
                break
            continue
        columns = [column.strip() for column in stripped.strip("|").split("|")]
        if len(columns) < 6 or columns[0] == "Run":
            continue
        fraction = progress_fraction_from_marker(columns[2])
        if fraction is not None:
            samples.append(fraction)
    return samples


def duration_from_manifest_entry(entry: object) -> Optional[float]:
    if not isinstance(entry, dict):
        return None
    duration = maybe_float(entry.get("duration_s"), default=-1.0)
    if duration > 0:
        return duration
    started_at = entry.get("started_at")
    ended_at = entry.get("ended_at")
    try:
        if started_at is None or ended_at is None:
            return None
        computed = float(ended_at) - float(started_at)
    except (TypeError, ValueError):
        return None
    return computed if computed > 0 else None


def successful_duration_by_slot(manifest: Dict[str, object]) -> Dict[tuple[str, str], float]:
    entries = manifest.get("finished", [])
    if not isinstance(entries, list):
        return {}
    durations: Dict[tuple[str, str], float] = {}
    ended_at_by_slot: Dict[tuple[str, str], float] = {}
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        slot = manifest_shard_slot(entry.get("spec"))
        duration = duration_from_manifest_entry(entry)
        if slot is None or duration is None:
            continue
        ended_at = maybe_float(entry.get("ended_at"), default=0.0)
        if slot not in durations or ended_at >= ended_at_by_slot.get(slot, 0.0):
            durations[slot] = duration
            ended_at_by_slot[slot] = ended_at
    return durations


def estimate_file_completion(
    manifest: Dict[str, object],
    successful: Sequence[Dict[str, object]],
    unresolved: Sequence[Dict[str, object]],
    total_shards: int,
    jobs_limit: int,
    now_ts: float,
) -> Dict[str, object]:
    if total_shards <= 0:
        return {}
    duration_map = successful_duration_by_slot(manifest)
    if not duration_map:
        return {}
    durations = list(duration_map.values())
    if not durations:
        return {}
    avg_duration = sum(durations) / len(durations)
    if avg_duration <= 0:
        return {}

    active_equivalent = 0.0
    for item in unresolved:
        fraction = progress_fraction_from_marker(str(item.get("progress", "")))
        if fraction is None:
            mtime = maybe_float(item.get("mtime"), default=0.0)
            if mtime > 0:
                fraction = max(0.0, min(0.99, (now_ts - mtime) / avg_duration))
        if fraction is not None:
            active_equivalent += max(0.0, min(0.99, fraction))

    jobs_limit = max(1, jobs_limit)
    effective_completed = min(total_shards, len(successful) + active_equivalent)
    remaining_equivalent = max(0.0, total_shards - effective_completed)
    estimated_total_seconds = max(0, int(round(avg_duration * total_shards / jobs_limit)))
    remaining_seconds = max(0, int(round(avg_duration * remaining_equivalent / jobs_limit)))
    return {
        "remaining_seconds": remaining_seconds,
        "estimated_total_seconds": estimated_total_seconds,
        "effective_completed": effective_completed,
        "average_shard_duration_seconds": avg_duration,
        "source": "manifest+progress",
    }


def extract_status_summary(manifest: Dict[str, object]) -> str:
    extract_status = manifest.get("extract_status", {})
    if not isinstance(extract_status, dict) or not extract_status:
        return "n/a"
    ordered_keys = ("raw_powheg", "yoda_merge", "yoda_nlo", "results", "diagnostic", "spin_diagnostic")
    parts = []
    for key in ordered_keys:
        if key in extract_status:
            parts.append(f"{key}={extract_status[key]}")
    for key in sorted(extract_status):
        if key not in ordered_keys:
            parts.append(f"{key}={extract_status[key]}")
    return ", ".join(parts) if parts else "n/a"


def list_screen_sessions() -> List[Dict[str, str]]:
    proc = subprocess.run(["screen", "-ls"], text=True, capture_output=True)
    if proc.returncode not in (0, 1):
        raise RuntimeError(proc.stderr.strip() or proc.stdout.strip() or "screen -ls failed")
    sessions: List[Dict[str, str]] = []
    for line in proc.stdout.splitlines():
        match = SCREEN_SESSION_RE.match(line)
        if not match:
            continue
        full_name = match.group("name")
        short_name = full_name.split(".", 1)[1] if "." in full_name else full_name
        sessions.append(
            {
                "full_name": full_name,
                "short_name": short_name,
                "state": match.group("state"),
            }
        )
    return sessions


def resolve_screen_session(name: str, sessions: Sequence[Dict[str, str]]) -> Optional[Dict[str, str]]:
    exact_matches = [
        session
        for session in sessions
        if name in (session["full_name"], session["short_name"])
    ]
    if len(exact_matches) == 1:
        return exact_matches[0]
    if len(exact_matches) > 1:
        raise ValueError(f"Screen session name {name!r} is ambiguous.")

    partial_matches = [
        session
        for session in sessions
        if (
            name in session["full_name"]
            or session["full_name"] in name
            or name in session["short_name"]
            or session["short_name"] in name
        )
    ]
    if len(partial_matches) == 1:
        return partial_matches[0]
    if len(partial_matches) > 1:
        raise ValueError(f"Screen session fragment {name!r} matches multiple sessions.")
    return None


def capture_screen_lines(session_name: str, window: Optional[str] = None) -> List[str]:
    with tempfile.NamedTemporaryFile(prefix="screen-hardcopy-", suffix=".txt", delete=False) as handle:
        temp_path = Path(handle.name)
    try:
        command = ["screen", "-S", session_name]
        if window:
            command.extend(["-p", window])
        command.extend(["-X", "hardcopy", "-h", str(temp_path)])
        proc = subprocess.run(
            command,
            text=True,
            capture_output=True,
        )
        if proc.returncode != 0:
            message = proc.stderr.strip() or proc.stdout.strip() or "screen hardcopy failed"
            raise RuntimeError(message)
        text = temp_path.read_text(errors="replace")
    finally:
        try:
            temp_path.unlink(missing_ok=True)
        except OSError:
            pass
    return [line.rstrip() for line in text.splitlines() if line.strip()]


def capture_screen_tail(session_name: str, tail_lines: int, window: Optional[str] = None) -> List[str]:
    lines = capture_screen_lines(session_name, window=window)
    if tail_lines <= 0:
        return []
    return lines[-tail_lines:]


def extract_screen_summary(lines: Sequence[str]) -> Dict[str, object]:
    summary: Dict[str, object] = {}
    candidate: Dict[str, object] = {}
    for line in lines:
        stripped = line.strip()
        tag_match = SCREEN_TAG_RE.match(stripped)
        if tag_match:
            if {"shards_line", "elapsed", "tag"}.issubset(candidate):
                summary = candidate
                break
            candidate = {"tag": tag_match.group("tag")}
            continue

        elapsed_match = SCREEN_ELAPSED_RE.match(stripped)
        if elapsed_match:
            if candidate:
                candidate["elapsed"] = elapsed_match.group("elapsed")
            continue

        shards_match = SCREEN_SHARDS_RE.match(stripped)
        if shards_match:
            if candidate:
                candidate["shards_line"] = stripped
                candidate["completed"] = int(shards_match.group("completed"))
                candidate["total"] = int(shards_match.group("total"))
                candidate["running"] = int(shards_match.group("running"))
                candidate["pending"] = int(shards_match.group("pending"))
                candidate["failed"] = int(shards_match.group("failed"))
                if {"shards_line", "elapsed", "tag"}.issubset(candidate):
                    summary = candidate
                    break
            continue

    if not summary and {"shards_line", "elapsed", "tag"}.issubset(candidate):
        summary = candidate
    active_progress_samples = extract_active_progress_samples(lines)
    active_equivalent = 0.0
    running = maybe_int(summary.get("running"), default=0)
    if active_progress_samples and running > 0:
        average_fraction = sum(active_progress_samples) / len(active_progress_samples)
        active_equivalent = average_fraction * running
        summary["active_progress_sample_count"] = len(active_progress_samples)
        summary["active_progress_average"] = average_fraction
    elif running > 0 and maybe_int(summary.get("completed"), default=0) > 0:
        active_equivalent = SCREEN_RUNNING_EQUIVALENT_FRACTION * running
    eta = estimate_completion(
        parse_elapsed_seconds(str(summary.get("elapsed", ""))),
        maybe_int(summary.get("completed"), default=0),
        maybe_int(summary.get("total"), default=0),
        active_equivalent,
    )
    if eta:
        eta["source"] = "screen-summary"
        summary["eta"] = eta
    return summary


def collect_screen_snapshot(name: str, tail_lines: int, window: Optional[str] = None) -> Dict[str, object]:
    # First try the session name exactly as provided. GNU screen accepts a
    # unique identifier or fragment here, which is more robust than parsing
    # `screen -ls` output up front.
    try:
        lines = capture_screen_lines(name, window=window)
        return {
            "requested": name,
            "status": "available",
            "session": name,
            "window": window or "",
            "summary": extract_screen_summary(lines),
            "tail": lines[-tail_lines:] if tail_lines > 0 else [],
        }
    except Exception as direct_exc:
        direct_error = str(direct_exc)

    sessions = list_screen_sessions()
    if not sessions:
        return {"requested": name, "status": "missing", "session": "", "window": window or "", "tail": [], "tail_error": direct_error}
    session = resolve_screen_session(name, sessions)
    if session is None:
        return {"requested": name, "status": "missing", "session": "", "window": window or "", "tail": [], "tail_error": direct_error}

    snapshot: Dict[str, object] = {
        "requested": name,
        "status": re.sub(r"[^a-z0-9]+", "-", str(session["state"]).strip().lower()).strip("-") or "unknown",
        "session": str(session["full_name"]),
        "window": window or "",
        "tail": [],
        "summary": {},
    }
    try:
        lines = capture_screen_lines(str(session["full_name"]), window=window)
        snapshot["summary"] = extract_screen_summary(lines)
        snapshot["tail"] = lines[-tail_lines:] if tail_lines > 0 else []
    except Exception as exc:
        snapshot["tail_error"] = str(exc)
    return snapshot


def logical_group_counts(items: Sequence[Dict[str, object]], total_by_run: Optional[Dict[str, int]] = None) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    for item in items:
        stem = str(item["stem"])
        counts[stem] = counts.get(stem, 0) + 1
    if total_by_run is not None:
        for stem in total_by_run:
            counts.setdefault(stem, 0)
    return counts


def collect_snapshot(args: argparse.Namespace) -> Dict[str, object]:
    now_ts = time.time()
    base_dir = args.base_dir.resolve()
    campaign_dir = base_dir / "campaigns" / args.tag
    manifest = load_manifest(campaign_dir)
    analysis_variant, include_lo, diagnostics_enabled, analysis_variants_by_setup = resolve_analysis_configuration(args)
    requested_setups = normalize_campaign_setups(args.setup)
    jobs_limit = max(1, maybe_int(manifest.get("jobs_limit"), args.jobs))
    target_shards = load_target_shards_from_manifest(manifest)
    if target_shards:
        jobs = []
        shards = target_shards
    else:
        manifest_event_counts = manifest.get("event_counts", {}) if isinstance(manifest.get("event_counts"), dict) else {}
        lo_events = maybe_int(manifest_event_counts.get("LO"), args.lo_events)
        posnlo_events = maybe_int(manifest_event_counts.get("POSNLO"), args.posnlo_events)
        negnlo_events = maybe_int(manifest_event_counts.get("NEGNLO"), args.negnlo_events)
        requested_shards = maybe_int(manifest.get("shards_per_logical_run"), args.shards)
        seed_base = maybe_int(manifest.get("seed_base"), args.seed_base)
        jobs = build_jobs(
            lo_events,
            posnlo_events,
            negnlo_events,
            analysis_variant=analysis_variant,
            analysis_variants_by_setup=analysis_variants_by_setup,
            setups=requested_setups,
            scale_variations=bool(args.scale_variations),
            campaign_tag=args.tag,
            include_lo=include_lo,
        )
        shards = build_shards(jobs, args.tag, requested_shards, jobs_limit, seed_base)
    expected_prefixes = [f"{spec['job']['stem']}-S{spec['seed']}-{spec['tag']}" for spec in shards]
    artifact_index = build_artifact_index(base_dir, expected_prefixes)
    launcher_dir = campaign_dir / "launcher-logs"
    launcher_names = build_launcher_log_index(launcher_dir)

    total_by_run: Dict[str, int] = {}
    for spec in shards:
        stem = str(spec["job"]["stem"])
        total_by_run[stem] = total_by_run.get(stem, 0) + 1

    successful: List[Dict[str, object]] = []
    unresolved: List[Dict[str, object]] = []
    missing: List[Dict[str, object]] = []
    launcher_mtimes: List[float] = []

    for spec in shards:
        stem = str(spec["job"]["stem"])
        prefix = f"{stem}-S{spec['seed']}-{spec['tag']}"
        launcher_name = f"{stem}-{spec['tag']}.launcher.log"
        launcher_path = launcher_dir / launcher_name if launcher_name in launcher_names else None
        artifacts = artifact_index.get(prefix, {"out_files": [], "log_files": [], "yoda_files": []})
        out_files = artifacts["out_files"]
        log_files = artifacts["log_files"]
        yoda_files = artifacts["yoda_files"]
        item = {
            "stem": stem,
            "tag": str(spec["tag"]),
            "seed": maybe_int(spec["seed"]),
            "events": maybe_int(spec["events"]),
            "launcher_log": str(launcher_path) if launcher_path is not None else "",
            "mtime": launcher_path.stat().st_mtime if launcher_path is not None and launcher_path.exists() else 0.0,
            "out_files": out_files,
            "log_files": log_files,
            "yoda_files": yoda_files,
            "progress": latest_progress_marker(launcher_path) if launcher_path is not None else "-",
        }
        if item["mtime"]:
            launcher_mtimes.append(float(item["mtime"]))
        if out_files or yoda_files:
            successful.append(item)
        elif launcher_path is not None or log_files:
            unresolved.append(item)
        else:
            missing.append(item)

    unresolved.sort(key=lambda item: (-float(item["mtime"]), str(item["stem"]), str(item["tag"])))
    unresolved_preview = []
    for item in unresolved[: max(0, args.max_listed)]:
        preview_item = dict(item)
        unresolved_preview.append(preview_item)

    successful_by_run = logical_group_counts(successful, total_by_run)
    unresolved_by_run = logical_group_counts(unresolved, total_by_run)
    missing_by_run = logical_group_counts(missing, total_by_run)
    logical_complete = sum(1 for stem, total in total_by_run.items() if successful_by_run.get(stem, 0) == total)
    logical_started = sum(1 for stem in total_by_run if successful_by_run.get(stem, 0) or unresolved_by_run.get(stem, 0))

    failed_entries = manifest.get("failed", [])
    if not isinstance(failed_entries, list):
        failed_entries = []

    total_shards = len(shards)
    complete_shards = len(successful)
    unresolved_shards = len(unresolved)
    missing_shards = len(missing)
    completion_pct = 0.0 if total_shards <= 0 else 100.0 * complete_shards / total_shards

    if total_shards > 0 and complete_shards == total_shards and unresolved_shards == 0:
        state = "complete"
    elif unresolved_shards > 0:
        state = "running"
    elif complete_shards > 0:
        state = "partial"
    else:
        state = "not-started"

    if failed_entries:
        state = "failed"

    screen_snapshot = (
        collect_screen_snapshot(args.screen_session, args.screen_tail_lines, window=args.screen_window)
        if args.screen_session
        else None
    )
    file_eta = estimate_file_completion(
        manifest,
        successful,
        unresolved,
        total_shards,
        jobs_limit,
        now_ts,
    )
    screen_eta = {}
    screen_total = 0
    if isinstance(screen_snapshot, dict):
        summary = screen_snapshot.get("summary")
        if isinstance(summary, dict):
            screen_total = maybe_int(summary.get("total"), default=0)
            screen_eta = summary.get("eta", {}) if isinstance(summary.get("eta"), dict) else {}
    screen_matches_target = screen_total in (0, total_shards)
    eta = file_eta or (screen_eta if screen_matches_target else {})

    return {
        "state": state,
        "host": socket.gethostname(),
        "timestamp": datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %Z"),
        "tag": args.tag,
        "base_dir": str(base_dir),
        "campaign_dir": str(campaign_dir),
        "setups": requested_setups,
        "analysis_variant": analysis_variant,
        "diagnostics_enabled": diagnostics_enabled,
        "manifest_present": bool(manifest),
        "extract_status_summary": extract_status_summary(manifest),
        "logical_runs_total": len(total_by_run),
        "logical_runs_started": logical_started,
        "logical_runs_complete": logical_complete,
        "shards_total": total_shards,
        "shards_complete": complete_shards,
        "shards_unresolved": unresolved_shards,
        "shards_missing": missing_shards,
        "completion_pct": completion_pct,
        "eta": eta,
        "file_eta": file_eta,
        "screen_eta": screen_eta,
        "screen_matches_target": screen_matches_target,
        "screen_total": screen_total,
        "unresolved_preview": unresolved_preview,
        "failed_count": len(failed_entries),
        "screen": screen_snapshot,
    }


def render_message(snapshot: Dict[str, object]) -> str:
    lines = [
        f"DISPOL status: {snapshot['state']}",
        f"Host: {snapshot['host']}",
        f"Time: {snapshot['timestamp']}",
        f"Tag: {snapshot['tag']}",
        f"Setups: {', '.join(snapshot['setups']) if snapshot['setups'] else 'default'}",
        f"Base dir: {snapshot['base_dir']}",
        f"Manifest: {'present' if snapshot['manifest_present'] else 'missing'}",
        (
            f"Logical runs: {snapshot['logical_runs_complete']}/{snapshot['logical_runs_total']} complete, "
            f"{snapshot['logical_runs_started']} started"
        ),
        (
            f"Shards: {snapshot['shards_complete']}/{snapshot['shards_total']} complete "
            f"({snapshot['completion_pct']:.1f}%), "
            f"{snapshot['shards_unresolved']} launched/unresolved, "
            f"{snapshot['shards_missing']} not yet launched"
        ),
    ]
    if snapshot["manifest_present"]:
        lines.append(f"Postprocess: {snapshot['extract_status_summary']}")
    if snapshot["failed_count"]:
        lines.append(f"Manifest failed shards: {snapshot['failed_count']}")
    eta = snapshot.get("eta")
    if isinstance(eta, dict) and eta:
        lines.append(
            "ETA: "
            f"{format_duration(maybe_int(eta.get('remaining_seconds'), 0))} remaining "
            f"(total ~{format_duration(maybe_int(eta.get('estimated_total_seconds'), 0))})"
        )
    elif snapshot.get("screen_matches_target") is False:
        lines.append(
            f"ETA: unavailable (screen session reports {snapshot.get('screen_total', 0)} total shards, "
            f"but the manifest target is {snapshot['shards_total']})"
        )
    screen = snapshot.get("screen")
    if isinstance(screen, dict):
        lines.append(
            f"Screen: {screen.get('status', 'unknown')}"
            + (
                f" ({screen.get('session')}{':' + str(screen.get('window')) if screen.get('window') else ''})"
                if screen.get("session")
                else ""
            )
        )
        summary = screen.get("summary")
        if isinstance(summary, dict) and summary:
            if summary.get("tag"):
                lines.append(f"Screen tag: {summary['tag']}")
            if summary.get("elapsed"):
                lines.append(f"Screen elapsed: {summary['elapsed']}")
            if summary.get("shards_line"):
                lines.append(f"Screen progress: {summary['shards_line']}")
            screen_eta = snapshot.get("screen_eta")
            if snapshot.get("screen_matches_target") is False:
                lines.append(
                    f"Screen target mismatch: screen shows {snapshot.get('screen_total', 0)} total shards, "
                    f"manifest target is {snapshot['shards_total']}"
                )
            if isinstance(screen_eta, dict) and screen_eta and not (isinstance(eta, dict) and eta):
                lines.append(
                    "Screen ETA (session-only): "
                    f"{format_duration(maybe_int(screen_eta.get('remaining_seconds'), 0))} remaining "
                    f"(total ~{format_duration(maybe_int(screen_eta.get('estimated_total_seconds'), 0))})"
                )
        tail_error = screen.get("tail_error")
        tail = screen.get("tail", [])
        if tail_error:
            lines.append(f"Screen tail: unavailable ({tail_error})")
        elif tail and not (isinstance(summary, dict) and summary.get("shards_line")):
            lines.append("")
            lines.append("Screen tail:")
            for line in tail:
                lines.append(f"> {line}")

    unresolved_preview = snapshot["unresolved_preview"]
    if unresolved_preview:
        lines.append("")
        lines.append("Recent unresolved shards:")
        for item in unresolved_preview:
            lines.append(
                f"- {item['stem']} [{item['tag']}] seed={item['seed']} progress={item['progress']}"
            )
    return "\n".join(lines)


def default_state_file(tag: str) -> Path:
    return Path("/tmp") / f"dispol-telegram-status-{tag}.json"


def fingerprint_message(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def load_state(path: Path) -> Dict[str, object]:
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def save_state(path: Path, message_hash: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "last_message_hash": message_hash,
        "updated_at": datetime.now().astimezone().isoformat(),
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True))


def send_telegram_message(token: str, chat_id: str, text: str, timeout_seconds: float) -> None:
    payload = urllib.parse.urlencode(
        {
            "chat_id": chat_id,
            "text": text,
            "disable_web_page_preview": "true",
        }
    ).encode("utf-8")
    request = urllib.request.Request(
        f"https://api.telegram.org/bot{token}/sendMessage",
        data=payload,
        method="POST",
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        body = response.read().decode("utf-8", errors="replace")
    parsed = json.loads(body)
    if not parsed.get("ok"):
        raise RuntimeError(f"Telegram API returned ok=false: {body}")


def run_once(args: argparse.Namespace) -> bool:
    snapshot = collect_snapshot(args)
    message = render_message(snapshot)
    message_hash = fingerprint_message(message)
    state_path = args.state_file
    if args.only_on_change:
        state = load_state(state_path)
        if state.get("last_message_hash") == message_hash:
            print(f"[skip] no status change for tag {args.tag}")
            return False

    if args.dry_run:
        print(message)
    else:
        if not args.telegram_bot_token:
            raise ValueError("Missing Telegram bot token. Use --telegram-bot-token or TELEGRAM_BOT_TOKEN.")
        if not args.telegram_chat_id:
            raise ValueError("Missing Telegram chat id. Use --telegram-chat-id or TELEGRAM_CHAT_ID.")
        send_telegram_message(args.telegram_bot_token, args.telegram_chat_id, message, args.timeout_seconds)
        print(f"[sent] Telegram update for {args.tag} at {snapshot['timestamp']}")

    if args.only_on_change:
        save_state(state_path, message_hash)
    return True


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=DEFAULT_BASE_DIR,
        help="Directory containing DISPOL cards, outputs, and campaigns/.",
    )
    parser.add_argument("-t", "--tag", required=True, help="Campaign tag to monitor.")
    parser.add_argument(
        "--setup",
        action="append",
        choices=SUPPORTED_SETUPS,
        default=[],
        help="Setup(s) used for the campaign. May be given multiple times.",
    )
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1, help="Maximum concurrent jobs used for the campaign.")
    parser.add_argument("--shards", type=int, default=0, help="Shards per logical run used for the campaign.")
    parser.add_argument("--seed-base", type=int, default=100000, help="Base RNG seed used for the campaign.")
    parser.add_argument("--lo-events", type=int, default=1_000_000_000, help="Events per LO run.")
    parser.add_argument("--posnlo-events", type=int, default=1_000_000_000, help="Events per POSNLO run.")
    parser.add_argument("--negnlo-events", type=int, default=10_000_000, help="Events per NEGNLO run.")
    parser.add_argument("--rivet", action="store_true", help="Campaign uses --rivet.")
    parser.add_argument("--rivetfo", action="store_true", help="Campaign uses --rivetfo.")
    parser.add_argument("--raw-powheg", action="store_true", help="Campaign uses --raw-powheg.")
    parser.add_argument("--include-lo", action="store_true", default=None, help="Campaign uses --include-lo.")
    parser.add_argument("--diagnostics", action="store_true", default=None, help="Campaign uses --diagnostics.")
    parser.add_argument("--scale-variations", action="store_true", default=None, help="Campaign uses --scale-variations.")
    parser.add_argument(
        "--telegram-bot-token",
        default=os.environ.get("TELEGRAM_BOT_TOKEN", ""),
        help="Telegram bot token. Defaults to TELEGRAM_BOT_TOKEN.",
    )
    parser.add_argument(
        "--telegram-chat-id",
        default=os.environ.get("TELEGRAM_CHAT_ID", ""),
        help="Telegram chat id. Defaults to TELEGRAM_CHAT_ID.",
    )
    parser.add_argument("--repeat-seconds", type=int, default=0, help="If > 0, keep sending updates at this interval.")
    parser.add_argument("--max-listed", type=int, default=5, help="Maximum unresolved shards to list in the message.")
    parser.add_argument("--screen-session", help="Optional GNU screen session name or unique fragment to inspect.")
    parser.add_argument("--screen-window", help="Optional GNU screen window index or title to inspect within the selected session.")
    parser.add_argument("--screen-tail-lines", type=int, default=6, help="Number of lines to capture from the screen scrollback tail.")
    parser.add_argument("--timeout-seconds", type=float, default=20.0, help="HTTP timeout for Telegram API calls.")
    parser.add_argument("--only-on-change", action="store_true", help="Skip sending when the rendered status message is unchanged.")
    parser.add_argument(
        "--state-file",
        type=Path,
        default=None,
        help="Optional external state file used by --only-on-change. Defaults to /tmp.",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print the message instead of sending it.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if args.state_file is None:
        args.state_file = default_state_file(args.tag)

    if args.repeat_seconds <= 0:
        run_once(args)
        return 0

    while True:
        try:
            run_once(args)
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            print(f"[error] {exc}", file=sys.stderr)
        time.sleep(args.repeat_seconds)


if __name__ == "__main__":
    raise SystemExit(main())
