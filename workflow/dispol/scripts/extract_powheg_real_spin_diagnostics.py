#!/usr/bin/env python3.10
"""
Extract dedicated POWHEG real-emission spin-vertex diagnostics from Herwig logs.

The script parses the log lines emitted by MENeutralCurrentDIS when the runtime
switch DiagnosePOWHEGRealSpinVertex is enabled:

  POWHEG_SPIN_EVENT
  POWHEG_SPIN_LEG

It summarizes, per selected log file, whether the realised 2->3 POWHEG state
has spinInfo(), productionVertex(), and sensible rho matrices on the incoming
and outgoing real-emission legs.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from prettytable import PrettyTable
except ImportError:
    PrettyTable = None


RUN_RE = re.compile(
    r"^DIS-POL-(?P<kind>POWHEG|LO)_(?P<hel>PP|PM|MP|MM|00)"
    r"(?:-(?P<piece>POSNLO|NEGNLO))?-(?P<ew>ALL|GAMMA|Z|SPINVAL|SPINCOMP)"
    r"(?:-(?P<analysis>RIVETPS-NOSPIN-UNPOL|RIVETPS-NOSPIN|RIVETPS-SPIN|RIVETFOFIXED|RIVETFO|RIVET))?"
    r"(?:-(?P<variant>.+))?$"
)
SEEDED_VARIANT_RE = re.compile(r"^S\d+-(?P<tag>.+)$")
SETUP_ORDER = {"GAMMA": 0, "Z": 1, "ALL": 2, "SPINVAL": 3, "SPINCOMP": 4}
HELICITY_ORDER = ("PP", "PM", "MP", "MM", "00")


@dataclass(frozen=True)
class EventSummary:
    event: int
    process: str
    enabled: bool
    parsed: bool
    all_spin: bool
    all_vertex: bool
    all_hard_vertex: bool
    all_sensible: bool


@dataclass(frozen=True)
class LegSummary:
    event: int
    process: str
    role: str
    incoming: bool
    pid: int
    has_spin: bool
    has_vertex: bool
    vertex_hard: bool
    production_location: int
    source: str
    trace: float
    trace_im: float
    antiherm: float
    min_diag: float
    max_diag: float
    max_diag_im: float
    d0: float
    d1: float
    d2: float
    finite: bool
    normalized: bool
    hermitian: bool
    sensible: bool


@dataclass(frozen=True)
class SelectedLog:
    path: Path
    paths: List[Path]
    run_name: str
    setup: str
    piece: str
    helicity: str
    analysis: str
    variant: str
    event_rows: List[EventSummary]
    leg_rows: List[LegSummary]


def normalize_variant_tag(variant: str) -> str:
    match = SEEDED_VARIANT_RE.match(variant)
    return match.group("tag") if match else variant


def variant_matches_tag(variant: str, preferred_tag: str) -> bool:
    normalized = normalize_variant_tag(variant)
    return bool(preferred_tag) and (
        normalized == preferred_tag or normalized.startswith(f"{preferred_tag}-s")
    )


def parse_run_name(run_name: str) -> Optional[Tuple[str, str, str, str, str]]:
    match = RUN_RE.match(run_name)
    if not match:
        return None
    kind = match.group("kind")
    piece = match.group("piece") or ("LO" if kind == "LO" else "")
    return (
        match.group("ew"),
        piece,
        match.group("hel"),
        match.group("analysis") or "",
        match.group("variant") or "",
    )


def parse_key_value_fields(payload: str) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for token in payload.strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        out[key] = value
    return out


def bool_field(fields: Dict[str, str], key: str) -> bool:
    return fields.get(key, "0") == "1"


def float_field(fields: Dict[str, str], key: str) -> float:
    try:
        return float(fields.get(key, "nan"))
    except ValueError:
        return float("nan")


def int_field(fields: Dict[str, str], key: str) -> int:
    try:
        return int(fields.get(key, "-1"))
    except ValueError:
        return -1


def parse_log(path: Path) -> Optional[SelectedLog]:
    event_rows: List[EventSummary] = []
    leg_rows: List[LegSummary] = []

    for line in path.read_text().splitlines():
        if line.startswith("POWHEG_SPIN_EVENT"):
            fields = parse_key_value_fields(line[len("POWHEG_SPIN_EVENT"):])
            event_rows.append(
                EventSummary(
                    event=int_field(fields, "event"),
                    process=fields.get("proc", ""),
                    enabled=bool_field(fields, "enabled"),
                    parsed=bool_field(fields, "parsed"),
                    all_spin=bool_field(fields, "allSpin"),
                    all_vertex=bool_field(fields, "allVertex"),
                    all_hard_vertex=bool_field(fields, "allHardVertex"),
                    all_sensible=bool_field(fields, "allSensible"),
                )
            )
        elif line.startswith("POWHEG_SPIN_LEG"):
            fields = parse_key_value_fields(line[len("POWHEG_SPIN_LEG"):])
            leg_rows.append(
                LegSummary(
                    event=int_field(fields, "event"),
                    process=fields.get("proc", ""),
                    role=fields.get("role", ""),
                    incoming=bool_field(fields, "incoming"),
                    pid=int_field(fields, "id"),
                    has_spin=bool_field(fields, "hasSpin"),
                    has_vertex=bool_field(fields, "hasVertex"),
                    vertex_hard=bool_field(fields, "vertexHard"),
                    production_location=int_field(fields, "loc"),
                    source=fields.get("source", ""),
                    trace=float_field(fields, "trace"),
                    trace_im=float_field(fields, "traceIm"),
                    antiherm=float_field(fields, "antiHerm"),
                    min_diag=float_field(fields, "minDiag"),
                    max_diag=float_field(fields, "maxDiag"),
                    max_diag_im=float_field(fields, "maxDiagIm"),
                    d0=float_field(fields, "d0"),
                    d1=float_field(fields, "d1"),
                    d2=float_field(fields, "d2"),
                    finite=bool_field(fields, "finite"),
                    normalized=bool_field(fields, "normalized"),
                    hermitian=bool_field(fields, "hermitian"),
                    sensible=bool_field(fields, "sensible"),
                )
            )

    if not event_rows and not leg_rows:
        return None

    run_name = path.stem
    parsed = parse_run_name(run_name)
    if parsed is None:
        return None
    setup, piece, helicity, analysis, variant = parsed

    return SelectedLog(
        path=path,
        paths=[path],
        run_name=run_name,
        setup=setup,
        piece=piece,
        helicity=helicity,
        analysis=analysis,
        variant=variant,
        event_rows=event_rows,
        leg_rows=leg_rows,
    )


def choose_better(current: SelectedLog, candidate: SelectedLog, preferred_tag: str, strict_tag: bool) -> bool:
    current_tag = normalize_variant_tag(current.variant)
    candidate_tag = normalize_variant_tag(candidate.variant)

    if strict_tag:
        return False

    current_matches = variant_matches_tag(current.variant, preferred_tag)
    candidate_matches = variant_matches_tag(candidate.variant, preferred_tag)
    if candidate_matches != current_matches:
        return candidate_matches
    if (current_tag == "") != (candidate_tag == ""):
        return candidate_tag != ""
    return candidate.path.stat().st_mtime > current.path.stat().st_mtime


def select_logs(base_dir: Path, preferred_tag: str, strict_tag: bool) -> List[SelectedLog]:
    selected: Dict[Tuple[str, str, str, str], SelectedLog] = {}
    for path in sorted(base_dir.rglob("DIS-POL-*.log")):
        parsed = parse_log(path)
        if parsed is None:
            continue
        key = (parsed.setup, parsed.piece, parsed.helicity, parsed.analysis)
        matches = variant_matches_tag(parsed.variant, preferred_tag)
        if strict_tag and preferred_tag and not matches:
            continue
        if key not in selected:
            selected[key] = parsed
            continue
        if choose_better(selected[key], parsed, preferred_tag, strict_tag):
            selected[key] = parsed
    return sorted(
        selected.values(),
        key=lambda item: (
            SETUP_ORDER.get(item.setup, 99),
            item.piece,
            HELICITY_ORDER.index(item.helicity) if item.helicity in HELICITY_ORDER else 99,
            item.analysis,
            item.variant,
        ),
    )


def render_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    if PrettyTable is not None:
        table = PrettyTable()
        table.field_names = list(headers)
        table.align = "l"
        for row in rows:
            table.add_row(list(row))
        return table.get_string()

    widths = [len(str(header)) for header in headers]
    for row in rows:
        for idx, cell in enumerate(row):
            widths[idx] = max(widths[idx], len(str(cell)))

    def fmt(row: Sequence[object]) -> str:
        return " | ".join(str(cell).ljust(widths[idx]) for idx, cell in enumerate(row))

    border = "-+-".join("-" * width for width in widths)
    parts = [fmt(headers), border]
    parts.extend(fmt(row) for row in rows)
    return "\n".join(parts)


def format_float(value: float) -> str:
    if not math.isfinite(value):
        return "nan"
    return f"{value:.6g}"


def build_summary(logs: Sequence[SelectedLog]) -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    event_summary: List[Dict[str, object]] = []
    leg_summary: List[Dict[str, object]] = []

    for log in logs:
        parsed_events = [row for row in log.event_rows if row.parsed]
        event_summary.append(
            {
                "run": log.run_name,
                "setup": log.setup,
                "piece": log.piece,
                "helicity": log.helicity,
                "analysis": log.analysis,
                "variant": log.variant,
                "events": len(log.event_rows),
                "parsed": len(parsed_events),
                "all_spin": sum(1 for row in parsed_events if row.all_spin),
                "all_vertex": sum(1 for row in parsed_events if row.all_vertex),
                "all_hard_vertex": sum(1 for row in parsed_events if row.all_hard_vertex),
                "all_sensible": sum(1 for row in parsed_events if row.all_sensible),
                "path": str(log.path),
            }
        )

        by_role: Dict[str, List[LegSummary]] = {}
        for row in log.leg_rows:
            by_role.setdefault(row.role, []).append(row)
        for role, role_rows in sorted(by_role.items()):
            finite_trace = [abs(row.trace - 1.0) for row in role_rows if math.isfinite(row.trace)]
            finite_antiherm = [row.antiherm for row in role_rows if math.isfinite(row.antiherm)]
            finite_min_diag = [row.min_diag for row in role_rows if math.isfinite(row.min_diag)]
            finite_max_diag = [row.max_diag for row in role_rows if math.isfinite(row.max_diag)]
            leg_summary.append(
                {
                    "run": log.run_name,
                    "setup": log.setup,
                    "piece": log.piece,
                    "helicity": log.helicity,
                    "analysis": log.analysis,
                    "variant": log.variant,
                    "role": role,
                    "entries": len(role_rows),
                    "has_spin": sum(1 for row in role_rows if row.has_spin),
                    "has_vertex": sum(1 for row in role_rows if row.has_vertex),
                    "vertex_hard": sum(1 for row in role_rows if row.vertex_hard),
                    "sensible": sum(1 for row in role_rows if row.sensible),
                    "max_trace_dev": max(finite_trace) if finite_trace else float("nan"),
                    "max_antiherm": max(finite_antiherm) if finite_antiherm else float("nan"),
                    "min_diag": min(finite_min_diag) if finite_min_diag else float("nan"),
                    "max_diag": max(finite_max_diag) if finite_max_diag else float("nan"),
                    "path": str(log.path),
                }
            )

    return event_summary, leg_summary


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=Path.cwd())
    parser.add_argument("-t", "--tag", default="", help="Preferred tag, e.g. testing13")
    parser.add_argument("--strict-tag", action="store_true", help="Only use logs matching the requested tag")
    parser.add_argument("--json-out", type=Path)
    parser.add_argument("--csv-out", type=Path)
    args = parser.parse_args(argv)

    logs = select_logs(args.base_dir, args.tag, args.strict_tag)
    event_summary, leg_summary = build_summary(logs)

    selected_rows = [
        [item.setup, item.piece, item.helicity, item.analysis or "-", item.run_name, item.variant or "-", item.path.name]
        for item in logs
    ]
    print("Selected spin diagnostic logs:")
    print(render_table(["Setup", "Piece", "Hel", "Analysis", "Run", "Variant", "Log file"], selected_rows))
    print()

    event_rows = [
        [
            row["setup"],
            row["piece"],
            row["helicity"],
            row["analysis"] or "-",
            row["events"],
            f'{row["all_spin"]}/{row["parsed"]}',
            f'{row["all_vertex"]}/{row["parsed"]}',
            f'{row["all_hard_vertex"]}/{row["parsed"]}',
            f'{row["all_sensible"]}/{row["parsed"]}',
        ]
        for row in event_summary
    ]
    print("Per-run event summary:")
    print(render_table(
        ["Setup", "Piece", "Hel", "Analysis", "Events", "allSpin", "allVertex", "allHard", "allSensible"],
        event_rows,
    ))
    print()

    leg_rows = [
        [
            row["setup"],
            row["piece"],
            row["helicity"],
            row["analysis"] or "-",
            row["role"],
            row["entries"],
            f'{row["has_spin"]}/{row["entries"]}',
            f'{row["has_vertex"]}/{row["entries"]}',
            f'{row["vertex_hard"]}/{row["entries"]}',
            f'{row["sensible"]}/{row["entries"]}',
            format_float(row["max_trace_dev"]),
            format_float(row["max_antiherm"]),
            format_float(row["min_diag"]),
            format_float(row["max_diag"]),
        ]
        for row in leg_summary
    ]
    print("Per-leg summary:")
    print(render_table(
        ["Setup", "Piece", "Hel", "Analysis", "Role", "Entries", "hasSpin", "hasVertex", "hard", "sensible",
         "max|tr-1|", "maxAntiHerm", "minDiag", "maxDiag"],
        leg_rows,
    ))

    payload = {
        "selected_logs": [
            {
                "run": item.run_name,
                "setup": item.setup,
                "piece": item.piece,
                "helicity": item.helicity,
                "analysis": item.analysis,
                "variant": item.variant,
                "path": str(item.path),
            }
            for item in logs
        ],
        "event_summary": event_summary,
        "leg_summary": leg_summary,
    }

    if args.json_out is not None:
        args.json_out.write_text(json.dumps(payload, indent=2, sort_keys=True))

    if args.csv_out is not None:
        with args.csv_out.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "section", "run", "setup", "piece", "helicity", "analysis", "variant", "role",
                    "events", "parsed", "all_spin", "all_vertex", "all_hard_vertex", "all_sensible",
                    "entries", "has_spin", "has_vertex", "vertex_hard", "sensible",
                    "max_trace_dev", "max_antiherm", "min_diag", "max_diag", "path",
                ],
            )
            writer.writeheader()
            for row in event_summary:
                writer.writerow({"section": "event", "role": "", **row})
            for row in leg_summary:
                writer.writerow({"section": "leg", "events": "", "parsed": "",
                                 "all_spin": "", "all_vertex": "", "all_hard_vertex": "",
                                 "all_sensible": "", **row})

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
