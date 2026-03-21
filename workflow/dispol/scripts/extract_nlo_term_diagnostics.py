#!/usr/bin/env python3.10
"""
Extract component-level NLO diagnostics from DIS POWHEG log files.

The script looks for the diagnostic lines emitted by DISBase::NLOWeight():

  NLO_TERM_CUM
  NLO_TERM_REAL
  NLO_TERM_APOL
  NLO_TERM_BORN
  NLO_TERM_COEFF

and combines them with the matching POSNLO/NEGNLO run cross sections from the
corresponding .out files. The resulting component tables are intended as a
diagnostic aid for comparing the polarized/unpolarized NLO correction pieces
between GAMMA, Z, and ALL, and for isolating the ALL - GAMMA - Z interference.

The per-component pb values are reconstructed as

  sigma_component(run) ~= sigma_run(.out) * F_component(run)

where F_component is the cumulative fraction of that specific POSNLO or NEGNLO
run cross section reported in the diagnostic logs.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, TypeVar

try:
    from prettytable import PrettyTable
except ImportError:
    PrettyTable = None


RUN_RE = re.compile(
    r"^DIS-POL-(?P<kind>POWHEG|LO)_(?P<hel>PP|PM|MP|MM|00)"
    r"(?:-(?P<piece>POSNLO|NEGNLO))?-(?P<ew>ALL|GAMMA|Z)"
    r"(?:-(?P<analysis>RIVETFOFIXED|RIVETFO|RIVET))?"
    r"(?:-(?P<variant>.+))?$"
)
SEEDED_VARIANT_RE = re.compile(r"^S\d+-(?P<tag>.+)$")
TOTAL_RE = re.compile(
    r"Total\s+\(from generated events\):\s+(\S+)\s+(\S+)\s+([0-9.]+)\((\d+)\)e([+-]?\d+)"
)
UNIT_RE = re.compile(r"Cross-section\s*\(([^)]+)\)")

SETUP_ORDER = {"GAMMA": 0, "Z": 1, "ALL": 2}
PRIMARY_COMPONENTS = (
    "F_virt",
    "F_cq_even",
    "F_cq_odd",
    "F_cg_even",
    "F_cg_odd",
    "F_rq",
    "F_rg",
    "F_total",
)
REAL_SPLIT_COMPONENTS = (
    "F_rq_even",
    "F_rq_odd",
    "F_rg_even",
    "F_rg_odd",
    "F_rq_chk",
    "F_rg_chk",
)
HELICITY_ORDER = ("PP", "PM", "MP", "MM", "00")
T = TypeVar("T")


@dataclass(frozen=True)
class Measurement:
    value_pb: float
    error_pb: float


@dataclass(frozen=True)
class OutRun:
    path: Path
    paths: List[Path]
    setup: str
    piece: str
    helicity: str
    analysis: str
    variant: str
    measurement: Measurement


@dataclass(frozen=True)
class DiagnosticRun:
    path: Path
    paths: List[Path]
    run_name: str
    setup: str
    piece: str
    helicity: str
    analysis: str
    variant: str
    n_events: Optional[int]
    values: Dict[str, float]


def nb_to_pb(value: float) -> float:
    return value * 1.0e3


def value_with_paren_error(mantissa_text: str, err_digits_text: str, exponent_text: str) -> Tuple[float, float]:
    mantissa = float(mantissa_text)
    exponent = int(exponent_text)
    value = mantissa * (10.0 ** exponent)
    decimals = len(mantissa_text.split(".", 1)[1]) if "." in mantissa_text else 0
    error = int(err_digits_text) * (10.0 ** (exponent - decimals))
    return value, error


def parse_measurement(text: str) -> Optional[Measurement]:
    unit_match = UNIT_RE.search(text)
    unit = unit_match.group(1).strip().lower() if unit_match else "nb"
    matches = list(TOTAL_RE.finditer(text))
    if not matches:
        return None

    value_raw, error_raw = value_with_paren_error(matches[-1].group(3), matches[-1].group(4), matches[-1].group(5))
    if unit == "pb":
        return Measurement(value_raw, error_raw)
    return Measurement(nb_to_pb(value_raw), nb_to_pb(error_raw))


def parse_run_name(run_name: str) -> Optional[Tuple[str, str, str, str, str]]:
    match = RUN_RE.match(run_name)
    if not match:
        return None
    kind = match.group("kind")
    helicity = match.group("hel")
    piece = match.group("piece") or ("LO" if kind == "LO" else "")
    setup = match.group("ew")
    analysis = match.group("analysis") or ""
    variant = match.group("variant") or ""
    return setup, piece, helicity, analysis, variant


def is_shard_variant(variant: str, preferred_tag: str) -> bool:
    normalized = normalize_variant_tag(variant)
    return bool(preferred_tag) and normalized.startswith(f"{preferred_tag}-s")


def variant_matches_tag(variant: str, preferred_tag: str) -> bool:
    normalized = normalize_variant_tag(variant)
    return bool(preferred_tag) and (normalized == preferred_tag or is_shard_variant(variant, preferred_tag))


def normalize_variant_tag(variant: str) -> str:
    match = SEEDED_VARIANT_RE.match(variant)
    if match:
        return match.group("tag")
    return variant


def parse_key_value_fields(payload: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    for token in payload.strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = value
    return fields


def parse_diagnostic_log(path: Path) -> Optional[DiagnosticRun]:
    last_cum: Optional[Dict[str, str]] = None
    last_real: Optional[Dict[str, str]] = None
    last_apol: Optional[Dict[str, str]] = None
    last_born: Optional[Dict[str, str]] = None
    last_coeff: Optional[Dict[str, str]] = None
    for line in path.read_text().splitlines():
        if line.startswith("NLO_TERM_CUM"):
            last_cum = parse_key_value_fields(line[len("NLO_TERM_CUM"):])
        elif line.startswith("NLO_TERM_REAL"):
            last_real = parse_key_value_fields(line[len("NLO_TERM_REAL"):])
        elif line.startswith("NLO_TERM_APOL"):
            last_apol = parse_key_value_fields(line[len("NLO_TERM_APOL"):])
        elif line.startswith("NLO_TERM_BORN"):
            last_born = parse_key_value_fields(line[len("NLO_TERM_BORN"):])
        elif line.startswith("NLO_TERM_COEFF"):
            last_coeff = parse_key_value_fields(line[len("NLO_TERM_COEFF"):])

    if last_cum is None and last_real is None and last_apol is None and last_born is None and last_coeff is None:
        return None

    source = last_cum or last_real or last_apol or last_born or last_coeff
    assert source is not None
    run_name = source.get("run", path.stem)
    parsed = parse_run_name(run_name)
    if parsed is None:
        return None
    setup, piece, helicity, analysis, variant = parsed

    values: Dict[str, float] = {}
    n_events: Optional[int] = None
    for payload in (last_cum, last_real, last_apol, last_born, last_coeff):
        if payload is None:
            continue
        if "n" in payload:
            try:
                n_events = int(payload["n"])
            except ValueError:
                pass
        for key, raw in payload.items():
            if key in ("run", "n"):
                continue
            try:
                values[key] = float(raw)
            except ValueError:
                continue

    return DiagnosticRun(
        path=path,
        paths=[path],
        run_name=run_name,
        setup=setup,
        piece=piece,
        helicity=helicity,
        analysis=analysis,
        variant=variant,
        n_events=n_events,
        values=values,
    )


def parse_out(path: Path) -> Optional[OutRun]:
    if path.suffix != ".out":
        return None
    parsed = parse_run_name(path.stem)
    if parsed is None:
        return None
    setup, piece, helicity, analysis, variant = parsed
    measurement = parse_measurement(path.read_text())
    if measurement is None:
        return None
    return OutRun(
        path=path,
        paths=[path],
        setup=setup,
        piece=piece,
        helicity=helicity,
        analysis=analysis,
        variant=variant,
        measurement=measurement,
    )


def run_preference(variant: str, path: Path, preferred_tag: str) -> Tuple[int, int, float, str]:
    preferred_penalty = 0 if variant_matches_tag(variant, preferred_tag) else 1
    default_penalty = 0 if normalize_variant_tag(variant) == "" else 1
    if not preferred_tag:
        preferred_penalty = 0
    return (preferred_penalty, default_penalty, -path.stat().st_mtime, path.name)


def combine_measurements(measurements: Sequence[Measurement]) -> Measurement:
    if not measurements:
        raise ValueError("combine_measurements() requires at least one measurement")
    positive = [m for m in measurements if m.error_pb > 0.0]
    if not positive:
        value = sum(m.value_pb for m in measurements) / len(measurements)
        return Measurement(value, 0.0)
    weights = [1.0 / (m.error_pb ** 2) for m in positive]
    total_weight = sum(weights)
    value = sum(w * m.value_pb for w, m in zip(weights, positive)) / total_weight
    error = (1.0 / total_weight) ** 0.5
    return Measurement(value, error)


def select_candidates_by_tag(
    candidates: Sequence[T],
    preferred_tag: str,
    strict_tag: bool,
    get_variant: Callable[[T], str],
) -> List[T]:
    items = list(candidates)
    if preferred_tag:
        shard_matches = [item for item in items if is_shard_variant(get_variant(item), preferred_tag)]
        exact_matches = [item for item in items if normalize_variant_tag(get_variant(item)) == preferred_tag]
        if shard_matches:
            return sorted(shard_matches, key=lambda item: get_variant(item))
        if exact_matches:
            return [sorted(exact_matches, key=lambda item: get_variant(item))[0]]
        if strict_tag:
            return []
    return items


def aggregate_out_group(candidates: Sequence[OutRun], preferred_tag: str) -> Optional[OutRun]:
    chosen = select_candidates_by_tag(candidates, preferred_tag, False, lambda item: item.variant)
    if not chosen:
        return None
    if preferred_tag:
        shard_matches = [item for item in chosen if is_shard_variant(item.variant, preferred_tag)]
        exact_matches = [item for item in chosen if normalize_variant_tag(item.variant) == preferred_tag]
        if shard_matches:
            chosen = shard_matches
        elif exact_matches:
            chosen = [exact_matches[0]]
        else:
            chosen = sorted(chosen, key=lambda item: run_preference(item.variant, item.path, preferred_tag))[:1]
    else:
        chosen = sorted(chosen, key=lambda item: run_preference(item.variant, item.path, preferred_tag))[:1]

    measurement = combine_measurements([item.measurement for item in chosen])
    paths = sorted((item.path for item in chosen), key=lambda item: item.name)
    return OutRun(
        path=paths[0],
        paths=paths,
        setup=chosen[0].setup,
        piece=chosen[0].piece,
        helicity=chosen[0].helicity,
        variant=preferred_tag or chosen[0].variant,
        measurement=measurement,
    )


def aggregate_diag_group(
    candidates: Sequence[DiagnosticRun],
    out_group: Optional[Sequence[OutRun]],
    preferred_tag: str,
    strict_tag: bool,
) -> Optional[DiagnosticRun]:
    chosen = select_candidates_by_tag(candidates, preferred_tag, strict_tag, lambda item: item.variant)
    if not chosen:
        return None
    if not preferred_tag:
        chosen = sorted(chosen, key=lambda item: run_preference(item.variant, item.path, preferred_tag))[:1]

    out_by_variant = {item.variant: item for item in (out_group or [])}
    weighted_entries: List[Tuple[float, DiagnosticRun]] = []
    for item in chosen:
        out = out_by_variant.get(item.variant)
        if out is not None and out.measurement.error_pb > 0.0:
            weight = 1.0 / (out.measurement.error_pb ** 2)
        else:
            weight = 1.0
        weighted_entries.append((weight, item))

    total_weight = sum(weight for weight, _ in weighted_entries)
    values: Dict[str, float] = {}
    all_keys = sorted({key for _, item in weighted_entries for key in item.values})
    for key in all_keys:
        accum = 0.0
        norm = 0.0
        for weight, item in weighted_entries:
            if key not in item.values:
                continue
            accum += weight * item.values[key]
            norm += weight
        if norm > 0.0:
            values[key] = accum / norm

    n_events = sum(item.n_events or 0 for _, item in weighted_entries) or None
    paths = sorted((item.path for _, item in weighted_entries), key=lambda item: item.name)
    run_name = chosen[0].run_name if len(chosen) == 1 else f"{chosen[0].setup}:{chosen[0].piece}:{chosen[0].helicity} ({len(chosen)} shards)"
    return DiagnosticRun(
        path=paths[0],
        paths=paths,
        run_name=run_name,
        setup=chosen[0].setup,
        piece=chosen[0].piece,
        helicity=chosen[0].helicity,
        variant=preferred_tag or chosen[0].variant,
        n_events=n_events,
        values=values,
    )


def select_best_runs(
    entries: Iterable[DiagnosticRun],
    out_entries: Iterable[OutRun],
    preferred_tag: str,
    strict_tag: bool,
) -> Tuple[Dict[Tuple[str, str, str], DiagnosticRun], Dict[Tuple[str, str, str], List[Path]]]:
    grouped: Dict[Tuple[str, str, str], List[DiagnosticRun]] = {}
    for entry in entries:
        grouped.setdefault((entry.setup, entry.piece, entry.helicity), []).append(entry)
    out_grouped: Dict[Tuple[str, str, str], List[OutRun]] = {}
    for entry in out_entries:
        out_grouped.setdefault((entry.setup, entry.piece, entry.helicity), []).append(entry)

    chosen: Dict[Tuple[str, str, str], DiagnosticRun] = {}
    duplicates: Dict[Tuple[str, str, str], List[Path]] = {}
    for key, candidates in grouped.items():
        aggregated = aggregate_diag_group(candidates, out_grouped.get(key), preferred_tag, strict_tag)
        if aggregated is None:
            continue
        chosen[key] = aggregated
        if len(candidates) > 1:
            duplicates[key] = [cand.path for cand in sorted(candidates, key=lambda item: run_preference(item.variant, item.path, preferred_tag))]
    return chosen, duplicates


def select_best_lo(
    entries: Iterable[OutRun],
    preferred_tag: str,
    strict_tag: bool,
) -> Tuple[Dict[Tuple[str, str, str], OutRun], Dict[Tuple[str, str, str], List[Path]]]:
    grouped: Dict[Tuple[str, str, str], List[OutRun]] = {}
    for entry in entries:
        grouped.setdefault((entry.setup, entry.piece, entry.helicity), []).append(entry)

    chosen: Dict[Tuple[str, str, str], OutRun] = {}
    duplicates: Dict[Tuple[str, str, str], List[Path]] = {}
    for key, candidates in grouped.items():
        selected = select_candidates_by_tag(candidates, preferred_tag, strict_tag, lambda item: item.variant)
        if not selected:
            continue
        aggregated = aggregate_out_group(selected, preferred_tag)
        if aggregated is None:
            continue
        chosen[key] = aggregated
        if len(candidates) > 1:
            duplicates[key] = [cand.path for cand in sorted(candidates, key=lambda item: run_preference(item.variant, item.path, preferred_tag))]
    return chosen, duplicates


def render_table(headers: Sequence[str], rows: Sequence[Sequence[str]], aligns: Optional[Sequence[str]] = None) -> List[str]:
    if not rows:
        return []
    alignments = list(aligns) if aligns is not None else ["l"] * len(headers)
    rows_as_text = [[str(cell) for cell in row] for row in rows]

    if PrettyTable is not None:
        table = PrettyTable()
        table.field_names = list(headers)
        for header, align in zip(headers, alignments):
            table.align[header] = "r" if align == "r" else "l"
        for row in rows_as_text:
            table.add_row(row)
        return table.get_string().splitlines()

    widths = []
    for i, header in enumerate(headers):
        widths.append(max(len(header), max((len(row[i]) for row in rows_as_text), default=0)))

    def format_row(row: Sequence[str], is_header: bool = False) -> str:
        out = []
        for i, cell in enumerate(row):
            width = widths[i]
            if is_header:
                out.append(cell.ljust(width))
            elif alignments[i] == "r":
                out.append(cell.rjust(width))
            else:
                out.append(cell.ljust(width))
        return "| " + " | ".join(out) + " |"

    sep = "+-" + "-+-".join("-" * width for width in widths) + "-+"
    lines = [sep, format_row(headers, True), sep]
    lines.extend(format_row(row) for row in rows_as_text)
    lines.append(sep)
    return lines


def fmt_pb(value: Optional[float]) -> str:
    if value is None:
        return "missing"
    return f"{value:+.6f} pb"


def fmt_unitless(value: Optional[float]) -> str:
    if value is None:
        return "missing"
    return f"{value:+.6f}"


def fmt_measurement(measurement: Measurement) -> str:
    return f"{measurement.value_pb:.6f} +- {measurement.error_pb:.6f} pb"


def combine_observables(setup: str, per_hel: Dict[str, float]) -> Dict[str, Optional[float]]:
    out = {"unpol_pb": None, "pol_pb": None}
    if setup == "GAMMA":
        if "00" in per_hel:
            out["unpol_pb"] = per_hel["00"]
        if "PP" in per_hel and "PM" in per_hel:
            out["pol_pb"] = 0.5 * (per_hel["PP"] - per_hel["PM"])
        return out

    needed = ("PP", "PM", "MP", "MM")
    if all(hel in per_hel for hel in needed):
        out["unpol_pb"] = 0.25 * (per_hel["PP"] + per_hel["PM"] + per_hel["MP"] + per_hel["MM"])
        out["pol_pb"] = 0.25 * (per_hel["PP"] + per_hel["MM"] - per_hel["PM"] - per_hel["MP"])
    return out


def component_estimates_for_piece(
    setup: str,
    piece: str,
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
    selected_out: Dict[Tuple[str, str, str], OutRun],
    components: Sequence[str],
) -> Dict[str, Dict[str, Optional[float]]]:
    per_component: Dict[str, Dict[str, Optional[float]]] = {}
    for component in components:
        per_hel: Dict[str, float] = {}
        for helicity in HELICITY_ORDER:
            diag = selected_diag.get((setup, piece, helicity))
            out = selected_out.get((setup, piece, helicity))
            if diag is None or out is None:
                continue
            if component not in diag.values:
                continue
            per_hel[helicity] = out.measurement.value_pb * diag.values[component]
        per_component[component] = combine_observables(setup, per_hel)
    return per_component


def subtract_estimates(
    first: Dict[str, Dict[str, Optional[float]]],
    second: Dict[str, Dict[str, Optional[float]]],
    components: Sequence[str],
) -> Dict[str, Dict[str, Optional[float]]]:
    result: Dict[str, Dict[str, Optional[float]]] = {}
    for component in components:
        out = {"unpol_pb": None, "pol_pb": None}
        for observable in ("unpol_pb", "pol_pb"):
            a = first.get(component, {}).get(observable)
            b = second.get(component, {}).get(observable)
            if a is None or b is None:
                continue
            out[observable] = a - b
        result[component] = out
    return result


def compute_all_estimates(
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
    selected_out: Dict[Tuple[str, str, str], OutRun],
    components: Sequence[str],
) -> Dict[str, Dict[str, Dict[str, Dict[str, Optional[float]]]]]:
    estimates: Dict[str, Dict[str, Dict[str, Dict[str, Optional[float]]]]] = {}
    for setup in sorted({key[0] for key in selected_diag}, key=lambda s: (SETUP_ORDER.get(s, 99), s)):
        pos = component_estimates_for_piece(setup, "POSNLO", selected_diag, selected_out, components)
        neg = component_estimates_for_piece(setup, "NEGNLO", selected_diag, selected_out, components)
        nlo = subtract_estimates(pos, neg, components)
        estimates[setup] = {"POSNLO": pos, "NEGNLO": neg, "NLO": nlo}
    return estimates


def interference_estimates(
    estimates: Dict[str, Dict[str, Dict[str, Dict[str, Optional[float]]]]],
    components: Sequence[str],
) -> Dict[str, Dict[str, Optional[float]]]:
    if not all(setup in estimates for setup in ("ALL", "GAMMA", "Z")):
        return {}
    out: Dict[str, Dict[str, Optional[float]]] = {}
    for component in components:
        values = {}
        for observable in ("unpol_pb", "pol_pb"):
            a = estimates["ALL"]["NLO"].get(component, {}).get(observable)
            g = estimates["GAMMA"]["NLO"].get(component, {}).get(observable)
            z = estimates["Z"]["NLO"].get(component, {}).get(observable)
            if a is None or g is None or z is None:
                values[observable] = None
            else:
                values[observable] = a - g - z
        out[component] = values
    return out


def build_report(
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
    diag_duplicates: Dict[Tuple[str, str, str], List[Path]],
    selected_out: Dict[Tuple[str, str, str], OutRun],
    out_duplicates: Dict[Tuple[str, str, str], List[Path]],
    preferred_tag: str,
    strict_tag: bool,
) -> str:
    lines: List[str] = []
    if not selected_diag:
        return "No NLO_TERM diagnostics found.\n"

    if preferred_tag:
        mode = "strict" if strict_tag else "fallback-enabled"
        lines.append(f"Tag selection: {preferred_tag} ({mode})")
        lines.append("")

    if diag_duplicates:
        lines.append("Diagnostic duplicates resolved by preference:")
        dup_rows = []
        for key in sorted(diag_duplicates, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
            setup, piece, helicity = key
            dup_rows.append([setup, piece, helicity, ", ".join(path.name for path in diag_duplicates[key])])
        lines.extend(render_table(["Setup", "Piece", "Helicity", "Candidates"], dup_rows))
        lines.append("")

    if out_duplicates:
        lines.append(".out duplicates resolved by preference:")
        dup_rows = []
        for key in sorted(out_duplicates, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
            setup, piece, helicity = key
            dup_rows.append([setup, piece, helicity, ", ".join(path.name for path in out_duplicates[key])])
        lines.extend(render_table(["Setup", "Piece", "Helicity", "Candidates"], dup_rows))
        lines.append("")

    run_rows = []
    for key in sorted(selected_diag, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
        entry = selected_diag[key]
        run_rows.append(
            [
                entry.setup,
                entry.piece,
                entry.helicity,
                entry.run_name,
                str(entry.n_events) if entry.n_events is not None else "n/a",
                f"{entry.values.get('F_total', float('nan')):.6f}" if "F_total" in entry.values else "missing",
                entry.path.name,
            ]
        )
    lines.append("Selected diagnostic runs:")
    lines.extend(render_table(["Setup", "Piece", "Hel", "Run", "n", "F_total", "Log file"], run_rows, aligns=("l", "l", "l", "l", "r", "r", "l")))
    lines.append("")

    if any("A_born0" in entry.values for entry in selected_diag.values()):
        a_rows = []
        for key in sorted(selected_diag, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
            entry = selected_diag[key]
            a_rows.append(
                [
                    entry.setup,
                    entry.piece,
                    entry.helicity,
                    fmt_unitless(entry.values.get("A_born0")),
                    fmt_unitless(entry.values.get("A_born_even")),
                    fmt_unitless(entry.values.get("A_born_odd")),
                    fmt_unitless(entry.values.get("A_qmap0")),
                    fmt_unitless(entry.values.get("A_qmap_even")),
                    fmt_unitless(entry.values.get("A_qmap_odd")),
                ]
            )
        lines.append("Weighted A_pol decomposition (Born and QCDC map):")
        lines.extend(
            render_table(
                ["Setup", "Piece", "Hel", "A_born0", "A_b_even", "A_b_odd", "A_q0", "A_q_even", "A_q_odd"],
                a_rows,
                aligns=("l", "l", "l", "r", "r", "r", "r", "r", "r"),
            )
        )
        lines.append("")

        g_rows = []
        for key in sorted(selected_diag, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
            entry = selected_diag[key]
            g_rows.append(
                [
                    entry.setup,
                    entry.piece,
                    entry.helicity,
                    fmt_unitless(entry.values.get("A_gR20")),
                    fmt_unitless(entry.values.get("A_gR2_even")),
                    fmt_unitless(entry.values.get("A_gR2_odd")),
                    fmt_unitless(entry.values.get("A_gR30")),
                    fmt_unitless(entry.values.get("A_gR3_even")),
                    fmt_unitless(entry.values.get("A_gR3_odd")),
                ]
            )
        lines.append("Weighted A_pol decomposition (BGF maps):")
        lines.extend(
            render_table(
                ["Setup", "Piece", "Hel", "A_gR20", "A_gR2_even", "A_gR2_odd", "A_gR30", "A_gR3_even", "A_gR3_odd"],
                g_rows,
                aligns=("l", "l", "l", "r", "r", "r", "r", "r", "r"),
            )
        )
        lines.append("")

    if any("B_closure" in entry.values for entry in selected_diag.values()):
        born_rows = []
        for key in sorted(selected_diag, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
            entry = selected_diag[key]
            born_rows.append(
                [
                    entry.setup,
                    entry.piece,
                    entry.helicity,
                    fmt_unitless(entry.values.get("B_me")),
                    fmt_unitless(entry.values.get("B_pred")),
                    fmt_unitless(entry.values.get("B_closure")),
                    fmt_unitless(entry.values.get("M_me")),
                    fmt_unitless(entry.values.get("M_pred")),
                    fmt_unitless(entry.values.get("M_closure")),
                    fmt_unitless(entry.values.get("ME_self")),
                ]
            )
        lines.append("Born/QCDC closure diagnostics:")
        lines.extend(
            render_table(
                ["Setup", "Piece", "Hel", "B_me", "B_pred", "B_close", "M_me", "M_pred", "M_close", "ME_self"],
                born_rows,
                aligns=("l", "l", "l", "r", "r", "r", "r", "r", "r", "r"),
            )
        )
        lines.append("")

    if any("CPl_pred" in entry.values for entry in selected_diag.values()):
        coeff_rows = []
        for key in sorted(selected_diag, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
            entry = selected_diag[key]
            coeff_rows.append(
                [
                    entry.setup,
                    entry.piece,
                    entry.helicity,
                    fmt_unitless(entry.values.get("C00_scale")),
                    fmt_unitless(entry.values.get("CPl_me")),
                    fmt_unitless(entry.values.get("CPl_pred")),
                    fmt_unitless(entry.values.get("CPl_close")),
                    fmt_unitless(entry.values.get("CPq_me")),
                    fmt_unitless(entry.values.get("CPq_pred")),
                    fmt_unitless(entry.values.get("CPq_close")),
                    fmt_unitless(entry.values.get("CPlPq_me")),
                    fmt_unitless(entry.values.get("CPlPq_pred")),
                    fmt_unitless(entry.values.get("CPlPq_close")),
                ]
            )
        lines.append("Direct me2 Born-basis coefficients:")
        lines.extend(
            render_table(
                ["Setup", "Piece", "Hel", "C00_scale", "CPl_me", "CPl_pred", "CPl_close",
                 "CPq_me", "CPq_pred", "CPq_close", "CPlPq_me", "CPlPq_pred", "CPlPq_close"],
                coeff_rows,
                aligns=("l", "l", "l", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r"),
            )
        )
        lines.append("")

    out_rows = []
    for key in sorted(selected_out, key=lambda item: (SETUP_ORDER.get(item[0], 99), item[0], item[1], item[2])):
        out = selected_out[key]
        if out.piece == "LO":
            continue
        out_rows.append([out.setup, out.piece, out.helicity, fmt_measurement(out.measurement), out.path.name])
    if out_rows:
        lines.append("Selected run cross sections:")
        lines.extend(render_table(["Setup", "Piece", "Helicity", "sigma_run", "Out file"], out_rows, aligns=("l", "l", "l", "r", "l")))
        lines.append("")

    primary_estimates = compute_all_estimates(selected_diag, selected_out, PRIMARY_COMPONENTS)
    split_estimates = compute_all_estimates(selected_diag, selected_out, REAL_SPLIT_COMPONENTS)
    for setup in sorted(primary_estimates, key=lambda s: (SETUP_ORDER.get(s, 99), s)):
        lines.append("=" * 92)
        lines.append(f"Setup: {setup}")
        lines.append("=" * 92)
        for piece in ("POSNLO", "NEGNLO", "NLO"):
            rows = []
            for component in PRIMARY_COMPONENTS:
                vals = primary_estimates[setup][piece][component]
                rows.append([component, fmt_pb(vals.get("unpol_pb")), fmt_pb(vals.get("pol_pb"))])
            lines.append(f"{piece} component estimates:")
            lines.extend(render_table(["Component", "Unpolarized", "Polarized"], rows, aligns=("l", "r", "r")))

            split_rows = []
            for component in REAL_SPLIT_COMPONENTS:
                vals = split_estimates[setup][piece][component]
                split_rows.append([component, fmt_pb(vals.get("unpol_pb")), fmt_pb(vals.get("pol_pb"))])
            lines.append(f"{piece} real-kernel split:")
            lines.extend(render_table(["Component", "Unpolarized", "Polarized"], split_rows, aligns=("l", "r", "r")))
        lines.append("")

    primary_interference = interference_estimates(primary_estimates, PRIMARY_COMPONENTS)
    split_interference = interference_estimates(split_estimates, REAL_SPLIT_COMPONENTS)
    if primary_interference:
        lines.append("=" * 92)
        lines.append("Interference: ALL - GAMMA - Z (NLO)")
        lines.append("=" * 92)
        rows = []
        for component in PRIMARY_COMPONENTS:
            vals = primary_interference[component]
            rows.append([component, fmt_pb(vals.get("unpol_pb")), fmt_pb(vals.get("pol_pb"))])
        lines.extend(render_table(["Component", "Unpolarized", "Polarized"], rows, aligns=("l", "r", "r")))
        lines.append("Real-kernel split:")
        split_rows = []
        for component in REAL_SPLIT_COMPONENTS:
            vals = split_interference[component]
            split_rows.append([component, fmt_pb(vals.get("unpol_pb")), fmt_pb(vals.get("pol_pb"))])
        lines.extend(render_table(["Component", "Unpolarized", "Polarized"], split_rows, aligns=("l", "r", "r")))
        lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def json_payload(
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
    diag_duplicates: Dict[Tuple[str, str, str], List[Path]],
    selected_out: Dict[Tuple[str, str, str], OutRun],
    out_duplicates: Dict[Tuple[str, str, str], List[Path]],
) -> Dict[str, object]:
    primary_estimates = compute_all_estimates(selected_diag, selected_out, PRIMARY_COMPONENTS)
    split_estimates = compute_all_estimates(selected_diag, selected_out, REAL_SPLIT_COMPONENTS)
    return {
        "selected_runs": [
            {
                "setup": entry.setup,
                "piece": entry.piece,
                "helicity": entry.helicity,
                "run_name": entry.run_name,
                "variant": entry.variant,
                "n_events": entry.n_events,
                "path": str(entry.path),
                "values": entry.values,
            }
            for _, entry in sorted(selected_diag.items(), key=lambda item: (SETUP_ORDER.get(item[0][0], 99), item[0][0], item[0][1], item[0][2]))
        ],
        "diag_duplicates": {
            f"{setup}:{piece}:{helicity}": [str(path) for path in paths]
            for (setup, piece, helicity), paths in diag_duplicates.items()
        },
        "selected_out": [
            {
                "setup": out.setup,
                "piece": out.piece,
                "helicity": out.helicity,
                "variant": out.variant,
                "path": str(out.path),
                "measurement": {"value_pb": out.measurement.value_pb, "error_pb": out.measurement.error_pb},
            }
            for _, out in sorted(selected_out.items(), key=lambda item: (SETUP_ORDER.get(item[0][0], 99), item[0][0], item[0][1], item[0][2]))
        ],
        "out_duplicates": {
            f"{setup}:{piece}:{helicity}": [str(path) for path in paths]
            for (setup, piece, helicity), paths in out_duplicates.items()
        },
        "primary_estimates": primary_estimates,
        "real_split_estimates": split_estimates,
        "primary_interference_nlo": interference_estimates(primary_estimates, PRIMARY_COMPONENTS),
        "real_split_interference_nlo": interference_estimates(split_estimates, REAL_SPLIT_COMPONENTS),
    }


def csv_rows(
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
    selected_out: Dict[Tuple[str, str, str], OutRun],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    primary_estimates = compute_all_estimates(selected_diag, selected_out, PRIMARY_COMPONENTS)
    split_estimates = compute_all_estimates(selected_diag, selected_out, REAL_SPLIT_COMPONENTS)

    for _, entry in sorted(selected_diag.items(), key=lambda item: (SETUP_ORDER.get(item[0][0], 99), item[0][0], item[0][1], item[0][2])):
        row = {
            "row_type": "selected_run",
            "setup": entry.setup,
            "piece": entry.piece,
            "helicity": entry.helicity,
            "component": "",
            "observable": "",
            "value_pb": "",
            "run_name": entry.run_name,
            "variant": entry.variant,
            "n_events": entry.n_events,
            "path": str(entry.path),
        }
        row.update(entry.values)
        rows.append(row)

    for _, out in sorted(selected_out.items(), key=lambda item: (SETUP_ORDER.get(item[0][0], 99), item[0][0], item[0][1], item[0][2])):
        rows.append(
            {
                "row_type": "selected_out",
                "setup": out.setup,
                "piece": out.piece,
                "helicity": out.helicity,
                "component": "",
                "observable": "sigma_run_pb",
                "value_pb": out.measurement.value_pb,
                "error_pb": out.measurement.error_pb,
                "run_name": "",
                "variant": out.variant,
                "n_events": "",
                "path": str(out.path),
            }
        )

    for name, estimate_map in (("primary_estimate", primary_estimates), ("real_split_estimate", split_estimates)):
        for setup, piece_map in estimate_map.items():
            for piece, comp_map in piece_map.items():
                for component, values in comp_map.items():
                    for observable in ("unpol_pb", "pol_pb"):
                        rows.append(
                            {
                                "row_type": name,
                                "setup": setup,
                                "piece": piece,
                                "helicity": "",
                                "component": component,
                                "observable": observable,
                                "value_pb": values.get(observable),
                                "run_name": "",
                                "variant": "",
                                "n_events": "",
                                "path": "",
                            }
                        )

    for name, comp_map in (
        ("primary_interference_nlo", interference_estimates(primary_estimates, PRIMARY_COMPONENTS)),
        ("real_split_interference_nlo", interference_estimates(split_estimates, REAL_SPLIT_COMPONENTS)),
    ):
        for component, values in comp_map.items():
            for observable in ("unpol_pb", "pol_pb"):
                rows.append(
                    {
                        "row_type": name,
                        "setup": "INTERFERENCE",
                        "piece": "NLO",
                        "helicity": "",
                        "component": component,
                        "observable": observable,
                        "value_pb": values.get(observable),
                        "run_name": "",
                        "variant": "",
                        "n_events": "",
                        "path": "",
                    }
                )

    return rows


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--base-dir",
        default=str(Path(__file__).resolve().parent),
        help="Directory to search recursively for POWHEG .log and matching .out files.",
    )
    parser.add_argument(
        "-t",
        "--tag",
        default="",
        help="Prefer runs with this variant tag (for example: testing).",
    )
    parser.add_argument(
        "--strict-tag",
        action="store_true",
        help="Require the requested tag exactly and do not fall back to untagged or other tagged runs.",
    )
    parser.add_argument(
        "--rivet",
        action="store_true",
        help="Resolve and report only -RIVET logs/outputs.",
    )
    parser.add_argument(
        "--rivetfo",
        action="store_true",
        help="Resolve and report only -RIVETFO logs/outputs.",
    )
    parser.add_argument(
        "--rivetfofixed",
        action="store_true",
        help="Resolve and report only -RIVETFOFIXED logs/outputs.",
    )
    parser.add_argument("--json-out", help="Optional path to write a structured JSON summary.")
    parser.add_argument("--csv-out", help="Optional path to write a flattened CSV summary.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if int(args.rivet) + int(args.rivetfo) + int(args.rivetfofixed) > 1:
        raise SystemExit("Use at most one of --rivet, --rivetfo, and --rivetfofixed.")
    base_dir = Path(args.base_dir).resolve()

    diag_entries = [entry for path in sorted(base_dir.rglob("DIS-POL-POWHEG*.log")) if (entry := parse_diagnostic_log(path)) is not None]
    out_entries = [entry for path in sorted(base_dir.rglob("DIS-POL-*.out")) if (entry := parse_out(path)) is not None]
    requested_analysis = (
        "RIVETFOFIXED" if args.rivetfofixed else
        ("RIVETFO" if args.rivetfo else ("RIVET" if args.rivet else ""))
    )
    diag_entries = [entry for entry in diag_entries if entry.analysis == requested_analysis]
    out_entries = [entry for entry in out_entries if entry.analysis == requested_analysis]

    selected_out, out_duplicates = select_best_lo(out_entries, args.tag, args.strict_tag)
    selected_diag, diag_duplicates = select_best_runs(diag_entries, out_entries, args.tag, args.strict_tag)

    report = build_report(selected_diag, diag_duplicates, selected_out, out_duplicates, args.tag, args.strict_tag)
    print(report, end="")

    if args.json_out:
        payload = json_payload(selected_diag, diag_duplicates, selected_out, out_duplicates)
        out = Path(args.json_out).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"Wrote JSON: {out}")

    if args.csv_out:
        rows = csv_rows(selected_diag, selected_out)
        out = Path(args.csv_out).resolve()
        out.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = sorted({key for row in rows for key in row})
        with out.open("w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow(row)
        print(f"Wrote CSV: {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
