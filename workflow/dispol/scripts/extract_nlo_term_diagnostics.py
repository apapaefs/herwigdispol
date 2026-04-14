#!/usr/bin/env python3.10
"""
Extract component-level NLO diagnostics from DIS POWHEG log files.

The script looks for the diagnostic lines emitted by DISBase::NLOWeight():

  NLO_TERM_SIGN
  NLO_TERM_CUM
  NLO_TERM_REAL
  NLO_TERM_ZSPIN
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
import math
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
# Support both plain shard variants like
#   S700560-plain48-z-termdiag-s001
# and prefixed variants like
#   TERMDIAG-SP040-S710000-plain49-gamma-termdiag-power-s001
SEEDED_VARIANT_RE = re.compile(r"(?:^|-)S\d+-(?P<tag>.+)$")
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
RAW_PROBE_COMPONENTS = (
    "F_cq_odd_raw_probe",
    "F_cq_odd_raw_shift",
    "F_rq_odd_raw_probe",
    "F_rq_odd_raw_shift",
)
Z_SHADOW_COMPONENTS = (
    "cq_odd_raw_probe_pb",
    "cq_odd_raw_shift_pb",
    "rq_odd_raw_probe_pb",
    "rq_odd_raw_shift_pb",
    "q_odd_raw_shift_pb",
)
Z_SPIN_COMPONENTS = (
    "virt_pb",
    "cq_even_pb",
    "cq_odd_pb",
    "cg_even_pb",
    "cg_odd_pb",
    "rq_even_pb",
    "rq_odd_pb",
    "rg_even_pb",
    "rg_odd_pb",
    "q_even_pb",
    "q_odd_pb",
    "g_even_pb",
    "g_odd_pb",
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
    generated_events: Optional[int]
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


def parse_generated_event_count(text: str) -> Optional[int]:
    matches = list(TOTAL_RE.finditer(text))
    if not matches:
        return None
    try:
        return int(matches[-1].group(1))
    except ValueError:
        return None


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
    match = SEEDED_VARIANT_RE.search(variant)
    if match:
        return match.group("tag")
    return variant


def parse_filter_values(values: Optional[Sequence[str]], uppercase: bool = False) -> Optional[set[str]]:
    if not values:
        return None
    out: set[str] = set()
    for value in values:
        for item in value.split(","):
            token = item.strip()
            if not token:
                continue
            out.add(token.upper() if uppercase else token)
    return out or None


def name_matches_filter(path: Path, run_name: str, needle: Optional[str]) -> bool:
    if not needle:
        return True
    needle_lower = needle.lower()
    haystacks = (
        path.name.lower(),
        path.stem.lower(),
        run_name.lower(),
    )
    return any(needle_lower in haystack for haystack in haystacks)


def entry_matches_filters(
    entry: object,
    setup_filter: Optional[set[str]],
    piece_filter: Optional[set[str]],
    helicity_filter: Optional[set[str]],
    name_filter: Optional[str],
) -> bool:
    setup = getattr(entry, "setup")
    piece = getattr(entry, "piece")
    helicity = getattr(entry, "helicity")
    path = getattr(entry, "path")
    run_name = getattr(entry, "run_name", path.stem)
    if setup_filter is not None and setup not in setup_filter:
        return False
    if piece_filter is not None and piece not in piece_filter:
        return False
    if helicity_filter is not None and helicity not in helicity_filter:
        return False
    return name_matches_filter(path, run_name, name_filter)


def path_matches_prefilters(
    path: Path,
    requested_analysis: str,
    setup_filter: Optional[set[str]],
    piece_filter: Optional[set[str]],
    helicity_filter: Optional[set[str]],
    name_filter: Optional[str],
) -> bool:
    if name_filter and not name_matches_filter(path, path.stem, name_filter):
        return False
    parsed = parse_run_name(path.stem)
    if parsed is None:
        return False
    setup, piece, helicity, analysis, _variant = parsed
    if analysis != requested_analysis:
        return False
    if setup_filter is not None and setup not in setup_filter:
        return False
    if piece_filter is not None and piece not in piece_filter:
        return False
    if helicity_filter is not None and helicity not in helicity_filter:
        return False
    return True


def parse_key_value_fields(payload: str) -> Dict[str, str]:
    fields: Dict[str, str] = {}
    for token in payload.strip().split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = value
    return fields


def parse_diagnostic_log(path: Path) -> Optional[DiagnosticRun]:
    last_sign: Optional[Dict[str, str]] = None
    last_cum: Optional[Dict[str, str]] = None
    last_real: Optional[Dict[str, str]] = None
    last_zspin: Optional[Dict[str, str]] = None
    last_apol: Optional[Dict[str, str]] = None
    last_born: Optional[Dict[str, str]] = None
    last_coeff: Optional[Dict[str, str]] = None
    for line in path.read_text().splitlines():
        if line.startswith("NLO_TERM_SIGN"):
            last_sign = parse_key_value_fields(line[len("NLO_TERM_SIGN"):])
        elif line.startswith("NLO_TERM_CUM"):
            last_cum = parse_key_value_fields(line[len("NLO_TERM_CUM"):])
        elif line.startswith("NLO_TERM_REAL"):
            last_real = parse_key_value_fields(line[len("NLO_TERM_REAL"):])
        elif line.startswith("NLO_TERM_ZSPIN"):
            last_zspin = parse_key_value_fields(line[len("NLO_TERM_ZSPIN"):])
        elif line.startswith("NLO_TERM_APOL"):
            last_apol = parse_key_value_fields(line[len("NLO_TERM_APOL"):])
        elif line.startswith("NLO_TERM_BORN"):
            last_born = parse_key_value_fields(line[len("NLO_TERM_BORN"):])
        elif line.startswith("NLO_TERM_COEFF"):
            last_coeff = parse_key_value_fields(line[len("NLO_TERM_COEFF"):])

    if (
        last_sign is None
        and last_cum is None
        and last_real is None
        and last_zspin is None
        and last_apol is None
        and last_born is None
        and last_coeff is None
    ):
        return None

    source = last_sign or last_cum or last_real or last_zspin or last_apol or last_born or last_coeff
    assert source is not None
    run_name = source.get("run", path.stem)
    parsed = parse_run_name(run_name)
    if parsed is None:
        return None
    setup, piece, helicity, analysis, variant = parsed

    values: Dict[str, float] = {}
    n_events: Optional[int] = None
    for payload in (last_sign, last_cum, last_real, last_zspin, last_apol, last_born, last_coeff):
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
    text = path.read_text()
    measurement = parse_measurement(text)
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
        generated_events=parse_generated_event_count(text),
        measurement=measurement,
    )


def run_preference(variant: str, path: Path, preferred_tag: str) -> Tuple[int, int, float, str]:
    preferred_penalty = 0 if variant_matches_tag(variant, preferred_tag) else 1
    default_penalty = 0 if normalize_variant_tag(variant) == "" else 1
    if not preferred_tag:
        preferred_penalty = 0
    return (preferred_penalty, default_penalty, -path.stat().st_mtime, path.name)


def shard_mean_weights(generated_events: Optional[Sequence[Optional[int]]], count: int) -> List[float]:
    if generated_events is None:
        return [1.0] * count
    events = list(generated_events)
    if len(events) != count:
        return [1.0] * count
    weights: List[float] = []
    for raw in events:
        if raw is None:
            return [1.0] * count
        try:
            weight = float(raw)
        except (TypeError, ValueError):
            return [1.0] * count
        if weight <= 0.0:
            return [1.0] * count
        weights.append(weight)
    return weights


def combine_measurements(
    measurements: Sequence[Measurement],
    generated_events: Optional[Sequence[Optional[int]]] = None,
) -> Measurement:
    if not measurements:
        raise ValueError("combine_measurements() requires at least one measurement")

    weights = shard_mean_weights(generated_events, len(measurements))
    total_weight = sum(weights)
    value = sum(w * m.value_pb for w, m in zip(weights, measurements)) / total_weight
    error = math.sqrt(sum((w * m.error_pb) ** 2 for w, m in zip(weights, measurements))) / total_weight
    return Measurement(value, error)


def diagnostic_shard_weight(item: DiagnosticRun, out: Optional[OutRun]) -> float:
    if item.n_events is not None and item.n_events > 0:
        return float(item.n_events)
    if out is not None and out.generated_events is not None and out.generated_events > 0:
        return float(out.generated_events)
    return 1.0


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

    measurement = combine_measurements(
        [item.measurement for item in chosen],
        [item.generated_events for item in chosen],
    )
    paths = sorted((item.path for item in chosen), key=lambda item: item.name)
    return OutRun(
        path=paths[0],
        paths=paths,
        setup=chosen[0].setup,
        piece=chosen[0].piece,
        helicity=chosen[0].helicity,
        analysis=chosen[0].analysis,
        variant=preferred_tag or chosen[0].variant,
        generated_events=sum(item.generated_events or 0 for item in chosen) or None,
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
        weight = diagnostic_shard_weight(item, out)
        weighted_entries.append((weight, item))
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
        analysis=chosen[0].analysis,
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


def sum_optional(*values: Optional[float]) -> Optional[float]:
    if any(value is None for value in values):
        return None
    return sum(value for value in values if value is not None)


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


def z_spin_component_estimates(
    primary_estimates: Dict[str, Dict[str, Dict[str, Dict[str, Optional[float]]]]],
    split_estimates: Dict[str, Dict[str, Dict[str, Dict[str, Optional[float]]]]],
) -> Dict[str, Dict[str, Optional[float]]]:
    if "Z" not in primary_estimates or "Z" not in split_estimates:
        return {}

    component_map: Dict[str, Dict[str, Optional[float]]] = {}
    for piece in ("POSNLO", "NEGNLO", "NLO"):
        primary_piece = primary_estimates["Z"][piece]
        split_piece = split_estimates["Z"][piece]
        piece_out: Dict[str, Optional[float]] = {
            "virt_pb": primary_piece["F_virt"].get("pol_pb"),
            "cq_even_pb": primary_piece["F_cq_even"].get("pol_pb"),
            "cq_odd_pb": primary_piece["F_cq_odd"].get("pol_pb"),
            "cg_even_pb": primary_piece["F_cg_even"].get("pol_pb"),
            "cg_odd_pb": primary_piece["F_cg_odd"].get("pol_pb"),
            "rq_even_pb": split_piece["F_rq_even"].get("pol_pb"),
            "rq_odd_pb": split_piece["F_rq_odd"].get("pol_pb"),
            "rg_even_pb": split_piece["F_rg_even"].get("pol_pb"),
            "rg_odd_pb": split_piece["F_rg_odd"].get("pol_pb"),
        }
        piece_out["q_even_pb"] = sum_optional(piece_out["cq_even_pb"], piece_out["rq_even_pb"])
        piece_out["q_odd_pb"] = sum_optional(piece_out["cq_odd_pb"], piece_out["rq_odd_pb"])
        piece_out["g_even_pb"] = sum_optional(piece_out["cg_even_pb"], piece_out["rg_even_pb"])
        piece_out["g_odd_pb"] = sum_optional(piece_out["cg_odd_pb"], piece_out["rg_odd_pb"])
        component_map[piece] = piece_out
    return component_map


def z_shadow_component_estimates(
    raw_probe_estimates: Dict[str, Dict[str, Dict[str, Dict[str, Optional[float]]]]],
) -> Dict[str, Dict[str, Optional[float]]]:
    if "Z" not in raw_probe_estimates:
        return {}

    component_map: Dict[str, Dict[str, Optional[float]]] = {}
    for piece in ("POSNLO", "NEGNLO", "NLO"):
        split_piece = raw_probe_estimates["Z"][piece]
        piece_out: Dict[str, Optional[float]] = {
            "cq_odd_raw_probe_pb": split_piece["F_cq_odd_raw_probe"].get("pol_pb"),
            "cq_odd_raw_shift_pb": split_piece["F_cq_odd_raw_shift"].get("pol_pb"),
            "rq_odd_raw_probe_pb": split_piece["F_rq_odd_raw_probe"].get("pol_pb"),
            "rq_odd_raw_shift_pb": split_piece["F_rq_odd_raw_shift"].get("pol_pb"),
        }
        piece_out["q_odd_raw_shift_pb"] = sum_optional(
            piece_out["cq_odd_raw_shift_pb"], piece_out["rq_odd_raw_shift_pb"]
        )
        component_map[piece] = piece_out
    return component_map


def dominant_clamp_sector(clip_q_absw_frac: float, clip_qm_absw_frac: float, clip_gm_absw_frac: float) -> str:
    sector_values = {
        "q_born": clip_q_absw_frac,
        "q_mapped": clip_qm_absw_frac,
        "g_mapped": clip_gm_absw_frac,
    }
    ordered = sorted(sector_values.items(), key=lambda item: item[1], reverse=True)
    top_name, top_value = ordered[0]
    second_value = ordered[1][1]
    if top_value <= 1e-12:
        return "none"
    if second_value > 1e-12 and abs(top_value - second_value) <= 0.05 * top_value:
        return "mixed"
    return top_name


def aggregate_z_clamp_summary(
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
) -> Dict[str, Dict[str, object]]:
    summaries: Dict[str, Dict[str, object]] = {}
    for piece in ("POSNLO", "NEGNLO"):
        n_eval_total = 0.0
        w_abs_total = 0.0
        clip_q_event_num = 0.0
        clip_q_absw_num = 0.0
        clip_qm_event_num = 0.0
        clip_qm_absw_num = 0.0
        clip_gm_event_num = 0.0
        clip_gm_absw_num = 0.0
        pq_raw_sum = pq_sum = pq_clip_sum = 0.0
        pqm_raw_sum = pqm_sum = pqm_clip_sum = 0.0
        pgm_raw_sum = pgm_sum = pgm_clip_sum = 0.0
        aq_raw_sum = aq_sum = aq_shift_sum = 0.0
        max_abs_pq_raw = max_abs_pq_clip = 0.0
        max_abs_pqm_raw = max_abs_pqm_clip = 0.0
        max_abs_pgm_raw = max_abs_pgm_clip = 0.0
        max_abs_aq_shift = 0.0

        for helicity in ("PP", "PM", "MP", "MM"):
            entry = selected_diag.get(("Z", piece, helicity))
            if entry is None or "w_abs" not in entry.values or "n_eval" not in entry.values:
                continue
            weight = entry.values["w_abs"]
            n_eval = entry.values["n_eval"]
            w_abs_total += weight
            n_eval_total += n_eval

            clip_q_event_num += entry.values.get("clipQ_event_frac", 0.0) * n_eval
            clip_q_absw_num += entry.values.get("clipQ_absw_frac", 0.0) * weight
            clip_qm_event_num += entry.values.get("clipQm_event_frac", 0.0) * n_eval
            clip_qm_absw_num += entry.values.get("clipQm_absw_frac", 0.0) * weight
            clip_gm_event_num += entry.values.get("clipGm_event_frac", 0.0) * n_eval
            clip_gm_absw_num += entry.values.get("clipGm_absw_frac", 0.0) * weight

            pq_raw_sum += entry.values.get("Pq_raw_absw_mean", 0.0) * weight
            pq_sum += entry.values.get("Pq_absw_mean", 0.0) * weight
            pq_clip_sum += entry.values.get("Pq_clip_absw_mean", 0.0) * weight
            pqm_raw_sum += entry.values.get("Pq_m_raw_absw_mean", 0.0) * weight
            pqm_sum += entry.values.get("Pq_m_absw_mean", 0.0) * weight
            pqm_clip_sum += entry.values.get("Pq_m_clip_absw_mean", 0.0) * weight
            pgm_raw_sum += entry.values.get("Pg_m_raw_absw_mean", 0.0) * weight
            pgm_sum += entry.values.get("Pg_m_absw_mean", 0.0) * weight
            pgm_clip_sum += entry.values.get("Pg_m_clip_absw_mean", 0.0) * weight
            aq_raw_sum += entry.values.get("A_qmap_raw_absw_mean", 0.0) * weight
            aq_sum += entry.values.get("A_qmap_absw_mean", 0.0) * weight
            aq_shift_sum += entry.values.get("A_qmap_shift_absw_mean", 0.0) * weight

            max_abs_pq_raw = max(max_abs_pq_raw, entry.values.get("maxAbsPq_raw", 0.0))
            max_abs_pq_clip = max(max_abs_pq_clip, entry.values.get("maxAbsPq_clip", 0.0))
            max_abs_pqm_raw = max(max_abs_pqm_raw, entry.values.get("maxAbsPq_m_raw", 0.0))
            max_abs_pqm_clip = max(max_abs_pqm_clip, entry.values.get("maxAbsPq_m_clip", 0.0))
            max_abs_pgm_raw = max(max_abs_pgm_raw, entry.values.get("maxAbsPg_m_raw", 0.0))
            max_abs_pgm_clip = max(max_abs_pgm_clip, entry.values.get("maxAbsPg_m_clip", 0.0))
            max_abs_aq_shift = max(max_abs_aq_shift, entry.values.get("maxAbsA_qmap_raw_shift", 0.0))

        if n_eval_total <= 0.0 or w_abs_total <= 0.0:
            continue

        clip_q_event_frac = clip_q_event_num / n_eval_total
        clip_q_absw_frac = clip_q_absw_num / w_abs_total
        clip_qm_event_frac = clip_qm_event_num / n_eval_total
        clip_qm_absw_frac = clip_qm_absw_num / w_abs_total
        clip_gm_event_frac = clip_gm_event_num / n_eval_total
        clip_gm_absw_frac = clip_gm_absw_num / w_abs_total

        summaries[piece] = {
            "n_eval": n_eval_total,
            "w_abs": w_abs_total,
            "clipQ_event_frac": clip_q_event_frac,
            "clipQ_absw_frac": clip_q_absw_frac,
            "clipQm_event_frac": clip_qm_event_frac,
            "clipQm_absw_frac": clip_qm_absw_frac,
            "clipGm_event_frac": clip_gm_event_frac,
            "clipGm_absw_frac": clip_gm_absw_frac,
            "Pq_raw_absw_mean": pq_raw_sum / w_abs_total,
            "Pq_absw_mean": pq_sum / w_abs_total,
            "Pq_clip_absw_mean": pq_clip_sum / w_abs_total,
            "Pq_m_raw_absw_mean": pqm_raw_sum / w_abs_total,
            "Pq_m_absw_mean": pqm_sum / w_abs_total,
            "Pq_m_clip_absw_mean": pqm_clip_sum / w_abs_total,
            "Pg_m_raw_absw_mean": pgm_raw_sum / w_abs_total,
            "Pg_m_absw_mean": pgm_sum / w_abs_total,
            "Pg_m_clip_absw_mean": pgm_clip_sum / w_abs_total,
            "A_qmap_raw_absw_mean": aq_raw_sum / w_abs_total,
            "A_qmap_absw_mean": aq_sum / w_abs_total,
            "A_qmap_shift_absw_mean": aq_shift_sum / w_abs_total,
            "maxAbsPq_raw": max_abs_pq_raw,
            "maxAbsPq_clip": max_abs_pq_clip,
            "maxAbsPq_m_raw": max_abs_pqm_raw,
            "maxAbsPq_m_clip": max_abs_pqm_clip,
            "maxAbsPg_m_raw": max_abs_pgm_raw,
            "maxAbsPg_m_clip": max_abs_pgm_clip,
            "maxAbsA_qmap_raw_shift": max_abs_aq_shift,
            "dominant_clamp_sector": dominant_clamp_sector(
                clip_q_absw_frac, clip_qm_absw_frac, clip_gm_absw_frac
            ),
        }
    return summaries


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
    z_spin_report: bool,
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
    raw_probe_estimates = compute_all_estimates(selected_diag, selected_out, RAW_PROBE_COMPONENTS)
    z_spin_estimates = z_spin_component_estimates(primary_estimates, split_estimates)
    z_shadow_estimates = z_shadow_component_estimates(raw_probe_estimates)
    z_clamp_summary = aggregate_z_clamp_summary(selected_diag)

    if z_spin_report and z_spin_estimates:
        lines.append("=" * 92)
        lines.append("Z sigma_LL component summary")
        lines.append("=" * 92)
        rows = []
        for component in Z_SPIN_COMPONENTS:
            rows.append(
                [
                    component,
                    fmt_pb(z_spin_estimates.get("POSNLO", {}).get(component)),
                    fmt_pb(z_spin_estimates.get("NEGNLO", {}).get(component)),
                    fmt_pb(z_spin_estimates.get("NLO", {}).get(component)),
                ]
            )
        lines.extend(
            render_table(
                ["Component", "POSNLO sigma_LL", "NEGNLO sigma_LL", "NLO sigma_LL"],
                rows,
                aligns=("l", "r", "r", "r"),
            )
        )
        lines.append("")

    if z_spin_report and z_shadow_estimates:
        lines.append("=" * 92)
        lines.append("Z raw-vs-clamped odd-sector shadow summary")
        lines.append("=" * 92)
        rows = []
        for component in Z_SHADOW_COMPONENTS:
            rows.append(
                [
                    component,
                    fmt_pb(z_shadow_estimates.get("POSNLO", {}).get(component)),
                    fmt_pb(z_shadow_estimates.get("NEGNLO", {}).get(component)),
                    fmt_pb(z_shadow_estimates.get("NLO", {}).get(component)),
                ]
            )
        lines.extend(
            render_table(
                ["Component", "POSNLO sigma_LL", "NEGNLO sigma_LL", "NLO sigma_LL"],
                rows,
                aligns=("l", "r", "r", "r"),
            )
        )
        lines.append("")

    if z_spin_report and z_clamp_summary:
        lines.append("=" * 92)
        lines.append("Z clamp summary")
        lines.append("=" * 92)
        frac_rows = []
        mean_rows = []
        max_rows = []
        for piece in ("POSNLO", "NEGNLO"):
            summary = z_clamp_summary.get(piece)
            if summary is None:
                continue
            frac_rows.append(
                [
                    piece,
                    f"{summary['n_eval']:.0f}",
                    f"{summary['w_abs']:.6e}",
                    fmt_unitless(summary["clipQ_event_frac"]),
                    fmt_unitless(summary["clipQ_absw_frac"]),
                    fmt_unitless(summary["clipQm_event_frac"]),
                    fmt_unitless(summary["clipQm_absw_frac"]),
                    fmt_unitless(summary["clipGm_event_frac"]),
                    fmt_unitless(summary["clipGm_absw_frac"]),
                    str(summary["dominant_clamp_sector"]),
                ]
            )
            mean_rows.append(
                [
                    piece,
                    fmt_unitless(summary["Pq_raw_absw_mean"]),
                    fmt_unitless(summary["Pq_absw_mean"]),
                    fmt_unitless(summary["Pq_clip_absw_mean"]),
                    fmt_unitless(summary["Pq_m_raw_absw_mean"]),
                    fmt_unitless(summary["Pq_m_absw_mean"]),
                    fmt_unitless(summary["Pq_m_clip_absw_mean"]),
                    fmt_unitless(summary["Pg_m_raw_absw_mean"]),
                    fmt_unitless(summary["Pg_m_absw_mean"]),
                    fmt_unitless(summary["Pg_m_clip_absw_mean"]),
                    fmt_unitless(summary.get("A_qmap_raw_absw_mean")),
                    fmt_unitless(summary.get("A_qmap_absw_mean")),
                    fmt_unitless(summary.get("A_qmap_shift_absw_mean")),
                ]
            )
            max_rows.append(
                [
                    piece,
                    fmt_unitless(summary["maxAbsPq_raw"]),
                    fmt_unitless(summary["maxAbsPq_clip"]),
                    fmt_unitless(summary["maxAbsPq_m_raw"]),
                    fmt_unitless(summary["maxAbsPq_m_clip"]),
                    fmt_unitless(summary["maxAbsPg_m_raw"]),
                    fmt_unitless(summary["maxAbsPg_m_clip"]),
                    fmt_unitless(summary.get("maxAbsA_qmap_raw_shift")),
                ]
            )
        lines.append("Clamp activity:")
        lines.extend(
            render_table(
                [
                    "Piece",
                    "n_eval",
                    "w_abs",
                    "clipQ_evt",
                    "clipQ_absw",
                    "clipQm_evt",
                    "clipQm_absw",
                    "clipGm_evt",
                    "clipGm_absw",
                    "dominant",
                ],
                frac_rows,
                aligns=("l", "r", "r", "r", "r", "r", "r", "r", "r", "l"),
            )
        )
        lines.append("Raw/clamped abs-weighted means:")
        lines.extend(
            render_table(
                [
                    "Piece",
                    "Pq_raw",
                    "Pq",
                    "clipQ",
                    "Pq_m_raw",
                    "Pq_m",
                    "clipQm",
                    "Pg_m_raw",
                    "Pg_m",
                    "clipGm",
                    "A_q_raw",
                    "A_q",
                    "A_q_shift",
                ],
                mean_rows,
                aligns=("l", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r", "r"),
            )
        )
        lines.append("Max absolute raw/clip values:")
        lines.extend(
            render_table(
                [
                    "Piece",
                    "max|Pq_raw|",
                    "max|clipQ|",
                    "max|Pq_m_raw|",
                    "max|clipQm|",
                    "max|Pg_m_raw|",
                    "max|clipGm|",
                    "max|dA_q|",
                ],
                max_rows,
                aligns=("l", "r", "r", "r", "r", "r", "r", "r"),
            )
        )
        lines.append("")

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
    raw_probe_estimates = compute_all_estimates(selected_diag, selected_out, RAW_PROBE_COMPONENTS)
    z_spin_estimates = z_spin_component_estimates(primary_estimates, split_estimates)
    z_shadow_estimates = z_shadow_component_estimates(raw_probe_estimates)
    z_clamp_summary = aggregate_z_clamp_summary(selected_diag)
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
        "raw_probe_estimates": raw_probe_estimates,
        "primary_interference_nlo": interference_estimates(primary_estimates, PRIMARY_COMPONENTS),
        "real_split_interference_nlo": interference_estimates(split_estimates, REAL_SPLIT_COMPONENTS),
        "z_spin_components": z_spin_estimates,
        "z_shadow_components": z_shadow_estimates,
        "z_clamp_summary": z_clamp_summary,
    }


def csv_rows(
    selected_diag: Dict[Tuple[str, str, str], DiagnosticRun],
    selected_out: Dict[Tuple[str, str, str], OutRun],
) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    primary_estimates = compute_all_estimates(selected_diag, selected_out, PRIMARY_COMPONENTS)
    split_estimates = compute_all_estimates(selected_diag, selected_out, REAL_SPLIT_COMPONENTS)
    raw_probe_estimates = compute_all_estimates(selected_diag, selected_out, RAW_PROBE_COMPONENTS)
    z_spin_estimates = z_spin_component_estimates(primary_estimates, split_estimates)
    z_shadow_estimates = z_shadow_component_estimates(raw_probe_estimates)
    z_clamp_summary = aggregate_z_clamp_summary(selected_diag)

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

    for name, estimate_map in (
        ("primary_estimate", primary_estimates),
        ("real_split_estimate", split_estimates),
        ("raw_probe_estimate", raw_probe_estimates),
    ):
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

    for piece, component_map in z_spin_estimates.items():
        for component, value in component_map.items():
            rows.append(
                {
                    "row_type": "z_spin_component",
                    "setup": "Z",
                    "piece": piece,
                    "helicity": "SIGMA_LL",
                    "component": component,
                    "observable": "pol_pb",
                    "value_pb": value,
                    "run_name": "",
                    "variant": "",
                    "n_events": "",
                    "path": "",
                }
            )

    for piece, component_map in z_shadow_estimates.items():
        for component, value in component_map.items():
            rows.append(
                {
                    "row_type": "z_shadow_component",
                    "setup": "Z",
                    "piece": piece,
                    "helicity": "SIGMA_LL",
                    "component": component,
                    "observable": "pol_pb",
                    "value_pb": value,
                    "run_name": "",
                    "variant": "",
                    "n_events": "",
                    "path": "",
                }
            )

    for piece, summary in z_clamp_summary.items():
        row = {
            "row_type": "z_clamp_summary",
            "setup": "Z",
            "piece": piece,
            "helicity": "",
            "component": "",
            "observable": "",
            "value_pb": "",
            "run_name": "",
            "variant": "",
            "n_events": summary.get("n_eval"),
            "path": "",
        }
        row.update(summary)
        rows.append(row)

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
    parser.add_argument(
        "--z-spin-report",
        action="store_true",
        help="Append a focused Z-only sigma_LL component and clamp summary.",
    )
    parser.add_argument(
        "--setup",
        action="append",
        help="Restrict to one or more setups (for example: Z or GAMMA,Z).",
    )
    parser.add_argument(
        "--piece",
        action="append",
        help="Restrict to one or more pieces (for example: POSNLO,NEGNLO).",
    )
    parser.add_argument(
        "--helicity",
        action="append",
        help="Restrict to one or more helicities (for example: PP,PM,MP,MM or 00).",
    )
    parser.add_argument(
        "--name-filter",
        help="Keep only runs whose file name or run name contains this substring.",
    )
    parser.add_argument("--json-out", help="Optional path to write a structured JSON summary.")
    parser.add_argument("--csv-out", help="Optional path to write a flattened CSV summary.")
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if int(args.rivet) + int(args.rivetfo) + int(args.rivetfofixed) > 1:
        raise SystemExit("Use at most one of --rivet, --rivetfo, and --rivetfofixed.")
    base_dir = Path(args.base_dir).resolve()

    requested_analysis = (
        "RIVETFOFIXED" if args.rivetfofixed else
        ("RIVETFO" if args.rivetfo else ("RIVET" if args.rivet else ""))
    )
    setup_filter = parse_filter_values(args.setup, uppercase=True)
    piece_filter = parse_filter_values(args.piece, uppercase=True)
    helicity_filter = parse_filter_values(args.helicity, uppercase=True)
    diag_entries = [
        entry
        for path in sorted(base_dir.rglob("DIS-POL-POWHEG*.log"))
        if path_matches_prefilters(
            path,
            requested_analysis,
            setup_filter,
            piece_filter,
            helicity_filter,
            args.name_filter,
        )
        if (entry := parse_diagnostic_log(path)) is not None
    ]
    out_entries = [
        entry
        for path in sorted(base_dir.rglob("DIS-POL-*.out"))
        if path_matches_prefilters(
            path,
            requested_analysis,
            setup_filter,
            piece_filter,
            helicity_filter,
            args.name_filter,
        )
        if (entry := parse_out(path)) is not None
    ]

    selected_out, out_duplicates = select_best_lo(out_entries, args.tag, args.strict_tag)
    selected_diag, diag_duplicates = select_best_runs(diag_entries, out_entries, args.tag, args.strict_tag)

    report = build_report(
        selected_diag,
        diag_duplicates,
        selected_out,
        out_duplicates,
        args.tag,
        args.strict_tag,
        args.z_spin_report,
    )
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
