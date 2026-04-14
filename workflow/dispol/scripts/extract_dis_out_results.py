#!/usr/bin/env python3.10
"""
Extract and compare DIS polarized cross sections from Herwig .out files.

The default run list matches the filenames provided in the note discussion:

  - LO runs:        DIS-POL-LO_{PP,PM,MP,MM,00}-{GAMMA,Z,ALL}[-RIVET|-RIVETFO].out
  - NLO runs:       DIS-POL-POWHEG_{PP,PM,MP,MM,00}-{POSNLO,NEGNLO}-{GAMMA,Z,ALL}[-RIVET|-RIVETFO][-test].out

For NLO, the physical cross section is

    sigma_NLO = sigma_POSNLO - sigma_NEGNLO

for each helicity configuration. The script computes, whenever the required
inputs are available,

    sigma_pol = (PP - PM) / 2
    sigma_avg = (PP + PM) / 2

and, when all four polarized helicity states are available,

    sigma_0  = (PP + PM + MP + MM) / 4
    sigma_l  = (PP + PM - MP - MM) / 4
    sigma_p  = (PP - PM + MP - MM) / 4
    sigma_LL = (PP + MM - PM - MP) / 4

for LO, POSNLO, NEGNLO, and the total NLO result. It also reports internal
consistency checks such as sigma_avg - sigma_00 and K-factors with respect to
LO when both orders are present for the same electroweak setup.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    from prettytable import PrettyTable
except ImportError:
    PrettyTable = None


DEFAULT_RUNS = [
    "DIS-POL-POWHEG_PP-POSNLO-ALL.out",
    "DIS-POL-POWHEG_00-POSNLO-ALL.out",
    "DIS-POL-POWHEG_PM-POSNLO-ALL.out",
    "DIS-POL-POWHEG_PP-NEGNLO-ALL.out",
    "DIS-POL-POWHEG_00-NEGNLO-ALL.out",
    "DIS-POL-POWHEG_PM-NEGNLO-ALL.out",
    "DIS-POL-POWHEG_MP-POSNLO-ALL.out",
    "DIS-POL-POWHEG_MM-POSNLO-ALL.out",
    "DIS-POL-POWHEG_MP-NEGNLO-ALL.out",
    "DIS-POL-POWHEG_MM-NEGNLO-ALL.out",
    "DIS-POL-POWHEG_PP-POSNLO-Z.out",
    "DIS-POL-POWHEG_00-POSNLO-Z.out",
    "DIS-POL-POWHEG_PM-POSNLO-Z.out",
    "DIS-POL-POWHEG_MM-POSNLO-Z.out",
    "DIS-POL-POWHEG_MP-POSNLO-Z.out",
    "DIS-POL-POWHEG_PP-NEGNLO-Z.out",
    "DIS-POL-POWHEG_00-NEGNLO-Z.out",
    "DIS-POL-POWHEG_PM-NEGNLO-Z.out",
    "DIS-POL-POWHEG_MM-NEGNLO-Z.out",
    "DIS-POL-POWHEG_MP-NEGNLO-Z.out",
    "DIS-POL-POWHEG_00-NEGNLO-GAMMA.out",
    "DIS-POL-POWHEG_PP-NEGNLO-GAMMA.out",
    "DIS-POL-POWHEG_PM-NEGNLO-GAMMA.out",
    "DIS-POL-POWHEG_00-POSNLO-GAMMA.out",
    "DIS-POL-POWHEG_PP-POSNLO-GAMMA.out",
    "DIS-POL-POWHEG_PM-POSNLO-GAMMA.out",
    "DIS-POL-LO_PM-GAMMA.out",
    "DIS-POL-LO_00-GAMMA.out",
    "DIS-POL-LO_PP-GAMMA.out",
    "DIS-POL-LO_00-ALL.out",
    "DIS-POL-LO_PP-ALL.out",
    "DIS-POL-LO_PM-ALL.out",
    "DIS-POL-LO_MP-ALL.out",
    "DIS-POL-LO_MM-ALL.out",
    "DIS-POL-LO_00-Z.out",
    "DIS-POL-LO_PP-Z.out",
    "DIS-POL-LO_PM-Z.out",
    "DIS-POL-LO_MP-Z.out",
    "DIS-POL-LO_MM-Z.out",
]

RIVETFOFIXED_DEFAULT_RUNS = [
    "DIS-POL-POWHEG_PP-POSNLO-ALL.out",
    "DIS-POL-POWHEG_00-POSNLO-ALL.out",
    "DIS-POL-POWHEG_PM-POSNLO-ALL.out",
    "DIS-POL-POWHEG_MP-POSNLO-ALL.out",
    "DIS-POL-POWHEG_MM-POSNLO-ALL.out",
    "DIS-POL-POWHEG_PP-POSNLO-Z.out",
    "DIS-POL-POWHEG_00-POSNLO-Z.out",
    "DIS-POL-POWHEG_PM-POSNLO-Z.out",
    "DIS-POL-POWHEG_MM-POSNLO-Z.out",
    "DIS-POL-POWHEG_MP-POSNLO-Z.out",
    "DIS-POL-POWHEG_00-POSNLO-GAMMA.out",
    "DIS-POL-POWHEG_PP-POSNLO-GAMMA.out",
    "DIS-POL-POWHEG_PM-POSNLO-GAMMA.out",
]

CC_DEFAULT_RUNS = [
    "DIS-POL-POWHEG_PP-POSNLO-CC.out",
    "DIS-POL-POWHEG_00-POSNLO-CC.out",
    "DIS-POL-POWHEG_PM-POSNLO-CC.out",
    "DIS-POL-POWHEG_PP-NEGNLO-CC.out",
    "DIS-POL-POWHEG_00-NEGNLO-CC.out",
    "DIS-POL-POWHEG_PM-NEGNLO-CC.out",
    "DIS-POL-POWHEG_MP-POSNLO-CC.out",
    "DIS-POL-POWHEG_MM-POSNLO-CC.out",
    "DIS-POL-POWHEG_MP-NEGNLO-CC.out",
    "DIS-POL-POWHEG_MM-NEGNLO-CC.out",
    "DIS-POL-LO_00-CC.out",
    "DIS-POL-LO_PP-CC.out",
    "DIS-POL-LO_PM-CC.out",
    "DIS-POL-LO_MP-CC.out",
    "DIS-POL-LO_MM-CC.out",
]

SETUP_ORDER = {
    "GAMMA": 0,
    "Z": 1,
    "ALL": 2,
}
DEFAULT_SETUPS = tuple(SETUP_ORDER)
SUPPORTED_SETUPS = DEFAULT_SETUPS + ("CC",)


TOTAL_RE = re.compile(
    r"Total\s+\(from generated events\):\s+(\S+)\s+(\S+)\s+([0-9.]+)\((\d+)\)e([+-]?\d+)"
)
UNIT_RE = re.compile(r"Cross-section\s*\(([^)]+)\)")

LO_NAME_RE = re.compile(
    r"^DIS-POL-LO_(?P<hel>PP|PM|MP|MM|00)-(?P<ew>ALL|GAMMA|Z|CC)"
    r"(?:-(?P<analysis>RIVETFOFIXED|RIVETFO|RIVET))?(?:-(?P<variant>[^.]+))?\.out$"
)
NLO_NAME_RE = re.compile(
    r"^DIS-POL-POWHEG_(?P<hel>PP|PM|MP|MM|00)-(?P<part>POSNLO|NEGNLO)-(?P<ew>ALL|GAMMA|Z|CC)"
    r"(?:-(?P<analysis>RIVETFOFIXED|RIVETFO|RIVET))?(?:-(?P<variant>[^.]+))?\.out$"
)
# Support both plain shard variants like
#   S700560-plain48-z-termdiag-s001
# and prefixed variants like
#   TERMDIAG-SP040-S710000-plain49-gamma-termdiag-power-s001
SEEDED_VARIANT_RE = re.compile(r"(?:^|-)S\d+-(?P<tag>.+)$")


@dataclass(frozen=True)
class Measurement:
    value_pb: float
    error_pb: float


def placeholder_measurement() -> Measurement:
    return Measurement(float("nan"), float("nan"))


def has_reference_measurement(refs: Dict[str, Measurement], stage: str) -> bool:
    meas = refs.get(stage)
    return meas is not None and math.isfinite(meas.value_pb) and math.isfinite(meas.error_pb)


def is_finite_measurement(meas: Optional[Measurement]) -> bool:
    return meas is not None and math.isfinite(meas.value_pb) and math.isfinite(meas.error_pb)


def emit_progress(enabled: bool, message: str) -> None:
    if enabled:
        print(f"[extract-dis] {message}", file=sys.stderr, flush=True)


HELICITY_ORDER = ("PP", "PM", "MP", "MM", "00")
THREE_STATE_TAGS = ("PP", "PM", "00")
FOUR_STATE_TAGS = ("PP", "PM", "MP", "MM")
STAGE_ORDER = ("LO", "POSNLO", "NEGNLO", "NLO")
COMPARISON_ORDER = ("LO_unpol", "NLO_unpol", "LO_pol", "NLO_pol")
INTERFERENCE_COMPARISON_ORDER = ("LO_unpol", "NLO_unpol", "LO_pol", "NLO_pol")

DISPLAY_LABELS = {
    "avg": "(PP+PM)/2",
    "pol": "(PP-PM)/2",
    "avg_minus_00": "(PP+PM)/2 - 00",
    "sigma0": "sigma0 = (PP+PM+MP+MM)/4",
    "sigma_l": "sigma_l = (PP+PM-MP-MM)/4",
    "sigma_p": "sigma_p = (PP-PM+MP-MM)/4",
    "sigma_ll": "sigma_LL = (PP+MM-PM-MP)/4",
    "sigma0_minus_00": "sigma0 - 00",
    "00": "sigma_00",
}

DELTA_OBSERVABLES = {"avg_minus_00", "sigma0_minus_00"}


# POLDIS reference cross sections for the current validation cut Q^2 > 49 GeV^2.
POLDIS_POL_REFS: Dict[str, Dict[str, Measurement]] = {
    "ALL": {
        "LO": Measurement(50.140712, 0.003077),
        "NLO": Measurement(48.288810, 0.003227),
        "NNLO": Measurement(48.158641, 0.003761),
    },
    "GAMMA": {
        "LO": Measurement(45.426932, 0.002804),
        "NLO": Measurement(43.669509, 0.003004),
        "NNLO": Measurement(43.544651, 0.003565),
    },
    "Z": {
        "LO": Measurement(0.054656, 0.000010),
        "NLO": Measurement(0.050150, 0.000010),
        "NNLO": Measurement(0.049715, 0.000010),
    },
    "CC": {
        "LO": Measurement(3.157534, 0.001497),
        "NLO": Measurement(3.035988, 0.001456),
        "NNLO": Measurement(3.021219, 0.001494),
    },
}

POLDIS_UNPOL_REFS: Dict[str, Dict[str, Measurement]] = {
    "ALL": {
        "LO": Measurement(2965.620920, 0.200533),
        "NLO": Measurement(2768.302391, 0.227469),
        "NNLO": Measurement(2689.313720, 0.299187),
    },
    "GAMMA": {
        "LO": Measurement(2952.498947, 0.200487),
        "NLO": Measurement(2756.008931, 0.227363),
        "NNLO": Measurement(2677.170133, 0.299018),
    },
    "Z": {
        "LO": Measurement(1.210471, 0.000119),
        "NLO": Measurement(1.153555, 0.000118),
        "NNLO": Measurement(1.139722, 0.000123),
    },
    "CC": {
        "LO": Measurement(9.431494, 0.003275),
        "NLO": Measurement(9.038495, 0.003235),
        "NNLO": Measurement(8.956922, 0.003349),
    },
}


def clone_reference_map(
    refs: Dict[str, Dict[str, Measurement]],
) -> Dict[str, Dict[str, Measurement]]:
    return {
        setup: {
            order: Measurement(meas.value_pb, meas.error_pb)
            for order, meas in orders.items()
        }
        for setup, orders in refs.items()
    }


def merge_reference_maps(
    default_refs: Dict[str, Dict[str, Measurement]],
    override_refs: Dict[str, Dict[str, Measurement]],
) -> Dict[str, Dict[str, Measurement]]:
    merged = clone_reference_map(default_refs)
    for setup, orders in override_refs.items():
        target = merged.setdefault(setup, {})
        for order, meas in orders.items():
            target[order] = Measurement(meas.value_pb, meas.error_pb)
    return merged


def parse_reference_measurement(payload: object, context: str) -> Measurement:
    if not isinstance(payload, dict):
        raise ValueError(f"Invalid {context}: expected object, got {type(payload).__name__}")
    try:
        return Measurement(float(payload["value_pb"]), float(payload["error_pb"]))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid {context}: expected value_pb/error_pb numbers") from exc


def load_poldis_reference_overrides(
    path: Path,
) -> Tuple[Dict[str, Dict[str, Measurement]], Dict[str, Dict[str, Measurement]]]:
    payload = json.loads(path.read_text())
    raw_setups = payload.get("setups")
    if not isinstance(raw_setups, dict):
        raise ValueError(f"Invalid POLDIS reference JSON {path}: missing top-level 'setups' object")

    unpol_refs: Dict[str, Dict[str, Measurement]] = {}
    pol_refs: Dict[str, Dict[str, Measurement]] = {}
    for setup, setup_payload in raw_setups.items():
        if not isinstance(setup, str) or not isinstance(setup_payload, dict):
            continue
        raw_unpol = setup_payload.get("unpolarized", {})
        raw_pol = setup_payload.get("polarized", {})
        if isinstance(raw_unpol, dict):
            unpol_orders: Dict[str, Measurement] = {}
            for order in ("LO", "NLO", "NNLO"):
                if order in raw_unpol:
                    unpol_orders[order] = parse_reference_measurement(
                        raw_unpol[order], f"{path}:{setup}:unpolarized:{order}"
                    )
            if unpol_orders:
                unpol_refs[setup] = unpol_orders
        if isinstance(raw_pol, dict):
            pol_orders: Dict[str, Measurement] = {}
            for order in ("LO", "NLO", "NNLO"):
                if order in raw_pol:
                    pol_orders[order] = parse_reference_measurement(
                        raw_pol[order], f"{path}:{setup}:polarized:{order}"
                    )
            if pol_orders:
                pol_refs[setup] = pol_orders
    return unpol_refs, pol_refs


def resolve_poldis_reference_maps(
    override_json: Optional[Path],
) -> Tuple[Dict[str, Dict[str, Measurement]], Dict[str, Dict[str, Measurement]], Optional[str]]:
    if override_json is None:
        return clone_reference_map(POLDIS_UNPOL_REFS), clone_reference_map(POLDIS_POL_REFS), None
    unpol_override, pol_override = load_poldis_reference_overrides(override_json)
    return (
        merge_reference_maps(POLDIS_UNPOL_REFS, unpol_override),
        merge_reference_maps(POLDIS_POL_REFS, pol_override),
        str(override_json),
    )


@dataclass
class RunSpec:
    requested_name: str
    path: Optional[Path]
    paths: List[Path]
    exists: bool
    helicity: str
    order: str
    ew: str
    piece: str
    analysis: str
    variant: str
    measurement: Optional[Measurement]
    unit_in_file: Optional[str]


def nb_to_pb(value: float) -> float:
    return value * 1.0e3


def value_with_paren_error(mantissa_text: str, err_digits_text: str, exponent_text: str) -> Tuple[float, float]:
    mantissa = float(mantissa_text)
    exponent = int(exponent_text)
    value = mantissa * (10.0 ** exponent)

    if "." in mantissa_text:
        decimals = len(mantissa_text.split(".", 1)[1])
    else:
        decimals = 0
    error = int(err_digits_text) * (10.0 ** (exponent - decimals))
    return value, error


def parse_measurement(text: str) -> Tuple[Optional[Measurement], Optional[str]]:
    unit_match = UNIT_RE.search(text)
    unit = unit_match.group(1).strip() if unit_match else None

    matches = list(TOTAL_RE.finditer(text))
    if not matches:
        return None, unit

    last = matches[-1]
    value_raw, error_raw = value_with_paren_error(last.group(3), last.group(4), last.group(5))

    if unit is None or unit.lower() == "nb":
        return Measurement(nb_to_pb(value_raw), nb_to_pb(error_raw)), unit or "nb"
    if unit.lower() == "pb":
        return Measurement(value_raw, error_raw), unit

    raise ValueError(f"Unsupported cross-section unit {unit!r}")


def parse_generated_event_count(text: str) -> Optional[int]:
    matches = list(TOTAL_RE.finditer(text))
    if not matches:
        return None
    try:
        return int(matches[-1].group(1))
    except ValueError:
        return None


def apply_analysis_suffix(name: str, analysis_variant: str) -> str:
    if not analysis_variant:
        return name
    suffix = f"-{analysis_variant}.out"
    if name.endswith(suffix):
        return name
    if name.endswith(".out"):
        return name[:-4] + suffix
    return name


def parse_requested_name(name: str) -> Tuple[str, str, str, str, str, str]:
    m = LO_NAME_RE.match(name)
    if m:
        return m.group("hel"), "LO", m.group("ew"), "LO", m.group("analysis") or "", m.group("variant") or ""

    m = NLO_NAME_RE.match(name)
    if m:
        return m.group("hel"), "NLO", m.group("ew"), m.group("part"), m.group("analysis") or "", m.group("variant") or ""

    raise ValueError(f"Unrecognized run filename pattern: {name}")


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


def candidate_sort_key(
    path: Path,
    requested_variant: str,
    actual_variant: str,
    exact_basename: str,
    preferred_tag: str,
) -> Tuple[int, int, int, int, str]:
    exact_name_penalty = 0 if path.name == exact_basename else 1
    preferred_penalty = 0 if variant_matches_tag(actual_variant, preferred_tag) else 1
    exact_variant_penalty = 0 if normalize_variant_tag(actual_variant) == requested_variant else 1
    empty_variant_penalty = 0 if normalize_variant_tag(actual_variant) == "" else 1
    if not preferred_tag:
        preferred_penalty = 0
    return (preferred_penalty, exact_name_penalty, exact_variant_penalty, empty_variant_penalty, path.name)


def find_equivalent_runs(base_dir: Path, requested_name: str, preferred_tag: str) -> List[Path]:
    helicity, order, ew, piece, analysis, variant = parse_requested_name(requested_name)
    candidates: List[Tuple[Tuple[int, int, int, int, str], Path]] = []
    for path in base_dir.rglob("DIS-POL-*.out"):
        try:
            cand_helicity, cand_order, cand_ew, cand_piece, cand_analysis, cand_variant = parse_requested_name(path.name)
        except ValueError:
            continue
        if (cand_helicity, cand_order, cand_ew, cand_piece, cand_analysis) != (helicity, order, ew, piece, analysis):
            continue
        candidates.append((candidate_sort_key(path, variant, cand_variant, requested_name, preferred_tag), path))

    if not candidates:
        return []
    candidates.sort(key=lambda item: item[0])
    paths = [path for _, path in candidates]
    if preferred_tag:
        shard_paths = []
        exact_paths = []
        for path in paths:
            try:
                _, _, _, _, _, cand_variant = parse_requested_name(path.name)
            except ValueError:
                continue
            if normalize_variant_tag(cand_variant) == preferred_tag:
                exact_paths.append(path)
            elif is_shard_variant(cand_variant, preferred_tag):
                shard_paths.append(path)
        if shard_paths:
            return sorted(shard_paths, key=lambda item: item.name)
        if exact_paths:
            return [exact_paths[0]]
    return [paths[0]]


def resolve_runs(base_dir: Path, requested_name: str, preferred_tag: str, strict_tag: bool = False) -> List[Path]:
    if preferred_tag:
        preferred = find_equivalent_runs(base_dir, requested_name, preferred_tag)
        if preferred:
            return preferred
        if strict_tag:
            return []

    direct = base_dir / requested_name
    if direct.exists():
        return [direct]

    matches = sorted(base_dir.rglob(requested_name))
    if matches:
        return [matches[0]]

    return find_equivalent_runs(base_dir, requested_name, preferred_tag)


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
        if not math.isfinite(weight) or weight <= 0.0:
            return [1.0] * count
        weights.append(weight)
    return weights


def combine_measurements(
    measurements: Sequence[Measurement],
    generated_events: Optional[Sequence[Optional[int]]] = None,
) -> Measurement:
    if not measurements:
        raise ValueError("combine_measurements() requires at least one measurement")

    # Campaign shards are independent samples of the same run configuration, so
    # merge them as a mean over shard event counts rather than by inverse-variance.
    weights = shard_mean_weights(generated_events, len(measurements))
    total_weight = sum(weights)
    value = sum(w * m.value_pb for w, m in zip(weights, measurements)) / total_weight
    error = math.sqrt(sum((w * m.error_pb) ** 2 for w, m in zip(weights, measurements))) / total_weight
    return Measurement(value, error)


def load_run(base_dir: Path, requested_name: str, preferred_tag: str, strict_tag: bool = False) -> RunSpec:
    helicity, order, ew, piece, analysis, variant = parse_requested_name(requested_name)
    paths = resolve_runs(base_dir, requested_name, preferred_tag, strict_tag=strict_tag)
    if not paths:
        return RunSpec(
            requested_name=requested_name,
            path=None,
            paths=[],
            exists=False,
            helicity=helicity,
            order=order,
            ew=ew,
            piece=piece,
            analysis=analysis,
            variant=variant,
            measurement=None,
            unit_in_file=None,
        )

    measurements: List[Measurement] = []
    generated_events: List[Optional[int]] = []
    units: List[str] = []
    for path in paths:
        text = path.read_text()
        measurement, unit = parse_measurement(text)
        if measurement is not None:
            measurements.append(measurement)
            generated_events.append(parse_generated_event_count(text))
        if unit is not None:
            units.append(unit)
    measurement = combine_measurements(measurements, generated_events) if measurements else None
    unit = units[0] if units else None
    return RunSpec(
        requested_name=requested_name,
        path=paths[0],
        paths=paths,
        exists=True,
        helicity=helicity,
        order=order,
        ew=ew,
        piece=piece,
        analysis=analysis,
        variant=variant,
        measurement=measurement,
        unit_in_file=unit,
    )


def add(a: Measurement, b: Measurement) -> Measurement:
    return Measurement(a.value_pb + b.value_pb, math.sqrt(a.error_pb ** 2 + b.error_pb ** 2))


def sub(a: Measurement, b: Measurement) -> Measurement:
    return Measurement(a.value_pb - b.value_pb, math.sqrt(a.error_pb ** 2 + b.error_pb ** 2))


def half_diff(a: Measurement, b: Measurement) -> Measurement:
    diff = sub(a, b)
    return Measurement(diff.value_pb / 2.0, diff.error_pb / 2.0)


def half_sum(a: Measurement, b: Measurement) -> Measurement:
    total = add(a, b)
    return Measurement(total.value_pb / 2.0, total.error_pb / 2.0)


def linear_combo(terms: Sequence[Tuple[float, Measurement]]) -> Measurement:
    value = sum(coeff * meas.value_pb for coeff, meas in terms)
    error = math.sqrt(sum((coeff * meas.error_pb) ** 2 for coeff, meas in terms))
    return Measurement(value, error)


def ratio(num: Measurement, den: Measurement) -> Optional[Measurement]:
    if abs(den.value_pb) < 1.0e-30:
        return None
    value = num.value_pb / den.value_pb
    # Linear error propagation for r = n/d:
    #   sigma_r^2 = (sigma_n/d)^2 + (n*sigma_d/d^2)^2
    # This stays well-defined even when the numerator is close to zero.
    error = math.sqrt(
        (num.error_pb / den.value_pb) ** 2
        + ((num.value_pb * den.error_pb) / (den.value_pb ** 2)) ** 2
    )
    return Measurement(value, error)


def decimals_from_error(error: float, sigfigs: int = 3) -> int:
    if error <= 0.0:
        return 6
    exponent = math.floor(math.log10(error))
    return max(0, sigfigs - 1 - exponent)


def fmt_xsec(m: Optional[Measurement]) -> str:
    if m is None:
        return "missing"
    if not is_finite_measurement(m):
        return "pending"
    decimals = decimals_from_error(abs(m.error_pb))
    return f"{m.value_pb:.{decimals}f} +- {m.error_pb:.{decimals}f} pb"


def fmt_delta(m: Optional[Measurement]) -> str:
    if m is None:
        return "missing"
    if not is_finite_measurement(m):
        return "pending"
    decimals = decimals_from_error(abs(m.error_pb))
    return f"{m.value_pb:+.{decimals}f} +- {m.error_pb:.{decimals}f} pb"


def fmt_ratio(m: Optional[Measurement]) -> str:
    if m is None:
        return "missing"
    if not is_finite_measurement(m):
        return "pending"
    return f"{m.value_pb:.6f} +- {m.error_pb:.6f}"


def fmt_pull(data: Measurement, ref: Measurement) -> str:
    diff = sub(data, ref)
    if diff.error_pb <= 0.0:
        return "n/a"
    return f"{abs(diff.value_pb) / diff.error_pb:.2f} sigma"


def pull_sigma(data: Measurement, ref: Measurement) -> Optional[float]:
    diff = sub(data, ref)
    if diff.error_pb <= 0.0:
        return None
    return abs(diff.value_pb) / diff.error_pb


def has_tags(runs: Dict[str, Measurement], tags: Sequence[str]) -> bool:
    return all(tag in runs for tag in tags)


def setup_sort_key(setup: str) -> Tuple[int, str]:
    return (SETUP_ORDER.get(setup, 99), setup)


def format_file_timestamp(path: Path) -> str:
    stamp = dt.datetime.fromtimestamp(path.stat().st_mtime).astimezone()
    return stamp.strftime("%Y-%m-%d %H:%M:%S %Z")


def measurement_from_dict(data: Dict[str, float]) -> Measurement:
    return Measurement(data["value_pb"], data["error_pb"])


def display_label(observable: str) -> str:
    return DISPLAY_LABELS.get(observable, observable)


def format_stage_observable(observable: str, measurement: Measurement) -> str:
    if observable in DELTA_OBSERVABLES:
        return fmt_delta(measurement)
    return fmt_xsec(measurement)


def render_table(
    headers: Sequence[str],
    rows: Sequence[Sequence[str]],
    aligns: Optional[Sequence[str]] = None,
) -> List[str]:
    if not rows:
        return []

    string_rows = [[str(cell) for cell in row] for row in rows]
    alignments = list(aligns) if aligns is not None else ["l"] * len(headers)

    if PrettyTable is not None:
        table = PrettyTable()
        table.field_names = list(headers)
        for header, align in zip(headers, alignments):
            table.align[header] = "r" if align == "r" else "l"
        for row in string_rows:
            table.add_row(row)
        return table.get_string().splitlines()

    widths = []
    for idx, header in enumerate(headers):
        cell_width = max((len(row[idx]) for row in string_rows), default=0)
        widths.append(max(len(header), cell_width))

    def format_row(row: Sequence[str], is_header: bool = False) -> str:
        cells: List[str] = []
        for idx, cell in enumerate(row):
            align = alignments[idx] if idx < len(alignments) else "l"
            width = widths[idx]
            if is_header:
                cells.append(cell.ljust(width))
            elif align == "r":
                cells.append(cell.rjust(width))
            else:
                cells.append(cell.ljust(width))
        return "| " + " | ".join(cells) + " |"

    separator = "+-" + "-+-".join("-" * width for width in widths) + "-+"
    lines = [separator, format_row(headers, is_header=True), separator]
    for row in string_rows:
        lines.append(format_row(row))
    lines.append(separator)
    return lines


def comparison_title(comp_name: str, comp: Dict[str, object]) -> str:
    if "title" in comp:
        return str(comp["title"])
    label = comp.get("label", "")
    if comp_name == "LO_unpol":
        return f"LO unpolarized ({label})"
    if comp_name == "NLO_unpol":
        return f"NLO unpolarized ({label})"
    if comp_name == "LO_pol":
        return f"LO polarized ({label})"
    if comp_name == "NLO_pol":
        return f"NLO polarized ({label})"
    return str(label or comp_name)


def comparison_cells(comp_name: str, comp: Dict[str, object]) -> List[str]:
    herwig = measurement_from_dict(comp["herwig"])
    poldis = measurement_from_dict(comp["poldis"])
    delta = measurement_from_dict(comp["delta"])
    ratio_payload = comp.get("ratio")
    ratio_text = "missing"
    if ratio_payload is not None:
        ratio_text = fmt_ratio(measurement_from_dict(ratio_payload))

    if comp_name.startswith("K_"):
        herwig_text = fmt_ratio(herwig)
        poldis_text = fmt_ratio(poldis)
        delta_text = f"{delta.value_pb:+.6f} +- {delta.error_pb:.6f}"
    else:
        herwig_text = fmt_xsec(herwig)
        poldis_text = fmt_xsec(poldis)
        delta_text = fmt_delta(delta)

    pull = comp.get("pull_sigma")
    pull_text = "n/a" if pull is None else f"{pull:.2f} sigma"
    return [comparison_title(comp_name, comp), herwig_text, poldis_text, delta_text, ratio_text, pull_text]


def select_pol_observable(
    setup: str, stage: Dict[str, object]
) -> Optional[Tuple[str, str, Measurement]]:
    derived = stage.get("derived", {})
    if setup in ("ALL", "Z", "CC"):
        if "sigma_ll" not in derived:
            return None
        meas = derived["sigma_ll"]
        return (
            "sigma_ll",
            "sigma_LL = (PP+MM-PM-MP)/4",
            Measurement(meas["value_pb"], meas["error_pb"]),
        )
    if "pol" in derived:
        meas = derived["pol"]
        return (
            "pol",
            "(PP-PM)/2",
            Measurement(meas["value_pb"], meas["error_pb"]),
        )
    return None


def select_unpol_observable(
    setup: str, stage: Dict[str, object]
) -> Optional[Tuple[str, str, Measurement]]:
    derived = stage.get("derived", {})
    helicity = stage.get("helicity", {})
    if setup in ("ALL", "Z", "CC") and "sigma0" in derived:
        meas = derived["sigma0"]
        return (
            "sigma0",
            "sigma0 = (PP+PM+MP+MM)/4",
            Measurement(meas["value_pb"], meas["error_pb"]),
        )
    if "00" in helicity:
        meas = helicity["00"]
        return (
            "00",
            "sigma_00",
            Measurement(meas["value_pb"], meas["error_pb"]),
        )
    if "sigma0" in derived:
        meas = derived["sigma0"]
        return (
            "sigma0",
            "sigma0 = (PP+PM+MP+MM)/4",
            Measurement(meas["value_pb"], meas["error_pb"]),
        )
    return None


def summarise_stage(label: str, runs: Dict[str, Measurement]) -> List[str]:
    lines: List[str] = [f"  {label}:"]
    for hel in HELICITY_ORDER:
        if hel in runs:
            lines.append(f"    {hel}: {fmt_xsec(runs.get(hel))}")

    any_derived = False
    if has_tags(runs, ("PP", "PM")):
        pol = half_diff(runs["PP"], runs["PM"])
        avg = half_sum(runs["PP"], runs["PM"])
        lines.append(f"    (PP-PM)/2: {fmt_xsec(pol)}")
        lines.append(f"    (PP+PM)/2: {fmt_xsec(avg)}")
        if "00" in runs:
            lines.append(f"    (PP+PM)/2 - 00: {fmt_delta(sub(avg, runs['00']))}")
        any_derived = True

    if has_tags(runs, FOUR_STATE_TAGS):
        sigma0 = linear_combo([(0.25, runs["PP"]), (0.25, runs["PM"]), (0.25, runs["MP"]), (0.25, runs["MM"])])
        sigma_l = linear_combo([(0.25, runs["PP"]), (0.25, runs["PM"]), (-0.25, runs["MP"]), (-0.25, runs["MM"])])
        sigma_p = linear_combo([(0.25, runs["PP"]), (-0.25, runs["PM"]), (0.25, runs["MP"]), (-0.25, runs["MM"])])
        sigma_ll = linear_combo([(0.25, runs["PP"]), (-0.25, runs["PM"]), (-0.25, runs["MP"]), (0.25, runs["MM"])])
        lines.append(f"    sigma0 = (PP+PM+MP+MM)/4: {fmt_xsec(sigma0)}")
        lines.append(f"    sigma_l = (PP+PM-MP-MM)/4: {fmt_xsec(sigma_l)}")
        lines.append(f"    sigma_p = (PP-PM+MP-MM)/4: {fmt_xsec(sigma_p)}")
        lines.append(f"    sigma_LL = (PP+MM-PM-MP)/4: {fmt_xsec(sigma_ll)}")
        if "00" in runs:
            lines.append(f"    sigma0 - 00: {fmt_delta(sub(sigma0, runs['00']))}")
        any_derived = True

    if not any_derived:
        lines.append("    Derived combinations: incomplete helicity set")
    return lines


def compare_lines(label: str, data: Measurement, ref: Measurement, unit: str = "pb") -> List[str]:
    diff = sub(data, ref)
    rat = ratio(data, ref)
    if unit == "pb":
        data_text = fmt_xsec(data)
        ref_text = fmt_xsec(ref)
        diff_text = fmt_delta(diff)
    else:
        data_text = fmt_ratio(data)
        ref_text = fmt_ratio(ref)
        diff_text = f"{diff.value_pb:+.6f} +- {diff.error_pb:.6f}"
    return [
        f"    {label}:",
        f"      Herwig = {data_text}",
        f"      POLDIS = {ref_text}",
        f"      Delta  = {diff_text}",
        f"      Ratio  = {fmt_ratio(rat)}",
        f"      Pull   = {fmt_pull(data, ref)}",
    ]


def measurement_to_dict(m: Measurement) -> Dict[str, float]:
    return {
        "value_pb": m.value_pb,
        "error_pb": m.error_pb,
    }


def stage_payload(runs: Dict[str, Measurement]) -> Dict[str, object]:
    payload: Dict[str, object] = {
        "complete_helicity_set": has_tags(runs, THREE_STATE_TAGS),
        "complete_four_helicity_set": has_tags(runs, FOUR_STATE_TAGS),
        "helicity": {hel: measurement_to_dict(runs[hel]) for hel in HELICITY_ORDER if hel in runs},
    }
    derived: Dict[str, object] = {}
    if has_tags(runs, ("PP", "PM")):
        pol = half_diff(runs["PP"], runs["PM"])
        avg = half_sum(runs["PP"], runs["PM"])
        derived["pol"] = measurement_to_dict(pol)
        derived["avg"] = measurement_to_dict(avg)
        if "00" in runs:
            derived["avg_minus_00"] = measurement_to_dict(sub(avg, runs["00"]))
    if has_tags(runs, FOUR_STATE_TAGS):
        sigma0 = linear_combo([(0.25, runs["PP"]), (0.25, runs["PM"]), (0.25, runs["MP"]), (0.25, runs["MM"])])
        sigma_l = linear_combo([(0.25, runs["PP"]), (0.25, runs["PM"]), (-0.25, runs["MP"]), (-0.25, runs["MM"])])
        sigma_p = linear_combo([(0.25, runs["PP"]), (-0.25, runs["PM"]), (0.25, runs["MP"]), (-0.25, runs["MM"])])
        sigma_ll = linear_combo([(0.25, runs["PP"]), (-0.25, runs["PM"]), (-0.25, runs["MP"]), (0.25, runs["MM"])])
        derived["sigma0"] = measurement_to_dict(sigma0)
        derived["sigma_l"] = measurement_to_dict(sigma_l)
        derived["sigma_p"] = measurement_to_dict(sigma_p)
        derived["sigma_ll"] = measurement_to_dict(sigma_ll)
        if "00" in runs:
            derived["sigma0_minus_00"] = measurement_to_dict(sub(sigma0, runs["00"]))
    if derived:
        payload["derived"] = derived
    return payload


def comparison_payload(data: Measurement, ref: Measurement) -> Dict[str, object]:
    diff = sub(data, ref)
    rat = ratio(data, ref)
    return {
        "herwig": measurement_to_dict(data),
        "poldis": measurement_to_dict(ref),
        "delta": measurement_to_dict(diff),
        "ratio": measurement_to_dict(rat) if rat is not None else None,
        "pull_sigma": pull_sigma(data, ref),
    }


def group_runs(specs: Iterable[RunSpec]) -> Dict[str, Dict[str, Dict[str, RunSpec]]]:
    grouped: Dict[str, Dict[str, Dict[str, RunSpec]]] = {}
    for spec in specs:
        grouped.setdefault(spec.ew, {}).setdefault(spec.piece, {})[spec.helicity] = spec
    return grouped


def existing_measurements(specs: Dict[str, RunSpec]) -> Dict[str, Measurement]:
    out: Dict[str, Measurement] = {}
    for hel, spec in specs.items():
        if spec.exists and spec.measurement is not None:
            out[hel] = spec.measurement
    return out


def selected_stage_observable(
    setups_payload: Dict[str, object],
    setup: str,
    stage_name: str,
    selector,
) -> Optional[Tuple[str, str, Measurement]]:
    setup_payload = setups_payload.get(setup)
    if setup_payload is None:
        return None
    stage_payload = setup_payload.get("stages", {}).get(stage_name)
    if stage_payload is None:
        return None
    return selector(setup, stage_payload)


def build_interference_comparisons(
    setups_payload: Dict[str, object],
    unpol_refs_map: Dict[str, Dict[str, Measurement]],
    pol_refs_map: Dict[str, Dict[str, Measurement]],
) -> Dict[str, object]:
    required_setups = ("ALL", "GAMMA", "Z")
    if not all(setup in setups_payload for setup in required_setups):
        return {}

    comparisons: Dict[str, object] = {}
    configs = (
        ("LO_unpol", "LO", select_unpol_observable, unpol_refs_map, "LO unpolarized interference"),
        ("NLO_unpol", "NLO", select_unpol_observable, unpol_refs_map, "NLO unpolarized interference"),
        ("LO_pol", "LO", select_pol_observable, pol_refs_map, "LO polarized interference"),
        ("NLO_pol", "NLO", select_pol_observable, pol_refs_map, "NLO polarized interference"),
    )

    for comp_name, stage_name, selector, refs, title in configs:
        all_choice = selected_stage_observable(setups_payload, "ALL", stage_name, selector)
        gamma_choice = selected_stage_observable(setups_payload, "GAMMA", stage_name, selector)
        z_choice = selected_stage_observable(setups_payload, "Z", stage_name, selector)
        if all_choice is None or gamma_choice is None or z_choice is None:
            continue

        herwig = sub(sub(all_choice[2], gamma_choice[2]), z_choice[2])
        poldis = sub(sub(refs["ALL"][stage_name], refs["GAMMA"][stage_name]), refs["Z"][stage_name])
        comp = comparison_payload(herwig, poldis)
        comp["observable"] = all_choice[0]
        comp["label"] = f"{all_choice[1]} interference = ALL - GAMMA - Z"
        comp["title"] = title
        comparisons[comp_name] = comp

    return comparisons


def build_summary(
    specs: Sequence[RunSpec],
    unpol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    pol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    ref_source: Optional[str] = None,
) -> Dict[str, object]:
    if unpol_refs_map is None:
        unpol_refs_map = clone_reference_map(POLDIS_UNPOL_REFS)
    if pol_refs_map is None:
        pol_refs_map = clone_reference_map(POLDIS_POL_REFS)
    grouped = group_runs(specs)
    summary: Dict[str, object] = {
        "missing_files": [spec.requested_name for spec in specs if not spec.exists],
        "unparsed_files": [
            {
                "requested_name": spec.requested_name,
                "path": str(spec.path) if spec.path is not None else None,
                "paths": [str(path) for path in spec.paths],
            }
            for spec in specs
            if spec.exists and spec.measurement is None
        ],
        "setups": {},
        "interference": {"comparisons": {}},
        "poldis_reference_source": ref_source,
    }

    setups_payload: Dict[str, object] = {}
    summary["setups"] = setups_payload

    for setup in sorted(grouped, key=setup_sort_key):
        setup_payload: Dict[str, object] = {
            "stages": {},
            "poldis_unpolarized_refs": {
                order: measurement_to_dict(meas)
                for order, meas in unpol_refs_map.get(setup, {}).items()
            },
            "poldis_polarized_refs": {
                order: measurement_to_dict(meas)
                for order, meas in pol_refs_map.get(setup, {}).items()
            },
            "comparisons": {},
        }
        setups_payload[setup] = setup_payload
        pieces = grouped[setup]
        stage_map: Dict[str, object] = setup_payload["stages"]  # type: ignore[assignment]

        lo_runs = existing_measurements(pieces.get("LO", {}))
        if lo_runs:
            stage_map["LO"] = stage_payload(lo_runs)

        pos_runs = existing_measurements(pieces.get("POSNLO", {}))
        neg_runs = existing_measurements(pieces.get("NEGNLO", {}))
        if pos_runs:
            stage_map["POSNLO"] = stage_payload(pos_runs)

        if neg_runs:
            stage_map["NEGNLO"] = stage_payload(neg_runs)

        common_nlo_hels = [hel for hel in HELICITY_ORDER if hel in pos_runs and hel in neg_runs]
        if common_nlo_hels:
            nlo_runs = {hel: sub(pos_runs[hel], neg_runs[hel]) for hel in common_nlo_hels}
            stage_map["NLO"] = stage_payload(nlo_runs)

        pol_refs = pol_refs_map.get(setup)
        if pol_refs:
            comparisons: Dict[str, object] = setup_payload["comparisons"]  # type: ignore[assignment]
            if "LO" in stage_map and has_reference_measurement(pol_refs, "LO"):
                lo_choice = select_pol_observable(setup, stage_map["LO"])  # type: ignore[arg-type]
                if lo_choice is not None:
                    comp = comparison_payload(lo_choice[2], pol_refs["LO"])
                    comp["observable"] = lo_choice[0]
                    comp["label"] = lo_choice[1]
                    comparisons["LO_pol"] = comp
            if "NLO" in stage_map and has_reference_measurement(pol_refs, "NLO"):
                nlo_choice = select_pol_observable(setup, stage_map["NLO"])  # type: ignore[arg-type]
                if nlo_choice is not None:
                    comp = comparison_payload(nlo_choice[2], pol_refs["NLO"])
                    comp["observable"] = nlo_choice[0]
                    comp["label"] = nlo_choice[1]
                    comparisons["NLO_pol"] = comp

        unpol_refs = unpol_refs_map.get(setup)
        if unpol_refs:
            comparisons = setup_payload["comparisons"]  # type: ignore[assignment]
            if "LO" in stage_map and has_reference_measurement(unpol_refs, "LO"):
                lo_choice = select_unpol_observable(setup, stage_map["LO"])  # type: ignore[arg-type]
                if lo_choice is not None:
                    comp = comparison_payload(lo_choice[2], unpol_refs["LO"])
                    comp["observable"] = lo_choice[0]
                    comp["label"] = lo_choice[1]
                    comparisons["LO_unpol"] = comp
            if "NLO" in stage_map and has_reference_measurement(unpol_refs, "NLO"):
                nlo_choice = select_unpol_observable(setup, stage_map["NLO"])  # type: ignore[arg-type]
                if nlo_choice is not None:
                    comp = comparison_payload(nlo_choice[2], unpol_refs["NLO"])
                    comp["observable"] = nlo_choice[0]
                    comp["label"] = nlo_choice[1]
                    comparisons["NLO_unpol"] = comp

    summary["interference"]["comparisons"] = build_interference_comparisons(
        setups_payload,
        unpol_refs_map,
        pol_refs_map,
    )

    return summary


def build_report(
    specs: Sequence[RunSpec],
    unpol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    pol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    ref_source: Optional[str] = None,
) -> str:
    summary = build_summary(
        specs,
        unpol_refs_map=unpol_refs_map,
        pol_refs_map=pol_refs_map,
        ref_source=ref_source,
    )
    lines: List[str] = []

    ref_source_text = summary.get("poldis_reference_source")
    if ref_source_text:
        lines.append(f"POLDIS reference overrides: {ref_source_text}")
        lines.append("")

    used_rows = []
    for spec in specs:
        if not spec.exists or spec.measurement is None or spec.path is None:
            continue
        if len(spec.paths) == 1:
            resolved_name = spec.path.name
            modified = format_file_timestamp(spec.path)
            resolved_path = str(spec.path)
        else:
            resolved_name = f"{len(spec.paths)} shards"
            modified = max(format_file_timestamp(path) for path in spec.paths)
            resolved_path = "; ".join(str(path) for path in spec.paths[:3])
            if len(spec.paths) > 3:
                resolved_path += "; ..."
        used_rows.append(
            [
                spec.requested_name,
                resolved_name,
                modified,
                resolved_path,
            ]
        )
    if used_rows:
        lines.append("Resolved input .out files used:")
        lines.extend(
            render_table(
                ["Requested slot", "Resolved file", "Modified", "Path"],
                used_rows,
                aligns=("l", "l", "l", "l"),
            )
        )
        lines.append("")

    missing = summary["missing_files"]
    if missing:
        lines.append("Missing requested files:")
        lines.extend(render_table(["Requested file"], [[name] for name in missing]))
        lines.append("")

    unparsed = summary["unparsed_files"]
    if unparsed:
        lines.append("Files found but no cross section could be extracted:")
        rows = [[item["requested_name"], item["path"] or ""] for item in unparsed]
        lines.extend(render_table(["Requested file", "Resolved path"], rows))
        lines.append("")

    stage_summaries = summary["setups"]
    for setup in sorted(stage_summaries, key=setup_sort_key):
        lines.append("=" * 88)
        lines.append(f"Setup: {setup}")
        lines.append("=" * 88)

        setup_summary = stage_summaries[setup]
        stages = setup_summary.get("stages", {})

        stage_rows: List[List[str]] = []
        for stage in STAGE_ORDER:
            if stage not in stages:
                continue
            stage_payload = stages[stage]
            for hel in HELICITY_ORDER:
                if hel in stage_payload.get("helicity", {}):
                    stage_rows.append([stage, hel, fmt_xsec(measurement_from_dict(stage_payload["helicity"][hel]))])
            for observable in (
                "avg",
                "pol",
                "avg_minus_00",
                "sigma0",
                "sigma_l",
                "sigma_p",
                "sigma_ll",
                "sigma0_minus_00",
            ):
                if observable in stage_payload.get("derived", {}):
                    stage_rows.append(
                        [
                            stage,
                            display_label(observable),
                            format_stage_observable(
                                observable, measurement_from_dict(stage_payload["derived"][observable])
                            ),
                        ]
                    )
        if stage_rows:
            lines.append("Herwig stages (NLO = POSNLO - NEGNLO):")
            lines.extend(render_table(["Stage", "Observable", "Value"], stage_rows, aligns=("l", "l", "r")))
        else:
            lines.append("Herwig stages: no files found")

        ref_rows: List[List[str]] = []
        unpol_refs = POLDIS_UNPOL_REFS.get(setup, {})
        pol_refs = POLDIS_POL_REFS.get(setup, {})
        for order in ("LO", "NLO", "NNLO"):
            if order not in unpol_refs and order not in pol_refs:
                continue
            ref_rows.append(
                [
                    order,
                    fmt_xsec(unpol_refs[order]) if order in unpol_refs else "missing",
                    fmt_xsec(pol_refs[order]) if order in pol_refs else "missing",
                ]
            )
        if ref_rows:
            lines.append("POLDIS references:")
            lines.extend(
                render_table(
                    ["Order", "Unpolarized", "Polarized"],
                    ref_rows,
                    aligns=("l", "r", "r"),
                )
            )

        comparison_rows: List[List[str]] = []
        for comp_name in COMPARISON_ORDER:
            comp = setup_summary.get("comparisons", {}).get(comp_name)
            if comp is None:
                continue
            comparison_rows.append(comparison_cells(comp_name, comp))
        if comparison_rows:
            lines.append("Herwig vs POLDIS:")
            lines.extend(
                render_table(
                    ["Comparison", "Herwig", "POLDIS", "Delta", "Ratio", "Pull"],
                    comparison_rows,
                    aligns=("l", "r", "r", "r", "r", "r"),
                )
            )

        lines.append("")

    interference = summary.get("interference", {}).get("comparisons", {})
    if interference:
        lines.append("=" * 88)
        lines.append("Interference: ALL - GAMMA - Z")
        lines.append("=" * 88)
        interference_rows = []
        for comp_name in INTERFERENCE_COMPARISON_ORDER:
            comp = interference.get(comp_name)
            if comp is None:
                continue
            interference_rows.append(comparison_cells(comp_name, comp))
        if interference_rows:
            lines.extend(
                render_table(
                    ["Comparison", "Herwig", "POLDIS", "Delta", "Ratio", "Pull"],
                    interference_rows,
                    aligns=("l", "r", "r", "r", "r", "r"),
                )
            )
            lines.append("")

    return "\n".join(lines).rstrip() + "\n"


def write_json(
    path: Path,
    specs: Sequence[RunSpec],
    unpol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    pol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    ref_source: Optional[str] = None,
) -> None:
    payload = build_summary(
        specs,
        unpol_refs_map=unpol_refs_map,
        pol_refs_map=pol_refs_map,
        ref_source=ref_source,
    )
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def csv_rows_from_summary(summary: Dict[str, object]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    setups = summary["setups"]
    for setup in sorted(setups, key=setup_sort_key):
        setup_summary = setups[setup]
        for stage_name, stage_payload in setup_summary.get("stages", {}).items():
            for hel, meas in stage_payload.get("helicity", {}).items():
                rows.append(
                    {
                        "row_type": "herwig_stage",
                        "setup": setup,
                        "stage": stage_name,
                        "observable": hel,
                        "value_pb": meas["value_pb"],
                        "error_pb": meas["error_pb"],
                    }
                )
            for obs, meas in stage_payload.get("derived", {}).items():
                rows.append(
                    {
                        "row_type": "herwig_stage",
                        "setup": setup,
                        "stage": stage_name,
                        "observable": obs,
                        "value_pb": meas["value_pb"],
                        "error_pb": meas["error_pb"],
                    }
                )

        for ref_name, meas in setup_summary.get("poldis_unpolarized_refs", {}).items():
            rows.append(
                {
                    "row_type": "poldis_reference",
                    "setup": setup,
                    "stage": ref_name,
                    "observable": "unpol",
                    "value_pb": meas["value_pb"],
                    "error_pb": meas["error_pb"],
                }
            )

        for ref_name, meas in setup_summary.get("poldis_polarized_refs", {}).items():
            rows.append(
                {
                    "row_type": "poldis_reference",
                    "setup": setup,
                    "stage": ref_name,
                    "observable": "pol",
                    "value_pb": meas["value_pb"],
                    "error_pb": meas["error_pb"],
                }
            )

        for comp_name, comp in setup_summary.get("comparisons", {}).items():
            row = {
                "row_type": "comparison",
                "setup": setup,
                "stage": comp_name,
                "observable": comp.get("observable", "pol"),
                "label": comp.get("label"),
                "herwig_value": comp["herwig"]["value_pb"],
                "herwig_error": comp["herwig"]["error_pb"],
                "poldis_value": comp["poldis"]["value_pb"],
                "poldis_error": comp["poldis"]["error_pb"],
                "delta_value": comp["delta"]["value_pb"],
                "delta_error": comp["delta"]["error_pb"],
                "pull_sigma": comp["pull_sigma"],
            }
            if comp["ratio"] is not None:
                row["ratio_value"] = comp["ratio"]["value_pb"]
                row["ratio_error"] = comp["ratio"]["error_pb"]
            rows.append(row)

    for comp_name, comp in summary.get("interference", {}).get("comparisons", {}).items():
        row = {
            "row_type": "interference_comparison",
            "setup": "INTERFERENCE",
            "stage": comp_name,
            "observable": comp.get("observable", "pol"),
            "label": comp.get("label"),
            "herwig_value": comp["herwig"]["value_pb"],
            "herwig_error": comp["herwig"]["error_pb"],
            "poldis_value": comp["poldis"]["value_pb"],
            "poldis_error": comp["poldis"]["error_pb"],
            "delta_value": comp["delta"]["value_pb"],
            "delta_error": comp["delta"]["error_pb"],
            "pull_sigma": comp["pull_sigma"],
        }
        if comp["ratio"] is not None:
            row["ratio_value"] = comp["ratio"]["value_pb"]
            row["ratio_error"] = comp["ratio"]["error_pb"]
        rows.append(row)

    return rows


def write_csv(
    path: Path,
    specs: Sequence[RunSpec],
    unpol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    pol_refs_map: Optional[Dict[str, Dict[str, Measurement]]] = None,
    ref_source: Optional[str] = None,
) -> None:
    rows = csv_rows_from_summary(
        build_summary(
            specs,
            unpol_refs_map=unpol_refs_map,
            pol_refs_map=pol_refs_map,
            ref_source=ref_source,
        )
    )
    fieldnames = [
        "row_type",
        "setup",
        "stage",
        "observable",
        "label",
        "value_pb",
        "error_pb",
        "herwig_value",
        "herwig_error",
        "poldis_value",
        "poldis_error",
        "delta_value",
        "delta_error",
        "ratio_value",
        "ratio_error",
        "pull_sigma",
    ]
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def normalize_requested_setups(setups: Sequence[str]) -> List[str]:
    if not setups:
        return list(DEFAULT_SETUPS)
    seen = set()
    ordered: List[str] = []
    for setup in setups:
        if setup not in seen:
            ordered.append(setup)
            seen.add(setup)
    return ordered


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "runs",
        nargs="*",
        default=DEFAULT_RUNS,
        help="Run filenames to parse. Defaults to the full standard DIS run list.",
    )
    parser.add_argument(
        "--base-dir",
        default=str(Path(__file__).resolve().parent),
        help="Directory to search for the requested files. Recursive search is used if needed.",
    )
    parser.add_argument(
        "-t",
        "--tag",
        default="",
        help="Prefer runs with this variant tag (for example: testing2).",
    )
    parser.add_argument(
        "--setup",
        action="append",
        choices=SUPPORTED_SETUPS,
        help="Restrict extraction to one or more electroweak setups.",
    )
    parser.add_argument(
        "--rivet",
        action="store_true",
        help="Resolve and report only -RIVET output files.",
    )
    parser.add_argument(
        "--rivetfo",
        action="store_true",
        help="Resolve and report only -RIVETFO output files.",
    )
    parser.add_argument(
        "--rivetfofixed",
        action="store_true",
        help="Resolve and report only -RIVETFOFIXED output files.",
    )
    parser.add_argument(
        "--json-out",
        help="Optional path to write the structured summary as JSON.",
    )
    parser.add_argument(
        "--csv-out",
        help="Optional path to write a flattened summary table as CSV.",
    )
    parser.add_argument(
        "--poldis-refs-json",
        help="Optional JSON file overriding the built-in POLDIS reference totals per setup/order.",
    )
    parser.add_argument(
        "--strict-tag",
        action="store_true",
        help="When -t/--tag is used, do not fall back to untagged or differently tagged files.",
    )
    progress_group = parser.add_mutually_exclusive_group()
    progress_group.add_argument(
        "--progress",
        dest="progress",
        action="store_true",
        help="Print progress messages to stderr while resolving and combining runs.",
    )
    progress_group.add_argument(
        "--no-progress",
        dest="progress",
        action="store_false",
        help="Suppress progress messages.",
    )
    parser.set_defaults(progress=None)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    if int(args.rivet) + int(args.rivetfo) + int(args.rivetfofixed) > 1:
        raise SystemExit("Use at most one of --rivet, --rivetfo, and --rivetfofixed.")
    base_dir = Path(args.base_dir).resolve()
    progress_enabled = sys.stderr.isatty() if args.progress is None else args.progress
    requested_analysis = (
        "RIVETFOFIXED" if args.rivetfofixed else
        ("RIVETFO" if args.rivetfo else ("RIVET" if args.rivet else ""))
    )
    override_path = Path(args.poldis_refs_json).resolve() if args.poldis_refs_json else None
    unpol_refs_map, pol_refs_map, ref_source = resolve_poldis_reference_maps(override_path)
    selected_setups = set(normalize_requested_setups(args.setup or ()))
    requested_base_runs = (
        RIVETFOFIXED_DEFAULT_RUNS
        if requested_analysis == "RIVETFOFIXED" and args.runs == DEFAULT_RUNS
        else args.runs
    )
    if args.runs == DEFAULT_RUNS and requested_analysis != "RIVETFOFIXED" and "CC" in selected_setups:
        requested_base_runs = list(requested_base_runs) + CC_DEFAULT_RUNS
    requested_runs = [apply_analysis_suffix(name, requested_analysis) for name in requested_base_runs]
    if selected_setups:
        requested_runs = [
            name for name in requested_runs
            if parse_requested_name(name)[2] in selected_setups
        ]
    emit_progress(progress_enabled, f"Resolving {len(requested_runs)} requested slots under {base_dir}")
    specs: List[RunSpec] = []
    for index, name in enumerate(requested_runs, start=1):
        emit_progress(progress_enabled, f"[{index}/{len(requested_runs)}] Resolving {name}")
        spec = load_run(base_dir, name, args.tag, strict_tag=args.strict_tag)
        specs.append(spec)
        if not spec.exists:
            emit_progress(progress_enabled, f"[{index}/{len(requested_runs)}] {name}: missing")
        elif spec.measurement is None:
            emit_progress(
                progress_enabled,
                f"[{index}/{len(requested_runs)}] {name}: resolved {len(spec.paths)} file(s), no cross section parsed",
            )
        else:
            emit_progress(
                progress_enabled,
                f"[{index}/{len(requested_runs)}] {name}: resolved {len(spec.paths)} file(s), {fmt_xsec(spec.measurement)}",
            )
    emit_progress(progress_enabled, "Building report")
    sys.stdout.write(
        build_report(
            specs,
            unpol_refs_map=unpol_refs_map,
            pol_refs_map=pol_refs_map,
            ref_source=ref_source,
        )
    )
    if args.json_out:
        json_path = Path(args.json_out).resolve()
        json_path.parent.mkdir(parents=True, exist_ok=True)
        emit_progress(progress_enabled, f"Writing JSON: {json_path}")
        write_json(
            json_path,
            specs,
            unpol_refs_map=unpol_refs_map,
            pol_refs_map=pol_refs_map,
            ref_source=ref_source,
        )
        sys.stdout.write(f"Wrote JSON: {json_path}\n")
    if args.csv_out:
        csv_path = Path(args.csv_out).resolve()
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        emit_progress(progress_enabled, f"Writing CSV: {csv_path}")
        write_csv(
            csv_path,
            specs,
            unpol_refs_map=unpol_refs_map,
            pol_refs_map=pol_refs_map,
            ref_source=ref_source,
        )
        sys.stdout.write(f"Wrote CSV: {csv_path}\n")
    emit_progress(progress_enabled, "Done")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
