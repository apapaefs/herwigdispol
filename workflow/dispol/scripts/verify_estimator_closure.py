#!/usr/bin/env python3
"""Verify the pure-Z estimator-level closure audit."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path


FIELD_RE = re.compile(r'([A-Za-z0-9_]+)=(".*?"|\S+)')
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-10


@dataclass(frozen=True)
class Issue:
    line_number: int
    marker: str
    key: tuple[str, str, str, str]
    issue: str


def is_close(lhs: float, rhs: float, *, abs_tol: float = ABS_TOL,
             rel_tol: float = REL_TOL) -> bool:
    return math.isclose(lhs, rhs, abs_tol=abs_tol, rel_tol=rel_tol)


def parse_fields(line: str, marker: str) -> dict[str, str] | None:
    pos = line.find(marker)
    if pos < 0:
        return None
    data: dict[str, str] = {}
    for key, raw in FIELD_RE.findall(line[pos:]):
        if raw.startswith('"') and raw.endswith('"'):
            data[key] = raw[1:-1]
        else:
            data[key] = raw
    return data


def float_field(record: dict[str, str], key: str) -> float:
    return float(record[key])


def trace_key(record: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        record.get("run", "NA"),
        record.get("audit_id", "0"),
        record.get("source", "NA"),
        record.get("term", "NA"),
    )


def closure_residual(w00: float, wpp: float, wpm: float, wmp: float, wmm: float) -> float:
    return 0.25 * (wpp + wpm + wmp + wmm) - w00


def closure_relative(residual: float, w00: float, wpp: float, wpm: float,
                     wmp: float, wmm: float) -> float:
    scale = max(abs(w00), abs(0.25 * (wpp + wpm + wmp + wmm)), 1.0e-30)
    return residual / scale


def parse_log(path: Path) -> tuple[
    list[tuple[int, dict[str, str]]],
    list[tuple[int, dict[str, str]]],
    list[tuple[int, dict[str, str]]],
]:
    weighted: list[tuple[int, dict[str, str]]] = []
    estimator: list[tuple[int, dict[str, str]]] = []
    totals: list[tuple[int, dict[str, str]]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = parse_fields(line, "DIS_WEIGHTED_TRACE")
            if record:
                weighted.append((line_number, record))
                continue
            record = parse_fields(line, "DIS_ESTIMATOR_TRACE")
            if record:
                estimator.append((line_number, record))
                continue
            record = parse_fields(line, "DIS_ESTIMATOR_TOTAL")
            if record:
                totals.append((line_number, record))
    return weighted, estimator, totals


def make_weighted_lookup(
    weighted_records: list[tuple[int, dict[str, str]]],
) -> dict[tuple[str, str, str, str], dict[str, dict[str, str]]]:
    lookup: dict[tuple[str, str, str, str], dict[str, dict[str, str]]] = {}
    for _, record in weighted_records:
        key = trace_key(record)
        lookup.setdefault(key, {})
        lookup[key][record.get("layer", "")] = record
    return lookup


def make_estimator_lookup(
    estimator_records: list[tuple[int, dict[str, str]]],
) -> dict[tuple[str, str, str, str], dict[str, dict[str, str]]]:
    lookup: dict[tuple[str, str, str, str], dict[str, dict[str, str]]] = {}
    for _, record in estimator_records:
        key = trace_key(record)
        lookup.setdefault(key, {})
        lookup[key][record.get("layer", "")] = record
    return lookup


def check_trace_closure(record: dict[str, str]) -> list[str]:
    w00 = float_field(record, "w00")
    wpp = float_field(record, "wPP")
    wpm = float_field(record, "wPM")
    wmp = float_field(record, "wMP")
    wmm = float_field(record, "wMM")
    residual = float_field(record, "closure_residual")
    relative = float_field(record, "closure_relative")
    exp_residual = closure_residual(w00, wpp, wpm, wmp, wmm)
    exp_relative = closure_relative(exp_residual, w00, wpp, wpm, wmp, wmm)
    issues: list[str] = []
    if not is_close(residual, exp_residual):
        issues.append(
            f"closure residual mismatch: actual={residual:.16e}, expected={exp_residual:.16e}"
        )
    if not is_close(relative, exp_relative):
        issues.append(
            f"closure relative mismatch: actual={relative:.16e}, expected={exp_relative:.16e}"
        )
    return issues


def check_bin_evaluate(record: dict[str, str],
                       weighted_upstream: dict[str, dict[str, str]] | None) -> list[str]:
    issues = check_trace_closure(record)
    if not weighted_upstream or "sampler" not in weighted_upstream:
        issues.append("missing matching DIS_WEIGHTED_TRACE sampler record")
        return issues
    upstream = weighted_upstream["sampler"]
    factor = float_field(record, "remap_factor")
    for suffix in ("00", "PP", "PM", "MP", "MM"):
        actual = float_field(record, f"w{suffix}")
        expected = float_field(upstream, f"sampler_{suffix.lower()}") * factor
        if not is_close(actual, expected):
            issues.append(
                f"bin_evaluate propagation mismatch for {suffix}: "
                f"actual={actual:.16e}, expected={expected:.16e}"
            )
    return issues


def check_general_generate(record: dict[str, str],
                           estimator_upstream: dict[str, dict[str, str]] | None) -> list[str]:
    issues = check_trace_closure(record)
    if not estimator_upstream or "bin_evaluate" not in estimator_upstream:
        issues.append("missing matching bin_evaluate record")
        return issues
    upstream = estimator_upstream["bin_evaluate"]
    ref = float_field(record, "reference_weight")
    if abs(ref) <= 1.0e-30:
        issues.append("reference_weight is too small for propagation test")
        return issues
    factor = float_field(record, "selection_factor") / ref
    for suffix in ("00", "PP", "PM", "MP", "MM"):
        actual = float_field(record, f"w{suffix}")
        expected = float_field(upstream, f"w{suffix}") * factor
        if not is_close(actual, expected):
            issues.append(
                f"general_generate propagation mismatch for {suffix}: "
                f"actual={actual:.16e}, expected={expected:.16e}"
            )
    return issues


def check_eventhandler_select(record: dict[str, str],
                              estimator_upstream: dict[str, dict[str, str]] | None) -> list[str]:
    issues = check_trace_closure(record)
    if not estimator_upstream or "general_generate" not in estimator_upstream:
        issues.append("missing matching general_generate record")
        return issues
    upstream = estimator_upstream["general_generate"]
    pre_weight = float_field(record, "pre_weight")
    if abs(pre_weight) <= 1.0e-30:
        issues.append("pre_weight is too small for propagation test")
        return issues
    factor = 1.0 / pre_weight
    for suffix in ("00", "PP", "PM", "MP", "MM"):
        actual = float_field(record, f"w{suffix}")
        expected = float_field(upstream, f"w{suffix}") * factor
        if not is_close(actual, expected):
            issues.append(
                f"eventhandler_select propagation mismatch for {suffix}: "
                f"actual={actual:.16e}, expected={expected:.16e}"
            )
    return issues


def check_total(record: dict[str, str]) -> list[str]:
    w00 = float_field(record, "total_00")
    wpp = float_field(record, "total_pp")
    wpm = float_field(record, "total_pm")
    wmp = float_field(record, "total_mp")
    wmm = float_field(record, "total_mm")
    residual = float_field(record, "closure_residual")
    relative = float_field(record, "closure_relative")
    exp_residual = closure_residual(w00, wpp, wpm, wmp, wmm)
    exp_relative = closure_relative(exp_residual, w00, wpp, wpm, wmp, wmm)
    issues: list[str] = []
    if not is_close(residual, exp_residual):
        issues.append(
            f"total residual mismatch: actual={residual:.16e}, expected={exp_residual:.16e}"
        )
    if not is_close(relative, exp_relative):
        issues.append(
            f"total relative mismatch: actual={relative:.16e}, expected={exp_relative:.16e}"
        )
    return issues


def run_trace_checks(
    weighted_records: list[tuple[int, dict[str, str]]],
    estimator_records: list[tuple[int, dict[str, str]]],
    *,
    layer_filter: str,
    source_filter: str,
) -> tuple[list[Issue], Counter]:
    weighted_lookup = make_weighted_lookup(weighted_records)
    estimator_lookup = make_estimator_lookup(estimator_records)
    issues: list[Issue] = []
    counts = Counter()

    for line_number, record in estimator_records:
        if record.get("mode", "run") != "run":
            continue
        if layer_filter != "all" and record.get("layer") != layer_filter:
            continue
        if source_filter != "all" and record.get("source") != source_filter:
            continue
        key = trace_key(record)
        counts[(record["layer"], record["source"], record.get("sign", "NA"),
                record.get("family", "NA"))] += 1
        if record["layer"] == "bin_evaluate":
            local = check_bin_evaluate(record, weighted_lookup.get(key))
        elif record["layer"] == "general_generate":
            local = check_general_generate(record, estimator_lookup.get(key))
        elif record["layer"] == "eventhandler_select":
            local = check_eventhandler_select(record, estimator_lookup.get(key))
        else:
            local = [f"unsupported estimator layer {record.get('layer', 'NA')}"]
        for issue in local:
            issues.append(Issue(line_number, "DIS_ESTIMATOR_TRACE", key, issue))

    return issues, counts


def run_total_checks(
    total_records: list[tuple[int, dict[str, str]]],
    *,
    estimator_filter: str,
    source_filter: str,
) -> tuple[list[Issue], Counter]:
    issues: list[Issue] = []
    counts = Counter()
    for line_number, record in total_records:
        estimator = record.get("estimator", "")
        if estimator_filter != "all" and estimator != estimator_filter:
            continue
        source = record.get("key", "")
        if source_filter != "all" and f"source={source_filter}" not in source:
            continue
        counts[(estimator, source)] += 1
        for issue in check_total(record):
            issues.append(Issue(line_number, "DIS_ESTIMATOR_TOTAL", ("NA", "0", estimator, source), issue))
    return issues, counts


def make_trace_record(**fields: str | float | int) -> dict[str, str]:
    return {key: str(value) for key, value in fields.items()}


def synthetic_trace_set(*, bad_remap: bool = False, bad_selection: bool = False) -> tuple[
    list[tuple[int, dict[str, str]]],
    list[tuple[int, dict[str, str]]],
]:
    weighted = [(
        1,
        make_trace_record(
            run="testrun", audit_id=1, source="me2", term="LO", layer="sampler",
            sampler_00=2.0, sampler_pp=3.0, sampler_pm=1.0, sampler_mp=1.0,
            sampler_mm=3.0, sampler_closure_residual=0.0, sampler_closure_rel=0.0,
        ),
    )]
    remap00 = 5.0 if bad_remap else 4.0
    bin_eval = (
        2,
        make_trace_record(
            run="testrun", audit_id=1, source="me2", term="LO", layer="bin_evaluate",
            mode="run", sign="quark", family="down", remap_factor=2.0,
            reference_weight=4.0, bias=1.0, global_max_weight=2.0,
            selection_factor=2.0, pre_weight=1.0, eventhandler_weight=0.0,
            actual_weight=4.0, w00=remap00, wPP=6.0, wPM=2.0, wMP=2.0, wMM=6.0,
            closure_residual=0.0, closure_relative=0.0,
        ),
    )
    gen_factor = 0.25 if bad_selection else 0.5
    general = (
        3,
        make_trace_record(
            run="testrun", audit_id=1, source="me2", term="LO", layer="general_generate",
            mode="run", sign="quark", family="down", remap_factor=1.0,
            reference_weight=4.0, bias=1.0, global_max_weight=2.0,
            selection_factor=2.0, pre_weight=1.0, eventhandler_weight=0.0,
            actual_weight=2.0, w00=1.0,
            wPP=3.0, wPM=1.0, wMP=1.0, wMM=3.0,
            closure_residual=0.0, closure_relative=0.0,
        ),
    )
    # Fix possibly string-valued w00 after construction.
    general[1]["w00"] = str(float(bin_eval[1]["w00"]) * gen_factor)
    general[1]["wPP"] = str(6.0 * gen_factor)
    general[1]["wPM"] = str(2.0 * gen_factor)
    general[1]["wMP"] = str(2.0 * gen_factor)
    general[1]["wMM"] = str(6.0 * gen_factor)
    general[1]["actual_weight"] = str(4.0 * gen_factor)
    event = (
        4,
        make_trace_record(
            run="testrun", audit_id=1, source="me2", term="LO", layer="eventhandler_select",
            mode="run", sign="quark", family="down", remap_factor=1.0,
            reference_weight=1.0, bias=1.0, global_max_weight=1.0,
            selection_factor=1.0, pre_weight=2.0, eventhandler_weight=1.0,
            actual_weight=1.0, w00=1.0, wPP=1.5, wPM=0.5, wMP=0.5, wMM=1.5,
            closure_residual=0.0, closure_relative=0.0,
        ),
    )
    return weighted, [bin_eval, general, event]


def synthetic_totals(*, bad_veto: bool = False) -> list[tuple[int, dict[str, str]]]:
    event_total_00 = 1.0 if not bad_veto else 1.2
    return [(
        10,
        make_trace_record(
            run="testrun", estimator="eventhandler_generated",
            key="source=me2,term=LO", total_00=event_total_00,
            total_pp=1.5, total_pm=0.5, total_mp=0.5, total_mm=1.5,
            closure_residual=0.0,
            closure_relative=0.0,
            actual_xsec_nb=1.0, actual_xsec_err_nb=0.1,
        ),
    )]


def run_self_test() -> None:
    weighted, estimator = synthetic_trace_set()
    issues, _ = run_trace_checks(weighted, estimator, layer_filter="all", source_filter="all")
    if issues:
        raise SystemExit(f"self-test baseline failed: {issues[0].issue}")

    totals = synthetic_totals()
    total_issues, _ = run_total_checks(totals, estimator_filter="all", source_filter="all")
    if total_issues:
        raise SystemExit(f"self-test total baseline failed: {total_issues[0].issue}")

    weighted_bad, estimator_bad = synthetic_trace_set(bad_remap=True)
    bad_issues, _ = run_trace_checks(weighted_bad, estimator_bad, layer_filter="bin_evaluate", source_filter="all")
    if not bad_issues:
        raise SystemExit("negative control failed: omitted remap factor was not detected")

    weighted_bad, estimator_bad = synthetic_trace_set(bad_selection=True)
    bad_issues, _ = run_trace_checks(weighted_bad, estimator_bad, layer_filter="general_generate", source_filter="all")
    if not bad_issues:
        raise SystemExit("negative control failed: omitted selection compensation was not detected")

    total_issues, _ = run_total_checks(synthetic_totals(bad_veto=True), estimator_filter="all", source_filter="all")
    if not total_issues:
        raise SystemExit("negative control failed: omitted veto correction was not detected")

    print("Self-test passed.")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-log", type=Path, help="Path to the Herwig .log file or filtered estimator log.")
    parser.add_argument("--layer-filter", choices=["bin_evaluate", "general_generate", "eventhandler_select", "all"], default="all")
    parser.add_argument("--estimator-filter", choices=["sampler_run", "sampler_combined", "eventhandler_generated", "all"], default="all")
    parser.add_argument("--source-filter", choices=["me2", "NLOWeight", "all"], default="all")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        run_self_test()
        return 0

    if not args.from_log:
        parser.error("--from-log is required unless --self-test is used")

    weighted, estimator, totals = parse_log(args.from_log)
    if not estimator:
        print(f"Trace file: {args.from_log}")
        print("No DIS_ESTIMATOR_TRACE records found. The log predates this audit or the new binary was not used.")
        return 1

    trace_issues, trace_counts = run_trace_checks(
        weighted, estimator,
        layer_filter=args.layer_filter,
        source_filter=args.source_filter,
    )
    total_issues, total_counts = run_total_checks(
        totals,
        estimator_filter=args.estimator_filter,
        source_filter=args.source_filter,
    )

    print(f"Trace file: {args.from_log}")
    checked_trace_records = sum(trace_counts.values())
    print(f"Estimator trace records checked: {checked_trace_records}")
    if trace_counts:
        print("Trace class counts:")
        for (layer, source, sign, family), count in sorted(trace_counts.items()):
            print(f"- layer={layer} source={source} sign={sign} family={family} n={count}")
    checked_totals = sum(total_counts.values())
    print(f"Estimator totals checked: {checked_totals}")
    if total_counts:
        print("Total class counts:")
        for (estimator_name, key), count in sorted(total_counts.items()):
            print(f"- estimator={estimator_name} key={key} n={count}")

    issues = trace_issues + total_issues
    if issues:
        first = issues[0]
        print(
            f"First issue at line {first.line_number} [{first.marker}] "
            f"key={first.key}: {first.issue}"
        )
        print(f"Total issues: {len(issues)}")
        return 1

    print(
        "Trace result: every checked estimator layer satisfies closure, "
        "the logged propagation factors are consistent, and the checked totals close."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
