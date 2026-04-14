#!/usr/bin/env python3
"""Verify the pure-Z weighted-integrand closure audit."""

from __future__ import annotations

import argparse
import math
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


TRACE_RE = re.compile(r'([A-Za-z0-9_]+)=(".*?"|\S+)')
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-10


@dataclass(frozen=True)
class TraceIssue:
    line_number: int
    key: tuple[str, str, str, str]
    layer: str
    mode: str
    issue: str


def is_close(lhs: float, rhs: float, *, abs_tol: float = ABS_TOL,
             rel_tol: float = REL_TOL) -> bool:
    return math.isclose(lhs, rhs, abs_tol=abs_tol, rel_tol=rel_tol)


def parse_trace_line(line: str) -> dict[str, str] | None:
    marker = line.find("DIS_WEIGHTED_TRACE")
    if marker < 0:
        return None
    line = line[marker:]
    data: dict[str, str] = {}
    for key, raw_value in TRACE_RE.findall(line):
        if raw_value.startswith('"') and raw_value.endswith('"'):
            data[key] = raw_value[1:-1]
        else:
            data[key] = raw_value
    return data


def closure_prefix(layer: str) -> str:
    if layer == "parton":
        return "parton"
    if layer.startswith("xcomb"):
        return "xcomb"
    if layer == "sampler":
        return "sampler"
    raise ValueError(f"Unsupported layer {layer}")


def record_key(record: dict[str, str]) -> tuple[str, str, str, str]:
    return (
        record.get("run", "NA"),
        record.get("audit_id", "0"),
        record.get("source", "NA"),
        record.get("term", "NA"),
    )


def float_field(record: dict[str, str], key: str) -> float:
    return float(record[key])


def check_closure(record: dict[str, str]) -> list[str]:
    layer = record["layer"]
    prefix = closure_prefix(layer)
    w00 = float_field(record, f"{prefix}_00")
    wpp = float_field(record, f"{prefix}_pp")
    wpm = float_field(record, f"{prefix}_pm")
    wmp = float_field(record, f"{prefix}_mp")
    wmm = float_field(record, f"{prefix}_mm")
    residual = float_field(record, f"{prefix}_closure_residual")
    rel = float_field(record, f"{prefix}_closure_rel")

    exp_residual = 0.25 * (wpp + wpm + wmp + wmm) - w00
    scale = max(abs(w00), abs(0.25 * (wpp + wpm + wmp + wmm)), 1.0e-30)
    exp_rel = exp_residual / scale

    issues: list[str] = []
    if not is_close(residual, exp_residual):
        issues.append(
            f"{layer} residual mismatch: actual={residual:.16e}, expected={exp_residual:.16e}"
        )
    if not is_close(rel, exp_rel):
        issues.append(
            f"{layer} relative residual mismatch: actual={rel:.16e}, expected={exp_rel:.16e}"
        )
    return issues


def check_xcomb_head(record: dict[str, str],
                     parton_run: dict[str, str] | None) -> list[str]:
    issues = check_closure(record)
    if parton_run is None:
        issues.append("missing matching parton run record for xcomb_head")
        return issues

    factor = (
        float_field(record, "pdf_weight") *
        float_field(record, "cut_weight") *
        float_field(record, "ckkw_factor") *
        float_field(record, "me_reweight_factor")
    )
    for suffix in ("00", "pp", "pm", "mp", "mm"):
        actual = float_field(record, f"xcomb_{suffix}")
        expected = float_field(parton_run, f"parton_{suffix}") * factor
        if not is_close(actual, expected):
            issues.append(
                f"xcomb_head propagation mismatch for {suffix}: "
                f"actual={actual:.16e}, expected={expected:.16e}"
            )
    return issues


def check_xcomb_group_total(record: dict[str, str]) -> list[str]:
    issues = check_closure(record)
    total = float_field(record, "group_total_xsec")
    head = float_field(record, "head_xsec")
    dep = float_field(record, "dep_xsec")
    expected = head + dep
    if not is_close(total, expected, abs_tol=1.0e-10, rel_tol=1.0e-8):
        issues.append(
            f"xcomb_group_total scalar mismatch: actual={total:.16e}, "
            f"expected={expected:.16e}"
        )
    return issues


def check_sampler(record: dict[str, str],
                  xcomb_run: dict[str, str] | None) -> list[str]:
    issues = check_closure(record)
    if xcomb_run is None:
        issues.append("missing matching xcomb record for sampler")
        return issues

    factor = float_field(record, "lumi_jac") * float_field(record, "lumi_value")
    xcomb_prefix = "xcomb"
    for suffix in ("00", "pp", "pm", "mp", "mm"):
        actual = float_field(record, f"sampler_{suffix}")
        expected = float_field(xcomb_run, f"{xcomb_prefix}_{suffix}") * factor
        if not is_close(actual, expected):
            issues.append(
                f"sampler propagation mismatch for {suffix}: "
                f"actual={actual:.16e}, expected={expected:.16e}"
            )
    return issues


def parse_records(path: Path) -> tuple[int, list[tuple[int, dict[str, str]]]]:
    trace_lines = 0
    records: list[tuple[int, dict[str, str]]] = []
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = parse_trace_line(line)
            if not record:
                continue
            trace_lines += 1
            records.append((line_number, record))
    return trace_lines, records


def filter_records(records: list[tuple[int, dict[str, str]]], *,
                   layer_filter: str,
                   source_filter: str,
                   raw_or_run: str) -> list[tuple[int, dict[str, str]]]:
    filtered: list[tuple[int, dict[str, str]]] = []
    for line_number, record in records:
        layer = record.get("layer", "")
        mode = record.get("mode", "run")
        source = record.get("source", "")
        if layer_filter != "all" and layer != layer_filter:
            continue
        if source_filter != "all" and source != source_filter:
            continue
        if raw_or_run != "all" and mode != raw_or_run:
            continue
        filtered.append((line_number, record))
    return filtered


def build_lookup(records: list[tuple[int, dict[str, str]]]) -> dict[tuple[str, str, str, str], dict[tuple[str, str], dict[str, str]]]:
    lookup: dict[tuple[str, str, str, str], dict[tuple[str, str], dict[str, str]]] = {}
    for _, record in records:
        key = record_key(record)
        lookup.setdefault(key, {})
        lookup[key][(record["layer"], record.get("mode", "run"))] = record
    return lookup


def choose_xcomb_record(group: dict[tuple[str, str], dict[str, str]]) -> dict[str, str] | None:
    if ("xcomb_group_total", "run") in group:
        return group[("xcomb_group_total", "run")]
    if ("xcomb_head", "run") in group:
        return group[("xcomb_head", "run")]
    return None


def run_checks(records: list[tuple[int, dict[str, str]]], *,
               layer_filter: str,
               source_filter: str,
               raw_or_run: str) -> tuple[list[TraceIssue], Counter]:
    filtered = filter_records(
        records,
        layer_filter=layer_filter,
        source_filter=source_filter,
        raw_or_run=raw_or_run,
    )
    lookup = build_lookup(records)
    issues: list[TraceIssue] = []
    counts = Counter()

    for line_number, record in filtered:
        key = record_key(record)
        counts[(record["layer"], record["source"], record.get("sign", "NA"),
                record.get("family", "NA"))] += 1
        group = lookup.get(key, {})

        layer = record["layer"]
        local_issues: list[str]
        if layer == "parton":
            local_issues = check_closure(record)
        elif layer == "xcomb_head":
            local_issues = check_xcomb_head(record, group.get(("parton", "run")))
        elif layer == "xcomb_group_total":
            local_issues = check_xcomb_group_total(record)
        elif layer == "sampler":
            local_issues = check_sampler(record, choose_xcomb_record(group))
        else:
            local_issues = [f"unsupported layer {layer}"]

        for issue in local_issues:
            issues.append(
                TraceIssue(
                    line_number=line_number,
                    key=key,
                    layer=layer,
                    mode=record.get("mode", "run"),
                    issue=issue,
                )
            )

    return issues, counts


def synthetic_records(*, omit_pdf: bool = False,
                      omit_lumi: bool = False) -> list[tuple[int, dict[str, str]]]:
    parton = {
        "run": "testrun",
        "audit_id": "1",
        "source": "me2",
        "term": "LO",
        "layer": "parton",
        "mode": "run",
        "qid": "1",
        "sign": "quark",
        "family": "down",
        "parton_c00": "1.0",
        "parton_cpl": "0.0",
        "parton_cpp": "0.5",
        "parton_cll": "0.5",
        "parton_00": "1.0",
        "parton_pp": "2.0",
        "parton_pm": "0.0",
        "parton_mp": "2.0",
        "parton_mm": "0.0",
        "parton_closure_residual": "0.0",
        "parton_closure_rel": "0.0",
        "nlo_wgt_raw": "0.0",
        "nlo_wgt_run": "0.0",
        "sigma_born": "1.0",
        "sigma_full_raw": "0.0",
        "sigma_full_run": "0.0",
    }
    factor_x = 6.0
    factor_s = 10.0
    x_factor = 3.0 if omit_pdf else factor_x
    s_factor = 1.0 if omit_lumi else factor_s
    xcomb = {
        "run": "testrun",
        "audit_id": "1",
        "source": "me2",
        "term": "LO",
        "layer": "xcomb_head",
        "mode": "run",
        "qid": "1",
        "sign": "quark",
        "family": "down",
        "pdf_weight": "2.0",
        "cut_weight": "3.0",
        "ckkw_factor": "1.0",
        "me_reweight_factor": "1.0",
        "head_xsec": f"{2.0 * factor_x:.16e}",
        "dep_xsec": "0.0",
        "group_total_xsec": f"{2.0 * factor_x:.16e}",
        "xcomb_00": f"{1.0 * x_factor:.16e}",
        "xcomb_pp": f"{2.0 * x_factor:.16e}",
        "xcomb_pm": "0.0",
        "xcomb_mp": f"{2.0 * x_factor:.16e}",
        "xcomb_mm": "0.0",
        "xcomb_closure_residual": "0.0",
        "xcomb_closure_rel": "0.0",
    }
    sampler = {
        "run": "testrun",
        "audit_id": "1",
        "source": "me2",
        "term": "LO",
        "layer": "sampler",
        "mode": "run",
        "qid": "1",
        "sign": "quark",
        "family": "down",
        "lumi_jac": "2.0",
        "lumi_value": "5.0",
        "sampler_integrand": f"{1.0 * factor_x * factor_s:.16e}",
        "sampler_00": f"{1.0 * factor_x * s_factor:.16e}",
        "sampler_pp": f"{2.0 * factor_x * s_factor:.16e}",
        "sampler_pm": "0.0",
        "sampler_mp": f"{2.0 * factor_x * s_factor:.16e}",
        "sampler_mm": "0.0",
        "sampler_closure_residual": "0.0",
        "sampler_closure_rel": "0.0",
    }
    return [(1, parton), (2, xcomb), (3, sampler)]


def run_self_test() -> int:
    good_issues, _ = run_checks(
        synthetic_records(),
        layer_filter="all",
        source_filter="all",
        raw_or_run="all",
    )
    if good_issues:
        print("Self-test failed on the nominal synthetic sample.", file=sys.stderr)
        for issue in good_issues:
            print(f"- {issue.issue}", file=sys.stderr)
        return 1

    bad_pdf_issues, _ = run_checks(
        synthetic_records(omit_pdf=True),
        layer_filter="all",
        source_filter="all",
        raw_or_run="all",
    )
    if not bad_pdf_issues:
        print("Self-test failed: omitting the PDF factor did not trigger a failure.", file=sys.stderr)
        return 1

    bad_lumi_issues, _ = run_checks(
        synthetic_records(omit_lumi=True),
        layer_filter="all",
        source_filter="all",
        raw_or_run="all",
    )
    if not bad_lumi_issues:
        print("Self-test failed: omitting the luminosity factor did not trigger a failure.", file=sys.stderr)
        return 1

    print("Self-test passed: nominal synthetic records validate, and the PDF/luminosity negative controls fail.")
    return 0


def render_counts(counts: Counter) -> None:
    print("Class counts:")
    for key, count in sorted(counts.items()):
        print(f"- layer={key[0]} source={key[1]} sign={key[2]} family={key[3]} n={count}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-log", help="Path to a Herwig .log file or extracted grep text.")
    parser.add_argument("--layer-filter",
                        choices=["parton", "xcomb_head", "xcomb_group_total", "sampler", "all"],
                        default="all")
    parser.add_argument("--source-filter",
                        choices=["me2", "NLOWeight", "all"],
                        default="all")
    parser.add_argument("--raw-or-run",
                        choices=["raw", "run", "all"],
                        default="all")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)

    if args.self_test:
        return run_self_test()

    if not args.from_log:
        parser.error("--from-log is required unless --self-test is used")

    trace_path = Path(args.from_log)
    if not trace_path.exists():
        print(f"Trace file not found: {trace_path}", file=sys.stderr)
        return 1

    trace_lines, records = parse_records(trace_path)
    if trace_lines == 0:
        print("No DIS_WEIGHTED_TRACE records were found.", file=sys.stderr)
        return 1

    issues, counts = run_checks(
        records,
        layer_filter=args.layer_filter,
        source_filter=args.source_filter,
        raw_or_run=args.raw_or_run,
    )
    filtered = filter_records(
        records,
        layer_filter=args.layer_filter,
        source_filter=args.source_filter,
        raw_or_run=args.raw_or_run,
    )
    if not filtered:
        print(
            "DIS_WEIGHTED_TRACE records were found, but none match the requested "
            "layer/source/mode filters.",
            file=sys.stderr,
        )
        return 1

    print(f"Trace file: {trace_path}")
    print(f"Weighted trace records checked: {len(filtered)}")
    render_counts(counts)

    if issues:
        print("\nWeighted-integrand mismatches:")
        for issue in issues[:20]:
            run, audit_id, source, term = issue.key
            print(
                f"- line={issue.line_number} run={run} audit_id={audit_id} "
                f"source={source} term={term} layer={issue.layer} "
                f"mode={issue.mode}: {issue.issue}"
            )
        if len(issues) > 20:
            print(f"- ... {len(issues) - 20} more mismatches omitted")
        return 1

    print(
        "\nTrace result: every checked weighted layer satisfies its logged closure, "
        "and the xcomb/sampler propagation matches the recorded PDF and luminosity factors."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
