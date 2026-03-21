#!/usr/bin/env python3.10
"""Summarize POWHEG raw-competition CSV diagnostics."""

from __future__ import annotations

import argparse
import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence


SUMMARY_CSV_GLOB = "*.rawpowheg.summary.csv"

NUMERIC_FIELDS = (
    "fallback",
    "Q2",
    "xB",
    "y",
    "s",
    "comp_trials",
    "comp_rejectXP",
    "comp_rejectVeto",
    "comp_weightNeg",
    "comp_weightHigh",
    "comp_xp",
    "comp_zp",
    "comp_xMapped",
    "comp_xT",
    "comp_xTMin",
    "comp_pT",
    "comp_phase",
    "comp_pdfRatio",
    "comp_alphaRatio",
    "comp_meAvg",
    "comp_wgt",
    "comp_pdfScale",
    "comp_alphaScale",
    "bgf_trials",
    "bgf_rejectXP",
    "bgf_rejectVeto",
    "bgf_weightNeg",
    "bgf_weightHigh",
    "bgf_xp",
    "bgf_zp",
    "bgf_xMapped",
    "bgf_xT",
    "bgf_xTMin",
    "bgf_pT",
    "bgf_phase",
    "bgf_pdfRatio",
    "bgf_alphaRatio",
    "bgf_meAvg",
    "bgf_wgt",
    "bgf_pdfScale",
    "bgf_alphaScale",
)


@dataclass(frozen=True)
class SummaryRow:
    source_log: str
    winner: str
    comp_status: str
    bgf_status: str
    fallback: int
    q2: float
    xb: float
    comp_trials: float
    comp_reject_xp: float
    comp_reject_veto: float
    comp_wgt: float
    comp_phase: float
    comp_pdf_ratio: float
    comp_me_avg: float
    comp_xp: float
    comp_xt: float
    comp_xtmin: float
    bgf_trials: float
    bgf_reject_xp: float
    bgf_reject_veto: float
    bgf_wgt: float
    bgf_phase: float
    bgf_pdf_ratio: float
    bgf_me_avg: float
    bgf_xp: float
    bgf_xt: float
    bgf_xtmin: float


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Summarize raw POWHEG competition summary CSV files from files and/or directories."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        help="Input summary CSV files or directories to scan recursively.",
    )
    parser.add_argument(
        "--xb-bin-width",
        type=float,
        default=0.05,
        help="Width of the xB bins used in the printed summary (default: 0.05).",
    )
    parser.add_argument(
        "--xb-max",
        type=float,
        default=0.35,
        help="Maximum xB edge to print (default: 0.35).",
    )
    parser.add_argument(
        "--min-bin-count",
        type=int,
        default=1,
        help="Only print xB bins with at least this many rows (default: 1).",
    )
    return parser.parse_args()


def as_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value == "":
        return 0.0
    return float(value)


def resolve_csv_inputs(inputs: Sequence[str]) -> list[str]:
    resolved: list[str] = []
    seen: set[str] = set()

    for raw_input in inputs:
        path = Path(raw_input).expanduser()
        if not path.exists():
            raise SystemExit(f"Input path does not exist: {path}")
        candidates = [path]
        if path.is_dir():
            candidates = sorted(
                candidate for candidate in path.rglob(SUMMARY_CSV_GLOB) if candidate.is_file()
            )
        for candidate in candidates:
            resolved_path = str(candidate.resolve())
            if resolved_path in seen:
                continue
            seen.add(resolved_path)
            resolved.append(resolved_path)
    return resolved


def iter_rows(paths: Sequence[str]) -> Iterable[SummaryRow]:
    for csv_path in paths:
        path = Path(csv_path)
        with path.open("rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle)
            for raw in reader:
                for key in NUMERIC_FIELDS:
                    raw.setdefault(key, "0")
                yield SummaryRow(
                    source_log=raw.get("source_log", ""),
                    winner=raw.get("winner", ""),
                    comp_status=raw.get("comp_status", ""),
                    bgf_status=raw.get("bgf_status", ""),
                    fallback=int(as_float(raw, "fallback")),
                    q2=as_float(raw, "Q2"),
                    xb=as_float(raw, "xB"),
                    comp_trials=as_float(raw, "comp_trials"),
                    comp_reject_xp=as_float(raw, "comp_rejectXP"),
                    comp_reject_veto=as_float(raw, "comp_rejectVeto"),
                    comp_wgt=as_float(raw, "comp_wgt"),
                    comp_phase=as_float(raw, "comp_phase"),
                    comp_pdf_ratio=as_float(raw, "comp_pdfRatio"),
                    comp_me_avg=as_float(raw, "comp_meAvg"),
                    comp_xp=as_float(raw, "comp_xp"),
                    comp_xt=as_float(raw, "comp_xT"),
                    comp_xtmin=as_float(raw, "comp_xTMin"),
                    bgf_trials=as_float(raw, "bgf_trials"),
                    bgf_reject_xp=as_float(raw, "bgf_rejectXP"),
                    bgf_reject_veto=as_float(raw, "bgf_rejectVeto"),
                    bgf_wgt=as_float(raw, "bgf_wgt"),
                    bgf_phase=as_float(raw, "bgf_phase"),
                    bgf_pdf_ratio=as_float(raw, "bgf_pdfRatio"),
                    bgf_me_avg=as_float(raw, "bgf_meAvg"),
                    bgf_xp=as_float(raw, "bgf_xp"),
                    bgf_xt=as_float(raw, "bgf_xT"),
                    bgf_xtmin=as_float(raw, "bgf_xTMin"),
                )


def mean(rows: Sequence[SummaryRow], getter) -> float:
    if not rows:
        return 0.0
    return sum(getter(row) for row in rows) / len(rows)


def ratio_string(numerator: int, denominator: int) -> str:
    if denominator == 0:
        return "0/0 (0.0%)"
    return f"{numerator}/{denominator} ({100.0 * numerator / denominator:.1f}%)"


def xt_ratio(xt: float, xtmin: float) -> float:
    if xtmin <= 0.0:
        return 0.0
    return xt / xtmin


def print_global_summary(rows: Sequence[SummaryRow], inputs: Sequence[str], paths: Sequence[str]) -> None:
    print("Inputs:")
    for path in inputs:
        print(f"  {path}")
    print(f"Resolved CSV files: {len(paths)}")
    if len(paths) <= 20:
        for path in paths:
            print(f"  {path}")
    else:
        for path in paths[:10]:
            print(f"  {path}")
        print("  ...")
        for path in paths[-10:]:
            print(f"  {path}")
    print(f"Rows: {len(rows)}")
    print(f"Winner counts: {dict(Counter(row.winner for row in rows))}")
    print(f"Fallback counts: {dict(Counter(row.fallback for row in rows))}")
    print(f"Compton status: {dict(Counter(row.comp_status for row in rows))}")
    print(f"BGF status: {dict(Counter(row.bgf_status for row in rows))}")
    print()

    fallback_rows = [row for row in rows if row.fallback == 1]
    print("Averages:")
    print(
        "  "
        f"xB={mean(rows, lambda row: row.xb):.5f} "
        f"Q2={mean(rows, lambda row: row.q2):.2f} "
        f"comp_rejectXP={mean(rows, lambda row: row.comp_reject_xp):.2f} "
        f"comp_rejectVeto={mean(rows, lambda row: row.comp_reject_veto):.2f} "
        f"bgf_rejectXP={mean(rows, lambda row: row.bgf_reject_xp):.2f} "
        f"bgf_rejectVeto={mean(rows, lambda row: row.bgf_reject_veto):.2f}"
    )
    if fallback_rows:
        print(
            "  "
            f"fallback rows={len(fallback_rows)} "
            f"mean comp_wgt={mean(fallback_rows, lambda row: row.comp_wgt):.5e} "
            f"mean bgf_wgt={mean(fallback_rows, lambda row: row.bgf_wgt):.5e} "
            f"mean comp_xT/xTMin={mean(fallback_rows, lambda row: xt_ratio(row.comp_xt, row.comp_xtmin)):.3f} "
            f"mean bgf_xT/xTMin={mean(fallback_rows, lambda row: xt_ratio(row.bgf_xt, row.bgf_xtmin)):.3f}"
        )
    print()


def print_bin_summary(rows: Sequence[SummaryRow], width: float, xb_max: float, min_bin_count: int) -> None:
    print("xB-binned summary:")
    edge = 0.0
    while edge < xb_max:
        next_edge = edge + width
        bin_rows = [row for row in rows if edge <= row.xb < next_edge]
        edge = next_edge
        if len(bin_rows) < min_bin_count:
            continue
        winner_counts = Counter(row.winner for row in bin_rows)
        comp_acc = sum(row.comp_status == "accepted" for row in bin_rows)
        bgf_acc = sum(row.bgf_status == "accepted" for row in bin_rows)
        print(
            f"  xB {next_edge - width:.2f}-{next_edge:.2f} "
            f"n={len(bin_rows)} "
            f"winners={dict(winner_counts)} "
            f"fallback={ratio_string(sum(row.fallback for row in bin_rows), len(bin_rows))}"
        )
        print(
            "    "
            f"comp accepted={ratio_string(comp_acc, len(bin_rows))} "
            f"mean veto={mean(bin_rows, lambda row: row.comp_reject_veto):.2f} "
            f"mean wgt={mean(bin_rows, lambda row: row.comp_wgt):.5f} "
            f"mean xp={mean(bin_rows, lambda row: row.comp_xp):.4f} "
            f"mean phase={mean(bin_rows, lambda row: row.comp_phase):.5f}"
        )
        print(
            "    "
            f"bgf accepted={ratio_string(bgf_acc, len(bin_rows))} "
            f"mean veto={mean(bin_rows, lambda row: row.bgf_reject_veto):.2f} "
            f"mean wgt={mean(bin_rows, lambda row: row.bgf_wgt):.5f} "
            f"mean xp={mean(bin_rows, lambda row: row.bgf_xp):.4f} "
            f"mean phase={mean(bin_rows, lambda row: row.bgf_phase):.5f}"
        )
    print()


def print_channel_subsets(rows: Sequence[SummaryRow]) -> None:
    labels = (
        ("winner=Compton", [row for row in rows if row.winner == "Compton"]),
        ("winner=BGF", [row for row in rows if row.winner == "BGF"]),
        ("winner=None", [row for row in rows if row.winner == "None"]),
    )
    print("Winner subsets:")
    for label, subset in labels:
        if not subset:
            continue
        print(
            "  "
            f"{label}: n={len(subset)} "
            f"mean xB={mean(subset, lambda row: row.xb):.5f} "
            f"mean Q2={mean(subset, lambda row: row.q2):.2f} "
            f"mean comp_veto={mean(subset, lambda row: row.comp_reject_veto):.2f} "
            f"mean bgf_veto={mean(subset, lambda row: row.bgf_reject_veto):.2f}"
        )
    print()


def main() -> int:
    args = parse_args()
    csv_paths = resolve_csv_inputs(args.inputs)
    rows = list(iter_rows(csv_paths))
    if not rows:
        raise SystemExit("No rows found in the supplied summary CSV inputs.")
    print_global_summary(rows, args.inputs, csv_paths)
    print_channel_subsets(rows)
    print_bin_summary(rows, args.xb_bin_width, args.xb_max, args.min_bin_count)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
