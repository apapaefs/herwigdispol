#!/usr/bin/env python3
"""Summarize low-Q2 POWHEG raw diagnostics from Herwig DIS logs."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable


NUMERIC_FIELDS = {
    "Q2",
    "xB",
    "y",
    "s",
    "fallback",
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
}


def parse_summary_line(line: str) -> Dict[str, object] | None:
    if not line.startswith("POWHEG_RAW_SUMMARY"):
        return None
    record: Dict[str, object] = {}
    for token in line.strip().split()[1:]:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        if key in NUMERIC_FIELDS:
            try:
                record[key] = float(value)
            except ValueError:
                record[key] = value
        else:
            record[key] = value
    return record


def iter_summaries(paths: Iterable[Path]) -> Iterable[Dict[str, object]]:
    for path in paths:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                record = parse_summary_line(line)
                if record is not None:
                    record["_source"] = str(path)
                    yield record


def update_channel_stats(stats: Dict[str, float], prefix: str, record: Dict[str, object]) -> None:
    for field in ("trials", "rejectXP", "rejectVeto", "weightNeg", "weightHigh"):
        key = f"{prefix}_{field}"
        value = float(record.get(key, 0.0))
        stats[field] += value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", help="Herwig log files to inspect")
    parser.add_argument("--q2-max", type=float, default=100.0,
                        help="Upper Q^2 cut in GeV^2 for the focused sample (default: 100)")
    parser.add_argument("--q2-min", type=float, default=49.0,
                        help="Lower Q^2 cut in GeV^2 for the focused sample (default: 49)")
    parser.add_argument("--show-sources", action="store_true",
                        help="Also print how many filtered events came from each log")
    args = parser.parse_args()

    paths = [Path(p) for p in args.logs]
    records = list(iter_summaries(paths))
    if not records:
        print("No POWHEG_RAW_SUMMARY lines found.")
        return 1

    filtered = [
        r for r in records
        if args.q2_min <= float(r.get("Q2", -1.0)) <= args.q2_max
    ]

    print(f"Total POWHEG summaries: {len(records)}")
    print(f"Filtered {args.q2_min:g} <= Q2 <= {args.q2_max:g} GeV^2: {len(filtered)}")
    if not filtered:
        return 0

    winners = Counter(str(r.get("winner", "UNKNOWN")) for r in filtered)
    fallbacks = sum(int(float(r.get("fallback", 0.0))) for r in filtered)
    comp_status = Counter(str(r.get("comp_status", "UNKNOWN")) for r in filtered)
    bgf_status = Counter(str(r.get("bgf_status", "UNKNOWN")) for r in filtered)
    source_counts = Counter(str(r.get("_source", "UNKNOWN")) for r in filtered)

    comp_totals = defaultdict(float)
    bgf_totals = defaultdict(float)
    for record in filtered:
        update_channel_stats(comp_totals, "comp", record)
        update_channel_stats(bgf_totals, "bgf", record)

    print("\nWinner counts:")
    for winner, count in winners.items():
        print(f"  {winner}: {count}")

    print(f"\nFallback count: {fallbacks}")

    print("\nCompton status counts:")
    for status, count in comp_status.items():
        print(f"  {status}: {count}")

    print("\nBGF status counts:")
    for status, count in bgf_status.items():
        print(f"  {status}: {count}")

    print("\nCompton aggregate diagnostics:")
    for key in ("trials", "rejectXP", "rejectVeto", "weightNeg", "weightHigh"):
        total = comp_totals[key]
        mean = total / len(filtered)
        print(f"  {key}: total={total:.0f} mean/event={mean:.3f}")

    print("\nBGF aggregate diagnostics:")
    for key in ("trials", "rejectXP", "rejectVeto", "weightNeg", "weightHigh"):
        total = bgf_totals[key]
        mean = total / len(filtered)
        print(f"  {key}: total={total:.0f} mean/event={mean:.3f}")

    if args.show_sources:
        print("\nFiltered events by source log:")
        for source, count in source_counts.items():
            print(f"  {source}: {count}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
