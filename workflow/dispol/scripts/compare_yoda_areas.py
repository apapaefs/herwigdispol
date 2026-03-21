#!/usr/bin/env python3.10
"""
compare_yoda_areas.py

Compare histogram integrals between a Herwig analyzed YODA and a POLDIS
reference YODA.

The script reports, for each requested observable:
- the true area: sum(bin_value * bin_width)
- the propagated area error
- the raw sum of bin heights/values

This is useful for diagnosing whether a mismatch is a real normalization
difference or a per-bin versus density convention issue.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

import yoda


DEFAULT_LABELS = ["Q2", "XBj", "Pt", "Mjj", "Eta", "Zeta", "pT1", "pT2", "pT2OverpT1", "pTAsym"]


def _resolve_path(aos: dict, label: str, ref: bool) -> str:
    prefixes = ["/REF/MC_DIS_BREIT", "/MC_DIS_BREIT"] if ref else ["/MC_DIS_BREIT", "/RAW/MC_DIS_BREIT"]
    candidates = [f"{prefix}/{label}" for prefix in prefixes]
    for candidate in candidates:
        if candidate in aos:
            return candidate
    preview = ", ".join(sorted(aos.keys())[:8])
    raise KeyError(f"Could not find {label}. Tried {candidates}. Available examples: {preview}")


def _bin_edges(obj) -> Sequence[float]:
    if hasattr(obj, "xEdges"):
        return obj.xEdges()
    raise RuntimeError(f"Object of type {type(obj).__name__} does not expose xEdges()")


def _bin_value_and_error(obj, index: int) -> Tuple[float, float]:
    bin_obj = obj.bin(index)
    if hasattr(bin_obj, "val"):
        value = float(bin_obj.val())
        error = 0.0
        if hasattr(bin_obj, "errAvg"):
            try:
                error = float(bin_obj.errAvg("stats"))
            except Exception:
                try:
                    error = float(bin_obj.errAvg())
                except Exception:
                    error = 0.0
        elif hasattr(bin_obj, "errMinus") and hasattr(bin_obj, "errPlus"):
            try:
                error = 0.5 * (abs(float(bin_obj.errMinus())) + abs(float(bin_obj.errPlus())))
            except Exception:
                error = 0.0
        return value, error

    if hasattr(bin_obj, "height"):
        value = float(bin_obj.height())
        error = float(bin_obj.heightErr()) if hasattr(bin_obj, "heightErr") else 0.0
        return value, error

    raise RuntimeError(f"Unsupported YODA bin type: {type(bin_obj).__name__}")


def histogram_summary(aos: dict, label: str, ref: bool) -> dict:
    path = _resolve_path(aos, label, ref=ref)
    obj = aos[path]
    edges = _bin_edges(obj)
    area = 0.0
    area_err2 = 0.0
    sum_heights = 0.0
    skipped = 0
    nbins = len(edges) - 1
    for index in range(nbins):
        value, error = _bin_value_and_error(obj, index)
        if not math.isfinite(value):
            skipped += 1
            continue
        if not math.isfinite(error):
            error = 0.0
        width = float(edges[index + 1] - edges[index])
        area += value * width
        area_err2 += (error * width) ** 2
        sum_heights += value
    return {
        "path": path,
        "nbins": nbins,
        "area": area,
        "area_err": math.sqrt(area_err2),
        "sum_heights": sum_heights,
        "skipped": skipped,
    }


def fmt_pm(value: float, error: float, precision: int = 6) -> str:
    return f"{value:.{precision}f} +- {error:.{precision}f}"


def compare_rows(herwig_aos: dict, poldis_aos: dict, labels: Iterable[str]) -> List[List[str]]:
    rows: List[List[str]] = []
    for label in labels:
        h = histogram_summary(herwig_aos, label, ref=False)
        p = histogram_summary(poldis_aos, label, ref=True)
        ratio = float("nan")
        ratio_err = float("nan")
        if p["area"] != 0.0:
            ratio = h["area"] / p["area"]
            ratio_err = math.sqrt(
                (h["area_err"] / p["area"]) ** 2
                + ((h["area"] * p["area_err"]) / (p["area"] ** 2)) ** 2
            )
        rows.append(
            [
                label,
                str(h["nbins"]),
                fmt_pm(h["area"], h["area_err"]),
                fmt_pm(p["area"], p["area_err"]),
                fmt_pm(h["sum_heights"], 0.0),
                fmt_pm(p["sum_heights"], 0.0),
                f'{h["skipped"]}/{p["skipped"]}',
                f"{ratio:.6f} +- {ratio_err:.6f}" if math.isfinite(ratio) else "nan",
            ]
        )
    return rows


def render_table(headers: Sequence[str], rows: Sequence[Sequence[str]]) -> str:
    try:
        from prettytable import PrettyTable

        table = PrettyTable()
        table.field_names = list(headers)
        for row in rows:
            table.add_row(list(row))
        return table.get_string()
    except Exception:
        widths = [len(h) for h in headers]
        for row in rows:
            for i, cell in enumerate(row):
                widths[i] = max(widths[i], len(str(cell)))

        def fmt_row(row: Sequence[str]) -> str:
            return " | ".join(str(cell).ljust(widths[i]) for i, cell in enumerate(row))

        sep = "-+-".join("-" * w for w in widths)
        out = [fmt_row(headers), sep]
        out.extend(fmt_row(row) for row in rows)
        return "\n".join(out)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare areas of Herwig and POLDIS YODA histograms.")
    parser.add_argument("herwig_yoda", help="Analyzed Herwig YODA file")
    parser.add_argument("poldis_yoda", help="POLDIS reference YODA file")
    parser.add_argument(
        "--labels",
        default=",".join(DEFAULT_LABELS),
        help="Comma-separated list of labels to compare (default: Q2,XBj,Pt,Mjj,Eta,Zeta,pT1,pT2,pT2OverpT1,pTAsym)",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    herwig_path = Path(args.herwig_yoda)
    poldis_path = Path(args.poldis_yoda)
    herwig_aos = yoda.read(str(herwig_path))
    poldis_aos = yoda.read(str(poldis_path))
    labels = [label.strip() for label in args.labels.split(",") if label.strip()]

    rows = compare_rows(herwig_aos, poldis_aos, labels)
    print(f"Herwig: {herwig_path}")
    print(f"POLDIS: {poldis_path}")
    print(render_table(
        ["Observable", "Bins", "Herwig area [pb]", "POLDIS area [pb]", "Herwig sum(y)", "POLDIS sum(y)", "Skipped H/P", "Area ratio"],
        rows,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
