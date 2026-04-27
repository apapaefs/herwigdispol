#!/usr/bin/env python3
"""Summarize low-Q2 fixed-order DIS diagnostics from NLO_AUDIT_* and NLO_TERM_* logs."""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


PREFIXES = {
    "NLO_AUDIT_OBJ",
    "NLO_AUDIT_PDF",
    "NLO_AUDIT_TERM",
    "NLO_AUDIT_NC",
    "NLO_TERM_UNIFORM_COMPARE",
    "NLO_TERM_RAWDELTA_COMPARE",
    "NLO_TERM_SIGN",
    "NLO_TERM_CUM",
    "NLO_TERM_REAL",
    "NLO_TERM_APOL",
    "NLO_TERM_BORN",
    "NLO_TERM_COEFF",
}


def parse_value(value: str):
    if value in {"true", "false"}:
        return value == "true"
    try:
        num = float(value)
    except ValueError:
        return value
    if math.isfinite(num) and num.is_integer():
        return int(num)
    return num


def parse_line(line: str):
    line = line.strip()
    if not line:
        return None
    parts = line.split()
    prefix = parts[0]
    if prefix not in PREFIXES:
        return None
    fields: Dict[str, object] = {}
    for token in parts[1:]:
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        fields[key] = parse_value(value)
    return prefix, fields


def iter_records(paths: Iterable[Path]):
    for path in paths:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                parsed = parse_line(line)
                if parsed is None:
                    continue
                prefix, fields = parsed
                fields["_source"] = str(path)
                yield prefix, fields


def numeric(value) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def mean(values: List[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def collect_values(records: List[dict], key: str) -> List[float]:
    out: List[float] = []
    for record in records:
        value = numeric(record.get(key))
        if value is not None and math.isfinite(value):
            out.append(value)
    return out


def format_stat(records: List[dict], key: str) -> str:
    values = collect_values(records, key)
    if not values:
        return "n/a"
    return f"{mean(values):.6g} [{min(values):.6g}, {max(values):.6g}]"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("logs", nargs="+", help="Herwig log files to inspect")
    parser.add_argument("--q2-min", type=float, default=49.0,
                        help="Lower Q^2 cut in GeV^2 (default: 49)")
    parser.add_argument("--q2-max", type=float, default=100.0,
                        help="Upper Q^2 cut in GeV^2 (default: 100)")
    parser.add_argument("--limit", type=int, default=20,
                        help="Maximum number of low-Q2 point rows to print (default: 20)")
    parser.add_argument("--show-points", action="store_true",
                        help="Print a compact low-Q2 point table")
    args = parser.parse_args()

    paths = [Path(p) for p in args.logs]

    audit_points: Dict[Tuple[str, int], Dict[str, dict]] = defaultdict(dict)
    term_snapshots: Dict[Tuple[str, int], Dict[str, dict]] = defaultdict(dict)

    for prefix, fields in iter_records(paths):
        run = str(fields.get("run", "UNKNOWN"))
        n_value = fields.get("n")
        if not isinstance(n_value, int):
            n_float = numeric(n_value)
            if n_float is None:
                continue
            n_value = int(round(n_float))
        key = (run, n_value)

        if prefix.startswith("NLO_AUDIT_"):
            audit_points[key][prefix] = fields
        else:
            term_snapshots[key][prefix] = fields

    low_q2_keys: List[Tuple[str, int]] = []
    for key, record_set in audit_points.items():
        anchor = record_set.get("NLO_AUDIT_NC") or record_set.get("NLO_AUDIT_PDF")
        if not anchor:
            continue
        q2 = numeric(anchor.get("Q2"))
        if q2 is None:
            continue
        if args.q2_min <= q2 <= args.q2_max:
            low_q2_keys.append(key)

    print(f"Audit points found: {len(audit_points)}")
    print(f"Low-Q2 audit points ({args.q2_min:g} <= Q2 <= {args.q2_max:g} GeV^2): {len(low_q2_keys)}")
    if not low_q2_keys:
        return 0

    low_nc = [audit_points[k].get("NLO_AUDIT_NC", {}) for k in low_q2_keys if "NLO_AUDIT_NC" in audit_points[k]]
    low_pdf = [audit_points[k].get("NLO_AUDIT_PDF", {}) for k in low_q2_keys if "NLO_AUDIT_PDF" in audit_points[k]]
    low_term = [audit_points[k].get("NLO_AUDIT_TERM", {}) for k in low_q2_keys if "NLO_AUDIT_TERM" in audit_points[k]]
    low_cum = [term_snapshots[k].get("NLO_TERM_CUM", {}) for k in low_q2_keys if "NLO_TERM_CUM" in term_snapshots[k]]
    low_real = [term_snapshots[k].get("NLO_TERM_REAL", {}) for k in low_q2_keys if "NLO_TERM_REAL" in term_snapshots[k]]

    runs = Counter(k[0] for k in low_q2_keys)
    helicities = Counter(str(audit_points[k].get("NLO_AUDIT_NC", audit_points[k].get("NLO_AUDIT_PDF", {})).get("hel", "UNKNOWN"))
                         for k in low_q2_keys)
    contribs = Counter(str(audit_points[k].get("NLO_AUDIT_NC", audit_points[k].get("NLO_AUDIT_PDF", {})).get("contrib", "UNKNOWN"))
                       for k in low_q2_keys)
    channels = Counter(str(audit_points[k].get("NLO_AUDIT_NC", {}).get("channel", "UNKNOWN"))
                       for k in low_q2_keys if "NLO_AUDIT_NC" in audit_points[k])

    print("\nCounts by run:")
    for run, count in runs.items():
        print(f"  {run}: {count}")

    print("\nCounts by helicity:")
    for hel, count in helicities.items():
        print(f"  {hel}: {count}")

    print("\nCounts by contribution:")
    for contrib, count in contribs.items():
        print(f"  {contrib}: {count}")

    if channels:
        print("\nCounts by NC channel:")
        for channel, count in channels.items():
            print(f"  {channel}: {count}")

    print("\nPoint-level low-Q2 means [min, max]:")
    point_keys = [
        ("Q2", low_nc or low_pdf),
        ("xB", low_nc or low_pdf),
        ("xp", low_nc or low_pdf),
        ("y", low_nc or low_pdf),
        ("Pq", low_pdf),
        ("Pq_m", low_pdf),
        ("Pg_m", low_pdf),
        ("qUnpol", low_nc),
        ("qPol", low_nc),
        ("gUnpol", low_nc),
        ("gPol", low_nc),
        ("qOddResponse", low_nc),
        ("gOddResponse", low_nc),
        ("mappedDenomRatio", low_nc),
        ("qcdcDenRatio", low_nc),
        ("a_born", low_term),
        ("a_q_mapped", low_term),
        ("a_g_R2", low_term),
        ("a_g_R3", low_term),
        ("virt", low_term),
        ("collq_even", low_term),
        ("collq_odd", low_term),
        ("collg_even", low_term),
        ("collg_odd", low_term),
        ("realq", low_term),
        ("realg_even", low_term),
        ("realg_odd", low_term),
        ("realg", low_term),
        ("wgt", low_term),
        ("wgt_old", low_term),
        ("wgt_new", low_term),
    ]
    for key, bucket in point_keys:
        print(f"  {key}: {format_stat(bucket, key)}")

    if low_cum:
        print("\nMatching cumulative snapshot means [min, max]:")
        for key in ("F_virt", "F_cq_even", "F_cq_odd", "F_cg_even", "F_cg_odd", "F_rq", "F_rg", "F_total"):
            print(f"  {key}: {format_stat(low_cum, key)}")

    if low_real:
        print("\nMatching cumulative real-piece means [min, max]:")
        for key in ("F_rq_even", "F_rq_odd", "F_rg_even", "F_rg_odd",
                    "F_rq_chk", "F_rg_chk", "F_rq_even_nc_probe", "F_rq_spin_nc_probe"):
            print(f"  {key}: {format_stat(low_real, key)}")

    if args.show_points:
        print("\nLow-Q2 point table:")
        print("run\tn\tQ2\txB\txp\ty\tchannel\tgUnpol\tgPol\tmappedDen\tqcdcDen\twgt\trealg_even\trealg_odd")
        shown = 0
        for key in sorted(low_q2_keys):
            nc = audit_points[key].get("NLO_AUDIT_NC", {})
            term = audit_points[key].get("NLO_AUDIT_TERM", {})
            anchor = nc or audit_points[key].get("NLO_AUDIT_PDF", {})
            print(
                f"{key[0]}\t{key[1]}\t"
                f"{anchor.get('Q2', 'n/a')}\t{anchor.get('xB', 'n/a')}\t{anchor.get('xp', 'n/a')}\t{anchor.get('y', 'n/a')}\t"
                f"{nc.get('channel', 'n/a')}\t{nc.get('gUnpol', 'n/a')}\t{nc.get('gPol', 'n/a')}\t"
                f"{nc.get('mappedDenomRatio', 'n/a')}\t{nc.get('qcdcDenRatio', 'n/a')}\t"
                f"{term.get('wgt', 'n/a')}\t{term.get('realg_even', 'n/a')}\t{term.get('realg_odd', 'n/a')}"
            )
            shown += 1
            if shown >= args.limit:
                break

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
