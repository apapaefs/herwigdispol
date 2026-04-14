#!/usr/bin/env python3
"""Verify the pure-Z antiquark sign convention in Herwig and POLDIS."""

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
DEFAULT_THETA = 0.490182090162431
DEFAULT_Y = 0.4
DEFAULT_PL = 1.0
DEFAULT_PZ = 1.0
DEFAULT_DELTA_OVER_F = 0.3
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-10


@dataclass(frozen=True)
class AnalyticRow:
    qid: int
    family: str
    qorder: bool
    herwig_pdf_qid: int
    poldis_pdf_qid: int
    eta_q: float
    cv_l: float
    ca_l: float
    cv_q: float
    ca_q: float
    herwig_signed_cvca: float
    poldis_signed_cvca: float
    z_even_norm: float
    herwig_z_signed_axial_norm: float
    poldis_z_signed_axial_norm: float
    poldis_g2_pdf_norm: float
    poldis_g4_pdf_norm: float
    herwig_g4_pdf_norm: float
    poldis_total_pdf_norm: float
    herwig_pq_coeff_norm: float
    herwig_plpq_coeff_norm: float


@dataclass(frozen=True)
class TraceCheck:
    record_index: int
    qid: int
    sign: str
    family: str
    issues: list[str]


def is_close(lhs: float, rhs: float, *, abs_tol: float = ABS_TOL,
             rel_tol: float = REL_TOL) -> bool:
    return math.isclose(lhs, rhs, abs_tol=abs_tol, rel_tol=rel_tol)


def sign_label(qid: int) -> str:
    return "antiquark" if qid < 0 else "quark"


def family_label(qid: int) -> str:
    return "up" if abs(qid) % 2 == 0 else "down"


def parse_trace_line(line: str) -> dict[str, str] | None:
    marker = line.find("DIS_BOOKKEEPING_TRACE")
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


def signed_flavors(explicit: list[int] | None) -> list[int]:
    if explicit:
        return explicit
    values = []
    for abs_qid in range(1, 6):
        values.extend([abs_qid, -abs_qid])
    return values


def family_couplings(theta: float, abs_qid: int) -> tuple[float, float]:
    sw = math.sin(theta)
    cw = math.cos(theta)
    if abs_qid % 2 == 0:
        t3 = 0.5
        charge = 2.0 / 3.0
    else:
        t3 = -0.5
        charge = -1.0 / 3.0
    cv = (t3 - 2.0 * charge * sw * sw) / (2.0 * sw * cw)
    ca = t3 / (2.0 * sw * cw)
    return cv, ca


def electron_couplings(theta: float) -> tuple[float, float]:
    sw = math.sin(theta)
    cw = math.cos(theta)
    t3 = -0.5
    charge = -1.0
    cv = (t3 - 2.0 * charge * sw * sw) / (2.0 * sw * cw)
    ca = t3 / (2.0 * sw * cw)
    return cv, ca


def y_plus(y: float) -> float:
    return 1.0 + (1.0 - y) ** 2


def y_minus(y: float) -> float:
    return 1.0 - (1.0 - y) ** 2


def effective_pdf_qid(qid: int, *, flip_antiquark_pdf_sign: bool) -> int:
    if qid < 0 and flip_antiquark_pdf_sign:
        return -qid
    return qid


def effective_eta_q(qid: int, *, flip_antiquark_coupling_sign: bool) -> float:
    if qid < 0 and flip_antiquark_coupling_sign:
        return 1.0
    return -1.0 if qid < 0 else 1.0


def analytic_row(qid: int, theta: float, y: float, pl: float, pz: float,
                 delta_over_f: float, *,
                 flip_antiquark_pdf_sign: bool = False,
                 flip_antiquark_coupling_sign: bool = False) -> AnalyticRow:
    if qid == 0 or abs(qid) > 5:
        raise ValueError(f"Unsupported flavour id {qid}")
    if not (0.0 < y < 1.0):
        raise ValueError("y must satisfy 0 < y < 1")

    cv_l, ca_l = electron_couplings(theta)
    cv_q, ca_q = family_couplings(theta, abs(qid))
    eta_q = -1.0 if qid < 0 else 1.0
    eta_q_poldis = effective_eta_q(
        qid, flip_antiquark_coupling_sign=flip_antiquark_coupling_sign
    )
    qid_poldis = effective_pdf_qid(
        qid, flip_antiquark_pdf_sign=flip_antiquark_pdf_sign
    )

    even_norm = (cv_l * cv_l + ca_l * ca_l) * (cv_q * cv_q + ca_q * ca_q)
    herwig_signed_cvca = eta_q * cv_q * ca_q
    poldis_signed_cvca = eta_q_poldis * cv_q * ca_q
    herwig_z_signed_axial_norm = 4.0 * cv_l * ca_l * herwig_signed_cvca
    poldis_z_signed_axial_norm = 4.0 * cv_l * ca_l * poldis_signed_cvca

    yp = y_plus(y)
    ym = y_minus(y)
    effective_delta_over_f = (
        -delta_over_f if (qid < 0 and flip_antiquark_pdf_sign) else delta_over_f
    )
    poldis_g2_pdf_norm = pz * effective_delta_over_f * ym * even_norm
    poldis_g4_pdf_norm = pz * effective_delta_over_f * yp * poldis_z_signed_axial_norm
    herwig_g4_pdf_norm = pz * delta_over_f * yp * herwig_z_signed_axial_norm
    poldis_total_pdf_norm = poldis_g2_pdf_norm + poldis_g4_pdf_norm

    ell = 2.0 / y - 1.0
    eta_l = 1.0
    d_q = eta_q * -2.0 * ca_q * cv_q * (ca_l * ca_l + cv_l * cv_l)
    d_lq = eta_l * eta_q * 4.0 * ca_l * ca_q * cv_l * cv_q
    n_q = eta_l * -4.0 * ca_l * cv_l * (ca_q * ca_q + cv_q * cv_q)
    n_lq = 2.0 * (ca_q * ca_q + cv_q * cv_q) * (ca_l * ca_l + cv_l * cv_l)
    herwig_pq_coeff_norm = (
        (1.0 + ell * ell) * (d_q + pl * d_lq) +
        ell * (n_q + pl * n_lq)
    )
    herwig_plpq_coeff_norm = (1.0 + ell * ell) * d_lq + ell * n_lq

    return AnalyticRow(
        qid=qid,
        family=family_label(qid),
        qorder=(qid > 0),
        herwig_pdf_qid=qid,
        poldis_pdf_qid=qid_poldis,
        eta_q=eta_q,
        cv_l=cv_l,
        ca_l=ca_l,
        cv_q=cv_q,
        ca_q=ca_q,
        herwig_signed_cvca=herwig_signed_cvca,
        poldis_signed_cvca=poldis_signed_cvca,
        z_even_norm=even_norm,
        herwig_z_signed_axial_norm=herwig_z_signed_axial_norm,
        poldis_z_signed_axial_norm=poldis_z_signed_axial_norm,
        poldis_g2_pdf_norm=poldis_g2_pdf_norm,
        poldis_g4_pdf_norm=poldis_g4_pdf_norm,
        herwig_g4_pdf_norm=herwig_g4_pdf_norm,
        poldis_total_pdf_norm=poldis_total_pdf_norm,
        herwig_pq_coeff_norm=herwig_pq_coeff_norm,
        herwig_plpq_coeff_norm=herwig_plpq_coeff_norm,
    )


def analytic_issues(row: AnalyticRow) -> list[str]:
    issues: list[str] = []
    if row.poldis_pdf_qid != row.herwig_pdf_qid:
        issues.append(
            f"PDF signed flavour mismatch: POLDIS uses {row.poldis_pdf_qid}, "
            f"Herwig uses {row.herwig_pdf_qid}"
        )
    if not is_close(row.poldis_signed_cvca, row.herwig_signed_cvca):
        issues.append(
            f"signed CV*CA mismatch: POLDIS={row.poldis_signed_cvca:.16e}, "
            f"Herwig={row.herwig_signed_cvca:.16e}"
        )
    if not is_close(row.poldis_z_signed_axial_norm, row.herwig_z_signed_axial_norm):
        issues.append(
            "normalized pure-Z axial carrier mismatch: "
            f"POLDIS={row.poldis_z_signed_axial_norm:.16e}, "
            f"Herwig={row.herwig_z_signed_axial_norm:.16e}"
        )
    if not is_close(row.poldis_g4_pdf_norm, row.herwig_g4_pdf_norm):
        issues.append(
            "PDF-weighted normalized G4 mismatch: "
            f"POLDIS={row.poldis_g4_pdf_norm:.16e}, "
            f"Herwig={row.herwig_g4_pdf_norm:.16e}"
        )
    return issues


def format_table(headers: list[str], rows: Iterable[list[str]]) -> str:
    table_rows = [headers, *rows]
    widths = [max(len(str(row[idx])) for row in table_rows) for idx in range(len(headers))]
    rendered = []
    for row_index, row in enumerate(table_rows):
        cells = [str(value).ljust(widths[idx]) for idx, value in enumerate(row)]
        rendered.append(" | ".join(cells))
        if row_index == 0:
            rendered.append("-+-".join("-" * width for width in widths))
    return "\n".join(rendered)


def render_analytic(args: argparse.Namespace) -> int:
    rows = [
        analytic_row(
            qid, args.theta, args.y, args.pl, args.pz, args.delta_over_f,
            flip_antiquark_pdf_sign=args.flip_antiquark_pdf_sign,
            flip_antiquark_coupling_sign=args.flip_antiquark_coupling_sign,
        )
        for qid in signed_flavors(args.flavor)
    ]
    checks = [(row, analytic_issues(row)) for row in rows]

    headers = [
        "qid",
        "family",
        "qorder",
        "pdf_qid(P)",
        "pdf_qid(H)",
        "signed_cvca(P)",
        "signed_cvca(H)",
        "z_g4_norm(P)",
        "z_g4_norm(H)",
        "poldis_g4_pdf",
        "herwig_g4_pdf",
        "poldis_total_pdf",
        "herwig_pq_norm",
        "herwig_plpq_norm",
        "status",
    ]
    body = []
    for row, issues in checks:
        body.append([
            row.qid,
            row.family,
            int(row.qorder),
            row.poldis_pdf_qid,
            row.herwig_pdf_qid,
            f"{row.poldis_signed_cvca:.8e}",
            f"{row.herwig_signed_cvca:.8e}",
            f"{row.poldis_z_signed_axial_norm:.8e}",
            f"{row.herwig_z_signed_axial_norm:.8e}",
            f"{row.poldis_g4_pdf_norm:.8e}",
            f"{row.herwig_g4_pdf_norm:.8e}",
            f"{row.poldis_total_pdf_norm:.8e}",
            f"{row.herwig_pq_coeff_norm:.8e}",
            f"{row.herwig_plpq_coeff_norm:.8e}",
            "OK" if not issues else "FAIL",
        ])

    print(
        "# Pure-Z antiquark sign audit\n"
        f"# theta={args.theta:.15f} y={args.y:.6f} Pl={args.pl:.6f} "
        f"Pz={args.pz:.6f} delta_over_f={args.delta_over_f:.6f}\n"
        "# P = POLDIS-style signed-flavour convention, "
        "H = Herwig signed-flavour convention\n"
    )
    print(format_table(headers, body))

    failures = [(row.qid, issues) for row, issues in checks if issues]
    if failures:
        print("\nAnalytic mismatches:")
        for qid, issues in failures:
            for issue in issues:
                print(f"- qid={qid}: {issue}")
        return 1

    print(
        "\nAnalytic result: the signed antiquark PDF label and the pure-Z axial "
        "carrier agree for all requested flavours."
    )
    return 0


def check_trace(path: Path, *, flip_antiquark_pdf_sign: bool = False,
                flip_antiquark_coupling_sign: bool = False) -> tuple[int, int, list[TraceCheck]]:
    checks: list[TraceCheck] = []
    trace_lines = 0
    total = 0
    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = parse_trace_line(line)
            if not record:
                continue
            trace_lines += 1
            if record.get("nc_has") != "1":
                continue
            if record.get("nc_channel") != "Z":
                continue

            total += 1
            qid = int(record["qid"])
            sign = sign_label(qid)
            family = family_label(qid)
            issues: list[str] = []

            qorder = bool(int(record["qorder"]))
            if qorder != (qid > 0):
                issues.append(
                    f"qorder={int(qorder)} but qid={qid} implies qorder={int(qid > 0)}"
                )

            expected_pdf_qid = effective_pdf_qid(
                qid, flip_antiquark_pdf_sign=flip_antiquark_pdf_sign
            )
            current_pdf_qid = int(record["cur_pdf_qid"])
            if current_pdf_qid != expected_pdf_qid:
                issues.append(
                    f"cur_pdf_qid={current_pdf_qid} but expected {expected_pdf_qid}"
                )

            expected_eta_q = -1.0 if qid < 0 else 1.0
            actual_eta_q = float(record["nc_etaQ"])
            if not is_close(actual_eta_q, expected_eta_q):
                issues.append(
                    f"nc_etaQ={actual_eta_q:.16e} but expected {expected_eta_q:.16e}"
                )

            cv_q = float(record["nc_CVq"])
            ca_q = float(record["nc_CAq"])
            expected_signed_cvca = (
                effective_eta_q(
                    qid, flip_antiquark_coupling_sign=flip_antiquark_coupling_sign
                ) * cv_q * ca_q
            )
            actual_signed_cvca = float(record["nc_signed_cvca"])
            if not is_close(actual_signed_cvca, expected_signed_cvca):
                issues.append(
                    "nc_signed_cvca mismatch: "
                    f"actual={actual_signed_cvca:.16e}, "
                    f"expected={expected_signed_cvca:.16e}"
                )

            axial_unsigned = float(record["z_axial_coupling"])
            expected_signed_axial = (
                effective_eta_q(
                    qid, flip_antiquark_coupling_sign=flip_antiquark_coupling_sign
                ) * axial_unsigned
            )
            actual_signed_axial = float(record["z_signed_axial_coupling"])
            if not is_close(actual_signed_axial, expected_signed_axial):
                issues.append(
                    "z_signed_axial_coupling mismatch: "
                    f"actual={actual_signed_axial:.16e}, "
                    f"expected={expected_signed_axial:.16e}"
                )

            checks.append(
                TraceCheck(
                    record_index=line_number,
                    qid=qid,
                    sign=sign,
                    family=family,
                    issues=issues,
                )
            )
    return trace_lines, total, checks


def render_trace(args: argparse.Namespace) -> int:
    trace_lines, total, checks = check_trace(
        Path(args.from_trace),
        flip_antiquark_pdf_sign=args.flip_antiquark_pdf_sign,
        flip_antiquark_coupling_sign=args.flip_antiquark_coupling_sign,
    )
    if total == 0:
        if trace_lines == 0:
            print("No DIS_BOOKKEEPING_TRACE records were found.", file=sys.stderr)
        else:
            print(
                "DIS_BOOKKEEPING_TRACE records were found, but none contain the new "
                "nc_* / z_* fields. Rerun with the rebuilt audit binary.",
                file=sys.stderr,
            )
        return 1

    grouped = Counter((check.sign, check.family) for check in checks)
    bad = [check for check in checks if check.issues]

    print(f"Trace file: {args.from_trace}")
    print(f"Pure-Z trace records checked: {total}")
    print("Class counts:")
    for key, count in sorted(grouped.items()):
        print(f"- sign={key[0]} family={key[1]} n={count}")

    if bad:
        print("\nTrace mismatches:")
        for check in bad[:20]:
            for issue in check.issues:
                print(
                    f"- line={check.record_index} qid={check.qid} "
                    f"sign={check.sign} family={check.family}: {issue}"
                )
        if len(bad) > 20:
            print(f"- ... {len(bad) - 20} more mismatching records omitted")
        return 1

    print(
        "\nTrace result: every pure-Z record uses the expected signed PDF flavour, "
        "qorder, and pure-Z axial carrier sign."
    )
    return 0


def run_self_test() -> int:
    base_kwargs = dict(
        theta=DEFAULT_THETA,
        y=DEFAULT_Y,
        pl=DEFAULT_PL,
        pz=DEFAULT_PZ,
        delta_over_f=DEFAULT_DELTA_OVER_F,
    )
    rows = [analytic_row(qid, **base_kwargs) for qid in signed_flavors(None)]
    if any(analytic_issues(row) for row in rows):
        print("Self-test failed: the good analytic baseline mismatched.", file=sys.stderr)
        return 1

    bad_pdf_rows = [
        analytic_row(
            qid, **base_kwargs, flip_antiquark_pdf_sign=True
        )
        for qid in signed_flavors(None)
    ]
    if not any(analytic_issues(row) for row in bad_pdf_rows if row.qid < 0):
        print("Self-test failed: bad antiquark PDF sign did not trigger a mismatch.", file=sys.stderr)
        return 1

    bad_coupling_rows = [
        analytic_row(
            qid, **base_kwargs, flip_antiquark_coupling_sign=True
        )
        for qid in signed_flavors(None)
    ]
    if not any(analytic_issues(row) for row in bad_coupling_rows if row.qid < 0):
        print("Self-test failed: bad antiquark coupling sign did not trigger a mismatch.", file=sys.stderr)
        return 1

    print("Self-test passed.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify the pure-Z antiquark sign convention used by Herwig and POLDIS."
    )
    parser.add_argument("--theta", type=float, default=DEFAULT_THETA,
                        help="Electroweak angle theta_W in radians.")
    parser.add_argument("--pl", type=float, default=DEFAULT_PL,
                        help="Lepton longitudinal polarisation used for contextual coefficients.")
    parser.add_argument("--pz", type=float, default=DEFAULT_PZ,
                        help="Hadron longitudinal polarisation used for contextual coefficients.")
    parser.add_argument("--y", type=float, default=DEFAULT_Y,
                        help="DIS inelasticity y used for contextual coefficients.")
    parser.add_argument("--flavor", type=int, nargs="*",
                        help="Signed incoming flavour ids to audit. Default: all +-1..+-5.")
    parser.add_argument("--delta-over-f", dest="delta_over_f", type=float,
                        default=DEFAULT_DELTA_OVER_F,
                        help="Contextual polarized-PDF ratio Delta f / f.")
    parser.add_argument("--from-trace",
                        help="Path to a file containing DIS_BOOKKEEPING_TRACE lines.")
    parser.add_argument("--self-test", action="store_true",
                        help="Run built-in positive and negative controls.")
    parser.add_argument("--flip-antiquark-pdf-sign", action="store_true",
                        help="Intentional bad control: treat antiquark PDF labels/signs incorrectly.")
    parser.add_argument("--flip-antiquark-coupling-sign", action="store_true",
                        help="Intentional bad control: treat antiquark axial carriers incorrectly.")
    return parser


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    if args.self_test:
        return run_self_test()
    if args.from_trace:
        return render_trace(args)
    return render_analytic(args)


if __name__ == "__main__":
    raise SystemExit(main())
