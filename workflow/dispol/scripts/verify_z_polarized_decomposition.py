#!/usr/bin/env python3
"""Verify the local pure-Z polarized coefficient decomposition."""

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
DEFAULT_DELTA_OVER_F = 0.3
ABS_TOL = 1.0e-12
REL_TOL = 1.0e-10


@dataclass(frozen=True)
class AnalyticRow:
    qid: int
    family: str
    qorder: bool
    z_poldis_g2_coeff: float
    z_poldis_g4_coeff: float
    z_poldis_pq_coeff: float
    z_poldis_plpq_coeff: float
    z_pq_coeff: float
    z_plpq_coeff: float
    z_pq_residual: float
    z_plpq_residual: float
    weighted_poldis_pq: float
    weighted_herwig_pq: float


@dataclass(frozen=True)
class TraceCheck:
    record_index: int
    source: str
    qid: int
    sign: str
    family: str
    issues: list[str]


def is_close(lhs: float, rhs: float, *, abs_tol: float = ABS_TOL,
             rel_tol: float = REL_TOL) -> bool:
    return math.isclose(lhs, rhs, abs_tol=abs_tol, rel_tol=rel_tol)


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
    values: list[int] = []
    for abs_qid in range(1, 6):
        values.extend([abs_qid, -abs_qid])
    return values


def sign_label(qid: int) -> str:
    return "antiquark" if qid < 0 else "quark"


def family_label(qid: int) -> str:
    return "up" if abs(qid) % 2 == 0 else "down"


def electron_couplings(theta: float) -> tuple[float, float]:
    sw = math.sin(theta)
    cw = math.cos(theta)
    t3 = -0.5
    charge = -1.0
    cv = (t3 - 2.0 * charge * sw * sw) / (2.0 * sw * cw)
    ca = t3 / (2.0 * sw * cw)
    return cv, ca


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


def eta_q(qid: int) -> float:
    return -1.0 if qid < 0 else 1.0


def y_plus(y: float) -> float:
    return 1.0 + (1.0 - y) ** 2


def y_minus(y: float) -> float:
    return 1.0 - (1.0 - y) ** 2


def z_propagator_prefactor(*, qid: int, cv_l: float, ca_l: float,
                           cv_q: float, ca_q: float, pure_z_even: float) -> float:
    lepton_even = cv_l * cv_l + ca_l * ca_l
    quark_even = cv_q * cv_q + ca_q * ca_q
    denom = lepton_even * quark_even
    if abs(denom) <= 1.0e-30:
        raise ValueError(f"Degenerate coupling denominator for qid={qid}")
    return pure_z_even / denom


def poldis_like_coefficients(*, z_sq: float, cv_l: float, ca_l: float,
                             cv_q: float, ca_q: float, eta_l: float,
                             eta_q_value: float, pl: float, y: float,
                             swap_y_factors: bool = False,
                             flip_signed_axial: bool = False) -> tuple[float, float, float, float]:
    lepton_even = cv_l * cv_l + ca_l * ca_l
    lepton_axial = cv_l * ca_l
    quark_even = cv_q * cv_q + ca_q * ca_q
    signed_axial = eta_q_value * cv_q * ca_q
    if flip_signed_axial:
        signed_axial *= -1.0

    inv_y2 = 1.0 / (y * y)
    z_g2 = (
        2.0 * z_sq * quark_even *
        (pl * lepton_even - 2.0 * eta_l * lepton_axial) *
        inv_y2
    )
    z_g4 = (
        -4.0 * z_sq * signed_axial *
        (lepton_even - 2.0 * eta_l * pl * lepton_axial) *
        inv_y2
    )

    yp = y_plus(y)
    ym = y_minus(y)
    if swap_y_factors:
        yp, ym = ym, yp

    z_pq = ym * z_g2 + yp * z_g4
    z_plpq = (
        ym * (2.0 * z_sq * quark_even * lepton_even * inv_y2) +
        yp * (8.0 * z_sq * eta_l * signed_axial * lepton_axial * inv_y2)
    )
    return z_g2, z_g4, z_pq, z_plpq


def herwig_coefficients(*, z_sq: float, cv_l: float, ca_l: float,
                        cv_q: float, ca_q: float, eta_l: float,
                        eta_q_value: float, pl: float, y: float) -> tuple[float, float]:
    ell = 2.0 / y - 1.0
    d_q = eta_q_value * -2.0 * ca_q * cv_q * (ca_l * ca_l + cv_l * cv_l) * z_sq
    d_lq = eta_l * eta_q_value * 4.0 * ca_l * ca_q * cv_l * cv_q * z_sq
    n_q = eta_l * -4.0 * ca_l * cv_l * (ca_q * ca_q + cv_q * cv_q) * z_sq
    n_lq = 2.0 * (ca_q * ca_q + cv_q * cv_q) * (ca_l * ca_l + cv_l * cv_l) * z_sq

    z_pq = (1.0 + ell * ell) * (d_q + pl * d_lq) + ell * (n_q + pl * n_lq)
    z_plpq = (1.0 + ell * ell) * d_lq + ell * n_lq
    return z_pq, z_plpq


def analytic_row(qid: int, theta: float, pl: float, y: float, delta_over_f: float, *,
                 swap_y_factors: bool = False,
                 flip_signed_axial: bool = False) -> AnalyticRow:
    if qid == 0 or abs(qid) > 5:
        raise ValueError(f"Unsupported flavour id {qid}")
    if not (0.0 < y < 1.0):
        raise ValueError("y must satisfy 0 < y < 1")

    cv_l, ca_l = electron_couplings(theta)
    cv_q, ca_q = family_couplings(theta, abs(qid))
    eta_l = 1.0
    eta_q_value = eta_q(qid)
    sw = math.sin(theta)
    cw = math.cos(theta)
    k = 1.0 / (sw * sw * cw * cw)
    z_sq = k * k

    z_poldis_g2, z_poldis_g4, z_poldis_pq, z_poldis_plpq = poldis_like_coefficients(
        z_sq=z_sq,
        cv_l=cv_l,
        ca_l=ca_l,
        cv_q=cv_q,
        ca_q=ca_q,
        eta_l=eta_l,
        eta_q_value=eta_q_value,
        pl=pl,
        y=y,
        swap_y_factors=swap_y_factors,
        flip_signed_axial=flip_signed_axial,
    )
    z_pq, z_plpq = herwig_coefficients(
        z_sq=z_sq,
        cv_l=cv_l,
        ca_l=ca_l,
        cv_q=cv_q,
        ca_q=ca_q,
        eta_l=eta_l,
        eta_q_value=eta_q_value,
        pl=pl,
        y=y,
    )

    return AnalyticRow(
        qid=qid,
        family=family_label(qid),
        qorder=(qid > 0),
        z_poldis_g2_coeff=z_poldis_g2,
        z_poldis_g4_coeff=z_poldis_g4,
        z_poldis_pq_coeff=z_poldis_pq,
        z_poldis_plpq_coeff=z_poldis_plpq,
        z_pq_coeff=z_pq,
        z_plpq_coeff=z_plpq,
        z_pq_residual=z_pq - z_poldis_pq,
        z_plpq_residual=z_plpq - z_poldis_plpq,
        weighted_poldis_pq=delta_over_f * z_poldis_pq,
        weighted_herwig_pq=delta_over_f * z_pq,
    )


def analytic_issues(row: AnalyticRow) -> list[str]:
    issues: list[str] = []
    if not is_close(row.z_poldis_pq_coeff, row.z_pq_coeff):
        issues.append(
            f"z_pq mismatch: POLDIS={row.z_poldis_pq_coeff:.16e}, "
            f"Herwig={row.z_pq_coeff:.16e}"
        )
    if not is_close(row.z_poldis_plpq_coeff, row.z_plpq_coeff):
        issues.append(
            f"z_plpq mismatch: POLDIS={row.z_poldis_plpq_coeff:.16e}, "
            f"Herwig={row.z_plpq_coeff:.16e}"
        )
    if not is_close(row.z_pq_residual, 0.0):
        issues.append(f"non-zero z_pq residual: {row.z_pq_residual:.16e}")
    if not is_close(row.z_plpq_residual, 0.0):
        issues.append(f"non-zero z_plpq residual: {row.z_plpq_residual:.16e}")
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
            qid,
            args.theta,
            args.pl,
            args.y,
            args.delta_over_f,
            swap_y_factors=args.swap_y_factors,
            flip_signed_axial=args.flip_signed_axial,
        )
        for qid in signed_flavors(args.flavor)
    ]
    checks = [(row, analytic_issues(row)) for row in rows]

    headers = [
        "qid",
        "family",
        "qorder",
        "z_poldis_g2",
        "z_poldis_g4",
        "z_poldis_pq",
        "z_poldis_plpq",
        "z_pq",
        "z_plpq",
        "weighted_poldis",
        "weighted_herwig",
        "status",
    ]
    body = []
    for row, issues in checks:
        body.append([
            row.qid,
            row.family,
            int(row.qorder),
            f"{row.z_poldis_g2_coeff:.8e}",
            f"{row.z_poldis_g4_coeff:.8e}",
            f"{row.z_poldis_pq_coeff:.8e}",
            f"{row.z_poldis_plpq_coeff:.8e}",
            f"{row.z_pq_coeff:.8e}",
            f"{row.z_plpq_coeff:.8e}",
            f"{row.weighted_poldis_pq:.8e}",
            f"{row.weighted_herwig_pq:.8e}",
            "OK" if not issues else "FAIL",
        ])

    print(
        "# Pure-Z polarized decomposition audit\n"
        f"# theta={args.theta:.15f} y={args.y:.6f} Pl={args.pl:.6f} "
        f"delta_over_f={args.delta_over_f:.6f}\n"
        "# analytic mode uses the common electroweak pure-Z prefactor only; "
        "trace mode keeps the full event-local zSq\n"
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
        "\nAnalytic result: the POLDIS-like Y_- G2 + Y_+ G4 reconstruction "
        "matches Herwig's local pure-Z coefficients for all requested flavours."
    )
    return 0


def check_trace(path: Path, *, source_filter: str) -> tuple[int, int, list[TraceCheck]]:
    checks: list[TraceCheck] = []
    trace_lines = 0
    total = 0
    required = {
        "source",
        "qid",
        "y",
        "rho_Pl_me",
        "nc_etaL",
        "nc_etaQ",
        "nc_CVl",
        "nc_CAl",
        "nc_CVq",
        "nc_CAq",
        "z_even_coupling",
        "z_pq_coeff",
        "z_plpq_coeff",
        "z_poldis_g2_coeff",
        "z_poldis_g4_coeff",
        "z_poldis_pq_coeff",
        "z_poldis_plpq_coeff",
        "z_pq_residual",
        "z_plpq_residual",
    }

    with path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = parse_trace_line(line)
            if not record:
                continue
            trace_lines += 1
            if record.get("nc_has") != "1" or record.get("nc_channel") != "Z":
                continue
            if source_filter != "all" and record.get("source") != source_filter:
                continue
            if not required.issubset(record.keys()):
                continue

            total += 1
            qid = int(record["qid"])
            source = record["source"]
            sign = sign_label(qid)
            family = family_label(qid)
            issues: list[str] = []

            y = float(record["y"])
            pl = float(record["rho_Pl_me"])
            eta_l_value = float(record["nc_etaL"])
            eta_q_value = float(record["nc_etaQ"])
            cv_l = float(record["nc_CVl"])
            ca_l = float(record["nc_CAl"])
            cv_q = float(record["nc_CVq"])
            ca_q = float(record["nc_CAq"])
            pure_z_even = float(record["z_even_coupling"])
            z_sq = z_propagator_prefactor(
                qid=qid,
                cv_l=cv_l,
                ca_l=ca_l,
                cv_q=cv_q,
                ca_q=ca_q,
                pure_z_even=pure_z_even,
            )

            exp_g2, exp_g4, exp_pq, exp_plpq = poldis_like_coefficients(
                z_sq=z_sq,
                cv_l=cv_l,
                ca_l=ca_l,
                cv_q=cv_q,
                ca_q=ca_q,
                eta_l=eta_l_value,
                eta_q_value=eta_q_value,
                pl=pl,
                y=y,
            )
            act_g2 = float(record["z_poldis_g2_coeff"])
            act_g4 = float(record["z_poldis_g4_coeff"])
            act_pq = float(record["z_poldis_pq_coeff"])
            act_plpq = float(record["z_poldis_plpq_coeff"])
            herwig_pq = float(record["z_pq_coeff"])
            herwig_plpq = float(record["z_plpq_coeff"])
            pq_residual = float(record["z_pq_residual"])
            plpq_residual = float(record["z_plpq_residual"])

            if not is_close(act_g2, exp_g2):
                issues.append(
                    f"z_poldis_g2_coeff mismatch: actual={act_g2:.16e}, expected={exp_g2:.16e}"
                )
            if not is_close(act_g4, exp_g4):
                issues.append(
                    f"z_poldis_g4_coeff mismatch: actual={act_g4:.16e}, expected={exp_g4:.16e}"
                )
            if not is_close(act_pq, exp_pq):
                issues.append(
                    f"z_poldis_pq_coeff mismatch: actual={act_pq:.16e}, expected={exp_pq:.16e}"
                )
            if not is_close(act_plpq, exp_plpq):
                issues.append(
                    f"z_poldis_plpq_coeff mismatch: actual={act_plpq:.16e}, expected={exp_plpq:.16e}"
                )
            if not is_close(herwig_pq, act_pq):
                issues.append(
                    f"z_pq_coeff mismatch: Herwig={herwig_pq:.16e}, POLDIS={act_pq:.16e}"
                )
            if not is_close(herwig_plpq, act_plpq):
                issues.append(
                    f"z_plpq_coeff mismatch: Herwig={herwig_plpq:.16e}, POLDIS={act_plpq:.16e}"
                )
            if not is_close(pq_residual, herwig_pq - act_pq):
                issues.append(
                    f"z_pq_residual mismatch: actual={pq_residual:.16e}, "
                    f"expected={(herwig_pq - act_pq):.16e}"
                )
            if not is_close(plpq_residual, herwig_plpq - act_plpq):
                issues.append(
                    f"z_plpq_residual mismatch: actual={plpq_residual:.16e}, "
                    f"expected={(herwig_plpq - act_plpq):.16e}"
                )

            checks.append(
                TraceCheck(
                    record_index=line_number,
                    source=source,
                    qid=qid,
                    sign=sign,
                    family=family,
                    issues=issues,
                )
            )
    return trace_lines, total, checks


def render_trace(args: argparse.Namespace) -> int:
    trace_path = Path(args.from_trace)
    if not trace_path.exists():
        print(f"Trace file not found: {trace_path}", file=sys.stderr)
        return 1
    trace_lines, total, checks = check_trace(
        trace_path,
        source_filter=args.source_filter,
    )
    if total == 0:
        if trace_lines == 0:
            print("No DIS_BOOKKEEPING_TRACE records were found.", file=sys.stderr)
        else:
            print(
                "DIS_BOOKKEEPING_TRACE records were found, but none contain the new "
                "pure-Z decomposition fields. Rerun with the rebuilt audit binary.",
                file=sys.stderr,
            )
        return 1

    grouped = Counter((check.source, check.sign, check.family) for check in checks)
    bad = [check for check in checks if check.issues]

    print(f"Trace file: {args.from_trace}")
    print(f"Pure-Z trace records checked: {total}")
    print("Class counts:")
    for key, count in sorted(grouped.items()):
        print(f"- source={key[0]} sign={key[1]} family={key[2]} n={count}")

    if bad:
        print("\nTrace mismatches:")
        for check in bad[:20]:
            for issue in check.issues:
                print(
                    f"- line={check.record_index} source={check.source} qid={check.qid} "
                    f"sign={check.sign} family={check.family}: {issue}"
                )
        if len(bad) > 20:
            print(f"- ... {len(bad) - 20} more mismatching records omitted")
        return 1

    print(
        "\nTrace result: every pure-Z record satisfies the local "
        "Y_- G2 + Y_+ G4 decomposition and the logged residuals vanish."
    )
    return 0


def run_self_test() -> int:
    base_kwargs = dict(
        theta=DEFAULT_THETA,
        pl=DEFAULT_PL,
        y=DEFAULT_Y,
        delta_over_f=DEFAULT_DELTA_OVER_F,
    )

    good_rows = [analytic_row(qid, **base_kwargs) for qid in signed_flavors(None)]
    if any(analytic_issues(row) for row in good_rows):
        print("Self-test failed: the good analytic baseline mismatched.", file=sys.stderr)
        return 1

    swapped_rows = [
        analytic_row(qid, **base_kwargs, swap_y_factors=True)
        for qid in signed_flavors(None)
    ]
    if not any(analytic_issues(row) for row in swapped_rows):
        print("Self-test failed: swapping Y_+ and Y_- did not trigger a mismatch.", file=sys.stderr)
        return 1

    flipped_rows = [
        analytic_row(qid, **base_kwargs, flip_signed_axial=True)
        for qid in signed_flavors(None)
    ]
    if not any(analytic_issues(row) for row in flipped_rows):
        print("Self-test failed: flipping the signed axial carrier did not trigger a mismatch.", file=sys.stderr)
        return 1

    print("Self-test passed.")
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Verify the local pure-Z polarized coefficient decomposition."
    )
    parser.add_argument("--theta", type=float, default=DEFAULT_THETA,
                        help="Electroweak angle theta_W in radians.")
    parser.add_argument("--pl", type=float, default=DEFAULT_PL,
                        help="Lepton longitudinal polarisation.")
    parser.add_argument("--y", type=float, default=DEFAULT_Y,
                        help="DIS inelasticity y.")
    parser.add_argument("--delta-over-f", dest="delta_over_f", type=float,
                        default=DEFAULT_DELTA_OVER_F,
                        help="Contextual polarized-PDF ratio Delta f / f.")
    parser.add_argument("--flavor", type=int, nargs="*",
                        help="Signed incoming flavour ids to audit. Default: all +-1..+-5.")
    parser.add_argument("--from-trace",
                        help="Path to a file containing DIS_BOOKKEEPING_TRACE lines.")
    parser.add_argument("--source-filter", choices=["me2", "NLOWeight", "all"],
                        default="all",
                        help="Restrict trace checks to one logged source.")
    parser.add_argument("--self-test", action="store_true",
                        help="Run built-in positive and negative controls.")
    parser.add_argument("--swap-y-factors", action="store_true",
                        help="Intentional bad control: swap Y_+ and Y_- in the POLDIS reconstruction.")
    parser.add_argument("--flip-signed-axial", action="store_true",
                        help="Intentional bad control: flip the signed axial carrier in the POLDIS reconstruction.")
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
