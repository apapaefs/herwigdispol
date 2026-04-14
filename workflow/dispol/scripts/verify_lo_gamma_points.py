#!/usr/bin/env python3.10
"""Verify sampled LO_GAMMA_POINT diagnostics against the standalone LO calculator."""

from __future__ import annotations

import argparse
import math
from pathlib import Path
from typing import Dict, Iterable, List

from lo_dis_nc_common import (
    BeamInputs,
    DEFAULT_ALPHA_EM,
    DEFAULT_BEAM_EA_GEV,
    DEFAULT_BEAM_EB_GEV,
    DEFAULT_DIFF_PDF,
    DEFAULT_GAMMA_Z,
    DEFAULT_LEPTON_CHARGE,
    DEFAULT_MODE_HERWIG_LO,
    DEFAULT_MW,
    DEFAULT_MZ,
    DEFAULT_SIN2_THETA_W,
    DEFAULT_SUM_PDF,
    EWInputs,
    import_lhapdf,
    plain_lo_gamma_hwmebase_jacobian,
    point_audit_observables,
)


NUMERIC_FIELDS = {
    "sample",
    "qid",
    "x1",
    "x2",
    "xB",
    "Q2",
    "y",
    "mu2",
    "sHat",
    "tHat",
    "uHat",
    "Pl",
    "Pq",
    "pdf_sum",
    "me2",
    "jacobian",
    "sigma_hat_nb",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--from-log", required=True, help="Path to a Herwig .log file containing LO_GAMMA_POINT lines.")
    parser.add_argument("--sum-pdf", default=DEFAULT_SUM_PDF, help="Unpolarized LHAPDF set.")
    parser.add_argument("--diff-pdf", default=DEFAULT_DIFF_PDF, help="Polarized-difference LHAPDF set.")
    parser.add_argument("--beam-ea", type=float, default=DEFAULT_BEAM_EA_GEV, help="Lepton beam energy in GeV.")
    parser.add_argument("--beam-eb", type=float, default=DEFAULT_BEAM_EB_GEV, help="Hadron beam energy in GeV.")
    parser.add_argument("--alpha-em", type=float, default=DEFAULT_ALPHA_EM, help="Fixed alpha_EM input.")
    parser.add_argument("--sin2-theta-w", type=float, default=DEFAULT_SIN2_THETA_W, help="sin^2(theta_W).")
    parser.add_argument("--mw", type=float, default=DEFAULT_MW, help="W mass in GeV.")
    parser.add_argument("--mz", type=float, default=DEFAULT_MZ, help="Z mass in GeV.")
    parser.add_argument("--gamma-z", type=float, default=DEFAULT_GAMMA_Z, help="Z width in GeV.")
    parser.add_argument("--lepton-charge", type=float, default=DEFAULT_LEPTON_CHARGE, help="Incoming lepton charge.")
    parser.add_argument("--pdf-rtol", type=float, default=1.0e-10, help="Relative tolerance for pdf_sum.")
    parser.add_argument("--pdf-atol", type=float, default=1.0e-12, help="Absolute tolerance for pdf_sum.")
    parser.add_argument("--jac-rtol", type=float, default=1.0e-10, help="Relative tolerance for the plain LO Jacobian.")
    parser.add_argument("--jac-atol", type=float, default=1.0e-12, help="Absolute tolerance for the plain LO Jacobian.")
    parser.add_argument("--sigma-rtol", type=float, default=1.0e-10, help="Relative tolerance for sigma_hat_nb.")
    parser.add_argument("--sigma-atol", type=float, default=1.0e-12, help="Absolute tolerance for sigma_hat_nb.")
    parser.add_argument(
        "--min-q2-cut",
        type=float,
        default=49.0,
        help="Plain-run SimpleDISCut MinQ2 in GeV^2 used in the HwMEBase angular map.",
    )
    return parser.parse_args()


def parse_lo_gamma_point_line(line: str) -> Dict[str, object] | None:
    if "LO_GAMMA_POINT" not in line:
        return None
    payload = line.split("LO_GAMMA_POINT", 1)[1].strip()
    fields: Dict[str, object] = {}
    for token in payload.split():
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        if key in NUMERIC_FIELDS:
            numeric = float(value)
            if key in {"sample", "qid"}:
                fields[key] = int(round(numeric))
            else:
                fields[key] = numeric
        else:
            fields[key] = value
    return fields


def load_points(path: Path) -> List[Dict[str, object]]:
    points: List[Dict[str, object]] = []
    for line in path.read_text().splitlines():
        parsed = parse_lo_gamma_point_line(line)
        if parsed is not None:
            points.append(parsed)
    return points


def within_tolerance(expected: float, observed: float, atol: float, rtol: float) -> bool:
    diff = abs(expected - observed)
    scale = max(abs(expected), abs(observed), 1.0)
    return diff <= atol or diff <= rtol * scale


def summarize_delta(expected: float, observed: float) -> str:
    diff = observed - expected
    rel = diff / expected if abs(expected) > 0.0 else float("inf")
    return f"expected={expected:.16e} observed={observed:.16e} diff={diff:.16e} rel={rel:.16e}"


def iter_failures(
    points: Iterable[Dict[str, object]],
    ew_inputs: EWInputs,
    beams: BeamInputs,
    sum_pdf,
    diff_pdf,
    pdf_atol: float,
    pdf_rtol: float,
    jac_atol: float,
    jac_rtol: float,
    sigma_atol: float,
    sigma_rtol: float,
    min_q2_cut: float,
):
    for point in points:
        expected_jacobian = plain_lo_gamma_hwmebase_jacobian(
            q2=float(point["Q2"]),
            s_hat=float(point["sHat"]),
            q2_min_cut=min_q2_cut,
        )
        observed_jac = float(point["jacobian"])
        expected_jac = float(expected_jacobian["jacobian"])
        if not within_tolerance(expected_jac, observed_jac, jac_atol, jac_rtol):
            yield {
                "kind": "jacobian",
                "sample": int(point["sample"]),
                "qid": int(point["qid"]),
                "details": summarize_delta(expected_jac, observed_jac),
            }
            continue

        expected = point_audit_observables(
            q2=float(point["Q2"]),
            y=float(point["y"]),
            qid=int(point["qid"]),
            ew_inputs=ew_inputs,
            beams=beams,
            sum_pdf=sum_pdf,
            diff_pdf=diff_pdf,
            channel="GAMMA",
            mode=DEFAULT_MODE_HERWIG_LO,
            pdf_q2=float(point["mu2"]),
            jacobian=expected_jac,
            s_hat=float(point["sHat"]),
            x_b=float(point["xB"]),
        )

        observed_pdf = float(point["pdf_sum"])
        expected_pdf = float(expected["pdf_sum"])
        if not within_tolerance(expected_pdf, observed_pdf, pdf_atol, pdf_rtol):
            yield {
                "kind": "pdf_sum",
                "sample": int(point["sample"]),
                "qid": int(point["qid"]),
                "details": summarize_delta(expected_pdf, observed_pdf),
            }
            continue

        observed_sigma = float(point["sigma_hat_nb"])
        expected_sigma = float(expected["sigma_hat_nb"])
        if not within_tolerance(expected_sigma, observed_sigma, sigma_atol, sigma_rtol):
            yield {
                "kind": "sigma_hat_nb",
                "sample": int(point["sample"]),
                "qid": int(point["qid"]),
                "details": summarize_delta(expected_sigma, observed_sigma),
            }


def main() -> int:
    args = parse_args()
    log_path = Path(args.from_log)
    if not log_path.exists():
        raise SystemExit(f"Log file not found: {log_path}")

    points = load_points(log_path)
    if not points:
        raise SystemExit(f"No LO_GAMMA_POINT lines found in {log_path}")

    ew_inputs = EWInputs(
        sin2_theta_w=args.sin2_theta_w,
        alpha_em=args.alpha_em,
        mw=args.mw,
        mz=args.mz,
        gamma_z=args.gamma_z,
        lepton_charge=args.lepton_charge,
    )
    beams = BeamInputs(args.beam_ea, args.beam_eb)
    lhapdf = import_lhapdf()
    sum_pdf = lhapdf.mkPDF(args.sum_pdf, 0)
    diff_pdf = lhapdf.mkPDF(args.diff_pdf, 0)

    failures = list(
        iter_failures(
            points,
            ew_inputs,
            beams,
            sum_pdf,
            diff_pdf,
            args.pdf_atol,
            args.pdf_rtol,
            args.jac_atol,
            args.jac_rtol,
            args.sigma_atol,
            args.sigma_rtol,
            args.min_q2_cut,
        )
    )

    print(f"Checked {len(points)} LO_GAMMA_POINT records from {log_path}.")
    if failures:
        first = failures[0]
        print(
            f"FAIL sample={first['sample']} qid={first['qid']} kind={first['kind']} {first['details']}"
        )
        print(f"Total failing records: {len(failures)}")
        return 1

    print("PASS pdf_sum, jacobian, and sigma_hat_nb agree with the standalone LO calculator.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
