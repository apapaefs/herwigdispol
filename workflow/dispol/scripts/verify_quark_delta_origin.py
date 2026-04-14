#!/usr/bin/env python3.10
"""
Verify that the quark delta coefficient matches exactly, and that the old
`-pi^2*CF/3` leftover was only a pointwise endpoint artifact of the split
verifier.

This script is diagnostic-only. It reuses the existing symbolic verifiers,
checks the direct Herwig-versus-POLDIS delta coefficient equality, and reports
the inherited pointwise endpoint proxy separately so it is not mistaken for a
physics mismatch.
"""

from __future__ import annotations

import argparse
import json
import math
from typing import Any, Dict, List, Sequence

import sympy as sp

import compare_polarized_nlo_terms as audit
import verify_qcdc_c4_presence as qcdc_verify
import verify_quark_nlo_split as split_verify


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--show-pieces",
        action="store_true",
        help="Print the explicit symbolic delta and artifact pieces.",
    )
    parser.add_argument(
        "--numeric-sanity",
        action="store_true",
        help="Evaluate the exact delta match and pointwise artifact with common Born weighting.",
    )
    parser.add_argument(
        "--channel",
        choices=("GAMMA", "Z", "ALL"),
        default=None,
        help="Restrict numeric sanity to one channel.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit a machine-readable summary.",
    )
    parser.add_argument("--q2", type=float, default=100.0, help="Control-point Q^2 in GeV^2.")
    parser.add_argument("--y", type=float, default=0.4, help="Control-point y value.")
    parser.add_argument("--xp", type=float, default=0.506313, help="Control-point xp value.")
    parser.add_argument(
        "--diff-pdf",
        default=audit.DEFAULT_DIFF_PDF,
        help="Polarized PDF set used for common Born-weighted sanity checks.",
    )
    parser.add_argument(
        "--sin2-theta-w",
        type=float,
        default=0.221639970483740179,
        help="sin^2(theta_W) used in the existing comparison workflow.",
    )
    parser.add_argument("--mz", type=float, default=91.1876, help="Z mass in GeV.")
    parser.add_argument("--gamma-z", type=float, default=2.4952, help="Z width in GeV.")
    parser.add_argument(
        "--s-hadronic",
        type=float,
        default=audit.DEFAULT_S,
        help="ep hadronic s in GeV^2.",
    )
    return parser.parse_args(argv)


def poldis_delta_target() -> Dict[str, Any]:
    cf, x_b = sp.symbols("CF xB", positive=True)
    logmu = sp.symbols("logmu", real=True)
    l = sp.log(1 - x_b)

    c_delta_raw = cf * (
        -2 * sp.pi**2 / 3
        - 9
        + 3 * logmu
        - 3 * l
        + 2 * l**2
        + 4 * l * logmu
    )
    c_delta_target = sp.simplify(sp.Rational(1, 2) * c_delta_raw.subs({logmu: 0}))

    return {
        "symbols": {"CF": cf, "xB": x_b, "logmu": logmu},
        "C_delta_target": c_delta_target,
    }


def herwig_virtual_delta(targets: Dict[str, Any]) -> Dict[str, Any]:
    cf = targets["symbols"]["CF"]
    x_b = targets["symbols"]["xB"]
    l = sp.log(1 - x_b)

    virt_delta_eff = sp.simplify(
        cf * (-sp.Rational(9, 2) - sp.pi**2 / 3 - sp.Rational(3, 2) * l + l**2)
    )
    direct_delta_diff = sp.simplify(virt_delta_eff - targets["C_delta_target"])

    return {
        "virt_delta_eff": virt_delta_eff,
        "direct_delta_diff": direct_delta_diff,
    }


def pointwise_endpoint_artifact(split_results: Dict[str, Any]) -> Dict[str, Any]:
    cf = split_results["targets"]["symbols"]["CF"]
    z = split_results["targets"]["symbols"]["z"]

    endpoint_kernel = sp.simplify(2 * cf * sp.log(z) / (1 - z))
    endpoint_integral = sp.simplify(-sp.pi**2 * cf / 3)
    pointwise_proxy = sp.simplify(split_results["pointwise_endpoint_proxy"]["C2_total_born_proxy_diff"])

    return {
        "symbols": {"CF": cf, "z": z},
        "endpoint_kernel": endpoint_kernel,
        "endpoint_integral": endpoint_integral,
        "pointwise_proxy": pointwise_proxy,
        "proxy_minus_integral": sp.simplify(pointwise_proxy - endpoint_integral),
    }


def evaluate_symbolic_result() -> Dict[str, Any]:
    split_results = split_verify.evaluate_symbolic_basis()
    targets = poldis_delta_target()
    virt_info = herwig_virtual_delta(targets)
    artifact = pointwise_endpoint_artifact(split_results)

    direct_ok = sp.simplify(virt_info["direct_delta_diff"]) == 0
    artifact_ok = sp.simplify(artifact["proxy_minus_integral"]) == 0

    if direct_ok and artifact_ok:
        verdict = "exact_delta_match_with_pointwise_artifact"
    elif direct_ok:
        verdict = "exact_delta_match_with_different_proxy"
    else:
        verdict = "different_delta_formula"

    return {
        "targets": targets,
        "virt_info": virt_info,
        "artifact": artifact,
        "verdict": verdict,
    }


def numeric_sanity(args: argparse.Namespace, symbolic: Dict[str, Any]) -> List[Dict[str, Any]]:
    channels = [args.channel] if args.channel else ["GAMMA", "Z", "ALL"]
    lhapdf = audit.import_lhapdf()
    diff_pdf = lhapdf.mkPDF(args.diff_pdf, 0)
    ew_inputs = audit.EWInputs(
        sin2_theta_w=args.sin2_theta_w,
        mz=args.mz,
        gamma_z=args.gamma_z,
        lepton_charge=-1.0,
    )

    x_b = args.q2 / (args.y * args.s_hadronic)
    if not (0.0 < x_b < 1.0):
        raise SystemExit(f"x_B={x_b:.6f} lies outside (0,1) for the chosen numeric sanity point.")

    poldis_alpha_s = diff_pdf.alphasQ2(args.q2)
    alpha_prefactor = poldis_alpha_s / (4.0 * math.pi)

    cf = 4.0 / 3.0
    direct_coeff = float(
        symbolic["virt_info"]["direct_delta_diff"].subs(
            {symbolic["targets"]["symbols"]["CF"]: cf, symbolic["targets"]["symbols"]["xB"]: x_b}
        ).evalf()
    )
    proxy_coeff = float(
        symbolic["artifact"]["pointwise_proxy"].subs(
            {symbolic["artifact"]["symbols"]["CF"]: cf}
        ).evalf()
    )

    rows: List[Dict[str, Any]] = []
    for channel in channels:
        for flavor in audit.FLAVORS:
            poldis_terms = audit.poldis_coefficients_single_flavor(
                args.q2,
                x_b,
                args.y,
                args.xp,
                channel,
                flavor,
                ew_inputs,
                diff_pdf,
                poldis_alpha_s,
            )
            common_abs_factor = alpha_prefactor * poldis_terms.born_abs

            rows.append(
                {
                    "channel": channel,
                    "flavor": flavor.label,
                    "directDeltaMismatchAbs": common_abs_factor * direct_coeff,
                    "pointwiseProxyAbs": common_abs_factor * proxy_coeff,
                }
            )
    return rows


def print_text_report(symbolic: Dict[str, Any], rows: List[Dict[str, Any]] | None, show_pieces: bool) -> None:
    print("Analytic quark-delta verifier")
    print("  caveat: direct delta agreement is tested here")
    print("  caveat: the old -pi^2*CF/3 residual is reported only as a pointwise endpoint artifact")

    print("\nVerdict")
    print(f"  verdict = {symbolic['verdict']}")

    print("\nDirect delta comparison")
    print(f"  C_delta_target    = {qcdc_verify.expr_string(symbolic['targets']['C_delta_target'])}")
    print(f"  virt_delta_eff    = {qcdc_verify.expr_string(symbolic['virt_info']['virt_delta_eff'])}")
    print(f"  direct_delta_diff = {qcdc_verify.expr_string(symbolic['virt_info']['direct_delta_diff'])}")

    print("\nPointwise endpoint artifact")
    print(f"  endpoint_kernel   = {qcdc_verify.expr_string(symbolic['artifact']['endpoint_kernel'])}")
    print(f"  endpoint_integral = {qcdc_verify.expr_string(symbolic['artifact']['endpoint_integral'])}")
    print(f"  pointwise_proxy   = {qcdc_verify.expr_string(symbolic['artifact']['pointwise_proxy'])}")
    print(f"  proxy-minus-int   = {qcdc_verify.expr_string(symbolic['artifact']['proxy_minus_integral'])}")

    if show_pieces:
        print("\nInterpretation")
        print("  The split verifier keeps the virtual term as born_abs * vdelta in a local-z algebra.")
        print("  That reproduces the z-dependent structure correctly but cannot be interpreted as a")
        print("  genuine delta(1-z) extraction. The inherited -pi^2*CF/3 is therefore a proxy artifact.")

    if rows is not None:
        print("\nNumeric sanity")
        print("  channel  flavor  directDeltaMismatchAbs  pointwiseProxyAbs")
        print("  -------  ------  ----------------------  -----------------")
        for row in rows:
            print(
                f"  {row['channel']:<7}  {row['flavor']:<6}  "
                f"{row['directDeltaMismatchAbs']:<22.8e}  {row['pointwiseProxyAbs']:<17.8e}"
            )


def build_json_report(symbolic: Dict[str, Any], rows: List[Dict[str, Any]] | None) -> Dict[str, Any]:
    out = {
        "verdict": symbolic["verdict"],
        "c_delta_target": qcdc_verify.expr_string(symbolic["targets"]["C_delta_target"]),
        "virt_delta_eff": qcdc_verify.expr_string(symbolic["virt_info"]["virt_delta_eff"]),
        "direct_delta_diff": qcdc_verify.expr_string(symbolic["virt_info"]["direct_delta_diff"]),
        "endpoint_kernel": qcdc_verify.expr_string(symbolic["artifact"]["endpoint_kernel"]),
        "endpoint_integral": qcdc_verify.expr_string(symbolic["artifact"]["endpoint_integral"]),
        "pointwise_proxy": qcdc_verify.expr_string(symbolic["artifact"]["pointwise_proxy"]),
        "proxy_minus_integral": qcdc_verify.expr_string(symbolic["artifact"]["proxy_minus_integral"]),
    }
    if rows is not None:
        out["numeric_sanity"] = rows
    return out


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    symbolic = evaluate_symbolic_result()
    rows = numeric_sanity(args, symbolic) if args.numeric_sanity else None

    if args.json:
        print(json.dumps(build_json_report(symbolic, rows), indent=2, sort_keys=True))
    else:
        print_text_report(symbolic, rows, show_pieces=args.show_pieces)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
