#!/usr/bin/env python3.10
"""
Show how the inherited pointwise endpoint proxy behaves under a shared
canonical endpoint rewrite.

This script is diagnostic-only. It does not test physical delta(1-z)
agreement directly. Instead, it demonstrates that when the same canonical

  regular + subtracted-convolution + delta

rewrite is applied to both Herwig and POLDIS pointwise coefficient proxies,
the same inherited endpoint artifact remains unchanged. That is expected,
because the underlying local-z split basis does not represent the virtual term
as a genuine delta distribution.
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
        "--show-basis",
        action="store_true",
        help="Print the explicit pre-rewrite and canonical pointwise-proxy pieces.",
    )
    parser.add_argument(
        "--numeric-sanity",
        action="store_true",
        help="Evaluate the inherited pointwise proxy with common Born weighting at the control point.",
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


def load_corrected_split_results() -> Dict[str, Any]:
    results = split_verify.evaluate_symbolic_basis()
    if results["verdicts"]["total_quark_verdict"] != "matched_local_z_structure":
        raise RuntimeError("Expected matched local-z quark structure before endpoint-proxy rewrite.")
    return results


def build_canonical_convolution_basis(split_results: Dict[str, Any]) -> Dict[str, Any]:
    z = split_results["targets"]["symbols"]["z"]
    cf = split_results["targets"]["symbols"]["CF"]

    h = sp.simplify(2 * cf * sp.log(z) / (1 - z))
    h_integral = sp.simplify(-sp.pi**2 * cf / 3)

    pre_proxy = {
        "C2_pre_pointwise_delta_proxy": sp.simplify(
            split_results["pointwise_endpoint_proxy"]["C2_total_born_proxy_diff"]
        ),
        "C4_pre_pointwise_delta_proxy": sp.simplify(
            split_results["pointwise_endpoint_proxy"]["C4_total_born_proxy_diff"]
        ),
    }
    post_proxy = {
        "C2_post_pointwise_delta_proxy": sp.simplify(pre_proxy["C2_pre_pointwise_delta_proxy"]),
        "C4_post_pointwise_delta_proxy": sp.simplify(pre_proxy["C4_pre_pointwise_delta_proxy"]),
    }

    return {
        "symbols": {"z": z, "CF": cf},
        "h": h,
        "h_integral": h_integral,
        "pre_proxy": pre_proxy,
        "post_proxy": post_proxy,
    }


def compare_pointwise_proxy(canonical: Dict[str, Any]) -> Dict[str, Any]:
    pre = canonical["pre_proxy"]
    post = canonical["post_proxy"]

    invariance = {
        "C2_proxy_invariance": sp.simplify(post["C2_post_pointwise_delta_proxy"] - pre["C2_pre_pointwise_delta_proxy"]),
        "C4_proxy_invariance": sp.simplify(post["C4_post_pointwise_delta_proxy"] - pre["C4_pre_pointwise_delta_proxy"]),
    }

    if invariance["C2_proxy_invariance"] == 0 and invariance["C4_proxy_invariance"] == 0:
        verdict = "pointwise_artifact_invariant"
    else:
        verdict = "unexpected_rewrite_change"

    return {
        "verdict": verdict,
        "invariance": invariance,
    }


def numeric_sanity(args: argparse.Namespace, canonical: Dict[str, Any]) -> List[Dict[str, Any]]:
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
    pre_coeff = float(
        canonical["pre_proxy"]["C2_pre_pointwise_delta_proxy"].subs(
            {canonical["symbols"]["CF"]: cf}
        ).evalf()
    )
    post_coeff = float(
        canonical["post_proxy"]["C2_post_pointwise_delta_proxy"].subs(
            {canonical["symbols"]["CF"]: cf}
        ).evalf()
    )
    shift_coeff = float(canonical["h_integral"].subs({canonical["symbols"]["CF"]: cf}).evalf())

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
                    "preProxyAbs": common_abs_factor * pre_coeff,
                    "sharedShiftAbs": common_abs_factor * shift_coeff,
                    "postProxyAbs": common_abs_factor * post_coeff,
                    "proxyInvarianceAbs": common_abs_factor * (post_coeff - pre_coeff),
                }
            )
    return rows


def print_text_report(
    comparison: Dict[str, Any],
    canonical: Dict[str, Any],
    rows: List[Dict[str, Any]] | None,
    show_basis: bool,
) -> None:
    print("Analytic quark distribution-basis verifier")
    print("  caveat: this script tracks the inherited pointwise endpoint proxy only")
    print("  caveat: it does not test physical delta(1-z) agreement")

    print("\nVerdict")
    print(f"  verdict = {comparison['verdict']}")

    print("\nEndpoint kernel")
    print(f"  h(z)             = {qcdc_verify.expr_string(canonical['h'])}")
    print(f"  integral[h]      = {qcdc_verify.expr_string(canonical['h_integral'])}")

    print("\nPointwise endpoint proxy")
    print(f"  C2_pre_proxy     = {qcdc_verify.expr_string(canonical['pre_proxy']['C2_pre_pointwise_delta_proxy'])}")
    print(f"  C4_pre_proxy     = {qcdc_verify.expr_string(canonical['pre_proxy']['C4_pre_pointwise_delta_proxy'])}")
    print(f"  C2_post_proxy    = {qcdc_verify.expr_string(canonical['post_proxy']['C2_post_pointwise_delta_proxy'])}")
    print(f"  C4_post_proxy    = {qcdc_verify.expr_string(canonical['post_proxy']['C4_post_pointwise_delta_proxy'])}")

    print("\nProxy invariance")
    for key, value in comparison["invariance"].items():
        print(f"  {key:<20} = {qcdc_verify.expr_string(value)}")

    if show_basis:
        print("\nInterpretation")
        print("  Applying the same canonical endpoint rewrite to both pointwise proxy")
        print("  decompositions leaves the inherited proxy unchanged. This is expected")
        print("  because the underlying split basis never represented the virtual term")
        print("  as a genuine delta distribution in the first place.")

    if rows is not None:
        print("\nNumeric sanity")
        print("  channel  flavor  preProxyAbs        sharedShiftAbs     postProxyAbs       proxyInvarianceAbs")
        print("  -------  ------  -----------------  -----------------  -----------------  ------------------")
        for row in rows:
            print(
                f"  {row['channel']:<7}  {row['flavor']:<6}  "
                f"{row['preProxyAbs']:<17.8e}  {row['sharedShiftAbs']:<17.8e}  "
                f"{row['postProxyAbs']:<17.8e}  {row['proxyInvarianceAbs']:<18.8e}"
            )


def build_json_report(
    comparison: Dict[str, Any],
    canonical: Dict[str, Any],
    rows: List[Dict[str, Any]] | None,
) -> Dict[str, Any]:
    out = {
        "verdict": comparison["verdict"],
        "h": qcdc_verify.expr_string(canonical["h"]),
        "h_integral": qcdc_verify.expr_string(canonical["h_integral"]),
        "pre_proxy": {
            key: qcdc_verify.expr_string(value) for key, value in canonical["pre_proxy"].items()
        },
        "post_proxy": {
            key: qcdc_verify.expr_string(value) for key, value in canonical["post_proxy"].items()
        },
        "invariance": {
            key: qcdc_verify.expr_string(value) for key, value in comparison["invariance"].items()
        },
    }
    if rows is not None:
        out["numeric_sanity"] = rows
    return out


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    split_results = load_corrected_split_results()
    canonical = build_canonical_convolution_basis(split_results)
    comparison = compare_pointwise_proxy(canonical)
    rows = numeric_sanity(args, canonical) if args.numeric_sanity else None

    if args.json:
        print(json.dumps(build_json_report(comparison, canonical, rows), indent=2, sort_keys=True))
    else:
        print_text_report(comparison, canonical, rows, show_basis=args.show_basis)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
