#!/usr/bin/env python3.10
"""
Verify the local-z quark NLO coefficient structure in the Herwig/POLDIS
comparison basis.

This script is diagnostic-only. It reuses the regular-piece reduction already
established in verify_qcdc_c4_presence.py, derives the Herwig singular quark
sector in the same LL absolute basis, and checks whether the z-dependent
regular plus singular pieces match in the full quark NLO coefficient.

Important caveat: this verifier keeps the virtual term as a pointwise factor
`born_abs * vdelta`, so it is not a faithful delta(1-z) distribution
decomposition. Any residual inferred from the pointwise "born-side" proxy must
be treated as an endpoint-localization artifact, not a physics mismatch.
"""

from __future__ import annotations

import argparse
import json
from typing import Any, Dict, List, Sequence

import sympy as sp

import compare_polarized_nlo_terms as audit
import verify_qcdc_c4_presence as qcdc_verify


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--show-basis",
        action="store_true",
        help="Print the explicit singular, regular, and total coefficient bases.",
    )
    parser.add_argument(
        "--numeric-sanity",
        action="store_true",
        help="Evaluate flavor-resolved absolute regular/singular/total deltas at the control point.",
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
        "--sum-pdf",
        default=audit.DEFAULT_SUM_PDF,
        help="Unpolarized PDF set used for Herwig numeric sanity checks.",
    )
    parser.add_argument(
        "--diff-pdf",
        default=audit.DEFAULT_DIFF_PDF,
        help="Polarized PDF set used for POLDIS and Herwig numeric sanity checks.",
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
        "--input-alpha-s",
        type=float,
        default=audit.DEFAULT_INPUT_ALPHA_S,
        help="Herwig Matchbox NLOAlphaS input value at the input scale.",
    )
    parser.add_argument(
        "--input-scale",
        type=float,
        default=audit.DEFAULT_INPUT_SCALE,
        help="Herwig Matchbox NLOAlphaS input scale in GeV.",
    )
    parser.add_argument(
        "--s-hadronic",
        type=float,
        default=audit.DEFAULT_S,
        help="ep hadronic s in GeV^2.",
    )
    return parser.parse_args(argv)


def poldis_quark_targets() -> Dict[str, Any]:
    regular = qcdc_verify.poldis_target()
    z = regular["symbols"]["z"]
    cf = regular["symbols"]["CF"]
    xb = sp.symbols("xB", positive=True)
    logmu = regular["symbols"]["logmu"]

    c_sing_raw = cf * (4 * sp.log(1 - z) / (1 - z) - sp.Rational(3) / (1 - z) + 4 * logmu / (1 - z))
    c_delta_raw = cf * (
        -2 * sp.pi**2 / 3
        - 9
        + 3 * logmu
        - 3 * sp.log(1 - xb)
        + 2 * sp.log(1 - xb) ** 2
        + 4 * sp.log(1 - xb) * logmu
    )

    herwig_norm = sp.Rational(1, 2)
    c_sing_target = sp.simplify(herwig_norm * c_sing_raw)
    c_delta_target = sp.simplify(herwig_norm * c_delta_raw)

    return {
        "symbols": {"z": z, "CF": cf, "xB": xb, "logmu": logmu},
        "regular": regular,
        "singular_raw": {
            "C2_sing": c_sing_raw,
            "C4_sing": c_sing_raw,
            "C2_delta": c_delta_raw,
            "C4_delta": c_delta_raw,
        },
        "herwig_normalized": {
            "C2_sing": c_sing_target,
            "C4_sing": c_sing_target,
            "C2_delta": c_delta_target,
            "C4_delta": c_delta_target,
        },
    }


def herwig_quark_singular_basis() -> Dict[str, Any]:
    z, ell = sp.symbols("z ell", positive=True)
    cf = sp.symbols("CF", positive=True)
    q0, q, dq0, dq = sp.symbols("q0 q dq0 dq", nonzero=True)
    d0, dl, dq_coeff, dlq = sp.symbols("D0 Dl Dq Dlq")
    n0, nl, nq, nlq = sp.symbols("N0 Nl Nq Nlq")
    vec, ax = sp.symbols("vec ax")
    k1, k2, vdelta = sp.symbols("K1 K2 Vdelta")

    pref = sp.Rational(2) / (1 + ell) ** 2
    ym = 4 * ell / (1 + ell) ** 2
    yp = 2 * (1 + ell**2) / (1 + ell) ** 2

    helicity_exprs: Dict[str, sp.Expr] = {}
    for label, pl, pz in audit.HELICITY_COMBINATIONS:
        p_q = pz * dq0 / q0
        d_even = d0 + pl * dl
        n_even = n0 + pl * nl
        d_spin_coeff = dq_coeff + pl * dlq
        n_spin_coeff = nq + pl * nlq
        born_abs = pref * (
            (1 + ell**2) * (q0 * d_even + pz * dq0 * d_spin_coeff)
            + ell * (q0 * n_even + pz * dq0 * n_spin_coeff)
        )

        coll_unpol = cf / z * (q / q0 * k1 + (q / q0 - z) * k2)
        coll_pol = cf / z * (dq / dq0 * k1 + (dq / dq0 - z) * k2)

        q_u_num = (1 + ell**2) * d_even + ell * n_even
        q_odd_num = (1 + ell**2) * d_spin_coeff + ell * n_spin_coeff

        coll_abs = pref * (q0 * q_u_num * coll_unpol + q0 * p_q * q_odd_num * coll_pol)
        virt_abs = born_abs * vdelta
        helicity_exprs[label] = sp.expand(coll_abs + virt_abs)

    ll_expr = sp.simplify(
        sp.Rational(1, 4)
        * (
            helicity_exprs["PP"]
            + helicity_exprs["MM"]
            - helicity_exprs["PM"]
            - helicity_exprs["MP"]
        )
    )
    ll_substituted = sp.simplify(ll_expr.subs({dlq: ax, nlq: 2 * vec}))

    common = ym * vec + yp * ax
    scalar = sp.simplify(ll_substituted / common)

    dq_coeff_linear = sp.simplify(sp.diff(scalar, dq))
    dq0_coeff_linear = sp.simplify(sp.diff(scalar, dq0))

    c_sing_eff = sp.simplify(z * dq_coeff_linear)
    c_delta_eff = sp.simplify(dq0_coeff_linear + c_sing_eff)

    return {
        "symbols": {
            "z": z,
            "ell": ell,
            "CF": cf,
            "xB": sp.symbols("xB", positive=True),
            "q0": q0,
            "q": q,
            "dq0": dq0,
            "dq": dq,
            "vec": vec,
            "ax": ax,
            "K1": k1,
            "K2": k2,
            "Vdelta": vdelta,
        },
        "pref": pref,
        "ym": ym,
        "yp": yp,
        "helicity_exprs": helicity_exprs,
        "ll_expr_raw": ll_expr,
        "ll_expr_substituted": ll_substituted,
        "scalar": scalar,
        "C2_sing_eff": sp.simplify(c_sing_eff),
        "C4_sing_eff": sp.simplify(c_sing_eff),
        "C2_delta_eff": sp.simplify(c_delta_eff),
        "C4_delta_eff": sp.simplify(c_delta_eff),
        "delta_common_from_singular": sp.simplify(c_sing_eff - cf * k2),
    }


def evaluate_symbolic_basis() -> Dict[str, Any]:
    targets = poldis_quark_targets()
    regular_basis = qcdc_verify.ll_projection_and_basis_change()
    regular_analysis = qcdc_verify.verdict_and_sanity(targets["regular"], regular_basis)
    singular_basis = herwig_quark_singular_basis()

    z = targets["symbols"]["z"]
    cf = targets["symbols"]["CF"]
    xb = targets["symbols"]["xB"]
    logmu = targets["symbols"]["logmu"]
    k1 = singular_basis["symbols"]["K1"]
    k2 = singular_basis["symbols"]["K2"]
    vdelta = singular_basis["symbols"]["Vdelta"]

    substitutions = {
        logmu: 0,
        k1: 1 - z - 2 * sp.log(z) / (1 - z) - (1 + z) * sp.log((1 - z) / z),
        k2: 2 * sp.log(1 - z) / (1 - z) - sp.Rational(3, 2) / (1 - z),
        vdelta: cf
        * (
            -2 * sp.pi**2 / 3
            - sp.Rational(9, 2)
            - sp.Rational(3, 2) * sp.log(1 - xb)
            + sp.log(1 - xb) ** 2
        ),
    }

    c2_sing_eff = sp.simplify(singular_basis["C2_sing_eff"].subs(substitutions))
    c4_sing_eff = sp.simplify(singular_basis["C4_sing_eff"].subs(substitutions))
    c2_delta_eff = sp.simplify(singular_basis["C2_delta_eff"].subs(substitutions))
    c4_delta_eff = sp.simplify(singular_basis["C4_delta_eff"].subs(substitutions))

    c2_sing_target = sp.simplify(targets["herwig_normalized"]["C2_sing"].subs({logmu: 0}))
    c4_sing_target = sp.simplify(targets["herwig_normalized"]["C4_sing"].subs({logmu: 0}))
    c2_delta_target = sp.simplify(targets["herwig_normalized"]["C2_delta"].subs({logmu: 0}))
    c4_delta_target = sp.simplify(targets["herwig_normalized"]["C4_delta"].subs({logmu: 0}))

    c2_sing_diff = sp.simplify(c2_sing_eff - c2_sing_target)
    c4_sing_diff = sp.simplify(c4_sing_eff - c4_sing_target)
    c2_delta_diff = sp.simplify(c2_delta_eff - c2_delta_target)
    c4_delta_diff = sp.simplify(c4_delta_eff - c4_delta_target)

    regular_diffs = regular_analysis["coeff_differences"]

    c2_total_u_eff = sp.simplify(regular_basis["C2_eff"] + c2_sing_eff)
    c4_total_u_eff = sp.simplify(regular_basis["C4_eff"] + c4_sing_eff)
    cl_total_eff = sp.simplify(regular_basis["CL_eff"])
    c2_total_born_eff = sp.simplify(c2_delta_eff - c2_sing_eff)
    c4_total_born_eff = sp.simplify(c4_delta_eff - c4_sing_eff)

    c2_total_u_target = sp.simplify(targets["regular"]["herwig_normalized"]["C2_reg"].subs({logmu: 0}) + c2_sing_target)
    c4_total_u_target = sp.simplify(targets["regular"]["herwig_normalized"]["C4_reg"].subs({logmu: 0}) + c4_sing_target)
    cl_total_target = sp.simplify(targets["regular"]["herwig_normalized"]["CL"])
    c2_total_born_target = sp.simplify(c2_delta_target - c2_sing_target)
    c4_total_born_target = sp.simplify(c4_delta_target - c4_sing_target)

    total_diffs = {
        "C2_total_u_diff": sp.simplify(c2_total_u_eff - c2_total_u_target),
        "C4_total_u_diff": sp.simplify(c4_total_u_eff - c4_total_u_target),
        "CL_total_diff": sp.simplify(cl_total_eff - cl_total_target),
        "C2_total_born_diff": sp.simplify(c2_total_born_eff - c2_total_born_target),
        "C4_total_born_diff": sp.simplify(c4_total_born_eff - c4_total_born_target),
    }

    has_regular = any(diff != 0 for diff in regular_diffs.values())
    has_singular_split = any(diff != 0 for diff in [c2_sing_diff, c4_sing_diff])
    has_local_total = any(diff != 0 for diff in [total_diffs["C2_total_u_diff"], total_diffs["C4_total_u_diff"], total_diffs["CL_total_diff"]])

    regular_verdict = "different_split" if has_regular else "match"
    singular_verdict = "different_split" if has_singular_split else "match"
    total_verdict = "matched_local_z_structure" if not has_local_total else "different_local_z_structure"

    return {
        "targets": targets,
        "regular_basis": regular_basis,
        "regular_analysis": regular_analysis,
        "singular_basis": singular_basis,
        "singular_effective": {
            "C2_sing_eff": c2_sing_eff,
            "C4_sing_eff": c4_sing_eff,
            "C2_delta_eff": c2_delta_eff,
            "C4_delta_eff": c4_delta_eff,
        },
        "singular_targets": {
            "C2_sing_target": c2_sing_target,
            "C4_sing_target": c4_sing_target,
            "C2_delta_target": c2_delta_target,
            "C4_delta_target": c4_delta_target,
        },
        "singular_diffs": {
            "C2_sing_diff": c2_sing_diff,
            "C4_sing_diff": c4_sing_diff,
            "C2_delta_diff": c2_delta_diff,
            "C4_delta_diff": c4_delta_diff,
        },
        "total_effective": {
            "C2_total_u_eff": c2_total_u_eff,
            "C4_total_u_eff": c4_total_u_eff,
            "CL_total_eff": cl_total_eff,
            "C2_total_born_eff": c2_total_born_eff,
            "C4_total_born_eff": c4_total_born_eff,
        },
        "total_targets": {
            "C2_total_u_target": c2_total_u_target,
            "C4_total_u_target": c4_total_u_target,
            "CL_total_target": cl_total_target,
            "C2_total_born_target": c2_total_born_target,
            "C4_total_born_target": c4_total_born_target,
        },
        "total_diffs": total_diffs,
        "pointwise_endpoint_proxy": {
            "C2_delta_proxy_diff": c2_delta_diff,
            "C4_delta_proxy_diff": c4_delta_diff,
            "C2_total_born_proxy_diff": total_diffs["C2_total_born_diff"],
            "C4_total_born_proxy_diff": total_diffs["C4_total_born_diff"],
        },
        "verdicts": {
            "regular_verdict": regular_verdict,
            "singular_verdict": singular_verdict,
            "total_quark_verdict": total_verdict,
        },
        "endpoint_localization_status": "not_tested_in_pointwise_basis",
    }


def numeric_sanity(args: argparse.Namespace) -> List[Dict[str, Any]]:
    channels = [args.channel] if args.channel else ["GAMMA", "Z", "ALL"]
    lhapdf = audit.import_lhapdf()
    sum_pdf = lhapdf.mkPDF(args.sum_pdf, 0)
    diff_pdf = lhapdf.mkPDF(args.diff_pdf, 0)
    ew_inputs = audit.EWInputs(
        sin2_theta_w=args.sin2_theta_w,
        mz=args.mz,
        gamma_z=args.gamma_z,
        lepton_charge=-1.0,
    )
    lambda_squared = audit.solve_matchbox_lambda_squared(args.input_alpha_s, args.input_scale)
    quark_masses2 = audit.default_matchbox_quark_masses2()

    x_b = args.q2 / (args.y * args.s_hadronic)
    if not (0.0 < x_b < 1.0):
        raise SystemExit(f"x_B={x_b:.6f} lies outside (0,1) for the chosen numeric sanity point.")

    herwig_alpha_s = audit.matchbox_running_alpha(args.q2, lambda_squared, quark_masses2)
    poldis_alpha_s = diff_pdf.alphasQ2(args.q2)

    rows: List[Dict[str, Any]] = []
    for channel in channels:
        for flavor in audit.FLAVORS:
            herwig_terms = audit.aggregate_herwig_terms(
                flavor,
                args.q2,
                args.y,
                args.xp,
                args.s_hadronic,
                channel,
                ew_inputs,
                sum_pdf,
                diff_pdf,
                herwig_alpha_s,
            )
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

            delta_regular_abs = herwig_terms.qcdc_mapped_abs - poldis_terms.qcdc_mapped_abs
            delta_singular_abs = herwig_terms.quark_collinear_abs - poldis_terms.quark_collinear_abs
            delta_total_abs = delta_regular_abs + delta_singular_abs

            rows.append(
                {
                    "channel": channel,
                    "flavor": flavor.label,
                    "delta_regular_abs": delta_regular_abs,
                    "delta_singular_abs": delta_singular_abs,
                    "delta_total_abs": delta_total_abs,
                }
            )
    return rows


def print_text_report(results: Dict[str, Any], rows: List[Dict[str, Any]] | None, show_basis: bool) -> None:
    regular_diffs = results["regular_analysis"]["coeff_differences"]
    singular_diffs = results["singular_diffs"]
    total_diffs = results["total_diffs"]
    proxy = results["pointwise_endpoint_proxy"]

    print("Analytic quark-NLO split verifier")
    print("  caveat: symbolic proof is in the unclamped, non-endpoint analytic region")
    print("  caveat: endpoint-localized delta(1-z) terms are not extracted physically here")
    print("  established baseline: c4-c2 is already present in the Herwig regular kernel")

    print("\nVerdicts")
    print(f"  regular_verdict     = {results['verdicts']['regular_verdict']}")
    print(f"  singular_verdict    = {results['verdicts']['singular_verdict']}")
    print(f"  total_quark_verdict = {results['verdicts']['total_quark_verdict']}")

    print("\nRegular leftovers")
    print(f"  C2_reg diff = {qcdc_verify.expr_string(regular_diffs['C2_eff_minus_target'])}")
    print(f"  C4_reg diff = {qcdc_verify.expr_string(regular_diffs['C4_eff_minus_target'])}")
    print(f"  CL diff     = {qcdc_verify.expr_string(regular_diffs['CL_eff_minus_target'])}")

    print("\nSingular split leftovers")
    print(f"  C2_sing diff  = {qcdc_verify.expr_string(singular_diffs['C2_sing_diff'])}")
    print(f"  C4_sing diff  = {qcdc_verify.expr_string(singular_diffs['C4_sing_diff'])}")

    print("\nMatched local-z totals")
    print(f"  C2_total(dq/z) diff   = {qcdc_verify.expr_string(total_diffs['C2_total_u_diff'])}")
    print(f"  C4_total(dq/z) diff   = {qcdc_verify.expr_string(total_diffs['C4_total_u_diff'])}")
    print(f"  CL_total diff         = {qcdc_verify.expr_string(total_diffs['CL_total_diff'])}")

    print("\nEndpoint localization")
    print(f"  status                = {results['endpoint_localization_status']}")
    print("  note                  = pointwise born-side residuals from this script are endpoint artifacts")

    if show_basis:
        print("\nSingular basis")
        print(f"  C2_sing_eff   = {qcdc_verify.expr_string(results['singular_effective']['C2_sing_eff'])}")
        print(f"  C2_sing_tgt   = {qcdc_verify.expr_string(results['singular_targets']['C2_sing_target'])}")
        print(f"  C2_delta_eff  = {qcdc_verify.expr_string(results['singular_effective']['C2_delta_eff'])}")
        print(f"  C2_delta_tgt  = {qcdc_verify.expr_string(results['singular_targets']['C2_delta_target'])}")
        print(f"  C4_sing_eff   = {qcdc_verify.expr_string(results['singular_effective']['C4_sing_eff'])}")
        print(f"  C4_delta_eff  = {qcdc_verify.expr_string(results['singular_effective']['C4_delta_eff'])}")

        print("\nRegular basis")
        print(f"  C2_reg_eff    = {qcdc_verify.expr_string(results['regular_basis']['C2_eff'])}")
        print(f"  C4_reg_eff    = {qcdc_verify.expr_string(results['regular_basis']['C4_eff'])}")
        print(f"  CL_eff        = {qcdc_verify.expr_string(results['regular_basis']['CL_eff'])}")

        print("\nTotal basis")
        print(f"  C2_total(dq/z)_eff   = {qcdc_verify.expr_string(results['total_effective']['C2_total_u_eff'])}")
        print(f"  C4_total(dq/z)_eff   = {qcdc_verify.expr_string(results['total_effective']['C4_total_u_eff'])}")
        print(f"  CL_total_eff         = {qcdc_verify.expr_string(results['total_effective']['CL_total_eff'])}")

        print("\nPointwise endpoint proxy (not a physical delta mismatch)")
        print(f"  C2_delta_proxy_diff       = {qcdc_verify.expr_string(proxy['C2_delta_proxy_diff'])}")
        print(f"  C4_delta_proxy_diff       = {qcdc_verify.expr_string(proxy['C4_delta_proxy_diff'])}")
        print(f"  C2_total_born_proxy_diff  = {qcdc_verify.expr_string(proxy['C2_total_born_proxy_diff'])}")
        print(f"  C4_total_born_proxy_diff  = {qcdc_verify.expr_string(proxy['C4_total_born_proxy_diff'])}")

    if rows is not None:
        print("\nNumeric sanity")
        print("  channel  flavor  deltaRegularAbs     deltaSingularAbs   deltaTotalAbs")
        print("  -------  ------  ------------------  -----------------  ------------------")
        for row in rows:
            print(
                f"  {row['channel']:<7}  {row['flavor']:<6}  "
                f"{row['delta_regular_abs']:<18.8e}  {row['delta_singular_abs']:<17.8e}  "
                f"{row['delta_total_abs']:<18.8e}"
            )


def build_json_report(results: Dict[str, Any], rows: List[Dict[str, Any]] | None) -> Dict[str, Any]:
    out = {
        "verdicts": results["verdicts"],
        "endpoint_localization_status": results["endpoint_localization_status"],
        "regular_diffs": {key: qcdc_verify.expr_string(value) for key, value in results["regular_analysis"]["coeff_differences"].items()},
        "singular_split_diffs": {
            key: qcdc_verify.expr_string(value)
            for key, value in results["singular_diffs"].items()
            if key in ("C2_sing_diff", "C4_sing_diff")
        },
        "local_total_diffs": {
            key: qcdc_verify.expr_string(value)
            for key, value in results["total_diffs"].items()
            if key in ("C2_total_u_diff", "C4_total_u_diff", "CL_total_diff")
        },
        "pointwise_endpoint_proxy": {
            key: qcdc_verify.expr_string(value)
            for key, value in results["pointwise_endpoint_proxy"].items()
        },
    }
    if rows is not None:
        out["numeric_sanity"] = rows
    return out


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    results = evaluate_symbolic_basis()
    rows = numeric_sanity(args) if args.numeric_sanity else None

    if args.json:
        print(json.dumps(build_json_report(results, rows), indent=2, sort_keys=True))
    else:
        print_text_report(results, rows, show_basis=args.show_basis)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
