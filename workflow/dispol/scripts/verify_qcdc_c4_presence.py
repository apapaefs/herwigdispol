#!/usr/bin/env python3.10
"""
Analytically verify whether the polarized NC quark-regular g4-g2 difference is
already encoded in Herwig's QCDC `realq` term.

This script is intentionally diagnostic-only. It does not change the production
code. It derives the POLDIS regular-coefficient target, reduces the current
Herwig `realq` expression to the LL absolute-cross-section basis, extracts the
effective coefficient functions, and reports:

  - whether the full quark-regular basis matches,
  - whether the specific c4-c2 difference is already present,
  - and, optionally, numerical sanity checks against the existing
    compare_polarized_nlo_terms.py diagnostic helpers.
"""

from __future__ import annotations

import argparse
import json
import math
from typing import Any, Dict, Iterable, List, Sequence

import sympy as sp

import compare_polarized_nlo_terms as audit


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--show-herwig-basis",
        action="store_true",
        help="Print the reduced Herwig LL basis before comparison.",
    )
    parser.add_argument(
        "--numeric-sanity",
        action="store_true",
        help="Evaluate the symbolic results numerically at the established control point.",
    )
    parser.add_argument(
        "--channel",
        choices=("GAMMA", "Z", "ALL"),
        default=None,
        help="Restrict numeric sanity checks to one channel.",
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
        help="Polarized PDF set used for numeric sanity checks.",
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


def expr_string(expr: sp.Expr) -> str:
    return sp.sstr(sp.factor(sp.simplify(expr)))


def poldis_target() -> Dict[str, Any]:
    z = sp.symbols("z", positive=True)
    cf = sp.symbols("CF", positive=True)
    logmu = sp.symbols("logmu", real=True)

    c2_reg = cf * (
        -2 * (1 + z) * sp.log(1 - z)
        - 2 * (1 + z**2) / (1 - z) * sp.log(z)
        + 4
        + 2 * z
        - 2 * (1 + z) * logmu
    )
    c4_reg = cf * (
        -2 * (1 + z) * sp.log(1 - z)
        - 2 * (1 + z**2) / (1 - z) * sp.log(z)
        + 6
        + 4 * z
        - 2 * (1 + z) * logmu
    )
    cl_reg = 4 * cf * z

    # compare_polarized_nlo_terms.py mirrors the POLDIS absolute piece with an
    # overall alpha_s/(4 pi) prefactor, while Herwig's realq kernel carries
    # alpha_s/(2 pi). For coefficient matching to Herwig's local regular term,
    # use the same convention scaled to the Herwig normalization.
    herwig_norm = sp.Rational(1, 2)
    c2_target = sp.simplify(herwig_norm * c2_reg)
    c4_target = sp.simplify(herwig_norm * c4_reg)
    cl_target = sp.simplify(herwig_norm * cl_reg)

    return {
        "symbols": {"z": z, "CF": cf, "logmu": logmu},
        "poldis_raw": {
            "C2_reg": c2_reg,
            "C4_reg": c4_reg,
            "CL": cl_reg,
            "C4_minus_C2": sp.simplify(c4_reg - c2_reg),
        },
        "herwig_normalized": {
            "C2_reg": c2_target,
            "C4_reg": c4_target,
            "CL": cl_target,
            "C4_minus_C2": sp.simplify(c4_target - c2_target),
        },
    }


def herwig_symbolic_reduction() -> Dict[str, Any]:
    z, ell = sp.symbols("z ell", positive=True)
    cf = sp.symbols("CF", positive=True)
    q, lo = sp.symbols("q lo", nonzero=True)
    db, nb, dm, nm = sp.symbols("D_born N_born D_mapped N_mapped", nonzero=True)

    pref = sp.Rational(2) / (1 + ell) ** 2
    q_ratio = q / lo
    qcdc_den_ratio = dm / db
    a_born = nb / db
    a_q_mapped = nm / dm
    e_kernel = 2 + 2 * ell**2 - z + 3 * z * ell**2
    a_kernel = 2 * z + 1

    born_abs = lo * pref * db * (1 + a_born * ell + ell**2)
    realq = (
        qcdc_den_ratio
        * cf
        / z
        / (1 + a_born * ell + ell**2)
        * q_ratio
        * (e_kernel + a_q_mapped * ell * a_kernel)
    )
    reduced = sp.simplify(born_abs * realq)

    return {
        "symbols": {
            "z": z,
            "ell": ell,
            "CF": cf,
            "q": q,
            "lo": lo,
            "D_born": db,
            "N_born": nb,
            "D_mapped": dm,
            "N_mapped": nm,
        },
        "pref": pref,
        "E": e_kernel,
        "A": a_kernel,
        "direct_realg_abs": born_abs * realq,
        "reduced_realg_abs": reduced,
    }


def ll_projection_and_basis_change() -> Dict[str, Any]:
    z, ell = sp.symbols("z ell", positive=True)
    cf = sp.symbols("CF", positive=True)
    q, dq = sp.symbols("q dq", nonzero=True)
    d0, dl, dq_coeff, dlq = sp.symbols("D0 Dl Dq Dlq")
    n0, nl, nq, nlq = sp.symbols("N0 Nl Nq Nlq")
    vec, ax = sp.symbols("vec ax")

    pref = sp.Rational(2) / (1 + ell) ** 2
    e_kernel = 2 + 2 * ell**2 - z + 3 * z * ell**2
    a_kernel = 2 * z + 1

    helicity_exprs: Dict[str, sp.Expr] = {}
    for label, pl, pz in audit.HELICITY_COMBINATIONS:
        pq_m = pz * dq / q
        d_mapped = d0 + pl * dl + pq_m * dq_coeff + pl * pq_m * dlq
        n_mapped = n0 + pl * nl + pq_m * nq + pl * pq_m * nlq
        helicity_exprs[label] = sp.expand(pref * cf / z * q * (d_mapped * e_kernel + ell * n_mapped * a_kernel))

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
    common = dq / z
    kernel = sp.simplify(ll_substituted / common)

    ym = 4 * ell / (1 + ell) ** 2
    yp = 2 * (1 + ell**2) / (1 + ell) ** 2
    y_sq = 4 / (1 + ell) ** 2

    vec_coeff = sp.simplify(sp.diff(kernel, vec))
    ax_coeff = sp.simplify(sp.diff(kernel, ax))

    c2_eff = sp.simplify(vec_coeff / ym)

    c4_eff, cl_eff = sp.symbols("C4_eff CL_eff")
    numerator_equation = sp.expand(
        (1 + ell) ** 2 * ax_coeff - (2 * (1 + ell**2) * c4_eff - 4 * cl_eff)
    )
    poly = sp.Poly(numerator_equation, ell)
    equations = [sp.Eq(coeff, 0) for coeff in poly.all_coeffs()]
    solution = sp.solve(equations, (c4_eff, cl_eff), dict=True)
    if not solution:
        raise RuntimeError("Could not solve for effective C4/CL coefficients.")
    solved = solution[0]

    return {
        "symbols": {
            "z": z,
            "ell": ell,
            "CF": cf,
            "q": q,
            "dq": dq,
            "vec": vec,
            "ax": ax,
        },
        "pref": pref,
        "ym": ym,
        "yp": yp,
        "y_sq": y_sq,
        "helicity_exprs": helicity_exprs,
        "ll_expr_raw": ll_expr,
        "ll_expr_substituted": ll_substituted,
        "kernel": kernel,
        "vec_coeff": vec_coeff,
        "ax_coeff": ax_coeff,
        "C2_eff": sp.simplify(c2_eff),
        "C4_eff": sp.simplify(solved[c4_eff]),
        "CL_eff": sp.simplify(solved[cl_eff]),
        "identities_used": {
            "Dlq": ax,
            "Nlq": 2 * vec,
        },
    }


def verdict_and_sanity(
    poldis: Dict[str, Any],
    herwig: Dict[str, Any],
) -> Dict[str, Any]:
    z = poldis["symbols"]["z"]
    logmu = poldis["symbols"]["logmu"]

    c2_target = poldis["herwig_normalized"]["C2_reg"]
    c4_target = poldis["herwig_normalized"]["C4_reg"]
    cl_target = poldis["herwig_normalized"]["CL"]

    c2_eff = herwig["C2_eff"]
    c4_eff = herwig["C4_eff"]
    cl_eff = herwig["CL_eff"]

    subs_logmu_zero = {logmu: 0}
    c2_diff = sp.simplify((c2_eff - c2_target).subs(subs_logmu_zero))
    c4_diff = sp.simplify((c4_eff - c4_target).subs(subs_logmu_zero))
    cl_diff = sp.simplify((cl_eff - cl_target).subs(subs_logmu_zero))

    c4_minus_c2_eff = sp.simplify(c4_eff - c2_eff)
    c4_minus_c2_target = sp.simplify((c4_target - c2_target).subs(subs_logmu_zero))
    c4_minus_c2_diff = sp.simplify(c4_minus_c2_eff - c4_minus_c2_target)

    if c2_diff == 0 and c4_diff == 0 and cl_diff == 0:
        verdict = "already_there"
    elif c2_diff == 0 and cl_diff == 0 and sp.simplify(c4_diff + c4_minus_c2_target) == 0:
        verdict = "missing_pure_c4_minus_c2"
    else:
        verdict = "different_leftover"

    return {
        "verdict": verdict,
        "c4_difference_verdict": "present" if c4_minus_c2_diff == 0 else "absent_or_reshuffled",
        "coeff_differences": {
            "C2_eff_minus_target": c2_diff,
            "C4_eff_minus_target": c4_diff,
            "CL_eff_minus_target": cl_diff,
        },
        "c4_difference": {
            "effective": c4_minus_c2_eff,
            "target": c4_minus_c2_target,
            "difference": c4_minus_c2_diff,
        },
    }


def numeric_sanity(
    args: argparse.Namespace,
    analysis: Dict[str, Any],
) -> List[Dict[str, Any]]:
    channels = [args.channel] if args.channel else ["GAMMA", "Z", "ALL"]
    lhapdf = audit.import_lhapdf()
    diff_pdf = lhapdf.mkPDF(args.diff_pdf, 0)
    ew_inputs = audit.EWInputs(
        sin2_theta_w=args.sin2_theta_w,
        mz=args.mz,
        gamma_z=args.gamma_z,
        lepton_charge=-1.0,
    )

    z = args.xp
    x_b = args.q2 / (args.y * args.s_hadronic)
    if not (0.0 < x_b < 1.0):
        raise SystemExit(f"x_B={x_b:.6f} lies outside (0,1) for the chosen numeric sanity point.")

    cf = 4.0 / 3.0
    yp = 1.0 + (1.0 - args.y) ** 2
    poldis_alpha_s = diff_pdf.alphasQ2(args.q2)
    a2pi = poldis_alpha_s / (2.0 * math.pi)

    rows: List[Dict[str, Any]] = []
    for channel in channels:
        for flavor in audit.FLAVORS:
            pieces = audit.poldis_channel_pieces(flavor, args.q2, ew_inputs, channel)
            poldis_terms = audit.poldis_coefficients_single_flavor(
                args.q2, x_b, args.y, z, channel, flavor, ew_inputs, diff_pdf, poldis_alpha_s
            )
            u = x_b / z
            fu = diff_pdf.xfxQ2(flavor.lhapdf_id, u, args.q2)
            c4_minus_c2_raw = 2.0 * cf * (1.0 + z)
            formula_abs = (x_b * (a2pi / 2.0) * yp * ((c4_minus_c2_raw * fu * pieces.ax / z))) / (x_b * x_b)

            coeff = audit.herwig_channel_coefficients(flavor, args.q2, ew_inputs, channel)
            channel_identity = {
                "Dlq_minus_ax": coeff.Dlq - pieces.ax,
                "Nlq_minus_2vec": coeff.Nlq - 2.0 * pieces.vec,
            }
            rows.append(
                {
                    "channel": channel,
                    "flavor": flavor.label,
                    "formula_abs": formula_abs,
                    "compare_script_abs": poldis_terms.inputs["qcdc_c4_minus_c2_abs"],
                    "difference": formula_abs - poldis_terms.inputs["qcdc_c4_minus_c2_abs"],
                    "ax": pieces.ax,
                    "vec": pieces.vec,
                    "Dlq_minus_ax": channel_identity["Dlq_minus_ax"],
                    "Nlq_minus_2vec": channel_identity["Nlq_minus_2vec"],
                }
            )
    return rows


def print_text_report(
    poldis: Dict[str, Any],
    herwig_reduction: Dict[str, Any],
    herwig_basis: Dict[str, Any],
    analysis: Dict[str, Any],
    numeric_rows: List[Dict[str, Any]] | None,
    show_herwig_basis: bool,
) -> None:
    print("Analytic QCDC c4-c2 verifier")
    print("  caveat: LL derivation is done in the unclamped analytic region")
    print("  identities used for basis change: Dlq = ax, Nlq = 2*vec")

    print("\nPOLDIS regular quark target")
    print(f"  C2_reg(raw)   = {expr_string(poldis['poldis_raw']['C2_reg'])}")
    print(f"  C4_reg(raw)   = {expr_string(poldis['poldis_raw']['C4_reg'])}")
    print(f"  CL(raw)       = {expr_string(poldis['poldis_raw']['CL'])}")
    print(f"  C4-C2(raw)    = {expr_string(poldis['poldis_raw']['C4_minus_C2'])}")
    print(f"  C2_target     = {expr_string(poldis['herwig_normalized']['C2_reg'].subs({poldis['symbols']['logmu']: 0}))}")
    print(f"  C4_target     = {expr_string(poldis['herwig_normalized']['C4_reg'].subs({poldis['symbols']['logmu']: 0}))}")
    print(f"  CL_target     = {expr_string(poldis['herwig_normalized']['CL'])}")
    print(f"  C4-C2(target) = {expr_string(poldis['herwig_normalized']['C4_minus_C2'].subs({poldis['symbols']['logmu']: 0}))}")

    if show_herwig_basis:
        print("\nHerwig reduction")
        print(f"  reduced realq_abs = {expr_string(herwig_reduction['reduced_realg_abs'])}")
        print(f"  LL kernel         = {expr_string(herwig_basis['kernel'])}")
        print(f"  vec coefficient   = {expr_string(herwig_basis['vec_coeff'])}")
        print(f"  ax coefficient    = {expr_string(herwig_basis['ax_coeff'])}")

    print("\nHerwig effective coefficients")
    print(f"  C2_eff = {expr_string(herwig_basis['C2_eff'])}")
    print(f"  C4_eff = {expr_string(herwig_basis['C4_eff'])}")
    print(f"  CL_eff = {expr_string(herwig_basis['CL_eff'])}")

    print("\nVerdict")
    print(f"  full_basis_verdict = {analysis['verdict']}")
    print(f"  c4_difference      = {analysis['c4_difference_verdict']}")
    print(f"  C2_eff-target      = {expr_string(analysis['coeff_differences']['C2_eff_minus_target'])}")
    print(f"  C4_eff-target      = {expr_string(analysis['coeff_differences']['C4_eff_minus_target'])}")
    print(f"  CL_eff-target      = {expr_string(analysis['coeff_differences']['CL_eff_minus_target'])}")
    print(f"  (C4-C2)_eff-target = {expr_string(analysis['c4_difference']['difference'])}")

    if numeric_rows is not None:
        print("\nNumeric sanity")
        print(
            "  channel  flavor  formulaAbs          compareScriptAbs     diff                 "
            "Dlq-ax             Nlq-2vec"
        )
        print(
            "  -------  ------  ------------------  ------------------  ------------------  "
            "-----------------  -----------------"
        )
        for row in numeric_rows:
            print(
                f"  {row['channel']:<7}  {row['flavor']:<6}  "
                f"{row['formula_abs']:<18.8e}  {row['compare_script_abs']:<18.8e}  "
                f"{row['difference']:<18.8e}  {row['Dlq_minus_ax']:<17.8e}  "
                f"{row['Nlq_minus_2vec']:<17.8e}"
            )


def build_json_report(
    poldis: Dict[str, Any],
    herwig_reduction: Dict[str, Any],
    herwig_basis: Dict[str, Any],
    analysis: Dict[str, Any],
    numeric_rows: List[Dict[str, Any]] | None,
) -> Dict[str, Any]:
    logmu = poldis["symbols"]["logmu"]
    report = {
        "verdict": analysis["verdict"],
        "c4_difference_verdict": analysis["c4_difference_verdict"],
        "poldis_target": {
            "C2_reg_raw": expr_string(poldis["poldis_raw"]["C2_reg"]),
            "C4_reg_raw": expr_string(poldis["poldis_raw"]["C4_reg"]),
            "CL_raw": expr_string(poldis["poldis_raw"]["CL"]),
            "C4_minus_C2_raw": expr_string(poldis["poldis_raw"]["C4_minus_C2"]),
            "C2_target": expr_string(poldis["herwig_normalized"]["C2_reg"].subs({logmu: 0})),
            "C4_target": expr_string(poldis["herwig_normalized"]["C4_reg"].subs({logmu: 0})),
            "CL_target": expr_string(poldis["herwig_normalized"]["CL"]),
            "C4_minus_C2_target": expr_string(poldis["herwig_normalized"]["C4_minus_C2"].subs({logmu: 0})),
        },
        "herwig_basis": {
            "reduced_realq_abs": expr_string(herwig_reduction["reduced_realg_abs"]),
            "ll_kernel": expr_string(herwig_basis["kernel"]),
            "C2_eff": expr_string(herwig_basis["C2_eff"]),
            "C4_eff": expr_string(herwig_basis["C4_eff"]),
            "CL_eff": expr_string(herwig_basis["CL_eff"]),
        },
        "differences": {
            "c2_diff": expr_string(analysis["coeff_differences"]["C2_eff_minus_target"]),
            "c4_diff": expr_string(analysis["coeff_differences"]["C4_eff_minus_target"]),
            "cl_diff": expr_string(analysis["coeff_differences"]["CL_eff_minus_target"]),
            "c4_minus_c2_diff": expr_string(analysis["c4_difference"]["difference"]),
        },
    }
    if numeric_rows is not None:
        report["numeric_sanity"] = numeric_rows
    return report


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    poldis = poldis_target()
    herwig_reduction = herwig_symbolic_reduction()
    herwig_basis = ll_projection_and_basis_change()
    analysis = verdict_and_sanity(poldis, herwig_basis)
    numeric_rows = numeric_sanity(args, analysis) if args.numeric_sanity else None

    if args.json:
        print(
            json.dumps(
                build_json_report(poldis, herwig_reduction, herwig_basis, analysis, numeric_rows),
                indent=2,
                sort_keys=True,
            )
        )
    else:
        print_text_report(
            poldis,
            herwig_reduction,
            herwig_basis,
            analysis,
            numeric_rows,
            show_herwig_basis=args.show_herwig_basis,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
