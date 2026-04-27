#!/usr/bin/env python3
"""Symbolic audit of the local NC BGF kernel: POLDIS vs Herwig.

This script compares the *local* boson-gluon-fusion real-emission object used
by POLDIS to the branchwise `R2/R3` construction implemented in Herwig.

It does not attempt a full phase-space or jet-analysis validation. Instead, it
answers one narrower physics question:

  Can the exact local gluon-channel tensor structures in
  `POLDIS-public/poldis.f::MATTHR` + `COMBINEGLUON`
  be rewritten exactly in the same branchwise form used by Herwig's mapped
  neutral-current BGF kernel?

The answer should be "yes" if the local matrix-element construction is
conceptually correct. If so, any remaining discrepancy with POLDIS must come
from somewhere other than the local BGF tensor algebra itself.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass

import sympy as sp


@dataclass(frozen=True)
class LocalKinematicBasis:
    """POLDIS local gluon-channel tensor basis."""

    k36: sp.Symbol  # (p3.p6)^2
    k37: sp.Symbol  # (p3.p7)^2
    k27: sp.Symbol  # (p2.p7)^2
    k26: sp.Symbol  # (p2.p6)^2


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--show-expanded",
        action="store_true",
        help="Print expanded intermediate symbolic expressions.",
    )
    return parser.parse_args()


def local_basis() -> LocalKinematicBasis:
    return LocalKinematicBasis(
        k36=sp.Symbol("k36"),
        k37=sp.Symbol("k37"),
        k27=sp.Symbol("k27"),
        k26=sp.Symbol("k26"),
    )


def poldis_local_structures(basis: LocalKinematicBasis) -> dict[str, sp.Expr]:
    """Exact local gluon-channel structures from POLDIS MATTHR.

    The naming follows the code in:
      - poldis.f lines 975-1007
      - poldis_gamma_nlo_blocks.f lines 7893-7901

    `vec` is the VVVV-like coupling combination (`CS` in COMBINEGLUON),
    `ax` is the AVAV-like coupling combination (`CSAX` in COMBINEGLUON).
    """

    vec, ax = sp.symbols("vec ax")

    gq_unpol = basis.k36 + basis.k37 + basis.k27 + basis.k26
    gqz_unpol = basis.k36 - basis.k37 + basis.k27 - basis.k26
    gq_pol = -basis.k36 + basis.k37 + basis.k27 - basis.k26
    gqz_pol = -basis.k36 - basis.k37 + basis.k27 + basis.k26

    sigma0 = sp.expand(vec * gq_unpol + ax * gqz_unpol)
    sigma_ll = sp.expand(vec * gq_pol + ax * gqz_pol)

    return {
        "vec": vec,
        "ax": ax,
        "gq_unpol": gq_unpol,
        "gqz_unpol": gqz_unpol,
        "gq_pol": gq_pol,
        "gqz_pol": gqz_pol,
        "sigma0": sigma0,
        "sigma_ll": sigma_ll,
    }


def herwig_branchwise_reconstruction(
    basis: LocalKinematicBasis, vec: sp.Expr, ax: sp.Expr
) -> dict[str, sp.Expr]:
    """Rewrite the local object in Herwig's branchwise R2/R3 language.

    The exact POLDIS object can be grouped as

      R2(q,+Pg)  = k27 (vec+ax) + k37 (vec-ax)
      R3(qc,-Pg) = k36 (vec+ax) + k26 (vec-ax)

    Equivalently:

      R2 = even2 * Dq   + odd2 * Nq
      R3 = even3 * Dqcc + odd3 * Nqcc

    with Dqcc = Dq = vec and Nqcc = -Nq = -ax.
    """

    even2 = basis.k27 + basis.k37
    odd2 = basis.k27 - basis.k37
    even3 = basis.k36 + basis.k26
    # For the charge-conjugate R3 branch the odd kernel is written with the
    # sign convention that matches N(q^c) = -N(q).
    odd3 = basis.k26 - basis.k36

    d_q = vec
    n_q = ax
    d_q_cc = vec
    n_q_cc = -ax

    r2_lr = sp.expand(basis.k27 * (vec + ax) + basis.k37 * (vec - ax))
    r3_lr = sp.expand(basis.k36 * (vec + ax) + basis.k26 * (vec - ax))

    r2_dn = sp.expand(even2 * d_q + odd2 * n_q)
    r3_dn = sp.expand(even3 * d_q_cc + odd3 * n_q_cc)

    sigma0_branch = sp.expand(r2_lr + r3_lr)
    sigma_ll_branch = sp.expand(r2_lr - r3_lr)

    return {
        "even2": even2,
        "odd2": odd2,
        "even3": even3,
        "odd3": odd3,
        "d_q": d_q,
        "n_q": n_q,
        "d_q_cc": d_q_cc,
        "n_q_cc": n_q_cc,
        "r2_lr": r2_lr,
        "r3_lr": r3_lr,
        "r2_dn": r2_dn,
        "r3_dn": r3_dn,
        "sigma0_branch": sigma0_branch,
        "sigma_ll_branch": sigma_ll_branch,
    }


def full_nc_coupling_structures() -> dict[str, sp.Expr]:
    """Exact coupling combinations from COMBINEGLUON.

    These are the symbolic versions of the coefficients multiplying the local
    gluon-channel tensor structures in POLDIS.
    """

    eq, cvq, caq, cve, cae, propz, intgz = sp.symbols(
        "eq cvq caq cve cae propz intgz"
    )

    vec = (
        eq**2
        + propz * (cvq**2 + caq**2) * (cve**2 + cae**2)
        + 2 * eq * propz * intgz * cvq * cve
    )
    ax = (
        propz * (4 * cvq * caq * cve * cae)
        + 2 * eq * propz * intgz * caq * cae
    )

    cc_map = {eq: -eq, cvq: -cvq, caq: caq}
    vec_cc = sp.expand(vec.subs(cc_map))
    ax_cc = sp.expand(ax.subs(cc_map))

    return {
        "eq": eq,
        "cvq": cvq,
        "caq": caq,
        "cve": cve,
        "cae": cae,
        "propz": propz,
        "intgz": intgz,
        "vec": sp.expand(vec),
        "ax": sp.expand(ax),
        "vec_cc": vec_cc,
        "ax_cc": ax_cc,
    }


def photon_limit(vec: sp.Expr, ax: sp.Expr) -> dict[str, sp.Expr]:
    """Photon-only limit of the COMBINEGLUON couplings."""
    eq = sp.Symbol("eq")
    photon_subs = {
        sp.Symbol("propz"): 0,
        sp.Symbol("intgz"): 0,
        sp.Symbol("cvq"): 0,
        sp.Symbol("caq"): 0,
        sp.Symbol("cve"): 0,
        sp.Symbol("cae"): 0,
    }
    return {
        "vec_gamma": sp.expand(vec.subs(photon_subs)),
        "ax_gamma": sp.expand(ax.subs(photon_subs)),
        "expected_vec_gamma": eq**2,
        "expected_ax_gamma": sp.Integer(0),
    }


def simplify_zero(label: str, expr: sp.Expr) -> tuple[str, sp.Expr]:
    return label, sp.factor(sp.simplify(expr))


def main() -> int:
    args = parse_args()

    basis = local_basis()
    poldis = poldis_local_structures(basis)
    herwig = herwig_branchwise_reconstruction(basis, poldis["vec"], poldis["ax"])
    couplings = full_nc_coupling_structures()
    gamma = photon_limit(couplings["vec"], couplings["ax"])

    checks = [
        simplify_zero("R2 left/right vs D/N form", herwig["r2_lr"] - herwig["r2_dn"]),
        simplify_zero("R3 left/right vs D/N form", herwig["r3_lr"] - herwig["r3_dn"]),
        simplify_zero("sigma0 exact local decomposition", poldis["sigma0"] - herwig["sigma0_branch"]),
        simplify_zero(
            "sigmaLL exact local decomposition",
            poldis["sigma_ll"] - herwig["sigma_ll_branch"],
        ),
        simplify_zero("charge conjugation leaves vec unchanged", couplings["vec_cc"] - couplings["vec"]),
        simplify_zero("charge conjugation flips ax", couplings["ax_cc"] + couplings["ax"]),
        simplify_zero("photon-limit vec", gamma["vec_gamma"] - gamma["expected_vec_gamma"]),
        simplify_zero("photon-limit ax", gamma["ax_gamma"] - gamma["expected_ax_gamma"]),
    ]

    print("Local NC BGF symbolic audit: POLDIS vs Herwig branchwise R2/R3\n")

    print("POLDIS local structures from MATTHR + COMBINEGLUON:")
    print(f"  sigma0  = {sp.sstr(poldis['sigma0'])}")
    print(f"  sigmaLL = {sp.sstr(poldis['sigma_ll'])}\n")

    print("Herwig-style branchwise reconstruction:")
    print(f"  R2(q,+Pg)   = {sp.sstr(herwig['r2_lr'])}")
    print(f"  R3(qc,-Pg)  = {sp.sstr(herwig['r3_lr'])}")
    print(f"  sigma0  = R2 + R3 = {sp.sstr(herwig['sigma0_branch'])}")
    print(f"  sigmaLL = R2 - R3 = {sp.sstr(herwig['sigma_ll_branch'])}\n")

    print("Full-NC COMBINEGLUON couplings:")
    print(f"  vec(q)   = {sp.sstr(couplings['vec'])}")
    print(f"  ax(q)    = {sp.sstr(couplings['ax'])}")
    print(f"  vec(qc)  = {sp.sstr(couplings['vec_cc'])}")
    print(f"  ax(qc)   = {sp.sstr(couplings['ax_cc'])}\n")

    if args.show_expanded:
        print("Expanded branch kernels:")
        print(f"  even2 = {sp.sstr(herwig['even2'])}")
        print(f"  odd2  = {sp.sstr(herwig['odd2'])}")
        print(f"  even3 = {sp.sstr(herwig['even3'])}")
        print(f"  odd3  = {sp.sstr(herwig['odd3'])}\n")

    print("Symbolic checks:")
    failures = 0
    for label, result in checks:
        ok = result == 0
        status = "OK" if ok else "FAIL"
        print(f"  [{status}] {label}: {sp.sstr(result)}")
        if not ok:
            failures += 1

    print()
    if failures:
        print("Conclusion: the local gluon-channel objects are not identical under the tested symbolic rewriting.")
        return 1

    print(
        "Conclusion: the exact local POLDIS gluon-channel tensor structures can be "
        "rewritten exactly as the same branchwise R2/R3 local object used by Herwig."
    )
    print(
        "This validates the local BGF matrix-element construction itself; any "
        "remaining Herwig/POLDIS discrepancy must come from somewhere else."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
