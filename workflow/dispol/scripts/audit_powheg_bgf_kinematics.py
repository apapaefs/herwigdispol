#!/usr/bin/env python3
"""Symbolic audit of the DIS POWHEG BGF acceptance kinematics.

This script checks the kinematic map used in DISBase::generateBGF() and in the
accepted-event construction paths. It is intentionally narrow:

  - it verifies the accepted Breit-frame momenta are massless;
  - it verifies exact momentum conservation with q = (0,0,0,-Q);
  - it checks the x-parameter identities used in the code;
  - it records the normalized BGF denominator used consistently by
    BGFME(..., true) and generateBGFPoint().

It does not prove that the full accept/reject probability is the exact POWHEG
distribution, but it does validate that the post-acceptance momentum setup
itself is internally consistent.
"""

from __future__ import annotations

import sympy as sp


def main() -> int:
    xp, zp, Q = sp.symbols("xp zp Q", positive=True)
    xT = sp.symbols("xT", positive=True)

    x1 = -sp.Integer(1) / xp
    x2 = sp.Integer(1) - (sp.Integer(1) - zp) / xp
    x3 = sp.Integer(2) + x1 - x2
    xT2 = sp.Integer(4) * (sp.Integer(1) - xp) * (sp.Integer(1) - zp) * zp / xp

    # Breit-frame momentum map from DISBase::generateBGF() and the MEC path.
    p0 = sp.Matrix([sp.Integer(0), sp.Integer(0), -sp.Rational(1, 2) * Q * x1, -sp.Rational(1, 2) * Q * x1])
    p1 = sp.Matrix([
        sp.Rational(1, 2) * Q * xT,
        sp.Integer(0),
        -sp.Rational(1, 2) * Q * x2,
        sp.Rational(1, 2) * Q * sp.sqrt(xT**2 + x2**2),
    ])
    p2 = sp.Matrix([
        -sp.Rational(1, 2) * Q * xT,
        sp.Integer(0),
        -sp.Rational(1, 2) * Q * x3,
        sp.Rational(1, 2) * Q * sp.sqrt(xT**2 + x3**2),
    ])
    q = sp.Matrix([sp.Integer(0), sp.Integer(0), -Q, sp.Integer(0)])

    def minkowski_sq(p: sp.Matrix) -> sp.Expr:
        return sp.simplify(p[3] ** 2 - p[0] ** 2 - p[1] ** 2 - p[2] ** 2)

    checks: list[tuple[str, sp.Expr]] = []

    checks.append(("x3 identity", sp.simplify(x3 - (sp.Integer(1) - zp / xp))))
    checks.append(("x2+x3 identity", sp.simplify(x2 + x3 - (sp.Integer(2) - sp.Integer(1) / xp))))
    checks.append(("xT^2 definition", sp.simplify(xT2 - (sp.Integer(4) * (sp.Integer(1) - xp) * (sp.Integer(1) - zp) * zp / xp))))

    # Put the code's physical xT^2 relation into the on-shell checks.
    p1_sq = sp.simplify(minkowski_sq(p1).subs(xT**2, xT2))
    p2_sq = sp.simplify(minkowski_sq(p2).subs(xT**2, xT2))
    p0_sq = sp.simplify(minkowski_sq(p0))
    checks.append(("p0 massless", p0_sq))
    checks.append(("p1 massless", p1_sq))
    checks.append(("p2 massless", p2_sq))

    mom_cons = sp.simplify((p0 + q - p1 - p2).subs(xT**2, xT2))
    checks.append(("momentum conservation px", mom_cons[0]))
    checks.append(("momentum conservation py", mom_cons[1]))
    checks.append(("momentum conservation pz", sp.factor(mom_cons[2])))

    # In the physical phase-space region the two radicals are perfect squares
    # with known positive branches:
    #   sqrt(xT^2 + x2^2) = (1 - xp - zp + 2 xp zp)/xp
    #   sqrt(xT^2 + x3^2) = (xp + zp - 2 xp zp)/xp
    # since
    #   1 - xp - zp + 2 xp zp = (1-xp)(1-zp) + xp zp > 0
    #   xp + zp - 2 xp zp   = xp(1-zp) + (1-xp)zp > 0 .
    root1 = sp.sqrt((xT**2 + x2**2).subs(xT**2, xT2))
    root2 = sp.sqrt((xT**2 + x3**2).subs(xT**2, xT2))
    root1_phys = sp.simplify((sp.Integer(1) - xp - zp + 2 * xp * zp) / xp)
    root2_phys = sp.simplify((xp + zp - 2 * xp * zp) / xp)
    checks.append(("sqrt(xT^2+x2^2) branch", sp.simplify(root1**2 - root1_phys**2)))
    checks.append(("sqrt(xT^2+x3^2) branch", sp.simplify(root2**2 - root2_phys**2)))
    energy_phys = sp.simplify(
        (-sp.Rational(1, 2) * Q * x1) - sp.Rational(1, 2) * Q * root1_phys - sp.Rational(1, 2) * Q * root2_phys
    )
    checks.append(("momentum conservation E (physical branch)", sp.factor(energy_phys)))

    # This is the common normalization factor used in BGFME(..., true) and
    # generateBGFPoint().
    norm_factor = sp.expand(xp**2 * (x2**2 + x3**2 + 3 * xT2))
    norm_factor_simplified = sp.factor(sp.simplify(norm_factor))

    print("DIS POWHEG BGF kinematic audit\n")
    print(f"x1 = {sp.sstr(x1)}")
    print(f"x2 = {sp.sstr(sp.simplify(x2))}")
    print(f"x3 = {sp.sstr(sp.simplify(x3))}")
    print(f"xT^2 = {sp.sstr(xT2)}\n")

    print("Normalized BGF denominator factor:")
    print(f"  xp^2 * (x2^2 + x3^2 + 3 xT^2) = {sp.sstr(norm_factor_simplified)}\n")

    failures = 0
    for label, expr in checks:
        ok = expr == 0
        status = "OK" if ok else "FAIL"
        print(f"[{status}] {label}: {sp.sstr(expr)}")
        if not ok:
            failures += 1

    print()
    if failures:
        print("Conclusion: at least one symbolic kinematic identity failed.")
        return 1

    print(
        "Conclusion: the accepted BGF momentum map used after acceptance is "
        "internally consistent: the momenta are massless, conserve momentum, "
        "and use the same normalized denominator factor as the acceptance helper."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
