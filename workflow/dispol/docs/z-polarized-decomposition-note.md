# Pure-Z Polarized Decomposition Audit

This note defines the second local pure-`Z` audit, after the antiquark sign convention itself was verified.

The question is narrower now:

Do Herwig's local pure-`Z` polarized coefficients

- `z_pq_coeff`
- `z_plpq_coeff`

match the same coefficient reconstructed from the POLDIS `Y_- G2 + Y_+ G4` basis?

The target is a **local dimensionless coefficient identity**. The production cross section is not being retested here.

## Scope and normalization

We work in pure `Z` only:

- Herwig: `_gammaZ == 2`
- POLDIS: `EWFLAG == 2`

The trace-level C++ audit keeps the full event-local pure-`Z` factor `zSq`, including the propagator dependence already present in `ncCoefficients(...)`.

The offline analytic script does not take `Q^2` as input, so its analytic table checks the same identity with the common pure-`Z` prefactor normalized to the electroweak piece only. This is sufficient for the algebraic equivalence check because the same `zSq` multiplies both sides.

The global factors that are not part of Herwig's stored local coefficient are stripped:

- `W * 2*pi*alpha^2 / (y Q^4) * NRM/NEV * X`

The signed flavour dependence is not stripped.

## Herwig basis

In `/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/MENeutralCurrentDIS.cc`, the pure-`Z` coefficient multiplying the incoming hadron-side quark polarization is

`z_pq_coeff = (1 + ell^2) * (Dq + Pl * Dlq) + ell * (Nq + Pl * Nlq)`

and the `Pl * Pq` coefficient is

`z_plpq_coeff = (1 + ell^2) * Dlq + ell * Nlq`

with

- `Dq  = -2 * etaQ * CAq * CVq * (CAl^2 + CVl^2) * zSq`
- `Dlq =  4 * etaL * etaQ * CAl * CAq * CVl * CVq * zSq`
- `Nq  = -4 * etaL * CAl * CVl * (CAq^2 + CVq^2) * zSq`
- `Nlq =  2 * (CAq^2 + CVq^2) * (CAl^2 + CVl^2) * zSq`

## POLDIS basis

In `/Users/apapaefs/Projects/HerwigPol/POLDIS/POLDIS-public/poldis.f`, the polarized inclusive pure-`Z` weight is built as

`Y_- * G2 + Y_+ * G4`

with

- `G2LO = (DELTA1 + SIGMA1) * X`
- `G4LO = (DELTAW1 + SIGMAW1) * X`

and for a fixed signed flavour `I`

- `COMBPDFVEC ~ F1(I) * (CVZ(I)^2 + CAZ(I)^2) * (CVZE^2 + CAZE^2)`
- `COMBPDFAX  ~ F1(I) * 4 * CVZ(I) * CAZ(I) * CVZE * CAZE`

For the current DIS setup, the signed antiquark dependence enters through

- `F1(I)` with signed `I`
- `CVZ(-I) = -CVZ(I)`
- `CAZ(-I) =  CAZ(I)`

## Local coefficient mapping

Define

- `L_even = CVl^2 + CAl^2`
- `L_ax   = CVl * CAl`
- `Q_even = CVq^2 + CAq^2`
- `Q_ax^s = etaQ * CVq * CAq`

Then the stripped POLDIS-like pure-`Z` coefficients in the same local normalization as Herwig are

`z_poldis_g2_coeff = 2 * zSq * Q_even * (Pl * L_even - 2 * etaL * L_ax) / y^2`

`z_poldis_g4_coeff = -4 * zSq * Q_ax^s * (L_even - 2 * etaL * Pl * L_ax) / y^2`

and therefore

`z_poldis_pq_coeff = Y_- * z_poldis_g2_coeff + Y_+ * z_poldis_g4_coeff`

Using

- `Y_+ / y^2 = (1 + ell^2) / 2`
- `Y_- / y^2 = ell`

this becomes

`z_poldis_pq_coeff`

`= ell * 2 * zSq * Q_even * (Pl * L_even - 2 * etaL * L_ax)`

`  + (1 + ell^2)/2 * (-4) * zSq * Q_ax^s * (L_even - 2 * etaL * Pl * L_ax)`

which is algebraically identical to

`(1 + ell^2) * (Dq + Pl * Dlq) + ell * (Nq + Pl * Nlq)`

after substituting the Herwig `D/N` coefficients above.

The corresponding `Pl * Pq` coefficient from the POLDIS basis is

`z_poldis_plpq_coeff`

`= Y_- * [2 * zSq * Q_even * L_even / y^2]`

`  + Y_+ * [8 * etaL * zSq * Q_ax^s * L_ax / y^2]`

which is algebraically identical to

`(1 + ell^2) * Dlq + ell * Nlq`

## Code-level identity to check

The C++ audit logs the following pairs:

- `z_poldis_pq_coeff` versus `z_pq_coeff`
- `z_poldis_plpq_coeff` versus `z_plpq_coeff`

and the residuals

- `z_pq_residual = z_pq_coeff - z_poldis_pq_coeff`
- `z_plpq_residual = z_plpq_coeff - z_poldis_plpq_coeff`

The exact local identities that should hold event by event are

- `z_poldis_pq_coeff == z_pq_coeff`
- `z_poldis_plpq_coeff == z_plpq_coeff`

If these fail for traced pure-`Z` events, the remaining issue is in the local coefficient assembly itself. If they pass, the next discrepancy must live downstream of this local decomposition.

## Offline checker

`/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/verify_z_polarized_decomposition.py`

checks the same identity in two ways:

The trace replay uses `DIS_BOOKKEEPING_TRACE` lines taken from the Herwig
`.log` files. The `.out` files remain the run cross-section summaries and are
not the source of these local bookkeeping diagnostics.

1. analytic flavour table for `qid = +/- 1..5`
2. trace replay from `DIS_BOOKKEEPING_TRACE` lines

The trace replay is the authoritative check for the full event-local propagator-dependent coefficient.
