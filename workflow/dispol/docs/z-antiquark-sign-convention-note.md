# Pure-Z Antiquark Sign Convention Audit

This note isolates one narrow question:

For pure `Z` (`Herwig _gammaZ == 2`, `POLDIS EWFLAG == 2`), do Herwig and POLDIS use the same sign convention for the incoming antiquark polarization ratio `Delta qbar / qbar`?

The answer should be checked at the level where the antiquark sign can actually enter:

1. the signed flavour label passed to the PDFs
2. the quark-line ordering in the helicity amplitude
3. the pure-`Z` axial carrier that is odd under `q <-> qbar`

This note does not reopen the mirror-slot question. In the traced DIS path, the incoming ordering is lepton-first, parton-second, and the sampled DIS traces are non-mirrored.

## POLDIS

In `/Users/apapaefs/Projects/HerwigPol/POLDIS/POLDIS-public/poldis.f`, the pure-`Z` fermion couplings are defined by flavour sign:

- `CVZ(-I) = -CVZ(I)`
- `CAZ(-I) =  CAZ(I)`

This is the key antiquark rule. The vector coupling is odd, the axial coupling is even.

The pure-`Z` axial carrier then enters in two places:

1. `COMBPDFAX`:

   `F(I) * PROPZ * (4 * CVZ(I) * CAZ(I) * CVZE * CAZE)`

2. `COMBINEQI`:

   `CSAX * PROPZ * (4 * CVZ(I) * CAZ(I) * CVZE * CAZE)`

So for a signed flavour `I`, the POLDIS antiquark sign lives in:

- the signed PDF slot `F(I)`
- the sign of `CVZ(I) * CAZ(I)`

For antiquarks, `CVZ(I) * CAZ(I)` flips sign because only `CVZ` flips.

At LO, the polarized pure-`Z` inclusive weight is assembled from:

- `G2LO = (DELTA1 + SIGMA1) * X`
- `G4LO = (DELTAW1 + SIGMAW1) * X`
- `WWW(0) ~ Y_- * G2LO + Y_+ * G4LO`

with `Y_+ = 1 + (1-y)^2` and `Y_- = 1 - (1-y)^2`.

The antiquark sign test therefore sits in the `G4` / `COMBPDFAX` piece. The `G2` / `COMBPDFVEC` piece is charge-even and is not the sign-sensitive part of this audit.

## Herwig

In `/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/MENeutralCurrentDIS.cc`, the incoming quark line is ordered by the sign of the actual signed incoming flavour:

- `qorder = (qin->id() > 0)` in `bornMEForMomenta()`
- the trace now also records `qid` and `qorder`

The signed flavour entering the DIS PDFs is the same signed flavour that labels the incoming parton:

- `mePartonData()[1]` is the real incoming signed flavour
- `DISBase.cc` uses that exact signed id in the `xfx()` calls for
  - `loPDF`
  - `dloPDF`
  - `qPDF`
  - `dqPDF`

For the pure-`Z` audit payload, Herwig stores:

- `qid`
- `qorder`
- `etaQ = sign(qid)`
- family couplings `CVq`, `CAq`
- the signed axial carrier `etaQ * CVq * CAq`

This is exactly the Herwig analogue of the POLDIS object `CVZ(I) * CAZ(I)`.

Because Herwig keeps the quark/antiquark sign in `etaQ` rather than in `CVq` itself, the mapping is:

`CVZ(I) * CAZ(I)  <->  etaQ * CVq * CAq`

for the same signed incoming flavour `I`.

## Analytic Mapping

Let `f > 0` denote the flavour family, and let `I = +/- f` be the signed incoming flavour.

POLDIS:

- `CVZ(+f) =  CV_f`
- `CVZ(-f) = -CV_f`
- `CAZ(+/-f) = CA_f`

So:

- `CVZ(+f) * CAZ(+f) =  CV_f * CA_f`
- `CVZ(-f) * CAZ(-f) = -CV_f * CA_f`

Herwig:

- `CVq = CV_f`
- `CAq = CA_f`
- `etaQ = +1` for `q`, `-1` for `qbar`

So:

- `etaQ * CVq * CAq = +CV_f * CA_f` for `q`
- `etaQ * CVq * CAq = -CV_f * CA_f` for `qbar`

Therefore the signed pure-`Z` axial carrier is analytically identical in the two codes.

## Decision Table

The table below states where the antiquark sign is supposed to come from.

| `qid` | family | `qorder` in Herwig | PDF flavour id used | sign source in POLDIS | sign source in Herwig |
| --- | --- | --- | --- | --- | --- |
| `+1` | down | `1` | `+1` | `CVZ(+1) * CAZ(+1)` | `etaQ=+1`, `CVq*CAq` |
| `-1` | down | `0` | `-1` | `CVZ(-1) * CAZ(-1)` and `F(-1)` | `etaQ=-1`, `CVq*CAq`, PDF id `-1` |
| `+2` | up | `1` | `+2` | `CVZ(+2) * CAZ(+2)` | `etaQ=+1`, `CVq*CAq` |
| `-2` | up | `0` | `-2` | `CVZ(-2) * CAZ(-2)` and `F(-2)` | `etaQ=-1`, `CVq*CAq`, PDF id `-2` |
| `+3` | down | `1` | `+3` | `CVZ(+3) * CAZ(+3)` | `etaQ=+1`, `CVq*CAq` |
| `-3` | down | `0` | `-3` | `CVZ(-3) * CAZ(-3)` and `F(-3)` | `etaQ=-1`, `CVq*CAq`, PDF id `-3` |
| `+4` | up | `1` | `+4` | `CVZ(+4) * CAZ(+4)` | `etaQ=+1`, `CVq*CAq` |
| `-4` | up | `0` | `-4` | `CVZ(-4) * CAZ(-4)` and `F(-4)` | `etaQ=-1`, `CVq*CAq`, PDF id `-4` |
| `+5` | down | `1` | `+5` | `CVZ(+5) * CAZ(+5)` | `etaQ=+1`, `CVq*CAq` |
| `-5` | down | `0` | `-5` | `CVZ(-5) * CAZ(-5)` and `F(-5)` | `etaQ=-1`, `CVq*CAq`, PDF id `-5` |

## What The Checker Verifies

`/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/verify_z_antiquark_sign.py` checks:

The trace input for this checker comes from the Herwig `.log` files, where
`DIS_BOOKKEEPING_TRACE` and `DIS_BOOKKEEPING_SUMMARY` are written. The
corresponding `.out` files remain the run-summary cross-section artifacts and
are not the source of these bookkeeping diagnostics.

1. analytic flavour-by-flavour agreement of
   - signed PDF id
   - signed axial carrier `CVZ(I) * CAZ(I)` vs `etaQ * CVq * CAq`
   - normalized pure-`Z` `G4` coefficient
2. trace-level agreement of
   - `qid`
   - `qorder`
   - `cur_pdf_qid`
   - `nc_signed_cvca`
   - `z_signed_axial_coupling`

The script also includes two negative controls:

- flip the antiquark PDF sign/slot
- flip the antiquark coupling sign

Both must fail for antiquark rows, otherwise the audit would be tautological.

## Scope Boundary

This note proves only the antiquark sign convention. It does not prove that the full Herwig pure-`Z` polarized observable equals the full POLDIS observable point by point, because the full observable mixes charge-even and charge-odd pieces and depends on the exact lepton-polarization convention used in each code.

The sign-sensitive statement proved here is narrower and cleaner:

For pure `Z`, the antiquark sign is supposed to enter through the same signed flavour id and the same odd axial carrier in both POLDIS and Herwig.
