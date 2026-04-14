# Conservative Polarized-NLO Representation Refactor

## Summary

Refactor the polarized NLO term assembly in `DISBase.cc` into one canonical, manifestly finite-looking basis without changing production behavior by default.

This is a **behavior-preserving phase-1 refactor**, not a physics rewrite. The goal is to remove the current mixed representation
- `dqRatio = dqPDF / dloPDF`
- `dgRatio = dgPDF / dloPDF`
- `Pq = Pz * dloPDF / loPDF`
- indirect assembly such as `Pq * collq_over_born_pol`

from the term assembly while preserving:
- the same PDF object sourcing
- the same clamps and floors
- the same endpoint guards
- the same default results

The refactored path should be introduced behind a switch and audited against the legacy path until old and new agree pointwise to floating-point precision.

## Motivation

The intended finite representation is already documented in `HerwigPolCodex.tex`:
- polarized terms should be proportional to `P_z Delta f_m / f_q`
- there should be no separate division by `Delta f_q(x_B)` in the final finite expression

The live code in `DISBase.cc` still uses a hybrid ratio-based representation plus guards:
- `Pq = Pz * dloPDF / loPDF`
- `deltaqOverLo = Pz * dqPDF / loPDF`
- `deltagOverLo = Pz * dgPDF / loPDF`
- `dqRatio = dqPDF / dloPDF`
- `dgRatio = dgPDF / dloPDF`
- `ratioFloor`
- `minDlo`
- clamping of `Pq`, `Pq_m`, `Pg_m`

A broader “finite-slope/probe” rewrite was attempted before `plain39` and rolled back after `plain38`, because it caused a broad regression in the shared quark-collinear sector. The production fix kept only the NC-odd repair:
- `collq_odd = qOddResponse * Pq * collq_over_born_pol`

So the safe next step is not another broad numerical rewrite. It is a conservative algebra-preserving refactor that makes the representation uniform while intentionally keeping the existing numerics.

## Files and Existing Context

Primary implementation files:
- `/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/DISBase.h`
- `/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/DISBase.cc`

Reference and validation context:
- `/Users/apapaefs/Projects/HerwigPol/HerwigPolCodex.tex`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/plain39-recovery-note.md`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/compare_polarized_nlo_terms.py`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/polarized_nlo_audit_note.md`

Key live code locations in `DISBase.cc`:
- unpolarized and polarized PDF setup around the current `loPDF`, `dqPDF`, `dgPDF`, `dloPDF`, `Pq`, `deltaqOverLo`, `deltagOverLo`
- polarized effective mapped polarizations `Pq_m`, `Pg_m`
- `collg_over_born_unpol`, `collg_over_born_pol`
- `collq_over_born_unpol`, `collq_over_born_pol`
- `collq_odd = qOddResponse * Pq * collq_over_born_pol`
- `realq`
- `realg_over_born_unpol`, `realg_over_born_pol`
- existing NLO diagnostic dumps and audit logging

## Implementation Goal

Introduce a canonical term basis for the polarized NLO assembly:
- `qRatio = qPDF / loPDF`
- `gRatio = gPDF / loPDF`
- `deltaqOverLo = Pz * dqPDF / loPDF`
- `deltagOverLo = Pz * dgPDF / loPDF`
- clamped `Pq`, `Pq_m`, `Pg_m`

Then express the polarized NLO pieces in that basis, without changing the numerical protections or production default behavior.

## Non-Goals

This patch must not:
- change the PDF object sourcing logic
- change `hadron_->pdf()` versus extractor-side sum PDF usage
- change `xp_` sampling or Jacobian logic
- change NC odd/even response formulas
- change unpolarized NLO terms
- change Born or LO behavior
- remove clamps or floors
- change default physics results

Those are separate follow-up questions.

## Public Interface Change

Add a new runtime switch on `DISBase`:
- `UseUniformPolarizedNLORepresentation`

Behavior:
- `No`:
  - keep the current production assembly
- `Yes`:
  - use the uniform canonical assembly

Default:
- `No`

Persist this switch through `DISBase` serialization.

## Detailed Implementation Steps

### 1. Add the new switch

In `DISBase`:
- add a private boolean member
- default it to `false`
- expose it as an interface switch named `UseUniformPolarizedNLORepresentation`
- document that it enables a behavior-preserving canonical refactor of the polarized NLO representation
- include it in `persistentOutput()` and `persistentInput()`

### 2. Add a small internal term bundle/helper

Create a helper struct or local bundle in `DISBase` that carries the canonical inputs:
- `loPDF`
- `qPDF`
- `gPDF`
- `dqPDF`
- `dgPDF`
- `dloPDF`
- `qRatio`
- `gRatio`
- `deltaqOverLo`
- `deltagOverLo`
- `Pq`
- `Pq_m`
- `Pg_m`
- `hasDiffPdf`
- `ratioFloor`
- `minDlo`

The helper should be internal only. No public API beyond the switch is needed.

### 3. Preserve current protection logic exactly

Do not change phase-1 numerics here.

Keep exactly as-is:
- `ratioFloor`
- `minDlo`
- the `hasDiffPdf` gating
- the conditions for zeroing terms
- clamping of `Pq`, `Pq_m`, `Pg_m`
- endpoint handling for `xp -> 1`
- any existing raw/clamped diagnostic support

This is crucial. The point is to change representation, not protection behavior.

### 4. Rewrite the polarized term assembly in the canonical basis

#### 4a. Quark collinear odd term

Replace the indirect assembly
- `collq_odd = qOddResponse * Pq * collq_over_born_pol`

with the explicit canonical form equivalent to the current math:
- define the same `k1` and `k2` kernel pieces as in the validation script
- assemble `collq_odd` as
  - `qOddResponse * [CF/xp * deltaqOverLo * k1 + CF/xp * (deltaqOverLo - Pq*xp) * k2]`

Important:
- phase 1 must preserve the current guards
- if the current code would zero the polarized quark collinear contribution because `abs(dloPDF) <= minDlo`, the new code must also zero it

#### 4b. Gluon collinear odd term

Keep the same current math, but express it directly in canonical form:
- `collg_odd` should use `deltagOverLo` rather than a mixed detour through `dgRatio` unless the legacy path is selected

The result must be algebraically identical to the current protected form.

#### 4c. Gluon real odd term

Likewise keep the same math, but express `realg_odd` directly in `deltagOverLo`.

#### 4d. Real-quark term

Do not change the math of `realq` in this phase.

Only refactor it into clearer subpieces:
- `qcdc_even`
- `qcdc_aq`

This is for readability and diagnostics only. The total `realq` must remain numerically unchanged.

### 5. Keep both paths live for A/B comparison

Implement the refactored path as an alternative branch:
- legacy path when `UseUniformPolarizedNLORepresentation No`
- uniform path when `UseUniformPolarizedNLORepresentation Yes`

For validation, also add a diagnostic mode that can compute both paths side-by-side even if only one is used for the final weight.

### 6. Add old-vs-new diagnostics

When:
- `DISDiagnostics` is on, and
- NLO term diagnostics are enabled

log the first few old/new differences for:
- `collq_odd`
- `collg_odd`
- `realg_odd`
- total `wgt`

Recommended diagnostic content:
- run name
- helicity label
- contribution label
- `xB`, `xp`, `Q2`, `mu2`
- old value
- new value
- absolute difference
- relative difference where stable

Use a dedicated diagnostic prefix so these lines are easy to grep separately from the existing `NLO_TERM_*` output.

### 7. Do not remove the legacy code in the same patch

Even if the new path matches exactly:
- keep the legacy path in place for this patch
- keep the switch defaulted to `No`

Any cleanup and deletion of the legacy representation should be a later patch after validation.

## Acceptance Criteria

The patch is acceptable only if all of the following hold.

### Build and load
- `HwMEDIS` compiles
- Herwig repository loading succeeds
- no new plugin-loading or persistency issues appear

### Default preservation
With `UseUniformPolarizedNLORepresentation No`:
- all current results remain unchanged within statistics
- the native DIS-window fix remains intact
- the finite-width spacelike-`Z` fix remains intact

### Pointwise equivalence
With diagnostics enabled:
- old and new `collq_odd` agree to floating-point precision
- old and new `collg_odd` agree to floating-point precision
- old and new `realg_odd` agree to floating-point precision
- total old/new `wgt` agreement is at numerical-noise level

### No new pathologies
With the new representation enabled:
- no increase in `NLO_TERM_PATHO`
- no new endpoint instabilities
- no new `NaN` or overflow behavior
- no regression of the `plain39` recovery pattern

## Validation Runs

### TERMDIAG pointwise checks
Run polarized TERMDIAG controls for:
- `GAMMA` POSNLO
- `GAMMA` NEGNLO
- `Z` POSNLO
- `Z` NEGNLO

Use:
- old path
- new path
- old-vs-new diagnostics enabled

Goal:
- confirm term-level equivalence before looking at campaign totals

### Small campaign A/B checks
With `UseUniformPolarizedNLORepresentation Yes`:
- rerun a small polarized NLO control set
- compare to switch-off baseline

Minimum control set:
- broad `GAMMA` polarized NLO
- broad `Z` polarized NLO
- broad `ALL` polarized NLO

Goal:
- confirm no statistical shift beyond expected noise

### Preserve known good fixes
Confirm the patch does not spoil:
- `plain44` conclusions for broad-window LO `GAMMA`
- `plain46/47` conclusions for finite-width `Z`
- current LO closure in polarized channels

## Risks

### Main risk
The representation cleanup may accidentally change the numerical order of operations in the shared quark-collinear sector and replay the `plain38` regression.

Mitigation:
- keep guards unchanged
- keep legacy path intact
- require pointwise term agreement before any campaign interpretation

### Secondary risk
The current code’s exact cancellation may depend on the existing mixed representation in ways that are not obvious from algebra alone.

Mitigation:
- compare old and new term-by-term
- do not combine this patch with any change to PDF sourcing, `xp_`, or clamp policy

## Follow-Up Work Explicitly Deferred

Do not mix these into the same patch:
- changing `hadron_->pdf()` versus extractor-side sum PDF usage
- removing `ratioFloor` or `minDlo`
- changing clamp behavior for `Pq`, `Pq_m`, `Pg_m`
- changing `xp_` parameterization
- changing the polarized coefficient-function convention
- removing the legacy representation after validation

## Working Assumptions

- The remaining polarized NLO issue is small enough that any safe refactor must begin with strict behavior preservation.
- The intended finite form in `HerwigPolCodex.tex` is the right conceptual target, but phase 1 should only move the code into that basis where it is algebraically identical to the current implementation.
- The previous broad finite rewrite failed because it changed more than representation. This patch should not repeat that.
