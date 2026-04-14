# Pure-Z Weighted Integrand Audit

This note defines the next pure-`Z` audit after the local coefficient checks.

The local algebra is already validated:

- pure-`Z` couplings are correct
- antiquark sign conventions are correct
- local helicity closure is exact
- the local `Y_- G2 + Y_+ G4` decomposition is exact in both `me2()` and `NLOWeight()`

So the remaining question is no longer the local formula. It is:

At what **weighted** layer does

`W00 = (WPP + WPM + WMP + WMM) / 4`

first fail?

## Weighted layers

The audit traces three layers:

- `parton`
  This is the partonic `dSigHatDR` layer inside the DIS matrix element.
- `xcomb_head` or `xcomb_group_total`
  This is the cross section after PDF weights, cuts, and any reweight factors in
  `StandardXComb` / `StdXCombGroup`.
- `sampler`
  This is the full sampled integrand after multiplying by the luminosity
  Jacobian and luminosity-function value in `StandardEventHandler`.

Each layer is written in the longitudinal basis

`W = C00 + Pl * CPl + Pp * CPp + Pl * Pp * CLL`

with the explicit helicity reconstructions

- `W00`
- `WPP`
- `WPM`
- `WMP`
- `WMM`

and the closure residual

`Delta = (WPP + WPM + WMP + WMM) / 4 - W00`

## NLO raw vs accepted

For the NLO correction there are two parton-level channels:

- `mode=raw`
  Uses the signed correction from `wgt`
- `mode=run`
  Uses the accepted contribution actually generated in the run:
  - `POSNLO`: `max(0, wgt)`
  - `NEGNLO`: `-max(0, -wgt)`

Only `mode=run` is propagated to the `xcomb_*` and `sampler` layers, because
those are the quantities that actually enter the generated cross section.

## What the trace means

The audit writes:

- `DIS_WEIGHTED_TRACE`
- `DIS_WEIGHTED_SUMMARY`

to the Herwig `.log` files.

The `.out` files remain the run cross-section summaries and are not the source
of these weighted-integrand diagnostics.

The intended interpretation is:

- if `parton` already fails, the problem is still inside the DIS weighted ME assembly
- if `parton` passes but `xcomb_head` fails, the problem is in PDF/cut/reweight application
- if `xcomb_head` passes but `xcomb_group_total` fails, the problem is in grouped dependent accumulation
- if `xcomb_group_total` passes but `sampler` fails, the problem is in luminosity folding
- if all layers pass, the remaining suspect is downstream extraction or reporting

## Offline checker

`/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/verify_weighted_integrand_closure.py`

replays `DIS_WEIGHTED_TRACE` records from the `.log` files and checks:

1. per-layer closure:
   `Delta == (WPP + WPM + WMP + WMM)/4 - W00`
2. `xcomb_head` propagation from the traced `parton` run record using
   `pdf_weight * cut_weight * ckkw_factor * me_reweight_factor`
3. `sampler` propagation from the latest available `xcomb_*` record using
   `lumi_jac * lumi_value`

The first failing layer is the point where the pure-`Z` closure should be
considered broken.
