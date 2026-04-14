# Pure-Z Estimator-Level Closure Audit

This audit starts only after the already-validated local pure-`Z` theory checks:

- couplings/signs are correct
- `Delta qbar / qbar` sign conventions are correct
- local helicity closure holds
- local `Y_- G2 + Y_+ G4` decomposition holds
- weighted integrand closure already holds through:
  - `parton`
  - `xcomb_head`
  - `sampler`

So the target here is narrower:

- find the first estimator/reporting layer where
  `sigma0 = (PP + PM + MP + MM) / 4`
  stops matching `00`

## Layers

The audited weighted layers are:

- `sampler_prebin`
  This is the already-validated `StandardEventHandler::dSigDR()` integrand from the weighted audit.
- `bin_evaluate`
  `sampler_prebin * remap_factor`
- `general_generate`
  For weighted runs:
  `bin_evaluate * reference_weight / (bias * global_max_weight) / reference_weight`
  so the net factor is `1 / (bias * global_max_weight)`.
- `eventhandler_select`
  The exact selected weight passed to `StandardEventHandler::select()` and then `XSecStat::select()`.

The cumulative estimators are:

- `sampler_run`
  The attempted-event estimator built from the same weights that feed
  `GeneralSampler::theSumWeights`, `theSumWeights2`, and `theAttempts`.
- `sampler_combined`
  The add-up-samplers path built from the same per-bin averages that feed
  `BinSampler::integratedXSec()`.
  This is emitted only when the sampler is actually using add-up mode.
- `eventhandler_generated`
  The generated-event estimator mirrored from the same
  `select / reject / reweight` bookkeeping as `XSecStat`.

## Longitudinal Basis

At every linear layer we track the same decomposition

`W = C00 + Pl * CPl + Pp * CPp + Pl * Pp * CLL`

with reconstructed helicity weights

- `W00`
- `WPP`
- `WPM`
- `WMP`
- `WMM`

and the exact closure target

`W00 = (WPP + WPM + WMP + WMM) / 4`

The logged residual is always

`closure_residual = (WPP + WPM + WMP + WMM) / 4 - W00`

with a relative normalization built from the larger of `|W00|` and the
helicity average.

## Weighted-Run Assumption

The current DIS debugging runs are treated as weighted-event runs.

If an unweighted or almost-unweighted nonlinear path is encountered, the
estimator trace is logged with

- `mode=skipped_nonlinear`

and should not be interpreted as a linear closure audit.

## Output Location

All estimator diagnostics live in the Herwig `.log` file, not the `.out` file.

- `.log`:
  `DIS_ESTIMATOR_TRACE`
  `DIS_ESTIMATOR_SUMMARY`
  `DIS_ESTIMATOR_TOTAL`
- `.out`:
  run cross-section summary only

## Decision Rule

The audit targets are:

- pointwise:
  - `bin_evaluate`
  - `general_generate`
  - `eventhandler_select`
- cumulative:
  - `sampler_run`
  - `sampler_combined`
  - `eventhandler_generated`

The first failing estimator layer wins.
