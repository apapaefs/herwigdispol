# Hard-process spin transfer to the angular-ordered shower

`QTildeShowerHandler:HardProcessSpin` is independent of `SpinCorrelations`.
It defaults to `Yes`. The off setting is an **LHE-like information-loss
control**, not an assertion of equivalence to every Les Houches setup.

| HardProcessSpin | SpinCorrelations | Treatment |
|---|---|---|
| Yes | Yes | Existing polarized hard-to-shower transfer and shower correlations |
| No | Yes | No hard-spin transfer or beam-polarized ISR conditioning; shower-generated correlations retained |
| Yes | No | Existing shower-spin-off control; polarized ISR conditioning is retained |
| No | No | Neither hard-spin transfer nor shower-generated spin correlations |

```text
set /Herwig/Shower/ShowerHandler:HardProcessSpin No
set /Herwig/Shower/ShowerHandler:SpinCorrelations Yes
```

## Scope and implementation

The first supported off-mode use is native LO `pp -> jj`, with incoming and
outgoing partons restricted to u,d,s,c,b and gluons, without MPI, hard resonance
decays, matrix-element shower corrections, matching or merging. Unsupported
configurations fail at shower initialization. DIS off mode is deliberately
unsupported. Other shower handlers keep their existing policy.

The base handler's virtual policy hook gates only the initialization of shower
beam polarization and longitudinal/transverse difference-PDF references.
Resetting precedes this gate every event. The hard polarized extractor and
matrix element are not modified. Consequently the hard cross section, event
weight, kinematics, flavours and colour assignment remain polarized.

`ShowerTree` detaches the entire inherited `SpinInfo` from its own hard-parton
copies, not from the original hard event. It does not just zero off-diagonal
density-matrix elements. Existing physical-state initialization then gives
quarks I2/2 and massless gluons diag(1/2,0,1/2). A rejected shower is reset to
clean progenitors; nothing clears the newly generated correlations during an
accepted shower.

The QTilde persistent class version is 1. Version-0 objects with the immediately
preceding supported field layout load with `HardProcessSpin Yes`. This is not a
promise to read arbitrary historical fork layouts. In particular, an older
Odysseus version-0 fork lacks the canonical `_powhegEmissionMode` field while
using the same version number; its saved payload is not a supported predecessor.
New campaigns must regenerate
their `.run` files and record the new setting explicitly. Rebuild all affected
Herwig components and verify runtime library loading; no ThePEG source change
is required. Keep the previous production installation untouched.

## Validation and limits

`../scripts/validation/run_hard_spin_validation.py` compiles a validation-only
plugin and creates a fresh output directory. It supports native full/off runs,
old-runtime default regression, a polarized LO DIS regression, and ordinary LHE
replay of captured unweighted 510 GeV hard events. It is never loaded by campaign
cards. Fixed-input fixtures check physical seeds, detached tensors, retry reset,
both-beam terminal matrices and PDF veto behavior. A PDF wrapper counts hard and
shower calls without altering the PDF value.

The counter distinguishes hard-extractor stack frames after an event veto:
the global current-shower pointer can outlive a rejected event, and that pointer
alone is not evidence that backward evolution called a polarized PDF.

The LHE export uses the standard four-parton record, `SPINUP=9`, positive unit
weights and captured colour connections/SCALUP. The replay uses an ordinary
extractor, `ReweightPDF No`, identical shower PDFs, explicit partner/recoil/scale
settings, MPI off, hadronization and decays. The native off mode is not enabled
on the imported LHE reference: the ordinary reader already lacks the hard QCD
tensor. See the [LHE format specification](https://arxiv.org/abs/hep-ph/0609017).

For validation, `WeightOption VarWeight` avoids the ordinary handler's automatic
`skipEvents()` resampling. Positive unit event weights are still required by the
audit. Do not pair accepted events by ordinal if an event was vetoed: the strict
ordered comparison intentionally refuses to certify such a sample. Use
`compare_hard_spin_validation.py --by-identity` for finite-file replay: it checks
unique hard inputs, treats a vetoed input as a zero-rate event, and retains all
exported hard inputs in the denominator. The exported hard cross section times
this rate includes the observed shower-veto acceptance. More than 0.1% missing
inputs fails closure. Duplicate or unknown inputs always fail.
`inspect_lhe_identity.py` diagnoses skipped/duplicated identities but cannot
approve closure. Unit tests cover the veto-aware estimator and covariance.

`compare_hard_spin_validation.py` checks each hard-event identity before computing
paired differences. It includes shared-hard-event covariance and reports the
achieved precision; small p-values or failure to resolve a difference are not
used as an equivalence criterion. Its predeclared rate gate is a simultaneous
95% upper bound of 5% relative on the baseline and third-jet rates above 2 and
3 GeV. Rare rates are reported, not silently accepted. This certifies only the
listed rate tests, not untested angular observables or generic LHE generators.
Do not promote the runtime or start the campaign pilot while interface closure
is unresolved. Production families use independent samples and errors, unlike
this paired validation.

Run results, exact library fingerprints, supported persistence checks and release
status belong in the project validation record, not inferred from source tests.
