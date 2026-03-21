# Git Bootstrap and Spin-Only POWHEG Real-Emission Vertex

## Summary
It is possible to build exact `2->3` helicity amplitudes and use them only to
attach spin-density information to the POWHEG real-emission event record while
leaving the validated subtraction, hardest-emission generation, and NLO weights
unchanged.

The implementation point is after
[`DISBase::generateHardest()`](../../../src/herwig/MatrixElement/DIS/DISBase.cc#L865)
has already selected the emission and created the fresh `2->3` particles. At
that point a `HardVertex` built from an exact `2->3`
[`ProductionMatrixElement`](../../../src/herwig/MatrixElement/ProductionMatrixElement.h)
can be attached purely for spin propagation. No subtraction term, no POWHEG
acceptance kernel, and no NLO normalization formula needs to be modified.

## Implementation changes

### 1. Git bootstrap
- initialize the curated repository in this clone
- add a repository-level `.gitignore` covering generated objects, LaTeX build
  files, Herwig run artifacts, shard outputs, and campaign summaries
- make an initial baseline commit before any spin-only changes
- develop on branch `codex/powheg-real-spin-vertex`

### 2. Spin-only feature flag
- add a DIS runtime switch `UsePOWHEGRealSpinVertex`
- default value: `false`
- with the switch off, behavior must remain identical to the validated state
- with the switch on, only the event-record spin setup of the real-emission
  particles is changed

### 3. Freeze subtraction and generation
Do not modify:
- `DISBase::ComptonME()`
- `DISBase::BGFME()`
- `DISBase::generateCompton()`
- `DISBase::generateBGF()`
- `DISBase::NLOWeight()`
- the validated neutral-current coefficient algebra already fixed in
  `MENeutralCurrentDIS`

The exact `2->3` amplitudes are for spin-density propagation only.

### 4. Add exact `2->3` helicity amplitudes
Implement spin-only helpers in `MENeutralCurrentDIS` for:
- QCDC: `e q -> e q g`
- BGF: `e g -> e q qbar`

Requirements:
- return a full `ProductionMatrixElement` for the realized `2->3` kinematics
- use the corrected neutral-current couplings and crossing conventions already
  present in `MENeutralCurrentDIS`
- use ThePEG helicity wavefunction classes:
  - `SpinorWaveFunction`
  - `SpinorBarWaveFunction`
  - `VectorWaveFunction`
- use the realized particle ordering from the `RealEmissionProcess`
- keep normalization consistent enough for spin-density extraction; exact
  overall normalization matching to the scalar POWHEG kernels is not required

### 5. Attach a real-emission `HardVertex`
Add a hook in `DISBase`, default no-op:
- `constructRealEmissionSpinVertex(RealEmissionProcessPtr proc, bool isCompton)`

Override it in `MENeutralCurrentDIS` and call it at the end of
`generateHardest()` after:
- the copied lepton particles are in `incoming()/outgoing()`
- `newin`, `newout`, and `emitted` are created
- color flow has been assigned

Inside the hook:
- ensure all five `2->3` particles have `spinInfo()` using the explicit ThePEG
  `constructSpinInfo(...)` helpers
- build a `HardVertex`
- set incoming rho matrices on the real-emission incoming lepton and mapped
  incoming parton
- attach `productionVertex(...)` on all five real-emission legs:
  - incoming lepton copy
  - mapped incoming parton
  - outgoing lepton copy
  - outgoing colored leg
  - emitted colored leg
- do not alter the underlying-Born `HardVertex`

### 6. Approximate fallback
If the exact `2->3` amplitudes block, fallback v1 may:
- assign direct diagonal rho matrices to the fresh POWHEG particles
- avoid attaching a fake `productionVertex()`
- propagate only line polarizations already known from the validated code
  (`Pq_m`, `Pg_m`)

This fallback is only a contingency and is not the target implementation.

## Test plan

### Git/bootstrap checks
- confirm `.gitignore` excludes generated Herwig outputs and campaign artifacts
- confirm the initial commit exists before spin-only edits
- confirm work continues on `codex/powheg-real-spin-vertex`

### Physics invariance checks
With `UsePOWHEGRealSpinVertex = false`:
- reproduce current `testing11` / `testing12` validation numbers within
  statistics

With `UsePOWHEGRealSpinVertex = true`:
- `extract_dis_out_results.py` must show no statistically significant movement
  in the `GAMMA`, `Z`, `ALL`, and interference validation tables
- `extract_nlo_term_diagnostics.py` closures and component estimates must also
  remain statistically unchanged

This is the acceptance criterion that proves subtraction and event generation
remain untouched.

### Event-record spin checks
For POWHEG real-emission events:
- all five `2->3` particles must have `spinInfo()`
- all five must have a `productionVertex()`
- rho matrices derived from the real-emission vertex must be Hermitian,
  positive, and unit trace
- incoming rho matrices on the real-emission state must match the polarization
  used to generate the event

### Shower checks
- confirm the shower sees nontrivial rho matrices on the real-emission
  particles
- compare one spin-sensitive shower observable with feature off versus on
- hard cross sections must not move; only shower-level spin propagation may
  change

## Assumptions and defaults
- exact `2->3` amplitudes are used only for `HardVertex` / spin-density
  construction
- the scalar POWHEG real-emission machinery remains the sole source of
  subtraction and emission generation
- generated outputs, campaign summaries, and seeded shard files remain ignored
  by git
