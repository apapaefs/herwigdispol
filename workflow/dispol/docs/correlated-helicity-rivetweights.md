# RIVETWEIGHTS Correlated Helicity Combination

This note describes the `RIVETWEIGHTS` correlated helicity-weight mode used for
neutral-current polarized DIS validation in Herwig/Rivet.

The short version is:

- generate only the unpolarized `00` event kinematics;
- attach the four helicity weights `PP`, `PM`, `MP`, and `MM` to that same
  event;
- fill the unpolarized and polarized Rivet histograms from those correlated
  weights;
- combine POSNLO and NEGNLO from already-normalized histogram bins, not from
  generic raw-accumulator YODA arithmetic.

This was introduced to reduce cancellation noise in selected polarized dijet
observables, especially low-`DQ2`, where independent helicity runs can subtract
large statistically unrelated samples.

## Motivation

The standard independent-helicity workflow estimates the polarized combination
from separate event samples:

```text
sigma    = (PP + PM + MP + MM) / 4
delta_LL = (PP + MM - PM - MP) / 4
```

For inclusive or pre-cut observables this usually behaves well. For selected
dijet observables, however, the second-jet cut can leave a relatively small
accepted population. Then the polarized combination can become a delicate
difference of large, independently fluctuating numbers.

`RIVETWEIGHTS` keeps the same helicity formulas, but evaluates all four
helicity weights on the same accepted `00` event kinematics. The event-by-event
correlation makes the cancellation much less noisy.

## Event-Level Weight Definition

For each accepted `00` POWHEG event, the DIS matrix element computes signed
target densities for the four longitudinal beam-helicity choices:

```text
PP = (Pl, Ph) = (+1, +1)
PM = (Pl, Ph) = (+1, -1)
MP = (Pl, Ph) = (-1, +1)
MM = (Pl, Ph) = (-1, -1)
```

Here `Pl` is the lepton polarization and `Ph` is the hadron polarization. The
hadron polarization enters the partonic matrix element through the polarized
PDF ratio, schematically

```text
Pq(x) = Ph * Delta f_q(x) / f_q(x).
```

The target density includes both:

- the Born matrix-element ratio for the target helicity state;
- the signed NLO density for the target helicity state.

In code this is represented schematically as

```text
rho_target = BornME2(target) * NLOWeightRaw(target)
rho_sample = BornME2(00)     * NLOWeightRaw(00)
```

The named optional event weight attached to the event is

```text
w_target = event.weight() * rho_target / abs(rho_sample)
```

The absolute value in the denominator is important because POSNLO and NEGNLO
samples are generated as separate sign pieces. The target density itself stays
signed, so the resulting `NEGNLO` optional weights can be negative.

The Herwig optional weights are:

```text
HERWIG_DIS_PP
HERWIG_DIS_PM
HERWIG_DIS_MP
HERWIG_DIS_MM
HERWIG_DIS_SIGMA
HERWIG_DIS_DELTA_LL
```

with

```text
HERWIG_DIS_SIGMA    = (PP + PM + MP + MM) / 4
HERWIG_DIS_DELTA_LL = (PP + MM - PM - MP) / 4
```

For pure photon `GAMMA`, the same four-weight representation is still available,
although only two independent lepton-hadron helicity structures are physically
needed. Keeping the same weight basis keeps the Rivet path uniform across
`GAMMA`, `Z`, and `ALL`.

## Propagation Path

The relevant implementation files are:

- `src/herwig/MatrixElement/DIS/DISBase.cc`
- `src/herwig/MatrixElement/DIS/DISBase.h`
- `src/thepeg/Config/HepMCHelper.h`
- `analyses/rivet/dis/MC_DIS_BREIT.cc`
- `workflow/dispol/scripts/run_validation_campaign.py`

The Herwig switch is:

```text
set /Herwig/MatrixElements/PowhegMEDISNCPol:GenerateRivetWeights Yes
```

The Rivet analysis option is:

```text
MC_DIS_BREIT:JETINPUT=MEPARTONS:RIVETWEIGHTS=YES
```

The event record path is:

```text
Herwig Event::optionalWeights()
  -> ThePEG HepMC3 weights with named GenRunInfo entries
  -> Rivet event weights
  -> MC_DIS_BREIT named-weight lookup
```

For HepMC3, `HepMCHelper.h` must write both the numeric weight vector and the
matching `GenRunInfo` weight names. Without names, Rivet can see a multiweight
vector but `MC_DIS_BREIT` cannot reliably identify `HERWIG_DIS_SIGMA` and
`HERWIG_DIS_DELTA_LL`.

The expected log signature is:

```text
RIVETWEIGHTS summary: events=<N>, used=<N>, missing=0
```

## MC_DIS_BREIT Filling Semantics

When `RIVETWEIGHTS=YES`, `MC_DIS_BREIT` fills:

- ordinary unpolarized histograms with `HERWIG_DIS_SIGMA`;
- polarized `D...` histograms with `HERWIG_DIS_DELTA_LL`.

The affected selected histograms include:

```text
Q2, XBj, Pt, Mjj, Eta, Zeta, pT1, pT2, pT2OverpT1, pTAsym
DQ2, DXBj, DPt, DMjj, DEta, DZeta, DpT1, DpT2, DpT2OverpT1, DpTAsym
```

The affected pre-cut histograms include:

```text
Q2PreCut, XBjPreCut, YPreCut, pT1PreCut, pT2PreCut
DQ2PreCut, DXBjPreCut, DYPreCut, DpT1PreCut, DpT2PreCut
```

The analysis does not reinterpret `TOP2PARTONS` helicity runs. `RIVETWEIGHTS`
is a separate `00`-only workflow and normally uses:

```text
JETINPUT=MEPARTONS
```

## Showered Parton-Level RIVETWEIGHTS

The campaign driver also exposes a separate `--rivetweights-shower` mode. It
keeps the same `00`-only correlated helicity weights and the same POSNLO/NEGNLO
normalized-bin combination as `--rivetweights`, but it leaves the parton shower
active.

This is intentionally separate from ordinary `--rivetweights`. The ordinary
mode continues to use the no-shower card shape and `JETINPUT=MEPARTONS`.
The showered sibling uses:

```text
MC_DIS_BREIT:JETINPUT=FULL:RIVETWEIGHTS=YES
```

The generated showered cards remove active no-shower controls such as
`ShowerHandler:LimitEmissions HardOnly`, while keeping hadronization and decays
off:

```text
set EventGenerator:EventHandler:HadronizationHandler NULL
set /Herwig/EventHandlers/EventHandler:DecayHandler NULL
```

Scale-envelope plots for `--rivetweights-shower` are rendered with blue Herwig
lines and blue scale-variation bands.

## MC_DIS_PS / SPINCOMP / SPINHAD

The same named-weight treatment is also wired into `MC_DIS_PS` for the
`main.tex` shower-spin plots. The analysis option is:

```text
MC_DIS_PS:JETINPUT=FULL:RIVETWEIGHTS=YES
```

In this mode `MC_DIS_PS` fills ordinary shower-sensitive histograms with
`HERWIG_DIS_SIGMA` and matching `D...` histograms with
`HERWIG_DIS_DELTA_LL`. This includes:

```text
NJets, pT3, pT3OverpT1, SumPtExtra, Phi3, PhiCurrentHemi, Broadening
DNJets, DpT3, DpT3OverpT1, DSumPtExtra, DPhi3, DPhiCurrentHemi, DBroadening
```

The campaign driver maps:

```text
--rivetweights --setup SPINCOMP --setup SPINHAD
```

onto 00-only weighted PS families:

```text
RIVETPSWEIGHTS-SPIN
RIVETPSWEIGHTS-NOSPIN
RIVETPSWEIGHTS-NOSPIN-UNPOL
```

These correspond to the same plot labels as the standard PS families:
`Full`, `Born-Only`, and `None`. The analyzer builds `ALL...` asymmetries and
`ALLBroadeningCumulative` directly from the correlated sigma/delta histograms
after POSNLO and NEGNLO have been added at the normalized-bin level.

Physics caveat: for shower-sensitive observables this is a correlated
hard-event helicity weighting of one showered event. It removes the
statistically independent PP/PM/MP/MM cancellation noise, but it is not a
separate helicity-conditioned shower resampling.

## POSNLO and NEGNLO Combination

This is the main practical caveat.

In the ordinary independent-helicity workflow, the campaign machinery often
builds physical NLO YODAs as:

```text
NLO = POSNLO - NEGNLO
```

That is not the right final operation for `RIVETWEIGHTS`.

In `RIVETWEIGHTS`, the optional weights already carry the signed target density.
The physical NLO result is therefore the additive combination of the signed
POSNLO and NEGNLO histograms:

```text
NLO = POSNLO + NEGNLO
```

There is a second subtlety: POSNLO and NEGNLO have different run
normalizations. Generic `yodamerge` arithmetic combines raw histogram
accumulators and then applies one merged scale. That destroys the separate
`ScaledBy` factors and can produce an incorrectly normalized result.

The observed failure mode in `rivetweights06_short` was:

- individual POSNLO `YPreCut` first bin: about `4.77e3 pb`, as expected;
- individual NEGNLO `YPreCut` first bin: about `-24 pb`, as expected;
- bad generic merged result: about `2.4e2 pb`, far too small.

The correct campaign-script behavior is to:

1. read the already-normalized POSNLO `Estimate1D` objects;
2. read the already-normalized NEGNLO `Estimate1D` objects;
3. add bin values directly;
4. combine bin errors in quadrature;
5. write the analyzed Herwig YODA from those normalized bins.

This is implemented in `build_rivetweights_objects()` via the normalized-bin
combiner rather than by using generic `yodamerge` for the final POS/NEG
analysis input.

## Campaign Usage

A short `RIVETWEIGHTS` campaign is run with:

```bash
python3 DISPOL/scripts/run_validation_campaign.py full \
  --base-dir DISPOL \
  -t rivetweightsXX_short \
  --jobs <ncores> \
  --shards <nshards> \
  --keep-going \
  --scale-variations \
  --poldis-refs-campaign <reference-tag> \
  --poldis-error-mode scale \
  --poldis skip \
  --rivetweights \
  --setup ALL \
  --force-prepare
```

The showered parton-level sibling is run the same way, replacing the flag:

```bash
python3 DISPOL/scripts/run_validation_campaign.py full \
  --base-dir DISPOL \
  -t rivetweights_showerXX_short \
  --jobs <ncores> \
  --shards <nshards> \
  --keep-going \
  --scale-variations \
  --poldis-refs-campaign <reference-tag> \
  --poldis-error-mode scale \
  --poldis skip \
  --rivetweights-shower \
  --setup ALL \
  --force-prepare
```

To repair plots for an existing campaign after a postprocessing-only fix, rerun
only:

```bash
python3 DISPOL/scripts/run_validation_campaign.py analyze-herwig \
  --base-dir DISPOL \
  -t rivetweights06_short \
  --setup ALL \
  --rivetweights \
  --scale-variations

python3 DISPOL/scripts/run_validation_campaign.py rivetplot \
  --base-dir DISPOL \
  -t rivetweights06_short \
  --setup ALL \
  --rivetweights \
  --scale-variations \
  --poldis-refs-campaign <reference-tag> \
  --poldis-error-mode scale
```

`rivetplot` does not take `--poldis skip`; use `--poldis-refs-campaign` to point
at an existing reference campaign.

## Validation Checks

Useful quick checks:

- logs contain `RIVETWEIGHTS summary: events=<N>, used=<N>, missing=0`;
- `YPreCut` has a first bin near the POLDIS scale, around `4.7e3 pb` in the
  `rivetweights06_short` ALL setup;
- `Q2PreCut`, `DQ2PreCut`, `Q2`, and `DQ2` exist in the analyzed YODA;
- selected `DQ2` is filled and should show reduced independent-helicity
  cancellation noise relative to the standard `RIVETFO` workflow.

The integrated unpolarized total in the Herwig text summary remains a separate
check on the generator-level normalization. The YODA plots check the accepted
Rivet histogram normalization.

## Current Scope and Limitations

The current `RIVETWEIGHTS` mode is:

- NC-only: `GAMMA`, `Z`, `ALL`, and the NC shower-study setups
  `SPINCOMP`/`SPINHAD`;
- `00`-only at the event-generation level;
- POSNLO/NEGNLO only, with LO skipped by default;
- a correlated fixed-order/Rivet weighting mode, not an exact alternative
  POWHEG Sudakov reweight;
- compatible with the existing scale-variation machinery, provided the
  normalized POS/NEG combination is used for each variation.

The `--rivetweights-shower` sibling is also NC-only, `00`-only, and NLO-only,
but it is a showered parton-level comparison rather than a no-shower
ME-parton comparison. It should not be mixed with `--fixed-order-powheg-no-shower`.

Charged-current support would require deciding the CC helicity basis and its
corresponding `sigma`/`delta` convention before extending the mode.
