# DISPOL Run Area

`DISPOL/` is the public run and validation area for the polarized DIS study.
It collects the inputs needed to generate Herwig DIS samples, analyze them with
Rivet, and compare the resulting observables to fixed-order POLDIS references.

## Contents

- `cards/`: Herwig input cards for LO, POWHEG, helicity, electroweak-channel,
  Rivet, fixed-order-like Rivet, parton-shower, and HepMC-weight checks.
- `analyses/rivet/dis/`: Rivet analysis sources for `MC_DIS_BREIT` and
  `MC_DIS_PS`.
- `scripts/`: validation and conversion drivers, including
  `run_validation_campaign.py`.
- `docs/`: workflow notes and audit notes.  Some notes retain historical
  diagnostic windows; the current standard cards use `100 < Q^2 < 2500 GeV^2`.

Generated campaign products are written under `DISPOL/campaigns/` and are
ignored by Git.

## Physics Conventions

The DIS cards generate physical cross sections for specified lepton and proton
helicities.  The conventional unpolarized and longitudinal double-spin
combinations are built as

```text
sigma    = (PP + PM + MP + MM) / 4
Delta_LL = (PP + MM - PM - MP) / 4
```

The `RIVETWEIGHTS` cards evaluate the four helicity configurations on the same
phase-space point and expose them as named event weights.  In that mode,
`MC_DIS_BREIT` and `MC_DIS_PS` fill ordinary histograms with
`HERWIG_DIS_SIGMA` and polarized `D...` histograms with
`HERWIG_DIS_DELTA_LL`.

Use `--rivetweights-shower` for the separate parton-shower-on,
hadronization-off sibling of the ordinary `--rivetweights` campaign mode.

## Prerequisites

Use a shell where the modified Herwig installation is active and where these
tools can be found:

- `Herwig`
- `rivet`, `rivet-build`, `rivet-mkhtml`, and YODA tools
- Python with `yoda`
- LHAPDF with the PDF sets referenced by the chosen cards
- the POLDIS compiler/runtime setup if POLDIS reference generation is requested

Build the Rivet analyses from the repository root:

```bash
rivet-build RivetMC_DIS.so DISPOL/analyses/rivet/dis/MC_DIS_BREIT.cc DISPOL/analyses/rivet/dis/MC_DIS_PS.cc
export RIVET_ANALYSIS_PATH="$PWD:${RIVET_ANALYSIS_PATH}"
```

## Driver Usage

The main driver is:

```bash
python3 DISPOL/scripts/run_validation_campaign.py --help
```

Pass `--base-dir DISPOL` in normal use.  Common setup names are `GAMMA`, `Z`,
`ALL`, `CC`, `SPINCOMP`, `SPINHAD`, and `SPINVAL`; availability depends on the
selected run mode and card set.

Create a small Herwig/RivetFO campaign:

```bash
python3 DISPOL/scripts/run_validation_campaign.py campaign \
  --base-dir DISPOL \
  -t test-dis \
  --setup ALL \
  --rivetfo \
  --jobs 4 \
  --shards 1
```

Analyze the output:

```bash
python3 DISPOL/scripts/run_validation_campaign.py postprocess --base-dir DISPOL -t test-dis --setup ALL --rivetfo
python3 DISPOL/scripts/run_validation_campaign.py analyze-herwig --base-dir DISPOL -t test-dis --setup ALL --rivetfo
python3 DISPOL/scripts/run_validation_campaign.py rivetplot --base-dir DISPOL -t test-dis --setup ALL --rivetfo
```

Run the correlated named-weight path:

```bash
python3 DISPOL/scripts/run_validation_campaign.py full \
  --base-dir DISPOL \
  -t test-weights \
  --setup ALL \
  --rivetweights \
  --jobs 4 \
  --shards 1
```

## Direct HepMC Weight Check

The card below enables the HepMC writer and named helicity weights:

```bash
Herwig read DISPOL/cards/DIS-POL-POWHEG_00-POSNLO-ALL-HEPMCWEIGHTS.in
Herwig run DIS-POL-POWHEG_00-POSNLO-ALL-HEPMCWEIGHTS.run -N 100 -t hepmcweights
```

The resulting HepMC file should advertise:

```text
HERWIG_DIS_PP
HERWIG_DIS_PM
HERWIG_DIS_MP
HERWIG_DIS_MM
HERWIG_DIS_SIGMA
HERWIG_DIS_DELTA_LL
```

The same weights are consumed by the custom Rivet analyses when
`RIVETWEIGHTS=YES` is present in the analysis option string.
