# Polarized DIS in Herwig 7

This repository contains the source code and validation machinery for our
implementation of longitudinally polarized deep-inelastic scattering (DIS) in
Herwig 7. It keeps together the modified Herwig and ThePEG source trees, the
POLDIS fixed-order reference code, the DIS Rivet analyses, the validation
scripts, and the paper source used in the study.

The preprint is in preparation. Until the arXiv entry and final citation are
available, please cite the relevant Herwig, ThePEG, Rivet, YODA, LHAPDF,
FastJet, and POLDIS references, and mention this repository if you use the code
or validation workflow.

## Status

This is research code, not a general-purpose Herwig distribution. It contains
the pieces needed to rebuild the polarized-DIS implementation and to repeat the
validation runs from source. Generated campaign output and build products are
not versioned.

The workflows kept here are the DIS NLO, POWHEG, correlated Rivet-weight, and
fixed-order/no-shower comparison paths used in the paper. Regenerate `.run`
files from the checked-in cards before starting new validation runs.

## Contents

| Path | Contents |
| --- | --- |
| `src/herwig/` | Modified Herwig 7.3.0 source tree |
| `src/thepeg/` | Modified ThePEG 2.3.0 source tree |
| `poldis/` | POLDIS source tree and DIS drivers |
| `analyses/rivet/dis/` | DIS Rivet analyses and metadata |
| `workflow/dispol/cards/` | Herwig input cards and source templates |
| `workflow/dispol/scripts/` | Campaign, analysis, conversion, and plotting tools |
| `workflow/dispol/docs/` | Workflow and validation notes |
| `docs/paper/` | Manuscript source and bibliography |

Generated libraries, `.run`, `.log`, `.out`, `.yoda`, `.csv`, `.html`, `.top`,
plot output, and rendered paper files are not tracked.

## Physics And Workflows

We implement longitudinally polarized lepton-hadron DIS with neutral-current and
charged-current exchange, NLO QCD corrections, POWHEG hardest-emission
generation, spin-density propagation, and correlated helicity weights for Rivet
analyses.

The main workflows are:

- `RIVETWEIGHTS`: correlated neutral-current helicity weights for low-noise
  differential validation.
- `RIVETFO`: POWHEG hardest-emission comparison mode with fixed-order scale and
  coupling choices.
- `RealOnly`: retained POWHEG comparison mode for isolating real-emission
  behavior.
- `FixedOrderNoShower`: QTilde mode that inserts the POWHEG real emission while
  suppressing subsequent shower evolution.

Older fixed-weight and dump-style diagnostic paths are not part of the public
workflow.

## Prerequisites

You will need, at least:

- a C++ compiler and `make`
- Autotools-compatible build tools
- Python 3
- LHAPDF and `lhapdf-config`
- Rivet and YODA, including the YODA Python module
- FastJet
- a Fortran compiler such as `gfortran` for POLDIS
- a LaTeX toolchain if rebuilding the manuscript

The exact paths are installation-dependent. The Herwig and ThePEG snapshots keep
their upstream `configure --help` output, which is still the best place to check
available build options.

## Build ThePEG And Herwig

The repository stores source trees only. One clean way to build ThePEG is:

```bash
mkdir -p build/thepeg
cd build/thepeg
../../src/thepeg/configure --prefix="$PWD/install" \
  --with-lhapdf=/path/to/lhapdf \
  --with-fastjet=/path/to/fastjet \
  --with-rivet=/path/to/rivet
make -j8
make install
```

Then build Herwig against that ThePEG installation:

```bash
cd /path/to/herwigdispol
mkdir -p build/herwig
cd build/herwig
../../src/herwig/configure --prefix="$PWD/install" \
  --with-thepeg="$PWD/../thepeg/install"
make -j8
make install
```

## Build POLDIS

The POLDIS tree includes wrappers for the DIS validation drivers:

```bash
cd poldis
./compile_dijet
./compile_dijet_rivetplots
./compile_singlejet
```

The wrappers use `gfortran` by default and honour `FC` if another Fortran
compiler is preferred. LHAPDF flags come from `lhapdf-config`.

## Build The Rivet Analyses

The DIS Rivet analyses are in `analyses/rivet/dis/`:

```bash
mkdir -p build/rivet
rivet-buildplugin build/rivet/RivetMC_DIS_BREIT.so analyses/rivet/dis/MC_DIS_BREIT.cc
rivet-buildplugin build/rivet/RivetMC_DIS_PS.so analyses/rivet/dis/MC_DIS_PS.cc
export RIVET_ANALYSIS_PATH="$PWD/build/rivet:$PWD/analyses/rivet/dis${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}"
```

The `.plot` and `.info` files are tracked. The compiled plugins are not.

## Run A Validation Workflow

Use `workflow/dispol/` as the workflow base directory. For a correlated
helicity-weight run:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py full \
  --base-dir workflow/dispol \
  -t <tag> \
  --setup ALL \
  --rivetweights \
  --jobs <ncores> \
  --shards <nshards>
```

For a fixed-order comparison run:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py full \
  --base-dir workflow/dispol \
  -t <tag> \
  --setup ALL \
  --rivetfo \
  --jobs <ncores> \
  --shards <nshards>
```

Useful follow-up stages for an existing campaign:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py postprocess --base-dir workflow/dispol -t <tag> --setup ALL
python3 workflow/dispol/scripts/run_validation_campaign.py analyze-herwig --base-dir workflow/dispol -t <tag> --setup ALL
python3 workflow/dispol/scripts/run_validation_campaign.py poldis-top --base-dir workflow/dispol --order NLO
python3 workflow/dispol/scripts/run_validation_campaign.py rivetplot --base-dir workflow/dispol -t <tag> --setup ALL
```

## Paper Reproduction Campaigns

The paper plots and tables were produced with large Tiresias campaigns. The
blocks below are the launch commands with tags and event counts exposed as shell
variables, so that they can be used either for smoke tests or for
production-statistics reruns.

Common paths:

```bash
BASE=${BASE:-/home/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL}
PYTHON=${PYTHON:-python3.10}
YODAMERGE=${YODAMERGE:-/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/bin/yodamerge}
```

Neutral-current integrated cross sections, used for the NC validation tables:

```bash
NC_TAG=${NC_TAG:-paper_nc_totals}
NC_JOBS=${NC_JOBS:-192}
NC_SHARDS=${NC_SHARDS:-200}
NC_SEED_BASE=${NC_SEED_BASE:-500000}
NC_LO_EVENTS=${NC_LO_EVENTS:-1280000000}
NC_POSNLO_EVENTS=${NC_POSNLO_EVENTS:-1280000000}
NC_NEGNLO_EVENTS=${NC_NEGNLO_EVENTS:-12800000}

"$PYTHON" "$BASE/run_validation_campaign.py" \
  campaign \
  --base-dir "$BASE" \
  -t "$NC_TAG" \
  --jobs "$NC_JOBS" \
  --shards "$NC_SHARDS" \
  --seed-base "$NC_SEED_BASE" \
  --lo-events "$NC_LO_EVENTS" \
  --posnlo-events "$NC_POSNLO_EVENTS" \
  --negnlo-events "$NC_NEGNLO_EVENTS" \
  --progress-interval 2 \
  --max-listed 40 \
  --force-prepare \
  --keep-going \
  --yoda-merge-tool "$YODAMERGE"
```

Charged-current integrated cross sections:

```bash
CC_TAG=${CC_TAG:-paper_cc_totals}
CC_JOBS=${CC_JOBS:-192}
CC_SHARDS=${CC_SHARDS:-200}
CC_SEED_BASE=${CC_SEED_BASE:-500000}
CC_LO_EVENTS=${CC_LO_EVENTS:-1280000000}
CC_POSNLO_EVENTS=${CC_POSNLO_EVENTS:-1280000000}
CC_NEGNLO_EVENTS=${CC_NEGNLO_EVENTS:-12800000}
CC_POLDIS_EVENTS=${CC_POLDIS_EVENTS:-200000000}

"$PYTHON" "$BASE/run_validation_campaign.py" \
  campaign \
  --base-dir "$BASE" \
  -t "$CC_TAG" \
  --setup CC \
  --pdf-profile nnpdf_paired \
  --poldis run \
  --poldis-events "$CC_POLDIS_EVENTS" \
  --poldis-jobs 4 \
  --poldis-variant-jobs 1 \
  --jobs "$CC_JOBS" \
  --shards "$CC_SHARDS" \
  --seed-base "$CC_SEED_BASE" \
  --lo-events "$CC_LO_EVENTS" \
  --posnlo-events "$CC_POSNLO_EVENTS" \
  --negnlo-events "$CC_NEGNLO_EVENTS" \
  --progress-interval 2 \
  --max-listed 40 \
  --force-prepare \
  --keep-going \
  --yoda-merge-tool "$YODAMERGE"
```

Correlated-weight no-shower validation, used for the parton-level differential
plots:

```bash
NOSHOWER_TAG=${NOSHOWER_TAG:-paper_noshower}
NOSHOWER_POLDIS_REFS_TAG=${NOSHOWER_POLDIS_REFS_TAG:-rivetweights_noshower02}
NOSHOWER_JOBS=${NOSHOWER_JOBS:-192}
NOSHOWER_SHARDS=${NOSHOWER_SHARDS:-200}
NOSHOWER_SEED_BASE=${NOSHOWER_SEED_BASE:-920000}
NOSHOWER_POSNLO_EVENTS=${NOSHOWER_POSNLO_EVENTS:-40000000}
NOSHOWER_NEGNLO_EVENTS=${NOSHOWER_NEGNLO_EVENTS:-2000000}

"$PYTHON" "$BASE/run_validation_campaign.py" \
  full \
  --base-dir "$BASE" \
  -t "$NOSHOWER_TAG" \
  --jobs "$NOSHOWER_JOBS" \
  --shards "$NOSHOWER_SHARDS" \
  --seed-base "$NOSHOWER_SEED_BASE" \
  --posnlo-events "$NOSHOWER_POSNLO_EVENTS" \
  --negnlo-events "$NOSHOWER_NEGNLO_EVENTS" \
  --progress-interval 2 \
  --max-listed 20 \
  --keep-going \
  --scale-variations \
  --poldis skip \
  --poldis-refs-campaign "$NOSHOWER_POLDIS_REFS_TAG" \
  --poldis-error-mode scale \
  --rivetweights \
  --setup ALL \
  --fixed-order-powheg-no-shower \
  --force-prepare \
  --yoda-merge-tool "$YODAMERGE"
```

Shower-spin comparison, used for the spin-comparison plots with and without
hadronization:

```bash
SPIN_TAG=${SPIN_TAG:-paper_spinweights}
SPIN_JOBS=${SPIN_JOBS:-192}
SPIN_SHARDS=${SPIN_SHARDS:-400}
SPIN_SEED_BASE=${SPIN_SEED_BASE:-400000}
SPIN_POSNLO_EVENTS=${SPIN_POSNLO_EVENTS:-320000000}
SPIN_NEGNLO_EVENTS=${SPIN_NEGNLO_EVENTS:-32000000}

"$PYTHON" "$BASE/run_validation_campaign.py" \
  full \
  --base-dir "$BASE" \
  -t "$SPIN_TAG" \
  --jobs "$SPIN_JOBS" \
  --shards "$SPIN_SHARDS" \
  --seed-base "$SPIN_SEED_BASE" \
  --posnlo-events "$SPIN_POSNLO_EVENTS" \
  --negnlo-events "$SPIN_NEGNLO_EVENTS" \
  --progress-interval 1 \
  --max-listed 40 \
  --keep-going \
  --poldis skip \
  --rivetweights \
  --setup SPINCOMP \
  --setup SPINHAD \
  --force-prepare \
  --yoda-merge-tool "$YODAMERGE"
```

The scripts write campaign output under:

- `workflow/dispol/campaigns/<tag>/`
- transient `.run`, `.log`, `.out`, `.yoda`, `.csv`, `.html`, and `.tex`
  products under `workflow/dispol/`

More workflow detail is in `workflow/dispol/docs/validation-workflow.md` and
`workflow/dispol/docs/correlated-helicity-rivetweights.md`.

## Tests

Run the workflow tests with:

```bash
python3 -m unittest discover workflow/dispol/scripts/tests
```

For source-level checks, also verify that the active DIS cards and scripts do
not regenerate retired validation-only settings, then rebuild the touched Herwig
modules in the target environment.

## Paper Source

The manuscript source lives in:

- `docs/paper/HerwigPolCodex.tex`
- `docs/paper/biblio.bib`

Rendered PDFs, auxiliary LaTeX products, and generated figures are excluded from
git. Citation metadata will be added once the preprint is available.

## Acknowledgments

We thank Daniel de Florian for useful discussions and for providing the POLDIS
code used for validation. AP acknowledges support from the US Department of
Energy, Office of Science, Office of Nuclear Physics under Award Number
DE-SC0025728.

## License

The modified Herwig 7 and ThePEG 2.3.0 source snapshots in this repository are
distributed under the GNU General Public License, version 3 (GPLv3). A top-level
copy is provided in `COPYING`, with component-level copies in:

- `src/herwig/COPYING`
- `src/thepeg/COPYING`

Additional component-specific notices and exceptions remain in the relevant
subdirectories.
