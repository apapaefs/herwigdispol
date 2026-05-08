# Polarized DIS in Herwig 7

This repository contains the release-facing source and reproducibility material
for a polarized deep-inelastic scattering (DIS) implementation in Herwig 7. It
collects the modified Herwig and ThePEG source snapshots, the curated POLDIS
fixed-order reference code, DIS Rivet analyses, validation workflow scripts, and
paper source used for the associated polarized-DIS study.

The associated pre-print is under preparation. Until the pre-print identifier
and citation information are available, please cite the relevant upstream
Herwig, ThePEG, Rivet, YODA, LHAPDF, FastJet, and POLDIS references used by your
study, and mention this repository when using the code or validation workflow.

## Status

This is a research software release snapshot, not a general-purpose Herwig
distribution. The repository is intended to make the DIS implementation and
validation workflow reproducible from source, while excluding generated
campaign products and build artifacts.

The release-facing workflow keeps the validated DIS NLO, POWHEG, correlated
Rivet-weight, and fixed-order/no-shower comparison paths. Generated `.run`
files should be regenerated from the checked-in source cards before validation
runs.

## Contents

| Path | Contents |
| --- | --- |
| `src/herwig/` | Modified Herwig 7.3.0 source snapshot |
| `src/thepeg/` | Modified ThePEG 2.3.0 source snapshot |
| `poldis/` | Curated POLDIS source tree and DIS drivers |
| `analyses/rivet/dis/` | DIS Rivet analyses and metadata |
| `workflow/dispol/cards/` | Herwig input cards and source templates |
| `workflow/dispol/scripts/` | Campaign, analysis, conversion, and plotting tools |
| `workflow/dispol/docs/` | Workflow and validation notes |
| `docs/paper/` | Manuscript source and bibliography |

Generated products such as compiled libraries, `.run`, `.log`, `.out`, `.yoda`,
`.csv`, `.html`, `.top`, plot outputs, and rendered paper files are intentionally
not tracked.

## Physics And Workflow Scope

The implementation supports longitudinally polarized lepton-hadron DIS with
neutral-current and charged-current exchange, NLO QCD corrections, POWHEG
hardest-emission generation, spin-density propagation, and correlated helicity
weights for Rivet analyses.

The main release workflows are:

- `RIVETWEIGHTS`: correlated neutral-current helicity weights for low-noise
  differential validation.
- `RIVETFO`: POWHEG hardest-emission comparison mode with fixed-order scale and
  coupling choices.
- `RealOnly`: retained POWHEG comparison mode for isolating real-emission
  behavior.
- `FixedOrderNoShower`: QTilde mode that inserts the POWHEG real emission while
  suppressing subsequent shower evolution.

Historical source paths for validation-only fixed-weight or dump-style studies
are not part of the release-facing workflow.

## Prerequisites

Typical local builds require:

- a C++ compiler and `make`
- Autotools-compatible build tools
- Python 3
- LHAPDF and `lhapdf-config`
- Rivet and YODA, including the YODA Python module
- FastJet
- a Fortran compiler such as `gfortran` for POLDIS
- a LaTeX toolchain if rebuilding the manuscript

Exact dependency locations are environment-dependent. The bundled Herwig and
ThePEG snapshots retain their upstream `configure --help` output for available
build options.

## Build ThePEG And Herwig

The repository stores source snapshots only. A clean build can be kept under
`build/`:

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

Build Herwig against the local ThePEG installation:

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

The curated POLDIS tree includes wrappers for the supported validation drivers:

```bash
cd poldis
./compile_dijet
./compile_dijet_rivetplots
./compile_singlejet
```

The wrappers use `gfortran` by default and honor `FC` if another Fortran
compiler is preferred. LHAPDF flags are resolved through `lhapdf-config`.

## Build The Rivet Analyses

The DIS Rivet analysis sources live in `analyses/rivet/dis/`:

```bash
mkdir -p build/rivet
rivet-buildplugin build/rivet/RivetMC_DIS_BREIT.so analyses/rivet/dis/MC_DIS_BREIT.cc
rivet-buildplugin build/rivet/RivetMC_DIS_PS.so analyses/rivet/dis/MC_DIS_PS.cc
export RIVET_ANALYSIS_PATH="$PWD/build/rivet:$PWD/analyses/rivet/dis${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}"
```

The `.plot` and `.info` files are tracked; compiled Rivet plugins are not.

## Run A Validation Workflow

Use `workflow/dispol/` as the workflow base directory. For example, a correlated
helicity-weight run can be launched with:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py full \
  --base-dir workflow/dispol \
  -t <tag> \
  --setup ALL \
  --rivetweights \
  --jobs <ncores> \
  --shards <nshards>
```

A fixed-order comparison run can be launched with:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py full \
  --base-dir workflow/dispol \
  -t <tag> \
  --setup ALL \
  --rivetfo \
  --jobs <ncores> \
  --shards <nshards>
```

Useful follow-up stages for an existing campaign are:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py postprocess --base-dir workflow/dispol -t <tag> --setup ALL
python3 workflow/dispol/scripts/run_validation_campaign.py analyze-herwig --base-dir workflow/dispol -t <tag> --setup ALL
python3 workflow/dispol/scripts/run_validation_campaign.py poldis-top --base-dir workflow/dispol --order NLO
python3 workflow/dispol/scripts/run_validation_campaign.py rivetplot --base-dir workflow/dispol -t <tag> --setup ALL
```

## Paper Reproduction Campaigns

The production campaigns used for the paper are large Tiresias-scale runs. The
commands below give generic versions of the launch commands, with tags and
event counts exposed as shell variables so they can be adjusted for a smoke
test or a fresh production run.

Common Tiresias-style paths:

```bash
BASE=${BASE:-/home/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL}
PYTHON=${PYTHON:-python3.10}
YODAMERGE=${YODAMERGE:-/home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/bin/yodamerge}
```

Neutral-current integrated cross sections, corresponding to the NC validation
tables:

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
validation plots:

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

Shower-spin comparison, used for the no-hadronization and hadronization spin
comparison plots:

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

Generated outputs are expected under:

- `workflow/dispol/campaigns/<tag>/`
- transient `.run`, `.log`, `.out`, `.yoda`, `.csv`, `.html`, and `.tex`
  products under `workflow/dispol/`

Additional workflow details are documented in
`workflow/dispol/docs/validation-workflow.md` and
`workflow/dispol/docs/correlated-helicity-rivetweights.md`.

## Tests

The workflow test suite can be run with:

```bash
python3 -m unittest discover workflow/dispol/scripts/tests
```

For source-level release checks, also verify that the active DIS source cards and
scripts do not regenerate retired validation-only settings, and rebuild the
touched Herwig modules in the target environment.

## Paper Source

The manuscript source lives in:

- `docs/paper/HerwigPolCodex.tex`
- `docs/paper/biblio.bib`

Rendered PDFs, auxiliary LaTeX products, and generated figures are intentionally
excluded from git. The associated pre-print is under preparation; citation
metadata will be added once it is available.

## Acknowledgments

We would like to thank Daniel de Florian for useful discussions and for
providing the POLDIS code used for validation. AP acknowledges support from the
US Department of Energy, Office of Science, Office of Nuclear Physics under
Award Number DE-SC0025728.

## License

The modified Herwig 7 and ThePEG 2.3.0 source snapshots bundled in this
repository are distributed under the terms of the GNU General Public License,
version 3 (GPLv3). A top-level copy of the GPLv3 text is provided in `COPYING`,
and component-level copies are retained in:

- `src/herwig/COPYING`
- `src/thepeg/COPYING`

Additional component-specific notices and exceptions remain in the relevant
subdirectories where applicable.
