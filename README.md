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
