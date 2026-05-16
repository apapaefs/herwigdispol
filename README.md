# Polarized Deep-Inelastic Scattering with Spin Correlations in Herwig 7

This repository is the research software and reproducibility companion for the
HerwigPol polarized deep-inelastic scattering implementation developed for
Herwig 7. It brings together the modified Herwig and ThePEG source snapshots,
the curated POLDIS fixed-order reference code, the custom Rivet analyses, the
DIS validation workflow, and the paper source in a single formal repository
layout.

The repository is intended to preserve the source-level ingredients needed to
rebuild and re-run the validated DIS studies. It therefore tracks code, input
cards, workflow drivers, and technical notes, while intentionally excluding
generated artifacts such as build products, campaign outputs, merged YODA
files, plots, and rendered paper outputs.

## Overview

The software collected here supports the implementation and validation of
polarized deep-inelastic scattering at next-to-leading order with spin
correlations in Herwig 7. In addition to the event-generator source changes,
the repository preserves the fixed-order comparison inputs, workflow
orchestration scripts, and analysis components required for end-to-end studies
of integrated and differential DIS observables.

The repository should therefore be understood as a curated scientific software
snapshot rather than as a general-purpose Herwig distribution. Its purpose is
to provide a reproducible basis for the specific DIS programme documented in
the accompanying paper and technical notes.

## Scope

This repository contains the source material needed to reproduce the DIS
studies based on:

- the modified ThePEG 2.3.0 snapshot in `src/thepeg/`
- the modified Herwig 7.3.0 snapshot in `src/herwig/`
- the curated POLDIS source tree and DIS-specific drivers in `poldis/`
- the DIS Rivet analyses `MC_DIS_BREIT` and `MC_DIS_PS`
- the DIS workflow cards, scripts, and notes in `workflow/dispol/`
- the paper source in `docs/paper/`

It does not attempt to archive unrelated side studies from the larger working
area, such as nearby HERMES, RHIC, or WJET materials.

## Provenance

The original clean source starting points for this line of work were:

- `~/Projects/HerwigPol/HwPol`
- `~/Projects/HerwigPol/ThePEGPol`

The source trees committed here are not pristine upstream checkouts. They are
curated snapshots copied from the active working trees used for the DIS
studies:

- `HerwigSource/Herwig-7.3.0` -> `src/herwig/`
- `HerwigSource/ThePEG-2.3.0` -> `src/thepeg/`

Generated libraries, build directories, tarballs, campaign outputs, and other
derived artifacts from the working area have been intentionally omitted.

## Repository Layout

- `src/herwig/`: modified Herwig 7.3.0 source snapshot
- `src/thepeg/`: modified ThePEG 2.3.0 source snapshot
- `poldis/`: curated POLDIS source tree and DIS-specific drivers
- `analyses/rivet/dis/`: DIS Rivet analyses, metadata, and helper wrappers
- `workflow/dispol/cards/`: checked-in DIS Herwig cards and templates
- `workflow/dispol/scripts/`: DIS workflow drivers and post-processing tools
- `workflow/dispol/docs/`: workflow notes and technical investigations
- `docs/paper/`: paper source, including `main.tex` and `biblio.bib`

## External Prerequisites

Typical local prerequisites include:

- a C++ toolchain and `make`
- Autotools-compatible build tooling for the bundled Herwig and ThePEG snapshots
- Python 3 with the `yoda` Python package available to the chosen interpreter
- LHAPDF with `lhapdf-config`
- Rivet and YODA command-line tools
- FastJet
- a Fortran compiler such as `gfortran` for POLDIS
- a LaTeX toolchain for building the paper draft

## Building ThePEG and Herwig

The repository stores source snapshots only. A clean approach is to build in a
separate `build/` area inside the clone.

Build ThePEG:

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

Notes:

- the exact external dependency locations depend on the local environment
- `src/thepeg/configure --help` and `src/herwig/configure --help` list the
  accepted options for each installation
- build outputs under `build/` are ignored by git

## Building POLDIS

The curated POLDIS directory contains wrapper scripts for the supported driver
variants:

```bash
cd poldis
./compile_dijet
./compile_dijet_rivetplots
./compile_singlejet
```

Useful details:

- the wrappers use `gfortran` by default and honor `FC` if another compiler is
  preferred
- LHAPDF flags are resolved through `lhapdf-config`
- the canonical DIS comparison driver is `poldis/user_dijet_rivetplots.f`
- the generated `poldis/poldis.x` executable is intentionally not tracked

## Building the Rivet Analyses

The DIS analysis sources live in `analyses/rivet/dis/`. A clean approach is to
build the generated plugins in a disposable build directory:

```bash
mkdir -p build/rivet
rivet-buildplugin build/rivet/RivetMC_DIS_BREIT.so analyses/rivet/dis/MC_DIS_BREIT.cc
rivet-buildplugin build/rivet/RivetMC_DIS_PS.so analyses/rivet/dis/MC_DIS_PS.cc
export RIVET_ANALYSIS_PATH="$PWD/build/rivet:$PWD/analyses/rivet/dis${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}"
```

The `.plot` and `.info` metadata files are tracked in the source tree, whereas
the compiled plugins are not.

## Workflow Drivers

The public workflow entry points live in `workflow/dispol/scripts/`. The main
drivers are:

- `run_validation_campaign.py`: campaign orchestration, post-processing, POLDIS
  conversion, and Rivet plotting
- `herwig_fixed_order_nlo.py`: standalone Python fixed-order neutral-current
  DIS NLO validator using Herwig's analytic components and signed HepMC3/Rivet
  callback events
- `analyze_DIS_polarized.py`: primary DIS YODA analysis step
- `analyze-DIS-polarized.py`: retained legacy name for compatibility
- `poldis_top_to_yoda.py`: conversion of POLDIS `.top` outputs into Rivet-style
  YODA references
- `poldis_top_to_yoda.sh`: lightweight wrapper around the converter
- `powheg_raw_momenta_to_yoda.py`: conversion of raw POWHEG momenta diagnostics
- `extract_dis_out_results.py`: extraction of text summaries from Herwig outputs
- `extract_nlo_term_diagnostics.py`: derivation of NLO term diagnostics
- `extract_powheg_real_spin_diagnostics.py`: summary of POWHEG real-emission
  spin diagnostics
- `recover_campaign_manifest.py`: reconstruction of a missing campaign manifest
  from shard artifacts and `launcher-logs/`
- `compare_nlo_gamma.py`: comparison of NLO gamma-channel outputs
- `compare_yoda_areas.py`: comparison of YODA integral behavior
- `analyze_raw_powheg_summary_csv.py`: inspection of raw POWHEG summary CSVs
- `parse_nlo_cum.py`: parsing of cumulative NLO summaries
- `rivet_mkhtml_safe.py`: safer HTML plot wrapper for Rivet
- `rivet_scale_plot_postprocess.py`: post-processing of scale-variation plot
  outputs

## Standalone Fixed-Order NLO Validator

The standalone validator in
`workflow/dispol/scripts/herwig_fixed_order_nlo.py` is an independent Python
implementation of the neutral-current DIS fixed-order NLO ingredients used in
the Herwig implementation. It is intended for component-level validation
against fixed-order references, especially POLDIS, without invoking the Herwig
POWHEG event-carrier path.

The implementation ports the Herwig choices for:

- beam energies, DIS cuts, helicity conventions, and neutral-current
  `GAMMA`, `Z`, and `ALL` setups
- Matchbox-style NLO running `alpha_s(MZ)=0.118`
- neutral-current electroweak coefficients and polarized `A_pol` factors
- LHAPDF unpolarized and polarized-difference PDF ratios
- the Herwig `NLOWeightRaw` virtual, collinear, QCDC, and BGF terms
- local QCDC/BGF real-emission kernels and Born-projected counterevents

It writes signed HepMC3 ASCII event streams and can run the same
`MC_DIS_BREIT:JETINPUT=MEPARTONS` Rivet analysis used for the Herwig and POLDIS
comparisons. The correlated polarized mode writes named helicity weights
(`HERWIG_DIS_PP`, `HERWIG_DIS_PM`, `HERWIG_DIS_MP`, `HERWIG_DIS_MM`,
`HERWIG_DIS_SIGMA`, and `HERWIG_DIS_DELTA_LL`) and derives the corresponding
`D*` and `ALL*` YODA objects after Rivet. Scale variations follow the Herwig
campaign convention with `ScaleFactorDown = 0.5`, nominal `1.0`, and
`ScaleFactorUp = 2.0`, using `mu^2 = Q^2 * factor^2`.

The validator does not use code from POLDIS or POWHEG-BOX, does not call back
into Herwig C++, and does not include showering, hadronization, decays, POWHEG
Sudakov generation, or POWHEG-selected real-emission carrier events. Its role is
to validate the fixed-order NLO components separately from the Herwig
POWHEG/no-shower event-carrier machinery.

A small local run, assuming Rivet/YODA/LHAPDF are on the environment path, is:

```bash
python3 workflow/dispol/scripts/herwig_fixed_order_nlo.py standalone-campaign \
  -t standalone_fo_q2gt100_smoke \
  --setup ALL \
  --window validation \
  --events 2000 \
  --shards 4 \
  --jobs 2 \
  --poldis-reference workflow/dispol/campaigns/rivetweights_noshowerMac03/analysis/_rivetplot_inputs/reference.scale.reference_all.sanitized.yoda.gz
```

For a higher-statistics validation run:

```bash
python3 workflow/dispol/scripts/herwig_fixed_order_nlo.py standalone-campaign \
  -t standalone_fo_q2gt100_highstat \
  --setup ALL \
  --window validation \
  --events 500000 \
  --shards 50 \
  --jobs 10 \
  --shard-monitor-interval 20 \
  --poldis-reference workflow/dispol/campaigns/rivetweights_noshowerMac03/analysis/_rivetplot_inputs/reference.scale.reference_all.sanitized.yoda.gz
```

During parallel runs the script prints shard progress and refreshes
`workflow/dispol/campaigns/<tag>/work/shard_status.tsv`, which records the
variation, shard id, event count, seed, state, return code, elapsed time, raw
YODA path, and stdout/stderr logs for each shard. Generated outputs are written
under `workflow/dispol/campaigns/<tag>/` and remain excluded from git.

## End-to-End DIS Workflow

The checked-in input cards live in `workflow/dispol/cards/`. The repository
workflow convention is to use `workflow/dispol/` as the `--base-dir`, which
keeps generated outputs under one subtree while leaving the curated cards in
`workflow/dispol/cards/`.

Typical steps from the repository root are:

1. Launch a campaign:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign \
  --base-dir workflow/dispol \
  -t testing13 \
  --rivet \
  --jobs 32 \
  --shards 0
```

2. Rebuild merged outputs and summaries for an existing tag:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py postprocess \
  --base-dir workflow/dispol \
  -t testing13 \
  --rivet
```

3. Build analyzed DIS YODAs:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py analyze-herwig \
  --base-dir workflow/dispol \
  -t testing13
```

4. Convert POLDIS `.top` files into a reference YODA:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py poldis-top \
  --base-dir workflow/dispol \
  --order NLO
```

5. Produce Rivet comparison plots:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py rivetplot \
  --base-dir workflow/dispol \
  -t testing13 \
  --setup ALL
```

6. Run the full chain:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py full \
  --base-dir workflow/dispol \
  -t testing13 \
  --rivet \
  --jobs 32 \
  --shards 0
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

Operational notes:

- `SPINCOMP` and `SPINHAD` shower cards are configured with low-verbosity
  `EventGenerator` settings so large parton-shower campaigns do not flood the
  `.log` files
- if a campaign crashes before
  `workflow/dispol/campaigns/<tag>/manifest.json` is written, it can be
  recovered with `workflow/dispol/scripts/recover_campaign_manifest.py`, then
  resumed with
  `workflow/dispol/scripts/run_validation_campaign.py full --rerun-failed-random-seed`

Additional workflow details are documented in
`workflow/dispol/docs/validation-workflow.md`.

## Paper Source

The paper draft lives in:

- `docs/paper/main.tex`
- `docs/paper/biblio.bib`

Only the TeX and bibliography sources are versioned. Rendered figure PDFs and
other LaTeX build products are intentionally excluded, so figures must be
regenerated locally or supplied separately before producing a final PDF.

The paper associated with this repository is titled:

- `Polarized Deep-Inelastic Scattering with Spin Correlations in Herwig 7`

## License

The modified Herwig 7 and ThePEG 2.3.0 source snapshots bundled in this
repository are distributed under the terms of the GNU General Public License,
version 3 (GPLv3). A top-level copy of the GPLv3 text is provided in
`COPYING`, and the corresponding component-level license files are retained in:

- `src/herwig/COPYING`
- `src/thepeg/COPYING`

Additional component-specific notices and exceptions remain in the relevant
subdirectories, where applicable.

## Acknowledgments

We would like to thank Daniel de Florian for useful discussions and for
providing the POLDIS code used for validation. AP acknowledges support from
the U.S. Department of Energy, Office of Science, Office of Nuclear Physics,
under Award Number DE-SC0025728.

## What Is Intentionally Excluded

The repository is source-first. The main excluded categories are:

- compiled libraries and in-tree build products
- campaign outputs and shard artifacts
- `.run`, `.log`, `.out`, `.yoda`, `.csv`, `.html`, and `.top` generated files
- compiled Rivet plugins
- rendered paper artifacts and figure PDFs
- stale duplicate files when a newer canonical source exists elsewhere in the
  curated layout

## Migration from the Older Working Layout

- `HwPolNotesNew/DISPOL` -> `workflow/dispol`
- DIS pieces from `HwPolNotesNew/analyses` -> `analyses/rivet/dis`
- `POLDIS/POLDIS-public` plus retained DIS helpers -> `poldis`
- `main.tex` and `biblio.bib` -> `docs/paper`

The canonical public entry points in this repository are the repo-relative
paths under `workflow/dispol/scripts/`, especially:

- `workflow/dispol/scripts/run_validation_campaign.py`
- `workflow/dispol/scripts/analyze_DIS_polarized.py`
- `workflow/dispol/scripts/poldis_top_to_yoda.py`
