# herwigdispol

`herwigdispol` is a curated source repository for the polarized deep-inelastic
scattering program developed around modified Herwig and ThePEG source
snapshots, the POLDIS fixed-order reference code, custom Rivet analyses, the
DIS validation workflow, and the current paper draft.

The repository is intended to support reproducibility from source. It preserves
the code, input cards, workflow drivers, and technical notes needed to rebuild
and re-run the DIS studies, while intentionally excluding generated artifacts
such as build products, campaign outputs, plots, merged YODA files, and
rendered paper outputs.

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
- `docs/paper/`: paper source, including `HerwigPolCodex.tex` and `biblio.bib`

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

- `docs/paper/HerwigPolCodex.tex`
- `docs/paper/biblio.bib`

Only the TeX and bibliography sources are versioned. Rendered figure PDFs and
other LaTeX build products are intentionally excluded, so figures must be
regenerated locally or supplied separately before producing a final PDF.

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
- `HerwigPolCodex.tex` and `biblio.bib` -> `docs/paper`

The canonical public entry points in this repository are the repo-relative
paths under `workflow/dispol/scripts/`, especially:

- `workflow/dispol/scripts/run_validation_campaign.py`
- `workflow/dispol/scripts/analyze_DIS_polarized.py`
- `workflow/dispol/scripts/poldis_top_to_yoda.py`
