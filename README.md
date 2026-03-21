# herwigdispol

This repository is a curated, DIS-focused basis checkout for the polarized DIS
work built around modified Herwig and ThePEG source snapshots, the POLDIS
fixed-order reference code, the custom Rivet analyses, the DIS workflow
drivers, and the current paper draft.

The aim is reproducibility from source, not archival storage of generated
outputs. Build products, campaign outputs, YODAs, plots, and rendered paper
artifacts are intentionally excluded from version control.

## Scope

This repository keeps the source needed to reproduce the DIS studies around:

- the modified ThePEG 2.3.0 snapshot in `src/thepeg/`
- the modified Herwig 7.3.0 snapshot in `src/herwig/`
- the curated POLDIS source and DIS driver in `poldis/`
- the DIS Rivet analyses `MC_DIS_BREIT` and `MC_DIS_PS`
- the DIS workflow cards, scripts, and notes in `workflow/dispol/`
- the paper source in `docs/paper/`

It intentionally does not try to preserve nearby HERMES, RHIC, WJET, or other
side-study material from the larger working area.

## Provenance

The original clean source starting points for this work were:

- `~/Projects/HerwigPol/HwPol`
- `~/Projects/HerwigPol/ThePEGPol`

The source trees committed here are not those pristine starting points. They
are curated snapshots copied from the active working trees that were used for
the DIS studies:

- `HerwigSource/Herwig-7.3.0` -> `src/herwig/`
- `HerwigSource/ThePEG-2.3.0` -> `src/thepeg/`

Generated libraries, build directories, tarballs, and campaign artifacts from
the working area were intentionally left out.

## Repository Layout

- `src/herwig/`: modified Herwig 7.3.0 source snapshot
- `src/thepeg/`: modified ThePEG 2.3.0 source snapshot
- `poldis/`: curated POLDIS source tree and DIS-specific drivers
- `analyses/rivet/dis/`: DIS Rivet analyses, metadata, and small helper wrappers
- `workflow/dispol/cards/`: checked-in DIS Herwig cards and templates
- `workflow/dispol/scripts/`: DIS workflow drivers and post-processing helpers
- `workflow/dispol/docs/`: workflow and technical notes
- `docs/paper/`: `HerwigPolCodex.tex` and `biblio.bib`

## External Prerequisites

You will typically need:

- a C++ toolchain and `make`
- Autotools-compatible build tooling for the bundled Herwig/ThePEG snapshots
- Python 3 with the `yoda` Python package available to the chosen interpreter
- LHAPDF with `lhapdf-config`
- Rivet and YODA command-line tools
- FastJet
- a Fortran compiler such as `gfortran` for POLDIS
- a LaTeX toolchain if you want to build the paper draft

## Building ThePEG And Herwig

The repository only stores source snapshots. A clean approach is to build in a
separate `build/` area inside the clone:

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

- the exact external dependency locations depend on your local setup
- `src/thepeg/configure --help` and `src/herwig/configure --help` list the
  accepted options for your environment
- build outputs under `build/` are ignored by git

## Building POLDIS

The curated POLDIS directory contains portable wrapper scripts that compile the
three supported driver variants in place:

```bash
cd poldis
./compile_dijet
./compile_dijet_rivetplots
./compile_singlejet
```

Useful details:

- the wrappers use `gfortran` by default and honor `FC` if you want another
  compiler
- they resolve LHAPDF flags via `lhapdf-config`
- the canonical DIS comparison driver is `poldis/user_dijet_rivetplots.f`
- `poldis/poldis.x` is generated locally and ignored by git

## Building The Rivet Analyses

The DIS analysis sources live in `analyses/rivet/dis/`. A clean way to build
them is to place the generated plugins in a throwaway build directory:

```bash
mkdir -p build/rivet
rivet-buildplugin build/rivet/RivetMC_DIS_BREIT.so analyses/rivet/dis/MC_DIS_BREIT.cc
rivet-buildplugin build/rivet/RivetMC_DIS_PS.so analyses/rivet/dis/MC_DIS_PS.cc
export RIVET_ANALYSIS_PATH="$PWD/build/rivet:$PWD/analyses/rivet/dis${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}"
```

The `.plot` and `.info` files are kept in the source tree, while the compiled
plugins are ignored by git.

## Workflow Drivers

The public workflow entrypoints live in `workflow/dispol/scripts/`. The main
ones are:

- `run_validation_campaign.py`: campaign orchestration, post-processing, POLDIS
  conversion, and Rivet plotting
- `analyze_DIS_polarized.py`: main DIS YODA analysis step
- `analyze-DIS-polarized.py`: retained legacy name for compatibility
- `poldis_top_to_yoda.py`: convert POLDIS `.top` outputs into Rivet-style YODA
  references
- `poldis_top_to_yoda.sh`: lightweight wrapper around the converter
- `powheg_raw_momenta_to_yoda.py`: convert raw POWHEG momenta diagnostics
- `extract_dis_out_results.py`: harvest text summaries from Herwig outputs
- `extract_nlo_term_diagnostics.py`: derive NLO term diagnostics
- `extract_powheg_real_spin_diagnostics.py`: summarize POWHEG real-emission spin
  diagnostics
- `compare_nlo_gamma.py`: compare NLO gamma-channel outputs
- `compare_yoda_areas.py`: compare YODA integral behavior
- `analyze_raw_powheg_summary_csv.py`: inspect raw POWHEG summary CSVs
- `parse_nlo_cum.py`: parse cumulative NLO summaries
- `rivet_mkhtml_safe.py`: safer HTML plot wrapper for Rivet
- `rivet_scale_plot_postprocess.py`: post-process scale-variation plot outputs

## End-To-End DIS Workflow

The checked-in input cards live in `workflow/dispol/cards/`. The repo-native
workflow convention is to use `workflow/dispol/` as the `--base-dir`, which
keeps generated outputs under one subtree while leaving the cards in
`workflow/dispol/cards/`.

Typical steps from the repository root:

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

6. Or run the full chain:

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

The more detailed workflow note is in
`workflow/dispol/docs/validation-workflow.md`.

## Paper Source

The paper draft lives in:

- `docs/paper/HerwigPolCodex.tex`
- `docs/paper/biblio.bib`

Only the TeX and bibliography sources are committed. Rendered figure PDFs and
other LaTeX build products are intentionally excluded, so figures must be
regenerated locally or supplied separately before producing a final PDF.

## What Is Intentionally Excluded

The repository is source-first. The main excluded categories are:

- compiled libraries and in-tree build products
- campaign outputs and shard artifacts
- `.run`, `.log`, `.out`, `.yoda`, `.csv`, `.html`, and `.top` generated files
- compiled Rivet plugins
- rendered paper artifacts and figure PDFs
- stale duplicate files when a newer canonical source exists elsewhere in the
  curated layout

## Migration From The Older Working Layout

- `HwPolNotesNew/DISPOL` -> `workflow/dispol`
- DIS pieces from `HwPolNotesNew/analyses` -> `analyses/rivet/dis`
- `POLDIS/POLDIS-public` plus retained DIS helpers -> `poldis`
- `HerwigPolCodex.tex` and `biblio.bib` -> `docs/paper`

The canonical public entrypoints in this repository are the repo-relative paths
under `workflow/dispol/scripts/`, especially:

- `workflow/dispol/scripts/run_validation_campaign.py`
- `workflow/dispol/scripts/analyze_DIS_polarized.py`
- `workflow/dispol/scripts/poldis_top_to_yoda.py`
