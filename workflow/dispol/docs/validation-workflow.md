# DISPOL Validation Workflow

This note records the repo-native commands for the curated DIS workflow in this
repository. All examples assume the repository root is the current working
directory and use:

- `workflow/dispol/cards/` for checked-in Herwig input cards and templates
- `workflow/dispol/scripts/` for the Python and shell workflow drivers
- `analyses/rivet/dis/` for the custom Rivet analyses
- `poldis/` for the fixed-order POLDIS source and build wrappers

The main entrypoint is:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py --help
```

## POWHEG comparison modes

The DIS POWHEG comparison machinery uses:

- `UseFixedOrderAlphaSInPOWHEGEmission`
- `UseQ2ScaleInPOWHEGEmission`
- `POWHEGEmissionComparisonMode`
- `POWHEGEmissionComparisonMaxAttempts`

Common settings:

```text
# native Herwig POWHEG
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission No
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission No
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode Default

# direct POLDIS-comparison mode
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode Default

# diagnostic no-Born-fallback mode
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode RealOnly
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMaxAttempts 100
```

## 1. Launch a validation campaign

Standard campaign:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign \
  --base-dir workflow/dispol \
  -t testing13 \
  --jobs 32 \
  --shards 0 \
  --seed-base 100000 \
  --lo-events 1000000 \
  --posnlo-events 1000000 \
  --negnlo-events 100000
```

Rebuild `.run` files after card edits:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign \
  --base-dir workflow/dispol \
  -t testing13 \
  --force-prepare \
  --jobs 32 \
  --shards 0
```

Restrict the run to a specific setup:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign \
  --base-dir workflow/dispol \
  -t testing13 \
  --setup ALL \
  --jobs 32 \
  --shards 0
```

Rivet-enabled POWHEG campaign:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign \
  --base-dir workflow/dispol \
  -t rivet01 \
  --rivet \
  --jobs 32 \
  --shards 0
```

Direct fixed-order diagnostic family:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign \
  --base-dir workflow/dispol \
  -t rivetfofixed01 \
  --rivetfofixed \
  --jobs 32 \
  --shards 0 \
  --posnlo-events 1000000
```

Notes:

- `--shards 0` auto-chooses the shard count from `--jobs`
- merged run-level YODAs and physical `NLO = POSNLO - NEGNLO` outputs are built
  automatically during campaign post-processing
- `--rivetfo` selects the scale-matched POWHEG fixed-order comparison cards
- `--rivetfofixed` runs the direct fixed-order event-record diagnostic family
  and only launches the meaningful `POSNLO` jobs

## 2. Rebuild merged YODAs and summaries

Use this when shard outputs already exist and you want to rebuild merged YODAs,
analysis products, and summaries:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py postprocess \
  --base-dir workflow/dispol \
  -t rivet01 \
  --rivet \
  --yoda-merge-tool auto
```

Preview only:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py postprocess \
  --base-dir workflow/dispol \
  -t rivet01 \
  --rivet \
  --dry-run
```

## 3. Build analyzed DIS YODAs

Default setup:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py analyze-herwig \
  --base-dir workflow/dispol \
  -t testing13
```

Explicit setup selection:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py analyze-herwig \
  --base-dir workflow/dispol \
  -t testing13 \
  --setup GAMMA \
  --setup Z \
  --setup ALL
```

Outputs are written under:

- `workflow/dispol/campaigns/<tag>/analysis/`

## 4. Convert POLDIS `.top` references into YODA

Default local files:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py poldis-top \
  --base-dir workflow/dispol \
  --order NLO
```

Explicit files and output:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py poldis-top \
  --base-dir workflow/dispol \
  --unpol workflow/dispol/dijets_unpol.top \
  --pol workflow/dispol/dijets_pol.top \
  --order NLO \
  --out workflow/dispol/POLDIS_MC_DIS_BREIT_ref.yoda.gz
```

Accepted reference orders are `LO`, `NLO`, and `NNLO`.

## 5. Produce Rivet comparison plots

Once analyzed Herwig YODAs and a POLDIS reference are available:

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py rivetplot \
  --base-dir workflow/dispol \
  -t testing13 \
  --setup ALL
```

The helper wrapper in `analyses/rivet/dis/rivet-mkhtml-hw-vs-poldis-example.sh`
can also be used directly:

```bash
analyses/rivet/dis/rivet-mkhtml-hw-vs-poldis-example.sh \
  workflow/dispol/campaigns/testing13/analysis/Herwig_MC_DIS_BREIT_ALL_NLO_polarized.yoda.gz \
  workflow/dispol/campaigns/testing13/refs/POLDIS_MC_DIS_BREIT_ref.yoda.gz \
  workflow/dispol/campaigns/testing13/plots_mc_vs_poldis
```

## 6. Run the full chain

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py full \
  --base-dir workflow/dispol \
  -t testing13 \
  --rivet \
  --jobs 32 \
  --shards 0
```

## Output locations

Generated products are intentionally kept out of version control. The main
locations are:

- `workflow/dispol/campaigns/<tag>/`
- `workflow/dispol/campaigns/<tag>/yoda/`
- `workflow/dispol/campaigns/<tag>/analysis/`
- `workflow/dispol/campaigns/<tag>/refs/`
- `workflow/dispol/campaigns/<tag>/plots_mc_vs_poldis_<setup>_<order>/`
- transient `.run`, `.log`, `.out`, `.yoda`, `.csv`, `.html`, and `.tex` files
  under `workflow/dispol/`

POLDIS `.top` files are generated inputs to the conversion step and are not
tracked in this curated repository.

## Help

```bash
python3 workflow/dispol/scripts/run_validation_campaign.py campaign --help
python3 workflow/dispol/scripts/run_validation_campaign.py postprocess --help
python3 workflow/dispol/scripts/run_validation_campaign.py analyze-herwig --help
python3 workflow/dispol/scripts/run_validation_campaign.py poldis-top --help
python3 workflow/dispol/scripts/run_validation_campaign.py rivetplot --help
python3 workflow/dispol/scripts/run_validation_campaign.py full --help
```
