# DISPOL Validation Workflow

This note records the main commands for the integrated DISPOL workflow driven by:

- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py`

## POWHEG emission comparison modes

The DIS POWHEG emission now has:

- two scale/source switches:
  - `UseFixedOrderAlphaSInPOWHEGEmission`
  - `UseQ2ScaleInPOWHEGEmission`
- one comparison mode:
  - `POWHEGEmissionComparisonMode`
- one retry parameter for the non-default comparison mode:
  - `POWHEGEmissionComparisonMaxAttempts`

Native Herwig POWHEG remains the default:

- `POWHEGEmissionComparisonMode Default`

`UseQ2ScaleInPOWHEGEmission` controls both:

- the POWHEG/MEC emission `alpha_s` scale
- the POWHEG/MEC emission-PDF numerator scale

`UseFixedOrderAlphaSInPOWHEGEmission` means:

- prefer the LHAPDF-set `alpha_s` from the beam PDF when available
- otherwise fall back to Herwig's fixed-order/model `alpha_s`

`POWHEGEmissionComparisonMode` means:

- `Default`: exact native Herwig POWHEG behavior
- `RealOnly`: native POWHEG generation, but retry and veto instead of falling back to Born

`POWHEGEmissionComparisonMaxAttempts` defaults to `100` and is only used when `POWHEGEmissionComparisonMode = RealOnly`.

Current/native POWHEG:

```text
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission No
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission No
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode Default
```

Direct POLDIS-comparison mode:

```text
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode Default
```

No-Born-fallback native POWHEG diagnostic:

```text
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode RealOnly
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMaxAttempts 100
```

## 1. Run a Herwig validation campaign

Standard campaign:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py campaign \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13 \
  --jobs 96 \
  --shards 0 \
  --seed-base 100000 \
  --lo-events 1000000000 \
  --posnlo-events 1000000000 \
  --negnlo-events 10000000
```

If you changed any `.in` cards and need the `.run` files rebuilt, add `--force-prepare`:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py campaign \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13 \
  --force-prepare \
  --jobs 96 \
  --shards 0 \
  --seed-base 100000 \
  --lo-events 1000000000 \
  --posnlo-events 1000000000 \
  --negnlo-events 10000000
```

To run only the `ALL` setup in the campaign stage:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py campaign \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13 \
  --setup ALL \
  --jobs 96 \
  --shards 0 \
  --seed-base 100000 \
  --lo-events 1000000000 \
  --posnlo-events 1000000000 \
  --negnlo-events 10000000
```

Rivet-enabled campaign:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py campaign \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t rivet01 \
  --rivet \
  --jobs 96 \
  --shards 0 \
  --seed-base 100000 \
  --lo-events 1000000000 \
  --posnlo-events 1000000000 \
  --negnlo-events 10000000
```

Notes:
- `--shards 0` auto-chooses shard count from `--jobs`
- shard outputs are merged automatically
- physical `NLO = POSNLO - NEGNLO` YODAs are built automatically
- for Rivet campaigns, merged YODAs are rebuilt with the custom `MC_DIS_BREIT`
  analysis path set in `RIVET_ANALYSIS_PATH`
- current `RIVETFO` cards use `MC_DIS_BREIT:JETINPUT=TOP2PARTONS`
- `--rivetfofixed` is available as a parallel fixed-order comparison family
  using `-RIVETFOFIXED.in` cards and the direct
  `MENeutralCurrentDISFixedOrder` matrix element with the cascade disabled
  In the current workflow, this is a POS-only direct-real diagnostic path:
  it runs only the meaningful `POSNLO` jobs and skips `LO`/`NEGNLO`.

Direct fixed-order Rivet campaign:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py campaign \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t rivetfofixed01 \
  --rivetfofixed \
  --jobs 96 \
  --shards 0 \
  --seed-base 100000 \
  --posnlo-events 1000000000
```

## 2. Rebuild merged YODAs and summaries for an existing tag

Use this when the raw shard outputs already exist and you want to:
- rebuild merged run-level YODAs
- rebuild physical `NLO = POSNLO - NEGNLO` YODAs
- rerun the text/CSV/JSON summaries
- rerun the Herwig DIS YODA analysis

Typical Rivet postprocess command:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py postprocess \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t rivet2 \
  --rivet \
  --yoda-merge-tool auto
```

Preview only:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py postprocess \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t rivet2 \
  --rivet \
  --dry-run
```

## 3. Analyze Herwig NLO YODAs into DIS polarized observables

Default setup is `ALL`:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py analyze-herwig \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13
```

Explicit setup:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py analyze-herwig \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13 \
  --setup GAMMA \
  --setup Z \
  --setup ALL
```

This writes analyzed YODAs under:

- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/analysis/`

The analyzer accepts:

- `--rivet` for the standard POWHEG+Rivet cards
- `--rivetfo` for the scale-matched POWHEG fixed-order-comparison cards
- `--rivetfofixed` for the direct fixed-order event-record cards

## 4. Convert POLDIS `.top` files into reference YODA

Default local `.top` files:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py poldis-top \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  --order NLO
```

Accepted POLDIS reference orders:
- `--order LO`
- `--order NLO`
- `--order NNLO`

Explicit files and output:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py poldis-top \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  --unpol /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/dijets_unpol.top \
  --pol /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/dijets_pol.top \
  --order NLO \
  --out /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/POLDIS_MC_DIS_BREIT_ref.yoda.gz
```

Alternative order selections:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py poldis-top \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  --order LO

python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py poldis-top \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  --order NNLO
```

## 5. Run everything in one go

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py full \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13 \
  --rivet \
  --jobs 96 \
  --shards 0 \
  --seed-base 100000 \
  --lo-events 1000000000 \
  --posnlo-events 1000000000 \
  --negnlo-events 10000000 \
  --setup ALL \
  --order NLO
```

Example full Rivet campaign command used for `rivet13`:

```bash
python /home/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py \
  --base-dir /home/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t rivet13 \
  --jobs 200 \
  --shards 0 \
  --seed-base 200000 \
  --lo-events 100000 \
  --posnlo-events 100000 \
  --negnlo-events 10000 \
  --progress-interval 2 \
  --max-listed 50 \
  --rivet \
  --yoda-merge-tool /home/apapaefs/Projects/Herwig/Herwig-pol-full-python3-rivet4/bin/yodamerge
```

This populates:

- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/results.*`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/diagnostic.*`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/yoda/`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/analysis/`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/refs/`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/plots_mc_vs_poldis_<setup>_<order>/`
- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/manifest.json`

## 6. Build Rivet comparison plots directly

Default `ALL` NLO comparison against the chosen POLDIS order:

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py rivetplot \
  --base-dir /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL \
  -t testing13 \
  --setup ALL \
  --order NLO
```

This uses the analyzed Herwig YODA from:

- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/analysis/Herwig_MC_DIS_BREIT_<setup>_NLO_polarized.yoda.gz`

and the POLDIS reference from:

- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/<tag>/refs/POLDIS_MC_DIS_BREIT_ref.yoda.gz`

It runs a command of the form:

```bash
rivet-mkhtml <Herwig analyzed yoda> <POLDIS ref yoda> --reflabel="POLDIS NLO" --ratiolabel "MC/POLDIS" -o <campaign plot dir>
```

## 7. Useful help commands

```bash
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py --help
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py campaign --help
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py postprocess --help
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py analyze-herwig --help
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py poldis-top --help
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py rivetplot --help
python3.10 /Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/run_validation_campaign.py full --help
```
