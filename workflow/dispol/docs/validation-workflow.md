# DISPOL Validation Workflow

This note records the generic validation workflow driven by
`DISPOL/scripts/run_validation_campaign.py`.  Run the commands from the
repository root and pass `--base-dir DISPOL` unless you intentionally maintain a
separate run area.

## POWHEG Emission Comparison Modes

The DIS POWHEG implementation exposes the following switches for validation
against fixed-order expectations:

- `UseFixedOrderAlphaSInPOWHEGEmission`
- `UseQ2ScaleInPOWHEGEmission`
- `POWHEGEmissionComparisonMode`
- `POWHEGEmissionComparisonMaxAttempts`

Native Herwig POWHEG generation is the default:

```text
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission No
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission No
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode Default
```

The fixed-order comparison cards use the same `Q^2` scale and PDF/alpha_s
treatment as the POLDIS-style reference comparison:

```text
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseFixedOrderAlphaSInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:UseQ2ScaleInPOWHEGEmission Yes
set /Herwig/MatrixElements/PowhegMEDISNCPol:POWHEGEmissionComparisonMode Default
```

The non-default `RealOnly` mode is a diagnostic path that retries the POWHEG
emission and vetoes failed trials rather than falling back to the underlying
Born configuration.

## Campaigns

Create a campaign with the active DIS cards:

```bash
python3 DISPOL/scripts/run_validation_campaign.py campaign \
  --base-dir DISPOL \
  -t <tag> \
  --setup ALL \
  --rivetfo \
  --jobs <ncores> \
  --shards <nshards>
```

Add `--force-prepare` after card edits to rebuild `.run` files.  `--shards 0`
lets the driver choose a shard count from `--jobs`.

Common setup names include `GAMMA`, `Z`, `ALL`, `CC`, `SPINCOMP`, `SPINHAD`,
and `SPINVAL`; the exact set depends on whether the campaign is LO, POWHEG,
Rivet, `RIVETFO`, parton-shower, or `RIVETWEIGHTS`.

## Postprocessing And Analysis

For an existing campaign, rebuild merged YODA files and summaries:

```bash
python3 DISPOL/scripts/run_validation_campaign.py postprocess \
  --base-dir DISPOL \
  -t <tag> \
  --setup ALL \
  --rivetfo
```

Analyze Herwig YODAs into the polarized DIS observable basis:

```bash
python3 DISPOL/scripts/run_validation_campaign.py analyze-herwig \
  --base-dir DISPOL \
  -t <tag> \
  --setup ALL \
  --rivetfo
```

Build Rivet comparison plots:

```bash
python3 DISPOL/scripts/run_validation_campaign.py rivetplot \
  --base-dir DISPOL \
  -t <tag> \
  --setup ALL \
  --rivetfo
```

Campaign products are written under `DISPOL/campaigns/<tag>/`, including raw
YODAs, merged physical `NLO = POSNLO - NEGNLO` YODAs, analyzed Herwig YODAs,
reference YODAs, plot directories, and `manifest.json`.

## POLDIS References

The `poldis-top` subcommand converts POLDIS Topdrawer output into the YODA
reference format expected by the plotting stage:

```bash
python3 DISPOL/scripts/run_validation_campaign.py poldis-top \
  --base-dir DISPOL \
  --order NLO
```

Use `--unpol`, `--pol`, and `--out` to override the default input and output
paths.  Accepted reference orders are `LO`, `NLO`, and `NNLO`.

## Correlated Helicity Weights

The `RIVETWEIGHTS` path evaluates all four helicity configurations on the same
phase-space point and writes the combinations as named weights.  It is useful
for cancellation-sensitive observables:

```bash
python3 DISPOL/scripts/run_validation_campaign.py full \
  --base-dir DISPOL \
  -t <tag> \
  --setup ALL \
  --rivetweights \
  --jobs <ncores> \
  --shards <nshards>
```

Use `--rivetweights-shower` for the showered parton-level sibling of the same
workflow. It keeps hadronization off, runs `MC_DIS_BREIT` on the full showered
partonic final state, and leaves ordinary `--rivetweights` unchanged.

See `DISPOL/docs/correlated-helicity-rivetweights.md` for the named-weight
semantics used by `MC_DIS_BREIT` and `MC_DIS_PS`.

## Help

Each driver stage has its own help output:

```bash
python3 DISPOL/scripts/run_validation_campaign.py --help
python3 DISPOL/scripts/run_validation_campaign.py campaign --help
python3 DISPOL/scripts/run_validation_campaign.py postprocess --help
python3 DISPOL/scripts/run_validation_campaign.py analyze-herwig --help
python3 DISPOL/scripts/run_validation_campaign.py poldis-top --help
python3 DISPOL/scripts/run_validation_campaign.py rivetplot --help
python3 DISPOL/scripts/run_validation_campaign.py full --help
```
