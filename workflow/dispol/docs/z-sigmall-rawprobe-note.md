# Z sigma_LL Raw-Probe Note

## Current status

The broad `plain38` regression is gone. The `plain39` code state restored the
good `GAMMA`/`ALL` behavior while keeping the NC-odd algebraic repair in
`DISBase.cc`.

The only notable remaining production discrepancy is the pure-`Z`,
spin-dependent observable:

- Herwig `Z` NLO `sigma_LL = 0.0491429 +- 0.0000397 pb`
- POLDIS `Z` NLO `sigma_LL = 0.0501390 +- 0.0000100 pb`
- delta `= -0.0009961 +- 0.0000410 pb`
- ratio `= 0.980133 +- 0.000816`

For context, the corresponding `plain39` results are otherwise good:

- `GAMMA` NLO unpolarized: Herwig `2755.663 pb`, POLDIS `2756.009 pb`
- `GAMMA` NLO polarized: Herwig `43.595 pb`, POLDIS `43.66951 pb`
- `ALL` NLO unpolarized ratio: `0.999822`
- `ALL` NLO polarized ratio: `0.994784`

## Benchmark hygiene

`results_testing15.csv` is now only a historical snapshot from `2026-03-10`.
It predates the `2026-04-07` `plain39` recovery and should not be used as the
current pure-`Z` benchmark.

Use instead:

- `/Users/apapaefs/trs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/plain39/results.txt`
- `/Users/apapaefs/trs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/plain39/results.csv`

This matters because the qualitative pattern changed:

- `testing15`: pure-`Z` `sigma0 - 00` was nearly zero and the Z normalization
  sat high
- `plain39`: pure-`Z` `sigma0 - 00` is clearly nonzero and both `sigma0` and
  `sigma_LL` sit low

## What was tested

Extra Z-only diagnostics were added to distinguish:

1. the `sigma_LL` contribution breakdown by term
2. raw-vs-clamped polarization effects for
   - `Pq`
   - `Pq_m`
   - `Pg_m`
3. shadow raw-probe odd-sector terms:
   - `collq_odd_rawProbe`
   - `realq_odd_rawProbe`

Relevant reports:

- `plain39` campaign:
  `/Users/apapaefs/trs/Projects/HerwigPol/HwPolNotesNew/DISPOL/campaigns/plain39/results.txt`
- Z TERMDIAG raw-probe report:
  `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/z_spin_report_rawprobe.txt`

## Main result

The raw-probe diagnostic rules out the leading clamp hypothesis.

The odd-sector decomposition in the Z-only TERMDIAG report is:

- `q_odd_pb = +0.004908 pb`
- `g_odd_pb = -0.005187 pb`

so the polarized Z NLO correction is controlled by a strong odd-sector
cancellation.

The key raw-probe result is:

- `q_odd_raw_shift_pb = -0.000510 pb`

with the split

- `cq_odd_raw_shift_pb = -0.000511 pb`
- `rq_odd_raw_shift_pb = +0.000001 pb`

Interpretation:

- replacing the clamped odd quark sector with the raw one would make
  `Z sigma_LL` smaller, not larger
- the shift is almost entirely in the quark collinear odd term `cq_odd`
- the mapped real-quark odd term `rq_odd` is essentially unaffected at the pb
  level relevant for the production discrepancy

So the clamp is not suppressing the positive odd contribution in the way needed
to explain the low `Z sigma_LL`. If anything, the clamp is slightly helping.

## Important diagnostic details

The report still shows large local raw/clamped changes in the mapped analyzing
power:

- POSNLO:
  - `clipQm_absw_frac = 0.027863`
  - `A_q_raw = 1.591789`
  - `A_q = 1.054444`
  - `A_q_shift = 0.537345`
- NEGNLO:
  - `clipQ_absw_frac = 0.236897`
  - `clipQm_absw_frac = 0.160088`
  - `A_q_raw = 1.300163`
  - `A_q = -0.222577`
  - `A_q_shift = 1.522739`

But these dramatic local changes do not translate into a large observable shift
in `rq_odd`. So large `A_q` excursions are not, by themselves, evidence for the
remaining production-level mismatch.

## What is ruled out

The following are no longer good explanations for the remaining pure-`Z`
discrepancy:

- missing `c4-c2`
- missing quark regular/singular coefficient-function terms
- the broad endpoint rewrite from `plain38`
- mapped quark odd-sector clamp suppression as the dominant cause

## Best current interpretation

The remaining discrepancy is probably not in the odd-sector clamp machinery.
The most plausible remaining targets are now:

1. a small pure-`Z` Born/LO axial normalization or convention mismatch that is
   already visible at LO
2. a broader pure-`Z` axial-sector normalization issue outside the odd clamp
3. an even-sector or analysis-level mismatch rather than an odd-sector one

## Recommended next debug step

Use `plain39`, not `results_testing15.csv`, as the comparison baseline.

The `2026-04-08` pure-`Z` LO A/B test was carried out with matched `plain40z`
and `plain41z` campaigns. The only intended physics change was enabling the
temporary LO rho-reconstruction toggle in `plain41z`.

Result:

1. LO `sigma0`: `1.2084085 -> 1.2084087 pb`
2. LO `sigma_LL`: `0.0542678 -> 0.0542712 pb`
3. LO `sigma0 - 00`: unchanged at `+0.001121 pb`

Conclusion:

1. the hadron-side longitudinal rho reconstruction is not the cause of the
   remaining pure-`Z` LO mismatch
2. the temporary LO rho toggle is not useful and should not be pursued further
3. the remaining suspects are still pure-`Z` normalization/convention issues or
   an analysis-level `sigma0` versus `00` closure problem

## One-line takeaway

The new raw-probe report shows that unclamping the odd quark sector would lower
`Z sigma_LL` by about `5.1e-4 pb`, so the clamp is not the source of the
remaining `plain39` pure-`Z` polarized discrepancy.
