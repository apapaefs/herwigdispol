# Plain39 Recovery Note

## Summary

`plain38` introduced a broad regression in the shared quark-collinear sector.
The dominant symptom was already visible in pure `GAMMA`, so the failure was
not specific to the NC odd term.

The production fix kept the NC-odd algebraic repair in
`DISBase.cc`:

```text
collq_odd = qOddResponse * Pq * collq_over_born_pol
```

and rolled back the broad finite-slope/probe rewrite of the shared quark
collinear kernels. The unpolarized and polarized `collq` formulas now use the
original pointwise expressions again.

## Campaign outcome

Using the `plain39` campaign results:

- `GAMMA` NLO agreement returned to the `plain37` level.
- `Z` and `ALL` also returned to the pre-`plain38` baseline.
- The large `plain38` POSNLO collapse disappeared.

Representative values:

- `GAMMA` NLO unpolarized: Herwig `2755.663 pb`, POLDIS `2756.009 pb`
- `GAMMA` NLO polarized: Herwig `43.595 pb`, POLDIS `43.66951 pb`
- `ALL` NLO unpolarized ratio: `0.999822`
- `ALL` NLO polarized ratio: `0.994784`

## Term diagnostics

Using `DIS-POL-POWHEG-ALL-TERMDIAG-plain39-check.txt`:

- `NLO_TERM_PATHO` count: `0`
- `NLO_TERM_ENDPOINT` count: `0`
- sampled `NLO_TERM_SIGN` lines all have `n_nan=0`

The mixed-helicity `NEGNLO` `VCG_neg` values are ordinary again:

- `PM-NEGNLO`: about `-4057`
- `MP-NEGNLO`: about `-4669`
- `PP-NEGNLO`: about `-4572`
- `MM-NEGNLO`: about `-5264`

So the rollback restored the good campaign agreement without reintroducing the
old endpoint catastrophe.
