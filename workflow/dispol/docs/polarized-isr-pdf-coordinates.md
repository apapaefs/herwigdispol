# Polarized QTilde ISR PDF coordinates

The six polarized PDF queries in `SudakovFormFactor.cc` use `xfx` with ordinary
momentum fractions. Four calls in `PDFVetoRatio` evaluate the longitudinal and
transverse difference PDFs at the child's `x` and parent's `x/z`; two calls in
`calculateHMatrix` evaluate the incoming density at its supplied fraction.

`PDFBase::xfl` takes `log(1/x)`. Passing ordinary `x` through that interface to
LHAPDF evaluates the numerator at `exp(-x)` while the unpolarized denominator
is still evaluated at `x`. Because both interfaces return momentum-weighted
PDFs, this also introduces the wrong momentum prefactor. Using `xfx` for both
numerator and denominator restores the intended signed difference-to-sum ratio.
The extractor's existing `xfl(..., li())` calls use logarithmic fractions
correctly.

The correction changes only six method names. Scales, normalization, clipping,
failure guards, interfaces and persistence are unchanged. Correct coordinates
do not guarantee positivity for an arbitrary pair of polarized/unpolarized
PDF sets. High-x support failures need separate investigation.

## Regressions

Run from this workflow directory (`DISPOL` or `workflow/dispol`):

```bash
python3 scripts/tests/test_polarized_qtilde_sudakov.py -v
```

This includes a source-contract check for all six calls alongside the existing
normalization/scale checks. Numerical coverage is supplied by a compiled
fixture loaded into an actual `Herwig read` process:

```bash
python3 scripts/tests/run_polarized_qtilde_pdf_coordinates.py \
  --prefix "$HERWIG_TEST_PREFIX" \
  --cxx "$HERWIG_TEST_CXX" \
  --library-dir "$HERWIG_TEST_SHOWER_DIR" \
  --output "$HERWIG_TEST_OUTPUT" \
  --real-pdfs
```

Set these variables to the intended installation, ABI-compatible C++ compiler,
isolated directory containing `HwShower.so`, and a fresh diagnostic output
directory. Load the runtime's normal module first. Omitting `--library-dir`
tests the prefix's installed shower. Omitting `--real-pdfs` runs the synthetic
tests without requiring the two NNPDF sets. The runner supports macOS dyld and
glibc loader traces; loader verification must succeed for a passing result.
The macOS implementation has been exercised; the Linux branch needs validation
on the target host.

The fixture calls the real `PDFVetoRatio` and `calculateHMatrix`. Its PDF
subclass overrides only `xfx`, so an accidental `xfl` call goes through the
actual base-class conversion. It records fractions, flavours and scales;
tests constant and x/flavour-dependent signed ratios, q-to-qg and g-to-qqbar,
complex and real transverse coherence, terminal/parent densities, PDF-scale
freezing, and disabled/zero-polarization controls. A scoped current shower
handler supplies the scale policy, and tested acceptances stay below `PDFmax`.

With `--real-pdfs`, member 0 of `NNPDF40_nlo_pch_as_01180` and
`NNPDFpol20_nlo_as_01180` reproduces the saved September anti-u state. The
correct old-density contraction is `0.483881310072547` and the normalized
q-to-qg spin factor is `0.995714491077483`.

Require failure against an old library and success against a corrected library
built with the same toolchain and dependencies. `result.json` records the
fixture hash, executable hash, compile/run commands, loaded shower/ThePEG/DIS
paths and hashes, check count, failures and loader verification. `build.log`,
`stdout.log` and `loader.log` supply the underlying evidence. A successful
`Herwig read` exit alone does not establish that the checks passed.

## Integration and deployment

The 24 October local paired build used QTilde snapshot
`3079f722aa9100043cdf5effa173bd169b0c9a78`, changing only the six PDF calls.
It passes all 410 native assertions with the correction; the old build fails
105. The Python suite passes all 10 tests. Newer runtime lineages require
their own compatible build and loader checks.
Seventeen 100-event pp/DIS runs completed, with byte-identical old/corrected
genuine-00 and zero-polarization control files. DIS logs also show legacy
ME-correction warnings and recoverable remnant failures; these smoke runs
establish execution and bypass behavior, not production-level DIS validation.

The change is present in the active and curated source trees. The production
installation has not been replaced. Before a runtime rollout, map its actual
executable, build source and loaded shower, port the fix to that maintained
lineage, and repeat the compiled and bounded event checks. The original full
failure-seed replay and observable-level comparisons remain separate follow-up
work. A corrected shower can change random decisions before the old failing
event; the saved-state fixture is the deterministic regression.
