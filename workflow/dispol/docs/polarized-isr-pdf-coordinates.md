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
Both macOS and native Linux paths have been exercised. On Linux an isolated
`--library-dir` also creates a private repository with that directory prepended.
ThePEG restores search paths from its saved repository after processing CLI
`-L`, so `-L` alone can silently select the installed shower. The runner checks
actual library initialization in `LD_DEBUG=libs`, and rejects wrong or multiple
shower libraries. The private repository changes no installed repository.

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
`stdout.log` and `loader.log` supply the underlying evidence; Linux isolated
tests also record the repository preparation command and `repository-build.log`.
A successful
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

On 24 October the Odysseus `herwig/pol` default was separately built and
validated using its native GCC 15.2 toolchain and existing QTilde objects.
Only the six-call translation unit differs between the paired builds. The old
library fails 105 of 410 checks, while the corrected library passes all 410,
with both loader selections verified. Twelve 100-event pp/DIS checks complete;
the pp ordinary-00, pp zero-polarization and DIS ordinary-00 old/new control
pairs are byte-identical. These DIS checks use `HardEmission None`; some logs
contain recoverable remnant or momentum-adjustment exceptions.

The corrected Odysseus library and source were then installed, with backups
of both original files. A further 410 checks pass through the normal installed
runtime without isolated-library arguments. The installed shower SHA-256 is
`1712ff3717af163d6f7bf847740f145d3de34daa1b7f3f74cee8308294de86f8`.
The executable, core Herwig, ThePEG and DIS libraries have unchanged hashes.
No HERMES production campaign was started. Its v5 inputs match the current
configuration, and a separate corrected-run command has been dry-run checked.

The local macOS production installation and other runtime lineages have not
been replaced. Map each intended runtime's actual executable, build source and
loaded shower, then repeat compatible compiled and bounded event checks.
The original full
failure-seed replay and observable-level comparisons remain separate follow-up
work. A corrected shower can change random decisions before the old failing
event; the saved-state fixture is the deterministic regression.
