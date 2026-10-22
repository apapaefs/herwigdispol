# Native gluon-fusion photon plus jet

`MEPP2GammaJetLoop` is a separate native `HwMEBase` process for
`g g -> gamma g`. The source is in
`src/herwig/MatrixElement/Hadron/`; `Makefile.am` registers it
in `HwMEHadron`, and `src/defaults/MatrixElements.in` creates
`/Herwig/MatrixElements/MEGammaJetLoop` when the default repository is rebuilt.
The existing `MEGammaJet` process is unchanged.

## Interfaces and normalization

- `AmplitudeMode Flat`: all 16 reduced transverse amplitudes are one. This
  diagnostic is defined in the local helicity basis and is not a physical
  Lorentz-invariant scattering amplitude. Longitudinal entries are zero.
- `AmplitudeMode Massless5`: massless u,d,s,c,b, with the top loop omitted.
- `AmplitudeMode Massive` (default): massless u,d,s, and particle-data masses
  for c,b,t. The u,d,s contributions cancel coherently. No incoherent
  flavour sum or mass-dependent damping factor is used.
- `ScalePreFactor` multiplies `2*s*t*u/(s*s+t*t+u*u)`, as in `MEGammaJet`.

Let `h` denote the one-fermion reduced light-by-light amplitudes, with
`h(++++)=h(-+++)=1` in the massless limit. With `Tr(Ta Tb)=delta_ab/2`,
the sum of the two reversed quark-box colour traces is `d_abc/2`;
the corresponding photon trace sum is `2*Nc`. Thus

```
M_abc = 2 sqrt(alphaEM(0) alphaS(scale)^3) d_abc sum_q Q_q h(m_q).
sum_abc d_abc^2 / (8*8) = 5/24.
```

`me2()` applies `5/24` once. `ProductionMatrixElement::average(rho1,rho2)`
performs the incoming-spin contraction once; there is no further factor of
one quarter and no identical-final-state factor. The unpolarized flat result
is `(10/3)*alphaEM*alphaS^3`. Incoming gluon density matrices are projected
onto indices 0 and 2 and normalized, including their off-diagonal entries.
This handles both the generic spin-one `I/3` initialization and the physical
densities from `PolarizedPartonExtractor`.

Two oriented three-gluon colour flows are selected with equal probabilities.
The exact hard colour sum is independent of that shower colour-flow choice.

## Analytic sources, branches and spin phases

The massless functions reproduce the epsilon-zero coefficients of
[Bargiela, Chakraborty and Gambuti, arXiv:2212.14069](https://arxiv.org/abs/2212.14069),
including the complex crossed two-minus amplitudes. Their `H` coefficients
are twice `h` in the normalization above. The test includes independent
values extracted from all eight one-loop ancillary files at `x=0.3`.

The massive functions implement section 4.1 of
[Ajjath, Chaubey and Shao, arXiv:2312.16966](https://arxiv.org/abs/2312.16966).
The photon-only fermion-box kinematics are adapted using the colour and
coupling derivation above. The masters are `f2/sqrt(x*(x-4))`, `f4` and
`f6/sqrt(x*y*(x*y-4*(x+y)))`; the last equals `m^4*D0/2`.
The Appendix A regions I/II and physical `+i0` prescription determine the
logarithm/dilogarithm branches. Box symmetry removes the need for a separate
region III implementation. The zero-mass limit is explicit.

The helicity coefficients alone do not specify ThePEG's polarization phases.
Instead of assigning phases to helicity labels, the implementation reconstructs
the paper's gauge-invariant tensor (section 2) and contracts the actual
ThePEG wavefunctions. In four dimensions only `A1+A2+A3` is physical;
the implementation chooses `Ai=AS/3` and derives the matching B coefficients
from the five helicity functions. This replaces the initially proposed
spinor-overlap conversion with a covariant contraction. The same builder is
used for `me2()` and `constructVertex()`, including event-record basis vectors.
One shared `HardVertex` is attached to all four particles.

The external order is `(g,g,gamma,g)`. The tensor uses all incoming momenta
and ACS invariants `s=s12`, `t=s23`, `u=s13`; this differs from the BCG naming
`t=s13`. Physical outgoing momenta are negated explicitly.

For `max(|s|,|t|,|u|)/m^2 < 0.5`, an 18-term convergent expansion of the
same scalar integrals removes analytically vanishing constant/linear terms,
preserving heavy-quark decoupling. At exactly `x=4`, the bubble's singular
part multiplies zero and its finite contribution is used. If every
`|s_ij|/m^2 > 1e12`, the massless asymptote avoids subtracting large logarithmic
terms. This last branch is a small-mass approximation outside the PHENIX
charm/bottom/top validation domain. The local long-double dilogarithm uses
reflection, inversion and a convergent series, following Herwig's existing
math routine strategy without changing shared math code.

## Build and use

A normal source build, following the repository's top-level instructions,
includes this process in `HwMEHadron` and creates the default object. Rebuild
the default repository and regenerate run files when installing the change.
The checked-in `Makefile.in` includes the test target and was regenerated
with Automake 1.16.1, matching the source snapshot.
If a Git checkout triggers regeneration with a missing `aclocal-1.16`, run
`autoreconf -fi` in `src/herwig` with the available Autotools before configuring.

For a bounded smoke with the rebuilt installation, run from the repository
root (the named output directory must be new):

```bash
mkdir -p build
mkdir build/gammajetloop-smoke
cp workflow/gammajetloop/*.in build/gammajetloop-smoke/
cd build/gammajetloop-smoke
Herwig read Massive.in
Herwig run Massive.run -N100 -s123459 -q
```

The cards create a separate `GammaJetLoopValidation` object, so they also
work when the default `MEGammaJetLoop` object already exists. They require
NNPDF40_nlo_pch_as_01180; the polarized examples additionally require
NNPDFpol20_nlo_as_01180 and this repository's polarized ThePEG extractor.

### Isolated Linux validation

Load the compatible polarized Herwig environment, then run:

```bash
source /etc/profile.d/modules.sh
module load herwig/pol
bash workflow/gammajetloop/build-isolated.sh \
  "$PWD/src/herwig" /new/isolated/build
```

The script refuses an existing build directory, compiles the native class and
amplitude test, records SHA256 hashes and creates `herwig-loop`. The installed
default repository pins the original Hadron library, so isolated validation
uses `HwGammaJetLoopValidation.so` and an explicit `LD_PRELOAD` in that launcher.
The normal source build uses the `HwMEHadron` registration instead.
No active installation is overwritten.

On Odysseus the tested source snapshot and build are under:

```
/home/apapaefs/Projects/HerwigPol/validation/gammajetloop-20261021/
```

The verified launcher is `final-build-v3/herwig-loop` under that directory.
Its plugin SHA256 is
`f9de82d60e019ccbb3c067d23ab5f2fdf010ac4a8a22c430d96d102b535c0f34`;
the dynamic-loader trace confirms this exact isolated library is initialized.

Copy the `.in` files from this directory into a new working directory and use
the isolated launcher for both reading and running. For example:

```bash
/new/isolated/build/herwig-loop read Massive.in
/new/isolated/build/herwig-loop run Massive.run -N100 -s123459 -q
```

All shower cards explicitly use QCD and spin correlations. QED partner
handling for the neutral hard state is outside this change. The standalone
validation cards turn MPI off; the opt-in phenomenology cards keep the
PHENIX common-card MPI setting. `MassiveParton.in` disables shower,
hadronization and decays. `Fiducial.in` applies the photon window directly.

## Focused tests and independent comparison

The [native helicity audit](native-helicity-validation.md) compares the actual
complex amplitudes stored in Herwig's `HardVertex` with MadLoop, using each
generator's actual polarization vectors. Its 22 October validation covers
68 events, all 1,088 transverse amplitudes and 17,408 interference products,
including polarized beams and bottom/top thresholds. All checks pass without
fitting helicity phases or modifying the production ME. The older point-probe
comparison below instead exercises the analytic helper with MadLoop vectors.

The Hadron `make check` target runs `GammaJetLoopAmplitudesTest`. It can also
be compiled without Herwig dependencies using C++11 and
`GammaJetLoopAmplitudes.cc`. It covers published massless and high-precision
massive reference values, heavy masses, thresholds, Ward identities, gauge
shifts, rotations/boosts, parity, incoming exchange and helicity sums.

`test_hard_vertex.cc` compiles against ThePEG plus the existing
`ProductionMatrixElement.cc` and `HardVertex.cc` sources. It checks the flat
normalization, the four-pure-helicity average, and direct density contractions
with off-diagonal incoming states, including Hermiticity, positivity and trace.

After loading the compatible Herwig/ThePEG environment, run both tests and
build the point probe from the repository root:

```bash
bash workflow/gammajetloop/check-source.sh \
  "$PWD/src/herwig" "$PWD/build/gammajetloop-tests"
```

The script refuses an existing output directory. The scalar test and probe
use C++11; the HardVertex test uses C++17 to match the validated installed
Herwig/ThePEG headers. The optional `check_integrals.py` requires `mpmath` and
prints high-precision scalar references; `--quadrature` adds slow numerical
Feynman-parameter integrations below threshold.

`madloop_driver.f` replaces `check_sa.f` **only in a diagnostic copy** of a
MadGraph 3.5.17 loop-induced export. It reports stable squared results, then
reevaluates without stability rotations to export the complex colour
amplitudes in the original frame. Its eight-component loop-vector convention
is important. `compare_madloop.py` compares amplitudes, all interference
products and the absolute squared result; it allows only one common phase.
The MG colour-flow coefficient is twice the `d_abc` coefficient.

Use the `point-probe` built by `check-source.sh`:

```bash
python3 workflow/gammajetloop/compare_madloop.py \
  --probe build/gammajetloop-tests/point-probe /path/to/massive/point-*.txt
python3 workflow/gammajetloop/compare_madloop.py \
  --probe build/gammajetloop-tests/point-probe --masses 0 0 -1 /path/to/massless/point-*.txt
```

Generate the massive reference with `import model loop_sm-c_mass` and
`generate g g > a g QCD=3 QED=1 [noborn=QCD]`, then
`output DIRECTORY -nojpeg`. Set the parameter-card masses to mc=1.27,
mb=4.18 and mt=172.69 GeV, with u,d,s massless. In the diagnostic export's
`SubProcesses/PV0_0_1_gg_ag/` directory, copy `madloop_driver.f` to
`check_sa.f`, run `make check`, and provide `E cos(theta) phi` on stdin:

```bash
./check > point-10.txt <<EOF
10 0.4 0.7
EOF
```

The eight massive benchmark inputs are `(10,0.4,0.7)`, `(4.17999,0.4,0.2)`,
`(4.18001,0.4,0.2)`, `(30,-0.6,1.2)`, `(100,0.99,0.3)`,
`(172.6899,0.3,0.8)`, `(172.6901,0.3,0.8)` and `(250,-0.4,1.7)`;
energies are in GeV and angles in radians. The checker exports its actual
couplings and polarizations, which the comparison uses directly.

For the exact massless reference use `loop_sm-no_b_mass` and
`generate g g > a g / t QCD=3 QED=1 [noborn=QCD]`, followed by ordinary
loop-induced `output DIRECTORY -nojpeg`. The standalone exporter in this MG
version failed. The generated placeholder `MadLoopParams.dat` must be populated
with valid numerical settings before using the checker. The original
`gg_massless5` campaign includes a massive top and is a different definition.

## Validated results, 21 October 2026

- Publication snapshot: the complete Hadron module builds from the shipped
  build templates and its `make check` passes on Odysseus. The standalone
  amplitude and HardVertex tests pass. All 11 example cards load, and Massive
  completes 100 QCD-showered events through both the isolated launcher and a
  fresh default repository using the rebuilt native `HwMEHadron.so.11`.
  The latter uses no validation preload; a loader trace verifies the new
  library path. These checks are isolated from the active installation.
- Eight massive MadLoop points: maximum relative squared-amplitude difference
  `8.5e-11`; maximum amplitude/interference error relative to the amplitude
  norm `1.1e-7`. Points bracket bottom/top thresholds and include forward
  scattering and nonzero azimuths.
- Two exact massless-five-flavour points: amplitude/interference errors below
  `3.3e-12`.
- Flat, Massless5 and Massive: 100 QCD-showered events each, including persisted
  run-file reloads. Four physical beam-helicity configurations and the
  polarized extractor with zero polarization also pass 100 events each.
- Combined tree/loop QCD sample: 3,000 events, including three loop events.
  This is a plumbing check, not a precise loop-fraction estimate.
- The delivered phenomenology Massive card also passes 100 events with MPI
  and PHENIX Rivet: eight inclusive-photon and seven isolated-photon entries.
  All three opt-in phenomenology cards pass the read/setup stage.
- Generation cuts: 20,000 native hard events give **935 +/- 4 pb**, versus
  **938.95 +/- 0.26175 pb** from the audited MadGraph campaign.
- PHENIX photon window: 5,000 native hard events give **18.0 +/- 0.1 pb**,
  versus **18.3189 +/- 0.4107 pb**. Native values/errors are rounded as printed
  by Herwig. These are parton-level comparisons, not isolated showered rates.

The synchronized setup is pp at 510 GeV, NNPDF40_nlo_pch_as_01180 member 0,
the native Herwig hard alphaS/alphaEM and scale, photon pT above 4 GeV and
|eta| below 3. The fiducial window is 6<pT<30 GeV and |eta|<0.25. Quark masses
are 1.27, 4.18 and 172.69 GeV. No production A_LL result is inferred from the
small beam-helicity smokes. The contribution is LO loop-induced
alphaEM*alphaS^3, not a complete NNLO prompt-photon prediction.

Generated evidence is stored locally under the ignored `evidence/` directory
and in the isolated Odysseus directory; it is not part of the source checkout.
Existing campaigns, the active Herwig installation and shared shower code
were not modified.
