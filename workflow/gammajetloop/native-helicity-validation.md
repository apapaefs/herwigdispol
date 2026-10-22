# Actual native helicity amplitudes versus MadLoop

On 22 October 2026, the amplitudes stored by the native `MEPP2GammaJetLoop`
in Herwig's `HardVertex` passed an independent MadLoop comparison at 68 Born
phase-space points. This tests the published production implementation at
`386e522e638cd8c186c0690646ea517756ca7f81`. No production ME change was needed.

The sample contains **1,088 complex transverse amplitudes** and **17,408
complex interference products**. It includes massive and exactly massless
five-flavour loops, physical polarized extractors, and bottom/top threshold
brackets. The four-event development pilot and one loader-verification event
are excluded from these counts.

## What is compared

`NativeHelicityAudit.cc` is an analysis plugin, compiled without either the ME
implementation or its analytic helper. It reads the existing shared native
`HardVertex`, all 81 amplitude entries, the four particles' stored vector
bases, incoming densities, outgoing vertex densities, momenta, masses and
couplings. Showering, hadronization and decays are disabled for this audit.
The observer also reevaluates the native `me2()` and compares it with the
colour-weighted contraction of the stored vertex. `lastME2()` is zero in this
HwMEBase workflow; it is recorded but is not used as a weight-validation gate.

`madloop_native_driver.f` evaluates those same positive-energy momenta in
MadGraph 3.5.17. It sets alphaS, alphaEM and the scale from each native event.
It exports both complex colour-flow coefficients and the **actual external
wavefunctions used by MadLoop**, including the momenta after its numerical
phase-space repair. Stable squared results are recorded separately from the
original-frame double-precision complex amplitude evaluation.

`compare_native_helicities.py` uses the Minkowski polarization overlaps

```
C_leg(h,m) = -eps_Herwig(h) . conj(eps_MadLoop(m))
M_MadLoop_in_Herwig_basis(h) = sum_m product_leg C_leg(h,m) * J_MadLoop(m)/2.
```

The overlaps determine helicity permutations and phases independently of the
amplitudes. Gauge differences proportional to each momentum are checked.
The two MadLoop colour coefficients must agree; each is twice Herwig's
coefficient of `d_abc`. Colour averaging is `5/24`, and the actual incoming
densities supply the spin contraction once. There is no additional symmetry
factor. The overall phase is fixed to **+1**, with no fitted phase or rescaling.
Every amplitude and every interference product is compared. Direct MadLoop
contractions also predict both native outgoing spin-density matrices.

The relative amplitude metric is
`abs(H-M)/max(abs(M), 1e-6*max_helicity_abs(M))`. The absolute-error floor was
inactive for all 1,088 amplitudes: the smallest amplitude was 0.00523 times
the largest at its event. Interference errors are divided by the square of
the largest amplitude at that event. Density errors are absolute entry errors.

## Results

| Sample | Points | Maximum relative amplitude error | Maximum scaled interference error | Required tolerance |
|---|---:|---:|---:|---:|
| Massive, ordinary extractor | 20 | 5.724e-7 | 5.472e-8 | 1e-6 |
| Massless5, top explicitly excluded | 16 | 1.252e-7 | 3.442e-8 | 1e-6 |
| Massive, polarized PP | 8 | 7.200e-8 | 2.949e-9 | 1e-6 |
| Massive, polarized PM | 8 | 3.329e-8 | 2.977e-9 | 1e-6 |
| Below bottom threshold | 4 | 4.230e-8 | 3.601e-9 | 1e-5 |
| Above bottom threshold | 4 | 1.788e-8 | 3.072e-9 | 1e-5 |
| Below top threshold | 4 | 2.583e-6 | 1.494e-8 | 1e-5 |
| Above top threshold | 4 | 8.212e-7 | 7.550e-9 | 1e-5 |

Across all samples, the maximum amplitude error divided by the largest
amplitude at the event is `5.472e-8`. The unpolarized squared result agrees to
`2.934e-9` relative, and the actual incoming-density-weighted result to
`3.513e-9`. Native `me2()` versus stored-vertex contraction agrees to
`7.395e-14`. All longitudinal amplitude entries are exactly zero.

The two outgoing density matrices agree with direct MadLoop contractions to
`8.037e-9`. They are Hermitian, have unit trace within `2.23e-16`, and their
transverse blocks have minimum eigenvalue `0.3613`. Nonzero off-diagonal
entries reach `0.1387`; this is not only a diagonal-density check. Including
the unused longitudinal state adds a zero eigenvalue.

As a negative control, multiplying one exported helicity in the pilot by
`exp(0.1 i)` leaves the helicity-summed rate unchanged but produces amplitude
and interference errors of `0.09996`: the comparison correctly fails.
Independent Python tests cover known basis phases, a helicity permutation,
gauge shifts, complex incoming coherence, an analytically known outgoing
density, and rejection of nonfinite input.

### A diagnostic phase artifact that was resolved

An initial exporter reconstructed MadLoop polarization vectors from the input
momenta. Native beam-transverse components can be tiny nonzero floating-point
numbers; MadLoop can round them to zero while restoring on-shell kinematics.
The helicity convention changes discontinuously on the beam axis, so those
reconstructed vectors did not describe the amplitudes actually evaluated.
Exporting `W(5:8,leg)` from MadLoop's generated wavefunction common block
resolved the discrepancy. The maximum momentum repair in the final sample
is only `6.322e-14` relative. This was an exporter issue, not a correction to
the Herwig ME or a fitted helicity phase.

## Reproduction

Use a compatible polarized Herwig/ThePEG environment with the published native
ME built into `HwMEHadron` and its default repository regenerated. This is
separate from the older preload-based helper validation. Set the following
paths for that installation; the isolated Odysseus paths are given below.

```bash
source /etc/profile.d/modules.sh
module load herwig/pol
task_source_repo=/absolute/path/to/herwigdispol
task_native_repo=/absolute/path/to/NativeDefaults.rpo
task_native_lib=/absolute/path/to/native-build/MatrixElement/Hadron/.libs
task_audit=/absolute/path/to/a/new/audit
bash "$task_source_repo/workflow/gammajetloop/build-helicity-audit.sh" \
  "$task_source_repo/src/herwig" "$task_audit/plugin"
mkdir "$task_audit/run"
cd "$task_audit/run"
cat > NativeAudit.in <<'EOF'
read NativeCommon.in
set /Herwig/MatrixElements/GammaJetLoopValidation:AmplitudeMode Massive
read NativeHelicityAudit.in
set /Herwig/Samplers/Sampler:BinSampler:InitialPoints 1000
saverun NativeAudit /Herwig/Generators/EventGenerator
EOF
Herwig read NativeAudit.in --repo="$task_native_repo" \
  -L "$task_native_lib" -L "$task_audit/plugin" \
  -I "$task_source_repo/workflow/gammajetloop"
Herwig run NativeAudit.run -N20 -s2042210 -q \
  -L "$task_native_lib" -L "$task_audit/plugin"
```

The output is `NativeAudit-S2042210.helicities`. Verify the executable with
`command -v Herwig` and the loaded Hadron library with a separate one-event
`LD_DEBUG=libs Herwig run ...` using a different seed. The repository file
selects the built native library; `-L` alone does not prove which library was
loaded. The source must match the native library used to generate the events.

Make a **new diagnostic copy** of each complete MadLoop export. Use
`loop_sm-c_mass` with `g g > a g QCD=3 QED=1 [noborn=QCD]` for Massive,
with mc=1.27, mb=4.18 and mt=172.69 GeV. For Massless5 use
`loop_sm-no_b_mass` and `g g > a g / t QCD=3 QED=1 [noborn=QCD]`,
with mc=mb=0. The older PHENIX `gg_massless5` campaign includes a top loop
and is unsuitable here. The [README](README.md) describes export generation.

Copy `madloop_native_driver.f` to `check_sa.f` in the diagnostic
`SubProcesses/PV0_0_1_gg_ag/` directory (`PV0_0_1_gg_ag_no_t/` for Massless5),
then run `make check`. The driver's common-block names and dimensions match
these MG 3.5.17 exports; check the generated Fortran declarations before using
it with another export/version. Populate `SubProcesses/MadLoop5_resources/`
`MadLoopParams.dat` with valid numerical values, not unexpanded placeholders.
The tested settings include `MLReductionLib=6|7|1`, `CTModeRun=-1`,
`MLStabThres=1e-7`, `NRotations_DP=1`, `ImprovePSPoint=2`,
`HelicityFilterLevel=1` and `DoubleCheckHelicityFilter=.false.`.

```bash
# Required when the diagnostic binary links the MG reduction libraries:
export LD_LIBRARY_PATH=/absolute/path/to/MG5_aMC_v3_5_17/HEPTools/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
python3 "$task_source_repo/workflow/gammajetloop/compare_native_helicities.py" \
  NativeAudit-S2042210.helicities --mode Massive \
  --madloop /absolute/path/to/diagnostic/SubProcesses/PV0_0_1_gg_ag \
  --output "$task_audit/comparison"
python3 "$task_source_repo/workflow/gammajetloop/test_native_helicity_comparison.py"
```

The comparison writes all 16 complex pairs per event, errors, provenance hashes
and raw MadLoop responses. A failed gate exits nonzero. `--reuse` reuses raw
responses only when their momenta, couplings, executable, parameter-card and
MadLoop-settings fingerprints match. It still repeats the numerical comparison.

For the full suite, the base settings are pp at 510 GeV,
NNPDF40_nlo_pch_as_01180, native alphaS and alphaEM(0), the native scale,
and particle-data masses 1.27/4.18/172.69 GeV. Use the following independent
run directories and overrides:

| Case | Events | Seed | Extra settings |
|---|---:|---:|---|
| massive | 20 | 2042210 | Base card above |
| massless5 | 16 | 2042211 | `AmplitudeMode Massless5` |
| polarized-pp | 8 | 2042212 | `PolarizedCommon.in`, beam polarizations +1,+1 |
| polarized-pm | 8 | 2042213 | `PolarizedCommon.in`, beam polarizations +1,-1 |
| bottom-below | 4 | 2042214 | 8.34 < MHat/GeV < 8.359, photon pT > 4 GeV |
| bottom-above | 4 | 2042215 | 8.361 < MHat/GeV < 8.38, photon pT > 4 GeV |
| top-below | 4 | 2042216 | 345.30 < MHat/GeV < 345.37, photon pT > 50 GeV |
| top-above | 4 | 2042217 | 345.39 < MHat/GeV < 345.46, photon pT > 50 GeV |

Set the windows through `/Herwig/Cuts/Cuts:MHatMin`, `MHatMax` and
`/Herwig/Cuts/PhotonKtCut:MinKT`, before `saverun`. For polarized cards, set
`/Herwig/Partons/PPPolarizedExtractor:FirstLongitudinalPolarization` and
`SecondLongitudinalPolarization`; `PolarizedCommon.in` uses
NNPDFpol20_nlo_as_01180. Pass `--tolerance 1e-5` for the four threshold cases.
The ordinary comparison tolerance is `1e-6`.

## Runtime and evidence

Odysseus uses `/home/apapaefs/Projects/Herwig/Herwig-Pol/bin/Herwig`.
The exact native build is under
`/home/apapaefs/Projects/HerwigPol/validation/gammajetloop-publication-20261021/`:

- Repository: `native-repository/NativeDefaults.rpo`.
- Native Hadron library: `native-build/MatrixElement/Hadron/.libs/HwMEHadron.so.11`.
- Library SHA256: `2fa9aac5df4c9fb22055a82e1df63d42e0e53e5d10df83d64e86e51ced5aa16b`.

The audit root is
`/home/apapaefs/Projects/HerwigPol/validation/gammajetloop-helicity-20261022/`.
It contains `runs/`, `comparisons/`, `cases.json`, `comparison-status.json`,
the isolated `audit-build-v2/NativeHelicityAudit.so`, and both diagnostic
MadLoop copies. The observer hash is
`18b56876f87ffbb1409c00d4f8eb0f9e2c080dbb9abe9ff3bd5264f90ecea3a7`.
The dynamic-loader trace confirms both library paths without a validation
preload. The active installation and completed campaigns are unchanged.

Local evidence is in the ignored
`workflow/gammajetloop/evidence/native-helicities-20261022/` directory.
The stored JSON reports include the full complex numbers, not only magnitudes.
These are pointwise hard-process tests. They do not constitute an exhaustive
phase-space proof, a new shower validation, or a production A_LL prediction;
the earlier QCD-shower and integrated-rate checks are documented separately.
