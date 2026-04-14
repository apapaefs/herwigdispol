# LO `GAMMA` Cut-Surface Audit

## Summary

The remaining LO `GAMMA` normalization mismatch is not coming from the local
DIS hard process. The local `LO_GAMMA_POINT` audit validated the Herwig
pointwise PDF lookup, hard `sigma_hat`, and local `HwMEBase` Jacobian against
the standalone LO calculator, and the `.out` extractor reproduces the raw
Herwig shard totals exactly. The surviving discrepancy is cut-window dependent.
The leading code-level suspect is still the non-equivalent ThePEG cut
treatment in `SimpleDISCut`, especially the fact that `MinQ2` shapes the phase
space while `MaxQ2`, `Maxy`, and `W2` are applied later as post-cuts. The
later cut-relaxation tests show, however, that the effect is broader than a
simple one-sided leak at only the upper cut boundaries.

## What Is Already Ruled Out

- The local LO `GAMMA` hard process is validated. The `LO_GAMMA_POINT` audit
  showed point-by-point agreement between Herwig and the standalone LO
  calculator for:
  - `pdf_sum`
  - local `sigma_hat_nb`
  - the local `HwMEBase` Jacobian
- The `.out` extraction and recombination are validated. The raw shard totals
  recombine exactly to the campaign-reported numbers for at least:
  - LO `00-GAMMA`
  - LO `00-Z`
- The `.log` integrated cross section and the `.out` attempted/generated
  cross sections agree for the same LO `GAMMA` shard. So the bias is already
  present before final `.out` formatting.
- `BinSampler:InitialPoints` is a null test. Matched-seed runs with
  `InitialPoints = 1e6` and `InitialPoints = 1e7` gave identical LO
  `00-GAMMA` shard outputs and identical combined means.
- `MonacoSampler` is also a null test. It gave the same LO `00-GAMMA`
  combined result as the default sampler.
- The quark `HardProcessMass = 0` test did not explain the original
  discrepancy. It created an artificial reshuffling effect and moved the
  plain-run cross sections in the wrong way.

## Code Path

The relevant cut/application path is:

- first hard-point cut:
  [`HwMEBase::generateKinematics()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/HwMEBase.cc#L197)
- integration-time second cut:
  [`StandardXComb::willPassCuts()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/Handlers/StandardXComb.cc#L261)
- later subprocess-level cut:
  [`StandardXComb::construct()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/Handlers/StandardXComb.cc#L794)
- subprocess construction:
  [`Tree2toNDiagram::construct()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/MatrixElement/Tree2toNDiagram.cc#L44)
- reshuffling trigger:
  [`StandardXComb::checkReshufflingNeeds()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/Handlers/StandardXComb.cc#L675)
- DIS cut implementation:
  [`SimpleDISCut::minTij()` / `passCuts()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/Cuts/SimpleDISCut.cc#L53)

The key `SimpleDISCut` detail is:

- `MinQ2` enters the phase-space generation through `minTij()`
- `MaxQ2`, `Miny`, `Maxy`, and `W2` are enforced later in `passCuts()`
- `passCuts()` computes:
  - `Q2 = -(p_in - p_out)^2`
  - `x = min(1, sqrt(currentSHat / SMax) * exp(+/- currentYHat))`
  - `y = Q2 / (SMax * x)`
  - `W2 = (1 - x) Q2 / x`
- the inequalities are strict: `>` / `<`, not inclusive

So the upper cut surfaces are treated as post-cuts rather than generation
variables. This differs from POLDIS, which generates directly inside the
`Q^2,y` window.

## Momentum Changes Between Cut Stages

There are two distinct "second cut" questions, and they matter differently.

For the integrated LO cross section, the relevant second cut is
[`StandardXComb::willPassCuts()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/Handlers/StandardXComb.cc#L261),
not the later subprocess-level cut in `construct()`.

In that integration-time path:

- `HwMEBase::generateKinematics()` applies the first cut on the generated
  `meMomenta()`
- `willPassCuts()` reapplies the cut on the same `meMomenta()`
- `willPassCuts()` may boost the outgoing momenta to the incoming-parton CM,
  but for LO DIS they are already in that frame
- both passes use the same `currentSHat/currentYHat`

So for the integrated LO cross section, the first and second cut applications
should agree up to negligible roundoff, and there is no evidence for a hidden
double-veto from different momenta in that path.

In the later subprocess/event-construction path, yes, the kinematics can still
change.

- The subprocess is rebuilt in
  [`Tree2toNDiagram::construct()`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/MatrixElement/Tree2toNDiagram.cc#L44)
  from `xc.meMomenta()`.
- If `StandardXComb::needsReshuffling()` is true, then `xc.reshuffle(pout)` is
  applied before the second cut.
- `needsReshuffling()` turns on when an outgoing particle has
  `mass() != hardProcessMass()`.

For the cleaned plain runs, the large reshuffling effect should be absent:

- the `HardProcessMass` / nominal-mass mismatch was removed with
  `AdjustNominalMass`
- showers are off in the plain cards
- the subprocess is rebuilt from `xc.meMomenta()`, not from
  `HwMEBase::rescaledMomenta()`

`DISBase::generateKinematics()` itself does not look suspicious for LO. For
`contrib_ == 0`, it mostly delegates to `HwMEBase::generateKinematics()` and
then caches `xB_` and `q2_`. The extra DIS remapping in that function is
NLO-only.

However, even with reshuffling off, there is still a CM/lab/CM round-trip
during subprocess construction. So small floating-point differences in the
outgoing lepton momentum remain possible in the later event record path.

This matters because `SimpleDISCut::passCuts()` mixes:

- `Q^2` from the outgoing lepton momentum
- `x` from `currentSHat/currentYHat`

At fixed `x`, even a small change in the outgoing lepton momentum shifts
`Q^2`, `y`, and `W^2`. But the current code tracing now suggests that this is
not the main explanation for the already-low integrated LO `.log` cross
section.

## Numerical Evidence

### Full window

LO `GAMMA`, plain-style window `[49,2500] x [0.2,0.6]`:

- standalone: `2952.415038984 pb`
- Herwig control run: `2951.260 +- 0.200 pb`

### Interior window

LO `GAMMA`, interior window `[100,1000] x [0.3,0.5]`:

- standalone: `543.4648337017 pb`
- Herwig: `543.470 +- 0.040 pb`

This is essentially perfect agreement.

### One-edge-at-a-time scans

Raise only `Q2Min` to `100`:

- standalone `[100,2500] x [0.2,0.6]`: `1188.320306944 pb`
- Herwig: `1187.681 +- 0.070427267447 pb`

Raise only `yMin` to `0.3`:

- standalone `[49,2500] x [0.3,0.6]`: `1870.561467278 pb`
- Herwig: `1869.190 +- 0.100 pb`

Lower only `Q2Max` to `1000`:

- standalone `[49,1000] x [0.2,0.6]`: `2920.705097461 pb`
- Herwig: `2920.530 +- 0.200 pb`

Lower only `yMax` to `0.5`:

- standalone `[49,2500] x [0.2,0.5]`: `2466.582880574 pb`
- Herwig: `2465.720 +- 0.200 pb`

### High-`Q^2` strip checks

The high-`Q^2` strip itself agrees well:

- `[1000,2500] x [0.2,0.6]`
  - standalone: `31.70993273651 pb`
  - Herwig: `31.7018 +- 0.0020 pb`
- `[1000,2500] x [0.2,0.5]`
  - standalone: `25.48671642419 pb`
  - Herwig: `25.4822 +- 0.0020 pb`
- `[1000,2500] x [0.5,0.6]`
  - standalone: `6.223216501124 pb`
  - Herwig: `6.22292 +- 0.00050 pb`

### Split-window non-additivity

The sharper evidence is that Herwig sub-runs do not recombine to the direct
full-window run.

`Q^2` split:

- Herwig `[49,1000] x [0.2,0.6]`: `2920.530 pb`
- Herwig `[1000,2500] x [0.2,0.6]`: `31.7018 pb`
- sum: `2952.2318 pb`
- direct Herwig full window `[49,2500] x [0.2,0.6]`: `2951.260 pb`

`y` split:

- standalone `[49,2500] x [0.2,0.5]`: `2466.582880574 pb`
- standalone `[49,2500] x [0.5,0.6]`: `485.8321570761 pb`
- standalone sum: `2952.4150376501 pb`

- Herwig `[49,2500] x [0.2,0.5]`: `2465.720 pb`
- Herwig `[49,2500] x [0.5,0.6]`: `484.948 +- 0.040 pb`
- Herwig sum: `2950.668 pb`

So the sub-runs do not recombine cleanly to the direct full-window behavior.

### Small cut-relaxation tests

These tests were meant to check whether the original deficit disappears when
the nominal cut surfaces are relaxed slightly.

Widen only the upper post-cut surfaces:

- standalone `[49,2501] x [0.2,0.601]`: `2956.821689873 pb`
- Herwig: `2955.770 +- 0.200 pb`
- gap: `-1.0517 pb`

Lower only the low-`y` boundary:

- standalone `[49,2500] x [0.199,0.6]`: `2965.502438869 pb`
- Herwig: `2964.350 +- 0.200 pb`
- gap: `-1.1524 pb`

Widen all three tested boundaries:

- standalone `[49,2501] x [0.199,0.601]`: `2969.909093042 pb`
- Herwig: `2968.910 +- 0.200 pb`
- gap: `-0.9991 pb`

These runs track the expected increase in phase space, but the baseline
`~1 pb` low-side offset survives almost unchanged. In particular:

- widening `MaxQ2` by `1 GeV^2` contributes almost nothing on the standalone
  side
- widening `Maxy` by `0.001` contributes a few pb
- lowering `Miny` by `0.001` contributes about `13 pb`

So the relaxation tests do not support the idea that the full discrepancy is
caused by a tiny loss at a single sharp boundary.

### `2 x 2` tiling of the full window

The full window was tiled into:

- `[49,1000] x [0.2,0.5]`
- `[49,1000] x [0.5,0.6]`
- `[1000,2500] x [0.2,0.5]`
- `[1000,2500] x [0.5,0.6]`

Herwig results:

- `[49,1000] x [0.2,0.5]`: `2440.670 +- 0.200 pb`
- `[49,1000] x [0.5,0.6]`: `479.549 +- 0.040 pb`
- `[1000,2500] x [0.2,0.5]`: `25.4822 +- 0.0020 pb`
- `[1000,2500] x [0.5,0.6]`: `6.22292 +- 0.00050 pb`

Standalone references:

- `[49,1000] x [0.2,0.5]`: `2441.096147136 pb`
- `[49,1000] x [0.5,0.6]`: `479.6089441879 pb`
- `[1000,2500] x [0.2,0.5]`: `25.48671642419 pb`
- `[1000,2500] x [0.5,0.6]`: `6.223216501124 pb`

So the tiles are individually close to the standalone values.

Two important consistency checks:

- low-`Q^2` slab from tiles:
  - `2440.670 + 479.549 = 2920.219 pb`
  - direct Herwig `[49,1000] x [0.2,0.6]`: `2920.530 pb`
  - compatible within the quoted run errors
- high-`Q^2` slab from tiles:
  - `25.4822 + 6.22292 = 31.70512 pb`
  - direct Herwig `[1000,2500] x [0.2,0.6]`: `31.7018 pb`
  - also compatible within the quoted run errors

But the broad high-`y` strip is not stable:

- tile sum for `[49,2500] x [0.5,0.6]`:
  - `479.549 + 6.22292 = 485.77192 pb`
- direct Herwig `[49,2500] x [0.5,0.6]`: `484.948 +- 0.040 pb`
- standalone `[49,2500] x [0.5,0.6]`: `485.8321570761 pb`

So the same high-`y` physical strip is close to the standalone result when it
is reconstructed from narrow `Q^2` tiles, but low when it is run directly over
the full `Q^2` span.

The full `2 x 2` tile sum is:

- Herwig tile sum:
  - `2440.670 + 479.549 + 25.4822 + 6.22292 = 2951.92412 pb`
- direct Herwig full window:
  - `2951.260 pb`
- standalone full window:
  - `2952.415038984 pb`

So the tile sum is substantially closer to the standalone full-window result
than the direct broad-window Herwig run.

## Current Interpretation

- `W^2` does not look like the active culprit in the standard plain runs.
  The plain cards set `MinW2 = 0`, `MaxW2 = 1E99`, so the explicit `W^2` cut
  is effectively inactive.
- In the standard plain window, `W^2 = y S - Q^2` stays comfortably positive
  throughout the physical region.
- The high-`Q^2` strip itself agrees well, so the discrepancy is not simply
  located in the high-`Q^2` corner.
- `DISBase::generateKinematics()` does not look suspicious for LO.
- The local LO DIS hard process is fine.
- For the integrated LO cross section, the first and second cut applications
  appear to use the same phase-space point in practice, so a hidden
  momentum-drift double-veto in that path is disfavored.

The most plausible remaining issue is the non-symmetric cut treatment in the
ThePEG cut path:

- `MinQ2` shapes the phase space through `minTij()`
- `MaxQ2`, `Maxy`, and `W2` are applied later as post-cuts
- `SimpleDISCut` recomputes `Q^2`, `y`, and `W^2` from the outgoing lepton
  momentum while holding `x` fixed from `currentSHat/currentYHat`
- but the LO deficit now looks less like a point-by-point cut mismatch and more
  like a broad-window integration or representation issue in the common
  `StandardXComb`/cuts path

But the later cut-relaxation tests show that the simplest "upper-boundary
leakage" explanation is not sufficient:

- widening only `MaxQ2` and `Maxy` does not remove the deficit
- lowering only `Miny` also does not remove the deficit
- the runs respond to real extra phase space about as expected, yet the
  underlying `~1 pb` offset remains
- the `2 x 2` tiles are individually much more stable than some of the direct
  broader windows

So the working interpretation is now:

- the local LO hard physics is validated
- the remaining plain-run LO `GAMMA` discrepancy is strongly cut-window
  dependent
- the issue is broader than a tiny one-edge leak and is more likely a
  full-window cut/integration representation problem
- the direct broad-window runs are the unstable objects; the narrow `Q^2` tiles
  behave much better
- the leading suspect remains the non-equivalent `SimpleDISCut` treatment of
  generated versus post-cut variables, now viewed through the broader
  `StandardXComb` integration/cuts path rather than a local matrix-element or
  PDF problem

## Campaign Outcome: `plain43` -> `plain44`

`plain44` enabled the native DIS-window generation path in the production
plain cards. The broad-window `GAMMA` problem that motivated this note is
therefore largely resolved at campaign level.

### `GAMMA`

The main broad-window LO discrepancy essentially disappeared.

- `plain43` LO unpolarized:
  - Herwig `2950.945 +- 0.212 pb`
  - POLDIS `2952.499 +- 0.200 pb`
  - delta `-1.554 +- 0.292 pb` (`5.32 sigma`)
- `plain44` LO unpolarized:
  - Herwig `2952.340 +- 0.141 pb`
  - POLDIS `2952.499 +- 0.200 pb`
  - delta `-0.159 +- 0.245 pb` (`0.65 sigma`)

So the native-window change moved the LO `GAMMA` direct broad-window result by
about `+1.395 pb`, which is exactly the direction suggested by the earlier
tiling studies.

The NLO `GAMMA` total also moved:

- `plain43` NLO unpolarized:
  - delta `-0.361 +- 0.311 pb`
- `plain44` NLO unpolarized:
  - delta `+0.836 +- 0.268 pb`

So the broad-window NLO `GAMMA` comparison remains statistically acceptable,
but it now sits slightly high instead of slightly low.

### `ALL`

`ALL` also improved materially.

- `plain43` LO unpolarized `sigma0`:
  - Herwig `2964.379 +- 0.106 pb`
  - POLDIS `2965.621 +- 0.201 pb`
  - delta `-1.242 +- 0.227 pb` (`5.48 sigma`)
- `plain44` LO unpolarized `sigma0`:
  - Herwig `2965.2188 +- 0.0707 pb`
  - POLDIS `2965.621 +- 0.201 pb`
  - delta `-0.402 +- 0.213 pb` (`1.89 sigma`)

The polarized/interference sectors also became much cleaner:

- `plain44` NLO polarized interference:
  - Herwig `4.564 +- 0.123 pb`
  - POLDIS `4.56915 +- 0.00441 pb`
  - delta `-0.005 +- 0.123 pb` (`0.04 sigma`)

### `Z`

`Z` improved in some internal consistency measures, but a residual systematic
offset remains.

- `plain43` LO unpolarized `sigma0`:
  - delta `-0.001302 +- 0.000131 pb`
- `plain44` LO unpolarized `sigma0`:
  - delta `+0.000692 +- 0.000122 pb`

- `plain43` NLO unpolarized `sigma0`:
  - delta `-0.001336 +- 0.000127 pb`
- `plain44` NLO unpolarized `sigma0`:
  - delta `+0.000768 +- 0.000121 pb`

So the sign flips, but the discrepancy remains at the few-`10^-4 pb` level and
is still statistically significant in `Z`.

The polarized `Z` LO quantity is now essentially perfect:

- `plain44` LO polarized `sigma_LL`:
  - Herwig `0.0546527 +- 0.0000280 pb`
  - POLDIS `0.0546560 +- 0.0000100 pb`
  - delta `-0.0000033 +- 0.0000298 pb` (`0.11 sigma`)

The `Z` NLO polarized quantity improved relative to `plain43`, but still shows
a residual offset:

- `plain43` NLO polarized `sigma_LL`:
  - delta `-0.0010055 +- 0.0000472 pb` (`21.28 sigma`)
- `plain44` NLO polarized `sigma_LL`:
  - delta `-0.0002440 +- 0.0000293 pb` (`8.33 sigma`)

### Interpretation after `plain44`

- The native DIS-window generation change fixed the main broad-window LO
  `GAMMA` failure that motivated this investigation.
- The same change materially improved `ALL`.
- The `Z` sector is not fully cured; the residual discrepancy is now much
  smaller in absolute size, but still statistically significant because the
  `Z` numbers are so precise.
- Internal helicity consistency is cleaner in `plain44`, for example:
  - `GAMMA` LO `(PP+PM)/2 - 00 = +0.047 +- 0.173 pb`
  - `Z` LO `sigma0 - 00 = +0.0002120 +- 0.0000684 pb`
  - `Z` NLO `sigma0 - 00 = -0.0000108 +- 0.0000597 pb`

So the broad cut-window pathology documented above appears to have been the
dominant source of the old `GAMMA` plain-run mismatch. The remaining work after
`plain44` is no longer to explain the original `~1.5 pb` LO `GAMMA` deficit,
but to understand the smaller residual `Z`-sector offsets.

It is important not to conflate the residual `plain44` offsets with the
original pathology.

The original problem had a distinctive pattern:

- broad-window non-additivity
- a systematic low bias
- strong `GAMMA` sensitivity
- and a direct cure when the full DIS window was moved into generation

The remaining `plain44` residuals look different:

- they are much smaller in absolute size
- the sign flips relative to `plain43` in several channels
- `GAMMA` NLO is now slightly high rather than low
- `Z` also shifts from low to high in the unpolarized totals
- there is not yet evidence that the same broad-window non-additivity survives

So the post-`plain44` residuals should be treated as a different class of
effect. The dominant broad-window generation/post-cut mismatch appears to have
been fixed.

## Post-`plain44` Next Steps

The next investigation target is no longer the old broad-window LO DIS failure.
It is the smaller NLO/`Z` residual structure that remains after the generation
fix.

Recommended order:

- check whether the native-generation `GAMMA` NLO and `Z` residuals still show
  any measurable broad-window non-additivity
- test sensitivity to the existing NLO `xp_` sampling map through
  `SamplingPower`
- only if those residuals are clearly sensitive to the NLO radiative map,
  consider a phase-2 reparameterization of `xp_`

In other words:

- the native Born-window fix appears to have solved the original dominant issue
- the next phase is a smaller NLO-level/electroweak-structure study, not a
  continuation of the same broad-window LO debugging

## Campaign Outcome: `plain44` -> `plain45` (`SamplingPower = 0.2`)

`plain45` kept the native DIS-window generation from `plain44` and changed only
the NLO sampling map through `SamplingPower = 0.2`.

The most important structural result is:

- LO is unchanged relative to `plain44`, as expected
- only the NLO pieces move
- so `SamplingPower` is probing the NLO `xp_` map rather than disturbing the
  native Born-window fix

### `GAMMA`

For `GAMMA`, lowering `SamplingPower` improves the NLO unpolarized comparison
but worsens the NLO polarized one.

- `plain44` NLO unpolarized `sigma_00`:
  - delta `+0.836 +- 0.268 pb` (`3.12 sigma`)
- `plain45` NLO unpolarized `sigma_00`:
  - delta `+0.683 +- 0.268 pb` (`2.55 sigma`)

- `plain44` NLO polarized `(PP-PM)/2`:
  - delta `-0.125 +- 0.101 pb` (`1.25 sigma`)
- `plain45` NLO polarized `(PP-PM)/2`:
  - delta `-0.220 +- 0.100 pb` (`2.19 sigma`)

The LO `GAMMA` success remains intact:

- `plain45` LO unpolarized `sigma_00`:
  - Herwig `2952.340 +- 0.141 pb`
  - POLDIS `2952.499 +- 0.200 pb`
  - delta `-0.159 +- 0.245 pb` (`0.65 sigma`)

### `Z`

For `Z`, the `SamplingPower = 0.2` change gives only modest improvement.

- `plain44` NLO unpolarized `sigma0`:
  - delta `+0.000768 +- 0.000121 pb` (`6.34 sigma`)
- `plain45` NLO unpolarized `sigma0`:
  - delta `+0.000758 +- 0.000121 pb` (`6.25 sigma`)

- `plain44` NLO polarized `sigma_LL`:
  - delta `-0.0002440 +- 0.0000293 pb` (`8.33 sigma`)
- `plain45` NLO polarized `sigma_LL`:
  - delta `-0.0002162 +- 0.0000297 pb` (`7.27 sigma`)

So `Z` moves in the expected direction for both NLO comparisons, but the
residuals remain statistically significant.

### `ALL`

For `ALL`, the change is mixed in the same way:

- `plain44` NLO unpolarized `sigma0`:
  - delta `+0.637 +- 0.238 pb` (`2.67 sigma`)
- `plain45` NLO unpolarized `sigma0`:
  - delta `+0.684 +- 0.238 pb` (`2.87 sigma`)

- `plain44` NLO polarized `sigma_LL`:
  - delta `-0.1303 +- 0.0711 pb` (`1.83 sigma`)
- `plain45` NLO polarized `sigma_LL`:
  - delta `-0.0972 +- 0.0709 pb` (`1.37 sigma`)

The `ALL` interference sector also shifts:

- `plain44` NLO polarized interference:
  - delta `-0.005 +- 0.123 pb` (`0.04 sigma`)
- `plain45` NLO polarized interference:
  - delta `+0.123 +- 0.123 pb` (`1.00 sigma`)

### Interpretation after `plain45`

`plain45` is useful even though it is not a clean across-the-board
improvement.

The main conclusion is:

- the remaining post-`plain44` NLO residuals are genuinely sensitive to the
  `xp_` sampling map
- but `SamplingPower = 0.2` is not a universal improvement over the old
  setting

So the `plain45` result strengthens the case that `xp_` is a live part of the
remaining NLO discrepancy story, but it does not identify `0.2` as a new
production default.

The practical reading is:

- the native Born-window fix from `plain44` remains the right production
  baseline
- the next `xp_` study should be a small map scan, not another immediate
  production switch

## Campaign Outcome: `plain45` -> `plain46` (finite-width spacelike `Z`)

`plain46` kept the native DIS-window generation from `plain44`, removed the
temporary `SamplingPower = 0.2` scan setting from `plain45`, and turned on the
finite-width spacelike `Z` propagator option in
`MENeutralCurrentDIS`.

The main effect is exactly where expected: the unpolarized `Z` normalization
shifts onto the POLDIS-like convention, while `GAMMA` stays unchanged.

### `GAMMA`

As expected, the `GAMMA` channel is essentially unchanged relative to
`plain44`:

- LO unpolarized `sigma_00`:
  - `plain46`: `2952.340 +- 0.141 pb`
  - POLDIS: `2952.499 +- 0.200 pb`
  - delta: `-0.159 +- 0.245 pb` (`0.65 sigma`)
- NLO unpolarized `sigma_00`:
  - `plain46`: `2756.845 +- 0.142 pb`
  - POLDIS: `2756.009 +- 0.227 pb`
  - delta: `+0.836 +- 0.268 pb` (`3.12 sigma`)
- LO polarized `(PP-PM)/2`:
  - `plain46`: `45.352 +- 0.100 pb`
  - POLDIS: `45.42693 +- 0.00280 pb`
  - delta: `-0.074 +- 0.100 pb` (`0.74 sigma`)
- NLO polarized `(PP-PM)/2`:
  - `plain46`: `43.544 +- 0.100 pb`
  - POLDIS: `43.66951 +- 0.00300 pb`
  - delta: `-0.125 +- 0.101 pb` (`1.25 sigma`)

So the finite-width spacelike-`Z` switch is not introducing any visible
cross-talk into the `GAMMA` control.

### `Z`

This is the clearest campaign-scale confirmation of the spacelike-`Z`
propagator interpretation.

Unpolarized `Z`:

- LO unpolarized `sigma0`:
  - `plain45`: `1.2111630 +- 0.0000280 pb`
  - `plain46`: `1.2103864 +- 0.0000280 pb`
  - POLDIS: `1.210471 +- 0.000119 pb`
  - `plain46` delta: `-0.000085 +- 0.000122 pb` (`0.69 sigma`)
- NLO unpolarized `sigma0`:
  - `plain45`: `1.1543131 +- 0.0000275 pb`
  - `plain46`: `1.1535932 +- 0.0000275 pb`
  - POLDIS: `1.153555 +- 0.000118 pb`
  - `plain46` delta: `+0.000038 +- 0.000121 pb` (`0.32 sigma`)

So the large `plain44/plain45` unpolarized `Z` discrepancy was indeed mostly a
propagator-convention issue rather than a broad-window integration failure.

The residual polarized/NLO `Z` offset remains:

- LO polarized `sigma_LL`:
  - `plain46`: `0.0546244 +- 0.0000280 pb`
  - POLDIS: `0.0546560 +- 0.0000100 pb`
  - delta: `-0.0000316 +- 0.0000297 pb` (`1.06 sigma`)
- NLO polarized `sigma_LL`:
  - `plain46`: `0.0498882 +- 0.0000275 pb`
  - POLDIS: `0.0501500 +- 0.0000100 pb`
  - delta: `-0.0002618 +- 0.0000293 pb` (`8.93 sigma`)

This is important because the finite-width propagator fixes the unpolarized
normalization very cleanly but does **not** fix the polarized NLO `Z`
residual. In fact, the polarized NLO `Z` comparison gets slightly worse than
in `plain45`.

Internal helicity consistency also separates into two effects:

- LO `sigma0 - 00` in `plain46`:
  - `+0.0001859 +- 0.0000683 pb`
- NLO `sigma0 - 00` in `plain46`:
  - `-0.0000165 +- 0.0000596 pb`

So the finite-width propagator solves the POLDIS-facing unpolarized `Z`
normalization issue, while the small broad-window LO `sigma0 - 00` effect
remains a separate question.

### `ALL`

The full `ALL` channel moves only mildly, as expected for a change that mainly
affects the small `Z` contribution:

- LO unpolarized `sigma0`:
  - `plain46`: `2965.2112 +- 0.0707 pb`
  - POLDIS: `2965.621 +- 0.201 pb`
  - delta: `-0.410 +- 0.213 pb` (`1.93 sigma`)
- NLO unpolarized `sigma0`:
  - `plain46`: `2768.9284 +- 0.0711 pb`
  - POLDIS: `2768.302 +- 0.227 pb`
  - delta: `+0.626 +- 0.238 pb` (`2.63 sigma`)
- LO polarized `sigma_LL`:
  - `plain46`: `50.2263 +- 0.0707 pb`
  - POLDIS: `50.1408 +- 0.0025 pb`
  - delta: `+0.0855 +- 0.0708 pb` (`1.21 sigma`)
- NLO polarized `sigma_LL`:
  - `plain46`: `48.1671 +- 0.0711 pb`
  - POLDIS: `48.2888 +- 0.0025 pb`
  - delta: `-0.1217 +- 0.0711 pb` (`1.71 sigma`)

### Interpretation after `plain46`

`plain46` resolves the outstanding question from the focused `Z` studies:

- the spacelike-`Z` propagator convention is the dominant source of the
  unpolarized `Z` mismatch against POLDIS-like references
- this is true not only in the focused interior-window study, but also in the
  full broad campaign results at both LO and NLO
- the remaining `Z` polarized NLO offset is therefore a different issue
- the small broad-window LO `sigma0 - 00` effect also remains logically
  separate from the propagator-convention fix

## Campaign Outcome: `plain46` -> `plain47` (4x statistics)

`plain47` repeated the finite-width spacelike-`Z` production setup with
substantially higher statistics. The main role of this campaign was to test
whether the `plain46` conclusions were statistically stable.

They are.

### `Z`

The finite-width spacelike-`Z` conclusion is now very solid.

Unpolarized `Z` remains in excellent agreement with the POLDIS-like reference:

- LO unpolarized `sigma0`:
  - `plain47`: `1.2103521 +- 0.0000165 pb`
  - POLDIS: `1.210471 +- 0.000119 pb`
  - delta: `-0.000119 +- 0.000120 pb` (`0.99 sigma`)
- NLO unpolarized `sigma0`:
  - `plain47`: `1.1535935 +- 0.0000146 pb`
  - POLDIS: `1.153555 +- 0.000118 pb`
  - delta: `+0.000039 +- 0.000119 pb` (`0.32 sigma`)

LO polarized `Z` is also now clearly fine:

- LO polarized `sigma_LL`:
  - `plain47`: `0.0546674 +- 0.0000165 pb`
  - POLDIS: `0.0546560 +- 0.0000100 pb`
  - delta: `+0.0000114 +- 0.0000193 pb` (`0.59 sigma`)

The earlier broad-window LO `sigma0 - 00` concern weakens substantially:

- `plain46` LO `sigma0 - 00`: `+0.0001859 +- 0.0000683 pb`
- `plain47` LO `sigma0 - 00`: `+0.0000396 +- 0.0000327 pb`

So with more statistics there is no longer strong evidence for a broad-window
LO `Z` `sigma0 - 00` problem.

The remaining `Z` issue is now much sharper:

- NLO polarized `sigma_LL`:
  - `plain47`: `0.0499144 +- 0.0000146 pb`
  - POLDIS: `0.0501500 +- 0.0000100 pb`
  - delta: `-0.0002356 +- 0.0000177 pb` (`13.29 sigma`)

### `GAMMA` and `ALL`

The higher-statistics `plain47` campaign makes it easier to see that the
remaining residual is not really `Z`-specific.

For `GAMMA`:

- LO polarized `(PP-PM)/2`:
  - `plain47`: `45.3933 +- 0.0500 pb`
  - POLDIS: `45.42693 +- 0.00280 pb`
  - delta: `-0.0336 +- 0.0501 pb` (`0.67 sigma`)
- NLO polarized `(PP-PM)/2`:
  - `plain47`: `43.4797 +- 0.0464 pb`
  - POLDIS: `43.66951 +- 0.00300 pb`
  - delta: `-0.1898 +- 0.0465 pb` (`4.08 sigma`)

For `ALL`:

- LO polarized `sigma_LL`:
  - `plain47`: `50.1571 +- 0.0353 pb`
  - POLDIS: `50.14071 +- 0.00308 pb`
  - delta: `+0.0164 +- 0.0354 pb` (`0.46 sigma`)
- NLO polarized `sigma_LL`:
  - `plain47`: `48.0844 +- 0.0330 pb`
  - POLDIS: `48.28881 +- 0.00323 pb`
  - delta: `-0.2045 +- 0.0332 pb` (`6.17 sigma`)

The relative ratios are strikingly similar:

- `GAMMA` NLO polarized ratio: `0.995654`
- `Z` NLO polarized ratio: `0.995303`
- `ALL` NLO polarized ratio: `0.995766`

The LO polarized sector is therefore fully consistent with closure at the
current precision, while the NLO polarized sector shows a common low bias.

At the full-cross-section level, the three polarized NLO channels are
consistent with a common fractional deficit of about `0.46%`:

- `GAMMA`: `0.435 +- 0.106 %`
- `Z`: `0.470 +- 0.035 %`
- `ALL`: `0.423 +- 0.069 %`

The weighted mean is about `0.458%`, with a good common-fit quality.

This is a stronger and more accurate summary than saying that the NLO
corrections `delta sigma = sigma_NLO - sigma_LO` are all off by the same
percentage. That `delta sigma` statement is not supported equally well across
the three channels, because their NLO K-factors differ substantially, with `Z`
having a much larger relative NLO correction than `GAMMA` or `ALL`.

So after `plain47`, the remaining live residual looks much more like a common
polarized NLO `sigma_NLO`-level low bias than a `Z`-only problem.

### Interpretation after `plain47`

`plain47` sharpens the overall picture considerably:

- the native Born-window fix remains validated
- the finite-width spacelike-`Z` propagator fix remains validated
- the old LO `Z` normalization concern is no longer live
- the old LO broad-window `Z` `sigma0 - 00` concern weakens strongly with more
  statistics
- the dominant remaining discrepancy is now a common polarized NLO low bias
  across `GAMMA`, `Z`, and `ALL`
- the best compact description of that pattern is a common polarized
  `sigma_NLO` fractional deficit of about `0.46%`, not a common fractional
  error on `delta sigma`

This strongly suggests that the next debugging target should be the shared
polarized NLO machinery in `DISBase`, rather than a `Z`-specific propagator or
LO helicity-combination issue.

The leading physics-level candidate is now a polarized NLO coefficient-function
or scheme mismatch that affects the polarized NLO path uniformly while leaving
LO and unpolarized NLO largely intact. That interpretation is consistent with
the `plain47` pattern and with the old POLDIS comments about scheme choices in
the polarized splitting/coefficient functions.

However, the Herwig DIS source does not currently expose an obvious explicit
`Larin`, `HVBM`, `gamma5`, or finite-renormalization switch in the matrix
element code, so this should still be treated as a leading hypothesis rather
than a demonstrated source-level cause. The main alternative remains that the
common effect is generated somewhere in the shared polarized NLO machinery of
`DISBase`, for example in the `xp_`-mapped polarized quark/gluon terms or the
polarized PDF-ratio construction.

## Focused `GAMMA` NLO `SamplingPower` Studies: `plain48` and `plain49`

The first focused `GAMMA` NLO follow-up was the broad-window
`SamplingPower` scan in `plain48`. With the available scan statistics, the
apparent `SamplingPower` dependence looked suggestive:

- `0.2`: `43.4497 +- 0.1001 pb`
- `0.4`: `43.6502 +- 0.1581 pb`
- `0.6`: `43.4797 +- 0.0464 pb`
- `0.8`: `43.4010 +- 0.1581 pb`

At that stage, the most suspicious feature was not the total shift itself but
the fact that the small `NEGNLO` polarized contribution moved noticeably across
the scan, suggesting sensitivity to the `xp_`-dependent subtraction balance.

To reduce the noise of the earlier single-shard TERMDIAG comparison, a
dedicated multi-shard TERMDIAG study was then run in `plain49` with:

- `SamplingPower = 0.4`
- `SamplingPower = 0.6`
- `PP/PM`
- `POSNLO/NEGNLO`
- `10 x 10^6` events per logical run

After fixing the TERMDIAG tag selector so that prefixed variants such as
`TERMDIAG-SP040-S<seed>-plain49-...` were resolved correctly, the cleaned
multi-shard result is:

- `0.4`: `NLO_pol = 43.3511 +- 0.2236 pb`
- `0.6`: `NLO_pol = 43.30075 +- 0.22361 pb`

So the residual difference between `0.4` and `0.6` is only about:

- `-0.05035 +- 0.31624 pb`

which is statistically negligible.

The closure diagnostics in the same TERMDIAG study are exact for both powers:

- `B_closure = 1`
- `M_closure = 1`
- `ME_self = 1`

So the local `GAMMA` POSNLO/NEGNLO algebra is still behaving as designed.

What does move with `SamplingPower` is the internal term decomposition. In the
cleaned `plain49` comparison, the NLO polarized component tables differ by
`O(0.1-1 pb)` in several places, for example:

- `F_virt`: `37.8791 -> 36.3066 pb`
- `F_rq`: `3.4220 -> 4.1403 pb`
- `F_rq_even`: `-5.8943 -> -4.7697 pb`
- `F_cg_odd`: `0.3184 -> 0.6824 pb`

but those shifts mostly cancel in the final NLO polarized total.

So the correct interpretation after `plain49` is:

- `SamplingPower` does affect how the polarized NLO subtraction pieces are
  sampled and partitioned statistically
- but at the current `10 x 10^6` level there is no evidence that changing
  `SamplingPower` from `0.6` to `0.4` produces a meaningful shift in the
  physical broad-window `GAMMA` polarized NLO cross section
- therefore a simple retuning of `SamplingPower` is **not** supported as the
  resolution of the remaining polarized NLO discrepancy

This weakens the earlier idea that one could simply choose a better
`SamplingPower` value to restore agreement with POLDIS. The better summary is:

- the total polarized NLO result is reasonably stable between `0.4` and `0.6`
- the decomposition into virtual/collinear/real pieces is not especially
  stable under that change
- so `SamplingPower` is probing the numerical organization of the subtraction
  terms more than it is acting as a clean physics-fixing knob

## Focused LO `Z` Studies: `plain46-zlo-int` and `plain47-zlo-int-fw`

Two focused LO `Z` studies were run to disentangle:

- the broad-window `sigma0 - 00` behavior
- the spacelike-`Z` propagator convention difference between POLDIS-like and
  Herwig-like calculations

### Interior-window `Z` check without finite-width toggle

The first focused study used the interior window
`[100,1000] x [0.3,0.5]` with the legacy helicity-path spacelike `Z`
propagator and compared it to the broad `plain45` control.

Broad control (`plain45`):

- `00 = 1.2109510 +- 0.0000624 pb`
- `sigma0 = 1.2111630 +- 0.0000280 pb`
- `sigma0 - 00 = +0.0002120 +- 0.0000684 pb`

Interior run (`plain46-zlo-int`):

- `00 = 0.3515960 +- 0.0000200 pb`
- `sigma0 = 0.35157475 +- 0.00000791 pb`
- `sigma0 - 00 = -0.00002125 +- 0.00002151 pb`

Standalone interior references:

- POLDIS-like unpolarized: `0.351364695273 pb`
- Herwig-like unpolarized: `0.351600889855 pb`

Key comparisons:

- `00 - herwig_like = -0.00000489 +- 0.00002000 pb`
- `sigma0 - herwig_like = -0.00002614 +- 0.00000791 pb`
- `sigma0 - poldis_like = +0.00021005 +- 0.00000791 pb`

So in the interior window:

- `sigma0 - 00` is consistent with zero
- `00` and `sigma0` both sit close to the Herwig-like standalone
- the broad-window `Z` LO `sigma0 - 00` mismatch does not survive as a
  generic interior-window effect

This means the broad-window `Z` issue is not simply a universal
`sigma0`/helicity-combination bug in the LO `Z` matrix element.

### Interior-window `Z` check with finite-width spacelike `Z`

The second focused study repeated the same interior-window setup, but turned on
the finite-width spacelike `Z` option in the helicity-amplitude path.

Interior run with finite-width propagator (`plain47-zlo-int-fw`):

- `00 = 0.3513650 +- 0.0000200 pb`
- `sigma0 = 0.35133875 +- 0.00000791 pb`
- `sigma0 - 00 = -0.00002625 +- 0.00002151 pb`

Standalone interior references:

- POLDIS-like unpolarized: `0.351364695273 pb`
- Herwig-like unpolarized: `0.351600889855 pb`

Key comparisons:

- `00 - herwig_like = -0.00023589 +- 0.00002000 pb`
- `sigma0 - herwig_like = -0.00026214 +- 0.00000791 pb`
- `sigma0 - poldis_like = -0.00002595 +- 0.00000791 pb`

This is exactly the expected pattern for a propagator-convention shift:

- both `00` and `sigma0` move downward by about `2.3e-4 pb`
- the shift is approximately the Herwig-like vs POLDIS-like spacelike-`Z`
  propagator gap in that window
- `sigma0 - 00` remains tiny and consistent with zero

So the focused finite-width test supports the interpretation that:

- the spacelike-`Z` propagator convention is the dominant source of the LO `Z`
  normalization difference against the POLDIS-like standalone
- but it is not the source of the broad-window `sigma0 - 00` effect

### Current status after the focused `Z` studies

The focused studies sharpen the diagnosis:

- the broad-window `Z` LO `sigma0 - 00` mismatch is not a generic interior
  helicity-combination failure
- the LO `Z` normalization difference against POLDIS-like references is largely
  a spacelike-`Z` propagator-convention issue
- the broad finite-width `Z` run is still needed to see whether turning on the
  finite-width propagator in the full broad window shifts `sigma0` onto the
  POLDIS-like normalization while leaving the separate broad-window `00`
  behavior unchanged

## NLO `alpha_s` Convention Mismatch for Polarized `\Delta\sigma`

There is a separate, smaller convention mismatch in the polarized NLO setup
that should be documented explicitly.

In POLDIS, the polarized NLO path uses the PDF-provided running coupling from
the active polarized PDF set:

- [`structure.f`](/Users/apapaefs/Projects/HerwigPol/POLDIS/POLDIS-public/structure.f#L277):
  `ALPHAS = alphasPDF(MUR)`
- [`poldis.f`](/Users/apapaefs/Projects/HerwigPol/POLDIS/POLDIS-public/poldis.f#L8162):
  `ALPHAS = alphasPDF(MUR)`

So for direct polarized `\Delta\sigma` calculations, POLDIS is using the
polarized PDF set's own `alpha_s`.

On the Herwig side, the fixed-order DIS NLO weight uses the common model-side
running coupling:

- [`DISBase.cc`](/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/DISBase.cc#L2636):
  `double aS = SM().alphaS(mu2);`

In the DISPOL cards this is normalized through the unpolarized Herwig-side
setup with `alpha_s(M_Z) = 0.118`.

The script
[`compare_pdf_alphas.py`](/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/compare_pdf_alphas.py)
was used to quantify the size of this mismatch at `M_Z`:

- `BDSSV24-NNLO` (polarized PDF): `alpha_s(M_Z) = 0.11916294`
- `PDF4LHC15_nnlo_100_pdfas` (unpolarized PDF): `alpha_s(M_Z) = 0.11800230`
- Herwig Matchbox `NLOAlphaS` with input `alpha_s(M_Z)=0.118`:
  `0.11800008`

So at `M_Z` the Herwig-side `alpha_s` used in the DIS cards is lower than the
polarized BDSSV24 value by about:

- `-0.00116286` in absolute terms
- approximately `-0.98%` relative to `BDSSV24-NNLO`

This matters only for the NLO pieces:

- it cannot explain any LO discrepancy
- it cannot explain the original broad-window LO `GAMMA` pathology
- it does not apply to genuine unpolarized `00` observables in the same way,
  since the mismatch is specifically about the polarized `\Delta\sigma`
  convention

It is therefore best classified as a subleading NLO convention effect rather
than the dominant residual bug.

The sign is also informative. In the polarized channels of interest, the NLO
correction is negative. A smaller Herwig-side `alpha_s` would make that
negative correction smaller in magnitude and therefore push the Herwig NLO
result upward. The remaining polarized NLO residuals after `plain44` are
generally low rather than high, so this `alpha_s` mismatch cannot be the main
explanation for their sign.

So the practical interpretation is:

- this `alpha_s` source mismatch is real and should be kept in mind for the
  remaining polarized NLO comparisons
- but it is too small and has the wrong sign to explain the dominant residual
  polarized NLO differences by itself

## Broad `GAMMA` Raw-Delta Branch Test: `plain53`

After the conservative polarized-NLO refactor was validated point-by-point,
the uniform effective-delta path was promoted to the default production path
(`Branch A`), while the raw finite-delta path was kept as an explicit
experimental option (`Branch B`).

The key implementation-level result is:

- `Branch A` reproduces the legacy polarized odd-sector assembly exactly at
  the TERMDIAG level
- `Branch B` changes the odd-sector kernels locally, as intended, by replacing
  the effective clamped combinations
  - `deltaqOverLo_eff = Pq * dqRatio`
  - `deltagOverLo_eff = Pq * dgRatio`
  with the raw finite-form ratios
  - `deltaqOverLo_raw = Pz * dqPDF / loPDF`
  - `deltagOverLo_raw = Pz * dgPDF / loPDF`

The broad `GAMMA` totals from `plain53` show that this raw-delta change has a
real but modest integrated effect.

### Hybrid PDF setup

- Herwig `Branch A`: `43.4190 +- 0.0707 pb`
- Herwig `Branch B`: `43.5306 +- 0.0707 pb`
- shift `B - A`: `+0.1116 +- 0.1000 pb`
- POLDIS reference: `43.661055 +- 0.009493 pb`

So in the documented hybrid setup:

- `Branch A - POLDIS = -0.2421 +- 0.0713 pb`
- `Branch B - POLDIS = -0.1305 +- 0.0713 pb`

This means the raw-delta option moves Herwig in the right direction for
`GAMMA`, but only by about `0.11 pb`.

### `NNPDF40 + NNPDFpol20` profile

- Herwig `Branch A`: `48.2603 +- 0.0707 pb`
- Herwig `Branch B`: `48.3113 +- 0.0707 pb`
- shift `B - A`: `+0.0510 +- 0.1000 pb`
- POLDIS reference: `48.622141 +- 0.010599 pb`

So in the cleaner `NNPDF40_nlo_pch_as_01180 + NNPDFpol20_nlo_as_01180`
profile:

- `Branch A - POLDIS = -0.3618 +- 0.0715 pb`
- `Branch B - POLDIS = -0.3108 +- 0.0715 pb`

Here the raw-delta option still helps, but much less.

### Practical conclusion

The `plain53` study does not support promoting the raw-delta option to the
production default.

The best summary is:

- `Branch A` should remain the default production path
- `Branch B` is a useful experimental probe of the odd-sector clamp
  sensitivity
- the clamp/effective-delta choice is part of the residual polarized NLO
  story, but it is not the whole story
- the size of the effect is profile-dependent, which weakens the case that it
  is a single universal physics fix

## Plausible Residual Systematics After `plain53`

From the cumulative diagnostics so far, the most plausible residual systematic
effects in the polarized broad-window NLO comparison are:

- the polarized and unpolarized PDFs are not fitted simultaneously
  - in the documented hybrid setup this is explicit:
    `PDF4LHC15_nnlo_100_pdfas` for the denominator and `BDSSV24-NNLO` for the
    polarized difference
  - even the `NNPDF40_nlo_pch_as_01180 + NNPDFpol20_nlo_as_01180` study is a
    cleaner family choice, not a truly joint simultaneous fit
  - so raw `Delta f / f` should not be over-interpreted as a strictly bounded
    event-level polarization variable

- the polarized NLO `alpha_s` convention differs between Herwig and POLDIS
  - POLDIS uses the active PDF-set `alpha_s`
  - Herwig uses the common model-side running coupling
  - this is subleading and has the wrong sign to explain the main deficit by
    itself, but it is still a real convention-level systematic

- the effective-vs-raw odd-sector representation is itself a small systematic
  at the current comparison level
  - `plain53` shows a broad `GAMMA` shift of roughly
    - `+0.11 pb` in the hybrid setup
    - `+0.05 pb` in the `NNPDF40 + NNPDFpol20` setup
  - so the treatment of clamped versus raw polarized deltas is not negligible
    compared to the remaining residual

- the mapped-spin side still mixes a physical polarization interpretation with
  fixed-order PDF ratios
  - `Pq`, `Pq_m`, and `Pg_m` are clamped before they feed the rho-matrix /
    `A_pol()` machinery
  - the raw finite-form odd kernels instead want the unclamped combinations
    `Pz Delta f / f`
  - `plain53` shows that this distinction has a measurable integrated effect,
    even if it is not by itself the full explanation

- PDF-profile dependence itself is a real comparison systematic
  - changing from the hybrid setup to the `NNPDF40 + NNPDFpol20` profile moves
    both the Herwig and POLDIS polarized `GAMMA` normalizations substantially
  - so any residual at the few-`1e-1 pb` level should be interpreted in the
    context of the chosen PDF profile, not as a profile-independent number

What looks less plausible now as a dominant residual systematic:

- the broad-window cut surface itself
  - the focused interior-window `GAMMA` comparison showed that the polarized
    NLO deficit survives away from the broad boundary surfaces

- `SamplingPower`
  - the multi-shard `plain49` study showed that it changes the term
    decomposition more than the final `NLO_pol` total
  - so it is not supported as the main fix

- a simple beam-PDF versus extractor-PDF object mismatch
  - the live audit did not support that hypothesis
