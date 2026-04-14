# BGF Odd-Sector Timeline Note

This note summarizes how the polarized-DIS NLO diagnosis evolved in the
repository's markdown notes, with a focus on the BGF odd sector and the
finite `realg` remainder.

## Why this note exists

By April 14, 2026, the current code-vs-note reading suggests a sharper
hypothesis than the older markdown history:

- the repository had already spent significant time auditing the shared
  polarized-NLO odd sector and the BGF real-emission term;
- but the earlier notes mostly treated this as a broad "odd-sector /
  projector / residual polarized-NLO machinery" problem;
- none of the earlier markdown notes found in this pass appears to state the
  narrower conclusion that the BGF polarized real odd coefficient in the
  single-channel `NLOWeight()` implementation should be read as a post-split
  `-(2 x_p - 1)` object rather than a pre-split `-2 (2 x_p - 1)` one.

So the goal here is to separate what was already established from what looks
newer.

## Chronology

### April 5 to April 8, 2026: the finite BGF remainder was already recognized as special

The strongest earlier summary is in:

- [CODEX_HANDOFF.md](/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/CODEX_HANDOFF.md#L65)

That note records that:

- the finite `realg` term had already been reworked away from the old full
  mapped `A_pol` structure;
- the replacement used the appropriate gluon projectors;
- the change improved `ALL` and interference behavior without disturbing
  `GAMMA`.

This is important context. It means the repository had already converged on
the idea that the BGF finite remainder should be handled in a more structured
gluon-projector basis than the older mapped-`A_pol` implementation.

What it did **not** yet say was that the normalization of the surviving odd
BGF remainder itself was off by a factor of two in the active single-channel
implementation.

### April 6, 2026: the live audit still framed the issue as a common polarized-NLO bias

The next anchor is:

- [polarized_nlo_audit_note.md](/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/polarized_nlo_audit_note.md#L1)

That note shows the audit mindset at the time:

- the persistent problem was a common polarized-NLO bias after `plain36`;
- the live target was still the shared `DISBase::NLOWeight()` machinery;
- the pointwise comparison explicitly included the Herwig BGF mapped real
  piece.

So by April 6, the BGF real term was already on the audit surface, but it was
still being treated as one ingredient inside a broader common-NLO mismatch.

### April 7, 2026: the repository rejected the "missing POLDIS G4NLOg" explanation

The clearest note on this point is:

- [code_audit_April7_mismatch.md](/Users/apapaefs/Projects/HerwigPol/code_audit_April7_mismatch.md#L153)

That note established an important negative result:

- the nonzero Herwig `bgf_split_odd` term should **not** be interpreted as
  evidence that POLDIS is missing an inclusive gluon contribution to `g4`;
- the local Herwig odd BGF term was instead interpreted as a local
  helicity-odd real-emission contribution that can survive pointwise while
  still canceling from the inclusive parity-violating sector after the full
  charge-conjugate sum and integration.

The same note also records where the suspicion sat at that stage:

- the strongest concern remained the even quark-collinear and QCDC pieces in
  the shared DIS NLO machinery;
- the mismatch was already visible in pure `GAMMA`, so the repo was not yet
  centering the diagnosis on a pure-`Z` or weak-only issue.

This note is very relevant because it shows that the repo had already thought
carefully about the meaning of a nonzero odd BGF term, and had **not**
concluded that the problem was simply "POLDIS forgot a gluon term."

However, this April 7 note still did **not** formulate the newer
pre-split/post-split normalization issue.

### April 11, 2026: `realg_odd` was isolated as a first-class diagnostic target

The implementation-planning view is captured in:

- [uniform-polarized-nlo-refactor-plan.md](/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/uniform-polarized-nlo-refactor-plan.md#L176)

That plan explicitly singled out:

- `collg_odd`
- `realg_odd`

as objects to be rewritten in a cleaner canonical form.

But the plan is equally explicit that this phase was supposed to preserve the
current math exactly. The aim was:

- clearer bookkeeping;
- direct expression in `deltagOverLo`;
- old-vs-new diagnostic agreement to floating-point precision.

So by April 11, the repository clearly knew `realg_odd` was important enough
to isolate and diagnose, but the workflow was still about representation
cleanup and protected A/B validation, not about an already-accepted change to
the BGF odd normalization itself.

### April 12, 2026: the odd-sector branch tests showed a real effect, but not a full solution

The most informative summary is in:

- [lo-gamma-cut-surface-note.md](/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/lo-gamma-cut-surface-note.md#L1080)

That note documents the `Branch A` vs `Branch B` odd-sector tests:

- `Branch A` kept the effective-delta odd-sector default;
- `Branch B` used raw finite-form odd kernels;
- the change moved Herwig in the right direction for broad `GAMMA`, including
  in the cleaner NNPDF profile;
- but the integrated shift was modest and profile-dependent.

The practical conclusion in that note is important:

- odd-sector changes are part of the residual story;
- but they are not the whole story;
- the result weakens the case for a single universal odd-sector fix.

This is probably the closest earlier markdown statement to the current debate.
It shows that the repository had already seen that the odd-sector structure
matters numerically, yet still judged the observed effect too small and too
profile-dependent to close the full polarized-NLO gap by itself.

### April 14, 2026: a newer, sharper normalization reading emerged

The newer reading is not something I found stated explicitly in an older
markdown note. It comes from comparing:

- the split derivation in [HerwigPolCodex.tex](/Users/apapaefs/Projects/HerwigPol/HerwigPolCodex.tex#L648),
- the BGF finite-remainder discussion in
  [HerwigPolCodex.tex](/Users/apapaefs/Projects/HerwigPol/HerwigPolCodex.tex#L975),
- and the active single-channel `NLOWeight()` assembly in
  [DISBase.cc](/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/DISBase.cc#L3141).

The sharpened point is:

- the codex writes a factor `2` for the **combined pre-split** odd BGF
  remainder;
- but it then says that after the underlying-Born split, each fixed channel
  receives half, so the single-channel kernel is
  `K_{g,real}^{(q)} = -(2 x_p - 1)`;
- the current `NLOWeight()` object is already a fixed underlying-Born
  single-channel implementation.

That is a narrower claim than the older markdown history.

In other words:

- the older notes had already identified the odd BGF sector as live and
  numerically relevant;
- the newer interpretation is that the active implementation may have been
  carrying a pre-split coefficient into a post-split single-channel object.

## Bottom line

The markdown history already tells a coherent story:

1. the odd BGF term is real and physically nontrivial;
2. it should not be confused with a missing inclusive `G4NLOg` in POLDIS;
3. the BGF finite remainder already needed a more careful projector-based
   treatment than the older mapped-`A_pol` form;
4. odd-sector modifications do move the integrated polarized-NLO result, but
   earlier branch tests suggested they were not obviously a complete fix.

What appears newer than the earlier markdown record is the **specific**
normalization hypothesis:

- the factor `2` belongs to the combined pre-split `R_2 + R_3` object;
- the active `NLOWeight()` implementation should instead use the post-split
  single-channel normalization.

That sharper interpretation should therefore be treated as a new refinement of
the older odd-sector audit history, not as a contradiction of it.
