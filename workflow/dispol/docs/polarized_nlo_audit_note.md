# Polarized NLO Audit Note

This note records the live-code audit targets for the persistent common
polarized NLO bias seen after `plain36`.

## Current runtime targets

- Live Herwig NLO path:
  - `DISBase::NLOWeight()` in
    `/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/DIS/DISBase.cc`
  - `HwMEBase::getBeamPolarization()` in
    `/Users/apapaefs/Projects/HerwigPol/HerwigSource/Herwig-7.3.0/MatrixElement/HwMEBase.cc`
  - `PolarizedPartonExtractor::getRhoMatrices()` in
    `/Users/apapaefs/Projects/HerwigPol/HerwigSource/ThePEG-2.3.0/PDF/PolarizedPartonExtractor.cc`

- New runtime diagnostics:
  - `NLO_AUDIT_OBJ`
  - `NLO_AUDIT_PDF`
  - `NLO_AUDIT_TERM`

These are switch-gated in `DISBase` and only sample accepted POSNLO/NEGNLO
events.

## Main plumbing question

The live code still has a potentially dangerous split between the objects used
for polarized and unpolarized hadron-side information:

- `HwMEBase::getBeamPolarization()` follows `lastExtractor()` and the polarized
  extractor path.
- `DISBase::NLOWeight()` rebuilds the polarized hadron-side ratios from
  `lastExtractor()->longitudinalDifferencePDF().second`.
- But the unpolarized denominator there is still `hadron_->pdf()`.

The new runtime audit adds the extractor-side hadron PDF from
`lastXComb().partonBinInstances().second->pdf()` so the log can tell us whether
the live event is using the same unpolarized PDF object in both places.

## Current Tiresias `.run` state

Representative plain-analysis runs on Tiresias still serialize the proton beam
and polarized extractor through `HardLOPDF`:

- `/Users/apapaefs/trs/Projects/HerwigPol/HwPolNotesNew/DISPOL/DIS-POL-POWHEG_PP-POSNLO-ALL.run`
- `/Users/apapaefs/trs/Projects/HerwigPol/HwPolNotesNew/DISPOL/DIS-POL-POWHEG_PP-NEGNLO-ALL.run`

The same runs also serialize:

- `QCD/RunningAlphaS -> /Herwig/Couplings/NLOAlphaS`
- `HardLOPDF` object entries carrying `PDF4LHC15_nnlo_100_pdfas`

So the live serialized state appears to be:

- routing still through `HardLOPDF`
- but the `HardLOPDF` object itself now points at the intended
  `PDF4LHC15_nnlo_100_pdfas`

The new runtime audit is meant to verify whether the event-time beam object,
the current `XComb` hadron bin, and the polarized extractor are all consistent
with that serialized picture.

## Standalone kernel comparison

The companion script

- `/Users/apapaefs/Projects/HerwigPol/HwPolNotesNew/DISPOL/compare_polarized_nlo_terms.py`

compares:

- Herwig Born polarization inputs `Pq`, `Pq_m`, `Pg_m`, `dqRatio`, `dgRatio`
- Herwig quark singular-like piece `(virtual correction + quark collinear)`
- Herwig gluon collinear piece
- Herwig `QCDC` mapped real piece
- Herwig `BGF` mapped real piece

against the closest flavor-resolved polarized POLDIS inclusive pieces built
from the `G2/G4/GL` coefficient functions.

This comparison is deliberately pointwise and pre-integration. It is not a
strict local-subtraction identity check, because POLDIS exposes inclusive
coefficient functions rather than Herwig’s separate virtual/collinear/real
subtraction terms.
