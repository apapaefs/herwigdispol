# PDF-derived spin-density positivity

`PolarizedPartonExtractor` constructs an incoming spin-half density matrix
from longitudinal and transverse difference-PDF ratios. Independent central
polarized and unpolarized fits can yield a finite polarization vector outside
the unit Bloch ball at rare large-x phase-space points. Such a vector is not a
positive semidefinite density matrix and can make a later, otherwise positive,
polarized backward-ISR contraction negative.

The maintained prescription projects the complete vector
`(pL, Re(pT), Im(pT))` radially onto the unit Bloch ball before the matrix is
constructed. It leaves physical vectors unchanged, preserves the direction of
an overshoot, consumes no random number, and rejects nonfinite inputs. It is
not an event or branching veto, and the strict material-negative ISR guard
remains active.

This is an interim density-construction prescription rather than a correction
to the PDFs or a derivation from NLO factorization. The exact STAR audit found
91 projected accepted hard points in 1.5 million events, no projected accepted
entry in the retained inclusive bins at pT >= 13.1 GeV or in any dijet bin,
and no resolved displacement in a 50,000-event same-seed counterfactual. The
full evidence and limitations are recorded in HerwigPol commit
`ad0c5486ad137e0e34c6274959c1fbf8323c5937`.
