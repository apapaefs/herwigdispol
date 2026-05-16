# Herwig And Standalone Fixed-Order NLO Correspondence

This note records the component-by-component correspondence between the neutral-current
DIS NLO pieces used by Herwig and the standalone Python fixed-order validator. It is
intended to support the validation statement that the standalone calculator uses the
same analytic Herwig ingredients, while making clear where the event carrier is
deliberately different.

Scope: neutral-current DIS with `GAMMA`, `Z`, and `ALL` exchange, longitudinal
helicity states `00`, `PP`, `PM`, `MP`, and `MM`, fixed-order NLO weights, and
`MC_DIS_BREIT:JETINPUT=MEPARTONS` Rivet analysis. The standalone path is Python-only:
it does not call a C++ bridge, does not select POWHEG emissions, and does not reuse
POLDIS or POWHEG-BOX code.

## Correspondence Table

| NLO ingredient | Herwig source | Standalone Python source | Correspondence statement |
| --- | --- | --- | --- |
| DIS run setup: beams, cuts, helicities, PDFs, NC setups | Herwig cards and matrix-element state feeding `DISBase` and `MENeutralCurrentDIS` | `DISPOL/scripts/herwig_fo/config.py`; CLI setup in `DISPOL/scripts/herwig_fixed_order_nlo.py` | Same validation scope is encoded in the Python run config: `18 x 275 GeV`, `100 < Q2 < 2500 GeV^2` for the validation window, `0.2 < y < 0.6`, paired unpolarized/polarized NNPDF profiles, and NC `GAMMA`, `Z`, `ALL`. |
| Scale choice and `alpha_s` | `DISBase::scale()` and `SM().alphaS(mu2)` in `DISBase::NLOWeightRaw()` (`HerwigSource/Herwig-7.3.0/MatrixElement/DIS/DISBase.cc:1859`) | `MatchboxAlphaS` in `DISPOL/scripts/herwig_fo/alphas.py`; `scale_factor` propagation through `DISPoint.mu2` | Python uses the same intended validation convention `mu2 = Q2 * scale_factor^2` and `alpha_s(MZ)=0.118`. This matches the standalone fixed-order/POLDIS comparison convention, not necessarily Herwig's generator-native POWHEG real-emission scale used in shower/no-shower POWHEG campaigns. |
| Born DIS kinematics and normalization | `HwMEBase::dSigHatDR()` and DIS Born state used by `DISBase::dSigHatDR()` (`DISBase.cc:1851`) plus `MENeutralCurrentDIS::helicityME()` (`MENeutralCurrentDIS.cc:366`) | `DISPOL/scripts/herwig_fo/born.py:17`, `born_weight()` at `born.py:70`, `born_event()` at `born.py:91` | Python reconstructs the same Born-projected DIS point, including the `x f(x)` luminosity, DIS Jacobians, helicity density dependence, and HepMC Born callback event. |
| Neutral-current gamma/Z/interference coefficients | `MENeutralCurrentDIS::ncCoefficients()` (`MENeutralCurrentDIS.cc:876`) and `MENeutralCurrentDIS::helicityME()` (`MENeutralCurrentDIS.cc:366`) | `nc_coefficients()` in `DISPOL/scripts/herwig_fo/ew.py:93` | The Python EW module ports the same NC coefficient structure for photon, Z, and interference terms. |
| Born analyzing power `A_pol` | `MENeutralCurrentDIS::A_pol()` (`MENeutralCurrentDIS.cc:761`) | `a_pol()` in `DISPOL/scripts/herwig_fo/ew.py:136` | Same helicity analyzing-power logic is used to build polarized Born and projected NLO factors. |
| Born helicity combinations | `DISBase::generateRivetWeights()` builds `HERWIG_DIS_PP`, `HERWIG_DIS_PM`, `HERWIG_DIS_MP`, `HERWIG_DIS_MM`, `HERWIG_DIS_SIGMA`, and `HERWIG_DIS_DELTA_LL` (`DISBase.cc:2199`) | `correlated_helicity_weights()` in `DISPOL/scripts/herwig_fixed_order_nlo.py:774`; correlated event generation at `herwig_fixed_order_nlo.py:1037` | Both use `sigma = (PP + PM + MP + MM)/4` and `Delta_LL = (PP + MM - PM - MP)/4`. The standalone default is a correlated named-weight stream so numerator and denominator are evaluated at the same phase-space points. |
| PDF ratios and effective parton polarizations | `DISBase::NLOWeightRaw()` PDF block, including `loPDF`, `qPDF`, `gPDF`, `dqPDF`, `dgPDF`, `Pq`, `Pq_m`, and `Pg_m` (`DISBase.cc:1859`) | `ratios_for_flavor()` in `DISPOL/scripts/herwig_fo/pdfs.py:141` | Python uses unpolarized PDFs for flux ratios and polarized difference PDFs only through effective longitudinal parton polarizations and `Delta f/f` ratios, matching the Herwig polarized-NLO representation. |
| `x_p` mapping and Jacobian | Production power-law map feeding `xp_` and `jac_`, used by `NLOWeightRaw()` (`DISBase.cc:1779`, `DISBase.cc:1859`) | `herwig_xp_power_map()` in `DISPOL/scripts/herwig_fo/nlo.py:28` | Same `x_p` power map and Jacobian convention are used when forming the inclusive NLO ratio and the real/counterevent weights. |
| Virtual correction | `virt` term in `DISBase::NLOWeightRaw()` (`DISBase.cc:1859`) | `nlo_terms()` virtual component in `DISPOL/scripts/herwig_fo/nlo.py:40` | Python ports the Herwig virtual-plus-subtracted-Born ratio, including the logarithms and division by the sampling Jacobian. |
| Quark collinear counterterm | `collq_k1`, `collq_k2`, `collq_even`, and `collq_odd` in `DISBase::NLOWeightRaw()` (`DISBase.cc:1859`) | `nlo_terms()` quark collinear component in `DISPOL/scripts/herwig_fo/nlo.py:112` | Same quark PDF-ratio structure, same `Pqq` kernel pieces, and same polarized odd/even NC response split. |
| Gluon collinear counterterm | `collg_over_born_unpol`, `collg_even`, and `collg_odd` in `DISBase::NLOWeightRaw()` (`DISBase.cc:1859`) | `nlo_terms()` gluon collinear component in `DISPOL/scripts/herwig_fo/nlo.py:98` | Same gluon PDF-ratio structure and polarized gluon response are used. |
| Collinear projector and NC response factors | `MENeutralCurrentDIS::collinearBlendWeights()` (`MENeutralCurrentDIS.cc:783`) and `MENeutralCurrentDIS::neutralCurrentResponse()` (`MENeutralCurrentDIS.cc:943`) | `collinear_blend_weights()` in `DISPOL/scripts/herwig_fo/ew.py:174`; `neutral_current_response()` in `ew.py:243` | The Python code uses the same Herwig projector logic separating unpolarized/even and polarized/odd NC response pieces. |
| Integrated QCDC finite real piece in `NLOWeightRaw` | `realq_prefactor`, `realq_even`, and `realq_odd` in `DISBase::NLOWeightRaw()` (`DISBase.cc:1859`) | `nlo_terms()` QCDC integrated component in `DISPOL/scripts/herwig_fo/nlo.py:142` | The inclusive NLO seed uses the same Herwig integrated QCDC finite term. The standalone real callback then adds and subtracts local QCDC events to redistribute this contribution over Rivet observables. |
| Integrated BGF finite real piece in `NLOWeightRaw` | `realg_over_born_unpol`, `realg_even`, and `realg_odd` in `DISBase::NLOWeightRaw()` (`DISBase.cc:1859`) | `nlo_terms()` BGF integrated component in `DISPOL/scripts/herwig_fo/nlo.py:147` | The inclusive NLO seed uses the same Herwig integrated BGF finite term. The standalone local BGF real/counterevent pair is used for differential fixed-order fills. |
| QCDC mapped denominator and real-emission denominator | `MENeutralCurrentDIS::qcdcMappedDenominatorRatio()` (`MENeutralCurrentDIS.cc:823`) and `realEmissionDenominatorFactor()` (`MENeutralCurrentDIS.cc:840`) | `qcdc_mapped_denominator_ratio()` in `DISPOL/scripts/herwig_fo/ew.py:205`; `real_emission_denominator_factor()` in `ew.py:224` | Same mapped polarized denominator factors are used for the local QCDC real kernel. This is the piece that fixed the polarized `DpT1PreCut` depletion in the standalone run. |
| QCDC local real-emission matrix element | `DISBase::ComptonME()` (`DISBase.cc:1393`) | `compton_azimuthal_coefficients()` in `DISPOL/scripts/herwig_fo/real.py:142`; `qcdc_real_weight_ratio()` in `real.py:251` | The standalone uses the same analytic `c0`, `c1`, `c2` azimuthal kernel, including mapped leading term, mapped analyzing power, and denominator ratio. Unlike Herwig's POWHEG/MEC event carrier, it writes this as a signed fixed-order real event plus a Born-projected counterevent. |
| BGF local real-emission matrix element | `DISBase::BGFME()` (`DISBase.cc:1453`) | `bgf_azimuthal_coefficients()` in `DISPOL/scripts/herwig_fo/real.py:184`; `bgf_real_weight_ratio()` in `real.py:264` | The standalone uses the same analytic BGF azimuthal kernel and the same NC spin projectors. Its fixed-order callback implements the endpoint partition needed for POLDIS-style differential closure. |
| Real-emission Breit-frame momenta | `DISBase::generateCompton()` and BGF generation flow (`DISBase.cc:1552`) | `_breit_two_body_momenta()` in `DISPOL/scripts/herwig_fo/real.py:69`; `real_event()` in `real.py:350` | Python builds the two-parton QCDC/BGF final state in the Breit frame and maps it to the lab with a DIS tetrad. The event record is closed with the mapped incoming parton remnant, so `MC_DIS_BREIT:JETINPUT=MEPARTONS` sees the hard partons directly. |
| Inclusive seed plus local real/counterevent fixed-order callback | Herwig does not expose this as its ordinary POWHEG carrier. `NLOWeightRaw()` supplies the inclusive ratio, while `ComptonME()`/`BGFME()` supply local real kernels. | `generate_events()` in `DISPOL/scripts/herwig_fixed_order_nlo.py:970`; `generate_correlated_events()` in `herwig_fixed_order_nlo.py:1037` | This is the main deliberate architectural difference. The standalone constructs fixed-order Rivet callbacks: Born-projected inclusive seed, local real event, and Born-projected counterevent. This is not the POWHEG-selected event path. |
| POS/NEG handling | `DISBase::NLOWeight()` clips `NLOWeightRaw()` into positive and negative contribution streams (`DISBase.cc:2151`) | Standalone HepMC events carry signed weights directly through `write_hepmc3()` in `DISPOL/scripts/herwig_fo/hepmc.py:65` | The physics density is the same signed NLO density, but the bookkeeping differs. Herwig campaign production normally separates POSNLO/NEGNLO streams; standalone writes signed callback events in one fixed-order stream. |
| Named HepMC/Rivet weights | `DISBase::generateRivetWeights()` (`DISBase.cc:2199`) and ThePEG HepMC named-weight propagation | `write_hepmc3()` in `DISPOL/scripts/herwig_fo/hepmc.py:65`; correlated bundles in `herwig_fixed_order_nlo.py:1037` | The standalone writes `HERWIG_DIS_*` named weights in HepMC3 and runs Rivet with `RIVETWEIGHTS=YES`, matching the Herwig correlated-helicity analysis contract. |
| Rivet analysis mode | Herwig campaign runs `MC_DIS_BREIT:JETINPUT=MEPARTONS:RIVETWEIGHTS=YES` for correlated no-shower validation | `run_rivet()` in `DISPOL/scripts/herwig_fo/rivet.py:24` | Same Rivet analysis and ME-parton input contract are used for the fixed-order comparison plots. |
| Scale variations | Herwig `ScaleFactor` variations in DIS validation campaigns; current POWHEG campaigns use Herwig-native real-emission central scales | `--scale-factor` and `--scale-variations` in `herwig_fixed_order_nlo.py`, with `mu2 = Q2 * scale_factor^2` | Standalone scale variations are fixed-order `Q2`-central variations. They are appropriate for POLDIS/Python fixed-order envelopes, but should not be described as identical to Herwig's native POWHEG emission-scale bands unless the Herwig campaign is explicitly configured to the same central scale. |

## Important Non-Identical Pieces

The standalone validator is not a clone of the Herwig event generator path. It is a
fixed-order callback calculator built from Herwig's analytic components.

- The Herwig POWHEG/MEC path selects real emissions through generation kernels,
  positive/negative contribution streams, veto/Sudakov machinery, and event-record
  construction. The standalone path writes signed fixed-order callback events:
  inclusive Born-projected seed, local QCDC/BGF real event, and Born-projected
  counterevent.
- The local QCDC and BGF matrix-element kernels are the same Herwig formulas, but
  the standalone uses them as physical fixed-order weights for Rivet fills rather
  than as POWHEG/MEC emission-proposal kernels.
- The standalone scale-variation convention is `mu2 = Q2 * scale_factor^2`.
  Herwig POWHEG validation campaigns may instead use generator-native emission
  scales for the real-emission alpha_s/PDF factors, depending on the card settings.
- No parton shower, hadronization, decays, Sudakov suppression, or POWHEG-selected
  hardest-emission path is included in the standalone validation.

Within those boundaries, the standalone calculation is a Herwig-component
fixed-order validator: the Born coefficients, PDF ratios, helicity projectors,
virtual term, collinear terms, integrated finite real pieces, and local QCDC/BGF
real-emission kernels are matched component by component to the active Herwig DIS
implementation.
