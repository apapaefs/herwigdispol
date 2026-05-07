// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISBase class.
//

#include "DISBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Repository/BaseRepository.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/PDF/LHAPDF6.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include <numeric>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <set>
#include "Herwig/Shower/RealEmissionProcess.h"
#include "ThePEG/PDF/PolarizedPartonExtractor.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

namespace {

// Evaluate the compact azimuthal kernel c0 + c1 cos(phi) + c2 cos^2(phi).
double azimuthalKernelValue(double c0, double c1, double c2, double cphi) {
  return c0 + c1 * cphi + c2 * cphi * cphi;
}

// Find a safe rejection-sampling envelope for the quadratic azimuthal kernel.
double azimuthalKernelMaximum(double c0, double c1, double c2) {
  double phimax = std::max(azimuthalKernelValue(c0, c1, c2, -1.0),
                           azimuthalKernelValue(c0, c1, c2,  1.0));
  if (std::abs(c2) > 1e-30) {
    const double cstar = -0.5 * c1 / c2;
    if (cstar > -1.0 && cstar < 1.0) {
      phimax = std::max(phimax, azimuthalKernelValue(c0, c1, c2, cstar));
    }
  }
  return phimax;
}

// Compute the PDF ratio used by POWHEG/MEC emission generation. The signed
// ratio is optional because the veto algorithm must sample with a positive
// probability but some callers still need the physical sign.
bool powhegEmissionPDFRatioForGeneration(tcPDFPtr pdf, tcBeamPtr beam,
                                         tcPDPtr numeratorParton,
                                         tcPDPtr denominatorParton,
                                         Energy2 numeratorScale,
                                         Energy2 denominatorScale,
                                         double numeratorX,
                                         double denominatorX,
                                         double & generationRatio,
                                         double * signedRatio = 0) {
  generationRatio = 0.0;
  if (signedRatio) *signedRatio = 0.0;
  if (!pdf || !beam || !numeratorParton || !denominatorParton) return false;

  const double numerator =
    pdf->xfx(beam, numeratorParton, numeratorScale, numeratorX);
  const double denominator =
    pdf->xfx(beam, denominatorParton, denominatorScale, denominatorX);
  if (!std::isfinite(numerator) || !std::isfinite(denominator) ||
      std::abs(denominator) <= 1e-30) {
    return false;
  }

  const double ratio = numerator / denominator;
  if (!std::isfinite(ratio)) return false;

  if (signedRatio) *signedRatio = ratio;
  // The veto algorithm needs a bona fide probability. Keep the signed ratio
  // available to callers that need it, but generate with its magnitude.
  generationRatio = std::abs(ratio);
  return true;
}

bool powhegAnalyzingPowerIsSafe(double value) {
  // For physical longitudinal spin states the DIS analyzing power should stay
  // inside |A_pol| <= 2. Values outside that range are a strong sign that the
  // denominator reconstruction has become numerically delicate.
  return std::isfinite(value) && std::abs(value) <= 2.0 + 1e-10;
}

bool powhegPositiveRatioIsSafe(double numerator,
                               double denominator,
                               double & ratio) {
  ratio = 1.0;
  if (!std::isfinite(numerator) || !std::isfinite(denominator)) return false;

  // The analytic helicity kernels factor through denominator pieces that are
  // expected to stay positive. When that response turns tiny or changes sign,
  // fall back to the stable Born kernel rather than feeding a signed or wildly
  // amplified ratio into the veto probability.
  const double floor = 1e-12 * std::max(1.0, std::abs(denominator));
  if (std::abs(denominator) <= floor) return false;
  if (numerator <= 0.0) return false;

  const double candidate = numerator / denominator;
  if (!std::isfinite(candidate) || candidate <= 0.0) return false;

  ratio = candidate;
  return true;
}

// Parse scalar cut-interface values returned as strings by ThePEG interfaces.
bool parseDoubleValue(const std::string & text, double & value) {
  std::istringstream is(text);
  is >> value;
  return bool(is);
}

// Extract the active switch option name from a ThePEG switch-interface string.
bool parseSwitchOptionName(const std::string & text, std::string & option) {
  const std::string::size_type begin = text.find('[');
  const std::string::size_type end = text.find(']', begin);
  if (begin == std::string::npos || end == std::string::npos || end <= begin + 1)
    return false;
  option = text.substr(begin + 1, end - begin - 1);
  return true;
}

}


// Set the release defaults for the shared DIS NLO, MEC, POWHEG and
// correlated-weight controls.
DISBase::DISBase()  : initial_(6.), final_(3.),
		      procProb_(0.35),
		      comptonInt_(0.), bgfInt_(0.),
		      comptonWeight_(50.), BGFWeight_(150.), 
		      pTmin_(0.1*GeV), 
		      scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.),
		      contrib_(0),
		      useFixedOrderAlphaSInPOWHEGEmission_(false),
		      useQ2ScaleInPOWHEGEmission_(false),
		      generateRivetWeights_(false),
		      useNativeDISWindowGeneration_(false),
		      powhegEmissionComparisonMode_(POWHEGEmissionComparisonModeDefault),
		      powhegEmissionComparisonMaxAttempts_(100),
		      usePOWHEGRealSpinVertex_(false),
		      leptonPolarization_(0.0),
		      comptonRawXP_(0.0), comptonRawZP_(0.0),
		      bgfRawXP_(0.0), bgfRawZP_(0.0),
		      xpSamplingRandom_(0.0), xpSamplingRho_(0.0),
		      xpSamplingRhomin_(0.0),
		      nativeDISWindow_(),
		      power_(0.1)
{}

DISBase::~DISBase() {}

// Persist only the current release-facing fields. Retired validation-only
// fields were intentionally dropped with class version 14.
void DISBase::persistentOutput(PersistentOStream & os) const {
  os << comptonInt_ << bgfInt_ << procProb_ << initial_ << final_ << alpha_
     << ounit(pTmin_,GeV) << comptonWeight_ << BGFWeight_ << gluon_
     << ounit(muF_,GeV) << scaleFact_ << scaleOpt_ << contrib_
     << useFixedOrderAlphaSInPOWHEGEmission_
     << useQ2ScaleInPOWHEGEmission_
     << useNativeDISWindowGeneration_
     << powhegEmissionComparisonMode_
     << powhegEmissionComparisonMaxAttempts_
     << usePOWHEGRealSpinVertex_
     << generateRivetWeights_ << power_;
}

// Read the current persistent layout. Older .run files are rejected explicitly
// so they cannot be silently misread after the release interface cleanup.
void DISBase::persistentInput(PersistentIStream & is, int version) {
  if(version != 14) {
    throw Exception()
      << "DISBase persistent input version " << version
      << " is no longer supported after the release interface cleanup. "
      << "Regenerate .run files from the source cards."
      << Exception::runerror;
  }

  is >> comptonInt_ >> bgfInt_ >> procProb_  >> initial_ >> final_ >> alpha_
     >> iunit(pTmin_,GeV) >> comptonWeight_ >> BGFWeight_ >> gluon_
     >> iunit(muF_,GeV) >> scaleFact_ >> scaleOpt_ >> contrib_
     >> useFixedOrderAlphaSInPOWHEGEmission_
     >> useQ2ScaleInPOWHEGEmission_
     >> useNativeDISWindowGeneration_
     >> powhegEmissionComparisonMode_
     >> powhegEmissionComparisonMaxAttempts_
     >> usePOWHEGRealSpinVertex_
     >> generateRivetWeights_ >> power_;
  if(powhegEmissionComparisonMode_ >
     POWHEGEmissionComparisonModeRealOnly) {
    powhegEmissionComparisonMode_ = POWHEGEmissionComparisonModeDefault;
  }
  if(powhegEmissionComparisonMaxAttempts_ == 0) {
    powhegEmissionComparisonMaxAttempts_ = 100;
  }
  nativeDISWindow_ = NativeDISWindowDefinition();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<DISBase,HwMEBase>
  describeHerwigDISBase("Herwig::DISBase", "HwMEDIS.so", 14);

// Extract longitudinal polarisation for spin-1/2: P = (rho++ - rho--)/(rho++ + rho--).
// ThePEG convention (PolarizedPartonExtractor): index 0 = NEGATIVE helicity,
//   index 1 (imax) = POSITIVE helicity.
//   rho(0,0) = (1-pL)/2,  rho(1,1) = (1+pL)/2.
double DISBase::longPol(const ThePEG::RhoDMatrix& rho) const {
  const double rhoMM = std::max(0.0, std::real(rho(0,0)));  // index 0 = minus
  const double rhoPP = std::max(0.0, std::real(rho(1,1)));  // index 1 = plus
  const double norm  = (rhoPP + rhoMM) > 0.0 ? (rhoPP + rhoMM) : 1.0;
  return (rhoPP - rhoMM) / norm;
}

// Build a diagonal longitudinal spin-density matrix for the particle type used
// by the DIS spin paths.
RhoDMatrix DISBase::longitudinalRhoMatrix(tcPDPtr data, double pol) const {
  if (!data) return RhoDMatrix();

  pol = std::max(-1.0, std::min(1.0, pol));
  RhoDMatrix rho(data->iSpin());

  if (data->iSpin() == PDT::Spin1Half) {
    const unsigned int imax = rho.iSpin() - 1;
    rho(0,0) = 0.5 * (1.0 - pol);
    rho(imax,imax) = 0.5 * (1.0 + pol);
    rho(0,imax) = 0.0;
    rho(imax,0) = 0.0;
  }
  else if (data->iSpin() == PDT::Spin1) {
    rho(0,0) = 0.5 * (1.0 - pol);
    rho(1,1) = 0.0;
    rho(2,2) = 0.5 * (1.0 + pol);
    rho(0,1) = rho(1,0) = 0.0;
    rho(0,2) = rho(2,0) = 0.0;
    rho(1,2) = rho(2,1) = 0.0;
  }

  return rho;
}

// Convert a beam-level longitudinal polarization and polarized-PDF ratio into
// an event-local incoming-parton rho matrix at the requested mapped x.
DISBase::MappedIncomingSpinDensity
DISBase::mappedIncomingSpinDensity(tcPDPtr parton,
	                                   double x,
	                                   Energy2 mu2,
	                                   const char * clampSource) const {
  MappedIncomingSpinDensity out{
    longitudinalRhoMatrix(parton, 0.0),
    0.0
  };
  (void) clampSource;

  if (!parton) return out;
  if (!hadron_) return out;
  if (x <= 0.0 || x >= 1.0) return out;

  const double Pz = getBeamPolarization(false).z();
  if (std::abs(Pz) <= 1e-12) return out;

  ThePEG::tcPDFPtr sumPdf = hadron_->pdf();
  if (!sumPdf) return out;

  ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr ppe =
    ThePEG::dynamic_ptr_cast<
      ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr>(lastExtractor());
  if (!ppe) return out;

  ThePEG::tcPDFPtr diffPdf = ppe->longitudinalDifferencePDF().second;
  if (!diffPdf) return out;

  const double sum = sumPdf->xfx(hadron_, parton, mu2, x) / x;
  const double diff = diffPdf->xfx(hadron_, parton, mu2, x) / x;

  if (!std::isfinite(sum) || !std::isfinite(diff)) return out;
  if (std::abs(sum) <= 1e-30) return out;

  const double rawPolarization = Pz * diff / sum;
  if (!std::isfinite(rawPolarization)) return out;

  out.clampedPolarization =
    std::max(-1.0, std::min(1.0, rawPolarization));
  out.rho = longitudinalRhoMatrix(parton, out.clampedPolarization);
  return out;
}

// Rebuild the incoming hadron-side rho matrix from the polarized extractor for
// NLO points, where reconstructed parton bins may not retain the correct
// polarized PDF information.
pair<RhoDMatrix,RhoDMatrix> DISBase::correctedLongitudinalRhoMatrices() const {
  pair<RhoDMatrix,RhoDMatrix> rho = getRhoMatrices();

  // Leave LO and unpolarized runs on the generic rho-matrix path.
  // At NLO the reconstructed parton bins can lose the polarized PDF
  // information on the hadron leg, so rebuild that entry from Pz*Delta q/q.
  if (contrib_ == 0 || !hadron_ || xB_ <= 0.0) return rho;

  const double Pz = getBeamPolarization(false).z();
  if (std::abs(Pz) <= 1e-12) return rho;

  ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr ppe =
    ThePEG::dynamic_ptr_cast<
      ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr>(lastExtractor());
  if (!ppe) return rho;

  ThePEG::tcPDFPtr diffPdf = ppe->longitudinalDifferencePDF().second;
  ThePEG::tcPDFPtr sumPdf = hadron_->pdf();
  if (!diffPdf || !sumPdf || !mePartonData()[1]) return rho;

  const Energy2 mu2 = scale();
  const double loPDF = sumPdf->xfx(hadron_, mePartonData()[1], mu2, xB_) / xB_;
  if (std::abs(loPDF) <= 1e-30) return rho;

  const double dloPDF = diffPdf->xfx(hadron_, mePartonData()[1], mu2, xB_) / xB_;
  double Pq = Pz * dloPDF / loPDF;
  const double clampedPq = std::max(-1.0, std::min(1.0, Pq));
  Pq = clampedPq;

  // DIS always has a spin-1/2 incoming quark here.
  rho.second = RhoDMatrix(mePartonData()[1]->iSpin());
  const unsigned int imax = rho.second.iSpin() - 1;
  rho.second(0,0) = 0.5 * (1.0 - Pq);
  rho.second(imax,imax) = 0.5 * (1.0 + Pq);
  rho.second(0,imax) = 0.0;
  rho.second(imax,0) = 0.0;

  return rho;
}

// Default collinear projector: the legacy scalar blend used by photon-like
// channels. Neutral/charged current implementations override where needed.
DISBase::CollinearBlendWeights
DISBase::collinearBlendWeights(tcPDPtr lin, tcPDPtr lout,
                               tcPDPtr qin, tcPDPtr qout,
                               Energy2 scale, double Pl, double Pq,
                               double ell) const {
  const double a = A_pol(lin, lout, qin, qout, scale, Pl, Pq);
  const double denom = 1. + a * ell + sqr(ell);
  const double f = (std::abs(denom) > 1e-30) ? a * ell / denom : 0.0;
  return {1.0 - f, f, 1.0 - f, f};
}

// Default mapped-denominator correction for QCDC real emission.
double DISBase::qcdcMappedDenominatorRatio(tcPDPtr, tcPDPtr,
                                           tcPDPtr, tcPDPtr,
                                           Energy2,
                                           double,
                                           double,
                                           double) const {
  return 1.0;
}

// Default real-emission denominator factor used when no current-specific
// parity decomposition is required.
double DISBase::realEmissionDenominatorFactor(tcPDPtr, tcPDPtr,
                                              tcPDPtr, tcPDPtr,
                                              Energy2,
                                              double,
                                              double) const {
  return 1.0;
}

// The base class keeps real-emission spin kernels on the Born polarization
// unless a current-specific implementation opts into mapped polarizations.
bool DISBase::useMappedPolarizedEmissionKernel() const {
  return false;
}

// Default neutral-current response hook. Returning false tells the common NLO
// assembly to use the photon-like fallback projectors.
bool DISBase::neutralCurrentResponse(tcPDPtr, tcPDPtr,
                                     tcPDPtr, tcPDPtr,
                                     Energy2,
                                     double,
                                     double,
                                     double,
                                     double,
                                     NeutralCurrentResponse &) const {
  return false;
}

// Optional spin-only POWHEG real-emission vertex hook. Current-specific classes
// attach an exact 2->3 HardVertex when the interface switch enables it.
void DISBase::constructRealEmissionSpinVertex(RealEmissionProcessPtr, bool) const {}

// Convenience wrapper for callers that only need the mapped longitudinal
// polarization rather than the full rho matrix.
double DISBase::mappedIncomingLongitudinalPolarization(tcPDPtr parton,
                                                       double x,
                                                       Energy2 mu2) const {
  return mappedIncomingSpinDensity(parton, x, mu2, "powheg-map").clampedPolarization;
}

// Locate the unique SimpleDISCut window that matches the current charged- or
// neutral-current matrix element, and cache it for Born generation.
bool DISBase::resolveNativeDISWindow() const {
  if (nativeDISWindow_.resolutionAttempted) return nativeDISWindow_.available;

  nativeDISWindow_ = NativeDISWindowDefinition();
  nativeDISWindow_.resolutionAttempted = true;

  const Cuts::TwoCutVector & twoCuts = lastCuts().twoCuts();
  NativeDISWindowDefinition matchedWindow;
  unsigned int matchedCount = 0;

  for (Cuts::TwoCutVector::const_iterator it = twoCuts.begin();
       it != twoCuts.end(); ++it) {
    IBPtr cutObject = *it;
    if (!cutObject) continue;

    const InterfaceBase * minQ2If = BaseRepository::FindInterface(cutObject, "MinQ2");
    const InterfaceBase * maxQ2If = BaseRepository::FindInterface(cutObject, "MaxQ2");
    const InterfaceBase * minYIf = BaseRepository::FindInterface(cutObject, "Miny");
    const InterfaceBase * maxYIf = BaseRepository::FindInterface(cutObject, "Maxy");
    const InterfaceBase * minW2If = BaseRepository::FindInterface(cutObject, "MinW2");
    const InterfaceBase * maxW2If = BaseRepository::FindInterface(cutObject, "MaxW2");
    const InterfaceBase * currentIf = BaseRepository::FindInterface(cutObject, "Current");

    if (!minQ2If || !maxQ2If || !minYIf || !maxYIf || !minW2If || !maxW2If || !currentIf)
      continue;

    double minQ2 = 0.0, maxQ2 = 0.0, minY = 0.0, maxY = 0.0, minW2 = 0.0, maxW2 = 0.0;
    if (!parseDoubleValue(minQ2If->exec(*cutObject, "get", ""), minQ2) ||
        !parseDoubleValue(maxQ2If->exec(*cutObject, "get", ""), maxQ2) ||
        !parseDoubleValue(minYIf->exec(*cutObject, "get", ""), minY) ||
        !parseDoubleValue(maxYIf->exec(*cutObject, "get", ""), maxY) ||
        !parseDoubleValue(minW2If->exec(*cutObject, "get", ""), minW2) ||
        !parseDoubleValue(maxW2If->exec(*cutObject, "get", ""), maxW2)) {
      continue;
    }

    std::string currentOption;
    if (!parseSwitchOptionName(currentIf->exec(*cutObject, "get", ""), currentOption))
      continue;

    const bool chargedCurrent = (currentOption == "Charged");
    const bool neutralCurrent = (currentOption == "Neutral");
    if (!chargedCurrent && !neutralCurrent) continue;
    if (chargedCurrent != usesChargedCurrentDISWindow()) continue;

    matchedWindow.available = true;
    matchedWindow.minQ2 = minQ2 * GeV2;
    matchedWindow.maxQ2 = maxQ2 * GeV2;
    matchedWindow.minY = minY;
    matchedWindow.maxY = maxY;
    matchedWindow.minW2 = minW2 * GeV2;
    matchedWindow.maxW2 = maxW2 * GeV2;
    ++matchedCount;
  }

  if (matchedCount != 1) {
    std::ostringstream msg;
    msg << "DISBase::generateKinematics falling back to legacy Born generation because ";
    if (matchedCount == 0) msg << "no matching DIS window cut";
    else msg << matchedCount << " matching DIS window cuts";
    msg << " were found for the current "
        << (usesChargedCurrentDISWindow() ? "charged-current" : "neutral-current")
        << " DIS process.\n";
    generator()->logWarning(Exception(msg.str(), Exception::warning));
    return false;
  }

  nativeDISWindow_.available = true;
  nativeDISWindow_.minQ2 = matchedWindow.minQ2;
  nativeDISWindow_.maxQ2 = matchedWindow.maxQ2;
  nativeDISWindow_.minY = matchedWindow.minY;
  nativeDISWindow_.maxY = matchedWindow.maxY;
  nativeDISWindow_.minW2 = matchedWindow.minW2;
  nativeDISWindow_.maxW2 = matchedWindow.maxW2;
  return true;
}

// Identify which incoming Born particle is the hadron-side beam and return the
// Bjorken x attached to that beam.
bool DISBase::determineBornHadronAndXB(tcBeamPtr & hadron, double & xB) const {
  hadron = tcBeamPtr();
  xB = 0.0;

  if (lastParticles().first && lastParticles().first->dataPtr() &&
      HadronMatcher::Check(*lastParticles().first->dataPtr())) {
    hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
    xB = lastX1();
    return true;
  }

  if (lastParticles().second && lastParticles().second->dataPtr() &&
      HadronMatcher::Check(*lastParticles().second->dataPtr())) {
    hadron = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
    xB = lastX2();
    return true;
  }

  return false;
}

// Convert the active DIS window in Q2, y, and W2 into tighter Born cos(theta)
// limits for the already-sampled xB point.
bool DISBase::tightenBornCosThetaWithNativeDISWindow(const HwMEBase::TwoToTwoKinematicsSetup & setup,
                                                     double xB,
                                                     double & ctmin,
                                                     double & ctmax) const {
  if (!resolveNativeDISWindow()) return true;
  if (xB <= 0.0 || xB >= 1.0) return false;
  if (setup.pq <= ZERO) return false;

  const Energy2 sBeam = lastCuts().SMax();
  if (sBeam <= ZERO) return false;

  const auto q2FromCth = [&](double cth) -> Energy2 {
    return setup.e0e2 - setup.m22 - setup.pq * cth;
  };
  const auto cthFromQ2 = [&](Energy2 q2) -> double {
    return (setup.e0e2 - setup.m22 - q2) / setup.pq;
  };
  const auto applyMinQ2Bound = [&](Energy2 q2min) {
    if (q2min > q2FromCth(ctmax)) ctmax = min(ctmax, cthFromQ2(q2min));
  };
  const auto applyMaxQ2Bound = [&](Energy2 q2max) {
    if (q2max < q2FromCth(ctmin)) ctmin = max(ctmin, cthFromQ2(q2max));
  };

  applyMinQ2Bound(nativeDISWindow_.minQ2);
  applyMaxQ2Bound(nativeDISWindow_.maxQ2);

  const Energy2 yMinQ2 = xB * sBeam * nativeDISWindow_.minY;
  const Energy2 yMaxQ2 = xB * sBeam * nativeDISWindow_.maxY;
  applyMinQ2Bound(yMinQ2);
  applyMaxQ2Bound(yMaxQ2);

  const double wFactor = xB / (1.0 - xB);
  applyMinQ2Bound(wFactor * nativeDISWindow_.minW2);
  applyMaxQ2Bound(wFactor * nativeDISWindow_.maxW2);

  return ctmin < ctmax;
}

// Register the public DIS interfaces retained for production generation,
// fixed-order comparisons, real-emission spin vertices, and correlated Rivet
// helicity weights.
void DISBase::Init() {
  
  static ClassDocumentation<DISBase> documentation
    ("The DISBase class provides the base class for the "
     "implementation of DIS type processes including the "
     "hard corrections in either the old-fashioned matrix "
     "element correction of POWHEG approaches");

  static Parameter<DISBase,double> interfaceProcessProbability
    ("ProcessProbability",
     "The probabilty of the QCD compton process for the process selection",
     &DISBase::procProb_, 0.3, 0.0, 1.,
     false, false, Interface::limited);

  static Reference<DISBase,ShowerAlpha> interfaceCoupling
    ("Coupling",
     "Pointer to the object to calculate the coupling for the correction",
     &DISBase::alpha_, false, false, true, true, false);
  
  static Parameter<DISBase,Energy> interfacepTMin
    ("pTMin",
     "The minimum pT",
     &DISBase::pTmin_, GeV, 1.*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<DISBase,double> interfaceComptonWeight
    ("ComptonWeight",
     "Weight for the overestimate ofthe compton channel",
     &DISBase::comptonWeight_, 50.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<DISBase,double> interfaceBGFWeight
    ("BGFWeight",
     "Weight for the overestimate of the BGF channel",
     &DISBase::BGFWeight_, 100.0, 0.0, 1000.0,
     false, false, Interface::limited);

  static Switch<DISBase,unsigned int> interfaceContribution
    ("Contribution",
     "Which contributions to the cross section to include",
     &DISBase::contrib_, 0, false, false);
  static SwitchOption interfaceContributionLeadingOrder
    (interfaceContribution,
     "LeadingOrder",
     "Just generate the leading order cross section",
     0);
  static SwitchOption interfaceContributionPositiveNLO
    (interfaceContribution,
     "PositiveNLO",
     "Generate the positive contribution to the full NLO cross section",
     1);
  static SwitchOption interfaceContributionNegativeNLO
    (interfaceContribution,
     "NegativeNLO",
     "Generate the negative contribution to the full NLO cross section",
     2);

  static Switch<DISBase,unsigned int> interfaceScaleOption
    ("ScaleOption",
     "Option for the choice of factorization (and renormalization) scale",
     &DISBase::scaleOpt_, 1, false, false);
  static SwitchOption interfaceDynamic
    (interfaceScaleOption,
     "Dynamic",
     "Dynamic factorization scale equal to the current sqrt(sHat())",
     1);
  static SwitchOption interfaceFixed
    (interfaceScaleOption,
     "Fixed",
     "Use a fixed factorization scale set with FactorizationScaleValue",
     2);

  static Parameter<DISBase,Energy> interfaceFactorizationScale
    ("FactorizationScale",
     "Value to use in the event of a fixed factorization scale",
     &DISBase::muF_, GeV, 100.0*GeV, 1.0*GeV, 500.0*GeV,
     true, false, Interface::limited);

  static Parameter<DISBase,double> interfaceScaleFactor
    ("ScaleFactor",
     "The factor used before Q2 if using a running scale",
     &DISBase::scaleFact_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<DISBase,double> interfaceSamplingPower
    ("SamplingPower",
     "Power for the sampling of xp",
     &DISBase::power_, 0.6, 0.0, 1.,
     false, false, Interface::limited);

  static Switch<DISBase,bool> interfaceUseFixedOrderAlphaSInPOWHEGEmission
    ("UseFixedOrderAlphaSInPOWHEGEmission",
     "Use a fixed-order alphaS source instead of the shower alphaS object for POWHEG/MEC emission weights, preferring LHAPDF alphaS when the beam PDF provides it.",
     &DISBase::useFixedOrderAlphaSInPOWHEGEmission_, false, false, false);
  static SwitchOption interfaceUseFixedOrderAlphaSInPOWHEGEmissionNo
    (interfaceUseFixedOrderAlphaSInPOWHEGEmission,
     "No",
     "Use the current shower-based alphaS implementation for POWHEG/MEC emission weights.",
     false);
  static SwitchOption interfaceUseFixedOrderAlphaSInPOWHEGEmissionYes
    (interfaceUseFixedOrderAlphaSInPOWHEGEmission,
     "Yes",
     "Use a fixed-order alphaS source for POWHEG/MEC emission weights, preferring LHAPDF alphaS when available from the beam PDF.",
     true);

  static Switch<DISBase,bool> interfaceUseQ2ScaleInPOWHEGEmission
    ("UseQ2ScaleInPOWHEGEmission",
     "Use Q2 instead of the native POWHEG emission scale for the POWHEG/MEC alphaS evaluation and emission-PDF numerator scale.",
     &DISBase::useQ2ScaleInPOWHEGEmission_, false, false, false);
  static SwitchOption interfaceUseQ2ScaleInPOWHEGEmissionNo
    (interfaceUseQ2ScaleInPOWHEGEmission,
     "No",
     "Use the current POWHEG emission scales.",
     false);
  static SwitchOption interfaceUseQ2ScaleInPOWHEGEmissionYes
    (interfaceUseQ2ScaleInPOWHEGEmission,
     "Yes",
     "Use Q2 as the alphaS scale and emission-PDF numerator scale in the POWHEG/MEC emission weights.",
     true);

  static Switch<DISBase,bool> interfaceGenerateRivetWeights
    ("GenerateRivetWeights",
     "Attach correlated NC helicity weights for Rivet analyses.",
     &DISBase::generateRivetWeights_, false, false, false);
  static SwitchOption interfaceGenerateRivetWeightsNo
    (interfaceGenerateRivetWeights,
     "No",
     "Do not attach correlated Rivet helicity weights.",
     false);
  static SwitchOption interfaceGenerateRivetWeightsYes
    (interfaceGenerateRivetWeights,
     "Yes",
     "Attach HERWIG_DIS_PP/PM/MP/MM/SIGMA/DELTA_LL optional weights.",
     true);

  static Switch<DISBase,bool> interfaceUseNativeDISWindowGeneration
    ("UseNativeDISWindowGeneration",
     "Tighten the underlying Born cos(theta) generation to the active SimpleDISCut DIS window before applying the final passCuts() safety check.",
     &DISBase::useNativeDISWindowGeneration_, false, false, false);
  static SwitchOption interfaceUseNativeDISWindowGenerationNo
    (interfaceUseNativeDISWindowGeneration,
     "No",
     "Keep the legacy HwMEBase Born generation and rely on passCuts() for the full DIS window.",
     false);
  static SwitchOption interfaceUseNativeDISWindowGenerationYes
    (interfaceUseNativeDISWindowGeneration,
     "Yes",
     "Tighten the Born angular generation directly to the DIS window when a unique matching SimpleDISCut is available.",
     true);

  static Switch<DISBase,unsigned int> interfacePOWHEGEmissionComparisonMode
    ("POWHEGEmissionComparisonMode",
     "Optional comparison-mode ladder for the POWHEG hardest emission. The default keeps the native Herwig POWHEG behaviour exactly.",
     &DISBase::powhegEmissionComparisonMode_, 0, false, false);
  static SwitchOption interfacePOWHEGEmissionComparisonModeDefault
    (interfacePOWHEGEmissionComparisonMode,
     "Default",
     "Use the native Herwig POWHEG hardest-emission generation.",
     DISBase::POWHEGEmissionComparisonModeDefault);
  static SwitchOption interfacePOWHEGEmissionComparisonModeRealOnly
    (interfacePOWHEGEmissionComparisonMode,
     "RealOnly",
     "Use native POWHEG emission generation, but veto instead of falling back to Born when no real emission is accepted.",
     DISBase::POWHEGEmissionComparisonModeRealOnly);

  static Parameter<DISBase,unsigned long> interfacePOWHEGEmissionComparisonMaxAttempts
    ("POWHEGEmissionComparisonMaxAttempts",
     "Maximum number of real-emission attempts in the non-default POWHEG comparison modes.",
     &DISBase::powhegEmissionComparisonMaxAttempts_, 100, 1, 100000,
     false, false, Interface::limited);

  static Switch<DISBase,bool> interfaceUsePOWHEGRealSpinVertex
    ("UsePOWHEGRealSpinVertex",
     "Attach an exact spin-only HardVertex to POWHEG real-emission states.",
     &DISBase::usePOWHEGRealSpinVertex_, false, false, false);
  static SwitchOption interfaceUsePOWHEGRealSpinVertexNo
    (interfaceUsePOWHEGRealSpinVertex,
     "No",
     "Do not attach a real-emission spin vertex.",
     false);
  static SwitchOption interfaceUsePOWHEGRealSpinVertexYes
    (interfaceUsePOWHEGRealSpinVertex,
     "Yes",
     "Attach a spin-only HardVertex to the POWHEG real-emission state.",
     true);

}

// Complete setup after the repository phase. The analytic Compton/BGF
// integrals are overestimate normalizations used by the emission samplers.
void DISBase::doinit() {
  HwMEBase::doinit();
  nativeDISWindow_ = NativeDISWindowDefinition();
  if ((hasMECorrection() || hasPOWHEGCorrection() != No) && !alpha_) {
    throw InitException()
      << "DISBase requires a non-null Coupling reference when MEC or POWHEG"
      << " corrections are enabled."
      << Exception::abortnow;
  }
  // Analytic channel integrals provide the normalization for rejection
  // sampling in the legacy MEC path.
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  comptonInt_ = 2.*(-21./20.-6./(5.*r5)*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
  bgfInt_ = 121./9.-56./r5*ath;
  gluon_ = getParticleData(ParticleID::g);
}

// Choose the alphaS scale for accepted real-emission generation.
Energy2 DISBase::powhegEmissionAlphaSScale(Energy2 q2, double xT) const {
  return useQ2ScaleInPOWHEGEmission_ ? q2 : 0.25*q2*sqr(xT);
}

// Choose the PDF numerator scale for emission PDF ratios.
Energy2 DISBase::powhegEmissionPDFScale(Energy2 q2, Energy2 mappedScale) const {
  return useQ2ScaleInPOWHEGEmission_ ? q2 : mappedScale;
}

// Fixed-order alphaS source used by release comparison modes when requested.
double DISBase::powhegEmissionFixedOrderAlphaSValue(Energy2 scale) const {
  if (pdf_) {
    const ThePEG::LHAPDF *lhaPdf = dynamic_cast<const ThePEG::LHAPDF *>(pdf_.operator->());
    if (lhaPdf) return lhaPdf->alphaS(scale);
  }
  return SM().alphaS(scale);
}

// AlphaS value entering the emission weight.
double DISBase::powhegEmissionAlphaSValue(Energy2 scale) const {
  return useFixedOrderAlphaSInPOWHEGEmission_ ? powhegEmissionFixedOrderAlphaSValue(scale)
                                              : alpha_->value(scale);
}

// Coupling overestimate used to generate trial emissions.
double DISBase::powhegEmissionAlphaSOverestimate(Energy2 referenceScale) const {
  return useFixedOrderAlphaSInPOWHEGEmission_ ? powhegEmissionFixedOrderAlphaSValue(referenceScale) :
    alpha_->overestimateValue();
}

// Ratio of the emission coupling to its overestimate for veto acceptance.
double DISBase::powhegEmissionAlphaSRatio(Energy2 scale, Energy2 referenceScale) const {
  if(!useFixedOrderAlphaSInPOWHEGEmission_) return alpha_->ratio(scale);
  const double alphaRef = powhegEmissionFixedOrderAlphaSValue(referenceScale);
  return alphaRef > 0.0 ? powhegEmissionFixedOrderAlphaSValue(scale)/alphaRef : 0.0;
}

// Extract the Born DIS legs, PDFs, Breit-frame rotation, and polarization state
// needed by the native POWHEG hardest-emission samplers.
void DISBase::initializePOWHEGEmissionState(RealEmissionProcessPtr born,
					    PPtr quark[2], PPtr lepton[2],
					    PPtr & hadron,
					    unsigned int & iqIn,
					    unsigned int & iqOut) {
  hadron = PPtr();
  iqIn = iqOut = 0;
  for(unsigned int ix=0; ix<2; ++ix) {
    quark[ix] = PPtr();
    lepton[ix] = PPtr();
  }

  for(unsigned int ix=0; ix<born->bornIncoming().size(); ++ix) {
    if(QuarkMatcher::Check(born->bornIncoming()[ix]->data())) {
      iqIn=ix;
      hadron = born->hadrons()[ix];
      quark [0] = born->bornIncoming()[ix];
      beam_ = dynamic_ptr_cast<tcBeamPtr>(hadron->dataPtr());
      hadron_ = beam_;
      xB_ = quark[0]->momentum().rho()/hadron->momentum().rho();
    }
    else if(LeptonMatcher::Check(born->bornIncoming()[ix]->data())) {
      lepton[0] = born->bornIncoming()[ix];
      leptons_[0] = lepton[0]->dataPtr();
    }
  }
  pdf_=beam_->pdf();
  assert(beam_&&pdf_&&quark[0]&&lepton[0]);

  for(unsigned int ix=0; ix<born->bornOutgoing().size(); ++ix) {
    if(QuarkMatcher::Check(born->bornOutgoing()[ix]->data())) {
      iqOut=ix;
      quark [1] = born->bornOutgoing()[ix];
    }
    else if(LeptonMatcher::Check(born->bornOutgoing()[ix]->data())) {
      lepton[1] = born->bornOutgoing()[ix];
      leptons_[1] = lepton[1]->dataPtr();
    }
  }
  assert(quark[1]&&lepton[1]);

  for(unsigned int ix=0; ix<2; ++ix) partons_[ix] = quark[ix]->dataPtr();

  q_ =lepton[0]->momentum()-lepton[1]->momentum();
  q2_ = -q_.m2();
  double yB =
    (                   q_*quark[0]->momentum())/
    (lepton[0]->momentum()*quark[0]->momentum());
  l_ = 2./yB-1.;

  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleRho();
  Lorentz5Momentum pb     = quark[0]->momentum();
  Axis axis(q_.vect().unit());
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot_ = LorentzRotation();
  if(axis.perp2()>1e-20) {
    rot_.setRotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
    rot_.rotateX(Constants::pi);
  }
  if(abs(1.-q_.e()/q_.vect().mag())>1e-6) rot_.boostZ( q_.e()/q_.vect().mag());
  pb *= rot_;
  if(pb.perp2()/GeV2>1e-20) {
    Boost trans = -1./pb.e()*pb.vect();
    trans.setZ(0.);
    rot_.boost(trans);
  }
  Lorentz5Momentum pl    = rot_*lepton[0]->momentum();
  rot_.rotateZ(-atan2(pl.y(),pl.x()));
  pl_[0]=rot_*lepton[0]->momentum();
  pl_[1]=rot_*lepton[1]->momentum();
  pq_[0]=rot_* quark[0]->momentum();
  pq_[1]=rot_* quark[1]->momentum();
  q_ *= rot_;

  leptonPolarization_ = 0.0;
  if(lepton[0]) leptonPolarization_ = longPol(lepton[0]->spinInfo()->rhoMatrix());
}

// Generate the native POWHEG hardest QCD emission, choosing the larger-pT
// candidate between QCDC/Compton and BGF and rebuilding the event record.
RealEmissionProcessPtr DISBase::generateNativePOWHEGHardest(RealEmissionProcessPtr born,
								    bool allowBornFallback) {
  PPtr quark[2],lepton[2];
  PPtr hadron;
  unsigned int iqIn(0),iqOut(0);
  initializePOWHEGEmissionState(born, quark, lepton, hadron, iqIn, iqOut);

  generateCompton();
  generateBGF();
  if(pTCompton_<ZERO&&pTBGF_<ZERO) {
    if(allowBornFallback) {
      born->pT()[ShowerInteraction::QCD] = pTmin_;
      return born;
    }
    return RealEmissionProcessPtr();
  }

  bool isCompton=pTCompton_>pTBGF_;
  bool isQuark = quark[0]->colourLine();
  bool FSR = false;
  if(iqIn==1)
    born->incoming().push_back(born->bornIncoming()[0]->dataPtr()->
			       produceParticle(born->bornIncoming()[0]->momentum()));
  if(iqOut==1)
    born->outgoing().push_back(born->bornOutgoing()[0]->dataPtr()->
			       produceParticle(born->bornOutgoing()[0]->momentum()));
  PPtr newout,newin,emitted;
  if(isCompton) {
    rot_.invert();
    for(unsigned int ix=0;ix<ComptonMomenta_.size();++ix) {
      ComptonMomenta_[ix].transform(rot_);
    }
    newout  = partons_[1]->produceParticle(ComptonMomenta_[1]);
    emitted = gluon_     ->produceParticle(ComptonMomenta_[2]);
    newin   = partons_[0]->produceParticle(ComptonMomenta_[0]);
    emitted->incomingColour(newin,!isQuark);
    emitted->colourConnect(newout,!isQuark);
    FSR = !ComptonISFS_;
    born->pT()[ShowerInteraction::QCD] = pTCompton_;
  }
  else {
    rot_.invert();
    for(unsigned int ix=0;ix<BGFMomenta_.size();++ix) {
      BGFMomenta_[ix].transform(rot_);
    }
    newin   = gluon_                   ->produceParticle(BGFMomenta_[0]);
    emitted = quark[0]->dataPtr()->CC()->produceParticle(BGFMomenta_[2]);
    newout  = quark[1]->dataPtr()      ->produceParticle(BGFMomenta_[1]);
    emitted->incomingColour(newin, isQuark);
    newout ->incomingColour(newin,!isQuark);
	    FSR = false;
	    born->pT()[ShowerInteraction::QCD] = pTBGF_;
	  }
	  double x = newin->momentum().rho()/hadron->momentum().rho();
	  if(born->incoming().size()==0)
	    born->x(make_pair(x,1.));
  else
    born->x(make_pair(1.,x));
  if(FSR) {
    born->emitter(born->outgoing().size()+2);
    born->spectator(born->incoming().size());
  }
  else {
    born->emitter(born->incoming().size());
    born->spectator(born->outgoing().size()+2);
  }
  born->emitted(4);
  born->incoming().push_back(newin );
  born->outgoing().push_back(newout);
  if(iqIn==0)
    born->incoming().push_back(born->bornIncoming()[1]->dataPtr()->
			       produceParticle(born->bornIncoming()[1]->momentum()));
  if(iqOut==0)
    born->outgoing().push_back(born->bornOutgoing()[1]->dataPtr()->
			       produceParticle(born->bornOutgoing()[1]->momentum()));
  born->outgoing().push_back(emitted);
  born->interaction(ShowerInteraction::QCD);
  if(usePOWHEGRealSpinVertex_) {
    constructRealEmissionSpinVertex(born,isCompton);
  }
  return born;
}

// RealOnly comparison mode: repeatedly request a genuine real emission and
// veto the event rather than falling back to Born kinematics.
RealEmissionProcessPtr
DISBase::generateComparisonModePOWHEGHardest(RealEmissionProcessPtr born) {
  const unsigned long maxAttempts =
    std::max<unsigned long>(1UL, powhegEmissionComparisonMaxAttempts_);

  for(unsigned long attempt=0; attempt<maxAttempts; ++attempt) {
    RealEmissionProcessPtr real = generateNativePOWHEGHardest(new_ptr(*born), false);
    if(real) return real;
  }

  throw Veto();
}

// Prepare the Born kinematic state used by the older matrix-element correction
// path. This mirrors the POWHEG setup but does not generate a real emission.
void DISBase::initializeMECorrection(RealEmissionProcessPtr born, double & initial,
				     double & final) {
  initial = initial_;
  final   = final_;
  leptonPolarization_ = 0.0;
  // incoming particles
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
      if(QuarkMatcher::Check(born->bornIncoming()[ix]->data())) {
      partons_[0] = born->bornIncoming()[ix]->dataPtr();
      pq_[0] = born->bornIncoming()[ix]->momentum();
      tcBeamPtr beam = dynamic_ptr_cast<tcBeamPtr>(born->hadrons()[ix]->dataPtr());
      beam_ = beam;
      hadron_ = beam;
      if (beam) {
        xB_ = pq_[0].rho()/born->hadrons()[ix]->momentum().rho();
        pdf_ = beam->pdf();
      }
    }
    else if(LeptonMatcher::Check(born->bornIncoming()[ix]->data())) {
      leptons_[0] = born->bornIncoming()[ix]->dataPtr();
      pl_[0] = born->bornIncoming()[ix]->momentum();
      leptonPolarization_ =
        longPol(born->bornIncoming()[ix]->spinInfo()->rhoMatrix());
    }
  }
  // outgoing particles
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(QuarkMatcher::Check(born->bornOutgoing()[ix]->data())) {
      partons_[1] = born->bornOutgoing()[ix]->dataPtr();
      pq_[1] = born->bornOutgoing()[ix]->momentum();
    }
    else if(LeptonMatcher::Check(born->bornOutgoing()[ix]->data())) {
      leptons_[1] = born->bornOutgoing()[ix]->dataPtr();
      pl_[1] = born->bornOutgoing()[ix]->momentum();
    }
  }
  // The Breit-frame DIS invariants drive both the MEC kernels and the
  // longitudinal spin factors.
  q_ =pl_[0]-pl_[1];
  q2_ = -q_.m2();
  double  yB = (q_*pq_[0])/(pl_[0]*pq_[0]); 
  l_ = 2./yB-1.;
}

// Apply the hard matrix-element correction by sampling one local DIS real
// emission and accepting it against the exact Compton/BGF kernel.
RealEmissionProcessPtr DISBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  static const double eps=1e-6;
  leptonPolarization_ = 0.0;
  // find the incoming and outgoing quarks and leptons
  PPtr quark[2],lepton[2];
  PPtr hadron;
  unsigned int iqIn(0),iqOut(0);
  // incoming particles
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
    if(QuarkMatcher::Check(born->bornIncoming()[ix]->data())) {
      iqIn=ix;
      quark[0] = born->bornIncoming()[ix];
      hadron   = born->hadrons()[ix];     
      beam_    = dynamic_ptr_cast<tcBeamPtr>(hadron->dataPtr());
      hadron_  = beam_;
      xB_ = quark[0]->momentum().rho()/hadron->momentum().rho();
    }
    else if(LeptonMatcher::Check(born->bornIncoming()[ix]->data())) {
      lepton[0] = born->bornIncoming()[ix];
      leptonPolarization_ =
        longPol(born->bornIncoming()[ix]->spinInfo()->rhoMatrix());
    }
  }
  pdf_ = beam_->pdf();
  assert(beam_&&pdf_&&quark[0]&&lepton[0]);
  // outgoing particles
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(QuarkMatcher::Check(born->bornOutgoing()[ix]->data())) {
      iqOut=ix;
      quark [1] = born->bornOutgoing()[ix];
    }
    else if(LeptonMatcher::Check(born->bornOutgoing()[ix]->data())) {
      lepton[1] = born->bornOutgoing()[ix];
    }
  }
  // momentum fraction
  assert(quark[1]&&lepton[1]);
  // calculate the matrix element
  AzimuthalKernelCoefficients azicoeff;
  // select the type of process
  bool BGF = UseRandom::rnd()>procProb_;
  double xp,zp,wgt,x1,x2,x3,xperp;
  // generate a QCD compton process
  if(!BGF) {
    wgt = generateComptonPoint(xp,zp);
    if(xp<eps) return RealEmissionProcessPtr();
    // common pieces
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    Energy2 pdfScale = powhegEmissionPDFScale(q2_, scale);
    Energy2 alphaScale = useQ2ScaleInPOWHEGEmission_ ? q2_ : scale;
    wgt *= 2./3./Constants::pi*powhegEmissionAlphaSValue(alphaScale)/procProb_;
    // PDF piece
    double pdfRatioGen = 0.0;
    if (!powhegEmissionPDFRatioForGeneration(pdf_, beam_,
                                             quark[0]->dataPtr(),
                                             quark[0]->dataPtr(),
                                             pdfScale, q2_,
                                             xB_/xp, xB_,
                                             pdfRatioGen)) {
      return RealEmissionProcessPtr();
    }
    wgt *= pdfRatioGen;
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = ComptonME(xp,x2,xperp,true);
  }
  // generate a BGF process
  else {
    wgt = generateBGFPoint(xp,zp);
    if(xp<eps) return RealEmissionProcessPtr();
    // common pieces 
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1);
    Energy2 pdfScale = powhegEmissionPDFScale(q2_, scale);
    Energy2 alphaScale = useQ2ScaleInPOWHEGEmission_ ? q2_ : scale;
    wgt *= 0.25/Constants::pi*powhegEmissionAlphaSValue(alphaScale)/(1.-procProb_);
    // PDF piece
    double pdfRatioGen = 0.0;
    if (!powhegEmissionPDFRatioForGeneration(pdf_, beam_,
                                             gluon_,
                                             quark[0]->dataPtr(),
                                             pdfScale, q2_,
                                             xB_/xp, xB_,
                                             pdfRatioGen)) {
      return RealEmissionProcessPtr();
    }
    wgt *= pdfRatioGen;
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,true);
  }
  // compute the azimuthal average of the weight
  wgt *= azicoeff.average();
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return RealEmissionProcessPtr();
  // if generate generate phi
  unsigned int itry(0);
  double phimax = azimuthalKernelMaximum(azicoeff.c0, azicoeff.c1, azicoeff.c2);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azimuthalKernelValue(azicoeff.c0, azicoeff.c1, azicoeff.c2, cphi);
    if (phiwgt > phimax * (1.0 + 1e-12)) {
      ostringstream wstring;
      wstring << "DISBase::applyHardMatrixElementCorrection() "
              << "Azimuthal envelope undershot kernel value"
              << " phimax = " << phimax
              << " phiwgt = " << phiwgt
              << " cphi = " << cphi << "\n";
      generator()->logWarning(Exception(wstring.str(), Exception::warning));
    }
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::applyHardMatrixElementCorrection() to"
				  << " generate phi" << Exception::eventerror;
  // construct lorentz transform from lab to breit frame
  Lorentz5Momentum phadron =  hadron->momentum();
  phadron.setMass(0.*GeV);
  phadron.rescaleEnergy();
  Lorentz5Momentum pcmf = phadron+0.5/xB_*q_;
  pcmf.rescaleMass();
  LorentzRotation rot(-pcmf.boostVector());
  Lorentz5Momentum pbeam = rot*phadron;
  Axis axis(pbeam.vect().unit());
  double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
  rot.rotate(-acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
  Lorentz5Momentum pl    = rot*pl_[0];
  rot.rotateZ(-atan2(pl.y(),pl.x()));
  pl_[0] *= rot;
  pl_[1] *= rot;
  pq_[0] *= rot;
  pq_[1] *= rot;
  // compute the new incoming and outgoing momenta
  Energy Q(sqrt(q2_));
  Lorentz5Momentum p1 = Lorentz5Momentum( 0.5*Q*xperp*cos(phi), 0.5*Q*xperp*sin(phi),
					  -0.5*Q*x2,0.*GeV,0.*GeV);
  p1.rescaleEnergy();
  Lorentz5Momentum p2 = Lorentz5Momentum(-0.5*Q*xperp*cos(phi),-0.5*Q*xperp*sin(phi),
					 -0.5*Q*x3,0.*GeV,0.*GeV);
  p2.rescaleEnergy();
  Lorentz5Momentum pin(0.*GeV,0.*GeV,-0.5*x1*Q,-0.5*x1*Q,0.*GeV);
  // we need the Lorentz transform back to the lab
  rot.invert();
  // transform the momenta to lab frame
  pin *= rot;
  p1  *= rot;
  p2  *= rot;
  // test to ensure outgoing particles can be put on-shell
  if(!BGF) {
    if(p1.e()<quark[1]->dataPtr()->constituentMass())  return RealEmissionProcessPtr();
    if(p2.e()<gluon_              ->constituentMass()) return RealEmissionProcessPtr();
  }
  else {
    if(p1.e()<quark[1]->dataPtr()      ->constituentMass()) return RealEmissionProcessPtr();
    if(p2.e()<quark[0]->dataPtr()->CC()->constituentMass()) return RealEmissionProcessPtr();
  }
  // create the new particles and real emission process
  bool isQuark = quark[0]->colourLine();
  bool FSR = false;
  // incoming lepton if first
  if(iqIn==1)
    born->incoming().push_back(born->bornIncoming()[0]->dataPtr()->
			       produceParticle(born->bornIncoming()[0]->momentum()));
  // outgoing lepton if first
  if(iqOut==1)
    born->outgoing().push_back(born->bornOutgoing()[0]->dataPtr()->
			       produceParticle(born->bornOutgoing()[0]->momentum()));
  PPtr newin,newout,emitted;
  // radiating system
  if(!BGF) {
    newin   = quark[0]->dataPtr()->produceParticle(pin);
    emitted = gluon_              ->produceParticle(p2 );
    newout  = quark[1]->dataPtr()->produceParticle(p1 );
    emitted->incomingColour(newin,!isQuark);
    emitted->colourConnect(newout,!isQuark);
    FSR = xp>zp;
  }
  else {
    newin   = gluon_                   ->produceParticle(pin);
    emitted = quark[0]->dataPtr()->CC()->produceParticle(p2 );
    newout  = quark[1]->dataPtr()      ->produceParticle(p1 );
    emitted->incomingColour(newin, isQuark);
	    newout ->incomingColour(newin,!isQuark);
	    FSR = false;
	  }
	  // set x
	  double x(xB_/xp);
  if(born->incoming().size()==0)
    born->x(make_pair(x,1.));
  else
    born->x(make_pair(1.,x));
  if(FSR) {
    born->emitter(born->outgoing().size()+2);
    born->spectator(born->incoming().size());
  }
  else {
    born->emitter(born->incoming().size());
    born->spectator(born->outgoing().size()+2);
  }
  born->emitted(4);
  // radiating particles
  born->incoming().push_back(newin );
  born->outgoing().push_back(newout);
  // incoming lepton if second
  if(iqIn==0)
    born->incoming().push_back(born->bornIncoming()[1]->dataPtr()->
			       produceParticle(born->bornIncoming()[1]->momentum()));
  // outgoing lepton if second
  if(iqOut==0)
    born->outgoing().push_back(born->bornOutgoing()[1]->dataPtr()->
			       produceParticle(born->bornOutgoing()[1]->momentum()));
  // radiated particle
  born->outgoing().push_back(emitted);
  born->interaction(ShowerInteraction::QCD);
  born->pT()[ShowerInteraction::QCD] = 0.5*Q*xperp;
  if(usePOWHEGRealSpinVertex_) {
    constructRealEmissionSpinVertex(born,!BGF);
  }
  return born;
}

// Veto soft shower branchings with the DIS matrix-element correction weight so
// the shower reproduces the hard-emission limit in both ISR and FSR regions.
bool DISBase::softMatrixElementVeto(PPtr parent,
				    PPtr progenitor,
				    const bool & fs,
				    const Energy & highestpT,
				    const vector<tcPDPtr> & ids,
				    const double & z,
				    const Energy & scale,
				    const Energy & pT) {
  bool veto = !UseRandom::rndbool(fs ? 1./final_ : 1./initial_);
  // check if me correction should be applied
  long id[2]={progenitor->id(),parent->id()};
  if(id[0]!=id[1]||id[1]==ParticleID::g) return veto; 
  // check if hardest so far
  if(pT<highestpT) return veto;
  double kappa(sqr(scale)/q2_);
  double zk((1.-z)*kappa);
  // final-state
  double wgt(0.);
  if(fs) {
    double zp=z,xp=1./(1.+z*zk);
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x2 = 1.-(1.-zp)/xp;
    AzimuthalKernelCoefficients azicoeff = ComptonME(xp,x2,xperp,false);
    wgt = azicoeff.average()*xp/(1.+sqr(z))/final_;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for FSR in DISBase::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  else {
    double xp = 2.*z/(1.+zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double zp = 0.5* (1.-zk+sqrt(sqr(1.+zk)-4.*z*zk));
    double xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    double x1 = -1./xp, x2 = 1.-(1.-zp)/xp, x3 = 2.+x1-x2;
    // compton
    if(ids[0]->id()!=ParticleID::g) {
      AzimuthalKernelCoefficients azicoeff = ComptonME(xp,x2,xperp,false);
      wgt = azicoeff.average()*xp*(1.-z)/(1.-xp)/(1.+sqr(z))/
	(1.-zp+xp-2.*xp*(1.-zp));
    }
    // BGF
    else {
      AzimuthalKernelCoefficients azicoeff = BGFME(xp,x2,x3,xperp,true);
      wgt = azicoeff.average()*xp/(1.-zp+xp-2.*xp*(1.-zp))/(sqr(z)+sqr(1.-z));
    }
    wgt /=initial_;
    if(wgt<.0||wgt>1.) {
      ostringstream wstring;
      wstring << "Soft ME correction weight too large or "
	      << "negative for ISR in DISBase::"
	      << "softMatrixElementVeto() soft weight " 
	      << " xp = " << xp << " zp = " << zp
	      << " weight = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(), 
					 Exception::warning) );
    }
  }
  // return whether or not vetoed
  return !UseRandom::rndbool(wgt);
}

// Sample the xp,zp variables for a QCDC/Compton real-emission point and return
// the analytic integral normalization used by the MEC path.
double DISBase::generateComptonPoint(double &xp, double & zp) {
  static const double maxwgt = 1.;
  double wgt;
  do {
    xp  = UseRandom::rnd();
    double zpmin = xp, zpmax = 1./(1.+xp*(1.-xp));
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    if(UseRandom::rndbool()) swap(xp,zp);
    double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp,x2=1.-(1.-zp)/xp;
    wgt *= 2.*(1.+sqr(xp)*(sqr(x2)+1.5*xperp2))/(1.-xp)/(1.-zp);
    if(wgt>maxwgt) {
      ostringstream wstring;
      wstring << "DISBase::generateComptonPoint "
	      << "Weight greater than maximum "
	      << "wgt = " << wgt << " maxwgt = 1\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return comptonInt_;
}

// Sample the xp,zp variables for a BGF real-emission point and return the
// analytic integral normalization used by the MEC path.
double DISBase::generateBGFPoint(double &xp, double & zp) {
  static const double maxwgt = 25.;
  double wgt;
  do {
    xp = UseRandom::rnd();
    double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
    zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
    wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
    wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
    if(wgt>maxwgt) {
      ostringstream wstring;
      wstring << "DISBase::generateBGFPoint "
	      << "Weight greater than maximum "
	      << "wgt = " << wgt << " maxwgt = 1\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
  }
  while(wgt<UseRandom::rnd()*maxwgt);
  return bgfInt_;
}

// Return the normalized azimuthal kernel for the QCDC/Compton channel. The
// kernel is expressed as c0+c1*cos(phi)+c2*cos^2(phi) for efficient sampling.
DISBase::AzimuthalKernelCoefficients
DISBase::ComptonME(double xp, double x2, double xperp,
		   bool norm) const {
  AzimuthalKernelCoefficients output;
  double cos2 =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2 = xperp/sqrt(sqr(x2)+sqr(xperp));
  double root = sqrt(std::max(0.0, sqr(l_)-1.));
  const double PqBorn =
    mappedIncomingLongitudinalPolarization(partons_[0], xB_, q2_);
  const double aZero = A_pol(leptons_[0], leptons_[1],
                             partons_[0], partons_[1],
                             q2_, leptonPolarization_, 0.0);
  double aBorn = A_pol(leptons_[0], leptons_[1],
                       partons_[0], partons_[1],
                       q2_, leptonPolarization_, PqBorn);
  if (!powhegAnalyzingPowerIsSafe(aBorn)) aBorn = aZero;

  double aMapped = aBorn;
  double mappedDenRatio = 1.0;
  if (useMappedPolarizedEmissionKernel() && xp > 0.0) {
    const double PqMapped =
      mappedIncomingLongitudinalPolarization(partons_[0], xB_/xp, q2_);
    const double dBorn =
      realEmissionDenominatorFactor(leptons_[0], leptons_[1],
                                    partons_[0], partons_[1],
                                    q2_, leptonPolarization_, PqBorn);
    const double dMapped =
      realEmissionDenominatorFactor(leptons_[0], leptons_[1],
                                    partons_[0], partons_[1],
                                    q2_, leptonPolarization_, PqMapped);
    const double aCandidate = A_pol(leptons_[0], leptons_[1],
                                    partons_[0], partons_[1],
                                    q2_, leptonPolarization_, PqMapped);
    double candidateRatio = 1.0;

    // The mapped kernel is preferable when it is well behaved, but the veto
    // step should not be driven by a helicity decomposition sitting on a tiny
    // or sign-flipped denominator. In that corner we revert to the stable Born
    // kernel instead of manufacturing a spurious signed acceptance weight.
    if (powhegAnalyzingPowerIsSafe(aCandidate) &&
        powhegPositiveRatioIsSafe(dMapped, dBorn, candidateRatio)) {
      aMapped = aCandidate;
      mappedDenRatio = candidateRatio;
    }
  }

  double lo(1.+aBorn*l_+sqr(l_));
  if (std::abs(lo) <= 1e-30) lo = (lo < 0.0 ? -1.0 : 1.0) * 1e-30;
  double denom = norm ? 1.+sqr(xp)*(sqr(x2)+1.5*sqr(xperp)) : 1.;
  double fact  = mappedDenRatio*sqr(xp)*(sqr(x2)+sqr(xperp))/lo;
  output.c0 = (1. + fact*(sqr(cos2)+aMapped*cos2*l_+sqr(l_)))/denom;
  output.c1 = fact*(-aMapped*cos2*root*sin2-2.*l_*root*sin2)/denom;
  output.c2 = fact*(sqr(root)*sqr(sin2))/denom;
  return output;
}

// Return the normalized azimuthal kernel for the BGF channel. The two final
// quark lines get separate mapped spin projectors in full polarized DIS.
DISBase::AzimuthalKernelCoefficients
DISBase::BGFME(double xp, double x2, double x3,
	       double xperp, bool norm) const {
  AzimuthalKernelCoefficients output;
  double cos2  =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2  = xperp/sqrt(sqr(x2)+sqr(xperp));
  double fact2 = sqr(xp)*(sqr(x2)+sqr(xperp));
  double cos3  =   x3 /sqrt(sqr(x3)+sqr(xperp));
  double sin3  = xperp/sqrt(sqr(x3)+sqr(xperp));
  double fact3 = sqr(xp)*(sqr(x3)+sqr(xperp));
  double root = sqrt(std::max(0.0, sqr(l_)-1.));

  const double PqBorn =
    mappedIncomingLongitudinalPolarization(partons_[0], xB_, q2_);
  const double aZero = A_pol(leptons_[0], leptons_[1],
                             partons_[0], partons_[1],
                             q2_, leptonPolarization_, 0.0);
  double aBorn = A_pol(leptons_[0], leptons_[1],
                       partons_[0], partons_[1],
                       q2_, leptonPolarization_, PqBorn);
  if (!powhegAnalyzingPowerIsSafe(aBorn)) aBorn = aZero;

  double r2DenRatio = 1.0;
  double r3DenRatio = 1.0;
  double aR2 = aBorn;
  double aR3 = aBorn;

  if (useMappedPolarizedEmissionKernel() && xp > 0.0) {
    const double PgMapped =
      mappedIncomingLongitudinalPolarization(gluon_, xB_/xp, q2_);
    tcPDPtr qbarIn;
    if (partons_[0]) qbarIn = partons_[0]->CC();
    tcPDPtr qbarOut;
    if (partons_[1]) qbarOut = partons_[1]->CC();
    if (!qbarIn) qbarIn = partons_[0];
    if (!qbarOut) qbarOut = partons_[1];

    const double dBorn =
      realEmissionDenominatorFactor(leptons_[0], leptons_[1],
                                    partons_[0], partons_[1],
                                    q2_, leptonPolarization_, PqBorn);
    const double dR2 =
      realEmissionDenominatorFactor(leptons_[0], leptons_[1],
                                    partons_[0], partons_[1],
                                    q2_, leptonPolarization_, PgMapped);
    const double dR3 =
      realEmissionDenominatorFactor(leptons_[0], leptons_[1],
                                    qbarIn, qbarOut,
                                    q2_, leptonPolarization_, -PgMapped);
    const double aR2Candidate = A_pol(leptons_[0], leptons_[1],
                                      partons_[0], partons_[1],
                                      q2_, leptonPolarization_, PgMapped);
    const double aR3Candidate = A_pol(leptons_[0], leptons_[1],
                                      qbarIn, qbarOut,
                                      q2_, leptonPolarization_, -PgMapped);
    double candidateRatio = 1.0;
    if (powhegAnalyzingPowerIsSafe(aR2Candidate) &&
        powhegPositiveRatioIsSafe(dR2, dBorn, candidateRatio)) {
      aR2 = aR2Candidate;
      r2DenRatio = candidateRatio;
    }
    if (powhegAnalyzingPowerIsSafe(aR3Candidate) &&
        powhegPositiveRatioIsSafe(dR3, dBorn, candidateRatio)) {
      aR3 = aR3Candidate;
      r3DenRatio = candidateRatio;
    }
  }

  output.c0 = fact2*r2DenRatio*(sqr(cos2)+aR2*cos2*l_+sqr(l_)) +
              fact3*r3DenRatio*(sqr(cos3)-aR3*cos3*l_+sqr(l_));
  output.c1 = - fact2*r2DenRatio*(aR2*cos2*root*sin2+2.*l_*root*sin2)
              - fact3*r3DenRatio*(aR3*cos3*root*sin3-2.*l_*root*sin3);
  output.c2 = fact2*r2DenRatio*(sqr(root)*sqr(sin2)) +
              fact3*r3DenRatio*(sqr(root)*sqr(sin3));
  double lo(1.+aBorn*l_+sqr(l_));
  if (std::abs(lo) <= 1e-30) lo = (lo < 0.0 ? -1.0 : 1.0) * 1e-30;
  double denom = norm ? sqr(xp)*(sqr(x3)+sqr(x2)+3.*sqr(xperp))*lo : lo;
  output.c0 /= denom;
  output.c1 /= denom;
  output.c2 /= denom;
  return output;
}

// Entry point used by the POWHEG shower machinery. Non-QCD interactions are
// ignored; comparison modes replace Born fallback with an explicit veto.
RealEmissionProcessPtr DISBase::generateHardest(RealEmissionProcessPtr born,
						ShowerInteraction inter) {
  // check if generating QCD radiation
  if(inter!=ShowerInteraction::QCD && inter!=ShowerInteraction::QEDQCD &&
     inter!=ShowerInteraction::ALL)
    return RealEmissionProcessPtr();
  if(powhegEmissionComparisonMode_ != POWHEGEmissionComparisonModeDefault) {
    return generateComparisonModePOWHEGHardest(born);
  }
  return generateNativePOWHEGHardest(born, true);
}

// Generate a candidate QCDC/Compton POWHEG emission and cache its Breit-frame
// momenta. A negative pT marks that no emission above the current threshold was
// accepted.
void DISBase::generateCompton() {
  comptonRawXP_ = 0.0;
  comptonRawZP_ = 0.0;
  // maximum value of the xT
  double xT = sqrt((1.-xB_)/xB_);
  double xTMin = 2.*pTmin_/sqrt(q2_);
  double zp;
  Energy2 alphaRefScale = powhegEmissionAlphaSScale(q2_,xTMin);
  // prefactor
  double a = powhegEmissionAlphaSOverestimate(alphaRefScale)*comptonWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  AzimuthalKernelCoefficients azicoeff;
  while(true) {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // zp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_||xp>1.) {
      if(xT<=xTMin) break;
      continue;
    }
    // phase-space piece of the weight
    wgt = 8.*(1.-xp)*zp/comptonWeight_;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    Energy2 pdfScale = powhegEmissionPDFScale(q2_, scale);
    double pdfRatioGen = 0.0;
    if (!powhegEmissionPDFRatioForGeneration(pdf_, beam_,
                                             partons_[0], partons_[0],
                                             pdfScale, q2_,
                                             xB_/xp, xB_,
                                             pdfRatioGen)) {
      if(xT<=xTMin) break;
      continue;
    }
    wgt *= pdfRatioGen;
    // me piece of the weight
    double x2 = 1.-(1.-zp)/xp;
    azicoeff = ComptonME(xp,x2,xT,false);
    Energy2 alphaScale = powhegEmissionAlphaSScale(q2_,xT);
    wgt *= 4./3.*powhegEmissionAlphaSRatio(alphaScale,alphaRefScale)
      * azicoeff.average();
    if(wgt>1.||wgt<0.) {
      ostringstream wstring;
      wstring << "DISBase::generateCompton() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
    if(xT<=xTMin) break;
    if(UseRandom::rnd()>wgt) {
      continue;
    }
    break;
  }
  if(xT<=xTMin) {
    pTCompton_=-GeV;
    return;
  }
  comptonRawXP_ = xp;
  comptonRawZP_ = zp;
  // generate phi
  unsigned int itry(0);
  double phimax = azimuthalKernelMaximum(azicoeff.c0, azicoeff.c1, azicoeff.c2);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azimuthalKernelValue(azicoeff.c0, azicoeff.c1, azicoeff.c2, cphi);
    if (phiwgt > phimax * (1.0 + 1e-12)) {
      ostringstream wstring;
      wstring << "DISBase::generateCompton() "
              << "Azimuthal envelope undershot kernel value"
              << " phimax = " << phimax
              << " phiwgt = " << phiwgt
              << " cphi = " << cphi << "\n";
      generator()->logWarning(Exception(wstring.str(), Exception::warning));
    }
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::generateCompton() to"
				  << " generate phi" << Exception::eventerror;
  // momenta for the configuration
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTCompton_ = 0.5*Q*xT;
  ComptonMomenta_.resize(3);
  ComptonMomenta_[0] = p0;
  ComptonMomenta_[1] = p1;
  ComptonMomenta_[2] = p2;
  ComptonISFS_ = zp>xp;
}

// Generate a candidate BGF POWHEG emission and cache its Breit-frame momenta.
// The lower pT bound starts from the accepted Compton candidate so the final
// winner is the hardest of the two channels.
void DISBase::generateBGF() {
  bgfRawXP_ = 0.0;
  bgfRawZP_ = 0.0;
  // maximum value of the xT
  double xT = (1.-xB_)/xB_;
  double xTMin = 2.*max(pTmin_,pTCompton_)/sqrt(q2_);
  double zp;
  Energy2 alphaRefScale = powhegEmissionAlphaSScale(q2_,xTMin);
  // prefactor
  double a = powhegEmissionAlphaSOverestimate(alphaRefScale)*BGFWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  AzimuthalKernelCoefficients azicoeff;
  while(true) {
    wgt = 0.;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // zp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    // check allowed
    if(xp<xB_||xp>1.) {
      if(xT<=xTMin) break;
      continue;
    }
    // phase-space piece of the weight
    wgt = 8.*sqr(1.-xp)*zp/BGFWeight_;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    Energy2 pdfScale = powhegEmissionPDFScale(q2_, scale);
    double pdfRatioGen = 0.0;
    if (!powhegEmissionPDFRatioForGeneration(pdf_, beam_,
                                             gluon_, partons_[0],
                                             pdfScale, q2_,
                                             xB_/xp, xB_,
                                             pdfRatioGen)) {
      if(xT<=xTMin) break;
      continue;
    }
    wgt *= pdfRatioGen;
    // me piece of the weight
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    azicoeff = BGFME(xp,x2,x3,xT,false);
    Energy2 alphaScale = powhegEmissionAlphaSScale(q2_,xT);
    wgt *= 0.5*powhegEmissionAlphaSRatio(alphaScale,alphaRefScale)
      * azicoeff.average();
    if(wgt>1.||wgt<0.) {
      ostringstream wstring;
      wstring << "DISBase::generateBGF() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
    if(xT<=xTMin) break;
    if(UseRandom::rnd()>wgt) {
      continue;
    }
    break;
  }
  if(xT<=xTMin) {
    pTBGF_=-GeV;
    return;
  }
  bgfRawXP_ = xp;
  bgfRawZP_ = zp;
  // generate phi
  unsigned int itry(0);
  double phimax = azimuthalKernelMaximum(azicoeff.c0, azicoeff.c1, azicoeff.c2);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azimuthalKernelValue(azicoeff.c0, azicoeff.c1, azicoeff.c2, cphi);
    if (phiwgt > phimax * (1.0 + 1e-12)) {
      ostringstream wstring;
      wstring << "DISBase::generateBGF() "
              << "Azimuthal envelope undershot kernel value"
              << " phimax = " << phimax
              << " phiwgt = " << phiwgt
              << " cphi = " << cphi << "\n";
      generator()->logWarning(Exception(wstring.str(), Exception::warning));
    }
    ++itry;
  }
  while (phimax*UseRandom::rnd() > phiwgt && itry<200);
  if(itry==200) throw Exception() << "Too many tries in DISMECorrection"
				  << "::generateBGF() to"
				  << " generate phi" << Exception::eventerror;
  // momenta for the configuration
  Energy Q(sqrt(q2_));
  double x1 = -1./xp;
  double x2 = 1.-(1.-zp)/xp;
  double x3 = 2.+x1-x2;
  Lorentz5Momentum p1( 0.5*Q*xT*cos(phi),  0.5*Q*xT*sin(phi),
		       -0.5*Q*x2, 0.5*Q*sqrt(sqr(xT)+sqr(x2)));
  Lorentz5Momentum p2(-0.5*Q*xT*cos(phi), -0.5*Q*xT*sin(phi),
		      -0.5*Q*x3, 0.5*Q*sqrt(sqr(xT)+sqr(x3)));
  Lorentz5Momentum p0(ZERO,ZERO,-0.5*Q*x1,-0.5*Q*x1);
  pTBGF_=0.5*Q*xT;
  BGFMomenta_.resize(3);
  BGFMomenta_[0]=p0;
  BGFMomenta_[1]=p1;
  BGFMomenta_[2]=p2;
}

// Add one integration dimension for the NLO xp variable when POS/NEG NLO
// contributions are being generated.
int DISBase::nDim() const {
  return HwMEBase::nDim() + (contrib_>0 ? 1 : 0 );
}

// Generate the Born DIS kinematics and, for NLO contributions, sample xp with
// the production power-law mapping used by NLOWeightRaw().
bool DISBase::generateKinematics(const double * r) {
  xpSamplingRandom_ = 0.0;
  xpSamplingRho_ = 0.0;
  xpSamplingRhomin_ = 0.0;
  hadron_ = tcBeamPtr();
  beam_ = tcBeamPtr();
  pdf_ = tcPDFPtr();
  xB_ = 0.0;
  q2_ = ZERO;

  HwMEBase::TwoToTwoKinematicsSetup setup;
  if(!setupTwoToTwoKinematics(r, setup)) return false;

  // Intentionally determine the hadron side before sampling cos(theta) so the
  // native DIS window can tighten the Born angular range from the fixed xB of
  // the current incoming-parton configuration.
  tcBeamPtr bornHadron;
  double bornXB = 0.0;
  const bool haveBornHadron = determineBornHadronAndXB(bornHadron, bornXB);

  const double legacyCtmin = setup.ctmin;
  const double legacyCtmax = setup.ctmax;
  double nativeCtmin = legacyCtmin;
  double nativeCtmax = legacyCtmax;

  if (useNativeDISWindowGeneration_ && resolveNativeDISWindow()) {
    if (!haveBornHadron) return false;
    if (!tightenBornCosThetaWithNativeDISWindow(setup, bornXB,
                                                nativeCtmin, nativeCtmax)) {
      return false;
    }

  }

  const double cth = getCosTheta(nativeCtmin, nativeCtmax, r[0]);
  bool cutsRejected = false;
  if(!finishTwoToTwoKinematics(setup, cth, &cutsRejected)) {
    (void) cutsRejected;
    return false;
  }

  q2_ = -(meMomenta()[0]-meMomenta()[2]).m2();
  if (haveBornHadron) {
    hadron_ = bornHadron;
    xB_ = bornXB;
  }
  beam_ = hadron_;
  if (beam_) pdf_ = beam_->pdf();
  if(contrib_!=0) {
    // xp
    int ndim=nDim();
    double rhomin = pow(1.-xB_,1.-power_); 
    double rho = r[ndim-1]*rhomin;
    xpSamplingRandom_ = r[ndim-1];
    xpSamplingRho_ = rho;
    xpSamplingRhomin_ = rhomin;
    xp_ = 1.-pow(rho,1./(1.-power_));
    jac_ = rhomin/(1.-power_)*pow(1.-xp_,power_);
    jacobian(jacobian()*jac_);
  }
  return true; 
}

// Factorization/renormalization scale used by the Born and NLO weights.
Energy2 DISBase::scale() const {
  return scaleOpt_ == 1 ? 
    -sqr(scaleFact_)*tHat() : sqr(scaleFact_*muF_);
}

// Return the signed-order differential cross section. POS/NEG NLO samples use
// the clipped NLO/Born ratio from NLOWeight().
CrossSection DISBase::dSigHatDR() const {
  const CrossSection sigmaHat = HwMEBase::dSigHatDR();
  return contrib_ == 0 ? sigmaHat : NLOWeight()*sigmaHat;
}

// Assemble the validated NLO/Born ratio for arbitrary longitudinal beam
// polarizations. The common form keeps unpolarized PDFs in the flux and uses
// polarized PDFs only through effective parton polarizations and delta ratios.
double DISBase::NLOWeightRaw(double overrideLeptonPolarization,
                             double overrideHadronPolarization,
                             bool overridePolarizations) const {

  // If only leading order is required return 1:
  if (contrib_ == 0) return 1.0;

  // scale and prefactors
  Energy2 mu2(scale());
  double aS = SM().alphaS(mu2);
  double CFfact = 4./3.*aS/Constants::twopi;
  double TRfact = 1./2.*aS/Constants::twopi;

  // LO + dipole subtracted virtual + collinear quark bit with LO pdf
  double virt = 1. + CFfact * (-4.5 - 1./3.*sqr(Constants::pi)
                              + 1.5*log(q2_/mu2/(1.-xB_))
                              + 2.*log(1.-xB_)*log(q2_/mu2)
                              + sqr(log(1.-xB_)));
  virt /= jac_;

  // The lepton rho matrix is reliable at NLO; rebuild the hadron-side
  // polarization below from Pz*Delta q(xB)/q(xB).
  double Pl = 0.0, Pq = 0.0;
  if (overridePolarizations) {
    Pl = overrideLeptonPolarization;
  }
  else {
    const std::pair<RhoDMatrix,RhoDMatrix> rho = getRhoMatrices();
    Pl = longPol(rho.first);
  }


  // --------------------------------------------------------------------------
  // PDFs and effective polarisations
  //
  // For the *physical cross section with polarised beams* σ(Pℓ,Pp), keep the
  // unpolarised PDFs (sum PDFs) in the flux factors. Use the difference PDF Δf
  // only to build Δf/f polarisations that enter analysing powers / density matrices.
  // --------------------------------------------------------------------------

  // Unpolarised (sum) PDF: f(x,Q^2) = f^+(x,Q^2) + f^-(x,Q^2)
  ThePEG::tcPDFPtr sumPdf = hadron_->pdf();
  ThePEG::tcPDFPtr extractorSumPdf;
  {
    const auto & pbins = lastXComb().partonBinInstances();
    if (pbins.second) extractorSumPdf = pbins.second->pdf();
  }

  // Longitudinal hadron beam polarisation (convention: false = proton beam)
  const double Pz =
    overridePolarizations ? overrideHadronPolarization :
    getBeamPolarization(false).z();

  ThePEG::tcPDFPtr diffPdf;
  {
    // Use the same extractor path as getBeamPolarization().
    // eventHandler()->partonExtractor() is not guaranteed to be the
    // polarized extractor in these NLO runs.
    ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr ppe =
      ThePEG::dynamic_ptr_cast<
        ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr>(lastExtractor());

    if (ppe) {
      // proton is treated as the "second" beam
      diffPdf = ppe->longitudinalDifferencePDF().second;
    }
  }
  

  
  // --- Unpolarised PDFs for flux ratios ---
  double loPDF = sumPdf->xfx(hadron_, mePartonData()[1], mu2, xB_) / xB_;
  if (loPDF == 0.0) return 0.0;

  tcPDPtr gluon = getParticleData(ParticleID::g);

  // PDFs evaluated at the mapped momentum fraction x = xB/xp_
  double gPDF = sumPdf->xfx(hadron_, gluon,            mu2, xB_/xp_) * xp_ / xB_;
  double qPDF = sumPdf->xfx(hadron_, mePartonData()[1], mu2, xB_/xp_) * xp_ / xB_;

  // --- Difference (polarised) PDFs ---
  // Needed for collinear counterterms and real-emission PDF ratios.
  // dqPDF = Δq(xB/xp), dgPDF = Δg(xB/xp), dloPDF = Δq(xB)
  double dqPDF = 0.0, dgPDF = 0.0, dloPDF = 0.0;
  bool hasDiffPdf = (diffPdf && std::abs(Pz) > 1e-12);

  if (hasDiffPdf) {
    dqPDF  = diffPdf->xfx(hadron_, mePartonData()[1], mu2, xB_/xp_) * xp_ / xB_;
    dgPDF  = diffPdf->xfx(hadron_, gluon,             mu2, xB_/xp_) * xp_ / xB_;
    dloPDF = diffPdf->xfx(hadron_, mePartonData()[1], mu2, xB_)      / xB_;
  }

  // Born-side parton polarization at xB, matched to the corrected rho matrix
  // used by the DIS Born ME.
  if (hasDiffPdf && std::abs(loPDF) > 1e-30) {
    Pq = Pz * dloPDF / loPDF;
    Pq = std::max(-1.0, std::min(1.0, Pq));
  }

  // Protect polarized ratios against tiny denominators near PDF nodes.
  const double ratioFloor = 1e-12;
  const double minDlo = std::max(ratioFloor, 1e-4 * std::abs(loPDF));
  const bool hasStableRawFiniteDelta =
    hasDiffPdf && std::abs(loPDF) > ratioFloor;
  const double deltaqOverLo_raw =
    hasStableRawFiniteDelta ? Pz * dqPDF / loPDF : 0.0;
  const double deltagOverLo_raw =
    hasStableRawFiniteDelta ? Pz * dgPDF / loPDF : 0.0;
  const double deltaqOverLo = deltaqOverLo_raw;
  const double deltagOverLo = deltagOverLo_raw;

  // --- Effective parton polarisations at x = xB/xp_ ---
  double Pq_m = Pq;
  double Pg_m = 0.0;

  if (hasDiffPdf) {
    if (std::abs(qPDF) > ratioFloor) Pq_m = Pz * (dqPDF / qPDF);
    else Pq_m = 0.0;

    if (std::abs(gPDF) > ratioFloor) Pg_m = Pz * (dgPDF / gPDF);
    else Pg_m = 0.0;
  }

  Pq_m = std::max(-1.0, std::min(1.0, Pq_m));
  Pg_m = std::max(-1.0, std::min(1.0, Pg_m));

  // Canonical ratios used by the uniform polarized-NLO representation.
  const double qRatio = qPDF / loPDF;
  const double gRatio = gPDF / loPDF;
  struct UniformPolarizedNLOInputs {
    double qRatio;
    double gRatio;
    double deltaqOverLo;
    double deltagOverLo;
    double Pq;
    double Pq_m;
    double Pg_m;
    double ratioFloor;
    double minDlo;
    bool hasDiffPdf;
    bool hasStableDiffRatio;
  } uniformInputs = {
    qRatio,
    gRatio,
    deltaqOverLo,
    deltagOverLo,
    Pq,
    Pq_m,
    Pg_m,
    ratioFloor,
    minDlo,
    hasDiffPdf,
    hasDiffPdf && std::abs(dloPDF) > minDlo
  };

  // Calculate lepton kinematic variables (needed for the Born projectors and
  // real emission).
  Lorentz5Momentum q = meMomenta()[0]-meMomenta()[2];
  double  yB = (q*meMomenta()[1])/(meMomenta()[0]*meMomenta()[1]);
  double l = 2./yB-1.;

  // Born analysing power at xB using the PDF-based hadron polarization Pq.
  double a_born = A_pol(mePartonData()[0], mePartonData()[2],
                        mePartonData()[1], mePartonData()[3],
                        q2_, Pl, Pq);
  const CollinearBlendWeights blend =
    collinearBlendWeights(mePartonData()[0], mePartonData()[2],
                          mePartonData()[1], mePartonData()[3],
                          q2_, Pl, Pq, l);
  NeutralCurrentResponse ncResponse;
  const bool hasNCResponse =
    neutralCurrentResponse(mePartonData()[0], mePartonData()[2],
                           mePartonData()[1], mePartonData()[3],
                           q2_, Pl, Pq, Pq_m, l, ncResponse);
  const double qOddResponse = hasNCResponse ? ncResponse.qOddResponse : 0.0;
  const double gOddResponse = hasNCResponse ? ncResponse.gOddResponse : 0.0;
  const double qOddWeight =
    hasNCResponse ? qOddResponse :
    ((std::abs(Pq) > ratioFloor) ? (blend.qPolarized / Pq) : 0.0);
  const double gOddWeight =
    hasNCResponse ? gOddResponse :
    ((std::abs(Pq) > ratioFloor) ? (blend.gPolarized / Pq) : 0.0);

  // Collinear counterterms:
  // The quark and gluon channels use projector weights derived from the Born
  // angular structure. In the photon limit these reduce to the old scalar
  // f_pol blend. For full NC the gluon split is still exact for massless light
  // quarks: the parity-odd spin-even term cancels between the charge-conjugate
  // R2/R3 channels, and only the parity-even spin-independent and
  // parity-odd spin-dependent gluon structures survive.

  double logRatio = log((1.-xp_)*q2_/xp_/mu2);

  // --- gluon collinear ---
  // These are divided by loPDF later, so we compute collX/loPDF for each part.
  //
  // unpolarised piece / loPDF:
  //   TR/xp * (gPDF/loPDF) * [2xp(1-xp) + (xp²+(1-xp)²)*log(...)]
  double collg_over_born_unpol =
    TRfact/xp_*(gPDF/loPDF)*(2.*xp_*(1.-xp_) + (sqr(xp_)+sqr(1.-xp_))*logRatio);

  double dqRatio = 0.0;
  double dgRatio = 0.0;
  if (uniformInputs.hasStableDiffRatio) {
    dqRatio = dqPDF / dloPDF;
    dgRatio = dgPDF / dloPDF;
  }
  const double deltaqOverLo_eff =
    uniformInputs.hasStableDiffRatio ? Pq * dqRatio : 0.0;
  const double deltagOverLo_eff =
    uniformInputs.hasStableDiffRatio ? Pq * dgRatio : 0.0;

  const double collg_even =
    blend.gUnpolarized * collg_over_born_unpol;
  const double collg_odd =
    (hasNCResponse || uniformInputs.hasStableDiffRatio)
    ? gOddWeight * TRfact/xp_ * deltagOverLo_eff *
      (2.*(1.-xp_) + (2.*xp_-1.)*logRatio)
    : 0.0;

  // --- quark collinear ---
  // Same Pqq kernel, but different PDF ratios for polarised part.
  // unpolarised piece / loPDF:
  const double collq_k1 =
    1.-xp_-2./(1.-xp_)*log(xp_)
    -(1.+xp_)*log((1.-xp_)/xp_*q2_/mu2);
  const double collq_k2 =
    2./(1.-xp_)*log(q2_*(1.-xp_)/mu2)-1.5/(1.-xp_);
  double collq_over_born_unpol =
    CFfact/xp_*uniformInputs.qRatio*collq_k1
    + CFfact/xp_*(uniformInputs.qRatio-xp_)*collq_k2;

  const double collq_even =
    blend.qUnpolarized * collq_over_born_unpol;
  const double collq_odd =
    uniformInputs.hasStableDiffRatio
    ? qOddWeight *
      (CFfact/xp_ * deltaqOverLo_eff * collq_k1
       + CFfact/xp_ * (deltaqOverLo_eff - uniformInputs.Pq * xp_) * collq_k2)
    : 0.0;

  // Real-emission kernels are written as sigma_real / sigma_Born.
  // Keep a_born in the denominator so the Born factor cancels exactly, and use
  // mapped analyzing powers at x = xB/xp_ for the spin dependence.

  // Mapped analysing powers for the real-emission kernels
  double a_q_mapped = A_pol(mePartonData()[0], mePartonData()[2],
                            mePartonData()[1], mePartonData()[3],
                            q2_, Pl, Pq_m);

  // q -> qg term (QCDC): denominator a_born, kernel a_q_mapped
  const double qcdcDenRatio =
    qcdcMappedDenominatorRatio(mePartonData()[0], mePartonData()[2],
                               mePartonData()[1], mePartonData()[3],
                               q2_, Pl, Pq, Pq_m);
  const double bornFactor = 1. + a_born*l + sqr(l);
  const double realq_prefactor =
    qcdcDenRatio * CFfact/xp_/bornFactor * uniformInputs.qRatio;
  const double realq_even =
    realq_prefactor * (2.+2.*sqr(l)-xp_+3.*xp_*sqr(l));
  const double realq_odd =
    realq_prefactor * a_q_mapped*l*(2.*xp_+1.);
  double realq = realq_even + realq_odd;

  // g -> q qbar term (BGF): use the same exact gluon projector split as the
  // collinear counterterm. For massless light quarks the integrated F3-like
  // gluon term cancels between the charge-conjugate R2/R3 channels, so the
  // finite remainder keeps only the F2/FL-like unpolarized channel and the
  // G2-like spin-dependent channel.
  const double realg_over_born_unpol =
    -TRfact/xp_ * gRatio *
    ((1.+sqr(l)+2.*(1.-3.*sqr(l))*xp_*(1.-xp_)) / (1.+sqr(l)));

  const double realg_even =
    blend.gUnpolarized * realg_over_born_unpol;
  const double realg_odd =
    (hasNCResponse || uniformInputs.hasStableDiffRatio)
    ? gOddWeight * (-1.0 * TRfact/xp_ * deltagOverLo_eff * (2.*xp_ - 1.))
    : 0.0;
  const double collq_over_born = collq_even + collq_odd;
  const double collg_over_born = collg_even + collg_odd;
  const double realg = realg_even + realg_odd;

  // Full NLO/Born ratio for the validated uniform polarized-NLO assembly.
  double wgt = virt + collq_over_born + collg_over_born + realq + realg;

  return wgt;
}

// Convert the signed NLO/Born ratio into the positive or negative event stream
// selected by the Contribution interface.
double DISBase::NLOWeight() const {
  const double wgt = NLOWeightRaw(0.0, 0.0, false);
  return contrib_ == 1 ? max(0., wgt) : max(0., -wgt);
}

// ThePEG hook for optional HepMC weights; DIS currently provides the correlated
// Rivet helicity-weight set when enabled.
std::map<std::string,double> DISBase::generateOptionalWeights() {
  return generateRivetWeights();
}

// Build the Born-side parton polarization for a requested hadron polarization
// using the same polarized-PDF ratio as the NLO weight.
double DISBase::rivetWeightBornPartonPolarization(double hadronPolarization) const {
  if (!hadron_ || mePartonData().size() < 2 || xB_ <= 0.0) return 0.0;
  if (std::abs(hadronPolarization) <= 1e-12) return 0.0;

  ThePEG::tcPDFPtr sumPdf = hadron_->pdf();
  if (!sumPdf) return 0.0;

  ThePEG::tcPDFPtr diffPdf;
  ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr ppe =
    ThePEG::dynamic_ptr_cast<
      ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr>(lastExtractor());
  if (ppe) diffPdf = ppe->longitudinalDifferencePDF().second;
  if (!diffPdf) return 0.0;

  const Energy2 mu2(scale());
  const double loPDF =
    sumPdf->xfx(hadron_, mePartonData()[1], mu2, xB_) / xB_;
  if (!std::isfinite(loPDF) || std::abs(loPDF) <= 1e-30) return 0.0;

  const double dloPDF =
    diffPdf->xfx(hadron_, mePartonData()[1], mu2, xB_) / xB_;
  if (!std::isfinite(dloPDF)) return 0.0;

  const double polarization = hadronPolarization * dloPDF / loPDF;
  if (!std::isfinite(polarization)) return 0.0;
  return std::max(-1.0, std::min(1.0, polarization));
}

// Default Born matrix element for correlated weights. Neutral/charged current
// classes override when explicit helicity polarizations are required.
double DISBase::rivetWeightBornME2(double, double) const {
  return me2();
}

// Produce correlated PP/PM/MP/MM/SIGMA/DELTA_LL weights on the sampled event so
// Rivet can fill helicity combinations without independent statistical noise.
std::map<std::string,double> DISBase::generateRivetWeights() const {
  std::map<std::string,double> weights;
  if (!generateRivetWeights_) return weights;

  weights["HERWIG_DIS_PP"] = 0.0;
  weights["HERWIG_DIS_PM"] = 0.0;
  weights["HERWIG_DIS_MP"] = 0.0;
  weights["HERWIG_DIS_MM"] = 0.0;
  weights["HERWIG_DIS_SIGMA"] = 0.0;
  weights["HERWIG_DIS_DELTA_LL"] = 0.0;

  const tcEventPtr event = generator()->currentEvent();
  const double eventWeight = event ? event->weight() : 1.0;
  if (!std::isfinite(eventWeight)) return weights;

  const double sampledPq = rivetWeightBornPartonPolarization(0.0);
  const double sampledBorn = rivetWeightBornME2(0.0, sampledPq);
  const double sampledRaw = NLOWeightRaw(0.0, 0.0, true);
  const double sampledDensity = sampledBorn * sampledRaw;
  const double sampledAbs = std::abs(sampledDensity);
  if (!std::isfinite(sampledAbs) || sampledAbs <= 1e-30) {
    return weights;
  }
  const double sampledDenominator = sampledAbs;
  if (!std::isfinite(sampledDenominator) || sampledDenominator <= 1e-30) return weights;

  struct HelicityTarget {
    const char * name;
    double leptonPol;
    double hadronPol;
  };
  static const HelicityTarget targets[] = {
    { "HERWIG_DIS_PP",  1.0,  1.0 },
    { "HERWIG_DIS_PM",  1.0, -1.0 },
    { "HERWIG_DIS_MP", -1.0,  1.0 },
    { "HERWIG_DIS_MM", -1.0, -1.0 }
  };

  for (const HelicityTarget & target : targets) {
    const double partonPol =
      rivetWeightBornPartonPolarization(target.hadronPol);
    const double born = rivetWeightBornME2(target.leptonPol, partonPol);
    const double raw = NLOWeightRaw(target.leptonPol, target.hadronPol, true);
    const double targetWeight = eventWeight * born * raw / sampledDenominator;
    weights[target.name] = std::isfinite(targetWeight) ? targetWeight : 0.0;
  }

  const double pp = weights["HERWIG_DIS_PP"];
  const double pm = weights["HERWIG_DIS_PM"];
  const double mp = weights["HERWIG_DIS_MP"];
  const double mm = weights["HERWIG_DIS_MM"];
  weights["HERWIG_DIS_SIGMA"] = 0.25 * (pp + pm + mp + mm);
  weights["HERWIG_DIS_DELTA_LL"] = 0.25 * (pp + mm - pm - mp);
  return weights;
}
