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
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.h"
#include "ThePEG/PDF/LHAPDF6.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include <numeric>
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <set>
#include "Herwig/Shower/RealEmissionProcess.h"
#include "ThePEG/PDF/PolarizedPartonExtractor.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

// namespace {
// using namespace Herwig;
// using namespace ThePEG::Helicity;
// 
// void debuggingMatrixElement(bool BGF,const Lorentz5Momentum & pin,
// 			    const Lorentz5Momentum & p1,
// 			    const Lorentz5Momentum & p2,
// 			    tcPDPtr gluon,
// 			    const Lorentz5Momentum & pl1,
// 			    const Lorentz5Momentum & pl2,
// 			    const Lorentz5Momentum & pq1,
// 			    const Lorentz5Momentum & pq2,
// 			    tcPDPtr lepton1,tcPDPtr lepton2,
// 			    tcPDPtr quark1 ,tcPDPtr quark2,
// 			    Energy2 Q2,double phi, double x2, double x3, 
// 			    double xperp, double zp, double xp,
// 			    const vector<double> & azicoeff, 
// 			    bool normalize) {
//   tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>
//     (CurrentGenerator::current().standardModel());
//   assert(hwsm);
//   vector<AbstractFFVVertexPtr> weakVertex;
//   vector<PDPtr> bosons;
//   AbstractFFVVertexPtr strongVertex = hwsm->vertexFFG();
//   if(lepton1->id()==lepton2->id()) {
//     weakVertex.push_back(hwsm->vertexFFZ());
//     bosons.push_back(hwsm->getParticleData(ParticleID::Z0));
//     weakVertex.push_back(hwsm->vertexFFP());
//     bosons.push_back(hwsm->getParticleData(ParticleID::gamma));
//   }
//   else {
//     weakVertex.push_back(hwsm->vertexFFW());
//     bosons.push_back(hwsm->getParticleData(ParticleID::Wplus));
//   }
//   if(!BGF) {
//     SpinorWaveFunction    l1,q1,qp1;
//     SpinorBarWaveFunction l2,q2,qp2;
//     VectorWaveFunction    gl(p2,gluon,outgoing);
//     if(lepton1->id()>0) {
//       l1  = SpinorWaveFunction   (pl1,lepton1,incoming);
//       l2  = SpinorBarWaveFunction(pl2,lepton2,outgoing);
//     }
//     else {
//       l1  = SpinorWaveFunction   (pl2,lepton2,outgoing);
//       l2  = SpinorBarWaveFunction(pl1,lepton1,incoming);
//     }
//     if(quark1->id()>0) {
//       q1  = SpinorWaveFunction   (pq1,quark1,incoming);
//       q2  = SpinorBarWaveFunction(pq2,quark2,outgoing);
//       qp1 = SpinorWaveFunction   (pin,quark1,incoming);
//       qp2 = SpinorBarWaveFunction(p1 ,quark2,outgoing);
//     }
//     else {
//       q1  = SpinorWaveFunction   (pq2,quark2,outgoing);
//       q2  = SpinorBarWaveFunction(pq1,quark1,incoming);
//       qp1 = SpinorWaveFunction   (p1 ,quark2,outgoing);
//       qp2 = SpinorBarWaveFunction(pin,quark1,incoming);
//     }
//     double lome(0.),realme(0.);
//     for(unsigned int lhel1=0;lhel1<2;++lhel1) {
//       l1.reset(lhel1);
//       for(unsigned int lhel2=0;lhel2<2;++lhel2) { 
// 	l2.reset(lhel2);
// 	for(unsigned int qhel1=0;qhel1<2;++qhel1) {
// 	  q1.reset(qhel1);
// 	  qp1.reset(qhel1);
// 	  for(unsigned int qhel2=0;qhel2<2;++qhel2) {
// 	    q2.reset(qhel2);
// 	    qp2.reset(qhel2);
// 	    // leading order matrix element 
// 	    Complex diagLO(0.);
// 	    for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 	      VectorWaveFunction inter = 
// 		weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 	      diagLO += weakVertex[ix]->evaluate(Q2,q1,q2,inter);
// 	    }
// 	    lome   += norm(diagLO);
// 	    // real emission matrix element
// 	    for(unsigned int ghel=0;ghel<2;++ghel) {
//  	      gl.reset(2*ghel);
// 	      Complex diagReal(0.);
// 	      for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 		VectorWaveFunction inter = 
// 		  weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 		SpinorWaveFunction off1 = 
// 		  strongVertex->evaluate(Q2,5,qp1.particle(),qp1,gl);
// 		Complex diag1 = weakVertex[ix]->evaluate(Q2,off1,qp2,inter);
// 		SpinorBarWaveFunction off2 = 
// 		  strongVertex->evaluate(Q2,5,qp2.particle(),qp2,gl);
// 		Complex diag2 = weakVertex[ix]->evaluate(Q2,qp1,off2,inter);
// 		diagReal += diag1+diag2;
// 	      }
// 	      realme += norm(diagReal);
// 	    }
// 	  }
// 	}
//       }
//     }
//     double test1 = realme/lome/hwsm->alphaS(Q2)*Q2*UnitRemoval::InvE2;
//     double cphi(cos(phi));
//     double test2;
//     if(normalize) {
//       test2 = 8.*Constants::pi/(1.-xp)/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi))*
// 	(1.+sqr(xp)*(sqr(x2)+1.5*sqr(xperp)));
//     }
//     else {
//       test2 = 8.*Constants::pi/(1.-xp)/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi));
//     }
//     cerr << "testing RATIO A  " << test1/test2 << "\n";
//   }
//   else {
//     SpinorWaveFunction    l1,q1,qp1;
//     SpinorBarWaveFunction l2,q2,qp2;
//     VectorWaveFunction    gl(pin,gluon,incoming);
//     if(lepton1->id()>0) {
//       l1  = SpinorWaveFunction   (pl1,lepton1,incoming);
//       l2  = SpinorBarWaveFunction(pl2,lepton2,outgoing);
//     }
//     else {
//       l1  = SpinorWaveFunction   (pl2,lepton2,outgoing);
//       l2  = SpinorBarWaveFunction(pl1,lepton1,incoming);
//     }
//     if(quark1->id()>0) {
//       q1  = SpinorWaveFunction   (pq1,quark1      ,incoming);
//       q2  = SpinorBarWaveFunction(pq2,quark2      ,outgoing);
//       qp2 = SpinorBarWaveFunction(p1    ,quark2      ,outgoing);
//       qp1 = SpinorWaveFunction   (p2    ,quark1->CC(),outgoing);
//     }
//     else {
//       q1  = SpinorWaveFunction   (pq2,quark2      ,outgoing);
//       q2  = SpinorBarWaveFunction(pq1,quark1      ,incoming);
//       qp2 = SpinorBarWaveFunction(p2    ,quark1->CC(),outgoing);
//       qp1 = SpinorWaveFunction   (p1    ,quark2      ,outgoing);
//     }
//     double lome(0.),realme(0.);
//     for(unsigned int lhel1=0;lhel1<2;++lhel1) {
//       l1.reset(lhel1);
//       for(unsigned int lhel2=0;lhel2<2;++lhel2) { 
// 	l2.reset(lhel2);
// 	for(unsigned int qhel1=0;qhel1<2;++qhel1) {
// 	  q1.reset(qhel1);
// 	  qp1.reset(qhel1);
// 	  for(unsigned int qhel2=0;qhel2<2;++qhel2) {
// 	    q2.reset(qhel2);
// 	    qp2.reset(qhel2);
// 	    // leading order matrix element 
// 	    Complex diagLO(0.);
// 	    for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 	      VectorWaveFunction inter = 
// 		weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 	      diagLO += weakVertex[ix]->evaluate(Q2,q1,q2,inter);
// 	    }
// 	    lome   += norm(diagLO);
// 	    // real emission matrix element
// 	    for(unsigned int ghel=0;ghel<2;++ghel) {
//   	      gl.reset(2*ghel);
// 	      Complex diagReal(0.);
// 	      for(unsigned int ix=0;ix<weakVertex.size();++ix) {
// 		VectorWaveFunction inter = 
// 		  weakVertex[ix]->evaluate(Q2,3,bosons[ix],l1,l2);
// 		SpinorWaveFunction off1 = 
// 		  strongVertex->evaluate(Q2,5,qp1.particle(),qp1,gl);
// 		Complex diag1 = weakVertex[ix]->evaluate(Q2,off1,qp2,inter);
// 		SpinorBarWaveFunction off2 = 
// 		  strongVertex->evaluate(Q2,5,qp2.particle(),qp2,gl);
// 		Complex diag2 = weakVertex[ix]->evaluate(Q2,qp1,off2,inter);
// 		diagReal += diag1+diag2;
// 	      }
// 	      realme += norm(diagReal);
// 	    }
// 	  }
// 	}
//       }
//     }
//     double test1 = realme/lome/hwsm->alphaS(Q2)*Q2*UnitRemoval::InvE2;
//     double cphi(cos(phi));
//     double test2;
//     if(normalize) {
//       test2 = 8.*Constants::pi/zp/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi))*
// 	sqr(xp)*(sqr(x3)+sqr(x2)+3.*sqr(xperp));
//     }
//     else {
//       test2 = 8.*Constants::pi/zp/(1.-zp)*
// 	(azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi));
//     }
//     cerr << "testing RATIO B " << test1/test2 << "\n";
//   }
// }
// 
// }
 
DISBase::POWHEGRawChannelDiagnostics::POWHEGRawChannelDiagnostics()
  : status("below_ptmin"),
    trials(0), rejectXP(0), rejectVeto(0), weightNeg(0), weightHigh(0),
    xp(0.0), zp(0.0), xMapped(0.0), xT(0.0), xTMin(0.0), pT(0.0),
    phase(0.0), pdfRatio(0.0), alphaRatio(0.0), meAvg(0.0), wgt(0.0),
    pdfScale(0.0), alphaScale(0.0)
{}

DISBase::POWHEGRawEventDiagnostics::POWHEGRawEventDiagnostics()
  : q2(0.0), xB(0.0), y(0.0), s(0.0), winner("None"), fallback(0),
    compton(), bgf()
{}

DISBase::DISBase()  : initial_(6.), final_(3.),
		      procProb_(0.35),
		      comptonInt_(0.), bgfInt_(0.),
		      comptonWeight_(50.), BGFWeight_(150.), 
		      pTmin_(0.1*GeV), 
		      scaleOpt_(1),  muF_(100.*GeV), scaleFact_(1.),
		      contrib_(0),
		      useFixedOrderAlphaSInPOWHEGEmission_(false),
		      useQ2ScaleInPOWHEGEmission_(false),
		      disDiagnostics_(false),
		      dumpNLOTermDiagnostics_(false),
		      powhegEmissionComparisonMode_(POWHEGEmissionComparisonModeDefault),
		      powhegEmissionComparisonMaxAttempts_(100),
		      usePOWHEGRealSpinVertex_(false),
		      dumpPOWHEGRawMomenta_(false),
		      powhegRawMomentaDumpMax_(0),
		      comptonRawXP_(0.0), comptonRawZP_(0.0),
		      bgfRawXP_(0.0), bgfRawZP_(0.0),
		      powhegRawDiagnostics_(),
		      diagnosePOWHEGRealSpinVertex_(false),
		      powhegRealSpinDiagMax_(20), powhegRealSpinDiagCount_(0),
		      powhegRawMomentaDumpCount_(0),
		      power_(0.1)
{}

DISBase::~DISBase() {}

void DISBase::persistentOutput(PersistentOStream & os) const {
  os << comptonInt_ << bgfInt_ << procProb_ << initial_ << final_ << alpha_
     << ounit(pTmin_,GeV) << comptonWeight_ << BGFWeight_ << gluon_
     << ounit(muF_,GeV) << scaleFact_ << scaleOpt_ << contrib_
     << useFixedOrderAlphaSInPOWHEGEmission_
     << useQ2ScaleInPOWHEGEmission_
     << disDiagnostics_
     << dumpNLOTermDiagnostics_
     << powhegEmissionComparisonMode_
     << powhegEmissionComparisonMaxAttempts_
     << usePOWHEGRealSpinVertex_
     << dumpPOWHEGRawMomenta_ << powhegRawMomentaDumpMax_
     << diagnosePOWHEGRealSpinVertex_
     << powhegRealSpinDiagMax_ << power_;
}

void DISBase::persistentInput(PersistentIStream & is, int version) {
  is >> comptonInt_ >> bgfInt_ >> procProb_  >> initial_ >> final_ >> alpha_
     >> iunit(pTmin_,GeV) >> comptonWeight_ >> BGFWeight_ >> gluon_
     >> iunit(muF_,GeV) >> scaleFact_ >> scaleOpt_ >> contrib_;
  if(version == 0) {
    useFixedOrderAlphaSInPOWHEGEmission_ = false;
    useQ2ScaleInPOWHEGEmission_ = false;
    disDiagnostics_ = false;
    dumpNLOTermDiagnostics_ = false;
    powhegEmissionComparisonMode_ = POWHEGEmissionComparisonModeDefault;
    powhegEmissionComparisonMaxAttempts_ = 100;
    diagnosePOWHEGRealSpinVertex_ = false;
    powhegRealSpinDiagMax_ = 20;
    is >> usePOWHEGRealSpinVertex_ >> power_;
  }
  else if(version == 1) {
    bool legacyFixedOrderSampling = false;
    is >> useFixedOrderAlphaSInPOWHEGEmission_
       >> useQ2ScaleInPOWHEGEmission_
       >> legacyFixedOrderSampling;
    disDiagnostics_ = false;
    dumpNLOTermDiagnostics_ = false;
    powhegEmissionComparisonMode_ = POWHEGEmissionComparisonModeDefault;
    powhegEmissionComparisonMaxAttempts_ = 100;
    is >> usePOWHEGRealSpinVertex_ >> diagnosePOWHEGRealSpinVertex_
       >> powhegRealSpinDiagMax_ >> power_;
  }
  else if(version == 2) {
    bool legacyFixedOrderSampling = false;
    is >> useFixedOrderAlphaSInPOWHEGEmission_
       >> useQ2ScaleInPOWHEGEmission_
       >> legacyFixedOrderSampling
       ;
    disDiagnostics_ = false;
    dumpNLOTermDiagnostics_ = false;
    is >> powhegEmissionComparisonMode_
       >> powhegEmissionComparisonMaxAttempts_
       >> usePOWHEGRealSpinVertex_;
    dumpPOWHEGRawMomenta_ = false;
    powhegRawMomentaDumpMax_ = 0;
    is >> diagnosePOWHEGRealSpinVertex_
       >> powhegRealSpinDiagMax_ >> power_;
  }
  else if(version == 3) {
    disDiagnostics_ = false;
    dumpNLOTermDiagnostics_ = false;
    is >> useFixedOrderAlphaSInPOWHEGEmission_
       >> useQ2ScaleInPOWHEGEmission_
       >> powhegEmissionComparisonMode_
       >> powhegEmissionComparisonMaxAttempts_
       >> usePOWHEGRealSpinVertex_;
    dumpPOWHEGRawMomenta_ = false;
    powhegRawMomentaDumpMax_ = 0;
    is >> diagnosePOWHEGRealSpinVertex_
       >> powhegRealSpinDiagMax_ >> power_;
  }
  else if(version == 4) {
    disDiagnostics_ = false;
    dumpNLOTermDiagnostics_ = false;
    is >> useFixedOrderAlphaSInPOWHEGEmission_
       >> useQ2ScaleInPOWHEGEmission_
       >> powhegEmissionComparisonMode_
       >> powhegEmissionComparisonMaxAttempts_
       >> usePOWHEGRealSpinVertex_
       >> dumpPOWHEGRawMomenta_ >> powhegRawMomentaDumpMax_
       >> diagnosePOWHEGRealSpinVertex_
       >> powhegRealSpinDiagMax_ >> power_;
  }
  else if(version == 5) {
    dumpNLOTermDiagnostics_ = false;
    is >> useFixedOrderAlphaSInPOWHEGEmission_
       >> useQ2ScaleInPOWHEGEmission_
       >> disDiagnostics_
       >> powhegEmissionComparisonMode_
       >> powhegEmissionComparisonMaxAttempts_
       >> usePOWHEGRealSpinVertex_
       >> dumpPOWHEGRawMomenta_ >> powhegRawMomentaDumpMax_
       >> diagnosePOWHEGRealSpinVertex_
       >> powhegRealSpinDiagMax_ >> power_;
  }
  else {
    is >> useFixedOrderAlphaSInPOWHEGEmission_
       >> useQ2ScaleInPOWHEGEmission_
       >> disDiagnostics_
       >> dumpNLOTermDiagnostics_
       >> powhegEmissionComparisonMode_
       >> powhegEmissionComparisonMaxAttempts_
       >> usePOWHEGRealSpinVertex_
       >> dumpPOWHEGRawMomenta_ >> powhegRawMomentaDumpMax_
       >> diagnosePOWHEGRealSpinVertex_
       >> powhegRealSpinDiagMax_ >> power_;
  }
  if(powhegEmissionComparisonMode_ >
     POWHEGEmissionComparisonModeRealOnly) {
    powhegEmissionComparisonMode_ = POWHEGEmissionComparisonModeDefault;
  }
  if(powhegEmissionComparisonMaxAttempts_ == 0) {
    powhegEmissionComparisonMaxAttempts_ = 100;
  }
  powhegRealSpinDiagCount_ = 0;
  powhegRawMomentaDumpCount_ = 0;
  resetPOWHEGRawDiagnostics();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<DISBase,HwMEBase>
  describeHerwigDISBase("Herwig::DISBase", "HwMEDIS.so", 6);

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
  Pq = std::max(-1.0, std::min(1.0, Pq));

  // DIS always has a spin-1/2 incoming quark here.
  rho.second = RhoDMatrix(mePartonData()[1]->iSpin());
  const unsigned int imax = rho.second.iSpin() - 1;
  rho.second(0,0) = 0.5 * (1.0 - Pq);
  rho.second(imax,imax) = 0.5 * (1.0 + Pq);
  rho.second(0,imax) = 0.0;
  rho.second(imax,0) = 0.0;

  return rho;
}

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

double DISBase::qcdcMappedDenominatorRatio(tcPDPtr, tcPDPtr,
                                           tcPDPtr, tcPDPtr,
                                           Energy2,
                                           double,
                                           double,
                                           double) const {
  return 1.0;
}

bool DISBase::bornClosureDiagnostics(double,
                                     double,
                                     double,
                                     double,
                                     BornClosureDiagnostics &) const {
  return false;
}

void DISBase::constructRealEmissionSpinVertex(RealEmissionProcessPtr, bool) const {}

void DISBase::diagnoseRealEmissionSpinState(RealEmissionProcessPtr, bool) const {}

bool DISBase::nextPOWHEGRealSpinDiagnosticSlot(unsigned long & index) const {
  if (!disDiagnostics_) return false;
  if (!diagnosePOWHEGRealSpinVertex_) return false;
  if (powhegRealSpinDiagMax_ != 0 &&
      powhegRealSpinDiagCount_ >= powhegRealSpinDiagMax_) return false;
  ++powhegRealSpinDiagCount_;
  index = powhegRealSpinDiagCount_;
  return true;
}

bool DISBase::beginPOWHEGRawDiagnosticEvent(unsigned long & index) const {
  if(!disDiagnostics_) return false;
  if(!dumpPOWHEGRawMomenta_) return false;
  if(powhegRawMomentaDumpMax_ != 0 &&
     powhegRawMomentaDumpCount_ >= powhegRawMomentaDumpMax_) return false;
  ++powhegRawMomentaDumpCount_;
  index = powhegRawMomentaDumpCount_;
  return true;
}

void DISBase::resetPOWHEGRawDiagnostics() {
  powhegRawDiagnostics_ = POWHEGRawEventDiagnostics();
}

double DISBase::currentPOWHEGRawBeamS() const {
  if(lastParticles().first && lastParticles().second) {
    return (lastParticles().first->momentum() +
            lastParticles().second->momentum()).m2()/GeV2;
  }
  return 0.0;
}

void DISBase::dumpPOWHEGRawSummary(unsigned long eventIndex) const {
  if(!disDiagnostics_ || !dumpPOWHEGRawMomenta_) return;

  const POWHEGRawChannelDiagnostics & comp = powhegRawDiagnostics_.compton;
  const POWHEGRawChannelDiagnostics & bgf  = powhegRawDiagnostics_.bgf;

  std::ostringstream line;
  line << std::scientific << std::setprecision(16)
       << "POWHEG_RAW_SUMMARY"
       << " event=" << eventIndex
       << " winner=" << powhegRawDiagnostics_.winner
       << " fallback=" << powhegRawDiagnostics_.fallback
       << " Q2=" << powhegRawDiagnostics_.q2
       << " xB=" << powhegRawDiagnostics_.xB
       << " y=" << powhegRawDiagnostics_.y
       << " s=" << powhegRawDiagnostics_.s
       << " comp_status=" << comp.status
       << " comp_trials=" << comp.trials
       << " comp_rejectXP=" << comp.rejectXP
       << " comp_rejectVeto=" << comp.rejectVeto
       << " comp_weightNeg=" << comp.weightNeg
       << " comp_weightHigh=" << comp.weightHigh
       << " comp_xp=" << comp.xp
       << " comp_zp=" << comp.zp
       << " comp_xMapped=" << comp.xMapped
       << " comp_xT=" << comp.xT
       << " comp_xTMin=" << comp.xTMin
       << " comp_pT=" << comp.pT
       << " comp_phase=" << comp.phase
       << " comp_pdfRatio=" << comp.pdfRatio
       << " comp_alphaRatio=" << comp.alphaRatio
       << " comp_meAvg=" << comp.meAvg
       << " comp_wgt=" << comp.wgt
       << " comp_pdfScale=" << comp.pdfScale
       << " comp_alphaScale=" << comp.alphaScale
       << " bgf_status=" << bgf.status
       << " bgf_trials=" << bgf.trials
       << " bgf_rejectXP=" << bgf.rejectXP
       << " bgf_rejectVeto=" << bgf.rejectVeto
       << " bgf_weightNeg=" << bgf.weightNeg
       << " bgf_weightHigh=" << bgf.weightHigh
       << " bgf_xp=" << bgf.xp
       << " bgf_zp=" << bgf.zp
       << " bgf_xMapped=" << bgf.xMapped
       << " bgf_xT=" << bgf.xT
       << " bgf_xTMin=" << bgf.xTMin
       << " bgf_pT=" << bgf.pT
       << " bgf_phase=" << bgf.phase
       << " bgf_pdfRatio=" << bgf.pdfRatio
       << " bgf_alphaRatio=" << bgf.alphaRatio
       << " bgf_meAvg=" << bgf.meAvg
       << " bgf_wgt=" << bgf.wgt
       << " bgf_pdfScale=" << bgf.pdfScale
       << " bgf_alphaScale=" << bgf.alphaScale
       << "\n";
  generator()->log() << line.str();
}

double DISBase::mappedIncomingLongitudinalPolarization(tcPDPtr parton,
                                                       double x,
                                                       Energy2 mu2) const {
  if (!parton || !hadron_ || x <= 0.0 || x >= 1.0) return 0.0;

  const double Pz = getBeamPolarization(false).z();
  if (std::abs(Pz) <= 1e-12) return 0.0;

  ThePEG::tcPDFPtr sumPdf = hadron_->pdf();
  if (!sumPdf) return 0.0;

  ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr ppe =
    ThePEG::dynamic_ptr_cast<
      ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr>(lastExtractor());
  if (!ppe) return 0.0;

  ThePEG::tcPDFPtr diffPdf = ppe->longitudinalDifferencePDF().second;
  if (!diffPdf) return 0.0;

  const double sum = sumPdf->xfx(hadron_, parton, mu2, x) / x;
  if (std::abs(sum) <= 1e-30) return 0.0;

  const double diff = diffPdf->xfx(hadron_, parton, mu2, x) / x;
  double pol = Pz * diff / sum;
  return std::max(-1.0, std::min(1.0, pol));
}

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

  static Switch<DISBase,bool> interfaceDISDiagnostics
    ("DISDiagnostics",
     "Master switch for optional DIS/POWHEG diagnostic log output. Specific diagnostic switches stay silent unless this is On.",
     &DISBase::disDiagnostics_, false, false, false);
  static SwitchOption interfaceDISDiagnosticsNo
    (interfaceDISDiagnostics,
     "Off",
     "Keep DIS/POWHEG diagnostic logging disabled.",
     false);
  static SwitchOption interfaceDISDiagnosticsYes
    (interfaceDISDiagnostics,
     "On",
     "Enable DIS/POWHEG diagnostic logging controlled by the specific diagnostic switches.",
     true);

  static Switch<DISBase,bool> interfaceDumpNLOTermDiagnostics
    ("DumpNLOTermDiagnostics",
     "Emit the periodic NLO_TERM_* summary diagnostics. Requires DISDiagnostics to be On.",
     &DISBase::dumpNLOTermDiagnostics_, false, false, false);
  static SwitchOption interfaceDumpNLOTermDiagnosticsNo
    (interfaceDumpNLOTermDiagnostics,
     "No",
     "Do not emit the periodic NLO_TERM_* summary diagnostics.",
     false);
  static SwitchOption interfaceDumpNLOTermDiagnosticsYes
    (interfaceDumpNLOTermDiagnostics,
     "Yes",
     "Emit the periodic NLO_TERM_* summary diagnostics.",
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

  static Switch<DISBase,bool> interfaceDumpPOWHEGRawMomenta
    ("DumpPOWHEGRawMomenta",
     "Dump the winning raw POWHEG hadronic pair before qtilde reconstruction changes it. Requires DISDiagnostics to be On.",
     &DISBase::dumpPOWHEGRawMomenta_, false, false, false);
  static SwitchOption interfaceDumpPOWHEGRawMomentaNo
    (interfaceDumpPOWHEGRawMomenta,
     "No",
     "Do not emit raw POWHEG momentum diagnostics.",
     false);
  static SwitchOption interfaceDumpPOWHEGRawMomentaYes
    (interfaceDumpPOWHEGRawMomenta,
     "Yes",
     "Emit structured raw POWHEG momentum diagnostics for the winning hardest emission.",
     true);

  static Parameter<DISBase,unsigned long> interfacePOWHEGRawMomentaDumpMax
    ("POWHEGRawMomentaDumpMax",
     "Maximum number of raw POWHEG momentum diagnostics to emit. Zero means unlimited. Active only when DISDiagnostics is On.",
     &DISBase::powhegRawMomentaDumpMax_, 0, 0, 1000000000,
     false, false, Interface::lowerlim);

  static Switch<DISBase,bool> interfaceDiagnosePOWHEGRealSpinVertex
    ("DiagnosePOWHEGRealSpinVertex",
     "Emit a dedicated diagnostic for the realised POWHEG 2->3 spin state. Requires DISDiagnostics to be On.",
     &DISBase::diagnosePOWHEGRealSpinVertex_, false, false, false);
  static SwitchOption interfaceDiagnosePOWHEGRealSpinVertexNo
    (interfaceDiagnosePOWHEGRealSpinVertex,
     "No",
     "Do not emit the real-emission spin-state diagnostic.",
     false);
  static SwitchOption interfaceDiagnosePOWHEGRealSpinVertexYes
    (interfaceDiagnosePOWHEGRealSpinVertex,
     "Yes",
     "Emit a dedicated diagnostic for the realised POWHEG 2->3 spin state.",
     true);

  static Parameter<DISBase,unsigned long> interfacePOWHEGRealSpinDiagMax
    ("POWHEGRealSpinDiagMax",
     "Maximum number of POWHEG real-emission spin diagnostics to emit. Zero means unlimited. Active only when DISDiagnostics is On.",
     &DISBase::powhegRealSpinDiagMax_, 20, 0, 1000000,
     false, false, Interface::limited);
}

void DISBase::doinit() {
  HwMEBase::doinit();
  powhegRawMomentaDumpCount_ = 0;
  resetPOWHEGRawDiagnostics();
  if ((hasMECorrection() || hasPOWHEGCorrection() != No) && !alpha_) {
    throw InitException()
      << "DISBase requires a non-null Coupling reference when MEC or POWHEG"
      << " corrections are enabled."
      << Exception::abortnow;
  }
  // integrals of me over phase space
  double r5=sqrt(5.),darg((r5-1.)/(r5+1.)),ath(0.5*log((1.+1./r5)/(1.-1./r5)));
  comptonInt_ = 2.*(-21./20.-6./(5.*r5)*ath+sqr(Constants::pi)/3.
		    -2.*Math::ReLi2(1.-darg)-2.*Math::ReLi2(1.-1./darg));
  bgfInt_ = 121./9.-56./r5*ath;
  // extract the gluon ParticleData objects
  gluon_ = getParticleData(ParticleID::g);
}

Energy2 DISBase::powhegEmissionAlphaSScale(Energy2 q2, double xT) const {
  return useQ2ScaleInPOWHEGEmission_ ? q2 : 0.25*q2*sqr(xT);
}

Energy2 DISBase::powhegEmissionPDFScale(Energy2 q2, Energy2 mappedScale) const {
  return useQ2ScaleInPOWHEGEmission_ ? q2 : mappedScale;
}

double DISBase::powhegEmissionFixedOrderAlphaSValue(Energy2 scale) const {
  if (pdf_) {
    const ThePEG::LHAPDF *lhaPdf = dynamic_cast<const ThePEG::LHAPDF *>(pdf_.operator->());
    if (lhaPdf) return lhaPdf->alphaS(scale);
  }
  return SM().alphaS(scale);
}

double DISBase::powhegEmissionAlphaSValue(Energy2 scale) const {
  return useFixedOrderAlphaSInPOWHEGEmission_ ? powhegEmissionFixedOrderAlphaSValue(scale)
                                              : alpha_->value(scale);
}

double DISBase::powhegEmissionAlphaSOverestimate(Energy2 referenceScale) const {
  return useFixedOrderAlphaSInPOWHEGEmission_ ? powhegEmissionFixedOrderAlphaSValue(referenceScale) :
    alpha_->overestimateValue();
}

double DISBase::powhegEmissionAlphaSRatio(Energy2 scale, Energy2 referenceScale) const {
  if(!useFixedOrderAlphaSInPOWHEGEmission_) return alpha_->ratio(scale);
  const double alphaRef = powhegEmissionFixedOrderAlphaSValue(referenceScale);
  return alphaRef > 0.0 ? powhegEmissionFixedOrderAlphaSValue(scale)/alphaRef : 0.0;
}

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

  double Pl = 0.0, Pq = 0.0;
  if (lepton[0]) Pl = longPol(lepton[0]->spinInfo()->rhoMatrix());
  if (quark [0]) Pq = longPol(quark [0]->spinInfo()->rhoMatrix());
  acoeff_ = A_pol(lepton[0]->dataPtr(), lepton[1]->dataPtr(),
                  quark [0]->dataPtr(),  quark [1]->dataPtr(),
                  q2_, Pl, Pq);
}

void DISBase::dumpPOWHEGRawMomenta(bool isCompton, unsigned long eventIndex) const {
  if(!disDiagnostics_ || !dumpPOWHEGRawMomenta_) return;
  const vector<Lorentz5Momentum> & rawMom =
    isCompton ? ComptonMomenta_ : BGFMomenta_;
  if(rawMom.size() < 3) return;

  LorentzRotation toLab = rot_;
  toLab.invert();

  Lorentz5Momentum pOut1Lab = rawMom[1];
  Lorentz5Momentum pOut2Lab = rawMom[2];
  pOut1Lab.transform(toLab);
  pOut2Lab.transform(toLab);

  const double y = 2.0 / (l_ + 1.0);
  const double chosenPT = (isCompton ? pTCompton_ : pTBGF_) / GeV;
  const double chosenXP = isCompton ? comptonRawXP_ : bgfRawXP_;
  const double chosenZP = isCompton ? comptonRawZP_ : bgfRawZP_;
  std::ostringstream line;
  line << std::scientific << std::setprecision(16)
       << "POWHEG_RAW_MOM"
       << " event=" << eventIndex
       << " proc=" << (isCompton ? "Compton" : "BGF")
       << " contrib=" << contrib_
       << " Q2=" << q2_/GeV2
       << " xB=" << xB_
       << " y=" << y
       << " Ee=" << meMomenta()[0].e()/GeV
       << " s=" << currentPOWHEGRawBeamS()
       << " xp=" << chosenXP
       << " zp=" << chosenZP
       << " pT=" << chosenPT
       << " p1Bx=" << rawMom[1].x()/GeV
       << " p1By=" << rawMom[1].y()/GeV
       << " p1Bz=" << rawMom[1].z()/GeV
       << " p1BE=" << rawMom[1].e()/GeV
       << " p2Bx=" << rawMom[2].x()/GeV
       << " p2By=" << rawMom[2].y()/GeV
       << " p2Bz=" << rawMom[2].z()/GeV
       << " p2BE=" << rawMom[2].e()/GeV
       << " p1Lx=" << pOut1Lab.x()/GeV
       << " p1Ly=" << pOut1Lab.y()/GeV
       << " p1Lz=" << pOut1Lab.z()/GeV
       << " p1LE=" << pOut1Lab.e()/GeV
       << " p2Lx=" << pOut2Lab.x()/GeV
       << " p2Ly=" << pOut2Lab.y()/GeV
       << " p2Lz=" << pOut2Lab.z()/GeV
       << " p2LE=" << pOut2Lab.e()/GeV
       << "\n";
  generator()->log() << line.str();
}

RealEmissionProcessPtr DISBase::generateNativePOWHEGHardest(RealEmissionProcessPtr born,
							    bool allowBornFallback) {
  PPtr quark[2],lepton[2];
  PPtr hadron;
  unsigned int iqIn(0),iqOut(0);
  resetPOWHEGRawDiagnostics();
  initializePOWHEGEmissionState(born, quark, lepton, hadron, iqIn, iqOut);
  powhegRawDiagnostics_.q2 = q2_/GeV2;
  powhegRawDiagnostics_.xB = xB_;
  powhegRawDiagnostics_.y = 2.0 / (l_ + 1.0);
  powhegRawDiagnostics_.s = currentPOWHEGRawBeamS();

  generateCompton();
  generateBGF();
  unsigned long rawEventIndex = 0;
  const bool dumpRawEvent = beginPOWHEGRawDiagnosticEvent(rawEventIndex);
  if(pTCompton_<ZERO&&pTBGF_<ZERO) {
    if(allowBornFallback) {
      powhegRawDiagnostics_.fallback = 1;
      if(dumpRawEvent) dumpPOWHEGRawSummary(rawEventIndex);
      born->pT()[ShowerInteraction::QCD] = pTmin_;
      return born;
    }
    if(dumpRawEvent) dumpPOWHEGRawSummary(rawEventIndex);
    return RealEmissionProcessPtr();
  }

  bool isCompton=pTCompton_>pTBGF_;
  powhegRawDiagnostics_.winner = isCompton ? "Compton" : "BGF";
  if(dumpRawEvent) {
    dumpPOWHEGRawSummary(rawEventIndex);
    dumpPOWHEGRawMomenta(isCompton, rawEventIndex);
  }
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
    if(diagnosePOWHEGRealSpinVertex_) {
      diagnoseRealEmissionSpinState(born,isCompton);
    }
  }
  return born;
}

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

void DISBase::initializeMECorrection(RealEmissionProcessPtr born, double & initial,
				     double & final) {
  initial = initial_;
  final   = final_;
  // incoming particles
  for(unsigned int ix=0;ix<born->bornIncoming().size();++ix) {
      if(QuarkMatcher::Check(born->bornIncoming()[ix]->data())) {
      partons_[0] = born->bornIncoming()[ix]->dataPtr();
      pq_[0] = born->bornIncoming()[ix]->momentum();
    }
    else if(LeptonMatcher::Check(born->bornIncoming()[ix]->data())) {
      leptons_[0] = born->bornIncoming()[ix]->dataPtr();
      pl_[0] = born->bornIncoming()[ix]->momentum();
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
  // extract the born variables
  q_ =pl_[0]-pl_[1];
  q2_ = -q_.m2();
  double  yB = (q_*pq_[0])/(pl_[0]*pq_[0]); 
  l_ = 2./yB-1.;

  // calculate the A coefficient for the correlations, using beam polarisations
  double Pl = 0.0, Pq = 0.0;
  // read rho from the born incoming legs (where available)
  for (unsigned int ix=0; ix<born->bornIncoming().size(); ++ix) {
    if (LeptonMatcher::Check(born->bornIncoming()[ix]->data())) {
      Pl = longPol(born->bornIncoming()[ix]->spinInfo()->rhoMatrix());
    } else if (QuarkMatcher::Check(born->bornIncoming()[ix]->data())) {
      Pq = longPol(born->bornIncoming()[ix]->spinInfo()->rhoMatrix());
    }
  }
  acoeff_ = A_pol(leptons_[0], leptons_[1],
                  partons_[0],  partons_[1],
                  q2_, Pl, Pq);
}

RealEmissionProcessPtr DISBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  static const double eps=1e-6;
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
      xB_ = quark[0]->momentum().rho()/hadron->momentum().rho();
    }
    else if(LeptonMatcher::Check(born->bornIncoming()[ix]->data())) {
      lepton[0] = born->bornIncoming()[ix];
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
  vector<double> azicoeff;
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
    wgt *= pdf_->xfx(beam_,quark[0]->dataPtr(),pdfScale,xB_/xp)/
           pdf_->xfx(beam_,quark[0]->dataPtr(),q2_  ,xB_);
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
    wgt *= pdf_->xfx(beam_,gluon_              ,pdfScale,xB_/xp)/
           pdf_->xfx(beam_,quark[0]->dataPtr(),q2_  ,xB_);
    // other bits
    xperp = sqrt(4.*(1.-xp)*(1.-zp)*zp/xp);
    x1 = -1./xp;
    x2 = 1.-(1.-zp)/xp;
    x3 = 2.+x1-x2;
    // matrix element pieces
    azicoeff = BGFME(xp,x2,x3,xperp,true);
  }
  // compute the azimuthal average of the weight
  wgt *= (azicoeff[0]+0.5*azicoeff[2]);
  // decide whether or not to accept the weight
  if(UseRandom::rnd()>wgt) return RealEmissionProcessPtr();
  // if generate generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
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
//   debuggingMatrixElement(BGF,pin,p1,p2,gluon_,pl_[0],pl_[1],pq_[0],pq_[1],
// 			 lepton[0]->dataPtr(),lepton[1]->dataPtr(),
// 			 quark [0]->dataPtr(),quark [1]->dataPtr(),
// 			 q2_,phi,x2,x3,xperp,zp,xp,azicoeff,true);
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
    if(diagnosePOWHEGRealSpinVertex_) {
      diagnoseRealEmissionSpinState(born,!BGF);
    }
  }
  return born;
}

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
    vector<double> azicoeff = ComptonME(xp,x2,xperp,false);
    wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.+sqr(z))/final_;
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
      vector<double> azicoeff = ComptonME(xp,x2,xperp,false);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp*(1.-z)/(1.-xp)/(1.+sqr(z))/
	(1.-zp+xp-2.*xp*(1.-zp));
    }
    // BGF
    else {
      vector<double> azicoeff = BGFME(xp,x2,x3,xperp,true);
      wgt = (azicoeff[0]+0.5*azicoeff[2])*xp/(1.-zp+xp-2.*xp*(1.-zp))/(sqr(z)+sqr(1.-z));
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
//   static const double maxwgt = 2.,npow=0.34,ac=1.0;
//   double wgt;
//   do {
//     double rho = UseRandom::rnd();
//     xp = 1.-pow(rho,1./(1.-npow));
//     wgt = (sqr(xp)+ac+sqr(1.-xp));
//     if(wgt>1.+ac) cerr << "testing violates BGF maxA " << wgt << "\n";
//   }
//   while(wgt<UseRandom::rnd()*(1.+ac));
//   double xpwgt = -((6.-5.*npow+sqr(npow))*ac-3.*npow+sqr(npow)+4) 
//     /(sqr(npow)*(npow-6.)+11.*npow-6.);
//   xpwgt *= pow(1.-xp,npow)/wgt;
//   double xp2(sqr(xp)),lxp(log(xp)),xp4(sqr(xp2)),lxp1(log(1.-xp));
//   double zpwgt = (2.*xp4*(lxp+lxp1-3.)+4.*xp*xp2*(3.-lxp-lxp1)
// 		  +xp2*(-13.+lxp+lxp1)+xp*(+7.+lxp+lxp1)-lxp-lxp1-1.)/(1.+xp-xp2);
//   do {
//     double zpmax = 1./(1.+xp*(1.-xp)), zpmin = 1.-zpmax;
//     zp = 1.-pow((1.-zpmin)/(1.-zpmax),UseRandom::rnd())*(1.-zpmax);
//     wgt = log((1.-zpmin)/(1.-zpmax))*(1.-zp);
//     double x1 = -1./xp;
//     double x2 = 1.-(1.-zp)/xp;
//     double x3 = 2.+x1-x2;
//     double xperp2 = 4.*(1.-xp)*(1.-zp)*zp/xp;
//     wgt *= sqr(xp)/(1.-zp)*(sqr(x3)+sqr(x2)+3.*xperp2);
//     if(wgt>maxwgt*zpwgt) cerr << "testing violates BGF maxB " << wgt/xpwgt << "\n";
//   }
//   while(wgt<UseRandom::rnd()*maxwgt);
//   return zpwgt*xpwgt;
}

vector<double> DISBase::ComptonME(double xp, double x2, double xperp,
				  bool norm) {
  vector<double> output(3,0.);
  double cos2 =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2 = xperp/sqrt(sqr(x2)+sqr(xperp));
  double root = sqrt(sqr(l_)-1.);
  output[0] = sqr(cos2)+acoeff_*cos2*l_+sqr(l_);
  output[1] = -acoeff_*cos2*root*sin2-2.*l_*root*sin2;
  output[2] = sqr(root)*sqr(sin2);
  double lo(1+acoeff_*l_+sqr(l_));
  double denom = norm ? 1.+sqr(xp)*(sqr(x2)+1.5*sqr(xperp)) : 1.;
  double fact  = sqr(xp)*(sqr(x2)+sqr(xperp))/lo;
  for(unsigned int ix=0;ix<output.size();++ix) 
    output[ix] = ((ix==0 ? 1. : 0.) +fact*output[ix])/denom;
  return output;
}

vector<double> DISBase::BGFME(double xp, double x2, double x3, 
			      double xperp, bool norm) {
  vector<double> output(3,0.);
  double cos2  =   x2 /sqrt(sqr(x2)+sqr(xperp));
  double sin2  = xperp/sqrt(sqr(x2)+sqr(xperp));
  double fact2 = sqr(xp)*(sqr(x2)+sqr(xperp));
  double cos3  =   x3 /sqrt(sqr(x3)+sqr(xperp));
  double sin3  = xperp/sqrt(sqr(x3)+sqr(xperp));
  double fact3 = sqr(xp)*(sqr(x3)+sqr(xperp));
  double root = sqrt(sqr(l_)-1.);
  output[0] = fact2*(sqr(cos2)+acoeff_*cos2*l_+sqr(l_)) +
              fact3*(sqr(cos3)-acoeff_*cos3*l_+sqr(l_));
  output[1] = - fact2*(acoeff_*cos2*root*sin2+2.*l_*root*sin2)
              - fact3*(acoeff_*cos3*root*sin3-2.*l_*root*sin3);
  output[2] = fact2*(sqr(root)*sqr(sin2)) +
              fact3*(sqr(root)*sqr(sin3));
  double lo(1+acoeff_*l_+sqr(l_));
  double denom = norm ? sqr(xp)*(sqr(x3)+sqr(x2)+3.*sqr(xperp))*lo : lo;
  for(unsigned int ix=0;ix<output.size();++ix) output[ix] /= denom;
  return output;
}

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

void DISBase::generateCompton() {
  POWHEGRawChannelDiagnostics & raw = powhegRawDiagnostics_.compton;
  comptonRawXP_ = 0.0;
  comptonRawZP_ = 0.0;
  // maximum value of the xT
  double xT = sqrt((1.-xB_)/xB_);
  double xTMin = 2.*pTmin_/sqrt(q2_);
  raw.xTMin = xTMin;
  double zp;
  Energy2 alphaRefScale = powhegEmissionAlphaSScale(q2_,xTMin);
  // prefactor
  double a = powhegEmissionAlphaSOverestimate(alphaRefScale)*comptonWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  vector<double> azicoeff;
  while(true) {
    wgt = 0.;
    raw.phase = 0.0;
    raw.pdfRatio = 0.0;
    raw.alphaRatio = 0.0;
    raw.meAvg = 0.0;
    raw.wgt = 0.0;
    raw.xMapped = 0.0;
    raw.pdfScale = 0.0;
    raw.alphaScale = 0.0;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // zp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    ++raw.trials;
    raw.xp = xp;
    raw.zp = zp;
    raw.xT = xT;
    raw.pT = 0.5*sqrt(q2_)*xT/GeV;
    // check allowed
    if(xp<xB_||xp>1.) {
      ++raw.rejectXP;
      if(xT<=xTMin) break;
      continue;
    }
    raw.xMapped = xB_/xp;
    // phase-space piece of the weight
    raw.phase = 8.*(1.-xp)*zp/comptonWeight_;
    wgt = raw.phase;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    Energy2 pdfScale = powhegEmissionPDFScale(q2_, scale);
    raw.pdfScale = pdfScale/GeV2;
    raw.pdfRatio = pdf_->xfx(beam_,partons_[0],pdfScale,xB_/xp)/
                   pdf_->xfx(beam_,partons_[0],q2_  ,xB_);
    wgt *= raw.pdfRatio;
    // me piece of the weight
    double x2 = 1.-(1.-zp)/xp;
    azicoeff = ComptonME(xp,x2,xT,false);
    Energy2 alphaScale = powhegEmissionAlphaSScale(q2_,xT);
    raw.alphaScale = alphaScale/GeV2;
    raw.alphaRatio = 4./3.*powhegEmissionAlphaSRatio(alphaScale,alphaRefScale);
    raw.meAvg = azicoeff[0]+0.5*azicoeff[2];
    wgt *= raw.alphaRatio*raw.meAvg;
    raw.wgt = wgt;
    if(wgt>1.||wgt<0.) {
      if(wgt>1.) ++raw.weightHigh;
      if(wgt<0.) ++raw.weightNeg;
      ostringstream wstring;
      wstring << "DISBase::generateCompton() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
    if(xT<=xTMin) break;
    if(UseRandom::rnd()>wgt) {
      ++raw.rejectVeto;
      continue;
    }
    break;
  }
  if(xT<=xTMin) {
    pTCompton_=-GeV;
    return;
  }
  raw.status = "accepted";
  comptonRawXP_ = xp;
  comptonRawZP_ = zp;
  // generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
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
//   debuggingMatrixElement(false,p0,p1,p2,gluon_,pl_[0],pl_[1],pq_[0],pq_[1],
// 			 leptons_[0],leptons_[1],
// 			 partons_[0],partons_[1],
// 			 q2_,phi,x2,x3,xT,zp,xp,azicoeff,false);
}

void DISBase::generateBGF() {
  POWHEGRawChannelDiagnostics & raw = powhegRawDiagnostics_.bgf;
  bgfRawXP_ = 0.0;
  bgfRawZP_ = 0.0;
  // maximum value of the xT
  double xT = (1.-xB_)/xB_;
  double xTMin = 2.*max(pTmin_,pTCompton_)/sqrt(q2_);
  raw.xTMin = xTMin;
  double zp;
  Energy2 alphaRefScale = powhegEmissionAlphaSScale(q2_,xTMin);
  // prefactor
  double a = powhegEmissionAlphaSOverestimate(alphaRefScale)*BGFWeight_/Constants::twopi;
  // loop to generate kinematics
  double wgt(0.),xp(0.);
  vector<double> azicoeff;
  while(true) {
    wgt = 0.;
    raw.phase = 0.0;
    raw.pdfRatio = 0.0;
    raw.alphaRatio = 0.0;
    raw.meAvg = 0.0;
    raw.wgt = 0.0;
    raw.xMapped = 0.0;
    raw.pdfScale = 0.0;
    raw.alphaScale = 0.0;
    // intergration variables dxT/xT^3
    xT *= 1./sqrt(1.-2.*log(UseRandom::rnd())/a*sqr(xT));
    // zp
    zp = UseRandom::rnd();
    xp = 1./(1.+0.25*sqr(xT)/zp/(1.-zp));
    ++raw.trials;
    raw.xp = xp;
    raw.zp = zp;
    raw.xT = xT;
    raw.pT = 0.5*sqrt(q2_)*xT/GeV;
    // check allowed
    if(xp<xB_||xp>1.) {
      ++raw.rejectXP;
      if(xT<=xTMin) break;
      continue;
    }
    raw.xMapped = xB_/xp;
    // phase-space piece of the weight
    raw.phase = 8.*sqr(1.-xp)*zp/BGFWeight_;
    wgt = raw.phase;
    // PDF piece of the weight
    Energy2 scale = q2_*((1.-xp)*(1-zp)*zp/xp+1.);
    Energy2 pdfScale = powhegEmissionPDFScale(q2_, scale);
    raw.pdfScale = pdfScale/GeV2;
    raw.pdfRatio = pdf_->xfx(beam_,gluon_     ,pdfScale,xB_/xp)/
                   pdf_->xfx(beam_,partons_[0],q2_  ,xB_);
    wgt *= raw.pdfRatio;
    // me piece of the weight
    double x1 = -1./xp;
    double x2 = 1.-(1.-zp)/xp;
    double x3 = 2.+x1-x2;
    azicoeff = BGFME(xp,x2,x3,xT,false);
    Energy2 alphaScale = powhegEmissionAlphaSScale(q2_,xT);
    raw.alphaScale = alphaScale/GeV2;
    raw.alphaRatio = 0.5*powhegEmissionAlphaSRatio(alphaScale,alphaRefScale);
    raw.meAvg = azicoeff[0]+0.5*azicoeff[2];
    wgt *= raw.alphaRatio*raw.meAvg;
    raw.wgt = wgt;
    if(wgt>1.||wgt<0.) {
      if(wgt>1.) ++raw.weightHigh;
      if(wgt<0.) ++raw.weightNeg;
      ostringstream wstring;
      wstring << "DISBase::generateBGF() "
	      << "Weight greater than one or less than zero"
	      << "wgt = " << wgt << "\n";
      generator()->logWarning( Exception(wstring.str(),
					 Exception::warning) );
    }
    if(xT<=xTMin) break;
    if(UseRandom::rnd()>wgt) {
      ++raw.rejectVeto;
      continue;
    }
    break;
  }
  if(xT<=xTMin) {
    pTBGF_=-GeV;
    return;
  }
  raw.status = "accepted";
  bgfRawXP_ = xp;
  bgfRawZP_ = zp;
  // generate phi
  unsigned int itry(0);
  double phimax = std::accumulate(azicoeff.begin(),azicoeff.end(),0.);
  double phiwgt,phi;
  do {
    phi = UseRandom::rnd()*Constants::twopi;
    double cphi(cos(phi));
    phiwgt = azicoeff[0]+azicoeff[1]*cphi+azicoeff[2]*sqr(cphi);
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
//   debuggingMatrixElement(true,p0,p1,p2,gluon_,pl_[0],pl_[1],pq_[0],pq_[1],
// 			 leptons_[0],leptons_[1],
// 			 partons_[0],partons_[1],
// 			 q2_,phi,x2,x3,xT,zp,xp,azicoeff,false);
}

int DISBase::nDim() const {
  return HwMEBase::nDim() + (contrib_>0 ? 1 : 0 );
}

bool DISBase::generateKinematics(const double * r) {
  // Born kinematics
  if(!HwMEBase::generateKinematics(r)) return false;
  if(contrib_!=0) {
    // hadron and momentum fraction
    if(HadronMatcher::Check(*lastParticles().first->dataPtr())) {
      hadron_ = dynamic_ptr_cast<tcBeamPtr>(lastParticles().first->dataPtr());
      xB_ = lastX1();
    }
    else {
      hadron_ = dynamic_ptr_cast<tcBeamPtr>(lastParticles().second->dataPtr());
      xB_ = lastX2();
    }
    // Q2
    q2_ = -(meMomenta()[0]-meMomenta()[2]).m2();
    // xp
    int ndim=nDim();
    double rhomin = pow(1.-xB_,1.-power_); 
    double rho = r[ndim-1]*rhomin;
    xp_ = 1.-pow(rho,1./(1.-power_));
    jac_ = rhomin/(1.-power_)*pow(1.-xp_,power_);
    jacobian(jacobian()*jac_);
  }
  return true; 
}

Energy2 DISBase::scale() const {
  return scaleOpt_ == 1 ? 
    -sqr(scaleFact_)*tHat() : sqr(scaleFact_*muF_);
}

CrossSection DISBase::dSigHatDR() const {
  return NLOWeight()*HwMEBase::dSigHatDR();
}

double DISBase::NLOWeight() const {

  // If only leading order is required return 1:
  if (contrib_ == 0) return 1.;

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
  {
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

  // Longitudinal hadron beam polarisation (convention: false = proton beam)
  const double Pz = getBeamPolarization(false).z();

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

#if 0
  // Diagnostic: print one representative point per quark flavour.
  {
    static std::set<int> seen_flavors;
    int qid = mePartonData()[1]->id();
    if (hasDiffPdf && seen_flavors.find(qid) == seen_flavors.end()) {
      seen_flavors.insert(qid);
      generator()->log() << "NLO_DIAG_FLAVOR"
        << " quark=" << mePartonData()[1]->PDGName()
        << " pid=" << qid
        << " xB=" << xB_ << " mu2=" << mu2/GeV2
        << " dloPDF=" << dloPDF
        << " loPDF=" << loPDF
        << " Pq=" << Pq
        << " dqPDF=" << dqPDF
        << " dgPDF=" << dgPDF
        << "\n";
    }
  }
#endif


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

#if 0
  // Diagnostic: compare the extractor paths used by the rho matrices and PDFs.
  {
    static bool first_diag = true;
    if (first_diag) {
      first_diag = false;

      auto le = lastExtractor();
      bool le_is_ppe = bool(ThePEG::dynamic_ptr_cast<
        ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::tptr>(le));

      auto pe = generator()->eventHandler()->partonExtractor();
      bool pe_is_ppe = bool(ThePEG::dynamic_ptr_cast<
        ThePEG::Ptr<ThePEG::PolarizedPartonExtractor>::pointer>(pe));

      bool same_obj = (le.operator->() == pe.operator->());

      generator()->log() << "NLO_DIAG_EXTRACTOR"
        << " lastExtractor_is_PPE=" << le_is_ppe
        << " eventHandler_PE_is_PPE=" << pe_is_ppe
        << " same_object=" << same_obj
        << " diffPdf=" << (diffPdf ? "OK" : "NULL")
        << " Pz=" << Pz
        << " hasDiffPdf=" << hasDiffPdf
        << "\n";
    }
  }
#endif

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

  // Collinear counterterms:
  // The quark and gluon channels use projector weights derived from the exact
  // Born angular structure. In the photon limit these reduce to the old scalar
  // f_pol blend. For full NC the quark and gluon projectors differ because the
  // parity-odd spin-even terms have no gluon kernel at this order.

  double logRatio = log((1.-xp_)*q2_/xp_/mu2);

  // --- gluon collinear ---
  // These are divided by loPDF later, so we compute collX/loPDF for each part.
  //
  // unpolarised piece / loPDF:
  //   TR/xp * (gPDF/loPDF) * [2xp(1-xp) + (xp²+(1-xp)²)*log(...)]
  double collg_over_born_unpol =
    TRfact/xp_*(gPDF/loPDF)*(2.*xp_*(1.-xp_) + (sqr(xp_)+sqr(1.-xp_))*logRatio);

  // polarised piece / (Pz*dloPDF):
  //   TR/xp * (Pz*Δg(x)/(Pz*Δq(xB))) * [2(1-xp) + (2xp-1)*log(...)]
  double collg_over_born_pol = 0.0;
  double dqRatio = 0.0;
  double dgRatio = 0.0;
  if (hasDiffPdf && std::abs(dloPDF) > minDlo) {
    dqRatio = dqPDF / dloPDF;
    dgRatio = dgPDF / dloPDF;
    collg_over_born_pol =
      TRfact/xp_*dgRatio*(2.*(1.-xp_) + (2.*xp_-1.)*logRatio);
  }

  double collg_over_born =
    blend.gUnpolarized * collg_over_born_unpol +
    blend.gPolarized   * collg_over_born_pol;

  // --- quark collinear ---
  // Same Pqq kernel, but different PDF ratios for polarised part.
  //
  // unpolarised piece / loPDF:
  double collq_over_born_unpol =
    CFfact/xp_*(qPDF/loPDF)*(1.-xp_-2./(1.-xp_)*log(xp_)-(1.+xp_)*log((1.-xp_)/xp_*q2_/mu2))+
    CFfact/xp_*(qPDF/loPDF-xp_)*(2./(1.-xp_)*log(q2_*(1.-xp_)/mu2)-1.5/(1.-xp_));

  // polarised piece / (Pz*dloPDF):
  double collq_over_born_pol = 0.0;
  if (hasDiffPdf && std::abs(dloPDF) > minDlo) {
    collq_over_born_pol =
      CFfact/xp_*dqRatio*(1.-xp_-2./(1.-xp_)*log(xp_)-(1.+xp_)*log((1.-xp_)/xp_*q2_/mu2))+
      CFfact/xp_*(dqRatio-xp_)*(2./(1.-xp_)*log(q2_*(1.-xp_)/mu2)-1.5/(1.-xp_));
  }

  double collq_over_born =
    blend.qUnpolarized * collq_over_born_unpol +
    blend.qPolarized   * collq_over_born_pol;

  // Real-emission kernels are written as sigma_real / sigma_Born.
  // Keep a_born in the denominator so the Born factor cancels exactly, and use
  // mapped analyzing powers at x = xB/xp_ for the spin dependence.

  // Mapped analysing powers for the real-emission kernels
  double a_q_mapped = A_pol(mePartonData()[0], mePartonData()[2],
                            mePartonData()[1], mePartonData()[3],
                            q2_, Pl, Pq_m);

  // BGF mapped analysing powers. These are kept for diagnostics and for
  // monitoring the exact R2/R3 crossing structure of A_pol(); the integrated
  // finite remainder below is assembled with the exact gluon projectors.

  // R2: Born parton (quark or antiquark as selected for the Born diagram)
  double a_g_R2 = A_pol(mePartonData()[0], mePartonData()[2],
                        mePartonData()[1], mePartonData()[3],
                        q2_, Pl, Pg_m);

  // R3: charge-conjugate of Born parton with opposite helicity (-Pg_m)
  tcPDPtr qbar_in  = mePartonData()[1]->CC();
  tcPDPtr qbar_out = mePartonData()[3]->CC();
  if (!qbar_in)  qbar_in  = tcPDPtr(mePartonData()[1]);
  if (!qbar_out) qbar_out = tcPDPtr(mePartonData()[3]);
  double a_g_R3 = A_pol(mePartonData()[0], mePartonData()[2],
                        qbar_in, qbar_out,
                        q2_, Pl, -Pg_m);

  // Unpolarised PDF ratios (no f_pol blending — see comment above)
  double qRatio = qPDF / loPDF;
  double gRatio = gPDF / loPDF;

  // q -> qg term (QCDC): denominator a_born, kernel a_q_mapped
  const double qcdcDenRatio =
    qcdcMappedDenominatorRatio(mePartonData()[0], mePartonData()[2],
                               mePartonData()[1], mePartonData()[3],
                               q2_, Pl, Pq, Pq_m);
  double realq = qcdcDenRatio * CFfact/xp_/(1.+a_born*l+sqr(l))*qRatio*
    (2.+2.*sqr(l)-xp_+3.*xp_*sqr(l)+a_q_mapped*l*(2.*xp_+1.));

  // g -> q qbar term (BGF): use the same exact gluon projectors as the
  // collinear counterterm. At this order the finite gluon remainder only
  // contributes to the F2/FL-like unpolarized channel and the G2-like
  // spin-dependent channel.
  const double realg_over_born_unpol =
    -TRfact/xp_ * gRatio *
    ((1.+sqr(l)+2.*(1.-3.*sqr(l))*xp_*(1.-xp_)) / (1.+sqr(l)));

  double realg_over_born_pol = 0.0;
  if (hasDiffPdf && std::abs(dloPDF) > minDlo) {
    realg_over_born_pol =
      -2.0 * TRfact/xp_ * dgRatio * (2.*xp_ - 1.);
  }

  double realg =
    blend.gUnpolarized * realg_over_born_unpol +
    blend.gPolarized   * realg_over_born_pol;

  // Full NLO/Born ratio.
  double wgt = virt + (collq_over_born + collg_over_born + realq + realg);

  struct APolDecomp {
    double actual;
    double at_zero;
    double even;
    double odd;
  };

  auto apolDecomp = [&](tcPDPtr lin, tcPDPtr lout,
                        tcPDPtr qin, tcPDPtr qout,
                        double P, double actual) {
    APolDecomp out{actual, A_pol(lin, lout, qin, qout, q2_, Pl, 0.0), 0.0, 0.0};
    if (std::abs(P) <= 1e-15) {
      out.even = out.at_zero;
      return out;
    }
    double flipped = A_pol(lin, lout, qin, qout, q2_, Pl, -P);
    out.even = 0.5 * (actual + flipped);
    out.odd  = 0.5 * (actual - flipped);
    return out;
  };

#if 0
  // Diagnostic dump for selected NLO points.
  {
    static unsigned long nlo_diag_counter_ = 0;
    ++nlo_diag_counter_;
    if (nlo_diag_counter_ <= 50 || nlo_diag_counter_ % 50000 == 0) {
      double realq_x_Born = CFfact/xp_*qRatio*
        (2.+2.*sqr(l)-xp_+3.*xp_*sqr(l)+a_q_mapped*l*(2.*xp_+1.));
      double realg_x_Born = -TRfact/xp_*gRatio*
        ((1.+sqr(l)+2.*(1.-3.*sqr(l))*xp_*(1.-xp_))
         +2.*l*(sqr(xp_)*a_g_R2+sqr(1.-xp_)*a_g_R3));

      double f2g_old = gPDF/xp_*TRfact*((sqr(1.-xp_)+sqr(xp_))*log((1.-xp_)/xp_)+
                                          8.*xp_*(1.-xp_)-1.);
      double f2q_old =
        loPDF/jac_*(1.+CFfact*(-1.5*log(1.-xB_)+sqr(log(1.-xB_))
                               -sqr(Constants::pi)/3.-4.5))
        +qPDF            *CFfact/xp_*(3.+2.*xp_-(1.+xp_)*log(1.-xp_)
                                      -(1.+sqr(xp_))/(1.-xp_)*log(xp_))
        +(qPDF-xp_*loPDF)*CFfact/xp_*(2.*log(1.-xp_)/(1.-xp_)-1.5/(1.-xp_));
      double wgt_old = (f2g_old+f2q_old)/loPDF;

      double Pq_from_pdf = (std::abs(Pz) > 1e-12 && loPDF != 0.0 && dloPDF != 0.0)
                           ? Pz * dloPDF / loPDF : 0.0;

      generator()->log() << "NLO_DIAG " << nlo_diag_counter_
        << " xB=" << xB_ << " xp=" << xp_ << " y=" << yB
        << " Q2=" << q2_/GeV2
        << " l=" << l
        << " Pl=" << Pl << " Pq=" << Pq
        << " Pz=" << Pz
        << " Pq_from_pdf=" << Pq_from_pdf
        << " Pq_m=" << Pq_m << " Pg_m=" << Pg_m
        << " f_pol_q=" << blend.qPolarized
        << " f_pol_g=" << blend.gPolarized
        << " a_born=" << a_born
        << " a_q_m=" << a_q_mapped
        << " a_gR2=" << a_g_R2 << " a_gR3=" << a_g_R3
        << "\n";
      generator()->log() << "NLO_DIAG_WGT"
        << " virt=" << virt
        << " collq_u=" << collq_over_born_unpol
        << " collq_p=" << collq_over_born_pol
        << " collq=" << collq_over_born
        << " collg_u=" << collg_over_born_unpol
        << " collg_p=" << collg_over_born_pol
        << " collg=" << collg_over_born
        << " realq=" << realq
        << " realg=" << realg
        << " wgt=" << wgt
        << " wgt_old=" << wgt_old
        << " contrib=" << contrib_
        << " jac=" << jac_
        << "\n";
      generator()->log() << "NLO_DIAG_PDFS"
        << " loPDF=" << loPDF
        << " qPDF=" << qPDF << " gPDF=" << gPDF
        << " dloPDF=" << dloPDF << " dqPDF=" << dqPDF << " dgPDF=" << dgPDF
        << " qR=" << qRatio << " gR=" << gRatio
        << "\n";
      generator()->log() << "NLO_DIAG_XBORN"
        << " rqxB=" << realq_x_Born
        << " rgxB=" << realg_x_Born
        << " rqxB_chk=" << realq*(1.+a_born*l+sqr(l))
        << " rgxB_chk=" << realg*(1.+a_born*l+sqr(l))
        << "\n";
    }
  }
#endif

  // Running component diagnostics for the NLO correction. The F_* entries are
  // cumulative fractions of the physical POSNLO/NEGNLO run cross section for
  // this specific run/helicity. They can be combined with the matching .out
  // cross sections to reconstruct the component contributions to sigma_LL and
  // to the ALL-GAMMA-Z interference.
  {
    double w_Born = me2() * jacobian() / jac_;
    double dSigmaBorn = jac_ * w_Born;
    bool keep = (contrib_ == 1) ? (wgt > 0.0) : (wgt < 0.0);
    double runSign = (contrib_ == 1) ? 1.0 : -1.0;

    static double sRun = 0.0;
    static double sVirt = 0.0;
    static double sCqEven = 0.0, sCqOdd = 0.0;
    static double sCgEven = 0.0, sCgOdd = 0.0;
    static double sRq = 0.0, sRg = 0.0;
    static double sRqEven = 0.0, sRqOdd = 0.0;
    static double sRgEven = 0.0, sRgOdd = 0.0;
    static double sABorn = 0.0, sABorn0 = 0.0, sABornEven = 0.0, sABornOdd = 0.0;
    static double sAQMap = 0.0, sAQMap0 = 0.0, sAQMapEven = 0.0, sAQMapOdd = 0.0;
    static double sAGR2 = 0.0, sAGR20 = 0.0, sAGR2Even = 0.0, sAGR2Odd = 0.0;
    static double sAGR3 = 0.0, sAGR30 = 0.0, sAGR3Even = 0.0, sAGR3Odd = 0.0;
    static double sBornMe = 0.0, sBornPred = 0.0, sBornClosure = 0.0;
    static double sMapMe = 0.0, sMapPred = 0.0, sMapClosure = 0.0;
    static double sMeSelf = 0.0;
    static double sCoeffScale = 0.0;
    static double sCoeffPlMe = 0.0, sCoeffPlPred = 0.0;
    static double sCoeffPqMe = 0.0, sCoeffPqPred = 0.0;
    static double sCoeffPlPqMe = 0.0, sCoeffPlPqPred = 0.0;
    static unsigned long born_n = 0;
    static unsigned long accepted_n = 0;
    static const unsigned long diagPeriod = 500000;

    if (keep) {
      ++accepted_n;
      double dSigmaRun = runSign * wgt * dSigmaBorn;
      sRun    += dSigmaRun;
      sVirt   += runSign * virt * dSigmaBorn;
      sCqEven += runSign * blend.qUnpolarized * collq_over_born_unpol * dSigmaBorn;
      sCqOdd  += runSign * blend.qPolarized   * collq_over_born_pol   * dSigmaBorn;
      sCgEven += runSign * blend.gUnpolarized * collg_over_born_unpol * dSigmaBorn;
      sCgOdd  += runSign * blend.gPolarized   * collg_over_born_pol   * dSigmaBorn;
      sRq     += runSign * realq * dSigmaBorn;
      sRg     += runSign * realg * dSigmaBorn;

      double born_d = 1. + a_born*l + sqr(l);
      if (std::abs(born_d) > 1e-30) {
        double rq_pref = CFfact/xp_/born_d * qRatio;
        sRqEven += runSign * rq_pref * (2.+2.*sqr(l)-xp_+3.*xp_*sqr(l)) * dSigmaBorn;
        sRqOdd  += runSign * rq_pref * a_q_mapped*l*(2.*xp_+1.) * dSigmaBorn;

        sRgEven += runSign * blend.gUnpolarized * realg_over_born_unpol * dSigmaBorn;
        sRgOdd  += runSign * blend.gPolarized   * realg_over_born_pol   * dSigmaBorn;
      }

      const APolDecomp bornParts = apolDecomp(mePartonData()[0], mePartonData()[2],
                                              mePartonData()[1], mePartonData()[3],
                                              Pq, a_born);
      const APolDecomp qMapParts = apolDecomp(mePartonData()[0], mePartonData()[2],
                                              mePartonData()[1], mePartonData()[3],
                                              Pq_m, a_q_mapped);
      const APolDecomp gR2Parts = apolDecomp(mePartonData()[0], mePartonData()[2],
                                             mePartonData()[1], mePartonData()[3],
                                             Pg_m, a_g_R2);
      const APolDecomp gR3Parts = apolDecomp(mePartonData()[0], mePartonData()[2],
                                             qbar_in, qbar_out,
                                             -Pg_m, a_g_R3);

      sABorn     += dSigmaRun * bornParts.actual;
      sABorn0    += dSigmaRun * bornParts.at_zero;
      sABornEven += dSigmaRun * bornParts.even;
      sABornOdd  += dSigmaRun * bornParts.odd;

      sAQMap     += dSigmaRun * qMapParts.actual;
      sAQMap0    += dSigmaRun * qMapParts.at_zero;
      sAQMapEven += dSigmaRun * qMapParts.even;
      sAQMapOdd  += dSigmaRun * qMapParts.odd;

      sAGR2     += dSigmaRun * gR2Parts.actual;
      sAGR20    += dSigmaRun * gR2Parts.at_zero;
      sAGR2Even += dSigmaRun * gR2Parts.even;
      sAGR2Odd  += dSigmaRun * gR2Parts.odd;

      sAGR3     += dSigmaRun * gR3Parts.actual;
      sAGR30    += dSigmaRun * gR3Parts.at_zero;
      sAGR3Even += dSigmaRun * gR3Parts.even;
      sAGR3Odd  += dSigmaRun * gR3Parts.odd;

      BornClosureDiagnostics bornDiag;
      if (bornClosureDiagnostics(Pl, Pq, Pq_m, l, bornDiag)
          && std::abs(bornDiag.me2Zero) > 1e-30
          && std::abs(bornDiag.sigmaZero) > 1e-30) {
        ++born_n;
        const double bornMe = bornDiag.me2Born / bornDiag.me2Zero;
        const double bornPred = bornDiag.sigmaBorn / bornDiag.sigmaZero;
        const double mapMe = bornDiag.me2Mapped / bornDiag.me2Zero;
        const double mapPred = bornDiag.sigmaMapped / bornDiag.sigmaZero;

        sBornMe += bornMe;
        sBornPred += bornPred;
        if (std::abs(bornPred) > 1e-30) sBornClosure += bornMe / bornPred;

        sMapMe += mapMe;
        sMapPred += mapPred;
        if (std::abs(mapPred) > 1e-30) sMapClosure += mapMe / mapPred;

        if (std::abs(bornDiag.me2Born) > 1e-30) sMeSelf += me2() / bornDiag.me2Born;

        sCoeffScale += bornDiag.coeffScale;
        sCoeffPlMe += bornDiag.coeffPlMe;
        sCoeffPlPred += bornDiag.coeffPlPred;
        sCoeffPqMe += bornDiag.coeffPqMe;
        sCoeffPqPred += bornDiag.coeffPqPred;
        sCoeffPlPqMe += bornDiag.coeffPlPqMe;
        sCoeffPlPqPred += bornDiag.coeffPlPqPred;
      }
    }

    if (disDiagnostics_ && dumpNLOTermDiagnostics_ && accepted_n > 0 &&
        accepted_n % diagPeriod == 0 && std::abs(sRun) > 1e-30) {
      double F_total = (sVirt + sCqEven + sCqOdd + sCgEven + sCgOdd + sRq + sRg) / sRun;
      const std::string runName = generator()->runName();

      generator()->log() << "NLO_TERM_CUM"
        << " run=" << runName
        << " n=" << accepted_n
        << " F_virt=" << sVirt/sRun
        << " F_cq_even=" << sCqEven/sRun
        << " F_cq_odd=" << sCqOdd/sRun
        << " F_cg_even=" << sCgEven/sRun
        << " F_cg_odd=" << sCgOdd/sRun
        << " F_rq=" << sRq/sRun
        << " F_rg=" << sRg/sRun
        << " F_total=" << F_total
        << "\n";

      generator()->log() << "NLO_TERM_REAL"
        << " run=" << runName
        << " n=" << accepted_n
        << " F_rq_even=" << sRqEven/sRun
        << " F_rq_odd=" << sRqOdd/sRun
        << " F_rg_even=" << sRgEven/sRun
        << " F_rg_odd=" << sRgOdd/sRun
        << " F_rq_chk=" << (sRqEven+sRqOdd)/sRun
        << " F_rg_chk=" << (sRgEven+sRgOdd)/sRun
        << "\n";

      generator()->log() << "NLO_TERM_APOL"
        << " run=" << runName
        << " n=" << accepted_n
        << " A_born=" << sABorn/sRun
        << " A_born0=" << sABorn0/sRun
        << " A_born_even=" << sABornEven/sRun
        << " A_born_odd=" << sABornOdd/sRun
        << " A_qmap=" << sAQMap/sRun
        << " A_qmap0=" << sAQMap0/sRun
        << " A_qmap_even=" << sAQMapEven/sRun
        << " A_qmap_odd=" << sAQMapOdd/sRun
        << " A_gR2=" << sAGR2/sRun
        << " A_gR20=" << sAGR20/sRun
        << " A_gR2_even=" << sAGR2Even/sRun
        << " A_gR2_odd=" << sAGR2Odd/sRun
        << " A_gR3=" << sAGR3/sRun
        << " A_gR30=" << sAGR30/sRun
        << " A_gR3_even=" << sAGR3Even/sRun
        << " A_gR3_odd=" << sAGR3Odd/sRun
        << "\n";

      if (born_n > 0) {
        generator()->log() << "NLO_TERM_BORN"
          << " run=" << runName
          << " n=" << born_n
          << " B_me=" << sBornMe/born_n
          << " B_pred=" << sBornPred/born_n
          << " B_closure=" << sBornClosure/born_n
          << " M_me=" << sMapMe/born_n
          << " M_pred=" << sMapPred/born_n
          << " M_closure=" << sMapClosure/born_n
          << " ME_self=" << sMeSelf/born_n
          << "\n";

        generator()->log() << "NLO_TERM_COEFF"
          << " run=" << runName
          << " n=" << born_n
          << " C00_scale=" << sCoeffScale/born_n
          << " CPl_me=" << sCoeffPlMe/born_n
          << " CPl_pred=" << sCoeffPlPred/born_n
          << " CPl_close=" << ((std::abs(sCoeffPlPred) > 1e-30) ? (sCoeffPlMe / sCoeffPlPred) : 0.0)
          << " CPq_me=" << sCoeffPqMe/born_n
          << " CPq_pred=" << sCoeffPqPred/born_n
          << " CPq_close=" << ((std::abs(sCoeffPqPred) > 1e-30) ? (sCoeffPqMe / sCoeffPqPred) : 0.0)
          << " CPlPq_me=" << sCoeffPlPqMe/born_n
          << " CPlPq_pred=" << sCoeffPlPqPred/born_n
          << " CPlPq_close=" << ((std::abs(sCoeffPlPqPred) > 1e-30) ? (sCoeffPlPqMe / sCoeffPlPqPred) : 0.0)
          << "\n";
      }
    }
  }

  //   double f2g = gPDF/xp_*TRfact*((sqr(1-xp_)+sqr(xp_))*log((1-xp_)/xp_)+
  // 				8*xp_*(1.-xp_)-1.);
  //   double f2q =
  //     loPDF/jac_*(1.+CFfact*(-1.5*log(1.-xB_)+sqr(log(1.-xB_))
  // 			   -sqr(Constants::pi)/3.-4.5))
  //     +qPDF            *CFfact/xp_*(3.+2.*xp_-(1.+xp_)*log(1.-xp_)
  // 				  -(1.+sqr(xp_))/(1.-xp_)*log(xp_))
  //     +(qPDF-xp_*loPDF)*CFfact/xp_*(2.*log(1.-xp_)/(1.-xp_)-1.5/(1.-xp_));
  //   double wgt = (f2g+f2q)/loPDF;

  return contrib_ == 1 ? max(0., wgt) : max(0., -wgt);
}
