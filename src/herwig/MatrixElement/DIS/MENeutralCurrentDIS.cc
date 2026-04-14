// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MENeutralCurrentDIS class.
//

#include <cmath>
#include <limits>
#include <string>
#include "MENeutralCurrentDIS.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/StandardMatchers.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Handlers/StandardXComb.h"
// #include "ThePEG/PDF/PolarizedBeamParticleData.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

namespace {

void ensureRealEmissionSpinInfo(const PPtr & part, bool incomingLeg) {
  if (!part || part->spinInfo()) return;

  const tcPDPtr data = part->dataPtr();
  if (!data) return;

  if (data->iSpin() == PDT::Spin1Half) {
    if (incomingLeg) {
      if (part->id() > 0) {
        vector<SpinorWaveFunction> waves;
        SpinorWaveFunction::calculateWaveFunctions(waves, part, incoming);
        SpinorWaveFunction::constructSpinInfo(waves, part, incoming, false);
      } else {
        vector<SpinorBarWaveFunction> waves;
        SpinorBarWaveFunction::calculateWaveFunctions(waves, part, incoming);
        SpinorBarWaveFunction::constructSpinInfo(waves, part, incoming, false);
      }
    } else {
      if (part->id() > 0) {
        vector<SpinorBarWaveFunction> waves;
        SpinorBarWaveFunction::calculateWaveFunctions(waves, part, outgoing);
        SpinorBarWaveFunction::constructSpinInfo(waves, part, outgoing, true);
      } else {
        vector<SpinorWaveFunction> waves;
        SpinorWaveFunction::calculateWaveFunctions(waves, part, outgoing);
        SpinorWaveFunction::constructSpinInfo(waves, part, outgoing, true);
      }
    }
  } else if (data->iSpin() == PDT::Spin1) {
    const bool massless =
      part->id() == ParticleID::g || part->id() == ParticleID::gamma ||
      part->id() == ParticleID::darkg;
    vector<LorentzPolarizationVector> waves;
    VectorWaveFunction::calculateWaveFunctions(waves, part,
                                               incomingLeg ? incoming : outgoing,
                                               massless, vector_phase);
    VectorWaveFunction::constructSpinInfo(waves, part,
                                          incomingLeg ? incoming : outgoing,
                                          !incomingLeg, massless);
  }
}

RhoDMatrix longitudinalRhoMatrix(tcPDPtr data, double pol) {
  pol = std::max(-1.0, std::min(1.0, pol));
  RhoDMatrix rho(data->iSpin());
  if (data->iSpin() == PDT::Spin1Half) {
    const unsigned int imax = rho.iSpin() - 1;
    rho(0,0) = 0.5 * (1.0 - pol);
    rho(imax,imax) = 0.5 * (1.0 + pol);
    rho(0,imax) = 0.0;
    rho(imax,0) = 0.0;
  } else if (data->iSpin() == PDT::Spin1) {
    rho(0,0) = 0.5 * (1.0 - pol);
    rho(1,1) = 0.0;
    rho(2,2) = 0.5 * (1.0 + pol);
    rho(0,1) = rho(1,0) = 0.0;
    rho(0,2) = rho(2,0) = 0.0;
    rho(1,2) = rho(2,1) = 0.0;
  }
  return rho;
}

struct FermionLineWaves {
  vector<SpinorWaveFunction> fermion;
  vector<SpinorBarWaveFunction> antifermion;
  bool incomingIsFermion;
};

FermionLineWaves buildFermionLineWaves(const PPtr & incomingPart,
                                       const PPtr & outgoingPart) {
  FermionLineWaves line;
  line.incomingIsFermion = incomingPart->id() > 0;

  SpinorWaveFunction fermionWF;
  SpinorBarWaveFunction antifermionWF;

  if (line.incomingIsFermion) {
    fermionWF = SpinorWaveFunction(incomingPart->momentum(),
                                   incomingPart->dataPtr(), incoming);
    antifermionWF = SpinorBarWaveFunction(outgoingPart->momentum(),
                                          outgoingPart->dataPtr(), outgoing);
  } else {
    fermionWF = SpinorWaveFunction(outgoingPart->momentum(),
                                   outgoingPart->dataPtr(), outgoing);
    antifermionWF = SpinorBarWaveFunction(incomingPart->momentum(),
                                          incomingPart->dataPtr(), incoming);
  }

  for (unsigned int ih = 0; ih < 2; ++ih) {
    fermionWF.reset(ih);
    antifermionWF.reset(ih);
    line.fermion.push_back(fermionWF);
    line.antifermion.push_back(antifermionWF);
  }

  return line;
}

vector<VectorWaveFunction> buildMasslessVectorWaves(const PPtr & part,
                                                    ThePEG::Helicity::Direction dir) {
  vector<VectorWaveFunction> waves;
  VectorWaveFunction proto(part->momentum(), part->dataPtr(), dir);
  proto.reset(0);
  waves.push_back(proto);
  proto.reset(2);
  waves.push_back(proto);
  return waves;
}

struct RealEmissionLegs {
  PPtr lin;
  PPtr pin;
  PPtr lout;
  PPtr out1;
  PPtr out2;
  std::string process;
  std::string out1Role;
  std::string out2Role;
};

bool collectRealEmissionLegs(const RealEmissionProcessPtr & proc,
                             bool isCompton,
                             RealEmissionLegs & legs) {
  if (!proc) return false;

  vector<PPtr> coloured;
  for (const auto & part : proc->incoming()) {
    if (ThePEG::LeptonMatcher::Check(part->data())) legs.lin = part;
    else legs.pin = part;
  }
  for (const auto & part : proc->outgoing()) {
    if (ThePEG::LeptonMatcher::Check(part->data())) legs.lout = part;
    else coloured.push_back(part);
  }
  if (!legs.lin || !legs.pin || !legs.lout) return false;

  if (isCompton) {
    PPtr qout, gout;
    for (const auto & part : coloured) {
      if (part->id() == ParticleID::g) gout = part;
      else qout = part;
    }
    if (!qout || !gout) return false;
    legs.out1 = qout;
    legs.out2 = gout;
    legs.process = "QCDC";
    legs.out1Role = "qout";
    legs.out2Role = "gout";
  } else {
    PPtr qout, qbout;
    for (const auto & part : coloured) {
      if (part->id() > 0) qout = part;
      else qbout = part;
    }
    if (!qout || !qbout) return false;
    legs.out1 = qout;
    legs.out2 = qbout;
    legs.process = "BGF";
    legs.out1Role = "qout";
    legs.out2Role = "qbout";
  }
  return true;
}

struct MatrixDiagnostics {
  double traceRe;
  double traceIm;
  double maxAntiHerm;
  double minDiag;
  double maxDiag;
  double maxDiagIm;
  double d0;
  double d1;
  double d2;
  bool finite;
  bool normalized;
  bool hermitian;
  bool sensible;
};

MatrixDiagnostics inspectMatrix(const RhoDMatrix & rho) {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  MatrixDiagnostics out{nan, nan, nan, nan, nan, nan, nan, nan, nan,
                        true, false, false, false};
  const unsigned int dim = std::abs(int(rho.iSpin()));
  if (dim == 0) {
    out.finite = false;
    return out;
  }

  Complex trace = 0.;
  out.maxAntiHerm = 0.0;
  out.minDiag = std::numeric_limits<double>::infinity();
  out.maxDiag = -std::numeric_limits<double>::infinity();
  out.maxDiagIm = 0.0;

  for (unsigned int i = 0; i < dim; ++i) {
    const Complex diag = rho(i, i);
    if (!std::isfinite(diag.real()) || !std::isfinite(diag.imag())) out.finite = false;
    trace += diag;
    out.minDiag = std::min(out.minDiag, diag.real());
    out.maxDiag = std::max(out.maxDiag, diag.real());
    out.maxDiagIm = std::max(out.maxDiagIm, std::abs(diag.imag()));
    if (i == 0) out.d0 = diag.real();
    else if (i == 1) out.d1 = diag.real();
    else if (i == 2) out.d2 = diag.real();
    for (unsigned int j = i + 1; j < dim; ++j) {
      const Complex aij = rho(i, j);
      const Complex aji = rho(j, i);
      if (!std::isfinite(aij.real()) || !std::isfinite(aij.imag()) ||
          !std::isfinite(aji.real()) || !std::isfinite(aji.imag())) {
        out.finite = false;
      }
      out.maxAntiHerm = std::max(out.maxAntiHerm, std::abs(aij - conj(aji)));
    }
  }

  out.traceRe = trace.real();
  out.traceIm = trace.imag();
  out.normalized = out.finite && std::abs(out.traceRe - 1.0) < 1e-8 &&
                   std::abs(out.traceIm) < 1e-8;
  out.hermitian = out.finite && out.maxAntiHerm < 1e-8 && out.maxDiagIm < 1e-8;
  out.sensible = out.normalized && out.hermitian &&
                 out.minDiag > -1e-8 && out.maxDiag < 1.0 + 1e-8;
  return out;
}

}

MENeutralCurrentDIS::MENeutralCurrentDIS() 
  : _minflavour(1), _maxflavour(5), _gammaZ(0),
    _useFiniteWidthSpacelikeZPropagator(false),
    _sinW(0.), _cosW(0.), _mz2(ZERO) {
  vector<unsigned int> mopt(2,0);
  mopt[0] = 1;
  massOption(mopt);
}

void MENeutralCurrentDIS::doinit() {
  DISBase::doinit();
  _z0    = getParticleData(ThePEG::ParticleID::Z0);
  _gamma = getParticleData(ThePEG::ParticleID::gamma);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() 
    << "Must be the Herwig StandardModel class in "
    << "MENeutralCurrentDIS::doinit" << Exception::abortnow;
  // vertices
  _theFFZVertex = hwsm->vertexFFZ();
  _theFFPVertex = hwsm->vertexFFP();
  _theFFGVertex = hwsm->vertexFFG();
  // electroweak parameters
  _sinW = generator()->standardModel()->sin2ThetaW();
  _cosW = sqrt(1.-_sinW);
  _sinW = sqrt(_sinW);
  _mz2 = sqr(_z0->mass());
#if 0
  // Debug: print coupling values during initialization.
  double alphaEM_val = generator()->standardModel()->alphaEM();
  generator()->log() << "MENeutralCurrentDIS: alphaEM (fixed) = " << alphaEM_val << "\n";
  double alphaS_val = generator()->standardModel()->alphaS();
  generator()->log() << "MENeutralCurrentDIS: alphaS (fixed) = " << alphaS_val << "\n";
#endif
}

void MENeutralCurrentDIS::getDiagrams() const {
  // which intermediates to include
  bool gamma = _gammaZ==0 || _gammaZ==1;
  bool Z0    = _gammaZ==0 || _gammaZ==2;
  // create the diagrams
  for(int ix=11;ix<=14;++ix) {
    for(unsigned int iz=0;iz<2;++iz) {
      tPDPtr lep = getParticleData(ix);
      if(iz==1) lep = lep->CC();
      for(int iy=_minflavour;iy<=_maxflavour;++iy) {
	tPDPtr quark = getParticleData(iy);
	// lepton quark scattering via gamma and Z
	if(gamma) add(new_ptr((Tree2toNDiagram(3), lep, _gamma, quark,
			       1, lep, 2, quark, -1)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), lep, _z0   , quark,
			       1, lep, 2, quark, -2)));
	// lepton antiquark scattering via gamma and Z
	quark = quark->CC();
	if(gamma) add(new_ptr((Tree2toNDiagram(3), lep, _gamma, quark,
			       1, lep, 2, quark, -3)));
	if(Z0)    add(new_ptr((Tree2toNDiagram(3), lep, _z0   , quark,
			       1, lep, 2, quark, -4)));
      }
    }
  }
}

unsigned int MENeutralCurrentDIS::orderInAlphaS() const {
  return 0;
}

unsigned int MENeutralCurrentDIS::orderInAlphaEW() const {
  return 2;
}

Selector<const ColourLines *>
MENeutralCurrentDIS::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c1("3 5");
  static ColourLines c2("-3 -5");
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 || diag->id() == -2 )
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MENeutralCurrentDIS::persistentOutput(PersistentOStream & os) const {
  os << _minflavour << _maxflavour << _gammaZ << _theFFZVertex << _theFFPVertex
     << _theFFGVertex
     << _gamma << _z0 << _sinW << _cosW << ounit(_mz2,GeV2)
     << _useFiniteWidthSpacelikeZPropagator;
}

void MENeutralCurrentDIS::persistentInput(PersistentIStream & is, int version) {
  is >> _minflavour >> _maxflavour >> _gammaZ >> _theFFZVertex >> _theFFPVertex
     >> _theFFGVertex
     >> _gamma >> _z0 >> _sinW >> _cosW >> iunit(_mz2,GeV2) ;
  if(version == 0) _useFiniteWidthSpacelikeZPropagator = false;
  else is >> _useFiniteWidthSpacelikeZPropagator;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MENeutralCurrentDIS,DISBase>
describeHerwigMENeutralCurrentDIS("Herwig::MENeutralCurrentDIS", "HwMEDIS.so", 1);

void MENeutralCurrentDIS::Init() {

  static ClassDocumentation<MENeutralCurrentDIS> documentation
    ("The MENeutralCurrentDIS class implements the matrix elements for leading-order "
     "neutral current deep inelastic scattering.");

  static Parameter<MENeutralCurrentDIS,int> interfaceMaxFlavour
    ("MaxFlavour",
     "The highest incoming quark flavour this matrix element is allowed to handle",
     &MENeutralCurrentDIS::_maxflavour, 5, 1, 5,
     false, false, Interface::limited);

  static Parameter<MENeutralCurrentDIS,int> interfaceMinFlavour
    ("MinFlavour",
     "The lightest incoming quark flavour this matrix element is allowed to handle",
     &MENeutralCurrentDIS::_minflavour, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<MENeutralCurrentDIS,unsigned int> interfaceGammaZ
    ("GammaZ",
     "Which terms to include",
     &MENeutralCurrentDIS::_gammaZ, 0, false, false);
  static SwitchOption interfaceGammaZAll
    (interfaceGammaZ,
     "All",
     "Include both gamma and Z terms",
     0);
  static SwitchOption interfaceGammaZGamma
    (interfaceGammaZ,
     "Gamma",
     "Only include the photon",
     1);
  static SwitchOption interfaceGammaZZ
    (interfaceGammaZ,
     "Z",
     "Only include the Z",
     2);

  static Switch<MENeutralCurrentDIS,bool>
    interfaceUseFiniteWidthSpacelikeZPropagator
    ("UseFiniteWidthSpacelikeZPropagator",
     "Whether to keep a finite width for spacelike Z exchange in the "
     "helicity-amplitude path. This affects the Born and real-emission DIS "
     "neutral-current amplitude construction only. The default No preserves "
     "the legacy zero-width spacelike Z behavior, while Yes uses the finite-"
     "width ThePEG propagator option iopt=7.",
     &MENeutralCurrentDIS::_useFiniteWidthSpacelikeZPropagator,
     false, false, false);
  static SwitchOption interfaceUseFiniteWidthSpacelikeZPropagatorYes
    (interfaceUseFiniteWidthSpacelikeZPropagator,
     "Yes",
     "Use the finite-width spacelike Z propagator in the helicity-amplitude path.",
     true);
  static SwitchOption interfaceUseFiniteWidthSpacelikeZPropagatorNo
    (interfaceUseFiniteWidthSpacelikeZPropagator,
     "No",
     "Use the legacy zero-width spacelike Z propagator in the helicity-amplitude path.",
     false);
}

Selector<MEBase::DiagramIndex>
MENeutralCurrentDIS::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if      ( diags[i]->id() == -1 || diags[i]->id() == -3 ) sel.insert(meInfo()[0], i);
    else if ( diags[i]->id() == -2 || diags[i]->id() == -4 ) sel.insert(meInfo()[1], i);
  }
  return sel;
}

double MENeutralCurrentDIS::helicityME(const pair<RhoDMatrix,RhoDMatrix> & rhoin ,
                                       vector<SpinorWaveFunction>    & f1,
				       vector<SpinorWaveFunction>    & f2,
				       vector<SpinorBarWaveFunction> & a1,
				       vector<SpinorBarWaveFunction> & a2,
				       bool lorder, bool qorder,
				       bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew (PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement gamma (PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  ProductionMatrixElement Zboson(PDT::Spin1Half,PDT::Spin1Half,
				 PDT::Spin1Half,PDT::Spin1Half);
  // which intermediates to include
  bool gam = _gammaZ==0 || _gammaZ==1;
  bool Z0  = _gammaZ==0 || _gammaZ==2;
  // declare the variables we need
  VectorWaveFunction inter[2];
  Complex diag1,diag2;
  // sum over helicities to get the matrix element
  unsigned int hel[4];
  unsigned int lhel1,lhel2,qhel1,qhel2;
  for(lhel1=0;lhel1<2;++lhel1) {
    for(lhel2=0;lhel2<2;++lhel2) {
      // intermediate for photon
      if(gam) inter[0]=_theFFPVertex->evaluate(mb2,1,_gamma,f1[lhel1],a1[lhel2]);
      // intermediate for Z
      if(Z0)    inter[1]=_theFFZVertex->evaluate(mb2,zHelicityPropagatorOption(),
                                                 _z0,f1[lhel1],a1[lhel2]);
      for(qhel1=0;qhel1<2;++qhel1) {
	for(qhel2=0;qhel2<2;++qhel2) {
	  hel[0] = lhel1;
	  hel[1] = qhel1;
	  hel[2] = lhel2;
	  hel[3] = qhel2;
	  if(!lorder) swap(hel[0],hel[2]);
	  if(!qorder) swap(hel[1],hel[3]);
	  // first the photon exchange diagram
	  diag1 = gam ?
	    _theFFPVertex->evaluate(mb2,f2[qhel1],a2[qhel2],inter[0]) : 0.;
	  // then the Z exchange diagram
	  diag2 = Z0 ?
	    _theFFZVertex->evaluate(mb2,f2[qhel1],a2[qhel2],inter[1]) : 0.;
	  // add up squares of individual terms
	  gamma (hel[0],hel[1],hel[2],hel[3]) = diag1;
	  Zboson(hel[0],hel[1],hel[2],hel[3]) = diag2;
	  // the full thing including interference
	  menew(hel[0],hel[1],hel[2],hel[3]) = diag1+diag2;
	}
      }
    }
  }
  // spin and colour factor (=1)
  double me[3]={menew .average(rhoin.first,rhoin.second),
                gamma .average(rhoin.first,rhoin.second),
                Zboson.average(rhoin.first,rhoin.second)};
  // info for diagram selection
  DVector save;
  save.push_back(me[1]);
  save.push_back(me[2]);
  meInfo(save);
  if(calc) _me.reset(menew);
  // analytic expression for testing
//   double test = 8.*sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(mb2))*
//     sqr(double(mePartonData()[1]->iCharge())/3.)/sqr(tHat())
//     *(sqr(sHat())+sqr(uHat())+4.*sqr(mePartonData()[0]->mass())*tHat())/4.;
//   cerr << "testing me " << me[0]/test << "\n";
  return me[0];
}

double MENeutralCurrentDIS::me2ForPolarizations(double Pl, double Pq) const {
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  bool lorder,qorder;
  SpinorWaveFunction    l1,q1;
  SpinorBarWaveFunction l2,q2;
  // lepton wave functions
  if(mePartonData()[0]->id()>0) {
    lorder=true;
    l1 = SpinorWaveFunction   (meMomenta()[0],mePartonData()[0],incoming);
    l2 = SpinorBarWaveFunction(meMomenta()[2],mePartonData()[2],outgoing);
  }
  else {
    lorder=false;
    l1 = SpinorWaveFunction   (meMomenta()[2],mePartonData()[2],outgoing);
    l2 = SpinorBarWaveFunction(meMomenta()[0],mePartonData()[0],incoming);
  }
  // quark wave functions
  if(mePartonData()[1]->id()>0) {
    qorder = true;
    q1 = SpinorWaveFunction   (meMomenta()[1],mePartonData()[1],incoming);
    q2 = SpinorBarWaveFunction(meMomenta()[3],mePartonData()[3],outgoing);
  }
  else {
    qorder = false;
    q1 = SpinorWaveFunction   (meMomenta()[3],mePartonData()[3],outgoing);
    q2 = SpinorBarWaveFunction(meMomenta()[1],mePartonData()[1],incoming);
  }
  // wavefunctions for various helicities
  for(unsigned int ix=0;ix<2;++ix) {
    l1.reset(ix); f1.push_back(l1);
    l2.reset(ix); a1.push_back(l2);
    q1.reset(ix); f2.push_back(q1);
    q2.reset(ix); a2.push_back(q2);
  }
  pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
  Pl = std::max(-1.0, std::min(1.0, Pl));
  Pq = std::max(-1.0, std::min(1.0, Pq));

  rhoin.first = RhoDMatrix(mePartonData()[0]->iSpin());
  const unsigned int lmax = rhoin.first.iSpin() - 1;
  rhoin.first(0,0) = 0.5 * (1.0 - Pl);
  rhoin.first(lmax,lmax) = 0.5 * (1.0 + Pl);
  rhoin.first(0,lmax) = 0.0;
  rhoin.first(lmax,0) = 0.0;

  rhoin.second = RhoDMatrix(mePartonData()[1]->iSpin());
  const unsigned int imax = rhoin.second.iSpin() - 1;
  rhoin.second(0,0) = 0.5 * (1.0 - Pq);
  rhoin.second(imax,imax) = 0.5 * (1.0 + Pq);
  rhoin.second(0,imax) = 0.0;
  rhoin.second(imax,0) = 0.0;
  return helicityME(rhoin,f1,f2,a1,a2,lorder,qorder,false);
}

double MENeutralCurrentDIS::me2ForPartonPolarization(double Pq) const {
  const std::pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
  return me2ForPolarizations(longPol(rhoin.first), Pq);
}

double MENeutralCurrentDIS::me2() const {
  const std::pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
  return me2ForPolarizations(longPol(rhoin.first), longPol(rhoin.second));
}

bool MENeutralCurrentDIS::pureLOGammaPointAuditChannel() const {
  return gammaZOption() == 1;
}

void MENeutralCurrentDIS::constructVertex(tSubProPtr sub) {
  // extract the particles in the hard process
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
  // sort out the ordering
  unsigned int order[4]={0,1,2,3};
  bool lorder(true),qorder(true);
  if(abs(hard[0]->id())<6) swap(hard[0],hard[1]);
  if(abs(hard[2]->id())<6) swap(hard[2],hard[3]);
  if(hard[0]->id()<0) {
    swap(order[0],order[2]);
    lorder = false;
  }
  if(hard[1]->id()<0) {
    swap(order[1],order[3]);
    qorder = false;
  }
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  SpinorWaveFunction   (f1,hard[order[0]], lorder ? incoming : outgoing, !lorder,true);
  SpinorWaveFunction   (f2,hard[order[1]], qorder ? incoming : outgoing, !qorder,true);
  SpinorBarWaveFunction(a1,hard[order[2]], lorder ? outgoing : incoming,  lorder,true);
  SpinorBarWaveFunction(a2,hard[order[3]], qorder ? outgoing : incoming,  qorder,true);
  pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
  helicityME(rhoin,f1,f2,a1,a2,lorder,qorder,true);
  // construct the vertex
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  // set the matrix element for the vertex
  hardvertex->ME(_me);
  // set the pointers and to and from the vertex
  hard[order[0]]->spinInfo()->rhoMatrix(rhoin.first );
  hard[order[1]]->spinInfo()->rhoMatrix(rhoin.second);
  for(unsigned int ix=0;ix<4;++ix)
    hard[ix]->spinInfo()->productionVertex(hardvertex);
}

void MENeutralCurrentDIS::constructRealEmissionSpinVertex(RealEmissionProcessPtr proc,
                                                          bool isCompton) const {
  RealEmissionLegs legs;
  if (!collectRealEmissionLegs(proc, isCompton, legs)) return;

  ProductionMatrixElement prodme;
  const Energy2 q2 = -(legs.lin->momentum() - legs.lout->momentum()).m2();

  if (isCompton) {
    prodme = qcdcRealEmissionME(legs.lin, legs.pin, legs.lout,
                                legs.out1, legs.out2, q2);
  } else {
    prodme = bgfRealEmissionME(legs.lin, legs.pin, legs.lout,
                               legs.out1, legs.out2, q2);
  }

  ensureRealEmissionSpinInfo(legs.lin, true);
  ensureRealEmissionSpinInfo(legs.pin, true);
  ensureRealEmissionSpinInfo(legs.lout, false);
  ensureRealEmissionSpinInfo(legs.out1, false);
  ensureRealEmissionSpinInfo(legs.out2, false);

  if (!legs.lin->spinInfo() || !legs.pin->spinInfo() || !legs.lout->spinInfo() ||
      !legs.out1->spinInfo() || !legs.out2->spinInfo()) return;

  RhoDMatrix leptonRho;
  for (const auto & part : proc->bornIncoming()) {
    if (ThePEG::LeptonMatcher::Check(part->data()) && part->spinInfo()) {
      leptonRho = part->spinInfo()->rhoMatrix();
      break;
    }
  }
  if (leptonRho.iSpin() == 0) {
    leptonRho = correctedLongitudinalRhoMatrices().first;
  }

  double xMapped = 0.0;
  if (proc->incoming().size() == 2) {
    xMapped = (proc->incoming()[0] == legs.pin) ? proc->x().first : proc->x().second;
  }
  const double mappedPol =
    mappedIncomingLongitudinalPolarization(legs.pin->dataPtr(), xMapped, q2);

  legs.lin->spinInfo()->rhoMatrix(leptonRho);
  legs.pin->spinInfo()->rhoMatrix(longitudinalRhoMatrix(legs.pin->dataPtr(), mappedPol));

  HardVertexPtr hardvertex = new_ptr(HardVertex());
  hardvertex->ME(prodme);

  legs.lin->spinInfo()->productionVertex(hardvertex);
  legs.pin->spinInfo()->productionVertex(hardvertex);
  legs.lout->spinInfo()->productionVertex(hardvertex);
  legs.out1->spinInfo()->productionVertex(hardvertex);
  legs.out2->spinInfo()->productionVertex(hardvertex);
}

void MENeutralCurrentDIS::diagnoseRealEmissionSpinState(RealEmissionProcessPtr proc,
                                                        bool isCompton) const {
  unsigned long diagIndex = 0;
  if (!nextPOWHEGRealSpinDiagnosticSlot(diagIndex)) return;

  RealEmissionLegs legs;
  if (!collectRealEmissionLegs(proc, isCompton, legs)) {
    generator()->log() << "POWHEG_SPIN_EVENT"
                       << " event=" << diagIndex
                       << " proc=" << (isCompton ? "QCDC" : "BGF")
                       << " enabled=" << (usePOWHEGRealSpinVertex() ? 1 : 0)
                       << " parsed=0\n";
    return;
  }

  struct LegView {
    const char * role;
    PPtr part;
    bool incoming;
  };
  const LegView views[] = {
    {"lin",  legs.lin,  true},
    {"pin",  legs.pin,  true},
    {"lout", legs.lout, false},
    {legs.out1Role.c_str(), legs.out1, false},
    {legs.out2Role.c_str(), legs.out2, false}
  };

  bool allSpin = true;
  bool allVertex = true;
  bool allHardVertex = true;
  bool allSensible = true;

  for (const auto & view : views) {
    const tSpinPtr spin = view.part ? view.part->spinInfo() : tSpinPtr();
    const bool hasSpin = bool(spin);
    const tcVertexPtr vertex = hasSpin ? spin->productionVertex() : tcVertexPtr();
    const bool hasVertex = bool(vertex);
    const bool vertexIsHard = bool(ThePEG::dynamic_ptr_cast<tcHardVertexPtr>(vertex));
    const int prodLoc = hasSpin ? spin->productionLocation() : -1;

    allSpin = allSpin && hasSpin;
    allVertex = allVertex && hasVertex;
    allHardVertex = allHardVertex && vertexIsHard;

    RhoDMatrix rho;
    if (hasSpin) {
      if (view.incoming) rho = spin->rhoMatrix();
      else if (hasVertex && prodLoc >= 0) rho = vertex->getRhoMatrix(prodLoc, true);
    }
    MatrixDiagnostics md = inspectMatrix(rho);
    allSensible = allSensible && md.sensible;

    generator()->log() << "POWHEG_SPIN_LEG"
                       << " event=" << diagIndex
                       << " proc=" << legs.process
                       << " role=" << view.role
                       << " incoming=" << (view.incoming ? 1 : 0)
                       << " id=" << (view.part ? view.part->id() : 0)
                       << " hasSpin=" << (hasSpin ? 1 : 0)
                       << " hasVertex=" << (hasVertex ? 1 : 0)
                       << " vertexHard=" << (vertexIsHard ? 1 : 0)
                       << " loc=" << prodLoc
                       << " source=" << (view.incoming ? "stored" : "vertex")
                       << " trace=" << md.traceRe
                       << " traceIm=" << md.traceIm
                       << " antiHerm=" << md.maxAntiHerm
                       << " minDiag=" << md.minDiag
                       << " maxDiag=" << md.maxDiag
                       << " maxDiagIm=" << md.maxDiagIm
                       << " d0=" << md.d0
                       << " d1=" << md.d1
                       << " d2=" << md.d2
                       << " finite=" << (md.finite ? 1 : 0)
                       << " normalized=" << (md.normalized ? 1 : 0)
                       << " hermitian=" << (md.hermitian ? 1 : 0)
                       << " sensible=" << (md.sensible ? 1 : 0)
                       << "\n";
  }

  generator()->log() << "POWHEG_SPIN_EVENT"
                     << " event=" << diagIndex
                     << " proc=" << legs.process
                     << " enabled=" << (usePOWHEGRealSpinVertex() ? 1 : 0)
                     << " parsed=1"
                     << " allSpin=" << (allSpin ? 1 : 0)
                     << " allVertex=" << (allVertex ? 1 : 0)
                     << " allHardVertex=" << (allHardVertex ? 1 : 0)
                     << " allSensible=" << (allSensible ? 1 : 0)
                     << "\n";
}

ProductionMatrixElement MENeutralCurrentDIS::qcdcRealEmissionME(PPtr lin, PPtr qin,
                                                                PPtr lout, PPtr qout,
                                                                PPtr gout,
                                                                Energy2 q2) const {
  const FermionLineWaves leptonLine = buildFermionLineWaves(lin, lout);
  const FermionLineWaves quarkLine = buildFermionLineWaves(qin, qout);
  const vector<VectorWaveFunction> gluonWaves =
    buildMasslessVectorWaves(gout, outgoing);

  ProductionMatrixElement prodme(PDT::Spin1Half, PDT::Spin1Half,
                                 PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1);

  for (unsigned int lhelF = 0; lhelF < 2; ++lhelF) {
    for (unsigned int lhelA = 0; lhelA < 2; ++lhelA) {
      const unsigned int helInL = leptonLine.incomingIsFermion ? lhelF : lhelA;
      const unsigned int helOutL = leptonLine.incomingIsFermion ? lhelA : lhelF;
      for (unsigned int qhelF = 0; qhelF < 2; ++qhelF) {
        for (unsigned int qhelA = 0; qhelA < 2; ++qhelA) {
          const unsigned int helInQ = quarkLine.incomingIsFermion ? qhelF : qhelA;
          const unsigned int helOutQ = quarkLine.incomingIsFermion ? qhelA : qhelF;

          for (unsigned int ghel = 0; ghel < gluonWaves.size(); ++ghel) {
            Complex amp = 0.;
            if (_gammaZ == 0 || _gammaZ == 1) {
              VectorWaveFunction inter =
                _theFFPVertex->evaluate(q2, 1, _gamma,
                                        leptonLine.fermion[lhelF],
                                        leptonLine.antifermion[lhelA]);
              SpinorWaveFunction off1 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.fermion[qhelF].particle()->CC(),
                                        quarkLine.fermion[qhelF],
                                        gluonWaves[ghel]);
              SpinorBarWaveFunction off2 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.antifermion[qhelA].particle()->CC(),
                                        quarkLine.antifermion[qhelA],
                                        gluonWaves[ghel]);
              amp += _theFFPVertex->evaluate(q2, off1,
                                             quarkLine.antifermion[qhelA], inter);
              amp += _theFFPVertex->evaluate(q2, quarkLine.fermion[qhelF],
                                             off2, inter);
            }
            if (_gammaZ == 0 || _gammaZ == 2) {
              VectorWaveFunction inter =
                _theFFZVertex->evaluate(q2, zHelicityPropagatorOption(), _z0,
                                        leptonLine.fermion[lhelF],
                                        leptonLine.antifermion[lhelA]);
              SpinorWaveFunction off1 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.fermion[qhelF].particle()->CC(),
                                        quarkLine.fermion[qhelF],
                                        gluonWaves[ghel]);
              SpinorBarWaveFunction off2 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.antifermion[qhelA].particle()->CC(),
                                        quarkLine.antifermion[qhelA],
                                        gluonWaves[ghel]);
              amp += _theFFZVertex->evaluate(q2, off1,
                                             quarkLine.antifermion[qhelA], inter);
              amp += _theFFZVertex->evaluate(q2, quarkLine.fermion[qhelF],
                                             off2, inter);
            }
            prodme(helInL, helInQ, helOutL, helOutQ, 2 * ghel) = amp;
          }
        }
      }
    }
  }

  return prodme;
}

ProductionMatrixElement MENeutralCurrentDIS::bgfRealEmissionME(PPtr lin, PPtr gin,
                                                               PPtr lout, PPtr qout,
                                                               PPtr qbout,
                                                               Energy2 q2) const {
  const FermionLineWaves leptonLine = buildFermionLineWaves(lin, lout);
  FermionLineWaves quarkLine;
  quarkLine.incomingIsFermion = false;
  SpinorWaveFunction qWF(qout->momentum(), qout->dataPtr(), outgoing);
  SpinorBarWaveFunction qbWF(qbout->momentum(), qbout->dataPtr(), outgoing);
  for (unsigned int ih = 0; ih < 2; ++ih) {
    qWF.reset(ih);
    qbWF.reset(ih);
    quarkLine.fermion.push_back(qWF);
    quarkLine.antifermion.push_back(qbWF);
  }
  const vector<VectorWaveFunction> gluonWaves =
    buildMasslessVectorWaves(gin, incoming);

  ProductionMatrixElement prodme(PDT::Spin1Half, PDT::Spin1,
                                 PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1Half);

  for (unsigned int lhelF = 0; lhelF < 2; ++lhelF) {
    for (unsigned int lhelA = 0; lhelA < 2; ++lhelA) {
      const unsigned int helInL = leptonLine.incomingIsFermion ? lhelF : lhelA;
      const unsigned int helOutL = leptonLine.incomingIsFermion ? lhelA : lhelF;
      for (unsigned int qhelF = 0; qhelF < 2; ++qhelF) {
        for (unsigned int qhelA = 0; qhelA < 2; ++qhelA) {
          for (unsigned int ghel = 0; ghel < gluonWaves.size(); ++ghel) {
            Complex amp = 0.;
            if (_gammaZ == 0 || _gammaZ == 1) {
              VectorWaveFunction inter =
                _theFFPVertex->evaluate(q2, 1, _gamma,
                                        leptonLine.fermion[lhelF],
                                        leptonLine.antifermion[lhelA]);
              SpinorWaveFunction off1 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.fermion[qhelF].particle()->CC(),
                                        quarkLine.fermion[qhelF],
                                        gluonWaves[ghel]);
              SpinorBarWaveFunction off2 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.antifermion[qhelA].particle()->CC(),
                                        quarkLine.antifermion[qhelA],
                                        gluonWaves[ghel]);
              amp += _theFFPVertex->evaluate(q2, off1,
                                             quarkLine.antifermion[qhelA], inter);
              amp += _theFFPVertex->evaluate(q2, quarkLine.fermion[qhelF],
                                             off2, inter);
            }
            if (_gammaZ == 0 || _gammaZ == 2) {
              VectorWaveFunction inter =
                _theFFZVertex->evaluate(q2, zHelicityPropagatorOption(), _z0,
                                        leptonLine.fermion[lhelF],
                                        leptonLine.antifermion[lhelA]);
              SpinorWaveFunction off1 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.fermion[qhelF].particle()->CC(),
                                        quarkLine.fermion[qhelF],
                                        gluonWaves[ghel]);
              SpinorBarWaveFunction off2 =
                _theFFGVertex->evaluate(q2, 5,
                                        quarkLine.antifermion[qhelA].particle()->CC(),
                                        quarkLine.antifermion[qhelA],
                                        gluonWaves[ghel]);
              amp += _theFFZVertex->evaluate(q2, off1,
                                             quarkLine.antifermion[qhelA], inter);
              amp += _theFFZVertex->evaluate(q2, quarkLine.fermion[qhelF],
                                             off2, inter);
            }
            prodme(helInL, 2 * ghel, helOutL, qhelF, qhelA) = amp;
          }
        }
      }
    }
  }

  return prodme;
}

double MENeutralCurrentDIS::A(tcPDPtr lin, tcPDPtr,
			      tcPDPtr qin, tcPDPtr, Energy2 q2) const {
  // photon only 
  if(_gammaZ==1) return 0.;
  const NCCoefficients coeff = ncCoefficients(lin, qin, q2);
  if (coeff.D0 == 0.0) return 0.0;
  return coeff.N0 / coeff.D0;
}

double MENeutralCurrentDIS::A_pol(tcPDPtr lin, tcPDPtr,
                                  tcPDPtr qin, tcPDPtr,
                                  Energy2 q2, double Pl, double Pq) const {

  // Photon-only: pure vector exchange gives A_gamma = 2*Pl*Pq.
  // There is no extra q/qbar sign flip in this branch; the R2/R3 difference
  // in BGF is handled by passing +/-Pg_m from DISBase.
  if (_gammaZ == 1) {
    return 2.0 * Pl * Pq;
  }

  const NCCoefficients coeff = ncCoefficients(lin, qin, q2);
  double num = coeff.N0 + Pl * coeff.Nl + Pq * coeff.Nq + Pl * Pq * coeff.Nlq;
  double den = coeff.D0 + Pl * coeff.Dl + Pq * coeff.Dq + Pl * Pq * coeff.Dlq;

  if (den == 0.0) return 0.0;
  return num / den;
}

DISBase::CollinearBlendWeights
MENeutralCurrentDIS::collinearBlendWeights(tcPDPtr lin, tcPDPtr,
                                           tcPDPtr qin, tcPDPtr,
                                           Energy2 q2, double Pl, double Pq,
                                           double ell) const {
  // Photon-only keeps the legacy scalar blend, which is exact in that limit.
  if (_gammaZ == 1) {
    return DISBase::collinearBlendWeights(lin, tcPDPtr(), qin, tcPDPtr(),
                                          q2, Pl, Pq, ell);
  }
  const NCCoefficients coeff = ncCoefficients(lin, qin, q2);
  const double D_even = coeff.D0 + Pl * coeff.Dl;
  const double D_spin = Pq * (coeff.Dq + Pl * coeff.Dlq);
  const double N_even = coeff.N0 + Pl * coeff.Nl;
  const double N_spin = Pq * (coeff.Nq + Pl * coeff.Nlq);

  const double sigmaBorn =
    (1.0 + sqr(ell)) * (D_even + D_spin) +
    ell * (N_even + N_spin);

  if (std::abs(sigmaBorn) <= 1e-30) {
    return {1.0, 0.0, 1.0, 0.0};
  }

  const double qUnpolarized =
    ((1.0 + sqr(ell)) * D_even + ell * N_even) / sigmaBorn;
  const double qPolarized =
    ((1.0 + sqr(ell)) * D_spin + ell * N_spin) / sigmaBorn;

  // At this order the gluon collinear counterterm only contributes to the
  // F2-like unpolarized channel and the G2-like spin-dependent channel.
  const double gUnpolarized =
    ((1.0 + sqr(ell)) * D_even) / sigmaBorn;
  const double gPolarized =
    (ell * N_spin) / sigmaBorn;

  return {qUnpolarized, qPolarized, gUnpolarized, gPolarized};
}

double MENeutralCurrentDIS::qcdcMappedDenominatorRatio(tcPDPtr lin, tcPDPtr,
                                                       tcPDPtr qin, tcPDPtr,
                                                       Energy2 q2, double Pl,
                                                       double PqBorn,
                                                       double PqMapped) const {
  if (_gammaZ == 1) return 1.0;
  const NCCoefficients coeff = ncCoefficients(lin, qin, q2);
  const double dBorn =
    coeff.D0 + Pl * coeff.Dl + PqBorn * coeff.Dq + Pl * PqBorn * coeff.Dlq;
  if (std::abs(dBorn) <= 1e-30) return 1.0;

  const double dMapped =
    coeff.D0 + Pl * coeff.Dl + PqMapped * coeff.Dq + Pl * PqMapped * coeff.Dlq;
  return dMapped / dBorn;
}

double MENeutralCurrentDIS::sigmaBornFactor(tcPDPtr lin, tcPDPtr qin, Energy2 q2,
                                            double Pl, double Pq,
                                            double ell) const {
  if (_gammaZ == 1) {
    return 1.0 + 2.0 * Pl * Pq * ell + sqr(ell);
  }
  const NCCoefficients coeff = ncCoefficients(lin, qin, q2);
  const double D_even = coeff.D0 + Pl * coeff.Dl;
  const double D_spin = Pq * (coeff.Dq + Pl * coeff.Dlq);
  const double N_even = coeff.N0 + Pl * coeff.Nl;
  const double N_spin = Pq * (coeff.Nq + Pl * coeff.Nlq);

  return (1.0 + sqr(ell)) * (D_even + D_spin) + ell * (N_even + N_spin);
}

MENeutralCurrentDIS::NCCoefficients
MENeutralCurrentDIS::ncCoefficients(tcPDPtr lin, tcPDPtr qin, Energy2 q2) const {
  NCCoefficients out{};
  const double etaL = (lin->id() < 0) ? -1.0 : 1.0;
  const double etaQ = (qin->id() < 0) ? -1.0 : 1.0;

  double Ql, Qq, CVl, CAl, CVq, CAq;

  if (std::abs(lin->id()) % 2 == 0) {
    Ql  = (_gammaZ == 0) ? generator()->standardModel()->enu() : 0.0;
    CVl = 0.25 * generator()->standardModel()->vnu();
    CAl = 0.25 * generator()->standardModel()->anu();
  } else {
    Ql  = (_gammaZ == 0) ? generator()->standardModel()->ee()  : 0.0;
    CVl = 0.25 * generator()->standardModel()->ve();
    CAl = 0.25 * generator()->standardModel()->ae();
  }

  if (std::abs(qin->id()) % 2 == 0) {
    Qq  = (_gammaZ == 0) ? generator()->standardModel()->eu() : 0.0;
    CVq = 0.25 * generator()->standardModel()->vu();
    CAq = 0.25 * generator()->standardModel()->au();
  } else {
    Qq  = (_gammaZ == 0) ? generator()->standardModel()->ed() : 0.0;
    CVq = 0.25 * generator()->standardModel()->vd();
    CAq = 0.25 * generator()->standardModel()->ad();
  }

  const double Ql2  = Ql * Ql;
  const double Qq2  = Qq * Qq;
  const double CVl2 = CVl * CVl;
  const double CAl2 = CAl * CAl;
  const double CVq2 = CVq * CVq;
  const double CAq2 = CAq * CAq;
  const double CVlCVq = CVl * CVq;
  const double CAlCAq = CAl * CAq;
  const double k = 1.0 / sqr(_sinW * _cosW);

  double zInt = 0.0;
  double zSq = 0.0;
  if (_gammaZ == 0 || _gammaZ == 2) {
    const Energy2 gammaZ = _z0->mass() * _z0->width();
    const auto den = sqr(q2 + _mz2) + sqr(gammaZ);
    const double rInt = double(q2 * (q2 + _mz2) / den);
    const double rSq = double(sqr(q2) / den);
    zInt = k * rInt;
    zSq = sqr(k) * rSq;
  }

  out.N0  = etaL * etaQ * 4.0 * CAlCAq * (zInt * Ql * Qq + 2.0 * zSq * CVlCVq);
  out.Nl  = etaQ * -4.0 * CAq * (zSq * CVq * (CAl2 + CVl2) + zInt * CVl * Ql * Qq);
  out.Nq  = etaL * -4.0 * CAl * (zSq * CVl * (CAq2 + CVq2) + zInt * CVq * Ql * Qq);
  out.Nlq = 2.0 * (zSq * (CAq2 + CVq2) * (CAl2 + CVl2)
                + 2.0 * zInt * CVlCVq * Ql * Qq
                + Ql2 * Qq2);

  out.D0  = Ql2 * Qq2
         + 2.0 * zInt * Ql * Qq * CVlCVq
         + zSq * (CVl2 + CAl2) * (CVq2 + CAq2);
  out.Dl  = etaL * -2.0 * CAl * (zSq * CVl * (CAq2 + CVq2) + zInt * CVq * Ql * Qq);
  out.Dq  = etaQ * -2.0 * CAq * (zSq * CVq * (CAl2 + CVl2) + zInt * CVl * Ql * Qq);
  out.Dlq = etaL * etaQ * 2.0 * CAlCAq * (2.0 * zSq * CVlCVq + zInt * Ql * Qq);

  return out;
}

bool MENeutralCurrentDIS::bornClosureDiagnostics(double Pl,
                                                 double PqBorn,
                                                 double PqMapped,
                                                 double ell,
                                                 BornClosureDiagnostics &out) const {
  if (!mePartonData()[0] || !mePartonData()[1]) return false;
  const Energy2 q2 = -(meMomenta()[0] - meMomenta()[2]).m2();
  out = BornClosureDiagnostics{};

  out.sigmaBorn = sigmaBornFactor(mePartonData()[0], mePartonData()[1], q2,
                                  Pl, PqBorn, ell);
  out.sigmaZero = sigmaBornFactor(mePartonData()[0], mePartonData()[1], q2,
                                  Pl, 0.0, ell);
  out.sigmaMapped = sigmaBornFactor(mePartonData()[0], mePartonData()[1], q2,
                                    Pl, PqMapped, ell);
  out.me2Born = me2ForPartonPolarization(PqBorn);
  out.me2Zero = me2ForPartonPolarization(0.0);
  out.me2Mapped = me2ForPartonPolarization(PqMapped);

  const double me00 = me2ForPolarizations(0.0, 0.0);
  const double me10 = me2ForPolarizations(1.0, 0.0);
  const double me01 = me2ForPolarizations(0.0, 1.0);
  const double me11 = me2ForPolarizations(1.0, 1.0);

  const NCCoefficients coeff = ncCoefficients(mePartonData()[0], mePartonData()[1], q2);
  const double pred00 = (1.0 + sqr(ell)) * coeff.D0 + ell * coeff.N0;
  const double pred10 = pred00 + (1.0 + sqr(ell)) * coeff.Dl + ell * coeff.Nl;
  const double pred01 = pred00 + (1.0 + sqr(ell)) * coeff.Dq + ell * coeff.Nq;
  const double pred11 = pred00
                      + (1.0 + sqr(ell)) * (coeff.Dl + coeff.Dq + coeff.Dlq)
                      + ell * (coeff.Nl + coeff.Nq + coeff.Nlq);

  if (std::abs(me00) > 1e-30 && std::abs(pred00) > 1e-30) {
    out.coeffScale = me00 / pred00;
    out.coeffPlMe = (me10 - me00) / me00;
    out.coeffPqMe = (me01 - me00) / me00;
    out.coeffPlPqMe = (me11 - me10 - me01 + me00) / me00;

    out.coeffPlPred = (pred10 - pred00) / pred00;
    out.coeffPqPred = (pred01 - pred00) / pred00;
    out.coeffPlPqPred = (pred11 - pred10 - pred01 + pred00) / pred00;
  }
  return true;
}
