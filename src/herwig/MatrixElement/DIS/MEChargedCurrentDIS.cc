// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEChargedCurrentDIS class.
//

#include <cmath>
#include <limits>
#include <string>
#include "MEChargedCurrentDIS.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
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

MEChargedCurrentDIS::MEChargedCurrentDIS() 
  : _maxflavour(5), _massopt(0) {
  vector<unsigned int> mopt(2,1);
  mopt[1] = _massopt;
  massOption(mopt);
}

void MEChargedCurrentDIS::doinit() {
  DISBase::doinit();
  _wp = getParticleData(ThePEG::ParticleID::Wplus );
  _wm = getParticleData(ThePEG::ParticleID::Wminus);
  // cast the SM pointer to the Herwig SM pointer
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() 
    << "Must be the Herwig StandardModel class in "
    << "MEChargedCurrentDIS::doinit" << Exception::abortnow;
  // vertices
  _theFFWVertex = hwsm->vertexFFW();
}

AbstractFFVVertexPtr MEChargedCurrentDIS::gluonVertex() const {
  tcHwSMPtr hwsm = ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) {
    throw Exception()
      << "Must be the Herwig StandardModel class in "
      << "MEChargedCurrentDIS::gluonVertex" << Exception::abortnow;
  }
  return hwsm->vertexFFG();
}


void MEChargedCurrentDIS::getDiagrams() const {
  // possible quarks
  typedef std::vector<pair<long,long> > Pairvector;
  Pairvector quarkpair;
  quarkpair.reserve(6);
  // don't even think of putting 'break' in here!
  switch(_maxflavour) {
  case 6:
    quarkpair.push_back(make_pair(ParticleID::s, ParticleID::t));
    quarkpair.push_back(make_pair(ParticleID::d, ParticleID::t));
    quarkpair.push_back(make_pair(ParticleID::b, ParticleID::t));
    [[fallthrough]];
  case 5:
    quarkpair.push_back(make_pair(ParticleID::b, ParticleID::c));
    quarkpair.push_back(make_pair(ParticleID::b, ParticleID::u));
    [[fallthrough]];
  case 4:
    quarkpair.push_back(make_pair(ParticleID::s, ParticleID::c));
    quarkpair.push_back(make_pair(ParticleID::d, ParticleID::c));
    [[fallthrough]];
  case 3:
    quarkpair.push_back(make_pair(ParticleID::s, ParticleID::u));
    [[fallthrough]];
  case 2:
    quarkpair.push_back(make_pair(ParticleID::d, ParticleID::u));
    [[fallthrough]];
  default:
    ;
  }
  // create the diagrams
  for(int il1=11;il1<=14;++il1) {
    int il2 = il1%2==0 ? il1-1 : il1+1;
    for(unsigned int iz=0;iz<2;++iz) {
      tcPDPtr lepin  = iz==1 ? getParticleData(il1) : getParticleData(-il1);
      tcPDPtr lepout = iz==1 ? getParticleData(il2) : getParticleData(-il2);
      tcPDPtr inter  = lepin->iCharge()-lepout->iCharge()==3 ? _wp : _wm;
      for(unsigned int iq=0;iq<quarkpair.size();++iq) {
	tcPDPtr first  = getParticleData(quarkpair[iq].first );
	tcPDPtr second = getParticleData(quarkpair[iq].second);
	if(inter==_wp) {
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, first       , 
		       1, lepout, 2, second     , -1)));
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, second->CC(), 
		       1, lepout, 2, first->CC(), -2)));
	}
	else {
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, second     , 
	  	       1, lepout, 2, first       , -1)));
	  add(new_ptr((Tree2toNDiagram(3), lepin, inter, first->CC(), 
	  	       1, lepout, 2, second->CC(), -2)));
	}
      }
    }
  }
}

unsigned int MEChargedCurrentDIS::orderInAlphaS() const {
  return 0;
}

unsigned int MEChargedCurrentDIS::orderInAlphaEW() const {
  return 2;
}

Selector<const ColourLines *>
MEChargedCurrentDIS::colourGeometries(tcDiagPtr diag) const {
  static ColourLines c1("3 5");
  static ColourLines c2("-3 -5");
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 )
    sel.insert(1.0, &c1);
  else
    sel.insert(1.0, &c2);
  return sel;
}

void MEChargedCurrentDIS::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _maxflavour << _wp << _wm << _massopt;
}

void MEChargedCurrentDIS::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _maxflavour >> _wp >> _wm >> _massopt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEChargedCurrentDIS,DISBase>
describeHerwigMEChargedCurrentDIS("Herwig::MEChargedCurrentDIS", "HwMEDIS.so");

void MEChargedCurrentDIS::Init() {

  static ClassDocumentation<MEChargedCurrentDIS> documentation
    ("The MEChargedCurrentDIS class implements the matrix elements "
     "for leading-order charged current deep inelastic scattering");

  static Parameter<MEChargedCurrentDIS,unsigned int> interfaceMaxFlavour
    ( "MaxFlavour",
      "The heaviest incoming quark flavour this matrix element is allowed to handle "
      "(if applicable).",
      &MEChargedCurrentDIS::_maxflavour, 5, 2, 6, false, false, true);

  static Switch<MEChargedCurrentDIS,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the mass of the outgoing quarks",
     &MEChargedCurrentDIS::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionMassless
    (interfaceMassOption,
     "Massless",
     "Treat the outgoing quarks as massless",
     0);
  static SwitchOption interfaceMassOptionMassive
    (interfaceMassOption,
     "Massive",
     "Treat the outgoing quarks as massive",
     1);

}

Selector<MEBase::DiagramIndex>
MEChargedCurrentDIS::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1., i);
  return sel;
}

double MEChargedCurrentDIS::helicityME(const pair<RhoDMatrix,RhoDMatrix> & rhoin,
				       vector<SpinorWaveFunction>    & f1,
				       vector<SpinorWaveFunction>    & f2,
				       vector<SpinorBarWaveFunction> & a1,
				       vector<SpinorBarWaveFunction> & a2,
				       bool lorder, bool qorder, bool calc) const {
  // scale
  Energy2 mb2(scale());
  // matrix element to be stored
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  // pick a W boson
  tcPDPtr ipart = (mePartonData()[0]->iCharge()-mePartonData()[1]->iCharge())==3 ?
    _wp : _wm;
  // declare the variables we need
  VectorWaveFunction inter;
  double me(0.);
  Complex diag;
  // sum over helicities to get the matrix element
  unsigned int hel[4];
  unsigned int lhel1,lhel2,qhel1,qhel2;
  for(lhel1=0;lhel1<2;++lhel1) {
    for(lhel2=0;lhel2<2;++lhel2) {
      // intermediate W
      inter = _theFFWVertex->evaluate(mb2,3,ipart,f1[lhel1],a1[lhel2]);
      for(qhel1=0;qhel1<2;++qhel1) {
	for(qhel2=0;qhel2<2;++qhel2) {
	  hel[0] = lhel1;
	  hel[1] = qhel1;
	  hel[2] = lhel2;
	  hel[3] = qhel2;
	  if(!lorder) swap(hel[0],hel[2]);
	  if(!qorder) swap(hel[1],hel[3]);
	  diag = _theFFWVertex->evaluate(mb2,f2[qhel1],a2[qhel2],inter);
	  menew(hel[0],hel[1],hel[2],hel[3]) = diag;
	}
      }
    }
  }
  // Average over the incoming spin density matrices.
  me = menew.average(rhoin.first,rhoin.second);
  // tcPolarizedBeamPDPtr beam[2] = 
  //   {dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[0]),
  //    dynamic_ptr_cast<tcPolarizedBeamPDPtr>(mePartonData()[1])};
  // if( beam[0] || beam[1] ) {
  //   RhoDMatrix rho[2] = {beam[0] ? beam[0]->rhoMatrix() : RhoDMatrix(mePartonData()[0]->iSpin()),
  //       		 beam[1] ? beam[1]->rhoMatrix() : RhoDMatrix(mePartonData()[1]->iSpin())};
  //   me = menew.average(rho[0],rho[1]);
  // }
  if(calc) _me.reset(menew);
  return me;
}

double MEChargedCurrentDIS::me2() const {
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
  const pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
  return helicityME(rhoin,f1,f2,a1,a2,lorder,qorder,false);
}

void MEChargedCurrentDIS::constructVertex(tSubProPtr sub) {
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
  const pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
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

void MEChargedCurrentDIS::constructRealEmissionSpinVertex(RealEmissionProcessPtr proc,
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

void MEChargedCurrentDIS::diagnoseRealEmissionSpinState(RealEmissionProcessPtr proc,
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

ProductionMatrixElement MEChargedCurrentDIS::qcdcRealEmissionME(PPtr lin, PPtr qin,
                                                                PPtr lout, PPtr qout,
                                                                PPtr gout,
                                                                Energy2 q2) const {
  const AbstractFFVVertexPtr ffgVertex = gluonVertex();
  const FermionLineWaves leptonLine = buildFermionLineWaves(lin, lout);
  const FermionLineWaves quarkLine = buildFermionLineWaves(qin, qout);
  const vector<VectorWaveFunction> gluonWaves =
    buildMasslessVectorWaves(gout, outgoing);
  const tcPDPtr mediator = (lin->dataPtr()->iCharge() - lout->dataPtr()->iCharge() == 3) ?
    _wp : _wm;

  ProductionMatrixElement prodme(PDT::Spin1Half, PDT::Spin1Half,
                                 PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1);

  for (unsigned int lhelF = 0; lhelF < 2; ++lhelF) {
    for (unsigned int lhelA = 0; lhelA < 2; ++lhelA) {
      const unsigned int helInL = leptonLine.incomingIsFermion ? lhelF : lhelA;
      const unsigned int helOutL = leptonLine.incomingIsFermion ? lhelA : lhelF;
      VectorWaveFunction inter =
        _theFFWVertex->evaluate(q2, 1, mediator,
                                leptonLine.fermion[lhelF],
                                leptonLine.antifermion[lhelA]);

      for (unsigned int qhelF = 0; qhelF < 2; ++qhelF) {
        for (unsigned int qhelA = 0; qhelA < 2; ++qhelA) {
          const unsigned int helInQ = quarkLine.incomingIsFermion ? qhelF : qhelA;
          const unsigned int helOutQ = quarkLine.incomingIsFermion ? qhelA : qhelF;

          for (unsigned int ghel = 0; ghel < gluonWaves.size(); ++ghel) {
            const SpinorWaveFunction off1 =
              ffgVertex->evaluate(q2, 5,
                                  quarkLine.fermion[qhelF].particle()->CC(),
                                  quarkLine.fermion[qhelF],
                                  gluonWaves[ghel]);
            const SpinorBarWaveFunction off2 =
              ffgVertex->evaluate(q2, 5,
                                  quarkLine.antifermion[qhelA].particle()->CC(),
                                  quarkLine.antifermion[qhelA],
                                  gluonWaves[ghel]);

            Complex amp = 0.;
            amp += _theFFWVertex->evaluate(q2, off1,
                                           quarkLine.antifermion[qhelA], inter);
            amp += _theFFWVertex->evaluate(q2, quarkLine.fermion[qhelF],
                                           off2, inter);
            prodme(helInL, helInQ, helOutL, helOutQ, 2 * ghel) = amp;
          }
        }
      }
    }
  }

  return prodme;
}

ProductionMatrixElement MEChargedCurrentDIS::bgfRealEmissionME(PPtr lin, PPtr gin,
                                                               PPtr lout, PPtr qout,
                                                               PPtr qbout,
                                                               Energy2 q2) const {
  const AbstractFFVVertexPtr ffgVertex = gluonVertex();
  const FermionLineWaves leptonLine = buildFermionLineWaves(lin, lout);
  const vector<VectorWaveFunction> gluonWaves =
    buildMasslessVectorWaves(gin, incoming);
  const tcPDPtr mediator = (lin->dataPtr()->iCharge() - lout->dataPtr()->iCharge() == 3) ?
    _wp : _wm;

  vector<SpinorWaveFunction> qbWaves;
  vector<SpinorBarWaveFunction> qWaves;
  SpinorWaveFunction qbWF(qbout->momentum(), qbout->dataPtr(), outgoing);
  SpinorBarWaveFunction qWF(qout->momentum(), qout->dataPtr(), outgoing);
  for (unsigned int ih = 0; ih < 2; ++ih) {
    qbWF.reset(ih);
    qWF.reset(ih);
    qbWaves.push_back(qbWF);
    qWaves.push_back(qWF);
  }

  ProductionMatrixElement prodme(PDT::Spin1Half, PDT::Spin1,
                                 PDT::Spin1Half, PDT::Spin1Half, PDT::Spin1Half);

  for (unsigned int lhelF = 0; lhelF < 2; ++lhelF) {
    for (unsigned int lhelA = 0; lhelA < 2; ++lhelA) {
      const unsigned int helInL = leptonLine.incomingIsFermion ? lhelF : lhelA;
      const unsigned int helOutL = leptonLine.incomingIsFermion ? lhelA : lhelF;
      VectorWaveFunction inter =
        _theFFWVertex->evaluate(q2, 1, mediator,
                                leptonLine.fermion[lhelF],
                                leptonLine.antifermion[lhelA]);

      for (unsigned int qhelQ = 0; qhelQ < 2; ++qhelQ) {
        for (unsigned int qhelQB = 0; qhelQB < 2; ++qhelQB) {
          for (unsigned int ghel = 0; ghel < gluonWaves.size(); ++ghel) {
            const SpinorWaveFunction off1 =
              ffgVertex->evaluate(q2, 5,
                                  qbWaves[qhelQB].particle()->CC(),
                                  qbWaves[qhelQB],
                                  gluonWaves[ghel]);
            const SpinorBarWaveFunction off2 =
              ffgVertex->evaluate(q2, 5,
                                  qWaves[qhelQ].particle()->CC(),
                                  qWaves[qhelQ],
                                  gluonWaves[ghel]);

            Complex amp = 0.;
            amp += _theFFWVertex->evaluate(q2, off1, qWaves[qhelQ], inter);
            amp += _theFFWVertex->evaluate(q2, qbWaves[qhelQB], off2, inter);
            prodme(helInL, 2 * ghel, helOutL, qhelQ, qhelQB) = amp;
          }
        }
      }
    }
  }

  return prodme;
}

double MEChargedCurrentDIS::A(tcPDPtr lin, tcPDPtr,
			      tcPDPtr qin, tcPDPtr, Energy2) const {
  double output = 2.;
  if(qin->id()<0) output *= -1.;
  if(lin->id()<0) output *= -1;
  return output;
}

double MEChargedCurrentDIS::A_pol(tcPDPtr lin, tcPDPtr lout,
                                  tcPDPtr qin, tcPDPtr qout,
                                  Energy2 scale, double, double) const {
  return A(lin,lout,qin,qout,scale);
}

DISBase::CollinearBlendWeights
MEChargedCurrentDIS::collinearBlendWeights(tcPDPtr, tcPDPtr,
                                           tcPDPtr qin, tcPDPtr,
                                           Energy2, double, double Pq,
                                           double) const {
  const double denom = ccHadronSpinFactor(qin, Pq);
  if (std::abs(denom) <= 1e-30) {
    return {1.0, 0.0, 1.0, 0.0};
  }

  const double etaQ = (qin->id() < 0) ? -1.0 : 1.0;
  const double unpolarized = 1.0 / denom;
  const double polarized = (-etaQ * Pq) / denom;
  return {unpolarized, polarized, unpolarized, polarized};
}

double MEChargedCurrentDIS::qcdcMappedDenominatorRatio(tcPDPtr, tcPDPtr,
                                                       tcPDPtr qin, tcPDPtr,
                                                       Energy2, double,
                                                       double PqBorn,
                                                       double PqMapped) const {
  const double born = ccHadronSpinFactor(qin, PqBorn);
  if (std::abs(born) <= 1e-30) return 1.0;
  const double mapped = ccHadronSpinFactor(qin, PqMapped);
  return mapped / born;
}

double MEChargedCurrentDIS::realEmissionDenominatorFactor(tcPDPtr, tcPDPtr,
                                                          tcPDPtr qin, tcPDPtr,
                                                          Energy2, double,
                                                          double Pq) const {
  return ccHadronSpinFactor(qin, Pq);
}

bool MEChargedCurrentDIS::useMappedPolarizedEmissionKernel() const {
  return true;
}

double MEChargedCurrentDIS::ccHadronSpinFactor(tcPDPtr qin, double Pq) const {
  const double etaQ = (qin->id() < 0) ? -1.0 : 1.0;
  return 1.0 - etaQ * Pq;
}
