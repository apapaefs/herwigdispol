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

using namespace Herwig;
using namespace ThePEG::Helicity;

namespace {

// Ensure realised POWHEG 2->3 legs carry spin information before the
// spin-correlation HardVertex is attached.
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

struct FermionLineWaves {
  vector<SpinorWaveFunction> fermion;
  vector<SpinorBarWaveFunction> antifermion;
  bool incomingIsFermion;
};

// Build helicity bases for a fermion line, independent of whether the incoming
// particle is the fermion or antifermion end of the line.
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

// Construct the two physical helicity states for an external massless vector.
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

// Extract the lepton, incoming coloured parton, and outgoing coloured legs from
// a realised POWHEG event in the ordering expected by the exact spin ME.
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

}

// Configure charged-current defaults. The outgoing-quark mass option is kept
// configurable for heavy-flavour channels such as top production.
MEChargedCurrentDIS::MEChargedCurrentDIS() 
  : _maxflavour(5), _massopt(0),
    _useFiniteWidthSpacelikeWPropagator(false) {
  vector<unsigned int> mopt(2,1);
  mopt[1] = _massopt;
  massOption(mopt);
}

// Resolve the charged-current particle data and the W vertex after repository
// setup, when the StandardModel instance is available.
void MEChargedCurrentDIS::doinit() {
  DISBase::doinit();
  _wp = getParticleData(ThePEG::ParticleID::Wplus );
  _wm = getParticleData(ThePEG::ParticleID::Wminus);
  tcHwSMPtr hwsm=ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() 
    << "Must be the Herwig StandardModel class in "
    << "MEChargedCurrentDIS::doinit" << Exception::abortnow;
  _theFFWVertex = hwsm->vertexFFW();
}

// Lazily access the QCD quark-gluon vertex used only by exact 2->3 spin MEs.
AbstractFFVVertexPtr MEChargedCurrentDIS::gluonVertex() const {
  tcHwSMPtr hwsm = ThePEG::dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if (!hwsm) {
    throw Exception()
      << "Must be the Herwig StandardModel class in "
      << "MEChargedCurrentDIS::gluonVertex" << Exception::abortnow;
  }
  return hwsm->vertexFFG();
}


// Register charged-current quark transitions allowed by MaxFlavour, including
// charge-conjugate channels for W+ and W- exchange.
void MEChargedCurrentDIS::getDiagrams() const {
  typedef std::vector<pair<long,long> > Pairvector;
  Pairvector quarkpair;
  quarkpair.reserve(6);
  // Intentional fall-through accumulates all lighter quark doublets.
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

// The Born charged-current hard process has no explicit powers of alphaS.
unsigned int MEChargedCurrentDIS::orderInAlphaS() const {
  return 0;
}

// Charged-current DIS is an electroweak 2->2 process at Born level.
unsigned int MEChargedCurrentDIS::orderInAlphaEW() const {
  return 2;
}

// Provide the colour flow for quark or antiquark charged-current scattering.
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

// Persist the charged-current vertex, flavour bound, W data, mass option, and
// finite-width W option used by the helicity-amplitude path.
void MEChargedCurrentDIS::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _maxflavour << _wp << _wm << _massopt
     << _useFiniteWidthSpacelikeWPropagator;
}

// Read the charged-current persistent state, defaulting the finite-width W
// switch for files written before that release field existed.
void MEChargedCurrentDIS::persistentInput(PersistentIStream & is, int version) {
  is >> _theFFWVertex >> _maxflavour >> _wp >> _wm >> _massopt;
  if(version == 0) _useFiniteWidthSpacelikeWPropagator = false;
  else is >> _useFiniteWidthSpacelikeWPropagator;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEChargedCurrentDIS,DISBase>
describeHerwigMEChargedCurrentDIS("Herwig::MEChargedCurrentDIS", "HwMEDIS.so", 1);

// Register user-facing charged-current controls.
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

  static Switch<MEChargedCurrentDIS,bool>
    interfaceUseFiniteWidthSpacelikeWPropagator
    ("UseFiniteWidthSpacelikeWPropagator",
     "Whether to keep a finite width for spacelike W exchange in the "
     "helicity-amplitude path. This affects the Born and real-emission DIS "
     "charged-current amplitude construction only. The default No preserves "
     "the legacy zero-width spacelike W behavior, while Yes uses the finite-"
     "width ThePEG propagator option iopt=7.",
     &MEChargedCurrentDIS::_useFiniteWidthSpacelikeWPropagator,
     false, false, false);
  static SwitchOption interfaceUseFiniteWidthSpacelikeWPropagatorYes
    (interfaceUseFiniteWidthSpacelikeWPropagator,
     "Yes",
     "Use the finite-width spacelike W propagator in the helicity-amplitude path.",
     true);
  static SwitchOption interfaceUseFiniteWidthSpacelikeWPropagatorNo
    (interfaceUseFiniteWidthSpacelikeWPropagator,
     "No",
     "Use the legacy zero-width spacelike W propagator in the helicity-amplitude path.",
     false);

}

// All charged-current diagrams have equal selector weight once kinematics pass.
Selector<MEBase::DiagramIndex>
MEChargedCurrentDIS::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) sel.insert(1., i);
  return sel;
}

// Build the W-exchange helicity amplitude for the current Born point, using
// the supplied incoming spin-density matrices for longitudinal polarization.
double MEChargedCurrentDIS::helicityME(const pair<RhoDMatrix,RhoDMatrix> & rhoin,
				       vector<SpinorWaveFunction>    & f1,
				       vector<SpinorWaveFunction>    & f2,
				       vector<SpinorBarWaveFunction> & a1,
				       vector<SpinorBarWaveFunction> & a2,
				       bool lorder, bool qorder, bool calc) const {
  Energy2 mb2(scale());
  ProductionMatrixElement menew(PDT::Spin1Half,PDT::Spin1Half,
				PDT::Spin1Half,PDT::Spin1Half);
  tcPDPtr ipart = (mePartonData()[0]->iCharge()-mePartonData()[1]->iCharge())==3 ?
    _wp : _wm;
  VectorWaveFunction inter;
  double me(0.);
  Complex diag;
  unsigned int hel[4];
  unsigned int lhel1,lhel2,qhel1,qhel2;
  for(lhel1=0;lhel1<2;++lhel1) {
    for(lhel2=0;lhel2<2;++lhel2) {
      inter = _theFFWVertex->evaluate(mb2,wHelicityPropagatorOption(3),
                                      ipart,f1[lhel1],a1[lhel2]);
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
  me = menew.average(rhoin.first,rhoin.second);
  if(calc) _me.reset(menew);
  return me;
}

// Standard Born ME entry point using the corrected event-local rho matrices.
double MEChargedCurrentDIS::me2() const {
  const pair<RhoDMatrix,RhoDMatrix> rhoin = correctedLongitudinalRhoMatrices();
  return me2ForPolarizations(longPol(rhoin.first), longPol(rhoin.second));
}

// Attach the Born spin-correlation vertex to the hard subprocess after
// normalizing the event-record ordering to lepton/quark fermion lines.
void MEChargedCurrentDIS::constructVertex(tSubProPtr sub) {
  ParticleVector hard;
  hard.push_back(sub->incoming().first);
  hard.push_back(sub->incoming().second);
  hard.push_back(sub->outgoing()[0]);
  hard.push_back(sub->outgoing()[1]);
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
  HardVertexPtr hardvertex=new_ptr(HardVertex());
  hardvertex->ME(_me);
  // Store the rho matrices used in the amplitude and point all hard legs at
  // the shared production vertex.
  hard[order[0]]->spinInfo()->rhoMatrix(rhoin.first );
  hard[order[1]]->spinInfo()->rhoMatrix(rhoin.second);
  for(unsigned int ix=0;ix<4;++ix)
    hard[ix]->spinInfo()->productionVertex(hardvertex);
}

// Build the exact spin-only HardVertex for a realised POWHEG real emission.
// Generation and weighting are already complete; this hook only records spin
// correlations on the accepted 2->3 event.
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
  const MappedIncomingSpinDensity mappedSpin =
    mappedIncomingSpinDensity(legs.pin->dataPtr(), xMapped, q2, "powheg-vertex");

  legs.lin->spinInfo()->rhoMatrix(leptonRho);
  legs.pin->spinInfo()->rhoMatrix(mappedSpin.rho);

  HardVertexPtr hardvertex = new_ptr(HardVertex());
  hardvertex->ME(prodme);

  legs.lin->spinInfo()->productionVertex(hardvertex);
  legs.pin->spinInfo()->productionVertex(hardvertex);
  legs.lout->spinInfo()->productionVertex(hardvertex);
  legs.out1->spinInfo()->productionVertex(hardvertex);
  legs.out2->spinInfo()->productionVertex(hardvertex);
}

// Exact W-exchange QCDC real-emission helicity amplitude for the POWHEG
// spin-correlation vertex, including both quark-line gluon-emission orderings.
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
        _theFFWVertex->evaluate(q2, wHelicityPropagatorOption(1), mediator,
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

// Exact W-exchange BGF real-emission helicity amplitude for the POWHEG
// spin-correlation vertex, with the incoming gluon in physical helicity states.
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
        _theFFWVertex->evaluate(q2, wHelicityPropagatorOption(1), mediator,
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

// Charged-current analyzing coefficient from the chiral W coupling and the
// lepton/quark charge-conjugation signs.
double MEChargedCurrentDIS::A(tcPDPtr lin, tcPDPtr,
			      tcPDPtr qin, tcPDPtr, Energy2) const {
  double output = 2.;
  if(qin->id()<0) output *= -1.;
  if(lin->id()<0) output *= -1;
  return output;
}

// The charged-current analyzing coefficient is independent of the incoming
// longitudinal polarization values once the chiral prefactor is handled.
double MEChargedCurrentDIS::A_pol(tcPDPtr lin, tcPDPtr lout,
                                  tcPDPtr qin, tcPDPtr qout,
                                  Energy2 scale, double, double) const {
  return A(lin,lout,qin,qout,scale);
}

// Split the charged-current hadron-side chiral prefactor into the pieces
// multiplying unpolarized and polarized NLO kernels.
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

// Correct QCDC denominators when the mapped incoming-parton polarization differs
// from the Born value.
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

// Return the charged-current chiral prefactor for the mapped real-emission
// denominator.
double MEChargedCurrentDIS::realEmissionDenominatorFactor(tcPDPtr, tcPDPtr,
                                                          tcPDPtr qin, tcPDPtr,
                                                          Energy2, double,
                                                          double Pq) const {
  return ccHadronSpinFactor(qin, Pq);
}

// Charged-current real-emission kernels follow the mapped x_B/x_p polarization
// because the chiral prefactor depends directly on the incoming parton spin.
bool MEChargedCurrentDIS::useMappedPolarizedEmissionKernel() const {
  return true;
}

// Hadron-side spin factor for V-A charged-current scattering.
double MEChargedCurrentDIS::ccHadronSpinFactor(tcPDPtr qin, double Pq) const {
  const double etaQ = (qin->id() < 0) ? -1.0 : 1.0;
  return 1.0 - etaQ * Pq;
}

// Re-evaluate the Born ME with explicit incoming lepton and parton
// polarizations for correlated weights and spin-factor checks.
double MEChargedCurrentDIS::me2ForPolarizations(double Pl, double Pq) const {
  vector<SpinorWaveFunction>    f1,f2;
  vector<SpinorBarWaveFunction> a1,a2;
  bool lorder,qorder;
  SpinorWaveFunction    l1,q1;
  SpinorBarWaveFunction l2,q2;
  // Normalize incoming/outgoing leptons and quarks to the fermion/antifermion
  // order expected by the helicity vertices.
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
  for(unsigned int ix=0;ix<2;++ix) {
    l1.reset(ix); f1.push_back(l1);
    l2.reset(ix); a1.push_back(l2);
    q1.reset(ix); f2.push_back(q1);
    q2.reset(ix); a2.push_back(q2);
  }
  const pair<RhoDMatrix,RhoDMatrix> rhoin =
    make_pair(longitudinalRhoMatrix(mePartonData()[0], Pl),
              longitudinalRhoMatrix(mePartonData()[1], Pq));
  return helicityME(rhoin,f1,f2,a1,a2,lorder,qorder,false);
}

double MEChargedCurrentDIS::rivetWeightBornME2(double Pl, double Pq) const {
  return me2ForPolarizations(Pl, Pq);
}
