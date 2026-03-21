#line 1 "./DMModel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMModel class.
//

#include "DMModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

DMModel::DMModel() : cDMmed_(0.), cSMmed_({1.0,1.0,1.0}) {}

IBPtr DMModel::clone() const {
  return new_ptr(*this);
}

IBPtr DMModel::fullclone() const {
  return new_ptr(*this);
}

void DMModel::persistentOutput(PersistentOStream & os) const {
  os << cDMmed_ << cSMmed_ << QQZpVertex_ << DMDMZpVertex_;
}

void DMModel::persistentInput(PersistentIStream & is, int) {
  is >> cDMmed_ >> cSMmed_ >> QQZpVertex_ >> DMDMZpVertex_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMModel,BSMModel>
describeHerwigDMModel("Herwig::DMModel", "HwDMModel.so");

void DMModel::Init() {

  static ClassDocumentation<DMModel> documentation
    ("The DMModel class is designed to implement a simple dark matter model"
     " with fermionic dark matter and a vector mediator, "
     "as described in  arXiv:1911.11147",
     "The DMModel class is designed to implement a simple dark matter model"
     " with fermionic dark matter and a vector mediator, "
     "as described in \\cite{Plehn:2019jeo}",
     "\\bibitem{Plehn:2019jeo}"
     "T.~Plehn, P.~Reimitz and P.~Richardson,"
     "%``Hadronic Footprint of GeV-Mass Dark Matter,''"
     "arXiv:1911.11147 [hep-ph]."
     "%%CITATION = ARXIV:1911.11147;%%");

  static Parameter<DMModel,double> interfacecDMmed
    ("cDMmed",
     "coupling of DM to dark mediator",
     &DMModel::cDMmed_, 1.0, 0., 10., false, false, Interface::limited);

  static ParVector<DMModel,double> interfacecSMmed
    ("cSMmed",
     "coupling of SM to dark mediator",
     &DMModel::cSMmed_, -1 , 1.0 , -10. , 10. , false, false, Interface::limited);

  static Reference<DMModel,AbstractFFVVertex> interfaceQQZpVertex
    ("Vertex/QQZpVertex",
     "The vertex coupling the quarks and the mediator",
     &DMModel::QQZpVertex_, false, false, true, false, false);

  static Reference<DMModel,AbstractFFVVertex> interfaceDMDMZpVertex
    ("Vertex/DMDMZpVertex",
     "The vertex coupling the DM to the mediator",
     &DMModel::DMDMZpVertex_, false, false, true, false, false);

}

void DMModel::doinit() {
  // additional vertices
  if(QQZpVertex_)   addVertex(QQZpVertex_);
  if(DMDMZpVertex_) addVertex(DMDMZpVertex_);
  BSMModel::doinit();
}
#line 1 "./DMMediatorQuarksVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMMediatorQuarksVertex class.
//

#include "DMMediatorQuarksVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DMModel.h"

using namespace Herwig;

DMMediatorQuarksVertex::DMMediatorQuarksVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr DMMediatorQuarksVertex::clone() const {
  return new_ptr(*this);
}

IBPtr DMMediatorQuarksVertex::fullclone() const {
  return new_ptr(*this);
}

void DMMediatorQuarksVertex::persistentOutput(PersistentOStream & os) const {
  os << cSMmed_;
}

void DMMediatorQuarksVertex::persistentInput(PersistentIStream & is, int) {
  is >> cSMmed_;
}

void DMMediatorQuarksVertex::setCoupling(Energy2 ,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  int iferm=abs(aa->id());
  assert(iferm>0 && iferm<4);
  norm(cSMmed_[iferm-1]);
  left(1.);
  right(1.);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMMediatorQuarksVertex,FFVVertex>
  describeHerwigDMMediatorQuarksVertex("Herwig::DMMediatorQuarksVertex", "HwDMModel.so");

void DMMediatorQuarksVertex::Init() {

  static ClassDocumentation<DMMediatorQuarksVertex> documentation
    ("The DMMediatorQuarksVertex class implements the coupling of the quarks to the mediator.");

}

void DMMediatorQuarksVertex::doinit() {
  cDMModelPtr model = dynamic_ptr_cast<cDMModelPtr>(generator()->standardModel());
  cSMmed_ = model->cSMmed();
  for(int ix=1;ix<4;++ix) {
    addToList(-ix, ix, 32);
  }
  FFVVertex::doinit();
}
#line 1 "./DMDMMediatorVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DMDMMediatorVertex class.
//

#include "DMDMMediatorVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DMModel.h"

using namespace Herwig;

DMDMMediatorVertex::DMDMMediatorVertex() : cDMmed_(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

IBPtr DMDMMediatorVertex::clone() const {
  return new_ptr(*this);
}

IBPtr DMDMMediatorVertex::fullclone() const {
  return new_ptr(*this);
}

void DMDMMediatorVertex::persistentOutput(PersistentOStream & os) const {
  os << cDMmed_;
}

void DMDMMediatorVertex::persistentInput(PersistentIStream & is, int) {
  is >> cDMmed_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DMDMMediatorVertex,FFVVertex>
describeHerwigDMDMMediatorVertex("Herwig::DMDMMediatorVertex", "HwDMModel.so");

void DMDMMediatorVertex::Init() {

  static ClassDocumentation<DMDMMediatorVertex> documentation
    ("The DMDMMediatorVertex class implements the couplnig of dark matter to the mediator.");

}

void DMDMMediatorVertex::setCoupling(Energy2 ,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  int iferm=abs(aa->id());
  assert(iferm==52);
  norm(cDMmed_);
  left(1.);
  right(1.);
}

void DMDMMediatorVertex::doinit() {
  cDMModelPtr model = dynamic_ptr_cast<cDMModelPtr>(generator()->standardModel());
  cDMmed_ = model->cDMmed();
  addToList(52, 52, 32);
  FFVVertex::doinit();
}
#line 1 "./MEDM2Mesons.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEDM2Mesons class.
//

#include "MEDM2Mesons.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Models/General/BSMModel.h"

using namespace Herwig;
typedef LorentzVector<complex<InvEnergy> > LorentzPolarizationVectorInvE;

MEDM2Mesons::MEDM2Mesons() : cDMmed_(0.), cSMmed_({0.0,0.0,0.0})
{}

Energy2 MEDM2Mesons::scale() const {
  return sHat();
}

unsigned int MEDM2Mesons::orderInAlphaS() const {
  return 0;
}

unsigned int MEDM2Mesons::orderInAlphaEW() const {
  return 0;
}

IBPtr MEDM2Mesons::clone() const {
  return new_ptr(*this);
}

IBPtr MEDM2Mesons::fullclone() const {
  return new_ptr(*this);
}

void MEDM2Mesons::doinit() {
  // make sure the current got initialised
  current_->init();
  // max energy
  Energy Emax = generator()->maximumCMEnergy();
  // loop over the modes
  int nmode=0;
  for(unsigned int imode=0;imode<current_->numberOfModes();++imode) {
     // get the external particles for this mode
     int iq(0),ia(0);
     tPDVector out = current_->particles(0,imode,iq,ia);
    current_->decayModeInfo(imode,iq,ia);
     if(iq==2&&ia==-2) continue;
     PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(incomingA_,out,1.,
						     incomingB_,Emax));
     PhaseSpaceChannel channel(mode);
     if(!current_->createMode(0,tcPDPtr(), FlavourInfo(), imode,mode,0,-1,channel,Emax)) continue;
     modeMap_[imode] = nmode;
     addMode(mode);
     ++nmode;
  }
  // cast the model
  Ptr<BSMModel>::ptr model = dynamic_ptr_cast<Ptr<BSMModel>::ptr>(generator()->standardModel());
  bool foundDM(false),foundU(false),foundD(false),foundS(false);
  // find the vertices we need and extract the couplings
  for(unsigned int ix = 0; ix < model->numberOfVertices(); ++ix ) {
    VertexBasePtr vertex = model->vertex(ix);
    if(vertex->getNpoint()!=3) continue;
    for(unsigned int iloc = 0;iloc < 3; ++iloc) {
      vector<long> ext = vertex->search(iloc, Mediator_->id());
      if(ext.empty()) continue;
      for(unsigned int ioff=0;ioff<ext.size();ioff+=3) {
	if(iloc!=2) assert(false);
	if((ext[ioff]==incomingA_->id() && ext[ioff+1]==incomingB_->id()) || 
	   (ext[ioff]==incomingB_->id() && ext[ioff+1]==incomingA_->id())) {
	  foundDM = true;
	  vertex->setCoupling(sqr(Emax),incomingA_,incomingB_,Mediator_);
	  cDMmed_ = vertex->norm();
	}
	else if(abs(ext[ioff])==1 && abs(ext[ioff+1])==1 &&  ext[ioff]==-ext[ioff+1]) {
	  foundD = true;
	  vertex->setCoupling(sqr(Emax),getParticleData(1),getParticleData(-1),Mediator_);
	  cSMmed_[0] = vertex->norm();
	}
	else if(abs(ext[ioff])==2 && abs(ext[ioff+1])==2 &&  ext[ioff]==-ext[ioff+1]) {
	  foundU = true;
	  vertex->setCoupling(sqr(Emax),getParticleData(2),getParticleData(-2),Mediator_);
	  cSMmed_[1] = vertex->norm();
	}
	else if(abs(ext[ioff])==3 && abs(ext[ioff+1])==3 &&  ext[ioff]==-ext[ioff+1]) {
	  foundS = true;
	  vertex->setCoupling(sqr(Emax),getParticleData(3),getParticleData(-3),Mediator_);
	  cSMmed_[2] = vertex->norm();
	}
      }
    }
  }
  if(!foundDM) {
    throw InitException() << "Cannot find DM coupling in MEDM2Mesons::doinit()";
  }
  if(!foundD) {
    throw InitException() << "Cannot find down quark coupling in MEDM2Mesons::doinit()";
  }
  if(!foundU) {
    throw InitException() << "Cannot find up quark coupling in MEDM2Mesons::doinit()";
  }
  if(!foundS) {
    throw InitException() << "Cannot find strange quark coupling in MEDM2Mesons::doinit()";
  }
  MEMultiChannel::doinit();
}

void MEDM2Mesons::persistentOutput(PersistentOStream & os) const {
  os << current_ << modeMap_ << incomingA_ << incomingB_ << Mediator_ << cDMmed_ << cSMmed_;
}

void MEDM2Mesons::persistentInput(PersistentIStream & is, int) {
  is >> current_ >> modeMap_ >> incomingA_ >> incomingB_ >> Mediator_ >> cDMmed_ >> cSMmed_;
}

//The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MEDM2Mesons,MEMultiChannel>
  describeHerwigMEDM2Mesons("Herwig::MEDM2Mesons", "Herwig.so");

void MEDM2Mesons::Init() {

  static ClassDocumentation<MEDM2Mesons> documentation
    ("The MEDM2Mesons class simulates the annhilation of"
     " DM particles to mesons at low energy");

  static Reference<MEDM2Mesons,WeakCurrent> interfaceWeakCurrent
    ("WeakCurrent",
     "The reference for the decay current to be used.",
     &MEDM2Mesons::current_, false, false, true, false, false);
  
  static Reference<MEDM2Mesons,ParticleData> interfaceIncomingA
    ("IncomingA",
     "First incoming particle",
     &MEDM2Mesons::incomingA_, false, false, true, false, false);
  
  static Reference<MEDM2Mesons,ParticleData> interfaceIncomingB
    ("IncomingB",
     "Second incoming particle",
     &MEDM2Mesons::incomingB_, false, false, true, false, false);

  static Reference<MEDM2Mesons,ParticleData> interfaceMediator_
    ("Mediator",
     "DM mediator",
     &MEDM2Mesons::Mediator_, false, false, true, false, false);

}


double MEDM2Mesons::me2(const int ichan) const {
  // compute the incoming current
  LorentzPolarizationVectorInvE lepton[2][2];
  if(incomingA_->iSpin()==PDT::Spin1Half && incomingB_->iSpin()==PDT::Spin1Half) {
    SpinorWaveFunction    em_in( meMomenta()[0],mePartonData()[0],incoming);
    SpinorBarWaveFunction ep_in( meMomenta()[1],mePartonData()[1],incoming);
    vector<SpinorWaveFunction> f1;
    vector<SpinorBarWaveFunction> a1;
    for(unsigned int ix=0;ix<2;++ix) {
      em_in.reset(ix);
      f1.push_back(em_in);
      ep_in.reset(ix);
      a1.push_back(ep_in);
    }
    // this should be coupling of DM to mediator/ mediator propagator
    complex<Energy> mmed = Mediator_->mass();
    complex<Energy2> mmed2 = sqr(mmed);
    complex<Energy> mwid = Mediator_->width();
    complex<Energy2> prop = sHat()-mmed2+Complex(0.,1.)*mmed*mwid;
    complex<InvEnergy2> pre = cDMmed_/prop;
    for(unsigned ix=0;ix<2;++ix) {
      for(unsigned iy=0;iy<2;++iy) {
	lepton[ix][iy]= pre*f1[ix].dimensionedWave().vectorCurrent(a1[iy].dimensionedWave());
      }
    }
  }
  // TODO think about other spins for the DM
  else
    assert(false);
  // work out the mapping for the hadron vector
  int nOut = int(meMomenta().size())-2;
  vector<unsigned int> constants(nOut+1);
  vector<PDT::Spin   > iSpin(nOut);
  vector<int> hadrons(nOut);
  int itemp(1);
  int ix(nOut);
  do {
    --ix;
    iSpin[ix]      = mePartonData()[ix+2]->iSpin();
    itemp         *= iSpin[ix];
    constants[ix]  = itemp;
    hadrons[ix]   = mePartonData()[ix+2]->id();
  }
  while(ix>0);
  constants[nOut] = 1;
  // calculate the matrix element
  me_.reset(ProductionMatrixElement(PDT::Spin1Half,PDT::Spin1Half,iSpin));
  // calculate the hadron current
  unsigned int imode = current_->decayMode(hadrons);
  Energy q = sqrt(sHat());
  vector<Lorentz5Momentum> momenta(meMomenta()   .begin()+2,   meMomenta().end());
  tPDVector out = mode(modeMap_.at(imode))->outgoing();
  if(ichan<0) iMode(modeMap_.at(imode));
  // get the hadronic currents for the I=1 and I=0 components
  vector<LorentzPolarizationVectorE> 
    hadronI0(current_->current(tcPDPtr(), FlavourInfo(IsoSpin::IZero, IsoSpin::I3Zero,Strangeness::Zero),
			       imode,ichan,q,out,momenta,DecayIntegrator::Calculate));
  vector<LorentzPolarizationVectorE> 
    hadronI1(current_->current(tcPDPtr(), FlavourInfo(IsoSpin::IOne, IsoSpin::I3Zero,Strangeness::Zero),
			       imode,ichan,q,out,momenta,DecayIntegrator::Calculate));
  vector<LorentzPolarizationVectorE> 
    hadronssbar(current_->current(tcPDPtr(), FlavourInfo(IsoSpin::IZero, IsoSpin::I3Zero,Strangeness::ssbar),
				  imode,ichan,q,out,momenta,DecayIntegrator::Calculate));
  // compute the matrix element
  vector<unsigned int> ihel(meMomenta().size());
  double output(0.);
  unsigned int hI0_size = hadronI0.size();
  unsigned int hI1_size = hadronI1.size();
  unsigned int hss_size = hadronssbar.size();
  unsigned int maxsize = max(max(hadronI0.size(),hadronI1.size()),hss_size);
  for(unsigned int hhel=0;hhel<maxsize;++hhel) {
    // map the index for the hadrons to a helicity state
    for(int ix=nOut;ix>0;--ix) {
      ihel[ix+1]=(hhel%constants[ix-1])/constants[ix];
    }
    // loop over the helicities of the incoming particles
    for(ihel[1]=0;ihel[1]<2;++ihel[1]){
      for(ihel[0]=0;ihel[0]<2;++ihel[0]) {
	Complex amp = 0.;
	// work on coefficients for the I1 and I0 bits
	if(hI0_size != 0 )
	  amp += Complex((cSMmed_[0]+cSMmed_[1])/sqrt(2.)*(lepton[ihel[0]][ihel[1]].dot(hadronI0[hhel])));
  	if(hI1_size !=0)
	  amp += Complex((cSMmed_[0]-cSMmed_[1])/sqrt(2.)*(lepton[ihel[0]][ihel[1]].dot(hadronI1[hhel])));
	if(hss_size !=0)
	  amp += Complex(cSMmed_[2]*                      (lepton[ihel[0]][ihel[1]].dot(hadronssbar[hhel])));
	me_(ihel)= amp;
  	output += std::norm(amp);
      }
    }
  }
  // symmetry factors
  map<long,int> ncount;
  double symmetry(1.);
  for(tPDPtr o : out) ncount[o->id()]+=1;
  for(map<long,int>::const_iterator it=ncount.begin();it!=ncount.end();++it) {
    symmetry *= it->second;
  }
  // prefactors
  output *= 0.25*sqr(pow(sqrt(sHat())/q,int(momenta.size()-2)));
  return output/symmetry;
}


void MEDM2Mesons::constructVertex(tSubProPtr) {
}
