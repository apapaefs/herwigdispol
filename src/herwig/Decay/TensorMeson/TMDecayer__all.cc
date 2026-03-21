#line 1 "./TensorMeson2PScalarDecayer.cc"
// -*- C++ -*-
//
// TensorMeson2PScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMeson2PScalarDecayer class.
//

#include "TensorMeson2PScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMeson2PScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMeson2PScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMeson2PScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMeson2PScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void TensorMeson2PScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV);
}

void TensorMeson2PScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMeson2PScalarDecayer,DecayIntegrator>
describeHerwigTensorMeson2PScalarDecayer("Herwig::TensorMeson2PScalarDecayer", "HwTMDecay.so");

void TensorMeson2PScalarDecayer::Init() {

  static ClassDocumentation<TensorMeson2PScalarDecayer> documentation
    ("The TensorMeson2PScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

  static Command<TensorMeson2PScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, scalars, coupling(1/GeV) and max weight for a decay",
     &TensorMeson2PScalarDecayer::setUpDecayMode, false);
  
  static Deleted<TensorMeson2PScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMeson2PScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in TensorMeson2PScalarDecayer have been deleted, please use SetUpDecayMode");
}

void TensorMeson2PScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

// matrix elememt for the process
double TensorMeson2PScalarDecayer::me2(const int,const Particle & part,
				       const tPDVector &,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin0,PDT::Spin0)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  // calculate the matrix element
  for(unsigned int ix=0;ix<5;++ix) {
    (*ME())(ix,0,0) = Complex(coupling_[imode()]/part.mass()*
			      ((tensors_[ix]*momenta[1])*momenta[0]));
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = Energy4(pow<4,1>(2*pcm))*sqr( coupling_[imode()]/part.mass())/120.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool TensorMeson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					       double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=7;
  return order;
}

void TensorMeson2PScalarDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMeson2PScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 0";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./TensorMesonVectorPScalarDecayer.cc"
// -*- C++ -*-
//
// TensorMesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorPScalarDecayer class.
//

#include "TensorMesonVectorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonVectorPScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMesonVectorPScalarDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void TensorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV2);
}

void TensorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonVectorPScalarDecayer,DecayIntegrator>
describeHerwigTensorMesonVectorPScalarDecayer("Herwig::TensorMesonVectorPScalarDecayer", "HwTMDecay.so");

void TensorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorPScalarDecayer> documentation
    ("The TensorMesonVectorPScalarDecayer class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static Command<TensorMesonVectorPScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, scalar, coupling(1/GeV2) and max weight for a decay",
     &TensorMesonVectorPScalarDecayer::setUpDecayMode, false);
  
  static Deleted<TensorMesonVectorPScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in TensorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorPScalarDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in TensorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorPScalarDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in TensorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorPScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in TensorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorPScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in TensorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");
}
void TensorMesonVectorPScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  VectorWaveFunction::constructSpinInfo(vectors_,decay[0],outgoing,true,
					decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double TensorMesonVectorPScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector &,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin0)));
  // check for photons
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  InvEnergy3 fact(coupling_[imode()]/part.mass());
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    for(unsigned int vhel=0;vhel<3;++vhel){
      if(vhel==1&&photon) (*ME())(inhel,vhel,0)=0.;
      else {
	LorentzVector<complex<InvEnergy> > vtemp=
	  fact*epsilon(momenta[0],vectors_[vhel],momenta[1]);
	(*ME())(inhel,vhel,0) = Complex((momenta[1]*tensors_[inhel]).dot(vtemp));
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = Energy4(pow<4,1>(2*pcm))*sqr( coupling_[imode()])/80.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool TensorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						   double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*sqr(dm.parent()->mass());
  mecode=8;
  return order;
}

void TensorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV2 << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMesonVectorPScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV2);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./TensorMesonVectorVectorDecayer.cc"
// -*- C++ -*-
//
// TensorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorVectorDecayer class.
//

#include "TensorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonVectorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMesonVectorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}


void TensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void TensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigTensorMesonVectorVectorDecayer("Herwig::TensorMesonVectorVectorDecayer", "HwTMDecay.so");

void TensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorVectorDecayer> documentation
    ("The TensorMesonVectorVectorDecayer class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static Command<TensorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(GeV) and max weight for a decay",
     &TensorMesonVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<TensorMesonVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<TensorMesonVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in TensorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

}
void TensorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,
					  decay[ix]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double TensorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  // check for photons
  bool photon[2] = {outgoing_[imode()].first ==ParticleID::gamma,
		    outgoing_[imode()].second==ParticleID::gamma};
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  Energy2 denom[2];
  for(unsigned int iy=0;iy<2;++iy) {
    denom[iy] = momenta[iy]*part.momentum()-momenta[iy].mass()*part.mass();
    vectors_[iy].resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      if(photon[iy] && ix==1) continue;
      vectors_[iy][ix] = HelicityFunctions::polarizationVector(-momenta[iy],ix,Helicity::outgoing);
    }
  }
  double fact(coupling_[imode()]/part.mass());
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    LorentzPolarizationVectorE v1 = tensors_[inhel].preDot (momenta[0]);
    LorentzPolarizationVectorE v2 = tensors_[inhel].postDot(momenta[1]);
    complex<Energy2> d0 = v1*momenta[1];
    for(unsigned int vhel1=0;vhel1<3;++vhel1) {
      LorentzPolarizationVector v3 = tensors_[inhel].postDot(vectors_[0][vhel1]);
      complex<InvEnergy> d1 = vectors_[0][vhel1]*part.momentum()/denom[0];
      complex<Energy> d4 = v2*vectors_[0][vhel1];
      for(unsigned int vhel2=0;vhel2<3;++vhel2) {
	complex<InvEnergy> d2 = vectors_[1][vhel2]*part.momentum()/denom[1];
	if ( (photon[0] && vhel1==1) || (photon[1] && vhel2==1))
	  (*ME())(inhel,vhel1,vhel2)=0.;
	else {
	  (*ME())(inhel,vhel1,vhel2)=fact*(v3*vectors_[1][vhel2] -(v1*vectors_[1][vhel2])*d1
					   -d4*d2+d0*d1*d2);
	}
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // double test = sqr(coupling_[imode()]/part.mass());
  // if(photon[0]&&photon[1]) test*=7./15.;
  // else if(photon[0]||photon[1]) test*=2./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  if(outgoing_[imode()].first == outgoing_[imode()].second) me *=0.5;
  return me;
}

bool TensorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						   double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
	imode=ix;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
	imode=ix;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
	imode=ix;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
	imode=ix;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]/dm.parent()->mass();
  mecode=0;
  return id1==outgoing_[imode].first&&id2==outgoing_[imode].second;
}

void TensorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]/GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g*GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./TensorMesonTensorScalarDecayer.cc"
// -*- C++ -*-
//
// TensorMesonTensorScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonTensorScalarDecayer class.
//

#include "TensorMesonTensorScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonTensorScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonTensorScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonTensorScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMesonTensorScalarDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void TensorMesonTensorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void TensorMesonTensorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonTensorScalarDecayer,DecayIntegrator>
describeHerwigTensorMesonTensorScalarDecayer("Herwig::TensorMesonTensorScalarDecayer", "HwTMDecay.so");

void TensorMesonTensorScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonTensorScalarDecayer> documentation
    ("The TensorMesonTensorScalarDecayer class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static Command<TensorMesonTensorScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, tensor, scalar, coupling(GeV) and max weight for a decay",
     &TensorMesonTensorScalarDecayer::setUpDecayMode, false);

}
void TensorMesonTensorScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_in_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  TensorWaveFunction::constructSpinInfo(tensors_out_,decay[0],outgoing,true,false);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double TensorMesonTensorScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin2,PDT::Spin0)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_in_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_out_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel);
    tensors_out_[ihel] = twave.wave();
  }
  double fact(coupling_[imode()]/part.mass());
  Energy2 denom = momenta[0]*part.momentum()-momenta[0].mass()*part.mass();
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    LorentzPolarizationVectorE v0 = tensors_in_[inhel].postDot(momenta[0]);
    Complex d0=(v0*momenta[0])/denom;
    for(unsigned int thel=0;thel<5;++thel) {
      LorentzPolarizationVectorE v1 = tensors_out_[thel].postDot(part.momentum());
      Complex d1=(v1*part.momentum())/denom;
      (*ME())(inhel,thel,0) = fact*(tensors_in_[inhel]*tensors_out_[thel]
				    -2.*(v0*v1)/denom +d0*d1);
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // double test =sqr( coupling_[imode()]/part.mass());
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool TensorMesonTensorScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						   double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]/dm.parent()->mass();
  mecode=6;
  return order;
}

void TensorMesonTensorScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]/GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMesonTensorScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g*GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./TensorMesonTensorPScalarDecayer.cc"
// -*- C++ -*-
//
// TensorMesonTensorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonTensorPScalarDecayer class.
//

#include "TensorMesonTensorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonTensorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonTensorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonTensorPScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMesonTensorPScalarDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void TensorMesonTensorPScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1./GeV);
}

void TensorMesonTensorPScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonTensorPScalarDecayer,DecayIntegrator>
describeHerwigTensorMesonTensorPScalarDecayer("Herwig::TensorMesonTensorPScalarDecayer", "HwTMDecay.so");

void TensorMesonTensorPScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonTensorPScalarDecayer> documentation
    ("The TensorMesonTensorPScalarDecayer class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static Command<TensorMesonTensorPScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, tensor, scalar, coupling(1/GeV) and max weight for a decay",
     &TensorMesonTensorPScalarDecayer::setUpDecayMode, false);

}
void TensorMesonTensorPScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_in_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  TensorWaveFunction::constructSpinInfo(tensors_out_,decay[0],outgoing,true,false);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double TensorMesonTensorPScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin2,PDT::Spin0)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_in_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_out_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel);
    tensors_out_[ihel] = twave.wave();
  }
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  Energy2 denom = momenta[0]*part.momentum()-momenta[0].mass()*part.mass();
  LorentzTensor<Energy2> t0 = epsilon(momenta[0],momenta[1]);
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    LorentzTensor<Energy2> t1 = tensors_in_[inhel].outerProduct(t0);
    LorentzVector<complex<Energy3> > v0 = t1.preDot(momenta[0]);
    for(unsigned int thel=0;thel<5;++thel) {
      LorentzPolarizationVectorE v1 = tensors_out_[thel].postDot(part.momentum());
      (*ME())(inhel,thel,0) = Complex(fact*(t1*tensors_out_[thel] -(v0*v1)/denom ));
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),momenta[1].mass());
  // double test = 0.5*sqr( coupling_[imode()]/part.mass()*pcm*part.mass());
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool TensorMesonTensorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						   double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=17;
  return order;
}

void TensorMesonTensorPScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMesonTensorPScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./TensorMesonVectorScalarDecayer.cc"
// -*- C++ -*-
//
// TensorMesonVectorScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonVectorScalarDecayer class.
//

#include "TensorMesonVectorScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonVectorScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonVectorScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonVectorScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMesonVectorScalarDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void TensorMesonVectorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_;
}

void TensorMesonVectorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonVectorScalarDecayer,DecayIntegrator>
describeHerwigTensorMesonVectorScalarDecayer("Herwig::TensorMesonVectorScalarDecayer", "HwTMDecay.so");

void TensorMesonVectorScalarDecayer::Init() {

  static ClassDocumentation<TensorMesonVectorScalarDecayer> documentation
    ("The TensorMesonVectorScalarDecayer class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static Command<TensorMesonVectorScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, tensor, scalar, coupling(GeV) and max weight for a decay",
     &TensorMesonVectorScalarDecayer::setUpDecayMode, false);

}
void TensorMesonVectorScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  VectorWaveFunction::constructSpinInfo(vectors_,decay[0],outgoing,true,
					decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double TensorMesonVectorScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin0)));
  // check for photons
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  InvEnergy fact(coupling_[imode()]/part.mass());
  Energy2 denom = momenta[0]*part.momentum()-momenta[0].mass()*part.mass();
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    LorentzPolarizationVectorE v0 = tensors_[inhel].postDot(momenta[0]-momenta[1]);
    Complex d0=(v0*momenta[0])/denom;
    for(unsigned int vhel=0;vhel<3;++vhel) {
      if(vhel==1&&photon) (*ME())(inhel,vhel,0)=0.;
      else {
	(*ME())(inhel,vhel,0) = Complex(fact*(v0*vectors_[vhel]-vectors_[vhel]*part.momentum()*d0));
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = sqr(coupling_[imode()]/part.mass())*sqr(pcm);
  // test *= photon ? 4./5. : 4./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool TensorMesonVectorScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						   double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode];
  mecode=0;
  return order;
}

void TensorMesonVectorScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix] << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string TensorMesonVectorScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./TensorMesonSpin3VectorDecayer.cc"
// -*- C++ -*-
//
// Spin3MesonTensorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorMesonSpin3VectorDecayer class.
//

#include "TensorMesonSpin3VectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void TensorMesonSpin3VectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void TensorMesonSpin3VectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters TensorMesonSpin3VectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int TensorMesonSpin3VectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void TensorMesonSpin3VectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void TensorMesonSpin3VectorDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TensorMesonSpin3VectorDecayer,DecayIntegrator>
describeHerwigTensorMesonSpin3VectorDecayer("Herwig::TensorMesonSpin3VectorDecayer", "HwTMDecay.so");

void TensorMesonSpin3VectorDecayer::Init() {

  static ClassDocumentation<TensorMesonSpin3VectorDecayer> documentation
    ("The TensorMesonSpin3VectorDecayer class is designed for the decay"
     " of a tensor meson to a tensor and pseudoscalar mesons.");

    static Command<TensorMesonSpin3VectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(GeV) and max weight for a decay",
     &TensorMesonSpin3VectorDecayer::setUpDecayMode, false);

}

string TensorMesonSpin3VectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 3";
  out.first = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy g = stof(stype)*GeV;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void TensorMesonSpin3VectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  Rank3TensorWaveFunction::constructSpinInfo(rank3_,decay[0],
					     outgoing,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_,decay[1],outgoing,true,
					decay[1]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double TensorMesonSpin3VectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin3,PDT::Spin1)));
  // check for photons
  bool photon(outgoing[1]->id()==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
  }
  rank3_.resize(7);
  Rank3TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  for(unsigned int ihel=0;ihel<7;++ihel) {
    twave.reset(ihel);
    rank3_[ihel] = twave.wave();
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Energy2 denom[2] = {part.momentum()*momenta[0]-part.mass()      *momenta[0].mass(),
   		      momenta[0]     *momenta[1]-momenta[0].mass()*momenta[1].mass()};
  double fact(coupling_[imode()]/part.mass());
  for(unsigned int ih1=0;ih1<7;++ih1) {
    for(unsigned int ih2=0;ih2<3;++ih2) {
      if(ih2==1 && photon) {
  	for(unsigned int ih0=0;ih0<5;++ih0)  (*ME())(ih0,ih1,ih2)=0.;
      }
      else {
	LorentzPolarizationVector v0  = vectors_[ih2] - momenta[1]*momenta[0].dot(vectors_[ih2])/denom[1];
	LorentzTensor<double>     t0  = rank3_[ih1].dot(v0,0);
	LorentzPolarizationVectorE v1 = t0.postDot(momenta[1]);
	Complex d1 = v1*momenta[1]/denom[0];
	for(unsigned int ih0=0;ih0<5;++ih0) {
	  LorentzPolarizationVectorE v2 = tensors_[ih0].postDot(momenta[0]);
	  Complex d2 = v2*momenta[0]/denom[0];
	  (*ME())(ih0,ih1,ih2) = fact*(t0*tensors_[ih0] - 2.*v1.dot(v2)/denom[0]+d1*d2);
	}
      }
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // double test = 7./5.*sqr(fact);
  // if(photon) test *=2./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool TensorMesonSpin3VectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						  double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]/dm.parent()->mass();
  mecode=22;
  return order;
}

void TensorMesonSpin3VectorDecayer::dataBaseOutput(ofstream & output,
						   bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]/GeV << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./PseudoTensorMesonTensorVectorDecayer.cc"
// -*- C++ -*-
//
// PseudoTensorMesonTensorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoTensorMesonTensorVectorDecayer class.
//

#include "PseudoTensorMesonTensorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoTensorMesonTensorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void PseudoTensorMesonTensorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters PseudoTensorMesonTensorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PseudoTensorMesonTensorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void PseudoTensorMesonTensorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_;
}

void PseudoTensorMesonTensorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoTensorMesonTensorVectorDecayer,DecayIntegrator>
describeHerwigPseudoTensorMesonTensorVectorDecayer("Herwig::PseudoTensorMesonTensorVectorDecayer", "HwTMDecay.so");

void PseudoTensorMesonTensorVectorDecayer::Init() {

  static ClassDocumentation<PseudoTensorMesonTensorVectorDecayer> documentation
    ("The PseudoTensorMesonTensorVectorDecayer class implements the"
     " decay of a tensor meson to a spin-1 particle and a pseduoscalar meson");

  static Command<PseudoTensorMesonTensorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, tensor, scalar, coupling(GeV) and max weight for a decay",
     &PseudoTensorMesonTensorVectorDecayer::setUpDecayMode, false);

}
void PseudoTensorMesonTensorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_in_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  TensorWaveFunction::constructSpinInfo(tensors_out_,decay[0],outgoing,true,false);
  // set up the spin information for the decay products
  VectorWaveFunction::constructSpinInfo(vectors_,decay[1],outgoing,true,
					decay[1]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double PseudoTensorMesonTensorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin2,PDT::Spin1)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_in_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false,tensor_phase);
  }
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_out_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel,tensor_phase);
    tensors_out_[ihel] = twave.wave();
  }
  bool photon(outgoing[1]->id()==ParticleID::gamma);
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,
							 Helicity::outgoing,vector_phase);
  }
  InvEnergy fact(coupling_[imode()]/part.mass());
  Energy2 denom = momenta[0]*part.momentum()-momenta[0].mass()*part.mass();
  // calculate the matrix element
  for(unsigned int ih2=0;ih2<3;++ih2) {
    if ( photon && ih2==1 ) {
      for(unsigned int ih0=0;ih0<5;++ih0)
	for(unsigned int ih1=0;ih1<5;++ih1) (*ME())(ih0,ih1,ih2)=0.;
    }
    else {
      LorentzTensor<Energy> eps = epsilon(momenta[1],vectors_[ih2]);
      for(unsigned int ih0=0;ih0<5;++ih0) {
	LorentzTensor<Energy> t0 = tensors_in_[ih0].outerProduct(eps);
	LorentzVector<complex<Energy2> > v0 = t0.preDot(momenta[0]);
	for(unsigned int ih1=0;ih1<5;++ih1) {
	  LorentzPolarizationVectorE v1 = tensors_out_[ih1].postDot(momenta[1]);
	  (*ME())(ih0,ih1,ih2) = Complex(fact*(t0*tensors_out_[ih1] - v0.dot(v1)/denom));
	}
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // const Energy & m0 = part.mass();
  // const Energy & m1 = momenta[0].mass();
  // const Energy & m2 = momenta[1].mass();
  // Energy2 m02(sqr(m0)),m12(sqr(m1)),m22(sqr(m2));
  // double test = sqr(fact)*(sqr(m02-m12)*(11*m02- 8*m0*m1 + 11*m12)
  // 			   - 2*(11*sqr(m02) - 52*m02*m12 + 11*sqr(m12))*m22 + 
  // 			   (11*m02 - 8*m0*m1 + 11*m12)*sqr(m22))/(120.*m02*m12);
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool PseudoTensorMesonTensorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
							 double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode];
  mecode=23;
  return order;
}

void PseudoTensorMesonTensorVectorDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix] << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PseudoTensorMesonTensorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./PseudoTensorMesonVectorVectorDecayer.cc"
// -*- C++ -*-
//
// PseudoTensorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoTensorMesonVectorVectorDecayer class.
//

#include "PseudoTensorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoTensorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()){
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void PseudoTensorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters PseudoTensorMesonVectorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PseudoTensorMesonVectorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].second&&id2   ==outgoing_[ix].first)||
	 (id2   ==outgoing_[ix].second&&id1   ==outgoing_[ix].first)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first)||
	 (id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void PseudoTensorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1./GeV);
}

void PseudoTensorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoTensorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigPseudoTensorMesonVectorVectorDecayer("Herwig::PseudoTensorMesonVectorVectorDecayer", "HwTMDecay.so");

void PseudoTensorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<PseudoTensorMesonVectorVectorDecayer> documentation
    ("The PseudoTensorMesonVectorVectorDecayer class implements the"
     " decay of a pseudotensor meson to two vector mesons");

  static Command<PseudoTensorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &PseudoTensorMesonVectorVectorDecayer::setUpDecayMode, false);

}
void PseudoTensorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,
					  decay[ix]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double PseudoTensorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  // check for photons
  bool photon[2] = {outgoing_[imode()].first ==ParticleID::gamma,
		    outgoing_[imode()].second==ParticleID::gamma};
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin2);
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  Energy2 denom[2];
  for(unsigned int iy=0;iy<2;++iy) {
    denom[iy] = momenta[iy]*part.momentum()-momenta[iy].mass()*part.mass();
    vectors_[iy].resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      if(photon[iy] && ix==1) continue;
      vectors_[iy][ix] = HelicityFunctions::polarizationVector(-momenta[iy],ix,Helicity::outgoing);
    }
  }
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),momenta[1].mass());
  InvEnergy2 c1 = 0.5*(sqr(part.mass())+sqr(momenta[0].mass())-sqr(momenta[1].mass()))/sqr(part.mass()*pcm);
  InvEnergy2 c2 = 0.5*(sqr(part.mass())-sqr(momenta[0].mass())+sqr(momenta[1].mass()))/sqr(part.mass()*pcm);
  // calculate the matrix element
  for(unsigned int inhel=0;inhel<5;++inhel) {
    LorentzPolarizationVectorE v0 = tensors_[inhel].postDot(momenta[0]-momenta[1]);
    LorentzVector<complex<Energy3> > v1 = epsilon(v0,momenta[0],momenta[1]);
    for(unsigned int vhel1=0;vhel1<3;++vhel1) {
      LorentzVector<complex<Energy2> > v3 = epsilon(part.momentum(),v0,vectors_[0][vhel1]);
      complex<InvEnergy> d1 = c1*vectors_[0][vhel1]*momenta[1];
      complex<Energy3> d3 =  vectors_[0][vhel1]*v1;
      for(unsigned int vhel2=0;vhel2<3;++vhel2) {
	if ( (photon[0] && vhel1==1) || (photon[1] && vhel2==1))
	  (*ME())(inhel,vhel1,vhel2)=0.;
	else {
	  complex<InvEnergy> d2 = c2*vectors_[1][vhel2]*momenta[0];
	  (*ME())(inhel,vhel1,vhel2)=Complex(fact*(v3.dot(vectors_[1][vhel2])
						   - d1*v1.dot(vectors_[1][vhel2]) - d2*d3));
	}
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // double test = 16./15.*sqr(coupling_[imode()]*pcm);
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool PseudoTensorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
							 double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); 
  bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].second&&id2==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2==outgoing_[ix].second&&id1==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].second&&id2bar==outgoing_[ix].first) {
	imode=ix;
	order=true;
      }
      if(id2bar==outgoing_[ix].second&&id1bar==outgoing_[ix].first) {
	imode=ix;
	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=9;
  return order;
}

void PseudoTensorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PseudoTensorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./Spin3Meson2PScalarDecayer.cc"
// -*- C++ -*-
//
// Spin3Meson2PScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Spin3Meson2PScalarDecayer class.
//

#include "Spin3Meson2PScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Spin3Meson2PScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void Spin3Meson2PScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters Spin3Meson2PScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int Spin3Meson2PScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void Spin3Meson2PScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV2);
}

void Spin3Meson2PScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Spin3Meson2PScalarDecayer,DecayIntegrator>
describeHerwigSpin3Meson2PScalarDecayer("Herwig::Spin3Meson2PScalarDecayer", "HwTMDecay.so");

void Spin3Meson2PScalarDecayer::Init() {

  static ClassDocumentation<Spin3Meson2PScalarDecayer> documentation
    ("The Spin3Meson2PScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

    static Command<Spin3Meson2PScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(1/GeV^2) and max weight for a decay",
     &Spin3Meson2PScalarDecayer::setUpDecayMode, false);

}

string Spin3Meson2PScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 3";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 0";
  out.first = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy2 g = stof(stype)/GeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void Spin3Meson2PScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  Rank3TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					     incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<decay.size();++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

// matrix elememt for the process
double Spin3Meson2PScalarDecayer::me2(const int,const Particle & part,
				      const tPDVector & ,
				      const vector<Lorentz5Momentum> & momenta,
				      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin3,PDT::Spin0,PDT::Spin0)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin3);
    Rank3TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  // calculate the matrix element
  Lorentz5Momentum pdiff=momenta[0]-momenta[1];
  for(unsigned int ix=0;ix<7;++ix) {
    (*ME())(ix,0,0) = Complex(coupling_[imode()]/part.mass()*
			      ((tensors_[ix].dot(pdiff,0)*pdiff)*pdiff));
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = Energy6(pow<6,1>(pcm))*sqr( coupling_[imode()]/part.mass())*128./35.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // cout << "testing masses " << part.mass()/GeV << " " << momenta[0].mass()/GeV << " " << momenta[1].mass()/GeV << "\n";
  // return the answer
  return output;
}

bool Spin3Meson2PScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					       double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*sqr(dm.parent()->mass());
  mecode=13;
  return order;
}

void Spin3Meson2PScalarDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*GeV2 << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./Spin3MesonVectorScalarDecayer.cc"
// -*- C++ -*-
//
// Spin3MesonVectorScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Spin3MesonVectorScalarDecayer class.
//

#include "Spin3MesonVectorScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Spin3MesonVectorScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void Spin3MesonVectorScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters Spin3MesonVectorScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int Spin3MesonVectorScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void Spin3MesonVectorScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV);
}

void Spin3MesonVectorScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Spin3MesonVectorScalarDecayer,DecayIntegrator>
describeHerwigSpin3MesonVectorScalarDecayer("Herwig::Spin3MesonVectorScalarDecayer", "HwTMDecay.so");

void Spin3MesonVectorScalarDecayer::Init() {

  static ClassDocumentation<Spin3MesonVectorScalarDecayer> documentation
    ("The Spin3MesonVectorScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

    static Command<Spin3MesonVectorScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(1/GeV) and max weight for a decay",
     &Spin3MesonVectorScalarDecayer::setUpDecayMode, false);

}

string Spin3MesonVectorScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 3";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  out.first = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy g = stof(stype)/GeV;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void Spin3MesonVectorScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  Rank3TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					     incoming,true,false);
  // set up the spin information for the decay products
  VectorWaveFunction::constructSpinInfo(vectors_,decay[0],outgoing,true,
					decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double Spin3MesonVectorScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin3,PDT::Spin1,PDT::Spin0)));
  // check for photons
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin3);
    Rank3TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Lorentz5Momentum pdiff=momenta[0]-momenta[1];
  Energy2 denom = part.momentum()*momenta[0]-part.mass()*momenta[0].mass();
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  for(unsigned int ih0=0;ih0<7;++ih0) {
    LorentzVector<complex<Energy2> > v0 = tensors_[ih0].dot(pdiff,0)*pdiff;
    complex<Energy3> d0 = v0*momenta[0];
    for(unsigned int ih1=0;ih1<3;++ih1) {
      if(photon && ih1==1) continue;
      (*ME())(ih0,ih1,0) = Complex(fact*(v0*vectors_[ih1] - (vectors_[ih1]*part.momentum())*d0/denom));
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = 32./15.*pow<4,1>(pcm)*sqr(fact);
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool Spin3MesonVectorScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						    double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*dm.parent()->mass();
  mecode=14;
  return order;
}

void Spin3MesonVectorScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*GeV << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./Spin3MesonVectorPScalarDecayer.cc"
// -*- C++ -*-
//
// Spin3MesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Spin3MesonVectorPScalarDecayer class.
//

#include "Spin3MesonVectorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Spin3MesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void Spin3MesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters Spin3MesonVectorPScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int Spin3MesonVectorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void Spin3MesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV/GeV2);
}

void Spin3MesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Spin3MesonVectorPScalarDecayer,DecayIntegrator>
describeHerwigSpin3MesonVectorPScalarDecayer("Herwig::Spin3MesonVectorPScalarDecayer", "HwTMDecay.so");

void Spin3MesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<Spin3MesonVectorPScalarDecayer> documentation
    ("The Spin3MesonVectorPScalarDecayer class is designed for the decay"
     " of a tensor meson to two (pseudo)-scalar mesons.");

    static Command<Spin3MesonVectorPScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(1/GeV^3) and max weight for a decay",
     &Spin3MesonVectorPScalarDecayer::setUpDecayMode, false);

}

string Spin3MesonVectorPScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 3";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1";
  out.first = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy3 g = stof(stype)/GeV/GeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void Spin3MesonVectorPScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  Rank3TensorWaveFunction::constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
					     incoming,true,false);
  // set up the spin information for the decay products
  VectorWaveFunction::constructSpinInfo(vectors_,decay[0],outgoing,true,
					decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double Spin3MesonVectorPScalarDecayer::me2(const int,const Particle & part,
				      const tPDVector & ,
				      const vector<Lorentz5Momentum> & momenta,
				      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin3,PDT::Spin1,PDT::Spin0)));
  // check for photons
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin3);
    Rank3TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Lorentz5Momentum pdiff=momenta[0]-momenta[1];
  InvEnergy4 fact(coupling_[imode()]/part.mass());
  for(unsigned int ih0=0;ih0<7;++ih0) {
    LorentzVector<complex<Energy4> > v0 = epsilon(part.momentum(),momenta[0],
						  tensors_[ih0].dot(pdiff,0)*pdiff);
    for(unsigned int ih1=0;ih1<3;++ih1) {
      if(photon && ih1==1) continue;
      (*ME())(ih0,ih1,0) = Complex(fact*(v0*vectors_[ih1]));
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = 128./105.*pow<6,1>(pcm)*sqr(fact*part.mass());
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool Spin3MesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					       double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*pow<3,1>(dm.parent()->mass());
  mecode=15;
  return order;
}

void Spin3MesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*GeV*GeV2 << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./Spin3MesonTensorPScalarDecayer.cc"
// -*- C++ -*-
//
// Spin3MesonTensorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Spin3MesonTensorPScalarDecayer class.
//

#include "Spin3MesonTensorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Spin3MesonTensorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void Spin3MesonTensorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters Spin3MesonTensorPScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int Spin3MesonTensorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void Spin3MesonTensorPScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV2);
}

void Spin3MesonTensorPScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Spin3MesonTensorPScalarDecayer,DecayIntegrator>
describeHerwigSpin3MesonTensorPScalarDecayer("Herwig::Spin3MesonTensorPScalarDecayer", "HwTMDecay.so");

void Spin3MesonTensorPScalarDecayer::Init() {

  static ClassDocumentation<Spin3MesonTensorPScalarDecayer> documentation
    ("The Spin3MesonTensorPScalarDecayer class is designed for the decay"
     " of a tensor meson to a tensor and pseudoscalar mesons.");

    static Command<Spin3MesonTensorPScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(1/GeV2) and max weight for a decay",
     &Spin3MesonTensorPScalarDecayer::setUpDecayMode, false);

}

string Spin3MesonTensorPScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 3";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 2";
  out.first = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy2 g = stof(stype)/GeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void Spin3MesonTensorPScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  Rank3TensorWaveFunction::constructSpinInfo(rank3_,const_ptr_cast<tPPtr>(&part),
					     incoming,true,false);
  // set up the spin information for the decay products
  TensorWaveFunction::constructSpinInfo(tensors_,decay[0],
					outgoing,true,false);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

// matrix elememt for the process
double Spin3MesonTensorPScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin3,PDT::Spin2,PDT::Spin0)));
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin3);
    Rank3TensorWaveFunction::
      calculateWaveFunctions(rank3_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  tensors_.resize(5);
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel,tensor_phase);
    tensors_[ihel] = twave.wave();
  }
  // calculate the matrix element
  Lorentz5Momentum pdiff=momenta[0]-momenta[1];
  Energy2 denom = part.momentum()*momenta[0]-part.mass()*momenta[0].mass();
  InvEnergy3 fact(coupling_[imode()]/part.mass());
  LorentzTensor<Energy2> eps = epsilon(part.momentum(),momenta[0]);
  for(unsigned int ih0=0;ih0<7;++ih0) {
    LorentzTensor<Energy> t0 = rank3_[ih0].dot(pdiff,0);
    LorentzTensor<Energy3> t1 = eps.outerProduct(t0);
    LorentzVector<complex<Energy4> > v0 = t1.postDot(momenta[0]);
    for(unsigned int ih1=0;ih1<5;++ih1) {
      LorentzVector<complex<Energy> > v1 = tensors_[ih1].postDot(part.momentum());
      (*ME())(ih0,ih1,0) = Complex(fact*(t1*tensors_[ih1] -v0*v1/denom ));
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  //  					     momenta[1].mass());
  // double test = 8./15.*pow<4,1>(pcm)*sqr(fact)*sqr(part.mass());
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool Spin3MesonTensorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						    double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]*sqr(dm.parent()->mass());
  mecode=16;
  return order;
}

void Spin3MesonTensorPScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*GeV2 << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./Spin3MesonTensorVectorDecayer.cc"
// -*- C++ -*-
//
// Spin3MesonTensorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Spin3MesonTensorVectorDecayer class.
//

#include "Spin3MesonTensorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/Rank3TensorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Spin3MesonTensorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void Spin3MesonTensorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters Spin3MesonTensorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt(0);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
     tPDPtr     in = getParticleData(incoming_[ix]);
     tPDVector out = {getParticleData(outgoing_[ix].first),
		      getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int Spin3MesonTensorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
						const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id(parent->id());
  int idbar = parent->CC() ? parent->CC()->id() : id;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  int imode(-1);
  unsigned int ix(0);
  cc=false;
  do {
    if(id   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(idbar==incoming_[ix]) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void Spin3MesonTensorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void Spin3MesonTensorVectorDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Spin3MesonTensorVectorDecayer,DecayIntegrator>
describeHerwigSpin3MesonTensorVectorDecayer("Herwig::Spin3MesonTensorVectorDecayer", "HwTMDecay.so");

void Spin3MesonTensorVectorDecayer::Init() {

  static ClassDocumentation<Spin3MesonTensorVectorDecayer> documentation
    ("The Spin3MesonTensorVectorDecayer class is designed for the decay"
     " of a tensor meson to a tensor and pseudoscalar mesons.");

    static Command<Spin3MesonTensorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling(GeV) and max weight for a decay",
     &Spin3MesonTensorVectorDecayer::setUpDecayMode, false);

}

string Spin3MesonTensorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 3";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin2)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 2";
  out.first = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy g = stof(stype)*GeV;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void Spin3MesonTensorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  Rank3TensorWaveFunction::constructSpinInfo(rank3_,const_ptr_cast<tPPtr>(&part),
					     incoming,true,false);
  // set up the spin information for the decay products
  TensorWaveFunction::constructSpinInfo(tensors_,decay[0],
					outgoing,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_,decay[1],outgoing,true,
					decay[1]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double Spin3MesonTensorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin3,PDT::Spin2,PDT::Spin1)));
  // check for photons
  bool photon(outgoing[1]->id()==ParticleID::gamma);
  // stuff for incoming particle
  if(meopt==Initialize) {
    rho_ = RhoDMatrix(PDT::Spin3);
    Rank3TensorWaveFunction::
      calculateWaveFunctions(rank3_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
  }
  tensors_.resize(5);
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel,tensor_phase);
    tensors_[ihel] = twave.wave();
  }
  vectors_.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Energy2 denom[2] = {part.momentum()*momenta[0]-part.mass()*momenta[0].mass(),
		      part.momentum()*momenta[1]-part.mass()*momenta[1].mass()};
  double fact(coupling_[imode()]/part.mass());
  for(unsigned int ih0=0;ih0<7;++ih0) {
    for(unsigned int ih2=0;ih2<3;++ih2) {
      if(ih2==1 && photon) {
	for(unsigned int ih1=0;ih1<5;++ih1)  (*ME())(ih0,ih1,ih2)=0.;
      }
      else {
	LorentzPolarizationVector v0  = vectors_[ih2] - momenta[1]*momenta[0].dot(vectors_[ih2])/denom[1];
	LorentzTensor<double>     t0  = rank3_[ih0].dot(v0,0);
	LorentzPolarizationVectorE v1 = t0.postDot(momenta[0]);
	Complex d1 = v1*momenta[0]/denom[0];
	for(unsigned int ih1=0;ih1<5;++ih1) {
	  LorentzPolarizationVectorE v2 = tensors_[ih1].postDot(momenta[1]);
	  Complex d2 = v2*momenta[1]/denom[0];
	  (*ME())(ih0,ih1,ih2) = fact*(t0*tensors_[ih1] - 2.*v1.dot(v2)/denom[0]+d1*d2);
	}
      }
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the answer
  // double test = sqr(fact);
  // if(photon) test *=2./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << test << " " << (output-test)/(output+test) << endl;
  // return the answer
  return output;
}

bool Spin3MesonTensorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						  double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  int idbar = dm.parent()->CC() ? dm.parent()->CC()->id() : id;
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());
  int id1bar = (**pit).CC() ? (**pit).CC()->id() : id1;
  ++pit;
  int id2((**pit).id());
  int id2bar = (**pit).CC() ? (**pit).CC()->id() : id2;
  unsigned int ix(0); bool order(false);
  do {
    if(id   ==incoming_[ix]) {
      if(id1==outgoing_[ix].first&&id2==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2==outgoing_[ix].first&&id1==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    if(idbar==incoming_[ix]&&imode<0) {
      if(id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second) {
  	imode=ix;
  	order=true;
      }
      if(id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second) {
  	imode=ix;
  	order=false;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  coupling=coupling_[imode]/dm.parent()->mass();
  mecode=21;
  return order;
}

void Spin3MesonTensorVectorDecayer::dataBaseOutput(ofstream & output,
						   bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]/GeV << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
