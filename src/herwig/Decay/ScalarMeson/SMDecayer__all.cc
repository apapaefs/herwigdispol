#line 1 "./EtaPiGammaGammaDecayer.cc"
// -*- C++ -*-
//
// EtaPiGammaGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiGammaGammaDecayer class.
//
#include "EtaPiGammaGammaDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void EtaPiGammaGammaDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _etamax  = mode(0)->maxWeight();
    _etapmax = mode(1)->maxWeight();
  }
}

EtaPiGammaGammaDecayer::EtaPiGammaGammaDecayer()
  : _grhoomega(12.924/GeV), _fpi(130.7*MeV),_rhomass(771.1*MeV),
    _rhowidth(149.2*MeV),_grho(_rhomass/_fpi),_mpi(ZERO),_rhoconst(0.),
    _localparameters(true),_ratiofpif8(1./1.3),_ratiofpif0(1./1.04),
    _theta(-Constants::pi/9.),_etamax(2.36858),_etapmax(0.006),
    _dconst(2), _econst(2) {
  // intermediates
  generateIntermediates(false);
}

void EtaPiGammaGammaDecayer::doinit() {
  DecayIntegrator::doinit();
  // set rho parameters if needed
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters) {
    _rhomass  = rho->mass();
    _rhowidth = rho->width();
  }
  // constant for the running rho width
  _mpi=getParticleData(ParticleID::pi0)->mass();
  Energy pcm =Kinematics::pstarTwoBodyDecay(_rhomass,_mpi,_mpi);
  _rhoconst=_rhomass*_rhomass*_rhowidth/(pcm*pcm*pcm);
  // set the prefactors
  double conv(sqrt(4.*Constants::pi*SM().alphaEM()));
  conv *=_fpi*_fpi*_grho/_rhomass/_rhomass;
  InvEnergy2 pre(2.*sqrt(3.)/9.*sqr(_grhoomega*conv));
  double fact[2];
  // constants for eta
  fact[0] = _ratiofpif8*cos(_theta)-sqrt(2.)*_ratiofpif0*sin(_theta);
  // constants for eta'
  fact[1] = _ratiofpif8*sin(_theta)+sqrt(2.)*_ratiofpif0*cos(_theta);
  for(unsigned int ix=0;ix<2;++ix) {
    _dconst[ix]=fact[ix]*pre;
    _econst[ix]=fact[ix]*pre;
  }
  // set up the phase space for the decays
  tPDPtr eta[2]={getParticleData(ParticleID::eta),
		 getParticleData(ParticleID::etaprime)};
  tPDVector out = {getParticleData(ParticleID::pi0),
		   getParticleData(ParticleID::gamma),
		   getParticleData(ParticleID::gamma)};
  for(unsigned int ix=0;ix<2;++ix) {
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(eta[ix],out,
				  (ix==0 ? _etamax : _etapmax)));
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rho,0,2,1,1,1,3));
    c1.weight(0.5);
    mode->addChannel(c1);
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rho,0,3,1,1,1,2));
    c2.weight(0.5);
    mode->addChannel(c2);
    addMode(mode);
  }
}

int EtaPiGammaGammaDecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  cc=false;
  int id;
  if(children.size()!=3) return -1;
  tPDVector::const_iterator pit = children.begin();
  unsigned int npi0(0),ngamma(0);
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::pi0)         ++npi0;
    else if(id==ParticleID::gamma)  ++ngamma;
  }
  if(!(npi0==1&&ngamma==2)) return -1;
  // number of the mode
  switch (parent->id()) {
  case ParticleID::eta     : return 0;
  case ParticleID::etaprime: return 1;
  default: return -1;
  }
}

void EtaPiGammaGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_grhoomega,1/GeV)<< ounit(_fpi,GeV)<< _grho 
     << ounit(_rhomass,GeV)<< ounit(_rhowidth,GeV)<< _localparameters 
     << _ratiofpif8 << _ratiofpif0 << _theta << _etamax << _etapmax 
     << _rhoconst << ounit(_mpi,GeV) << ounit(_dconst,1/GeV2) 
     << ounit(_econst,1/GeV2);
}

void EtaPiGammaGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_grhoomega,1/GeV) >> iunit(_fpi,GeV)>> _grho 
     >> iunit(_rhomass,GeV)>> iunit(_rhowidth,GeV)>> _localparameters 
     >> _ratiofpif8 >> _ratiofpif0 >> _theta >> _etamax >> _etapmax 
     >> _rhoconst >> iunit(_mpi,GeV) >> iunit(_dconst,1/GeV2) 
     >> iunit(_econst,1/GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPiGammaGammaDecayer,DecayIntegrator>
describeHerwigEtaPiGammaGammaDecayer("Herwig::EtaPiGammaGammaDecayer", "HwSMDecay.so");

void EtaPiGammaGammaDecayer::Init() {

  static ClassDocumentation<EtaPiGammaGammaDecayer> documentation
    ("The EtaPiGammaGammaDecayer class implements a VMD model for the"
     " decay of the eta or etaprime to a pion and two photons.",
     "The decays of $\\eta,\\eta'\\to\\pi^0\\gamma\\gamma$ were simulated using"
     " the matrix elements of \\cite{Holstein:2001bt}",
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Parameter<EtaPiGammaGammaDecayer,InvEnergy> interfacegrhoomega
    ("grhoomega",
     "The couping of the rho, omega and a pion",
     &EtaPiGammaGammaDecayer::_grhoomega, 1./GeV, 12.924/GeV, ZERO, 100./GeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant",
     &EtaPiGammaGammaDecayer::_fpi, MeV, 130.7*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfacegrho
    ("grho",
     "Rho decay constant",
     &EtaPiGammaGammaDecayer::_grho, 5.9, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho meson",
     &EtaPiGammaGammaDecayer::_rhomass, MeV, 771.1*MeV, 500.0*MeV, 1000.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho meson",
     &EtaPiGammaGammaDecayer::_rhowidth, MeV, 149.2*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceRatioFpiF8
    ("RatioFpiF8",
     "The ratio of the decay constant Fpi to F8",
     &EtaPiGammaGammaDecayer::_ratiofpif8, 1./1.3, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceRatioFpiF0
    ("RatioFpiF0",
     "The ratio of the decay constant Fpi to F0",
     &EtaPiGammaGammaDecayer::_ratiofpif0, 1./1.04, 0.0, 10.0,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceTheta
    ("Theta",
     "The eta etaprime mixing angle",
     &EtaPiGammaGammaDecayer::_theta, -Constants::pi/9., -Constants::pi, Constants::pi,
     false, false, true);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaMax
    ("EtaMax",
     "THe maximum weight for the eta decay",
     &EtaPiGammaGammaDecayer::_etamax, 1.35, -1.0e12, 1.0e12,
     false, false, false);

  static Parameter<EtaPiGammaGammaDecayer,double> interfaceEtaPrimeMax
    ("EtaPrimeMax",
     "THe maximum weight for the eta prime decay",
     &EtaPiGammaGammaDecayer::_etapmax, 0.006, -1.0e12, 1.0e12,
     false, false, false);

  static Switch<EtaPiGammaGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the parameters",
     &EtaPiGammaGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);
}

void EtaPiGammaGammaDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix+1],
					  outgoing,true,true);
}

double EtaPiGammaGammaDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    _vectors[ix].resize(3);
    for(unsigned int ihel=0;ihel<3;ihel+=2) {
      _vectors[ix][ihel] = HelicityFunctions::polarizationVector(-momenta[ix+1],ihel,Helicity::outgoing);
    }
  }
  // dot products we need
  Energy2 q1dotq2(momenta[1]*momenta[2]),
    pdotq1(part.momentum()*momenta[1]),
    pdotq2(part.momentum()*momenta[2]);
  complex<Energy> e1dotq2[3],e1dotp[3],e2dotq1[3],e2dotp[3];
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) {
      e1dotq2[ix]=ZERO;
      e1dotp[ix] =ZERO;
      e2dotq1[ix]=ZERO;
      e2dotp[ix] =ZERO;
    }
    else {
      e1dotq2[ix] =_vectors[0][ix]*momenta[2];
      e1dotp[ix]  =_vectors[0][ix]*part.momentum();
      e2dotq1[ix] =_vectors[1][ix]*momenta[1];
      e2dotp[ix]  =_vectors[1][ix]*part.momentum();
    }
  }
  // the momentum dependent pieces of the matrix element
  Complex ii(0.,1.);
  Energy2 mpi2(sqr(momenta[0].mass())),meta2(sqr(part.mass())),
    mrho2(sqr(_rhomass)),
    t(mpi2+2.*((momenta[0])*(momenta[1]))),
    u(mpi2+2.*((momenta[0])*(momenta[2])));
  Energy q(sqrt(t)),pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
  complex<Energy2> tgamma(ii*pcm*pcm*pcm*_rhoconst/q);
  q=sqrt(u);pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy2> ugamma(ii*pcm*pcm*pcm*_rhoconst/q);
  complex<InvEnergy2> prop1(1./(mrho2-t-tgamma)),prop2(1./(mrho2-u-ugamma));
  complex<InvEnergy2> Dfact(_dconst[imode()]*(prop1*(pdotq2-meta2)
					      +prop2*(pdotq1-meta2)));
  complex<InvEnergy4> Efact(_econst[imode()]*(prop1+prop2));
  Complex e1dote2;
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      if(ix==1||iy==1) (*ME())(0,0,ix,iy)=0.;
      else {
	e1dote2=_vectors[0][ix].dot(_vectors[1][iy]);
	(*ME())(0,0,ix,iy) = 
	  Complex(Dfact*complex<Energy2>(e1dote2*q1dotq2-
					 e1dotq2[ix]*e2dotq1[iy])
		  -Efact*complex<Energy4>(-e1dote2*pdotq1*pdotq2
					  -e1dotp[ix]*e2dotp[iy]*q1dotq2
					  +e1dotq2[ix]*e2dotp[iy]*pdotq1
					  +e1dotp[ix]*e2dotq1[iy]*pdotq2));
      }
    }
  }
  double me(ME()->contract(_rho).real());
  // test of the me
  // Energy M(part.mass());
  // Energy2 M2(M*M);
  // Energy2 s1(2.*(momenta[1]*momenta[2]));
  // Energy2 s2(M2-2.*(part.momentum()*momenta[1]));
  // Energy2 s3(M2-2.*(part.momentum()*momenta[2]));
  // cout << "testing the matrix element " << (
  //  2*(2*(Dfact*conj(Dfact)).real() + 2*(Dfact*conj(Efact)).real()*M2 
  //     + (Efact*conj(Efact)).real()*M2*M2)*
  //     s1*s1 - 2*(Efact*conj(Efact)).real()*M2*s1*(M2 - s2)*
  //  (M2 - s3) +(Efact*conj(Efact)).real()*(M2 - s2)*(M2 - s2)*
  //  (M2-s3)*(M2-s3))/8. - me << endl;
  return me;
}
 
double EtaPiGammaGammaDecayer::
threeBodyMatrixElement(const int imodeb, const Energy2 q2,const  Energy2 s3,
		       const Energy2 s2,const Energy2 s1,const Energy ,
		       const Energy ,const Energy ) const {
  // compute the prefactors
  Energy2 mrho2 = sqr(_rhomass);
  Energy q = sqrt(s3);
  Energy pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  Complex ii(0.,1.);
  complex<Energy2> tgamma(ii*pcm*pcm*pcm*_rhoconst/q);
  q = sqrt(s2);
  pcm = Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi);
  complex<Energy2> ugamma(ii*pcm*pcm*pcm*_rhoconst/q);
  complex<InvEnergy2> prop1(1./(mrho2-s3-tgamma)), prop2(1./(mrho2-s2-ugamma));
  Energy2 pdotq2(0.5*(q2-s3)), pdotq1(0.5*(q2-s2));
  complex<InvEnergy2> Dfact(_dconst[imodeb]*(prop1*(pdotq2-q2)+prop2*(pdotq1-q2)));
  complex<InvEnergy4> Efact(_econst[imodeb]*(prop1+prop2));
  InvEnergy4 D2 = (Dfact*conj(Dfact)).real();
  InvEnergy8 E2((Efact*conj(Efact)).real());
  InvEnergy6 ED((Efact*conj(Dfact)).real());
  return (2 * (2*D2 + 2*ED*q2 + E2*sqr(q2)) * sqr(s1)
	  - double(2*E2*q2*s1*(q2-s2)*(q2-s3))
	  + double(E2*sqr(q2-s2)*sqr(q2-s3))
	  )/8.;
}

WidthCalculatorBasePtr 
EtaPiGammaGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights; inweights.push_back(0.5); inweights.push_back(0.5);
  Energy mrho(getParticleData(ParticleID::rhoplus)->mass());
  Energy wrho(getParticleData(ParticleID::rhoplus)->width());
  vector<Energy> inmass;  inmass.push_back(mrho);  inmass.push_back(mrho);
  vector<Energy> inwidth; inwidth.push_back(wrho); inwidth.push_back(wrho);
  vector<int> intype; intype.push_back(1); intype.push_back(2);
  vector<double> inpow(2,0.0);
  return new_ptr(ThreeBodyAllOnCalculator<EtaPiGammaGammaDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,
		  imode,_mpi,ZERO,ZERO));
}

void EtaPiGammaGammaDecayer::dataBaseOutput(ofstream & output, 
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":grhoomega " << _grhoomega*GeV << "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV  << "\n";
  output << "newdef " << name() << ":grho " << _grho << "\n";
  output << "newdef " << name() << ":RhoMass " << _rhomass/MeV << "\n";
  output << "newdef " << name() << ":RhoWidth " << _rhowidth/MeV << "\n";
  output << "newdef " << name() << ":RatioFpiF8 " << _ratiofpif8 << "\n";
  output << "newdef " << name() << ":RatioFpiF0 " << _ratiofpif0 << "\n";
  output << "newdef " << name() << ":Theta " << _theta  << "\n";
  output << "newdef " << name() << ":EtaMax " << _etamax << "\n";
  output << "newdef " << name() << ":EtaPrimeMax " << _etapmax << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./EtaPiPiGammaDecayer.cc"
// -*- C++ -*-
//
// EtaPiPiGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiGammaDecayer class.
//

#include "EtaPiPiGammaDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "Herwig/Decay/FormFactors/OmnesFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


DescribeClass<EtaPiPiGammaDecayer,DecayIntegrator>
describeHerwigEtaPiPiGammaDecayer("Herwig::EtaPiPiGammaDecayer",
				  "HwSMDecay.so");
HERWIG_INTERPOLATOR_CLASSDESC(EtaPiPiGammaDecayer,double,Energy)



void EtaPiPiGammaDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<_maxweight.size();++ix)
      _maxweight[ix]=mode(ix)->maxWeight();
  }
}

EtaPiPiGammaDecayer::EtaPiPiGammaDecayer() 
  : _incoming(2), _coupling(2), _maxweight(2), _option(2) {
  // the pion decay constant
  _fpi=130.7*MeV;
  // the rho mass
  _mrho=0.7711*GeV;
  _rhowidth=0.1492*GeV;
  // the constants for the omnes function form
  _aconst=0.5/_mrho/_mrho;
  _cconst=1.0;
  // use local values of the parameters
  _localparameters=true;
  // the modes
  // eta decay
  _incoming[0] = 221; 
  _option[0] = 3; 
  _coupling[0] = 5.060e-3; 
  _maxweight[0] = 3.95072; 
  // eta' decay
  _incoming[1] = 331; 
  _option[1] = 3; 
  _coupling[1] = 4.278e-3; 
  _maxweight[1] = 3.53141; 
  _rhoconst=0.;
  _mpi=ZERO;
  // intermediates
  generateIntermediates(false);
}

void EtaPiPiGammaDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the consistence of the parameters
  unsigned int isize=_incoming.size();
  if(isize!=_coupling.size()||isize!=_option.size()||isize!=_maxweight.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "EtaPiPiGammaDecayer::doinit()" << Exception::abortnow;
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!_localparameters) {
    _mrho=rho->mass();
    _rhowidth=rho->width();
  }
  _mpi=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(_mrho,_mpi,_mpi));
  _rhoconst=sqr(_mrho)*_rhowidth/pow<3,1>(pcm);
  // set up the modes
  tPDVector out = {getParticleData(ParticleID::piplus),
		   getParticleData(ParticleID::piminus),
		   getParticleData(ParticleID::gamma)};
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<_coupling.size();++ix) {
    tPDPtr in = getParticleData(_incoming[ix]);
    mode = new_ptr(PhaseSpaceMode(in,out,_maxweight[ix]));
    mode->addChannel((PhaseSpaceChannel(mode),0,rho,0,3,1,1,1,2));
    addMode(mode);
  }
}

int EtaPiPiGammaDecayer::modeNumber(bool & cc,tcPDPtr parent,
				    const tPDVector & children) const {
  int imode(-1);
  // check number of external particles
  if(children.size()!=3){return imode;}
  // check the outgoing particles
  unsigned int npip(0),npim(0),ngamma(0);
  tPDVector::const_iterator pit = children.begin();
  int id;
  for(;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::gamma)   ++ngamma;
  }
  if(!(npip==1&&npim==1&&ngamma==1)) return imode;
  unsigned int ix(0);
  id=parent->id();
  do{if(id==_incoming[ix]){imode=ix;}++ix;}
  while(imode<0&&ix<_incoming.size());
  cc=false;
  return imode;
}

void EtaPiPiGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_fpi,MeV) << _incoming << _coupling << _maxweight << _option 
     << ounit(_aconst,1/MeV2) << _cconst <<ounit(_mrho,MeV) << ounit(_rhowidth,MeV) 
     << _rhoconst << ounit(_mpi,MeV) << _localparameters << omnesFunction_;
}

void EtaPiPiGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_fpi,MeV) >> _incoming >> _coupling >> _maxweight >> _option 
     >> iunit(_aconst,1/MeV2) >> _cconst >>iunit(_mrho,MeV) >> iunit(_rhowidth,MeV) 
     >> _rhoconst >> iunit(_mpi,MeV) >> _localparameters >> omnesFunction_;
}

void EtaPiPiGammaDecayer::Init() {

  static ClassDocumentation<EtaPiPiGammaDecayer> documentation
    ("The EtaPiPiGammaDecayer class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma",
     "The decays of $\\eta,\\eta'\\to\\pi^+\\pi^-\\gamma$ were simulated"
     " using the matrix elements from \\cite{Venugopal:1998fq,Holstein:2001bt}",
     "\\bibitem{Venugopal:1998fq} E.~P.~Venugopal and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 57} (1998) 4397 [arXiv:hep-ph/9710382].\n"
     "%%CITATION = PHRVA,D57,4397;%%\n"
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Reference<EtaPiPiGammaDecayer,OmnesFunction> interfaceOmnesFunction
    ("OmnesFunction",
     "Omnes function for the matrix element",
     &EtaPiPiGammaDecayer::omnesFunction_, false, false, true, false, false);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiGammaDecayer::_fpi, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static ParVector<EtaPiPiGammaDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code for the incoming particle",
     &EtaPiPiGammaDecayer::_incoming,
     0, 0, 0, -10000000, 10000000, false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decay mode",
     &EtaPiPiGammaDecayer::_coupling,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<EtaPiPiGammaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &EtaPiPiGammaDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &EtaPiPiGammaDecayer::_mrho, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &EtaPiPiGammaDecayer::_rhowidth, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<EtaPiPiGammaDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &EtaPiPiGammaDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Parameter<EtaPiPiGammaDecayer,double> interfaceOmnesC
    ("OmnesC",
     "The constant c for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_cconst, 1.0, -10., 10.,
     false, false, false);

  static Parameter<EtaPiPiGammaDecayer,InvEnergy2> interfaceOmnesA
    ("OmnesA",
     "The constant a for the Omnes form of the prefactor",
     &EtaPiPiGammaDecayer::_aconst, 1./GeV2, 0.8409082/GeV2, ZERO,
     10./GeV2,
     false, false, false);

  static ParVector<EtaPiPiGammaDecayer,int> interfaceOption
    ("Option",
     "The form of the prefactor 0 is a VMD model using M Gamma for the width term,"
     "1 is a VMD model using q Gamma for the width term,"
     "2. analytic form of the Omnes function,"
     "3. experimental form of the Omnes function.",
     &EtaPiPiGammaDecayer::_option,
     0, 0, 0, 0, 4, false, false, true);
}
void EtaPiPiGammaDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  VectorWaveFunction::constructSpinInfo(_vectors,decay[2],
					outgoing,true,true);
}

double EtaPiPiGammaDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
  }
  _vectors.resize(3);
  for(unsigned int ix=0;ix<3;ix+=2) {
    _vectors[ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  }
  // prefactor for the matrix element
  complex<InvEnergy3> pre(_coupling[imode()]*2.*sqrt(2.)/(_fpi*_fpi*_fpi));
  Lorentz5Momentum ppipi(momenta[0]+momenta[1]);
  ppipi.rescaleMass();
  Energy q(ppipi.mass());
  Energy2 q2(q*q);
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(_option[imode()]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(q2/(_mrho*_mrho-q2-ii*_mrho*pcm*pcm*pcm*_rhoconst/q2));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(_option[imode()]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(q2/(_mrho*_mrho-q2-ii*pcm*pcm*pcm*_rhoconst/q));
    fact=(1.+1.5*resfact);
  }
  // omnes function
  else if(_option[imode()]==2 || _option[imode()]==3) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*q2)/omnesFunction_->D(q2));
  }
  pre = pre*fact;
  LorentzPolarizationVector epstemp(pre*Helicity::epsilon(momenta[0],
							  momenta[1],
							  momenta[2]));
  // compute the matrix element
  vector<unsigned int> ispin(4,0);
  for(ispin[3]=0;ispin[3]<3;++ispin[3]) {
    if(ispin[3]==1) (*ME())(ispin)=0.;
    else            (*ME())(ispin)=epstemp.dot(_vectors[ispin[3]]);
  }
  // contract the whole thing
  return ME()->contract(_rho).real();
}

double EtaPiPiGammaDecayer::
threeBodyMatrixElement(const int imodeb,const Energy2 ,const  Energy2 s3,const 
		       Energy2 s2,const Energy2 s1,const Energy ,
		       const Energy ,const Energy ) const {
  complex<InvEnergy3> pre(_coupling[imodeb]*2.*sqrt(2.)/pow<3,1>(_fpi));
  Energy q(sqrt(s3));
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(_option[imodeb]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(s3/(_mrho*_mrho-s3-ii*_mrho*pcm*pcm*pcm*_rhoconst/s3));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(_option[imodeb]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,_mpi,_mpi));
    Complex resfact(s3/(_mrho*_mrho-s3-ii*pcm*pcm*pcm*_rhoconst/q));
    fact=(1.+1.5*resfact);
  }
  // omnes function
  else if(_option[imodeb]==2||_option[imodeb]==3) {
    fact=(1.-_cconst+_cconst*(1.+_aconst*s3)/omnesFunction_->D(s3));
  }
  pre =pre*fact;
  InvEnergy6 factor((pre*conj(pre)).real());
  Energy2 mpi2(_mpi*_mpi);
  return factor*((-mpi2*(-2*mpi2+s1+s2)*(-2*mpi2+s1+s2)+(mpi2-s1)*(mpi2-s2)*s3)/4.);
}

WidthCalculatorBasePtr 
EtaPiPiGammaDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int id(dm.parent()->id()),imode(1);
  if(id==ParticleID::eta){imode=0;}
  // construct the integrator
  vector<double> inweights(1,1.);
  vector<Energy> inmass(1,getParticleData(ParticleID::rho0)->mass());
  vector<Energy> inwidth(1,getParticleData(ParticleID::rho0)->width());
  vector<int> intype(1,1);
  vector<double> inpow(1,0.0);
  WidthCalculatorBasePtr 
    output(new_ptr(ThreeBodyAllOnCalculator<EtaPiPiGammaDecayer>
		   (inweights,intype,inmass,inwidth,inpow,*this,imode,_mpi,_mpi,ZERO)));
  return output;
}

void EtaPiPiGammaDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":fpi             " << _fpi/MeV         << "\n";
  output << "newdef " << name() << ":RhoMass         " << _mrho/MeV        << "\n";
  output << "newdef " << name() << ":RhoWidth        " << _rhowidth/MeV    << "\n";
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":OmnesC          " << _cconst          << "\n";
  output << "newdef " << name() << ":OmnesA          " << _aconst*GeV2     << "\n";
  for(unsigned int ix=0;ix<2;++ix) {
    output << "newdef " << name() << ":Incoming    " << ix << "  " 
	   << _incoming[ix]    << "\n";
    output << "newdef " << name() << ":Coupling    " << ix << "  " 
	   << _coupling[ix]    << "\n";
    output << "newdef " << name() << ":MaxWeight   " << ix << "  " 
	   << _maxweight[ix]   << "\n";
    output << "newdef " << name() << ":Option      " << ix << "  " 
	   << _option[ix]      << "\n";
  }
  omnesFunction_->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":OmnesFunction " << omnesFunction_->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./EtaPiPiFermionsDecayer.cc"
// -*- C++ -*-
//
// EtaPiPiFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiFermionsDecayer class.
//

#include "EtaPiPiFermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "Herwig/Decay/FormFactors/OmnesFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


DescribeClass<EtaPiPiFermionsDecayer,DecayIntegrator>
describeHerwigEtaPiPiFermionsDecayer("Herwig::EtaPiPiFermionsDecayer",
				  "HwSMDecay.so");
HERWIG_INTERPOLATOR_CLASSDESC(EtaPiPiFermionsDecayer,double,Energy)



void EtaPiPiFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
      maxWeight_[ix]=mode(ix)->maxWeight();
    }
  }
}

EtaPiPiFermionsDecayer::EtaPiPiFermionsDecayer() 
  : fPi_(130.7*MeV), mRho_(0.7711*GeV), rhoWidth_(0.1492*GeV) {
  // the constants for the omnes function form
  aConst_=0.5/mRho_/mRho_;
  cConst_=1.0;
  // use local values of the parameters
  localParameters_=true;
  // the modes
  mPi_=ZERO;
  // intermediates
  generateIntermediates(false);
}

void EtaPiPiFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=coupling_.size()||isize!=option_.size()||isize!=lepton_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "EtaPiPiFermionsDecayer::doinit()" << Exception::abortnow;
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!localParameters_) {
    mRho_=rho->mass();
    rhoWidth_=rho->width();
  }
  mPi_=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(mRho_,mPi_,mPi_));
  rhoConst_=sqr(mRho_)*rhoWidth_/pow<3,1>(pcm);
  // set up the modes
  tPDPtr gamma = getParticleData(ParticleID::gamma);
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<coupling_.size();++ix) {
    tPDVector out = {getParticleData(ParticleID::piplus),
		     getParticleData(ParticleID::piminus),
		     getParticleData( lepton_[ix]),
		     getParticleData(-lepton_[ix])};
    tPDPtr in = getParticleData(incoming_[ix]);
    mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,rho,0,gamma,1,1,1,2,2,3,2,4));
    newChannel.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel);
    addMode(mode);
  }
}

int EtaPiPiFermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
				    const tPDVector & children) const {
  int imode(-1);
  // check number of external particles
  if(children.size()!=4){return imode;}
  // check the outgoing particles
  unsigned int npip(0),npim(0),nl(0);
  int il(0);
  tPDVector::const_iterator pit = children.begin();
  int id;
  for(;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else {
      il = abs(id);
      ++nl;
    }
  }
  if(!(npip==1&&npim==1&&nl==2)) return imode;
  unsigned int ix(0);
  id=parent->id();
  do {
    if(id==incoming_[ix] && il==lepton_[ix]) imode=ix;
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void EtaPiPiFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(fPi_,MeV) << incoming_ << coupling_ << maxWeight_ << lepton_ << option_ 
     << ounit(aConst_,1/MeV2) << cConst_ <<ounit(mRho_,MeV) << ounit(rhoWidth_,MeV) 
     << rhoConst_ << ounit(mPi_,MeV) << localParameters_ << omnesFunction_;
}

void EtaPiPiFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fPi_,MeV) >> incoming_ >> coupling_ >> maxWeight_ >> lepton_ >> option_ 
     >> iunit(aConst_,1/MeV2) >> cConst_ >>iunit(mRho_,MeV) >> iunit(rhoWidth_,MeV) 
     >> rhoConst_ >> iunit(mPi_,MeV) >> localParameters_ >> omnesFunction_;
}

void EtaPiPiFermionsDecayer::Init() {

  static ClassDocumentation<EtaPiPiFermionsDecayer> documentation
    ("The EtaPiPiFermionsDecayer class is design for the decay of"
     " the eta and eta prime to pi+pi-gamma",
     "The decays of $\\eta,\\eta'\\to\\pi^+\\pi^-\\gamma$ were simulated"
     " using the matrix elements from \\cite{Venugopal:1998fq,Holstein:2001bt}",
     "\\bibitem{Venugopal:1998fq} E.~P.~Venugopal and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 57} (1998) 4397 [arXiv:hep-ph/9710382].\n"
     "%%CITATION = PHRVA,D57,4397;%%\n"
     "\\bibitem{Holstein:2001bt} B.~R.~Holstein,\n"
     " Phys.\\ Scripta {\\bf T99} (2002) 55 [arXiv:hep-ph/0112150].\n"
     "%%CITATION = PHSTB,T99,55;%%\n");

  static Reference<EtaPiPiFermionsDecayer,OmnesFunction> interfaceOmnesFunction
    ("OmnesFunction",
     "Omnes function for the matrix element",
     &EtaPiPiFermionsDecayer::omnesFunction_, false, false, true, false, false);

  static Parameter<EtaPiPiFermionsDecayer,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &EtaPiPiFermionsDecayer::fPi_, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static Command<EtaPiPiFermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, lepton, option for VMD, coupling and maxweight."
     "There are three options for for VMD, 0 is a VMD model using M Gamma for the width term,"
     " 1 is a VMD model using q Gamma for the width term,"
     "2. analytic form of the Omnes function,"
     "3. experimental form of the Omnes function.",
     &EtaPiPiFermionsDecayer::setUpDecayMode, false);

  static Parameter<EtaPiPiFermionsDecayer,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &EtaPiPiFermionsDecayer::mRho_, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<EtaPiPiFermionsDecayer,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &EtaPiPiFermionsDecayer::rhoWidth_, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<EtaPiPiFermionsDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &EtaPiPiFermionsDecayer::localParameters_, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static Parameter<EtaPiPiFermionsDecayer,double> interfaceOmnesC
    ("OmnesC",
     "The constant c for the Omnes form of the prefactor",
     &EtaPiPiFermionsDecayer::cConst_, 1.0, -10., 10.,
     false, false, false);

  static Parameter<EtaPiPiFermionsDecayer,InvEnergy2> interfaceOmnesA
    ("OmnesA",
     "The constant a for the Omnes form of the prefactor",
     &EtaPiPiFermionsDecayer::aConst_, 1./GeV2, 0.8409082/GeV2, ZERO,
     10./GeV2,
     false, false, false);

}

void EtaPiPiFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[2],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[3],outgoing,true);
}

double EtaPiPiFermionsDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,
					 PDT::Spin1Half,PDT::Spin1Half)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the spinors
  wavebar_.resize(2);
  wave_   .resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wavebar_[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[2],ihel,Helicity::outgoing);
    wave_   [ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[3],ihel,Helicity::outgoing);
  }
  Lorentz5Momentum pff(momenta[2]+momenta[3]);
  pff.rescaleMass();
  Energy2 mff2(pff.mass()*pff.mass());
  // prefactor for the matrix element
  complex<InvEnergy4> pre(part.mass()*coupling_[imode()]*2.*sqrt(2.)/(fPi_*fPi_*fPi_)*
			  sqrt(SM().alphaEM(mff2)*4.*Constants::pi)/mff2);
  Lorentz5Momentum ppipi(momenta[0]+momenta[1]);
  ppipi.rescaleMass();
  Energy q(ppipi.mass());
  Energy2 q2(q*q);
  Complex ii(0.,1.);
  // first VMD option
  Complex fact;
  if(option_[imode()]==0) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,mPi_,mPi_));
    Complex resfact(q2/(mRho_*mRho_-q2-ii*mRho_*pcm*pcm*pcm*rhoConst_/q2));
    fact=(1.+1.5*resfact);
  }
  // second VMD option
  else if(option_[imode()]==1) {
    Energy pcm(Kinematics::pstarTwoBodyDecay(q,mPi_,mPi_));
    Complex resfact(q2/(mRho_*mRho_-q2-ii*pcm*pcm*pcm*rhoConst_/q));
    fact=(1.+1.5*resfact);
  }
  // omnes function
  else if(option_[imode()]==2 || option_[imode()]==3) {
    fact=(1.-cConst_+cConst_*(1.+aConst_*q2)/omnesFunction_->D(q2));
  }
  pre = pre*fact;
  LorentzVector<complex<InvEnergy> >
    epstemp(pre*Helicity::epsilon(momenta[0],momenta[1],pff));
  // compute the matrix element
  vector<unsigned int> ispin(5,0);
  for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
    for(ispin[4]=0;ispin[4]<2;++ispin[4]) {
      LorentzPolarizationVectorE fcurrent = wave_[ispin[4]].vectorCurrent(wavebar_[ispin[3]]);
      (*ME())(ispin) = Complex(epstemp.dot(fcurrent));
    }
  }
  // contract the whole thing
  return ME()->contract(rho_).real();
}

void EtaPiPiFermionsDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":fpi             " << fPi_/MeV         << "\n";
  output << "newdef " << name() << ":RhoMass         " << mRho_/MeV        << "\n";
  output << "newdef " << name() << ":RhoWidth        " << rhoWidth_/MeV    << "\n";
  output << "newdef " << name() << ":LocalParameters " << localParameters_ << "\n";
  output << "newdef " << name() << ":OmnesC          " << cConst_          << "\n";
  output << "newdef " << name() << ":OmnesA          " << aConst_*GeV2     << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << lepton_[ix] << " " << option_[ix] << " "  
	   << coupling_[ix] << " " << maxWeight_[ix] << "\n";
  }
  omnesFunction_->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":OmnesFunction " << omnesFunction_->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string EtaPiPiFermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "Outgoing fermion with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Outgoing fermion with id " + std::to_string(out) + "does not have spin 1/2";
  // vmd option  
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDinclude = stoi(stype);
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
  lepton_.push_back(out);
  coupling_.push_back(g);
  option_.push_back(VMDinclude);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./EtaPiPiPiDecayer.cc"
// -*- C++ -*-
//
// EtaPiPiPiDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EtaPiPiPiDecayer class.
//
#include "EtaPiPiPiDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/PDT/OneOffShellCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void EtaPiPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void EtaPiPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()||isize!=prefactor_.size()||
     isize!=charged_.size()||isize!=a_.size()||
     isize!=b_.size()||isize!=c_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in EtaPiPiPiDecayer::doinit()"
  			  << Exception::runerror;
  // external particles for the modes
  tPDVector outneut(3),outcharged(3);
  outneut[0]    = getParticleData(ParticleID::pi0);
  outneut[1]    = getParticleData(ParticleID::pi0);
  outcharged[0] = getParticleData(ParticleID::piplus);
  outcharged[1] = getParticleData(ParticleID::piminus);
  tPDPtr rho(getParticleData(113));
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr incoming = getParticleData(incoming_[ix]);
    outneut[2]    = getParticleData(outgoing_[ix]);
    outcharged[2] = getParticleData(outgoing_[ix]);
    // the pi+pi- mode
    if(charged_[ix]) {
      mode = new_ptr(PhaseSpaceMode(incoming,outcharged,maxWeight_[ix]));
    }
    // the pi0pi0 mode
    else {
      mode = new_ptr(PhaseSpaceMode(incoming,outneut,maxWeight_[ix]));
    }
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,rho,0,3,1,1,1,2));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.0);
    mode->addChannel(newChannel);
    addMode(mode);
  }
  resetIntermediate(rho,600.*MeV,600.*MeV);
}

int EtaPiPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				 const tPDVector & children) const {
  if(children.size()!=3) return -1;
  unsigned int npi0(0),npip(0),npim(0); int id,iother(0);
  tPDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)           ++npip;
    else if(id==ParticleID::piminus)     ++npim;
    else if(id==ParticleID::pi0&&npi0<2) ++npi0;
    else iother=id;
  }
  bool charged;
  if(npim==1&&npip==1) {
    charged=true;
    if(npi0==1) iother=ParticleID::pi0;
  }
  else if(npi0==2) charged=false;
  else return -1;
  // find the mode
  id=parent->id();
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id==incoming_[ix]&&iother==outgoing_[ix]&&charged_[ix]==charged) 
      imode=ix;
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void EtaPiPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << charged_ << prefactor_ << a_ << b_ << c_  
     << maxWeight_;
}

void EtaPiPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> charged_ >> prefactor_ >> a_ >> b_ >> c_ 
     >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<EtaPiPiPiDecayer,DecayIntegrator>
describeHerwigEtaPiPiPiDecayer("Herwig::EtaPiPiPiDecayer", "HwSMDecay.so");

void EtaPiPiPiDecayer::Init() {

  static ClassDocumentation<EtaPiPiPiDecayer> documentation
    ("The EtaPiPiPiDecayer class performs the decay of a scalar meson to"
     " two pions and another meson using a simple paramterisation of the dalitz plot.",
     "The decay of eta to two pions follows \\cite{Beisert:2003zs,Gormley:1970qz,Tippens:2001fm}.",
     "%\\cite{Beisert:2003zs}\n"
     "\\bibitem{Beisert:2003zs}\n"
     "  N.~Beisert and B.~Borasoy,\n"
     "  %``Hadronic decays of eta and eta' with coupled channels,''\n"
     "  Nucl.\\ Phys.\\  A {\\bf 716}, 186 (2003)\n"
     "  [arXiv:hep-ph/0301058].\n"
     "  %%CITATION = NUPHA,A716,186;%%\n"
     "%\\cite{Gormley:1970qz}\n"
     "\\bibitem{Gormley:1970qz}\n"
     "  M.~Gormley, E.~Hyman, W.~Y.~Lee, T.~Nash, J.~Peoples, C.~Schultz and S.~Stein,\n"
     "   ``Experimental determination of the dalitz-plot distribution of the decays\n"
     "   eta $\\to$ pi+ pi- pi0 and eta $\\to$ pi+ pi- gamma, and the branching ratio\n"
     "  %eta $\\to$ pi+ pi- gamma/eta $\\to$ pi+,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 2}, 501 (1970).\n"
     "  %%CITATION = PHRVA,D2,501;%%\n"
     "%\\cite{Tippens:2001fm}\n"
     "\\bibitem{Tippens:2001fm}\n"
     "  W.~B.~Tippens {\\it et al.}  [Crystal Ball Collaboration],\n"
     "  %``Determination of the quadratic slope parameter in eta $\\to$ 3pi0 decay,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 87}, 192001 (2001).\n"
     "  %%CITATION = PRLTA,87,192001;%%\n"
     );

  static Command<EtaPiPiPiDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the decay mode (incoming, outgoing, charged/neutral pions, prefactor, a, b, c parameters and maximum weight",
     &EtaPiPiPiDecayer::setUpDecayMode, false);

  static Deleted<EtaPiPiPiDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<EtaPiPiPiDecayer> interfaceOutgoing
    ("Outgoing","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfaceCharged
    ("Charged","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfacePrefactor
    ("Prefactor","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfacea
    ("a","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfaceb
    ("b","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfacec
    ("c","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<EtaPiPiPiDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in EtaPiPiPiDecayer have been deleted, please use SetUpDecayMode");

}

void EtaPiPiPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}



double EtaPiPiPiDecayer::me2(const int,const Particle & part,
					const tPDVector &,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the matrix element
  // compute the variables we need
  Lorentz5Momentum ps(part.momentum()-momenta[2]);
  ps.rescaleMass();
  Lorentz5Momentum pu(part.momentum()-momenta[0]);
  pu.rescaleMass();
  Lorentz5Momentum pt(part.momentum()-momenta[1]);
  pt.rescaleMass();
  Energy2 s(ps.mass2()),u(pu.mass2()),t(pt.mass2());
  Energy m34(0.5*(momenta[0].mass()+momenta[1].mass()));
  Energy msum(momenta[2].mass()+2.*m34);
  Energy Q(part.mass()-msum);
  Energy2 Mmm2((part.mass()-momenta[2].mass())*(part.mass()-momenta[2].mass()));
  // compute the variables
  double x(0.5*sqrt(3.)*(u-t)/part.mass()/Q),x2(x*x);
  double y(0.5*msum/part.mass()*(Mmm2-s)/m34/Q-1),y2(y*y);
  double me(prefactor_[imode()]*(1+a_[imode()]*y+b_[imode()]*y2+c_[imode()]*x2));
  if(me<0.) me=0.;
  (*ME())(0,0,0,0)=sqrt(me);
  return me;
}

InvEnergy EtaPiPiPiDecayer::threeBodydGammads(const int imodeb, const Energy2 q2,
					   const  Energy2 s, const Energy m1,
					   const Energy m2, const Energy m3) const {
  Energy q(sqrt(q2)),m34(m1+m2),msum(m34+m3),Q(q-msum);
  Energy2 Mmm2((q-m3)*(q-m3)),m12(m1*m1),m22(m2*m2),m32(m3*m3);
  double y(0.5*msum/q*(Mmm2-s)/m34/Q-1),y2(y*y);
  InvEnergy2 xfact=0.5*sqrt(3.)/q/Q;
  Energy2 xc(q2+m12+m22+m32-s);
  Energy rs(sqrt(s)),e2star(0.5*(s-m12+m22)/rs),e3star(0.5*(q2-s-m32)/rs);
  Energy e2sm(sqrt(e2star*e2star-m22)),e3sm(sqrt(e3star*e3star-m32));
  Energy2 a(2*e2star*e3star+m22+m32),b(2*e2sm*e3sm);
  Energy2 output=2*b*(1+a_[imodeb]*y+b_[imodeb]*y2+c_[imodeb]*xfact*xfact*(xc*xc))
    +c_[imodeb]*(-8.*xfact*xfact*xc*a*b
		 +4.*2*b*(3.*a*a+b*b)/3.*xfact*xfact);
  using Constants::pi;
  return output*prefactor_[imodeb]/256./pi/pi/pi/q2/q;
}


WidthCalculatorBasePtr 
EtaPiPiPiDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  int idout(0),id,imode(-1);
  unsigned int npi0(0),ix(0);
  ParticleMSet::const_iterator pit(dm.products().begin());
  for( ;pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::pi0&&npi0<2)                            ++npi0;
    else if(id!=ParticleID::piplus&&id!=ParticleID::piminus) idout=id;
  }
  if(npi0==1) idout=ParticleID::pi0;
  bool charged(npi0<2);
  id=dm.parent()->id();
  do {
    if(id==incoming_[ix]&&idout==outgoing_[ix]&&charged_[ix]==charged) 
      imode=ix;
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  Energy mpi;
  if(charged){mpi=getParticleData(ParticleID::piplus)->mass();}
  else{mpi=getParticleData(ParticleID::pi0)->mass();}
  Energy m[3]={mpi,mpi,getParticleData(outgoing_[imode])->mass()};
  WidthCalculatorBasePtr 
    temp(new_ptr(ThreeBodyAllOn1IntegralCalculator<EtaPiPiPiDecayer>
  		 (1,-1000.*MeV,ZERO,0.0,*this,imode,m[0],m[1],m[2])));
  if(outgoing_[imode]==ParticleID::eta) {
    tcGenericMassGeneratorPtr test;
    tGenericMassGeneratorPtr massptr;
    if(getParticleData(outgoing_[imode])->massGenerator()) {
      test=dynamic_ptr_cast<tcGenericMassGeneratorPtr>
  	(getParticleData(outgoing_[imode])->massGenerator());
      massptr=const_ptr_cast<tGenericMassGeneratorPtr>(test);
    }
    if(massptr) {
      massptr->init();
      return new_ptr(OneOffShellCalculator(3,temp,massptr,ZERO));
    }
  }
  return temp;
} 
  
void EtaPiPiPiDecayer::dataBaseOutput(ofstream & output,
				      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix]   << " "
	   << outgoing_[ix]  << " " << charged_[ix]  << " " << prefactor_[ix]  << " "
	   << a_[ix]  << " " << b_[ix]  << " " << c_[ix]  << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string EtaPiPiPiDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 0";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  bool charge = stoi(stype);
  // get the couplings
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double a = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double b = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double c = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_ .push_back(in);
  outgoing_ .push_back(out);
  charged_  .push_back(charge);
  prefactor_.push_back(g);
  a_        .push_back(a);
  b_        .push_back(b);
  c_        .push_back(c);
  maxWeight_.push_back(wgt);
  return "";
}
#line 1 "./PScalar4FermionsDecayer.cc"
// -*- C++ -*-
//
// PScalar4FermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalar4FermionsDecayer class.
//

#include "PScalar4FermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalar4FermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      maxweight_[ix] = mode(ix)->maxWeight();
  }
}

PScalar4FermionsDecayer::PScalar4FermionsDecayer() {
  // intermediates
  generateIntermediates(false);
}

void PScalar4FermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!=maxweight_.size() || isize!=includeVMD_.size()|| isize!=VMDid_.size()    ||
     isize!=VMDmass_.size()  || isize!=VMDwidth_.size())
    throw InitException() << "Inconsistent parameters in PScalar4FermionsDecayer"
  			  << Exception::abortnow;
  // create the integration channels for each mode
  tPDPtr gamma=getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in = getParticleData(incoming_[ix]);
    tPDVector out={getParticleData( outgoing_[ix].first),
		   getParticleData(-outgoing_[ix].first),
		   getParticleData( outgoing_[ix].second),
		   getParticleData(-outgoing_[ix].second)};
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,gamma,0,gamma,1,1,1,2,2,3,2,4));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    newChannel.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel);
    PhaseSpaceChannel newChannel2((PhaseSpaceChannel(mode),0,gamma,0,gamma,1,3,1,2,2,1,2,4));
    newChannel2.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    newChannel2.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel2);
    addMode(mode);
  }
  // set up the values for the VMD factor if needed (copy the default mass and width 
  //                                                 into the array)
  for(unsigned ix=0;ix<isize;++ix) {
    if(includeVMD_[ix]==1) {
      VMDmass_[ix]=getParticleData(VMDid_[ix])->mass();
      VMDwidth_[ix]=getParticleData(VMDid_[ix])->width();
    }
  }
}

int PScalar4FermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  // must be four outgoing particles
  if(children.size()!=4) return -1;
  // get the id's of the outgoing particles
  int id[4]={0,0,0,0}; 
  bool done[4]={false,false,false,false}; 
  unsigned int ix(0),iy(0);
  // ids of the particles
  int id0(parent->id()),idtemp(-1),idl1(-1),idl2(-1),idt[2];
  tPDVector::const_iterator pit = children.begin();
  for ( ;pit!=children.end();++pit) {
    id[ix]=(**pit).id();
    done[ix]=false;
    ++ix;
  }
  // find the two lepton pairs
  // find the first fermion
  ix=0;
  do {
    if( id[ix]>0 && !done[ix] ) {
      done[ix]=true;
      idtemp=id[ix];
    }
    ++ix;
  }
  while(ix<4&&idtemp<0);
  if(idtemp<0) return -1;
  // find its antiparticle
  ix=0;
  do {
    if( id[ix]==-idtemp && !done[ix] ) {
      done[ix]=true;
      idl1=idtemp;
    }
    ++ix;
  } while( ix<4 && idl1<0 );
  if(idl1<0) return -1;
  // find the second particle antiparticle pair
  for(ix=0;ix<4;++ix) {
    if(!done[ix]) {
      idt[iy]=id[ix];
      ++iy;
    }
  }
  if(idt[0]==-idt[1]) idl2=abs(idt[0]);
  if(idl2<0) return -1;
  // loop over the modes and see if this is one of them
  ix=0;
  int imode(-1);
  do {
    if(incoming_[ix]==id0) {
      if((idl1==outgoing_[ix].first&&idl2==outgoing_[ix].second)||
	 (idl2==outgoing_[ix].first&&idl1==outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void PScalar4FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/MeV) 
     << incoming_ << outgoing_ << maxweight_ 
     << includeVMD_ << VMDid_ 
     << ounit(VMDmass_,MeV) << ounit(VMDwidth_,MeV);
}

void PScalar4FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/MeV) 
     >> incoming_ >> outgoing_ >> maxweight_ 
     >> includeVMD_ >> VMDid_ 
     >> iunit(VMDmass_,MeV) >> iunit(VMDwidth_,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalar4FermionsDecayer,DecayIntegrator>
describeHerwigPScalar4FermionsDecayer("Herwig::PScalar4FermionsDecayer", "HwSMDecay.so");

void PScalar4FermionsDecayer::Init() {

  static ClassDocumentation<PScalar4FermionsDecayer> documentation
    ("The PScalar4FermionsDecayer class is designed for the decay"
     " of a pseudosclar meson to four fermions. It is intended for the decay of"
     "the pion to two electron-positron pairs.");

  static Command<PScalar4FermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, 1st fermion, 2nd fermion, coupling(1/GeV), includeVMD,"
     " VMD id, VMD mass(GeV), VMD width(GeV) and max weight for a decay."
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalar4FermionsDecayer::setUpDecayMode, false);

  static Deleted<PScalar4FermionsDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceOutcoming1
    ("Outgoing1","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceOutcoming2
    ("Outgoing2","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceIncludeVMD
    ("IncludeVMD","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceVMDID
    ("VMDID","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceVMDmass
    ("VMDmass","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalar4FermionsDecayer> interfaceVMDwidth
    ("VMDwidth","The old methods of setting up a decay in PScalar4FermionsDecayer have been deleted, please use SetUpDecayMode");

}
void PScalar4FermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<2;++ix) {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_[ix],decay[2*ix  ],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave_[ix]   ,decay[2*ix+1],outgoing,true);
  }
}

double PScalar4FermionsDecayer::me2(const int,const Particle & part,
				    const tPDVector &,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
  					 PDT::Spin1Half,PDT::Spin1Half)));
  bool identical((outgoing_[imode()].first==outgoing_[imode()].second));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the spinors
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix].resize(2);
    wave_   [ix].resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel) {
      wavebar_[ix][ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[2*ix  ],ihel,Helicity::outgoing);
      wave_   [ix][ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[2*ix+1],ihel,Helicity::outgoing);
    }
  }
  // // momenta of the outgoing photons
  Lorentz5Momentum poff[4];
  poff[0]=momenta[0]+momenta[1];
  poff[0].rescaleMass();
  poff[1]=momenta[2]+momenta[3];
  poff[1].rescaleMass();
  if(identical) {
    poff[2]=momenta[2]+momenta[1];
    poff[2].rescaleMass();
    poff[3]=momenta[0]+momenta[3];
    poff[3].rescaleMass();
  }
  // compute the currents for the two leptonic decays
  LorentzPolarizationVectorE current[4][2][2];
  unsigned int it,ix,iy,iz;
  for(iz=0;iz<2;++iz) {
    it = iz==0 ? 1 : 0;
    for(ix=0;ix<2;++ix) {
      for(iy=0;iy<2;++iy) {
  	current[iz][iy][ix] = wave_[iz][ix].vectorCurrent(wavebar_[iz][iy]);
  	// the second two currents      
  	if(identical)
  	  current[iz+2][iy][ix] = wave_[it][ix].vectorCurrent(wavebar_[iz][iy]);
      }
    }
  }
  // invariants
  Energy2 m12(sqr(poff[0].mass()));
  Energy2 m34(sqr(poff[1].mass()));
  Energy2 m14(ZERO), m23(ZERO);
  complex<InvEnergy4> prop1(1./m12/m34),prop2(0./sqr(MeV2));
  Complex ii(0.,1.);
  if(identical) {
    m14=poff[2].mass()*poff[2].mass();
    m23=poff[3].mass()*poff[3].mass();
    prop2=1./m14/m23;
  }
  // the VMD factor if needed
  if(includeVMD_[imode()]>0) {
    Energy2 mrho2(VMDmass_[imode()]*VMDmass_[imode()]);
    Energy2 mwrho(VMDmass_[imode()]*VMDwidth_[imode()]);
    prop1 = prop1*(-mrho2+ii*mwrho)/(m12-mrho2+ii*mwrho)*
                  (-mrho2+ii*mwrho)/(m34-mrho2+ii*mwrho);
    if(identical) {
      prop2 = prop2*(-mrho2+ii*mwrho)/(m14-mrho2+ii*mwrho)*
  	            (-mrho2+ii*mwrho)/(m23-mrho2+ii*mwrho);
    }
  }
  // prefactor
  Complex pre(coupling_[imode()]*4.*Constants::pi
  	      *SM().alphaEM()*part.mass());
  Complex diag;
  // now compute the matrix element
  LorentzVector<complex<Energy3> > eps;
  vector<unsigned int> ispin(5,0);
  for(ispin[1]=0; ispin[1]<2;++ispin[1]) {
    for(ispin[2]=0;ispin[2]<2;++ispin[2]) {
      for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
  	for(ispin[4]=0;ispin[4]<2;++ispin[4]) {
  	  // the first diagram
  	  eps = epsilon(current[0][ispin[1]][ispin[2]],poff[1],
  			current[1][ispin[3]][ispin[4]]);
  	  diag = Complex(prop1*(eps*poff[0]));
  	  // exchanged diagram if identical particles
  	  //  (sign due normal ordering) 
  	  if(identical) {
  	    eps = epsilon(current[2][ispin[1]][ispin[4]],poff[3],
  			  current[3][ispin[3]][ispin[2]]);
  	    diag-= Complex(prop2*(eps*poff[2]));
  	  }
  	  (*ME())(ispin)=pre*diag;
  	}
      }
    }
  }
  double me=ME()->contract(rho_).real();
  if(identical) me *= 0.25;
  return me;
}

// output the setup info for the particle database
void PScalar4FermionsDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << includeVMD_[ix] << " "
	   << VMDid_[ix] << " " << VMDmass_[ix]/GeV << " "
	   << VMDwidth_[ix]/GeV << " " << maxweight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string PScalar4FermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // vmd option  
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDinclude = stoi(stype);
  // vmd id
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDid = stoi(stype);
  // vmd mass
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double VMDmass = stof(stype);
  // vmd width
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double VMDwidth = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  includeVMD_.push_back(VMDinclude);
  VMDid_.push_back(VMDid);
  VMDmass_.push_back(VMDmass*GeV);
  VMDwidth_.push_back(VMDwidth*GeV);
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./PScalarLeptonNeutrinoDecayer.cc"
// -*- C++ -*-
//
// PScalarLeptonNeutrinoDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarLeptonNeutrinoDecayer class.
//

#include "PScalarLeptonNeutrinoDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarLeptonNeutrinoDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  unsigned int iz(0),ix,iy;
  if(initialize()) {
    for(ix=0;ix<incoming_.size();++ix) {
      for(iy=0;iy<leptons_[ix];++iy) {
	maxWeight_[ix][iy] = mode(iz)->maxWeight();
	++iz;
      }
    }
  }
}

void PScalarLeptonNeutrinoDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize(incoming_.size());
  if(isize!=decayConstant_.size()||isize!=leptons_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in PScalarLeptonNeutrinoDecayer"
			  << Exception::abortnow;
  // create the integration channels
  tPDPtr nu[3]={getParticleData(ParticleID::nu_e),
		getParticleData(ParticleID::nu_mu),
		getParticleData(ParticleID::nu_tau)};
  tPDPtr nubar[3]={getParticleData(ParticleID::nu_ebar),
		   getParticleData(ParticleID::nu_mubar),
		   getParticleData(ParticleID::nu_taubar)};
  tPDPtr lep[3]={getParticleData(ParticleID::eminus),
		 getParticleData(ParticleID::muminus),
		 getParticleData(ParticleID::tauminus)};
  tPDPtr lepbar[3]={getParticleData(ParticleID::eplus),
		    getParticleData(ParticleID::muplus),
		    getParticleData(ParticleID::tauplus)};
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in = getParticleData(incoming_[ix]);
    for(unsigned int iy=0;iy<leptons_[ix];++iy) {
      tPDVector out(2);
      if(in->iCharge()>0) {
	out[0]=lepbar[iy];
	out[1]=nu[iy];
      }
      else {
	out[0]=lep[iy];
	out[1]=nubar[iy];
      }
      double wgt = maxWeight_[ix][iy];
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,wgt));
      addMode(mode);
    }
  }
}

int PScalarLeptonNeutrinoDecayer::modeNumber(bool & cc,tcPDPtr parent,
					     const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id()),id0bar(id0);
  if(parent->CC()){id0bar=-id0;}
  tPDVector::const_iterator pit = children.begin();
  int id;
  unsigned ilep(4);
  for(;pit!=children.end();++pit) {
    id=abs((**pit).id());
    if(id>=11&&id<=16&&id%2==0) ilep=(id-10)/2;
  }
  // find the channel we need
  bool found(false);
  int ichan(-1);
  unsigned int ix(0);
  do {
    if(id0   ==incoming_[ix]||id0bar==incoming_[ix]) {
      found=true;ichan+=ilep;
      cc=id0bar==incoming_[ix];
    }
    else {
      ichan+=leptons_[ix];
    }
    ++ix;
  }
  while (!found&&ix<incoming_.size());
  if(found) imode=ichan;
  return imode;
}


void PScalarLeptonNeutrinoDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << ounit(decayConstant_,GeV) << leptons_ << maxWeight_;
}

void PScalarLeptonNeutrinoDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> iunit(decayConstant_,GeV) >> leptons_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarLeptonNeutrinoDecayer,DecayIntegrator>
describeHerwigPScalarLeptonNeutrinoDecayer("Herwig::PScalarLeptonNeutrinoDecayer", "HwSMDecay.so");

void PScalarLeptonNeutrinoDecayer::Init() {

  static ClassDocumentation<PScalarLeptonNeutrinoDecayer> documentation
    ("The PScalarLeptonNeutrinoDecayer class is the base class for"
     " the implementation of leptonic decay of a pseudoscalar meson to a lepton"
     "and a neutrino.");

  static Command<PScalarLeptonNeutrinoDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, decay constant no of leptonic decays, max weights",
     &PScalarLeptonNeutrinoDecayer::setUpDecayMode, false);

  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceLeptons
    ("Leptons","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceMaxWeightElectron
    ("MaxWeightElectron","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceMaxWeightMuon
    ("MaxWeightMuon","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceMaxWeightTau
    ("MaxWeightTau","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<PScalarLeptonNeutrinoDecayer> interfaceDecayConstant
    ("DecayConstant","The old methods of setting up a decay in "
     "PScalarLeptonNeutrinoDecayer have been deleted, please use SetUpDecayMode");
}

void PScalarLeptonNeutrinoDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(0);
  for(unsigned ix=0;ix<decay.size();++ix) {
    long id=decay[ix]->id();
    if(id<=-11&&id>=-16)    ianti=ix;
    else if(id>=11&&id<=16) iferm=ix;
  }
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  // set up the spin information for the decay products
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double PScalarLeptonNeutrinoDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // work out which decay constant to use
  int icoup(0),id(abs(part.id()));
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    if(id==abs(incoming_[ix])) icoup=ix;
  }
  // find the particles
  unsigned int iferm(0),ianti(0);
  for(unsigned ix=0;ix<outgoing.size();++ix) {
    id=outgoing[ix]->id();
    if(id<=-11&&id>=-16)    ianti=ix;
    else if(id>=11&&id<=16) iferm=ix;
  }
  int idferm = outgoing[iferm]->id();
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate the spinors
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // the prefactor
  Energy premass =  idferm%2==0 ? momenta[ianti].mass() : momenta[iferm].mass();
  InvEnergy pre = premass * 2.*decayConstant_[icoup]*SM().fermiConstant()/part.mass();
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  for(ispin[ianti+1]=0;ispin[ianti+1]<2;++ispin[ianti+1]) {
    for(ispin[iferm+1]=0;ispin[iferm+1]<2;++ispin[iferm+1]) {
      (*ME())(ispin)= idferm%2==0 ? 
	Complex(pre*wave_[ispin[ianti+1]].rightScalar(wavebar_[ispin[iferm+1]])) :
	Complex(pre*wave_[ispin[ianti+1]].leftScalar( wavebar_[ispin[iferm+1]])) ;
    }
  }
  double me = 0.5*ME()->contract(rho_).real();
  // test of the matrix element
  // Energy mass = idferm%2==0 ? momenta[ianti].mass() : momenta[iferm].mass();
  // double test = 0.5*sqr(decayConstant_[icoup]*SM().fermiConstant()*2.*mass/part.mass())*
  //   (sqr(part.mass())-sqr(mass));
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << (me-test)/(me+test) << endl;
  return me;
}

void PScalarLeptonNeutrinoDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << decayConstant_[ix]/MeV << " " << leptons_[ix];
    for(unsigned int iy=0;iy<3;++iy) output << " " << maxWeight_[ix][iy];
    output << "\n";
  }
  if(header) 
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string PScalarLeptonNeutrinoDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // no of leptons
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int lepton = stoi(stype);
  // and the maximum weights
  array<double,3> wgt;
  for(unsigned int iy=0;iy<3;++iy) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    wgt[iy] = stof(stype);
  }
  // store the information
  incoming_.push_back(in);
  decayConstant_.push_back(g*MeV);
  leptons_.push_back(lepton);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./PScalarPScalarVectorDecayer.cc"
// -*- C++ -*-
//
// PScalarPScalarVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarPScalarVectorDecayer class.
//

#include "PScalarPScalarVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarPScalarVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void PScalarPScalarVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() || isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in PScalarPScalarVectorDecayer" 
  			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData( incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PScalarPScalarVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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


void PScalarPScalarVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxWeight_;
}

void PScalarPScalarVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarPScalarVectorDecayer,DecayIntegrator>
describeHerwigPScalarPScalarVectorDecayer("Herwig::PScalarPScalarVectorDecayer", "HwSMDecay.so");

void PScalarPScalarVectorDecayer::Init() {

  static ClassDocumentation<PScalarPScalarVectorDecayer> documentation
    ("The PScalarPScalarVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");
  
  static Command<PScalarPScalarVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling (dimensionless) and max weight for a decay",
     &PScalarPScalarVectorDecayer::setUpDecayMode, false);
  
  static Deleted<PScalarPScalarVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceOutcomingP
    ("OutgoingP","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceOutcomingV
    ("OutgoingV","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarPScalarVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalarPScalarVectorDecayer have been deleted, please use SetUpDecayMode");
}

void PScalarPScalarVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  VectorWaveFunction::constructSpinInfo(vectors_,decay[1],
					outgoing,true,false);
}

double PScalarPScalarVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  if(meopt==Initialize) {
     ScalarWaveFunction::
       calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  vectors_.resize(3);
  bool massless = outgoing[1]->id()==ParticleID::gamma;
  for(unsigned int ix=0;ix<3;++ix) {
    if(massless && ix==1) continue;
    vectors_[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // calculate the matrix element
  Lorentz5Momentum psum(part.momentum()+momenta[0]);
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(0,0,ix) = Complex(coupling_[imode()]/part.mass()*(vectors_[ix]*psum));
  }
  double me=ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = 4.*sqr(coupling_[imode()]*pcm/momenta[1].mass());
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return me;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarPScalarVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & mecode,
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
  unsigned int ix(0); bool order(true);
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
  coupling=coupling_[imode];
  mecode=10;
  return order;
}

// output the setup information for the particle database
void PScalarPScalarVectorDecayer::dataBaseOutput(ofstream & output,
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

string PScalarPScalarVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
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
#line 1 "./PScalarVectorFermionsDecayer.cc"
// -*- C++ -*-
//
// PScalarVectorFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorFermionsDecayer class.
//

#include "PScalarVectorFermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarVectorFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      maxweight_[ix]=mode(ix)->maxWeight();
  }
}

PScalarVectorFermionsDecayer::PScalarVectorFermionsDecayer() {
  // intermediates
  generateIntermediates(false);
}

void PScalarVectorFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!=maxweight_.size() || isize!=includeVMD_.size()||
     isize!=VMDid_.size()     || isize!=VMDmass_.size()   || isize!=VMDwidth_.size())
    throw InitException() << "Inconsistent parameters in PScalarVectorFermionsDecayer"
  			  << Exception::abortnow;
  // create the integration channel for each mode 
  tPDPtr gamma(getParticleData(ParticleID::gamma));
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second),
		     getParticleData(-outgoing_[ix].second)};
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    PhaseSpaceChannel newChannel((PhaseSpaceChannel(mode),0,gamma,0,1,1,2,1,3));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    mode->addChannel(newChannel);
    addMode(mode);
  }
  // set up the values for the VMD factor if needed (copy the default mass and width 
  //                                                 into the array)
  for(unsigned ix=0;ix<isize;++ix) {
    if(includeVMD_[ix]==1) {
      VMDmass_[ix]=getParticleData(VMDid_[ix])->mass();
      VMDwidth_[ix]=getParticleData(VMDid_[ix])->width();
    }
  }
}

int PScalarVectorFermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  int imode(-1);
  // must be three outgoing particles
  if(children.size()!=3) return imode;
  // ids of the particles
  int id0(parent->id()),idf[2]={0,0},idv(0);
  unsigned int nf(0);
  tPDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    if((**pit).iSpin()==PDT::Spin1) {
      idv=(**pit).id();
    }
    else {
      idf[nf]=(**pit).id();
      ++nf;
    }
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==idv)
      {if((idf[0]==outgoing_[ix].second&&idf[1]==-outgoing_[ix].second)||
	  (idf[1]==outgoing_[ix].second&&idf[0]==-outgoing_[ix].second)){imode=ix;}}
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  cc=false;
  return imode;
}

void PScalarVectorFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/MeV) 
     << incoming_ << outgoing_ << maxweight_
     << includeVMD_ << VMDid_ 
     << ounit(VMDmass_,MeV) << ounit(VMDwidth_,MeV);
}

void PScalarVectorFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/MeV) 
     >> incoming_ >> outgoing_ >> maxweight_
     >> includeVMD_ >> VMDid_ 
     >> iunit(VMDmass_,MeV) >> iunit(VMDwidth_,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarVectorFermionsDecayer,DecayIntegrator>
describeHerwigPScalarVectorFermionsDecayer("Herwig::PScalarVectorFermionsDecayer", "HwSMDecay.so");

void PScalarVectorFermionsDecayer::Init() {

  static ClassDocumentation<PScalarVectorFermionsDecayer> documentation
    ("The PScalarVectorFermionsDecayer class is designed"
     " for the decay of a pseudoscalar meson to a photon and a"
     "fermion-antifermion pair");

  static Command<PScalarVectorFermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, fermion, coupling(1/GeV), includeVMD,"
     " VMD id, VMD mass(GeV), VMD width(GeV) and max weight for a decay."
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &PScalarVectorFermionsDecayer::setUpDecayMode, false);
  
  static Deleted<PScalarVectorFermionsDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceOutcomingV
    ("OutgoingVector","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceOutcomingF
    ("OutgoingFermion","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceOutcomingA
    ("OutgoingAntiFermion","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceIncludeVMD
    ("IncludeVMD","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceVMDID
    ("VMDID","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceVMDmass
    ("VMDmass","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorFermionsDecayer> interfaceVMDwidth
    ("VMDwidth","The old methods of setting up a decay in PScalarVectorFermionsDecayer have been deleted, please use SetUpDecayMode");

}

void PScalarVectorFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  // set up the spin information for the decay products
  VectorWaveFunction::
    constructSpinInfo(vectors_,decay[0],outgoing,true,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[1],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[2],outgoing,true);
}

double PScalarVectorFermionsDecayer::me2(const int,const Particle & part,
					 const tPDVector &,
					 const vector<Lorentz5Momentum> & momenta,
					 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1Half,
					 PDT::Spin1Half)));
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // calculate polarization vector
  vectors_.resize(3);
  for(unsigned int ihel=0;ihel<3;ihel+=2) {
    vectors_[ihel] = HelicityFunctions::polarizationVector(-momenta[0],ihel,Helicity::outgoing);
  }
  // calculate the spinors
  wavebar_.resize(2);
  wave_   .resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wavebar_[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[1],ihel,Helicity::outgoing);
    wave_   [ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[2],ihel,Helicity::outgoing);
  }
  // now compute the matrix element
  Complex ii(0.,1.);
  Lorentz5Momentum pff(momenta[1]+momenta[2]);
  pff.rescaleMass();
  Energy2 mff2(pff.mass()*pff.mass());
  // compute the prefactor
  complex<InvEnergy3> pre(coupling_[imode()]/mff2);
  // the VMD factor
  if(includeVMD_[imode()]>0) {
    Energy2 mrho2=VMDmass_[imode()]*VMDmass_[imode()];
    Energy2 mwrho=VMDmass_[imode()]*VMDwidth_[imode()];
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  LorentzVector<complex<Energy3> > eps;
  LorentzVector<complex<Energy> > fcurrent;
  // compute the matrix element
  vector<unsigned int> ispin(4);ispin[0]=0;
  for(ispin[3]=0;ispin[3]<2;++ispin[3]) {
    for(ispin[2]=0;ispin[2]<2;++ispin[2]) {
      fcurrent = wave_[ispin[3]].vectorCurrent(wavebar_[ispin[2]]);
      // compute the current for this part
      eps = epsilon(momenta[0],pff,fcurrent);
      for(ispin[1]=0;ispin[1]<3;++ispin[1]) {
	(*ME())(ispin) = Complex(pre *vectors_[ispin[1]].dot(eps));
      }
    }	  
  }
  double me = ME()->contract(rho_).real();
  //   //code to test the matrix element against the analytic result
  //   Energy   m[4]={part.mass(),momenta[0].mass(),momenta[1].mass(),momenta[2].mass()};
  //   Energy2 m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  //   Lorentz5Momentum p12=momenta[0]+momenta[1];p12.rescaleMass();
  //   Energy2 m122(p12.mass2());
  //   Complex output( ((pre*conj(pre)).real()*(
  // 				 -2*m122*m122*mff2 - mff2*mff2*mff2 + 
  // 				 m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
  // 					m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
  // 				 2*m[2]*(m[2]*m2[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
  // 				 m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
  // 				 m2[0]*m2[0] +mff2*mff2*
  // 				 (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
  // 				 mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
  // 				       2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
  // 				       m2[0]*m2[0]) + 2*m122*
  // 				 (-mff2*mff2 - (m2[2] - m2[3])*
  // 				  (m2[1] - m2[0]) + 
  // 				  mff2*(m2[1] + m2[2] + m2[3] + 
  // 					m2[0])))));
  //   cout << "testing the matrix element " 
  //        << real(output) << " " << me << " " << test2 << endl;
  return me;
}

// method to return an object to calculate the 3 or higher body partial width
WidthCalculatorBasePtr 
PScalarVectorFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id()),idf[2]={0,0},idv(0);
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit) {
    if((**pit).iSpin()==PDT::Spin1){idv=(**pit).id();}
    else{idf[nf]=(**pit).id();++nf;}
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==idv) {
      if((idf[0]==outgoing_[ix].second&&idf[1]==-outgoing_[ix].second)||
	 (idf[1]==outgoing_[ix].second&&idf[0]==-outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  // get the masses we need
  Energy m[3]={getParticleData( outgoing_[imode].first )->mass(),
	       getParticleData( outgoing_[imode].second)->mass(),
	       getParticleData(-outgoing_[imode].second)->mass()};
  return 
    new_ptr(ThreeBodyAllOn1IntegralCalculator<PScalarVectorFermionsDecayer>
	    (3,-1000.*MeV,-0.9*MeV,-0.9,*this,imode,m[0],m[1],m[2]));
}

InvEnergy PScalarVectorFermionsDecayer::threeBodydGammads(const int imodeb,
							  const Energy2 q2, 
							  const Energy2 mff2, 
							  const Energy m1,
							  const Energy m2,
							  const Energy m3) const {
  // the masses of the external particles
  Energy q=sqrt(q2);
  Energy2 m12=m1*m1;
  Energy2 m22=m2*m2;
  Energy2 m32=m3*m3;
  // calculate the prefactor
  Complex ii(0.,1.);
  complex<InvEnergy3> pre = coupling_[imodeb] / mff2;
  // the VMD factor
  if(includeVMD_[imodeb]>0) {
    Energy2 mrho2=VMDmass_[imodeb]*VMDmass_[imodeb];
    Energy2 mwrho=VMDmass_[imodeb]*VMDwidth_[imodeb];
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  InvEnergy6 factor=real(pre*conj(pre));
  // compute the pieces from the integration limits
  Energy mff=sqrt(mff2);
  Energy e2star = 0.5*(mff2-m32+m22)/mff;
  Energy e1star = 0.5*(q2-mff2-m12)/mff;
  Energy e1sm = sqrt(e1star*e1star-m12);
  Energy e2sm = sqrt(e2star*e2star-m22);
  Energy2 a = 2*e1star*e2star+m12+m22;
  Energy2 b = 2*e1sm*e2sm;
  // term independent of s3
  Energy8 me = 2*b*(2*(m12*(mff2*mff2 + 4*mff2*m2*m3 -(m22 - m32)*(m22 - m32)) + 
		       2*m2*(m12 +m22)*m3*(-mff2 +m22 + q2))
		    +(m12 +m22)*(m12 +m22)*(-mff2 +m22 - 2*m2*m3 - m32)
		    -(mff2 +m22 + 2*m2*m3 - m32)*(-mff2 +m22 + q2)*(-mff2 +m22 + q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2 - (m22 - m32)*(m12 - q2) + 
		  mff2*(m12 + m22 + m32 + q2)));
  // quadratic term
  me+=-4.*mff2*b*(3.*a*a+b*b)/3.;

  // phase space factors
  using Constants::pi;
  return -factor * me/256./pi/pi/pi/q2/q;
}

// output the setup information for the particle database
void PScalarVectorFermionsDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*GeV << " " << includeVMD_[ix] << " "
	   << VMDid_[ix] << " " << VMDmass_[ix]/GeV << " " << VMDwidth_[ix]/GeV << " "
	   << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string PScalarVectorFermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
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
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  // vmd option  
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDinclude = stoi(stype);
  // vmd id
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int VMDid = stoi(stype);
  // vmd mass
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double VMDmass = stof(stype);
  // vmd width
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double VMDwidth = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  includeVMD_.push_back(VMDinclude);
  VMDid_.push_back(VMDid);
  VMDmass_.push_back(VMDmass*GeV);
  VMDwidth_.push_back(VMDwidth*GeV);
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./PScalarVectorVectorDecayer.cc"
// -*- C++ -*-
//
// PScalarVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PScalarVectorVectorDecayer class.
//
#include "PScalarVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PScalarVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      maxweight_[ix] = mode(ix)->maxWeight();
  }
}

void PScalarVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() || isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in PScalarVectorVectorDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  vector<double> wgt;
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  = getParticleData( incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int PScalarVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  unsigned int ix(0);
  int imode(-1);
  do {
    if(incoming_[ix]==id0) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==incoming_[ix]&&imode<0) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

void PScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/GeV) << incoming_ << outgoing_ << maxweight_;
}

void PScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/GeV) >> incoming_ >> outgoing_ >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PScalarVectorVectorDecayer,DecayIntegrator>
describeHerwigPScalarVectorVectorDecayer("Herwig::PScalarVectorVectorDecayer", "HwSMDecay.so");

void PScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<PScalarVectorVectorDecayer> documentation
    ("The PScalarVectorVectorDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");

  static Command<PScalarVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &PScalarVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<PScalarVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<PScalarVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in PScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");
}

void PScalarVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  bool photon[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->id()==ParticleID::gamma;
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,photon[ix]);
}

double PScalarVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  bool photon[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],ihel,Helicity::outgoing);
    }
  }
  // now compute the matrix element
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon[0] && ix==1) continue;
    for(unsigned int iy=0;iy<3;++iy) {
      if(photon[1] && iy==1) continue;
      (*ME())(0,ix,iy)=Complex(fact*epsilon(vectors_[0][ix],momenta[1],
					    vectors_[1][iy])*momenta[0]);
    }
  }
  double output = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 2.*sqr(fact*part.mass())*
  //   sqr(Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),momenta[1].mass()));
  // cerr << "testing the matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << output << " " << (output-test)/(output+test) << "\n";
  return output;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool PScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  bool order(true);
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
  while(imode<0&&ix<incoming_.size());
  coupling=coupling_[imode]*dm.parent()->mass();
  itype = 3;
  return order;
}

// output the setup info for the particle database
void PScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]*GeV << " " << maxweight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PScalarVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
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
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./ScalarMesonTensorScalarDecayer.cc"
// -*- C++ -*-
//
// ScalarMesonTensorScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonTensorScalarDecayer class.
//

#include "ScalarMesonTensorScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarMesonTensorScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) 
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void ScalarMesonTensorScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "ScalarMesonTensorScalarDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    PhaseSpaceModePtr mode;
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int ScalarMesonTensorScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					       const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id0   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==incoming_[ix]&&imode<0) {
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

void ScalarMesonTensorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/GeV) << incoming_ << outgoing_ << maxWeight_;
}

void ScalarMesonTensorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/GeV) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarMesonTensorScalarDecayer,DecayIntegrator>
describeHerwigScalarMesonTensorScalarDecayer("Herwig::ScalarMesonTensorScalarDecayer", "HwSMDecay.so");

void ScalarMesonTensorScalarDecayer::Init() {

  static ClassDocumentation<ScalarMesonTensorScalarDecayer> documentation
    ("The ScalarMesonTensorScalarDecayer class is designed for"
     " the decay of a pseduoscalar meson to two spin-1 particles.");
  
  static Command<ScalarMesonTensorScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, tensor, scalar, coupling (1/GeV) and max weight for a decay",
     &ScalarMesonTensorScalarDecayer::setUpDecayMode, false);
  
  static Deleted<ScalarMesonTensorScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in ScalarMesonTensorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarMesonTensorScalarDecayer> interfaceOutcomingT
    ("OutgoingT","The old methods of setting up a decay in ScalarMesonTensorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarMesonTensorScalarDecayer> interfaceOutcomingS
    ("OutgoingS","The old methods of setting up a decay in ScalarMesonTensorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarMesonTensorScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in ScalarMesonTensorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarMesonTensorScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in ScalarMesonTensorScalarDecayer have been deleted, please use SetUpDecayMode");
}

void ScalarMesonTensorScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  TensorWaveFunction::constructSpinInfo(tensors_,decay[0],
					outgoing,true,false);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

double ScalarMesonTensorScalarDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin2,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel);
    tensors_[ihel] = twave.wave();
  }
  // calculate the matrix element
  InvEnergy2 fact(coupling_[imode()]/part.mass());
  LorentzPolarizationVectorE vtemp;
  for(unsigned int ix=0;ix<5;++ix) {
    vtemp = tensors_[ix]*part.momentum(); 
    (*ME())(0,ix,0) = Complex(fact * momenta[1].dot(vtemp));
  }
  double me = ME()->contract(rho_).real();
  // // test of the matrix element
  // Energy pcm = Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					     momenta[1].mass());
  // double test = 2.*pow<4,1>(pcm)*sqr(coupling_[imode()]*part.mass())/
  //   3./pow<4,1>(momenta[0].mass());
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << (me-test)/(me+test) << "\n";
  // output the answer
  return me;
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarMesonTensorScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
						   double & coupling) const {
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
  int imode(-1);
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
  itype = 11;
  return order;
}

// output the setup information for the particle database
void ScalarMesonTensorScalarDecayer::dataBaseOutput(ofstream & output,
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

string ScalarMesonTensorScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
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
#line 1 "./ScalarScalarScalarDecayer.cc"
// -*- C++ -*-
//
// ScalarScalarScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarScalarScalarDecayer class.
//

#include "ScalarScalarScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarScalarScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix) ) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void ScalarScalarScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() || isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in ScalarScalarScalarDecayer" 
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

int ScalarScalarScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  if(children.size()!=2) return -1;
  int id0(parent->id());
  int id0bar = parent->CC() ? parent->CC()->id() : id0;
  int id1(children[0]->id());
  int id1bar = children[0]->CC() ? children[0]->CC()->id() : id1;
  int id2(children[1]->id());
  int id2bar = children[1]->CC() ? children[1]->CC()->id() : id2;
  // loop over the modes and see if this is one of them
  unsigned int ix(0);
  int imode(-1);
  do {
    if(id0   ==incoming_[ix]) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) {
	imode=ix;
	cc=false;
      }
    }
    if(id0bar==incoming_[ix]&&imode<0) {
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

void ScalarScalarScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,MeV) << incoming_ << outgoing_ << maxWeight_;
}

void ScalarScalarScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,MeV) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarScalarScalarDecayer,DecayIntegrator>
describeHerwigScalarScalarScalarDecayer("Herwig::ScalarScalarScalarDecayer", "HwSMDecay.so");

void ScalarScalarScalarDecayer::Init() {

  static ClassDocumentation<ScalarScalarScalarDecayer> documentation
    ("The ScalarScalarScalarDecayer class is designed for the"
     " decay of a scalar meson to two scalar mesons including off-shell effects");

  static Command<ScalarScalarScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing scalars, coupling(MeV) and max weight for a decay",
     &ScalarScalarScalarDecayer::setUpDecayMode, false);
  
  static Deleted<ScalarScalarScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarScalarScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in ScalarScalarScalarDecayer have been deleted, please use SetUpDecayMode");

}

void ScalarScalarScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double ScalarScalarScalarDecayer::me2(const int,const Particle & part,
				      const tPDVector &,
				      const vector<Lorentz5Momentum> &,
				      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  double fact(coupling_[imode()]/part.mass());
  (*ME())(0,0,0) = fact;
  return sqr(fact);
}

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarScalarScalarDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
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
  bool order(true);
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
  itype = 6;
  return order;
}

// output the setup information for the particle database
void ScalarScalarScalarDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second  << " "
	   << coupling_[ix]/MeV << " " << maxWeight_[ix]  << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string ScalarScalarScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
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
  coupling_.push_back(g*MeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./SemiLeptonicScalarDecayer.cc"
// -*- C++ -*-
//
// SemiLeptonicScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicScalarDecayer class.
//

#include "SemiLeptonicScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SemiLeptonicScalarDecayer::SemiLeptonicScalarDecayer() {
  // intermediates
  generateIntermediates(true);
}

void SemiLeptonicScalarDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _maxwgt.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _maxwgt.push_back(mode(ix)->maxWeight());
    }
  }
}

void SemiLeptonicScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // and the form factors
  _form->init();
  _modemap.clear();
  for(unsigned int ix=0;ix<_form->numberOfFactors();++ix) {
    // get the external particles for this mode
    int id0(0),id1(0);
    _form->particleID(ix,id0,id1);
    tPDPtr  in = getParticleData(id0);
    tPDPtr out = getParticleData(id1);
    _modemap.push_back(numberModes());
    if(!in || !out) continue;
    int Wcharge =(in->iCharge()-out->iCharge());
    Energy min = in->mass()+in->widthUpCut()
      -out->mass()+out->widthLoCut();
    for(unsigned int iy=0;iy<_current->numberOfModes();++iy) {
      int iq(0),ia(0);
      _current->decayModeInfo(iy,iq,ia);
      tPDVector outV = {out};
      tPDVector ptemp=_current->particles(Wcharge,iy,iq,ia);
      outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
      // create the mode
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,outV,1.));
      // create the first piece of the channel
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      // and the rest
      bool done = _current->createMode(Wcharge,tcPDPtr(),FlavourInfo(),iy,mode,1,0,channel,min);
      if(done) {
	// the maximum weight
	double maxweight = _maxwgt.size()>numberModes() ? _maxwgt[numberModes()] : 2.;
	mode->maxWeight(maxweight);
	addMode(mode);
      }
    }
  }
}

bool SemiLeptonicScalarDecayer::accept(tcPDPtr parent, 
				       const tPDVector & children) const {
  // find the non-lepton
  int imes(0),idtemp,idin(parent->id());
  vector<int> idother; bool dummy;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) imes=idtemp;
    else               idother.push_back(idtemp);
  }
  // check that the form factor exists
  if(_form->formFactorNumber(idin,imes,dummy)<0) return false;
  // and the current
  return _current->accept(idother);
}

int  SemiLeptonicScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					   const tPDVector & children) const {
  // find the ids of the particles for the decay current
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,imes(0),idin(parent->id());
  vector<int> idother;
  cc=false;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) imes=idtemp;
    else               idother.push_back(idtemp);
  }
  return _modemap[_form->formFactorNumber(idin,imes,cc)]
    +_current->decayMode(idother);  
}


void SemiLeptonicScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap;
}

void SemiLeptonicScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SemiLeptonicScalarDecayer,DecayIntegrator>
describeHerwigSemiLeptonicScalarDecayer("Herwig::SemiLeptonicScalarDecayer", "HwSMDecay.so");

void SemiLeptonicScalarDecayer::Init() {

  static ClassDocumentation<SemiLeptonicScalarDecayer> documentation
    ("The SemiLeptonicScalarDecayer class is designed for the"
    "semi-leptonic decay of a (pseudo)-scalar meson.");

  static Reference<SemiLeptonicScalarDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &SemiLeptonicScalarDecayer::_current, true, true, true, false, false);

  static Reference<SemiLeptonicScalarDecayer,ScalarFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor",
     &SemiLeptonicScalarDecayer::_form, true, true, true, false, false);

  static ParVector<SemiLeptonicScalarDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weights for the decays",
     &SemiLeptonicScalarDecayer::_maxwgt,
     0, 0, 0, 0, 100., false, false, true);

}

void SemiLeptonicScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin0)
    ScalarWaveFunction::
      constructSpinInfo(decay[0],outgoing,true);
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin1)
    VectorWaveFunction::
      constructSpinInfo(_vectors,decay[0],outgoing,true,false);
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin2)
    TensorWaveFunction::
      constructSpinInfo(_tensors,decay[0],outgoing,true,false);
  // and the stuff from the current
  _current->constructSpinInfo(ParticleVector(decay.begin()+1,decay.end()));
}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicScalarDecayer::me2(const int , const Particle & part,
				      const tPDVector & outgoing,
				      const vector<Lorentz5Momentum> & momenta,
				      MEOption meopt) const {
  // get the information on the form-factor
  int jspin(0),id0(part.id()),id1(outgoing[0]->id());
  bool cc(false);
  unsigned int iloc(_form->formFactorNumber(id0,id1,cc));
  int spect,iq,ia;
  _form->formFactorInfo(iloc,jspin,spect,iq,ia);
  if(!ME()) {
    if(jspin==0)
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
    else if(jspin==1)       
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
    else if(jspin==2)       
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half)));
  }
  // initialisation
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    // work out the mapping for the lepton vector
    _constants.resize(outgoing.size()+1);
    _ispin.resize(outgoing.size());
    _imes=0;
    unsigned int itemp(1);
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      _ispin[ix]=outgoing[ix]->iSpin();
      if(abs(outgoing[ix]->id())<=16) {
  	itemp*=_ispin[ix];
  	_constants[ix]=itemp;
      }
      else _imes=ix;
    }
    _constants[outgoing.size()]=1;
    _constants[_imes]=_constants[_imes+1];
  }
  // get the wavefunctions of the decay products
  switch(outgoing[0]->iSpin()) {
  case PDT::Spin0:
    break;
  case PDT::Spin1:
    _vectors.resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel)
      _vectors[ihel] = HelicityFunctions::polarizationVector(-momenta[0],ihel,
							     Helicity::outgoing);
    break;
  case PDT::Spin2:
    {
      TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
      _tensors.resize(5);
      for(unsigned int ihel=0;ihel<5;++ihel) {
	twave.reset(ihel);
	_tensors[ihel] = twave.wave();
      }
    }
    break;
  default:
    assert(false);
  }
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the hadronic current for the decay
  Complex ii(0.,1.);
  vector<LorentzPolarizationVectorE> hadron;
  if(jspin==0) {
    Complex fp,f0;
    _form->ScalarScalarFormFactor(q2,iloc,id0,id1,part.mass(),momenta[0].mass(),
  				  f0,fp);
    Complex pre((sqr(part.mass())-sqr(momenta[0].mass()))/q2*(f0-fp));
    hadron.push_back(fp*sum+(pre*q));
  }
  else if(jspin==1) {
    Complex A0,A1,A2,A3,V;
    complex<Energy> dot;
    Energy MP(part.mass()),MV(momenta[0].mass()),msum(MP+MV),mdiff(MP-MV);
    _form->ScalarVectorFormFactor(q2,iloc,id0,id1,MP,MV,A0,A1,A2,V);
    A3 = Complex(0.5/MV*(msum*A1-mdiff*A2));
    if(cc) V*=-1.;
    // compute the hadron currents
    for(unsigned int ix=0;ix<3;++ix) {
      // dot product
      dot = _vectors[ix]*part.momentum();
      // current
      hadron.push_back(-ii*msum*A1*_vectors[ix]
  		       +ii*A2/msum*dot*sum
  		       +2.*ii*MV/q2*(A3-A0)*dot*q
  		       +2.*V/msum*Helicity::epsilon(_vectors[ix],part.momentum(),
  						    momenta[0]));
    }
  }
  else if(jspin==2) {
    complex<InvEnergy2> h,bp,bm;
    complex<double> k;
    complex<Energy2> dot;
    _form->ScalarTensorFormFactor(q2,iloc,id0,id1,part.mass(),momenta[0].mass(),
  				  h,k,bp,bm);
    if(!cc) h*=-1.;
    LorentzPolarizationVectorE dotv;
    // compute the hadron currents
    for(unsigned int ix=0;ix<5;++ix) {
      dotv = _tensors[ix]*part.momentum();
      dot = dotv*part.momentum();
      hadron.push_back(ii*h*Helicity::epsilon(dotv,sum,q)
  		       -k*dotv-bp*dot*sum-bm*dot*q);
    }
  }
  Energy scale;
  int mode=(abs(outgoing[1]->id())-11)/2;
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(tcPDPtr(),FlavourInfo(),
			     mode,-1,scale,tPDVector(outgoing.begin()+1,outgoing.end()),
			     vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),
			     meopt));
  // compute the matrix element
  vector<unsigned int> ihel(outgoing.size()+1);
  for(unsigned int mhel=0;mhel<hadron.size();++mhel) {
    for(unsigned int lhel=0;lhel<lepton.size();++lhel) {
      // map the index for the leptons to a helicity state
      for(unsigned int ix=outgoing.size();ix>0;--ix) {
  	if(ix-1!=_imes) ihel[ix]=(lhel%_constants[ix-1])/_constants[ix];
      }
      // helicities of mesons
      ihel[0]=0;
      ihel[_imes+1]=mhel;
      (*ME())(ihel) = Complex(lepton[lhel].dot(hadron[mhel])*SM().fermiConstant());
    }
  }
  // store the matrix element
  double ckm(1.);
  if(iq<=6) {
    if(iq%2==0) ckm = SM().CKM(abs(iq)/2-1,(abs(ia)-1)/2);
    else        ckm = SM().CKM(abs(ia)/2-1,(abs(iq)-1)/2);
  }
  // return the answer
  return 0.5*(ME()->contract(_rho)).real()*ckm; 
}
 
// output the setup information for the particle database
void SemiLeptonicScalarDecayer::dataBaseOutput(ofstream & output,
					       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    output << "insert " << name() << ":MaximumWeight " << ix << " " 
	   << _maxwgt[ix] << "\n";
  }
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":Current " << _current->name() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":FormFactor " << _form->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./ScalarMesonFactorizedDecayer.cc"
// -*- C++ -*-
//
// ScalarMesonFactorizedDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarMesonFactorizedDecayer class.
//

#include "ScalarMesonFactorizedDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;


ScalarMesonFactorizedDecayer::ScalarMesonFactorizedDecayer() 
// default values of the couplings (taken from ZPC34, 103)
  : _a1b(1.10), _a2b(-0.24), _a1c(1.30), _a2c(-0.55) { 
  // intermediates
  generateIntermediates(true);
}

void ScalarMesonFactorizedDecayer::rebind(const TranslationMap & trans)
  {
  _ckm = trans.translate(_ckm);
  DecayIntegrator::rebind(trans);
}

IVector ScalarMesonFactorizedDecayer::getReferences() {
  IVector ret = DecayIntegrator::getReferences();
  ret.push_back(_ckm);
  return ret;
}

void ScalarMesonFactorizedDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the ckm object
  _ckm=dynamic_ptr_cast<Ptr<StandardCKM>::pointer>(SM().CKM());
  if(!_ckm) throw InitException() << "ScalarMesonFactorizedDecayer::doinit() "
  				  << "the CKM object must be the Herwig one"
  				  << Exception::runerror;
  // get the CKM matrix (unsquared for interference)
  Complex ckmmat[3][3];
  vector< vector<Complex > > CKM(_ckm->getUnsquaredMatrix(SM().families()));
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy){
      ckmmat[ix][iy]=CKM[ix][iy];
    }
  }
  // make sure the currents and form factors got initialised
  for(unsigned int ix=0;ix<_current.size();++ix)
    _current[ix]->init();
  for(unsigned int ix=0;ix<_form.size();++ix)
    _form[ix]->init();
  // find all the possible modes
  vector<unsigned int> tformmap[2],tcurrmap[2];
  vector<int> inquark,outquark,currq,curra;
  tPDVector incoming;
  vector<tPDVector> outgoing;
  // loop over the modes in the form factors and currents
  for(unsigned int iform=0;iform<_form.size();++iform) {
    for(unsigned int ix=0;ix<_form[iform]->numberOfFactors();++ix) {
      // information from the form-factor
      int id0,id1,jspin,spect,inq,outq;
      _form[iform]->particleID(ix,id0,id1);
      _form[iform]->formFactorInfo(ix,jspin,spect,inq,outq);
      // particles from the form factor
      tPDPtr in  = getParticleData(id0);
      tPDPtr out = getParticleData(id1);
      // charge of the decay products
      int Wcharge = in->iCharge()-out->iCharge();
      // max mass for the particles in the current
      Energy min = in->massMax()-out->massMin();
      for(unsigned int icurr=0;icurr<_current.size();++icurr) {
  	for(unsigned int iy=0;iy<_current[icurr]->numberOfModes();++iy) {
	  // get the particles from the current
	  int iq,ia;
	  _current[icurr]->decayModeInfo(iy,iq,ia);
	  tPDVector ptemp=_current[icurr]->particles(Wcharge,iy,iq,ia);
	  tPDVector outV = {out};
	  outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
	  // create the mode
	  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,outV,1.));
	  // create the first piece of the channel
	  PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1)); 	  
   	  Energy minb=ZERO;
   	  for(unsigned int iz=0;iz<ptemp.size();++iz)
	    minb += ptemp[iz]->massMin();
  	  // add this mode to the list
   	  if(outV.size()>1&&minb<min&&
  	     (Wcharge!=0||(Wcharge==0&&((inq>0&&inq%2!=iq%2)||
   					(inq<0&&abs(inq)%2!=abs(ia)%2))))) {
   	    tformmap[0].push_back(iform);tformmap[1].push_back(ix);
   	    tcurrmap[0].push_back(icurr);tcurrmap[1].push_back(iy);
 	    incoming.push_back( in );
 	    outgoing.push_back(outV);
 	    inquark.push_back(inq);outquark.push_back(outq);
 	    currq.push_back(iq);curra.push_back(ia);
 	  }
   	  // if the meson in the current is neutral try the CC mode
   	  if(Wcharge==0&&iq!=-ia&&((inq>0&&inq%2!=iq%2)||
   				   (inq<0&&abs(inq)%2!=abs(ia)%2))) {
   	    // get the particles from the current
   	    tPDVector ptemp=_current[icurr]->particles(Wcharge,iy,-ia,-iq);
	    outV = {out};
	    outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
   	    minb=ZERO;
   	    for(unsigned int iz=0;iz<ptemp.size();++iz)
	      minb+=ptemp[iz]->massMin();
   	    if(outV.size()>1&&minb<min) {
   	      tformmap[0].push_back(iform);tformmap[1].push_back(ix);
   	      tcurrmap[0].push_back(icurr);tcurrmap[1].push_back(iy);
	      incoming.push_back(in);
   	      outgoing.push_back(outV);
   	      inquark.push_back(inq);outquark.push_back(outq);
   	      currq.push_back(-ia);curra.push_back(-iq);
   	    }
   	  }
  	}
      }
    }
  }
  // loop over the modes and find the dupliciates
  static const double ort(sqrt(0.5));
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    while (true) {
      if(outgoing[ix].empty()) break;
      vector<bool> modecc;
      vector<unsigned int> modeloc;
      findModes(ix,incoming,outgoing,modeloc,modecc);
      // if more than two outgoing only allow one diagram
      if ( outgoing[ix].size()>2 && !modeloc.empty() ) break;
      // create the mode and set the particles as for the first instance
      PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],1.));
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      Energy min = incoming[ix]->massMax()-outgoing[ix][0]->massMin();
      int Wcharge = incoming[ix]->iCharge()-outgoing[ix][0]->iCharge();
      bool done = _current[tcurrmap[0][ix]]->
	createMode(Wcharge,tcPDPtr(),FlavourInfo(),tcurrmap[1][ix],mode,1,0,channel,min);
      if(!done) throw InitException() << "Failed to construct mode in "
   				      << "ScalarMesonFactorizedDecayer::doinit()." 
   				      << Exception::abortnow;
      // set the parameters for the additional modes
      vector<unsigned int> tformpart(1,0),ttform[2],ttcurr[2];
      ttform[0].push_back(tformmap[0][ix]);ttform[1].push_back(tformmap[1][ix]);
      ttcurr[0].push_back(tcurrmap[0][ix]);ttcurr[1].push_back(tcurrmap[1][ix]);
      int id    = outgoing[ix][0]->id();
      int idbar = outgoing[ix][0]->CC() ? outgoing[ix][0]->CC()->id() : id;
      for(unsigned int iy=0;iy<modeloc.size();++iy) {
   	ttform[0].push_back(tformmap[0][modeloc[iy]]);
   	ttform[1].push_back(tformmap[1][modeloc[iy]]);
   	ttcurr[0].push_back(tcurrmap[0][modeloc[iy]]);
   	ttcurr[1].push_back(tcurrmap[1][modeloc[iy]]);
   	unsigned int iz=0;
   	do {
   	  if(( modecc[iy]&&outgoing[modeloc[iy]][iz]->id()==idbar)||
   	     (!modecc[iy]&&outgoing[modeloc[iy]][iz]->id()==id))
   	    tformpart.push_back(iz);
   	  ++iz;
   	}
   	while(tformpart.size()!=iy+2&&iz<2);
      }
      // calculate ckm factors
      vector<Complex> tCKM;
      for(unsigned int iy=0;iy<ttcurr[0].size();++iy) {
  	// get the quarks involved in the process
	int iq,ia,inq,outq;
  	if(iy==0) {
   	  iq=currq[ix];ia=curra[ix];
   	  inq=inquark[ix];outq=outquark[ix];
	}
 	else {
 	  if(!modecc[iy-1]) {
 	    iq=currq[modeloc[iy-1]];ia=curra[modeloc[iy-1]];
 	    inq=inquark[modeloc[iy-1]];outq=outquark[modeloc[iy-1]];
 	  }
 	  else {
 	    ia=-currq[modeloc[iy-1]];iq=-curra[modeloc[iy-1]];
 	    inq=-inquark[modeloc[iy-1]];outq=-outquark[modeloc[iy-1]];
 	  }
 	}
	int id0,id1;
	_form[ttform[0][iy]]->particleID(ttform[1][iy],id0,id1);
  	int Wcharge = getParticleData(id0)->iCharge()-getParticleData(id1)->iCharge();
   	Complex ckm=1.;
   	if(Wcharge!=0) {
   	  if(abs(iq)%2==0)  ckm *= conj(ckmmat[abs(iq)/2-1][(abs(ia)-1)/2]);
   	  else              ckm *= conj(ckmmat[abs(ia)/2-1][(abs(iq)-1)/2]);
   	  if(abs(inq)%2==0) ckm *= ckmmat[abs(inq)/2-1][(abs(outq)-1)/2];
   	  else              ckm *= ckmmat[abs(outq)/2-1][(abs(inq)-1)/2];
   	  if(abs(inq)==5)   ckm*=_a1b;
   	  else              ckm*=_a1c;
   	}
   	else {
   	  if(inq>0) {
   	    if(abs(inq)%2==0)  ckm *= ckmmat[abs(inq)/2-1][(abs(iq)-1)/2];
   	    else               ckm *= ckmmat[abs(iq)/2-1][(abs(inq)-1)/2];
   	    if(abs(outq)%2==0) ckm *= conj(ckmmat[abs(outq)/2-1][(abs(ia)-1)/2]);
   	    else               ckm *= conj(ckmmat[abs(ia)/2-1][(abs(outq)-1)/2]);
   	  }
   	  else {
   	    if(abs(inq)%2==0)  ckm *= ckmmat[abs(inq)/2-1][(abs(ia)-1)/2];
   	    else               ckm *= ckmmat[abs(ia)/2-1][(abs(inq)-1)/2];
   	    if(abs(outq)%2==0) ckm *= conj(ckmmat[abs(outq)/2-1][(abs(iq)-1)/2]);
   	    else               ckm *= conj(ckmmat[abs(iq)/2-1][(abs(outq)-1)/2]);
   	  }
   	  if(abs(inq)==5) ckm*=_a2b;
   	  else            ckm*=_a2c;
   	}
   	if((abs(inq)%2==0&&inq<0)||(abs(inq)%2!=0&&inq>0)){ckm=conj(ckm);}
   	tCKM.push_back(ckm);
      }
      // special if the particles are idential add additional modes and 
      // identical particle factors
      if(outgoing[ix][0]->id()==outgoing[ix][1]->id()&&outgoing[ix].size()==2) {
  	unsigned int isize=ttcurr[0].size();
   	for(unsigned int iy=0;iy<isize;++iy) {
  	  ttcurr[0].push_back(ttcurr[0][iy]);ttcurr[1].push_back(ttcurr[1][iy]);
  	  ttform[0].push_back(ttform[0][iy]);ttform[1].push_back(ttform[1][iy]);
  	  if(tformpart[iy]==0){tformpart.push_back(1);}
  	  else{tformpart.push_back(0);}
  	  tCKM[iy]*=ort;tCKM.push_back(tCKM[iy]);
  	}
      }
      // add the parameters for the mode to the list
      _currentmapA.push_back(ttcurr[0]);_currentmapB.push_back(ttcurr[1]);
      _formmapA.push_back(ttform[0]);_formmapB.push_back(ttform[1]);
      _formpart.push_back(tformpart);
      _CKMfact.push_back(tCKM);
      // add the mode to the list
      double maxweight(0.);
      if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
      // the weights for the channels
      vector<double> channelwgts;
     if(_wgtloc.size()>numberModes()&&
 	 _wgtloc[numberModes()]+mode->channels().size()<=_weights.size()) {
       vector<double>::iterator start=_weights.begin()+_wgtloc[numberModes()];
       vector<double>::iterator end  = start+mode->channels().size();
 	channelwgts=vector<double>(start,end);
     }
     else {
 	channelwgts.resize(mode->channels().size(),1./(mode->channels().size()));
     }
     // don't need channels for two body decays
     if(outgoing[ix].size()==2) {
       channelwgts.clear();
       mode = new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],maxweight));
     }
     else {
       mode->maxWeight(maxweight);
       mode->setWeights(channelwgts);
     }
     addMode(mode);
     // resize the duplicate modes to remove them
     for(unsigned int iy=0;iy<modeloc.size();++iy) outgoing[modeloc[iy]] = tPDVector();
     break;
    }
  }
}

void ScalarMesonFactorizedDecayer::doinitrun() {
  unsigned int ix,iy;
  for(ix=0;ix<_current.size();++ix) _current[ix]->initrun();
  for(ix=0;ix<_form.size();++ix)    _form[ix]->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _weights.clear();
    _wgtloc.clear();
    _wgtmax.clear();
    for(ix=0;ix<numberModes();++ix) {
      _wgtmax.push_back(mode(ix)->maxWeight());
      _wgtloc.push_back(_weights.size());
      for(iy=0;iy<mode(ix)->channels().size();++iy) {
	_weights.push_back(mode(ix)->channels()[iy].weight());
      }
    }
  }
}

bool ScalarMesonFactorizedDecayer::accept(tcPDPtr parent,
					  const tPDVector & children) const {
  // N.B. this is a necessary but not sufficient test
  bool allowed(false),dummy;
  // find the ids of the particles
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  vector<int> ids,idcurr;
  int id(parent->id());
  for( ; pit!=pend;++pit) ids.push_back((**pit).id());
  // loop over the possible particles in the formfactor
  unsigned int ipart(0),iform,icurr,ix;
  do {
    idcurr.clear();
    for(ix=0;ix<ids.size();++ix){if(ix!=ipart){idcurr.push_back(ids[ix]);}}
    iform=0;
    do {
      // check if possible from the form factor
      if(_form[iform]->formFactorNumber(id,ids[ipart],dummy)>=0) {
	// check if possible from the current
	icurr=0;
	do {
	  allowed=_current[icurr]->accept(idcurr);
	  ++icurr;
	}
	while(!allowed&&icurr<_current.size());
      }
      ++iform;
    }
    while(!allowed&&iform<_form.size());
    ++ipart;
  }
  while(!allowed&&ipart<ids.size());
  return allowed;
}

int ScalarMesonFactorizedDecayer::modeNumber(bool & cc,tcPDPtr parent,
					     const tPDVector & children) const {
  int imode(-1);
  // id's of the particles and CC
  // of the parent
  int id0(parent->id()),id0bar(id0);
  if(parent->CC())  id0bar = parent->CC()->id();
  // of the products
  vector<int> ids,idbars;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ;pit!=pend;++pit) {
    ids.push_back((**pit).id());
    if((**pit).CC()) idbars.push_back((**pit).CC()->id());
    else             idbars.push_back(ids.back());
  }
  // loop over the modes
  cc=false;
  unsigned int ix=0;
  do {
    // particle mode
    if(id0==mode(ix)->incoming().first->id()&&
       ids.size()==mode(ix)->outgoing().size()) {
      unsigned int nfound=0;
      vector<bool> done(ids.size(),false);
      for(unsigned int iy=0;iy<ids.size();++iy) {
   	int idtemp = mode(ix)->outgoing()[iy]->id();
	unsigned int iz=0;
	bool found=false;
	do{
	  if(idtemp==ids[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<ids.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==ids.size()) {
	cc=false;
	imode=ix;
      }
    }
    // CC mode
    if(id0bar==mode(ix)->incoming().first->id()&&
       ids.size()==mode(ix)->outgoing().size()) {
      unsigned int nfound=0;
      vector<bool> done(ids.size(),false);
      for(unsigned int iy=0;iy<idbars.size();++iy) {
  	int idtemp=mode(ix)->outgoing()[iy]->id();
	unsigned int iz=0;
	bool found=false;
	do {
	  if(idtemp==idbars[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<idbars.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==idbars.size()) {
	cc=true;
	imode=ix;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<numberModes());
  if(imode<0) {
    string mode = parent->PDGName() + "->";
    for(unsigned int ix=0;ix<children.size();++ix) 
      mode += children[ix]->PDGName() +",";
    throw DecayIntegratorError() << "Unable to find the mode " << mode << " in " 
				  << name() 
				  << " ScalarMesonFactorizedDecayer::decay()" 
				  << Exception::abortnow;
  }
  return imode;
}


void ScalarMesonFactorizedDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _ckm 
     << _a1b << _a2b << _a1c << _a2c 
     << _currentmapA << _currentmapB 
     << _formmapA << _formmapB << _formpart << _wgtloc 
     << _wgtmax << _weights << _CKMfact ;
}

void ScalarMesonFactorizedDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _ckm 
     >> _a1b >> _a2b >> _a1c >> _a2c 
     >> _currentmapA >> _currentmapB 
     >> _formmapA >> _formmapB >> _formpart >> _wgtloc
     >> _wgtmax >> _weights >> _CKMfact;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarMesonFactorizedDecayer,DecayIntegrator>
describeHerwigScalarMesonFactorizedDecayer("Herwig::ScalarMesonFactorizedDecayer", "HwSMDecay.so");

void ScalarMesonFactorizedDecayer::Init() {

  static ClassDocumentation<ScalarMesonFactorizedDecayer> documentation
    ("The ScalarMesonFactorizedDecayer class is designed for the weak decay of"
     " scalar mesons using the factorization approximation.");

  static RefVector<ScalarMesonFactorizedDecayer,WeakCurrent> interfaceCurrents
    ("Currents",
     "A vector of references to the currents",
     &ScalarMesonFactorizedDecayer::_current, -1, false, false, true, false, false);

  static RefVector<ScalarMesonFactorizedDecayer,ScalarFormFactor> interfaceFormFactors
    ("FormFactors",
     "A vector of references to the form-factors",
     &ScalarMesonFactorizedDecayer::_form, -1, false, false, true, false, false);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea1Bottom
    ("a1Bottom",
     "The factorization paramter a_1 for decays of bottom baryons",
     &ScalarMesonFactorizedDecayer::_a1b, 1.1, -10.0, 10.0,
     false, false, true);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea2Bottom
    ("a2Bottom",
     "The factorization paramter a_2 for decays of bottom baryons",
     &ScalarMesonFactorizedDecayer::_a2b, -0.24, -10.0, 10.0,
     false, false, true);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea1Charm
    ("a1Charm",
     "The factorization paramter a_1 for decays of charm baryons",
     &ScalarMesonFactorizedDecayer::_a1c, 1.3, -10.0, 10.0,
     false, false, true);

  static Parameter<ScalarMesonFactorizedDecayer,double> interfacea2Charm
    ("a2Charm",
     "The factorization paramter a_2 for decays of charm baryons",
     &ScalarMesonFactorizedDecayer::_a2c, -0.55, -10.0, 10.0,
     false, false, true);

  static ParVector<ScalarMesonFactorizedDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &ScalarMesonFactorizedDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<ScalarMesonFactorizedDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &ScalarMesonFactorizedDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<ScalarMesonFactorizedDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &ScalarMesonFactorizedDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);
}
 
void ScalarMesonFactorizedDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
  					incoming,true);
  // get the wavefunctions of the decay products
  for(unsigned int ix=0;ix<decay.size();++ix) {
    switch(decay[ix]->dataPtr()->iSpin()) {
    case PDT::Spin0:
      ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
      break;
    case PDT::Spin1:
      VectorWaveFunction::constructSpinInfo(_vectors[ix],decay[ix],outgoing,
					    true,false);
      break;
    case PDT::Spin2:
      TensorWaveFunction::constructSpinInfo(_tensors[ix],decay[ix],outgoing,
					    true,false);
      break;
    default:
      assert(false);
    }
  }
}

double ScalarMesonFactorizedDecayer::me2(const int ichan, const Particle & part,
					 const tPDVector & outgoing,
					 const vector<Lorentz5Momentum> & momenta,
					 MEOption meopt) const {
  if(!ME()) {
    // create the matrix element
    vector<PDT::Spin> spin;
    for(unsigned int ix=0;ix<outgoing.size();++ix)
      spin.push_back(outgoing[ix]->iSpin());
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,spin)));
  }
  // initialisation
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _vectors.resize(outgoing.size());
    _tensors.resize(outgoing.size());
  }
  // get the wavefunctions of the decay products
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    switch(outgoing[ix]->iSpin()) {
    case PDT::Spin0:
      break;
    case PDT::Spin1:
      _vectors[ix].resize(3);
      for(unsigned int ihel=0;ihel<3;++ihel)
	_vectors[ix][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],ihel,
								   Helicity::outgoing);
      break;
    case PDT::Spin2:
      {
	TensorWaveFunction twave(momenta[ix],outgoing[ix],Helicity::outgoing);
	_tensors[ix].resize(5);
	for(unsigned int ihel=0;ihel<5;++ihel) {
	  twave.reset(ihel);
	  _tensors[ix][ihel] = twave.wave();
	}
      }
    break;
    default:
      assert(false);
    }
  }
  ME()->zero();
  // find the mode
  unsigned int mode(imode());
  int id0(part.id());
  Complex ii(0.,1.);
  vector<unsigned int> ihel(outgoing.size());
  // loop over the different diagrams
  Energy MP(part.mass()),scale;
  double pre;
  for(unsigned int iy=0;iy<_CKMfact[mode].size();++iy) {
    Energy MV = momenta[_formpart[mode][iy]].mass();
    int id1   = outgoing[_formpart[mode][iy]]->id();
    int id0t,id1t;
    _form[_formmapA[mode][iy]]->particleID(_formmapB[mode][iy],id0t,id1t);
    bool cc(id0t!=id0);
    // calculate the form-factor part
    vector<LorentzPolarizationVectorE> form;
    Lorentz5Momentum q   = part.momentum()-momenta[_formpart[mode][iy]];
    q.rescaleMass();
    Lorentz5Momentum sum = part.momentum()+momenta[_formpart[mode][iy]];
    sum.rescaleMass();
    Energy2 q2=q.mass2();
    if(outgoing[_formpart[mode][iy]]->iSpin()==1) {
      Complex fp,f0;
      _form[_formmapA[mode][iy]]->ScalarScalarFormFactor(q2,_formmapB[mode][iy],
   							 id0,id1,MP,MV,f0,fp);
      pre=(MP*MP-MV*MV)/q2;
      form.push_back(fp*sum+pre*(f0-fp)*q);
    }
    else if(outgoing[_formpart[mode][iy]]->iSpin()==3) {
      Energy msum  = MP+MV;
      Energy mdiff = MP-MV;
      Complex A0,A1,A2,V;
       _form[_formmapA[mode][iy]]->ScalarVectorFormFactor(q2,_formmapB[mode][iy],id0,
							  id1,MP,MV,A0,A1,A2,V);
       if(cc) V=-V;
       Complex A3 = 0.5/MV*(msum*A1-mdiff*A2);
       // compute the hadron currents
       for(unsigned int ix=0;ix<3;++ix) {
	 // dot product
	 complex<Energy> dot = _vectors[_formpart[mode][iy]][ix]*part.momentum();
	 // current
	 form.push_back(-ii*msum*A1*_vectors[_formpart[mode][iy]][ix]
			+ii*A2/msum*dot*sum
			+2.*ii*MV/q2*(A3-A0)*dot*q
			+2.*V/msum*epsilon(_vectors[_formpart[mode][iy]][ix],
					   part.momentum(),
					   momenta[_formpart[mode][iy]])); 
       }
    }
    else if(outgoing[_formpart[mode][iy]]->iSpin()==5) {
      Complex k;
      complex<InvEnergy2> h,bp,bm;
      _form[_formmapA[mode][iy]]->ScalarTensorFormFactor(q2,_formmapB[mode][iy],
  							 id0,id1,MP,MV,h,k,bp,bm);
      if(cc) h=-h;
      // compute the hadron currents
      for(unsigned int ix=0;ix<5;++ix) {
	LorentzPolarizationVectorE dotv =
	  _tensors[_formpart[mode][iy]][ix]*part.momentum();
   	complex<Energy2> dot = dotv*part.momentum();
   	form.push_back(ii*h*epsilon(dotv,sum,q)-k*dotv
   		       -bp*dot*sum-bm*dot*q);
      }
    }
    // find the particles for the current
    tPDVector cpart;
    vector<Lorentz5Momentum> cmom;
    for(unsigned int ix=0;ix<outgoing.size();++ix) {
      if(ix!=_formpart[mode][iy]) {
	cpart.push_back(outgoing[ix]);
	cmom .push_back(momenta[ix]);
      }
    }
    unsigned int ix=outgoing.size();
    vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
    int itemp(1);
    do {
      --ix;
      if(ix!=_formpart[mode][iy]) {
	itemp*=outgoing[ix]->iSpin();
	constants[ix]=itemp;
      }
    }
    while(ix!=0);
    constants[outgoing.size()]=1;
    if(_formpart[mode][iy]!=outgoing.size())
      constants[_formpart[mode][iy]]=constants[_formpart[mode][iy]+1];
    // calculate the current
    vector<LorentzPolarizationVectorE>
      curr=_current[_currentmapA[mode][iy]]->
      current(tcPDPtr(),FlavourInfo(),_currentmapB[mode][iy],ichan,scale,cpart,cmom,meopt);
    pre = (pow(part.mass()/scale,int(cpart.size()-2)));
    // loop over the helicities to calculate the matrix element
    ihel[0]=0;
    for(unsigned int chel=0;chel<curr.size();++chel) {
      for(ix=outgoing.size();ix>0;--ix) {
	if(ix!=_formpart[mode][iy]+1)
	  ihel[ix]=(chel%constants[ix-1])/constants[ix];
      }
      for(unsigned int fhel=0;fhel<form.size();++fhel) {
	ihel[_formpart[mode][iy]+1]=fhel;
	(*ME())(ihel) += Complex(pre*_CKMfact[mode][iy]*
				 form[fhel].dot(curr[chel])*SM().fermiConstant());
      }
    }
  }
  // perform the contraction
  return 0.5*(ME()->contract(_rho)).real();
}
  
void ScalarMesonFactorizedDecayer::findModes(unsigned int imode,
					     tPDVector & incoming,
					     vector<tPDVector> & outgoing,
					     vector<unsigned int> & loc,
					     vector<bool> & cc) {
  // get the id's for the mode
  // incoming
  int id_in    = incoming[imode]->id();
  int idbar_in = incoming[imode]->CC() ?
    incoming[imode]->CC()->id() : incoming[imode]->id();
  // outgoing
  vector<int> id_out,idbar_out;
  for(unsigned int ix=0;ix<outgoing[imode].size();++ix) {
    id_out.push_back(outgoing[imode][ix]->id());
    if(outgoing[imode][ix]->CC())
      idbar_out.push_back(outgoing[imode][ix]->CC()->id());
    else
      idbar_out.push_back(id_out[ix]);
  }
  // loop over the modes
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    if(ix==imode||outgoing[ix].empty()) continue;
    assert(!outgoing[ix].empty());
    assert(incoming[ix]);
    // the particle mode
    if(incoming[ix]->id()==id_in&&outgoing[ix].size()==id_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound = 0;
      for(unsigned int iy=0;iy<id_out.size();++iy) {
	int idtemp=outgoing[ix][iy]->id();
     	unsigned int iz(0);
     	bool found=false;
    	do {
	  if(idtemp==id_out[iz]&&!done[iz]) {
     	    done[iz]=true;
     	    found=true;
    	  }
    	  ++iz;
    	}
     	while(iz<id_out.size()&&!found);
     	if(found) ++nfound;
	if(nfound==id_out.size()) {
	  cc.push_back(false);
	  loc.push_back(ix);
        }
      }
    }
    // the charge conjugate mode
    if(incoming[ix]->id()==idbar_in&&outgoing[ix].size()==idbar_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound = 0;
      for(unsigned int iy=0;iy<idbar_out.size();++iy) {
    	int idtemp=outgoing[ix][iy]->id();
	unsigned int iz(0);
   	bool found=false;
   	do {
     	  if(idtemp==idbar_out[iz]&&!done[iz]) {
     	    done[iz]=true;
     	    found=true;
     	  }
     	  ++iz;
     	}
     	while(iz<idbar_out.size()&&!found);
    	if(found) ++nfound;
      }
      if(nfound==idbar_out.size()) {
	cc.push_back(false);
	loc.push_back(ix);
      }
    }
  }
}

void ScalarMesonFactorizedDecayer::dataBaseOutput(ofstream & output,
						  bool header) const {
  unsigned int ix;
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":a1Bottom "  << _a1b << "\n";
  output << "newdef " << name() << ":a2Bottom "  << _a2b << "\n";
  output << "newdef " << name() << ":a1Charm "   << _a1c << "\n";
  output << "newdef " << name() << ":a2Charm "   << _a2c << "\n";
  for(ix=0;ix<_current.size();++ix) {
    _current[ix]->dataBaseOutput(output,false,true);
    output << "insert " << name() << ":Currents " << ix << " " 
	   << _current[ix]->name() << " \n";
  }
  for(ix=0;ix<_form.size();++ix) {
    _form[ix]->dataBaseOutput(output,false,true);
    output << "insert " << name() << ":FormFactors " << ix << " " 
	   << _form[ix]->name() << " \n";
  }
  for(ix=0;ix<_wgtloc.size();++ix) {
    output << "insert " << name() << ":WeightLocation " << ix << " " 
	   << _wgtloc[ix] << "\n";
  }
  for(ix=0;ix<_wgtmax.size();++ix) {
    output << "insert " << name() << ":MaximumWeight "  << ix << " " 
	   << _wgtmax[ix] << "\n";
  }
  for(ix=0;ix<_weights.size();++ix) {
    output << "insert " << name() << ":Weights "        << ix << " " 
	   << _weights[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./ScalarVectorVectorDecayer.cc"
// -*- C++ -*-
//
// ScalarVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarVectorVectorDecayer class.
//

#include "ScalarVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void ScalarVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void ScalarVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "ScalarVectorVectorDecayerDecayer" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr in     =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first ),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode=new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int ScalarVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  cc = false;
  // check that at least some modes exist
  // must be two outgoing particles
  if(incoming_.size()==0||children.size()!=2) return -1;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  int ccid0(parent     ->CC() ? parent     ->CC()->id() : parent     ->id());
  int ccid1(children[0]->CC() ? children[0]->CC()->id() : children[0]->id());
  int ccid2(children[1]->CC() ? children[1]->CC()->id() : children[1]->id());
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  int imode(-1);
  do {
    if(incoming_[ix]==id0) {
      if((outgoing_[ix].first==id1&&outgoing_[ix].second==id2)||
	 (outgoing_[ix].first==id2&&outgoing_[ix].second==id1)) {
	imode=ix;
	break;
      }
    }
    if(incoming_[ix]==ccid0) {
      if((outgoing_[ix].first==ccid1&&outgoing_[ix].second==ccid2)||
	 (outgoing_[ix].first==ccid2&&outgoing_[ix].second==ccid1)) {
	imode=ix;
	cc=true;
	break;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

void ScalarVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,GeV) << incoming_ << outgoing_ << maxWeight_;
}

void ScalarVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,GeV) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarVectorVectorDecayer,DecayIntegrator>
describeHerwigScalarVectorVectorDecayer("Herwig::ScalarVectorVectorDecayer", "HwSMDecay.so");

void ScalarVectorVectorDecayer::Init() {

  static ClassDocumentation<ScalarVectorVectorDecayer> documentation
    ("The ScalarVectorVectorDecayer class is designed for"
     " the decay of a scalar meson to two spin-1 particles.");

  static Command<ScalarVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(1/GeV) and max weight for a decay",
     &ScalarVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<ScalarVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<ScalarVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in ScalarVectorVectorDecayer have been deleted, please use SetUpDecayMode");
}

void ScalarVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix],decay[ix],
					  outgoing,true,decay[ix]->id()==ParticleID::gamma);
}

double ScalarVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  bool photon[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],ihel,Helicity::outgoing);
    }
  }
  // now compute the matrix element
  double fact(coupling_[imode()]/part.mass());
  Energy2 p1p2(momenta[0]*momenta[1]);
  unsigned int ix,iy;
  for(ix=0;ix<3;++ix) {
    if(photon[0] && ix==1) continue;
    for(iy=0;iy<3;++iy) {
      if(photon[1] && iy==1) continue;
      (*ME())(0,ix,iy)=Complex(fact*(vectors_[0][ix].dot(vectors_[1][iy])-
				     (vectors_[1][iy]*momenta[0])*
				     (vectors_[0][ix]*momenta[1])/(p1p2-momenta[0].mass()*momenta[1].mass())));
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = sqr(coupling_[imode()]/part.mass())*
  //   (2.*sqr(pcm*part.mass())+3.*sqr(momenta[0].mass()*momenta[1].mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  return me;
}

// output the setup info for the particle database
void ScalarVectorVectorDecayer::dataBaseOutput(ofstream & output,
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

// specify the 1-2 matrix element to be used in the running width calculation
bool ScalarVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm, int & itype,
					      double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id());
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id());++pit;
  int id2((**pit).id());
  unsigned int ix(0);
  do {
    if(incoming_[ix]==id) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  coupling=coupling_[imode]/dm.parent()->mass();
  itype = 12;
  return id1==outgoing_[imode].first&&id2==outgoing_[imode].second;
}

string ScalarVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 0";
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
#line 1 "./Scalar2FermionsDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Scalar2FermionsDecayer class.
//

#include "Scalar2FermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void Scalar2FermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}


void Scalar2FermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "FermionDecayer::doiin() " << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

Scalar2FermionsDecayer::Scalar2FermionsDecayer() {
  // don't include intermediates
  generateIntermediates(false);
}

int Scalar2FermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(incoming_[ix]==id   ) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(incoming_[ix]==idbar) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

IBPtr Scalar2FermionsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr Scalar2FermionsDecayer::fullclone() const {
  return new_ptr(*this);
}

void Scalar2FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxweight_;
}

void Scalar2FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxweight_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<Scalar2FermionsDecayer,DecayIntegrator>
describeHerwigScalar2FermionsDecayer("Herwig::Scalar2FermionsDecayer", "HwSMDecay.so");

void Scalar2FermionsDecayer::Init() {

  static ClassDocumentation<Scalar2FermionsDecayer> documentation
    ("The Scalar2FermionsDecayer class implements the decay of a scalar meson "
     "to a fermion and antifermion.");
  
  static Command<Scalar2FermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, fermion, antifermion), coupling and max weight for a decay",
     &Scalar2FermionsDecayer::setUpDecayMode, false);

}

void Scalar2FermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  // set up the spin information for the decay products
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}


double Scalar2FermionsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // prefactor
  InvEnergy pre(coupling_[imode()]/part.mass());
  // now compute the ME
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      Complex temp = pre*wave_[ix].scalar(wavebar_[iy]);
      if(iferm>ianti) (*ME())(0,ix,iy)=temp;
      else            (*ME())(0,iy,ix)=temp;
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 2*sqr(coupling_[imode()])*(1.-sqr((momenta[0].mass()+momenta[1].mass())/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool Scalar2FermionsDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					   double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
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
  coupling=coupling_[imode];
  mecode=25;
  return order;
}

// output the setup information for the particle database
void Scalar2FermionsDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string Scalar2FermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
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
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./PseudoScalar2FermionsDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoScalar2FermionsDecayer class.
//

#include "PseudoScalar2FermionsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoScalar2FermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}


void PseudoScalar2FermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=coupling_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "FermionDecayer::doiin() " << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

PseudoScalar2FermionsDecayer::PseudoScalar2FermionsDecayer() {
  // don't include intermediates
  generateIntermediates(false);
}

int PseudoScalar2FermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(incoming_[ix]==id   ) {
      if((id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second)||
	 (id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second)) imode=ix;
    }
    if(incoming_[ix]==idbar) {
      if((id1bar==outgoing_[ix].first&&id2bar==outgoing_[ix].second)||
	 (id2bar==outgoing_[ix].first&&id1bar==outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  return imode;
}

IBPtr PseudoScalar2FermionsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr PseudoScalar2FermionsDecayer::fullclone() const {
  return new_ptr(*this);
}

void PseudoScalar2FermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxweight_;
}

void PseudoScalar2FermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxweight_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoScalar2FermionsDecayer,DecayIntegrator>
describeHerwigPseudoScalar2FermionsDecayer("Herwig::PseudoScalar2FermionsDecayer", "HwSMDecay.so");

void PseudoScalar2FermionsDecayer::Init() {

  static ClassDocumentation<PseudoScalar2FermionsDecayer> documentation
    ("The PseudoScalar2FermionsDecayer class implements the decay of a pseudoscalar meson "
     "to a fermion and antifermion.");
  
  static Command<PseudoScalar2FermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, fermion, antifermion), coupling and max weight for a decay",
     &PseudoScalar2FermionsDecayer::setUpDecayMode, false);

}

void PseudoScalar2FermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  // set up the spin information for the decay products
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}


double PseudoScalar2FermionsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    ScalarWaveFunction::calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // prefactor
  InvEnergy pre(coupling_[imode()]/part.mass());
  // now compute the ME
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      Complex temp = pre*wave_[ix].pseudoScalar(wavebar_[iy]);
      if(iferm>ianti) (*ME())(0,ix,iy)=temp;
      else            (*ME())(0,iy,ix)=temp;
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 2*sqr(coupling_[imode()])*(1.-sqr((momenta[0].mass()-momenta[1].mass())/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool PseudoScalar2FermionsDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
						double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()){idbar=dm.parent()->CC()->id();}
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()){id1bar=(**pit).CC()->id();}
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()){id2bar=(**pit).CC()->id();}
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
  coupling=coupling_[imode];
  mecode=24;
  return order;
}

// output the setup information for the particle database
void PseudoScalar2FermionsDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string PseudoScalar2FermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
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
  maxweight_.push_back(wgt);
  // success
  return "";
}
