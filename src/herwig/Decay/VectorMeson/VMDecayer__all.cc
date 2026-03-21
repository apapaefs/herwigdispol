#line 1 "./a1ThreePionCLEODecayer.cc"
// -*- C++ -*-
//
// a1ThreePionCLEODecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionCLEODecayer class.
//

#include "a1ThreePionCLEODecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
using Constants::pi;
 
a1ThreePionCLEODecayer::a1ThreePionCLEODecayer() 
  : _rhomass(2), _rhowidth(2), _f2mass(1.275*GeV), _f2width(0.185*GeV), 
    _pf2cc(ZERO), _pf200(ZERO), _f0mass(1.186*GeV), _f0width(0.350*GeV),
    _pf0cc(ZERO), _pf000(ZERO), _sigmamass(0.860*GeV), _sigmawidth(0.880*GeV), 
    _psigmacc(ZERO), _psigma00(ZERO), _mpi0(ZERO), _mpic(ZERO),
    _coupling(45.57/GeV), _rhomagP(2), _rhophaseP(2), _rhomagD(2),
    _rhophaseD(2), _f2mag(0.71/GeV2), _f2phase(0.56*pi), _f2coup(ZERO,ZERO),
    _f0mag(0.77), _f0phase(-0.54*pi), _f0coup(0.,0.), _sigmamag(2.10),
    _sigmaphase(0.23*pi), _sigmacoup(0.,0.), _localparameters(true),
    _zerowgts(9), _onewgts(9), _twowgts(9), _threewgts(12), _zeromax(13.0704),
    _onemax(6.91104), _twomax(6.94654), _threemax(6.40086) {
  // rho masses and widths
  _rhomass[0] = 0.7743*GeV; _rhowidth[0] = 0.1491*GeV;
  _rhomass[1] = 1.370 *GeV; _rhowidth[1] = 0.386 *GeV;
  // p-wave rho and rho prime
  _rhomagP[0] = 1.  ; _rhophaseP[0] = 0.;
  _rhomagP[1] = 0.12; _rhophaseP[1] = 0.99*pi;
  // d-wave rho and rho prime
  _rhomagD[0] = 0.37/GeV2; _rhophaseD[0] = -0.15*pi;
  _rhomagD[1] = 0.87/GeV2; _rhophaseD[1] =  0.53*pi;
  // set up the integration channels
  _zerowgts[0]  = 0.132162;_zerowgts[1]  = 0.116638;_zerowgts[2]  = 0.121088;
  _zerowgts[3]  = 0.10656 ;_zerowgts[4]  = 0.102577;_zerowgts[5]  = 0.101169;
  _zerowgts[6]  = 0.104587;_zerowgts[7]  = 0.104663;_zerowgts[8]  = 0.110557;
  _onewgts[0]   = 0.177017;_onewgts[1]   = 0.176011;_onewgts[2]   = 0.110129;
  _onewgts[3]   = 0.108023;_onewgts[4]   = 0.110553;_onewgts[5]   = 0.109976;
  _onewgts[6]   = 0.088634;_onewgts[7]   = 0.059104;_onewgts[8]   = 0.060553;
  _twowgts[0]   = 0.173357;_twowgts[1]   = 0.172283;_twowgts[2]   = 0.116031;
  _twowgts[3]   = 0.114642;_twowgts[4]   = 0.109058;_twowgts[5]   = 0.114073;
  _twowgts[6]   = 0.080946;_twowgts[7]   = 0.060135;_twowgts[8]   = 0.059477;
  _threewgts[0] = 0.125022;_threewgts[1] = 0.129911;_threewgts[2] = 0.074165;
  _threewgts[3] = 0.075813;_threewgts[4 ]= 0.071154;_threewgts[5 ]= 0.077730;
  _threewgts[6] = 0.082255;_threewgts[7 ]= 0.086761;_threewgts[8 ]= 0.067106;
  _threewgts[9] = 0.070171;_threewgts[10]= 0.070146;_threewgts[11]= 0.069767;
  // generation of intermediates
  generateIntermediates(true);
}
  
void a1ThreePionCLEODecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  // possible intermediate particles
  // the different rho resonances
  tPDPtr rhop[3] = {getParticleData(213),getParticleData(100213),
		    getParticleData(30213)};
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  tPDPtr rhom[3] = {getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // the f_2
  tPDPtr f2=getParticleData(225);
  // the f_0
  tPDPtr f0=getParticleData(10221);
  // set up the integration channels
  // decay mode a_10 -> pi0 pi0 pi0
  tPDPtr     in =  a10;
  tPDVector out = {pi0,pi0,pi0};
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,_zeromax));
  vector<PhaseSpaceChannel> channels;
  // there are six sigma channels
  for(unsigned int ix=0;ix<3;++ix) {
    tPDPtr temp;
    if(ix==0)      temp = sigma;
    else if(ix==1) temp = f2;
    else if(ix==2) temp = f0;
    if(temp) {
      channels.push_back((PhaseSpaceChannel(mode),0,temp,0,1,1,2,1,3));
      channels.push_back((PhaseSpaceChannel(mode),0,temp,0,2,1,1,1,3));
      channels.push_back((PhaseSpaceChannel(mode),0,temp,0,3,1,1,1,2));
    }
  }
  if(_zerowgts.size()!=channels.size()) 
    _zerowgts=vector<double>(channels.size(),
			     1./channels.size());
  for(unsigned int ix=0;ix<channels.size();++ix) {
    channels[ix].weight(_zerowgts[ix]);
    mode->addChannel(channels[ix]);
  }
  addMode(mode);
  // decay mode a_1+ -> pi+ pi0 pi0
  in   = a1p;
  out = {pi0,pi0,pip};
  mode = new_ptr(PhaseSpaceMode(in,out,_onemax));
  channels.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho+ channel
    channels.push_back((PhaseSpaceChannel(mode),0,rhop[ix],0,1,1,2,1,3));
    // second rho+ channel
    channels.push_back((PhaseSpaceChannel(mode),0,rhop[ix],0,2,1,1,1,3));
  }
  // the sigma channel
  if(sigma) {
    channels.push_back((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
  }
  //  the f_2  channel
  if(f2) {
    channels.push_back((PhaseSpaceChannel(mode),0,f2,0,3,1,1,1,2));
  }
  // the f_0 channel
  if(f0) {
    channels.push_back((PhaseSpaceChannel(mode),0,f0,0,3,1,1,1,2));
  }
  if(_onewgts.size()!=channels.size()) 
    _onewgts=vector<double>(channels.size(),1./channels.size());
  for(unsigned int ix=0;ix<channels.size();++ix) {
    channels[ix].weight(_onewgts[ix]);
    mode->addChannel(channels[ix]);
  }
  addMode(mode);
  // decay mode a_10 -> pi+ pi- pi0
  in = a10;
  out = {pip,pim,pi0};
  mode = new_ptr(PhaseSpaceMode(in,out,_twomax));
  channels.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho channel
    channels.push_back((PhaseSpaceChannel(mode),0,rhom[ix],0,1,1,2,1,3));
    // second channel
    channels.push_back((PhaseSpaceChannel(mode),0,rhom[ix],0,2,1,1,1,3));
  }
  // sigma channel
  if(sigma) {
    channels.push_back((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
  }
  // f_2 channel
  if(f2) {
    channels.push_back((PhaseSpaceChannel(mode),0,f2,0,3,1,1,1,2));
  }
  // f_0 channel
  if(f0) {
    channels.push_back((PhaseSpaceChannel(mode),0,f0,0,3,1,1,1,2));
  }
  if(_twowgts.size()!=channels.size()) 
    _twowgts=vector<double>(channels.size(),1./channels.size());
  for(unsigned int ix=0;ix<channels.size();++ix) {
    channels[ix].weight(_twowgts[ix]);
    mode->addChannel(channels[ix]);
  }
  addMode(mode);
  // decay mode a_1+ -> pi+ pi+ pi-
  in  = a1p;
  out = {pip,pip,pim};
  mode = new_ptr(PhaseSpaceMode(in,out,_threemax));
  channels.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    // the neutral rho channels
    if(!rho0[ix]) continue;
    // first channel
    channels.push_back((PhaseSpaceChannel(mode),0,rho0[ix],0,1,1,2,1,3));
    // interchanged channel
    channels.push_back((PhaseSpaceChannel(mode),0,rho0[ix],0,2,1,1,1,3));
  }
  // the sigma channels
  if(sigma) {
    channels.push_back((PhaseSpaceChannel(mode),0,sigma,0,1,1,2,1,3));
    channels.push_back((PhaseSpaceChannel(mode),0,sigma,0,2,1,1,1,3));
  }
  // the f_2 channels
  if(f2) {
    channels.push_back((PhaseSpaceChannel(mode),0,f2,0,1,1,2,1,3));
    channels.push_back((PhaseSpaceChannel(mode),0,f2,0,2,1,1,1,3));
  }
  // the f_0 channel
  if(f0) {
    channels.push_back((PhaseSpaceChannel(mode),0,f0,0,1,1,2,1,3));
    channels.push_back((PhaseSpaceChannel(mode),0,f0,0,2,1,1,1,3));
  }
  if(_threewgts.size()!=channels.size()) 
    _threewgts=vector<double>(channels.size(),1./channels.size());
  for(unsigned int ix=0;ix<channels.size();++ix) {
    channels[ix].weight(_threewgts[ix]);
    mode->addChannel(channels[ix]);
  }
  addMode(mode);
  // if using local parameters set the values in the phase space channels
  if(_localparameters) {
    for(unsigned int iy=0;iy<_rhomass.size();++iy) {
      resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
      resetIntermediate(rhop[iy],_rhomass[iy],_rhowidth[iy]);
      resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
    }
    resetIntermediate(sigma,_sigmamass,_sigmawidth);
    resetIntermediate(f2,_f2mass,_f2width);
    resetIntermediate(f0,_f0mass,_f0width);
    // make sure the rho array has enough masses
    if(_rhomass.size()<3) {
      for(unsigned int ix=_rhomass.size();ix<3;++ix) {
	_rhomass.push_back(rhop[ix]->mass());
	_rhowidth.push_back(rhop[ix]->width());
      }
    }
  }
  // set the local variables if needed
  else {
    // masses and widths for the particles
    _rhomass.resize(3);_rhowidth.resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      _rhomass[ix]=rhop[ix]->mass();
      _rhowidth[ix]=rhop[ix]->width();
    }
    if(f2) {
      _f2mass=f2->mass();
      _f2width=f2->width();
    }
    if(f0) {
      _f0mass=f0->mass();
      _f0width=f0->width();
    }
    if(sigma) {
      _sigmamass=sigma->mass();
      _sigmawidth=sigma->width();
    }
  }
  // parameters for the breit-wigners
  _mpic=pip->mass();
  _mpi0=pi0->mass();
  // momenta of the decay products for on-shell particles
  _psigmacc = Kinematics::pstarTwoBodyDecay(_sigmamass,_mpic,_mpic);
  _psigma00 = Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi0,_mpi0);
  _pf2cc    = Kinematics::pstarTwoBodyDecay(_f2mass   ,_mpic,_mpic);
  _pf200    = Kinematics::pstarTwoBodyDecay(_f2mass   ,_mpi0,_mpi0);
  _pf0cc    = Kinematics::pstarTwoBodyDecay(_f0mass   ,_mpic,_mpic);
  _pf000    = Kinematics::pstarTwoBodyDecay(_f0mass   ,_mpi0,_mpi0); 
  _prhocc.resize(3);_prhoc0.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    _prhocc[ix] = Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpic,_mpic);
    _prhoc0[ix] = Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpic,_mpi0);
  }
  // couplings for the different modes
  Complex ii(0.,1.);
  _rhocoupP.resize(_rhomagP.size());
  for(unsigned int ix=0;ix<_rhomagP.size();++ix)
    _rhocoupP[ix]=_rhomagP[ix]*(cos(_rhophaseP[ix])+ii*sin(_rhophaseP[ix]));
  _rhocoupD.resize(_rhomagD.size());
  for(unsigned int ix=0;ix<_rhomagD.size();++ix)
    _rhocoupD[ix]=_rhomagD[ix]*(cos(_rhophaseD[ix])+ii*sin(_rhophaseD[ix]));
  _f0coup=_f0mag*(cos(_f0phase)+ii*sin(_f0phase));
  _f2coup=_f2mag*(cos(_f2phase)+ii*sin(_f2phase));
  _sigmacoup=_sigmamag*(cos(_sigmaphase)+ii*sin(_sigmaphase));
}


inline void a1ThreePionCLEODecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    // get the weights for the different channels
    for(unsigned int ix=0;ix<_zerowgts.size();++ix)
      _zerowgts[ix]=mode(0)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_onewgts.size();++ix)
      _onewgts[ix]=mode(1)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_twowgts.size();++ix)
      _twowgts[ix]=mode(2)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_threewgts.size();++ix)
      _threewgts[ix]=mode(3)->channels()[ix].weight();
    // get the maximum weight
    _zeromax  = mode(0)->maxWeight();
    _onemax   = mode(1)->maxWeight();
    _twomax   = mode(2)->maxWeight();
    _threemax = mode(3)->maxWeight();
  }
}

int a1ThreePionCLEODecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  if(children.size()!=3) return -1;
  int id(parent->id());
  // check the pions
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,npi0(0),npiplus(0),npiminus(0);
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
  }
  int imode(-1);
  // a_1+ decay modes
  if(id==ParticleID::a_1plus) {
    cc=false;
    if(npiplus==1&&npi0==2)          imode=1;
    else if(npiplus==2&&npiminus==1) imode=3;
    }
  // a_1- modes
  else if(id==ParticleID::a_1minus) {
    cc=true;
    if(npiminus==1&&npi0==2)         imode=1;
    else if(npiminus==2&&npiplus==1) imode=3;
  }
  // a_0 modes
  else if(id==ParticleID::a_10) {
    cc=false;
    if(npiminus==1&&npiplus==1&&npi0==1) imode=2;
    else if(npi0==3)                     imode=0;
  }
  return imode;
}

void a1ThreePionCLEODecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) 
     << ounit(_prhocc,GeV) << ounit(_prhoc0,GeV) 
     << ounit(_f2mass,GeV) << ounit(_f2width,GeV) 
     << ounit(_pf2cc,GeV) 
     << ounit(_pf200,GeV) << ounit(_f0mass,GeV) 
     << ounit(_f0width,GeV) << ounit(_pf0cc,GeV) << ounit(_pf000,GeV) 
     << ounit(_sigmamass,GeV) << ounit(_sigmawidth,GeV)
     << ounit(_psigmacc,GeV) << ounit(_psigma00,GeV) 
     << ounit(_mpi0,GeV) << ounit(_mpic,GeV) 
     << ounit(_coupling,1/GeV) << _rhomagP << _rhophaseP 
     << _rhocoupP << ounit(_rhomagD,1/GeV2)<< _rhophaseD 
     << ounit(_rhocoupD,1/GeV2) << ounit(_f2mag,1/GeV2)
     << _f2phase << ounit(_f2coup,1/GeV2)
     << _f0mag << _f0phase << _f0coup 
     << _sigmamag << _sigmaphase << _sigmacoup 
     << _localparameters 
     << _zerowgts << _onewgts << _twowgts << _threewgts 
     << _zeromax << _onemax << _twomax << _threemax;
}
  
void a1ThreePionCLEODecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) 
     >> iunit(_prhocc,GeV) >> iunit(_prhoc0,GeV) 
     >> iunit(_f2mass,GeV) >> iunit(_f2width,GeV) 
     >> iunit(_pf2cc,GeV)
     >> iunit(_pf200,GeV) >> iunit(_f0mass,GeV) 
     >> iunit(_f0width,GeV) >> iunit(_pf0cc,GeV) >> iunit(_pf000,GeV) 
     >> iunit(_sigmamass,GeV) >> iunit(_sigmawidth,GeV) 
     >> iunit(_psigmacc,GeV) >> iunit(_psigma00,GeV) 
     >> iunit(_mpi0,GeV) >> iunit(_mpic,GeV) 
     >> iunit(_coupling,1/GeV) >> _rhomagP >> _rhophaseP
     >> _rhocoupP >> iunit(_rhomagD,1/GeV2) >> _rhophaseD 
     >> iunit(_rhocoupD,1/GeV2)>>iunit(_f2mag,1/GeV2) 
     >> _f2phase >> iunit(_f2coup,1/GeV2)
     >> _f0mag >> _f0phase >> _f0coup
     >> _sigmamag >> _sigmaphase >> _sigmacoup
     >> _localparameters 
     >> _zerowgts >> _onewgts >> _twowgts >> _threewgts
     >> _zeromax >> _onemax >> _twomax >> _threemax;
}
  
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<a1ThreePionCLEODecayer,DecayIntegrator>
describeHerwiga1ThreePionCLEODecayer("Herwig::a1ThreePionCLEODecayer", "HwVMDecay.so");
  
void a1ThreePionCLEODecayer::Init() {
  
  static ClassDocumentation<a1ThreePionCLEODecayer> documentation
    ("The a1ThreePionCLEODecayer class performs the decay of the "
     "a_1 to three pions using the model of CLEO",
     "The decay of a_1 to three pions was modelled after \\cite{Asner:1999kj}.",
     "%\\cite{Asner:1999kj}\n"
     "\\bibitem{Asner:1999kj}\n"
     "  D.~M.~Asner {\\it et al.}  [CLEO Collaboration],\n"
     "   ``Hadronic structure in the decay tau- --> nu/tau pi- pi0 pi0 and the  sign\n"
     "  %of the tau neutrino helicity,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 61}, 012002 (2000)\n"
     "  [arXiv:hep-ex/9902022].\n"
     "  %%CITATION = PHRVA,D61,012002;%%\n"
     );

  static ParVector<a1ThreePionCLEODecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionCLEODecayer::_rhomass,
     GeV, 0, ZERO, -10000*GeV, 10000*GeV, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionCLEODecayer::_rhowidth,
     GeV, 0, ZERO, -10000*GeV, 10000*GeV, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_2Mass
    ("f_2Mass",
     "The mass of the f_2 meson",
     &a1ThreePionCLEODecayer::_f2mass, GeV, 1.275*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_2Width
    ("f_2Width",
     "The width of the f_2 meson",
     &a1ThreePionCLEODecayer::_f2width, GeV, 0.185*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_0Mass
    ("f_0Mass",
     "The mass of the f_0 meson",
     &a1ThreePionCLEODecayer::_f0mass, GeV, 1.186*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacef_0Width
    ("f_0Width",
     "The width of the f_0 meson",
     &a1ThreePionCLEODecayer::_f0width, GeV, 0.350*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacesigmaMass
    ("sigmaMass",
     "The mass of the sigma meson",
     &a1ThreePionCLEODecayer::_sigmamass, GeV, 0.860*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,Energy> interfacesigmaWidth
    ("sigmaWidth",
     "The width of the sigma meson",
     &a1ThreePionCLEODecayer::_sigmawidth, GeV, 0.880*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionCLEODecayer::_coupling, 1./GeV, 45.57/GeV, -0./GeV, 1000./GeV,
     false, false, false);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhomagP
    ("RhoPWaveMagnitude",
     "The magnitude of the couplings for the p-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagP,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseP
    ("RhoPWavePhase",
     "The phase of the couplings for the p-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseP,
     0, 0, 0, -pi, pi, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,InvEnergy2> interfacerhomagD
    ("RhoDWaveMagnitude",
     "The magnitude of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhomagD,
     1/MeV2, 0, ZERO, ZERO, 10000/MeV2, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacerhophaseD
    ("RhoDWavePhase",
     "The phase of the couplings for the d-wave rho currents",
     &a1ThreePionCLEODecayer::_rhophaseD,
     0, 0, 0, -pi, pi, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Phase
    ("f0Phase",
     "The phase of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0phase, 0.54*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef2Phase
    ("f2Phase",
     "The phase of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2phase, 0.56*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaPhase
    ("sigmaPhase",
     "The phase of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmaphase, 0.23*pi, -pi, pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionCLEODecayer,double> interfacef0Magnitude
    ("f0Magnitude",
     "The magnitude of the f_0 scalar current",
     &a1ThreePionCLEODecayer::_f0mag, 0.77, 0.0, 10,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,InvEnergy2> interfacef2Magnitude
    ("f2Magnitude",
     "The magnitude of the f_2 tensor current",
     &a1ThreePionCLEODecayer::_f2mag, 1./GeV2, 0.71/GeV2, ZERO, 10./GeV2,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfacesigmaMagnitude
    ("sigmaMagnitude",
     "The magnitude of the sigma scalar current",
     &a1ThreePionCLEODecayer::_sigmamag, 2.1, 0.0, 10,
     false, false, true);

  static Switch<a1ThreePionCLEODecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1ThreePionCLEODecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particleData objects",
     false);

  static ParVector<a1ThreePionCLEODecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionCLEODecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceZeroMax
    ("ZeroMax",
     "The maximum weight for the integration fo the channel a_1^0->pi0pi0pi0",
     &a1ThreePionCLEODecayer::_zeromax, 0.0716349E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1ThreePionCLEODecayer::_onemax,1.23756E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1ThreePionCLEODecayer::_twomax,2.43819E3, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionCLEODecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1ThreePionCLEODecayer::_threemax, 1.38754E3, 0.0, 10000.0,
     false, false, true);
}

void a1ThreePionCLEODecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double a1ThreePionCLEODecayer::me2(const int ichan, const Particle & part,
				   const tPDVector & ,
				   const vector<Lorentz5Momentum> & momenta,
				   MEOption meopt) const {
  useMe();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  // momentum of the incoming particle
  Lorentz5Momentum Q=part.momentum();
  Energy2 q2=Q.mass2();
  // identify the mesons
  // calculate the invariants and form factors
  Lorentz5Momentum ps1=momenta[1]+momenta[2];
  Lorentz5Momentum ps2=momenta[0]+momenta[2];
  Lorentz5Momentum ps3=momenta[0]+momenta[1];
  ps1.rescaleMass();
  ps2.rescaleMass();
  ps3.rescaleMass();
  Energy2 s1=ps1.mass2(),s2=ps2.mass2(),s3=ps3.mass2();
  complex<InvEnergy> F1,F2,F3;
  formFactors(imode(),ichan,q2,s1,s2,s3,F1,F2,F3);
  // use the form-factors to compute the current
  LorentzPolarizationVector output=
    LorentzPolarizationVector(F1*momenta[1])-LorentzPolarizationVector(F1*momenta[2])
    - LorentzPolarizationVector(F2*momenta[2])+LorentzPolarizationVector(F2*momenta[0])
    + LorentzPolarizationVector(F3*momenta[0])-LorentzPolarizationVector(F3*momenta[1]);
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) 
    (*ME())(ix,0,0,0)=output.dot(_vectors[ix]);
  // answer
  double out = ME()->contract(_rho).real();
  // test of the answer
//   double test = threeBodyMatrixElement(imode(),sqr(part.mass()),
// 				       s3,s2,s1,momenta[0].mass(),momenta[1].mass(), 
// 				       momenta[2].mass());
//   if(ichan<0) cerr << "testing matrix element " << part.PDGName() << " -> "
//        << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
//        << outgoing[2]->PDGName() << out << " " << test << " " 
//        << (out-test)/(out+test) << "\n";  
  // return the answer
  return out;
}

// matrix element for the running a_1 width
double a1ThreePionCLEODecayer::
threeBodyMatrixElement(const int iopt,const Energy2 q2, const Energy2 s3,
		       const Energy2 s2,const Energy2 s1, const Energy m1, 
		       const Energy m2 ,const Energy m3) const {
  Energy2 m12=m1*m1,m22=m2*m2,m32=m3*m3;
  // calculate the form factors
  complex<InvEnergy> F1,F2,F3;
  formFactors(iopt,-1,q2,s1,s2,s3,F1,F2,F3);
  // analytic calculation of the matrix element
  double dot1=( F1*conj(F1)*(2.*m22+2.*m32-s1)+F2*conj(F2)*(2.*m12+2.*m32-s2)
		+F3*conj(F3)*(2.*m12+2.*m22-s3)-F1*conj(F2)*( s1+s2-s3-4.*m32)
		+F1*conj(F3)*( s1-s2+s3-4.*m22)-F2*conj(F3)*(-s1+s2+s3-4.*m12)).real();
  complex<Energy> dot2 = 0.5*(F1*(s3-m32-s2+m22)-F2*(s1-m12-s3+m32)+F3*(s2-m22-s1+m12));
  return (-dot1+(dot2*conj(dot2)).real()/q2)/3.;
}

WidthCalculatorBasePtr 
a1ThreePionCLEODecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit) {
    if(abs((**pit).id())==ParticleID::piplus) ++ncharged;
  }
  // integrator to perform the integral
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  vector<Energy> inmass(2,_rhomass[0]),inwidth(2,_rhowidth[0]);
  vector<double> inpow(2,0.0);
  Energy m[3];
  m[0] = ncharged<2 ? _mpi0 : _mpic;
  m[1] = m[0];
  m[2] = (ncharged==0||ncharged==2) ? _mpi0 : _mpic;
  return new_ptr(ThreeBodyAllOnCalculator<a1ThreePionCLEODecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,ncharged,m[0],m[1],m[2]));
}

// calculate the form factos
void a1ThreePionCLEODecayer::formFactors(int iopt,int ichan,
					 Energy2 q2,Energy2 s1,Energy2 s2,
					 Energy2 s3,
					 complex<InvEnergy> & FF1,
					 complex<InvEnergy> & FF2,
					 complex<InvEnergy> & FF3) const {
  Complex F1, F2, F3;
  InvEnergy fact = _coupling;
  // a_1^0 pi0 pi0 pi0 mode
  if(iopt==0) {
    fact*=1./sqrt(6.);
    // compute the breit wigners we need
    Complex sigbws1 = sigmaBreitWigner(s1,1);
    Complex sigbws2 = sigmaBreitWigner(s2,1);
    Complex sigbws3 = sigmaBreitWigner(s3,1);
    Complex f0bws1  = f0BreitWigner(s1,1);
    Complex f0bws2  = f0BreitWigner(s2,1);
    Complex f0bws3  = f0BreitWigner(s3,1);
    Complex f2bws1  = f2BreitWigner(s1,1);
    Complex f2bws2  = f2BreitWigner(s2,1);
    Complex f2bws3  = f2BreitWigner(s3,1);
    if(ichan<0) {
      // the scalar terms
      F1=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	-2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2=2./3.*(_sigmacoup*sigbws3+_f0coup*f0bws3)
	-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	+2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
      complex<Energy2> Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
      complex<Energy2> Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
      F1 += Complex(_f2coup*( 0.5*(s3-s2)*f2bws1-Dfact2+Dfact3));
      F2 += Complex(_f2coup*( 0.5*(s3-s1)*f2bws2-Dfact1+Dfact3));
      F3 += Complex(_f2coup*(-0.5*(s1-s2)*f2bws3-Dfact1+Dfact2));
    }
    else if(ichan==0) {
      F2=-2./3.*_sigmacoup*sigbws1;
      F3=-2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==1) {
      F1=-2./3.*_sigmacoup*sigbws2;
      F3=+2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==2) {
      F1= 2./3.*_sigmacoup*sigbws3;
      F2= 2./3.*_sigmacoup*sigbws3;
    }
    else if(ichan==3) {
      complex<Energy2> Dfact1 = 1./18.*(4.*_mpi0*_mpi0-s1)*(q2+s1-_mpi0*_mpi0)/s1*f2bws1;
      F1+=Complex(_f2coup*0.5*(s3-s2)*f2bws1);
      F2-=Complex(_f2coup*Dfact1); 
      F3-=Complex(_f2coup*Dfact1);
    }
    else if(ichan==4) {
      complex<Energy2> Dfact2 = 1./18.*(4.*_mpi0*_mpi0-s2)*(q2+s2-_mpi0*_mpi0)/s2*f2bws2;
      F2+=Complex(_f2coup*0.5*(s3-s1)*f2bws2);
      F1-=Complex(_f2coup*Dfact2);
      F3+=Complex(_f2coup*Dfact2);
    }
    else if(ichan==5) {
      complex<Energy2> Dfact3 = 1./18.*(4.*_mpi0*_mpi0-s3)*(q2-_mpi0*_mpi0+s3)/s3*f2bws3;
      F3+=Complex(-_f2coup*0.5*(s1-s2)*f2bws3);
      F1+=Complex(_f2coup*Dfact3);
      F2+=Complex(_f2coup*Dfact3);
    }
    else if(ichan==6) {
      F2=-2./3.*_f0coup*f0bws1;
      F3=-2./3.*_f0coup*f0bws1;
    }
    else if(ichan==7) {
      F1=-2./3.*_f0coup*f0bws2;
      F3=+2./3.*_f0coup*f0bws2;
    }
    else if(ichan==8) {
      F1= 2./3.*_f0coup*f0bws3;
      F2= 2./3.*_f0coup*f0bws3;
    }
  }
  // a_1^+ -> pi0 pi0 pi+
  else if(iopt==1) {
    fact *= 1./sqrt(2.);
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
    }
    f0bw  = f0BreitWigner(s3,1);
    sigbw = sigmaBreitWigner(s3,1);
    f2bw  = f2BreitWigner(s3,1);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2+=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3+=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;
      F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()) F1+=_rhocoupP[ires]*rhos1bw[ires];
      Energy2 Dfact2=-1./3.*((s3-_mpic*_mpic)-(s2-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F2+=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()) F2+=_rhocoupP[ires]*rhos2bw[ires];
      Energy2 Dfact1=-1./3.*((s3-_mpic*_mpic)-(s1-_mpi0*_mpi0));
      if(ires<_rhocoupD.size()) {
	F1+=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3-=Complex(_rhocoupD[ires]*Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpic*_mpic+s3)*(4.*_mpi0*_mpi0-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  // a_1^0 ->pi+pi-pi0
  else if(iopt==2) {
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bw,sigbw,f2bw;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,1);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,1);
    }
    f0bw  =f0BreitWigner(s3,0);
    sigbw =sigmaBreitWigner(s3,0);
    f2bw  =f2BreitWigner(s3,0);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1+=_rhocoupP[ix]*rhos1bw[ix];
	F2+=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
      Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1+=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2+=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3+=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      Complex scalar=2./3.*(_sigmacoup*sigbw+_f0coup*f0bw);
      F1+=scalar;
      F2+=scalar;
      // the tensor terms
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      if(ires<_rhocoupP.size()) F1+=_rhocoupP[ires]*rhos1bw[ires];
      Energy2 Dfact2=-1./3.*(s3-_mpi0*_mpi0-s2+_mpic*_mpic);
      if(ires<_rhocoupD.size()) {
	F2+=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      if(ires<_rhocoupP.size()) F2+=_rhocoupP[ires]*rhos2bw[ires];
      Energy2 Dfact1=-1./3.*(s3-_mpi0*_mpi0-s1+_mpic*_mpic);
      if(ires<_rhocoupD.size()) {
	F1+=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3-=Complex(_rhocoupD[ires]*-Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F1+=2./3.*_sigmacoup*sigbw;
      F2+=2./3.*_sigmacoup*sigbw;
    }
    else if(ichan==7) {
      Complex Dfact3=1./18./s3*_f2coup*(q2-_mpi0*_mpi0+s3)*(4.*_mpic*_mpic-s3)*f2bw;
      F1+=Dfact3;
      F2+=Dfact3;
      F3-=Complex(0.5*_f2coup*(s1-s2)*f2bw);
    }
    else if(ichan==8) {
      F1+=2./3.*_f0coup*f0bw;
      F2+=2./3.*_f0coup*f0bw;
    }
  }
  // a_1^+ -> pi+ pi+ pi- mode
  else {
    fact *= 1./sqrt(2.);
    // compute the breit wigners we need
    Complex rhos1bw[3],rhos2bw[3],f0bws1,sigbws1,f2bws1,f0bws2,sigbws2,f2bws2;
    for(unsigned int ix=0,N=max(_rhocoupP.size(),_rhocoupD.size());ix<N;++ix) {
      rhos1bw[ix]=rhoBreitWigner(ix,s1,0);
      rhos2bw[ix]=rhoBreitWigner(ix,s2,0);
    }
    f0bws1  =f0BreitWigner(s1,0);
    sigbws1 =sigmaBreitWigner(s1,0);
    f2bws1  =f2BreitWigner(s1,0);
    f0bws2  =f0BreitWigner(s2,0);
    sigbws2 =sigmaBreitWigner(s2,0);
    f2bws2  =f2BreitWigner(s2,0);
    if(ichan<0) {
      // the p-wave rho terms
      for(unsigned int ix=0;ix<_rhocoupP.size();++ix) {
	F1-=_rhocoupP[ix]*rhos1bw[ix];
	F2-=_rhocoupP[ix]*rhos2bw[ix];
      }
      // the D-wave rho terms
      Energy2 Dfact1=1./3.*(s1-s3);
      Energy2 Dfact2=1./3.*(s2-s3);
      for(unsigned int ix=0;ix<_rhocoupD.size();++ix) {
	F1-=Complex(Dfact1*_rhocoupD[ix]*rhos2bw[ix]);
	F2-=Complex(Dfact2*_rhocoupD[ix]*rhos1bw[ix]);
	F3-=Complex(_rhocoupD[ix]*(Dfact2*rhos1bw[ix]-Dfact1*rhos2bw[ix]));
      }
      // the scalar terms
      F1-=2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      F2-=2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1);
      F3+=-2./3.*(_sigmacoup*sigbws1+_f0coup*f0bws1)
	+2./3.*(_sigmacoup*sigbws2+_f0coup*f0bws2);
      // the tensor terms
      complex<Energy2> sfact1 
	= 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      complex<Energy2> sfact2 
	= 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1+=Complex(_f2coup*(0.5*(s3-s2)*f2bws1-sfact2));
      F2+=Complex(_f2coup*(0.5*(s3-s1)*f2bws2-sfact1));
      F3+=Complex(_f2coup*(-sfact1+sfact2));
    }
    else if(ichan%2==0&&ichan<=4) {
      unsigned int ires=ichan/2;
      Energy2 Dfact2=1./3.*(s2-s3);
      if(ires<_rhocoupP.size()) F1-=_rhocoupP[ires]*rhos1bw[ires];
      if(ires<_rhocoupD.size()){
	F2-=Complex(Dfact2*_rhocoupD[ires]*rhos1bw[ires]);
	F3-=Complex(_rhocoupD[ires]*Dfact2*rhos1bw[ires]);
      }
    }
    else if(ichan%2==1&&ichan<=5) {
      unsigned int ires=(ichan-1)/2;
      Energy2 Dfact1=1./3.*(s1-s3);
      if(ires<_rhocoupP.size()) F2-=_rhocoupP[ires]*rhos2bw[ires];
      if(ires<_rhocoupD.size()) {
	F1-=Complex(Dfact1*_rhocoupD[ires]*rhos2bw[ires]);
	F3+=Complex(_rhocoupD[ires]*Dfact1*rhos2bw[ires]);
      }
    }
    else if(ichan==6) {
      F2-=2./3.*_sigmacoup*sigbws1;
      F3-=2./3.*_sigmacoup*sigbws1;
    }
    else if(ichan==7) {
      F1-=2./3.*_sigmacoup*sigbws2;
      F3+=2./3.*_sigmacoup*sigbws2;
    }
    else if(ichan==8) {
      complex<Energy2> sfact1 = 1./18.*(4.*_mpic*_mpic-s1)*(q2+s1-_mpic*_mpic)/s1*f2bws1;
      F1+=Complex(_f2coup*0.5*(s3-s2)*f2bws1);
      F2-=Complex(_f2coup*sfact1);
      F3-=Complex(_f2coup*sfact1);
    }
    else if(ichan==9) {
      complex<Energy2> sfact2 = 1./18.*(4.*_mpic*_mpic-s2)*(q2+s2-_mpic*_mpic)/s2*f2bws2;
      F1-=Complex(_f2coup*sfact2);
      F2+=Complex(_f2coup*0.5*(s3-s1)*f2bws2);
      F3+=Complex(_f2coup*sfact2);
    }
    else if(ichan==10) {
      F2-=2./3.*_f0coup*f0bws1;
      F3-=2./3.*_f0coup*f0bws1;
    }
    else if(ichan==11) {
      F1-=2./3.*_f0coup*f0bws2;
      F3+=2./3.*_f0coup*f0bws2;
    }
  }
  FF1 = F1 * fact;
  FF2 = F2 * fact;
  FF3 = F3 * fact;
} 

// output the setup information for the particle database
void a1ThreePionCLEODecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // masses and widths of the intermediate particles
  output << "newdef " << name() << ":f_2Mass "    << _f2mass/GeV     << "\n";
  output << "newdef " << name() << ":f_2Width "   << _f2width/GeV    << "\n";
  output << "newdef " << name() << ":f_0Mass "    << _f0mass/GeV     << "\n";
  output << "newdef " << name() << ":f_0Width "   << _f0width/GeV    << "\n";
  output << "newdef " << name() << ":sigmaMass "  << _sigmamass/GeV  << "\n";
  output << "newdef " << name() << ":sigmaWidth " << _sigmawidth/GeV << "\n";
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
  }
  // couplings and phases for different channels
  output << "newdef " << name() << ":f0Phase " << _f0phase << "\n";
  output << "newdef " << name() << ":f2Phase " << _f2phase<< "\n";
  output << "newdef " << name() << ":sigmaPhase " << _sigmaphase<< "\n";
  output << "newdef " << name() << ":f0Magnitude " << _f0mag<< "\n";
  output << "newdef " << name() << ":f2Magnitude " << _f2mag*GeV2 << "\n";
  output << "newdef " << name() << ":sigmaMagnitude " << _sigmamag << "\n";
  output << "newdef " << name() << ":Coupling " << _coupling*GeV << "\n";
  for(unsigned int ix=0;ix<_rhomagP.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoPWaveMagnitude " << ix << " " 
		    << _rhomagP[ix] << "\n";
    else     output << "insert " << name() << ":RhoPWaveMagnitude " << ix << " " 
		    << _rhomagP[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhophaseP.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoPWavePhase " << ix << " " 
		    << _rhophaseP[ix] << "\n";
    else     output << "insert " << name() << ":RhoPWavePhase " << ix << " " 
		    << _rhophaseP[ix] << "\n";
  }  
  for(unsigned int ix=0;ix<_rhomagD.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoDWaveMagnitude " << ix << " " 
		    << _rhomagD[ix]*MeV2 << "\n";
    else     output << "insert " << name() << ":RhoDWaveMagnitude " << ix << " " 
		    << _rhomagD[ix]*MeV2 << "\n";
  }
  for(unsigned int ix=0;ix<_rhophaseD.size();++ix) {
    if(ix<2) output << "newdef    " << name() << ":RhoDWavePhase " << ix << " " 
		    << _rhophaseD[ix] << "\n";
    else     output << "insert " << name() << ":RhoDWavePhase " << ix << " " 
		    << _rhophaseD[ix] << "\n";
  }
  // use local values of the masses etc.
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  // integration weights for the different channels
  for(unsigned int ix=0;ix<_zerowgts.size();++ix) {
    output << "newdef " << name() << ":AllNeutralWeights " 
	   << ix << " " << _zerowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_onewgts.size();++ix) {
    output << "newdef " << name() << ":OneChargedWeights " 
	   << ix << " " << _onewgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_twowgts.size();++ix) {
    output << "newdef " << name() << ":TwoChargedWeights " 
	   << ix << " " << _twowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_threewgts.size();++ix) {
    output << "newdef " << name() << ":ThreeChargedWeights " 
	   << ix << " " << _threewgts[ix] << "\n";
  }
  // maximum weights for the different  channels
  output << "newdef " << name() << ":ZeroMax "  << _zeromax  << "\n";
  output << "newdef " << name() << ":OneMax "   << _onemax   << "\n";
  output << "newdef " << name() << ":TwoMax "   << _twomax   << "\n";
  output << "newdef " << name() << ":ThreeMax " << _threemax << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./a1SimpleDecayer.cc"
// -*- C++ -*-
//
// a1SimpleDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1SimpleDecayer class.
//

#include "a1SimpleDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/WidthCalculatorBase.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

a1SimpleDecayer::a1SimpleDecayer() 
  // rho masses, widths and weights
  : _rhomass ({0.773*GeV,1.370*GeV,1.750*GeV}),
    _rhowidth({0.145*GeV,0.510*GeV,0.120*GeV}),
    _rhowgts({1.0,-0.145,0.}),_localparameters(true), 
    _coupling(47.95/GeV),
    // integration weights
    _onemax(5.4474), _twomax(5.47784), _threemax(5.40185),
    _onewgts  ({0.235562,0.231098,0.131071,0.131135,0.135841,0.135294}), 
    _twowgts  ({0.236208,0.229481,0.131169,0.133604,0.132685,0.136854}),
    _threewgts({0.234259,0.233634,0.135922,0.129231,0.133949,0.133005}),
    _mpi(ZERO) {
  generateIntermediates(true);
}

void a1SimpleDecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  // the different rho resonances
  tPDPtr rhop[3] = {getParticleData(213),getParticleData(100213),
		    getParticleData(30213)};
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  tPDPtr rhom[3] = {getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // decay mode a_1+ -> pi+ pi0 pi0
  tPDPtr in = a1p;
  tPDVector out = {pi0,pi0,pip};
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,_onemax));
  unsigned int nrho(0);
  for(unsigned int ix=0;ix<3;++ix) if(rhop[ix]) ++nrho;
  if(_onewgts.size()!=2*nrho) _onewgts=vector<double>(2*nrho,0.5/nrho);
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho+ channel
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rhop[ix],0,2,1,1,1,3));
    c1.weight(_onewgts[2*ix]);
    mode->addChannel(c1);
    // second rho+ channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rhop[ix],0,1,1,2,1,3));
    c2.weight(_onewgts[2*ix+1]);
    mode->addChannel(c2);
  }
  addMode(mode);
  // decay mode a_10 -> pi+ pi- pi0
  in = a10;
  out = {pip,pim,pi0};
  mode = new_ptr(PhaseSpaceMode(in,out,_twomax));
  if(_twowgts.size()!=2*nrho) _twowgts=vector<double>(2*nrho,0.5/nrho);
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho channel
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rhop[ix],0,2,1,1,1,3));
    c1.weight(_twowgts[2*ix]);
    mode->addChannel(c1);
    // second channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rhom[ix],0,1,1,2,1,3));
    c2.weight(_twowgts[2*ix+1]);
    mode->addChannel(c2);
  }
  addMode(mode);
  // decay mode a_1+ -> pi+ pi+ pi-
  in = a1p;
  out = {pip,pip,pim};
  mode = new_ptr(PhaseSpaceMode(in,out,_threemax));
  nrho = 0;
  for(unsigned int ix=0;ix<3;++ix) if(rho0[ix]) ++nrho;
  if(_threewgts.size()!=2*nrho) _threewgts=vector<double>(2*nrho,0.5/nrho);
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rho0[ix]) continue;
    // the neutral rho channels
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rho0[ix],0,2,1,1,1,3));
    c1.weight(_threewgts[2*ix]);
    mode->addChannel(c1);
    // interchanged channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rho0[ix],0,1,1,2,1,3));
    c2.weight(_threewgts[2*ix+1]);
    mode->addChannel(c2);
  }
  addMode(mode);
  // if using local parameters set the values in the phase space channels
  if(_localparameters) {
    for(unsigned int iy=0;iy<_rhomass.size();++iy) {
      resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
      resetIntermediate(rhop[iy],_rhomass[iy],_rhowidth[iy]);
      resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
    }
    // make sure the rho array has enough masses
    if(_rhomass.size()<3) {
      for(unsigned int ix=_rhomass.size();ix<3;++ix) {
	_rhomass.push_back(rhop[ix]->mass());
	_rhowidth.push_back(rhop[ix]->width());
      }
    }
  }
  // set the local variables if needed
  else {
    // masses and widths for the particles
    _rhomass.resize(3);_rhowidth.resize(3);
    for(unsigned int ix=0;ix<3;++ix) {
      if(!rhop[ix]) continue;
      _rhomass[ix]=rhop[ix]->mass();
      _rhowidth[ix]=rhop[ix]->width();
    }
  }
  _mpi = pip->mass();
}

void a1SimpleDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    // get the weights for the different channels
    for(unsigned int ix=0;ix<_onewgts.size();++ix)
      _onewgts[ix]=mode(0)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_twowgts.size();++ix)
      _twowgts[ix]=mode(1)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_threewgts.size();++ix)
      _threewgts[ix]=mode(2)->channels()[ix].weight();
    // get the maximum weight
    _onemax   = mode(0)->maxWeight();
    _twomax   = mode(1)->maxWeight();
    _threemax = mode(2)->maxWeight();
  }
}

void a1SimpleDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) << _rhowgts 
     << _localparameters << ounit(_coupling,1./GeV) << _onemax
     << _twomax << _threemax << _onewgts << _twowgts << _threewgts
     << ounit(_mpi,GeV);
}

void a1SimpleDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) >> _rhowgts 
     >> _localparameters >> iunit(_coupling,1./GeV) >> _onemax
     >> _twomax >> _threemax >> _onewgts >> _twowgts >> _threewgts
     >> iunit(_mpi,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<a1SimpleDecayer,DecayIntegrator>
describeHerwiga1SimpleDecayer("Herwig::a1SimpleDecayer", "HwVMDecay.so");

void a1SimpleDecayer::Init() {

  static ClassDocumentation<a1SimpleDecayer> documentation
    ("The a1SimpleDecayer class implements a simple model for the decay of"
     " the a_1 to three pions based on the approach of Kuhn and Santanmaria,"
     " Z.Phys. C48, 445 (1990)",
     "The decays of the $a_1$ were modelled using the approach of "
     "\\cite{Kuhn:1990ad}.\n",
     "\\bibitem{Kuhn:1990ad} J.~H.~Kuhn and A.~Santamaria,\n"
     "Z.\\ Phys.\\  C {\\bf 48} (1990) 445.\n"
     "%%CITATION = ZEPYA,C48,445;%%\n");

  static Switch<a1SimpleDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1SimpleDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particleData objects",
     false);

  static Parameter<a1SimpleDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1SimpleDecayer::_coupling, 1./GeV, 47.95/GeV, ZERO, 100./GeV,
     false, false, Interface::limited);

  static ParVector<a1SimpleDecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonaces",
     &a1SimpleDecayer::_rhomass, MeV, 3, 775.*MeV, ZERO, 10000*MeV,
     false, false, Interface::limited);

  static ParVector<a1SimpleDecayer,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the different rho resonances",
     &a1SimpleDecayer::_rhowidth, MeV, 3, 141*MeV, 0.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);
  
  static ParVector<a1SimpleDecayer,double> interfaceRhoWeights
    ("RhoWeights",
     "Weight for the different rho resonances",
     &a1SimpleDecayer::_rhowgts, 3, 0.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<a1SimpleDecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1SimpleDecayer::_onemax, 5.57613, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1SimpleDecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1SimpleDecayer::_twomax, 5.61218, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1SimpleDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1SimpleDecayer::_threemax, 5.5384, 0.0, 10000.0,
     false, false, true);
  
  static ParVector<a1SimpleDecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1SimpleDecayer::_onewgts,
     0, 0, 0, 0., 1., false, false, true);

  static ParVector<a1SimpleDecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1SimpleDecayer::_twowgts,
     0, 0, 0, 0., 1., false, false, true);

  static ParVector<a1SimpleDecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1SimpleDecayer::_threewgts,
     0, 0, 0, 0., 1., false, false, true);

}

int a1SimpleDecayer::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  if(children.size()!=3) return -1;
  int id(parent->id());
  // check the pions
  int npi0(0),npiplus(0),npiminus(0);
  for(auto child : children) {
    int idtemp= child->id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
  }
  int imode(-1);
  // a_1+ decay modes
  if(id==ParticleID::a_1plus) {
    cc=false;
    if(npiplus==1&&npi0==2)          imode=0;
    else if(npiplus==2&&npiminus==1) imode=2;
  }
  // a_1- modes
  else if(id==ParticleID::a_1minus) {
    cc=true;
    if(npiminus==1&&npi0==2)         imode=0;
    else if(npiminus==2&&npiplus==1) imode=2;
  }
  // a_0 modes
  else if(id==ParticleID::a_10) {
    cc=false;
    if(npiminus==1&&npiplus==1&&npi0==1) imode=1;
  }
  return imode;
}

void a1SimpleDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double a1SimpleDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  Lorentz5Vector<complex<Energy> > current;
  Energy2 s1 = (momenta[1]+momenta[2]).m2();
  Energy2 s2 = (momenta[0]+momenta[2]).m2();
  if(ichan<0) {
    current = rhoFormFactor(s2,-1)*(momenta[0]-momenta[2])
    +rhoFormFactor(s1,-1)*(momenta[1]-momenta[2]);
  }
  else if(ichan<3) {
    current = 
      rhoFormFactor(s2,ichan)*(momenta[0]-momenta[2]);
  }
  else if(ichan<6) {
    current = 
      rhoFormFactor(s1,-1)*(momenta[1]-momenta[2]);
  }
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix)
    (*ME())(ix,0,0,0) = Complex(_coupling*current.dot(_vectors[ix]));
  // matrix element and identical particle factor
  double output=ME()->contract(_rho).real();
  if(imode()!=1) output*=0.5;
  // test the output
  // Energy2 s3 = (momenta[0]+momenta[1]).m2();
  // double test = threeBodyMatrixElement(imode(),sqr(part.mass()),
  // 				       s3,s2,s1,momenta[0].mass(),
  // 				       momenta[1].mass(), 
  // 				       momenta[2].mass());
  // if(ichan<0) cerr << "testing matrix element " << part.PDGName() << " -> "
  // 		   << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  // 		   << outgoing[2]->PDGName() << output << " " << test << " " 
  // 		   << (output-test)/(output+test) << "\n";  
  // return the answer
  return output;
}

double a1SimpleDecayer::
threeBodyMatrixElement(const int iopt,const Energy2 q2, const Energy2 s3,
		       const Energy2 s2,const Energy2 s1, const Energy m1, 
		       const Energy m2 ,const Energy m3) const {
  Energy2 v12  = (s2-2.*sqr(m1)-2.*sqr(m3))+0.25*sqr(s1-s3-sqr(m1)+sqr(m3))/q2;
  Energy2 v22  = (s1-2.*sqr(m2)-2.*sqr(m3))+0.25*sqr(s2-s3-sqr(m2)+sqr(m3))/q2;
  Energy2 v1v2 = (0.5*q2-s3-0.5*(3*sqr(m3)-sqr(m1)-sqr(m2)))
    +0.25*(s1-s3-sqr(m1)+sqr(m3))*(s2-s3-sqr(m2)+sqr(m3))/q2;
  Complex rho1=rhoFormFactor(s2,-1);
  Complex rho2=rhoFormFactor(s1,-1);
  double me = sqr(_coupling)*real(v12*rho1*conj(rho1)+v22*rho2*conj(rho2)
				  +2.*v1v2*rho1*conj(rho2))/3.;
  if(iopt!=1) me *= 0.5;
  return me;
}

WidthCalculatorBasePtr
a1SimpleDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit) {
    if(abs((**pit).id())==ParticleID::piplus) ++ncharged;
  }
  --ncharged;
  // integrator to perform the integral
  vector<double> inweights;inweights.push_back(0.5);inweights.push_back(0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  vector<Energy> inmass(2,_rhomass[0]),inwidth(2,_rhowidth[0]);
  vector<double> inpow(2,0.0);
  Energy mpi0=getParticleData(ParticleID::pi0)->mass();
  Energy mpic=getParticleData(ParticleID::piplus)->mass();
  Energy m[3];
  m[0] = ncharged<2 ? mpi0 : mpic;
  m[1] = m[0];
  m[2] = (ncharged==0||ncharged==2) ? mpi0 : mpic;
  return new_ptr(ThreeBodyAllOnCalculator<a1SimpleDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,ncharged,m[0],m[1],m[2]));
}


void a1SimpleDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Coupling " << _coupling*GeV << "\n";
  output << "newdef " << name() << ":OneMax   " <<   _onemax << "\n";
  output << "newdef " << name() << ":TwoMax   " <<   _twomax << "\n";
  output << "newdef " << name() << ":ThreeMax " << _threemax << "\n";
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    output << "newdef    " << name() << ":RhoMasses " << ix << " "
	   << _rhomass[ix]/MeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    output << "newdef    " << name() << ":RhoWidths " << ix << " "
	   << _rhowidth[ix]/MeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowgts.size();++ix) {
    output << "newdef    " << name() << ":RhoWeights " << ix << " "
	   << _rhowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_onewgts.size();++ix) {
    output << "newdef " << name() << ":OneChargedWeights " 
	   << ix << " " << _onewgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_twowgts.size();++ix) {
    output << "newdef " << name() << ":TwoChargedWeights " 
	   << ix << " " << _twowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_threewgts.size();++ix) {
    output << "newdef " << name() << ":ThreeChargedWeights " 
	   << ix << " " << _threewgts[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./a1ThreePionDecayer.cc"
// -*- C++ -*-
//
// a1ThreePionDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the a1ThreePionDecayer class.
//

#include "a1ThreePionDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

a1ThreePionDecayer::a1ThreePionDecayer() 
  : _rhomass(1,0.7761*GeV), _rhowidth(1,0.1445*GeV), _sigmamass(0.8*GeV),
    _sigmawidth(0.8*GeV), _psigma(ZERO), _mpi(ZERO), _mpi2(ZERO),
    _lambda2(1.2*GeV2), _a1mass2(1.23*1.23*GeV2),
    _zsigma(0.), _zmag(1.3998721), _zphase(0.43585036),
    _rhomag(1,1.), _rhophase(1,0.), _coupling(90.44), 
    _localparameters(true),
    _zerowgts ({0.339108,0.335601,0.325291}),
    _onewgts  ({0.19616 ,0.191408,0.12137 ,0.115498,0.12729 ,0.127183,0.12109 }),
    _twowgts  ({0.188163,0.192479,0.121658,0.12135 ,0.127298,0.124835,0.124217}),
    _threewgts({0.153071,0.165741,0.107509,0.10275 ,0.109738,0.11254 ,0.125344,0.123307}) ,
    _zeromax(19.144), _onemax(7.83592), 
    _twomax(6.64804), _threemax(6.66296) {
  // generation of intermediates
  generateIntermediates(true);
}

void a1ThreePionDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    // get the weights for the different channels
    for(unsigned int ix=0;ix<_zerowgts.size();++ix)
      _zerowgts[ix]=mode(0)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_onewgts.size();++ix)
      _onewgts[ix]=mode(1)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_twowgts.size();++ix)
      _twowgts[ix]=mode(2)->channels()[ix].weight();
    for(unsigned int ix=0;ix<_threewgts.size();++ix)
      _threewgts[ix]=mode(3)->channels()[ix].weight();
    // get the maximum weight
    _zeromax  = mode(0)->maxWeight();
    _onemax   = mode(1)->maxWeight();
    _twomax   = mode(2)->maxWeight();
    _threemax = mode(3)->maxWeight();
  }
}

void a1ThreePionDecayer::doinit() {
  DecayIntegrator::doinit();
  // particles we need for the external state
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  // possible intermediate particles
  // the different rho resonances
  tPDPtr rhop[3] = {getParticleData(213),getParticleData(100213),
		    getParticleData(30213)};
  tPDPtr rho0[3] = {getParticleData(113),getParticleData(100113),
		    getParticleData(30113)};
  tPDPtr rhom[3] = {getParticleData(-213),getParticleData(-100213),
		    getParticleData(-30213)};
  // the sigma
  tPDPtr sigma = getParticleData(9000221);
  // set up the phase space integration
  // decay mode a_0 -> pi0 pi0 pi0
  tPDPtr in = a10;
  tPDVector out={pi0,pi0,pi0};
  if(sigma) {
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,_zeromax));
    if(_zerowgts.size()!=3) _zerowgts=vector<double>(3,1./3.);
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,sigma,0,1,1,2,1,3));
    c1.weight(_zerowgts[0]);
    mode->addChannel(c1);
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,sigma,0,2,1,1,1,3));
    c2.weight(_zerowgts[1]);
    mode->addChannel(c2);
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
    c3.weight(_zerowgts[2]);
    mode->addChannel(c3);
    addMode(mode);
  }
  // decay mode a_1+ -> pi+ pi0 pi0
  in = a1p;
  out = {pi0,pi0,pip};
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,_onemax));
  unsigned int nrho(0);
  for(unsigned int ix=0;ix<3;++ix) if(rhop[ix]) ++nrho;
  if(_onewgts.size()!=2*nrho+1) _onewgts=vector<double>(2*nrho+1,1./(2.*nrho+1.));
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho+ channel
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rhop[ix],0,1,1,2,1,3));
    c1.weight(_onewgts[2*ix]);
    mode->addChannel(c1);
    // second rho+ channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rhop[ix],0,2,1,1,1,3));
    c2.weight(_onewgts[2*ix+1]);
    mode->addChannel(c2);
  }
  // the sigma channel
  if(sigma) {
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
    c3.weight(_onewgts.back());
    mode->addChannel(c3);
  }
  addMode(mode);
  // decay mode a_1 -> pi+ pi- pi0
  in = a10;
  out = {pip,pim,pi0};
  mode = new_ptr(PhaseSpaceMode(in,out,_twomax));
  if(_twowgts.size()!=2*nrho+1) _twowgts=vector<double>(2*nrho+1,1./(2.*nrho+1.));
  for(unsigned int ix=0;ix<3;++ix) {
    if(!rhop[ix]) continue;
    // first rho channel
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rhop[ix],0,1,1,2,1,3));
    c1.weight(_twowgts[2*ix]);
    mode->addChannel(c1);
    // second channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rhom[ix],0,2,1,1,1,3));
    c2.weight(_twowgts[2*ix+1]);
    mode->addChannel(c2);
  }
  // sigma channel
  if(sigma) {
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,3,1,1,1,2));
    c3.weight(_twowgts.back());
    mode->addChannel(c3);
  }
  addMode(mode);
  // decay mode a_1+ -> pi+ pi+ pi-
  in = a1p;
  out = {pip,pip,pim};
  mode = new_ptr(PhaseSpaceMode(in,out,_threemax));
  nrho = 0;
  for(unsigned int ix=0;ix<3;++ix) if(rho0[ix]) ++nrho;
  if(_threewgts.size()!=2*nrho+2) _threewgts=vector<double>(2*nrho+2,1./(2.*nrho+2.));
  for(unsigned int ix=0;ix<3;++ix) {
    // the neutral rho channels
    if(!rho0[ix]) continue;
    // the neutral rho channels
    PhaseSpaceChannel c1((PhaseSpaceChannel(mode),0,rho0[ix],0,1,1,2,1,3));
    c1.weight(_threewgts[2*ix]);
    mode->addChannel(c1);
    // interchanged channel
    PhaseSpaceChannel c2((PhaseSpaceChannel(mode),0,rho0[ix],0,2,1,1,1,3));
    c2.weight(_threewgts[2*ix+1]);
    mode->addChannel(c2);
  }
  // the sigma channels
  if(sigma) {
    PhaseSpaceChannel c3((PhaseSpaceChannel(mode),0,sigma,0,1,1,2,1,3));
    c3.weight(_threewgts[6]);
    mode->addChannel(c3);
    PhaseSpaceChannel c4((PhaseSpaceChannel(mode),0,sigma,0,2,1,1,1,3));
    c4.weight(_threewgts[7]);
    mode->addChannel(c4);
  }
  addMode(mode);
  // set up the parameters 
  _mpi=getParticleData(ParticleID::piplus)->mass();
  _mpi2=sqr(_mpi);
  if(_localparameters) {
    if(_rhomass.size()<_rhocoupling.size()) {
      unsigned int itemp=_rhomass.size();
      _rhomass.resize(_rhocoupling.size());
      _rhowidth.resize(_rhocoupling.size());
      for(unsigned int ix=itemp;ix<_rhocoupling.size();++ix) {
	_rhomass[ix]=rhop[ix]->mass();
	_rhowidth[ix]=rhop[ix]->width();
      }
      // reset the intermediates in the phase space integration if needed
      resetIntermediate(sigma,_sigmamass,_sigmawidth);
      for(unsigned int iy=0;iy<_rhocoupling.size();++iy) {
	resetIntermediate(rho0[iy],_rhomass[iy],_rhowidth[iy]);
	resetIntermediate(rhop[iy],_rhomass[iy],_rhowidth[iy]);
	resetIntermediate(rhom[iy],_rhomass[iy],_rhowidth[iy]);
      }
    }
  }
  else {
    _a1mass2=sqr(getParticleData(ParticleID::a_1plus)->mass());
    if(sigma) {
      _sigmamass=sigma->mass();
      _sigmawidth=sigma->width();
    }
    _rhomass.resize(_rhocoupling.size());
    _rhowidth.resize(_rhocoupling.size());
    for(unsigned int ix=0;ix<_rhocoupling.size();++ix) {
      _rhomass[ix]=rhop[ix]->mass();
      _rhowidth[ix]=rhop[ix]->width();
    }
  }
  // parameters for the resonances
  // for the sigma
  _psigma=Kinematics::pstarTwoBodyDecay(_sigmamass,_mpi,_mpi);
  // for the rho
  _prho.resize(_rhomass.size());_hm2.resize(_rhomass.size());
  _dhdq2m2.resize(_rhomass.size());_rhoD.resize(_rhomass.size());
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    _prho[ix]    = Kinematics::pstarTwoBodyDecay(_rhomass[ix],_mpi,_mpi);
    _hm2[ix]     = hFunction(_rhomass[ix]);
    _dhdq2m2[ix] = dhdq2Parameter(ix);
    _rhoD[ix]    = DParameter(ix);
  }
  // convert the magnitude and phase of z into a phase
  _zsigma = _zmag*Complex(cos(_zphase),sin(_zphase));
  // convert rho couplings
  for(unsigned int ix=0;ix<_rhomag.size();++ix) {
    _rhocoupling.push_back(_rhomag[ix]*Complex(cos(_rhophase[ix]),sin(_rhophase[ix])));
  }
}
  
int a1ThreePionDecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  if(children.size()!=3) return -1;
  int id(parent->id());
  // check the pions
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,npi0(0),npiplus(0),npiminus(0);
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
  }
  int imode(-1);
  // a_1+ decay modes
  if(id==ParticleID::a_1plus) {
    cc=false;
    if(npiplus==1&&npi0==2)          imode=1;
    else if(npiplus==2&&npiminus==1) imode=3;
  }
  // a_1- modes
  else if(id==ParticleID::a_1minus) {
    cc=true;
    if(npiminus==1&&npi0==2)         imode=1;
    else if(npiminus==2&&npiplus==1) imode=3;
  }
  // a_0 modes
  else if(id==ParticleID::a_10) {
    cc=false;
    if(npiminus==1&&npiplus==1&&npi0==1) imode=2;
    else if(npi0==3)                     imode=0;
  }
  return imode;
}
  
void a1ThreePionDecayer::persistentOutput(PersistentOStream & os) const {
   os << ounit(_rhomass,GeV) << ounit(_rhowidth,GeV) << ounit(_prho,GeV) 
      << ounit(_hm2,GeV2) << ounit(_rhoD,GeV2) << _dhdq2m2 <<  ounit(_sigmamass,GeV)
      << ounit(_sigmawidth,GeV) << ounit(_psigma,GeV) << ounit(_mpi,GeV)
      << ounit(_mpi2,GeV2) << ounit(_lambda2,GeV2) << ounit(_a1mass2,GeV2) << _zsigma  
      << _rhocoupling << _coupling << _localparameters << _zerowgts << _onewgts 
      << _twowgts << _threewgts << _zeromax << _zmag << _zphase
      << _onemax << _twomax << _threemax << _coupling << _rhomag << _rhophase;
}
  
void a1ThreePionDecayer::persistentInput(PersistentIStream & is, int) {
   is >> iunit(_rhomass,GeV) >> iunit(_rhowidth,GeV) >> iunit(_prho,GeV) 
      >> iunit(_hm2,GeV2) >> iunit(_rhoD,GeV2) >> _dhdq2m2 >>  iunit(_sigmamass,GeV)
      >> iunit(_sigmawidth,GeV) >> iunit(_psigma,GeV) >> iunit(_mpi,GeV) 
      >> iunit(_mpi2,GeV2) >> iunit(_lambda2,GeV2) >> iunit(_a1mass2,GeV2) >> _zsigma
      >> _rhocoupling >> _coupling >> _localparameters >> _zerowgts >> _onewgts 
      >> _twowgts >> _threewgts >> _zeromax >> _zmag >> _zphase
      >> _onemax >> _twomax >> _threemax >> _coupling >> _rhomag >> _rhophase;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<a1ThreePionDecayer,DecayIntegrator>
describeHerwiga1ThreePionDecayer("Herwig::a1ThreePionDecayer", "HwVMDecay.so");
  
void a1ThreePionDecayer::Init() {
    
  static ClassDocumentation<a1ThreePionDecayer> documentation
    ("The a1ThreePionDecayer class is designed to decay the a_1 "
     "resonance to three pions using a model based on that used in the modelling "
     "of tau->4 pions.",
     "The decay of the $a_1$ resonance to three pions uses a model based on"
     "tau to four pions, \\cite{Bondar:2002mw}.",
     "%\\cite{Bondar:2002mw}\n"
     "\\bibitem{Bondar:2002mw}\n"
     "  A.~E.~Bondar, S.~I.~Eidelman, A.~I.~Milstein, T.~Pierzchala, N.~I.~Root, Z.~Was and M.~Worek,\n"
     "   ``Novosibirsk hadronic currents for tau --> 4pi channels of tau decay\n"
     "  %library TAUOLA,''\n"
     "  Comput.\\ Phys.\\ Commun.\\  {\\bf 146}, 139 (2002)\n"
     "  [arXiv:hep-ph/0201149].\n"
     "  %%CITATION = CPHCB,146,139;%%\n"
     );

  static Switch<a1ThreePionDecayer,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the intermediate resonances masses and widths",
     &a1ThreePionDecayer::_localparameters, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use the local values",
     true);
  static SwitchOption interfaceLocalParametersDefault
    (interfaceLocalParameters,
     "ParticleData",
     "Use the values from the particleData objects",
     false);

  static Parameter<a1ThreePionDecayer,double> interfaceCoupling
    ("Coupling",
     "The overall coupling for the decay",
     &a1ThreePionDecayer::_coupling, 90.44, 0.0, 1000.0,
     false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhomass
    ("RhoMasses",
     "The masses of the different rho resonnaces",
     &a1ThreePionDecayer::_rhomass,
     GeV, 0, ZERO, ZERO, 10000*GeV, false, false, true);

  static ParVector<a1ThreePionDecayer,Energy> interfacerhowidth
    ("RhoWidths",
     "The widths of the different rho resonnaces",
     &a1ThreePionDecayer::_rhowidth,
     GeV, 0, ZERO, ZERO, 10000*GeV, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfaceRhoMagnitude
    ("RhoMagnitude",
     "The magnitude of the rho couplings",
     &a1ThreePionDecayer::_rhomag, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static ParVector<a1ThreePionDecayer,double> interfaceRhoPhase
    ("RhoPhase",
     "The phase of the rho coupling",
     &a1ThreePionDecayer::_rhophase, -1, 0., 0.0, 2.*Constants::pi,
     false, false, Interface::limited);

  static Parameter<a1ThreePionDecayer,Energy2> interfaceLambda2
    ("Lambda2",
     "The value of the mass scale squared to use in the form-factor",
     &a1ThreePionDecayer::_lambda2, GeV2, 1.2*GeV2, 0.0001*GeV2, 10.0*GeV2,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy2> interfacea1mass2
    ("a1mass2",
     "The local value of the square of the a_1 mass",
     &a1ThreePionDecayer::_a1mass2, GeV2, 1.5129*GeV2, 0.5*GeV2, 10.0*GeV2,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaMass
    ("SigmaMass",
     "The local value of the sigma mass",
     &a1ThreePionDecayer::_sigmamass, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<a1ThreePionDecayer,Energy> interfaceSigmaWidth
    ("SigmaWidth",
     "The local value of the sigma width",
     &a1ThreePionDecayer::_sigmawidth, GeV, 0.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacezerowgts
    ("AllNeutralWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi0pi0pi0",
     &a1ThreePionDecayer::_zerowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfaceonewgts
    ("OneChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi0pi0",
     &a1ThreePionDecayer::_onewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacetwowgts
    ("TwoChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^0->pi+pi-pi0",
     &a1ThreePionDecayer::_twowgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<a1ThreePionDecayer,double> interfacethreewgts
    ("ThreeChargedWeights",
     "The weights of the different channels to use for the integration of"
     " the decay a_1^+->pi+pi+pi-",
     &a1ThreePionDecayer::_threewgts,
     0, 0, 0, -10000, 10000, false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceZeroMax
    ("ZeroMax",
     "The maximum weight for the integration fo the channel a_1^0->pi0pi0pi0",
     &a1ThreePionDecayer::_zeromax, 107.793, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceOneMax
    ("OneMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi0pi0",
     &a1ThreePionDecayer::_onemax, 1088.96, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceTwoMax
    ("TwoMax",
     "The maximum weight for the integration fo the channel a_1^0->pi+pi-pi0",
     &a1ThreePionDecayer::_twomax, 1750.73, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "The maximum weight for the integration fo the channel a_1^+->pi+pi+pi-",
     &a1ThreePionDecayer::_threemax, 739.334, 0.0, 10000.0,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceSigmaMagnitude
    ("SigmaMagnitude",
     "magnitude of the relative sigma coupling",
     &a1ThreePionDecayer::_zmag, 1.3998721, 0.0, 10.0e20,
     false, false, true);

  static Parameter<a1ThreePionDecayer,double> interfaceSigmaPhase
    ("SigmaPhase",
     "phase of the relative sigma coupling",
     &a1ThreePionDecayer::_zphase, 0.43585036, 0.0, Constants::twopi,
     false, false, true);
}

void   a1ThreePionDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double a1ThreePionDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  useMe();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  // momentum of the incoming particle
  Lorentz5Momentum Q=part.momentum();
  // momenta of the intermediates
  Energy2 s1=(momenta[1]+momenta[2]).m2();
  Energy2 s2=(momenta[0]+momenta[2]).m2();
  Energy2 s3=(momenta[0]+momenta[1]).m2();
  Energy2 dot01=Q*momenta[0];
  Energy2 dot02=Q*momenta[1];
  Energy2 dot03=Q*momenta[2];
  // vector for the output
  LorentzVector<complex<Energy3> > output;
  // a_10 -> pi0pi0pi0
  if(imode()==0) {
    //the breit-wigners
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    Complex sig3=sigmaBreitWigner(s3);
    // compute the vector
    LorentzPolarizationVectorE tmpoutput;
    if(ichan<0) {
      tmpoutput= sig1*(momenta[0])+sig2*(momenta[1])
	+sig3*(momenta[2]);
    }
    else if(ichan==0) tmpoutput=sig1*(momenta[0]);
    else if(ichan==1) tmpoutput=sig2*(momenta[1]);
    else if(ichan==2) tmpoutput=sig3*(momenta[2]);
    // the coupling z and identical particle factor
    output = tmpoutput * _zsigma* 1./sqrt(6.) *Q.mass2();
  }
  // a_1+ -> pi0pi0pi+
  else if(imode()==1) {
    // scalar propagator
    Complex sig1 = sigmaBreitWigner(s3);
    // sigma terms
    if(ichan<0||ichan==6) 
      output = _zsigma*Q.mass2()*sig1*momenta[2];
    // the rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output +=rho1*(dot03*(momenta[1])-
		       dot02*(momenta[2]));
      }
      if(ichan<0||ichan==2*ix+1){
	output +=rho2*(dot03*(momenta[0])-
		       dot01*(momenta[2]));
      }
    }
    // the identical particle factor
    output *= 1./sqrt(2.);
  }
  // a_10->pi+pi-pi0
  else if(imode()==2) {
    // the sigma terms
    Complex sig1=sigmaBreitWigner(s3);
    if(ichan<0||ichan==6)
      output = _zsigma*Q.mass2()*sig1*momenta[2];
    // rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1=_rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2=_rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output+=rho1*(dot03*(momenta[1])
		      -dot02*(momenta[2]));
      }
      if(ichan<0||ichan==2*ix+1) {
	output+=rho2*(dot03*(momenta[0])
		      -dot01*(momenta[2]));
      }
    }
  }
  // a1+ -> pi+pi+pi-
  else if(imode()==3) {
    // the scalar propagators 
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    // sigma terms
    LorentzPolarizationVectorE tmpoutput;
    if(ichan<0||ichan==6) tmpoutput+=sig1*(momenta[0]);
    if(ichan<0||ichan==7) tmpoutput+=sig2*(momenta[1]);
    output = tmpoutput * _zsigma * Q.mass2();
    // rho terms
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      Complex rho1 = _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      Complex rho2 = _rhocoupling[ix]*rhoBreitWigner(s2,ix);
      if(ichan<0||ichan==2*ix) {
	output-=rho1*( dot03*(momenta[1])-
		       dot02*(momenta[2]));
      }
      if(ichan<0||ichan==2*ix+1) {
	output-=rho2*( dot03*(momenta[0])-
		       dot01*(momenta[2]));
      }
    }
    // the identical particle factor
    output *= 1./sqrt(2.);
  }
  // form-factor
  LorentzPolarizationVector outputFinal 
    = output * a1FormFactor(Q.mass2())*_coupling/(Q.mass()*sqr(_rhomass[0]));
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix)
    (*ME())(ix,0,0,0)=outputFinal.dot(_vectors[ix]);
  // return the answer
  double out = ME()->contract(_rho).real();
  // test of the answer
  // double test = threeBodyMatrixElement(imode(),sqr(part.mass()),s3,s2,s1,
  // 				       momenta[0].mass(),momenta[1].mass(), 
  // 				       momenta[2].mass());
  // if(ichan<0) cerr << "testing matrix element " << part.PDGName() << " -> "
  // 		   << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  // 		   << outgoing[2]->PDGName() << " " << out << " " << test << " " 
  // 		   << (out-test)/(out+test) << "\n";
  return out;
}

// matrix element for the running a_1 width
double a1ThreePionDecayer::
threeBodyMatrixElement(const int iopt,const Energy2 q2, const Energy2 s3,
		       const Energy2 s2,const Energy2 s1, const Energy m1, 
		       const Energy m2 ,const Energy m3) const {
  Energy6 meout(0.*pow<3,1>(GeV2));
  Energy2 m12(sqr(m1)),m22(sqr(m2)),m32(sqr(m3));
  Energy2 dot01(q2-s1+m12),dot02(q2-s2+m22),dot03(q2-s3+m32),
    dot12(s3-m12-m22),dot13(s2-m12-m32),dot23(s1-m22-m32);
  if(iopt==0) {
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    Complex sig3=sigmaBreitWigner(s3);
    Energy2 metemp = 
      real(0.25*sig1*conj(sig1)*lambda(q2,s1,m12)/q2+
	   0.25*sig2*conj(sig2)*lambda(q2,s2,m22)/q2+
	   0.25*sig3*conj(sig3)*lambda(q2,s3,m32)/q2+
 	   sig1*conj(sig2)*(-dot12+0.5*dot01*dot02/q2)+
 	   sig1*conj(sig3)*(-dot13+0.5*dot01*dot03/q2)+
 	   sig2*conj(sig3)*(-dot23+0.5*dot02*dot03/q2));
    meout = metemp*real(_zsigma*conj(_zsigma))/6.*sqr(q2);
  }
  else if(iopt==1||iopt==2) {
    // the sigma terms
    Complex sig=sigmaBreitWigner(s3);
    Complex rho1,rho2;
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      rho1 += _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      rho2 += _rhocoupling[ix]*rhoBreitWigner(s2,ix);
    }
    meout =
      0.25*lambda(q2,m32,s3)*q2*norm(_zsigma*sig)+
      0.25*norm(rho1)*(dot23*dot02*dot03-m32*sqr(dot02)-m22*sqr(dot03))+
      0.25*norm(rho2)*(dot13*dot01*dot03-m32*sqr(dot01)-m12*sqr(dot03))-
      0.5*real(_zsigma*sig*conj(rho1))*q2*(dot03*dot23-2.*m32*dot02)-
      0.5*real(_zsigma*sig*conj(rho2))*q2*(dot03*dot13-2.*m32*dot01)-
      0.25*real(rho1*conj(rho2))*(sqr(dot03)*dot12-dot03*dot02*dot13
				  -dot03*dot01*dot23+2.*m32*dot02*dot01);
    if(iopt==1) meout *= 0.5;
  }
  else if(iopt==3) {
    Complex sig1=sigmaBreitWigner(s1);
    Complex sig2=sigmaBreitWigner(s2);
    Complex rho1,rho2;
    for(int ix=0,N=_rhocoupling.size();ix<N;++ix) {
      rho1 += _rhocoupling[ix]*rhoBreitWigner(s1,ix);
      rho2 += _rhocoupling[ix]*rhoBreitWigner(s2,ix);
    }
    meout =
      0.25*lambda(q2,m12,s1)*q2*norm(_zsigma*sig1)+
      0.25*lambda(q2,m22,s2)*q2*norm(_zsigma*sig2)+
      0.25*norm(rho1)*(dot23*dot02*dot03-m32*sqr(dot02)-m22*sqr(dot03))+
      0.25*norm(rho2)*(dot13*dot01*dot03-m32*sqr(dot01)-m12*sqr(dot03))-
      0.25*real(rho1*conj(rho2))*(sqr(dot03)*dot12-dot03*dot02*dot13
				  -dot03*dot01*dot23+2.*m32*dot02*dot01)-
      real(_zsigma*sig1*conj(_zsigma*sig2))*q2*(q2*dot12-0.5*dot02*dot01)+
      0.5*real(_zsigma*sig1*conj(rho1))*q2*(dot03*dot12-dot02*dot13)+
      0.5*real(_zsigma*sig2*conj(rho1))*q2*(2.*dot03*m22-dot02*dot23)+
      0.5*real(_zsigma*sig1*conj(rho2))*q2*(2.*dot03*m12-dot01*dot13)+
      0.5*real(_zsigma*sig2*conj(rho2))*q2*(dot03*dot12-dot01*dot23);
    meout *= 0.5;
  }
  return meout*a1FormFactor(q2)*sqr(_coupling/sqr(_rhomass[0]))/q2/3.;
}

WidthCalculatorBasePtr 
a1ThreePionDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  ParticleMSet::const_iterator pit  = dm.products().begin();
  ParticleMSet::const_iterator pend = dm.products().end();
  int ncharged=0;
  for( ; pit!=pend;++pit) {
    if(abs((**pit).id())==ParticleID::piplus) ++ncharged;
  }
  // integrator to perform the integral
  vector<double> inweights(2,0.5);
  vector<int> intype;intype.push_back(2);intype.push_back(3);
  vector<Energy> inmass(2,_rhomass[0]),inwidth(2,_rhowidth[0]);
  vector<double> inpow(2,0.0);
  Energy m[3];
  Energy mpi0=getParticleData(ParticleID::pi0)->mass();
  Energy mpic=getParticleData(ParticleID::piplus)->mass();
  m[0] = ncharged<2 ? mpi0 : mpic;
  m[1] = m[0];
  m[2] = (ncharged==0||ncharged==2) ? mpi0 : mpic;
  return new_ptr(ThreeBodyAllOnCalculator<a1ThreePionDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,ncharged,m[0],m[1],m[2]));
}

// output the setup information for the particle database
void a1ThreePionDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":LocalParameters " << _localparameters << "\n";
  output << "newdef " << name() << ":Coupling " << _coupling     << "\n";
  output << "newdef " << name() << ":Lambda2 "  << _lambda2/GeV2 << "\n";
  output << "newdef " << name() << ":a1mass2 "  << _a1mass2/GeV2 << "\n";
  output << "newdef " << name() << ":SigmaMass "  << _sigmamass/GeV  << "\n";
  output << "newdef " << name() << ":SigmaWidth " << _sigmawidth/GeV << "\n";
  output << "newdef " << name() << ":SigmaMagnitude " << _zmag << "\n";
  output << "newdef " << name() << ":SigmaPhase " << _zphase << "\n";
  for(unsigned int ix=0;ix<_rhomag.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoMagnitude " << ix << " " 
		    << _rhomag[ix] << "\n";
    else     output << "insert " << name() << ":RhoMagnitude " << ix << " " 
		    << _rhomag[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhophase.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoPhase " << ix << " " 
		    << _rhophase[ix] << "\n";
    else     output << "insert " << name() << ":RhoPhase " << ix << " " 
		    << _rhophase[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_rhomass.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoMasses " << ix << " " 
		    << _rhomass[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<_rhowidth.size();++ix) {
    if(ix<1) output << "newdef    " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
    else     output << "insert " << name() << ":RhoWidths " << ix << " " 
		    << _rhowidth[ix]/GeV << "\n";
  }
  // integration weights for the different channels
  for(unsigned int ix=0;ix<_zerowgts.size();++ix) {
    output << "newdef " << name() << ":AllNeutralWeights " 
	   << ix << " " << _zerowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_onewgts.size();++ix) {
    output << "newdef " << name() << ":OneChargedWeights " 
	   << ix << " " << _onewgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_twowgts.size();++ix) {
    output << "newdef " << name() << ":TwoChargedWeights " 
	   << ix << " " << _twowgts[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_threewgts.size();++ix) {
    output << "newdef " << name() << ":ThreeChargedWeights " 
	   << ix << " " << _threewgts[ix] << "\n";
  }
  output << "newdef " << name() << ":ZeroMax "  << _zeromax  << "\n";
  output << "newdef " << name() << ":OneMax "   << _onemax   << "\n";
  output << "newdef " << name() << ":TwoMax "   << _twomax   << "\n";
  output << "newdef " << name() << ":ThreeMax " << _threemax << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./VectorMeson2FermionDecayer.cc"
// -*- C++ -*-
//
// VectorMeson2FermionDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2FermionDecayer class.
//

#include "VectorMeson2FermionDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2FermionDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2FermionDecayer::doinit() {
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

VectorMeson2FermionDecayer::VectorMeson2FermionDecayer() {
  // don't include intermediates
  generateIntermediates(false);
}

int VectorMeson2FermionDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void VectorMeson2FermionDecayer::persistentOutput(PersistentOStream & os) const {
  os << coupling_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2FermionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> coupling_ >> incoming_ >> outgoing_ >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMeson2FermionDecayer,DecayIntegrator>
describeHerwigVectorMeson2FermionDecayer("Herwig::VectorMeson2FermionDecayer", "HwVMDecay.so");

void VectorMeson2FermionDecayer::Init() {

  static ClassDocumentation<VectorMeson2FermionDecayer> documentation
    ("The VectorMeson2FermionDecayer class is designed for the decay "
     "of vector mesons to fermions. It is mainly used for the decay of vector mesons "
     "to electrons and muons.");
  
  static Command<VectorMeson2FermionDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, fermion, antifermion), coupling and max weight for a decay",
     &VectorMeson2FermionDecayer::setUpDecayMode, false);

  static Deleted<VectorMeson2FermionDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceOutcoming1
    ("OutgoingFermion","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceOutcoming2
    ("OutgoingAntiFermion","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2FermionDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMeson2FermionDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMeson2FermionDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double VectorMeson2FermionDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // prefactor
  InvEnergy pre(coupling_[imode()]/part.mass());
  // now compute the currents
  LorentzPolarizationVector temp;
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      temp = pre*wave_[ix].vectorCurrent(wavebar_[iy]);
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // double test = 4.*sqr(coupling_[imode()])/3.*
  //   (1.+2.*sqr(momenta[0].mass()/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMeson2FermionDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  mecode=2;
  return order;
}

// output the setup information for the particle database
void VectorMeson2FermionDecayer::dataBaseOutput(ofstream & output,
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

string VectorMeson2FermionDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
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
#line 1 "./VectorMeson2SpinHalfBaryonsDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2SpinHalfBaryonsDecayer class.
//

#include "VectorMeson2SpinHalfBaryonsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2SpinHalfBaryonsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2SpinHalfBaryonsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=gm_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!= maxweight_.size()||
     isize!= phi_.size() || isize!= ge_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "BaryonsDecayer::doiin() " << Exception::runerror;
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

int VectorMeson2SpinHalfBaryonsDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

IBPtr VectorMeson2SpinHalfBaryonsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VectorMeson2SpinHalfBaryonsDecayer::fullclone() const {
  return new_ptr(*this);
}

void VectorMeson2SpinHalfBaryonsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ge_ << gm_ << phi_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2SpinHalfBaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> ge_ >> gm_ >> phi_ >> incoming_ >> outgoing_ >> maxweight_;
}


DescribeClass<VectorMeson2SpinHalfBaryonsDecayer,DecayIntegrator>
describeHerwigVectorMeson2SpinHalfBaryonsDecayer("Herwig::VectorMeson2SpinHalfBaryonsDecayer", "HwVMDecay.so");

void VectorMeson2SpinHalfBaryonsDecayer::Init() {

  static ClassDocumentation<VectorMeson2SpinHalfBaryonsDecayer> documentation
    ("The VectorMeson2SpinHalfBaryonsDecayer class is designed for "
     "the decay of vector mesons to baryons.");
  
  static Command<VectorMeson2SpinHalfBaryonsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, baryon, antibaryon), GE, GM, phase and max weight for a decay",
     &VectorMeson2SpinHalfBaryonsDecayer::setUpDecayMode, false);
}

void VectorMeson2SpinHalfBaryonsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // outgoing fermion
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  // outgoing antifermion
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double VectorMeson2SpinHalfBaryonsDecayer::me2(const int,const Particle & part,
					       const tPDVector & outgoing,
					       const vector<Lorentz5Momentum> & momenta,
					       MEOption meopt) const {
  // initialze me
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[iferm],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[ianti],ix,Helicity::outgoing);
  }
  // coefficients
  Complex GM = gm_[imode()];
  Complex GE = ge_[imode()]*exp(Complex(0.,phi_[imode()]));
  LorentzPolarizationVector c2 = -2.*outgoing[0]->mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()))*
    (GM-GE)*(momenta[iferm]-momenta[ianti]);
  // now compute the currents
  LorentzPolarizationVector temp;
  //double mesum(0.);
  for(unsigned ix=0;ix<2;++ix) {
    for(unsigned iy=0;iy<2;++iy) {
      LorentzPolarizationVector temp = (GM*wave_[ix].vectorCurrent(wavebar_[iy])+c2*wave_[ix].scalar(wavebar_[iy]))/part.mass();
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
	//mesum += norm(vectors_[iz].dot(temp));
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // cerr << "testing decay " << part.PDGName() << " -> " << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << "\n";
  // cerr << "testing ME " << mesum/3. << " " << me << " " << 4./3.*(norm(GM)+2.*sqr(outgoing[0]->mass()/part.mass())*norm(GE)) << "\n";
  // cerr << "testing gamma " << mesum/3./8./Constants::pi*sqrt(sqr(part.mass())-4.*sqr(outgoing[0]->mass()))/MeV << "\n";
  // return the answer
  return me;
}

// output the setup information for the particle database
void VectorMeson2SpinHalfBaryonsDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << gm_[ix] << " " << ge_[ix] << " " << phi_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMeson2SpinHalfBaryonsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
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
  // get the couplings
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ge_.push_back(stof(stype));
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  gm_.push_back(stof(stype));
  // and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  phi_.push_back(stof(stype));
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./VectorMeson2SpinThreeHalfBaryonsDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2SpinThreeHalfBaryonsDecayer class.
//

#include "VectorMeson2SpinThreeHalfBaryonsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMeson2SpinThreeHalfBaryonsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters arew consistent
  unsigned int isize=gm_.size();
  if(isize!=incoming_.size()  || isize!=outgoing_.size()  ||
     isize!= maxweight_.size()||
     isize!= phi_.size() || isize!= ge_.size())
    throw InitException() << "Inconsistent parameters in VectorMeson2"
			   << "BaryonsDecayer::doiin() " << Exception::runerror;
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

int VectorMeson2SpinThreeHalfBaryonsDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

IBPtr VectorMeson2SpinThreeHalfBaryonsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VectorMeson2SpinThreeHalfBaryonsDecayer::fullclone() const {
  return new_ptr(*this);
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ge_ << gm_ << phi_ << incoming_ << outgoing_ << maxweight_;
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> ge_ >> gm_ >> phi_ >> incoming_ >> outgoing_ >> maxweight_;
}


DescribeClass<VectorMeson2SpinThreeHalfBaryonsDecayer,DecayIntegrator>
describeHerwigVectorMeson2SpinThreeHalfBaryonsDecayer("Herwig::VectorMeson2SpinThreeHalfBaryonsDecayer", "HwVMDecay.so");

void VectorMeson2SpinThreeHalfBaryonsDecayer::Init() {

  static ClassDocumentation<VectorMeson2SpinThreeHalfBaryonsDecayer> documentation
    ("The VectorMeson2SpinThreeHalfBaryonsDecayer class is designed for "
     "the decay of vector mesons to baryons.");
  
  static Command<VectorMeson2SpinThreeHalfBaryonsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, baryon, antibaryon), GE, GM, phase and max weight for a decay",
     &VectorMeson2SpinThreeHalfBaryonsDecayer::setUpDecayMode, false);
}

void VectorMeson2SpinThreeHalfBaryonsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=decay[iferm]->id()) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // outgoing fermion
  RSSpinorBarWaveFunction::
    constructSpinInfo(wave2bar_,decay[iferm],outgoing,true);
  // outgoing antifermion
  RSSpinorWaveFunction::
    constructSpinInfo(wave2_,decay[ianti],outgoing,true);
}

double VectorMeson2SpinThreeHalfBaryonsDecayer::me2(const int,const Particle & part,
				       const tPDVector & outgoing,
				       const vector<Lorentz5Momentum> & momenta,
				       MEOption meopt) const {
  // initialze me
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin3Half,PDT::Spin3Half)));
  // fermion and antifermion
  unsigned int iferm(0),ianti(1);
  if(outgoing_[imode()].first!=outgoing[iferm]->id()) swap(iferm,ianti);
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave2_.resize(4);
  wave2bar_.resize(4);
  wave_.resize(4);
  wavebar_.resize(4);
  RSSpinorBarWaveFunction swave(momenta[iferm],outgoing[iferm],Helicity::outgoing);
  RSSpinorWaveFunction    awave(momenta[ianti],outgoing[ianti],Helicity::outgoing);
  LorentzPolarizationVector vtemp = part.momentum()/part.mass();
  for(unsigned int ix=0;ix<4;++ix) {
    swave.reset(ix);
    awave.reset(ix);
    wave2bar_[ix] = swave.dimensionedWf();
    wavebar_ [ix] = wave2bar_[ix].dot(vtemp);
    wave2_   [ix] = awave.dimensionedWf();
    wave_    [ix] = wave2_[ix].dot(vtemp);
  }
  // coefficients
  Complex GM = gm_[imode()];
  Complex GE = ge_[imode()]*exp(Complex(0.,phi_[imode()]));
  LorentzPolarizationVector c2 = -2.*outgoing[0]->mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()))*
    (GM-GE)*(momenta[iferm]-momenta[ianti]);
  // now compute the currents
  for(unsigned ix=0;ix<4;++ix) {
    for(unsigned iy=0;iy<4;++iy) {
      // q(al)q(be) piece
      LorentzPolarizationVector temp2 = (GM*wave_[ix].vectorCurrent(wavebar_[iy])+c2*wave_[ix].scalar(wavebar_[iy]))*
	2.*part.mass()/(4.*sqr(outgoing[0]->mass())-sqr(part.mass()));
      // g(al)g(be) * GM-GE piece
      LorentzPolarizationVector temp3 = wave2_[ix].generalScalar(wave2bar_[iy],1.,1.)*c2/part.mass();
      // g(al)g(be) * gamma^mu
      LorentzPolarizationVector temp1(GM/part.mass()*(wave2bar_[iy](0,3)*wave2_[ix](0,0) + wave2bar_[iy](0,2)*wave2_[ix](0,1) - 
						      wave2bar_[iy](0,1)*wave2_[ix](0,2) - wave2bar_[iy](0,0)*wave2_[ix](0,3) + 
						      wave2bar_[iy](1,3)*wave2_[ix](1,0) + wave2bar_[iy](1,2)*wave2_[ix](1,1) - 
						      wave2bar_[iy](1,1)*wave2_[ix](1,2) - wave2bar_[iy](1,0)*wave2_[ix](1,3) + 
						      wave2bar_[iy](2,3)*wave2_[ix](2,0) + wave2bar_[iy](2,2)*wave2_[ix](2,1) - 
						      wave2bar_[iy](2,1)*wave2_[ix](2,2) - wave2bar_[iy](2,0)*wave2_[ix](2,3) - 
						      wave2bar_[iy](3,3)*wave2_[ix](3,0) - wave2bar_[iy](3,2)*wave2_[ix](3,1) + 
						      wave2bar_[iy](3,1)*wave2_[ix](3,2) + wave2bar_[iy](3,0)*wave2_[ix](3,3)),
				      Complex(0,1)*GM/part.mass()*(wave2bar_[iy](0,3)*wave2_[ix](0,0) - wave2bar_[iy](0,2)*wave2_[ix](0,1) - 
								   wave2bar_[iy](0,1)*wave2_[ix](0,2) + wave2bar_[iy](0,0)*wave2_[ix](0,3) + 
								   wave2bar_[iy](1,3)*wave2_[ix](1,0) - wave2bar_[iy](1,2)*wave2_[ix](1,1) - 
								   wave2bar_[iy](1,1)*wave2_[ix](1,2) + wave2bar_[iy](1,0)*wave2_[ix](1,3) + 
								   wave2bar_[iy](2,3)*wave2_[ix](2,0) - wave2bar_[iy](2,2)*wave2_[ix](2,1) - 
								   wave2bar_[iy](2,1)*wave2_[ix](2,2) + wave2bar_[iy](2,0)*wave2_[ix](2,3) - 
								   wave2bar_[iy](3,3)*wave2_[ix](3,0) + wave2bar_[iy](3,2)*wave2_[ix](3,1) + 
								   wave2bar_[iy](3,1)*wave2_[ix](3,2) - wave2bar_[iy](3,0)*wave2_[ix](3,3)),
				      GM/part.mass()*(wave2bar_[iy](0,2)*wave2_[ix](0,0) - wave2bar_[iy](0,3)*wave2_[ix](0,1) - 
						      wave2bar_[iy](0,0)*wave2_[ix](0,2) + wave2bar_[iy](0,1)*wave2_[ix](0,3) + 
						      wave2bar_[iy](1,2)*wave2_[ix](1,0) - wave2bar_[iy](1,3)*wave2_[ix](1,1) - 
						      wave2bar_[iy](1,0)*wave2_[ix](1,2) + wave2bar_[iy](1,1)*wave2_[ix](1,3) + 
						      wave2bar_[iy](2,2)*wave2_[ix](2,0) - wave2bar_[iy](2,3)*wave2_[ix](2,1) - 
						      wave2bar_[iy](2,0)*wave2_[ix](2,2) + wave2bar_[iy](2,1)*wave2_[ix](2,3) - 
						      wave2bar_[iy](3,2)*wave2_[ix](3,0) + wave2bar_[iy](3,3)*wave2_[ix](3,1) + 
						      wave2bar_[iy](3,0)*wave2_[ix](3,2) - wave2bar_[iy](3,1)*wave2_[ix](3,3)),
				      GM/part.mass()*(-wave2bar_[iy](0,2)*wave2_[ix](0,0) - wave2bar_[iy](0,3)*wave2_[ix](0,1) - 
						      wave2bar_[iy](0,0)*wave2_[ix](0,2) - wave2bar_[iy](0,1)*wave2_[ix](0,3) - 
						      wave2bar_[iy](1,2)*wave2_[ix](1,0) - wave2bar_[iy](1,3)*wave2_[ix](1,1) - 
						      wave2bar_[iy](1,0)*wave2_[ix](1,2) - wave2bar_[iy](1,1)*wave2_[ix](1,3) - 
						      wave2bar_[iy](2,2)*wave2_[ix](2,0) - wave2bar_[iy](2,3)*wave2_[ix](2,1) - 
						      wave2bar_[iy](2,0)*wave2_[ix](2,2) - wave2bar_[iy](2,1)*wave2_[ix](2,3) + 
						      wave2bar_[iy](3,2)*wave2_[ix](3,0) + wave2bar_[iy](3,3)*wave2_[ix](3,1) + 
						      wave2bar_[iy](3,0)*wave2_[ix](3,2) + wave2bar_[iy](3,1)*wave2_[ix](3,3)));
      LorentzPolarizationVector temp = temp1+temp2+temp3;
      for(unsigned int iz=0;iz<3;++iz) {
	if(iferm>ianti) (*ME())(iz,ix,iy)=vectors_[iz].dot(temp);
	else            (*ME())(iz,iy,ix)=vectors_[iz].dot(temp);
      }
    }
  }
  // double mesum = ME()->contract(RhoDMatrix(PDT::Spin1)).real();
  // generator()->log() << "testing decay " << part.PDGName() << " -> " << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << "\n";
  // generator()->log() << "testing ME " << mesum  << " " << me << " " << 1./3.*(40.*norm(GM)/9.+16.*sqr(outgoing[0]->mass()/part.mass())*norm(GE)) << "\n";
  // return the answer
  return ME()->contract(rho_).real();
}

// output the setup information for the particle database
void VectorMeson2SpinThreeHalfBaryonsDecayer::dataBaseOutput(ofstream & output,
						bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << gm_[ix] << " " << ge_[ix] << " " << phi_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMeson2SpinThreeHalfBaryonsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin3Half)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 1/2";
  // get the couplings
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ge_.push_back(stof(stype));
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  gm_.push_back(stof(stype));
  // and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  phi_.push_back(stof(stype));
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./VectorMeson2MesonDecayer.cc"
// -*- C++ -*-
//
// VectorMeson2MesonDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMeson2MesonDecayer class.
//

#include "VectorMeson2MesonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig; 
using namespace ThePEG::Helicity;
 
void VectorMeson2MesonDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix]= mode(ix)->maxWeight();
  }
}

void VectorMeson2MesonDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() ||
     isize!=coupling_.size()) {
    throw InitException() << "Inconsistent parameters in "
			  << "VectorMeson2MesonDecayer" << Exception::runerror;
  }
  // set up the integration channelsx
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr     in =  getParticleData( incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

VectorMeson2MesonDecayer::VectorMeson2MesonDecayer() {
  // don't generate intermediates
  generateIntermediates(false);
}
 
int VectorMeson2MesonDecayer::modeNumber(bool & cc,tcPDPtr parent,
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
    if(idbar==incoming_[ix]&&imode<0) {
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
  
void VectorMeson2MesonDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_;
}
  
void VectorMeson2MesonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_;
}
  
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMeson2MesonDecayer,DecayIntegrator>
describeHerwigVectorMeson2MesonDecayer("Herwig::VectorMeson2MesonDecayer", "HwVMDecay.so");

void VectorMeson2MesonDecayer::Init() {
  
  static ClassDocumentation<VectorMeson2MesonDecayer> documentation
    ("The VectorMeson2MesonDecayer class is designed to implement "
     "the decay of vector mesons to 2 scalar mesons via a current which is the "
     "difference of the momenta of the two scalars. The order of the scalar meson "
     "momenta does not matter as it only changes the sign of the matrix element.");

  static Command<VectorMeson2MesonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, pseudoscalars, coupling and max weight for a decay",
     &VectorMeson2MesonDecayer::setUpDecayMode, false);

  static Deleted<VectorMeson2MesonDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMeson2MesonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMeson2MesonDecayer have been deleted, please use SetUpDecayMode");
  
}
 
void VectorMeson2MesonDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}
 
double VectorMeson2MesonDecayer::me2(const int,const Particle & part,
				     const tPDVector & ,
				     const vector<Lorentz5Momentum> & momenta,
				     MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  // difference of the momenta
  Lorentz5Vector<double> pdiff = (momenta[0]-momenta[1]) 
    * coupling_[imode()]/part.mass();
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) 
    (*ME())(ix,0,0)=vectors_[ix].dot(pdiff);
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = 4.*sqr(coupling_[imode()]*pcm/part.mass())/3.;
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}
 
bool VectorMeson2MesonDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
					     double & coupling) const {
  int imode(-1);
  int id(dm.parent()->id()),idbar(id);
  if(dm.parent()->CC()) idbar=dm.parent()->CC()->id();
  ParticleMSet::const_iterator pit(dm.products().begin());
  int id1((**pit).id()),id1bar(id1);
  if((**pit).CC()) id1bar=(**pit).CC()->id();
  ++pit;
  int id2((**pit).id()),id2bar(id2);
  if((**pit).CC()) id2bar=(**pit).CC()->id();
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
  mecode=0;
  return order;
}

// output the setup information for the particle database
void VectorMeson2MesonDecayer::dataBaseOutput(ofstream & output,
					      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string VectorMeson2MesonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./VectorMesonPScalarFermionsDecayer.cc"
// -*- C++ -*-
//
// VectorMesonPScalarFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonPScalarFermionsDecayer class.
//
//  Author: Peter Richardson
//

#include "VectorMesonPScalarFermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonPScalarFermionsDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      maxweight_[ix] = mode(ix)->maxWeight();
      weight_[ix]    = mode(ix)->channels()[1].weight();
    }
  }
}

VectorMesonPScalarFermionsDecayer::VectorMesonPScalarFermionsDecayer() {
  // intermediates
  generateIntermediates(false);
}

void VectorMesonPScalarFermionsDecayer::doinit() {
  DecayIntegrator::doinit();
  // check the parameters are consistent
  unsigned int isize(coupling_.size());
  if(isize!=incoming_.size()  || isize!=outgoing_.size() ||
     isize!=maxweight_.size()|| isize!=includeVMD_.size()||
     isize!=VMDid_.size()     || isize!=VMDmass_.size()  || isize!=VMDwidth_.size()||
     isize!=weight_.size())
    throw InitException() << "Inconsistent parameters in VectorMesonPScalar"
			  << "FermionsDecayer" << Exception::abortnow;
  // create the integration channel for each mode
  tPDPtr gamma(getParticleData(ParticleID::gamma)),rho;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    rho=getParticleData(VMDid_[ix]);
    tPDPtr in     =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second),
		     getParticleData(-outgoing_[ix].second)};
    PhaseSpaceModePtr newmode = new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    // photon channel
    PhaseSpaceChannel newChannel ((PhaseSpaceChannel(newmode),0,gamma,0,1,1,2,1,3));
    newChannel.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,-1.1);
    newChannel.weight(1.-weight_[ix]);
    newmode->addChannel(newChannel);
    // vmd channel
    PhaseSpaceChannel newChannel2((PhaseSpaceChannel(newmode),0,rho,0,1,1,2,1,3));
    newChannel2.weight(weight_[ix]);
    newmode->addChannel(newChannel2);
    addMode(newmode);
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

int VectorMesonPScalarFermionsDecayer::modeNumber(bool & cc,tcPDPtr parent,
						  const tPDVector & children) const {
  int imode(-1);
  // must be three outgoing particles
  if(children.size()!=3){return imode;}
  // ids of the particles
  int id0(parent->id()),idf[2]={0,0},ids(0);
  unsigned int nf(0);
  tPDVector::const_iterator pit = children.begin();
  for( ;pit!=children.end();++pit) {
    if((**pit).iSpin()==PDT::Spin0) ids=(**pit).id();
    else {
      idf[nf]=(**pit).id();
      ++nf;
    }
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==ids) {
      if((idf[0]==outgoing_[ix].second&&idf[1]==-outgoing_[ix].second)||
	 (idf[1]==outgoing_[ix].second&&idf[0]==-outgoing_[ix].second)) imode=ix;
    }
    ++ix;
  }
  while(imode<0&&ix<incoming_.size());
  // perform the decay
  cc=false;
  return imode;
}

void VectorMesonPScalarFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1/GeV) << incoming_ << outgoing_
     << maxweight_ << weight_ << includeVMD_ << VMDid_ << ounit(VMDmass_,GeV) 
     << ounit(VMDwidth_,GeV);
}

void VectorMesonPScalarFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1/GeV) >> incoming_ >> outgoing_
     >> maxweight_ >> weight_ >> includeVMD_ >> VMDid_ >> iunit(VMDmass_,GeV) 
     >> iunit(VMDwidth_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonPScalarFermionsDecayer,DecayIntegrator>
describeHerwigVectorMesonPScalarFermionsDecayer("Herwig::VectorMesonPScalarFermionsDecayer", "HwVMDecay.so");

void VectorMesonPScalarFermionsDecayer::Init() {

  static ClassDocumentation<VectorMesonPScalarFermionsDecayer> documentation
    ("The VectorMesonPScalarFermionsDecayer class is designed to "
     "perform the decay of a vector meson to a pseudoscalar meson and a "
     "fermion-antifermion pair.");

  static Command<VectorMesonPScalarFermionsDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, fermion, coupling(1/GeV), includeVMD,"
     " VMD id, VMD mass(GeV), VMD width(GeV) and max weight and weight of second channel for a decay."
     "There are three options for 0 the VMD factor is not included, for 1 the factor "
     "is included using the default mass and width of the particle specified by"
     " VMDID, and for 2 the factor is included using the mass and width specified"
     " by VMDwidth and VMDmass.",
     &VectorMesonPScalarFermionsDecayer::setUpDecayMode, false);

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceOutcomingP
    ("OutgoingPseudoScalar","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceOutcomingF
    ("OutgoingFermion","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceOutcomingA
    ("OutgoingAntiFermion","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceWeight
    ("Weight","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceIncludeVMD
    ("IncludeVMD","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceVMDID
    ("VMDID","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceVMDmass
    ("VMDmass","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonPScalarFermionsDecayer> interfaceVMDwidth
    ("VMDwidth","The old methods of setting up a decay in VectorMesonPScalarFermionsDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMesonPScalarFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  ScalarWaveFunction::constructSpinInfo(decay[0],outgoing,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[1],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[2],outgoing,true);
}

double VectorMesonPScalarFermionsDecayer::me2(const int, const Particle & part,
					      const tPDVector & ,
					      const vector<Lorentz5Momentum> & momenta,
					      MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,
					 PDT::Spin1Half,PDT::Spin1Half)));
  // initialization
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  wave_.resize(2);
  wavebar_.resize(2);
  for(unsigned int ix=0;ix<2;++ix) {
    wavebar_[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[1],ix,Helicity::outgoing);
    wave_   [ix] = HelicityFunctions::dimensionedSpinor   (-momenta[2],ix,Helicity::outgoing);
  }
  // the factor for the off-shell photon
  Lorentz5Momentum pff(momenta[1]+momenta[2]);
  pff.rescaleMass();
  Energy2 mff2(pff.mass2());
  // prefactor
  complex<InvEnergy3> pre(coupling_[imode()]/mff2);
  Complex ii(0.,1.);
  // the VMD factor
  if(includeVMD_[imode()]>0) {
    Energy2 mrho2(VMDmass_[imode()]*VMDmass_[imode()]);
    Energy2 mwrho(VMDmass_[imode()]*VMDwidth_[imode()]);
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  // calculate the matrix element
  LorentzPolarizationVector temp;
  unsigned int ix,iy,iz;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      temp=pre*epsilon(part.momentum(),pff,
				    wave_[ix].vectorCurrent(wavebar_[iy]));
      for(iz=0;iz<3;++iz) 
	(*ME())(iz,0,iy,ix)=temp.dot(vectors_[iz]); 
    }
  }
  // code for the spin averaged me for testing only
  // Energy  m[4]={part.mass(),momenta[0].mass(),
  // 		momenta[1].mass(),momenta[2].mass()};
  // Energy2 m2[4]={m[0]*m[0],m[1]*m[1],m[2]*m[2],m[3]*m[3]};
  // Lorentz5Momentum p12=momenta[0]+momenta[1];p12.rescaleMass();
  // Energy2 m122(p12.mass2());
  // cout << "testing the matrix element " 
  //      << -1./3.*(pre*conj(pre)).real()*(-2*m122*m122*mff2 - mff2*mff2*mff2 + 
  // 		   m2[1]*(2*m2[2]*m2[3] - 2*m2[3]*m2[3] + 
  // 			  m2[1]*(m2[2] - 2*m[2]*m[3] - m2[3])) - 
  // 		   2*m[2]*(m2[2]*m[2] - 2*m2[1]*m[3] - m[2]*m2[3])*
  // 		   m2[0] - (m2[2] + 2*m[2]*m[3] - m2[3])*
  // 		   m2[0]*m2[0] + mff2*mff2*
  // 		   (2*m2[1] + (m[2] - m[3])*(m[2] - m[3]) + 2*m2[0]) - 
  // 		   mff2*(m2[1]*m2[1] + 2*m2[1]*m[2]*(m[2] - 2*m[3]) + 
  // 			 2*m2[2]*m2[3] - 2*(2*m[2] - m[3])*m[3]*m2[0] + 
  // 			 m2[0]*m2[0]) + 2*m122*
  // 		   (-mff2*mff2 - (m[2] - m[3])*(m[2] + m[3])*(m[1] - m[0])*
  // 		    (m[1] + m[0]) + mff2*
  // 		    (m2[1] + m2[2] + m2[3] + m2[0])))
  //      << endl;
  // return the answer
  return ME()->contract(rho_).real();
}

WidthCalculatorBasePtr 
VectorMesonPScalarFermionsDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // workout which mode we are doing
  int imode(-1);
  // ids of the particles
  int id0(dm.parent()->id());
  int idf[2] = {0,0};
  int ids(0);
  unsigned int nf(0);
  ParticleMSet::const_iterator pit = dm.products().begin();
  for( ;pit!=dm.products().end();++pit) {
    if((**pit).iSpin()==PDT::Spin0){ids=(**pit).id();}
    else{idf[nf]=(**pit).id();++nf;}
  }
  // loop over the modes and see if this is one of them
  unsigned int ix=0;
  do {
    if(incoming_[ix]==id0&&outgoing_[ix].first==ids) {
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
    new_ptr(ThreeBodyAllOn1IntegralCalculator<VectorMesonPScalarFermionsDecayer>
	    (3,-1000.*MeV,-0.8*MeV,-0.8,*this,imode,m[0],m[1],m[2]));
} 

InvEnergy VectorMesonPScalarFermionsDecayer::
threeBodydGammads(int imodeb, const Energy2 q2, const Energy2 mff2, 
		  const  Energy m1, const Energy m2, const  Energy m3) const {
  // the masses of the external particles
  Energy q(sqrt(q2));
  Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
  // prefactor
  complex<InvEnergy3> pre(coupling_[imodeb]/mff2);
  Complex ii(0.,1.);
  // the VMD factor
  if(includeVMD_[imodeb]>0) {
    Energy2 mrho2(VMDmass_[imodeb]*VMDmass_[imodeb]),
      mwrho(VMDmass_[imodeb]*VMDwidth_[imodeb]);
    pre = pre*(-mrho2+ii*mwrho)/(mff2-mrho2+ii*mwrho);
  }
  InvEnergy6 factor(real(pre*conj(pre)));
  // compute the pieces from the integration limits
  Energy mff(sqrt(mff2)),e2star(0.5*(mff2-m32+m22)/mff),e1star(0.5*(q2-mff2-m12)/mff),
    e1sm(sqrt(e1star*e1star-m12)),e2sm(sqrt(e2star*e2star-m22));
  Energy2 a(2*e1star*e2star+m12+m22),b(2*e1sm*e2sm);
  // term independent of s3
  Energy8 me = 2*b*(-mff2*mff2*mff2 +m12*m12*(m22 - 2*m2*m3 - m32) - 
		   2*m22*(m22 - m32)*q2 -(m22 + 2*m2*m3 - m32)*q2*q2 + 
		   mff2*mff2*(2*m12 +(m2-m3)*(m2-m3)+2*q2) + 2*m12*m3*
		   ((m22-m32)*m3 + 2*m2*q2) - 
		   mff2*(m12*m12 + 2*m12*m2*(m2 - 2*m3) + 2*m22*m32 - 
			 2*(2*m2 - m3)*m3*q2 + q2*q2));
  // linear term
  me+= 2.*a*b*(2*(-mff2*mff2-(m22-m32)*(m12-q2)+mff2*(m12+m22+m32+q2)));
  // quadratic term
  me+=2*b*(3.*a*a+b*b)/3.*(-2*mff2);
  // phase space factors
  using Constants::pi;
  return -factor * me/768./pi/pi/pi/q2/q;
}

// output the setup information for the particle database
void VectorMesonPScalarFermionsDecayer::dataBaseOutput(ofstream & output,
						       bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first  << " " << outgoing_[ix].second << "  " 
	   << coupling_[ix]*GeV << " " << includeVMD_[ix] << " "
	   << VMDid_[ix] << " " << VMDmass_[ix]/GeV    << " "
	   << VMDwidth_[ix]/GeV << " " << maxweight_[ix]  << " "
	   << weight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}


string VectorMesonPScalarFermionsDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt2 = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g/GeV);
  includeVMD_.push_back(VMDinclude);
  VMDid_.push_back(VMDid);
  VMDmass_.push_back(VMDmass*GeV);
  VMDwidth_.push_back(VMDwidth*GeV);
  maxweight_.push_back(wgt);
  weight_.push_back(wgt2);
  // success
  return "";
}
#line 1 "./VectorMesonVectorPScalarDecayer.cc"
// -*- C++ -*-
//
// VectorMesonVectorPScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorPScalarDecayer class.
//

#include "VectorMesonVectorPScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorPScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) {
      if(mode(ix)) maxweight_[ix] = mode(ix)->maxWeight();
    }
  }
}

void VectorMesonVectorPScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size()  || 
     isize!=maxweight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMesonVectorPScalarDecayer" << Exception::runerror;
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

VectorMesonVectorPScalarDecayer::VectorMesonVectorPScalarDecayer()  {
  // intermediates
  generateIntermediates(false);
}

int VectorMesonVectorPScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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


void VectorMesonVectorPScalarDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxweight_ << ounit(coupling_,1/MeV);
}

void VectorMesonVectorPScalarDecayer::persistentInput(PersistentIStream & is, int)  {
  is >> incoming_ >> outgoing_ >> maxweight_ >> iunit(coupling_,1/MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonVectorPScalarDecayer,DecayIntegrator>
describeHerwigVectorMesonVectorPScalarDecayer("Herwig::VectorMesonVectorPScalarDecayer", "HwVMDecay.so");

void VectorMesonVectorPScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorPScalarDecayer> documentation
    ("The VectorMesonVectorPScalarDecayer class is designed for the "
     "decay of a vector meson to another vector meson, or the photon, and a "
     "pseudoscalar meson.");

  static Command<VectorMesonVectorPScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, pesudovector), coupling(1/MeV) and max weight for a decay",
     &VectorMesonVectorPScalarDecayer::setUpDecayMode, false);

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceOutcomingVector
    ("OutgoingVector","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceOutcomingPScalar
    ("OutgoingPScalar","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorPScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonVectorPScalarDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMesonVectorPScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[0],
					outgoing,true,decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

double VectorMesonVectorPScalarDecayer::me2(const int,const Particle & part,
					    const tPDVector & ,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  // is the vector massless
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  
  // compute the matrix element
  LorentzPolarizationVector vtemp;
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1&&photon) {
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=0.;
    }
    else {
      vtemp=coupling_[imode()]/part.mass()*
	epsilon(part.momentum(),vectors_[1][ix],momenta[0]);
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=vectors_[0][iy].dot(vtemp);
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = sqr(coupling_[imode()]*pcm)*2./3.;
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMesonVectorPScalarDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
    if(id==incoming_[ix]) {
      if(id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second) {
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
  coupling = coupling_[imode]*dm.parent()->mass();  
  mecode = 1;
  return order;
}

// output the setup info for the particle database
void VectorMesonVectorPScalarDecayer::dataBaseOutput(ofstream & output,
						     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix]*MeV << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string VectorMesonVectorPScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  InvEnergy g = stof(stype)/MeV;
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
#line 1 "./VectorMesonVectorScalarDecayer.cc"
// -*- C++ -*-
//
// VectorMesonVectorScalarDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorScalarDecayer class.
//
#include "VectorMesonVectorScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorScalarDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix) 
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void VectorMesonVectorScalarDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "VectorMesonVectorScalarDecayer::doinit()" 
			  << Exception::abortnow;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

int VectorMesonVectorScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void VectorMesonVectorScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,1/GeV);
}

void VectorMesonVectorScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,1/GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonVectorScalarDecayer,DecayIntegrator>
describeHerwigVectorMesonVectorScalarDecayer("Herwig::VectorMesonVectorScalarDecayer", "HwVMDecay.so");

void VectorMesonVectorScalarDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorScalarDecayer> documentation
    ("The VectorMesonVectorScalarDecayer class is designed for the "
     "decay of a vector meson to a vector meson, or the photon, and a "
     "scalar meson.");

  static Command<VectorMesonVectorScalarDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vector, scalar, coupling(1/GeV) and max weight for a decay",
     &VectorMesonVectorScalarDecayer::setUpDecayMode, false);

  static Deleted<VectorMesonVectorScalarDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonVectorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorScalarDecayer> interfaceOutcoming1
    ("FirstOutgoing","The old methods of setting up a decay in VectorMesonVectorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorScalarDecayer> interfaceOutcoming2
    ("SecondOutgoing","The old methods of setting up a decay in VectorMesonVectorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorScalarDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonVectorScalarDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<VectorMesonVectorScalarDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonVectorScalarDecayer have been deleted, please use SetUpDecayMode");

}

void VectorMesonVectorScalarDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[0],
					outgoing,true,decay[0]->id()==ParticleID::gamma);
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}
  
double VectorMesonVectorScalarDecayer::me2(const int,const Particle & part,
					   const tPDVector &,
					   const vector<Lorentz5Momentum> & momenta,
					   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  // is the vector massless
  bool photon(outgoing_[imode()].first==ParticleID::gamma);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // compute the matrix element
  Energy2 p0dotpv(part.momentum()*momenta[0]);
  complex<Energy> epsdot(ZERO);
  InvEnergy2 pre(coupling_[imode()]/part.mass());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1&&photon) {
      for(unsigned int iy=0;iy<3;++iy) (*ME())(iy,ix,0)=0.;
    }
    else {
      epsdot=vectors_[1][ix]*part.momentum();
      for(unsigned int iy=0;iy<3;++iy) {
	(*ME())(iy,ix,0)=Complex(pre*vectors_[0][iy].dot(p0dotpv*vectors_[1][ix]-
							 epsdot*momenta[0]));
      }
    }
  }
  double me = ME()->contract(rho_).real(); 
  // test of the matrix element
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = sqr(coupling_[imode()])/3.*(2.*sqr(pcm)+3.*sqr(momenta[0].mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMesonVectorScalarDecayer::twoBodyMEcode(const DecayMode & dm,
						     int & mecode,
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
    if(id==incoming_[ix]) {
      if(id1   ==outgoing_[ix].first&&id2   ==outgoing_[ix].second) {
	imode=ix;
	order=true;
      }
      if(id2   ==outgoing_[ix].first&&id1   ==outgoing_[ix].second) {
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
  coupling = coupling_[imode]*dm.parent()->mass();  
  mecode = 4;
  return order;
}

void VectorMesonVectorScalarDecayer::dataBaseOutput(ofstream & output,
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
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

string VectorMesonVectorScalarDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  coupling_.push_back(g/GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

#line 1 "./VectorMesonVectorVectorDecayer.cc"
// -*- C++ -*-
//
// VectorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonVectorVectorDecayer class.
//

#include "VectorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void VectorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxweight_[ix]=mode(ix)->maxWeight();
  }
}

void VectorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()  ||
     isize!=maxweight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "VectorMesonVectorVectorDecayer" << Exception::runerror;
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

VectorMesonVectorVectorDecayer::VectorMesonVectorVectorDecayer() {
  // intermediates
  generateIntermediates(false);
}

int VectorMesonVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void VectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxweight_ << coupling_;
}

void VectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxweight_ >> coupling_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigVectorMesonVectorVectorDecayer("Herwig::VectorMesonVectorVectorDecayer", "HwVMDecay.so");

void VectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonVectorVectorDecayer> documentation
    ("The VectorMesonVectorVectorDecayer class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static Command<VectorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling and max weight for a decay",
     &VectorMesonVectorVectorDecayer::setUpDecayMode, false);
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceOutgoing1
    ("Outgoing1","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceOutgoing2
    ("Outgoing2","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<VectorMesonVectorVectorDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in VectorMesonVectorVectorDecayer have been deleted, please use SetUpDecayMode");

}

string VectorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  coupling_.push_back(g);
  maxweight_.push_back(wgt);
  // success
  return "";
}

void VectorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix+1],decay[ix],
					  outgoing,true,decay[ix]->id()==ParticleID::gamma);
}

double VectorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix) 
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix+1].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix+1][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],
								   ihel,Helicity::outgoing);
    }
  }
  // work out the dot products we need for the matrix element
  Energy2 p1p2((momenta[0])*(momenta[1]));
  complex<Energy> p1eps2[3],p2eps1[3];
  for(unsigned int ix=0;ix<3;++ix) {
    p1eps2[ix]=vectors_[2][ix]*(momenta[0]);
    p2eps1[ix]=vectors_[1][ix]*(momenta[1]);
  }
  // compute the matrix element
  Lorentz5Momentum pdiff(momenta[0]-momenta[1]);
  Energy2 m12(momenta[0].mass()*momenta[0].mass()),m22(momenta[1].mass()*momenta[1].mass());
  InvEnergy3 fact(2.*coupling_[imode()]/(part.mass()*part.mass()*part.mass()));
  LorentzPolarizationVector vtemp;
  for(unsigned int ipol1=0;ipol1<3;++ipol1) {
    for(unsigned int ipol2=0;ipol2<3;++ipol2) {
      Complex eps1eps2=vectors_[1][ipol1].dot(vectors_[2][ipol2]);
      vtemp=fact*(p1eps2[ipol2]*p2eps1[ipol1]*pdiff
		  +p1eps2[ipol2]*m22*vectors_[1][ipol1]
		  -p2eps1[ipol1]*m12*vectors_[2][ipol2]
		  +eps1eps2*(-p1p2*pdiff+m12*momenta[1]
			     -m22*momenta[0]));
      for(unsigned int inpol=0;inpol<3;++inpol) 
	(*ME())(inpol,ipol1,ipol2)=vectors_[0][inpol].dot(vtemp);
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element;
  // Energy pcm=Kinematics::pstarTwoBodyDecay(part.mass(),momenta[0].mass(),
  // 					   momenta[1].mass());
  // double test = 8./3.*sqr(coupling_[imode()]*pcm/part.mass())*
  //   (1.+sqr(momenta[0].mass()/part.mass())+sqr(momenta[1].mass()/part.mass()));
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool VectorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling = coupling_[imode]; 
  mecode = 5;
  return order; 
}

// output the setup information for the particle database
void VectorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header){output << "update decayers set parameters=\"";}
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
#line 1 "./VectorMesonTensorVectorDecayer.cc"
// -*- C++ -*-
//
// VectorMesonTensorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorMesonTensorVectorDecayer class.
//

#include "VectorMesonTensorVectorDecayer.h"
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

void VectorMesonTensorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void VectorMesonTensorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistence of the parameters
  unsigned int isize=incoming_.size();
  if(isize!=outgoing_.size() || isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters VectorMesonTensorVectorDecayer" 
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

int VectorMesonTensorVectorDecayer::modeNumber(bool & cc, tcPDPtr parent, 
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


void VectorMesonTensorVectorDecayer::persistentOutput(PersistentOStream & os) const  {
  os << incoming_ << outgoing_ << maxWeight_ << ounit(coupling_,GeV);
}

void VectorMesonTensorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> iunit(coupling_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorMesonTensorVectorDecayer,DecayIntegrator>
describeHerwigVectorMesonTensorVectorDecayer("Herwig::VectorMesonTensorVectorDecayer", "HwTMDecay.so");

void VectorMesonTensorVectorDecayer::Init() {

  static ClassDocumentation<VectorMesonTensorVectorDecayer> documentation
    ("The VectorMesonTensorVectorDecayer class performs the"
     " decay of a tensor meson to two scalar mesons.");

  static Command<VectorMesonTensorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, vectors, coupling(GeV) and max weight for a decay",
     &VectorMesonTensorVectorDecayer::setUpDecayMode, false);

}

void VectorMesonTensorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  TensorWaveFunction::constructSpinInfo(tensors_,decay[0],
					outgoing,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[1],
					outgoing,true,decay[1]->id()==ParticleID::gamma);
}

// matrix elememt for the process
double VectorMesonTensorVectorDecayer::me2(const int,const Particle & part,
					    const tPDVector & outgoing,
					    const vector<Lorentz5Momentum> & momenta,
					    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin2,PDT::Spin1)));
  // check for photons
  bool photon = outgoing[1]->id()==ParticleID::gamma;
  // // stuff for incoming particle
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  TensorWaveFunction twave(momenta[0],outgoing[0],Helicity::outgoing);
  tensors_.resize(5);
  for(unsigned int ihel=0;ihel<5;++ihel) {
    twave.reset(ihel);
    tensors_[ihel] = twave.wave();
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  Energy2 denom[2] = {momenta[0]*part.momentum()-momenta[0].mass()*   part   .mass(),
		      momenta[0]*momenta[1]     -momenta[0].mass()*momenta[1].mass()};
  double fact(coupling_[imode()]/part.mass());
  // calculate the matrix element
  for(unsigned int ih1=0;ih1<5;++ih1) {
    LorentzPolarizationVectorE v1 = tensors_[ih1].preDot (part.momentum());
    LorentzPolarizationVectorE v2 = tensors_[ih1].postDot(momenta[1]);
    complex<Energy2> d0 = v1*momenta[1];
    for(unsigned int ih0=0;ih0<3;++ih0) {
      LorentzPolarizationVector v3 = tensors_[ih1].postDot(vectors_[0][ih0]);
      complex<InvEnergy> d1 = vectors_[0][ih0]*momenta[0]/denom[0];
      complex<Energy> d4 = v2*vectors_[0][ih0];
      for(unsigned int ih2=0;ih2<3;++ih2) {
	complex<InvEnergy> d2 = vectors_[1][ih2]*momenta[0]/denom[1];
	if ( photon && ih2==1)
	  (*ME())(ih0,ih1,ih2)=0.;
	else {
	  (*ME())(ih0,ih1,ih2)=fact*(v3*vectors_[1][ih2] -(v1*vectors_[1][ih2])*d1
				     -d4*d2+d0*d1*d2);
	}
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the answer
  // double test = 5.*sqr(coupling_[imode()]/part.mass())/3.;
  // if(photon) test *=2./3.;
  // cout << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " " 
  //      << me << " " << test << " " << (me-test)/(me+test) << endl;
  // return the answer
  return me;
}

bool VectorMesonTensorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  mecode=19;
  return id1==outgoing_[imode].first&&id2==outgoing_[imode].second;
}

void VectorMesonTensorVectorDecayer::dataBaseOutput(ofstream & output,
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

string VectorMesonTensorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  coupling_.push_back(g*GeV);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./PseudoVectorMesonVectorVectorDecayer.cc"
// -*- C++ -*-
//
// PseudoVectorMesonVectorVectorDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PseudoVectorMesonVectorVectorDecayer class.
//

#include "PseudoVectorMesonVectorVectorDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PseudoVectorMesonVectorVectorDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<incoming_.size();++ix)
      if(mode(ix)) maxWeight_[ix]=mode(ix)->maxWeight();
  }
}

void PseudoVectorMesonVectorVectorDecayer::doinit() {
  DecayIntegrator::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()  ||
     isize!=maxWeight_.size() || isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "PseudoVectorMesonVectorVectorDecayer" << Exception::runerror;
  // set up the integration channels
  PhaseSpaceModePtr mode;
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
		     getParticleData(outgoing_[ix].second)};
    if(in&&out[0]&&out[1]) 
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    else
      mode=PhaseSpaceModePtr();
    addMode(mode);
  }
}

PseudoVectorMesonVectorVectorDecayer::PseudoVectorMesonVectorVectorDecayer() {
  // intermediates
  generateIntermediates(false);
}

int PseudoVectorMesonVectorVectorDecayer::modeNumber(bool & cc,tcPDPtr parent,
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

void PseudoVectorMesonVectorVectorDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_;
}

void PseudoVectorMesonVectorVectorDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PseudoVectorMesonVectorVectorDecayer,DecayIntegrator>
describeHerwigPseudoVectorMesonVectorVectorDecayer("Herwig::PseudoVectorMesonVectorVectorDecayer", "HwVMDecay.so");

void PseudoVectorMesonVectorVectorDecayer::Init() {

  static ClassDocumentation<PseudoVectorMesonVectorVectorDecayer> documentation
    ("The PseudoVectorMesonVectorVectorDecayer class is designed for the "
     "decay of a vector meson to two vector particles, either photons or other "
     "vector mesons.");

  static Command<PseudoVectorMesonVectorVectorDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles, coupling and max weight for a decay",
     &PseudoVectorMesonVectorVectorDecayer::setUpDecayMode, false);

}

string PseudoVectorMesonVectorVectorDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
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
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}

void PseudoVectorMesonVectorVectorDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(vectors_[ix+1],decay[ix],
					  outgoing,true,decay[ix]->id()==ParticleID::gamma);
}

double PseudoVectorMesonVectorVectorDecayer::me2(const int,const Particle & part,
					const tPDVector & outgoing,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix) 
    photon[ix] = outgoing[ix]->id()==ParticleID::gamma;
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    vectors_[ix+1].resize(3);
    for(unsigned int ihel=0;ihel<3;++ihel) {
      if(photon[ix] && ihel==1) continue;
      vectors_[ix+1][ihel] = HelicityFunctions::polarizationVector(-momenta[ix],
								   ihel,Helicity::outgoing);
    }
  }
  // work out the dot products we need for the matrix element
  complex<InvEnergy> fact(coupling_[imode()]/part.mass());
  Lorentz5Momentum pref = part.momentum();
  if(photon[0]) pref = momenta[0];
  else if(photon[1]) pref =momenta[1];
  for(unsigned int ih0=0;ih0<3;++ih0) {
    for(unsigned int ih1=0;ih1<3;++ih1) {
      if(photon[0] && ih1==1) continue;
      LorentzPolarizationVectorE       v0 = epsilon(pref,vectors_[0][ih0],vectors_[1][ih1]);
      for(unsigned int ih2=0;ih2<3;++ih2) {
	if(photon[1] && ih2==1) continue;
   	(*ME())(ih0,ih1,ih2)= Complex(fact*(v0*vectors_[2][ih2]));
      }
    }
  }
  double me = ME()->contract(rho_).real();
  // test of the matrix element;
  // double test;
  // if(photon[0]||photon[1]) {
  //   Energy2 mout2 = sqr(momenta[0].mass()+momenta[1].mass());
  //   test = sqr(coupling_[imode()])/6.*sqr(sqr(part.mass()) - mout2)*(sqr(part.mass()) + mout2)/(pow<4,1>(part.mass())*mout2);
  // }
  // else {
  //   Energy2 m02(sqr(part.mass())),m12(sqr(momenta[0].mass())),m22(sqr(momenta[1].mass()));
  //   test = sqr(coupling_[imode()])/6.*(sqr(m02)*(m12 + m22) + sqr(m12 - m22)*(m12 + m22) - 
  // 				       2*m02*(sqr(m12) - 4*m12*m22 + sqr(m22)))/(m02*m12*m22);
  // }
  // cerr << "testing matrix element for " << part.PDGName() << " -> " 
  //      << outgoing[0]->PDGName() << " " << outgoing[1]->PDGName() << " "
  //      << me << " " << test << " " << (me-test)/(me+test) << "\n";
  // return the answer
  return me;
}

bool PseudoVectorMesonVectorVectorDecayer::twoBodyMEcode(const DecayMode & dm,int & mecode,
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
  coupling = coupling_[imode]; 
  mecode = 20;
  return order; 
}

// output the setup information for the particle database
void PseudoVectorMesonVectorVectorDecayer::dataBaseOutput(ofstream & output,
						    bool header) const {
  if(header){output << "update decayers set parameters=\"";}
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " "
	   << coupling_[ix] << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./OniumToOniumPiPiDecayer.cc"
// -*- C++ -*-
//
// OniumToOniumPiPiDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OniumToOniumPiPiDecayer class.
//

#include "OniumToOniumPiPiDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void OniumToOniumPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
    if(initialize()) maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

void OniumToOniumPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // check consistency of the vectors
  unsigned int isize=incoming_.size();
  if(outgoing_.size()!=isize||maxWeight_.size()!=2*isize||
     coupling_.size()!=isize||cA_.size()      !=isize||
     cB_      .size()!=isize||cC_.size()      !=isize)
    throw InitException() << "Inconsistent size of the parameter vectors in "
			  << "OniumToOniumPiPiDecayer"
			  << Exception::runerror;
  // construct the decay channels
  tPDPtr pip(getParticleData(ParticleID::piplus ));
  tPDPtr pim(getParticleData(ParticleID::piminus));
  tPDPtr pi0(getParticleData(ParticleID::pi0    ));
  tPDPtr rho0(getParticleData(113)); 
  for(unsigned int ix=0;ix<isize;++ix) {
    tPDPtr     in = getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix]),pip,pim};
    for(unsigned int iy=0;iy<2;++iy) {
      // pi0 pi0
      if(iy==1) {
	out[1]=pi0;
	out[2]=pi0;
      }
      // construct the phase-space mode
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1,0,rho0,1,2,1,3));
      mode->addChannel(channel);
      // reset the resonance parameters
      mode->resetIntermediate(rho0,2*in->mass(),in->mass());
      // add the mode
      addMode(mode);
    }
  }
}

void OniumToOniumPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << maxWeight_ << coupling_
     << ounit(cA_,1./GeV2) << ounit(cB_,1./GeV2) << ounit(cC_,1./GeV2);
}

void OniumToOniumPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> maxWeight_ >> coupling_
     >> iunit(cA_,1./GeV2) >> iunit(cB_,1./GeV2) >> iunit(cC_,1./GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OniumToOniumPiPiDecayer,DecayIntegrator>
describeHerwigOniumToOniumPiPiDecayer("Herwig::OniumToOniumPiPiDecayer", "HwVMDecay.so");

void OniumToOniumPiPiDecayer::Init() {

  static ClassDocumentation<OniumToOniumPiPiDecayer> documentation
    ("The OniumToOniumPiPiDecayer class uses the matrix element of "
     "Brown and Cahn, PRL35, 1 (1975), for"
     " the decay of onium resonaces to lighter states and pion pairs."
     " The results of hep-ex/9909038 are used for psi'->psi and "
     " arXiv:0706.2317 for Upsilon(3S) and Upsilon(2S) decays."
     " The remaining parameters are choosen to approximately reproduce"
     " the distributions from hep-ex/0604031 and hep-ex/0508023.",
     "The decays of onium resonances to lighter states and pion pairs were modelled"
     " using the matrix element of \\cite{Brown:1975dz}. The results of "
     "\\cite{Bai:1999mj} are used for $\\psi'\\to\\psi$ and "
     "\\cite{Cronin-Hennessy:2007sj} for $\\Upsilon(3S)$ and $\\Upsilon(2S)$ decays."
     " The remaining parameters are choosen to approximately reproduce"
     " the distributions from \\cite{Aubert:2006bm} and \\cite{Adam:2005mr}.",
     "\\bibitem{Brown:1975dz} L.~S.~Brown and R.~N.~Cahn,"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 35} (1975) 1."
     "%%CITATION = PRLTA,35,1;%%\n"
     "\\bibitem{Bai:1999mj} J.~Z.~Bai {\\it et al.}  [BES Collaboration],"
     "Phys.\\ Rev.\\  D {\\bf 62} (2000) 032002 [arXiv:hep-ex/9909038]."
     "%%CITATION = PHRVA,D62,032002;%%\n"
     "\\bibitem{Cronin-Hennessy:2007sj} D.~Cronin-Hennessy{\\it et al.} "
     "[CLEO Collaboration], arXiv:0706.2317 [hep-ex]."
     "%%CITATION = ARXIV:0706.2317;%%\n"
     "\\bibitem{Aubert:2006bm} B.~Aubert {\\it et al.}  [BABAR Collaboration],"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 96} (2006) 232001 [arXiv:hep-ex/0604031]."
     "%%CITATION = PRLTA,96,232001;%%\n"
     "\\bibitem{Adam:2005mr} N.~E.~Adam {\\it et al.}  [CLEO Collaboration],"
     "Phys.\\ Rev.\\ Lett.\\  {\\bf 96} (2006) 082004 [arXiv:hep-ex/0508023]."
     "%%CITATION = PRLTA,96,082004;%%");

  static Command<OniumToOniumPiPiDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing, coupling, A, B, C (re and im parts, 1/MeV2) and max weights for a decay",
     &OniumToOniumPiPiDecayer::setUpDecayMode, false);

  static Deleted<OniumToOniumPiPiDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceOutgoing
    ("Outgoing","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceCoupling
    ("Coupling","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceReA
    ("ReA","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceReB
    ("ReB","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceReC
    ("ReC","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceImA
    ("ImA","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceImB
    ("ImB","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<OniumToOniumPiPiDecayer> interfaceImC
    ("ImC","The old methods of setting up a decay in OniumToOniumPiPiDecayer have been deleted, please use SetUpDecayMode");
}

int OniumToOniumPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  cc=false;
  int imode(-1);
  long idin(parent->id());
  if(children.size()!=3) return -1;
  unsigned int npip(0),npim(0),npi0(0);
  long idother(0),id;
  for(tPDVector::const_iterator pit=children.begin();
      pit!=children.end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::pi0)     ++npi0;
    else idother=id;
  }
  // check pi+ pi- or pi0 pi0 and outgoing state
  if(!((npip==1&&npim==1)||npi0==2)||idother==0) return -1;
  unsigned int ix=0;
  do {
    if(idin==incoming_[ix]&&idother==outgoing_[ix]) imode=ix;
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return npi0==2 ? 2*imode+1 : 2*imode;
}

void OniumToOniumPiPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[0],
					outgoing,true,false);
  for(unsigned int ix=1;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double OniumToOniumPiPiDecayer::me2(const int,const Particle & part,
					const tPDVector & ,
					const vector<Lorentz5Momentum> & momenta,
					MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
						const_ptr_cast<tPPtr>(&part),
						incoming,false);
  }
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[0],ix,Helicity::outgoing);
  }
  // compute the matrix element
  complex<InvEnergy2> A(cA_[imode()/2]),B(cB_[imode()/2]),C(cC_[imode()/2]);
  Energy2 q2  =(momenta[1]+momenta[2]).m2();
  Energy2 mpi2=sqr(momenta[1].mass());
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Complex dota = vectors_[0][ix].dot(vectors_[1][iy]);
      complex<Energy2> dotb = 
	(vectors_[0][ix]*momenta[1])*(vectors_[1][iy]*momenta[2])+
	(vectors_[0][ix]*momenta[2])*(vectors_[1][iy]*momenta[1]);
      (*ME())(ix,iy,0,0)= coupling_[imode()/2]*
	Complex(A*dota*(q2-2.*mpi2)+B*dota*momenta[1].e()*momenta[2].e()
	 +C*dotb);
    }
  }
  // matrix element
  double output=ME()->contract(rho_).real();
  if(imode()%2==1) output*=0.5;
  // test of the matrix element
  // Energy2 s1=(momenta[1]+momenta[2]).m2();
  // Energy2 s2=(momenta[0]+momenta[2]).m2();
  // Energy2 s3=(momenta[1]+momenta[0]).m2();
  // double test=threeBodyMatrixElement(imode(),sqr(part.mass()),
  // 				     s3,s2,s1,momenta[0].mass(),
  // 				     momenta[1].mass(),momenta[2].mass());
  // cerr << "testing " << output << " " << test << " " << (output-test)/(output+test) << "\n";
  // return the answer
  return output;
}

// output the setup information for the particle database
void OniumToOniumPiPiDecayer::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  // the rest of the parameters
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix] << " " << coupling_[ix] << " "
	   << cA_[ix].real()*MeV2 << " " << cA_[ix].imag()*MeV2 << " "
	   << cB_[ix].real()*MeV2 << " " << cB_[ix].imag()*MeV2 << " "
	   << cC_[ix].real()*MeV2 << " " << cC_[ix].imag()*MeV2 << " "
	   << maxWeight_[2*ix] << " " << maxWeight_[2*ix+1] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}

WidthCalculatorBasePtr OniumToOniumPiPiDecayer::
threeBodyMEIntegrator(const DecayMode & dm) const {
  int imode(-1);
  long idin(dm.parent()->id());
  unsigned int npip(0),npim(0),npi0(0);
  long idother(0),id;
  for(ParticleMSet::const_iterator pit=dm.products().begin();
      pit!=dm.products().end();++pit) {
    id=(**pit).id();
    if(id==ParticleID::piplus)       ++npip;
    else if(id==ParticleID::piminus) ++npim;
    else if(id==ParticleID::pi0)     ++npi0;
    else idother=id;
  }
  unsigned int ix=0;
  do {
    if(idin==incoming_[ix]&&idother==outgoing_[ix]) imode=ix;
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  imode = npi0==2 ? 2*imode+1 : 2*imode;
  // construct the integrator
  vector<double> inweights(1,1.);
  Energy scale=getParticleData(incoming_[ix-1])->mass();
  Energy m1=getParticleData(outgoing_[ix-1])->mass();
  Energy mpi = npi0==2 ? getParticleData(ParticleID::pi0)->mass() :
    getParticleData(ParticleID::piplus)->mass();
  vector<int> intype(1,3);
  vector<Energy> inmass (1,scale);
  vector<Energy> inwidth(1,scale);
  vector<double> inpow(1,0.0);
  return new_ptr(ThreeBodyAllOnCalculator<OniumToOniumPiPiDecayer>
		 (inweights,intype,inmass,inwidth,inpow,
		  *this,imode,m1,mpi,mpi));
}

double OniumToOniumPiPiDecayer::
threeBodyMatrixElement(const int imode, const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, 
		       const Energy m1, const Energy m2, const Energy m3) const {
  Energy q=sqrt(q2);
  Energy e2 = 0.5*(q2+sqr(m2)-s2)/q;
  Energy e3 = 0.5*(q2+sqr(m3)-s3)/q;
  Complex amp = cA_[imode/2]*(s1-sqr(m2)-sqr(m3))+cB_[imode/2]*e2*e3;
  Energy2 dot = 0.5*(q2+sqr(m1)-s1);
  double output=(2.+sqr(dot/q/m1))*real(amp*conj(amp))*sqr(coupling_[imode/2])/3.;
  if(imode%2==1) output*=0.5;
  return output;
}

string OniumToOniumPiPiDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 1";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  complex<InvEnergy2> coup[3];
  Complex ii(0.,1.);
  for(unsigned int ix=0;ix<3;++ix) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double re = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double im = stof(stype);
    coup[ix] = re/MeV2+im/MeV2*ii;
  }
  pair<double,double> wgt;
  wgt.first = stof(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  wgt.second = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  cA_.push_back(coup[0]);
  cB_.push_back(coup[1]);
  cC_.push_back(coup[2]);
  maxWeight_.push_back(wgt.first);
  maxWeight_.push_back(wgt.second);
  // success
  return "";
}
#line 1 "./f1RhoPiPiDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the f1RhoPiPiDecayer class.
//

#include "f1RhoPiPiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

f1RhoPiPiDecayer::f1RhoPiPiDecayer() :
  ga1RhoPi_(4.8*GeV), gf1a1Pi_(9.77/GeV), maxWeight_({1.,1.}) {
  generateIntermediates(true);
}

// normally not implemented but do it here to get rid of the a_1 which
// can't be on-shell
ParticleVector f1RhoPiPiDecayer::decay(const Particle & parent,
				       const tPDVector & children) const {
  ParticleVector output = DecayIntegrator::decay(parent,children);
  ParticleVector::iterator it =output.begin();
  for(;it!=output.end();++it) {
    long id = (**it).id();
    if(id==20113 || id==20213 || id==-20213)
      break;
  }
  if(it!=output.end()) {
    PPtr a1 = *it;
    output.erase(it);
    for(PPtr child : a1->children()) {
      output.push_back(child);
      a1->abandonChild(child);
    }
  }
  return output;
}


IBPtr f1RhoPiPiDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr f1RhoPiPiDecayer::fullclone() const {
  return new_ptr(*this);
}

void f1RhoPiPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(ga1RhoPi_,GeV) << ounit(gf1a1Pi_,1./GeV)
     << ounit(ma1_,GeV) << ounit(ga1_,GeV) << maxWeight_;
}

void f1RhoPiPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ga1RhoPi_,GeV) >> iunit(gf1a1Pi_,1./GeV)
     >> iunit(ma1_,GeV) >> iunit(ga1_,GeV) >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<f1RhoPiPiDecayer,DecayIntegrator>
describeHerwigf1RhoPiPiDecayer("Herwig::f1RhoPiPiDecayer", "HwVMDecay.so");

void f1RhoPiPiDecayer::Init() {

  static ClassDocumentation<f1RhoPiPiDecayer> documentation
    ("The f1RhoPiPiDecayer class implements a simple model for "
     "f1 -> pipirho via an intermediate a_1");

  static ParVector<f1RhoPiPiDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "Maximum weights for the decays",
     &f1RhoPiPiDecayer::maxWeight_, 2, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<f1RhoPiPiDecayer,Energy> interfacega1RhoPi
    ("ga1RhoPi",
     "Coupling of the a_1 to rho pi",
     &f1RhoPiPiDecayer::ga1RhoPi_, GeV, 4.8*GeV, 0.*GeV, 20.*GeV,
     false, false, Interface::limited);
  
  static Parameter<f1RhoPiPiDecayer,InvEnergy> interfacegf1a1Pi
    ("gf1a1Pi",
     "Coupling of f_1 to a_1 pi",
     &f1RhoPiPiDecayer::gf1a1Pi_, 1./GeV, 1.0/GeV, 0.0/GeV, 10.0/GeV,
     false, false, Interface::limited);

}

void f1RhoPiPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
  tPDPtr f1 = getParticleData(ParticleID::f_1);
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  tPDPtr rhop = getParticleData(ParticleID::rhoplus);
  tPDPtr rhom = getParticleData(ParticleID::rhominus);
  tPDPtr rho0 = getParticleData(ParticleID::rho0);
  // decay mode f_1 -> pi+ pi- rho0
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(f1,{pip,pim,rho0},maxWeight_[0]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,3));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,1,1,3));
  addMode(mode);
  // decay mode f_1 -> pi- pi0 rho+
  mode = new_ptr(PhaseSpaceMode(f1,{pim,pi0,rhop},maxWeight_[1]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,1,1,2,1,3));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,1,1,3));
  addMode(mode);
  ma1_ = a1p->mass();
  ga1_ = a1p->width(); 
}

void f1RhoPiPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<maxWeight_.size();++ix)
      maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

int f1RhoPiPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  if(children.size()!=3 || parent->id()!=ParticleID::f_1) return -1;
  // check the pions
  int npi0(0),npiplus(0),npiminus(0),idrho(0);
  for(auto child : children) {
    int idtemp= child->id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
    else                                 idrho=idtemp;
  }
  cc = false;
  // f_1 -> pi+pi-rho0 mode
  if(idrho==113 && npiplus==1 && npiminus==1)         return 0;
  // f_1 -> pi-pi0rho+ mode
  else if  (idrho== 213 && npiminus==1 && npi0==1)    return 1;
  else if  (idrho==-213 && npiplus ==1 && npi0==1) {
    cc=true;
    return 1;
  }
  // not found
  return -1;
}

void f1RhoPiPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
  					incoming,true,false);
  // set up the spin information for the decay products
  // pions
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
  // rho
  VectorWaveFunction::constructSpinInfo(vectors_[1],decay[2],outgoing,true,false);
}

double f1RhoPiPiDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin1)));
  useMe();
  // polarization vectors for incoming
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
  }
  // polarization vectors for rho
  vectors_[1].resize(3);
  for(unsigned int ix=0;ix<3;++ix)
    vectors_[1][ix] = HelicityFunctions::polarizationVector(-momenta[2],ix,Helicity::outgoing);
  // breit wigners
  Lorentz5Momentum pa1[2] = {part.momentum()-momenta[0],
			     part.momentum()-momenta[1]};
  for(unsigned int ix=0;ix<2;++ix) pa1[ix].rescaleMass();
  complex<InvEnergy2> bwa1[2] = {Resonance::BreitWignera1(pa1[0].mass2(),ma1_,ga1_)/sqr(ma1_),
				 Resonance::BreitWignera1(pa1[1].mass2(),ma1_,ga1_)/sqr(ma1_)};
  if(ichan>0) bwa1[ ichan == 0 ? 1 : 0 ] = ZERO;
  // compute the matrix element
  for(unsigned int ihel=0;ihel<3;++ihel) {
    LorentzVector<complex<Energy2> > pol[2] = {epsilon(vectors_[0][ihel],part.momentum(),pa1[0]),
					       epsilon(vectors_[0][ihel],part.momentum(),pa1[1])};
    for(unsigned int ohel=0;ohel<3;++ohel) {
      (*ME())(ihel,0,0,ohel) = Complex(gf1a1Pi_*ga1RhoPi_*(LorentzPolarizationVector(pol[0]*bwa1[0])-
							   LorentzPolarizationVector(pol[1]*bwa1[1])).dot(vectors_[1][ohel]));
    }
  } 
  // matrix element
  return ME()->contract(rho_).real();
}

void f1RhoPiPiDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ga1RhoPi " << ga1RhoPi_/GeV << "\n";
  output << "newdef " << name() << ":gf1a1Pi "  << gf1a1Pi_*GeV  << "\n";
  for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
    output << "newdef    " << name() << ":maxWeight " << ix << " "
  	   << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
#line 1 "./f1FourPiDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the f1FourPiDecayer class.
//

#include "f1FourPiDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"

using namespace Herwig;

f1FourPiDecayer::f1FourPiDecayer() :
  gRhoPiPi_(6.), ga1RhoPi_(4.8*GeV), gf1a1Pi_(9.77/GeV), maxWeight_({1.,1.}) {
  generateIntermediates(true);
}

// normally not implemented but do it here to get rid of the a_1 which
// can't be on-shell
ParticleVector f1FourPiDecayer::decay(const Particle & parent,
				       const tPDVector & children) const {
  ParticleVector output = DecayIntegrator::decay(parent,children);
  ParticleVector::iterator it =output.begin();
  for(;it!=output.end();++it) {
    long id = (**it).id();
    if(id==20113 || id==20213 || id==-20213)
      break;
  }
  if(it!=output.end()) {
    PPtr a1 = *it;
    output.erase(it);
    for(PPtr child : a1->children()) {
      output.push_back(child);
      a1->abandonChild(child);
    }
  }
  return output;
}

IBPtr f1FourPiDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr f1FourPiDecayer::fullclone() const {
  return new_ptr(*this);
}

void f1FourPiDecayer::persistentOutput(PersistentOStream & os) const {
  os << gRhoPiPi_ << ounit(ga1RhoPi_,GeV) << ounit(gf1a1Pi_,1./GeV)
     << ounit(ma1_,GeV) << ounit(ga1_,GeV)
     << ounit(mrho_,GeV) << ounit(grho_,GeV) << maxWeight_;
}

void f1FourPiDecayer::persistentInput(PersistentIStream & is, int) {
  is >> gRhoPiPi_ >> iunit(ga1RhoPi_,GeV) >> iunit(gf1a1Pi_,1./GeV)
     >> iunit(ma1_,GeV) >> iunit(ga1_,GeV)
     >> iunit(mrho_,GeV) >> iunit(grho_,GeV) >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<f1FourPiDecayer,DecayIntegrator>
describeHerwigf1FourPiDecayer("Herwig::f1FourPiDecayer", "HwVMDecay.so");

void f1FourPiDecayer::Init() {

  static ClassDocumentation<f1FourPiDecayer> documentation
    ("The f1FourPiDecayer class implements a simple model for "
     "f1 -> pipirho via an intermediate a_1");

  static ParVector<f1FourPiDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "Maximum weights for the decays",
     &f1FourPiDecayer::maxWeight_, 2, 1.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Parameter<f1FourPiDecayer,double> interfacegRhoPiPi
    ("gRhoPiPi",
     "The coupling of the rho to two pions",
     &f1FourPiDecayer::gRhoPiPi_, 6., 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<f1FourPiDecayer,Energy> interfacega1RhoPi
    ("ga1RhoPi",
     "Coupling of the a_1 to rho pi",
     &f1FourPiDecayer::ga1RhoPi_, GeV, 4.8*GeV, 0.*GeV, 20.*GeV,
     false, false, Interface::limited);
  
  static Parameter<f1FourPiDecayer,InvEnergy> interfacegf1a1Pi
    ("gf1a1Pi",
     "Coupling of f_1 to a_1 pi",
     &f1FourPiDecayer::gf1a1Pi_, 1./GeV, 1.0/GeV, 0.0/GeV, 10.0/GeV,
     false, false, Interface::limited);

}

void f1FourPiDecayer::doinit() {
  DecayIntegrator::doinit();
  // pointers to the particles we need as external particles
  tPDPtr f1 = getParticleData(ParticleID::f_1);
  tPDPtr a1p = getParticleData(ParticleID::a_1plus);
  tPDPtr a1m = getParticleData(ParticleID::a_1minus);
  tPDPtr a10 = getParticleData(ParticleID::a_10);
  tPDPtr pip = getParticleData(ParticleID::piplus);
  tPDPtr pim = getParticleData(ParticleID::piminus);
  tPDPtr pi0 = getParticleData(ParticleID::pi0);
  tPDPtr rhop = getParticleData(ParticleID::rhoplus);
  tPDPtr rhom = getParticleData(ParticleID::rhominus);
  tPDPtr rho0 = getParticleData(ParticleID::rho0);
  // decay mode f_1 -> pi+ pi- pi+ pi-
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(f1,{pip,pim,pip,pim},maxWeight_[0]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rho0,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,1,1,rho0,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,2,1,rho0,2,1,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,3,1,rho0,2,1,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rho0,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,1,1,rho0,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,4,1,rho0,2,1,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,3,1,rho0,2,1,2,2));
  addMode(mode);
  // decay mode f_1 -> pi+ pi0 pi- pi0
  mode = new_ptr(PhaseSpaceMode(f1,{pip,pi0,pim,pi0},maxWeight_[0]));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rhop,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,1,1,rhop,2,3,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rhop,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,1,1,rhop,2,3,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,4,1,rhom,2,1,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,3,1,rhom,2,1,2,2));
  mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,2,1,rhom,2,1,2,4));
  mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,3,1,rhom,2,1,2,4));
  addMode(mode);
  // masses of intermediates
  ma1_ = a1p->mass();
  ga1_ = a1p->width(); 
  mrho_ = rhop->mass();
  grho_ = rhop->width(); 
}

void f1FourPiDecayer::doinitrun() {
  DecayIntegrator::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<maxWeight_.size();++ix)
      maxWeight_[ix] = mode(ix)->maxWeight();
  }
}

int f1FourPiDecayer::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  if(children.size()!=4 || parent->id()!=ParticleID::f_1) return -1;
  // check the pions
  int npi0(0),npiplus(0),npiminus(0);
  for(auto child : children) {
    int idtemp= child->id();
    if(idtemp==ParticleID::piplus)       ++npiplus;
    else if(idtemp==ParticleID::piminus) ++npiminus;
    else if(idtemp==ParticleID::pi0)     ++npi0;
  }
  cc = false;
  // f_1 -> 2pi+ 2pi-mode
  if(npiplus==2 && npiminus==2)                    return 0;
  // f_1 -> pi+ pi- pi0 pi0 mode
  else if  (npiplus ==1 && npiminus==1 && npi0==2) return 1;
  // not found
  return -1;
}

void f1FourPiDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vector_,const_ptr_cast<tPPtr>(&part),
  					incoming,true,false);
  // set up the spin information for the pions
  for(unsigned int ix=0;ix<4;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double f1FourPiDecayer::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  // polarization vectors for incoming
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vector_,rho_,
  					       const_ptr_cast<tPPtr>(&part),
  					       incoming,false);
  }
  // breit wigners
  Lorentz5Momentum pa1[4] = {part.momentum()-momenta[0],part.momentum()-momenta[1],
			     part.momentum()-momenta[2],part.momentum()-momenta[3]};
  complex<InvEnergy2> bwa1[4];
  for(unsigned int ix=0;ix<4;++ix) {
    pa1[ix].rescaleMass();
    bwa1[ix] = Resonance::BreitWignera1(pa1[ix].mass2(),ma1_,ga1_)/sqr(ma1_);
  }
  // compute the matrix element (as a current and then contract)
  LorentzVector<complex<InvEnergy> > current;
  // decay mode f_1 -> pi+ pi- pi+pi-
  double sym(0.5);
  if(imode()==0) {
    sym=0.25;
    complex<InvEnergy2> bwrho[4] =
      {Resonance::BreitWignerPWave((momenta[0]+momenta[1]).m2(),mrho_,grho_,momenta[0].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[0]+momenta[3]).m2(),mrho_,grho_,momenta[0].mass(),momenta[3].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[1]).m2(),mrho_,grho_,momenta[2].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[3]).m2(),mrho_,grho_,momenta[2].mass(),momenta[3].mass())/sqr(mrho_)};
    if(ichan<=0) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rho0,2,3,2,4));
      current += bwa1[0]*bwrho[3]*epsilon(part.momentum(),pa1[0],momenta[2]-momenta[3]);
    }
    if(ichan<0||ichan==1) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,1,1,rho0,2,3,2,4));
      current -= bwa1[1]*bwrho[3]*epsilon(part.momentum(),pa1[1],momenta[2]-momenta[3]);
    }
    if(ichan<0||ichan==2) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,2,1,rho0,2,1,2,4));
      current += bwa1[2]*bwrho[1]*epsilon(part.momentum(),pa1[2],momenta[0]-momenta[3]);
    }
    if(ichan<0||ichan==3) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,2,1,3,1,rho0,2,1,2,4));
      current -= bwa1[1]*bwrho[1]*epsilon(part.momentum(),pa1[1],momenta[0]-momenta[3]);
    }
    if(ichan<0||ichan==4) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rho0,2,3,2,2));
      current += bwa1[0]*bwrho[2]*epsilon(part.momentum(),pa1[0],momenta[2]-momenta[1]);
    }
    if(ichan<0||ichan==5) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,1,1,rho0,2,3,2,2));
      current -= bwa1[3]*bwrho[2]*epsilon(part.momentum(),pa1[3],momenta[2]-momenta[1]);
    }
    if(ichan<0||ichan==6) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,3,1,4,1,rho0,2,1,2,2));
      current += bwa1[2]*bwrho[0]*epsilon(part.momentum(),pa1[2],momenta[0]-momenta[1]);
    }
    if(ichan<0||ichan==7) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,4,1,3,1,rho0,2,1,2,2));
      current -= bwa1[3]*bwrho[0]*epsilon(part.momentum(),pa1[3],momenta[0]-momenta[1]);
    }
  }
  // decay mode f_1 -> pi+ pi0 pi- pi0
  else {
    complex<InvEnergy2> bwrho[4] =
      {Resonance::BreitWignerPWave((momenta[0]+momenta[1]).m2(),mrho_,grho_,momenta[0].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[0]+momenta[3]).m2(),mrho_,grho_,momenta[0].mass(),momenta[3].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[1]).m2(),mrho_,grho_,momenta[2].mass(),momenta[1].mass())/sqr(mrho_),
       Resonance::BreitWignerPWave((momenta[2]+momenta[3]).m2(),mrho_,grho_,momenta[2].mass(),momenta[3].mass())/sqr(mrho_)};
    double f1 = (momenta[2].mass2()-momenta[3].mass2())/sqr(mrho_);
    double f2 = 1+f1;
    f1 = 1.-f1;
    if(ichan<=0) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,2,1,rhop,2,3,2,4));
      current += bwa1[0]*bwrho[3]*epsilon(part.momentum(),pa1[0],f1*momenta[2]-f2*momenta[3]);
    }
    if(ichan<0||ichan==1) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,1,1,rhop,2,3,2,4));
      current -= bwa1[1]*bwrho[3]*epsilon(part.momentum(),pa1[1],f1*momenta[2]-f2*momenta[3]);
    }
    if(ichan<0||ichan==2) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1m,0,1,1,4,1,rhop,2,3,2,2));
      current += bwa1[0]*bwrho[2]*epsilon(part.momentum(),pa1[0],f1*momenta[2]-f2*momenta[1]);
    }
    if(ichan<0||ichan==3) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,1,1,rhop,2,3,2,2));
      current -= bwa1[3]*bwrho[2]*epsilon(part.momentum(),pa1[3],f1*momenta[2]-f2*momenta[1]);
    }
    if(ichan<0||ichan==4) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,4,1,rhom,2,1,2,2));
      current += bwa1[2]*bwrho[0]*epsilon(part.momentum(),pa1[2],f1*momenta[0]-f2*momenta[1]);
    }
    if(ichan<0||ichan==5) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,4,1,3,1,rhom,2,1,2,2));
      current -= bwa1[3]*bwrho[0]*epsilon(part.momentum(),pa1[3],f1*momenta[0]-f2*momenta[1]);
    }
    if(ichan<0||ichan==6) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a1p,0,3,1,2,1,rhom,2,1,2,4));
      current += bwa1[2]*bwrho[1]*epsilon(part.momentum(),pa1[2],f1*momenta[0]-f2*momenta[3]);
    }
    if(ichan<0||ichan==7) {
      // mode->addChannel((PhaseSpaceChannel(mode),0,a10,0,2,1,3,1,rhom,2,1,2,4));
      current -= bwa1[1]*bwrho[1]*epsilon(part.momentum(),pa1[1],f1*momenta[0]-f2*momenta[3]);
    }
  }
  // contract the current
  double pre = gRhoPiPi_*gf1a1Pi_*ga1RhoPi_;
  for(unsigned int ihel=0;ihel<3;++ihel)
    (*ME())(ihel,0,0,0,0) = Complex(pre*current.dot(vector_[ihel])*part.mass());
  // matrix element
  return sym*ME()->contract(rho_).real();
}

void f1FourPiDecayer::dataBaseOutput(ofstream & output,
					    bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DecayIntegrator base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":gRhoPiPi " << gRhoPiPi_     << "\n";
  output << "newdef " << name() << ":ga1RhoPi " << ga1RhoPi_/GeV << "\n";
  output << "newdef " << name() << ":gf1a1Pi "  << gf1a1Pi_*GeV  << "\n";
  for(unsigned int ix=0;ix<maxWeight_.size();++ix) {
    output << "newdef    " << name() << ":maxWeight " << ix << " "
  	   << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
  		    << fullName() << "\";" << endl;
}
