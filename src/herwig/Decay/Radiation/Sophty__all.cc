#line 1 "./SOPHTY.cc"
// -*- C++ -*-
//
// SOPHTY.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SOPHTY class.
//

#include "SOPHTY.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "FFDipole.h"
#include "IFDipole.h"

using namespace Herwig;

void SOPHTY::persistentOutput(PersistentOStream & os) const {
  os << FFDipole_ << IFDipole_ << colouredOption_;
}

void SOPHTY::persistentInput(PersistentIStream & is, int) {
  is >> FFDipole_ >> IFDipole_ >> colouredOption_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SOPHTY,DecayRadiationGenerator>
describeHerwigSOPHTY("Herwig::SOPHTY", "HwSOPHTY.so");

void SOPHTY::Init() {
  
  static ClassDocumentation<SOPHTY> documentation
    ("The SOPHTY class implements photon radiation in decays",
     "QED in particle decays was generated using the approach described in "
     "\\cite{Hamilton:2006xz}.",
     "\\bibitem{Hamilton:2006xz} K.~Hamilton and P.~Richardson,"
     "JHEP 07 (2006) 010.");
  
  static Reference<SOPHTY,FFDipole> interfaceFFDipole
    ("FFDipole",
     "The final-final dipole",
     &SOPHTY::FFDipole_, false, false, true, false, false);
  
  static Reference<SOPHTY,IFDipole> interfaceIFDipole
    ("IFDipole",
     "_ifdipole",
     &SOPHTY::IFDipole_, false, false, true, false, false);

  static Switch<SOPHTY,unsigned int> interfaceColouredTreatment
    ("ColouredTreatment",
     "Option for the treatment of QED radiation in decays involving coloured particles.",
     &SOPHTY::colouredOption_, 0, false, false);
  static SwitchOption interfaceColouredTreatmentNone
    (interfaceColouredTreatment,
     "None",
     "Generate no QED radiation to avoid problems with the interplay"
     " of QCD and QED radiation",
     0);
  static SwitchOption interfaceColouredTreatmentRadiation
    (interfaceColouredTreatment,
     "Radiation",
     "Generate radiation from the coloured particles.",
     1);

}

ParticleVector SOPHTY::generatePhotons(const Particle & p,ParticleVector children,
				       tDecayIntegratorPtr decayer) {
  if ( children.size() != 2 ) return children;
  // if not generating radiation from coloured particles
  // return if there are any coloured particles
  if(colouredOption_==0) {
    bool coloured = p.dataPtr()->coloured();
    for(unsigned int ix=0;ix<children.size();++ix) {
      coloured |= children[ix]->dataPtr()->coloured();
    }
    if(coloured) return children;
  }
  useMe();
  // final-final dipole
  if(p.dataPtr()->iCharge()==0) {
    if(children[0]->dataPtr()->iCharge()!=0&&
       children[1]->dataPtr()->iCharge()!=0)
      return FFDipole_->generatePhotons(p,children,decayer);
    else
      return children;
  }
  // initial final dipole
  else {
    if((children[0]->dataPtr()->iCharge()==0&&
	children[1]->dataPtr()->iCharge()!=0)||
       (children[0]->dataPtr()->iCharge()!=0&&
	children[1]->dataPtr()->iCharge()==0))
      return IFDipole_->generatePhotons(p,children);
    else
      return children;
  }
}
#line 1 "./FFDipole.cc"
// -*- C++ -*-
//
// FFDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFDipole class.
//

#include "FFDipole.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "YFSFormFactors.h"
#include "Herwig/Decay/PhaseSpaceMode.h"
#include "Herwig/Decay/DecayIntegrator.h"

using namespace Herwig;

void FFDipole::persistentOutput(PersistentOStream & os) const {
  os << ounit(_emin,GeV) << ounit(_eminrest,GeV) << ounit(_eminlab,GeV) 
     << _maxwgt << _weightOutput
     << _mode << _maxtry << _energyopt << _betaopt << _dipoleopt;
}

void FFDipole::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_emin,GeV) >> iunit(_eminrest,GeV) >> iunit(_eminlab,GeV) 
     >> _maxwgt >> _weightOutput
     >> _mode >> _maxtry >> _energyopt >> _betaopt >> _dipoleopt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFDipole,Interfaced>
describeHerwigFFDipole("Herwig::FFDipole", "HwSOPHTY.so");

void FFDipole::Init() {

  static ClassDocumentation<FFDipole> documentation
    ("The FFDipole class implements the final-final dipole for the SOPTHY algorithm");

  static Switch<FFDipole,unsigned int> interfaceUnWeight
    ("UnWeight",
     "Control the type of unweighting to perform, only one should be used the"
     " other options are for debugging purposes.",
     &FFDipole::_mode, 1, false, false);
  static SwitchOption interfaceUnWeightNoUnweighting
    (interfaceUnWeight,
     "NoUnweighting",
     "Perform no unweighting",
     0);
  static SwitchOption interfaceUnWeightAllWeights
    (interfaceUnWeight,
     "AllWeights",
     "Include all the weights",
     1);
  static SwitchOption interfaceUnWeightNoJacobian
    (interfaceUnWeight,
     "NoJacobian",
     "Only include the dipole and YFS weights",
     2);
  static SwitchOption interfaceUnWeightDipole
    (interfaceUnWeight,
     "Dipole",
     "Only include the dipole weight",
     3);
  static SwitchOption interfaceUnWeightYFS
    (interfaceUnWeight,
     "YFS",
     "Only include the YFS weight",
     4);
  static SwitchOption interfaceUnWeightNLO
    (interfaceUnWeight,
     "NLO",
     "Weight to get the stict NLO rate",
     5);

  static Parameter<FFDipole,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight",
     &FFDipole::_maxtry, 500, 10, 100000,
     false, false, Interface::limited);

  static Parameter<FFDipole,Energy> interfaceMinimumEnergyBoosted
    ("MinimumEnergyBoosted",
     "The minimum energy of the photons in the boosted frame in which"
     " they are generated.",
     &FFDipole::_emin, MeV, 1.e-6*MeV, ZERO, 100.0*MeV,
     false, false, Interface::limited);

  static Parameter<FFDipole,Energy> interfaceMinimumEnergyRest
    ("MinimumEnergyRest",
     "The minimum energy of the photons in the rest frame of the decaying particle",
     &FFDipole::_eminrest, MeV, 100.0*MeV, 1.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<FFDipole,Energy> interfaceMinimumEnergyLab
    ("MinimumEnergyLab",
     "The minimum energy of the photons in the lab frame",
     &FFDipole::_eminlab, MeV, 100.0*MeV, 1.0*MeV, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<FFDipole,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for unweighting",
     &FFDipole::_maxwgt, 7.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Switch<FFDipole,unsigned int> interfaceEnergyCutOff
    ("EnergyCutOff",
     "The type of cut-off on the photon energy to apply",
     &FFDipole::_energyopt, 1, false, false);
  static SwitchOption interfaceEnergyCutOffBoostedFrame
    (interfaceEnergyCutOff,
     "BoostedFrame",
     "Only apply cut-off in boosted frame",
     0);
  static SwitchOption interfaceEnergyCutOffRestFrame
    (interfaceEnergyCutOff,
     "RestFrame",
     "Apply cut-off in rest frame",
     1);
  static SwitchOption interfaceEnergyCutOff2
    (interfaceEnergyCutOff,
     "LabFrame",
     "Apply cut-off in lab frame",
     2);

  static Switch<FFDipole,unsigned int> interfaceBetaOption
    ("BetaOption",
     "Option for the inclusive of the higher beta coefficients",
     &FFDipole::_betaopt, 4, false, false);
  static SwitchOption interfaceBetaOptionNone
    (interfaceBetaOption,
     "None",
     "No higher betas included",
     0);
  static SwitchOption interfaceBetaOptionCollinear
    (interfaceBetaOption,
     "Collinear",
     "Include the collinear approx",
     1);
  static SwitchOption interfaceBetaOptionCollinearVirtA
    (interfaceBetaOption,
     "CollinearVirtualA",
     "Include the collinear approx with virtual corrections",
     2);
  static SwitchOption interfaceBetaOptionCollinearVirtB
    (interfaceBetaOption,
     "CollinearVirtualB",
     "Include the collinear approx with virtual corrections",
     3);
  static SwitchOption interfaceBetaOptionExact
    (interfaceBetaOption,
     "Exact",
     "Include the exact higher order terms if available",
     4);

  static Switch<FFDipole,unsigned int> interfaceDipoleOption
    ("DipoleOption",
     "Option for generating the primary dipole distribution",
     &FFDipole::_dipoleopt, 0, false, false);
  static SwitchOption interfaceDipoleOptionNoMass
    (interfaceDipoleOption,
     "NoMass",
     "Don't include the mass terms in the primary distribution",
     0);
  static SwitchOption interfaceDipoleOptionMass
    (interfaceDipoleOption,
     "Mass",
     "Include the mass terms in the primary distribution",
     1);

  static Switch<FFDipole,bool> interfaceWeightOutput
    ("WeightOutput",
     "Whether or not to output the average weight for testing",
     &FFDipole::_weightOutput, false, false, false);
  static SwitchOption interfaceWeightOutputNo
    (interfaceWeightOutput,
     "No",
     "Don't output the average",
     false);
  static SwitchOption interfaceWeightOutputYes
    (interfaceWeightOutput,
     "Yes",
     "Output the average",
     true);

}

void FFDipole::printDebugInfo(const Particle & p,
			      const ParticleVector & children,
			      double wgt) const {
  generator()->log() << "Input masses " 
		     << p.mass()/GeV << " -> " 
		     << children[0]->mass()/GeV << " " 
		     << children[1]->mass()/GeV << '\n'; 
  generator()->log() << "Momenta\n";
  generator()->log() << "parent " << p.momentum()/GeV << '\n';
  for(unsigned int ix=0;ix<2;++ix)
    generator()->log() << "charged " << ix << " " 
		       << _qnewlab[ix]/GeV << " " 
		       << children[ix]->momentum()/GeV << '\n';
  for(unsigned int ix=0;ix<_multiplicity;++ix) {
    generator()->log() << "photons " << ix << " "
		       << "phocut " << _photcut[ix] << ' '
		       << _llab[ix]/GeV << '\n';
  }
  generator()->log() << "wgt         : " << wgt          << '\n';
  generator()->log() << "_mewgt      : " << _mewgt       << '\n';
  generator()->log() << "_jacobianwgt: " << _jacobianwgt << '\n';
  generator()->log() << "_yfswgt     : " << _yfswgt      << '\n';
  generator()->log() << "_dipolewgt  : " << _dipolewgt   << '\n';
  generator()->log() << "dipoleopt   : " << _dipoleopt   << '\n';
}

ParticleVector FFDipole::generatePhotons(const Particle & p,
					 ParticleVector children,
					 tDecayIntegratorPtr decayer) {
  _parent = const_ptr_cast<tPPtr>(&p);
  // set the decayer
  _decayer=decayer;
  // set parameters which won't change in the event loop
  // masses of the particles
  _m[0] = p.mass(); 
  _m[1] = children[0]->mass();
  _m[2] = children[1]->mass();
  // set the maximum photon energy (exact - no approximations here).
  _emax=(0.5*(_m[0]-sqr(_m[1]+_m[2])/_m[0]))*_m[0]/(_m[1]+_m[2]);
  // check masses non-zero
  for(unsigned int ix=0;ix<2;++ix) {
    if(children[ix]->mass()<1e-4*GeV) { 
      ostringstream message;
      message << "FFDipole::generatePhotons() trying to generate QED radiation from "
	      << children[ix]->dataPtr()->PDGName() << "\n with mass " << children[ix]->mass()/GeV
	      << "which is much smaller than the mass of the electron.\n"
	      << "This is probably due to reading events from a LHEF,\nskipping radiation in this case.\n";
      generator()->logWarning( Exception(message.str(), Exception::warning));
      return children;
    }
  }
  // momenta before radiation in lab
  for(unsigned int ix=0;ix<2;++ix)
    _qlab[ix]=children[ix]->momentum();
  // get the charges of the particles in units of the positron charge
  _charge=children[0]->dataPtr()->iCharge()*children[1]->dataPtr()->iCharge()/9.;
  // boost the momenta to the rest frame
  Boost boostv(-p.momentum().boostVector());
  // boost the particles to the parent rest frame
  // and set the initial momenta of the charged particles 
  // in the dipole rest frame: currently this is the same 
  // as the boson rest frame...
  for(unsigned int ix=0;ix<2;++ix) {
    children[ix]->deepBoost(boostv);
    _qdrf[ix]=children[ix]->momentum();
    _qprf[ix]=children[ix]->momentum();
  }
  _parent->boost(boostv);
  // perform the unweighting
  double wgt;
  unsigned int ntry(0);
  do {
    ++ntry;
    wgt = makePhotons(-boostv,children);

    // Error checks
    if ( std::isnan(wgt) ) {
      generator()->log() << "Infinite weight for decay " 
			 << p.PDGName() << " " 
			 << children[0]->PDGName() 
			 << " " << children[1]->PDGName()
			 << '\n';
      wgt = 0.0;
    }
    else if ( wgt < 0.0 && _mode != 5 ) {
      generator()->log() << "Negative weight for decay " 
			 << p.PDGName() << " " 
			 << children[0]->PDGName() 
			 << " " << children[1]->PDGName()
			 << "in FFDipole: Weight = " << wgt << '\n';
      if ( Debug::level ) 
	printDebugInfo(p,children,wgt);
    }
    else if ( wgt > _maxwgt ) {
      generator()->log() << "Weight "<< wgt<<" exceeds maximum for decay " 
			 << p.PDGName() << ' '
			 << children[0]->PDGName() 
			 << " " << children[1]->PDGName()
			 << " in FFDipole:\nresetting maximum weight.\n"
			 << "Old Maximum = " << _maxwgt;

      _maxwgt = min(1.1 * wgt, 10.0);
      generator()->log() << " New Maximum = " << wgt << '\n';
      if ( Debug::level && _mode!=5 ) 
	printDebugInfo(p,children,wgt);
    }
    // End of error checks

    _wgtsum += wgt;
    _wgtsq  += sqr(wgt);
    ++_nweight;

  }
  while ( wgt<(_maxwgt*UseRandom::rnd()) && ntry<_maxtry );
  if(ntry>=_maxtry) {
    generator()->log() << "FFDipole failed to generate QED radiation for the decay " 
		       << p.PDGName() << " -> " 
		       << children[0]->PDGName() << " "
		       << children[1]->PDGName() << '\n';
    _parent->boost(-boostv);
    for(unsigned int ix=0;ix<2;++ix)
      children[ix]->deepBoost(-boostv);
    return children;
  }
  // produce products after radiation if needed
  if(_multiplicity>0) {
    // change the momenta of the children, they are currently
    // in original rest frame 
    for(unsigned int ix=0;ix<2;++ix) {
      // unit vector along direction
      Boost br = children[ix]->momentum().vect().unit();
      // calculate the boost vector using expression accurate for beta->1
      double beta(sqrt((_qdrf[ix].e()+_m[ix+1])*(_qdrf[ix].e()-_m[ix+1]))/
		  _qdrf[ix].e());
      double ombeta(sqr(_m[ix+1]/_qdrf[ix].e())/(1.+beta));
      double betap(sqrt((_qnewdrf[ix].e()+_m[ix+1])*(_qnewdrf[ix].e()-_m[ix+1]))
		   /_qnewdrf[ix].e());
      double ombetap(sqr(_m[ix+1]/_qnewdrf[ix].e())/(1.+betap));
      // boost to get correct momentum in dipole rest frame
      double bv = -(ombetap-ombeta)/(beta*ombetap + ombeta);
      br *= bv;
      children[ix]->deepBoost(br);
      // boost to the parent rest frame
      Lorentz5Momentum pnew(_bigLdrf);
      pnew.setMass(_m[0]);
      pnew.rescaleEnergy();
      br = pnew.findBoostToCM();
      children[ix]->deepBoost(br);
      // boost back to the lab
      children[ix]->deepBoost(-boostv);
    }
    // add the photons to the event record
    tcPDPtr photon=getParticleData(ParticleID::gamma);
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      // add if not removed because energy too low
      if(!_photcut[ix]) {
	PPtr newphoton=new_ptr(Particle(photon));
	newphoton->set5Momentum(_llab[ix]);
	children.push_back(newphoton);
      }
    }
    _parent->boost(-boostv);

    //printDebugInfo(p, children, wgt);

    return children;
  }
  // otherwise just return the original particles
  else {
    for(unsigned int ix=0;ix<2;++ix)
      children[ix]->deepBoost(-boostv);
    _parent->boost(-boostv);
    return children;
  }
}

// member which generates the photons
double FFDipole::makePhotons(const Boost & boostv,
			     const ParticleVector & children) {
  // set the initial parameters
  // number of photons (zero)
  _multiplicity=0;
  // zero size of photon vectors
  _ldrf.clear();
  _lprf.clear();
  _llab.clear();
  // zero size of angle storage
  _sinphot.clear();
  _cosphot.clear();
  _photcut.clear();
  _photonwgt.clear();
  // zero total momenta of the photons
  _bigLdrf=Lorentz5Momentum();
  _bigLprf=Lorentz5Momentum();
  // set the initial values of the reweighting factors to one
  _dipolewgt   = 1.0;
  _yfswgt      = 1.0;
  _jacobianwgt = 1.0;
  _mewgt       = 1.0;
  // calculate the velocities of the charged particles (crude/overvalued)
  double beta1(sqrt((_qdrf[0].e()+_m[1])*(_qdrf[0].e()-_m[1]))/_qdrf[0].e());
  double beta2(sqrt((_qdrf[1].e()+_m[2])*(_qdrf[1].e()-_m[2]))/_qdrf[1].e());
  // calculate 1-beta to avoid numerical problems
  double ombeta1(sqr(_m[1]/_qdrf[0].e())/(1.+beta1));
  double ombeta2(sqr(_m[2]/_qdrf[1].e())/(1.+beta2));
  // calculate the average photon multiplicity
  double aver(YFSFormFactors::nbarFF(beta1,ombeta1,beta2,ombeta2,_charge,
				     _emax,_emin,_dipoleopt==1));
  // calculate the number of photons using the poisson
  _multiplicity = _mode !=5 ? UseRandom::rndPoisson(aver) : 1;
  // calculate the first part of the YFS factor
  // (N.B. crude form factor is just exp(-aver) to get a poisson)
  _yfswgt *= exp(aver);
  // if photons produced
  if(_multiplicity>0) {
    _photonwgt.resize(_multiplicity);
    // generate the photon momenta with respect to q1 
    // keeping track of the weight
    for(unsigned int ix=0;ix<_multiplicity;++ix)
      _photonwgt[ix] = photon(beta1,ombeta1,beta2,ombeta2);
    // rotate the photons so in dipole rest frame rather 
    // than angle measured w.r.t q1 first work out the rotation
    SpinOneLorentzRotation rotation;
    rotation.setRotateZ(-_qdrf[0].phi());
    rotation.rotateY(_qdrf[0].theta());
    rotation.rotateZ(_qdrf[0].phi());
    // rotate the total
    _bigLdrf *= rotation;
    // rotate the photons
    for(unsigned int ix=0;ix<_multiplicity;++ix)
      _ldrf[ix]*=rotation;
    // boost the momenta without any removal of low energy photons
    // resize arrays
    _photcut.resize(_multiplicity,false);
    _lprf.resize(_multiplicity);
    _llab.resize(_multiplicity);
    // perform the boost
    if(!boostMomenta(boostv)){return 0.;}
    // apply the cut on the photon energy if needed
    unsigned int nremoved(removePhotons());
    // redo the boost if we have removed photons
    if(nremoved!=0){if(!boostMomenta(boostv)){return 0.;}}
    // form factor part of the removal term to remove existing cut
    if(_energyopt!=0) _dipolewgt *= 
      YFSFormFactors::exponentialYFSFF(beta1,ombeta1,beta2,ombeta2,
				       _qdrf[0].e(),_qdrf[1].e(),
				       _m[1],_m[2],_m[0]*_m[0],
				       _charge,_emin);
    // calculate the new dipole weight
    // calculate velocities and 1-velocites
    beta1=sqrt((_qnewdrf[0].e()+_m[1])*(_qnewdrf[0].e()-_m[1]))/_qnewdrf[0].e();
    beta2=sqrt((_qnewdrf[1].e()+_m[2])*(_qnewdrf[1].e()-_m[2]))/_qnewdrf[1].e();
    ombeta1=sqr(_m[1]/_qnewdrf[0].e())/(1.+beta1);
    ombeta2=sqr(_m[2]/_qnewdrf[1].e())/(1.+beta2);
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(!_photcut[ix])
	_dipolewgt *= exactDipoleWeight(beta1,ombeta1,beta2,ombeta2,ix)/
	  _photonwgt[ix];
    }
    // calculate the weight for the photon removal
    Energy2 s((_qnewdrf[0]+_qnewdrf[1]).m2());
    // calculate the second part of the yfs form factor
    // this is different for the different photon removal options
    // option with no removal
    if(_energyopt==0) {
      _yfswgt *= 
	YFSFormFactors::exponentialYFSFF(beta1,ombeta1,beta2,ombeta2,
					 _qnewdrf[0].e(),_qnewdrf[1].e(),
					 _m[1],_m[2],s,_charge,_emin);
    }
    // weight for option with cut in the rest frame
    else if(_energyopt==1) {
      // yfs piece
      double nbeta1(sqrt( (_qnewprf[0].e()+_m[1])*(_qnewprf[0].e()-_m[1]))
		    /_qnewprf[0].e());
      double nbeta2(sqrt( (_qnewprf[1].e()+_m[2])*(_qnewprf[1].e()-_m[2]))
		    /_qnewprf[1].e());
      double nomb1 (sqr(_m[1]/_qnewprf[0].e())/(1.+nbeta1));
      double nomb2 (sqr(_m[2]/_qnewprf[1].e())/(1.+nbeta2));
      _yfswgt *= 
	YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
					 _qnewprf[0].e(),_qnewprf[1].e(),
					 _m[1],_m[2],s,_charge,_eminrest);
      // dipole piece
      // Find the momenta of the particles of original particles in new rest frame 
      Lorentz5Momentum pnew(_bigLdrf.x(),_bigLdrf.y(),
			    _bigLdrf.z(),_bigLdrf.e(),_m[0]);
      pnew.rescaleEnergy();
      SpinOneLorentzRotation boost(pnew.findBoostToCM());
      Lorentz5Momentum q1=boost*_qdrf[0];
      Lorentz5Momentum q2=boost*_qdrf[1];
      // use this to calculate the form factor
      nbeta1=sqrt( (q1.e()+_m[1])*(q1.e()-_m[1]))/q1.e();
      nbeta2=sqrt( (q2.e()+_m[2])*(q2.e()-_m[2]))/q2.e();
      nomb1 =sqr(_m[1]/q1.e())/(1.+nbeta1);
      nomb2 =sqr(_m[2]/q2.e())/(1.+nbeta2);
      _dipolewgt /=YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
						    q1.e(),q2.e(),
						    _m[1],_m[2],_m[0]*_m[0],
						    _charge,_eminrest);
    }
    // weight for option with cut in the rest frame
    else if(_energyopt==2) {
      // yfs piece
      double nbeta1(sqrt( (_qnewlab[0].e()+_m[1])*(_qnewlab[0].e()-_m[1]))
		    /_qnewlab[0].e());
      double nbeta2(sqrt( (_qnewlab[1].e()+_m[2])*(_qnewlab[1].e()-_m[2]))
		    /_qnewlab[1].e());
      double nomb1 (sqr(_m[1]/_qnewlab[0].e())/(1.+nbeta1));
      double nomb2 (sqr(_m[2]/_qnewlab[1].e())/(1.+nbeta2));
      _yfswgt *= 
	YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
					 _qnewlab[0].e(),_qnewlab[1].e(),
					 _m[1],_m[2],s,_charge,_eminlab);
      // dipole piece
      // Find the momenta of the particles of original particles in new rest frame 
      Lorentz5Momentum pnew(_bigLdrf.x(),_bigLdrf.y(),
			    _bigLdrf.z(),_bigLdrf.e(),_m[0]);
      pnew.rescaleEnergy();
      SpinOneLorentzRotation boost(pnew.findBoostToCM());
      Lorentz5Momentum q1=boost*_qdrf[0];
      Lorentz5Momentum q2=boost*_qdrf[1];
      // then boost to the lab
      boost.setBoost(boostv);
      q1 *=boost;
      q2 *=boost;
      // use this to calculate the form factor
      nbeta1=sqrt( (q1.e()+_m[1])*(q1.e()-_m[1]))
	/q1.e();
      nbeta2=sqrt( (q2.e()+_m[2])*(q2.e()-_m[2]))
	/q2.e();
      nomb1 =sqr(_m[1]/q1.e())/(1.+nbeta1);
      nomb2 =sqr(_m[2]/q2.e())/(1.+nbeta2);
      _dipolewgt /=YFSFormFactors::exponentialYFSFF(nbeta1,nomb1,nbeta2,nomb2,
						    q1.e(),q2.e(),_m[1],_m[2],
						    _m[0]*_m[0],_charge,_eminlab);
    }
    // Calculating jacobian weight
    _jacobianwgt = jacobianWeight();
    // Calculate the weight for the corrections
    _mewgt = meWeight(children);
  }
  // otherwise copy momenta
  else {
    for(unsigned int ix=0;ix<2;++ix) {
      _qnewdrf[ix]=_qdrf[ix];
      _qnewprf[ix]=_qprf[ix]; 
      _qnewlab[ix]=_qlab[ix]; 
    }
    _jacobianwgt = 1.0;
    _yfswgt*=YFSFormFactors::exponentialYFSFF(beta1,ombeta1,beta2,ombeta2,
					      _qdrf[0].e(),_qdrf[1].e(),
					      _m[1],_m[2],_m[0]*_m[0],
					      _charge,_emin);
    _dipolewgt   = 1.0;
  }
  double wgt;
  if(_mode!=5) {
    // virtual corrections
    _mewgt += virtualWeight(children);
    // calculate the weight depending on the option
    if(_mode==0)      wgt = _maxwgt;
    else if(_mode==1) wgt = _mewgt*_yfswgt*_jacobianwgt*_dipolewgt;
    else if(_mode==2) wgt = _jacobianwgt*_yfswgt*_dipolewgt;
    else if(_mode==3) wgt = _yfswgt*_dipolewgt;
    else              wgt = _yfswgt;
  }
  // special to test NLO results
  else {
    double beta1   = sqrt((_qdrf[0].e()+_m[1])*(_qdrf[0].e()-_m[1]))/_qdrf[0].e();
    double beta2   = sqrt((_qdrf[1].e()+_m[2])*(_qdrf[1].e()-_m[2]))/_qdrf[1].e();
    double ombeta1 = sqr(_m[1]/_qdrf[0].e())/(1.+beta1);
    double ombeta2 = sqr(_m[2]/_qdrf[1].e())/(1.+beta2);
    double yfs = YFSFormFactors::YFSFF(beta1,ombeta1,beta2,ombeta2,
				       _qdrf[0].e(),_qdrf[1].e(),
				       _m[1],_m[2],_m[0]*_m[0],
				       _charge,_emin);
    double nbar = YFSFormFactors::nbarFF(beta1,ombeta1,beta2,ombeta2,_charge,
					 _emax,_emin,_dipoleopt==1);
    wgt = 1.+virtualWeight(children)+yfs+nbar*_dipolewgt*_mewgt*_jacobianwgt;
  }
  return wgt;
}

double FFDipole::photon(double beta1,double ombeta1,
                        double beta2,double ombeta2) {
  // generate the polar angle
  double r1,r2,costh,sinth,opbc,ombc;
  // relative weights for the two terms
  double Pp(log((1+beta2)/ombeta2));
  double Pm(log((1+beta1)/ombeta1));
  Pp/=(Pp+Pm);
  // generate the angle
  double wgt=1.;
  do {
    r1=UseRandom::rnd();
    r2=UseRandom::rnd();
    // 1/(1+bc) branch
    if(r1<=Pp) {
      opbc  = pow(1.+beta2,r2)*pow(ombeta2,1.-r2);
      costh = -1./beta2*(1.-opbc);
      ombc  = 1.-beta1*costh;
      sinth = sqrt(opbc*(2.-opbc)-(1.+beta2)*ombeta2*sqr(costh));
    }
    // 1/(1-bc) branch
    else {
      ombc  = pow(1.+beta1,1.-r2)*pow(ombeta1,r2);
      costh = 1./beta1*(1.-ombc);
      opbc  = 1.+beta2*costh;
      sinth = sqrt(ombc*(2.-ombc)-(1.+beta1)*ombeta1*sqr(costh));
    }
    // wgt for rejection
    if(_dipoleopt==1)
      wgt = 1.-0.5/(1.+beta1*beta2)*(ombeta1*(1.+beta1)*opbc/ombc+
				     ombeta2*(1.+beta2)*ombc/opbc);
  }
  while(UseRandom::rnd()>wgt);
  // generate the polar angle randomly in -pi->+pi
  double phi(-pi+UseRandom::rnd()*2.*pi);
  // generate the ln(energy) uniformly in ln(_emin)->ln(_emax)
  Energy en(pow(_emax/_emin,UseRandom::rnd())*_emin);
  // calculate the weight (omit the pre and energy factors
  //                       which would cancel later anyway)
  if(_dipoleopt==0)
    wgt = 0.5*(1.+beta1*beta2)/opbc/ombc;
  else
    wgt = 0.25*(2.*(1.+beta1*beta2)/opbc/ombc
		-ombeta1*(1.+beta1)/sqr(ombc)
		-ombeta2*(1.+beta2)/sqr(opbc));
  // store the angles
  _cosphot.push_back(costh);
  _sinphot.push_back(sinth);
  // store the four vector for the photon
  _ldrf.push_back(Lorentz5Momentum(en*sinth*cos(phi),
				   en*sinth*sin(phi),
				   en*costh,en,
				   ZERO));
  // add the photon momentum to the total
  _bigLdrf+=_ldrf.back();
  // return the weight
  return wgt;
}

double FFDipole::meWeight(const ParticleVector & children) {
  if(_multiplicity==0) return 1.;
  // option which does nothing
  if(_betaopt==0) {
    return 1.;
  }
  // collinear approx
  else if(_betaopt <= 3) {
    return collinearWeight(children);
  }
  else if (_betaopt == 4 ) {
    if(_decayer&&_decayer->hasRealEmissionME()) {
      double outwgt=1.;
      // values of beta etc to evaluate the dipole
      double beta1(sqrt( (_qnewdrf[0].e()+_m[1])*(_qnewdrf[0].e()-_m[1]))/
 		   _qnewdrf[0].e());
      double beta2(sqrt( (_qnewdrf[1].e()+_m[2])*(_qnewdrf[1].e()-_m[2]))/
 		   _qnewdrf[1].e());
      double ombeta1(sqr(_m[1]/_qnewdrf[0].e())/(1.+beta1));
      double ombeta2(sqr(_m[2]/_qnewdrf[1].e())/(1.+beta2));
      // storage of the weights
      ParticleVector ptemp;
      for(unsigned int ix=0;ix<children.size();++ix)
	ptemp.push_back(new_ptr(Particle(children[ix]->dataPtr())));
      ptemp.push_back(new_ptr(Particle(getParticleData(ParticleID::gamma))));
      for(unsigned int i=0;i<_multiplicity;++i) {
	PPtr new_parent = new_ptr(Particle(*_parent));
	if(_photcut[i]) continue;
	// compute the angle terms
	// if cos is greater than zero use result accurate as cos->1
	double opbc,ombc;
	if(_cosphot[i]>0) {
	  opbc=1.+beta2*_cosphot[i];
	  ombc=ombeta1+beta1*sqr(_sinphot[i])/(1.+_cosphot[i]);
	}
	// if cos is less    than zero use result accurate as cos->-1
	else {
	  opbc=ombeta2+beta2*sqr(_sinphot[i])/(1.-_cosphot[i]);
	  ombc=1.-beta1*_cosphot[i];
	}
	// dipole factor for denominator
	double dipole = 2./opbc/ombc*(1.+beta1*beta2
				      -0.5*ombeta1*(1.+beta1)*opbc/ombc		 
				      -0.5*ombeta2*(1.+beta2)*ombc/opbc); 
	// energy and momentum of the photon
	Energy L0(_ldrf[i].e()),modL(_ldrf[i].rho());
	// 3-momenta of charged particles
	Energy modq(_qdrf[0].rho());
	// calculate the energy of the fermion pair
	Energy newE12(-L0+sqrt(sqr(_m[0])+sqr(modL)));
	// 3-momentum rescaling factor (NOT energy rescaling).
	double kappa(Kinematics::pstarTwoBodyDecay(newE12,_m[1],_m[2])/modq);
	// calculate the rescaled momenta
	Lorentz5Momentum porig[3];
	for(unsigned int ix=0;ix<2;++ix) {
	  porig[ix] = kappa*_qdrf[ix];
	  porig[ix].setMass(_m[ix+1]);
	  porig[ix].rescaleEnergy();
	}
	porig[2] = _ldrf[i];
	// calculate the momentum of the decaying particle in dipole rest frame
	Lorentz5Momentum pnew(_ldrf[i].x(),_ldrf[i].y(),
			      _ldrf[i].z(),_ldrf[i].e(),_m[0]);
	pnew.rescaleEnergy();
	// Find the momenta of the particles in the rest frame of the parent...
	// First get the boost from the parent particle
	Boost boost = pnew.findBoostToCM();
	LorentzRotation rot1(-boost, pnew.e()/pnew.mass());
	// check the photon energy
	Lorentz5Momentum ptest = _ldrf[i];
	ptest.boost(boost);
	if(_energyopt==1&&ptest.e()<_eminrest) continue;
	new_parent->transform(rot1);
	// rotation to put the emitter along the z axis
	// first particle emits
	unsigned int iemit = _cosphot[i]>0. ? 0 : 1;
	LorentzRotation rot2;
	rot2.setRotateZ(-porig[iemit].phi());
	rot2.rotateY(porig[iemit].theta());
	rot2.rotateZ(porig[iemit].phi());
	rot2.invert();
	// Boost the momenta of the charged particles
	for(unsigned int ix=0;ix<3;++ix) {
 	  porig[ix].transform(rot2);
	  ptemp[ix]->set5Momentum(porig[ix]);
	}
 	new_parent->transform(rot2);
	if(_cosphot[i]>0.) {
	  outwgt -= _decayer->
	    realEmissionME(_decayer->imode(),*new_parent,ptemp,
			   0,_cosphot[i],_sinphot[i],rot1,rot2)/
	    (_charge/sqr(_ldrf[i].e())*dipole);
	}
	else {
	  outwgt -= _decayer->
	    realEmissionME(_decayer->imode(),*new_parent,ptemp,
			   1,-_cosphot[i],_sinphot[i],rot1,rot2)/
	    (_charge/sqr(_ldrf[i].e())*dipole);
	}
	rot1.invert();
 	rot2.invert();
 	new_parent->transform(rot2);
 	new_parent->transform(rot1);
      }
      return outwgt;
    }
    else
      return collinearWeight(children);
  }
  return 1.;
}

double FFDipole::collinearWeight(const ParticleVector & children) {
  double outwgt=1.;
  // spins of the decay products
  PDT::Spin spin1(children[0]->dataPtr()->iSpin());
  PDT::Spin spin2(children[1]->dataPtr()->iSpin());
  // values of beta etc to evaluate the dipole
  double beta1(sqrt( (_qnewdrf[0].e()+_m[1])*(_qnewdrf[0].e()-_m[1]))/
	       _qnewdrf[0].e());
  double beta2(sqrt( (_qnewdrf[1].e()+_m[2])*(_qnewdrf[1].e()-_m[2]))/
	       _qnewdrf[1].e());
  double ombeta1(sqr(_m[1]/_qnewdrf[0].e())/(1.+beta1));
  double ombeta2(sqr(_m[2]/_qnewdrf[1].e())/(1.+beta2));
  // storage of the weights
  double twgt,dipole;
  double opbc,ombc;
  // compute the collinear approx
  for(unsigned int i=0;i<_multiplicity;++i) {
    if(_photcut[i]) continue;
    // compute the angle terms
    // if cos is greater than zero use result accurate as cos->1
    if(_cosphot[i]>0) {
      opbc=1.+beta2*_cosphot[i];
      ombc=ombeta1+beta1*sqr(_sinphot[i])/(1.+_cosphot[i]);
    }
    // if cos is less    than zero use result accurate as cos->-1
    else {
      opbc=ombeta2+beta2*sqr(_sinphot[i])/(1.-_cosphot[i]);
      ombc=1.-beta1*_cosphot[i];
    }
    // dipole factor for denominator
    dipole = 2.*(1.+beta1*beta2
		 -0.5*ombeta1*(1.+beta1)*opbc/ombc		 
		 -0.5*ombeta2*(1.+beta2)*ombc/opbc); 
    twgt=0.;
    // correction for the first particle
    double ratio(_ldrf[i].e()/_qnewdrf[0].e());
    if(spin1==PDT::Spin0)      twgt += 0.;
    else if(spin1==PDT::Spin1Half) 
      twgt += opbc*ratio/(1.+(1.+beta1*beta2)/ratio/opbc);
    else              
      twgt += 2.*sqr(opbc*ratio) *
	(+1./(1+beta1*beta2+_ldrf[i].e()/_qnewdrf[1].e()*ombc)
	 +(1.+beta1*beta2)/sqr(1.+beta1*beta2
			       +_ldrf[i].e()/_qnewdrf[0].e()*opbc));
    // correction for the second particle
    ratio =_ldrf[i].e()/_qnewdrf[1].e();
    if(spin2==PDT::Spin0)      twgt += 0.;
    else if(spin2==PDT::Spin1Half) 
      twgt += ombc*ratio/(1.+(1.+beta1*beta2)/ratio/ombc);
    else       
      twgt += 2.*sqr(ombc*ratio) *
	(1./(1. + beta1*beta2 + _ldrf[i].e()/_qnewdrf[0].e()*opbc)
	 + (1.+beta1*beta2) / sqr(1. + beta1*beta2
				  + _ldrf[i].e()/_qnewdrf[1].e()*ombc));
    twgt/=dipole;
    outwgt+=twgt;
  }
  return outwgt;
}

bool FFDipole::boostMomenta(const Boost & boostv) {
  // total energy  and momentum of photons
  Energy L0(_bigLdrf.e()),modL(_bigLdrf.rho());
  // 3-momenta of charged particles
  Energy modq(_qdrf[0].rho());
  // calculate the energy of the fermion pair
  Energy newE12(-L0+sqrt(_m[0]*_m[0]+modL*modL));
  // check this is allowed
  if(newE12<_m[1]+_m[2]){return false;}
  // 3-momentum rescaling factor (NOT energy rescaling).
  double kappa(Kinematics::pstarTwoBodyDecay(newE12,_m[1],_m[2])/modq);
  // calculate the rescaled momenta
  for(unsigned int ix=0;ix<2;++ix) {
    _qnewdrf[ix] = kappa*_qdrf[ix];
    _qnewdrf[ix].setMass(_m[ix+1]);
    _qnewdrf[ix].rescaleEnergy();
  }
  // calculate the momentum of the decaying particle in dipole rest frame
  Lorentz5Momentum pnew(_bigLdrf.x(),_bigLdrf.y(),
			_bigLdrf.z(),_bigLdrf.e(),_m[0]);
  pnew.rescaleEnergy();
  // Find the momenta of the particles in the rest frame 
  // of the parent...
  // First get the boost from the parent particle
  SpinOneLorentzRotation boost(pnew.findBoostToCM());
  // Boost the momenta of the charged particles
  for(unsigned int ix=0;ix<2;++ix) _qnewprf[ix]=boost*_qnewdrf[ix];
  // Boost the total photon momentum
  _bigLprf=boost*_bigLdrf;
  // Boost the individual photon momenta
  for(unsigned int ix=0;ix<_multiplicity;++ix){_lprf[ix]=boost*_ldrf[ix];}
  // Now boost from the parent rest frame to the lab frame
  boost.setBoost(boostv);
  // Boosting charged particles
  for(unsigned int ix=0;ix<2;++ix){_qnewlab[ix]=boost*_qnewprf[ix];}
  // Boosting total photon momentum
  _bigLlab=boost*_bigLprf;
  // Boosting individual photon momenta
  for(unsigned int ix=0;ix<_multiplicity;++ix){_llab[ix]=boost*_lprf[ix];}
  return true;
}

unsigned int FFDipole::removePhotons() {
  unsigned int nremoved(0);
  // apply the cut in the rest frame
  if(_energyopt==1) {
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(_lprf[ix].e()<_eminrest) {
	++nremoved;
	_photcut[ix]=true;
	_bigLdrf-=_ldrf[ix];
	_ldrf[ix]=Lorentz5Momentum();
      }
    }
  }
  // apply the cut in the lab frame
  else if(_energyopt==2) {
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(_llab[ix].e()<_eminlab) {
	++nremoved;
	_photcut[ix]=true;
	_bigLdrf-=_ldrf[ix];
	_ldrf[ix]=Lorentz5Momentum();
      }
    }
  }
  // correction factor for dipoles if needed
  if(_dipoleopt==0&&nremoved!=0) {
    // calculate the velocities of the charged particles (crude/overvalued)
    double beta1(sqrt((_qdrf[0].e()+_m[1])*(_qdrf[0].e()-_m[1]))/_qdrf[0].e());
    double beta2(sqrt((_qdrf[1].e()+_m[2])*(_qdrf[1].e()-_m[2]))/_qdrf[1].e());
    // calculate 1-beta to avoid numerical problems
    double ombeta1(sqr(_m[1]/_qdrf[0].e())/(1.+beta1));
    double ombeta2(sqr(_m[2]/_qdrf[1].e())/(1.+beta2));
    // calculate the weights
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      if(_photcut[ix]) _dipolewgt *= 
	exactDipoleWeight(beta1,ombeta1,beta2,ombeta2,ix)/_photonwgt[ix];
    }
  }
  // return number of remove photons
  return nremoved;
}

double FFDipole::virtualWeight(const ParticleVector & children) {
  double output = 0.;
  // Virtual corrections for beta_0:
  // These should be zero for the scalar case as there is no
  // collinear singularity going by the dipoles above...
  // Use mass of decaying particle...
  if(_betaopt==2) { 
    if((children[0]->dataPtr()->iSpin())==2&&
       (children[1]->dataPtr()->iSpin())==2
       ) {
      output += (1.0*YFSFormFactors::_alpha/pi)
	* log(sqr(_m[0]/_m[1]));
    }
  }
  // OR Use invariant mass of final state children...
  else if(_betaopt==3) { 
    if((children[0]->dataPtr()->iSpin())==2&&
       (children[1]->dataPtr()->iSpin())==2
       ) {
      output += (1.0*YFSFormFactors::_alpha/pi)
	* log((_qnewprf[0]+_qnewprf[1]).m2()/sqr(_m[1]));
    }
  }
  else if (_betaopt==4) {
    if(_decayer&&_decayer->hasOneLoopME()) {
      output += 
	_decayer->oneLoopVirtualME(_decayer->imode(),*_parent,
				   children);
    }
    else {
      output += (1.0*YFSFormFactors::_alpha/pi)
	* log(sqr(_m[0]/_m[1]));
    }
  }
  return output;
}

void FFDipole::dofinish() {
  Interfaced::dofinish();
  if(_weightOutput) {
    _wgtsum /= double(_nweight);
    _wgtsq  /= double(_nweight);
    _wgtsq = max(_wgtsq - sqr(_wgtsum),0.);
    _wgtsq /= double(_nweight);
    _wgtsq = sqrt(_wgtsq);
    generator()->log() << "The average weight for QED Radiation in " << fullName() 
		       << " was " << _wgtsum << " +/- " << _wgtsq << '\n';
  }
}
#line 1 "./IFDipole.cc"
// -*- C++ -*-
//
// IFDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFDipole class.
//

#include "IFDipole.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"

using namespace ThePEG::Helicity;

using namespace Herwig;

void IFDipole::persistentOutput(PersistentOStream & os) const {
  os << _alpha << ounit(_emin,GeV) << _maxwgt
     << _mode  << _maxtry << _energyopt  << _betaopt;
}

void IFDipole::persistentInput(PersistentIStream & is, int) {
  is >> _alpha >> iunit(_emin,GeV) >> _maxwgt
     >> _mode  >> _maxtry >> _energyopt  >> _betaopt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<IFDipole,Interfaced>
describeHerwigIFDipole("Herwig::IFDipole", "Herwig.so");

void IFDipole::Init() {
  static ClassDocumentation<IFDipole> documentation
    ("The IFDipole class implements the initial-final dipole for the SOPTHY algorithm");

  static Switch<IFDipole,unsigned int> interfaceUnWeight
    ("UnWeight",
     "Control the type of unweighting to perform, only one should be used the"
     " other options are for debugging purposes.",
     &IFDipole::_mode, 1, false, false);
  static SwitchOption interfaceUnWeightNoUnweighting
    (interfaceUnWeight,
     "NoUnweighting",
     "Perform no unweighting",
     0);
  static SwitchOption interfaceUnWeightAllWeights
    (interfaceUnWeight,
     "AllWeights",
     "Include all the weights",
     1);
  static SwitchOption interfaceUnWeightNoJacobian
    (interfaceUnWeight,
     "NoJacobian",
     "Only include the dipole and YFS weights",
     2);
  static SwitchOption interfaceUnWeightDipole
    (interfaceUnWeight,
     "Dipole",
     "Only include the dipole weight",
     3);
  static SwitchOption interfaceUnWeightYFS
    (interfaceUnWeight,
     "YFS",
     "Only include the YFS weight",
     4);

  static Parameter<IFDipole,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight",
     &IFDipole::_maxtry, 500, 10, 100000,
     false, false, Interface::limited);

  static Parameter<IFDipole,Energy> interfaceMinimumEnergyRest
    ("MinimumEnergyRest",
     "The minimum energy of the photons in the rest frame of the decaying particle",
     &IFDipole::_emin, MeV, 1.*MeV, ZERO, 10000.0*MeV,
     false, false, Interface::limited);

  static Parameter<IFDipole,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for unweighting",
     &IFDipole::_maxwgt, 2.0, 0.0, 100.0,
     false, false, Interface::limited);

  static Switch<IFDipole,unsigned int> interfaceEnergyCutOff
    ("EnergyCutOff",
     "The type of cut-off on the photon energy to apply",
     &IFDipole::_energyopt, 1, false, false);
  static SwitchOption interfaceEnergyCutOffRestFrame
    (interfaceEnergyCutOff,
     "RestFrame",
     "Apply cut-off in rest frame",
     1);
  static SwitchOption interfaceEnergyCutOff2
    (interfaceEnergyCutOff,
     "LabFrame",
     "Apply cut-off in lab frame",
     2);

  static Switch<IFDipole,unsigned int> interfaceBetaOption
    ("BetaOption",
     "Option for the inclusive of the higher beta coefficients",
     &IFDipole::_betaopt, 1, false, false);
  static SwitchOption interfaceBetaOptionNone
    (interfaceBetaOption,
     "None",
     "No higher betas included",
     0);
  static SwitchOption interfaceBetaOptionCollinear
    (interfaceBetaOption,
     "Collinear",
     "Include the collinear approx",
     1);
  static SwitchOption interfaceBetaOptionCollinearVirtA
    (interfaceBetaOption,
     "CollinearVirtualA",
     "Include the collinear approx with virtual corrections",
     2);
  static SwitchOption interfaceBetaOptionCollinearVirtB
    (interfaceBetaOption,
     "CollinearVirtualB",
     "Include the collinear approx with virtual corrections",
     3);
  static SwitchOption interfaceBetaOptionExact
    (interfaceBetaOption,
     "Exact",
     "Include the exact higher order terms if available",
     4);

}

ParticleVector IFDipole::generatePhotons(const Particle & p,ParticleVector children) {
  // set parameters which won't change in the event loop
  // masses of the particles
  _m[0] = p.mass(); 
  _m[1] = children[0]->mass();
  _m[2] = children[1]->mass();
  // momenta before radiation in lab
  for(unsigned int ix=0;ix<2;++ix){_qlab[ix]=children[ix]->momentum();}
  // get the charges of the particles in units of the positron charge
  // chrg1 is the charge of the parent and chrg2 is the charge of the
  // charged child. Also we create a map between the arguments of 
  // _q???[X] _m[X] etc so that 
  // _q???[_map[0]] and _m[_map[0]] are the momenta and masses of 
  // the charged child while
  // _q???[_map[1]] and _m[_map[1]] are the momenta and masses of 
  // the neutral child.
  _chrg1 = p.dataPtr()->iCharge()/3.0;
  if(children[1]->dataPtr()->iCharge()/3.0==0.0) {
    _chrg2  = children[0]->dataPtr()->iCharge()/3.0; 
    _map[0] = 0; _map[1] = 1;
  }
  else if(children[0]->dataPtr()->iCharge()/3.0==0.0) {
    _chrg2  = children[1]->dataPtr()->iCharge()/3.0; 
    _map[0] = 1; _map[1] = 0; 
  }
  // check the radiating particle is not massless
  // if(children[1]->mass()<
  if(children[_map[0]]->mass()<1e-4*GeV) { 
    ostringstream message;
    message << "IFDipole::generatePhotons() trying to generate QED radiation from "
	    << children[_map[0]]->dataPtr()->PDGName() << "\n with mass " << children[_map[0]]->mass()/GeV
	    << "which is much smaller than the mass of the electron.\n"
	    << "This is probably due to reading events from a LHEF,\nskipping radiation in this case.\n";
    generator()->logWarning( Exception(message.str(), Exception::warning));
    return children;
  } 
  // boost the momenta to the rest frame
  Boost boostv(p.momentum().boostVector());
  // boost the particles to the parent rest frame
  // and set the initial momenta of the charged particles 
  // in the dipole rest frame: currently this is the same 
  // as the boson rest frame...
  for(unsigned int ix=0;ix<2;++ix) {
    // KMH - 08/11/05 - This used to be boostv instead of -boostv
    // -boostv is the boost from the lab to the parent rest frame
    // whereas boostv goes the other way!!!
    children[ix]->deepBoost(-boostv);
    _qprf[ix]=children[ix]->momentum();
  }
  // perform the unweighting
  double wgt;
  unsigned int ntry(0);
  do {
    wgt =makePhotons(boostv,children);
    ++ntry;
    // Record warnings about large and weird weights in the .log file.
    if(wgt>_maxwgt||wgt<0.0||std::isnan(wgt)) {
      generator()->log() << "IFDipole.cc:\n";
      if(wgt>_maxwgt) {
	generator()->log() << "Weight exceeds maximum for decay!\n"; 
      } 
      if(wgt<0.0) {
	generator()->log() << "Weight is negative! \n"; 
      }
      if(std::isnan(wgt)) {
	generator()->log() << "Weight is NAN! \n";
	wgt = 0.;
      }
      generator()->log() << p.PDGName() << " " 
			 << children[0]->PDGName() << " " 
			 << children[1]->PDGName()
			 << endl 
			 << " Current Maximum = " << _maxwgt 
			 << endl
			 << " Current Weight  = " << wgt  
			 << endl;
      generator()->log() << "Photon Multiplicity      : " 
			 << _multiplicity                          << endl
			 << "Original Parent rest frame momenta: " << endl
			 << "charged child: " << ounit(_qprf[_map[0]],GeV) << endl
			 << "neutral child: " << ounit(_qprf[_map[1]],GeV) << endl
			 << "Parent rest frame momenta: "          << endl
			 << "charged child: " << ounit(_qnewprf[_map[0]],GeV)<< endl
			 << "neutral child: " << ounit(_qnewprf[_map[1]],GeV)<< endl
			 << "photons      : " << ounit(_bigLprf,GeV) << endl
			 << "Weights      : "                      << endl
			 << "_dipolewgt   : " << _dipolewgt        << endl
			 << "_yfswgt      : " << _yfswgt           << endl
			 << "_jacobianwgt : " << _jacobianwgt      << endl
			 << "_mewgt       : " << _mewgt            << endl;
      for(unsigned int ct=0;ct<_multiplicity;ct++) {
	generator()->log() << "_cosphot[" << ct << "]: " << _cosphot[ct] << endl;
	generator()->log() << "_sinphot[" << ct << "]: " << _sinphot[ct] << endl;
      }
      if(wgt>_maxwgt) {
	if(wgt<15.0) { 
	  generator()->log() << "Resetting maximum weight" 
			     << endl << " New Maximum = " << wgt  << endl;
	  _maxwgt=wgt; 
	} else {
	  generator()->log() << "Maximum weight set to limit (15)" << endl;
	  _maxwgt=15.0; 
	}
      }
    }
  } while (wgt<(_maxwgt*UseRandom::rnd()) && ntry<_maxtry);
  if(ntry>=_maxtry) {
    generator()->log() << "IFDipole Failed to generate QED radiation for the decay " 
		       << p.PDGName() << " -> " 
		       << children[0]->PDGName() << " "
		       << children[1]->PDGName() << endl;
    return children;
  }
  // produce products after radiation if needed
  if(_multiplicity>0) {
    // change the momenta of the children, they are currently
    // in parent rest frame
    for(unsigned int ix=0;ix<2;++ix) {
      LorentzRotation boost(solveBoost(_qnewprf[ix],children[ix]->momentum()));
      children[ix]->deepTransform(boost);
      // boost back to the lab
      // KMH - 08/11/05 - This used to be -boostv instead of boostv
      // -boostv is the boost from the lab to the parent rest frame
      // whereas boostv goes the other way!!!
      children[ix]->deepBoost(boostv);
    }
    // add the photons to the event record
    tcPDPtr photon=getParticleData(ParticleID::gamma);
    for(unsigned int ix=0;ix<_multiplicity;++ix) {
      PPtr newphoton=new_ptr(Particle(photon));
      newphoton->set5Momentum(_llab[ix]);
      children.push_back(newphoton);
    }
    return children;
  }
  // otherwise just return the orginial particles
  // boosted back to lab
  else {
    for(unsigned int ix=0;ix<children.size();++ix)
      children[ix]->deepBoost(boostv);
    return children;
  }
}

// member which generates the photons
double IFDipole::makePhotons(Boost boostv,ParticleVector children) {
  // set the initial parameters
  // number of photons (zero)
  _multiplicity=0;
  // zero size of photon vectors
  _lprf.clear();
  _llab.clear();
  // zero size of angle storage
  _sinphot.clear();
  _cosphot.clear();
  // zero total momenta of the photons
  _bigLprf=Lorentz5Momentum();
  // set the initial values of the reweighting factors to one
  _dipolewgt   = 1.0;
  _yfswgt      = 1.0;
  _jacobianwgt = 1.0;
  _mewgt       = 1.0;
  // set the maximum photon energy (exact - no approximations here).
  double boost_factor = 1.0;
  _emax=(0.5*(_m[0]-sqr(_m[1]+_m[2])/_m[0]))*boost_factor;
  // calculate the velocities of the children (crude/overvalued)
  double beta1(sqrt( (_qprf[_map[0]].e()+_m[_map[0]+1])
                    *(_qprf[_map[0]].e()-_m[_map[0]+1])
                   )
                   /_qprf[_map[0]].e());
  double beta2(sqrt( (_qprf[_map[1]].e()+_m[_map[1]+1]) 
                    *(_qprf[_map[1]].e()-_m[_map[1]+1])
                   )
                   /_qprf[_map[1]].e());
  // calculate 1-beta to avoid numerical problems
  double ombeta1(sqr(_m[_map[0]+1]/_qprf[_map[0]].e())/(1.+beta1));
  double ombeta2(sqr(_m[_map[1]+1]/_qprf[_map[1]].e())/(1.+beta2));
  // calculate the average photon multiplicity
  double aver(nbar(beta1,ombeta1));
  // calculate the number of photons using the poisson
  _multiplicity = UseRandom::rndPoisson(aver);
  // calculate the first part of the YFS factor
  _yfswgt/=crudeYFSFormFactor(beta1,ombeta1); 
  // generate the photon momenta with respect to q1 
  // keeping track of the weight
  double dipoles(1.);
  for(unsigned int ix=0;ix<_multiplicity;++ix)
  { dipoles *= photon(beta1,ombeta1); }
  // calculate contributions to the dipole weights so far
  _dipolewgt /=dipoles;

  // now do the momentum reshuffling
  Lorentz5Momentum pmom(ZERO,ZERO,ZERO,_m[0],_m[0]);
  if(_multiplicity>0) {
      // total energy  and momentum of photons
      Energy L0(_bigLprf.e()),modL(_bigLprf.rho());
      // squared invariant mass of final state fermions...
      Energy2 m122 = sqr(_m[0]-L0)-sqr(modL);
      if(m122<sqr(_m[1]+_m[2])) return 0.;
      // 3-momenta of charged particles
      Energy modq(_qprf[_map[0]].rho());
      // total photon momentum perpendicular to charged child...
      Energy LT(_bigLprf.perp());
      // kallen function...
      Energy4 kallen = ( m122 - sqr(_m[1]+_m[2]) )
	             * ( m122 - sqr(_m[1]-_m[2]) );
      // discriminant of rho...
      Energy4 droot = kallen-4.*sqr(_m[_map[0]+1]*LT);
      if(droot<ZERO) return 0.;
      double disc = (_m[0]-L0) *  sqrt(droot) / (2.*modq*(m122+LT*LT));
      // calculate the energy rescaling factor
      double rho  = disc-_bigLprf.z()
	          * (m122+sqr(_m[_map[0]+1])-sqr(_m[_map[1]+1]))
  	          / (2.*modq*(m122+LT*LT));
      // calculate the rescaled charged child momentum
      _qnewprf[_map[0]]=rho*_qprf[_map[0]];
      _qnewprf[_map[0]].setMass(_m[_map[0]+1]);
      _qnewprf[_map[0]].rescaleEnergy();
      // rotate the photons so in parent rest frame rather 
      // than angle measured w.r.t q1 first work out the rotation
      SpinOneLorentzRotation rotation;
      rotation.setRotateZ(-_qprf[_map[0]].phi());
      rotation.rotateY(_qprf[_map[0]].theta());
      rotation.rotateZ(_qprf[_map[0]].phi());
      // rotate the total
      _bigLprf*=rotation;
      // rotate the photons
      for(unsigned int ix=0;ix<_multiplicity;++ix){_lprf[ix]*=rotation;}
      // calculate the rescaled neutral child momentum
      _qnewprf[_map[1]]=pmom-_qnewprf[_map[0]]-_bigLprf;
      _qnewprf[_map[1]].setMass(_m[_map[1]+1]);
      _qnewprf[_map[1]].rescaleEnergy();
      // calculate the new dipole weight
      // Note this (weight) is Lorentz invariant
      // calculate velocities and 1-velocites
      beta1=sqrt( (_qnewprf[_map[0]].e()+_m[_map[0]+1])
                 *(_qnewprf[_map[0]].e()-_m[_map[0]+1]))
                /_qnewprf[_map[0]].e();
      beta2=sqrt( (_qnewprf[_map[1]].e()+_m[_map[1]+1])
                 *(_qnewprf[_map[1]].e()-_m[_map[1]+1]))
                /_qnewprf[_map[1]].e();
      ombeta1=sqr(_m[_map[0]+1]/_qnewprf[_map[0]].e())/(1.+beta1);
      ombeta2=sqr(_m[_map[1]+1]/_qnewprf[_map[1]].e())/(1.+beta2);
      for(unsigned int ix=0;ix<_multiplicity;++ix)
	{_dipolewgt*=exactDipoleWeight(beta1,ombeta1,ix);}
      // calculate the second part of the yfs form factor
      _yfswgt*=exactYFSFormFactor(beta1,ombeta1,beta2,ombeta2);
      // Now boost from the parent rest frame to the lab frame
      SpinOneLorentzRotation boost(boostv);
      // Boosting charged particles
      for(unsigned int ix=0;ix<2;++ix){_qnewlab[ix]=boost*_qnewprf[ix];}
      // Boosting total photon momentum
      _bigLlab=boost*_bigLprf;
      // Boosting individual photon momenta
      for(unsigned int ix=0;ix<_multiplicity;++ix)
        {_llab.push_back(boost*_lprf[ix]);}
      // Calculating jacobian weight
      _jacobianwgt = jacobianWeight();
      // Calculating beta^1  weight
      _mewgt = meWeight(children);
      // Apply phase space vetos...
      if(kallen<(4.*sqr(_m[_map[0]+1]*LT))||m122<sqr(_m[1]+_m[2])||rho<0.0) {  
//           generator()->log() << "Outside Phase Space" << endl;
//           generator()->log() << "Photon Multiplicity: " 
//                              << _multiplicity                          << endl
//                              << "Original Parent rest frame momenta: " << endl
//                              << "charged child: " << _qprf[_map[0]]    << endl
//                              << "neutral child: " << _qprf[_map[1]]    << endl
//                              << "rescaling    : " << rho               << endl
//                              << "Parent rest frame momenta: "          << endl
//                              << "charged child: " << _qnewprf[_map[0]] << endl
//                              << "neutral child: " << _qnewprf[_map[1]] << endl
//                              << "photons      : " << _bigLprf          << endl
//                              << endl;
	_dipolewgt   = 0.0 ;
	_yfswgt      = 0.0 ;
	_jacobianwgt = 0.0 ;
	_mewgt       = 0.0 ;
      }
      _qprf[_map[0]].rescaleEnergy();
      _qprf[_map[1]].rescaleEnergy();
      _qnewprf[_map[0]].rescaleEnergy();
      _qnewprf[_map[1]].rescaleEnergy();
      if( ((abs(_m[0]-_bigLprf.e()-_qnewprf[0].e()-_qnewprf[1].e())>0.00001*MeV)||
           (abs(      _bigLprf.x()+_qnewprf[0].x()+_qnewprf[1].x())>0.00001*MeV)||
           (abs(      _bigLprf.y()+_qnewprf[0].y()+_qnewprf[1].y())>0.00001*MeV)||
           (abs(      _bigLprf.z()+_qnewprf[0].z()+_qnewprf[1].z())>0.00001*MeV))
         &&(_dipolewgt*_jacobianwgt*_yfswgt*_mewgt>0.0)) {
	Lorentz5Momentum ptotal = _bigLprf+_qnewprf[0]+_qnewprf[1];
	ptotal.setE(ptotal.e()-_m[0]);
	generator()->log() 
	  <<   "Warning! Energy Not Conserved! tol = 0.00001 MeV"
	  << "\nwgt               = " << _dipolewgt*_yfswgt*_jacobianwgt*_mewgt
	  << "\nrho               = " << rho
	  << "\nmultiplicity      = " << _multiplicity
	  << "\n_qprf[_map[0]]    = " << _qprf[_map[0]]/GeV
	  << "\n_qprf[_map[1]]    = " << _qprf[_map[1]]/GeV
	  << "\n_qnewprf[_map[0]] = " << _qnewprf[_map[0]]/GeV << " " 
	  << _qnewprf[_map[0]].m()/GeV << " " << _m[_map[0]+1]/GeV
	  << "\n_qnewprf[_map[1]] = " << _qnewprf[_map[1]]/GeV << " " 
	  << _qnewprf[_map[1]].m()/GeV << " " << _m[_map[1]+1]/GeV
	  << "\n_bigLprf          = " << _bigLprf/GeV
	  << "\n_bigLprf.m2()     = " << _bigLprf.m2()/GeV2
	  << "\n_total out -in    = " << ptotal/GeV
	  << "\nRejecting Event.    " << "\n";
	_dipolewgt   = 0.0 ;
	_yfswgt      = 0.0 ;
	_jacobianwgt = 0.0 ;
	_mewgt       = 0.0 ;
      }
    }
  // otherwise copy momenta
  else
    { for(unsigned int ix=0;ix<2;++ix) {
        _qnewprf[ix]=_qprf[ix];
        _qnewlab[ix]=_qlab[ix]; 
      }
      _jacobianwgt = 1.0;
      // calculate the second part of the yfs form factor
      _yfswgt*=exactYFSFormFactor(beta1,ombeta1,beta2,ombeta2);
      _dipolewgt   = 1.0;
    }
// Virtual corrections for beta_0:
// These should be zero for the scalar case as there is no
// collinear singularity going by the dipoles above...
// Use mass of decaying particle...
  if(_betaopt==2) { 
    if((children[_map[0]]->dataPtr()->iSpin())==2) {
      _mewgt += (0.5*_alpha/pi) 
              * log(sqr(_m[0]
                       /_m[_map[0]+1])
                   );
    }
  }
// OR Use invariant mass of final state children...
  if(_betaopt==3) { 
    if((children[_map[0]]->dataPtr()->iSpin())==2) {
      _mewgt += (0.5*_alpha/pi)
              * log((_qnewprf[0]+_qnewprf[1]).m2()
                   /sqr(_m[_map[0]+1])
                   );
    }
  }
  // calculate the weight depending on the option
  double wgt;
  if(_mode==0){wgt=_maxwgt;}
  else if(_mode==1){wgt=_mewgt*_jacobianwgt*_yfswgt*_dipolewgt;}
  else if(_mode==2){wgt=_jacobianwgt*_yfswgt*_dipolewgt;}
  else if(_mode==3){wgt=_yfswgt*_dipolewgt;}
  else             {wgt=_yfswgt;}
  return wgt;
}

double IFDipole::photon(double beta1,double ombeta1)
{
  // generate the azimuthal angle randomly in -pi->+pi
  double phi(-pi+UseRandom::rnd()*2.*pi);
  // generate the polar angle
  double r(UseRandom::rnd());
  double costh,sinth,ombc;
  ombc  = pow(1.+beta1,1.-r)*pow(ombeta1,r);
  costh = 1./beta1*(1.-ombc);
  sinth = sqrt(ombc*(2.-ombc)-(1.+beta1)*ombeta1*sqr(costh));
  // generate the ln(energy) uniformly in ln(_emin)->ln(_emax)
  Energy energy   = pow(_emax/_emin,UseRandom::rnd())*_emin;
  // calculate the weight (omit the pre and energy factors
  //                       which would cancel later anyway)
  double wgt = 2./ombc;
  // store the angles
  _cosphot.push_back(costh);
  _sinphot.push_back(sinth);
  // store the four vector for the photon
  _lprf.push_back(Lorentz5Momentum(energy*sinth*cos(phi),energy*sinth*sin(phi),
				   energy*costh,energy,ZERO));
  // add the photon momentum to the total
  _bigLprf+=_lprf.back();
  // return the weight
  return wgt;
}

double IFDipole::meWeight(ParticleVector children)
{
  unsigned int spin = children[_map[0]]->dataPtr()->iSpin();
  double mewgt = 1.0;
  double beta1=sqrt( (_qnewprf[_map[0]].e()+_m[_map[0]+1])
                    *(_qnewprf[_map[0]].e()-_m[_map[0]+1]))
                   /_qnewprf[_map[0]].e();
  double ombeta1=sqr(_m[_map[0]+1]/_qnewprf[_map[0]].e())/(1.+beta1);
  // option which does nothing
  if(_betaopt==0){mewgt=1.;}
  // collinear approx
  else if(_betaopt==1||_betaopt==2||_betaopt==3)
    {
      double ombc;
      InvEnergy2 dipole;
      for(unsigned int i=0;i<_multiplicity;++i) {
	double opbc;
        if(_cosphot[i]<0.0)
          { opbc=ombeta1+beta1*sqr(_sinphot[i])/(1.-_cosphot[i]); }
        // if cos is greater than zero use result accurate as cos->-1
        else
          { opbc=1.+beta1*_cosphot[i]; }
        // if cos is greater than zero use result accurate as cos->1
        if(_cosphot[i]>0.0)
          { ombc=ombeta1+beta1*sqr(_sinphot[i])/(1.+_cosphot[i]); }
        // if cos is less    than zero use result accurate as cos->-1
        else
          { ombc=1.-beta1*_cosphot[i]; }
        if(((_qnewprf[_map[0]].z()>ZERO)&&(_qprf[_map[0]].z()<ZERO))||
           ((_qnewprf[_map[0]].z()<ZERO)&&(_qprf[_map[0]].z()>ZERO))) {
          dipole = sqr(beta1*_sinphot[i]/(opbc*_lprf[i].e()));
        } else {
          dipole = sqr(beta1*_sinphot[i]/(ombc*_lprf[i].e()));
	}
	// here "dipole" is the exact dipole function divided by alpha/4pi^2.

        if(spin==2) {
	  Energy magpi= sqrt( sqr(_qnewprf[_map[0]].x())
			      + sqr(_qnewprf[_map[0]].y())
			      + sqr(_qnewprf[_map[0]].z())
			      );

	  mewgt += sqr(_lprf[i].e())*_qnewprf[_map[0]].e()*ombc
	         / (sqr(magpi*_sinphot[i])*(_qnewprf[_map[0]].e()+_lprf[i].e()));
        }
        else if(spin==3) {
	  Energy2 pik  = _qnewprf[_map[0]].e()*_lprf[i].e()
	               - _qnewprf[_map[0]].x()*_lprf[i].x()
                       - _qnewprf[_map[0]].y()*_lprf[i].y()
	               - _qnewprf[_map[0]].z()*_lprf[i].z();

	  Energy2 pjk = _m[0]*_lprf[i].e();

	  Energy2 pipj = _m[0]*_qnewprf[_map[0]].e();

          mewgt += (2.*pjk*pipj/(pik*sqr(pipj+pjk))
		   +2.*pjk/(pik*(pipj+pik))
	           )/dipole;
        }
        else {
          mewgt = 1.0;
        }
      }
    }
  return mewgt;
}

double IFDipole::exactYFSFormFactor(double beta1,double ombeta1,
					   double beta2,double ombeta2) {
  double Y    = 0.0    ;
  double b    = beta1  ;
  double omb  = ombeta1;
  double c    = beta2  ;
  double omc  = ombeta2;
  double arg1 = -omc/(2.*c);
  double arg2 = -omb*omc/(2.*(b+c));
  double arg3 = 2.*b/(1.+b);
  if(_m[_map[1]+1]!=ZERO) {
    Y = _chrg1*_chrg2*(_alpha/(2.*pi))*(
         log(_m[0]*_m[_map[1]+1]/sqr(2.*_emin))
        +log(_m[_map[0]+1]*_m[_map[1]+1]/sqr(2.*_emin))
        -(1./b )*log((1.+b)/omb)*log(sqr(_m[_map[1]+1]/(2.*_emin)))
        -(1./b )*log(omb/(1.+b))
        -(0.5/b )*sqr(log(omb/(1.+b)))
        +((b+c  )/(b*omc))*log((b+c  )/(b*omc))
        -((c+b*c)/(b*omc))*log((c+b*c)/(b*omc))
        +((b+c  )/(b+b*c))*log((b+c  )/(b+b*c))
        -((c*omb)/(b+b*c))*log((c*omb)/(b+b*c))
        +(0.5/b)*(   sqr(log(  (b+c)/(b*omc)))-sqr(log((c+b*c)/(b*omc)))
       	           + sqr(log((c*omb)/(b+b*c)))-sqr(log((b+  c)/(b+b*c)))
	         )  
        +(2./b )*(   real(Math::Li2(arg1))
                   - real(Math::Li2(arg2))
                   - real(Math::Li2(arg3))
  	         )
        +(1./b )*log((b+c)/(b+b*c))*log((1.+c)/(2.*c))
        -(1./b )*log((c*omb)/(b*(1.+c)))*log((1.+b)*(1.+c)/(2.*(b+c)))
        -(1./b )*log((2.*c/b)*((b+c)/(omc*(1.+c))))*log((b+c)/(c*omb))
        );
  }
  else if(_m[_map[1]+1]==ZERO) {
    Y = _chrg1*_chrg2*(_alpha/(2.*pi))*(
         log(sqr(_m[0]/(2.*_emin)))
        +log(sqr(_m[_map[0]+1]/(2.*_emin)))
        -(1./b )*log((1.+b)/omb)
                *log((sqr(_m[0])-sqr(_m[_map[0]+1]))/sqr(2.*_emin))
        -0.5*log(omb*(1.+b)/sqr(2.*b))
        +((1.+b)/(2.*b))*log((1.+b)/(2.*b))
        -(   omb/(2.*b))*log(   omb/(2.*b))
        -(1./b )*log((1.-b)/(1.+b))
        +1.
        +(0.5/b)*sqr(log(   omb/(2.*b)))
        -(0.5/b)*sqr(log((1.+b)/(2.*b)))  
        -(0.5/b)*sqr(log((1.-b)/(1.+b)))
        -(2. /b)*real(Math::Li2(arg3))
        );
  }
  return exp(Y);
}

double IFDipole::jacobianWeight() {
  // calculate the velocities of the children (crude/overvalued)
  Energy mag1old  = sqrt( (_qprf[_map[0]].e()   +_m[_map[0]+1])
                         *(_qprf[_map[0]].e()   -_m[_map[0]+1])
                        );
  Energy mag1new  = sqrt( (_qnewprf[_map[0]].e()+_m[_map[0]+1])
                         *(_qnewprf[_map[0]].e()-_m[_map[0]+1])
			);
  Energy magL     = sqrt( sqr(_bigLprf.x())
			  + sqr(_bigLprf.y())
			  + sqr(_bigLprf.z())
			  );

// 14/12/05 - KMH - This was another mistake. This is supposed to be 
// the angel between _qnewprf[_map[0]] and _bigLprf instead of
// between _qnewprf[0] and _bigLprf. Stupid. Hopefully this weight
// is correct now.
//  double cos1L    = (_qnewprf[0].x()*_bigLprf.x()
//                    +_qnewprf[0].y()*_bigLprf.y()
//                    +_qnewprf[0].z()*_bigLprf.z()
//                    )
//                    /(mag1new*magL);
  double cos1L    = (_qnewprf[_map[0]].x()*_bigLprf.x()
                    +_qnewprf[_map[0]].y()*_bigLprf.y()
                    +_qnewprf[_map[0]].z()*_bigLprf.z()
                    )
                    /(mag1new*magL);

  return abs(  (_m[0]*sqr(mag1new)/mag1old)
             / (  mag1new*(_m[0]-_bigLprf.e())
                +_qnewprf[_map[0]].e()*magL*cos1L
	       )
            );
}

LorentzRotation IFDipole::solveBoost(const Lorentz5Momentum & q, 
				     const Lorentz5Momentum & p ) const {
  Energy modp = p.vect().mag();
  Energy modq = q.vect().mag();
  double betam = (p.e()*modp-q.e()*modq)/(sqr(modq)+sqr(modp)+p.mass2());
  Boost beta = -betam*q.vect().unit();
  ThreeVector<Energy2> ax = p.vect().cross( q.vect() ); 
  double delta = p.vect().angle( q.vect() );
  LorentzRotation R;
  using Constants::pi;
  if ( ax.mag2()/GeV2/MeV2 > 1e-16 ) {
    R.rotate( delta, unitVector(ax) ).boost( beta );
  } 
  else {
    if(p.mass()>ZERO) {
      R.boost(p.findBoostToCM(),p.e()/p.mass());
      R.boost(q.boostVector(),q.e()/q.mass());
    }
    else {
      if(modp>modq) beta = -betam*p.vect().unit();
      R.boost( beta );
    }
  } 
  return R;
}

void IFDipole::doinit() {
  Interfaced::doinit();
  // get the value fo alpha from the Standard Model object
  _alpha=generator()->standardModel()->alphaEM();
}
#line 1 "./YFSFormFactors.cc"
// -*- C++ -*-
//
// YFSFormFactors.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the YFSFormFactors class.
//

#include "YFSFormFactors.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include <cassert>

using namespace Herwig;

using Constants::pi;
using Herwig::Math::ReLi2;

const double  YFSFormFactors::_alpha=1./137.03599911;

const Energy  YFSFormFactors::_mgamma=1e-10*MeV;

const Energy2 YFSFormFactors::_tcut=1.e-11*GeV2;

const Energy YFSFormFactors::_ecut=1e-6*GeV;

double YFSFormFactors::ReBIF(Energy  m0      ,Energy  m1      , Energy2 t       ,
			     double  charge  ,bool    includegamma,
			     Energy  mgamma) {
  // mass squared for speed
  Energy2 m02(m0*m0),m12(m1*m1),nu(0.5*(m02+m12-t)),mprod(m0*m1);
  double Anu,vfinite;
  double output;
  // t>0
  if(t>_tcut) {
    // parameters
    Energy2 lambda(sqrt((nu-mprod)*(nu+mprod)));
    double eta(0.5*m12*t/lambda/(lambda+nu-m12)),zeta((lambda+nu)*eta/m12);
    // simple A functions for virtual piece
    InvEnergy2 A;
    if(lambda>1e-6*GeV2){A=(log((lambda+nu)/mprod)/lambda);}
    else{A=1./mprod;}
    double A1((m02-m12)/t*log(m0/m1)-2.*sqr(lambda)/t*A-2.);
    InvEnergy2 A3(A*log(2.*lambda/mprod)
		  +1./lambda*
		  (+0.25*(log((lambda+nu)/m02)+2.*log((lambda-nu+m02)/t  ))*
		   log((lambda+nu)/m02)
		   +0.25*(log((lambda+nu)/m12)-2.*log((lambda+nu-m12)/m12))*
		   log((lambda+nu)/m12)
		   +0.5*(log(eta)*log(1.+eta)-log(zeta)*log(1.+zeta))
		   +ReLi2(-eta)-ReLi2(-zeta)));
    Anu=nu*A;
    vfinite=0.5*A1-nu*A3;
  }
  // t==0
  else {
    // virtual part of the dipole
    Anu = (m02+m12)/(m02-m12)*log(m0/m1);
    vfinite=0.5*(Anu-1.);
  }
  if(includegamma){output=-_alpha*charge/pi*((Anu-1.)*log(sqr(mgamma)/mprod)+vfinite);}
  else            {output=-_alpha*charge/pi*((Anu-1.)*log(MeV2/mprod)+vfinite);}
  //  assert(isfinite(output));
  return output;
}

double YFSFormFactors::ReBFF(Energy m1,Energy m2,Energy2 s,double  charge,
			     bool    includegamma,Energy  mgamma) {
  // masses etc
  Energy2 m12(m1*m1),m22(m2*m2),mu(0.5*(s-m12-m22)),mprod(m1*m2);
  // parameters
  double ratio(m1*m2/mu),rho(sqrt((1.-ratio)*(1.+ratio)));
  Energy2 prod(mu*(1.+rho));
  // the finite piece
  double vfinite(mu*rho/s*log(prod/mprod)+0.5*(m12-m22)/s*log(m1/m2)
		 +1./rho*(pi*pi-0.5*log(prod/m12)*log(prod/m22)
			  -0.5*sqr(log((m12+prod)/(m22+prod)))
			  -ReLi2(2.*mu*rho/(m12+prod))
			  -ReLi2(2.*mu*rho/(m22+prod)))-1.);
  // the cut-off piece
  double Anu(log(prod/mprod)/rho),output;
  if(includegamma){output=-_alpha*charge/pi*((Anu-1.)*log(sqr(mgamma)/mprod)+vfinite);}
  else            {output=-_alpha*charge/pi*((Anu-1.)*log(MeV2/mprod)+vfinite);}
  //  assert(isfinite(output));
  return output;
}

double YFSFormFactors::BtildeIF(double  beta0   ,double  ombeta0 ,
				double  beta1   ,double  ombeta1 ,
				Energy  en0     ,Energy  en1     ,
				Energy  m0      ,Energy  m1      , 
				Energy2 t       ,double  charge  ,
				Energy  emin    ,bool    includegamma,
				Energy  mgamma) {
  // coefficient of the divergent piece
  Energy2 mprod(m0*m1),nu(0.5*(m0*m0+m1*m1-t));
  double Anu;
  if(nu-mprod>1e-12*GeV2) {
    Energy2 lambda(sqrt((nu-mprod)*(nu+mprod)));
    Anu=nu/lambda*log((lambda+nu)/mprod);
  }
  else
    {Anu=1.;}
  // finite piece
  double rfinite(-0.5*A4single(beta0,ombeta0)-0.5*A4single(beta1,ombeta1)
		 +nu*A4IF(beta0,ombeta0,beta1,ombeta1,en0,en1,m0,m1,t));
  //  assert(isfinite(rfinite));
  // return the answer
  double output; 
  if(includegamma) {
    output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/mgamma)+rfinite);
  }
  else {
    output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/MeV)+rfinite);
  }
  //  assert(isfinite(output));
  return output;
}

double YFSFormFactors::BtildeFF(double  beta1   ,double  ombeta1 ,
				       double  beta2   ,double  ombeta2 ,
				       Energy  en1     ,Energy  en2     ,
				       Energy  m1      ,Energy  m2      , 
				       Energy2 s       ,double  charge  ,
				       Energy  emin    ,bool    includegamma,
				       Energy  mgamma) {
  // masses etc
  Energy2 m12(m1*m1),m22(m2*m2),mu(0.5*(s-m12-m22)),mprod(m1*m2);
  // parameters
  double ratio(m1*m2/mu),rho(sqrt((1.-ratio)*(1.+ratio)));
  Energy2 prod(mu*(1.+rho));
  // finite piece
  double rfinite(-0.5*A4single(beta1,ombeta1)-0.5*A4single(beta2,ombeta2)
		 +mu*A4FFFull(en1,en2,beta1,beta2,m1,m2,s));
  double Anu(log(prod/mprod)/rho);
  // return the answer
  double output; 
  if(includegamma){output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/mgamma)+rfinite);}
  else            {output=-_alpha*charge/pi*((Anu-1.)*2.*log(2.*emin/MeV)+rfinite);}
  //  assert(isfinite(output));
  return output;
}

InvEnergy2 YFSFormFactors::A4FFFull(Energy  inen1  ,Energy inen2,
					   double  beta1,double beta2,
					   Energy   inm1  ,Energy inm2,Energy2 s    ) {
  Energy en1(inen1),en2(inen2),m1(inm1),m2(inm2);
  // order the particles so en1>en2
  if(inen1*beta1<inen2*beta2) {
    en1=inen2;
    en2=inen1;
    m1=inm2;
    m2=inm1;
  }
  Energy Delta(en1-en2);
  Energy Omega(en1+en2),delta(m1-m2),omega(m1+m2);
  Energy2 Q2(s-2.*(m1*m1+m2*m2));
  Energy root(sqrt(Delta*Delta+Q2));
  Energy eta[2]={sqrt((en2-m2)*(en2+m2)),sqrt((en1-m1)*(en1+m1))+root};
  if(0.5*(s-m1*m1-m2*m2)>en1*en2){eta[0]=-eta[0];}
  Energy2 root2(sqrt((Q2+omega*omega)*(Q2+delta*delta)));
  double Y[2];
  // various limits
  Energy y[4];
  y[0]=0.5*(root-Omega+(omega*delta+root2)/(root+Delta));
  y[1]=y[0]-root2/(root+Delta);
  y[2]=0.5*(root+Omega+(omega*delta+root2)/(root-Delta));
  y[3]=y[2]-root2/(root-Delta);
  // the Y function at both limits
  for(unsigned int ix=0;ix<2;++ix)
    {Y[ix]=Zij(eta[ix],y[0],y[3])+Zij(eta[ix],y[1],y[0])
	+Zij(eta[ix],y[2],y[1])-Zij(eta[ix],y[2],y[3])
	+0.5*Xijkl(eta[ix],y[0],y[1],y[2],y[3])*Xijkl(eta[ix],y[1],y[2],y[0],y[3]);}
  // the answer
  // the Z function at both limits
  double output(0.);
  if(abs(Delta)>_ecut) {
    output=log(abs((root-Delta)/(root+Delta)))*(+Xijkl(eta[1],y[0],y[3],y[1],y[2])
						-Xijkl(eta[0],y[0],y[3],y[1],y[2]));
  }
  return 1./root2*(output+Y[1]-Y[0]);
}

InvEnergy2 YFSFormFactors::A4IF(double  beta0   ,double  ombeta0 ,
				       double  beta1   ,double  ombeta1 ,
				       Energy  en0  ,Energy en1  , Energy  m0   ,Energy m1   ,
				       Energy2 t) {
  // this is the general function so pick the special case
  if(t>_tcut){
    // rest frame of decaying particle t!=0
    if(abs(en0-m0)<_ecut){return A4IFRest(m0,m1,beta1,ombeta1,en1);}
    // rest frame of decay product t!=0
    else if(abs(en1-m1)<_ecut){return A4IFRest(m1,m0,beta0,ombeta0,en0);}
    // general frame t!=0
    else
      {return A4IFFull(beta0,beta1,en0,en1,m0,m1,t);}
  }
  else {
    // rest frame of decaying particle t=0
    if(abs(en0-m0)<_ecut){return A4IFRestZero(m0,m1);}
    // rest frame of decay products t=0
    else if(abs(en1-m1)<_ecut){return A4IFRestZero(m1,m0);}
    // general frame t=0
    else{return A4IFZero(beta0,beta1,ombeta1,en0,en1,m0,m1);}
  }
}

InvEnergy2 YFSFormFactors::A4IFZero(double  beta0, double beta1, double ombeta1,
				    Energy  en0,
				    Energy en1  , Energy  m0   , Energy m1) {
  Energy  Delta = en0-en1;
  Energy2 mu2 = (m0-m1)*(m0+m1);
  long double z[2]={ beta1*en1/Delta, beta0*en0/Delta-1. };
  long double y[3],xi[3];
  y[0]=en1/Delta;
  y[1]=y[0]-0.5*mu2/sqr(Delta);
  y[2]=-y[0]+2.*m1*m1/mu2;
  for(unsigned int ix = 0; ix < 3; ++ix) {
    if ( ix == 0 ) xi[0]  = -ombeta1*y[0]  / (z[1] - y[0] );
    else           xi[ix] = (z[0] - y[ix]) / (z[1] - y[ix]);
  }
  long double U[2];
  for(unsigned int ix=0;ix<2;++ix) {
    // U[ix] = 0.5*sqr(log(abs((z[ix]-y[0])*(z[ix]-y[1])/(z[ix]-y[2]))))
    //   +log(abs(z[ix]-y[0]))*log(abs(z[ix]-y[0])/sqr(z[ix]-y[1]))
    //   +2.*ReLi2((y[1]-y[0])/(z[ix]-y[0]))
    //   +2.*ReLi2((y[2]-y[1])/(z[ix]-y[1]));
    const long double a = ix==0 ? -ombeta1*y[0] : z[ix]-y[0];
    const long double b = z[ix]-y[1];
    const long double c = z[ix]-y[2];
    const long double A = abs(a*b/c);
    const long double B = abs(a);
    const long double C = B/sqr(b);
    const long double D = (y[1]-y[0])/a;
    const long double E = (y[2]-y[1])/b;
    U[ix] = 0.5*sqr(log(A)) + log(B)*log(C) + 2.*ReLi2(D) + 2.*ReLi2(E);
  }
  return 1./mu2*(log(2.*sqr(Delta)/mu2)*log(abs(xi[1]*xi[2]/xi[0]))+U[1]-U[0]);
}

InvEnergy2 YFSFormFactors::A4IFRest(Energy m0   ,Energy m1, double beta1,
					   double ombeta1, Energy E1) {
  Energy  Mfact0 = m0-E1*ombeta1;
  Energy  Mfact1 = m0-E1*(1.+beta1);
  Energy2 Mfact2 = m0*E1*(1.+beta1)-m1*m1;
  Energy2 Mfact3 = m0*E1*ombeta1-m1*m1;
  Energy2 qprod(m0*E1*beta1);
  return 0.5/qprod*(+log(abs(Mfact0/Mfact1))*log(E1*(1.+beta1)/m0)
		    -2.*log(abs(2.*beta1*E1*Mfact0/m0/m1))*log(E1*(1.+beta1)/m1)
		    +2.*ReLi2(E1/m0*ombeta1)-2.*ReLi2(E1/m0*(1.+beta1))
		    +ReLi2(-0.5*Mfact1/beta1/E1)-ReLi2( 0.5*Mfact0/beta1/E1)
		    +ReLi2( 0.5*Mfact2/qprod   )-ReLi2(-0.5*Mfact3/qprod));
}

InvEnergy2 YFSFormFactors::A4IFFull(Velocity beta0,Velocity beta1,
					   Energy  en0  ,Energy en1  ,
					   Energy  m0   ,Energy m1   , Energy2 t) {
  Energy Delta(en0-en1),Omega(en0+en1),delta(m0-m1),omega(m0+m1);
  Energy  T(sqrt(sqr(Delta)-t)),V(Delta+T);
  Energy2 kappa(sqrt((sqr(omega)-t)*(sqr(delta)-t)));
  long double y[4]={-0.5/T*(T+Omega-(omega*delta+kappa)*V/t),
		    -0.5/T*(T+Omega-(omega*delta-kappa)*V/t),
		    -0.5/T*(T-Omega+(omega*delta+kappa)/V),
		    -0.5/T*(T-Omega+(omega*delta-kappa)/V)};
  long double z[2]={beta1*en1/T,beta0*en0/T-1.};
  double Y[2],lfact(log(abs(V*V/t)));
  for(unsigned int ix=0;ix<2;++ix) {
    Y[ix] = lfact*Xijkl(z[ix],y[0],y[3],y[1],y[2])
      +Zij(z[ix],y[0],y[3])
      +Zij(z[ix],y[1],y[0])
      +Zij(z[ix],y[2],y[1])
      -Zij(z[ix],y[2],y[3])
      +0.5*Xijkl(z[ix],y[0],y[1],y[2],y[3])*Xijkl(z[ix],y[1],y[2],y[0],y[3]);
  }
  return (Y[1]-Y[0])/kappa;
}
