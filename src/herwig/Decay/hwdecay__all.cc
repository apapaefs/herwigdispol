#line 1 "./DecayIntegrator.cc"
// -*- C++ -*-
//
// DecayIntegrator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayIntegrator class.
//

#include "DecayIntegrator.h"
#include "PhaseSpaceMode.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/WidthCalculatorBase.h"

using namespace Herwig;

void DecayIntegrator::persistentOutput(PersistentOStream & os) const {
  os << modes_ << nIter_ << nPoint_ << nTry_
     << photonGen_ << generateInter_ << ounit(eps_,GeV) << warnings_;
}

void DecayIntegrator::persistentInput(PersistentIStream & is, int) {
  is >> modes_ >> nIter_ >> nPoint_ >> nTry_
     >> photonGen_ >> generateInter_ >> iunit(eps_,GeV) >> warnings_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<DecayIntegrator,HwDecayerBase>
describeHerwigDecayIntegrator("Herwig::DecayIntegrator", "Herwig.so");

void DecayIntegrator::Init() {
  static ClassDocumentation<DecayIntegrator> documentation
    ("The DecayIntegrator class is a base decayer class "
     "including a multi-channel integrator.");

  static Parameter<DecayIntegrator,unsigned int> interfaceIteration
    ("Iteration",
     "Number of iterations for the initialization of the phase space",
     &DecayIntegrator::nIter_, 10, 0, 100,
     false, false, true);  
  
  static Parameter<DecayIntegrator,unsigned int> interfacePoints
    ("Points",
     "number of phase space points to generate in the initialisation.",
     &DecayIntegrator::nPoint_, 10000, 1, 1000000000,
     false, false, true);
  
  static Parameter<DecayIntegrator,unsigned int> interfaceNtry
    ("Ntry",
     "Number of attempts to generate the decay",
     &DecayIntegrator::nTry_, 500, 0, 100000,
     false, false, true);

  static Reference<DecayIntegrator,DecayRadiationGenerator> interfacePhotonGenerator
    ("PhotonGenerator",
     "Object responsible for generating photons in the decay.",
     &DecayIntegrator::photonGen_, false, false, true, true, false);
 
  static Switch<DecayIntegrator,bool> interfaceGenerateIntermediates
    ("GenerateIntermediates",
     "Whether or not to include intermediate particles in the output",
     &DecayIntegrator::generateInter_, false, false, false);
  static SwitchOption interfaceGenerateIntermediatesNoIntermediates
    (interfaceGenerateIntermediates,
     "No",
     "Don't include the intermediates",
     false);
  static SwitchOption interfaceGenerateIntermediatesIncludeIntermediates
    (interfaceGenerateIntermediates,
     "Yes",
     "include the intermediates",
     true);

   static Switch<DecayIntegrator, bool> InterfacePhaseSpaceWarning
     ("PhaseSpaceWarning",
      "Switch on/off text warnings in PhaseSpaceMode class",
      &DecayIntegrator::warnings_, false, false, false);
   static SwitchOption on
     (InterfacePhaseSpaceWarning,"on","turn on the warnings", true);
   static SwitchOption off
     (InterfacePhaseSpaceWarning,"off","turn off the warnings", false);

}

double DecayIntegrator::oneLoopVirtualME(unsigned int ,
					  const Particle &, 
					  const ParticleVector &) {
  throw Exception()
    << "DecayIntegrator::oneLoopVirtualME() called. This should"
    << " have been overidden in an inheriting class if it is used"
    << Exception::runerror;
}

InvEnergy2 DecayIntegrator::realEmissionME(unsigned int,
					    const Particle &, 
					    ParticleVector &,
					    unsigned int,
					    double, double, 
					    const LorentzRotation &,
					    const LorentzRotation &) {
  throw Exception()
    << "DecayIntegrator::realEmmisionME() called. This should"
    << " have been overidden in an inheriting class if it is used"
    << Exception::runerror;
}

ParticleVector DecayIntegrator::decay(const Particle & parent,
				      const tPDVector & children) const {
  // return empty vector if products heavier than parent
  Energy mout(ZERO);
  for(tPDPtr pd : children)  mout += pd->massMin();
  if(mout>parent.mass()) return ParticleVector();
  // generate the decay
  bool cc;
  iMode_ = modeNumber(cc,parent.dataPtr(),children);
  if(numberModes()==0) return ParticleVector();
  return modes_[iMode_]->generateDecay(parent,this,generateInter_,cc);
}

void DecayIntegrator::doinitrun() {
  HwDecayerBase::doinitrun();
  if ( initialize() && Debug::level > 1 ) 
    CurrentGenerator::current().log() << "Start of the initialisation for " 
				      << name() << "\n";
  for(unsigned int ix=0;ix<modes_.size();++ix) {
    if(!modes_[ix]) continue;
    modes_[ix]->initrun();
    iMode_=ix;
    modes_[ix]->initializePhaseSpace(initialize(),this);
  }
}

// output the information for the database
void DecayIntegrator::dataBaseOutput(ofstream & output,bool header) const {
  // header for MySQL
  if(header) output << "update decayers set parameters=\"";
  output << "newdef " << name() << ":Iteration " << nIter_ << "\n";
  output << "newdef " << name() << ":Ntry " << nTry_ << "\n";
  output << "newdef " << name() << ":Points " << nPoint_ << "\n";
  //if(_photongen){;}
  output << "newdef " << name() << ":GenerateIntermediates " << generateInter_ << " \n";
  // footer for MySQL
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";\n";
  }
}

// set the code for the partial width
void DecayIntegrator::setPartialWidth(const DecayMode & dm, int imode) {
  int ifound = findMode(dm);
  if(ifound>=0) modes_[ifound]->setPartialWidth(imode);
}

WidthCalculatorBasePtr 
DecayIntegrator::threeBodyMEIntegrator(const DecayMode &) const {
  return WidthCalculatorBasePtr();
}

int DecayIntegrator::findMode(const DecayMode & dm) {
  int imode(-1);
  vector<int> extid;
  bool found(false);
  int id;
  unsigned int ix(0),iy,N,iz,tmax,nmatched;
  if(modes_.size()==0) return -1;
  do {
    if(!modes_[ix]) {
      ++ix;
      continue;
    }
    tcPDPtr in = modes_[ix]->incoming().first;
    tcPDPtr cc = modes_[ix]->incoming().first->CC();
    tmax=1;if(!cc){++tmax;}
    for(iz=0;iz<tmax;++iz) {
      extid.clear();
      // check the parent
      if(dm.parent()!=in && dm.parent()!=cc) continue;
      if(dm.parent()->id()==in->id()&&iz==0) {
	for(iy=0,N=modes_[ix]->numberOfParticles();iy<N;++iy) {
	  extid.push_back(modes_[ix]->outgoing()[iy]->id());
	}
      }
      else if(dm.parent()->id()==in->id()&&iz==1) {
	for(iy=0,N=modes_[ix]->numberOfParticles();iy<N;++iy) {
	  tcPDPtr cc2=modes_[ix]->outgoing()[iy]->CC();
	  extid.push_back( cc2 ? cc2->id() : modes_[ix]->outgoing()[iy]->id());
	}
      }
      else if(cc&&dm.parent()->id()==cc->id()) {
	for(iy=0,N=modes_[ix]->numberOfParticles();iy<N;++iy) {
	  tcPDPtr cc2 = modes_[ix]->outgoing()[iy]->CC();
	  extid.push_back( cc2 ? cc2->id() : modes_[ix]->outgoing()[iy]->id());
	}
      }
      // if the parents match
      if(!extid.empty()) {
	vector<bool> matched(extid.size(),false);
	bool done;
	nmatched=0;
	ParticleMSet::const_iterator pit = dm.products().begin();
	do {
	  id=(**pit).id();
	  done=false;
	  iy=0;
	  do {
	    if(id==extid[iy]&&!matched[iy]) {
	      matched[iy]=true;
	      ++nmatched;
	      done=true;
	    }
	    ++iy;
	  }
	  while(iy<extid.size()&&!done);
	  ++pit;
	}
	while(pit!=dm.products().end());
	if(nmatched==extid.size()) {
	  imode=ix;
	  found=true;
	}
      }
    }
    ++ix;
  }
  while(!found&&ix<modes_.size());
  return imode;
}

// the matrix element to be integrated for the me
double DecayIntegrator::threeBodyMatrixElement(const int,const Energy2,
					       const Energy2,
					       const Energy2,const Energy2,
					       const Energy, const Energy, 
					       const Energy) const {
  throw Exception() 
    << "Calling the virtual DecayIntegrator::threeBodyMatrixElement"
    << "method. This must be overwritten in the classes "
    << "inheriting from DecayIntegrator where it is needed"
    << Exception::runerror;
}

// the differential three body decay rate with one integral performed
InvEnergy DecayIntegrator::threeBodydGammads(const int, const Energy2,
					     const Energy2,
					     const Energy, const Energy, 
					     const Energy) const {
  throw Exception() 
    << "Calling the virtual DecayIntegrator::threeBodydGammads()" 
    <<"method. This must be overwritten in the classes "
    << "inheriting from DecayIntegrator where it is needed"
    << Exception::runerror;
}

// generate the momenta for the decay
ParticleVector DecayIntegrator::generate(bool inter,bool cc,
					 const unsigned int & imode,
					 const Particle & inpart) const {
  iMode_=imode;
  return modes_[imode]->generateDecay(inpart,this,inter,cc);
}

void  DecayIntegrator::addMode(PhaseSpaceModePtr mode) const {
  modes_.push_back(mode);
  if(mode) mode->init();
}

ostream & Herwig::operator<<(ostream & os, const DecayIntegrator & decay) { 
  os << "The integrator has " << decay.modes_.size() << " modes"  << endl;
  for(unsigned int ix=0;ix<decay.modes_.size();++ix) {
    os << "Information on mode " << ix << endl;
    os << *(decay.modes_[ix]);
  }
  return os;
}
  
// reset the properities of all intermediates
void DecayIntegrator::resetIntermediate(tcPDPtr part, Energy mass, Energy width) {
  if(!part) return;
  for(unsigned int ix=0,N=modes_.size();ix<N;++ix) {
    modes_[ix]->resetIntermediate(part,mass,width);
  }
}

Energy DecayIntegrator::initializePhaseSpaceMode(unsigned int imode,bool init, bool onShell) const{
  tcPhaseSpaceModePtr cmodeptr=mode(imode);
  tPhaseSpaceModePtr modeptr = const_ptr_cast<tPhaseSpaceModePtr>(cmodeptr);
  modeptr->init();
  return modeptr->initializePhaseSpace(init,this,onShell);
}
#line 1 "./PhaseSpaceChannel.cc"
// -*- C++ -*-
//
// PhaseSpaceChannel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhaseSpaceChannel class.
//
// Author: Peter Richardson
// 

#include "PhaseSpaceChannel.h"
#include "PhaseSpaceMode.h"
#include "Herwig/Utilities/Kinematics.h"

/** 
 *  Constructor with incoming particles
 */
PhaseSpaceChannel::PhaseSpaceChannel(tPhaseSpaceModePtr inm, bool skip) : mode_(inm), weight_(1.),
									  initialized_(false),
									  skipFirst_(skip) {
  if(!inm->incoming().second)
     intermediates_.push_back(PhaseSpaceResonance(inm->incoming().first));
}

void PhaseSpaceChannel::init(tPhaseSpaceModePtr mode) {
  mode_=mode;
  if(initialized_) return;
  initialized_=true;
  // find the descendents
  for(PhaseSpaceResonance & res : intermediates_)
    findChildren(res,res.descendents);
  // ensure intermediates either have the width set, or
  // can't possibly be on-shell
  // first the maximum energy release
  Energy massmax = mode->eMax_;
  for(tcPDPtr part : mode->outgoing_)
    massmax -= mode->testOnShell_ ? part->mass() : part->massMin();
  for(PhaseSpaceResonance & res : intermediates_) {
    if(!res.particle || res.mWidth!=ZERO ||
       res.jacobian != PhaseSpaceResonance::BreitWigner) continue;
    Energy massmin(ZERO);
    for(const int & ext : res.descendents)
      massmin += mode->testOnShell_ ? 
	mode->outgoing_[ext-1]->mass() : mode->outgoing_[ext-1]->massMin();
    // check if can be on-shell
    Energy mass = sqrt(res.mass2);
    if(mass>=massmin&&mass<=massmax+massmin) {
      string modeout = mode->incoming_.first->PDGName() + " ";
      if(mode->incoming_.second)
	modeout += mode->incoming_.second->PDGName() + " ";
      for( tcPDPtr out : mode->outgoing_)
	modeout += out->PDGName() + " ";
      throw InitException() << "Width zero for " << res.particle->PDGName()
 			      << " in PhaseSpaceChannel::init() "
			    << modeout << Exception::runerror;
    }
  }
}

ostream & Herwig::operator<<(ostream & os, const PhaseSpaceChannel & channel) {
  // output of the external particles
  if(!channel.mode_->incoming().second) {
    os << "Channel for the decay of " << channel.mode_->incoming().first->PDGName() << " -> ";
  }
  else 
    os << "Channel for " << channel.mode_->incoming().first->PDGName()
       << ", " << channel.mode_->incoming().second->PDGName()
       << " -> ";
  for(tcPDPtr part : channel.mode_->outgoing())
    os << part->PDGName() << " ";
  os << endl;
  os << "Proceeds in following steps ";
  for(unsigned int ix=0;ix<channel.intermediates_.size();++ix) {
    os << channel.intermediates_[ix].particle->PDGName() << " -> ";
    if(channel.intermediates_[ix].children.first>0) {
      os << channel.mode_->outgoing()[channel.intermediates_[ix].children.first-1]->PDGName()  
  	 << "(" << channel.intermediates_[ix].children.first<< ") ";
    }
    else {
      os << channel.intermediates_[-channel.intermediates_[ix].children.first].particle->PDGName() 
     	 << "(" << channel.intermediates_[ix].children.first<< ") ";
    }
    if(channel.intermediates_[ix].children.second>0) {
      os << channel.mode_->outgoing()[channel.intermediates_[ix].children.second-1]->PDGName()  
  	 << "(" << channel.intermediates_[ix].children.second<< ") ";
    }
    else {
      os << channel.intermediates_[-channel.intermediates_[ix].children.second].particle->PDGName() 
     	 << "(" << channel.intermediates_[ix].children.second<< ") ";
    }
    os << endl;
  }
  return os;
}
void PhaseSpaceChannel::initrun(tPhaseSpaceModePtr mode) {
  mode_=mode;
  if(!mode_->testOnShell_) return;
  // ensure intermediates either have the width set, or
  // can't possibly be on-shell
  Energy massmax = mode_->incoming().first->massMax();
  for(tPDPtr out : mode_->outgoing()) 
    massmax -= out->massMin();
  for(unsigned int ix=1;ix<intermediates_.size();++ix) { 
    if(intermediates_[ix].mWidth==ZERO && intermediates_[ix].jacobian==PhaseSpaceResonance::BreitWigner) {
      Energy massmin(ZERO);
      for(const int & iloc : intermediates_[ix].descendents) 
	massmin += mode_->outgoing()[iloc-1]->massMin();
      // check if can be on-shell
      if(intermediates_[ix].mass2>=sqr(massmin)&&
	 intermediates_[ix].mass2<=sqr(massmax+massmin)) {
	string modeout = mode_->incoming().first->PDGName() + " -> ";
	for(tPDPtr out : mode_->outgoing()) 
	  modeout += out->PDGName() + " ";
	throw Exception() << "Width zero for " << intermediates_[ix].particle->PDGName()
			  << " in PhaseSpaceChannel::initrun() "
			  << modeout
			  << Exception::runerror;
      }
    }
  }
}

// generate the momenta of the external particles
vector<Lorentz5Momentum> 
PhaseSpaceChannel::generateMomenta(const Lorentz5Momentum & pin,
				   const vector<Energy> & massext) const {
  // storage of the momenta of the external particles
  vector<Lorentz5Momentum> pexternal(massext.size());
  // and the internal particles
  vector<Lorentz5Momentum> pinter(1,pin);
  pinter.resize(intermediates_.size());
  // masses of the intermediate particles
  vector<Energy> massint(intermediates_.size());
  massint[0]=pin.mass();
  // generate all the decays in the chain
  for(unsigned int ix=0;ix<intermediates_.size();++ix) { 
    int idau[2] = {abs(intermediates_[ix].children.first),
		   abs(intermediates_[ix].children.second)};
    // if both decay products off-shell
    if(intermediates_[ix].children.first<0&&intermediates_[ix].children.second<0) {
      double rnd  = mode_->rStack_.top();
      mode_->rStack_.pop();
      Energy lowerb[2];
       // lower limits on the masses of the two resonances
      for(unsigned int iy=0;iy<2;++iy) {
   	lowerb[iy]=ZERO;
   	bool massless=true;
   	for(const int & des : intermediates_[idau[iy]].descendents) {
   	  if(massext[des-1]!=ZERO) massless = false;
   	  lowerb[iy]+=massext[des-1];
  	}
   	if(massless) lowerb[iy] = mode_->epsilonPS(); 
      }
      double rnd2 = mode_->rStack_.top();
      mode_->rStack_.pop();
      // randomize the order
      if(rnd<0.5) {
	rnd*=2.;
  	// mass of the first resonance
  	Energy upper = massint[ix]-lowerb[1];
  	Energy lower = lowerb[0];
   	massint[idau[0]]=generateMass(intermediates_[idau[0]],lower,upper,rnd );
   	// mass of the second resonance
   	upper = massint[ix]-massint[idau[0]];
   	lower = lowerb[1];
   	massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper,rnd2);
      }
      else {
	rnd = 2.*rnd-1.;
  	// mass of the second resonance
  	Energy upper = massint[ix]-lowerb[0];
  	Energy lower = lowerb[1];
   	massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper,rnd2);
   	// mass of the first resonance
   	upper = massint[ix]-massint[idau[1]];
   	lower = lowerb[0];
   	massint[idau[0]]=generateMass(intermediates_[idau[0]],lower,upper,rnd );
      }
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massint[idau[1]],
   		   pinter[idau[0]],pinter[idau[1]]);
    }
    // only first off-shell
    else if(intermediates_[ix].children.first<0) {
      double rnd  = mode_->rStack_.top();
      mode_->rStack_.pop();
      // compute the limits of integration
      Energy upper = massint[ix]-massext[idau[1]-1];
      Energy lower = ZERO;
      bool massless=true;
      for(const int & des : intermediates_[idau[0]].descendents) {
   	if(massext[des-1]!=ZERO) massless = false;
   	lower+=massext[des-1];
      }
      if(massless) lower = mode_->epsilonPS();
      massint[idau[0]] = generateMass(intermediates_[idau[0]],lower,upper,rnd);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massint[idau[0]],massext[idau[1]-1], 
		   pinter[idau[0]],pexternal[idau[1]-1]);
    }
    // only second off-shell
    else if(intermediates_[ix].children.second<0) {
      double rnd  = mode_->rStack_.top();
      mode_->rStack_.pop();
      // compute the limits of integration
      Energy upper = massint[ix]-massext[idau[0]-1];
      Energy lower = ZERO;
      bool massless=true;
      for(const int & des : intermediates_[idau[1]].descendents) {	
 	if(massext[des-1]!=ZERO) massless = false;
 	lower+=massext[des-1];
      }
      if(massless) lower = mode_->epsilonPS();
      massint[idau[1]]=generateMass(intermediates_[idau[1]],lower,upper,rnd);
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massext[idau[0]-1],massint[idau[1]], 
   		   pexternal[idau[0]-1],pinter[idau[1]]);
    }
    // both on-shell
    else {
      // generate the momenta of the decay products
      twoBodyDecay(pinter[ix],massext[idau[0]-1],massext[idau[1]-1], 
   		   pexternal[idau[0]-1],pexternal[idau[1]-1]);
    }
  }
  // return the external momenta
  return pexternal;
}
void PhaseSpaceChannel::twoBodyDecay(const Lorentz5Momentum & p,
				     const Energy m1, const Energy m2,
				     Lorentz5Momentum & p1,
				     Lorentz5Momentum & p2 ) const {
  static const double eps=1e-6;
  double ctheta = 2.*mode_->rStack_.top()-1.;
  mode_->rStack_.pop();
  double phi = Constants::twopi*mode_->rStack_.top();
  mode_->rStack_.pop();
  Kinematics::generateAngles(ctheta,phi);
  Axis unitDir1=Kinematics::unitDirection(ctheta,phi);
  Momentum3 pstarVector;
  Energy min=p.mass();
  if ( min >= m1 + m2  &&  m1 >= ZERO  &&  m2 >= ZERO  ) {
    pstarVector = unitDir1 * Kinematics::pstarTwoBodyDecay(min,m1,m2);
  }
  else if( m1 >= ZERO  &&  m2 >= ZERO && (m1+m2-min)/(min+m1+m2)<eps) {
    pstarVector = Momentum3();
  }
  else {
    throw PhaseSpaceError() << "Two body decay cannot proceed "
				 << "p = " << p / GeV 
				 << " p.m() = " << min / GeV
				 << " -> " << m1/GeV 
				 << ' ' << m2/GeV << Exception::eventerror;
  }
  p1 = Lorentz5Momentum(m1, pstarVector);
  p2 = Lorentz5Momentum(m2,-pstarVector);
  // boost from CM to LAB
  Boost bv = p.boostVector();
  double gammarest = p.e()/p.mass();
  p1.boost( bv , gammarest );
  p2.boost( bv , gammarest );
}
double PhaseSpaceChannel::atanhelper(const PhaseSpaceResonance & res,
				     Energy limit) const {
  return atan2( sqr(limit) - res.mass2, res.mWidth );
}

// return the weight for a given resonance
InvEnergy2 PhaseSpaceChannel::massWeight(const PhaseSpaceResonance & res,
					 Energy moff, Energy lower,
					 Energy upper) const {
  InvEnergy2 wgt = ZERO;
  if(lower>upper) {
    string modestring = mode_->incoming().first->PDGName() + " -> ";
    for(tPDPtr part :mode_->outgoing()) modestring += part->PDGName() + " ";
    throw PhaseSpaceError() << "PhaseSpaceChannel::massWeight not allowed in"
			    << modestring << " "
			    << res.particle->PDGName() << "   " 
			    << moff/GeV << " " << lower/GeV << " " << upper/GeV
			    << Exception::eventerror;
  }
  // use a Breit-Wigner 
  if ( res.jacobian == PhaseSpaceResonance::BreitWigner ) {
    double rhomin  = atanhelper(res,lower);
    double rhomax  = atanhelper(res,upper) - rhomin;
    if ( rhomax != 0.0 ) {
      Energy2 moff2=moff*moff-res.mass2;
      wgt = res.mWidth/rhomax/(moff2*moff2+res.mWidth*res.mWidth);
    }
    else {
      wgt = 1./((sqr(upper)-sqr(lower))*sqr(sqr(moff)-res.mass2)/
 		(sqr(lower)-res.mass2)/(sqr(upper)-res.mass2));
    }
  }
  // power law
  else if(res.jacobian == PhaseSpaceResonance::Power) {
    double rhomin = pow(sqr(lower/MeV),res.power+1.);
    double rhomax = pow(sqr(upper/MeV),res.power+1.)-rhomin;
    wgt = (res.power+1.)/rhomax*pow(sqr(moff/MeV),res.power)
      /MeV/MeV;
  }
  else if(res.jacobian == PhaseSpaceResonance::OnShell ) {
    wgt = 1./Constants::pi/res.mWidth;
  } 
  else {
    throw PhaseSpaceError() << "Unknown type of Jacobian in " 
			    << "PhaseSpaceChannel::massWeight"
			    << Exception::eventerror;
  }
  return wgt;
}

Energy PhaseSpaceChannel::generateMass(const PhaseSpaceResonance & res,
				       Energy lower,Energy upper,
				       const double & rnd) const {
  static const Energy eps=1e-9*MeV;
  if(lower<eps) lower=eps;
  Energy mass=ZERO;
  if(lower>upper) {
    string modestring = mode_->incoming().first->PDGName() + " -> ";
    for(tPDPtr part :mode_->outgoing()) modestring += part->PDGName() + " ";
    throw PhaseSpaceError() << "PhaseSpaceChannel::generateMass"
			    << " not allowed in"
			    << modestring << " "
			    << res.particle->PDGName()
			    << " " << lower/GeV << " " << upper/GeV
			    << Exception::eventerror;
  }
  if(abs(lower-upper)/(lower+upper)>2e-10) {
    lower +=1e-10*(lower+upper);
    upper -=1e-10*(lower+upper);
  }
  else 
    return 0.5*(lower+upper);
  // use a Breit-Wigner
  if(res.jacobian==PhaseSpaceResonance::BreitWigner) {
    if(res.mWidth!=ZERO) {
      Energy2 lower2 = sqr(lower);
      Energy2 upper2 = sqr(upper);
      
      double rhomin = atan2((lower2 - res.mass2),res.mWidth);
      double rhomax = atan2((upper2 - res.mass2),res.mWidth)-rhomin;
      double rho = rhomin+rhomax*rnd;
      Energy2 mass2 = max(lower2,min(upper2,res.mass2+res.mWidth*tan(rho)));
      if(mass2<ZERO) mass2 = ZERO;
      mass = sqrt(mass2);
    }
    else {
      mass = sqrt(res.mass2+
  		  (sqr(lower)-res.mass2)*(sqr(upper)-res.mass2)/
  		  (sqr(lower)-res.mass2-rnd*(sqr(lower)-sqr(upper))));
    }
  }
  // use a power-law
  else if(res.jacobian == PhaseSpaceResonance::Power) {
    double rhomin = pow(sqr(lower/MeV),res.power+1.);
    double rhomax = pow(sqr(upper/MeV),res.power+1.)-rhomin;
    double rho = rhomin+rhomax*rnd;
    mass = pow(rho,0.5/(res.power+1.))*MeV;
  }
  else if(res.jacobian == PhaseSpaceResonance::OnShell) {
    mass = sqrt(res.mass2);
  } 
  else {
    throw PhaseSpaceError() << "Unknown type of Jacobian in " 
			    << "PhaseSpaceChannel::generateMass" 
			    << Exception::eventerror;
  }
  if(mass<lower+1e-10*(lower+upper))      mass=lower+1e-10*(lower+upper);
  else if(mass>upper-1e-10*(lower+upper)) mass=upper-1e-10*(lower+upper);
  return mass;
}

// generate the weight for this channel given a phase space configuration
double PhaseSpaceChannel::generateWeight(const vector<Lorentz5Momentum> & output) const {
  using Constants::pi;
  // include the prefactor due to the weight of the channel
  double wgt = weight_;
  // work out the masses of the intermediate particles
  vector<Energy> intmass;
  for(const PhaseSpaceResonance & res: intermediates_) {
    Lorentz5Momentum pinter;
    for(const int & des : res.descendents) pinter += output[des-1];
    pinter.rescaleMass();
    intmass.push_back( pinter.mass() );
  }
  Energy2 scale(sqr(intmass[0]));
  // calculate the terms for each of the decays
  for(unsigned int ix=0;ix<intermediates_.size();++ix) {
    int idau[2] = {abs(intermediates_[ix].children.first),
		   abs(intermediates_[ix].children.second)};
    // if both decay products off-shell
    if(intermediates_[ix].children.first<0&&intermediates_[ix].children.second<0) {
      // lower limits on the masses of the two resonances
      Energy lowerb[2];
      for(unsigned int iy=0;iy<2;++iy) {
       	lowerb[iy]=ZERO;
   	for(const int & des : intermediates_[idau[iy]].descendents) {
   	  lowerb[iy]+=output[des-1].mass();
  	}
      }
      // undo effect of randomising
      // weight for the first order
      // contribution of first resonance
      Energy upper = intmass[ix]-lowerb[1];
      Energy lower = lowerb[0];
      InvEnergy2 wgta=massWeight(intermediates_[idau[0]],
				 intmass[idau[0]],lower,upper);
      // contribution of second resonance
      upper = intmass[ix]-intmass[idau[0]];
      lower = lowerb[1];
      InvEnergy4 wgta2 = wgta*massWeight(intermediates_[idau[1]],
					 intmass[idau[1]],lower,upper);
      // weight for the second order
      upper = intmass[ix]-lowerb[0];
      lower = lowerb[1];
      InvEnergy2 wgtb=massWeight(intermediates_[idau[1]],
				 intmass[idau[1]],lower,upper);
      upper = intmass[ix]-intmass[idau[1]];
      lower = lowerb[0];
      InvEnergy4 wgtb2=wgtb*massWeight(intermediates_[idau[0]],
				       intmass[idau[0]],lower,upper);
      // weight factor
      wgt *=0.5*sqr(scale)*(wgta2+wgtb2);
      // factor for the kinematics
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[0]],
						 intmass[idau[1]]);
      if(pcm>ZERO)
      	wgt *= intmass[ix]*8.*pi*pi/pcm;
      else
      	wgt = 0.;
    }
    // only first off-shell
    else if(intermediates_[ix].children.first<0) {
      // compute the limits of integration
      Energy upper = intmass[ix]-output[idau[1]-1].mass();
      Energy lower = ZERO;
      for(const int & des : intermediates_[idau[0]].descendents) {
	lower += output[des-1].mass();
      }
      wgt *=scale*massWeight(intermediates_[idau[0]],intmass[idau[0]],lower,upper);
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[0]],
   					  output[idau[1]-1].mass());
      if(pcm>ZERO)
   	wgt *= intmass[ix]*8.*pi*pi/pcm;
      else
   	wgt = 0.;
    }
    // only second off-shell
    else if(intermediates_[ix].children.second<0) {
      // compute the limits of integration
      Energy upper = intmass[ix]-output[idau[0]-1].mass(); 
      Energy lower = ZERO;
      for(const int & des : intermediates_[idau[1]].descendents) {
	lower += output[des-1].mass();
      }
      wgt *=scale*massWeight(intermediates_[idau[1]],intmass[idau[1]],lower,upper);
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],intmass[idau[1]],
						 output[idau[0]-1].mass());
      if(pcm>ZERO)
      	wgt *=intmass[ix]*8.*pi*pi/pcm;
      else
      	wgt=0.;
    }
    // both on-shell
    else {
      Energy pcm = Kinematics::pstarTwoBodyDecay(intmass[ix],output[idau[1]-1].mass(),
						 output[idau[0]-1].mass());
      if(pcm>ZERO)
	wgt *=intmass[ix]*8.*pi*pi/pcm;
      else
	wgt = 0.;
    }
  }
  // finally the overall factor
  wgt /= pi;
  // return the answer
  return wgt;
}

// generate the final-state particles including the intermediate resonances
void PhaseSpaceChannel::generateIntermediates(bool cc, const Particle & in,
					      ParticleVector & out) {
  // create the particles
  // incoming particle
  ParticleVector external;
  external.push_back(const_ptr_cast<tPPtr>(&in));
  // outgoing
  for(unsigned int ix=0;ix<out.size();++ix)
    external.push_back(out[ix]);
  out.clear();
  // now create the intermediates
  ParticleVector resonance;
  resonance.push_back(external[0]);
  // Lorentz5Momentum pinter;
  for(unsigned ix=1;ix<intermediates_.size();++ix) {
    Lorentz5Momentum pinter;
    for(const int & des : intermediates_[ix].descendents)
      pinter += external[des]->momentum();
    pinter.rescaleMass();
    PPtr respart = (cc&&intermediates_[ix].particle->CC()) ? 
      intermediates_[ix].particle->CC()->produceParticle(pinter) : 
      intermediates_[ix].particle      ->produceParticle(pinter);
    resonance.push_back(respart);
  }
  // set up the mother daughter relations
  for(unsigned int ix=1;ix<intermediates_.size();++ix) {
    resonance[ix]->addChild( intermediates_[ix].children.first <0 ? 
   			     resonance[-intermediates_[ix].children.first ] :
			     external[intermediates_[ix].children.first   ]);
    
    resonance[ix]->addChild( intermediates_[ix].children.second<0 ? 
   			     resonance[-intermediates_[ix].children.second] :
			     external[intermediates_[ix].children.second  ]);
    if(resonance[ix]->dataPtr()->stable())
      resonance[ix]->setLifeLength(Lorentz5Distance());
  }
  // construct the output with the particles in the first step
  out.push_back( intermediates_[0].children.first >0 ?
		 external[intermediates_[0].children.first] :
		 resonance[-intermediates_[0].children.first ]);
  out.push_back( intermediates_[0].children.second>0 ?
		 external[intermediates_[0].children.second] :
		 resonance[-intermediates_[0].children.second]);
}
 
ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer PhaseSpaceChannel::createDiagram() const {
  assert(mode_->incoming().second);
  // create the diagram with incoming and s chnnel particle
  ThePEG::Ptr<ThePEG::Tree2toNDiagram>::pointer diag =
    new_ptr((Tree2toNDiagram(2), mode_->incoming().first,
	     mode_->incoming().second,1,intermediates_[0].particle));
  map<int,int> children;
  map<unsigned int,int> res;
  int ires=3;
  res[0]=3;
  // add the intermediates
  for(unsigned int ix=0;ix<intermediates_.size();++ix) {
    if(intermediates_[ix].children.first>0)
      children[intermediates_[ix].children.first]= res[ix];
    else {
      ires+=1;
      diag = new_ptr((*diag,res[ix],intermediates_[-intermediates_[ix].children.first].particle));
      res[-intermediates_[ix].children.first]=ires;
    }
    if(intermediates_[ix].children.second>0)
      children[intermediates_[ix].children.second]= res[ix];
    else {
      ires+=1;
      diag = new_ptr((*diag,res[ix],intermediates_[-intermediates_[ix].children.second].particle));
      res[-intermediates_[ix].children.second]=ires;
    }
  }
  // add the children in the corret order
  for(map<int,int>::const_iterator it=children.begin();it!=children.end();++it) {
    diag = new_ptr((*diag,it->second,mode_->outgoing()[it->first-1]));
  }
  return diag;
}

bool PhaseSpaceChannel::checkKinematics() {
  // recalculate the masses and widths of the resonances
  for(auto inter : intermediates_) {
    inter.mass2 = sqr(inter.particle->mass());
    inter.mWidth = inter.particle->mass()*inter.particle->width();
  }
  Energy massmax = mode_->incoming().first->massMax();
  for(tPDPtr out : mode_->outgoing()) 
    massmax -= out->massMin();
  for(unsigned int ix=1;ix<intermediates_.size();++ix) {
    Energy massmin(ZERO);
    for(const int & iloc : intermediates_[ix].descendents) 
      massmin += mode_->outgoing()[iloc-1]->massMin();
    if(intermediates_[ix].mass2>=sqr(massmin)&&
       intermediates_[ix].mass2<=sqr(massmax+massmin)&&
       intermediates_[ix].mWidth==ZERO)
      return false;
  }
  return true;
}
#line 1 "./PhaseSpaceMode.cc"
// -*- C++ -*-
//
// PhaseSpaceMode.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhaseSpaceMode class.
//

#include "PhaseSpaceMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void PhaseSpaceMode::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << channels_ << maxWeight_ << outgoing_ << outgoingCC_
     << partial_ << widthGen_ << massGen_ << BRsum_ << testOnShell_
     << ounit(eMax_,GeV) << nRand_;
}

void PhaseSpaceMode::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> channels_ >> maxWeight_ >> outgoing_ >> outgoingCC_
     >> partial_ >> widthGen_ >> massGen_ >> BRsum_ >> testOnShell_
     >> iunit(eMax_,GeV) >> nRand_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PhaseSpaceMode,Base>
describeHerwigPhaseSpaceMode("Herwig::PhaseSpaceMode", "libHerwig.so");

void PhaseSpaceMode::Init() {

  static ClassDocumentation<PhaseSpaceMode> documentation
    ("The PhaseSpaceMode classs contains a number of phase space"
     " channels for the integration of a particular decay mode");

}

ParticleVector
PhaseSpaceMode::generateDecay(const Particle & inpart,
			      tcDecayIntegratorPtr decayer,
			      bool intermediates,bool cc) {
  decayer->ME(DecayMEPtr());
  eps_=decayer->eps_;
  // compute the prefactor
  Energy prewid = (widthGen_&&partial_>=0) ?
    widthGen_->partialWidth(partial_,inpart.mass()) :
    incoming_.first->width();
  InvEnergy pre = prewid>ZERO ? 1./prewid : 1./MeV;
  // boosts to/from rest
  Boost bv =-inpart.momentum().boostVector();
  double gammarest = inpart.momentum().e()/inpart.momentum().mass();
  LorentzRotation boostToRest( bv,gammarest);
  LorentzRotation boostFromRest(-bv,gammarest);
  // construct a new particle which is at rest
  Particle inrest(inpart);
  inrest.transform(boostToRest);
  double wgt(0.);
  vector<Lorentz5Momentum> momenta(outgoing_.size());
  int ichan;
  unsigned int ncount(0);
  try {
    do {
      // phase-space pieces of the weight
      fillStack();
      wgt = pre*weight(ichan,inrest.momentum(),momenta);
      // matrix element piece
      wgt *= decayer->me2(-1,inrest,!cc ? outgoing_ : outgoingCC_,
			  momenta,
			  ncount ==0 ? DecayIntegrator::Initialize : DecayIntegrator::Calculate);
      ++ncount;
      if(wgt>maxWeight_) {
        if(decayer->warnings()) {
	CurrentGenerator::log() << "Resetting max weight for decay " 
				<< inrest.PDGName() << " -> ";
	for(tcPDPtr part : outgoing_)
	  CurrentGenerator::log() << "  " << part->PDGName();
	CurrentGenerator::log() << "  " << maxWeight_ << "  " << wgt 
				<< "  " << inrest.mass()/MeV << "\n";
      }
	maxWeight_=wgt;
      }
    }
    while(maxWeight_*UseRandom::rnd()>wgt&&ncount<decayer->nTry_);
  }
  catch (Veto) {
    // restore the incoming particle to its original state
    inrest.transform(boostFromRest);
    while(!rStack_.empty()) rStack_.pop();
    throw Veto();
  }
  if(ncount>=decayer->nTry_) {
    if(decayer->warnings()) {
    CurrentGenerator::log() << "The decay " << inrest.PDGName() << " -> ";
    for(tcPDPtr part : outgoing_)
      CurrentGenerator::log() << "  " << part->PDGName();
    CurrentGenerator::log() << "  " << maxWeight_ << " " << decayer->nTry_
			    << " is too inefficient for the particle "
			    << inpart << "vetoing the decay \n";
    }
    momenta.clear();
    throw Veto();
  }
  // construct the particles
  ParticleVector output(outgoing_.size());
  for(unsigned int ix=0;ix<outgoing_.size();++ix)
    output[ix] = (!cc ? outgoing_[ix] : outgoingCC_[ix] )->produceParticle(momenta[ix]);
  // set up the spin correlations
  decayer->constructSpinInfo(inrest,output);
  const_ptr_cast<tPPtr>(&inpart)->spinInfo(inrest.spinInfo());
  constructVertex(inpart,output,decayer);
  // return if intermediate particles not required
  if(channels_.empty()||!intermediates) {
    for(tPPtr part : output) part->transform(boostFromRest);
  }
  // find the intermediate particles
  else {
    // select the channel
    vector<double> mewgts(channels_.size(),0.0);
    double total=0.;
    for(unsigned int ix=0,N=channels_.size();ix<N;++ix) {
      mewgts[ix]=decayer->me2(ix,inrest,!cc ? outgoing_ : outgoingCC_,
			      momenta,DecayIntegrator::Calculate);
      total+=mewgts[ix];
    }
    // randomly pick a channel
    total *= UseRandom::rnd();
    int iChannel = -1;
    do {
      ++iChannel;
      total-=mewgts[iChannel];
    }
    while(iChannel<int(channels_.size()) && total>0.);
    iChannel_ = iChannel;
    // apply boost
    for(tPPtr part : output) part->transform(boostFromRest);
    channels_[iChannel].generateIntermediates(cc,inpart,output);
  }
  decayer->ME(DecayMEPtr());
  // return particles;
  return output;
}

// flat phase space generation and weight
Energy PhaseSpaceMode::flatPhaseSpace(const Lorentz5Momentum & in,
				      vector<Lorentz5Momentum> & momenta,
				      bool onShell) const {
  double wgt(1.);
  // masses of the particles
  vector<Energy> mass = externalMasses(in.mass(),wgt,onShell);
  // two body decay
  double ctheta  = 2.*rStack_.top() - 1.;
  rStack_.pop();
  double phi  = Constants::twopi*rStack_.top();
  rStack_.pop();
  if(! Kinematics::twoBodyDecay(in, mass[0], mass[1],
				ctheta, phi, momenta[0],momenta[1])) 
    throw Exception() << "Incoming mass - Outgoing mass negative in "
 		      << "PhaseSpaceMode::flatPhaseSpace()"
		      << Exception::eventerror;
  wgt *= Kinematics::pstarTwoBodyDecay(in.mass(),mass[0],mass[1])
    /8./Constants::pi/in.mass();
  return wgt*in.mass();
}

// generate the masses of the external particles
vector<Energy> PhaseSpaceMode::externalMasses(Energy inmass,double & wgt,
					      bool onShell) const {
  vector<Energy> mass(outgoing_.size());
  vector<int> notdone;
  Energy mlow(ZERO);
  // set masses of stable particles and limits 
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    // get the mass of the particle if can't use weight
    if(onShell) {
      mass[ix] = outgoing_[ix]->mass();
      if(massGen_[ix]) rStack_.pop();
    }
    else if(!massGen_[ix] || outgoing_[ix]->stable()) {
      if(massGen_[ix]) rStack_.pop();
      mass[ix] = outgoing_[ix]->generateMass();
      mlow += mass[ix];
    }
    else {
      mass[ix] = ZERO;
      notdone.push_back(ix);
      mlow+=max(outgoing_[ix]->mass()-outgoing_[ix]->widthLoCut(),ZERO);
    }
  }
  if(mlow>inmass) throw Veto();
  // now we need to generate the masses for the particles we haven't
  for( ;!notdone.empty();) {
    unsigned int iloc=long(UseRandom::rnd()*(notdone.size()-1)); 
    Energy low = max(outgoing_[notdone[iloc]]->mass()-outgoing_[notdone[iloc]]->widthLoCut(),ZERO);
    mlow-=low;
    double wgttemp;
    mass[notdone[iloc]]= massGen_[notdone[iloc]]->mass(wgttemp,*outgoing_[notdone[iloc]],low,inmass-mlow,rStack_.top());
    wgttemp /= BRsum_[notdone[iloc]];
    rStack_.pop();
    assert(mass[notdone[iloc]]>=low&&mass[notdone[iloc]]<=inmass-mlow);
    wgt   *= wgttemp;
    mlow  += mass[notdone[iloc]];
    notdone.erase(notdone.begin()+iloc);
  }
  return mass;
}

// construct the vertex for spin corrections
void PhaseSpaceMode::constructVertex(const Particle & inpart, 
				     const ParticleVector & decay,
				     tcDecayIntegratorPtr decayer) const {
  // construct the decay vertex
  VertexPtr vertex(new_ptr(DecayVertex()));
  DVertexPtr Dvertex(dynamic_ptr_cast<DVertexPtr>(vertex));
  // set the incoming particle for the decay vertex
  (inpart.spinInfo())->decayVertex(vertex);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    decay[ix]->spinInfo()->productionVertex(vertex);
  }
  // set the matrix element
  Dvertex->ME(decayer->ME());
  decayer->ME(DecayMEPtr());
}

// initialise the phase space
Energy PhaseSpaceMode::initializePhaseSpace(bool init, tcDecayIntegratorPtr decayer,
					    bool onShell) {
  decayer->ME(DecayMEPtr());
  eps_=decayer->eps_;
  Energy output(ZERO);
  // ensure that the weights add up to one
  if(!channels_.empty()) {
    double temp=0.;
    for(const PhaseSpaceChannel & channel : channels_) temp+= channel.weight();
    for(PhaseSpaceChannel & channel : channels_) channel.weight(channel.weight()/temp);
  }
  if(!init) return ZERO;
  ThePEG::PPtr inpart=incoming_.first->produceParticle();
  // now if using flat phase space
  maxWeight_=0.;
  if(channels_.empty()) {
    vector<Lorentz5Momentum> momenta(outgoing_.size());
    double wsum=0.,wsqsum=0.;
    Energy m0,mmin(ZERO);
    for(tcPDPtr part : outgoing_) mmin += part->massMin();
    for(unsigned int ix=0;ix<decayer->nPoint_;++ix) {
      // set the mass of the decaying particle
      m0 = !onShell ? inpart->dataPtr()->generateMass() : inpart->dataPtr()->mass();
      double wgt=0.;
      if(m0<=mmin) continue;
      inpart->set5Momentum(Lorentz5Momentum(m0));
      // compute the prefactor
      Energy prewid = (widthGen_&&partial_>=0) ?
	widthGen_->partialWidth(partial_,inpart->mass()) :
	incoming_.first->width();
      InvEnergy pre = prewid>ZERO ? 1./prewid : 1./MeV;
      // generate the weight for this point
      try {
	int dummy;
	// phase-space piece
	fillStack();
	wgt = pre*weight(dummy,inpart->momentum(),momenta,onShell);
	// matrix element piece
	wgt *= decayer->me2(-1,*inpart,outgoing_,momenta,DecayIntegrator::Initialize);
      }
      catch (Veto) {
	while(!rStack_.empty()) rStack_.pop();
	wgt=0.;
      }
      if(wgt>maxWeight_) maxWeight_ = wgt;
      wsum   += wgt;
      wsqsum += sqr(wgt);
    }
    wsum /= decayer->nPoint_;
    wsqsum=wsqsum/decayer->nPoint_-sqr(wsum);
    if(wsqsum<0.) wsqsum=0.;
    wsqsum=sqrt(wsqsum/decayer->nPoint_);
    Energy fact = (widthGen_&&partial_>=0) ? 
      widthGen_->partialWidth(partial_,inpart->nominalMass()) :
      inpart->dataPtr()->width();
    if(fact==ZERO) fact=MeV;
     // factor for the weight with spin correlations
    maxWeight_ *= inpart->dataPtr()->iSpin()==1 ? 1.1 : 1.6;
    if ( Debug::level > 1 ) {
      // ouptut the information on the initialisation
      CurrentGenerator::log() << "Initialized the phase space for the decay " 
  			      << incoming_.first->PDGName() << " -> ";
      for(tPDPtr part : outgoing_)
	CurrentGenerator::log() << part->PDGName() << " ";
      CurrentGenerator::log() << "\n";
      if(fact!=MeV) CurrentGenerator::log() << "The branching ratio is " << wsum 
 					    << " +/- " << wsqsum << "\n";
      CurrentGenerator::log() << "The partial width is " << wsum*fact/MeV 
 			      << " +/- " << wsqsum*fact/MeV << " MeV\n";
      CurrentGenerator::log() << "The partial width is " 
 			      << wsum*fact/6.58212E-22/MeV 
 			      << " +/- " << wsqsum*fact/6.58212E-22/MeV<< " s-1\n";
      CurrentGenerator::log() << "The maximum weight is " 
 			      << maxWeight_ << endl;
    }
    output=wsum*fact;
  }
  else {
    vector<Lorentz5Momentum> momenta(outgoing_.size());
    double totsum(0.),totsq(0.);
    for(unsigned int iy=0;iy<decayer->nIter_;++iy) {
      // zero the maximum weight
      maxWeight_=0.;
      vector<double> wsum(channels_.size(),0.),wsqsum(channels_.size(),0.);
      vector<int> nchan(channels_.size(),0);
      totsum = 0.;
      totsq  = 0.;
      Energy mmin(ZERO);
      for(tcPDPtr part : outgoing_) mmin += part->massMin();
      for(unsigned int ix=0;ix<decayer->nPoint_;++ix) {
     	Energy m0 = !onShell ? incoming_.first->generateMass() : incoming_.first->mass();
     	double wgt=0.; 
     	int ichan(-1);
     	if(m0>mmin) {
	  inpart->set5Momentum(Lorentz5Momentum(m0));
	  // compute the prefactor
	  Energy prewid= (widthGen_&&partial_>=0) ? 
	    widthGen_->partialWidth(partial_,inpart->mass()) :
	    inpart->dataPtr()->width();
	  InvEnergy pre = prewid>ZERO ? 1./prewid : 1./MeV;
	  // generate the weight for this point
	  try {
	    fillStack();
	    wgt = pre*weight(ichan,inpart->momentum(),momenta,onShell);
	    // matrix element piece
	    wgt *= decayer->me2(-1,*inpart,outgoing_,momenta,DecayIntegrator::Initialize);
	  }
	  catch (Veto) {
	    wgt=0.;
	  }
     	}
     	if(wgt>maxWeight_) maxWeight_=wgt;
    	if(ichan>=0) {
    	  wsum[ichan]   += wgt;
    	  wsqsum[ichan] += sqr(wgt);
    	  ++nchan[ichan];
    	}
    	totsum+=wgt;
    	totsq+=wgt*wgt;
      }
      totsum=totsum/decayer->nPoint_;
      totsq=totsq/decayer->nPoint_-sqr(totsum);
      if(totsq<0.) totsq=0.;
      totsq=sqrt(totsq/decayer->nPoint_);
      if ( Debug::level > 1 )
    	CurrentGenerator::log() << "The branching ratio is " << iy << " " 
    				<< totsum << " +/- " << totsq 
    				<< maxWeight_ << "\n";
      // compute the individual terms
      double total(0.);
      for(unsigned int ix=0;ix<channels_.size();++ix) {
     	if(nchan[ix]!=0) {
     	  wsum[ix]=wsum[ix]/nchan[ix];
     	  wsqsum[ix]=wsqsum[ix]/nchan[ix];
     	  if(wsqsum[ix]<0.) wsqsum[ix]=0.;
     	  wsqsum[ix]=sqrt(wsqsum[ix]/nchan[ix]);
     	}
     	else {
     	  wsum[ix]=0;
     	  wsqsum[ix]=0;
     	}
     	total+=sqrt(wsqsum[ix])*channels_[ix].weight();
      }
      if(total>0.) {
    	for(unsigned int ix=0;ix<channels_.size();++ix) {
	  channels_[ix].weight(sqrt(wsqsum[ix])*channels_[ix].weight()/total);
     	}
      }
    }
    // factor for the weight with spin correlations
    maxWeight_*= inpart->dataPtr()->iSpin()==1 ? 1.1 : 1.6;
    // output the information on the initialisation
    Energy fact = (widthGen_&&partial_>=0) ? 
      widthGen_->partialWidth(partial_,inpart->nominalMass()) :
      inpart->dataPtr()->width();
    if(fact==ZERO) fact=MeV;
    output=totsum*fact;
    if ( Debug::level > 1 ) {
      CurrentGenerator::log() << "Initialized the phase space for the decay " 
     			      << incoming_.first->PDGName() << " -> ";
      for(tcPDPtr part : outgoing_)
	CurrentGenerator::log() << part->PDGName() << " ";
      CurrentGenerator::log() << "\n";
      if(fact!=MeV) CurrentGenerator::log() << "The branching ratio is " << totsum 
     					    << " +/- " << totsq << "\n";
      CurrentGenerator::log() << "The partial width is " << totsum*fact/MeV 
     			      << " +/- " << totsq*fact/MeV << " MeV\n";
      CurrentGenerator::log() << "The partial width is " 
     			      << totsum*fact/6.58212E-22/MeV 
     			      << " +/- " << totsq*fact/6.58212E-22/MeV 
     			      << " s-1\n";
      CurrentGenerator::log() << "The maximum weight is " 
     			      << maxWeight_ << "\n";
      CurrentGenerator::log() << "The weights for the different phase" 
			      << " space channels are \n";
      for(unsigned int ix=0,N=channels_.size();ix<N;++ix) {
    	CurrentGenerator::log() << "Channel " << ix 
    				<< " had weight " << channels_[ix].weight()
    				<< "\n";
      }
    }
    CurrentGenerator::log() << flush;
  }
  decayer->ME(DecayMEPtr());
  return output;
}
 
// generate a phase-space point using multichannel phase space
Energy PhaseSpaceMode::channelPhaseSpace(int & ichan, const Lorentz5Momentum & in, 
					 vector<Lorentz5Momentum> & momenta,
					 bool onShell) const {
  double wgt(rStack_.top());
  rStack_.pop();
  // select a channel
  ichan=-1;
  do {
    ++ichan;
    wgt-=channels_[ichan].weight();
  }
  while(ichan<int(channels_.size())&&wgt>0.);
  // generate the momenta
  if(ichan==int(channels_.size())) {
    throw Exception() << "PhaseSpaceMode::channelPhaseSpace()"
		      << " failed to select a channel" 
		      << Exception::abortnow;
  }
  // generate the masses of the external particles
  double masswgt(1.);
  vector<Energy> mass(externalMasses(in.mass(),masswgt,onShell));
  momenta=channels_[ichan].generateMomenta(in,mass);
  // compute the denominator of the weight
  wgt=0.;
  for(const PhaseSpaceChannel & channel : channels_) 
    wgt += channel.generateWeight(momenta);
  // return the weight
  return wgt!=0. ? in.mass()*masswgt/wgt : ZERO;
}

void PhaseSpaceMode::init() {
  // get mass and width generators
  outgoingCC_.clear();
  for(tPDPtr part : outgoing_) {
    outgoingCC_.push_back(part->CC() ? part->CC() : part);
  }
  massGen_.resize(outgoing_.size());
  BRsum_.resize(outgoing_.size());
  widthGen_ = dynamic_ptr_cast<cGenericWidthGeneratorPtr>(incoming_.first->widthGenerator());
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    assert(outgoing_[ix]);
    massGen_[ix]= dynamic_ptr_cast<cGenericMassGeneratorPtr>(outgoing_[ix]->massGenerator());
    if(outgoing_[ix]->stable()) {
      BRsum_[ix] = 1.;
    }
    else {
      double total(0.);
      for(tDMPtr mode : outgoing_[ix]->decayModes()) {
	if(mode->on()) total += mode->brat();
      }
      BRsum_[ix] = total;
    }
  }
  // get max energy for decays
  if(!incoming_.second)
    eMax_ = testOnShell_ ? incoming_.first->mass() : incoming_.first->massMax();
  // work out how many random numbers we need
  nRand_ = 3*outgoing_.size()-4;
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    if(massGen_[ix]) ++nRand_;
  }
  if(channels_.empty()) return;
  ++nRand_;
  // ensure weights sum to one
  double sum(0.);
  for(const PhaseSpaceChannel & channel : channels_)
    sum+=channel.weight();
  for(PhaseSpaceChannel & channel : channels_)
    channel.weight(channel.weight()/sum);
}

void PhaseSpaceMode::initrun() {
  // update the mass and width generators
  if(incoming_.first->widthGenerator()!=widthGen_)
    widthGen_ = dynamic_ptr_cast<cGenericWidthGeneratorPtr>(incoming_.first->widthGenerator());
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    if(massGen_[ix]!=outgoing_[ix]->massGenerator())
      massGen_[ix] = dynamic_ptr_cast<cGenericMassGeneratorPtr>(outgoing_[ix]->massGenerator());
  }
  for(PhaseSpaceChannel & channel : channels_) channel.initrun(this);
  if(widthGen_) const_ptr_cast<GenericWidthGeneratorPtr>(widthGen_)->initrun();
  tcGenericWidthGeneratorPtr wtemp;
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    wtemp=
      dynamic_ptr_cast<tcGenericWidthGeneratorPtr>(outgoing_[ix]->widthGenerator());
    if(wtemp) const_ptr_cast<tGenericWidthGeneratorPtr>(wtemp)->initrun();
  }
  // ensure weights sum to one
  double sum(0.);
  for(const PhaseSpaceChannel & channel : channels_)
    sum+=channel.weight();
  for(PhaseSpaceChannel & channel : channels_)
    channel.weight(channel.weight()/sum);
  nRand_ = 3*outgoing_.size()-4;
  for(unsigned int ix=0;ix<outgoing_.size();++ix) {
    if(massGen_[ix]) ++nRand_;
  }
  if(channels_.empty()) return;
  ++nRand_;
}
#line 1 "./HwDecayerBase.cc"
// -*- C++ -*-
//
// HwDecayerBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayerBase class.
//

#include "HwDecayerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/ShowerHandler.h"

using namespace Herwig;

bool HwDecayerBase::accept(const DecayMode & dm) const {
  // get the primary products
  tPDVector products=dm.orderedProducts();
  // add products for which the decay mode is all ready specified
  if(!dm.cascadeProducts().empty()) {
    for(ModeMSet::const_iterator mit=dm.cascadeProducts().begin();
	mit!=dm.cascadeProducts().end();++mit) {
      products.push_back(const_ptr_cast<PDPtr>((**mit).parent()));
    }
  }
  // can this mode be handled ?
  return accept(dm.parent(),products);
}

ParticleVector HwDecayerBase::decay(const DecayMode & dm,
				    const Particle & p) const {
  // handling of the decay including the special features of the
  // DecayMode  
  // get the primary products
  tPDVector products=dm.orderedProducts();
  // add products for which the decay mode is all ready specified
  if(!dm.cascadeProducts().empty()) {
    for(ModeMSet::const_iterator mit=dm.cascadeProducts().begin();
	mit!=dm.cascadeProducts().end();++mit) {
      products.push_back(const_ptr_cast<PDPtr>((**mit).parent()));
    }
  }
  // perform the primary decay
  ParticleVector output=decay(p,products);
  // perform the secondary decays
  if(!dm.cascadeProducts().empty()) {
    unsigned int iloc=dm.orderedProducts().size();
    for(ModeMSet::const_iterator mit=dm.cascadeProducts().begin();
	mit!=dm.cascadeProducts().end();++mit) {
      if(!(*mit)->decayer()) 
	throw Exception() << "Decay mode " << (**mit).tag() 
			  << "does not have a decayer, can't perform"
			  << "decay in  HwDecayerBase::decay()"
			  << Exception::eventerror;
      ParticleVector children=(*mit)->decayer()->decay(**mit,*output[iloc]);
      for(unsigned int ix=0;ix<children.size();++ix) {
	output[iloc]->addChild(children[ix]);
      }
      ++iloc;
    }
  }
  return output;
}

void HwDecayerBase::persistentOutput(PersistentOStream & os) const {
  os << _initialize << _dbOutput;
}

void HwDecayerBase::persistentInput(PersistentIStream & is, int) {
  is >> _initialize >> _dbOutput;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<HwDecayerBase,Decayer>
describeHerwigHwDecayerBase("Herwig::HwDecayerBase", "Herwig.so");

void HwDecayerBase::Init() {

  static ClassDocumentation<HwDecayerBase> documentation
    ("The HwDecayerBase class is the base class for Decayers in Hw++.");

  static Switch<HwDecayerBase,bool> interfaceInitialize
    ("Initialize",
     "Initialization of the phase space calculation",
     &HwDecayerBase::_initialize, false, false, false);
  static SwitchOption interfaceInitializeon
    (interfaceInitialize,
     "Yes",
     "At initialisation find max weight and optimise the integration",
     true);
  static SwitchOption interfaceInitializeoff
    (interfaceInitialize,
     "No",
     "Use the maximum weight and channel weights supplied for the integration",
     false);

  static Switch<HwDecayerBase,bool> interfaceDatabaseOutput
    ("DatabaseOutput",
     "Whether to print the database information",
     &HwDecayerBase::_dbOutput, false, false, false);
  static SwitchOption interfaceDatabaseOutputYes
    (interfaceDatabaseOutput,
     "Yes",
     "Output information on the decayer initialization",
     true);
  static SwitchOption interfaceDatabaseOutputNo
    (interfaceDatabaseOutput,
     "No",
     "Do not output information about the decayer initialization",
     false);
}

void HwDecayerBase::dofinish() {
  Decayer::dofinish();
  if(initialize() && databaseOutput()) {
    string fname = CurrentGenerator::current().filename() + 
      string("-") + name() + string(".output");
    ofstream output(fname.c_str());
    dataBaseOutput(output,true);
  }
}

bool HwDecayerBase::softMatrixElementVeto(PPtr , PPtr,
					  const bool &,
					  const Energy & ,
					  const vector<tcPDPtr> & ,
					  const double & ,
					  const Energy & ,
					  const Energy & ) {
  return false;
}

RealEmissionProcessPtr HwDecayerBase::generateHardest(RealEmissionProcessPtr) {
  return RealEmissionProcessPtr();
}

void HwDecayerBase::initializeMECorrection(RealEmissionProcessPtr , double & ,
			    double & ) {
  assert(false);
}

RealEmissionProcessPtr HwDecayerBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr) {
  assert(false);
  return RealEmissionProcessPtr();
}

void HwDecayerBase::fixRho(RhoDMatrix & rho) const {
  if(ShowerHandler::currentHandlerIsSet() &&
     !ShowerHandler::currentHandler()->spinCorrelations())
    rho.reset();
}
#line 1 "./HwDecayHandler.cc"
// -*- C++ -*-
//
// HwDecayHandler.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwDecayHandler class.
//

#include "HwDecayHandler.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/Decayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DecayIntegrator.h"
#include "PhaseSpaceMode.h"
#include "ThePEG/PDT/MixedParticleData.h"
#include "Herwig/Utilities/EnumParticles.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

void HwDecayHandler::
handle(EventHandler &, const tPVector & tagged,
       const Hint &) {
  // First go through the tagged particles for unstable ones
  tPVector parents;
  for(int i = 0, N = tagged.size(); i<N; ++i) {
    if(tagged[i]) {
      // add to parents if not stable
      if(!tagged[i]->data().stable() &&
	 tagged[i]->data().id() != ParticleID::Remnant &&
	 _excluded.find( tagged[i]->dataPtr() ) == _excluded.end() ) {
	parents.push_back(tagged[i]);
      }
      // if stable and has spinInfo set the developed flag
      else {
	develop(tagged[i]);
      }
    }
  }
  // if nothing to be decayed return
  if(parents.empty()) return;
  useMe();
  // Create a new step, decay all particles and add their children to the step
  StepPtr newstep = _newstep ? newStep() : currentStep();
  for(int i = 0, N = parents.size(); i<N; ++i) {
    performDecay(newstep->find(parents[i]), *newstep);
  }
}

// perform decay method including modifications for spin correlations
// and for the decayer to specify intermediate decay products
void HwDecayHandler::performDecay(tPPtr parent, Step & s) const {
  long ntry = 0;
  tcSpinPtr hwspin;
  tcMixedParticleDataPtr 
    mixdata=dynamic_ptr_cast<tcMixedParticleDataPtr>(parent->dataPtr());
  if(mixdata) {
    pair<bool,Length> mixing = mixdata->generateLifeTime();
    develop(parent);
    parent->setLifeLength(Distance());
    PPtr newparent;
    if(mixing.first) {
      newparent = parent->dataPtr()->CC()->
	produceParticle(parent->momentum());
    }
    else {
      newparent = parent->dataPtr()      ->
	produceParticle(parent->momentum());
    }
    newparent->setLabVertex(parent->labDecayVertex());
    Lorentz5Distance lifeLength(mixing.second,
				parent->momentum().vect()*
				(mixing.second/parent->mass()));
    newparent->setLifeLength(lifeLength);
    s.addDecayProduct(parent, newparent);
    parent = newparent;
  }
  else if ( maxLifeTime() >= ZERO ) {
    if( ( lifeTimeOption() && parent->lifeLength().tau() > maxLifeTime())||
	(!lifeTimeOption() && parent->data().cTau()      > maxLifeTime()) ) {
      parent->setLifeLength(Distance());
      develop(parent);
      return;
    }
  }
  while ( true ) {
    // exit if fails
    if ( ++ntry >= maxLoop() ) 
      throw Exception() << "Too many tries " << maxLoop() << "to generate decay of "
			<< *parent << "in "
			<< "HwDecayHandler::performDecay" << Exception::eventerror;
    // select the decay mode
    tDMPtr   dm(parent->data().selectMode(*parent));
    // check we found a decay mode and it had a decayer
    if ( !dm ) {
      generator()->log() << *generator()->currentEvent() << "\n";
      generator()->log() << *parent << "\n";
      throw Exception() << "No DecayModes for " << parent->PDGName()
			<< " in HwDecayHandler::performDecay" 
			<< Exception::eventerror;
    }
    if ( !dm->decayer() ) throw Exception() << "No decayer for DecayMode of " 
					    << parent->PDGName()
					    << " in HwDecayHandler::performDecay" 
					    << Exception::eventerror;
    try {
      ParticleVector children = dm->decayer()->decay(*dm, *parent);
      if(children.empty()) continue;
      assert(parent->children().empty());
      // generate radiation in the decay
      tDecayIntegratorPtr hwdec=dynamic_ptr_cast<tDecayIntegratorPtr>(dm->decayer());
      if (hwdec && hwdec->canGeneratePhotons())
	children = hwdec->generatePhotons(*parent,children);
      // set up parent
      parent->decayMode(dm);
      // add children
      for ( int i = 0, N = children.size(); i < N; ++i ) {
	children[i]->setLabVertex(parent->labDecayVertex());
	if ( !s.addDecayProduct(parent, children[i]) ) 
	  throw Exception() << "Failed to add child " 
			    << children[i]->PDGName() 
			    << " in decay of " << parent->PDGName() 
			    << Exception::eventerror;
      }
      parent->scale(ZERO);
      // loop over the children
      for ( int i = 0, N = children.size(); i < N; ++i ) {
	// if the child has already been decayed add products to the record
	if(children[i]->decayed()) addDecayedParticle(children[i],s);
	// if not stable decay the child
	else if (!children[i]->data().stable() &&
		 _excluded.find( children[i]->dataPtr() ) == _excluded.end() ) {
	  performDecay(children[i], s);
	}
	// if stable and has spinInfo set up decay matrices etc.
	else {
	  develop(children[i]);
	}
      }
      // sort out the spinInfo for the parent after the decays
      if(parent->spinInfo()) parent->spinInfo()->develop();
      return;
    }
    catch (Veto) 
      {}
  }
}

// method to add an intermediate which has already been decayed to the event record
void HwDecayHandler::addDecayedParticle(tPPtr parent, Step & s) const {
  for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
    parent->children()[i]->setLabVertex(parent->labDecayVertex());
    s.addDecayProduct(parent->children()[i]);
  }
  parent->scale(ZERO);
  for ( int i = 0, N = parent->children().size(); i < N; ++i ) {
    if((parent->children()[i])->decayed()) {
      for(unsigned int ix=0;ix<(parent->children()[i])->children().size();++ix)
	addDecayedParticle(parent->children()[i],s);
    }
    else if ( ! parent->children()[i]->data().stable() &&
	      _excluded.find( parent->children()[i]->dataPtr() ) == _excluded.end() ) {
      performDecay(parent->children()[i], s);
    }
    else {
      develop(parent->children()[i]);
    }
  }
  return;
}

void HwDecayHandler::persistentOutput(PersistentOStream & os) const {
  os << _newstep << _excluded << _excludedVector;
}

void HwDecayHandler::persistentInput(PersistentIStream & is, int)  {
  is >> _newstep >> _excluded >> _excludedVector;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HwDecayHandler,DecayHandler>
describeHerwigHwDecayHandler("Herwig::HwDecayHandler", "Herwig.so");

void HwDecayHandler::Init() {

  static ClassDocumentation<HwDecayHandler> documentation
    ("This is the handler for decays in Herwig.",
     "Decays in Herwig include full spin correlations, based on \\cite{Richardson:2001df}.",
     "%\\cite{Richardson:2001df}\n"
     "\\bibitem{Richardson:2001df}\n"
     "  P.~Richardson,\n"
     "  ``Spin correlations in Monte Carlo simulations,''\n"
     "  JHEP {\\bf 0111}, 029 (2001)\n"
     "  [arXiv:hep-ph/0110108].\n"
     "  %%CITATION = JHEPA,0111,029;%%\n"
     );

  static Switch<HwDecayHandler,bool> interfaceNewStep
    ("NewStep",
     "Add the particles in a new step",
     &HwDecayHandler::_newstep, true, false, false);
  static SwitchOption interfaceNewStepNew
    (interfaceNewStep,
     "Yes",
     "Add particles in a new step",
     true);
  static SwitchOption interfaceNewStepCurrent
    (interfaceNewStep,
     "No",
     "Add them in the current step",
     false);

  static RefVector<HwDecayHandler,ParticleData> interfaceExcluded
    ("Excluded",
     "Particles which should not be decayed",
     &HwDecayHandler::_excludedVector, -1, false, false, true, false, false);

}

void HwDecayHandler::doinit() {
  DecayHandler::doinit();
  _excluded = set<tcPDPtr>(_excludedVector.begin(),_excludedVector.end());
}
#line 1 "./DecayVertex.cc"
// -*- C++ -*-
//
// DecayVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayVertex class.
//
//  Author: Peter Richardson
//

#include <ThePEG/EventRecord/SpinInfo.h>
#include "DecayVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using namespace ThePEG;

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<DecayVertex,HelicityVertex>
describeHerwigDecayVertex("Herwig::DecayVertex", "Herwig.so");
    
void DecayVertex::Init() {
  
  static ClassDocumentation<DecayVertex> documentation
    ("The DecayVertex class is the implementation of a "
     "vertex for a decay for the Herwig spin correlation algorithm");
  
}

// method to get the rho matrix for a given outgoing particle
RhoDMatrix DecayVertex::getRhoMatrix(int i,bool recursive) const {
  // get the rho matrix of the decaying particle
  RhoDMatrix input;
  tcSpinPtr inspin = incoming()[0];
  assert(inspin);
  if(recursive&&inspin->productionVertex()&&
     inspin->iSpin()!=PDT::Spin0) {
    input = inspin->productionVertex()->
      getRhoMatrix(inspin->productionLocation(),true);
    inspin->rhoMatrix() = input;
    inspin->needsUpdate();
  }
  else {
    input = inspin->rhoMatrix();
  }
  // get the D matrices for the outgoing particles
  vector<RhoDMatrix> rhoout(outgoing().size()-1);
  for(int ix=0,N=outgoing().size();ix<N;++ix) {
    if(ix<i)      
      rhoout[ix]   = outgoing()[ix]->DMatrix();
    else if(ix>i) 
      rhoout[ix-1] = outgoing()[ix]->DMatrix();
  }
  // calculate the spin density matrix
  return matrixElement_->calculateRhoMatrix(i,input,rhoout);
}

// method to get the D matrix for an incoming particle
RhoDMatrix DecayVertex::getDMatrix(int) const {
  tcSpinPtr inspin = incoming()[0];
  if(inspin->developed()==SpinInfo::Developed) 
    return inspin->DMatrix();
  // get the decay matrices for the outgoing particles
  vector<RhoDMatrix> Dout(outgoing().size());
  for(unsigned int ix=0,N=outgoing().size();ix<N;++ix) {
    tcSpinPtr hwspin = outgoing()[ix];
    if(hwspin->developed()!=SpinInfo::Developed) 
      hwspin->develop();
    Dout[ix] = hwspin->DMatrix();
  }
  // calculate the spin density matrix and return the answer
  return matrixElement_->calculateDMatrix(Dout);
}
#line 1 "./DecayMatrixElement.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayMatrixElement class.
//

#include "DecayMatrixElement.h"

using namespace Herwig;

DecayMatrixElement::~DecayMatrixElement() {}

#line 1 "./TwoBodyDecayMatrixElement.cc"
// -*- C++ -*-
//
// TwoBodyMatrixElement.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyMatrixElement class.
//
// Author: Peter Richardson
//

#include "TwoBodyDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG;
    
// calculate the decay matrix for this decay
RhoDMatrix TwoBodyDecayMatrixElement::
calculateDMatrix(const vector<RhoDMatrix> & rhoout) const {
  // rhomatrix to be returned
  RhoDMatrix output(inspin(), false);
  // loop over all helicity components of the matrix element
  for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
    for(int jhel1=0;jhel1<outspin()[0];++jhel1) {
      if(rhoout[0](ihel1,jhel1)==0.)continue;
      for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
        for(int jhel2=0;jhel2<outspin()[1];++jhel2) {
          if(rhoout[1](ihel2,jhel2)==0.)continue;
  	  for(int ihel0=0;ihel0<inspin();++ihel0) {
            for(int jhel0=0;jhel0<inspin();++jhel0) {
      
	      output(ihel0,jhel0) += 
		     matrixElement_[ihel0][ihel1][ihel2]*
		conj(matrixElement_[jhel0][jhel1][jhel2])*
		rhoout[0](ihel1,jhel1)*rhoout[1](ihel2,jhel2);
	    }
	  }
	}
      }
    }
  }
  // ensure unit trace for the matrix
  output.normalize();
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix TwoBodyDecayMatrixElement::
calculateRhoMatrix(int id,const RhoDMatrix & rhoin,
		   const vector<RhoDMatrix> & rhoout) const {
  // rhomatrix to be returned
  RhoDMatrix output(outspin()[id], false);
  // loop over all helicity components of the matrix element
  if(id==0) {
    for(int ihel0=0;ihel0<inspin();++ihel0) {
      for(int jhel0=0;jhel0<inspin();++jhel0) {
        if (rhoin(ihel0,jhel0)==0.)continue;
        for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
          for(int jhel2=0;jhel2<outspin()[1];++jhel2) {
            if(rhoout[0](ihel2,jhel2)==0.)continue;
            for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
              if(matrixElement_[ihel0][ihel1][ihel2]==0.)continue;
              for(int jhel1=0;jhel1<outspin()[0];++jhel1) {
                output(ihel1,jhel1) +=
                matrixElement_[ihel0][ihel1][ihel2]*
                conj(matrixElement_[jhel0][jhel1][jhel2])*
                rhoin(ihel0,jhel0)*
                rhoout[0](ihel2,jhel2);
              }
            }
          }
        }
      }
    }
  }
  else {
    for(int ihel0=0;ihel0<inspin();++ihel0) {
      for(int jhel0=0;jhel0<inspin();++jhel0) {
        if (rhoin(ihel0,jhel0)==0.)continue;
	for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
	  for(int jhel1=0;jhel1<outspin()[0];++jhel1) {
            if (rhoout[0](ihel1,jhel1)==0.)continue;
	    for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
	      for(int jhel2=0;jhel2<outspin()[1];++jhel2) {
		output(ihel2,jhel2) += 
		  matrixElement_[ihel0][ihel1][ihel2]*
		  conj(matrixElement_[jhel0][jhel1][jhel2])*
		  rhoin(ihel0,jhel0)*
		  rhoout[0](ihel1,jhel1);
	      }
	    }
	  }
	}
      }
    }
  }
  // return the answer
  output.normalize();
  return output;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex TwoBodyDecayMatrixElement::contract(const RhoDMatrix & in) const {
  Complex me=0.;
  for(int ihel0=0;ihel0<inspin();++ihel0) {
    for(int jhel0=0;jhel0<inspin();++jhel0) {
      for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
	  for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
	      me +=  matrixElement_[ihel0][ihel1][ihel2]*
		conj(matrixElement_[jhel0][ihel1][ihel2])*in(ihel0,jhel0);
	  }
      }
    }
  }
  return me;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex TwoBodyDecayMatrixElement::contract(const TwoBodyDecayMatrixElement & con, 
					    const RhoDMatrix & in) {
  Complex me=0.;
  for(int ihel0=0;ihel0<inspin();++ihel0) {
    for(int jhel0=0;jhel0<inspin();++jhel0) {
      for(int ihel1=0;ihel1<outspin()[0];++ihel1) {
	for(int ihel2=0;ihel2<outspin()[1];++ihel2) {
	  // compute the term
	  me += matrixElement_[ihel0][ihel1][ihel2]*
	    conj(con.matrixElement_[jhel0][ihel1][ihel2])*in(ihel0,jhel0);
	}
      }
    }
  }
  return me;
}
#line 1 "./GeneralDecayMatrixElement.cc"
// -*- C++ -*-
//
// GeneralDecayMatrixElement.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GeneralDecayMatrixElement class.
//
// Author: Peter Richardson
//

#include "GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG;
    
// calculate the decay matrix for this decay
RhoDMatrix GeneralDecayMatrixElement::calculateDMatrix(const vector<RhoDMatrix> & rhoout) const {
  // vectors for the helicities
  vector<int> ihel1(outspin().size()+1),ihel2(outspin().size()+1);
  // rhomatrix to be returned
  RhoDMatrix output(inspin(), false);
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  unsigned int ix,iy,iz;
  int ixa,iya;
  for(ix=0;ix<matrixElement_.size();++ix) {
    // map the vector index to the helicities
    for(ixa=outspin().size();ixa>=0;--ixa) 
      ihel1[ixa]=(ix%constants_[ixa])/constants_[ixa+1];
    // inner loop
    for(iy=0;iy<matrixElement_.size();++iy) {
      // map the vector index to the helicities	   
      for(iya=outspin().size();iya>=0;--iya)
	ihel2[iya]=(iy%constants_[iya])/constants_[iya+1];
      // matrix element piece
      temp=matrixElement_[ix]*conj(matrixElement_[iy]);
      // spin density matrices for the outgoing particles
      for(iz=0;iz<outspin().size();++iz)
	temp*=rhoout[iz](ihel1[iz+1],ihel2[iz+1]);
      output(ihel1[0],ihel2[0])+=temp;
    }
  }
  // ensure unit trace for the matrix
  output.normalize();
  // return the answer
  return output;
}

// calculate the rho matrix for a given outgoing particle
RhoDMatrix GeneralDecayMatrixElement::
calculateRhoMatrix(int id,const RhoDMatrix & rhoin,
		   const vector<RhoDMatrix> & rhoout) const {
  // vectors for the helicities
  vector<int> ihel1(outspin().size()+1),ihel2(outspin().size()+1);
  // rhomatrix to be returned
  RhoDMatrix output(outspin()[id], false);
  // loop over all helicity components of the matrix element
  // outer loop
  Complex temp;
  unsigned int ix,iy,iz;
  int ixa,iya;
  for(ix=0;ix<matrixElement_.size();++ix) {
    // map the vector index to the helicities
    for(ixa=outspin().size();ixa>=0;--ixa)
      ihel1[ixa]=(ix%constants_[ixa])/constants_[ixa+1];
    // inner loop
    for(iy=0;iy<matrixElement_.size();++iy) {
      // map the vector index to the helicities	   
      for(iya=outspin().size();iya>=0;--iya)
	ihel2[iya]=(iy%constants_[iya])/constants_[iya+1];
      // matrix element piece
      temp=matrixElement_[ix]*conj(matrixElement_[iy]);
      // spin denisty matrix for the incoming particle
      temp *= rhoin(ihel1[0],ihel2[0]);
      // spin density matrix for the outgoing particles
      for(iz=0;iz<outspin().size()-1;++iz) {
	if(int(iz)<id) temp*=rhoout[iz](ihel1[iz+1],ihel2[iz+1]);
	else           temp*=rhoout[iz](ihel1[iz+2],ihel2[iz+2]);
      }
      // add to the rho matrix
      output(ihel1[id+1],ihel2[id+1])+=temp;
    }
  }
  // return the answer
  output.normalize();
  return output;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex GeneralDecayMatrixElement::contract(const RhoDMatrix & in) const {
  unsigned int ispin(abs(int(inspin())));
  Complex me=0.;
  for(unsigned int ix=0;ix<constants_[1];++ix) {
    for(unsigned int inhel1=0;inhel1<ispin;++inhel1) {
      for(unsigned int inhel2=0;inhel2<ispin;++inhel2) {
	// compute the term
	me+=matrixElement_[inhel1*constants_[1]+ix]*
	  conj(matrixElement_[inhel2*constants_[1]+ix])*in(inhel1,inhel2);
      }
    }
  }
  return me;
}

// contract the matrix element with the rho matrix of the incoming particle
Complex GeneralDecayMatrixElement::contract(const GeneralDecayMatrixElement & con, 
					    const RhoDMatrix & in) {
  unsigned int ispin(abs(int(inspin())));
  Complex me=0.;
  unsigned int ix,inhel1,inhel2;
  for(ix=0;ix<constants_[1];++ix) {
    for(inhel1=0;inhel1<ispin;++inhel1) {
      for(inhel2=0;inhel2<ispin;++inhel2) {
	// compute the term
	me+=matrixElement_[inhel1*constants_[1]+ix]*
	  conj(con.matrixElement_[inhel2*constants_[1]+ix])*in(inhel1,inhel2);
      }
    }
  }
  return me;
}
#line 1 "./BranchingRatioReweighter.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BranchingRatioReweighter class.
//

#include "BranchingRatioReweighter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/EventRecord/StandardSelectors.h"
#include "ThePEG/Handlers/StandardEventHandler.h"
#include "Herwig/Utilities/EnumParticles.h"

using namespace Herwig;

BranchingRatioReweighter::BranchingRatioReweighter() {}

void BranchingRatioReweighter::
handle(EventHandler & eh,const tPVector & ,const Hint & ) {
  tEventPtr event = eh.currentEvent();
  // weight
  double weight = 1.;
  // get all the particles
  set<tcPPtr> particles;
  event->select(inserter(particles),ThePEG::AllSelector());
  for(set<tcPPtr>::const_iterator it=particles.begin();it!=particles.end();++it) {
    // skip stable
    if((**it).dataPtr()->stable()) continue;
    // skip remnant and clusters
    if((**it).id()==ParticleID::Remnant || 
       (**it).id()==ParticleID::Cluster) continue;
    // if spacelike skip
    if((**it).mass()<ZERO) continue;
    if(*it == event->incoming().first || 
       *it == event->incoming().second ) continue;
    // find unique particles
    bool unique = true;
    for(unsigned int ix=0;ix<(**it).children().size();++ix) {
      if((**it).children()[ix]->id()==(**it).id()) {
	unique = false;
	break;
      }
    }
    if(!unique) continue;
    weight *= (**it).dataPtr()->decaySelector().sum();
  }
  // do the reweighting
  if ( dynamic_cast<StandardEventHandler*>(&eh) ) {
    StandardEventHandler& seh = 
      dynamic_cast<StandardEventHandler&>(eh);
    seh.reweight(weight);
  }
}

IBPtr BranchingRatioReweighter::clone() const {
  return new_ptr(*this);
}

IBPtr BranchingRatioReweighter::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type description system in ThePEG.
DescribeNoPIOClass<BranchingRatioReweighter,StepHandler>
describeHerwigBranchingRatioReweighter("Herwig::BranchingRatioReweighter", "Herwig.so");

void BranchingRatioReweighter::Init() {

  static ClassDocumentation<BranchingRatioReweighter> documentation
    ("The BranchingRatioReweighter class reweights events if some"
     " decay modes are switched off");

}

#line 1 "./PerturbativeDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PerturbativeDecayer class.
//

#include "PerturbativeDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace Herwig;

void PerturbativeDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(pTmin_,GeV) << oenum(inter_) << alphaS_ << alphaEM_
     << useMEforT2_ << C_ << ymax_ << phaseOpt_;
}

void PerturbativeDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(pTmin_,GeV) >> ienum(inter_) >> alphaS_ >> alphaEM_
     >> useMEforT2_ >> C_ >> ymax_ >> phaseOpt_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<PerturbativeDecayer,DecayIntegrator>
describeHerwigPerturbativeDecayer("Herwig::PerturbativeDecayer",
				  "Herwig.so HwPerturbativeDecay.so");

void PerturbativeDecayer::Init() {

  static ClassDocumentation<PerturbativeDecayer> documentation
    ("The PerturbativeDecayer class is the mase class for "
     "perturbative decays in Herwig");

  static Parameter<PerturbativeDecayer,Energy> interfacepTmin
    ("pTmin",
     "Minimum transverse momentum from gluon radiation",
     &PerturbativeDecayer::pTmin_, GeV, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Switch<PerturbativeDecayer,ShowerInteraction> interfaceInteractions
    ("Interactions",
     "which interactions to include for the hard corrections",
     &PerturbativeDecayer::inter_, ShowerInteraction::QCD, false, false);
  static SwitchOption interfaceInteractionsQCD
    (interfaceInteractions,
     "QCD",
     "QCD Only",
     ShowerInteraction::QCD);
  static SwitchOption interfaceInteractionsQED
    (interfaceInteractions,
     "QED",
     "QED only",
     ShowerInteraction::QED);
  static SwitchOption interfaceInteractionsQCDandQED
    (interfaceInteractions,
     "QCDandQED",
     "Both QCD and QED",
     ShowerInteraction::QEDQCD);

  static Reference<PerturbativeDecayer,ShowerAlpha> interfaceAlphaS
    ("AlphaS",
     "Object for the coupling in the generation of hard QCD radiation",
     &PerturbativeDecayer::alphaS_, false, false, true, true, false);

  static Reference<PerturbativeDecayer,ShowerAlpha> interfaceAlphaEM
    ("AlphaEM",
     "Object for the coupling in the generation of hard QED radiation",
     &PerturbativeDecayer::alphaEM_, false, false, true, true, false);

  static Switch<PerturbativeDecayer,bool> interfaceUseMEForT2
    ("UseMEForT2",
     "Use the matrix element correction, if available to fill the T2"
     " region for the decay shower and don't fill using the shower",
     &PerturbativeDecayer::useMEforT2_, true, false, false);
  static SwitchOption interfaceUseMEForT2Shower
    (interfaceUseMEForT2,
     "Shower",
     "Use the shower to fill the T2 region",
     false);
  static SwitchOption interfaceUseMEForT2ME
    (interfaceUseMEForT2,
     "ME",
     "Use the Matrix element to fill the T2 region",
     true);

    static Parameter<PerturbativeDecayer,double> interfacePrefactor
    ("Prefactor",
     "The prefactor for the sampling of the powheg Sudakov",
     &PerturbativeDecayer::C_, 6.3, 0.0, 1e10,
     false, false, Interface::limited);

  static Parameter<PerturbativeDecayer,double> interfaceYMax
    ("YMax",
     "The maximum value for the rapidity",
     &PerturbativeDecayer::ymax_, 10., 0.0, 100.,
     false, false, Interface::limited);

  static Switch<PerturbativeDecayer,unsigned int> interfacePhaseSpaceOption
    ("PhaseSpaceOption",
     "Option for the phase-space sampling",
     &PerturbativeDecayer::phaseOpt_, 0, false, false);
  static SwitchOption interfacePhaseSpaceOptionFixedYLimits
    (interfacePhaseSpaceOption,
     "FixedYLimits",
     "Use a fixed limit for the rapidity",
     0);
  static SwitchOption interfacePhaseSpaceOptionVariableYLimits
    (interfacePhaseSpaceOption,
     "VariableYLimits",
     "Change limit for the rapidity with pT",
     1);

}

double PerturbativeDecayer::matrixElementRatio(const Particle & , 
					       const ParticleVector & ,
					       const ParticleVector & , 
					       MEOption ,
					       ShowerInteraction ) {
  throw Exception() << "Base class PerturbativeDecayer::matrixElementRatio() "
		    << "called, should have an implementation in the inheriting class"
		    << Exception::runerror;
  return 0.;
}

RealEmissionProcessPtr PerturbativeDecayer::generateHardest(RealEmissionProcessPtr born) {
  return getHardEvent(born,false,inter_);
}

RealEmissionProcessPtr PerturbativeDecayer::applyHardMatrixElementCorrection(RealEmissionProcessPtr born) {
  return getHardEvent(born,true,ShowerInteraction::QCD);
}

RealEmissionProcessPtr PerturbativeDecayer::getHardEvent(RealEmissionProcessPtr born,
							 bool inDeadZone,
							 ShowerInteraction inter) {
  // check one incoming
  assert(born->bornIncoming().size()==1);
  // search for coloured/charged particles
  bool colouredParticles=born->bornIncoming()[0]->dataPtr()->coloured();
  bool chargedParticles=born->bornIncoming()[0]->dataPtr()->charged();
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(born->bornOutgoing()[ix]->dataPtr()->coloured())
      colouredParticles=true;
    if(born->bornOutgoing()[ix]->dataPtr()->charged())
      chargedParticles=true;
  }
  // if no coloured/charged particles return
  if ( !colouredParticles && !chargedParticles ) return RealEmissionProcessPtr();
  if ( !colouredParticles && inter==ShowerInteraction::QCD ) return RealEmissionProcessPtr();
  if ( ! chargedParticles && inter==ShowerInteraction::QED ) return RealEmissionProcessPtr();
  // check exactly two outgoing particles
  if(born->bornOutgoing().size()!=2) return RealEmissionProcessPtr();
  // for decay b -> a c 
  // set progenitors
  PPtr cProgenitor = born->bornOutgoing()[0];
  PPtr aProgenitor = born->bornOutgoing()[1];
  // get the decaying particle
  PPtr bProgenitor = born->bornIncoming()[0];
  // identify which dipoles are required
  vector<DipoleType> dipoles;
  if(!identifyDipoles(dipoles,aProgenitor,bProgenitor,cProgenitor,inter)) {
    return RealEmissionProcessPtr();
  }
  Energy trialpT = pTmin_;
  LorentzRotation eventFrame;
  vector<Lorentz5Momentum> momenta;
  vector<Lorentz5Momentum> trialMomenta(4);
  PPtr finalEmitter, finalSpectator;
  PPtr trialEmitter, trialSpectator;
  DipoleType finalType(FFa,ShowerInteraction::QCD);
  for (int i=0; i<int(dipoles.size()); ++i) {
    if(dipoles[i].type==FFg) continue;
    // assign emitter and spectator based on current dipole
    if (dipoles[i].type==FFc || dipoles[i].type==IFc || dipoles[i].type==IFbc) {
      trialEmitter   = cProgenitor;
      trialSpectator = aProgenitor;
    }
    else if (dipoles[i].type==FFa || dipoles[i].type==IFa || dipoles[i].type==IFba) {
      trialEmitter   = aProgenitor;
      trialSpectator = cProgenitor;
    }
    
    // find rotation from lab to frame with the spectator along -z
    LorentzRotation trialEventFrame(bProgenitor->momentum().findBoostToCM());
    Lorentz5Momentum pspectator = (trialEventFrame*trialSpectator->momentum());
    trialEventFrame.rotateZ( -pspectator.phi() );
    trialEventFrame.rotateY( -pspectator.theta() - Constants::pi );
    // invert it
    trialEventFrame.invert();
    // try to generate an emission
    pT_ = pTmin_;
    vector<Lorentz5Momentum> trialMomenta 
      = hardMomenta(bProgenitor, trialEmitter, trialSpectator,
		    dipoles, i, inDeadZone);
    // select dipole which gives highest pT emission
    if(pT_>trialpT) {
      trialpT        = pT_;
      momenta        = trialMomenta;
      eventFrame     = trialEventFrame;
      finalEmitter   = trialEmitter;
      finalSpectator = trialSpectator;
      finalType      = dipoles[i];
      if (dipoles[i].type==FFc || dipoles[i].type==FFa ) {
      	if((momenta[3]+momenta[1]).m2()-momenta[1].m2()>
	   (momenta[3]+momenta[2]).m2()-momenta[2].m2()) {
      	  swap(finalEmitter,finalSpectator);
      	  swap(momenta[1],momenta[2]);
      	}
      }
    }
  }
  pT_ = trialpT;
  // if no emission return
  if(momenta.empty()) {
    if(inter==ShowerInteraction::ALL || inter==ShowerInteraction::QEDQCD || inter==ShowerInteraction::QCD)
      born->pT()[ShowerInteraction::QCD] = pTmin_;
    if(inter==ShowerInteraction::ALL || inter==ShowerInteraction::QEDQCD || inter==ShowerInteraction::QED)
      born->pT()[ShowerInteraction::QED] = pTmin_;
    return born;
  }
  // rotate momenta back to the lab
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    momenta[ix] *= eventFrame;
  }
 
  // set maximum pT for subsequent branchings
  if(inter==ShowerInteraction::ALL || inter==ShowerInteraction::QEDQCD || inter==ShowerInteraction::QCD)
    born->pT()[ShowerInteraction::QCD] = pT_;
  if(inter==ShowerInteraction::ALL || inter==ShowerInteraction::QEDQCD || inter==ShowerInteraction::QED)
    born->pT()[ShowerInteraction::QED] = pT_;

  // get ParticleData objects
  tcPDPtr b = bProgenitor   ->dataPtr();
  tcPDPtr e = finalEmitter  ->dataPtr();
  tcPDPtr s = finalSpectator->dataPtr();
  
  tcPDPtr boson  = getParticleData(finalType.interaction==ShowerInteraction::QCD ?
				   ParticleID::g : ParticleID::gamma);

  // create new ShowerParticles
  PPtr emitter   = e    ->produceParticle(momenta[1]);
  PPtr spectator = s    ->produceParticle(momenta[2]);
  PPtr gauge     = boson->produceParticle(momenta[3]);
  PPtr incoming  = b    ->produceParticle(bProgenitor->momentum());

  // insert the particles
  born->incoming().push_back(incoming);
  unsigned int iemit(0),ispect(0);
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    if(born->bornOutgoing()[ix]==finalEmitter) {
      born->outgoing().push_back(emitter);
      iemit = born->outgoing().size();
    }
    else if(born->bornOutgoing()[ix]==finalSpectator) {
      born->outgoing().push_back(spectator);
      ispect = born->outgoing().size();
    }
  }
  born->outgoing().push_back(gauge);
  if(!spectator->dataPtr()->coloured() ||
     (finalType.type != FFa && finalType.type!=FFc) ) ispect = 0;
  born->emitter(iemit);
  born->spectator(ispect);
  born->emitted(3);
  // boost if being use as ME correction
  if(inDeadZone) {
    if(finalType.type==IFa || finalType.type==IFba) {
      LorentzRotation trans(cProgenitor->momentum().findBoostToCM());
      trans.boost(spectator->momentum().boostVector());
      born->transformation(trans);
    }
    else if(finalType.type==IFc || finalType.type==IFbc) {
      LorentzRotation trans(bProgenitor->momentum().findBoostToCM());
      trans.boost(spectator->momentum().boostVector());
      born->transformation(trans);
    }
  }
  // set the interaction
  born->interaction(finalType.interaction);
  // set up colour lines
  getColourLines(born);
  // return the tree
  return born;
}

bool PerturbativeDecayer::identifyDipoles(vector<DipoleType>  & dipoles,
					  PPtr & aProgenitor,
					  PPtr & bProgenitor,
					  PPtr & cProgenitor,
					  ShowerInteraction inter) const {
  enhance_ = 1.;
  // identify any QCD dipoles
  if(inter==ShowerInteraction::QCD ||
     inter==ShowerInteraction::ALL || inter==ShowerInteraction::QEDQCD ) {
    PDT::Colour bColour = bProgenitor->dataPtr()->iColour();
    PDT::Colour cColour = cProgenitor->dataPtr()->iColour();
    PDT::Colour aColour = aProgenitor->dataPtr()->iColour();
    
    // decaying colour singlet
    if    (bColour==PDT::Colour0 ) {
      if ((cColour==PDT::Colour3    && aColour==PDT::Colour3bar) ||
	  (cColour==PDT::Colour3bar && aColour==PDT::Colour3)    ||
	  (cColour==PDT::Colour8    && aColour==PDT::Colour8)){
	dipoles.push_back(DipoleType(FFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc,ShowerInteraction::QCD));
	if(aProgenitor->id()==ParticleID::g &&
	   cProgenitor->id()==ParticleID::g ) {
	  enhance_ = 1.5;
	  dipoles.push_back(DipoleType(FFg,ShowerInteraction::QCD));
	}
      }
    }
    // decaying colour triplet
    else if (bColour==PDT::Colour3 ) {
      if (cColour==PDT::Colour3 && aColour==PDT::Colour0){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour0 && aColour==PDT::Colour3){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour8 && aColour==PDT::Colour3){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour3 && aColour==PDT::Colour8){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if(cColour==PDT::Colour3bar && aColour==PDT::Colour3bar) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
    }
    // decaying colour anti-triplet 
    else if (bColour==PDT::Colour3bar) {
      if ((cColour==PDT::Colour3bar && aColour==PDT::Colour0)){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
      }
      else if ((cColour==PDT::Colour0 && aColour==PDT::Colour3bar)){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));      
      }
      else if (cColour==PDT::Colour8 && aColour==PDT::Colour3bar){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour3bar && aColour==PDT::Colour8){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
      else if(cColour==PDT::Colour3 && aColour==PDT::Colour3) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFc ,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(FFa ,ShowerInteraction::QCD));
      }
    }
    // decaying colour octet
    else if (bColour==PDT::Colour8){
      if ((cColour==PDT::Colour3    && aColour==PDT::Colour3bar) ||
	  (cColour==PDT::Colour3bar && aColour==PDT::Colour3)){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour8 && aColour==PDT::Colour0){
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
      else if (cColour==PDT::Colour0 && aColour==PDT::Colour8){
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
      }
    }
    // decaying colour sextet
    else if(bColour==PDT::Colour6) {
      if (cColour==PDT::Colour3 && aColour==PDT::Colour3) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
    }
    // decaying colour antisextet
    else if(bColour==PDT::Colour6bar) {
      if (cColour==PDT::Colour3bar && aColour==PDT::Colour3bar) {
	dipoles.push_back(DipoleType(IFba,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFa,ShowerInteraction::QCD));
	dipoles.push_back(DipoleType(IFc,ShowerInteraction::QCD));
      }
    }
  }
  // QED dipoles
  if(inter==ShowerInteraction::ALL || inter==ShowerInteraction::QEDQCD ||
     inter==ShowerInteraction::QED) {
    const bool & bCharged = bProgenitor->dataPtr()->charged();
    const bool & cCharged = cProgenitor->dataPtr()->charged();
    const bool & aCharged = aProgenitor->dataPtr()->charged();
    // initial-final
    if(bCharged && aCharged) {
      dipoles.push_back(DipoleType(IFba,ShowerInteraction::QED));
      dipoles.push_back(DipoleType(IFa ,ShowerInteraction::QED));
    }
    if(bCharged && cCharged) {
      dipoles.push_back(DipoleType(IFbc,ShowerInteraction::QED));
      dipoles.push_back(DipoleType(IFc ,ShowerInteraction::QED));
    }
    // final-state
    if(aCharged && cCharged) {
      dipoles.push_back(DipoleType(FFa,ShowerInteraction::QED));
      dipoles.push_back(DipoleType(FFc,ShowerInteraction::QED));
    }
  }
  // check colour structure is allowed
  return !dipoles.empty();
}

vector<Lorentz5Momentum>  PerturbativeDecayer::hardMomenta(PPtr in, PPtr emitter, 
							   PPtr spectator, 
							   const vector<DipoleType>  &dipoles, 
							   int i, bool inDeadZone) {
  // get masses of the particles
  mb_  = in       ->momentum().mass();
  e_   = emitter  ->momentum().mass()/mb_;
  s_   = spectator->momentum().mass()/mb_;
  e2_  = sqr(e_);
  s2_  = sqr(s_);

  vector<Lorentz5Momentum> particleMomenta;
  Energy2 lambda = sqr(mb_)*sqrt(1.+sqr(s2_)+sqr(e2_)-2.*s2_-2.*e2_-2.*s2_*e2_);    

  // calculate A
  double pre = C_;
  // multiply by the colour factor of the dipole
  // ISR
  if (dipoles[i].type==IFba || dipoles[i].type==IFbc) {
    pre *= colourCoeff(in->dataPtr(),emitter->dataPtr(),spectator->dataPtr(),dipoles[i]);
  }
  // radiation from a/c with initial-final connection
  else if (dipoles[i].type==IFa || dipoles[i].type==IFc) {
    pre *= colourCoeff(emitter->dataPtr(),in->dataPtr(),spectator->dataPtr(),dipoles[i]);
  }
  // radiation from a/c with final-final connection
  else if (dipoles[i].type==FFa || dipoles[i].type==FFc) {
    pre *= colourCoeff(emitter->dataPtr(),spectator->dataPtr(),in->dataPtr(),dipoles[i]);
  }
  double A = 2.*abs(pre)/Constants::twopi;
  // factor due sampling choice
  if(phaseOpt_==0) A *= ymax_;
  // coupling factor
  if(dipoles[i].interaction==ShowerInteraction::QCD)
    A *= alphaS() ->overestimateValue();
  else
    A *= alphaEM()->overestimateValue();

  Energy pTmax = 0.5*mb_*(1.-sqr(s_+e_));

  // if no possible branching return
  if ( pTmax < pTmin_ ) return particleMomenta;
  // loop over the two regions
  for(unsigned int j=0;j<2;++j) {
    Energy pT=pTmax;
    vector<Lorentz5Momentum> momenta(4);
    while (pT >= pTmin_) {
      double ymax;
      // overestimate with flat y limit
      if(phaseOpt_==0) {
	pT *= pow(UseRandom::rnd(),(1./A));
	ymax=ymax_;
      }
      // pT sampling including tighter pT dependent y limit
      else {
	pT = 2.*pTmax*exp(-sqrt(-2.*log(UseRandom::rnd())/A+sqr(log(2.*pTmax/pT))));
	// choice of limit overestimate ln(2*pTmax/pT) (true limit acosh(pTmax/pT))
	ymax = log(2.*pTmax/pT);
      }
      if (pT < pTmin_) break;
      double phi = UseRandom::rnd()*Constants::twopi;
      double y   = ymax*(2.*UseRandom::rnd()-1.);
      double xs, xe, xe_z, xg;
      // check if the momenta are physical
      if (!calcMomenta(j, pT, y, phi, xg, xs, xe,
		       xe_z, momenta)) 
	continue;
      // check if point lies within phase space
      if (!psCheck(xg, xs)) continue;
      // check if point lies within the dead-zone (if required)
      if(inDeadZone && !inTotalDeadZone(xg,xs,dipoles,i)) continue;
      // decay products for 3 body decay
      PPtr inpart   = in        ->dataPtr()->produceParticle(momenta[0]);
      ParticleVector decay3;
      decay3.push_back(emitter  ->dataPtr()->produceParticle(momenta[1]));
      decay3.push_back(spectator->dataPtr()->produceParticle(momenta[2]));
      if(dipoles[i].interaction==ShowerInteraction::QCD)
	decay3.push_back(getParticleData(ParticleID::g    )->produceParticle(momenta[3]));
      else
	decay3.push_back(getParticleData(ParticleID::gamma)->produceParticle(momenta[3]));
      // decay products for 2 body decay
      Lorentz5Momentum p1(ZERO,ZERO, lambda/2./mb_,(mb_/2.)*(1.+e2_-s2_),mb_*e_);
      Lorentz5Momentum p2(ZERO,ZERO,-lambda/2./mb_,(mb_/2.)*(1.+s2_-e2_),mb_*s_);
      ParticleVector decay2;
      decay2.push_back(emitter  ->dataPtr()->produceParticle(p1));
      decay2.push_back(spectator->dataPtr()->produceParticle(p2));
      if (dipoles[i].type==FFa || dipoles[i].type==IFa || dipoles[i].type==IFba) {
	swap(decay2[0],decay2[1]);
	swap(decay3[0],decay3[1]);
      }
      // calculate matrix element ratio R/B
      double meRatio = matrixElementRatio(*inpart,decay2,decay3,Initialize,dipoles[i].interaction);
      // calculate dipole factor
      double dipoleSum(0.),numerator(0.);
      for (int k=0; k<int(dipoles.size()); ++k) {
	// skip dipoles which are not of the interaction being considered
	if(dipoles[k].interaction!=dipoles[i].interaction) continue;
	pair<double,double> dipole = calculateDipole(dipoles[k],*inpart,decay3);
	dipoleSum += abs(dipole.first);
	if (k==i) numerator = abs(dipole.second);
      }
      meRatio *= numerator/dipoleSum;
      // calculate jacobian
      Energy2 denom = (mb_-momenta[3].e())*momenta[2].vect().mag() -
	momenta[2].e()*momenta[3].z(); 
      InvEnergy2  J  = (momenta[2].vect().mag2())/(lambda*denom);
      // calculate weight
      double weight = enhance_*meRatio*fabs(sqr(pT)*J)/pre/Constants::twopi; 
      if(dipoles[i].interaction==ShowerInteraction::QCD)
	weight *= alphaS() ->ratio(pT*pT);
      else
	weight *= alphaEM()->ratio(pT*pT);
      // accept point if weight > R
      if (pT > pT_ && weight > UseRandom::rnd()) {
	particleMomenta=momenta;
	if (weight > 1.) {
	  generator()->log() << "WEIGHT PROBLEM " << fullName() << " " << weight << "\n";
	  generator()->log() << xe << " " << xs << " " << xg << "\n";
	  for(unsigned int ix=0;ix<particleMomenta.size();++ix)
	    generator()->log() << particleMomenta[ix]/GeV << "\n";
	}
	pT_ = pT;
	break;
      }
    }
  }
  return particleMomenta;
}

bool PerturbativeDecayer::calcMomenta(int j, Energy pT, double y, double phi,
					double& xg, double& xs, double& xe, double& xe_z,
					vector<Lorentz5Momentum>& particleMomenta) {
  // calculate xg
  xg = 2.*pT*cosh(y) / mb_;
  if (xg>(1. - sqr(e_ + s_)) || xg<0.) return false;
  // calculate the two values of zs
  double xT  = 2.*pT / mb_;
  double zg = 2.*pT*sinh(y) / mb_;
  double A = (sqr(xT) - 4. * xg + 4.);
  double B = 2. * zg * (s2_ - e2_ - xg + 1.);
  double det = -4. * (-sqr(s2_) + (2. * e2_ + sqr(xT) - 2. * xg + 2.) * s2_ - sqr(e2_ + xg - 1.)) * sqr(xg - 2.);
  if (det<0.) return false;
  double zs= j==0 ? (-B+sqrt(det))/A : (-B-sqrt(det))/A;
  // zs must be negative
  if(zs>0.) return false;
  xs = sqrt(sqr(zs)+4.*s2_);
  // check value of xs is physical
  if (xs>(1.+s2_-e2_) || xs<2.*s_) return false;
  // calculate xe
  xe = 2.-xs-xg;     
  // check value of xe is physical
  if (xe>(1.+e2_-s2_) || xe<2.*e_) return false;       
  // calculate xe_z
  xe_z = -zg-zs;
  // calculate 4 momenta
  particleMomenta[0].setE   ( mb_);
  particleMomenta[0].setX   ( ZERO);
  particleMomenta[0].setY   ( ZERO);
  particleMomenta[0].setZ   ( ZERO);
  particleMomenta[0].setMass( mb_);

  particleMomenta[1].setE   ( mb_*xe/2.);
  particleMomenta[1].setX   (-pT*cos(phi));
  particleMomenta[1].setY   (-pT*sin(phi));
  particleMomenta[1].setZ   ( mb_*xe_z/2.);
  particleMomenta[1].setMass( mb_*e_);

  particleMomenta[2].setE   ( mb_*xs/2.);
  particleMomenta[2].setX   ( ZERO);
  particleMomenta[2].setY   ( ZERO);
  particleMomenta[2].setZ   ( mb_*zs/2.);
  particleMomenta[2].setMass( mb_*s_);

  particleMomenta[3].setE   ( pT*cosh(y));
  particleMomenta[3].setX   ( pT*cos(phi));
  particleMomenta[3].setY   ( pT*sin(phi));
  particleMomenta[3].setZ   ( pT*sinh(y));
  particleMomenta[3].setMass( ZERO);

  return true;
}

bool PerturbativeDecayer::psCheck(const double xg, const double xs) {

  // check is point is in allowed region of phase space
  double xe_star = (1.-s2_+e2_-xg)/sqrt(1.-xg);
  double xg_star = xg/sqrt(1.-xg);

  if ((sqr(xe_star)-4.*e2_) < 1e-10) return false;
  double xs_max = (4.+4.*s2_-sqr(xe_star+xg_star)+ 
		   sqr(sqrt(sqr(xe_star)-4.*e2_)+xg_star))/ 4.;
  double xs_min = (4.+4.*s2_-sqr(xe_star+xg_star)+ 
		   sqr(sqrt(sqr(xe_star)-4.*e2_)-xg_star))/ 4.;

  if (xs < xs_min || xs > xs_max) return false;

  return true;
}

pair<double,double> PerturbativeDecayer::calculateDipole(const DipoleType & dipoleId,
							 const Particle & inpart,
							 const ParticleVector & decay3) {
  // calculate dipole for decay b->ac
  pair<double,double> dipole = make_pair(0.,0.);
  double x1 = 2.*decay3[0]->momentum().e()/mb_;
  double x2 = 2.*decay3[1]->momentum().e()/mb_;
  double xg = 2.*decay3[2]->momentum().e()/mb_;
  double mu12 = sqr(decay3[0]->mass()/mb_);
  double mu22 = sqr(decay3[1]->mass()/mb_);
  tcPDPtr part[3] = {inpart.dataPtr(),decay3[0]->dataPtr(),decay3[1]->dataPtr()};
  if(dipoleId.type==FFa || dipoleId.type == IFa || dipoleId.type == IFba) {
    swap(part[1],part[2]);
    swap(x1,x2);
    swap(mu12,mu22);
  }
  // radiation from b with initial-final connection 
  if (dipoleId.type==IFba || dipoleId.type==IFbc) {
    dipole.first  = -2./sqr(xg);
    dipole.first *= colourCoeff(part[0],part[1],part[2],dipoleId);
  }
  // radiation from a/c with initial-final connection
  else if (dipoleId.type==IFa || dipoleId.type==IFc) {
    double z  = 1. - xg/(1.-mu22+mu12);
    dipole.first = (-2.*mu12/sqr(1.-x2+mu22-mu12) + (1./(1.-x2+mu22-mu12))*
	      (2./(1.-z)-dipoleSpinFactor(part[1],z))); 
    dipole.first *= colourCoeff(part[1],part[0],part[2],dipoleId);
  }
  // radiation from a/c with final-final connection
  else if (dipoleId.type==FFa || dipoleId.type==FFc) {
    double z = 1. + ((x1-1.+mu22-mu12)/(x2-2.*mu22));
    double y = (1.-x2-mu12+mu22)/(1.-mu12-mu22);
    double vt = sqrt((1.-sqr(e_+s_))*(1.-sqr(e_-s_)))/(1.-mu12-mu22);
    double v  = sqrt(sqr(2.*mu22+(1.-mu12-mu22)*(1.-y))-4.*mu22)
      /(1.-y)/(1.-mu12-mu22);
    if(part[1]->iSpin()!=PDT::Spin1) {
      dipole.first = (1./(1.-x2+mu22-mu12))*
	((2./(1.-z*(1.-y)))-vt/v*(dipoleSpinFactor(part[1],z)+(2.*mu12/(1.+mu22-mu12-x2))));
    }
    else {
      dipole.first  = (1./(1.-x2+mu22-mu12))*
	(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))+(z*(1.-z)-2.)/v-vt/v*(2.*mu12/(1.+mu22-mu12-x2)));
      dipole.second = (1./(1.-x2+mu22-mu12))*
	(2./(1.-z*(1.-y))+(z*(1.-z)-2.)/v-vt/v*(2.*mu12/(1.+mu22-mu12-x2)));
    dipole.second   *= colourCoeff(part[1],part[2],part[0],dipoleId);
    }
    dipole.first *= colourCoeff(part[1],part[2],part[0],dipoleId);
  }
  // special for the case that all particles are gluons
  else if(dipoleId.type==FFg) {
    double z = (1.-x2)/xg;
    double y = 1.-xg;
    dipole.first = 1./(1.-xg)*(1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))+(z*(1.-z)-2.));
    dipole.first *= colourCoeff(part[1],part[2],part[0],dipoleId);
  }
  else
    assert(false);
  // coupling prefactors
  if(dipole.second==0.) dipole.second=dipole.first;
  dipole.first  *= 8.*Constants::pi;
  dipole.second *= 8.*Constants::pi;
  // return the answer
  return dipole;
}

double PerturbativeDecayer::dipoleSpinFactor(tcPDPtr part, double z){
  // calculate the spin dependent component of the dipole  
  if      (part->iSpin()==PDT::Spin0)
    return 2.;
  else if (part->iSpin()==PDT::Spin1Half)
    return (1. + z);
  else if (part->iSpin()==PDT::Spin1)
    return -(z*(1.-z) - 1./(1.-z) + 1./z -2.);
  return 0.;
}

namespace {

double colourCharge(PDT::Colour icol) {
  switch(icol) {
  case PDT::Colour0 :
    return 0.;
  case PDT::Colour3 : case PDT::Colour3bar :
    return 4./3.;
  case PDT::Colour8:
    return 3.;
  case PDT::Colour6 : case PDT::Colour6bar :
    return 10./3.;
  default :
    assert(false);
    return 0.;
  }
}
}

double PerturbativeDecayer::colourCoeff(tcPDPtr emitter,
					tcPDPtr spectator,
					tcPDPtr other,
					DipoleType dipole) {
  if(dipole.interaction==ShowerInteraction::QCD) {
    double emitterColour   = colourCharge(emitter  ->iColour());
    double spectatorColour = colourCharge(spectator->iColour());
    double otherColour     = colourCharge(other    ->iColour());
    double val = 0.5*(sqr(emitterColour)+sqr(spectatorColour)-sqr(otherColour))/emitterColour;
    return val;
  }
  else {
    double val = double(emitter->iCharge()*spectator->iCharge())/9.;
    // FF dipoles
    if(dipole.type==FFa || dipole.type == FFc) return -val;
    // IF dipoles
    else                                       return  val;
  }
}

void PerturbativeDecayer::getColourLines(RealEmissionProcessPtr real) {
  // extract the particles
  vector<tPPtr> branchingPart;
  branchingPart.push_back(real->incoming()[0]);
  for(unsigned int ix=0;ix<real->outgoing().size();++ix) {
    branchingPart.push_back(real->outgoing()[ix]);
  }
  vector<unsigned int> sing,trip,atrip,oct,sex,asex;
  for (size_t ib=0;ib<branchingPart.size()-1;++ib) {
    if     (branchingPart[ib]->dataPtr()->iColour()==PDT::Colour0   ) sing. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3   ) trip. push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour3bar) atrip.push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour8   ) oct.  push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour6   ) sex.  push_back(ib);
    else if(branchingPart[ib]->dataPtr()->iColour()==PDT::Colour6bar) asex. push_back(ib);
  }
  // decaying colour singlet
  if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour0) {
    // 0 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	branchingPart[atrip[0]]->colourConnect(branchingPart[   3   ]);
	branchingPart[    3   ]->colourConnect(branchingPart[trip[0]]);
      }
      else {
	branchingPart[atrip[0]]->colourConnect(branchingPart[trip[0]]);
      }
    }
    // 0 -> 8 8
    else if (oct.size()==2 ) {
      if(real->interaction()==ShowerInteraction::QCD) {
	bool col = UseRandom::rndbool();
	branchingPart[oct[0]]->colourConnect(branchingPart[   3  ],col);
	branchingPart[   3  ]->colourConnect(branchingPart[oct[1]],col);
	branchingPart[oct[1]]->colourConnect(branchingPart[oct[0]],col);
      }
      else {
	branchingPart[oct[0]]->colourConnect(branchingPart[oct[1]]);
	branchingPart[oct[1]]->colourConnect(branchingPart[oct[0]]);
      }
    }
    else 
      assert(real->interaction()==ShowerInteraction::QED);
  }
  // decaying colour triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3 ) {
    // 3 -> 3 0
    if (trip.size()==2 && sing.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	branchingPart[3]->incomingColour(branchingPart[trip[0]]);
	branchingPart[3]-> colourConnect(branchingPart[trip[1]]);
      }
      else {
	branchingPart[trip[1]]->incomingColour(branchingPart[trip[0]]);
      }
    }
    // 3 -> 3 8
    else if (trip.size()==2 && oct.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	// 8 emit incoming partner
	if(real->emitter()==oct[0]&&real->spectator()==0) {
	  branchingPart[  3   ]->incomingColour(branchingPart[trip[0]]);
	  branchingPart[  3   ]-> colourConnect(branchingPart[oct[0] ]);
	  branchingPart[oct[0]]-> colourConnect(branchingPart[trip[1]]);
	}
	// 8 emit final spectator or vice veras
	else {
	  branchingPart[oct[0]]->incomingColour(branchingPart[trip[0]]);
	  branchingPart[oct[0]]-> colourConnect(branchingPart[   3   ]);
	  branchingPart[   3  ]-> colourConnect(branchingPart[trip[1]]);
	}
      }
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[trip[0]]);
	branchingPart[oct[0]]-> colourConnect(branchingPart[trip[1]]);
      }
    }
    // 3  -> 3bar 3bar
    else if(trip.size() ==1 && atrip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	if(real->emitter()==atrip[0]) {
	  branchingPart[3]->colourConnect(branchingPart[atrip[0]],true);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[ trip[0]],false),
			       ColourLine::create(branchingPart[       3],true ),
			       ColourLine::create(branchingPart[atrip[1]],true)};
	  col[0]->setSinkNeighbours(col[1],col[2]);
	}
	else {
	  branchingPart[3]->colourConnect(branchingPart[atrip[1]],true);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[ trip[0]],false),
			       ColourLine::create(branchingPart[atrip[0]],true ),
			       ColourLine::create(branchingPart[       3],true)};
	  col[0]->setSinkNeighbours(col[1],col[2]);
	}
      }
      else {
	tColinePtr col[3] = {ColourLine::create(branchingPart[ trip[0]],false),
			     ColourLine::create(branchingPart[atrip[0]],true ),
			     ColourLine::create(branchingPart[atrip[1]],true)};
	col[0]->setSinkNeighbours(col[1],col[2]);
      }
    }
    else
      assert(false);
  }
  // decaying colour anti-triplet
  else if (branchingPart[0]->dataPtr()->iColour()==PDT::Colour3bar) {
    // 3bar -> 3bar 0
    if (atrip.size()==2 && sing.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	branchingPart[3]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[3]-> colourConnect(branchingPart[atrip[1]],true);
      }
      else {
	branchingPart[atrip[1]]->incomingColour(branchingPart[atrip[0]],true);
      }
    }
    // 3 -> 3 8
    else if (atrip.size()==2 && oct.size()==1){
      if(real->interaction()==ShowerInteraction::QCD) {
	// 8 emit incoming partner
	if(real->emitter()==oct[0]&&real->spectator()==0) {
	  branchingPart[   3  ]->incomingColour(branchingPart[atrip[0]],true);
	  branchingPart[   3  ]-> colourConnect(branchingPart[oct[0]  ],true);
	  branchingPart[oct[0]]-> colourConnect(branchingPart[atrip[1]],true);
	}
	// 8 emit final spectator or vice veras
	else {
	  if(real->interaction()==ShowerInteraction::QCD) {
	    branchingPart[oct[0]]->incomingColour(branchingPart[atrip[0]],true);
	    branchingPart[oct[0]]-> colourConnect(branchingPart[   3    ],true);
	    branchingPart[3]-> colourConnect(branchingPart[atrip[1]]     ,true);
	  }
	}
      }
      else {
	branchingPart[oct[0]]->incomingColour(branchingPart[atrip[0]],true);
	branchingPart[oct[0]]-> colourConnect(branchingPart[atrip[1]],true);
      }
    }
    // 3bar  -> 3 3 
    else if(atrip.size() ==1 && trip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	if(real->emitter()==trip[0]) {
	  branchingPart[3]->colourConnect(branchingPart[trip[0]],false);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[atrip[0]],true ),
			       ColourLine::create(branchingPart[       3],false),
			       ColourLine::create(branchingPart[ trip[1]],false)};
	  col[0]->setSourceNeighbours(col[1],col[2]);
	}
	else {
	  branchingPart[3]->colourConnect(branchingPart[trip[1]],false);
	  tColinePtr col[3] = {ColourLine::create(branchingPart[atrip[0]],true ),
			       ColourLine::create(branchingPart[ trip[0]],false),
			       ColourLine::create(branchingPart[       3],false)};
	  col[0]->setSourceNeighbours(col[1],col[2]);
	}
      }
      else {
	tColinePtr col[3] = {ColourLine::create(branchingPart[atrip[0]],true ),
			     ColourLine::create(branchingPart[ trip[0]],false),
			     ColourLine::create(branchingPart[ trip[1]],false)};
	col[0]->setSourceNeighbours(col[1],col[2]);
      }
    }
    else
      assert(false);
  }
  // decaying colour octet
  else if(branchingPart[0]->dataPtr()->iColour()==PDT::Colour8 ) {
    // 8 -> 3 3bar
    if (trip.size()==1 && atrip.size()==1) {
      if(real->interaction()==ShowerInteraction::QCD) {
	// 3 emits
	if(trip[0]==real->emitter()) {
	  branchingPart[3]       ->incomingColour(branchingPart[oct[0]] );
	  branchingPart[3]       -> colourConnect(branchingPart[trip[0]]);
	  branchingPart[atrip[0]]->incomingColour(branchingPart[oct[0]],true);
	}
	// 3bar emits
	else {
	  branchingPart[3]       ->incomingColour(branchingPart[oct[0]]  ,true);
	  branchingPart[3]       -> colourConnect(branchingPart[atrip[0]],true);
	  branchingPart[trip[0]]->incomingColour(branchingPart[oct[0]]  );
	}
      }
      else {
	branchingPart[trip[0]]->incomingColour(branchingPart[oct[0]] );
	branchingPart[atrip[0]]->incomingColour(branchingPart[oct[0]],true);
      }
    }
    // 8 -> 8 0 
    else if (sing.size()==1 && oct.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	bool col = UseRandom::rndbool();
	branchingPart[   3  ]->colourConnect (branchingPart[oct[1]], col);
	branchingPart[   3  ]->incomingColour(branchingPart[oct[0]], col);
	branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]],!col);
      }
      else {
	branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]]);
	branchingPart[oct[1]]->incomingColour(branchingPart[oct[0]],true);
      }
    }
    else
      assert(false);
  }
  // sextet
  else if(branchingPart[0]->dataPtr()->iColour() == PDT::Colour6) {
    if(trip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	if(trip[0]==real->emitter()) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[3]);
	  branchingPart[3]       -> colourConnect(branchingPart[trip[0]]);
	  cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[trip[1]]);
	}
	else {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[3]);
	  branchingPart[3]       -> colourConnect(branchingPart[trip[1]]);
	  cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[trip[0]]);
	}
      }
      else {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	for(unsigned int ix=0;ix<2;++ix) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->colourLine(cline);
	  cline->addColoured(branchingPart[trip[ix]]);
	}
      }
    }
    else
      assert(false);
  }
  // antisextet
  else if(branchingPart[0]->dataPtr()->iColour() == PDT::Colour6bar) {
    if(atrip.size()==2) {
      if(real->interaction()==ShowerInteraction::QCD) {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	if(atrip[0]==real->emitter()) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[3]);
	  branchingPart[3]->antiColourConnect(branchingPart[atrip[0]]);
	  cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[atrip[1]]);
	}
	else {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[3]);
	  branchingPart[3]->antiColourConnect(branchingPart[atrip[1]]);
	  cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addAntiColoured(branchingPart[trip[0]]);
	}
      }
      else {
	Ptr<MultiColour>::pointer parentColour = 
	  dynamic_ptr_cast<Ptr<MultiColour>::pointer>
	  (branchingPart[0]->colourInfo());
	for(unsigned int ix=0;ix<2;++ix) {
	  ColinePtr cline = new_ptr(ColourLine());
	  parentColour->antiColourLine(cline);
	  cline->addColoured(branchingPart[atrip[ix]],true);
	}
      }
    }
    else
      assert(false);
  }
  else
    assert(false);
}

PerturbativeDecayer::phaseSpaceRegion
PerturbativeDecayer::inInitialFinalDeadZone(double xg, double xa,
					    double a, double c) const {
  double lam    = sqrt(1.+a*a+c*c-2.*a-2.*c-2.*a*c);
  double kappab = 1.+0.5*(1.-a+c+lam);
  double kappac = kappab-1.+c;
  double kappa(0.);
  // check whether or not in the region for emission from c
  double r = 0.5;
  if(c!=0.) r += 0.5*c/(1.+a-xa);
  double pa = sqrt(sqr(xa)-4.*a);
  double z = ((2.-xa)*(1.-r)+r*pa-xg)/pa;
  if(z<1. && z>0.) {
    kappa = (1.+a-c-xa)/(z*(1.-z));
    if(kappa<kappac)
      return emissionFromC;
  }
  // check in region for emission from b (T1)
  double cq = sqr(1.+a-c)-4*a;
  double bq = -2.*kappab*(1.-a-c);
  double aq = sqr(kappab)-4.*a*(kappab-1);
  double dis = sqr(bq)-4.*aq*cq;
  z=1.-(-bq-sqrt(dis))/2./aq;
  double w = 1.-(1.-z)*(kappab-1.);
  double xgmax = (1.-z)*kappab;
  // possibly in T1 region
  if(xg<xgmax) {
    z = 1.-xg/kappab;
    kappa=kappab;
  }
  // possibly in T2 region
  else {
    aq = 4.*a;
    bq = -4.*a*(2.-xg);
    cq = sqr(1.+a-c-xg);
    dis = sqr(bq)-4.*aq*cq;
    z = (-bq-sqrt(dis))/2./aq;
    kappa = xg/(1.-z);
  }
  // compute limit on xa
  double u = 1.+a-c-(1.-z)*kappa;
  w = 1.-(1.-z)*(kappa-1.);
  double v = sqr(u)-4.*z*a*w;
  if(v<0. && v>-1e-10) v= 0.;
  v = sqrt(v);
  if(xa<0.5*((u+v)/w+(u-v)/z)) {
    if(xg<xgmax)
      return emissionFromA1;
    else if(useMEforT2_)
      return deadZone;
    else
      return emissionFromA2;
  }
  else
    return deadZone;
}

PerturbativeDecayer::phaseSpaceRegion
PerturbativeDecayer::inFinalFinalDeadZone(double xb, double xc,
					  double b, double c) const {
  // basic kinematics
  double lam = sqrt(1.+b*b+c*c-2.*b-2.*c-2.*b*c);
  // check whether or not in the region for emission from b
  double r = 0.5;
  if(b!=0.) r+=0.5*b/(1.+c-xc);
  double pc = sqrt(sqr(xc)-4.*c);
  double z = -((2.-xc)*r-r*pc-xb)/pc;
  if(z<1. and z>0.) {
    if((1.-b+c-xc)/(z*(1.-z))<0.5*(1.+b-c+lam)) return emissionFromB;
  }
  // check whether or not in the region for emission from c
  r = 0.5;
  if(c!=0.) r+=0.5*c/(1.+b-xb);
  double pb = sqrt(sqr(xb)-4.*b);
  z = -((2.-xb)*r-r*pb-xc)/pb;
  if(z<1. and z>0.) {
    if((1.-c+b-xb)/(z*(1.-z))<0.5*(1.-b+c+lam)) return emissionFromC;
  }
  return deadZone;
}

bool PerturbativeDecayer::inTotalDeadZone(double xg, double xs,
					  const vector<DipoleType>  & dipoles,
					  int i) {
  double xb,xc,b,c;
  if(dipoles[i].type==FFa || dipoles[i].type == IFa || dipoles[i].type == IFba) {
    xc = xs;
    xb = 2.-xg-xs;
    b = e2_;
    c = s2_;
  }
  else {
    xb = xs;
    xc = 2.-xg-xs;
    b = s2_;
    c = e2_;
  }
  for(unsigned int ix=0;ix<dipoles.size();++ix) {
    if(dipoles[ix].interaction!=dipoles[i].interaction)
      continue;
    // should also remove negative QED dipoles but shouldn't be an issue unless we
    // support QED ME corrections
    switch (dipoles[ix].type) {
    case FFa :
      if(inFinalFinalDeadZone(xb,xc,b,c)!=deadZone) return false;
      break;
    case FFc :
      if(inFinalFinalDeadZone(xc,xb,c,b)!=deadZone) return false;
      break;
    case IFa : case IFba:
      if(inInitialFinalDeadZone(xg,xc,c,b)!=deadZone) return false;
      break;
    case IFc : case IFbc:
      if(inInitialFinalDeadZone(xg,xb,b,c)!=deadZone) return false;
      break;
    case FFg:
      break;
    }
  }
  return true;
}
