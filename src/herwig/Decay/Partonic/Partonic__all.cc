#line 1 "./DarkQuarkoniumDecayer.cc"
// -*- C++ -*-
//
// DarkQuarkoniumDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "DarkQuarkoniumDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>
#include <cassert>

using namespace Herwig;

DarkQuarkoniumDecayer::DarkQuarkoniumDecayer() : MECode(0) {} 

IBPtr DarkQuarkoniumDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr DarkQuarkoniumDecayer::fullclone() const {
  return new_ptr(*this);
}

void DarkQuarkoniumDecayer::Init() {
  
  static ClassDocumentation<DarkQuarkoniumDecayer> documentation
    ("The DarkQuarkoniumDecayer performs partonic decays of quarkonium"
     " resonances");
  
  static Switch<DarkQuarkoniumDecayer,int> interfaceMECode
    ("MECode",
     "The code for the ME type to use in the decay",
     &DarkQuarkoniumDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Use a phase-space distribution",
     0);
  static SwitchOption interfaceMECodeOrePowell
    (interfaceMECode,
     "OrePowell",
     "Use the Ore-Powell matrix element",
     130);
  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DarkQuarkoniumDecayer,PartonicDecayerBase>
describeHerwigDarkQuarkoniumDecayer("Herwig::DarkQuarkoniumDecayer", "HwShower.so HwPartonicDecay.so");

bool DarkQuarkoniumDecayer::accept(tcPDPtr, const tPDVector & children) const {
  return (children.size() == 3 || children.size() == 2);
}

ParticleVector DarkQuarkoniumDecayer::decay(const Particle & p,
					const tPDVector & children) const {
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
  }
  assert(partons.size()==2 || partons.size()==3);
  Lorentz5Momentum products[3];
  Energy gluMass = getParticleData(ParticleID::g)->constituentMass();
  for(unsigned int i = 0; i<partons.size(); i++) {
    if(partons[i]->id() == ParticleID::g) products[i].setMass(gluMass);
    else                                  products[i].setMass(partons[i]->mass());
  }
  if(partons.size() == 3) {
    // 2 gluon or quark + dark pion decay
    // Ore & Powell orthopositronium matrix element
    if(MECode == 130) { 
      double x1, x2, x3, test;
      do {
	// if decay fails return empty vector.
	if (! Kinematics::threeBodyDecay(p.momentum(), products[0],
					 products[1],  products[2]))
	  return ParticleVector();
	x1 = 2.*(p.momentum()*products[0])/sqr(p.mass());
	x2 = 2.*(products[1]*products[2])/sqr(p.mass());
	x3 = 2. - x1 - x2;
	test = sqr(x1*(1.-x1)) + sqr(x2*(1.-x2)) + sqr(x3*(1.-x3));
	test /= sqr(x1*x2*x3);
      } 
      while(test < 2.*UseRandom::rnd());
    }
    else {
      if (! Kinematics::threeBodyDecay(p.momentum(), products[0], 
				       products[1],products[2]))
	return ParticleVector();
    }
    // test the momenta
    for(unsigned int i = 0; i<partons.size(); i++)
      partons[i]->set5Momentum(products[i]);
    // Now set colour connections
    partons[0]->colourNeighbour(partons[1]);
    if(abs(partons[0]->id()) == ParticleID::g) 
      partons[0]->colourNeighbour(partons[1]);
  } 
  // two decay children
  // 2 gluon or q-qbar decay
  else {
    double Theta, Phi;
    Kinematics::generateAngles(Theta, Phi);
    Energy p1 = partons[0]->mass();
    Energy p2 = partons[1]->mass();
    if(p1 == ZERO) p1 = gluMass;
    if(p2 == ZERO) p2 = gluMass;
    if (! Kinematics::twoBodyDecay(p.momentum(), p1, p2, Theta, Phi,
				   products[0], products[1]))
      return ParticleVector();
    for(unsigned int i = 0; i<partons.size(); i++)
      partons[i]->set5Momentum(products[i]);
    int first(0), second(1); 
    if(partons[0]->id() < 0) swap(first,second);
    partons[first]->antiColourNeighbour(partons[second]);
    if(abs(partons[first]->id()) == ParticleID::g) 
      partons[first]->colourNeighbour(partons[second]);
  }
  return partons;
}
   
void DarkQuarkoniumDecayer::persistentOutput(PersistentOStream &os) const { 
  os << MECode;
}

void DarkQuarkoniumDecayer::persistentInput(PersistentIStream &is, int) { 
  is >> MECode;
}

void DarkQuarkoniumDecayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./QuarkoniumDecayer.cc"
// -*- C++ -*-
//
// QuarkoniumDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "QuarkoniumDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>
#include <cassert>

using namespace Herwig;

QuarkoniumDecayer::QuarkoniumDecayer() : MECode(0) {} 

IBPtr QuarkoniumDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr QuarkoniumDecayer::fullclone() const {
  return new_ptr(*this);
}

void QuarkoniumDecayer::Init() {
  
  static ClassDocumentation<QuarkoniumDecayer> documentation
    ("The QuarkoniumDecayer performs partonic decays of quarkonium"
     " resonances");
  
  static Switch<QuarkoniumDecayer,int> interfaceMECode
    ("MECode",
     "The code for the ME type to use in the decay",
     &QuarkoniumDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Use a phase-space distribution",
     0);
  static SwitchOption interfaceMECodeOrePowell
    (interfaceMECode,
     "OrePowell",
     "Use the Ore-Powell matrix element",
     130);
  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QuarkoniumDecayer,PartonicDecayerBase>
describeHerwigQuarkoniumDecayer("Herwig::QuarkoniumDecayer", "HwShower.so HwPartonicDecay.so");

bool QuarkoniumDecayer::accept(tcPDPtr, const tPDVector & children) const {
  return (children.size() == 3 || children.size() == 2);
}

ParticleVector QuarkoniumDecayer::decay(const Particle & p,
					const tPDVector & children) const {
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
  }
  assert(partons.size()==2 || partons.size()==3);
  Lorentz5Momentum products[3];
  Energy gluMass = getParticleData(ParticleID::g)->constituentMass();
  for(unsigned int i = 0; i<partons.size(); i++) {
    if(partons[i]->id() == ParticleID::g) products[i].setMass(gluMass);
    else                                  products[i].setMass(partons[i]->mass());
  }
  if(partons.size() == 3) {
    // 3-gluon or 2-gluon + photon decay
    // Ore & Powell orthopositronium matrix element
    if(MECode == 130) { 
      double x1, x2, x3, test;
      do {
	// if decay fails return empty vector.
	if (! Kinematics::threeBodyDecay(p.momentum(), products[0],
					 products[1],  products[2]))
	  return ParticleVector();
	x1 = 2.*(p.momentum()*products[0])/sqr(p.mass());
	x2 = 2.*(products[1]*products[2])/sqr(p.mass());
	x3 = 2. - x1 - x2;
	test = sqr(x1*(1.-x1)) + sqr(x2*(1.-x2)) + sqr(x3*(1.-x3));
	test /= sqr(x1*x2*x3);
      } 
      while(test < 2.*UseRandom::rnd());
    }
    else {
      if (! Kinematics::threeBodyDecay(p.momentum(), products[0], 
				       products[1],products[2]))
	return ParticleVector();
    }
    // test the momenta
    for(unsigned int i = 0; i<partons.size(); i++)
      partons[i]->set5Momentum(products[i]);
    // Now set colour connections
    if(partons[2]->id() == ParticleID::g) {
      partons[0]->colourNeighbour(partons[1]);
      partons[1]->colourNeighbour(partons[2]);
      partons[2]->colourNeighbour(partons[0]);
    } 
    else {
      partons[0]->colourNeighbour(partons[1]);
      partons[0]->antiColourNeighbour(partons[1]);
    }
  } 
  // two decay children
  // 2 gluon or q-qbar decay
  else {
    double Theta, Phi;
    Kinematics::generateAngles(Theta, Phi);
    Energy p1 = partons[0]->mass();
    Energy p2 = partons[1]->mass();
    if(p1 == ZERO) p1 = gluMass;
    if(p2 == ZERO) p2 = gluMass;
    if (! Kinematics::twoBodyDecay(p.momentum(), p1, p2, Theta, Phi,
				   products[0], products[1]))
      return ParticleVector();
    for(unsigned int i = 0; i<partons.size(); i++)
      partons[i]->set5Momentum(products[i]);
    int first(0), second(1); 
    if(partons[0]->id() < 0) swap(first,second);
    partons[first]->antiColourNeighbour(partons[second]);
    if(abs(partons[first]->id()) == ParticleID::g) 
      partons[first]->colourNeighbour(partons[second]);
  }
  return partons;
}
   
void QuarkoniumDecayer::persistentOutput(PersistentOStream &os) const { 
  os << MECode;
}

void QuarkoniumDecayer::persistentInput(PersistentIStream &is, int) { 
  is >> MECode;
}

void QuarkoniumDecayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./HeavyDecayer.cc"
// -*- C++ -*-
//
// HeavyDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HeavyDecayer class.
//

#include "HeavyDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include <ThePEG/PDT/EnumParticles.h>
#include <ThePEG/PDT/DecayMode.h>
#include <ThePEG/Interface/ClassDocumentation.h>
#include <ThePEG/Interface/Switch.h>
#include "Herwig/Utilities/Kinematics.h"
#include <ThePEG/Interface/Parameter.h>
#include <ThePEG/Persistency/PersistentOStream.h>
#include <ThePEG/Persistency/PersistentIStream.h>
#include <ThePEG/Repository/UseRandom.h>
#include <ThePEG/Utilities/Maths.h>

using namespace Herwig;

HeavyDecayer::HeavyDecayer() : MECode(0) {} 

IBPtr HeavyDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr HeavyDecayer::fullclone() const {
  return new_ptr(*this);
}

void HeavyDecayer::Init() {

  static ClassDocumentation<HeavyDecayer> documentation
    ("Class to decay all particles in HERWIG by the algorithms used in HERWIG 6.4");
  
  static Switch<HeavyDecayer,int> interfaceMECode
    ("MECode",
     "The code for the ME type to use in the decay",
     &HeavyDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Use a phase-space distribution",
     0);
  static SwitchOption interfaceMECodeWeak
    (interfaceMECode,
     "Weak",
     "Use the weak V-A matrix element",
     100);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HeavyDecayer,PartonicDecayerBase>
describeHerwigHeavyDecayer("Herwig::HeavyDecayer", "HwShower.so HwPartonicDecay.so");

bool HeavyDecayer::accept(tcPDPtr parent, const tPDVector & children) const { 
  long id = parent->id();
  int flav1, flav2;
  if((id / 1000)%10) {
    flav1 = (id/1000)%10;
    flav2 = (id/10)%100;
  } 
  else {
    flav1 = id/100;  
    flav2 = (id/10)%10;
  }
  if(!flav1 || !flav2) return false;
  return children.size()==4;
}

ParticleVector HeavyDecayer::decay(const Particle & p,
				   const tPDVector & children) const {
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
  }
  Lorentz5Momentum products[4];
  for(int i = 0; i<4; i++) products[i].setMass(partons[i]->mass());
  // Fraction of momentum to spectator
  double xs = partons[3]->mass()/p.momentum().e();
  // Fraction of momentum to heavy quark (which is weakly decaying)
  double xb = 1.-xs;
  partons[3]->setMomentum(p.momentum()*xs);
  // Get the particle that is decaying
  long idSpec = partons[3]->id();
  long idQ = (p.id()/1000)%10;
  if(!idQ) idQ = (p.id()/100)%10;
  // Now the odd case of a B_c where the c decays, not the b
  if(idSpec == idQ) idQ = (p.id()/10)%10;
  PPtr inter = getParticleData(idQ)->produceParticle(p.momentum()*xb);
  // Three Body Decay
  // Free Massless (V-A)*(V-A) ME
  if(MECode == 100) {
    Energy Mw = getParticleData(ParticleID::Wplus)->mass();
    Energy GamW = getParticleData(ParticleID::Wplus)->width();
    Energy2 EMwSq = sqr(Mw);
    Energy4 GMwSq = sqr(Mw*GamW);
    Energy4 EmLim = GMwSq + sqr(EMwSq - sqr(inter->mass()-p.mass()));
    Energy4 EmTest;
    do {
      Kinematics::threeBodyDecay(inter->momentum(),products[0], products[1], 
				 products[2], &VAWt);
      Energy2 pw2 = (products[0]+products[1]).m2();
      EmTest = sqr(pw2-EMwSq);
    } 
    while((EmTest+GMwSq)*rnd() > EmLim);
  }
  else Kinematics::threeBodyDecay(inter->momentum(),products[0], products[1],
				  products[2]);
  // set the momenta of the products
  for(int i = 0; i<3; i++) partons[i]->setMomentum(products[i]);
  // Set up colour connections based on the diagram above and input order
  if(partons[0]->coloured()) {
    if(partons[0]->id() > 0) partons[0]->antiColourNeighbour(partons[1]);
    else                     partons[0]->colourNeighbour(    partons[1]);
  }
  if(partons[2]->coloured()) {
    if(partons[2]->id() > 0) partons[2]->antiColourNeighbour(partons[3]);
    else                     partons[2]->colourNeighbour(partons[3]);
  }
  return partons;
}
   
void HeavyDecayer::persistentOutput(PersistentOStream &os) const {
  os << MECode; 
}

void HeavyDecayer::persistentInput(PersistentIStream &is, int) {
  is >> MECode; 
}

double HeavyDecayer::VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3) {
  return (t1-t0)*(t0-t2)*t3; 
}

void HeavyDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./WeakPartonicDecayer.cc"
// -*- C++ -*-
//
// WeakPartonicDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakPartonicDecayer class.
//

#include "WeakPartonicDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ConstituentParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

WeakPartonicDecayer::WeakPartonicDecayer() : MECode(0), _radprob(0.0), _maxtry(300),
					     _threemax(3.), _fourmax(3.)
{}

IBPtr WeakPartonicDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr WeakPartonicDecayer::fullclone() const {
  return new_ptr(*this);
}

bool WeakPartonicDecayer::accept(tcPDPtr parent, const tPDVector & prod) const {
  // check we can find the flavours of the quarks in the decaying meson
  long id = parent->id();
  int flav1, flav2;
  if((id / 1000)%10) {
    flav1 = (id/1000)%10;
    flav2 = (id/10)%100;
  }
  else {
    flav1 = id/100;
    flav2 = (id/10)%10;
  }
  if(!flav1 || !flav2) return false;
  // if two decay products one must be in triplet and one antitriplet
  if(prod.size()==2) {
    if((prod[0]->iColour()==PDT::Colour3&&prod[1]->iColour()==PDT::Colour3bar)||
       (prod[0]->iColour()==PDT::Colour3bar&&prod[1]->iColour()==PDT::Colour3))
      return true;
  }
  else if(prod.size()==3) {
    if(((prod[0]->iColour()==PDT::Colour3   &&prod[2]->iColour()==PDT::Colour3bar)||
	(prod[0]->iColour()==PDT::Colour3bar&&prod[2]->iColour()==PDT::Colour3))
       &&prod[1]->iColour()==PDT::Colour8) return true;
  }
  else if(prod.size()==4) {
    // first two particles should be leptons or q qbar
    if((prod[0]->id()>=11&&prod[0]->id()<=16&&prod[1]->id()<=-11&&prod[1]->id()>=-16)||
       (prod[1]->id()>=11&&prod[1]->id()<=16&&prod[0]->id()<=-11&&prod[0]->id()>=-16)||
       (prod[0]->iColour()==PDT::Colour3    &&prod[1]->iColour()==PDT::Colour3bar   )||
       (prod[1]->iColour()==PDT::Colour3    &&prod[0]->iColour()==PDT::Colour3bar   )) {
      // third particle quark and fourth colour anti-triplet or
      // thrid particle antiquark and fourth colour triplet
      if((prod[2]->iColour()==PDT::Colour3bar&&prod[3]->iColour()==PDT::Colour3   )||
	 (prod[2]->iColour()==PDT::Colour3   &&prod[3]->iColour()==PDT::Colour3bar))
	return true;
    }
    else return false;
  }
  return false;
}

ParticleVector WeakPartonicDecayer::decay(const Particle & parent,
					  const tPDVector & children) const {
  static tcPDPtr gluon=getParticleData(ParticleID::g);
  // make the particles
  ParticleVector partons;
  for(unsigned int ix=0;ix<children.size();++ix) {
    partons.push_back(children[ix]->produceParticle());
    // these products have the mass but should have constituent mass
    partons[ix]->set5Momentum(Lorentz5Momentum(children[ix]->constituentMass()));
  }
  // 2-body decays
  if(partons.size()==2) {
    // no gluon if not select based on probability or if three body not allowed
    if(UseRandom::rnd()>_radprob||
       parent.mass()<gluon->constituentMass()+partons[0]->mass()+partons[1]->mass()) {
      double ctheta,phi;
      Lorentz5Momentum pout[2];
      for(unsigned int ix=0;ix<2;++ix) pout[ix].setMass(partons[ix]->mass());
      Kinematics::generateAngles(ctheta,phi);
      Kinematics::twoBodyDecay(parent.momentum(),pout[0].mass(),pout[1].mass(),
			       ctheta,phi,pout[0],pout[1]);
      for(unsigned int ix=0; ix<2;++ix) partons[ix]->setMomentum(pout[ix]);
      if(partons[0]->dataPtr()->iColour()==PDT::Colour3) {
	partons[0]->antiColourNeighbour(partons[1]);
      }
      else {
	partons[0]->    colourNeighbour(partons[1]);
      }
    }
    else {
      Lorentz5Momentum pout[3];
      for(unsigned int ix=0;ix<2;++ix) pout[ix].setMass(partons[ix]->mass());
      // add gluon
      partons.push_back(gluon->produceParticle());
      partons.back()->set5Momentum(gluon->constituentMass());
      // momentum of gluon
      pout[2] = Lorentz5Momentum(gluon->constituentMass());
      Kinematics::threeBodyDecay(parent.momentum(),pout[1],pout[0],pout[2]);
      for(unsigned int ix=0; ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
      if(partons[0]->dataPtr()->iColour()==PDT::Colour3) {
	partons[0]->antiColourNeighbour(partons[2]);
	partons[1]->    colourNeighbour(partons[2]);
      }
      else {
	partons[0]->    colourNeighbour(partons[2]);
	partons[1]->antiColourNeighbour(partons[2]);
      }
    }
  }
  // 3-body decays
  else if(partons.size()==3) {
    // set masses of products
    Lorentz5Momentum pout[3],pin(parent.momentum());
    for(unsigned int ix=0;ix<3;++ix) pout[ix].setMass(partons[ix]->mass());
    double xs(partons[2]->mass()/pin.mass()),xb(1.-xs);
    pout[2]=xs*pin;
    // Get the particle quark that is decaying
    long idQ, idSpec;
    idSpec = partons[2]->id();
    idQ = (parent.id()/1000)%10;
    if(!idQ) idQ = (parent.id()/100)%10;
    // Now the odd case of a B_c where the c decays, not the b
    if(idSpec == idQ) idQ = (parent.id()/10)%10;
    // momentum of the decaying quark
    PPtr inter = getParticleData(idQ)->produceParticle(parent.momentum()*xb);
    // two body decay of heavy quark
    double ctheta,phi;
    Kinematics::generateAngles(ctheta,phi);
    Kinematics::twoBodyDecay(inter->momentum(),pout[0].mass(),pout[1].mass(),
			     ctheta,phi,pout[0],pout[1]);
    // set the momenta of the decay products
    for(unsigned int ix=0; ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
    // make the colour connections
    // quark first
    if(partons[0]->data().iColour()==PDT::Colour3) {
      partons[0]->antiColourNeighbour(partons[1]);
      partons[1]->colourNeighbour(partons[0]);
      partons[1]->antiColourNeighbour(partons[2]);
      partons[2]->colourNeighbour(partons[1]);
    }
    // antiquark first
    else {
      partons[0]->colourNeighbour(partons[1]);
      partons[1]->antiColourNeighbour(partons[0]);
      partons[1]->colourNeighbour(partons[2]);
      partons[2]->antiColourNeighbour(partons[1]);
    }
  }
  // 4-body decays
  else if(partons.size()==4) {
    // swap 0 and 1 if needed
    if(partons[1]->dataPtr()->iColour()!=PDT::Colour0&&
       partons[1]->dataPtr()->iColour()!=partons[2]->dataPtr()->iColour())
      swap(partons[0],partons[1]);
    // get the momenta of the decaying quark and the spectator
    Lorentz5Momentum pin(parent.momentum());
    double xs(partons[3]->mass()/pin.mass()),xb(1.-xs);
    Lorentz5Momentum pspect(xs*pin),pdec(xb*pin);
    pspect.setMass(partons[3]->mass());
    pdec.rescaleMass();
    // Get the particle quark that is decaying
    long idSpec = partons[3]->id();
    long idQ = (abs(parent.id())/1000)%10;
    if(!idQ) idQ = (abs(parent.id())/100)%10;
    // special for doubly heavy baryons, diquark decays
    if (abs(parent.id())%2==0 && idQ>=4 &&
	((abs(parent.id())%1000)/100>=4 || (abs(parent.id())%100)/10>=4)) {
      long id2=max((abs(parent.id())%1000)/100, (abs(parent.id())%100)/10);
      if (id2>idQ)     idQ=id2*1000+idQ*100+1;
      else if(idQ<id2) idQ=id2*100+idQ*1000+1;
      else             idQ=id2*100+idQ*1000+3;
      if(parent.id()<0) idQ=-idQ;
    }
    else {
      // Now the odd case of a B_c where the c decays, not the b
      if(abs(idSpec) == idQ) idQ = (abs(parent.id())/10)%10;
      // change sign if spectator quark or antidiquark
      if((idSpec>0&&idSpec<6)||idSpec<-6) idQ = -idQ;
    }
    // check if W products coloured
    bool Wcol = partons[0]->coloured();
    // particle data object
    tcPDPtr dec = getParticleData(idQ);
    // spin density matrix for the decaying quark
    RhoDMatrix rhoin(PDT::Spin1Half,true);
    if(parent.dataPtr()->iSpin()!=PDT::Spin0 && parent.spinInfo()) {
      parent.spinInfo()->decay();
      RhoDMatrix rhoHadron = parent.spinInfo()->rhoMatrix();
      // particles with spin 0 diquark
      if(abs(parent.id())==5122 || abs(parent.id())==4122 ||
	 abs(parent.id())==5122 || abs(parent.id())==4122 ||
	 abs(parent.id())==5132 || abs(parent.id())==4132 ||
	 abs(parent.id())==5232 || abs(parent.id())==4232) {
	for(unsigned int ix=0;ix<2;++ix) rhoin(ix,ix) = rhoHadron(ix,ix);
      }
      // particles with spin 1 diquark
      else if(abs(parent.id())==5332 || abs(parent.id())==4332) {
	rhoin(0,0) = 2./3.*rhoHadron(1,1)+1./3.*rhoHadron(0,0);
	rhoin(1,1) = 2./3.*rhoHadron(0,0)+1./3.*rhoHadron(1,1);
      }
    }
    // momenta of the decay products
    vector<Lorentz5Momentum> pout(3,Lorentz5Momentum());
    for(unsigned int ix=0;ix<3;++ix) pout[ix].setMass(partons[ix]->mass());
    // charges of the exchanged boson and check if colour rearranged
    int c1 = dec                  ->iCharge()-partons[2]->dataPtr()->iCharge();
    int c2 = partons[0]->dataPtr()->iCharge()+partons[1]->dataPtr()->iCharge();
    bool rearranged = !(c1==c2&&abs(c1)==3);
    if(MECode==0) rearranged=false;
    if(rearranged) {
      int c3 = dec                  ->iCharge()-partons[1]->dataPtr()->iCharge();
      int c4 = partons[0]->dataPtr()->iCharge()+partons[2]->dataPtr()->iCharge();
      if(!(c3==c4&&abs(c3)==3)) {
	generator()->log() << "Unknown order for colour rearranged decay"
			   << " in WeakPartonicDecayer::decay()\n";
	generator()->log() << c1 << " " << c2 << " " << c3 << " " << c4 << "\n";
	generator()->log() << parent << "\n" << dec->PDGName() << "\n";
	for(unsigned int ix=0;ix<4;++ix) generator()->log() << *partons[ix] << "\n";
	throw Exception()  << "Unknown order for colour rearranged decay"
			   << " in WeakPartonicDecayer::decay() "
			   << Exception::runerror;
      }
      swap(pout[1]   ,pout[2]   );
      swap(partons[1],partons[2]);
    }
    // decide if three or four body using prob
    bool threeBody = UseRandom::rnd() > _radprob;
    // if four body not kinematically possible must be three body
    if(pdec.mass()<gluon->constituentMass()+pout[0].mass()+
       pout[1].mass()+pout[2].mass()) threeBody=true;
    // if code ==0 always three body
    if(MECode==0) threeBody=true;
    // three body decay
    if( threeBody ) {
      if(MECode==0) {
				Kinematics::threeBodyDecay(pdec,pout[1],pout[0],pout[2]);
				// set momenta of particles
				for(unsigned int ix=0;ix<pout.size();++ix)
					partons[ix]->setMomentum(pout[ix]);
      }
      else {
	// generate the kinematics
	double wgt(0.);
	Energy2 mb2max = sqr(pdec.mass()    - pout[2].mass());
	Energy2 mb2min = sqr(pout[0].mass() + pout[1].mass());
	unsigned int ntry = 0;
	do {
	  ++ntry;
	  Energy2 mb2 = (mb2max-mb2min)*UseRandom::rnd()+mb2min;
	  double CosAngle, AzmAngle;
	  // perform first decay
	  Lorentz5Momentum p01;
	  p01.setMass(sqrt(mb2));
	  Kinematics::generateAngles(CosAngle,AzmAngle);
	  Kinematics::twoBodyDecay(pdec,pout[2].mass(),p01.mass(),
				   CosAngle,AzmAngle,pout[2],p01);
	  // perform second decay
	  Kinematics::generateAngles(CosAngle,AzmAngle);
	  Kinematics::twoBodyDecay(p01,pout[0].mass(),pout[1].mass(),
				   CosAngle,AzmAngle,pout[0],pout[1]);
	  // kinematic piece of the weight
	  wgt =
	    Kinematics::pstarTwoBodyDecay(pdec.mass(),p01    .mass(),pout[2].mass())/pdec.mass()*
	    Kinematics::pstarTwoBodyDecay(p01 .mass(),pout[0].mass(),pout[1].mass())/p01.mass();
	  // piece to improve weight variation (not kinematics dependent)
	  // and integration over m23^2 (N.B. m23^2 fac on bottom for efficiency)
	  wgt *= pdec.mass()/Kinematics::pstarTwoBodyDecay(pdec.mass(),sqrt(mb2min),pout[2].mass())
	    /(mb2max-mb2min)*sqr(pdec.mass());
	  // set momenta of particles
	  for(unsigned int ix=0;ix<pout.size();++ix) partons[ix]->setMomentum(pout[ix]);
	  // matrix element piece
	  if(dec->iSpin()==2)
	    wgt *= threeBodyMatrixElement(dec,rhoin,pdec,partons);
	  else
	    wgt *= 16.*(pdec*pout[1])*(pout[0]*pout[2])/sqr(sqr(pdec.mass()));
	  // check doesn't violate max
	  if(wgt>_threemax) {
	    ostringstream message;
	    message << "Maximum weight for three-body decay "
		    << "violated in WeakPartonicDecayer::decay()"
		    << "Maximum = " << _threemax << " weight = " << wgt;
	    generator()->logWarning( Exception(message.str(),Exception::warning) );
	  }
	}
	while( wgt < _threemax*UseRandom::rnd() && ntry < _maxtry );
	if(ntry==_maxtry) throw Exception()
	  << "Too many attempts to generate three body kinematics in "
	  << "WeakPartonicDecayer::decay()" << Exception::eventerror;
      }
      partons[3]->setMomentum(pspect);
      // set up the colour connections
      if(rearranged) swap(partons[1],partons[2]);
      if(Wcol) {
	if(partons[0]->data().iColour()==PDT::Colour3)
	  partons[0]->antiColourNeighbour(partons[1]);
	else
	  partons[0]->    colourNeighbour(partons[1]);
      }
      if(partons[2]->data().iColour()==PDT::Colour3) {
	partons[2]->antiColourNeighbour(partons[3]);
      }
      else {
	partons[2]->    colourNeighbour(partons[3]);
      }
    }
    // four body decay
    else {
      // generate the extra gluon
      partons.push_back(gluon->produceParticle());
      partons.back()->set5Momentum(gluon->constituentMass());
      // momentum of gluon
      pout.push_back(Lorentz5Momentum(gluon->constituentMass()));
      // generate the kinematics
      Energy2 ms2min(sqr(pout[0].mass()+pout[1].mass()+pout[2].mass()));
      Energy2 ms2max(sqr(pdec.mass()-pout[3].mass()));
      double wgt(0.);
      unsigned int ntry=0;
      bool initial = true;
      do {
	++ntry;
	Energy2 ms2 = ms2min+UseRandom::rnd()*(ms2max-ms2min);
	Energy ms  = sqrt(ms2);
	// and the W
	Energy2 mb2max = sqr(ms            -pout[2].mass());
	Energy2 mb2min = sqr(pout[0].mass()+pout[1].mass());
	Energy2 mb2 = (mb2max-mb2min)*UseRandom::rnd()+mb2min;
	wgt = (mb2max-mb2min)/(ms2max-mb2min);
	// perform first decay
	Lorentz5Momentum ps;
	double CosAngle,AzmAngle;
	Kinematics::generateAngles(CosAngle,AzmAngle);
	Kinematics::twoBodyDecay(pdec,pout[3].mass(),ms,CosAngle,AzmAngle,pout[3],ps);
	// generate the kinematics
	// perform second decay
	Kinematics::generateAngles(CosAngle,AzmAngle);
	Lorentz5Momentum p01;
	p01.setMass(sqrt(mb2));
	Kinematics::twoBodyDecay(ps,pout[2].mass(),p01.mass(),
				 CosAngle,AzmAngle,pout[2],p01);
	// perform third decay
	Kinematics::generateAngles(CosAngle,AzmAngle);
	Kinematics::twoBodyDecay(p01,pout[0].mass(),pout[1].mass(),
				 CosAngle,AzmAngle,pout[0],pout[1]);
	// kinematic piece of the weight
	wgt *= 16.*
	  Kinematics::pstarTwoBodyDecay(pdec.mass(),pout[3].mass(),ms            )/pdec.mass()*
	  Kinematics::pstarTwoBodyDecay(ms         ,p01    .mass(),pout[2].mass())/ms*
	  Kinematics::pstarTwoBodyDecay(p01 .mass(),pout[0].mass(),pout[1].mass())/p01.mass();
	wgt *= fourBodyMatrixElement(pdec,pout[2],pout[0],pout[1],pout[3],Wcol,initial);
	// check doesn't violate max
	if(wgt>_threemax) {
	  ostringstream message;
	  message << "Maximum weight for four-body decay "
		  << "violated in WeakPartonicDecayer::decay()"
		  << "Maximum = " << _fourmax << " weight = " << wgt;
	  generator()->logWarning( Exception(message.str(),Exception::warning) );
	}
      }
      while ( wgt < _fourmax*UseRandom::rnd() && ntry < _maxtry );
      if(ntry==_maxtry) throw Exception()
	<< "Too many attempts to generate four body kinematics in "
	<< "WeakPartonicDecayer::decay()" << Exception::eventerror;
      // set momenta of particles
      for(unsigned int ix=0;ix<3;++ix) partons[ix]->setMomentum(pout[ix]);
      partons[3]->setMomentum(pspect);
      partons[4]->setMomentum(pout[3]);
      // special for tau leptons to get correlations
      threeBodyMatrixElement(dec,rhoin,pdec,partons);
      // set up the colour connections
      if(rearranged) swap(partons[1],partons[2]);
      // radiation from initial-state
      if(initial) {
	if(Wcol) {
	  if(partons[0]->data().iColour()==PDT::Colour3)
	    partons[0]->antiColourNeighbour(partons[1]);
	  else
	    partons[0]->    colourNeighbour(partons[1]);
	}
	if(partons[2]->data().iColour()==PDT::Colour3) {
	  partons[2]->antiColourNeighbour(partons[4]);
	  partons[3]->    colourNeighbour(partons[4]);
	}
	else {
	  partons[2]->    colourNeighbour(partons[4]);
	  partons[3]->antiColourNeighbour(partons[4]);
	}
      }
      // radiation from final-state
      else {
	if(partons[0]->data().iColour()==PDT::Colour3) {
	  partons[0]->antiColourNeighbour(partons[4]);
	  partons[1]->    colourNeighbour(partons[4]);
	}
	else {
	  partons[0]->    colourNeighbour(partons[4]);
	  partons[1]->antiColourNeighbour(partons[4]);
	}
	if(partons[2]->data().iColour()==PDT::Colour3) {
	  partons[2]->antiColourNeighbour(partons[3]);
	}
	else {
	  partons[2]->    colourNeighbour(partons[3]);
	}
      }
    }
  }
  return partons;
}

void WeakPartonicDecayer::persistentOutput(PersistentOStream & os) const {
  os << MECode << _radprob << _maxtry << _threemax << _fourmax;
}

void WeakPartonicDecayer::persistentInput(PersistentIStream & is, int) {
  is >> MECode >> _radprob >> _maxtry >> _threemax >> _fourmax;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<WeakPartonicDecayer,PartonicDecayerBase>
describeHerwigWeakPartonicDecayer("Herwig::WeakPartonicDecayer", "HwShower.so HwPartonicDecay.so");

void WeakPartonicDecayer::Init() {

  static ClassDocumentation<WeakPartonicDecayer> documentation
    ("The WeakPartonicDecayer class performs partonic decays of hadrons containing a "
     "heavy quark.");

  static Switch<WeakPartonicDecayer,int> interfaceMECode
    ("MECode",
     "The code for the type of matrix element to be used.",
     &WeakPartonicDecayer::MECode, 0, false, false);
  static SwitchOption interfaceMECodePhaseSpace
    (interfaceMECode,
     "PhaseSpace",
     "Phase space decays",
     0);
  static SwitchOption interfaceMECodeWeak
    (interfaceMECode,
     "Weak",
     "Weak matrix element",
     100);

  static Parameter<WeakPartonicDecayer,double> interfaceRadiationProbability
    ("RadiationProbability",
     "The probability that QCD radiation produces an extra q qbar pair",
     &WeakPartonicDecayer::_radprob, 0., 0.0, 1.,
     false, false, Interface::limited);

  static Parameter<WeakPartonicDecayer,unsigned int> interfaceMaxTry
    ("MaxTry",
     "The maximum number of attempts to generate the kinematics",
     &WeakPartonicDecayer::_maxtry, 300, 10, 1000,
     false, false, Interface::limited);

  static Parameter<WeakPartonicDecayer,double> interfaceThreeMax
    ("ThreeMax",
     "Maximum weight for sampling of three-body decays",
     &WeakPartonicDecayer::_threemax, 3.0, 1.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<WeakPartonicDecayer,double> interfaceFourMax
    ("FourMax",
     "Maximum weight for sampling of four-body decays",
     &WeakPartonicDecayer::_fourmax, 3.0, 1.0, 1000.0,
     false, false, Interface::limited);

}

double WeakPartonicDecayer::VAWt(Energy2 t0, Energy2 t1, Energy2 t2, InvEnergy4 t3) {
  return (t1-t0)*(t0-t2)*t3;
}

void WeakPartonicDecayer::dataBaseOutput(ofstream & output,
					 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":MECode " << MECode << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}


double WeakPartonicDecayer::
threeBodyMatrixElement(tcPDPtr dec, const RhoDMatrix & rhoin,
		       Lorentz5Momentum & pdec,
		       ParticleVector & partons) const {
  // spinors
  LorentzSpinor   <SqrtEnergy> w0[2],w2[2];
  LorentzSpinorBar<SqrtEnergy> w1[2],w3[2];
  // spinors for the decaying particle and first product
  if(dec->id()>0) {
    SpinorWaveFunction    win(pdec,dec,0,incoming);
    w0[0] = win.dimensionedWave();
    win.reset(1);
    w0[1] = win.dimensionedWave();
    SpinorBarWaveFunction wout(partons[2]->momentum(),
			       partons[2]->dataPtr(),0,outgoing);
    w1[0] = wout.dimensionedWave();
    wout.reset(1);
    w1[1] = wout.dimensionedWave();
  }
  else {
    SpinorBarWaveFunction win(pdec,dec,0,incoming);
    w1[0] = win.dimensionedWave();
    win.reset(1);
    w1[1] = win.dimensionedWave();
    SpinorWaveFunction wout(partons[2]->momentum(),
			    partons[2]->dataPtr(),0,outgoing);
    w0[0] = wout.dimensionedWave();
    wout.reset(1);
    w0[1] = wout.dimensionedWave();
  }
  // spinors for the W decay products
  bool lorder = true;
  if(partons[0]->id()<0) {
    SpinorWaveFunction    wout2(partons[0]->momentum(),
				partons[0]->dataPtr(),0,outgoing);
    SpinorBarWaveFunction wout3(partons[1]->momentum(),
				partons[1]->dataPtr(),0,outgoing);
    lorder = partons[0]->dataPtr()->charged();
    w2[0] = wout2.dimensionedWave();
    w3[0] = wout3.dimensionedWave();
    wout2.reset(1);
    wout3.reset(1);
    w2[1] = wout2.dimensionedWave();
    w3[1] = wout3.dimensionedWave();
  }
  else {
    SpinorWaveFunction    wout2(partons[1]->momentum(),
				partons[1]->dataPtr(),0,outgoing);
    SpinorBarWaveFunction wout3(partons[0]->momentum(),
				partons[0]->dataPtr(),0,outgoing);
    lorder = partons[1]->dataPtr()->charged();
    w2[0] = wout2.dimensionedWave();
    w3[0] = wout3.dimensionedWave();
    wout2.reset(1);
    wout3.reset(1);
    w2[1] = wout2.dimensionedWave();
    w3[1] = wout3.dimensionedWave();
  }
  bool tau = abs(partons[0]->id())==ParticleID::tauminus || abs(partons[1]->id())==ParticleID::tauminus;
  // calculate the currents
  LorentzPolarizationVectorE Jbc[2][2],Jdec[2][2];
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      if(dec->id()>0)
	Jbc [ix][iy] = w0[ix].leftCurrent(w1[iy]);
      else
	Jbc [ix][iy] = w0[iy].leftCurrent(w1[ix]);
      if(lorder)
	Jdec[ix][iy] = w2[ix].leftCurrent(w3[iy]);
      else
	Jdec[ix][iy] = w2[iy].leftCurrent(w3[ix]);
    }
  }
  // compute the matrix element
  Complex me[2][2][2][2];
  double total=0.;
  for(unsigned int i0=0;i0<2;++i0) {
    for(unsigned int i1=0;i1<2;++i1) {
      for(unsigned int i2=0;i2<2;++i2) {
	for(unsigned int i3=0;i3<2;++i3) {
	  me[i0][i1][i2][i3] = Jbc[i0][i1].dot(Jdec[i2][i3])/sqr(pdec.mass());
	  total += rhoin(i0,i0).real()*norm(me[i0][i1][i2][i3]);
	}
      }
    }
  }
  total *=2.;
  if(tau) {
    RhoDMatrix rho(PDT::Spin1Half);
    for(unsigned int it1=0;it1<2;++it1) {
      for(unsigned int it2=0;it2<2;++it2) {
	for(unsigned int i0=0;i0<2;++i0) {
	  for(unsigned int i1=0;i1<2;++i1) {
	    for(unsigned int i2=0;i2<2;++i2) {
	      rho(it1,it2) += me[i0][i1][it1][i2 ]*conj(me[i0][i1][it2][i2 ]);
	    }
	  }
	}
      }
    }
    // normalize matrix to unit trace
    rho.normalize();
    for(unsigned int ix=0;ix<2;++ix) {
      if(abs(partons[ix]->id())!=ParticleID::tauminus) continue;
      bool loc = partons[ix]->id() < 0;
      // create the spin info object
      FermionSpinPtr spin = new_ptr(FermionSpinInfo(partons[ix]->momentum(),true));
      // assign spinors
      for(unsigned int iy=0;iy<2;++iy) {
	spin->setBasisState(iy, loc ? w2[iy] : w3[iy].bar());
      }
      // assign rho
      spin->rhoMatrix() = rho;
      // assign spin info
      partons[ix]->spinInfo(spin);
    }
  }
  return total;
}

double WeakPartonicDecayer::
fourBodyMatrixElement(Lorentz5Momentum & p0,Lorentz5Momentum & p1,
		      Lorentz5Momentum & p2,Lorentz5Momentum & p3,
		      Lorentz5Momentum & pg, bool Wcol, bool & initial) const {
  Energy2 d01(p0*p1),d02(p0*p2),d03(p0*p3),d0g(p0*pg);
  Energy2 d12(p1*p2),d13(p1*p3),d1g(p1*pg);
  Energy2 d23(p2*p3),d2g(p2*pg),d3g(p3*pg);
  Energy2 m02(sqr(p0.mass())),m12(sqr(p1.mass())),m22(sqr(p2.mass())),
    m32(sqr(p3.mass()));
  Energy2 mei =
    +1./d0g/d1g  *( -d01*d12*d3g+d01*d03*d2g+2*d01*d03*d12 )
    +1./d0g      *( d12*d3g-d03*d12-d02*d03 )
    +1./d1g      *( d12*d13+d03*d2g+d03*d12 )
    +m12/sqr(d1g)*( -d03*d2g-d03*d12 )
    +m02/sqr(d0g)*(  d12*d3g-d03*d12 );
  Energy2 mef = !Wcol ? ZERO :
    +1./d2g/d3g  *( d0g*d12*d23+d03*d1g*d23+2*d03*d12*d23 )
    +1./d2g      *( d03*d1g+d03*d12-d02*d12 )
    +1./d3g      *( d0g*d12-d03*d13+d03*d12 )
    +m32/sqr(d3g)*( -d0g*d12-d03*d12 )
    +m22/sqr(d2g)*( -d03*d1g-d03*d12 );
  initial = mef/(mei+mef)<UseRandom::rnd();
  return 0.5*(mei+mef)/sqr(p0.mass()-p1.mass()-p2.mass()-p3.mass()-pg.mass());
}
#line 1 "./BtoSGammaDecayer.cc"
// -*- C++ -*-
//
// BtoSGammaDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaDecayer class.
//

#include "BtoSGammaDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr BtoSGammaDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr BtoSGammaDecayer::fullclone() const {
  return new_ptr(*this);
}

bool BtoSGammaDecayer::accept(tcPDPtr , const tPDVector & children) const {
  // should be three decay products
  if(children.size()!=3) return false;
  // photon should be last
  if(children[2]->id()!=ParticleID::gamma) return false;
  // strange should be first
  if(abs(children[0]->id())!=ParticleID::s) return false;
  // first and second should form a colour singlet
  if((children[0]->iColour()==PDT::Colour3&&
      children[1]->iColour()==PDT::Colour3bar)||
     (children[1]->iColour()==PDT::Colour3&&
      children[0]->iColour()==PDT::Colour3bar)) return true;
  else return false;
}

ParticleVector BtoSGammaDecayer::decay(const Particle & parent,
				       const tPDVector & prod) const {
  ParticleVector children;
  for(unsigned int ix=0;ix<prod.size();++ix) {
    children.push_back(prod[ix]->produceParticle());
  }
  // momenta of the decay products
  Lorentz5Momentum pout[3],phad;
  pout[0].setMass(children[0]->dataPtr()->constituentMass());
  pout[1].setMass(children[1]->dataPtr()->constituentMass());
  pout[2].setMass(ZERO);
  // first calculate the hadronic mass spectrum
  phad.setMass(_hadronicmass->hadronicMass(parent.mass(),pout[0].mass()+pout[1].mass()));
  // two body decay to hadronic cluster and photon
  double ctheta,phi;
  Kinematics::generateAngles(ctheta,phi);
  Kinematics::twoBodyDecay(parent.momentum(),pout[2].mass(),phad.mass(),
			   ctheta,phi,pout[2],phad);
  // two body decay of the cluster
  Kinematics::generateAngles(ctheta,phi);
  Kinematics::twoBodyDecay(phad,pout[0].mass(),pout[1].mass(),
			   ctheta,phi,pout[0],pout[1]);
  // set momenta of decay products
  for(unsigned int ix=0;ix<3;++ix){children[ix]->setMomentum(pout[ix]);}
  // make the colour connections
  // quark first
  if(children[0]->data().iColour()==PDT::Colour3) {
    children[0]->antiColourNeighbour(children[1]);
    children[1]->colourNeighbour(children[0]);
  }
  // antiquark first
  else {
    children[0]->colourNeighbour(children[1]);
    children[1]->antiColourNeighbour(children[0]);
  }
  return children;
}


void BtoSGammaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _hadronicmass;
}

void BtoSGammaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hadronicmass;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BtoSGammaDecayer,PartonicDecayerBase>
describeHerwigBtoSGammaDecayer("Herwig::BtoSGammaDecayer", "HwShower.so HwPartonicDecay.so");

void BtoSGammaDecayer::Init() {

  static ClassDocumentation<BtoSGammaDecayer> documentation
    ("The BtoSGammaDecayer class performs to the exclusive decay B to s gamma");

  static Reference<BtoSGammaDecayer,BtoSGammaHadronicMass> interfaceHadronicMass
    ("HadronicMass",
     "Pointer to the object computing the hadronic mass spectrum.",
     &BtoSGammaDecayer::_hadronicmass, false, false, true, false, false);

}

void BtoSGammaDecayer::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the PartonicDecayerBase base class
  PartonicDecayerBase::dataBaseOutput(output,false);
  _hadronicmass->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":HadronicMass " 
	 << _hadronicmass->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./PartonicDecayerBase.cc"
// -*- C++ -*-
//
// PartonicDecayerBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonicDecayerBase class.
//

#include "PartonicDecayerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Hadronization/CluHadConfig.h"
#include "Herwig/Shower/QTilde/Base/ShowerParticle.h"

using namespace Herwig;

PartonicDecayerBase::PartonicDecayerBase() : _exclusive(true), 
					     _partontries(100), _inter(false),
					     _shower(false)
{}

void PartonicDecayerBase::doinit() {
  HwDecayerBase::doinit();
  if(_shower) {
    if(!_partnerFinder) {
      throw InitException() << "PartonicDecayerBase::doinit() " << fullName()
			    << " must have the PartonFinder set if"
			    << " showers are being performecd in partonic decays"
			    << Exception::runerror;
    }
    if(!_splittingGenerator) {
      throw InitException() << "PartonicDecayerBase::doinit() " << fullName()
			    << " must have the SplittingGenerator set if"
			    << " showers are being performecd in partonic decays"
			    << Exception::runerror;
    }
    if(!_reconstructor) {
      throw InitException() << "PartonicDecayerBase::doinit() " << fullName()
			    << " must have the KinematicsReconstructor set if"
			    << " showers are being performecd in partonic decays"
			    << Exception::runerror;
    }

  }
}

ParticleVector PartonicDecayerBase::decay(const DecayMode & dm,
					  const Particle & p) const {
  // handling of the decay including the special features of the
  // DecayMode  
  // get the primary products
  tPDVector products=dm.orderedProducts();
  // add products for which the decay mode is all ready specified
  if(!dm.cascadeProducts().empty()) 
    throw Exception() << "PartonicDecayerBase::decay() cannot handle"
		      << " cascadeProducts() from the DecayMode "
		      << Exception::runerror;
  unsigned int ptry(0);
  bool hadronized(false);
  ParticleVector outpart;
  ParticleVector outhad;
  // copy of particle to act as parent so hadronization konws about it
  PPtr ptemp(new_ptr(Particle(p)));
  do {
    hadronized=false;
    // increment loop counter
    ++ptry;
    // perform the primary decay
    ParticleVector partons=decay(p,products);
    // must add partons are children so reshuffling vs leptons will work
    for(unsigned int ix=0;ix<partons.size();++ix) ptemp->addChild(partons[ix]);
    PVector currentlist = partons;
    // perform the shower if needed
    if(_shower) {
      currentlist = shower(p, partons);
      if(currentlist.empty()) continue;
    }
    // split the gluons
    _partonSplitter->split(currentlist);
    // form the clusters
    ClusterVector clusters = _clusterFinder->formClusters(currentlist);
    _clusterFinder->reduceToTwoComponents(clusters);
    tPVector finalHadrons = _clusterFissioner->fission(clusters,false);
    bool lightOK = _lightClusterDecayer->decay(clusters,finalHadrons);
    // abandon child here so always done
    for(unsigned int ix=0;ix<partons.size();++ix) ptemp->abandonChild(partons[ix]);
    // try again if can't reshuffle
    if(!lightOK) continue;
    // decay the remaining clusters
    _clusterDecayer->decay(clusters,finalHadrons);
    hadronized = !duplicateMode(p,finalHadrons);
    if(hadronized) {
      outhad  = ParticleVector(finalHadrons.begin(),finalHadrons.end());
      outpart = partons;
    }
  }
  while(!hadronized&&ptry<_partontries);
  if(ptry>=_partontries) return ParticleVector();
  // return decay products
  if(_inter) {
    return outpart;
  }
  else {
    for(unsigned int ix=0;ix<outhad.size();++ix) {
      tParticleVector parents=outhad[ix]->parents();
      for(unsigned int iy=0;iy<parents.size();++iy) {
	parents[iy]->abandonChild(outhad[ix]);
      }
    }
    for(unsigned int ix=0;ix<outpart.size();++ix) {
      if(!outpart[ix]->coloured()) outhad.push_back(outpart[ix]);
    }
    return outhad;
  }
}

void PartonicDecayerBase::persistentOutput(PersistentOStream & os) const {
  os << _partonSplitter << _clusterFinder << _clusterFissioner
    << _lightClusterDecayer << _clusterDecayer << _exclusive << _partontries
     << _inter << _shower << _splittingGenerator << _partnerFinder << _reconstructor;
}

void PartonicDecayerBase::persistentInput(PersistentIStream & is, int) {
  is >> _partonSplitter >> _clusterFinder >> _clusterFissioner
    >> _lightClusterDecayer >> _clusterDecayer >> _exclusive >> _partontries
     >> _inter >> _shower >> _splittingGenerator >> _partnerFinder >> _reconstructor;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<PartonicDecayerBase,HwDecayerBase>
describeHerwigPartonicDecayerBase("Herwig::PartonicDecayerBase", "HwShower.so HwPartonicDecay.so");

void PartonicDecayerBase::Init() {

  static ClassDocumentation<PartonicDecayerBase> documentation
    ("The PartonicDecayerBase class is the base class for partonic decays"
     " and handles the hadronization of the decays");

  static Reference<PartonicDecayerBase,PartonSplitter> 
    interfacePartonSplitter("PartonSplitter", 
		      "A reference to the PartonSplitter object", 
		      &Herwig::PartonicDecayerBase::_partonSplitter,
		      false, false, true, false);

  static Reference<PartonicDecayerBase,ClusterFinder> 
    interfaceClusterFinder("ClusterFinder", 
		      "A reference to the ClusterFinder object", 
		      &Herwig::PartonicDecayerBase::_clusterFinder,
		      false, false, true, false);

  static Reference<PartonicDecayerBase,ClusterFissioner> 
    interfaceClusterFissioner("ClusterFissioner", 
		      "A reference to the ClusterFissioner object", 
		      &Herwig::PartonicDecayerBase::_clusterFissioner,
		      false, false, true, false);

  static Reference<PartonicDecayerBase,LightClusterDecayer> 
    interfaceLightClusterDecayer("LightClusterDecayer", 
		    "A reference to the LightClusterDecayer object", 
		    &Herwig::PartonicDecayerBase::_lightClusterDecayer,
		    false, false, true, false);

  static Reference<PartonicDecayerBase,ClusterDecayer> 
    interfaceClusterDecayer("ClusterDecayer", 
		       "A reference to the ClusterDecayer object", 
		       &Herwig::PartonicDecayerBase::_clusterDecayer,
		       false, false, true, false);

  static Switch<PartonicDecayerBase,bool> interface_exclusive
    ("Exclusive",
     "Ensure that the hadrons produced in the partonic decays of bottom"
     " and charm baryons do not duplicate the inclusive modes.",
     &PartonicDecayerBase::_exclusive, true, false, false);
  static SwitchOption interface_exclusiveNoDuplication
    (interface_exclusive,
     "Yes",
     "Forbid duplication",
     true);
  static SwitchOption interface_exclusiveDuplication
    (interface_exclusive,
     "No",
     "Duplication allowed",
     false);
  
  static Switch<PartonicDecayerBase,bool> interfaceIntermediates
    ("Intermediates",
     "Whether or not to include the intermediate particles produced by the"
     " cluster alogorithm in the event record.",
     &PartonicDecayerBase::_inter, false, false, false);
  static SwitchOption interfaceIntermediatesIntermediates
    (interfaceIntermediates,
     "Yes",
     "Include the intermediates",
     true);
  static SwitchOption interfaceIntermediatesNoIntermediates
    (interfaceIntermediates,
     "No",
     "Don't include the intermediates.",
     false);

  static Parameter<PartonicDecayerBase,unsigned int> interfacePartonic_Tries
    ("Partonic_Tries",
     "Number of attempts to generator the hadronisation of the decay",
     &PartonicDecayerBase::_partontries, 100, 1, 1000,
     false, false, Interface::limited);
  
  static Switch<PartonicDecayerBase,bool> interfaceShower
    ("Shower",
     "Whether or not to perform the parton shower before the hadronization",
     &PartonicDecayerBase::_shower, false, false, false);
  static SwitchOption interfaceShowerYes
    (interfaceShower,
     "Yes",
     "Perform the shower",
     true);
  static SwitchOption interfaceShowerNo
    (interfaceShower,
     "No",
     "Don't perform the shower",
     false);

  static Reference<PartonicDecayerBase,PartnerFinder> interfacePartnerFinder
    ("PartnerFinder",
     "The parton finder",
     &PartonicDecayerBase::_partnerFinder, false, false, true, true, false);

  static Reference<PartonicDecayerBase,SplittingGenerator> interfaceSplittingGenerator
    ("SplittingGenerator",
     "The splitting generator",
     &PartonicDecayerBase::_splittingGenerator, false, false, true, true, false);
  
  static Reference<PartonicDecayerBase,KinematicsReconstructor> interfaceKinematicsReconstructor
    ("KinematicsReconstructor",
     "The kinematics reconstructor",
     &PartonicDecayerBase::_reconstructor, false, false, true, true, false);

}

bool PartonicDecayerBase::duplicateMode(const Particle & parent,
					const vector<tPPtr> & hadrons) const {
  // if not exclusive return
  if(!_exclusive) return false;
  // now find the hadrons and check them
  cParticleMSet hadronsb;
  bool found(false);
  for (unsigned ix = 0; ix < hadrons.size(); ++ix)
    hadronsb.insert(hadrons[ix]->dataPtr());
  // now check particle's decay modes 
  Selector<tDMPtr>::const_iterator modeptr 
    = parent.dataPtr()->decaySelector().begin();
  Selector<tDMPtr>::const_iterator end     
    = parent.dataPtr()->decaySelector().end();
  // check not a duplicate of a known mode
  for(;modeptr!=end;++modeptr) {
    tcDMPtr mode=(*modeptr).second;
    // check same number of products
    if(mode->products().size() != hadronsb.size()) continue;
    ParticleMSet::const_iterator dit;
    cParticleMSet::const_iterator pit;
    for(dit=mode->products().begin(), pit=hadronsb.begin();
	dit!=mode->products().end(); ++dit,++pit) {
      if((*dit)!=(*pit)) break;
    }
    if(dit != mode->products().end()) continue;
    found = true;
    break;
  }
  return found;
}

void PartonicDecayerBase::dataBaseOutput(ofstream & output,bool header) const {
  // header for MySQL
  if(header) output << "update decayers set parameters=\"";
  // parameters
  output << "newdef  " << name() << ":PartonSplitter " 
	 << _partonSplitter->name() << " \n";
  output << "newdef  " << name() << ":ClusterFinder " 
	 << _clusterFinder->name() << " \n";
  output << "newdef  " << name() << ":ClusterFissioner " 
	 << _clusterFissioner->name() << " \n";
  output << "newdef  " << name() << ":LightClusterDecayer " 
	 << _lightClusterDecayer->name() << " \n";
  output << "newdef  " << name() << ":ClusterDecayer " 
	 << _clusterDecayer->name() << " \n";
  output << "newdef  " << name() << ":Exclusive " <<  _exclusive<< " \n";
  output << "newdef  " << name() << ":Intermediates " << _inter << " \n";
  output << "newdef  " << name() << ":Shower " << _shower << " \n";
  if(_splittingGenerator)
    output << "newdef  " << name() << ":SplittingGenerator " << _splittingGenerator->fullName() << " \n";
  if(_partnerFinder)
    output << "newdef  " << name() << ":PartonFinder " << _partnerFinder->fullName() << " \n";
  if(_reconstructor)
    output << "newdef  " << name() << ":KinematicsReconstructor " << _reconstructor->fullName() << " \n";
  output << "newdef  " << name() << ":Shower " << _shower << " \n";
  output << "newdef  " << name() << ":Partonic_Tries " << _partontries << " \n";
  // footer for MySQL
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

PVector PartonicDecayerBase::shower(const Particle & parent,
				    ParticleVector & partons) const {
  try {
    // set up the shower
    // create the shower particles
    ShowerParticleVector particles;
    for(auto p : partons) {
      particles.push_back(new_ptr(ShowerParticle(*p,2,true)));
    }
    // set up the partners and evolution scales
    _partnerFinder->setInitialEvolutionScales(particles,true,ShowerInteraction::QCD,false);
    // do the shower
    PVector output;
    for(tShowerParticlePtr p : particles) {
      p->initializeFinalState();
      timeLikeShower(p,Branching(),true,output);
    }
    // check if in CMF frame
    Boost beta_cm = parent.momentum().findBoostToCM();
    bool gottaBoost(false);
    LorentzRotation trans = LorentzRotation();
    if(beta_cm.mag() > 1e-12) {
      gottaBoost = true;
      trans.boost(beta_cm);
    }
    bool radiated(false);
    vector<JetKinStruct> jetKinematics;
    for(unsigned int ix=0;ix<particles.size();++ix) {
      JetKinStruct tempJetKin;      
      tempJetKin.parent = particles[ix];
      if(gottaBoost) {
	_reconstructor->deepTransform(tempJetKin.parent,trans);
      }
      tempJetKin.p = particles[ix]->momentum();
      radiated |= _reconstructor->reconstructTimeLikeJet(particles[ix],particles[ix]);
      tempJetKin.q = particles[ix]->momentum();
      jetKinematics.push_back(tempJetKin);
    }
    // find the rescaling factor
    double k = 0.0;
    if(radiated) {
      k = _reconstructor->solveKfactor(parent.mass(), jetKinematics);
      // perform the rescaling and boosts
      for(const JetKinStruct & jet : jetKinematics) {
	LorentzRotation Trafo = _reconstructor->solveBoost(k, jet.q, jet.p);
	_reconstructor->deepTransform(jet.parent,Trafo);
      }
    }
    if(gottaBoost) {
      LorentzRotation rinv = trans.inverse();
      for(const JetKinStruct & jet : jetKinematics) {
	_reconstructor->deepTransform(jet.parent,rinv);
      }
    }
    // add shower particles as children
    for(unsigned int ix=0;ix<partons.size();++ix)
      partons[ix]->addChild(particles[ix]);
    // return final-state of shower
    return output;
  }
  catch (KinematicsReconstructionVeto & ) {
    return PVector();
  }
}

bool PartonicDecayerBase::timeLikeShower(tShowerParticlePtr particle,
					 Branching fb, bool first, PVector & output) const {
  // generate the emission
  if(!fb.kinematics) {
    fb=_splittingGenerator->chooseForwardBranching(*particle,1.,ShowerInteraction::QCD);
    if(fb.kinematics) fb.hard = false;
  }
  // no emission, return
  if(!fb.kinematics) {
    if(particle->spinInfo()) particle->spinInfo()->develop();
    output.push_back(particle);
    return false;
  }
  Branching fc[2] = {Branching(),Branching()};
  // has emitted
  // Assign the shower kinematics to the emitting particle.
  particle->showerKinematics(fb.kinematics);
  // create the children
  ShowerParticleVector children;
  children.reserve(2);
  for(unsigned int ix=0;ix<fb.ids.size()-1;++ix) {
    children.push_back(new_ptr(ShowerParticle(fb.ids[ix+1],true)));
    children[ix]->set5Momentum(Lorentz5Momentum(fb.ids[ix+1]->mass()));
  }
  // update the children
  particle->showerKinematics()->updateChildren(particle, children,1,fb.type);
  // select branchings for children
  for(unsigned int ix=0;ix<children.size()-1;++ix) {
    fc[ix] = _splittingGenerator->chooseForwardBranching(*children[ix],1.,ShowerInteraction::QCD);
  }
  // shower the children
  for(unsigned int ix=0;ix<children.size();++ix) {
    if(fc[ix].kinematics) {
      timeLikeShower(children[ix],fc[ix],false,output);
    }
    else {
      output.push_back(children[ix]);
    }
    if(children[ix]->spinInfo()) children[ix]->spinInfo()->develop();
  }
  particle->showerKinematics()->updateParent(particle, children,1,fb.type);
  // branching has happened
  if(first&&!children.empty())
    particle->showerKinematics()->resetChildren(particle,children);
  if(particle->spinInfo()) particle->spinInfo()->develop();
  return true;
}
