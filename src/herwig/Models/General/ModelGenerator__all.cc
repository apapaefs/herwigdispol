#line 1 "./ModelGenerator.cc"
// -*- C++ -*-
//
// ModelGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ModelGenerator class.
//

#include "ModelGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "BSMWidthGenerator.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "Herwig/Decay/DecayIntegrator.h"
#include "ThePEG/Repository/BaseRepository.h"

using namespace Herwig;

IBPtr ModelGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr ModelGenerator::fullclone() const {
  return new_ptr(*this);
}

void ModelGenerator::persistentOutput(PersistentOStream & os) const {
  os << hardProcessConstructors_ << _theDecayConstructor << particles_ 
     << offshell_ << Offsel_ << BRnorm_ << twoBodyOnly_ << howOffShell_
     << Npoints_ << Iorder_ << BWshape_ << brMin_ << decayOutput_;
}

void ModelGenerator::persistentInput(PersistentIStream & is, int) {
  is >> hardProcessConstructors_ >> _theDecayConstructor >> particles_
     >> offshell_ >> Offsel_ >> BRnorm_ >> twoBodyOnly_ >> howOffShell_
     >> Npoints_ >> Iorder_ >> BWshape_ >> brMin_ >> decayOutput_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<ModelGenerator,Interfaced>
describeThePEGModelGenerator("Herwig::ModelGenerator", "Herwig.so");


void ModelGenerator::Init() {

  static ClassDocumentation<ModelGenerator> documentation
    ("This class controls the the use of BSM physics.",
     "BSM physics was produced using the algorithm of "
     "\\cite{Gigg:2007cr,Gigg:2008yc}",
     "\\bibitem{Gigg:2007cr} M.~Gigg and P.~Richardson, \n"
     "Eur.\\ Phys.\\ J.\\  C {\\bf 51} (2007) 989.\n"
     "%%CITATION = EPHJA,C51,989;%%\n"
     " %\\cite{Gigg:2008yc}\n"
     "\\bibitem{Gigg:2008yc}\n"
     "  M.~A.~Gigg and P.~Richardson,\n"
     "  %``Simulation of Finite Width Effects in Physics Beyond the Standard Model,''\n"
     "  arXiv:0805.3037 [hep-ph].\n"
     "  %%CITATION = ARXIV:0805.3037;%%\n"
     );
 
  static RefVector<ModelGenerator,HardProcessConstructor> 
    interfaceHardProcessConstructors
    ("HardProcessConstructors",
     "The objects to construct hard processes",
     &ModelGenerator::hardProcessConstructors_, -1, 
     false, false, true, false, false);

  static Reference<ModelGenerator,Herwig::DecayConstructor> 
     interfaceDecayConstructor
     ("DecayConstructor",
      "Pointer to DecayConstructor helper class",
      &ModelGenerator::_theDecayConstructor, false, false, true, false);
  
  static RefVector<ModelGenerator,ThePEG::ParticleData> interfaceModelParticles
    ("DecayParticles",
     "ParticleData pointers to the particles requiring spin correlation "
     "decayers. If decay modes do not exist they will also be created.",
     &ModelGenerator::particles_, -1, false, false, true, false);
    
  static RefVector<ModelGenerator,ParticleData> interfaceOffshell
    ("Offshell",
     "The particles to treat as off-shell",
     &ModelGenerator::offshell_, -1, false, false, true, false);

  static Switch<ModelGenerator,int> interfaceWhichOffshell
    ("WhichOffshell",
     "A switch to determine which particles to create mass and width "
     "generators for.",
     &ModelGenerator::Offsel_, 0, false, false);
  static SwitchOption interfaceWhichOffshellSelected
    (interfaceWhichOffshell,
     "Selected",
     "Only create mass and width generators for the particles specified",
     0);
  static SwitchOption interfaceWhichOffshellAll
    (interfaceWhichOffshell,
     "All",
     "Treat all particles specified in the DecayParticles "
     "list as off-shell",
     1);
  
  static Switch<ModelGenerator,bool> interfaceBRNormalize
    ("BRNormalize",
     "Whether to normalize the partial widths to BR*total width for an "
     "on-shell particle",
     &ModelGenerator::BRnorm_, true, false, false);
  static SwitchOption interfaceBRNormalizeNormalize
    (interfaceBRNormalize,
     "Yes",
     "Normalize the partial widths",
     true);
  static SwitchOption interfaceBRNormalizeNoNormalize
    (interfaceBRNormalize,
     "No",
     "Do not normalize the partial widths",
     false);

  static Parameter<ModelGenerator,int> interfacePoints
    ("InterpolationPoints",
     "Number of points to use for interpolation tables when needed",
     &ModelGenerator::Npoints_, 10, 5, 1000,
     false, false, true);
  
  static Parameter<ModelGenerator,unsigned int> 
    interfaceInterpolationOrder
    ("InterpolationOrder", "The interpolation order for the tables",
     &ModelGenerator::Iorder_, 1, 1, 5,
     false, false, Interface::limited);

  static Switch<ModelGenerator,int> interfaceBreitWignerShape
    ("BreitWignerShape",
     "Controls the shape of the mass distribution generated",
     &ModelGenerator::BWshape_, 0, false, false);
  static SwitchOption interfaceBreitWignerShapeDefault
    (interfaceBreitWignerShape,
     "Default",
     "Running width with q in numerator and denominator width factor",
     0);
  static SwitchOption interfaceBreitWignerShapeFixedWidth
    (interfaceBreitWignerShape,
     "FixedWidth",
     "Use a fixed width",
     1);
  static SwitchOption interfaceBreitWignerShapeNoq
    (interfaceBreitWignerShape,
     "Noq",
     "Use M rather than q in the numerator and denominator width factor",
     2);
  static SwitchOption interfaceBreitWignerShapeNoNumerator
    (interfaceBreitWignerShape,
     "NoNumerator",
     "Neglect the numerator factors",
     3);

  static Switch<ModelGenerator,bool> interfaceTwoBodyOnly
    ("TwoBodyOnly",
     "Whether to use only two-body or all modes in the running width calculation",
     &ModelGenerator::twoBodyOnly_, false, false, false);
  static SwitchOption interfaceTwoBodyOnlyYes
    (interfaceTwoBodyOnly,
     "Yes",
     "Only use two-body modes",
     true);
  static SwitchOption interfaceTwoBodyOnlyNo
    (interfaceTwoBodyOnly,
     "No",
     "Use all modes",
     false);
  
  static Parameter<ModelGenerator,double> interfaceMinimumBR
    ("MinimumBR",
     "The minimum branching fraction to include",
     &ModelGenerator::brMin_, 1e-6, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<ModelGenerator,unsigned int> interfaceDecayOutput
    ("DecayOutput",
     "Option to control the output of the decay mode information",
     &ModelGenerator::decayOutput_, 1, false, false);
  static SwitchOption interfaceDecayOutputNone
    (interfaceDecayOutput,
     "None",
     "No output",
     0);
  static SwitchOption interfaceDecayOutputPlain
    (interfaceDecayOutput,
     "Plain",
     "Default plain text output",
     1);
  static SwitchOption interfaceDecayOutputSLHA
    (interfaceDecayOutput,
     "SLHA",
     "Output in the Susy Les Houches Accord format",
     2);

  static Parameter<ModelGenerator,double> interfaceMinimumWidthFraction
    ("MinimumWidthFraction",
     "Minimum fraction of the particle's mass the width can be"
     " for the off-shell treatment.",
     &ModelGenerator::minWidth_, 1e-6, 1e-15, 1.,
     false, false, Interface::limited);

  static Parameter<ModelGenerator,double> interfaceHowMuchOffShell
    ("HowMuchOffShell",
     "The multiple of the particle's width by which it is allowed to be off-shell",
     &ModelGenerator::howOffShell_, 5., 0.0, 100.,
     false, false, Interface::limited);

}

namespace {
  
  /// Helper function for sorting by mass
  inline bool massIsLess(tcPDPtr a, tcPDPtr b) {
    return a->mass() < b->mass();
  }

  // Helper function to find minimum possible mass of a particle
  inline Energy minimumMass(tcPDPtr parent) {
    Energy output(Constants::MaxEnergy);
    for(set<tDMPtr>::const_iterator dit = parent->decayModes().begin();
	dit != parent->decayModes().end(); ++dit) {
      Energy outMass(ZERO);
      for(unsigned int ix=0;ix<(**dit).orderedProducts().size();++ix) {
	outMass += (**dit).orderedProducts()[ix]->massMin();
      }
      output = min(output,outMass);
    }
    return output;
  }
}

void ModelGenerator::doinit() {
  useMe();
  Interfaced::doinit();
  // make sure the model is initialized
  Ptr<Herwig::StandardModel>::pointer model 
    = dynamic_ptr_cast<Ptr<Herwig::StandardModel>::pointer>(generator()->standardModel());
  model->init();
  // and the vertices
  for(size_t iv = 0; iv < model->numberOfVertices(); ++iv)
    model->vertex(iv)->init();
  // uniq and sort DecayParticles list by mass
  set<PDPtr> tmp(particles_.begin(),particles_.end());
  particles_.assign(tmp.begin(),tmp.end());
  sort(particles_.begin(),particles_.end(),massIsLess);
  //create decayers and decaymodes (if necessary)
  if( _theDecayConstructor ) {
    _theDecayConstructor->init();
    _theDecayConstructor->createDecayers(particles_,brMin_);
  }

  // write out decays with spin correlations
  ostream & os = CurrentGenerator::current().misc();
  ofstream ofs;
  if ( decayOutput_ > 1 ) {
    string filename 
      = CurrentGenerator::current().filename() + "-BR.spc";
    ofs.open(filename.c_str());
  }

  if(decayOutput_!=0) {
    if(decayOutput_==1) {
      os << "# The decay modes listed below will have spin\n"
	  << "# correlations included when they are generated.\n#\n#";
    }
    else {
      ofs << "#  Herwig decay tables in SUSY Les Houches accord format\n";
      ofs << "Block DCINFO                           # Program information\n";
      ofs << "1   Herwig          # Decay Calculator\n";
      ofs << "2   " << generator()->strategy()->versionstring() 
	  << "     # Version number\n";
    }
  }
  //create mass and width generators for the requested particles
  set<PDPtr> offShell;
  if( Offsel_ == 0 ) offShell = set<PDPtr>(offshell_.begin() ,offshell_.end() );
  else               offShell = set<PDPtr>(particles_.begin(),particles_.end());
  
  for(PDVector::iterator pit = particles_.begin(); 
      pit != particles_.end(); ++pit) {
    tPDPtr parent = *pit;
    // Check decays for ones where quarks cannot be put on constituent
    // mass-shell
    checkDecays(parent);
    parent->reset();
    // Now switch off the modes if needed
    for(DecaySet::const_iterator it=parent->decayModes().begin();
	it!=parent->decayModes().end();++it) {
      if( _theDecayConstructor->disableDecayMode((**it).tag()) ) {
	DMPtr mode = *it;
	generator()->preinitInterface(mode, "Active", "set", "No");
	if(mode->CC()) {
	  DMPtr CCmode = mode->CC();
	  generator()->preinitInterface(CCmode, "Active", "set", "No");
	}
      }
    }
    parent->update();
    if( parent->CC() ) parent->CC()->synchronize();
    if( parent->decaySelector().empty() ) {
      parent->stable(true);
      parent->width(ZERO);
      parent->widthCut(ZERO);
      parent->massGenerator(tGenericMassGeneratorPtr());
      parent->widthGenerator(tGenericWidthGeneratorPtr());
    }
    else {
      if(parent->mass()*minWidth_>parent->width()) {
	parent->massGenerator(tGenericMassGeneratorPtr());
	parent->widthGenerator(tGenericWidthGeneratorPtr());
      }
      else {
	if( offShell.find(*pit) != offShell.end() ) {
	  createWidthGenerator(*pit);
	}
	else {
	  parent->massGenerator(tGenericMassGeneratorPtr());
	  parent->widthGenerator(tGenericWidthGeneratorPtr());
	}
      }

    }
    // set the offshellness 
    Energy minMass = minimumMass(parent);
    Energy offShellNess = howOffShell_*parent->width();
    if(minMass>parent->mass()-offShellNess) {
      offShellNess = parent->mass()-minMass;
    }
    parent->widthCut(offShellNess);
    
    if( parent->massGenerator() ) {

      parent->massGenerator()->reset();
      if(decayOutput_==1)
	os << "# " <<parent->PDGName() << " will be considered off-shell.\n#\n";
    }
    if( parent->widthGenerator() ) parent->widthGenerator()->reset();
  }
  // loop again to initialise mass and width generators
  // switch off modes and write output
  for(PDVector::iterator pit = particles_.begin();
      pit != particles_.end(); ++pit) {
    tPDPtr parent = *pit;
    if(parent->widthGenerator())
      parent->widthGenerator()->init();
    if(parent->massGenerator())
      parent->massGenerator()->init();
    // output the modes if needed
    if( !parent->decaySelector().empty() ) {
      if ( decayOutput_ == 2 )
	writeDecayModes(ofs, parent);
      else
	writeDecayModes(os, parent);
    }
  }

  //Now construct hard processes given that we know which
  //objects have running widths
  for(unsigned int ix=0;ix<hardProcessConstructors_.size();++ix) {
    hardProcessConstructors_[ix]->init();
    hardProcessConstructors_[ix]->constructDiagrams();
  }
}

void ModelGenerator::checkDecays(PDPtr parent) {
  if( parent->stable() ) {
    if(parent->coloured())
      cerr << "Warning: No decays for coloured particle " << parent->PDGName() << "\n\n" 
	   << "have been calcluated in BSM model.\n"
	   << "This may cause problems in the hadronization phase.\n"
	   << "You may have forgotten to switch on the decay mode calculation using\n"
	   << "  set TwoBodyDC:CreateDecayModes Yes\n"
	   << "  set ThreeBodyDC:CreateDecayModes Yes\n"
	   << "  set WeakDecayConstructor:CreateDecayModes Yes\n"
	   << "or the decays of this particle are missing from your\n"
	   << "input spectrum and decay file in the SLHA format.\n\n";
    return;
  }
  DecaySet::iterator dit = parent->decayModes().begin();
  DecaySet::iterator dend = parent->decayModes().end();
  Energy oldwidth(parent->width()), newwidth(ZERO);
  bool rescalebrat(false);
  double brsum(0.);
  for(; dit != dend; ++dit ) {
    if( !(**dit).on() ) continue;
    Energy release((**dit).parent()->mass());
    tPDVector::const_iterator pit = (**dit).orderedProducts().begin();
    tPDVector::const_iterator pend =(**dit).orderedProducts().end();
    for( ; pit != pend; ++pit ) {
      release -= (**pit).constituentMass();
    }
    if( (**dit).brat() < brMin_ || release < ZERO ) {
      if( release < ZERO )
	cerr << "Warning: The shower cannot be generated using this decay " 
	     << (**dit).tag() << " because it is too close to threshold.\nIt "
	     << "will be switched off and the branching fractions of the "
	     << "remaining modes rescaled.\n";
      rescalebrat = true;
      generator()->preinitInterface(*dit, "Active", "set", "No");
      generator()->preinitInterface(*dit, "BranchingRatio", 
				    "set", "0.0");
      DecayIntegratorPtr decayer = dynamic_ptr_cast<DecayIntegratorPtr>((**dit).decayer());
      if(decayer) {
      	generator()->preinitInterface(decayer->fullName(), "Initialize", "set","0");
      }
    }
    else {
      brsum += (**dit).brat();
      newwidth += (**dit).brat()*oldwidth;
    }
  }
  // if no modes left set stable
  if(newwidth==ZERO) {
    parent->stable(true);
    parent->width(ZERO);
    parent->widthCut(ZERO);
    parent->massGenerator(tGenericMassGeneratorPtr());
    parent->widthGenerator(tGenericWidthGeneratorPtr());
  }
  // otherwise rescale if needed
  else if( ( rescalebrat || abs(brsum - 1.) > 1e-12 ) && !parent->decayModes().empty()) {
    dit = parent->decayModes().begin();
    dend = parent->decayModes().end();
    double factor = oldwidth/newwidth;
    brsum = 0.;
    for( ; dit != dend; ++dit ) {
      if( !(**dit).on() ) continue;
      double newbrat = ((**dit).brat())*factor;
      brsum += newbrat;
      ostringstream brf;
      brf << setprecision(13) << newbrat;
      generator()->preinitInterface(*dit, "BranchingRatio",
				    "set", brf.str());
    }
    parent->width(newwidth);
    if( newwidth > ZERO ) parent->cTau(hbarc/newwidth);
  }
}

namespace {
  struct DecayModeOrdering {
    bool operator()(tcDMPtr m1, tcDMPtr m2) const {
      if(m1->brat()!=m2->brat()) {
	return m1->brat()>m2->brat();
      }
      else {
	if(m1->products().size()==m2->products().size()) {
	  ParticleMSet::const_iterator it1=m1->products().begin();
	  ParticleMSet::const_iterator it2=m2->products().begin();
	  do {
	    if((**it1).id()!=(**it2).id()) {
	      return (**it1).id()>(**it2).id();
	    }
	    ++it1;
	    ++it2;
	  }
	  while(it1!=m1->products().end()&&
		it2!=m2->products().end());
	  assert(false);
	}
	else
	  return m1->products().size()<m2->products().size();
      }
      return false;
    }
  };
}

void ModelGenerator::writeDecayModes(ostream & os, tcPDPtr parent) const {
  if(decayOutput_==0) return;
  set<tcDMPtr,DecayModeOrdering> modes(parent->decayModes().begin(),
				       parent->decayModes().end());
  if(decayOutput_==1) {
    os << " Parent: " << parent->PDGName() << "  Mass (GeV): " 
       << parent->mass()/GeV << "  Total Width (GeV): " 
       << parent->width()/GeV << endl;
    os << std::left << std::setw(40) << '#' 
       << std::left << std::setw(20) << "Partial Width/GeV"
       << std::left << std::setw(20) << "BR" << "Yes/No\n";
    for(set<tcDMPtr,DecayModeOrdering>::iterator dit=modes.begin();
	dit!=modes.end();++dit)
      os << std::left << std::setw(40) << (**dit).tag() 
	 << std::left << std::setw(20) << (**dit).brat()*parent->width()/GeV 
	 << std::left << std::setw(20)  << (**dit).brat()
	 << ((**dit).on() ? "Yes" : "No" ) << '\n';
    os << "#\n#";
  }
  else if(decayOutput_==2) {
    os << "#    \t PDG \t Width\n";
    os << "DECAY\t" << parent->id() << "\t" << parent->width()/GeV << "\t # " << parent->PDGName() << "\n";
    for(set<tcDMPtr,DecayModeOrdering>::iterator dit=modes.begin();
	dit!=modes.end();++dit) {
      os << "\t" << std::left << std::setw(10) 
	 << (**dit).brat() << "\t" << (**dit).orderedProducts().size() 
	 << "\t";
      for(unsigned int ix=0;ix<(**dit).orderedProducts().size();++ix)
	os << std::right << std::setw(10)
	   << (**dit).orderedProducts()[ix]->id() ;
      for(unsigned int ix=(**dit).orderedProducts().size();ix<4;++ix)
	os << "\t";
      os << "# " << (**dit).tag() << "\n";
    }
  }
}

void ModelGenerator::createWidthGenerator(tPDPtr p) {
  string wn = p->fullName() + string("-WGen");
  string mn = p->fullName() + string("-MGen");
  GenericMassGeneratorPtr mgen = dynamic_ptr_cast<GenericMassGeneratorPtr>
    (generator()->preinitCreate("Herwig::GenericMassGenerator", mn));
  BSMWidthGeneratorPtr wgen = dynamic_ptr_cast<BSMWidthGeneratorPtr>
    (generator()->preinitCreate("Herwig::BSMWidthGenerator", wn));

  //set the particle interface
  mgen->particle(p);
  wgen->particle(p);

  //set the generator interfaces in the ParticleData object
  generator()->preinitInterface(p, "Mass_generator","set", mn);
  generator()->preinitInterface(p, "Width_generator","set", wn);
  //allow the branching fraction of this particle type to vary
  p->variableRatio(true);
  if( p->CC() ) p->CC()->variableRatio(true);
  
  //initialize the generators
  generator()->preinitInterface(mgen, "Initialize", "set", "Yes");
  generator()->preinitInterface(wgen, "Initialize", "set", "Yes");

  string norm = BRnorm_ ? "Yes" : "No";
  generator()->preinitInterface(wgen, "BRNormalize", "set", norm);
  string twob = twoBodyOnly_ ? "Yes" : "No";
  generator()->preinitInterface(wgen, "TwoBodyOnly", "set", twob);
  ostringstream os;
  os << Npoints_;
  generator()->preinitInterface(wgen, "Points", "set", os.str());
  os.str("");
  os << Iorder_;
  generator()->preinitInterface(wgen, "InterpolationOrder", "set",
				  os.str());
  os.str("");
  os << BWshape_;
  generator()->preinitInterface(mgen, "BreitWignerShape", "set", 
				  os.str());
}
#line 1 "./DecayConstructor.cc"
// -*- C++ -*-
//
// DecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayConstructor class.
//

#include "DecayConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/BaseRepository.h"
#include <iterator>

using namespace Herwig;
using namespace ThePEG;

IBPtr DecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr DecayConstructor::fullclone() const {
  return new_ptr(*this);
}
  
void DecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << NBodyDecayConstructors_ << QEDGenerator_;
}

void DecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> NBodyDecayConstructors_ >> QEDGenerator_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DecayConstructor,Interfaced>
describeHerwigDecayConstructor("Herwig::DecayConstructor", "Herwig.so");

void DecayConstructor::Init() {

  static ClassDocumentation<DecayConstructor> documentation
    ("There is no documentation for the TwoBodyDecayConstructor class");

  static RefVector<DecayConstructor,Herwig::NBodyDecayConstructorBase> 
    interfaceNBodyDecayConstructors
    ("NBodyDecayConstructors",
     "Vector of references to NBodyDecayConstructors",
     &DecayConstructor::NBodyDecayConstructors_, -1, false, false, true,
     false, false);

  static ParVector<DecayConstructor,string> interfaceDisableModes
    ("DisableModes",
     "A list of decay modes to disable",
     &DecayConstructor::_disableDMTags, -1, string(""), string(""), string(""),
     false, false, Interface::nolimits);

  static Reference<DecayConstructor,DecayRadiationGenerator> interfaceQEDGenerator
    ("QEDGenerator",
     "Object to generate QED radiation in particle decays",
     &DecayConstructor::QEDGenerator_, false, false, true, true, false);

}

/** A helper function for for_each to sort the decay mode tags into the 
 *  standard order.
 */
namespace {

  void adjustFSOrder(string & tag) {
    string::size_type sep = tag.find(">");
    string head = tag.substr(0, sep + 1);
    string products = tag.substr(sep + 1);
    OrderedParticles finalstate;
    bool loopbreak(true);
    while ( loopbreak ) {
      sep = products.find(",");
      string child;
      if( sep != string::npos ) {
	child = products.substr(0, sep);
	products = products.substr(sep + 1);
      }
      else {
	child = string(products.begin(), products.end() - 1);
	loopbreak = false;
      }
      PDPtr p = BaseRepository::GetObject<PDPtr>
	(string("/Herwig/Particles/" + child));
      if( p ) finalstate.insert(p);
    }
    if( finalstate.empty() ) return;
    tag = head;
    OrderedParticles::const_iterator iend = finalstate.end();
    OrderedParticles::size_type count(0), npr(finalstate.size());
    for( OrderedParticles::const_iterator it = finalstate.begin(); 
	 it != iend;  ++it ) {
      tag += (**it).name();
      if( ++count != npr ) tag += string(",");
    }
    tag += string(";");
   }
}

namespace {
  /// Helper function for sorting by number of outgoing lines
  inline bool orderNBodyConstructors(tNBodyDecayConstructorBasePtr a,
				     tNBodyDecayConstructorBasePtr b) {
    return a->numBodies() < b->numBodies();
  }
}

void DecayConstructor::doinit() {
  Interfaced::doinit();
  //Need to check that the stored decay mode tags have the
  //products in the standard order
  for_each( _disableDMTags.begin(), _disableDMTags.end(), adjustFSOrder );
  sort(NBodyDecayConstructors_.begin(), NBodyDecayConstructors_.end(),
       orderNBodyConstructors);
}

void DecayConstructor::createDecayers(const PDVector & particles,
				      double minBR) {
  _minBR = minBR;
  if ( particles.empty() || NBodyDecayConstructors_.empty() ) return;
  // turn the vector into a set to avoid duplicates
  set<PDPtr,NBodyDecayConstructorBase::MassOrdering> particleSet(particles.begin(),particles.end());
  // remove any antiparticles
  for(set<PDPtr>::iterator it=particleSet.begin();it!=particleSet.end();++it) {
    PDPtr cc = (**it).CC();
    if(!cc) continue;
    set<PDPtr>::iterator ic = particleSet.find(cc);
    if(ic!=particleSet.end()) particleSet.erase(ic);
  }
  // set the decay list in the NBodyDecayConstructors
  typedef vector<NBodyDecayConstructorBasePtr>::iterator NBDecayIterator;
  NBDecayIterator it =  NBodyDecayConstructors_.begin();
  NBDecayIterator iend = NBodyDecayConstructors_.end();
  for( ; it != iend; ++it ) {
    (**it).init();
    (**it).decayConstructor(this);
    (**it).DecayList(particleSet);
  }
}

bool DecayConstructor::disableDecayMode(string tag) const {
  if( _disableDMTags.empty() ) return false;
  vector<string>::const_iterator dit = _disableDMTags.begin();
  vector<string>::const_iterator dend = _disableDMTags.end();
  bool disable(false);
  for( ; dit != dend; ++dit ) {
    if( *dit == tag ) {
      disable = true;
      break;
    }
  }
  return disable;
}
#line 1 "./NBodyDecayConstructorBase.cc"
// -*- C++ -*-
//
// NBodyDecayConstructorBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NBodyDecayConstructorBase class.
//

#include "NBodyDecayConstructorBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "DecayConstructor.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig; 
using namespace ThePEG;

void NBodyDecayConstructorBase::persistentOutput(PersistentOStream & os ) const {
  os << init_ << iteration_ << points_ << info_ << decayConstructor_
     << createModes_ << minReleaseFraction_ << maxBoson_ << maxList_
     << removeOnShell_ << excludeEffective_ << includeTopOnShell_
     << excludedVerticesVector_ << excludedVerticesSet_ 
     << excludedParticlesVector_ << excludedParticlesSet_
     << removeFlavourChangingVertices_ << removeSmallVertices_
     << minVertexNorm_;
}

void NBodyDecayConstructorBase::persistentInput(PersistentIStream & is , int) {
  is >> init_ >> iteration_ >> points_ >> info_ >> decayConstructor_
     >> createModes_ >> minReleaseFraction_ >> maxBoson_ >> maxList_
     >> removeOnShell_ >> excludeEffective_ >> includeTopOnShell_
     >> excludedVerticesVector_ >> excludedVerticesSet_
     >> excludedParticlesVector_ >> excludedParticlesSet_
     >> removeFlavourChangingVertices_ >> removeSmallVertices_
     >> minVertexNorm_;
}

// Static variable needed for the type description system in ThePEG.
DescribeAbstractClass<NBodyDecayConstructorBase,Interfaced>
describeThePEGNBodyDecayConstructorBase("Herwig::NBodyDecayConstructorBase", "Herwig.so");

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<NBodyDecayConstructorBase,Interfaced>
describeHerwigNBodyDecayConstructorBase("Herwig::NBodyDecayConstructorBase", "Herwig.so");

void NBodyDecayConstructorBase::Init() {

  static ClassDocumentation<NBodyDecayConstructorBase> documentation
    ("The NBodyDecayConstructorBase class is the base class for the automatic"
     "construction of the decay modes");
  
  static Switch<NBodyDecayConstructorBase,bool> interfaceInitializeDecayers
    ("InitializeDecayers",
     "Initialize new decayers",
     &NBodyDecayConstructorBase::init_, true, false, false);
  static SwitchOption interfaceInitializeDecayersInitializeDecayersYes
    (interfaceInitializeDecayers,
     "Yes",
     "Initialize new decayers to find max weights",
     true);
  static SwitchOption interfaceInitializeDecayersNo
    (interfaceInitializeDecayers,
     "No",
     "Use supplied weights for integration",
     false);
  
  static Parameter<NBodyDecayConstructorBase,int> interfaceInitIteration
    ("InitIteration",
     "Number of iterations to optimise integration weights",
     &NBodyDecayConstructorBase::iteration_, 1, 0, 10,
     false, false, true);

  static Parameter<NBodyDecayConstructorBase,int> interfaceInitPoints
    ("InitPoints",
     "Number of points to generate when optimising integration",
     &NBodyDecayConstructorBase::points_, 1000, 100, 100000000,
     false, false, true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceOutputInfo
    ("OutputInfo",
     "Whether to output information about the decayers",
     &NBodyDecayConstructorBase::info_, false, false, false);
  static SwitchOption interfaceOutputInfoNo
    (interfaceOutputInfo,
     "No",
     "Do not output information regarding the created decayers",
     false);
  static SwitchOption interfaceOutputInfoYes
    (interfaceOutputInfo,
     "Yes",
     "Output information regarding the decayers",
     true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceCreateDecayModes
    ("CreateDecayModes",
     "Whether to create the ThePEG::DecayMode objects as well as the decayers",
     &NBodyDecayConstructorBase::createModes_, true, false, false);
  static SwitchOption interfaceCreateDecayModesYes
    (interfaceCreateDecayModes,
     "Yes",
     "Create the ThePEG::DecayMode objects",
     true);
  static SwitchOption interfaceCreateDecayModesNo
    (interfaceCreateDecayModes,
     "No",
     "Only create the Decayer objects",
     false);

  static Switch<NBodyDecayConstructorBase,unsigned int> interfaceRemoveOnShell
    ("RemoveOnShell",
     "Remove on-shell diagrams as should be treated as a sequence of 1->2 decays",
     &NBodyDecayConstructorBase::removeOnShell_, 1, false, false);
  static SwitchOption interfaceRemoveOnShellYes
    (interfaceRemoveOnShell,
     "Yes",
     "Remove the diagrams if neither the production of decay or the intermediate"
     " can happen",
     1);
  static SwitchOption interfaceRemoveOnShellNo
    (interfaceRemoveOnShell,
     "No",
     "Never remove the intermediate",
     0);
  static SwitchOption interfaceRemoveOnShellProduction
    (interfaceRemoveOnShell,
     "Production",
     "Remove the diagram if the on-shell production of the intermediate is allowed",
     2);

  static RefVector<NBodyDecayConstructorBase,VertexBase> interfaceExcludedVertices
    ("ExcludedVertices",
     "Vertices which are not included in the three-body decayers",
     &NBodyDecayConstructorBase::excludedVerticesVector_,
     -1, false, false, true, true, false);

  static RefVector<NBodyDecayConstructorBase,ParticleData> interfaceExcludedIntermediates
    ("ExcludedIntermediates",
     "Excluded intermediate particles",
     &NBodyDecayConstructorBase::excludedParticlesVector_,
     -1, false, false, true, true, false);

  static Switch<NBodyDecayConstructorBase,bool> interfaceExcludeEffectiveVertices
    ("ExcludeEffectiveVertices",
     "Exclude effectice vertices",
     &NBodyDecayConstructorBase::excludeEffective_, true, false, false);
  static SwitchOption interfaceExcludeEffectiveVerticesYes
    (interfaceExcludeEffectiveVertices,
     "Yes",
     "Exclude the effective vertices",
     true);
  static SwitchOption interfaceExcludeEffectiveVerticesNo
    (interfaceExcludeEffectiveVertices,
     "No",
     "Don't exclude the effective vertices",
     false);
  
  static Parameter<NBodyDecayConstructorBase,double> interfaceMinReleaseFraction
    ("MinReleaseFraction",
     "The minimum energy release for a three-body decay, as a "
     "fraction of the parent mass.",
     &NBodyDecayConstructorBase::minReleaseFraction_, 1e-3, 0.0, 1.0,
     false, false, Interface::limited);

  static Switch<NBodyDecayConstructorBase,unsigned int> interfaceMaximumGaugeBosons
    ("MaximumGaugeBosons",
     "Maximum number of electroweak gauge bosons"
     " to be produced as decay products",
     &NBodyDecayConstructorBase::maxBoson_, 1, false, false);
  static SwitchOption interfaceMaximumGaugeBosonsNone
    (interfaceMaximumGaugeBosons,
     "None",
     "Produce no W/Zs",
     0);
  static SwitchOption interfaceMaximumGaugeBosonsSingle
    (interfaceMaximumGaugeBosons,
     "Single",
     "Produce at most one W/Zs",
     1);
  static SwitchOption interfaceMaximumGaugeBosonsDouble
    (interfaceMaximumGaugeBosons,
     "Double",
     "Produce at most two W/Zs",
     2);
  static SwitchOption interfaceMaximumGaugeBosonsTriple
    (interfaceMaximumGaugeBosons,
     "Triple",
     "Produce at most three W/Zs",
     3);

  static Switch<NBodyDecayConstructorBase,unsigned int> interfaceMaximumNewParticles
    ("MaximumNewParticles",
     "Maximum number of particles from the list of "
     "decaying particles to be allowed as decay products",
     &NBodyDecayConstructorBase::maxList_, 0, false, false);
  static SwitchOption interfaceMaximumNewParticlesNone
    (interfaceMaximumNewParticles,
     "None",
     "No particles from the list",
     0);
  static SwitchOption interfaceMaximumNewParticlesOne
    (interfaceMaximumNewParticles,
     "One",
     "A single particle from the list",
     1);
  static SwitchOption interfaceMaximumNewParticlesTwo
    (interfaceMaximumNewParticles,
     "Two",
     "Two particles from the list",
     2);
  static SwitchOption interfaceMaximumNewParticlesThree
    (interfaceMaximumNewParticles,
     "Three",
     "Three particles from the list",
     3);
  static SwitchOption interfaceMaximumNewParticlesFour
    (interfaceMaximumNewParticles,
     "Four",
     "Four particles from the list",
     4);

  static Switch<NBodyDecayConstructorBase,bool> interfaceIncludeOnShellTop
    ("IncludeOnShellTop",
     "Include the on-shell diagrams involving t -> bW",
     &NBodyDecayConstructorBase::includeTopOnShell_, false, false, false);
  static SwitchOption interfaceIncludeOnShellTopYes
    (interfaceIncludeOnShellTop,
     "Yes",
     "Inlude them",
     true);
  static SwitchOption interfaceIncludeOnShellTopNo
    (interfaceIncludeOnShellTop,
     "No",
     "Don't include them",
     true);

  static Switch<NBodyDecayConstructorBase,bool> interfaceRemoveSmallVertices
    ("RemoveSmallVertices",
     "Remove vertices with norm() below minVertexNorm",
     &NBodyDecayConstructorBase::removeSmallVertices_, false, false, false);
  static SwitchOption interfaceRemoveSmallVerticesYes
    (interfaceRemoveSmallVertices,
     "Yes",
     "Remove them",
     true);
  static SwitchOption interfaceRemoveSmallVerticesNo
    (interfaceRemoveSmallVertices,
     "No",
     "Don't remove them",
     false);

  static Parameter<NBodyDecayConstructorBase,double> interfaceMinVertexNorm
    ("MinVertexNorm",
     "Minimum allowed value of the notm() of the vertex if removing small vertices",
     &NBodyDecayConstructorBase::minVertexNorm_, 1e-8, 1e-300, 1.,
     false, false, Interface::limited);

  static Switch<NBodyDecayConstructorBase,bool> interfaceRemoveFlavourChangingVertices
    ("RemoveFlavourChangingVertices",
     "Remove flavour changing interactions with the photon and gluon",
     &NBodyDecayConstructorBase::removeFlavourChangingVertices_, false, false, false);
  static SwitchOption interfaceRemoveFlavourChangingVerticesYes
    (interfaceRemoveFlavourChangingVertices,
     "Yes",
     "Remove them",
     true);
  static SwitchOption interfaceRemoveFlavourChangingVerticesNo
    (interfaceRemoveFlavourChangingVertices,
     "No",
     "Don't remove them",
     false);

}

void NBodyDecayConstructorBase::setBranchingRatio(tDMPtr dm, Energy pwidth) {
  // if zero width just set BR to zero
  if(pwidth==ZERO) {
    generator()->preinitInterface(dm, "BranchingRatio","set", "0.");
    generator()->preinitInterface(dm, "OnOff","set", "Off");
    return;
  }
  // Need width and branching ratios for all currently created decay modes
  PDPtr parent = const_ptr_cast<PDPtr>(dm->parent());
  DecaySet modes = parent->decayModes();
  unsigned int nmodes=0;
  for( auto dm : modes ) {
    if(dm->on()) ++nmodes;
  }
  if( nmodes==0 ) return;
  double dmbrat(0.);
  if( nmodes == 1 ) {
    parent->width(pwidth);
    if( pwidth > ZERO ) parent->cTau(hbarc/pwidth);
    dmbrat = 1.;
  }
  else {
    Energy currentwidth(parent->width());
    Energy newWidth(currentwidth + pwidth);
    parent->width(newWidth);
    if( newWidth > ZERO ) parent->cTau(hbarc/newWidth);
    //need to reweight current branching fractions if there are any
    double factor = newWidth > ZERO ? double(currentwidth/newWidth) : 0.;
    for(DecaySet::const_iterator dit = modes.begin(); 
	dit != modes.end(); ++dit) {
      if( **dit == *dm || !(**dit).on() ) continue; 
      double newbrat = (**dit).brat()*factor;
      ostringstream brf;
      brf << setprecision(13)<< newbrat;
      generator()->preinitInterface(*dit, "BranchingRatio",
				    "set", brf.str());
    }
    //set brat for current mode
    dmbrat = newWidth > ZERO ? double(pwidth/newWidth) : 0.;
  }
  ostringstream br;
  br << setprecision(13) << dmbrat;
  generator()->preinitInterface(dm, "BranchingRatio",
				"set", br.str());
}

void NBodyDecayConstructorBase::setDecayerInterfaces(string fullname) const {
  if( initialize() ) {
    ostringstream value;
    value << initialize();
    generator()->preinitInterface(fullname, "Initialize", "set",
				  value.str());
    value.str("");
    value << iteration();
    generator()->preinitInterface(fullname, "Iteration", "set",
				  value.str());
    value.str("");
    value << points();
    generator()->preinitInterface(fullname, "Points", "set",
				  value.str());
  }
  // QED stuff if needed
  if(decayConstructor()->QEDGenerator())
    generator()->preinitInterface(fullname, "PhotonGenerator", "set",
				  decayConstructor()->QEDGenerator()->fullName());
  string outputmodes;
  if( info() ) outputmodes = string("Output");
  else outputmodes = string("NoOutput");
  generator()->preinitInterface(fullname, "OutputModes", "set",
				outputmodes);
}

void NBodyDecayConstructorBase::doinit() {
  Interfaced::doinit();
  excludedVerticesSet_ = set<VertexBasePtr>(excludedVerticesVector_.begin(),
					    excludedVerticesVector_.end());
  excludedParticlesSet_ = set<PDPtr>(excludedParticlesVector_.begin(),
				     excludedParticlesVector_.end());
  if(removeOnShell_==0&&numBodies()>2) 
    generator()->log() << "Warning: Including diagrams with on-shell "
		       << "intermediates in " << numBodies() << "-body BSM decays, this"
		       << " can lead to double counting and is not"
		       << " recommended unless you really know what you are doing\n"
		       << "This can be switched off using\n set "
		       << fullName() << ":RemoveOnShell Yes\n"; 
}

namespace {

void constructIdenticalSwaps(unsigned int depth,
			     vector<vector<unsigned int> > identical,
			     vector<unsigned int> channelType,
			     list<vector<unsigned int> > & swaps) {
  if(depth==0) {
    unsigned int size = identical.size();
    for(unsigned ix=0;ix<size;++ix) {
      for(unsigned int iy=2;iy<identical[ix].size();++iy)
	identical.push_back(identical[ix]);
    }
  }
  if(depth+1!=identical.size()) {
    constructIdenticalSwaps(depth+1,identical,channelType,swaps);
  }
  else {
    swaps.push_back(channelType);
  }
  for(unsigned int ix=0;ix<identical[depth].size();++ix) {
    for(unsigned int iy=ix+1;iy<identical[depth].size();++iy) {
      vector<unsigned int> newType=channelType;
      swap(newType[identical[depth][ix]],newType[identical[depth][iy]]);
      // at bottom of chain
      if(depth+1==identical.size()) {
	swaps.push_back(newType);
      }
      else {
	constructIdenticalSwaps(depth+1,identical,newType,swaps);
      }
    }
  }
}

void identicalFromSameDecay(unsigned int & loc, const NBVertex & vertex,
			    vector<vector<unsigned int> > & sameDecay) {
  list<pair<PDPtr,NBVertex> >::const_iterator it = vertex.vertices.begin();
  while(it!=vertex.vertices.end()) {
    if(it->second.incoming) {
      identicalFromSameDecay(loc,it->second,sameDecay);
      ++it;
      continue;
    }
    ++loc;
    long id = it->first->id();
    ++it;
    if(it == vertex.vertices.end()) break;
    if(it->second.incoming) continue;
    if(it->first->id()!=id) continue;
    sameDecay.push_back(vector<unsigned int>());
    sameDecay.back().push_back(loc-1);
    while(it != vertex.vertices.end() 
	  && !it->second.incoming
	  && it->first->id()==id) {
      ++loc;
      ++it;
      sameDecay.back().push_back(loc-1);
    }
  };
}

}

void NBodyDecayConstructorBase::DecayList(const set<PDPtr,MassOrdering> & particles) {
  if( particles.empty() ) return;
  // cast the StandardModel to the Hw++ one to get the vertices
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  unsigned int nv(model->numberOfVertices());
  // loop over the particles and create the modes
  for(set<PDPtr,MassOrdering>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    // get the decaying particle
    tPDPtr parent = *ip;
    if ( Debug::level > 0 )
      Repository::cout() << "Constructing N-body decays for " 
			 << parent->PDGName() << '\n';
    // first create prototype 1->2 decays
    std::stack<PrototypeVertexPtr> prototypes;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model->vertex(iv);
      if(excluded(vertex)) continue;
      PrototypeVertex::createPrototypes(parent, vertex, prototypes,this);
    }
    // now expand them into potential decay modes
    list<vector<PrototypeVertexPtr> > modes;
    while(!prototypes.empty()) {
      // get the first prototype from the stack
      PrototypeVertexPtr proto = prototypes.top();
      prototypes.pop();
      // multiplcity too low
      if(proto->npart<numBodies()) {
	// loop over all vertices and expand
	for(unsigned int iv = 0; iv < nv; ++iv) {
	  VertexBasePtr vertex = model->vertex(iv);
	  if(excluded(vertex) ||
	     proto->npart+vertex->getNpoint()>numBodies()+2) continue;
	  PrototypeVertex::expandPrototypes(proto,vertex,prototypes,
					    excludedParticlesSet_,this);
	}
      }
      // multiplcity too high disgard
      else if(proto->npart>numBodies()) {
	continue;
      }
      // right multiplicity
      else {
	// check it's kinematical allowed for physical masses
	if( proto->incomingMass() < proto->outgoingMass() ) continue;
	// and for constituent masses
	Energy outgoingMass = proto->outgoingConstituentMass();
	if( proto->incomingMass() < proto->outgoingConstituentMass() ) continue;
	// remove processes which aren't kinematically allowed within
	// the release fraction
	if( proto->incomingMass() - outgoingMass <=
	    minReleaseFraction_ * proto->incomingMass() ) continue;
	// check the external particles
	if(!proto->checkExternal()) continue;
	// check if first piece on-shell
	bool onShell = proto->canBeOnShell(removeOnShell(),proto->incomingMass(),true);
	// special treatment for three-body top decays
	if(onShell) {
	  if(includeTopOnShell_ &&
	     abs(proto->incoming->id())==ParticleID::t && proto->npart==3) {
	    unsigned int nprimary=0;
	    bool foundW=false, foundQ=false;
	    for(OrderedVertices::const_iterator it = proto->outgoing.begin();
		it!=proto->outgoing.end();++it) {
	      if(abs(it->first->id())==ParticleID::Wplus)
		foundW = true;
	      if(abs(it->first->id())==ParticleID::b ||
		 abs(it->first->id())==ParticleID::s ||
		 abs(it->first->id())==ParticleID::d)
		foundQ = true;
	      ++nprimary;
	    }
	    if(!(foundW&&foundQ&&nprimary==2)) continue;
	  }
	  else continue;
	}
	// check if should be added to an existing decaymode
	bool added = false;
	for(list<vector<PrototypeVertexPtr> >::iterator it=modes.begin();
	    it!=modes.end();++it) {
	  // is the same ?
 	  if(!(*it)[0]->sameDecay(*proto)) continue;
	  // it is the same
	  added = true;
	  // check if it is a duplicate
	  bool already = false;
	  for(unsigned int iz = 0; iz < (*it).size(); ++iz) {
	    if( *(*it)[iz] == *proto) {
	      already = true;
	      break;
	    }
	  }
	  if(!already) (*it).push_back(proto);
	  break;
	}
	if(!added) modes.push_back(vector<PrototypeVertexPtr>(1,proto));
      }
    }
    // now look at the decay modes
    for(list<vector<PrototypeVertexPtr> >::iterator mit = modes.begin();
	mit!=modes.end();++mit) {
      // count the number of gauge bosons and particles from the list
      unsigned int nlist(0),nbos(0);
      for(OrderedParticles::const_iterator it=(*mit)[0]->outPart.begin();
	  it!=(*mit)[0]->outPart.end();++it) {
	if(abs((**it).id()) == ParticleID::Wplus ||
	   abs((**it).id()) == ParticleID::Z0) ++nbos;
	if(particles.find(*it)!=particles.end()) ++nlist;
	if((**it).CC() && particles.find((**it).CC())!=particles.end()) ++nlist;
      }
      // if too many ignore the mode
      if(nbos > maxBoson_ || nlist > maxList_) continue;
      // translate the prototypes into diagrams
      vector<NBDiagram> newDiagrams;
      double symfac(1.);
      bool possibleOnShell=false;
      for(unsigned int ix=0;ix<(*mit).size();++ix) {
	symfac = 1.;
	possibleOnShell |= (*mit)[ix]->possibleOnShell;
	NBDiagram templateDiagram = NBDiagram((*mit)[ix]);
	// extract the ordering
	vector<int> order(templateDiagram.channelType.size());
	for(unsigned int iz=0;iz<order.size();++iz) {
	  order[templateDiagram.channelType[iz]-1]=iz;
	}
	// find any identical particles
	vector<vector<unsigned int> > identical;
	OrderedParticles::const_iterator it=templateDiagram.outgoing.begin();
	unsigned int iloc=0;
	while(it!=templateDiagram.outgoing.end()) {
	  OrderedParticles::const_iterator lt = templateDiagram.outgoing.lower_bound(*it);
	  OrderedParticles::const_iterator ut = templateDiagram.outgoing.upper_bound(*it);
	  unsigned int nx=0;
	  for(OrderedParticles::const_iterator jt=lt;jt!=ut;++jt) {++nx;}
	  if(nx==1) {
	    ++it;
	    ++iloc;
	  }
	  else {
	    symfac *= factorial(nx);
	    identical.push_back(vector<unsigned int>());
	    for(OrderedParticles::const_iterator jt=lt;jt!=ut;++jt) {
	      identical.back().push_back(order[iloc]);
	      ++iloc;
	    }
	    it = ut;
	  }
	}
	// that's it if there outgoing are unqiue
	if(identical.empty()) {
	  newDiagrams.push_back(templateDiagram);
	  continue;
	}
	// find any identical particles which shouldn't be swapped
	unsigned int loc=0;
	vector<vector<unsigned int> > sameDecay;
	identicalFromSameDecay(loc,templateDiagram,sameDecay);
	// compute the swaps
	list<vector<unsigned int> > swaps;
	constructIdenticalSwaps(0,identical,templateDiagram.channelType,swaps);
	// special if identical from same decay
	if(!sameDecay.empty()) {
	  for(vector<vector<unsigned int> >::const_iterator st=sameDecay.begin();
	      st!=sameDecay.end();++st) {
	    for(list<vector<unsigned int> >::iterator it=swaps.begin();
		it!=swaps.end();++it) {
	      for(unsigned int iy=0;iy<(*st).size();++iy) {
		for(unsigned int iz=iy+1;iz<(*st).size();++iz) {
		  if((*it)[(*st)[iy]]>(*it)[(*st)[iz]])
		    swap((*it)[(*st)[iy]],(*it)[(*st)[iz]]);
		}
	      }
	    }
	  }
	}
	// remove any dupliciates
	for(list<vector<unsigned int> >::iterator it=swaps.begin();
	    it!=swaps.end();++it) {
	  for(list<vector<unsigned int> >::iterator jt=it;
	      jt!=swaps.end();++jt) {
	    if(it==jt) continue;
	    bool different=false;
	    for(unsigned int iz=0;iz<(*it).size();++iz) {
	      if((*it)[iz]!=(*jt)[iz]) {
		different = true;
		break;
	      }
	    }
	    if(!different) {
	      jt = swaps.erase(jt);
	      --jt;
	    }
	  }
	}
	// special for identical decay products
	if(templateDiagram.vertices.begin()->second.incoming) {
	  if(   templateDiagram.vertices.begin() ->first ==
		(++templateDiagram.vertices.begin())->first) {
	    if(*(   (*mit)[ix]->outgoing.begin() ->second) ==
	       *((++(*mit)[ix]->outgoing.begin())->second)) {
	      for(list<vector<unsigned int> >::iterator it=swaps.begin();
		  it!=swaps.end();++it) {
		for(list<vector<unsigned int> >::iterator jt=it;
		    jt!=swaps.end();++jt) {
		  if(it==jt) continue;
		  if((*it)[0]==(*jt)[2]&&(*it)[1]==(*jt)[3]) {
		    jt = swaps.erase(jt);
		    --jt;
		  }
		}
	      }
	    }
	  }
	}
	for(list<vector<unsigned int> >::iterator it=swaps.begin();
	    it!=swaps.end();++it) {
	  newDiagrams.push_back(templateDiagram);
	  newDiagrams.back().channelType = *it;
	}
      }
      // create the decay
      if( Debug::level > 1 ) {
	generator()->log() << "Mode: ";
	generator()->log() << (*mit)[0]->incoming->PDGName() << " -> ";
	for(OrderedParticles::const_iterator it=(*mit)[0]->outPart.begin();
	    it!=(*mit)[0]->outPart.end();++it)
	  generator()->log() << (**it).PDGName() << " ";
	generator()->log() << "\n";
	generator()->log() << "There are " << (*mit).size() << " diagrams\n";
	for(unsigned int iy=0;iy<newDiagrams.size();++iy) {
	  generator()->log() << "Diagram: " << iy << "\n";
	  generator()->log() << newDiagrams[iy] << "\n";
	}
      }
      createDecayMode(newDiagrams,possibleOnShell,symfac);
    }
  }
}

void NBodyDecayConstructorBase::createDecayMode(vector<NBDiagram> &,
						bool, double) {
  assert(false);
  throw Exception() << "In NBodyDecayConstructorBase::createDecayMode() which"
		    << " should have be overridden in an inheriting class "
		    << fullName()
		    << Exception::abortnow; 
}
#line 1 "./TwoBodyDecayConstructor.cc"
// -*- C++ -*-
//
// TwoBodyDecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoBodyDecayConstructor class.
//

#include "TwoBodyDecayConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "DecayConstructor.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractFFSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractSSTVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractRFSVertex.fh"
#include "ThePEG/Helicity/Vertex/AbstractRFVVertex.fh"
#include "VectorCurrentDecayConstructor.h"

using namespace Herwig;
using ThePEG::Helicity::VertexBasePtr;

IBPtr TwoBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr TwoBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void TwoBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << alphaQCD_ << alphaQED_ << oenum(inter_);
}

void TwoBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is  >> alphaQCD_ >> alphaQED_>> ienum(inter_);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoBodyDecayConstructor,NBodyDecayConstructorBase>
describeHerwigTwoBodyDecayConstructor("Herwig::TwoBodyDecayConstructor", "Herwig.so");

void TwoBodyDecayConstructor::Init() {

  static ClassDocumentation<TwoBodyDecayConstructor> documentation
    ("The TwoBodyDecayConstructor implements to creation of 2 body decaymodes "
     "and decayers that do not already exist for the given set of vertices.");
  
  static Reference<TwoBodyDecayConstructor,ShowerAlpha> interfaceShowerAlphaQCD
    ("AlphaQCD",
     "The coupling for QCD corrections",
     &TwoBodyDecayConstructor::alphaQCD_, false, false, true, false, false);
  
  static Reference<TwoBodyDecayConstructor,ShowerAlpha> interfaceShowerAlphaQED
    ("AlphaQED",
     "The coupling for QED corrections",
     &TwoBodyDecayConstructor::alphaQED_, false, false, true, false, false);
  
  static Switch<TwoBodyDecayConstructor,ShowerInteraction> interfaceInteractions
    ("Interactions",
     "which interactions to include for the hard corrections",
     &TwoBodyDecayConstructor::inter_, ShowerInteraction::QCD, false, false);
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
  static SwitchOption interfaceInteractionsQCDQED
    (interfaceInteractions,
     "QEDQCD",
     "Both QCD and QED",
     ShowerInteraction::QEDQCD);
  static SwitchOption interfaceInteractionsDARK
    (interfaceInteractions,
     "DARK",
     "DARK Only",
     ShowerInteraction::DARK);
 
  static SwitchOption interfaceInteractionsALL
   (interfaceInteractions,
    "ALL",
    "ALL",
    ShowerInteraction::ALL);

}

void TwoBodyDecayConstructor::DecayList(const set<PDPtr,MassOrdering> & particles) {
  // special for weak decays
  for(unsigned int ix=0;ix<decayConstructor()->decayConstructors().size();++ix) {
    Ptr<Herwig::VectorCurrentDecayConstructor>::pointer 
      weak = dynamic_ptr_cast<Ptr<Herwig::VectorCurrentDecayConstructor>::pointer >
      (decayConstructor()->decayConstructors()[ix]);
    if(!weak) continue;
    weakMassCut_ = max(weakMassCut_,weak->massCut());
  }
  if( particles.empty() ) return;
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  unsigned int nv(model->numberOfVertices());
  
  for(set<PDPtr,MassOrdering>::const_iterator ip=particles.begin();
      ip!=particles.end();++ip) {
    tPDPtr parent = *ip;
    if ( Debug::level > 0 )
      Repository::cout() << "Constructing 2-body decays for " 
			 << parent->PDGName() << '\n';
    multiset<TwoBodyDecay> decays;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      if(excluded(model->vertex(iv)) || 
	 model->vertex(iv)->getNpoint()>3) continue;
      for(unsigned int il = 0; il < 3; ++il) 
	createModes(parent, model->vertex(iv), il,decays);
    }
    if( !decays.empty() ) createDecayMode(decays);
  }
}

void TwoBodyDecayConstructor::
createModes(tPDPtr inpart, VertexBasePtr vertex,
	    unsigned int list, multiset<TwoBodyDecay> & modes) {
  if( !vertex->isIncoming(inpart) || vertex->getNpoint() != 3 )
    return;
  Energy m1(inpart->mass());
  tPDPtr ccpart = inpart->CC() ? inpart->CC() : inpart;
  long id = ccpart->id();
  tPDVector decaylist = vertex->search(list, ccpart);
  tPDVector::size_type nd = decaylist.size();
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //allowed on-shell decay?
    if( m1 <= pb->mass() + pc->mass() ) continue;
    // double counting with current decays?
    if(inpart->iSpin()==PDT::Spin1 && inpart->iCharge()==0 &&
       pb->id()==-pc->id() && abs(pb->id())<=3 && inpart->mass() <= weakMassCut_ ) {
      continue;
    }
    //vertices are defined with all particles incoming
    modes.insert( TwoBodyDecay(inpart,pb, pc, vertex) );
  }
} 

GeneralTwoBodyDecayerPtr
TwoBodyDecayConstructor::createDecayer(TwoBodyDecay decay,
				       vector<VertexBasePtr> vertices) {
  string name;
  using namespace Helicity::VertexType;
  PDT::Spin in   = decay.parent_->iSpin();
  // PDT::Spin out1 = decay.children_.first ->iSpin();
  PDT::Spin out2 = decay.children_.second->iSpin();
  switch(decay.vertex_->getName()) {
  case FFV :
    if(in == PDT::Spin1Half) {
      name = "FFVDecayer";
      if(out2==PDT::Spin1Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "VFFDecayer";
    }
    break;
  case FFS :
    if(in == PDT::Spin1Half) {
      name = "FFSDecayer";
      if(out2==PDT::Spin1Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "SFFDecayer";
    }
    break;
  case VVS :
    if(in == PDT::Spin1) {
      name = "VVSDecayer";
      if(out2==PDT::Spin1)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "SVVDecayer";
    }
    break;
  case VSS :
    if(in == PDT::Spin1) {
      name = "VSSDecayer";
    }
    else {
      name = "SSVDecayer";
      if(out2==PDT::Spin0)
	swap(decay.children_.first,decay.children_.second);
    }
    break;
  case VVT :
    name = in==PDT::Spin2 ? "TVVDecayer" : "Unknown";
    break;
  case FFT :
    name = in==PDT::Spin2 ? "TFFDecayer" : "Unknown";
    break;
  case SST :
    name = in==PDT::Spin2 ? "TSSDecayer" : "Unknown";
    break;
  case SSS :
    name = "SSSDecayer";
    break;
  case VVV :
    name = "VVVDecayer";
    break;
  case RFS :
    if(in==PDT::Spin1Half) {
      name = "FRSDecayer";
      if(out2==PDT::Spin3Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else if(in==PDT::Spin0) {
      name = "SRFDecayer";
      if(out2==PDT::Spin3Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else {
      name = "Unknown";
    }
    break;
  case RFV :
    if(in==PDT::Spin1Half) {
      name = "FRVDecayer";
      if(out2==PDT::Spin3Half)
	swap(decay.children_.first,decay.children_.second);
    }
    else
      name = "Unknown";
    break;
  default : Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      <<  "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName() << Exception::runerror;
  }
  if(name=="Unknown") 
    Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      <<  "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName() << Exception::runerror;
  ostringstream fullname;
  fullname << "/Herwig/Decays/" << name << "_" << decay.parent_->PDGName() 
	   << "_" << decay.children_.first ->PDGName() 
	   << "_" << decay.children_.second->PDGName();
  string classname = "Herwig::" + name;
  GeneralTwoBodyDecayerPtr decayer;
  decayer = dynamic_ptr_cast<GeneralTwoBodyDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
  if(!decayer) 
    Throw<NBodyDecayConstructorError>() 
      << "Error: Cannot assign " << decay.vertex_->fullName() << " to a decayer. " 
      << "Decay is " << decay.parent_->PDGName() << " -> "
      << decay.children_.first ->PDGName() << " " 
      << decay.children_.second->PDGName() << Exception::runerror;
  // set the strong coupling for radiation
  generator()->preinitInterface(decayer, "AlphaS" , "set", alphaQCD_->fullName());
  // set the EM     coupling for radiation
  generator()->preinitInterface(decayer, "AlphaEM", "set", alphaQED_->fullName());
  // set the type of interactions for the correction
  if(inter_==ShowerInteraction::QCD)
    generator()->preinitInterface(decayer, "Interactions", "set", "QCD");
  else if(inter_==ShowerInteraction::QED)
    generator()->preinitInterface(decayer, "Interactions", "set", "QED");
  else
    generator()->preinitInterface(decayer, "Interactions", "set", "QEDQCD");
  // get the vertices for radiation from the external legs
  map<ShowerInteraction,VertexBasePtr> inRad,fourRad;
  vector<map<ShowerInteraction,VertexBasePtr> > outRad(2);
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    inRad[inter] = radiationVertex(decay.parent_,inter);
    outRad[0][inter] = radiationVertex(decay.children_.first ,inter);
    outRad[1][inter] = radiationVertex(decay.children_.second,inter);
    // get any contributing 4 point vertices
    fourRad[inter]   = radiationVertex(decay.parent_,inter, decay.children_);
  }

  // set info on decay
  decayer->setDecayInfo(decay.parent_,decay.children_,vertices,
  			inRad,outRad,fourRad);
  // initialised the decayer
  setDecayerInterfaces(fullname.str());
  decayer->init();
  return decayer;
}

void TwoBodyDecayConstructor::
createDecayMode(multiset<TwoBodyDecay> & decays) {
  tPDPtr inpart = decays.begin()->parent_;
  for( multiset<TwoBodyDecay>::iterator dit = decays.begin();
       dit != decays.end(); ) {
    TwoBodyDecay mode = *dit;
    // get all the moees with the same in and outgoing particles
    pair<multiset<TwoBodyDecay>::iterator,
	 multiset<TwoBodyDecay>::iterator> range = decays.equal_range(mode);
    // construct the decay mode
    tPDPtr pb((mode).children_.first), pc((mode).children_.second);
    string tag = inpart->name() + "->" + pb->name() + "," + 
      pc->name() + ";";
    // Does it exist already ?
    tDMPtr dm = generator()->findDecayMode(tag);
    // find the vertices
    vector<VertexBasePtr> vertices;
    for ( multiset<TwoBodyDecay>::iterator dit2 = range.first;
	  dit2 != range.second; ++dit2) {
      vertices.push_back(dit2->vertex_);
    }
    dit=range.second;
    // now create DecayMode objects that do not already exist      
    if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
      tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
      if(ndm) {
	inpart->stable(false);
	GeneralTwoBodyDecayerPtr decayer=createDecayer(mode,vertices);
	if(!decayer) continue;
	generator()->preinitInterface(ndm, "Decayer", "set",
				      decayer->fullName());
	generator()->preinitInterface(ndm, "Active", "set", "Yes");
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	setBranchingRatio(ndm, width);
	if(width==ZERO || ndm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
      else
	Throw<NBodyDecayConstructorError>() 
	  << "TwoBodyDecayConstructor::createDecayMode - Needed to create "
	  << "new decaymode but one could not be created for the tag " 
	  << tag << Exception::warning;
    }
    else if( dm ) {
      if(dm->brat()<decayConstructor()->minimumBR()) {
	continue;
      }
      if((dm->decayer()->fullName()).find("Mambo") != string::npos) {
	inpart->stable(false);
	GeneralTwoBodyDecayerPtr decayer=createDecayer(mode,vertices);
	if(!decayer) continue;
	generator()->preinitInterface(dm, "Decayer", "set", 
				      decayer->fullName());
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	if(width/(dm->brat()*inpart->width())<1e-10) {
	  string message = "Herwig calculation of the partial width for the decay mode "
	    + inpart->PDGName() + " -> " + pb->PDGName() + " " + pc->PDGName()
	    + " is zero.\n This will cause problems with the calculation of"
	    + " spin correlations.\n It is probably due to inconsistent parameters"
	    + " and decay modes being passed to Herwig via the SLHA file.\n"
	    + " Zeroing the branching ratio for this mode.";
	  setBranchingRatio(dm,ZERO);
	  generator()->logWarning(NBodyDecayConstructorError(message,Exception::warning));
	}
      }
    }
  }
  // update CC mode if it exists
  if( inpart->CC() ) inpart->CC()->synchronize();
}

VertexBasePtr TwoBodyDecayConstructor::radiationVertex(tPDPtr particle,
						       ShowerInteraction inter,
						       tPDPair children) {
  tHwSMPtr model = dynamic_ptr_cast<tHwSMPtr>(generator()->standardModel());
  map<tPDPtr,VertexBasePtr>::iterator rit = radiationVertices_[inter].find(particle);
  tPDPtr cc = particle->CC() ? particle->CC() : particle;
  if(children==tPDPair() && rit!=radiationVertices_[inter].end()) return rit->second;
  unsigned int nv(model->numberOfVertices());
  long bosonID = inter==ShowerInteraction::QCD ? ParticleID::g : ParticleID::gamma;
  tPDPtr gluon = getParticleData(bosonID);
  // look for radiation vertices for incoming and outgoing particles
  for(unsigned int iv=0;iv<nv;++iv) {
    VertexBasePtr vertex = model->vertex(iv);
    // look for 3 point vertices
    if (children==tPDPair()){
      if( !vertex->isIncoming(particle) ||  vertex->getNpoint() != 3 ||
	  !vertex->isOutgoing(particle) || !vertex->isOutgoing(gluon)) continue;      
      for(unsigned int list=0;list<3;++list) {
	tPDVector decaylist = vertex->search(list, particle);
	for( tPDVector::size_type i = 0; i < decaylist.size(); i += 3 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
	  if( pb->id() == bosonID ) swap(pa, pb);
	  if( pc->id() == bosonID ) swap(pa, pc);
	  if( pb->id() != particle->id()) swap(pb, pc);
	  if( pa->id() != bosonID) continue;
	  if( pb       != particle)      continue;
	  if( pc       != cc)            continue;
	  radiationVertices_[inter][particle] = vertex; 
	  return vertex;
	}
      }
    }
    // look for 4 point vertex including a gluon
    else {           
      if( !vertex->isIncoming(particle)       ||  vertex->getNpoint()!=4              ||
      	  !vertex->isOutgoing(children.first) || !vertex->isOutgoing(children.second) || 
	  !vertex->isOutgoing(gluon)) continue;
      
      for(unsigned int list=0;list<4;++list) {
	tPDVector decaylist = vertex->search(list, particle);
	for( tPDVector::size_type i = 0; i < decaylist.size(); i += 4 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i+1]), pc(decaylist[i+2]), pd(decaylist[i+3]);
	  // order so that a = g, b = parent
	  if( pb->id() == bosonID ) swap(pa, pb);
	  if( pc->id() == bosonID ) swap(pa, pc);
	  if( pd->id() == bosonID ) swap(pa, pd);
	  if( pc->id() == particle->id()) swap(pb, pc);
	  if( pd->id() == particle->id()) swap(pb, pd);
	  if( pa->id() != bosonID)  continue;
	  if( pb->id() != particle->id()) continue;

	  if( !((abs(pd->id()) == abs(children. first->id()) &&
		 abs(pc->id()) == abs(children.second->id())) ||
		(abs(pc->id()) == abs(children. first->id()) &&
		 abs(pd->id()) == abs(children.second->id()))))
	    continue;

	  return vertex;
	}
      }
    }
  }
  return VertexBasePtr();
}


#line 1 "./TwoToTwoProcessConstructor.cc"
// -*- C++ -*-
//
// TwoToTwoProcessConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TwoToTwoProcessConstructor class.
//

#include "TwoToTwoProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include <sstream>

using std::stringstream;

using namespace Herwig;

TwoToTwoProcessConstructor::TwoToTwoProcessConstructor() : 
  Nout_(0), nv_(0), allDiagrams_(true),
  processOption_(0), scaleChoice_(0), scaleFactor_(1.) 
{}

IBPtr TwoToTwoProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr TwoToTwoProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void TwoToTwoProcessConstructor::doinit() {
  HardProcessConstructor::doinit();
  if((processOption_==2 || processOption_==4) &&
     outgoing_.size()!=2)
    throw InitException() 
      << "Exclusive processes require exactly"
      << " two outgoing particles but " << outgoing_.size()
      << "have been inserted in TwoToTwoProcessConstructor::doinit()." 
      << Exception::runerror;
  if(processOption_==4 && incoming_.size()!=2)
    throw InitException() 
      << "Exclusive processes require exactly"
      << " two incoming particles but " << incoming_.size()
      << "have been inserted in TwoToTwoProcessConstructor::doinit()." 
      << Exception::runerror;
  Nout_ = outgoing_.size();
  PDVector::size_type ninc = incoming_.size();
  // exit if nothing to do
  if(Nout_==0||ninc==0) return;
  //create vector of initial-state pairs
  for(PDVector::size_type i = 0; i < ninc; ++i) {
    for(PDVector::size_type j = 0; j < ninc; ++j) {
      tPDPair inc = make_pair(incoming_[i], incoming_[j]);
      
      if( (inc.first->iSpin() > inc.second->iSpin()) ||
	  (inc.first->iSpin() == inc.second->iSpin() &&
	   inc.first->id() < inc.second->id()) )
	swap(inc.first, inc.second);

      if( !HPC_helper::duplicateIncoming(inc,incPairs_) ) {
	incPairs_.push_back(inc);
      }
    }
  }
  // excluded vertices
  excludedVertexSet_ = 
    set<VertexBasePtr>(excludedVertexVector_.begin(),
		       excludedVertexVector_.end());
}


void TwoToTwoProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << vertices_ << incoming_ << outgoing_
     << allDiagrams_ << processOption_
     << scaleChoice_ << scaleFactor_ << excluded_ << excludedExternal_
     << excludedVertexVector_ << excludedVertexSet_;
}

void TwoToTwoProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> vertices_ >> incoming_ >> outgoing_
     >> allDiagrams_ >> processOption_
     >> scaleChoice_ >> scaleFactor_ >> excluded_ >> excludedExternal_
     >> excludedVertexVector_ >> excludedVertexSet_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TwoToTwoProcessConstructor,HardProcessConstructor>
describeHerwigTwoToTwoProcessConstructor("Herwig::TwoToTwoProcessConstructor", "Herwig.so");

void TwoToTwoProcessConstructor::Init() {

  static ClassDocumentation<TwoToTwoProcessConstructor> documentation
    ("TwoToTwoProcessConstructor constructs the possible diagrams for "
     "a process given the external particles");
 
  static RefVector<TwoToTwoProcessConstructor,ThePEG::ParticleData> interfaceIn
    ("Incoming",
     "Pointers to incoming particles",
     &TwoToTwoProcessConstructor::incoming_, -1, false, false, true, false);

  static RefVector<TwoToTwoProcessConstructor,ThePEG::ParticleData> interfaceOut
    ("Outgoing",
     "Pointers to incoming particles",
     &TwoToTwoProcessConstructor::outgoing_, -1, false, false, true, false);
 
  static Switch<TwoToTwoProcessConstructor,bool> interfaceIncludeAllDiagrams
    ("IncludeEW",
     "Switch to decide which diagrams to include in ME calc.",
     &TwoToTwoProcessConstructor::allDiagrams_, true, false, false);
  static SwitchOption interfaceIncludeAllDiagramsNo
    (interfaceIncludeAllDiagrams,
     "No",
     "Only include QCD diagrams",
     false);
  static SwitchOption interfaceIncludeAllDiagramsYes
   (interfaceIncludeAllDiagrams,
     "Yes",
    "Include EW+QCD.",
    true);

  static Switch<TwoToTwoProcessConstructor,unsigned int> interfaceProcesses
    ("Processes",
     "Whether to generate inclusive or exclusive processes",
     &TwoToTwoProcessConstructor::processOption_, 0, false, false);
  static SwitchOption interfaceProcessesSingleParticleInclusive
    (interfaceProcesses,
     "SingleParticleInclusive",
     "Require at least one particle from the list of outgoing particles"
     " in the hard process",
     0);
  static SwitchOption interfaceProcessesTwoParticleInclusive
    (interfaceProcesses,
     "TwoParticleInclusive",
     "Require that both the particles in the hard processes are in the"
     " list of outgoing particles",
     1);
  static SwitchOption interfaceProcessesExclusive
    (interfaceProcesses,
     "Exclusive",
     "Require that both the particles in the hard processes are in the"
     " list of outgoing particles in every hard process",
     2);
  static SwitchOption interfaceProcessesVeryExclusive
    (interfaceProcesses,
     "VeryExclusive",
     "Require that both the incoming and outgoing particles in the hard processes are in the"
     " list of outgoing particles in every hard process",
     4);

  static Switch<TwoToTwoProcessConstructor,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "&TwoToTwoProcessConstructor::scaleChoice_",
     &TwoToTwoProcessConstructor::scaleChoice_, 0, false, false);
  static SwitchOption interfaceScaleChoiceDefault
    (interfaceScaleChoice,
     "Default",
     "Use if sHat if intermediates all colour neutral, otherwise the transverse mass",
     0);
  static SwitchOption interfaceScaleChoicesHat
    (interfaceScaleChoice,
     "sHat",
     "Always use sHat",
     1);
  static SwitchOption interfaceScaleChoiceTransverseMass
    (interfaceScaleChoice,
     "TransverseMass",
     "Always use the transverse mass",
     2);
  static SwitchOption interfaceScaleChoiceGeometicMean
    (interfaceScaleChoice,
     "MaxMT",
     "Use the maximum of m^2+p_T^2 for the two particles",
     3);

  static Parameter<TwoToTwoProcessConstructor,double> interfaceScaleFactor
    ("ScaleFactor",
     "The prefactor used in the scale calculation. The scale used is"
     " that defined by scaleChoice multiplied by this prefactor",
     &TwoToTwoProcessConstructor::scaleFactor_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static RefVector<TwoToTwoProcessConstructor,ThePEG::ParticleData> interfaceExcluded
    ("Excluded",
     "Particles which are not allowed as intermediates",
     &TwoToTwoProcessConstructor::excluded_, -1, false, false, true, false, false);

  static RefVector<TwoToTwoProcessConstructor,ParticleData> interfaceExcludedExternal
    ("ExcludedExternal",
     "Particles which are not allowed as outgoing particles",
     &TwoToTwoProcessConstructor::excludedExternal_, -1,
     false, false, true, false, false);

  static RefVector<TwoToTwoProcessConstructor,VertexBase> interfaceExcludedVertices
    ("ExcludedVertices",
     "Vertices which are not included in the 2 -> 2 scatterings",
     &TwoToTwoProcessConstructor::excludedVertexVector_, -1, false, false, true, true, false);

}

void TwoToTwoProcessConstructor::constructDiagrams() {
  if(incPairs_.empty() || outgoing_.empty() || !subProcess() ) return;
  nv_ = model()->numberOfVertices();
  //make sure  vertices are initialised
  for(unsigned int ix = 0; ix < nv_; ++ix ) {
    VertexBasePtr vertex = model()->vertex(ix);
    if(excludedVertexSet_.find(vertex) != 
       excludedVertexSet_.end()) continue;
    vertices_.push_back(vertex);
  }
  nv_ = vertices_.size();
  //Create necessary diagrams
  vector<tcPDPair>::size_type is;
  PDVector::size_type os;
  for(is = 0; is < incPairs_.size(); ++is) {
    tPDPair ppi = incPairs_[is]; 
    for(os = 0; os < Nout_; ++os) { 
      long fs = outgoing_[os]->id();
      for(size_t iv = 0; iv < nv_; ++iv) {
	tVertexBasePtr vertexA = vertices_[iv];
	//This skips an effective vertex and the EW ones if 
	// we only want the strong diagrams
	if( !allDiagrams_ && vertexA->orderInGs() == 0 ) 
	  continue;

	if(vertexA->getNpoint() == 3) {
	  //scattering diagrams
	  createTChannels(ppi, fs, vertexA);
	  
	  //resonance diagrams
	  if( vertexA->isIncoming(ppi.first) &&  
	      vertexA->isIncoming(ppi.second) )
	    createSChannels(ppi, fs, vertexA);
	}
	else {
	  makeFourPointDiagrams(ppi.first->id(), ppi.second->id(),
				fs, vertexA);
	}
      }
    }
  }
  // need to find all of the diagrams that relate to the same process
  // first insert them into a map which uses the '<' operator 
  // to sort the diagrams 
  multiset<HPDiagram> grouped;
  HPDVector::iterator dit = processes_.begin();
  HPDVector::iterator dend = processes_.end();
  for( ; dit != dend; ++dit) {
    grouped.insert(*dit);
  }
  assert( processes_.size() == grouped.size() );
  processes_.clear();
  typedef multiset<HPDiagram>::const_iterator set_iter;
  set_iter it = grouped.begin(), iend = grouped.end();
  while( it != iend ) {
    pair<set_iter,set_iter> range = grouped.equal_range(*it);
    set_iter itb = range.first;
    HPDVector process;
    for( ; itb != range.second; ++itb ) {
      process.push_back(*itb);
    }
    // if inclusive enforce the exclusivity
    if(processOption_==2 || processOption_==4) {
      if(!((process[0].outgoing. first==outgoing_[0]->id()&&
	    process[0].outgoing.second==outgoing_[1]->id())||
	   (process[0].outgoing. first==outgoing_[1]->id()&&
	    process[0].outgoing.second==outgoing_[0]->id()))) {
	process.clear();
	it = range.second;
	continue;
      }
      if(processOption_==4) {
	if(!((process[0].incoming. first==incoming_[0]->id()&&
	      process[0].incoming.second==incoming_[1]->id())||
	     (process[0].incoming. first==incoming_[1]->id()&&
	      process[0].incoming.second==incoming_[0]->id()))) {
	  process.clear();
	  it = range.second;
	  continue;
	}
      }
    }
    // check no zero width s-channel intermediates
    for( dit=process.begin(); dit != process.end(); ++dit) {
      tPDPtr out1 = getParticleData(dit->outgoing.first );
      tPDPtr out2 = getParticleData(dit->outgoing.second);
      if(dit->channelType == HPDiagram::sChannel && 
	 dit->intermediate->width()==ZERO &&
	 dit->intermediate->mass() > out1->mass()+ out2->mass()) {
	tPDPtr in1 = getParticleData(dit->incoming.first );
	tPDPtr in2 = getParticleData(dit->incoming.second);
	throw Exception() << "Process with zero width resonant intermediates\n"
			  << dit->intermediate->PDGName() 
			  << " can be on-shell in the process "
			  << in1 ->PDGName() << " " <<  in2->PDGName() << " -> "
			  << out1->PDGName() << " " << out2->PDGName() 
			  << " but has zero width.\nEither set the width, enable "
			  << "calculation of its decays, and hence the width,\n"
			  << "or disable it as a potential intermediate using\n"
			  << "insert " << fullName() << ":Excluded 0 "
			  << dit->intermediate->fullName() << "\n---\n"
			  << Exception::runerror;
      }
    }    
    if(find(excludedExternal_.begin(),excludedExternal_.end(),
	    getParticleData(process[0].outgoing. first))!=excludedExternal_.end()) {
      process.clear();
      it = range.second;
      continue;
    }
    if(find(excludedExternal_.begin(),excludedExternal_.end(),
	    getParticleData(process[0].outgoing.second))!=excludedExternal_.end()) {
      process.clear();
      it = range.second;
      continue;
    }
    // finally if the process is allow assign the colour flows
    for(unsigned int ix=0;ix<process.size();++ix) assignToCF(process[ix]);
    // create the matrix element
    createMatrixElement(process);
    process.clear();
    it = range.second;
  }
}

void TwoToTwoProcessConstructor::
createSChannels(tcPDPair inpp, long fs, tVertexBasePtr vertex) {
  //Have 2 incoming particle and a vertex, find the possible offshell
  //particles
  pair<long,long> inc = make_pair(inpp.first->id(), inpp.second->id());
  tPDSet offshells = search(vertex, inpp.first->id(), incoming,
			   inpp.second->id(), incoming, outgoing);
  tPDSet::const_iterator it;
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    if(find(excluded_.begin(),excluded_.end(),*it)!=excluded_.end()) continue;
    for(size_t iv = 0; iv < nv_; ++iv) {
      tVertexBasePtr vertexB = vertices_[iv];
      if( vertexB->getNpoint() != 3) continue;
      if( !allDiagrams_ && vertexB->orderInGs() == 0 ) continue;
      
      tPDSet final;
      if( vertexB->isOutgoing(getParticleData(fs)) &&
	  vertexB->isIncoming(*it) )
	final = search(vertexB, (*it)->id(), incoming, fs,
		       outgoing, outgoing);
      //Now make diagrams
      if(!final.empty()) 
	makeDiagrams(inc, fs, final, *it, HPDiagram::sChannel,
		     make_pair(vertex, vertexB), make_pair(true,true));
    }
  }

}

void TwoToTwoProcessConstructor::
createTChannels(tPDPair inpp, long fs, tVertexBasePtr vertex) {
  pair<long,long> inc = make_pair(inpp.first->id(), inpp.second->id());
  //first try a with c
  tPDSet offshells = search(vertex, inpp.first->id(), incoming, fs,
			   outgoing, outgoing);
  tPDSet::const_iterator it;
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    if(find(excluded_.begin(),excluded_.end(),*it)!=excluded_.end()) continue;
     for(size_t iv = 0; iv < nv_; ++iv) {
       tVertexBasePtr vertexB = vertices_[iv];
       if( vertexB->getNpoint() != 3 ) continue;
       if( !allDiagrams_ && vertexB->orderInGs() == 0 ) continue;
       tPDSet final;
       if( vertexB->isIncoming(inpp.second) )
	 final = search(vertexB, inc.second, incoming, (*it)->id(),
			incoming, outgoing);
       if( !final.empty() )
	 makeDiagrams(inc, fs, final, *it, HPDiagram::tChannel,
		      make_pair(vertex, vertexB), make_pair(true,true));
     }
  }
  //now try b with c
  offshells = search(vertex, inpp.second->id(), incoming, fs,
			   outgoing, incoming);
  for(it = offshells.begin(); it != offshells.end(); ++it) {
    if(find(excluded_.begin(),excluded_.end(),*it)!=excluded_.end()) continue;
    for(size_t iv = 0; iv < nv_; ++iv) {
       tVertexBasePtr vertexB = vertices_[iv];
       if( vertexB->getNpoint() != 3 ) continue;
       if( !allDiagrams_ && vertexB->orderInGs() == 0 ) continue;

       tPDSet final;
       if( vertexB->isIncoming(inpp.first) )
	 final = search(vertexB, inc.first, incoming, (*it)->id(),
			outgoing, outgoing);
       if( !final.empty() )
	 makeDiagrams(inc, fs, final, *it, HPDiagram::tChannel,
		      make_pair(vertexB, vertex), make_pair(true, false));
    }
  }

}

void TwoToTwoProcessConstructor::makeFourPointDiagrams(long parta, long partb,
						   long partc, VBPtr vert) {
  if(processOption_>=1) {
    PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					getParticleData(partc));
    if(loc==outgoing_.end()) return;
  }
  tPDSet ext = search(vert, parta, incoming, partb,incoming, partc, outgoing);
  if( ext.empty() ) return;
  IDPair in(parta, partb);
  for(tPDSet::const_iterator iter=ext.begin(); iter!=ext.end();
      ++iter) {
    if(processOption_>=1) {
      PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					  *iter);
      if(loc==outgoing_.end()) continue;
    }
    HPDiagram nhp(in,make_pair(partc, (*iter)->id()));
    nhp.vertices = make_pair(vert, vert);
    nhp.channelType = HPDiagram::fourPoint;
    fixFSOrder(nhp);
    if(!checkOrder(nhp)) continue;
    if( !duplicate(nhp, processes_) ) processes_.push_back(nhp);
  }
}

void 
TwoToTwoProcessConstructor::makeDiagrams(IDPair in, long out1, const tPDSet & out2, 
				     PDPtr inter, HPDiagram::Channel chan, 
				     VBPair vertexpair, BPair cross) {
  if(processOption_>=1) {
    PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					getParticleData(out1));
    if(loc==outgoing_.end()) return;
  }
  for(tPDSet::const_iterator it = out2.begin(); it != out2.end(); ++it) {
    if(processOption_>=1) {
      PDVector::const_iterator loc = find(outgoing_.begin(),outgoing_.end(),
					  *it);
      if(loc==outgoing_.end()) continue;
    }
    HPDiagram nhp( in,make_pair(out1, (*it)->id()) );
    nhp.intermediate = inter;
    nhp.vertices = vertexpair;
    nhp.channelType = chan;
    nhp.ordered = cross;
    fixFSOrder(nhp);
    if(!checkOrder(nhp)) continue;
    if( !duplicate(nhp, processes_) ) processes_.push_back(nhp);
  }
}

set<tPDPtr> 
TwoToTwoProcessConstructor::search(VBPtr vertex, long part1, direction d1, 
			       long part2, direction d2, direction d3) {
  if(vertex->getNpoint() != 3) return tPDSet();
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  vector<long> ext;
  tPDSet third;
  for(unsigned int ix = 0;ix < 3; ++ix) {
    vector<long> pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0; ix < ext.size(); ix += 3) {
    long id0 = ext.at(ix);
    long id1 = ext.at(ix+1);
    long id2 = ext.at(ix+2);
    int pos;
    if((id0 == part1 && id1 == part2) ||
       (id0 == part2 && id1 == part1))
      pos = ix + 2;
    else if((id0 == part1 && id2 == part2) ||
	    (id0 == part2 && id2 == part1))
      pos = ix + 1;
    else if((id1 == part1 && id2 == part2) ||
	    (id1 == part2 && id2 == part1))
      pos = ix;
    else
      pos = -1;
    if(pos >= 0) {
      tPDPtr p = getParticleData(ext[pos]);
      if(d3 == incoming && p->CC()) p = p->CC();
      third.insert(p);
    }
  }
  
  return third;
}

set<tPDPtr>
TwoToTwoProcessConstructor::search(VBPtr vertex,
				   long part1, direction d1,
				   long part2, direction d2,
				   long part3, direction d3,
				   direction d4) {
  if(vertex->getNpoint() != 4) return tPDSet();
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  if(d3 == incoming && getParticleData(part3)->CC()) part3 = -part3;
  vector<long> ext;
  tPDSet fourth;
  for(unsigned int ix = 0;ix < 4; ++ix) {
    vector<long> pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0;ix < ext.size(); ix += 4) {
    long id0 = ext.at(ix); long id1 = ext.at(ix + 1);
    long id2 = ext.at(ix + 2); long id3 = ext.at(ix + 3);
    int pos;
    if((id0 == part1 && id1 == part2 && id2 == part3) ||
       (id0 == part1 && id1 == part3 && id2 == part2) ||
       (id0 == part2 && id1 == part1 && id2 == part3) ||
       (id0 == part2 && id1 == part3 && id2 == part1) ||
       (id0 == part3 && id1 == part1 && id2 == part2) ||
       (id0 == part3 && id1 == part2 && id2 == part1))
      pos = ix + 3;
    else  if((id0 == part1 && id1 == part2 && id3 == part3) ||
	     (id0 == part1 && id1 == part3 && id3 == part2) ||
	     (id0 == part2 && id1 == part1 && id3 == part3) ||
	     (id0 == part2 && id1 == part3 && id3 == part1) ||
	     (id0 == part3 && id1 == part1 && id3 == part2) ||
	     (id0 == part3 && id1 == part2 && id3 == part1))
      pos = ix + 2;
    else if((id0 == part1 && id2 == part2 && id3 == part3) ||
	    (id0 == part1 && id2 == part3 && id3 == part2) ||
	    (id0 == part2 && id2 == part1 && id3 == part3) ||
	    (id0 == part2 && id2 == part3 && id3 == part1) ||
	    (id0 == part3 && id2 == part1 && id3 == part2) ||
	    (id0 == part3 && id2 == part2 && id3 == part1))
      pos = ix + 1;
    else if((id1 == part1 && id2 == part2 && id3 == part3) ||
	    (id1 == part1 && id2 == part3 && id3 == part2) ||
	    (id1 == part2 && id2 == part1 && id3 == part3) ||
	    (id1 == part2 && id2 == part3 && id3 == part1) ||
	    (id1 == part3 && id2 == part1 && id3 == part2) ||
	    (id1 == part3 && id2 == part2 && id3 == part1))
      pos = ix;
    else 
      pos = -1;
    
    if(pos >= 0) {
      tPDPtr p = getParticleData(ext[pos]);
      if(d4 == incoming && p->CC()) 
	p = p->CC();
      fourth.insert(p);
    }
  } 
  return fourth;
}

void 
TwoToTwoProcessConstructor::createMatrixElement(const HPDVector & process) const {
  if ( process.empty() ) return;
  // external particles
  tcPDVector extpart(4);
  extpart[0] = getParticleData(process[0].incoming.first);
  extpart[1] = getParticleData(process[0].incoming.second);
  extpart[2] = getParticleData(process[0].outgoing.first);
  extpart[3] = getParticleData(process[0].outgoing.second);
  // create the object
  string objectname ("/Herwig/MatrixElements/");
  string classname = MEClassname(extpart, objectname);
  GeneralHardMEPtr matrixElement = dynamic_ptr_cast<GeneralHardMEPtr>
      (generator()->preinitCreate(classname, objectname));
  if( !matrixElement ) {
    std::stringstream message;
    message << "TwoToTwoProcessConstructor::createMatrixElement "
	    << "- No matrix element object could be created for "
	    << "the process " 
	    << extpart[0]->PDGName() << " " << extpart[0]->iSpin() << "," 
	    << extpart[1]->PDGName() << " " << extpart[1]->iSpin() << "->" 
	    << extpart[2]->PDGName() << " " << extpart[2]->iSpin() << "," 
	    << extpart[3]->PDGName() << " " << extpart[3]->iSpin() 
	    << ".  Constructed class name: \"" << classname << "\"";
    generator()->logWarning(TwoToTwoProcessConstructorError(message.str(),Exception::warning));
    return;
  }
  // choice for the scale
  unsigned int scale;
  if(scaleChoice_==0) {
    // check coloured initial and final state
    bool inColour  = ( extpart[0]->coloured() ||
		       extpart[1]->coloured());
    bool outColour = ( extpart[2]->coloured() ||
		       extpart[3]->coloured());
    if(inColour&&outColour) {
      bool coloured = false;
      for(unsigned int ix=0;ix<process.size();++ix) {
	if(process[ix].intermediate&&
	   process[ix].intermediate->coloured()) {
	  coloured = true;
	  break;
	}
      }
      scale = coloured ? 1 : 0;
    }
    else {
      scale = 0;
    } 
  }
  else {
    scale = scaleChoice_-1;
  }
  // set the information
  matrixElement->setProcessInfo(process, colourFlow(extpart),
				debug(), scale, scaleFactor_);
  // insert it
  generator()->preinitInterface(subProcess(), "MatrixElements", 
				subProcess()->MEs().size(),
				"insert", matrixElement->fullName()); 
}

string TwoToTwoProcessConstructor::MEClassname(const vector<tcPDPtr> & extpart, 
					   string & objname) const {
  string classname("Herwig::ME");
  for(vector<tcPDPtr>::size_type ix = 0; ix < extpart.size(); ++ix) {
    if(ix == 2) classname += "2";
    if(extpart[ix]->iSpin() == PDT::Spin0) classname += "s";
    else if(extpart[ix]->iSpin() == PDT::Spin1) classname += "v";
    else if(extpart[ix]->iSpin() == PDT::Spin1Half) classname += "f";
    else if(extpart[ix]->iSpin() == PDT::Spin3Half) classname += "r";
    else if(extpart[ix]->iSpin() == PDT::Spin2) classname += "t";
    else {
      std::stringstream message;
      message << "MEClassname() : Encountered an unknown spin for "
	      << extpart[ix]->PDGName() << " while constructing MatrixElement "
	      << "classname " << extpart[ix]->iSpin();
      generator()->logWarning(TwoToTwoProcessConstructorError(message.str(),Exception::warning));
    }
  }
  objname += "ME" + extpart[0]->PDGName() + extpart[1]->PDGName() + "2" 
    + extpart[2]->PDGName() + extpart[3]->PDGName();
  return classname;  
}
#line 1 "./HardProcessConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HardProcessConstructor class.
//

#include "HardProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void HardProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << debug_ << subProcess_ << model_;
}

void HardProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> debug_ >> subProcess_ >> model_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<HardProcessConstructor,Interfaced>
describeHerwigHardProcessConstructor("Herwig::HardProcessConstructor", "Herwig.so");

void HardProcessConstructor::Init() {

  static ClassDocumentation<HardProcessConstructor> documentation
    ("Base class for implementation of the automatic generation of hard processes");

  static Switch<HardProcessConstructor,bool> interfaceDebugME
    ("DebugME",
     "Print comparison with analytical ME",
     &HardProcessConstructor::debug_, false, false, false);
  static SwitchOption interfaceDebugMEYes
    (interfaceDebugME,
     "Yes",
     "Print the debug information",
     true);
  static SwitchOption interfaceDebugMENo
    (interfaceDebugME,
     "No",
     "Do not print the debug information",
     false);

}

void HardProcessConstructor::doinit() {
  Interfaced::doinit();
  EGPtr eg = generator();
  model_ = dynamic_ptr_cast<HwSMPtr>(eg->standardModel());
  if(!model_)
    throw InitException() << "HardProcessConstructor:: doinit() - "
			  << "The model pointer is null!"
			  << Exception::abortnow;
  if(!eg->eventHandler()) {
    throw
      InitException() << "HardProcessConstructor:: doinit() - "
		      << "The eventHandler pointer was null therefore "
		      << "could not get SubProcessHandler pointer "
		      << Exception::abortnow;
  }
  string subProcessName =
    eg->preinitInterface(eg->eventHandler(), "SubProcessHandlers", "get","");
  subProcess_ = eg->getObject<SubProcessHandler>(subProcessName);
  if(!subProcess_ && HardProcessConstructor::debug_) {
    ostringstream s;
    s << "HardProcessConstructor:: doinit() - "
      << "There was an error getting the SubProcessHandler "
      << "from the current event handler. ";
    generator()->logWarning( Exception(s.str(), Exception::warning) );
  }
}

GeneralHardME::ColourStructure HardProcessConstructor::
colourFlow(const tcPDVector & extpart) const {
  PDT::Colour ina = extpart[0]->iColour();
  PDT::Colour inb = extpart[1]->iColour();
  PDT::Colour outa = extpart[2]->iColour();
  PDT::Colour outb = extpart[3]->iColour();

  // incoming colour neutral
  if(ina == PDT::Colour0 && inb == PDT::Colour0) {
    if( outa == PDT::Colour0 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour11to11;
    }
    else if( outa == PDT::Colour3 && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour11to33bar;
    } 
    else if( outa == PDT::DarkColourFundamental && outb == PDT::DarkColourAntiFundamental) {
      return GeneralHardME::Colour11to33bar;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour11to88;
    }
    else
      assert(false);
  }
  // incoming 3 3
  else if(ina == PDT::Colour3 && inb == PDT::Colour3 ) {
    if( outa == PDT::Colour3 && outb == PDT::Colour3 ) {
      return GeneralHardME::Colour33to33;
    }
    else if( outa == PDT::Colour6 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour33to61;
    }
    else if( outa == PDT::Colour0 && outb == PDT::Colour6 ) {
      return GeneralHardME::Colour33to16;
    }
    else if ( outa == PDT::Colour0 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour33to13bar;
    }
    else if ( outb == PDT::Colour0 && outa == PDT::Colour3bar) {
      return GeneralHardME::Colour33to3bar1;
    }
    else if ( outa == PDT::Colour8 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour33to83bar;
    }
    else if ( outb == PDT::Colour8 && outa == PDT::Colour3bar) {
      return GeneralHardME::Colour33to3bar8;
    }
    else
      assert(false);
  }
  // incoming 3bar 3bar
  else if(ina == PDT::Colour3bar && inb == PDT::Colour3bar ) {
    if( outa == PDT::Colour3bar && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour3bar3barto3bar3bar;
    }
    else if( outa == PDT::Colour6bar && outb == PDT::Colour0) {
      return GeneralHardME::Colour3bar3barto6bar1;
    }
    else if ( outa == PDT::Colour0 && outb == PDT::Colour6bar ) {
      return GeneralHardME::Colour3bar3barto16bar;
    }
    else if ( outa == PDT::Colour0 && outb == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto13;
    }
    else if ( outb == PDT::Colour0 && outa == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto31;
    }
    else if ( outa == PDT::Colour8 && outb == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto83;
    }
    else if ( outb == PDT::Colour8 && outa == PDT::Colour3) {
      return GeneralHardME::Colour3bar3barto38;
    }
    else
      assert(false);
  }
  // incoming 3 3bar
  else if(ina == PDT::Colour3 && inb == PDT::Colour3bar ) {
    if( outa == PDT::Colour0 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour33barto11;
    }
    else if( outa == PDT::Colour3 && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour33barto33bar;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour33barto88;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour33barto81;
    }
    else if( outa == PDT::Colour0 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour33barto18;
    }
    else if( outa == PDT::Colour6 && outb == PDT::Colour6bar) {
      return GeneralHardME::Colour33barto66bar;
    }
    else if( outa == PDT::Colour6bar && outb == PDT::Colour6) {
      return GeneralHardME::Colour33barto6bar6;
    }
    else
      assert(false);
  }
  // incoming 88
  else if(ina == PDT::Colour8 && inb == PDT::Colour8 ) {
    if( outa == PDT::Colour0 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour88to11;
    }
    else if( outa == PDT::Colour3 && outb == PDT::Colour3bar ) {
      return GeneralHardME::Colour88to33bar;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour88to88;
    }
    else if( outa == PDT::Colour8 && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour88to81;
    }
    else if( outa == PDT::Colour0 && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour88to18;
    }
    else if( outa == PDT::Colour6 && outb == PDT::Colour6bar ) {
      return GeneralHardME::Colour88to66bar;
    }
    else
      assert(false);
  }
  // incoming 38
  else if(ina == PDT::Colour3 && inb == PDT::Colour8 ) {
    if(outa == PDT::Colour3 && outb == PDT::Colour0) {
      return GeneralHardME::Colour38to31;
    }
    else if(outa == PDT::Colour0 && outb == PDT::Colour3) {
      return GeneralHardME::Colour38to13;
    }
    else if(outa == PDT::Colour3 && outb == PDT::Colour8) {
      return GeneralHardME::Colour38to38;
    }
    else if(outa == PDT::Colour8 && outb == PDT::Colour3) {
      return GeneralHardME::Colour38to83;
    }
    else if(outa == PDT::Colour3bar && outb == PDT::Colour6){
      return GeneralHardME::Colour38to3bar6;
    }
    else if(outa == PDT::Colour6 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour38to63bar;
    }
    else if(outa == PDT::Colour3bar && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour38to3bar3bar;
    }
    else
      assert(false);
  }
  // incoming 3bar8
  else if(ina == PDT::Colour3bar && inb == PDT::Colour8 ) {
    if(outa == PDT::Colour3bar && outb == PDT::Colour0 ) {
      return GeneralHardME::Colour3bar8to3bar1;
    }
    else if(outa == PDT::Colour0 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour3bar8to13bar;
    }
    else if(outa == PDT::Colour3bar && outb == PDT::Colour8 ) {
      return GeneralHardME::Colour3bar8to3bar8;
    }
    else if(outa == PDT::Colour8 && outb == PDT::Colour3bar) {
      return GeneralHardME::Colour3bar8to83bar;
    }
    else if(outa == PDT::Colour3 && outb == PDT::Colour3) {
      return GeneralHardME::Colour3bar8to33;
    }
    else
      assert(false);
  }
  // unknown colour flow
  else
    assert(false);
  return GeneralHardME::UNDEFINED;
}


void HardProcessConstructor::fixFSOrder(HPDiagram & diag) {
  tcPDPtr psa = getParticleData(diag.incoming.first);
  tcPDPtr psb = getParticleData(diag.incoming.second);
  tcPDPtr psc = getParticleData(diag.outgoing.first);
  tcPDPtr psd = getParticleData(diag.outgoing.second);

  //fix a spin order
  if( psc->iSpin() < psd->iSpin() ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }

  if( psc->iSpin() == psd->iSpin() &&
      psc->id() < 0 && psd->id() > 0 ) {
    swap(diag.outgoing.first, diag.outgoing.second);
    if(diag.channelType == HPDiagram::tChannel) {
      diag.ordered.second = !diag.ordered.second;
    }
    return;
  }
}

void HardProcessConstructor::assignToCF(HPDiagram & diag) {
  if(diag.channelType == HPDiagram::tChannel) {
    if(diag.ordered.second) tChannelCF(diag);
    else                    uChannelCF(diag);
  }
  else if(diag.channelType == HPDiagram::sChannel) {
    sChannelCF(diag);
  }
  else if (diag.channelType == HPDiagram::fourPoint) {
    fourPointCF(diag);
  }
  else
    assert(false);
}

void HardProcessConstructor::tChannelCF(HPDiagram & diag) {
  tcPDPtr ia = getParticleData(diag.incoming.first );
  tcPDPtr ib = getParticleData(diag.incoming.second);
  tcPDPtr oa = getParticleData(diag.outgoing.first );
  tcPDPtr ob = getParticleData(diag.outgoing.second);
  PDT::Colour ina  = ia->iColour();
  PDT::Colour inb  = ib->iColour();
  PDT::Colour outa = oa->iColour();
  PDT::Colour outb = ob->iColour();
  vector<CFPair> cfv(1, make_pair(0, 1.));
  if(diag.intermediate->iColour() == PDT::Colour0) {
    if(ina==PDT::Colour0) {
      cfv[0] = make_pair(0, 1);
    }
    else if(ina==PDT::Colour3 || ina==PDT::Colour3bar) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(2, 1);
      }
    }
    else if(ina==PDT::Colour8) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(7, -1);
      }
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour8) {
    if(ina==PDT::Colour8&&outa==PDT::Colour8&&
       inb==PDT::Colour8&&outb==PDT::Colour8) {
      cfv[0]=make_pair(2,  2.);
      cfv.push_back(make_pair(3, -2.));
      cfv.push_back(make_pair(1, -2.));
      cfv.push_back(make_pair(4,  2.));
    }
    else if(ina==PDT::Colour8&&outa==PDT::Colour0&&
	    inb==PDT::Colour8&&outb==PDT::Colour8&&
	    (oa->iSpin()==PDT::Spin0||oa->iSpin()==PDT::Spin1Half||
	     oa->iSpin()==PDT::Spin3Half)) {
      cfv[0] = make_pair(0,-1);
    }
    else if(ina==PDT::Colour8&&outa==PDT::Colour8&&
	    inb==PDT::Colour8&&outb==PDT::Colour0&&
	    (ob->iSpin()==PDT::Spin0||ob->iSpin()==PDT::Spin1Half||
	     ob->iSpin()==PDT::Spin3Half)) {
      cfv[0] = make_pair(0,-1);
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour3 ||
	  diag.intermediate->iColour() == PDT::Colour3bar) {
    if(outa == PDT::Colour0 || outb == PDT::Colour0) {
      if( outa == PDT::Colour6    || outb == PDT::Colour6   ||
	  outa == PDT::Colour6bar || outb == PDT::Colour6bar) {
	cfv[0] = make_pair(0,0.5);
	cfv.push_back(make_pair(1,0.5));
      }
      else if ((ina==PDT::Colour3 && inb == PDT::Colour3 &&
		(outa == PDT::Colour3bar || outb == PDT::Colour3bar))||
	       (ina==PDT::Colour3bar && inb == PDT::Colour3bar &&
		(outa == PDT::Colour3 || outb == PDT::Colour3 ))) {
	cfv[0] = make_pair(0,1.);
      }
      else {
	cfv[0] = make_pair(0,1.);
      }
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(4,1.);
      cfv.push_back(make_pair(5,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour6bar) {
      cfv[0] = make_pair(4, 1.);
      for(unsigned int ix=5;ix<8;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour6 || outa ==PDT::Colour6bar ||
	    outb==PDT::Colour6 || outb ==PDT::Colour6bar ) {
      assert(false);
    }
    else if(ina==PDT::Colour3    && inb==PDT::Colour3    ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3bar)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3bar))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3bar)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3bar)) {
	cfv[0] = make_pair(1,1.);
      }
    }
    else if(ina==PDT::Colour3bar && inb==PDT::Colour3bar ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3)) {
	double sign = diag.intermediate->iSpin()==PDT::Spin0 ? -1. : 1.;
	cfv[0] = make_pair(1,sign);
      }
    }
    else if((ina==PDT::Colour3    && inb==PDT::Colour8) ||
	    (ina==PDT::Colour3bar && inb==PDT::Colour8) ||
	    (inb==PDT::Colour3    && ina==PDT::Colour8) ||
	    (inb==PDT::Colour3bar && ina==PDT::Colour8) ) {
      if((outa==PDT::Colour3    && outb==PDT::Colour3    ) ||
	 (outa==PDT::Colour3bar && outb==PDT::Colour3bar)) {
	cfv[0] = make_pair(1,1.);
      }
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour6 ||
	  diag.intermediate->iColour() == PDT::Colour6bar) {
    if(ina==PDT::Colour8 && inb==PDT::Colour8) {
      cfv[0] = make_pair(0, 1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
      for(unsigned int ix=4;ix<8;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour3bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(0,1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(4,1.);
      cfv.push_back(make_pair(5,1.));
    }
  }
  diag.colourFlow = cfv;
}

void HardProcessConstructor::uChannelCF(HPDiagram & diag) {
  tcPDPtr ia = getParticleData(diag.incoming.first );
  tcPDPtr ib = getParticleData(diag.incoming.second);
  tcPDPtr oa = getParticleData(diag.outgoing.first );
  tcPDPtr ob = getParticleData(diag.outgoing.second);
  PDT::Colour ina  = ia->iColour();
  PDT::Colour inb  = ib->iColour();
  PDT::Colour outa = oa->iColour();
  PDT::Colour outb = ob->iColour();
  PDT::Colour offshell = diag.intermediate->iColour();
  vector<CFPair> cfv(1, make_pair(1, 1.));
  if(offshell == PDT::Colour8) {
    if(outa == PDT::Colour0 &&
       outb == PDT::Colour0) {
      cfv[0].first = 0;
    }
    else if( outa != outb ) {
      if(outa == PDT::Colour0 ||
	 outb == PDT::Colour0) {
	cfv[0].first = 0;
      }
      else if(ina  == PDT::Colour3 && inb  == PDT::Colour8 &&
	      outb == PDT::Colour3 && outa == PDT::Colour8) {
	tPDPtr off = diag.intermediate;
	if(off->CC()) off=off->CC();
	if(off->iSpin()!=PDT::Spin1Half ||
	   diag.vertices.second->allowed(off->id(),diag.outgoing.first,diag.incoming.second)) {
	  cfv[0].first = 0;
	  cfv.push_back(make_pair(1, -1.));
	}
	else {
	  cfv[0].first = 1;
	  cfv.push_back(make_pair(0, -1.));
	}
      }
      else if(ina  == PDT::Colour3bar && inb  == PDT::Colour8 &&
	      outb == PDT::Colour3bar && outa == PDT::Colour8) {
	tPDPtr off = diag.intermediate;
	if(off->CC()) off=off->CC();
	if(off->iSpin()!=PDT::Spin1Half ||
	   diag.vertices.second->allowed(diag.outgoing.first,off->id(),diag.incoming.second)) {
	  cfv[0].first = 0;
	  cfv.push_back(make_pair(1, -1.));
	}
	else {
	  cfv[0].first = 1;
	  cfv.push_back(make_pair(0, -1.));
	}
      }
      else {
	cfv[0].first = 0;
	cfv.push_back(make_pair(1, -1.));
      }
    }
    else if(outa==PDT::Colour8&&ina==PDT::Colour8) {
      cfv[0]=make_pair(4, 2.);
      cfv.push_back(make_pair(5, -2.));
      cfv.push_back(make_pair(0, -2.));
      cfv.push_back(make_pair(2,  2.));
    }
  }
  else if(offshell == PDT::Colour3 || offshell == PDT::Colour3bar) {
    if( outa == PDT::Colour0 || outb == PDT::Colour0 ) {
      if( outa == PDT::Colour6    || outb == PDT::Colour6   ||
	  outa == PDT::Colour6bar || outb == PDT::Colour6bar) {
	cfv[0] = make_pair(0,0.5);
	cfv.push_back(make_pair(1,0.5));
      }
      else if ((ina==PDT::Colour3 && inb == PDT::Colour3 &&
		(outa == PDT::Colour3bar || outb == PDT::Colour3bar))||
	       (ina==PDT::Colour3bar && inb == PDT::Colour3bar &&
		(outa == PDT::Colour3 || outb == PDT::Colour3 ))) {
	double sign = diag.intermediate->iSpin()==PDT::Spin0 ? -1. : 1.;
	cfv[0] = make_pair(0,sign);
      }
      else {
	cfv[0] = make_pair(0,1.);
      }
    }
    else if(outa==PDT::Colour3bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(4,1.);
      cfv.push_back(make_pair(5,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(0,1.);
      for(int ix=0; ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour6bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(4,1.);
      for(int ix=5; ix<8;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(ina==PDT::Colour0 && inb==PDT::Colour0) {
      cfv[0] = make_pair(0,1.);
    }
    else if(ina==PDT::Colour3    && inb==PDT::Colour3    ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3bar)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3bar))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3bar)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3bar)) {
	double sign = diag.intermediate->iSpin()==PDT::Spin0 ? -1. : 1.;
	cfv[0] = make_pair(2,sign);
      }
    }
    else if(ina==PDT::Colour3bar && inb==PDT::Colour3bar ) {
      if((outa==PDT::Colour0 && outb==PDT::Colour3)||
	 (outb==PDT::Colour0 && outa==PDT::Colour3))
	cfv[0] = make_pair(0,1.);
      else if((outa==PDT::Colour8 && outb==PDT::Colour3)||
	      (outb==PDT::Colour8 && outa==PDT::Colour3)) {
	cfv[0] = make_pair(2,1.);
      }
    }
    else if(((ina==PDT::Colour3    && inb==PDT::Colour8) ||
	     (ina==PDT::Colour3bar && inb==PDT::Colour8) ||
	     (inb==PDT::Colour3    && ina==PDT::Colour8) ||
	     (inb==PDT::Colour3bar && ina==PDT::Colour8)) &&
	    ((outa==PDT::Colour3    && outb==PDT::Colour3    ) ||
	     (outa==PDT::Colour3bar && outb==PDT::Colour3bar))) {
      cfv[0] = make_pair(2, 1.);
    }
    else if(( ina==PDT::Colour3    &&  inb==PDT::Colour3bar &&
	      outa==PDT::Colour3    && outb==PDT::Colour3bar)) {
      cfv[0] = make_pair(2, 1.);
      cfv.push_back(make_pair(3,-1.));
    }
  }
  else if( offshell == PDT::Colour0 ) {
    if(ina==PDT::Colour0) {
      cfv[0] = make_pair(0, 1);
    }
    else if(ina==PDT::Colour3 || ina==PDT::Colour3bar) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || inb==PDT::Colour3bar) {
	cfv[0] = make_pair(3, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(2, 1);
      }
    }
    else if(ina==PDT::Colour8) {
      if( inb == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(inb==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(inb==PDT::Colour8) {
	cfv[0] = make_pair(8, -1);
      }
    }
  }
  else if(diag.intermediate->iColour() == PDT::Colour6 ||
	  diag.intermediate->iColour() == PDT::Colour6bar) {
    if(ina==PDT::Colour8 && inb==PDT::Colour8) {
      cfv[0] = make_pair(0, 1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
      for(unsigned int ix=8;ix<12;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
    else if(outa==PDT::Colour3bar && outb==PDT::Colour6) {
      cfv[0] = make_pair(4, 1.);
      cfv.push_back(make_pair(5,1.));
    }
    else if(outa==PDT::Colour6 && outb==PDT::Colour3bar) {
      cfv[0] = make_pair(0, 1.);
      for(unsigned int ix=1;ix<4;++ix)
	cfv.push_back(make_pair(ix,1.));
    }
  }
  diag.colourFlow = cfv;
}

void HardProcessConstructor::sChannelCF(HPDiagram & diag) {
  tcPDPtr pa = getParticleData(diag.incoming.first);
  tcPDPtr pb = getParticleData(diag.incoming.second);
  PDT::Colour ina = pa->iColour();
  PDT::Colour inb = pb->iColour();
  PDT::Colour offshell = diag.intermediate->iColour();
  tcPDPtr pc = getParticleData(diag.outgoing.first);
  tcPDPtr pd = getParticleData(diag.outgoing.second);
  PDT::Colour outa = pc->iColour();
  PDT::Colour outb = pd->iColour();
  vector<CFPair> cfv(1);
  if(offshell == PDT::Colour8) {
    if(ina  == PDT::Colour0 || inb  == PDT::Colour0 ||
       outa == PDT::Colour0 || outb == PDT::Colour0) {
      cfv[0] = make_pair(0, 1);
    }
    else {
      bool incol   = ina  == PDT::Colour8 && inb  == PDT::Colour8;
      bool outcol  = outa == PDT::Colour8 && outb == PDT::Colour8;
      bool intrip  = ina  == PDT::Colour3 && inb  == PDT::Colour3bar;
      bool outtrip = outa == PDT::Colour3 && outb == PDT::Colour3bar;
      bool outsex  = outa == PDT::Colour6 && outb == PDT::Colour6bar;
      bool outsexb = outa == PDT::Colour6bar && outb == PDT::Colour6;
      if(incol || outcol) {
	// Require an additional minus sign for a scalar/fermion
	// 33bar final state due to the way the vertex rules are defined.
	int prefact(1);
	if( ((pc->iSpin() == PDT::Spin1Half && pd->iSpin() == PDT::Spin1Half) ||
	     (pc->iSpin() == PDT::Spin0     && pd->iSpin() == PDT::Spin0    ) ||
	     (pc->iSpin() == PDT::Spin1     && pd->iSpin() == PDT::Spin1    )) &&
	    (outa        == PDT::Colour3   && outb        == PDT::Colour3bar) )
	  prefact = -1;
	if(incol && outcol) {
	  cfv[0] = make_pair(0, -2.);
	  cfv.push_back(make_pair(1,  2.));
	  cfv.push_back(make_pair(3,  2.));
	  cfv.push_back(make_pair(5,  -2.));
	}
	else if(incol && outsex) {
	  cfv[0].first = 4;
	  cfv[0].second =  prefact;
	  for(unsigned int ix=1;ix<4;++ix)
	    cfv.push_back(make_pair(4+ix, prefact));
	  for(unsigned int ix=0;ix<4;++ix)
	    cfv.push_back(make_pair(8+ix,-prefact));
	}
	else {
	  cfv[0].first = 0;
	  cfv[0].second = -prefact;
	  cfv.push_back(make_pair(1, prefact));
	}
      }
      else if( (  intrip && !outtrip ) ||
	       ( !intrip &&  outtrip ) ) {
	if(!outsex)
	  cfv[0] = make_pair(0, 1);
	else {
	  cfv[0] = make_pair(0, 1.);
	  for(unsigned int ix=0;ix<3;++ix)
	    cfv.push_back(make_pair(ix+1, 1.));
	}
      }
      else if((intrip && outsex) || (intrip && outsexb)) {
	cfv[0] = make_pair(0,1.);
	for(int ix=1; ix<4; ++ix)
	  cfv.push_back(make_pair(ix,1.));
      }
      else
	cfv[0] = make_pair(1, 1);
    }
  }
  else if(offshell == PDT::Colour0) {
    if( ina == PDT::Colour0 ) {
      cfv[0] = make_pair(0, 1);
    }
    else if(ina==PDT::Colour3 || ina==PDT::Colour3bar) {
      if( outa == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(outa==PDT::Colour3 || outa==PDT::Colour3bar) {
	cfv[0] = make_pair(3, 1);
      }
      else if(outa==PDT::Colour8) {
	cfv[0] = make_pair(2, 1);
      }
      else if(outa==PDT::Colour6 || outa==PDT::Colour6bar) {
	cfv[0] = make_pair(8, 1.);
	cfv.push_back(make_pair(9,1.));
      }
      else
	assert(false);
    }
    else if(ina==PDT::Colour8) {
      if( outa == PDT::Colour0 ) {
	cfv[0] = make_pair(0, 1);
      }
      else if(outa==PDT::Colour3 || outb==PDT::Colour3bar) {
	cfv[0] = make_pair(2, 1);
      }
      else if(outa==PDT::Colour8) {
	cfv[0] = make_pair(6, 1);
      }
    }
  }
  else if(offshell == PDT::Colour3 || offshell == PDT::Colour3bar) {
    if(outa == PDT::Colour6    || outa == PDT::Colour6bar ||
       outb == PDT::Colour6bar || outb == PDT::Colour6) {
      cfv[0] = make_pair(6, 1.);
      cfv.push_back(make_pair(7,1.));
    }
    else if((ina  == PDT::Colour3    && inb  == PDT::Colour3) ||
	    (ina  == PDT::Colour3bar && inb  == PDT::Colour3bar)) {
      if((outa == PDT::Colour3    && outb == PDT::Colour3   ) ||
	 (outa == PDT::Colour3bar && outb == PDT::Colour3bar)) {
	cfv[0]      = make_pair(2, 1.);
	cfv.push_back(make_pair(3,-1.));
      }
      else
	cfv[0] = make_pair(0,1.);
    }
    else if(((ina==PDT::Colour3    && inb==PDT::Colour8) ||
	     (ina==PDT::Colour3bar && inb==PDT::Colour8) ||
	     (inb==PDT::Colour3    && ina==PDT::Colour8) ||
	     (inb==PDT::Colour3bar && ina==PDT::Colour8) ) &&
	    ((outa==PDT::Colour3    && outb==PDT::Colour3    ) ||
	     (outa==PDT::Colour3bar && outb==PDT::Colour3bar))) {
      cfv[0] = make_pair(0,1.);
    }
    else {
      if(outa == PDT::Colour0 || outb == PDT::Colour0)
	cfv[0] = make_pair(0, 1);
      else
	cfv[0] = make_pair(1, 1);
    }
  }
  else if( offshell == PDT::Colour6 || offshell == PDT::Colour6bar) {
    if((ina  == PDT::Colour3    && inb  == PDT::Colour3    &&
	outa == PDT::Colour3    && outb == PDT::Colour3   ) ||
       (ina  == PDT::Colour3bar && inb  == PDT::Colour3bar &&
	outa == PDT::Colour3bar && outb == PDT::Colour3bar)) {
      cfv[0]      = make_pair(2,0.5);
      cfv.push_back(make_pair(3,0.5));
    }
    else if((ina  == PDT::Colour3    && inb  == PDT::Colour3    &&
	     ((outa == PDT::Colour6    && outb == PDT::Colour0)||
	      (outb == PDT::Colour6    && outa == PDT::Colour0))) ||
	    (ina  == PDT::Colour3bar && inb  == PDT::Colour3bar &&
	     ((outa == PDT::Colour6bar && outb == PDT::Colour0)||
	      (outb == PDT::Colour6bar && outa == PDT::Colour0)))) {
      cfv[0]      = make_pair(0,0.5);
      cfv.push_back(make_pair(1,0.5));
    }
    else
      assert(false);
  }
  else {
    if(outa == PDT::Colour0 || outb == PDT::Colour0)
      cfv[0] = make_pair(0, 1);
    else
      cfv[0] = make_pair(1, 1);
  }
  diag.colourFlow = cfv;
}

void HardProcessConstructor::fourPointCF(HPDiagram & diag) {
  using namespace ThePEG::Helicity;
  // count the colours
  unsigned int noct(0),ntri(0),nsng(0),nsex(0),nf(0);
  vector<tcPDPtr> particles;
  for(unsigned int ix=0;ix<4;++ix) {
    particles.push_back(getParticleData(diag.ids[ix]));
    PDT::Colour col = particles.back()->iColour();
    if(col==PDT::Colour0)                            ++nsng;
    else if(col==PDT::Colour3||col==PDT::Colour3bar) ++ntri;
    else if(col==PDT::Colour8)                       ++noct;
    else if(col==PDT::Colour6||col==PDT::Colour6bar) ++nsex;
    if(particles.back()->iSpin()==2) nf+=1;
  }
  if(nsng==4 || (ntri==2&&nsng==2) ||
     (noct==3            && nsng==1) ||
     (ntri==2 && noct==1 && nsng==1) ||
     (noct == 2 && nsng == 2) ) {
    vector<CFPair> cfv(1,make_pair(0,1));
    diag.colourFlow = cfv;
  }
  else if(noct==4) {
    // flows for SSVV, VVVV is handled in me class
    vector<CFPair> cfv(6);
    cfv[0] = make_pair(0, -2.);
    cfv[1] = make_pair(1, -2.);
    cfv[2] = make_pair(2, +4.);
    cfv[3] = make_pair(3, -2.);
    cfv[4] = make_pair(4, +4.);
    cfv[5] = make_pair(5, -2.);
    diag.colourFlow = cfv;
  }
  else if(ntri==2&&noct==2) {
    vector<CFPair> cfv(2);
    cfv[0] = make_pair(0, 1);
    cfv[1] = make_pair(1, 1);
    if(nf==2) cfv[1].second = -1.;
    diag.colourFlow = cfv;
  }
  else if(nsex==2&&noct==2) {
    vector<CFPair> cfv;
    for(unsigned int ix=0;ix<4;++ix)
      cfv.push_back(make_pair(ix  ,2.));
    for(unsigned int ix=0;ix<8;++ix)
      cfv.push_back(make_pair(4+ix,1.));
    diag.colourFlow = cfv;
  }
  else if(ntri==4) {
    // get the order from the vertex
    vector<long> temp;
    for(unsigned int ix=0;ix<4;++ix) {
      temp = diag.vertices.first->search(ix,diag.outgoing.first);
      if(!temp.empty()) break;
    }
    // compute the mapping
    vector<long> ids;
    ids.push_back( particles[0]->CC() ? -diag.incoming.first  : diag.incoming.first );
    ids.push_back( particles[1]->CC() ? -diag.incoming.second : diag.incoming.second);
    ids.push_back(  diag.outgoing.first );
    ids.push_back(  diag.outgoing.second);
    vector<unsigned int> order = {0,1,2,3};
    vector<bool> matched(4,false);
    for(unsigned int ix=0;ix<temp.size();++ix) {
      for(unsigned int iy=0;iy<ids.size();++iy) {
	if(matched[iy]) continue;
	if(temp[ix]==ids[iy]) {
	  matched[iy] = true;
	  order[ix]=iy;
	  break;
	}
      }
    }
    // 3 3 -> 3 3
    if((particles[0]->iColour()==PDT::Colour3 &&
	particles[1]->iColour()==PDT::Colour3) ||
       (particles[0]->iColour()==PDT::Colour3bar &&
	particles[1]->iColour()==PDT::Colour3bar) ) {
      if(diag.vertices.first->colourStructure()==ColourStructure::SU3I12I34) {
	if( (order[0]==0 && order[1]==2) || (order[2]==0 && order[3]==2) ||
	    (order[0]==2 && order[1]==0) || (order[2]==2 && order[3]==0))
      	  diag.colourFlow = vector<CFPair>(1,make_pair(2,1.));
      	else
      	  diag.colourFlow = vector<CFPair>(1,make_pair(3,1.));
      }
      else if(diag.vertices.first->colourStructure()==ColourStructure::SU3I14I23) {
      	if( (order[0]==0 && order[3]==2) || (order[1]==0 && order[2]==2) ||
      	    (order[0]==2 && order[3]==0) || (order[1]==2 && order[2]==0))
      	  diag.colourFlow = vector<CFPair>(1,make_pair(2,1.));
      	else
      	  diag.colourFlow = vector<CFPair>(1,make_pair(3,1.));
      }
      else if(diag.vertices.first->colourStructure()==ColourStructure::SU3T21T43) {
	if( (order[1]==0 && order[0]==2) || (order[3]==0 && order[2]==2) ||
	    (order[1]==2 && order[0]==0) || (order[3]==2 && order[2]==0))
	  diag.colourFlow = vector<CFPair>(1,make_pair(0,1.));
	else
	  diag.colourFlow = vector<CFPair>(1,make_pair(1,1.));
      }
      else if(diag.vertices.first->colourStructure()==ColourStructure::SU3T23T41) {
	if( (order[1]==0 && order[2]==2) || (order[3]==0 && order[0]==2) ||
	    (order[1]==2 && order[2]==0) || (order[3]==2 && order[0]==0))
	  diag.colourFlow = vector<CFPair>(1,make_pair(0,1.));
	else
	  diag.colourFlow = vector<CFPair>(1,make_pair(1,1.));
      }
      else
	assert(false);
    }
    else if((particles[0]->iColour()==PDT::Colour3 &&
	     particles[1]->iColour()==PDT::Colour3bar) ||
	    (particles[0]->iColour()==PDT::Colour3bar &&
	     particles[1]->iColour()==PDT::Colour3)) {
      if(diag.vertices.first->colourStructure()==ColourStructure::SU3I12I34) {
      	if( (order[0]==0 && order[1]==1) || (order[2]==0 && order[3]==0) ||
      	    (order[0]==1 && order[1]==0) || (order[2]==1 && order[3]==1))
      	  diag.colourFlow = vector<CFPair>(1,make_pair(3,1.));
      	else
      	  diag.colourFlow = vector<CFPair>(1,make_pair(2,1.));
      }
      else if(diag.vertices.first->colourStructure()==ColourStructure::SU3I14I23) {
      	if( (order[0]==0 && order[3]==1) || (order[0]==2 && order[3]==3) ||
      	    (order[0]==1 && order[3]==0) || (order[0]==3 && order[3]==2))
      	  diag.colourFlow = vector<CFPair>(1,make_pair(3,1.));
      	else
      	  diag.colourFlow = vector<CFPair>(1,make_pair(2,1.));
      }
      else if(diag.vertices.first->colourStructure()==ColourStructure::SU3T21T43) {
       	if( (order[1]==0 && order[0]==1) || (order[3]==0 && order[2]==1) ||
      	    (order[1]==1 && order[0]==0) || (order[3]==1 && order[2]==0))
       	  diag.colourFlow = vector<CFPair>(1,make_pair(1,1.));
      	else
      	  diag.colourFlow = vector<CFPair>(1,make_pair(0,1.));
      }
      else if(diag.vertices.first->colourStructure()==ColourStructure::SU3T23T41) {
	if( (order[1]==0 && order[2]==1) || (order[1]==1 && order[2]==0) ||
	    (order[1]==3 && order[2]==2) || (order[1]==2 && order[2]==3))
	  diag.colourFlow = vector<CFPair>(1,make_pair(1,1.));
	else
	  diag.colourFlow = vector<CFPair>(1,make_pair(2,1.));
      }
      else
	assert(false);
    }
    else {
      assert(false);
    }
  }
  else {
    assert(false);
  }
}

namespace {
  // Helper functor for find_if in duplicate function.
  class SameDiagramAs {
  public:
    SameDiagramAs(const HPDiagram & diag) : a(diag) {}
    bool operator()(const HPDiagram & b) const {
      return a == b;
    }
  private:
    HPDiagram a;
  };
}

bool HardProcessConstructor::duplicate(const HPDiagram & diag,
				       const HPDVector & group) const {
  //find if a duplicate diagram exists
  HPDVector::const_iterator it =
    find_if(group.begin(), group.end(), SameDiagramAs(diag));
  return it != group.end();
}

bool HardProcessConstructor::checkOrder(const HPDiagram & diag) const {
  for(map<string,pair<unsigned int,int> >::const_iterator it=model_->couplings().begin();
      it!=model_->couplings().end();++it) {
    int order=0;
    if(diag.vertices.first ) order += diag.vertices.first ->orderInCoupling(it->second.first);
    if(diag.vertices.second&&diag.vertices.first->getNpoint()==3)
      order += diag.vertices.second->orderInCoupling(it->second.first);
    if(order>it->second.second) return false;
  }
  return true;
}
#line 1 "./HiggsVectorBosonProcessConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsVectorBosonProcessConstructor class.
//

#include "HiggsVectorBosonProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/General/GeneralfftoVH.h"

using namespace Herwig;

HiggsVectorBosonProcessConstructor::HiggsVectorBosonProcessConstructor()
  : _type(true), _shapeOpt(1) {
}

IBPtr HiggsVectorBosonProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr HiggsVectorBosonProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void HiggsVectorBosonProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << _vector << _higgs << _type << _shapeOpt << _alpha;
}

void HiggsVectorBosonProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _vector >> _higgs >> _type >> _shapeOpt >> _alpha;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HiggsVectorBosonProcessConstructor,HardProcessConstructor>
describeHerwigHiggsVectorBosonProcessConstructor("Herwig::HiggsVectorBosonProcessConstructor", "Herwig.so");

void HiggsVectorBosonProcessConstructor::Init() {

  static ClassDocumentation<HiggsVectorBosonProcessConstructor> documentation
    ("The HiggsVectorBosonProcessConstructor class generates hard process for"
     " Higgs boson production in assoication with a vector boson in general models.");

  static RefVector<HiggsVectorBosonProcessConstructor,ParticleData> interfaceVectorBoson
    ("VectorBoson",
     "The possible outgoing vector bosons, must be W/Z",
     &HiggsVectorBosonProcessConstructor::_vector, -1, false, false, true, false, false);

  static RefVector<HiggsVectorBosonProcessConstructor,ParticleData> interfaceHiggsBoson
    ("HiggsBoson",
     "The possible Higgs bosons",
     &HiggsVectorBosonProcessConstructor::_higgs, -1, false, false, true, false, false);

  static Switch<HiggsVectorBosonProcessConstructor,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &HiggsVectorBosonProcessConstructor::_shapeOpt, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeYes
    (interfaceShapeOption,
     "OnShell",
     "Produce the Higgs on-shell",
     0);

  static Switch<HiggsVectorBosonProcessConstructor,bool> interfaceCollisionType
    ("CollisionType",
     "Type of collision",
     &HiggsVectorBosonProcessConstructor::_type, true, false, false);
  static SwitchOption interfaceCollisionTypeLepton
    (interfaceCollisionType,
     "Lepton",
     "Lepton-Lepton collisions",
     false);
  static SwitchOption interfaceCollisionTypeHadron
    (interfaceCollisionType,
     "Hadron",
     "Hadron-Hadron collisions",
     true);

  static Reference<HiggsVectorBosonProcessConstructor,ShowerAlpha> interfaceAlphaQCD
    ("AlphaQCD",
     "The strong coupling used in the shower for MME or POWHEG corrections.",
     &HiggsVectorBosonProcessConstructor::_alpha, false, false, true, false, false);

}

void HiggsVectorBosonProcessConstructor::constructDiagrams() {
  if(_vector.empty()||_higgs.empty() || !subProcess() ) return;
  // initialise the particles
  for(unsigned int ix=0;ix<_vector.size();++ix)
    _vector[ix]->init();
  for(unsigned int ix=0;ix<_higgs.size();++ix)
    _higgs[ix]->init();
  for(PDVector::const_iterator iv=_vector.begin();
      iv!=_vector.end();++iv) {
    // skip if combination not possible
    if(_type==false && (**iv).id()!=ParticleID::Z0)
      continue;
    else if(_type==true && (abs((**iv).id()) != ParticleID::Wplus &&
			    (**iv).id()      != ParticleID::Z0))
      continue;
    // loop over the possible Higgs bosons
    for(PDVector::const_iterator ih=_higgs.begin();
	ih!=_higgs.end();++ih) {
      // check higgs is neutral and scalar
      if((**ih).iCharge()!=0 || (**ih).coloured() ||
	 (**ih).iSpin()!=PDT::Spin0) continue;
      // find a suitable vertex
      for(unsigned int nv = 0; nv < model()->numberOfVertices(); ++nv ) {
	VertexBasePtr vertex = model()->vertex(nv);
	AbstractVVSVertexPtr svert =
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(vertex);
	if(!svert) continue;
	if(vertex->getNpoint() != 3) continue;
	if(!vertex->isIncoming(*iv)) continue;
	if(!vertex->isOutgoing(*iv)) continue;
	if(!vertex->isOutgoing(*ih)) continue;
	// create the MatrixElement object
	string objectname ("/Herwig/MatrixElements/");
	string classname("Herwig::GeneralfftoVH");
	if(_type) objectname += "MEPP2";
	else      objectname += "MEee2";
	objectname += (**iv).PDGName();
	objectname += (**ih).PDGName();
	GeneralfftoVHPtr matrixElement = dynamic_ptr_cast<GeneralfftoVHPtr>
	  (generator()->preinitCreate(classname, objectname));
  if(abs((**iv).id()) == ParticleID::Z0 && !vertex->allowed(23,23,(**ih).id()))
    continue;
  if(abs((**iv).id()) == ParticleID::Wplus && !vertex->allowed(-24,24,(**ih).id()))
    continue;
	if( !matrixElement )
	  throw Exception()
	    << "HiggsVectorBosonProcessConstructor::constructDiagrams() "
	    << " Failed to construct matrix element for "
	    << (**iv).PDGName() << " + "
	    << (**ih).PDGName() << " production"
	    << Exception::runerror;
	GeneralfftoVH::Process process = GeneralfftoVH::Lepton;
	if(_type) {
	  if((**iv).id()==ParticleID::Z0)
	    process = GeneralfftoVH::HadronZ;
	  else if((**iv).id()==ParticleID::Wplus)
	    process = GeneralfftoVH::HadronWplus;
	  else if((**iv).id()==ParticleID::Wminus)
	    process = GeneralfftoVH::HadronWminus;
	}
	// set the coupling
	generator()->preinitInterface(matrixElement, "Coupling",
				      "set", _alpha->fullName());
	// set the information
	matrixElement->setProcessInfo( process, *ih, svert,_shapeOpt);
	// insert it
	generator()->preinitInterface(subProcess(), "MatrixElements",
				      subProcess()->MEs().size(),
				      "insert", matrixElement->fullName());
      }
    }
  }
}
#line 1 "./HiggsVBFProcessConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HiggsVBFProcessConstructor class.
//

#include "HiggsVBFProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/General/GeneralfftoffH.h"

using namespace Herwig;

HiggsVBFProcessConstructor::HiggsVBFProcessConstructor()
  : _type(true), _shapeOpt(1), _intermediates(0) {
}

IBPtr HiggsVBFProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr HiggsVBFProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void HiggsVBFProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << _higgs << _type << _shapeOpt;
}

void HiggsVBFProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _higgs >> _type >> _shapeOpt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HiggsVBFProcessConstructor,HardProcessConstructor>
describeHerwigHiggsVBFProcessConstructor("Herwig::HiggsVBFProcessConstructor", "Herwig.so");

void HiggsVBFProcessConstructor::Init() {

  static ClassDocumentation<HiggsVBFProcessConstructor> documentation
    ("The HiggsVBFProcessConstructor class generates hard processes for"
     " Higgs boson production in association with a vector boson in general models.");

  static RefVector<HiggsVBFProcessConstructor,ParticleData> interfaceHiggsBoson
    ("HiggsBoson",
     "The possible Higgs bosons",
     &HiggsVBFProcessConstructor::_higgs, -1, false, false, true, false, false);

  static Switch<HiggsVBFProcessConstructor,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &HiggsVBFProcessConstructor::_shapeOpt, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeYes
    (interfaceShapeOption,
     "OnShell",
     "Produce the Higgs on-shell",
     0);

  static Switch<HiggsVBFProcessConstructor,unsigned int> interfaceProcess
    ("Process",
     "Which processes to include",
     &HiggsVBFProcessConstructor::_intermediates, 0, false, false);
  static SwitchOption interfaceProcessBoth
    (interfaceProcess,
     "Both",
     "Include both WW and ZZ processes",
     0);
  static SwitchOption interfaceProcessWW
    (interfaceProcess,
     "WW",
     "Only include WW processes",
     1);
  static SwitchOption interfaceProcessZZ
    (interfaceProcess,
     "ZZ",
     "Only include ZZ processes",
     2);

  static Switch<HiggsVBFProcessConstructor,bool> interfaceCollisionType
    ("CollisionType",
     "Type of collision",
     &HiggsVBFProcessConstructor::_type, true, false, false);
  static SwitchOption interfaceCollisionTypeLepton
    (interfaceCollisionType,
     "Lepton",
     "Lepton-Lepton collisions",
     false);
  static SwitchOption interfaceCollisionTypeHadron
    (interfaceCollisionType,
     "Hadron",
     "Hadron-Hadron collisions",
     true);

}

void HiggsVBFProcessConstructor::constructDiagrams() {
  if(_higgs.empty() || !subProcess() ) return;
  tPDPtr Wplus  = getParticleData(ParticleID::Wplus);
  tPDPtr Wminus = getParticleData(ParticleID::Wminus);
  tPDPtr Z0     = getParticleData(ParticleID::Z0);
  for(unsigned int ix=0;ix<_higgs.size();++ix)
    _higgs[ix]->init();
  for(unsigned int ix=0;ix<2;++ix) {
    if( ( ix == 0 && _intermediates == 2 ) ||
	( ix == 1 && _intermediates == 1 )) continue;
    // loop over the possible Higgs bosons
    for(PDVector::const_iterator ih=_higgs.begin();
	ih!=_higgs.end();++ih) {
      // check higgs is neutral and scalar
      if((**ih).iCharge()!=0 || (**ih).coloured() ||
	 (**ih).iSpin()!=PDT::Spin0) continue;
      // find a suitable vertex
      for(unsigned int nv = 0; nv < model()->numberOfVertices(); ++nv ) {
	VertexBasePtr vertex = model()->vertex(nv);
	AbstractVVSVertexPtr svert =
	  dynamic_ptr_cast<AbstractVVSVertexPtr>(vertex);
	if(!svert) continue;
	if(vertex->getNpoint() != 3) continue;
	// check outgoing higgs allowed
	if(!vertex->isOutgoing(*ih)) continue;
	// check incoming W+W- or ZZ allowed
	if(ix==0) {
	  if(!vertex->isIncoming(Wminus)||
	     !vertex->isIncoming(Wplus)) continue;
    if(!vertex->allowed(-24,24,(**ih).id())) continue;
	}
	else {
	  if(!vertex->isIncoming(Z0)) continue;
	  if(!vertex->allowed(23,23,(**ih).id())) continue;
	}
 	// create the MatrixElement object
 	string objectname ("/Herwig/MatrixElements/");
 	string classname("Herwig::GeneralfftoffH");
 	if(_type) objectname += "MEPP2";
 	else      objectname += "MEee2";
	string bos = ix==0 ? "W+W+" : "ZOZO";
	objectname += bos;
 	objectname += (**ih).PDGName();
	GeneralfftoffHPtr matrixElement = dynamic_ptr_cast<GeneralfftoffHPtr>
	  (generator()->preinitCreate(classname, objectname));
	if( !matrixElement )
	  throw Exception()
	    << "HiggsVBFProcessConstructor::constructDiagrams() "
	    << " Failed to construct matrix element for "
	    << bos  << " + "
	    << (**ih).PDGName() << " production"
	    << Exception::runerror;
	GeneralfftoffH::Process process = _type ?
	  GeneralfftoffH::Hadron : GeneralfftoffH::Lepton;
	// set the information
	matrixElement->setProcessInfo( process, *ih, svert,_shapeOpt,
				       ix+1 );
	// insert it
	generator()->preinitInterface(subProcess(), "MatrixElements",
				      subProcess()->MEs().size(),
				      "insert", matrixElement->fullName());
      }
    }
  }
}
#line 1 "./QQHiggsProcessConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QQHiggsProcessConstructor class.
//

#include "QQHiggsProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/General/GeneralQQHiggs.h"

using namespace Herwig;

QQHiggsProcessConstructor::QQHiggsProcessConstructor() 
  : _process(0), _quarkFlavour(0), _shapeOpt(1)
{}

IBPtr QQHiggsProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr QQHiggsProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void QQHiggsProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << _process << _quarkFlavour << _higgs << _shapeOpt;
}

void QQHiggsProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> _process >> _quarkFlavour >> _higgs >> _shapeOpt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<QQHiggsProcessConstructor,HardProcessConstructor>
describeHerwigQQHiggsProcessConstructor("Herwig::QQHiggsProcessConstructor", "Herwig.so");

void QQHiggsProcessConstructor::Init() {

  static ClassDocumentation<QQHiggsProcessConstructor> documentation
    ("The QQHiggsProcessConstructor class generates hard processes for the"
     " production of the Higgs boson in association with a heavy quark-antiquark"
     " pair in general models.");

  static RefVector<QQHiggsProcessConstructor,ParticleData> interfaceHiggsBoson
    ("HiggsBoson",
     "The possible Higgs bosons",
     &QQHiggsProcessConstructor::_higgs, -1, false, false, true, false, false);

  static Switch<QQHiggsProcessConstructor,unsigned int> interfaceShapeOption
    ("ShapeScheme",
     "Option for the treatment of the Higgs resonance shape",
     &QQHiggsProcessConstructor::_shapeOpt, 2, false, false);
  static SwitchOption interfaceStandardShapeFixed
    (interfaceShapeOption,
     "FixedBreitWigner",
     "Breit-Wigner s-channel resonanse",
     1);
  static SwitchOption interfaceStandardShapeRunning
    (interfaceShapeOption,
     "MassGenerator",
     "Use the mass generator to give the shape",
     2);
  static SwitchOption interfaceStandardShapeYes
    (interfaceShapeOption,
     "OnShell",
     "Produce the Higgs on-shell",
     0);

  static Switch<QQHiggsProcessConstructor,unsigned int> interfaceProcess
    ("Process",
     "Which subprocesses to include",
     &QQHiggsProcessConstructor::_process, 0, false, false);
  static SwitchOption interfaceProcessAll
    (interfaceProcess,
     "All",
     "Include all subprocesses",
     0);
  static SwitchOption interfaceProcess1
    (interfaceProcess,
     "gg",
     "Include only gg -> QQbarH processes",
     1);
  static SwitchOption interfaceProcessqbarqbarqbarqbar
    (interfaceProcess,
     "qqbar",
     "Include only qbar qbar -> QQbarH processes",
     2);

  static Switch<QQHiggsProcessConstructor,unsigned int> interfaceQuarkType
    ("QuarkType",
     "The type of quark",
     &QQHiggsProcessConstructor::_quarkFlavour, 6, false, false);
  static SwitchOption interfaceQuarkTypeBottom
    (interfaceQuarkType,
     "Bottom",
     "Produce bottom-antibottom",
     5);
  static SwitchOption interfaceQuarkTypeTop
    (interfaceQuarkType,
     "Top",
     "Produce top-antitop",
     6);
  static SwitchOption interfaceQuarkTypeBottomTop
    (interfaceQuarkType,
     "BottomandTop",
     "Produce bottom-antibottom and top-antitop",
     0);

}


void QQHiggsProcessConstructor::constructDiagrams() {
  if(_higgs.empty() || !subProcess() ) return;
  // initialize the Higgs bosons
  for(unsigned int ix=0;ix<_higgs.size();++ix)
    _higgs[ix]->init();
  long qmin = _quarkFlavour == 0 ? 5 : _quarkFlavour;
  long qmax = _quarkFlavour == 0 ? 6 : _quarkFlavour;
  for(long iq=qmin;iq<=qmax;++iq) {
    tPDPtr qk = getParticleData(iq);
    tPDPtr qb = qk->CC();
    // loop over the possible Higgs bosons
    for(PDVector::const_iterator ih=_higgs.begin();
	ih!=_higgs.end();++ih) {
      // check higgs is neutral and scalar
      if((**ih).iCharge()!=0 || (**ih).coloured() ||
	 (**ih).iSpin()!=PDT::Spin0) continue;
      // find a suitable vertex
      for(unsigned int nv = 0; nv < model()->numberOfVertices(); ++nv ) {
	VertexBasePtr vertex = model()->vertex(nv);
 	AbstractFFSVertexPtr svert = 
 	  dynamic_ptr_cast<AbstractFFSVertexPtr>(vertex);
 	if(!svert) continue;
	if(vertex->getNpoint() != 3) continue;
	// check q qb allowed
	if(!vertex->isOutgoing(qk)||
	   !vertex->isOutgoing(qb)) continue;
	// check outgoing higgs allowed
	if(!vertex->isOutgoing(*ih)) continue;
  	// create the MatrixElement object
  	string objectname ("/Herwig/MatrixElements/");
  	string classname("Herwig::GeneralQQHiggs");
  	objectname += "MEPP2";
	if(iq==5) objectname += "bbbar";
	else      objectname += "ttbar";
 	objectname += (**ih).PDGName();
	GeneralQQHiggsPtr matrixElement = dynamic_ptr_cast<GeneralQQHiggsPtr>
	  (generator()->preinitCreate(classname, objectname));
	if( !matrixElement )
	  throw Exception()
	    << "QQHiggsProcessConstructor::constructDiagrams() "
	    << " Failed to construct matrix element for "
	    << qk->PDGName() << " + " << qb->PDGName() << " + " 
	    << (**ih).PDGName() << " production"
	    << Exception::runerror;
	// set the information
	matrixElement->setProcessInfo( iq, *ih, svert,_shapeOpt,
				       _process );
	// insert it
	generator()->preinitInterface(subProcess(), "MatrixElements", 
				      subProcess()->MEs().size(),
				      "insert", matrixElement->fullName());
      }
    }
  }
}
#line 1 "./ThreeBodyDecayConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ThreeBodyDecayConstructor class.
//

#include "ThreeBodyDecayConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/General/GeneralThreeBodyDecayer.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/Throw.h"
#include "DecayConstructor.h"
#include "WeakCurrentDecayConstructor.h"

using namespace Herwig;

IBPtr ThreeBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ThreeBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void ThreeBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << interOpt_ << widthOpt_ << intOpt_ << relErr_
     << includeIntermediatePhotons_;
}

void ThreeBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is  >> interOpt_ >> widthOpt_ >> intOpt_ >> relErr_
      >> includeIntermediatePhotons_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ThreeBodyDecayConstructor,NBodyDecayConstructorBase>
describeHerwigThreeBodyDecayConstructor("Herwig::ThreeBodyDecayConstructor", "Herwig.so");

void ThreeBodyDecayConstructor::Init() {

  static ClassDocumentation<ThreeBodyDecayConstructor> documentation
    ("The ThreeBodyDecayConstructor class constructs the three body decay modes");

  static Switch<ThreeBodyDecayConstructor,bool> interfaceIncludeIntermediatePhotons
    ("IncludeIntermediatePhotons",
     "Whether or not ot allow intermediate photons",
     &ThreeBodyDecayConstructor::includeIntermediatePhotons_, false, false, false);
  static SwitchOption interfaceIncludeIntermediatePhotonsYes
    (interfaceIncludeIntermediatePhotons,
     "Yes",
     "Include them",
     true);
  static SwitchOption interfaceIncludeIntermediatePhotonsNo
    (interfaceIncludeIntermediatePhotons,
     "No",
     "Don't include them",
     false);

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &ThreeBodyDecayConstructor::widthOpt_, 1, false, false);
  static SwitchOption interfaceWidthOptionFixed
    (interfaceWidthOption,
     "Fixed",
     "Use fixed widths",
     1);
  static SwitchOption interfaceWidthOptionRunning
    (interfaceWidthOption,
     "Running",
     "Use running widths",
     2);
  static SwitchOption interfaceWidthOptionZero
    (interfaceWidthOption,
     "Zero",
     "Set the widths to zero",
     3);

  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceIntermediateOption
    ("IntermediateOption",
     "Option for the inclusion of intermediates in the event",
     &ThreeBodyDecayConstructor::interOpt_, 0, false, false);
  static SwitchOption interfaceIntermediateOptionAlways
    (interfaceIntermediateOption,
     "Always",
     "Always include the intermediates",
     1);
  static SwitchOption interfaceIntermediateOptionNever
    (interfaceIntermediateOption,
     "Never",
     "Never include the intermediates",
     2);
  static SwitchOption interfaceIntermediateOptionOnlyIfOnShell
    (interfaceIntermediateOption,
     "OnlyIfOnShell",
     "Only if there are on-shell diagrams",
     0);


  static Switch<ThreeBodyDecayConstructor,unsigned int> interfaceIntegrationOption
    ("IntegrationOption",
     "Option for the integration of the partial width",
     &ThreeBodyDecayConstructor::intOpt_, 1, false, false);
  static SwitchOption interfaceIntegrationOptionAllPoles
    (interfaceIntegrationOption,
     "AllPoles",
     "Include all potential poles",
     0);
  static SwitchOption interfaceIntegrationOptionShallowestPole
    (interfaceIntegrationOption,
     "ShallowestPole",
     "Only include the  shallowest pole",
     1);

  static Parameter<ThreeBodyDecayConstructor,double> interfaceRelativeError
    ("RelativeError",
     "The relative error for the GQ integration",
     &ThreeBodyDecayConstructor::relErr_, 1e-2, 1e-10, 1.,
     false, false, Interface::limited);

}

void ThreeBodyDecayConstructor::DecayList(const set<PDPtr,MassOrdering> & particles) {
  if( particles.empty() ) return;
  // special for weak decays
  for(unsigned int ix=0;ix<decayConstructor()->decayConstructors().size();++ix) {
    Ptr<Herwig::WeakCurrentDecayConstructor>::pointer 
      weak = dynamic_ptr_cast<Ptr<Herwig::WeakCurrentDecayConstructor>::pointer >
      (decayConstructor()->decayConstructors()[ix]);
    if(!weak) continue;
    weakMassCut_ = max(weakMassCut_,weak->massCut());
  }
  NBodyDecayConstructorBase::DecayList(particles);
}

GeneralThreeBodyDecayerPtr ThreeBodyDecayConstructor::
createDecayer(vector<TBDiagram> & diagrams, bool inter,
	      double symfac) const {
  if(diagrams.empty()) return GeneralThreeBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing           ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // get the object name
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, outgoing, objectname);
  if(classname=="") return GeneralThreeBodyDecayerPtr();
  // create the object
  GeneralThreeBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralThreeBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer and return if doesn't work
  vector<PDPtr> outVector(outgoing.begin(),outgoing.end());
  if(!decayer->setDecayInfo(incoming,outVector,diagrams,symfac))
    return GeneralThreeBodyDecayerPtr();
  // set decayer options from base class
  setDecayerInterfaces(objectname);
  // options for partial width integration
  ostringstream value;
  value << intOpt_;
  generator()->preinitInterface(objectname, "PartialWidthIntegration", "set",
				value.str());
  value.str("");
  value << relErr_;
  generator()->preinitInterface(objectname, "RelativeError", "set",
				value.str());
  // set the width option
  value.str("");
  value << widthOpt_;
  generator()->preinitInterface(objectname, "WidthOption", "set", value.str());
  // set the intermediates option
  value.str("");
  value << inter;
  generator()->preinitInterface(objectname, "GenerateIntermediates", "set", 
				value.str());
  // initialize the decayer
  decayer->init();
  // return the decayer
  return decayer;
}

string ThreeBodyDecayConstructor::
DecayerClassName(tcPDPtr incoming, const OrderedParticles & outgoing, 
		 string & objname) const {
  string classname("Herwig::");
  // spins of the outgoing particles
  unsigned int ns(0),nf(0),nv(0);
  objname += incoming->PDGName() + "2";
  for(OrderedParticles::const_iterator it=outgoing.begin();
      it!=outgoing.end();++it) {
    if     ((**it).iSpin()==PDT::Spin0    ) ++ns;
    else if((**it).iSpin()==PDT::Spin1Half) ++nf;
    else if((**it).iSpin()==PDT::Spin1    ) ++nv;
    objname += (**it).PDGName();
  }
  objname   += "Decayer";
  if(incoming->iSpin()==PDT::Spin0) {
    if(ns==1&&nf==2) classname += "StoSFFDecayer";
    else if(nf==2&&nv==1) classname += "StoFFVDecayer";
    else             classname  = "";
  }
  else if(incoming->iSpin()==PDT::Spin1Half) {
    if(nf==3) classname += "FtoFFFDecayer";
    else if(nf==1&&nv==2) classname += "FtoFVVDecayer";
    else      classname  = "";
  }
  else if(incoming->iSpin()==PDT::Spin1) {
    if(nf==2&&nv==1) classname += "VtoFFVDecayer";
    else classname = "";
  }
  else {
    classname="";
  }
  return classname;
}

void ThreeBodyDecayConstructor::
createDecayMode(vector<NBDiagram> & mode,
		bool possibleOnShell,
		double symfac) {
  // convert the diagrams from the general to the three body structure
  vector<TBDiagram> diagrams;
  for(unsigned int iy=0;iy<mode.size();++iy) {
    diagrams.push_back(TBDiagram(mode[iy]));
    // determine the type
    diagrams.back().channelType = diagrams.back().intermediate ?
      TBDiagram::Channel(mode[iy].channelType[0]-1) : TBDiagram::fourPoint;
    // remove weak processes simulated using the weak current
    if(weakMassCut_>ZERO && diagrams.back().intermediate &&
       abs(diagrams.back().intermediate->id())==ParticleID::Wplus) {
      Energy deltaM = 
    	getParticleData(diagrams.back().incoming)->mass() - 
    	getParticleData(diagrams.back().outgoing)->mass();
      if(deltaM<weakMassCut_) diagrams.pop_back();
    }
    // remove intermediate photons
    else if(!includeIntermediatePhotons_ && 
	    diagrams.back().intermediate &&
    	    abs(diagrams.back().intermediate->id())==ParticleID::gamma)
      diagrams.pop_back();
  }
  if(diagrams.empty()) return;
  // check if possible on-shell internal particles
  bool inter = interOpt_ == 1 || (interOpt_ == 0 && possibleOnShell);
  // incoming particle
  tPDPtr inpart = getParticleData(diagrams[0].incoming);
  // outgoing particles
  OrderedParticles outgoing;
  outgoing.insert(getParticleData(diagrams[0].outgoing));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.first ));
  outgoing.insert(getParticleData(diagrams[0].outgoingPair.second));
  // sort out ordering and labeling of diagrams
  vector<PDPtr> outVector(outgoing.begin(),outgoing.end());
  for(unsigned int ix=0;ix<diagrams.size();++ix) {
    if( ( diagrams[ix].channelType == TBDiagram::channel23 && 
	  outVector[1]->id() != diagrams[ix].outgoingPair.first) ||
	( diagrams[ix].channelType == TBDiagram::channel13 && 
	  outVector[0]->id() != diagrams[ix].outgoingPair.first) || 
	( diagrams[ix].channelType == TBDiagram::channel12 && 
	  outVector[0]->id() != diagrams[ix].outgoingPair.first) ) 
      swap(diagrams[ix].outgoingPair.first, diagrams[ix].outgoingPair.second);
  }
  // construct the tag for the decay mode
  string tag = inpart->name() + "->";
  unsigned int iprod=0;
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end(); ++it) {
    ++iprod;
    tag += (**it).name();
    if(iprod != 3) tag += ",";
  }
  tag += ";";
  tDMPtr dm = generator()->findDecayMode(tag);
  // create mode if needed
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    // create the decayer
    GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    OrderedParticles::const_iterator pit=outgoing.begin();
    tPDPtr pa = *pit; ++pit;
    tPDPtr pb = *pit; ++pit;
    tPDPtr pc = *pit;
    Energy width = 
      decayer->partialWidth(make_pair(inpart,inpart->mass()),
			    make_pair(pa,pa->mass()) , 
			    make_pair(pb,pb->mass()) , 
			    make_pair(pc,pc->mass()));
    if ( Debug::level > 1 )
      generator()->log() << "Partial width is: " << width / GeV << " GeV\n";
    if(width==ZERO) {
      if ( Debug::level > 1 )
	generator()->log() << "Partial width for " 
			   << tag << " zero so mode not created \n";
      generator()->preinitRemove(decayer);
      return;
    }
    tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
    if(!ndm)
      throw NBodyDecayConstructorError() 
	<< "ThreeBodyDecayConstructor::createDecayMode - Needed to create "
	<< "new decaymode but one could not be created for the tag " 
	<< tag << Exception::warning;
    generator()->preinitInterface(ndm, "Decayer", "set",
				  decayer->fullName());
    generator()->preinitInterface(ndm, "Active", "set", "Yes");
    if(!ndm->decayer()) {
      generator()->log() << "Can't set the decayer for " 
			 << tag << " so mode not created \n";
      return;
    }
    setBranchingRatio(ndm, width);
    if(ndm->brat()<decayConstructor()->minimumBR()) {
      generator()->preinitInterface(decayer->fullName(),
				    "Initialize", "set","0");
    }
    // incoming particle is now unstable
    inpart->stable(false);
    if(Debug::level > 1 ) {
      generator()->log() << "Calculated partial width for mode "
			 << tag << " is " << width/GeV << "\n";
    }
  }
  else if( dm ) {
    if(dm->brat()<decayConstructor()->minimumBR()) {
      return;
    }
    if((dm->decayer()->fullName()).find("Mambo") != string::npos) {
      // create the decayer
      GeneralThreeBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
      if(decayer) {
	generator()->preinitInterface(dm, "Decayer", "set", 
				      decayer->fullName());
	// check not zero
	OrderedParticles::const_iterator pit=outgoing.begin();
	tPDPtr pa = *pit; ++pit;
	tPDPtr pb = *pit; ++pit;
	tPDPtr pc = *pit;
	Energy width = 
	  decayer->partialWidth(make_pair(inpart,inpart->mass()),
				make_pair(pa,pa->mass()) , 
				make_pair(pb,pb->mass()) , 
				make_pair(pc,pc->mass()));
	if(width/(dm->brat()*inpart->width())<1e-10) {
	  string message = "Herwig calculation of the partial width for the decay mode "
	    + inpart->PDGName() + " -> " + pa->PDGName() + " " + pb->PDGName() + " " + pc->PDGName()
	    + " is zero.\n This will cause problems with the calculation of"
	    + " spin correlations.\n It is probably due to inconsistent parameters"
	    + " and decay modes being passed to Herwig via the SLHA file.\n"
	    + " Zeroing the branching ratio for this mode.";
	  setBranchingRatio(dm,ZERO);
	  generator()->logWarning(NBodyDecayConstructorError(message,Exception::warning));
	}
      }
      // incoming particle is now unstable
      inpart->stable(false);
    }
  }
  //update CC mode if it exists
  if( inpart->CC() ) inpart->CC()->synchronize();
}
#line 1 "./FourBodyDecayConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourBodyDecayConstructor class.
//

#include "FourBodyDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/General/GeneralFourBodyDecayer.h"
#include "DecayConstructor.h"

using namespace Herwig;

IBPtr FourBodyDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr FourBodyDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void FourBodyDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << interOpt_ << widthOpt_ << particles_;
}

void FourBodyDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> interOpt_ >> widthOpt_ >> particles_;
}

DescribeClass<FourBodyDecayConstructor,NBodyDecayConstructorBase>
describeFourBodyDecayConstructor("Herwig::FourBodyDecayConstructor","Herwig.so");

void FourBodyDecayConstructor::Init() {

  static ClassDocumentation<FourBodyDecayConstructor> documentation
    ("The FourBodyDecayConstructor class implements a small number"
     " of 4-body decays in general models");

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceWidthOption
    ("WidthOption",
     "Option for the treatment of the widths of the intermediates",
     &FourBodyDecayConstructor::widthOpt_, 1, false, false);
  static SwitchOption interfaceWidthOptionFixed
    (interfaceWidthOption,
     "Fixed",
     "Use fixed widths",
     1);
  static SwitchOption interfaceWidthOptionRunning
    (interfaceWidthOption,
     "Running",
     "Use running widths",
     2);
  static SwitchOption interfaceWidthOptionZero
    (interfaceWidthOption,
     "Zero",
     "Set the widths to zero",
     3);

  static Switch<FourBodyDecayConstructor,unsigned int> interfaceIntermediateOption
    ("IntermediateOption",
     "Option for the inclusion of intermediates in the event",
     &FourBodyDecayConstructor::interOpt_, 0, false, false);
  static SwitchOption interfaceIntermediateOptionAlways
    (interfaceIntermediateOption,
     "Always",
     "Always include the intermediates",
     1);
  static SwitchOption interfaceIntermediateOptionNever
    (interfaceIntermediateOption,
     "Never",
     "Never include the intermediates",
     2);
  static SwitchOption interfaceIntermediateOptionOnlyIfOnShell
    (interfaceIntermediateOption,
     "OnlyIfOnShell",
     "Only if there are on-shell diagrams",
     0);

  static RefVector<FourBodyDecayConstructor,ParticleData> interfaceParticles
    ("Particles",
     "Particles to override the choice in the DecayConstructor for 4-body decays,"
     " if empty the defaults from the DecayConstructor are used.",
     &FourBodyDecayConstructor::particles_, -1, false, false, true, true, false);

  static Switch<FourBodyDecayConstructor,bool> interfaceParticleType
    ("ParticleType",
     "Which types of particles to calculate four body decay modes for",
     &FourBodyDecayConstructor::particleType_, false, false, false);
  static SwitchOption interfaceParticleTypeStable
    (interfaceParticleType,
     "Stable",
     "Only calculate four-body decays in no 2/3 body modes",
     false);
  static SwitchOption interfaceParticleTypeAll
    (interfaceParticleType,
     "All",
     "Calculate 4-body modes for all particles",
     true);

}

void FourBodyDecayConstructor::DecayList(const set<PDPtr,MassOrdering> & particles) {
  if( particles.empty() ) return;
  set<PDPtr,MassOrdering> new_particles;
  for(set<PDPtr,MassOrdering>::const_iterator it=particles.begin();it!=particles.end();++it) {
    if(!particles_.empty() && find(particles_.begin(),particles_.end(),*it)==particles_.end()) continue;
    if(!(**it).stable()&&!particleType_) continue;
    new_particles.insert(*it);
  }
  if(!new_particles.empty())
    NBodyDecayConstructorBase::DecayList(new_particles);
}

void FourBodyDecayConstructor::
createDecayMode(vector<NBDiagram> & diagrams,
		bool possibleOnShell, double symfac) {
  // some basic checks for the modes we are interested in
  // only looking at scalars
  if(diagrams[0].incoming->iSpin()!=PDT::Spin0) return;
  // which decay to 4 fermions
  unsigned int nferm=0;
  for(OrderedParticles::const_iterator it=diagrams[0].outgoing.begin();
      it!=diagrams[0].outgoing.end();++it) {
    if((**it).iSpin()==PDT::Spin1Half) ++nferm;
  }
  if(nferm!=4) return;
  // check for on-shell intermediates
  bool inter = interOpt_ == 1 || (interOpt_ == 0 && possibleOnShell);
  // incoming particle  
  tPDPtr inpart = diagrams[0].incoming;
  // outgoing particles
  OrderedParticles outgoing=diagrams[0].outgoing;
  // incoming particle is now unstable
  inpart->stable(false);
  // construct the tag for the decay mode
  string tag = inpart->name() + "->";
  for(OrderedParticles::const_iterator it = outgoing.begin();
      it != outgoing.end(); ++it) {
    if(it!=outgoing.begin()) tag += ",";
    tag += (**it).name();
  }
  tag += ";";
  tDMPtr dm = generator()->findDecayMode(tag);
  // create mode if needed
  if( createDecayModes() && (!dm || inpart->id() == ParticleID::h0) ) {
    // create the decayer
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    // create the decay mode
    tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
    if(ndm) {
      string test = generator()->preinitInterface(ndm, "Decayer", "set",
						  decayer->fullName());
      generator()->preinitInterface(ndm, "Active", "set", "Yes");
      Energy width = 
	decayer->partialWidth(inpart,outgoing);
      setBranchingRatio(ndm, width);
    }
    else 
      throw NBodyDecayConstructorError() 
	<< "FourBodyDecayConstructor::createDecayMode - Needed to create "
	<< "new decaymode but one could not be created for the tag " 
	<< tag << Exception::warning;
  }
  // otherwise 
  else if (dm && (dm->decayer()->fullName()).find("Mambo") != string::npos) {
    // create the decayer
    GeneralFourBodyDecayerPtr decayer = createDecayer(diagrams,inter,symfac);
    if(!decayer) {
      if(Debug::level > 1 ) generator()->log() << "Can't create the decayer for " 
					       << tag << " so mode not created\n";
      return;
    }
    generator()->preinitInterface(dm, "Decayer", "set", 
				  decayer->fullName());
  }
  //update CC mode if it exists
  if( inpart->CC() )
    inpart->CC()->synchronize();
}

GeneralFourBodyDecayerPtr 
FourBodyDecayConstructor::createDecayer(vector<NBDiagram> & diagrams, 
					bool inter, double symfac) const {
  if(diagrams.empty()) return GeneralFourBodyDecayerPtr();
  // extract the external particles for the process
  PDPtr incoming = diagrams[0].incoming;
  // outgoing particles
  vector<PDPtr> outgoing(diagrams[0].outgoing.begin(),
			 diagrams[0].outgoing.end());
  // get the name for the object
  string objectname ("/Herwig/Decays/");
  string classname = DecayerClassName(incoming, diagrams[0].outgoing, objectname);
  if(classname=="") return GeneralFourBodyDecayerPtr();
  // create the object
  GeneralFourBodyDecayerPtr decayer = 
    dynamic_ptr_cast<GeneralFourBodyDecayerPtr>
    (generator()->preinitCreate(classname, objectname));
  // set up the decayer and return if doesn't work
  if(!decayer->setDecayInfo(incoming,outgoing,diagrams,symfac))
    return GeneralFourBodyDecayerPtr();
  // set decayer options from base class
  setDecayerInterfaces(objectname);
  // set the width option
  ostringstream value;
  value << widthOpt_;
  generator()->preinitInterface(objectname, "WidthOption", "set", value.str());
  // set the intermediates option
  ostringstream value2;
  value2 << inter;
  generator()->preinitInterface(objectname, "GenerateIntermediates", "set", 
 				value2.str());
  // initialize the decayer
  decayer->init();
  // return the decayer
  return decayer;
}

string  FourBodyDecayConstructor::DecayerClassName(tcPDPtr incoming,
						   const OrderedParticles & outgoing, 
						   string & objname) const {
  string classname("Herwig::");
  // spins of the outgoing particles
  unsigned int ns(0),nf(0),nv(0);
  objname += incoming->PDGName() + "2";
  for(OrderedParticles::const_iterator it=outgoing.begin();
      it!=outgoing.end();++it) {
    if     ((**it).iSpin()==PDT::Spin0    ) ++ns;
    else if((**it).iSpin()==PDT::Spin1Half) ++nf;
    else if((**it).iSpin()==PDT::Spin1    ) ++nv;
    objname += (**it).PDGName();
  }
  objname   += "Decayer";
  if(incoming->iSpin()==PDT::Spin0) {
    if(nf==4) classname += "StoFFFFDecayer";
    else      classname  = "";
  }
  else {
    classname="";
  }
  return classname;
}
#line 1 "./WeakCurrentDecayConstructor.cc"
// -*- C++ -*-
//
// WeakCurrentDecayConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WeakCurrentDecayConstructor class.
//

#include "WeakCurrentDecayConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/Vertex/AbstractFFVVertex.fh"
#include "DecayConstructor.h"

using namespace Herwig;
using ThePEG::Helicity::VertexBasePtr;

IBPtr WeakCurrentDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr WeakCurrentDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void WeakCurrentDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  model_ = dynamic_ptr_cast<Ptr<Herwig::StandardModel>::pointer>(generator()->standardModel());
  unsigned int isize=decayTags_.size();
  if(isize!=_norm .size()||isize!=_current.size())
    throw InitException() << "Invalid sizes for the decay mode vectors in "
			  << " WeakCurrentDecayConstructor " 
			  << decayTags_.size() << " " << _norm.size() << " " 
			  << _current.size() << Exception::runerror;
  // get the particles from the tags
  for(unsigned int ix=0;ix<decayTags_.size();++ix) {
    _current[ix]->init();
    particles_.push_back(vector<tPDPtr>());
    string tag=decayTags_[ix];
    do {
      string::size_type next = min(tag.find(','), tag.find(';'));
      particles_.back().push_back(generator()->findParticle(tag.substr(0,next)));
      if(!particles_.back().back()) 
	throw Exception() << "Failed to find particle " << tag.substr(0,next)
			  << " in DecayMode " << decayTags_[ix]
			  << " in WeakCurrentDecayConstructor::doinit()"
			  << Exception::runerror;
      if(tag[next]==';') break;
      tag = tag.substr(next+1);
    }
    while(true);
  }
}

void WeakCurrentDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << ounit(massCut_,GeV) << decayTags_ << particles_ << _norm << _current;
}

void WeakCurrentDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(massCut_,GeV) >> decayTags_ >> particles_ >> _norm >> _current;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<WeakCurrentDecayConstructor,NBodyDecayConstructorBase>
describeHerwigWeakCurrentDecayConstructor("Herwig::WeakCurrentDecayConstructor", "Herwig.so");

void WeakCurrentDecayConstructor::Init() {

  static ClassDocumentation<WeakCurrentDecayConstructor> documentation
    ("The WeakCurrentDecayConstructor class implemets the decay of BSM particles "
     "to low mass hadronic states using the Weak current");

  static ParVector<WeakCurrentDecayConstructor,string> interfaceDecayModes
    ("DecayModes",
     "The decays of the weak current",
     &WeakCurrentDecayConstructor::decayTags_, -1, "", "", "",
     false, false, Interface::nolimits);

  static ParVector<WeakCurrentDecayConstructor,double> interfaceNormalisation
    ("Normalisation",
     "The normalisation of the different modes",
     &WeakCurrentDecayConstructor::_norm, -1, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static RefVector<WeakCurrentDecayConstructor,WeakCurrent> interfaceCurrent
    ("Current",
     "The current for the decay mode",
     &WeakCurrentDecayConstructor::_current, -1, false, false, true, false, false);

  static Parameter<WeakCurrentDecayConstructor,Energy> interfaceMassCut
    ("MassCut",
     "The maximum mass difference for the decay",
     &WeakCurrentDecayConstructor::massCut_, GeV, 5.0*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void WeakCurrentDecayConstructor::DecayList(const set<PDPtr,MassOrdering> & part) {
  if( part.empty() ) return;
  unsigned int nv(model_->numberOfVertices());
  for(set<PDPtr,MassOrdering>::const_iterator ip=part.begin();ip!=part.end();++ip) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      for(unsigned int ilist = 0; ilist < 3; ++ilist) { 
	vector<TwoBodyDecay> decays =
	  createModes(*ip, model_->vertex(iv),ilist);
	if(!decays.empty()) createDecayMode(decays);
      }
    }
  }
}
  
vector<TwoBodyDecay> WeakCurrentDecayConstructor::createModes(const PDPtr inpart,
							      const VertexBasePtr vert,
							      unsigned int ilist) {
  int id = inpart->id();
  if( !vert->isIncoming(inpart) || vert->getNpoint() != 3 )
    return vector<TwoBodyDecay>();
  Energy m1(inpart->mass());
  vector<tPDPtr> decaylist;
  decaylist = vert->search(ilist,inpart);
  tPDVector::size_type nd = decaylist.size();
  vector<TwoBodyDecay> decays;
  for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
    tPDPtr pa(decaylist[i]), pb(decaylist.at(i + 1)), 
      pc(decaylist.at(i + 2));
    if( pb->id() == id ) swap(pa, pb);
    if( pc->id() == id ) swap(pa, pc);
    //One of the products must be a W
    Energy mp(ZERO);
    if( abs(pb->id()) == ParticleID::Wplus )
      mp = pc->mass();
    else if( abs(pc->id()) == ParticleID::Wplus )
      mp = pb->mass();
    else 
      continue;
    //allowed on-shell decay and passes mass cut
    if( m1 >= pb->mass() + pc->mass() ) continue;
    if( m1 < mp ) continue;
    if( m1 - mp >= massCut_ ) continue;
    //vertices are defined with all particles incoming
    if( pb->CC() ) pb = pb->CC();
    if( pc->CC() ) pc = pc->CC();
    decays.push_back( TwoBodyDecay(inpart,pb, pc, vert) );
    if(abs(decays.back().children_.second->id())!=ParticleID::Wplus)
      swap(decays.back().children_.first,decays.back().children_.second);
    assert(abs(decays.back().children_.second->id())==ParticleID::Wplus);
  }
  return decays;
}

GeneralCurrentDecayerPtr  WeakCurrentDecayConstructor::createDecayer(PDPtr in, PDPtr out1,
								     vector<tPDPtr> outCurrent,
								     VertexBasePtr vertex,
								     WeakCurrentPtr current) {
  string name;
  using namespace ThePEG::Helicity::VertexType;
  switch(vertex->getName()) {
  case FFV : 
    name = "FFVCurrentDecayer";
    break;
  default :
    ostringstream message;
    message << "Invalid vertex for decays of " << in->PDGName() << " -> " << out1->PDGName() 
	    << " via weak current " << vertex->fullName() << "\n";
    generator()->logWarning(NBodyDecayConstructorError(message.str(),
						       Exception::warning));
    return GeneralCurrentDecayerPtr();
  }
  ostringstream fullname;
  fullname << "/Herwig/Decays/" << name << "_" << in->PDGName() << "_"
	   << out1->PDGName();
  for(unsigned int ix=0;ix<outCurrent.size();++ix)
    fullname  << "_" << outCurrent[ix]->PDGName();
  string classname = "Herwig::" + name;
  GeneralCurrentDecayerPtr decayer = dynamic_ptr_cast<GeneralCurrentDecayerPtr>
    (generator()->preinitCreate(classname,fullname.str()));
  decayer->setDecayInfo(in,out1,outCurrent,vertex,current,massCut_);
  // set decayer options from base class
  setDecayerInterfaces(fullname.str());
  // initialize the decayer
  decayer->init();
  // return the decayer
  return decayer;
}

void WeakCurrentDecayConstructor::
createDecayMode(vector<TwoBodyDecay> & decays) {
  assert(!decays.empty());
  for(unsigned int ix = 0; ix < decays.size(); ++ix) {
    PDVector particles(3);
    particles[0] = decays[ix].parent_;
    particles[1] = decays[ix].children_.first ;
    bool Wplus=decays[ix].children_.second->id()==ParticleID::Wplus;
    for(unsigned int iy=0;iy<_current.size();++iy) {
      particles.resize(2);
      vector<tPDPtr> wprod=particles_[iy];
      int icharge=0;
      Energy msum = particles[0]->mass()-particles[1]->mass();
      for(unsigned int iz=0;iz<wprod.size();++iz) {
	icharge += wprod[iz]->iCharge();
	msum -=wprod[iz]->mass();
      }
      if(msum<=ZERO) continue;
      bool cc = (Wplus&&icharge==-3)||(!Wplus&&icharge==3);
      OrderedParticles outgoing;
      outgoing.insert(particles[1]);
      for(unsigned int iz=0;iz<wprod.size();++iz) {
 	if(cc&&wprod[iz]->CC())  wprod[iz]=wprod[iz]->CC();
	outgoing.insert(wprod[iz]);
      }
      // check outgoing particles initialised
      for(unsigned int iz=0;iz<wprod.size();++iz) wprod[iz]->init();
      // create the tag for the decay mode
      string tag = particles[0]->PDGName() + "->";
      OrderedParticles::const_iterator it = outgoing.begin();
      do {
	tag += (**it).name();
	++it;
	if(it!=outgoing.end()) tag +=",";
	else                   tag +=";";
      }
      while(it!=outgoing.end());
      // find the decay mode
      tDMPtr dm= generator()->findDecayMode(tag);
      if( !dm && createDecayModes() ) {
	// create the decayer
	GeneralCurrentDecayerPtr decayer = createDecayer(particles[0],particles[1],
							 wprod,decays[ix].vertex_,
							 _current[iy]);
	if(!decayer) continue;
	// calculate the width
	Energy pWidth = _norm[iy]*decayer->partialWidth(particles[0],particles[1],wprod);
	if(pWidth<=ZERO) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	  continue;
	}
	tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
	if(!ndm) throw NBodyDecayConstructorError() 
		   << "WeakCurrentDecayConstructor::createDecayMode - Needed to create "
		   << "new decaymode but one could not be created for the tag " 
		   << tag
		   << Exception::warning;
	generator()->preinitInterface(ndm, "Decayer", "set",
				      decayer->fullName());
	generator()->preinitInterface(ndm, "Active", "set", "Yes");
	setBranchingRatio(ndm, pWidth);
	particles[0]->stable(false);
	if(ndm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
      else if (dm) {
	// create the decayer
	GeneralCurrentDecayerPtr decayer = createDecayer(particles[0],particles[1],
							 wprod,decays[ix].vertex_,
							 _current[iy]);
	if(!decayer) continue;
	generator()->preinitInterface(dm, "Decayer", "set", decayer->fullName());
	particles[0]->stable(false);
	if(createDecayModes()) {
	  // calculate the width
	  Energy pWidth = _norm[iy]*decayer->partialWidth(particles[0],particles[1],wprod);
	  if(pWidth<=ZERO) {
	    generator()->preinitInterface(decayer->fullName(),
					  "Initialize", "set","0");
	    continue;
	  }
	  generator()->preinitInterface(dm, "Active", "set", "Yes");
	  particles[0]->width(particles[0]->width()*(1.-dm->brat()));
	  setBranchingRatio(dm, pWidth);
	}
	if(dm->brat()<decayConstructor()->minimumBR()) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	}
      }
    }
  }
}
#line 1 "./VectorCurrentDecayConstructor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorCurrentDecayConstructor class.
//

#include "VectorCurrentDecayConstructor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/General/VectorCurrentDecayer.h"


namespace {
struct ParticleOrdering {
  /**
   *  Operator for the ordering
   * @param p1 The first ParticleData object
   * @param p2 The second ParticleData object
   */
  bool operator() (tcPDPtr p1, tcPDPtr p2) const {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};
}

using namespace Herwig;

IBPtr VectorCurrentDecayConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr VectorCurrentDecayConstructor::fullclone() const {
  return new_ptr(*this);
}

void VectorCurrentDecayConstructor::persistentOutput(PersistentOStream & os) const {
  os << ounit(massCut_,GeV) << current_;
}

void VectorCurrentDecayConstructor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(massCut_,GeV) >> current_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorCurrentDecayConstructor,NBodyDecayConstructorBase>
  describeHerwigVectorCurrentDecayConstructor("Herwig::VectorCurrentDecayConstructor", "Herwig.so");

void VectorCurrentDecayConstructor::Init() {

  static ClassDocumentation<VectorCurrentDecayConstructor> documentation
    ("The VectorCurrentDecayConstructor class constructs the decays of low mass vector bosons"
     " to hadrons via the weak current");

  static RefVector<VectorCurrentDecayConstructor,WeakCurrent> interfaceCurrent
    ("Current",
     "The current for the decay mode",
     &VectorCurrentDecayConstructor::current_, -1, false, false, true, false, false);

  static Parameter<VectorCurrentDecayConstructor,Energy> interfaceMassCut
    ("MassCut",
     "The maximum mass difference for the decay",
     &VectorCurrentDecayConstructor::massCut_, GeV, 2.0*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

}

void VectorCurrentDecayConstructor::doinit() {
  NBodyDecayConstructorBase::doinit();
  model_ = dynamic_ptr_cast<Ptr<Herwig::StandardModel>::pointer>(generator()->standardModel());
}

void VectorCurrentDecayConstructor::DecayList(const set<PDPtr,MassOrdering> & particles) {
  if( particles.empty() ) return;
  unsigned int nv(model_->numberOfVertices());
  for(PDPtr part : particles) {
    if(part->iSpin()!=PDT::Spin1) continue;
    if(part->iCharge()!=0) continue;
    bool foundD(false),foundU(false),foundS(false);
    if(part->mass()>massCut_) continue;
    for(unsigned int iv = 0; iv < nv; ++iv) {
      VertexBasePtr vertex = model_->vertex(iv);
      if( !vertex->isIncoming(part) || vertex->getNpoint() != 3 ) continue;
      for(unsigned int iloc = 0;iloc < 3; ++iloc) {
	vector<long> ext = vertex->search(iloc, part->id());
	if(ext.empty()) continue;
	for(unsigned int ioff=0;ioff<ext.size();ioff+=3) {
	  if(iloc!=2) assert(false);
	  if(abs(ext[ioff])==1 && abs(ext[ioff+1])==1 &&  ext[ioff]==-ext[ioff+1]) {
	    foundD = true;
	  }
	  else if(abs(ext[ioff])==2 && abs(ext[ioff+1])==2 &&  ext[ioff]==-ext[ioff+1]) {
	    foundU = true;
	  }
	  else if(abs(ext[ioff])==3 && abs(ext[ioff+1])==3 &&  ext[ioff]==-ext[ioff+1]) {
	    foundS = true;
	  }
	}
      }
    }
    if(!foundD && !foundU && !foundS) continue;
    for(tWeakCurrentPtr current : current_) {
      current->init();
      for(unsigned int imode=0;imode<current->numberOfModes();++imode) {
	// get the external particles for this mode
	int iq(0),ia(0);
	tPDVector out = current->particles(0,imode,iq,ia);
	current->decayModeInfo(imode,iq,ia);
	if(iq==2&&ia==-2) continue;
	// order the particles
	bool skip=false;
	for(unsigned int ix=0;ix<out.size();++ix) {
	  if(!out[ix]) {
	    skip=true;
	    break;
	  }
	}
	if(skip) continue;
	multiset<tcPDPtr,ParticleOrdering> outgoing(out.begin(),out.end());
	Energy minMass(ZERO);
	string tag = part->PDGName() + "->";
	bool first=false;
	int charge(0);
	for(tcPDPtr part : outgoing) {
	  if(!first)
	    first=true;
	  else
	    tag+=",";
	  tag+=part->PDGName();
	  minMass+=part->mass();
	  charge+=part->iCharge();
	}
	tag+=";";
	if(minMass>part->mass()) continue;
	if(charge!=0) continue;
	// create the decayer
	ostringstream fullname;
	fullname << "/Herwig/Decays/DMMediator_" << part->PDGName();
	for(tcPDPtr part : out)
	  fullname  << "_" << part->PDGName();
	string classname = "Herwig::VectorCurrentDecayer";
	VectorCurrentDecayerPtr decayer = dynamic_ptr_cast<VectorCurrentDecayerPtr>
	  (generator()->preinitCreate(classname,fullname.str()));
	decayer->setDecayInfo(part,out,current);
	// // set decayer options from base class
	// setDecayerInterfaces(fullname.str());
	// initialize the decayer
	decayer->init();
	// calculate the width
	Energy pWidth = decayer->partialWidth(part,out);
	if(pWidth<=ZERO) {
	  generator()->preinitInterface(decayer->fullName(),
					"Initialize", "set","0");
	  continue;
	}
	tDMPtr ndm = generator()->preinitCreateDecayMode(tag);
	generator()->preinitInterface(ndm, "Decayer", "set", decayer->fullName());
	part->stable(false);
	generator()->preinitInterface(ndm, "Active", "set", "Yes");
	setBranchingRatio(ndm, pWidth);
      }
    }
  }
}
#line 1 "./ResonantProcessConstructor.cc"
// -*- C++ -*-
//
// ResonantProcessConstructor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ResonantProcessConstructor class.
//

#include "ResonantProcessConstructor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

using namespace Herwig;

IBPtr ResonantProcessConstructor::clone() const {
  return new_ptr(*this);
}

IBPtr ResonantProcessConstructor::fullclone() const {
  return new_ptr(*this);
}

void ResonantProcessConstructor::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << intermediates_ << outgoing_
     << processOption_ << scaleChoice_ << scaleFactor_;
}

void ResonantProcessConstructor::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> intermediates_ >> outgoing_
     >> processOption_ >> scaleChoice_ >> scaleFactor_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ResonantProcessConstructor,HardProcessConstructor>
describeHerwigResonantProcessConstructor("Herwig::ResonantProcessConstructor", "Herwig.so");

void ResonantProcessConstructor::Init() {

  static ClassDocumentation<ResonantProcessConstructor> documentation
    ("This class is designed solely to contruct resonant processes using"
     "a provided set of intermediate particles");

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceOffshell
    ("Intermediates",
     "A vector of offshell particles for resonant diagrams",
     &ResonantProcessConstructor::intermediates_, -1, false, false, true,
     false);

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceIncoming
    ("Incoming",
     "A vector of incoming particles for resonant diagrams",
     &ResonantProcessConstructor::incoming_, -1, false, false, true,
     false);

  static RefVector<ResonantProcessConstructor, ParticleData> interfaceOutgoing
    ("Outgoing",
     "A vector of outgoing particles for resonant diagrams",
     &ResonantProcessConstructor::outgoing_, -1, false, false, true,
     false);

  static Switch<ResonantProcessConstructor,unsigned int> interfaceProcesses
    ("Processes",
     "Whether to generate inclusive or exclusive processes",
     &ResonantProcessConstructor::processOption_, 0, false, false);
  static SwitchOption interfaceProcessesSingleParticleInclusive
    (interfaceProcesses,
     "SingleParticleInclusive",
     "Require at least one particle from the list of outgoing particles"
     " in the hard process",
     0);
  static SwitchOption interfaceProcessesTwoParticleInclusive
    (interfaceProcesses,
     "TwoParticleInclusive",
     "Require that both the particles in the hard processes are in the"
     " list of outgoing particles",
     1);
  static SwitchOption interfaceProcessesExclusive
    (interfaceProcesses,
     "Exclusive",
     "Require that both the particles in the hard processes are in the"
     " list of outgoing particles in every hard process",
     2);
  static SwitchOption interfaceProcessesInclusive
    (interfaceProcesses,
     "Inclusive",
     "Generate all modes which are allowed for the on-shell intermediate particle",
     3);
  static SwitchOption interfaceProcessesVeryExclusive
    (interfaceProcesses,
     "VeryExclusive",
     "Require that both the incoming and outgoing particles in the hard processes are in the"
     " list of outgoing particles in every hard process",
     4);

  static Switch<ResonantProcessConstructor,unsigned int> interfaceScaleChoice
    ("ScaleChoice",
     "&ResonantProcessConstructor::scaleChoice_",
     &ResonantProcessConstructor::scaleChoice_, 1, false, false);
  static SwitchOption interfaceScaleChoicesHat
    (interfaceScaleChoice,
     "sHat",
     "Always use sHat",
     1);
  static SwitchOption interfaceScaleChoiceTransverseMass
    (interfaceScaleChoice,
     "TransverseMass",
     "Always use the transverse mass",
     2);
  static SwitchOption interfaceScaleChoiceGeometicMean
    (interfaceScaleChoice,
     "MaxMT",
     "Use the maximum of m^2+p_T^2 for the two particles",
     3);

  static Parameter<ResonantProcessConstructor,double> interfaceScaleFactor
    ("ScaleFactor",
     "The prefactor used in the scale calculation. The scale used is"
     " sHat multiplied by this prefactor",
     &ResonantProcessConstructor::scaleFactor_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

}

void ResonantProcessConstructor::doinit() {
  HardProcessConstructor::doinit();
  if((processOption_==2 || processOption_==4) &&
     outgoing_.size()!=2)
    throw InitException()
      << "Exclusive processes require exactly"
      << " two outgoing particles but " << outgoing_.size()
      << "have been inserted in ResonantProcessConstructor::doinit()."
      << Exception::runerror;
  if(processOption_==4 && incoming_.size()!=2)
    throw InitException()
      << "Exclusive processes require exactly"
      << " two incoming particles but " << incoming_.size()
      << "have been inserted in ResonantProcessConstructor::doinit()."
      << Exception::runerror;
}

void ResonantProcessConstructor::constructDiagrams() {
  size_t ninc = incoming_.size() , ninter = intermediates_.size();
  if(ninc == 0 || ninter == 0  || !subProcess() ) return;
  // find the incoming particle pairs
  vector<tPDPair> incPairs;
  for(PDVector::size_type i = 0; i < ninc; ++i) {
    for(PDVector::size_type j = 0; j < ninc; ++j) {
      tPDPair inc = make_pair(incoming_[i], incoming_[j]);
      if( (inc.first->iSpin() > inc.second->iSpin()) ||
	  (inc.first->iSpin() == inc.second->iSpin() &&
	   inc.first->id() < inc.second->id()) )
	swap(inc.first, inc.second);
      if( !HPC_helper::duplicateIncoming(inc,incPairs) ) {
	incPairs.push_back(inc);
      }
    }
  }
  size_t nvertices = model()->numberOfVertices();
  //To construct resonant diagrams loop over the incoming particles, intermediates
  //and vertices to find allowed diagrams. Need to exclude the diagrams that have
  //the intermediate as an external particle as well
  for(vector<tcPDPair>::size_type is = 0; is < incPairs.size(); ++is) {
    tPDPair ppi = incPairs[is];
    for(vector<PDPtr>::size_type ik = 0; ik < ninter ; ++ik) {
      long ipart = intermediates_[ik]->id();
      for(size_t iv = 0; iv < nvertices; ++iv) {
	VBPtr vertex = model()->vertex(iv);
	if(vertex->getNpoint() > 3) continue;
	long part1 = ppi.first->CC() ? -ppi.first->id() : ppi.first->id();
	long part2 = ppi.second->CC() ? -ppi.second->id() : ppi.second->id();
	if(vertex->allowed(part1, part2, ipart) ||
	   vertex->allowed(part1, ipart, part2) ||
	   vertex->allowed(part2, part1, ipart) ||
	   vertex->allowed(part2, ipart, part1) ||
	   vertex->allowed(ipart, part1, part2) ||
	   vertex->allowed(ipart, part2, part1) ) {
	  constructVertex2(make_pair(ppi.first->id(), ppi.second->id()), vertex,
			   intermediates_[ik]);
	}
      }
    }
  }
  //Create matrix element for each process
  const HPDVector::size_type ndiags = diagrams_.size();
  vector<string> diagLib;
  for(HPDVector::size_type ix = 0; ix < ndiags; ++ix) {
    vector<tcPDPtr> extpart(4);
    extpart[0] = getParticleData(diagrams_[ix].incoming.first);
    extpart[1] = getParticleData(diagrams_[ix].incoming.second);
    extpart[2] = getParticleData(diagrams_[ix].outgoing.first);
    extpart[3] = getParticleData(diagrams_[ix].outgoing.second);
    string objectname("");
    string classname = MEClassname(extpart, diagrams_[ix].intermediate, objectname);
    if(std::find(diagLib.begin(), diagLib.end(), objectname) == diagLib.end()) {
      createMatrixElement(diagrams_[ix]);
      diagLib.push_back(objectname);
    }
    else {
      cerr<<"Warning: ResonantProcessConstructor::constructDiagrams - "
          <<"UFO model tried to produce douplicate diagrams for "
          << extpart[0]->PDGName() << " " << extpart[1]->PDGName() << "->"
          << extpart[2]->PDGName() << " " << extpart[3]->PDGName() << " "
          << "for the class: \"" << classname <<" "
          << "in object " << objectname << "\". "
          <<"neglectiting douplicat diagram.\n";
      continue;
    }
  }
}

void ResonantProcessConstructor::
constructVertex2(IDPair in, VertexBasePtr vertex,
		 PDPtr partc) {
  //We have the left hand part of the diagram, just need all the possibilities
  //for the RHS
  size_t nvertices = model()->numberOfVertices();
  if(processOption_!=3) {
    for(size_t io = 0; io < outgoing_.size(); ++io) {
      tcPDPtr outa = outgoing_[io];
      for(size_t iv = 0; iv < nvertices; ++iv) {
	VBPtr vertex2 = model()->vertex(iv);
	if(vertex2->getNpoint() > 3) continue;
	tPDSet outb = search(vertex2, partc->id(), incoming, outa->id(), outgoing,
			     outgoing);
	for(tPDSet::const_iterator ita = outb.begin(); ita != outb.end(); ++ita)
	  makeResonantDiagram(in, partc, outa->id(),(**ita).id(),
			      make_pair(vertex, vertex2));
      }
    }
  }
  else {
    long idRes = !partc->CC() ? partc->id() : partc->CC()->id();
    for(size_t iv = 0; iv < nvertices; ++iv) {
      VBPtr vertex2 = model()->vertex(iv);
      if(vertex2->getNpoint() > 3) continue;
      for(unsigned int ix = 0;ix < 3; ++ix) {
	vector<long> pdlist = vertex2->search(ix, idRes);
	for(unsigned int iy=0;iy<pdlist.size();iy+=3) {
	  long out1 = ix==0 ? pdlist.at(iy+1) : pdlist.at(iy  );
	  long out2 = ix==2 ? pdlist.at(iy+1) : pdlist.at(iy+2);
	  if(partc->mass() < getParticleData(out1)->mass() +
	     getParticleData(out2)->mass()) continue;
	  makeResonantDiagram(in, partc, out1, out2,
			      make_pair(vertex, vertex2));
	}
      }
    }
  }
}

void ResonantProcessConstructor::
makeResonantDiagram(IDPair in, PDPtr offshell, long outa, long outb,
		     VBPair vertpair) {
  assert(vertpair.first && vertpair.second);
  if( abs(outa) == abs(offshell->id()) ||
      abs(outb) == abs(offshell->id())) return;
  HPDiagram newdiag(in,make_pair(outa,outb));
  newdiag.intermediate = offshell;
  newdiag.vertices = vertpair;
  newdiag.channelType = HPDiagram::sChannel;
  fixFSOrder(newdiag);
  assignToCF(newdiag);
  if(duplicate(newdiag,diagrams_)) return;
  // if inclusive enforce the exclusivity
  if(processOption_==1) {
    PDVector::const_iterator loc =
      std::find(outgoing_.begin(),outgoing_.end(),
		getParticleData(newdiag.outgoing. first));
    if(loc==outgoing_.end()) return;
    loc =
      std::find(outgoing_.begin(),outgoing_.end(),
		getParticleData(newdiag.outgoing.second));
    if(loc==outgoing_.end()) return;
  }
  else if(processOption_==2 || processOption_==4 ) {
    if(!((newdiag.outgoing. first==outgoing_[0]->id()&&
	  newdiag.outgoing.second==outgoing_[1]->id())||
	 (newdiag.outgoing. first==outgoing_[1]->id()&&
	  newdiag.outgoing.second==outgoing_[0]->id())))
      return;
    if(processOption_==4) {
      if(!((newdiag.incoming. first==incoming_[0]->id()&&
	    newdiag.incoming.second==incoming_[1]->id())||
	   (newdiag.incoming. first==incoming_[1]->id()&&
	    newdiag.incoming.second==incoming_[0]->id())))
	return;
    }
  }
  // add to the list
  if(checkOrder(newdiag)) diagrams_.push_back(newdiag);
}

set<tPDPtr>
ResonantProcessConstructor::search(VBPtr vertex, long part1, direction d1,
				   long part2, direction d2, direction d3) {
  if(d1 == incoming && getParticleData(part1)->CC()) part1 = -part1;
  if(d2 == incoming && getParticleData(part2)->CC()) part2 = -part2;
  vector<long> ext;
  tPDSet third;
  for(unsigned int ix = 0;ix < 3; ++ix) {
    vector<long> pdlist = vertex->search(ix, part1);
    ext.insert(ext.end(), pdlist.begin(), pdlist.end());
  }
  for(unsigned int ix = 0; ix < ext.size(); ix += 3) {
    long id0 = ext.at(ix);
    long id1 = ext.at(ix+1);
    long id2 = ext.at(ix+2);
    int pos;
    if((id0 == part1 && id1 == part2) ||
       (id0 == part2 && id1 == part1))
      pos = ix + 2;
    else if((id0 == part1 && id2 == part2) ||
	    (id0 == part2 && id2 == part1))
      pos = ix + 1;
    else if((id1 == part1 && id2 == part2) ||
	    (id1 == part2 && id2 == part1))
      pos = ix;
    else
      pos = -1;
    if(pos >= 0) {
      tPDPtr p = getParticleData(ext[pos]);
      if(d3 == incoming && p->CC()) p = p->CC();
      third.insert(p);
    }
  }

  return third;
}

IDPair ResonantProcessConstructor::
find(long part, const vector<PDPtr> & out) const {
  vector<PDPtr>::size_type iloc(0);
  bool found(false);
  do {
    if(out[iloc]->id() == part) found = true;
    else ++iloc;
  }
  while(found == false && iloc < out.size());
  //found offshell
  IDPair outids;
  if(iloc == 0)
    outids = make_pair(out[1]->id(), out[2]->id());
  else if(iloc == 1)
    outids = make_pair(out[0]->id(), out[2]->id());
  else
    outids = make_pair(out[0]->id(), out[1]->id());
  return outids;
}

void ResonantProcessConstructor::
createMatrixElement(const HPDiagram & diag) const {
  vector<tcPDPtr> extpart(4);
  extpart[0] = getParticleData(diag.incoming.first);
  extpart[1] = getParticleData(diag.incoming.second);
  extpart[2] = getParticleData(diag.outgoing.first);
  extpart[3] = getParticleData(diag.outgoing.second);
  string objectname ("/Herwig/MatrixElements/");
  string classname = MEClassname(extpart, diag.intermediate, objectname);
  GeneralHardMEPtr matrixElement = dynamic_ptr_cast<GeneralHardMEPtr>
    (generator()->preinitCreate(classname, objectname));
  if( !matrixElement ) {
    throw RPConstructorError()
      << "ResonantProcessConstructor::createMatrixElement "
      << "- No matrix element object could be created for "
      << "the process "
      << extpart[0]->PDGName() << " " << extpart[0]->iSpin() << ","
      << extpart[1]->PDGName() << " " << extpart[1]->iSpin() << "->"
      << extpart[2]->PDGName() << " " << extpart[2]->iSpin() << ","
      << extpart[3]->PDGName() << " " << extpart[3]->iSpin()
      << ".  Constructed class name: \"" << classname << "\"\n"
      << Exception::warning;
    return;
  }
  matrixElement->setProcessInfo(HPDVector(1, diag),
				colourFlow(extpart), debug(),scaleChoice_-1,scaleFactor_);
  generator()->preinitInterface(subProcess(), "MatrixElements",
				subProcess()->MEs().size(),
				"insert", matrixElement->fullName());
}

string ResonantProcessConstructor::
MEClassname(const vector<tcPDPtr> & extpart, tcPDPtr inter,
	    string & objname) const {
 string classname("Herwig::ME");
  for(vector<tcPDPtr>::size_type ix = 0; ix < extpart.size(); ++ix) {
    if(ix == 2) classname += "2";
    if(extpart[ix]->iSpin() == PDT::Spin0) classname += "s";
    else if(extpart[ix]->iSpin() == PDT::Spin1) classname += "v";
    else if(extpart[ix]->iSpin() == PDT::Spin1Half) classname += "f";
    else if(extpart[ix]->iSpin() == PDT::Spin2) classname += "t";
    else
      throw RPConstructorError()
	<< "MEClassname() : Encountered an unknown spin for "
	<< extpart[ix]->PDGName() << " while constructing MatrixElement "
	<< "classname " << extpart[ix]->iSpin() << Exception::warning;
  }
  objname += "ME" + extpart[0]->PDGName() + extpart[1]->PDGName() + "2"
    + inter->PDGName() + "2"
    + extpart[2]->PDGName() + extpart[3]->PDGName();
  return classname;
}
#line 1 "./VVSLoopVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSLoopVertex class.
//

#include "VVSLoopVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;
using namespace ThePEG;
namespace LT = Looptools;

IBPtr VVSLoopVertex::clone() const {
  return new_ptr(*this);
}

IBPtr VVSLoopVertex::fullclone() const {
  return new_ptr(*this);
}

void VVSLoopVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(masses,GeV) << type << couplings << Npart_;
}

void VVSLoopVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(masses,GeV) >> type >> couplings >> Npart_;
}
void VVSLoopVertex::dofinish() {
  if(loopToolsInit_) Looptools::ltexi();
  GeneralVVSVertex::dofinish();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVSLoopVertex,Helicity::GeneralVVSVertex>
describeHerwigVVSLoopVertex("Herwig::VVSLoopVertex", "Herwig.so");

void VVSLoopVertex::Init() {

  static ClassDocumentation<VVSLoopVertex> documentation
    ("The VVSLoopVertex class calculates the tensor integral"
     " coefficients using Looptools.");

}

void VVSLoopVertex::setCoupling(Energy2, tcPDPtr p1, tcPDPtr,tcPDPtr) {
  if(!loopToolsInit_) {
    Looptools::ltini();
    loopToolsInit_ = true;
  }
  //Kinematic invariants
  double ps2 = invariant(0,0) / MeV2;
  double pv1s = invariant(1,1) / MeV2;
  double pv2s = invariant(2,2) / MeV2;
  Complex a(0.),b(0.),c(0.),d(0.),e(0.),f(0.);
  for(unsigned int i = 0; i< Npart_;++i) {
    double lmass = masses[i] / MeV;
    double mls = sqr(lmass);
    // left coupling for fermions or the total coupling
    // for scalars or vectors
    Complex lc = couplings[i].first;
    // double p1p2 = 0.5*(ps2-pv1s-pv2s);
    Complex C0A = LT::C0i(LT::cc0,pv2s,ps2,pv1s,mls,mls,mls);
    long theC = LT::Cget(pv2s,ps2,pv1s,    mls,mls,mls);
    // Complex C1A  = LT::Cval(LT::cc1,theC);
    // Complex C2A  = LT::Cval(LT::cc2,theC);
    // Complex C11A = LT::Cval(LT::cc11,theC);
    Complex C12A = LT::Cval(LT::cc12,theC);
    // Complex C22A = LT::Cval(LT::cc22,theC);
    // Complex C00A = LT::Cval(LT::cc00,theC);
    Complex C0B = LT::C0i(LT::cc0,pv1s,ps2,pv2s,mls,mls,mls);
    theC = LT::Cget(pv1s,ps2,pv2s,    mls,mls,mls);
    // Complex C1B  = LT::Cval(LT::cc1,theC);
    // Complex C2B  = LT::Cval(LT::cc2,theC);
    // Complex C11B = LT::Cval(LT::cc11,theC);
    Complex C12B = LT::Cval(LT::cc12,theC);
    // Complex C22B = LT::Cval(LT::cc22,theC);
    // Complex C00B = LT::Cval(LT::cc00,theC);
    if(type[i] == PDT::Spin1Half) {
      Complex lpr = lc + couplings[i].second;
      Complex loop = 2.*lpr*lmass*( - 4.*(C12A + C12B) + C0A  + C0B );
      a -= loop;
      d += loop;
      f +=  2.*(lc - couplings[i].second)*lmass*(C0A+C0B);
      // a +=  2.*lpr*lmass*(2.*C12A-C0A+2.*C12B-C0B+
      // 			  (1. - pv1s*(C22A+C11B)- pv2s*(C11A+C22B) + mls*(C0A+C0B) )/p1p2);
      // b += 8.*lpr*lmass*(2.*C22A + C2A  +2.*C11B + C1B);
      // c += 2.*lpr*lmass*(- 4.*(C12A + C12B) - 2.*(C2A + C1A + C2B + C1B) - C0A - C0B);
      // d +=  2.*lpr*lmass*( - 4.*(C12A + C12B) + C0A  + C0B );
      // e +=  4.*lpr*lmass*( 2.*(C11A + C22B) + C1A  + C2B);
    }
    else if(type[i] == PDT::Spin1) {
      //Complex B0W(LT::B0(ps2 ,mls,mls));
      double mr(sqr(p1->mass()/masses[i]));
      Complex loop = 2.*lc*( -(6.+mr)*(C12A+C12B) + 4.*(C0A +C0B));
      a -= loop;
      d += loop;
      // a += lc*(-8.*(C0A+C0B) +(6.+mr)*(2.*(C00A+C00B)-B0W)/p1p2);
      // b += lc*(6.+mr)*(2.*(C22A+C11B)+C2A+C1B);
      // c += 0.5*lc*(6.+mr)*(  - 4.*(C12A+C12B) - 2.*(C1A+C2A+C1B+C2B) - C0A - C0B );  
      // d += 2.*lc*( -(6.+mr)*(C12A+C12B) + 4.*(C0A +C0B)); 
      // e += lc*(6.+mr)*( 2.*(C11A+C22B) + C1A + C2B );
    }    
    else if(type[i] == PDT::Spin0) {
      Complex loop = 4.*lc*(C12A+C12B);
      a -= loop;
      d += loop;
      // a += -2.*lc*( 2.*(C00A+C00B) - LT::B0(ps2,mls,mls))/p1p2;
      // b += -2.*lc* ( 2.*(C22A + C11B) + C2A + C1B);
      // c +=    -lc*(  - 4.*(C12A + C12B) - 2.*(C2A + C1A + C2B + C1B) - C0A - C0B);
      // d +=  4.*lc*(C12A+C12B);
      // e +=    -lc*( 2.*(C11A + C22B) + C1A + C2B );
    }
    else {
      throw Helicity::HelicityConsistencyError() 
	<< "SVVLoopVertex::setCoupling - Incorrect particle in SVV loop. "
	<< "Spin: " << type[i]
	<< Exception::runerror;
    }
  }
  //Looptools defines integrals differently
  double fact = 1./16./sqr(Constants::pi);
  a00(fact*a);
  a11(fact*b);
  a12(fact*c);
  a21(fact*d);
  a22(fact*e);
  aEp(fact*f);
}
#line 1 "./GenericHGGVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericHGGVertex class.
//

#include "GenericHGGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

GenericHGGVertex::GenericHGGVertex() : setup_(false), q2Last_(ZERO), coupLast_(0.), idLast_(0) {
  orderInGs(2);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

IBPtr GenericHGGVertex::clone() const {
  return new_ptr(*this);
}

IBPtr GenericHGGVertex::fullclone() const {
  return new_ptr(*this);
}

void GenericHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << bosons_ << setup_ << vertices_ << model_;
}

void GenericHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> bosons_ >> setup_ >> vertices_ >> model_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GenericHGGVertex,VVSLoopVertex>
describeHerwigGenericHGGVertex("Herwig::GenericHGGVertex", "Herwig.so");

void GenericHGGVertex::Init() {

  static ClassDocumentation<GenericHGGVertex> documentation
    ("The GenericHGGVertex class implements the coupling of"
     " the Higgs bosons to gluons in a generic model.");

  static RefVector<GenericHGGVertex,ParticleData> interfaceBosons
    ("Bosons",
     "The Higgs bosons in the model.",
     &GenericHGGVertex::bosons_, -1, false, false, true, false, false);

}

void GenericHGGVertex::doinit() {
  for(unsigned int ix=0;ix<bosons_.size();++ix) {
    addToList(21,21,bosons_[ix]->id());
  }
  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void GenericHGGVertex::setCoupling (Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				    tcPDPtr part3) {
  if(!setup_) initializeVertex();
  assert(part1->id()==ParticleID::g && part2->id()==ParticleID::g);
  // find the particles in the loop
  map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(part3);
  // check there are some
  if(it==vertices_.end()) {
    norm(0.);
    return;
  }
  Looptools::clearcache();
  // overall factor
  if( q2 != q2Last_ || coupLast_ == 0. || part3->id() != idLast_ ) {
    q2Last_ = q2;
    idLast_ = part3->id();
    coupLast_ = sqr(strongCoupling(q2));
    // loop over the loop particles
    masses.clear();
    type.clear();
    couplings.clear();
    setNParticles(it->second.size());
    for(unsigned int ix=0;ix<it->second.size();++ix) {
      masses.push_back(model_->mass(q2,it->second[ix].particle));
      if(it->second[ix].particle->iSpin()==PDT::Spin0) {
	type.push_back(PDT::Spin0);
	it->second[ix].scalar->setCoupling(q2,part3,it->second[ix].particle,it->second[ix].particle->CC());
	couplings.push_back(make_pair(0.5*it->second[ix].scalar->norm(),0.5*it->second[ix].scalar->norm()));
      }
      else if(it->second[ix].particle->iSpin()==PDT::Spin1Half) {
	type.push_back(PDT::Spin1Half);
	assert(it->second[ix].fermion);
	it->second[ix].fermion->setCoupling(q2,it->second[ix].particle,it->second[ix].particle->CC(),part3);
	Complex coupling = it->second[ix].fermion->norm();
	Complex lc = it->second[ix].fermion->left ();
	Complex rc = it->second[ix].fermion->right();
	couplings.push_back(make_pair(0.5*coupling*lc,0.5*coupling*rc));
      }
      else
	assert(false);
    }
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
  }
  norm(coupLast_);
}

void GenericHGGVertex::initializeVertex() {
  // get the model
  model_ = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if(!model_) throw InitException();
  // loop over all the vertices
  unsigned int nv(model_->numberOfVertices());
  for(unsigned int ib=0;ib<bosons_.size();++ib) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      // 3-point vertex with boson as incoming
      if( model_->vertex(iv)->getNpoint()>3) continue;
      if( !model_->vertex(iv)->isIncoming(bosons_[ib])) continue;
      for(unsigned int il = 0; il < 3; ++il) { 
	tPDVector decaylist = model_->vertex(iv)->search(il, bosons_[ib]);
	tPDVector::size_type nd = decaylist.size();
	for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
	  if( pb->id() == bosons_[ib]->id() ) swap(pa, pb);
	  if( pc->id() == bosons_[ib]->id() ) swap(pa, pc);
	  // check coloured and particle antiparticle
	  if( pb->CC() != pc || !pb->coloured() || !pc->coloured()) 
	    continue;
	  // only know how to deal with triplets
	  assert(pb->iColour()==PDT::Colour3 || pb->iColour()==PDT::Colour3bar);
	  // scalar loop
	  if(pb->iSpin()==PDT::Spin0) {
	    SSSVertexPtr vertex = dynamic_ptr_cast<SSSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,vertex,FFSVertexPtr()));
	    }
	    else {
	      vertices_.insert(make_pair(bosons_[ib],vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,vertex,FFSVertexPtr()))));
	    }
	  }
	  // fermion loop
	  else if(pb->iSpin()==PDT::Spin1Half) {
	    FFSVertexPtr vertex = dynamic_ptr_cast<FFSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),vertex));
	    }
	    else {
	     vertices_.insert(make_pair(bosons_[ib],vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),vertex))));
	    }
	  }
	  else
	    assert(false);
	}
      }
    }
  }
  // set up now
  setup_ = true;
}
#line 1 "./GenericHPPVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericHPPVertex class.
//

#include "GenericHPPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

GenericHPPVertex::GenericHPPVertex() : setup_(false), q2Last_(ZERO), coupLast_(0.), idLast_(0) {
  orderInGs(0);
  orderInGem(3);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr GenericHPPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr GenericHPPVertex::fullclone() const {
  return new_ptr(*this);
}

void GenericHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << bosons_ << setup_ << vertices_ << model_;
}

void GenericHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> bosons_ >> setup_ >> vertices_ >> model_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GenericHPPVertex,VVSLoopVertex>
describeHerwigGenericHPPVertex("Herwig::GenericHPPVertex", "Herwig.so");

void GenericHPPVertex::Init() {

  static ClassDocumentation<GenericHPPVertex> documentation
    ("The GenericHPPVertex class implements the coupling of"
     " the Higgs bosons to gluons in a generic model.");

  static RefVector<GenericHPPVertex,ParticleData> interfaceBosons
    ("Bosons",
     "The Higgs bosons in the model.",
     &GenericHPPVertex::bosons_, -1, false, false, true, false, false);

}

void GenericHPPVertex::doinit() {
  for(unsigned int ix=0;ix<bosons_.size();++ix) {
    addToList(22,22,bosons_[ix]->id());
  }
  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void GenericHPPVertex::setCoupling (Energy2 q2, tcPDPtr part1, tcPDPtr part2,
				    tcPDPtr part3) {
  if(!setup_) initializeVertex();
  assert(part1->id()==ParticleID::gamma && part2->id()==ParticleID::gamma);
  // find the particles in the loop
  map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(part3);
  // check there are some
  if(it==vertices_.end()) {
    norm(0.);
    return;
  }
  Looptools::clearcache();
  // overall factor
  if( q2 != q2Last_ || coupLast_ == 0. || part3->id() != idLast_ ) {
    q2Last_ = q2;
    idLast_ = part3->id();
    coupLast_ = sqr(electroMagneticCoupling(q2));
    // loop over the loop particles
    masses.clear();
    type.clear();
    couplings.clear();
    setNParticles(it->second.size());
    for(unsigned int ix=0;ix<it->second.size();++ix) {
      masses.push_back(model_->mass(q2,it->second[ix].particle));
      // charge and colour factor
      double fact(1.);
      fact *= sqr(double(it->second[ix].particle->iCharge())/3.);
      if(it->second[ix].particle->iColour()==PDT::Colour3||
	 it->second[ix].particle->iColour()==PDT::Colour3bar)
	fact *=3.;
      else if(it->second[ix].particle->iColour()==PDT::Colour0)
	fact *= 1.;
      else {
	assert(false);
      }
      // spin-0
      if(it->second[ix].particle->iSpin()==PDT::Spin0) {
   	type.push_back(PDT::Spin0);
   	it->second[ix].scalar->setCoupling(q2,part3,it->second[ix].particle,it->second[ix].particle->CC());
  	couplings.push_back(make_pair(fact*it->second[ix].scalar->norm(),
				      fact*it->second[ix].scalar->norm()));
      }
      else if(it->second[ix].particle->iSpin()==PDT::Spin1Half) {
   	type.push_back(PDT::Spin1Half);
   	assert(it->second[ix].fermion);
   	it->second[ix].fermion->setCoupling(q2,it->second[ix].particle,it->second[ix].particle->CC(),part3);
   	Complex coupling = fact*it->second[ix].fermion->norm();
   	Complex lc = it->second[ix].fermion->left ();
   	Complex rc = it->second[ix].fermion->right();
   	couplings.push_back(make_pair(coupling*lc,coupling*rc));
      }
      else if(it->second[ix].particle->iSpin()==PDT::Spin1) {
   	type.push_back(PDT::Spin1);
   	assert(it->second[ix].vector);
   	it->second[ix].vector->setCoupling(q2,it->second[ix].particle,it->second[ix].particle->CC(),part3);
  	couplings.push_back(make_pair(fact*it->second[ix].vector->norm(),
				      fact*it->second[ix].vector->norm()));
      }
      else
   	assert(false);
    }
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
  }
  norm(coupLast_);
}

void GenericHPPVertex::initializeVertex() {
  // get the model
  model_ = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if(!model_) throw InitException();
  // loop over all the vertices
  unsigned int nv(model_->numberOfVertices());
  for(unsigned int ib=0;ib<bosons_.size();++ib) {
    for(unsigned int iv = 0; iv < nv; ++iv) {
      // 3-point vertex with boson as incoming
      if( model_->vertex(iv)->getNpoint()>3) continue;
      if( !model_->vertex(iv)->isIncoming(bosons_[ib])) continue;
      for(unsigned int il = 0; il < 3; ++il) { 
	tPDVector decaylist = model_->vertex(iv)->search(il, bosons_[ib]);
	tPDVector::size_type nd = decaylist.size();
	for( tPDVector::size_type i = 0; i < nd; i += 3 ) {
	  tPDPtr pa(decaylist[i]), pb(decaylist[i + 1]), pc(decaylist[i + 2]);
	  if( pb->id() == bosons_[ib]->id() ) swap(pa, pb);
	  if( pc->id() == bosons_[ib]->id() ) swap(pa, pc);
	  // check coloured and particle antiparticle
	  if( pb->CC() != pc || !pb->charged() || !pc->charged()) 
	    continue;
	  // // scalar loop
	  if(pb->iSpin()==PDT::Spin0) {
	    SSSVertexPtr vertex = dynamic_ptr_cast<SSSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,vertex,FFSVertexPtr(),VVSVertexPtr()));
	    }
	    else {
	      vertices_.insert(make_pair(bosons_[ib],
					 vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,vertex,
									   FFSVertexPtr(),VVSVertexPtr()))));
	    }
	  }
	  // fermion loop
	  else if(pb->iSpin()==PDT::Spin1Half) {
	    FFSVertexPtr vertex = dynamic_ptr_cast<FFSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),vertex,VVSVertexPtr()));
	    }
	    else {
	      vertices_.insert(make_pair(bosons_[ib],
					 vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),
									   vertex,VVSVertexPtr()))));
	    }
	  }
	  // vector loop
	  else if(pb->iSpin()==PDT::Spin1) {
	    VVSVertexPtr vertex = dynamic_ptr_cast<VVSVertexPtr>(model_->vertex(iv));
	    if(!vertex) continue;
	    map<cPDPtr,vector<Interaction> >::iterator it = vertices_.find(bosons_[ib]);
	    if(it!=vertices_.end()) {
	      it->second.push_back(Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),FFSVertexPtr(),vertex));
	    }
	    else {
	      vertices_.insert(make_pair(bosons_[ib],
					 vector<Interaction>(1,Interaction(pb->id()>0?pb:pc,SSSVertexPtr(),
									   FFSVertexPtr(),vertex))));
	    }
	  }
	  else
	    assert(false);
	}
      }
    }
  }
  // set up now
  setup_ = true;
}
#line 1 "./BSMWidthGenerator.cc"
// -*- C++ -*-
//
// BSMWidthGenerator.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BSMWidthGenerator class.
//

#include "BSMWidthGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/General/GeneralTwoBodyDecayer.h"

using namespace Herwig;

IBPtr BSMWidthGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr BSMWidthGenerator::fullclone() const {
  return new_ptr(*this);
}

void BSMWidthGenerator::persistentOutput(PersistentOStream & os) const {
  os << theModes;
}

void BSMWidthGenerator::persistentInput(PersistentIStream & is, int) {
  is >> theModes;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BSMWidthGenerator,GenericWidthGenerator>
describeHerwigBSMWidthGenerator("Herwig::BSMWidthGenerator", "Herwig.so");

void BSMWidthGenerator::Init() {

  static ClassDocumentation<BSMWidthGenerator> documentation
    ("A width generator for BSM particles.");

}

void BSMWidthGenerator::setupMode(tcDMPtr mode, tDecayIntegratorPtr decayer, 
				  unsigned int) {
  tcGeneralTwoBodyDecayerPtr dec = 
    dynamic_ptr_cast<tcGeneralTwoBodyDecayerPtr>(decayer);
  theModes.push_back(make_pair(mode, dec));
}

Energy BSMWidthGenerator::partial2BodyWidth(int iloc, Energy m0,
					    Energy m1, Energy m2) const {
  if( m0 < (m1 + m2) ) return ZERO;
  //need pointers to particles involved
  tcDMPtr dm = theModes[iloc].first;
  ParticleMSet::const_iterator pit = dm->products().begin();
  tcPDPtr parta = *pit;
  ++pit;
  tcPDPtr partb = *pit;
  int dummya(0);
  double dummyb(0.);
  if(theModes[iloc].second) {
    if( !theModes[iloc].second->twoBodyMEcode(*dm, dummya, dummyb) )
      swap(parta, partb);
    return theModes[iloc].second->partialWidth(make_pair(dm->parent(), m0),
					       make_pair(parta, m1),
					       make_pair(partb, m2));
  }
  else
    return theModes[iloc].first->brat()*theModes[iloc].first->parent()->width();
}
#line 1 "./PrototypeVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FourBodyDecayConstructor class.
//

#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/Exception.h"
#include "PrototypeVertex.h"
#include "NBodyDecayConstructorBase.h"

using namespace Herwig;

void PrototypeVertex::createPrototypes(tPDPtr inpart, VertexBasePtr vertex,
				       std::stack<PrototypeVertexPtr> & prototypes,
				       NBodyDecayConstructorBasePtr decayCon) {
  if(!vertex->isIncoming(inpart)) return;
  tPDPtr ccpart = inpart->CC() ? inpart->CC() : inpart;
  long id = ccpart->id();
  Energy2 scale = sqr(inpart->mass());
  for(unsigned int list=0;list<vertex->getNpoint();++list) {
    tPDVector decaylist = vertex->search(list, ccpart);
    tPDVector::size_type nd = decaylist.size();
    for( tPDVector::size_type i = 0; i < nd; i += vertex->getNpoint() ) {
      tPDVector pout(decaylist.begin()+i,
		     decaylist.begin()+i+vertex->getNpoint());
      OrderedVertices out;
      for(unsigned int ix=1;ix<pout.size();++ix) {
	if(pout[ix]->id() == id ) swap(pout[0], pout[ix]);
	out.insert(make_pair(pout[ix],PrototypeVertexPtr()));
      }
      // remove flavour changing vertices
      if(decayCon->removeFlavourChangingVertices()&&vertex->getNpoint()==3) {
	if((pout[1]->id()==ParticleID::gamma||
	    pout[1]->id()==ParticleID::g) &&
	   inpart!=pout[2]) continue;
	else if((pout[2]->id()==ParticleID::gamma||
		 pout[2]->id()==ParticleID::g) &&
		inpart!=pout[1]) continue;
      }
      // remove small vertices if needed
      if(decayCon->removeSmallVertices()) {
	if(vertex->getNpoint()==3)
	  vertex->setCoupling(scale,pout[0],pout[1],pout[2]);
	else if(vertex->getNpoint()==4)
	  vertex->setCoupling(scale,pout[0],pout[1],pout[2],pout[3]);
	if(abs(vertex->norm())<decayCon->minVertexNorm()) continue;
      }
      if(vertex->getNpoint()==3) {
	// remove radiation
	if((inpart==pout[1] && (pout[2]->id()==ParticleID::gamma||
				pout[2]->id()==ParticleID::g||
				pout[2]->id()==ParticleID::Z0)) 
	   ||
	   (inpart==pout[2] && (pout[1]->id()==ParticleID::gamma||
				pout[1]->id()==ParticleID::g||
				pout[1]->id()==ParticleID::Z0)))
	  continue;
      }
      prototypes.push(new_ptr(PrototypeVertex(inpart,out,vertex,
					      int(vertex->getNpoint())-1)));
    }
  }
}

PrototypeVertexPtr PrototypeVertex::replicateTree(PrototypeVertexPtr parent,
						  PrototypeVertexPtr oldChild,
						  PrototypeVertexPtr & newChild) {
  PrototypeVertexPtr newParent = 
    new_ptr(PrototypeVertex(parent->incoming,OrderedVertices(),
			    parent->vertex,parent->npart));
  for(OrderedVertices::const_iterator it = parent->outgoing.begin();
      it!=parent->outgoing.end();++it) {
    PrototypeVertexPtr child = it->second ? 
      replicateTree(it->second,oldChild,newChild) :
      PrototypeVertexPtr();
    newParent->outgoing.insert(make_pair(it->first,child));
    if(child) child->parent = newParent;
    if(it->second==oldChild) newChild=child;
  }
  if(oldChild==parent) newChild=newParent;
  return newParent;
}

void PrototypeVertex::expandPrototypes(PrototypeVertexPtr proto, VertexBasePtr vertex,
				       std::stack<PrototypeVertexPtr> & prototypes,
				       const set<PDPtr> & excluded,
				       NBodyDecayConstructorBasePtr decayCon) {
  for(OrderedVertices::const_iterator it = proto->outgoing.begin();
      it!=proto->outgoing.end();++it) {
    if(it->second) {
      expandPrototypes(it->second,vertex,prototypes,excluded,decayCon);
    }
    else {
      if(!vertex->isIncoming(it->first)) continue;
      if(excluded.find(it->first)!=excluded.end()) continue;
      if(it->first->CC() && 
	 excluded.find(it->first->CC())!=excluded.end()) continue;

      tcPDPtr ccpart = it->first->CC() ? it->first->CC() : it->first;
      long id = ccpart->id();

      PrototypeVertexPtr parent=proto;
      while(parent->parent) parent=parent->parent;
      Energy2 scale(sqr(parent->incoming->mass()));
      for(unsigned int il = 0; il < vertex->getNpoint(); ++il) {
	tPDVector decaylist = vertex->search(il,ccpart);
	tPDVector::size_type nd = decaylist.size();
	for( tPDVector::size_type i = 0; i < nd; i += vertex->getNpoint() ) {
	  tPDVector pout(decaylist.begin()+i,
			 decaylist.begin()+i+vertex->getNpoint());
	  OrderedVertices outgoing;
	  for(unsigned int iy=1;iy<pout.size();++iy) {
	    if(pout[iy]->id() == id ) swap(pout[0], pout[iy]);
	    outgoing.insert(make_pair(pout[iy],PrototypeVertexPtr()));
	  }
	  if(vertex->getNpoint()==3) {
	    if((it->first==pout[1] && (pout[2]->id()==ParticleID::gamma||
				       pout[2]->id()==ParticleID::g||
				       pout[2]->id()==ParticleID::Z0)) 
	       ||
	       (it->first==pout[2] && (pout[1]->id()==ParticleID::gamma||
				       pout[1]->id()==ParticleID::g||
				       pout[1]->id()==ParticleID::Z0)))
	      continue;
	    // remove weak decays of quarks other than top
	    if(StandardQCDPartonMatcher::Check(pout[0]->id()) &&
	       ((StandardQCDPartonMatcher::Check(pout[1]->id()) &&
		 abs(pout[2]->id())==ParticleID::Wplus)||
		(StandardQCDPartonMatcher::Check(pout[2]->id())&&
		 abs(pout[1]->id())==ParticleID::Wplus))) continue;
	    // remove weak decays of leptons
	    if((abs(pout[0]->id())>=11&&abs(pout[0]->id())<=16) &&
		(((abs(pout[1]->id())>=11&&abs(pout[1]->id())<=16) &&
		  abs(pout[2]->id())==ParticleID::Wplus)||
		((abs(pout[2]->id())>=11&&abs(pout[2]->id())<=16)&&
		 abs(pout[1]->id())==ParticleID::Wplus))) continue;
	    }
	  // remove flavour changing vertices
	  if(decayCon->removeFlavourChangingVertices()&&vertex->getNpoint()==3) {
	    if((pout[1]->id()==ParticleID::gamma||
		pout[1]->id()==ParticleID::g) &&
	       it->first!=pout[2]) continue;
	    else if((pout[2]->id()==ParticleID::gamma||
		     pout[2]->id()==ParticleID::g) &&
		    it->first!=pout[1]) continue;
	  }
	  // remove small vertices if needed
	  if(decayCon->removeSmallVertices()) {
	    if(vertex->getNpoint()==3)
	      vertex->setCoupling(scale,pout[0],pout[1],pout[2]);
	    else if(vertex->getNpoint()==4)
	      vertex->setCoupling(scale,pout[0],pout[1],pout[2],pout[3]);
	    if(abs(vertex->norm())<decayCon->minVertexNorm()) continue;
	  }
	  PrototypeVertexPtr newBranch = 
	    new_ptr(PrototypeVertex(it->first,
				    outgoing,vertex,int(vertex->getNpoint())-1));
	  PrototypeVertexPtr newChild;
	  PrototypeVertexPtr newVertex = replicateTree(parent,proto,newChild);
	  newBranch->parent = newChild;
	  OrderedVertices::iterator kt = newChild->outgoing.begin();
	  for(OrderedVertices::const_iterator jt = proto->outgoing.begin();
	      jt!=it;++jt,++kt) {;}
	  pair< tPDPtr, PrototypeVertexPtr > newPair = make_pair(kt->first,newBranch);
	  newChild->outgoing.erase(kt);
	  newChild->outgoing.insert(newPair);
	  newChild->incrementN(newBranch->npart-1);
	  prototypes.push(newVertex);
	}
      }
    }
  }
}

bool PrototypeVertex::canBeOnShell(unsigned int opt,Energy maxMass,bool first) {
  Energy outMass = outgoingMass();
  if(!first) { 
    bool in  = maxMass>incomingMass();
    bool out = incomingMass()>outMass;
    if(opt!=0) {
      if( in && ( out || opt==2 ) ) return true;
    }
    else if (incoming->width() == ZERO) {
      tPrototypeVertexPtr original = this;
      while(original->parent) {
	original=original->parent;
      };
      ostringstream name;
      name << original->incoming->PDGName() << " -> ";
      for(OrderedParticles::const_iterator it = original->outPart.begin();
	  it!= original->outPart.end();++it)
	name << (**it).PDGName() << " ";
      Throw<InitException>() 
	<< "Trying to include on-shell diagram in decay"
	<< name.str() << "including on-shell "
	<< incoming->PDGName() << " with zero width.\n"
	<< "You should make sure that the width for the intermediate is either"
	<< " read from an SLHA file or the intermediate is included in the "
	<< "DecayParticles list of the ModelGenerator.\n"
	<< Exception::runerror;
    }
  }
  else maxMass = incomingMass();
  // check the decay products
  for(OrderedVertices::const_iterator it = outgoing.begin();
      it!=outgoing.end();++it) {
    if(!it->second) continue;
    Energy newMax = maxMass-outMass+it->second->outgoingMass();
    if(it->second->canBeOnShell(opt,newMax,false)) {
      if(first) possibleOnShell = true;
      return true;
    }
  }
  return false;
}

bool PrototypeVertex::sameDecay(const PrototypeVertex & x) const {
  if(incoming != x.incoming) return false;
  if(outPart.empty())     setOutgoing();
  if(x.outPart.empty()) x.setOutgoing();
  OrderedParticles::const_iterator cit =   outPart.begin();
  OrderedParticles::const_iterator cjt = x.outPart.begin();
  if(x.npart!=npart) return false;
  for(;cit!=outPart.end();++cit,++cjt) {
    if(*cit!=*cjt) return false;
  }
  return true;
}

namespace {
  void constructTag(vector<unsigned int> & order,NBVertex & vertex,
		    const OrderedParticles & outgoing,
		    vector<bool> & matched) {
    for(list<pair<PDPtr,NBVertex> >::iterator it=vertex.vertices.begin();
	it!=vertex.vertices.end();++it) {
      // search down the tree
      if(it->second.vertex) {
	constructTag(order,it->second,outgoing,matched);
      }
      // identify this particle
      else {
	unsigned int ix=0;
	for(OrderedParticles::const_iterator pit=outgoing.begin();
	    pit!=outgoing.end();++pit) {
	  if(*pit==it->first&&!matched[ix]) {
	    matched[ix] = true;
	    order.push_back(ix+1);
	    break;
	  }
	  ++ix;
	}
      }
    }
  }
}

NBDiagram::NBDiagram(PrototypeVertexPtr proto) 
  : NBVertex(proto),
    colourFlow       (vector<CFPair>(1,make_pair(1,1.))),
    largeNcColourFlow(vector<CFPair>(1,make_pair(1,1.))) {
  if(!proto) return;
  // finally work out the channel and the ordering
  vector<bool> matched(outgoing.size(),false);
  for(list<pair<PDPtr,NBVertex> >::iterator it=vertices.begin();
      it!=vertices.end();++it) {
    // search down the tree
    if(it->second.vertex) {
      constructTag(channelType,it->second,outgoing,matched);
    }
    // identify this particle
    else {
      unsigned int ix=0;
      for(OrderedParticles::const_iterator pit=outgoing.begin();
	  pit!=outgoing.end();++pit) {
	if(*pit==it->first&&!matched[ix]) {
	  matched[ix] = true;
	  channelType.push_back(ix+1);
	  break;
	}
	++ix;
      }
    }
  }
}

NBVertex::NBVertex(PrototypeVertexPtr proto) { 
  if(!proto) return;
  incoming = proto->incoming;
  outgoing = proto->outPart;
  vertex   = proto->vertex;
  // create the vertices
  for(OrderedVertices::const_iterator it=proto->outgoing.begin();
      it!=proto->outgoing.end();++it) {
    vertices.push_back(make_pair(it->first,NBVertex(it->second)));
  }
  // now let's re-order so that branchings are at the end
  for(list<pair<PDPtr,NBVertex> >::iterator it=vertices.begin();
      it!=vertices.end();++it) {
    if(!it->second.incoming) continue;
    list<pair<PDPtr,NBVertex> >::iterator jt=it;
    for( ; jt!=vertices.end();++jt) {
      if(jt==it) continue;
      if(!jt->second.incoming) {
	break;
      }
    }
    if(jt!=vertices.end()) {
      list<pair<PDPtr,NBVertex> >::iterator kt = it;
      while(kt!=jt) {
	list<pair<PDPtr,NBVertex> >::iterator lt = kt;
	++lt;
	swap(*kt,*lt);
	kt=lt;
      }
    }
  }
}
#line 1 "./BSMModel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BSMModel class.
//

#include "BSMModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/MassGenerator.h"
#include "ThePEG/PDT/WidthGenerator.h"
#include "ThePEG/PDT/DecayMode.h"

#include "PrototypeVertex.h"

using namespace Herwig;

BSMModel::BSMModel() : decayFile_(), readDecays_(true),
		       topModesFromFile_(false),
		       tolerance_(1e-6), allowedToResetSMWidths_(false)
{}

void BSMModel::persistentOutput(PersistentOStream & os) const {
  os << decayFile_ << topModesFromFile_ << tolerance_ << allowedToResetSMWidths_;
}

void BSMModel::persistentInput(PersistentIStream & is, int) {
  is >> decayFile_ >> topModesFromFile_ >> tolerance_ >> allowedToResetSMWidths_;
}

DescribeAbstractClass<BSMModel,Herwig::StandardModel>
describeHerwigBSMModel("Herwig::BSMModel", "Herwig.so");

void BSMModel::Init() {

  static ClassDocumentation<BSMModel> documentation
    ("The BSMModel class provides a base class for BSM models including the"
     " features to read decays in the SLHA format");

  static Parameter<BSMModel,string> interfaceDecayFileName
    ("DecayFileName",
     "Name of the file from which to read decays in the SLHA format",
     &BSMModel::decayFile_, "",
     false, false);

  static Switch<BSMModel,bool> interfaceTopModes
    ("TopModes",
     "Whether ro use the Herwig SM top decays or those from the SLHA file",
     &BSMModel::topModesFromFile_, false, false, false);
  static SwitchOption interfaceTopModesFile
    (interfaceTopModes,
     "File",
     "Take the modes from the files",
     true);
  static SwitchOption interfaceTopModesHerwig
    (interfaceTopModes,
     "Herwig",
     "Use the SM ones", false);

  static Parameter<BSMModel,double> interfaceBRTolerance
    ("BRTolerance",
     "Tolerance for the sum of branching ratios to be difference from one.",
     &BSMModel::tolerance_, 1e-6, 1e-8, 0.01,
     false, false, Interface::limited);

  static Switch<BSMModel,bool> interfaceAllowedToResetSMWidths
    ("AllowedToResetSMWidths",
     "Whether or not the widths of SM particles can be reset via the SLHA file",
     &BSMModel::allowedToResetSMWidths_, false, false, false);
  static SwitchOption interfaceAllowedToResetSMWidthsNo
    (interfaceAllowedToResetSMWidths,
     "No",
     "Not allowed",
     false);
  static SwitchOption interfaceAllowedToResetSMWidthsYes
    (interfaceAllowedToResetSMWidths,
     "Yes",
     "Allowed",
     true);

}

void BSMModel::doinit() {
  StandardModel::doinit();
  // check if need to read decays
  if(decayFile()==""||!readDecays_) return;
  decayRead();
}

void BSMModel::decayRead() {
  // read decays
  CFileLineReader cfile;
  cfile.open(decayFile_);
  if( !cfile ) throw SetupException() 
		 << "BSMModel::doinit - An error occurred in opening the "
		 << "decay file \"" << decayFile_ << "\"."
		 << Exception::runerror;
  //Before reading the spectrum/decay files the SM higgs 
  //decay modes, mass and width generators need to be turned off.
  PDPtr h0 = getParticleData(ParticleID::h0);
  h0->widthGenerator(WidthGeneratorPtr());
  h0->massGenerator(MassGenPtr());
  h0->width(ZERO);
  h0->stable(true);
  DecaySet::const_iterator dit = h0->decayModes().begin();
  DecaySet::const_iterator dend = h0->decayModes().end();
  for( ; dit != dend; ++dit ) {
    generator()->preinitInterface(*dit, "BranchingRatio", "set", "0.");
    generator()->preinitInterface(*dit, "Active", "set", "No");
  }
  // if taking the top modes from the file
  // delete the SM stuff
  if(topModesFromFile_) {
    PDPtr top = getParticleData(ParticleID::t);
    top->widthGenerator(WidthGeneratorPtr());
    top->massGenerator(MassGenPtr());
    DecaySet::const_iterator dit = top->decayModes().begin();
    DecaySet::const_iterator dend = top->decayModes().end();
    for( ; dit != dend; ++dit ) {
      generator()->preinitInterface(*dit, "BranchingRatio", "set", "0.");
      generator()->preinitInterface(*dit, "Active", "set", "No");
    }
  }
  // read first line and check if this is a Les Houches event file
  cfile.readline();
  bool lesHouches = cfile.find("<LesHouchesEvents");
  bool reading = !lesHouches;
  if(lesHouches) cfile.readline();
  // function pointer for putting all characters to lower case.
  int (*pf)(int) = tolower;
  while (true) {
    string line = cfile.getline();
    // check for start of slha block in SLHA files
    if(lesHouches && !reading) {
      if(line.find("<slha")==0) reading = true;
      if(!cfile.readline()) break;
      continue;
    }
    // ignore comment lines
    if(line[0] == '#') {
      if(!cfile.readline()) break;
      continue;
    }
    // make everything lower case
    transform(line.begin(), line.end(), line.begin(), pf);
    // start of a block
    if(line.find("decay") == 0) {
      readDecay(cfile, line);
      if(!cfile) break;
      continue;
    }
    else if( lesHouches && line.find("</slha") == 0 ) {
      break;
    }
    if(!cfile.readline()) break;
  }
}

void BSMModel::readDecay(CFileLineReader & cfile, 
			 string decay) const{
  // extract parent PDG code and width
  long parent(0);
  Energy width(ZERO);
  istringstream iss(decay);
  string dummy;
  iss >> dummy >> parent >> iunit(width, GeV);
  PDPtr inpart = getBSMParticleData(parent);
  if(!topModesFromFile_&&abs(parent)==ParticleID::t) {
    cfile.readline();
    return;
  }
  if(!inpart)  throw SetupException() 
  		 << "BSMModel::readDecay() - "
  		 << "A ParticleData object with the PDG code "
  		 << parent << " does not exist. " 
  		 << Exception::runerror;
  // check this ain't a SM particle
  if(abs(parent)<=5||abs(parent)==23||abs(parent)==24||
     (abs(parent)>=11&&abs(parent)<=16)) {
    if(allowedToResetSMWidths_) {
      cerr << "BSMModel::readDecay() Resetting width of " 
	   << inpart->PDGName() << " using SLHA "
	   << "file,\nthis can affect parts of the Standard Model simulation and"
	   << " is strongly discouraged.\n";
      inpart->width(width);
      if( width > ZERO ) inpart->cTau(hbarc/width);
      inpart->widthCut(5.*width);
    }
    else {
      cerr << "BSMModel::readDecay() You have tried to Reset the width of " 
	   << inpart->PDGName() << " using an SLHA "
	   << "file,\nthis can affect parts of the Standard Model simulation and"
	   << " is not allowed by default.\n If you really want to be this stupid"
	   << " set AllowedToResetSMWidths to Yes for this model.\n";
    }
  }
  else {
    inpart->width(width);
    if( width > ZERO ) inpart->cTau(hbarc/width);
    inpart->widthCut(5.*width);
  }
  
  Energy inMass = inpart->mass();
  string prefix(inpart->name() + "->");
  double brsum(0.);
  unsigned int nmode = 0;
  while(cfile.readline()) {
    string line = cfile.getline();
    line = StringUtils::stripws(line);
    // skip comments
    if(line[0] == '#') continue;
    // reached the end
    if( line[0] == 'B' || line[0] == 'b' ||
	line[0] == 'D' || line[0] == 'd' ||
	line[0] == '<' ) {
      cfile.resetline();
      break;
    }
    // read the mode
    // get the branching ratio and no of decay products
    istringstream is(line);
    double brat(0.);
    unsigned int nda(0),npr(0);
    is >> brat >> nda;
    vector<tcPDPtr> products,bosons;
    Energy mout(ZERO),moutnoWZ(ZERO);
    string tag = prefix;
    multiset<tcPDPtr,ParticleOrdering> outgoing;
    int charge = -inpart->iCharge();
    while( true ) {
      long t;
      is >> t;
      if( is.fail() ) break; 
      if( t == abs(parent) ) {
  	throw SetupException() 
  	  << "An error occurred while read a decay of the " 
  	  << inpart->PDGName() << ". One of its products has the same PDG code "
  	  << "as the parent particle. Please check the SLHA file.\n"
  	  << Exception::runerror;
      }
      tcPDPtr p = getBSMParticleData(t);
      if( !p ) {
  	throw SetupException()
  	  << "BSMModel::readDecay() - An unknown PDG code has been encounterd "
  	  << "while reading a decay mode. ID: " << t
  	  << Exception::runerror;
      }
      charge += p->iCharge();
      ++npr;
      outgoing.insert(p);
      Energy mass =  p->mass();
      mout += mass;
      if(abs(p->id())==ParticleID::Wplus||p->id()==ParticleID::Z0) {
  	bosons.push_back(p);
      }
      else {
  	products.push_back(p);
  	moutnoWZ += mass;
      }
    }
    if( npr != nda ) {
      throw SetupException()
 	<< "BSMModel::readDecay - While reading a decay of the " 
 	<< inpart->PDGName() << " from an SLHA file, an inconsistency "
 	<< "between the number of decay products and the value in "
 	<< "the 'NDA' column was found. Please check if the spectrum "
 	<< "file is correct.\n"
 	<< Exception::warning;
    }

    // must be at least two decay products
    if(npr<=1) continue;

    // create the tag
    for(multiset<tcPDPtr,ParticleOrdering>::iterator it=outgoing.begin();
	it!=outgoing.end();++it)
      tag += (**it).name() + ",";
    tag.replace(tag.size() - 1, 1, ";");
    if(charge!=0) {
      cerr << "BSMModel::readDecay() "
	   << "Decay mode " << tag << " read from SLHA file does not conserve charge,"
	   << "\nare you really sure you want to do this?\n";
    }
    ++nmode;
    if(nmode==1) {
      if(abs(parent)<=5||abs(parent)==23||abs(parent)==24||
	 (abs(parent)>=11&&abs(parent)<=16)) {
	cerr << "BSMModel::readDecay() Resetting the decays of " 
	     << inpart->PDGName() << " using SLHA "
	     << "file,\nthis can affect parts of the Standard Model simulation,"
	     << " give unexpected results and"
	     << " is strongly discouraged.\n";
	cerr << "Switching off all the internal modes so only those from the SLHA file "
	     << "are used, this may have unintended consequences\n";
	for(DecaySet::iterator it=inpart->decayModes().begin();
	    it!=inpart->decayModes().end();++it) {
	  generator()->preinitInterface(*it, "Active", "set", "No");
	}
	if(inpart->CC()) {
	  for(DecaySet::iterator it=inpart->CC()->decayModes().begin();
	      it!=inpart->CC()->decayModes().end();++it) {
	    generator()->preinitInterface(*it, "Active", "set", "No");
	  }
	}
      }
    }
    // normal option
    if(mout<=inMass) {
      inpart->stable(false);
      brsum += brat;
      createDecayMode(tag, brat);
    }
    // no possible off-shell gauge bosons throw it away
    else if(bosons.empty() || bosons.size()>2 ||
	    moutnoWZ>=inMass) {
      cerr << "BSMModel::readDecay() "
	   << "The decay " << tag << " cannot proceed for on-shell "
	   << "particles, skipping it.\n";
    }
    else {
      Energy maxMass = inMass - moutnoWZ;
      string newTag = prefix;
      for(unsigned int ix=0;ix<products.size();++ix)
	newTag += products[ix]->name() + ",";
      if(bosons.size()==1) {
	cerr << "BSMModel::readDecay() "
	     << "The decay " << tag << " cannot proceed for on-shell\n"
	     << "particles, replacing gauge boson with its decay products\n";
	vector<pair<double,string> > modes = 
	  createWZDecayModes(newTag,brat,bosons[0],maxMass);
	for(unsigned int ix=0;ix<modes.size();++ix) {
	  modes[ix].second.replace(modes[ix].second.size() - 1, 1, ";");
	  createDecayMode(modes[ix].second,modes[ix].first);
	  brsum += modes[ix].first;
	}
      }
      else if(bosons.size()==2) {
	bool identical = bosons[0]->id()==bosons[1]->id();
	if(maxMass>bosons[0]->mass()&&maxMass>bosons[1]->mass()) {
	  cerr << "BSMModel::readDecay() "
	       << "The decay " << tag << " cannot proceed for on-shell\n"
	       << "particles, replacing one of the gauge bosons"
	       << " with its decay products\n";
	  unsigned int imax = identical ? 1 : 2;
	  if(imax==2) brat *= 0.5;
	  for(unsigned int ix=0;ix<imax;++ix) {
	    string newTag2 = newTag+bosons[ix]->name()+',';
	    unsigned int iother = ix==0 ? 1 : 0;
	    vector<pair<double,string> > modes = 
	      createWZDecayModes(newTag2,brat,bosons[iother],maxMass);
	    for(unsigned int ix=0;ix<modes.size();++ix) {
	      modes[ix].second.replace(modes[ix].second.size() - 1, 1, ";");
	      createDecayMode(modes[ix].second,modes[ix].first);
	      brsum += modes[ix].first;
	    }
	  }
	}
	else {
	  cerr << "BSMModel::readDecay() "
	       << "The decay " << tag << " cannot proceed for on-shell\n"
	       << "particles, and has too many off-shell gauge bosons,"
	       << " skipping it.\n";
	}
      }
      else {
	cerr << "BSMModel::readDecay() "
	     << "The decay " << tag << " cannot proceed for on-shell\n"
	     << "particles, and has too many outgoing gauge bosons skipping it.\n";
      }
    }
  }
  if( abs(brsum - 1.) > tolerance_ && nmode!=0 ) {
    cerr << "Warning: The total branching ratio for " << inpart->PDGName()
  	 << " from the spectrum file does not sum to 1.\nThe branching fractions"
  	 << " will be rescaled. Difference from 1 is " 
	 << abs(brsum - 1.) << "\n";
  }
  if(nmode>0) {
    inpart->update();
    inpart->reset();
    if(inpart->CC()) {
      inpart->CC()->update();
      inpart->CC()->reset();
    }
    if(inpart->massGenerator())
      inpart->massGenerator()->reset();
    if(inpart->widthGenerator())
      inpart->widthGenerator()->reset();
  }
}

void BSMModel::createDecayMode(string tag, double brat) const {
  tDMPtr dm = generator()->findDecayMode(tag);
  if(!dm) dm = generator()->preinitCreateDecayMode(tag);
  generator()->preinitInterface(dm, "Active", "set", "Yes");
  generator()->preinitInterface(dm, "Decayer", "set","/Herwig/Decays/Mambo");
  ostringstream brf;
  brf << setprecision(13)<< brat;
  generator()->preinitInterface(dm, "BranchingRatio","set", brf.str());
  if(dm->CC()) {
    generator()->preinitInterface(dm->CC(), "Active", "set", "Yes");
    generator()->preinitInterface(dm->CC(), "Decayer", "set","/Herwig/Decays/Mambo");
    generator()->preinitInterface(dm->CC(), "BranchingRatio","set", brf.str());
  }
}


vector<pair<double,string> > 
BSMModel::createWZDecayModes(string tag, double brat,
			     tcPDPtr boson, Energy maxMass) const {
  vector<pair<double,string> > modes;
  double sum(0.);
  for(DecaySet::const_iterator dit=boson->decayModes().begin();
      dit!=boson->decayModes().end();++dit) {
    tcDMPtr mode = *dit;
    if(!mode->on()) continue;
    string extra;
    Energy outMass(ZERO);
    for(ParticleMSet::const_iterator pit=mode->products().begin();
	pit!=mode->products().end();++pit) {
      extra += (**pit).name() + ",";
      outMass += (**pit).mass();
    }
    if(outMass<maxMass) {
      sum += mode->brat();
      modes.push_back(make_pair(mode->brat(),tag+extra));
    }
  }
  for(unsigned int ix=0;ix<modes.size();++ix) 
    modes[ix].first *= brat/sum;
  return modes;
}
