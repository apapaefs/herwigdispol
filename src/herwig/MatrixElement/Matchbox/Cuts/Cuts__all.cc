#line 1 "./IdentifiedParticleCut.cc"
// -*- C++ -*-
//
// IdentifiedParticleCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IdentifiedParticleCut class.
//

#include "IdentifiedParticleCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IdentifiedParticleCut::IdentifiedParticleCut() 
  : thePtMin(0.*GeV), thePtMax(Constants::MaxEnergy) {}

IBPtr IdentifiedParticleCut::clone() const {
  return new_ptr(*this);
}

IBPtr IdentifiedParticleCut::fullclone() const {
  return new_ptr(*this);
}

bool IdentifiedParticleCut::passCuts(tcCutsPtr parent,
				     tcPDPtr ptype, LorentzMomentum p) const {

  if ( !matcher()->check(*ptype) ||
       ( thePtMin == ZERO && thePtMax == Constants::MaxEnergy &&
	 theYRanges.empty() ) )
    return true;

  double weight = 1.0;

  if ( !parent->isInside<CutTypes::Momentum>(p.perp(),ptMin(),ptMax(),weight) ) {
    parent->lastCutWeight(0.0);
    return false;
  }

  double y = p.rapidity() + parent->currentYHat();
  for ( vector<pair<double,double> >::const_iterator dy = yRanges().begin();
	dy != yRanges().end(); ++dy ) {
    if ( !parent->isInside<CutTypes::Rapidity>(y,dy->first,dy->second,weight) ) {
      parent->lastCutWeight(0.0);
      return false;
    }
  }

  parent->lastCutWeight(weight);
  return true;

}

void IdentifiedParticleCut::describe() const {

  CurrentGenerator::log()
    << "IdentifiedParticleCut '" << name() << "' matching "
    << "'" << matcher()->name() << "'";
  CurrentGenerator::log() << " within:\n";

  CurrentGenerator::log() 
    << "pt  = " << ptMin()/GeV << " .. " << ptMax()/GeV << " GeV\n";

  for ( vector<pair<double,double> >::const_iterator r = yRanges().begin();
	r != yRanges().end(); ++r ) {
    CurrentGenerator::log() << "y   = " << r->first << " .. " << r->second << "\n";
  }

}

string IdentifiedParticleCut::doYRange(string in) {
  istringstream ins(in);
  double first, second;
  ins >> first >> second;
  if ( first > second )
    swap(first,second);
  theYRanges.push_back(make_pair(first,second));
  return "";
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IdentifiedParticleCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePtMin,GeV) << ounit(thePtMax,GeV)
     << theYRanges << theMatcher;
}

void IdentifiedParticleCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePtMin,GeV) >> iunit(thePtMax,GeV)
     >> theYRanges >> theMatcher;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IdentifiedParticleCut,OneCutBase>
  describeHerwigIdentifiedParticleCut("Herwig::IdentifiedParticleCut", "HwMatchboxCuts.so");

void IdentifiedParticleCut::Init() {

  static ClassDocumentation<IdentifiedParticleCut> documentation
    ("IdentifiedParticleCut implements cuts on single momenta.");

  static Parameter<IdentifiedParticleCut,Energy> interfacePtMin
    ("PtMin",
     "The minimum pt required.",
     &IdentifiedParticleCut::thePtMin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<IdentifiedParticleCut,Energy> interfacePtMax
    ("PtMax",
     "The maximum pt allowed.",
     &IdentifiedParticleCut::thePtMax, GeV, Constants::MaxEnergy, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Command<IdentifiedParticleCut> interfaceYRange
    ("YRange",
     "Insert a rapidity range.",
     &IdentifiedParticleCut::doYRange, false);

  static Reference<IdentifiedParticleCut,MatcherBase> interfaceMatcher
    ("Matcher",
     "A matcher for particles to cut on.",
     &IdentifiedParticleCut::theMatcher, false, false, true, false, false);

}

#line 1 "./MatchboxDeltaRCut.cc"
// -*- C++ -*-
//
// MatchboxDeltaRCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxDeltaRCut class.
//

#include "MatchboxDeltaRCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxDeltaRCut::MatchboxDeltaRCut() 
  : theDeltaRMin(0.0), theDeltaRMax(Constants::MaxRapidity), 
    theDeltaYMin(0.0), theDeltaYMax(Constants::MaxRapidity),
    theDeltaPhiMin(0.0), theDeltaPhiMax(2.0*Constants::pi) {}

IBPtr MatchboxDeltaRCut::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxDeltaRCut::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxDeltaRCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
  if ( !match ||
       (theDeltaRMin == 0.0 && theDeltaRMax == Constants::MaxRapidity && 
	theDeltaYMin == 0.0 && theDeltaYMax == Constants::MaxRapidity &&
	theDeltaPhiMin == 0.0 && theDeltaPhiMax == 2.0*Constants::pi) ) return true;
  if ( inci || incj ) return true;

  double weight = 1.0;

  double dY = abs(pi.rapidity() - pj.rapidity());
  double dPhi = abs(pi.phi() - pj.phi());
  if ( dPhi > Constants::pi ) dPhi = 2.0*Constants::pi - dPhi;
  double dR = sqrt(sqr(dY) + sqr(dPhi));
  if ( !parent->isInside<CutTypes::Rapidity>(dY,deltaYMin(),deltaYMax(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }
  if ( !parent->isInside<CutTypes::Azimuth>(dPhi,deltaPhiMin(),deltaPhiMax(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }
  if ( !parent->isInside<CutTypes::Rapidity>(dR,deltaRMin(),deltaRMax(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

void MatchboxDeltaRCut::describe() const {

  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching distances between: '"
    << theFirstMatcher->name() << "' and '"
    << theSecondMatcher->name() << "':\n"
    << "DeltaRMin = " << theDeltaRMin << " \n"
    << "DeltaRMax = " << theDeltaRMax << " \n"
    << "DeltaPhiMin = " << theDeltaPhiMin << " \n"
    << "DeltaPhiMax = " << theDeltaPhiMax << " \n"
    << "DeltaYMin = " << theDeltaYMin << " \n"
    << "DeltaYMax = " << theDeltaYMax << " \n\n";

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxDeltaRCut::persistentOutput(PersistentOStream & os) const {
  os << theDeltaYMin << theDeltaYMax 
     << theDeltaPhiMin << theDeltaPhiMax 
     << theDeltaRMin << theDeltaRMax 
     << theFirstMatcher << theSecondMatcher;
}

void MatchboxDeltaRCut::persistentInput(PersistentIStream & is, int) {
  is >> theDeltaYMin >> theDeltaYMax 
     >> theDeltaPhiMin >> theDeltaPhiMax 
     >> theDeltaRMin >> theDeltaRMax 
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxDeltaRCut,TwoCutBase>
  describeHerwigMatchboxDeltaRCut("Herwig::MatchboxDeltaRCut", "HwMatchboxCuts.so");

void MatchboxDeltaRCut::Init() {

  static ClassDocumentation<MatchboxDeltaRCut> documentation
    ("This class implements cuts on legoplot, rapidity and azimuthal separation, "
     "i.e. on the \\f$\\Delta R\\f$-measure and on \\f$\\Delta Y\\f$ and \\f$\\Delta \\phi\\f$. "
     "By default the cuts are only applied to coloured particles, but "
     "may optionally be applied to all particle types. ");

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaRMin
    ("DeltaRMin",
     "The minimum allowed for the legoplot distance "
     "\\f$\\Delta R_{ij}=\\sqrt{\\Delta \\phi_{ij}^2+\\Delta Y_{ij}^2}\\f$ ",
     &MatchboxDeltaRCut::theDeltaRMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaRMax
    ("DeltaRMax",
     "The maximum allowed for the legoplot distance "
     "\\f$\\Delta R_{ij}=\\sqrt{\\Delta \\phi_{ij}^2+\\Delta Y_{ij}^2}\\f$ ",
     &MatchboxDeltaRCut::theDeltaRMax, Constants::MaxRapidity, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaPhiMin
    ("DeltaPhiMin",
     "The minimum allowed for the azimuthal separation "
     "\\f$\\Delta \\phi_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaPhiMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaPhiMax
    ("DeltaPhiMax",
     "The maximum allowed for the azimuthal separation "
     "\\f$\\Delta \\phi_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaPhiMax, 2.0*Constants::pi, 0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaYMin
    ("DeltaYMin",
     "The minimum allowed for the rapidity separation "
     "\\f$\\Delta Y_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaYMin, 0.0, 0.0, 0,
     false, false, Interface::lowerlim);

  static Parameter<MatchboxDeltaRCut,double> interfaceDeltaYMax
    ("DeltaYMax",
     "The maximum allowed for the rapidity separation "
     "\\f$\\Delta Y_{ij}\\f$ ",
     &MatchboxDeltaRCut::theDeltaYMax, Constants::MaxRapidity, 0, 0,
     false, false, Interface::lowerlim);

  static Reference<MatchboxDeltaRCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "If non-null only particles matching this object will be affected "
     "by the cut. ",
     &MatchboxDeltaRCut::theFirstMatcher, true, false, true, true, false);
//      &MatchboxDeltaRCut::theFirstMatcher, false, false, true, false, false);

  static Reference<MatchboxDeltaRCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "If non-null only particles matching this object will be affected "
     "by the cut. ",
     &MatchboxDeltaRCut::theSecondMatcher, true, false, true, true, false);
//      &MatchboxDeltaRCut::theSecondMatcher, false, false, true, false, false);

}

#line 1 "./MissingPtCut.cc"
// -*- C++ -*-
//
// MissingPtCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MissingPtCut class.
//

#include "MissingPtCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MissingPtCut::MissingPtCut() 
  : thePtMissMin(0.*GeV), thePtMissMax(Constants::MaxEnergy) {}

IBPtr MissingPtCut::clone() const {
  return new_ptr(*this);
}

IBPtr MissingPtCut::fullclone() const {
  return new_ptr(*this);
}

bool MissingPtCut::passCuts(tcCutsPtr parent, const tcPDVector & ptype, 
                            const vector<LorentzMomentum> & p) const {

  if ( thePtMissMin == ZERO && thePtMissMax == Constants::MaxEnergy )
    return true;

  // Energy ptMissSum = 0.0*GeV;
  LorentzMomentum momentumMissSum;
  bool nonu = true;

  for ( int i = 0, N = ptype.size(); i < N; ++i ) {

    if ( invisibleParticles().size() == 0 ) {
      if ( matcher()->check(*ptype[i]) ) {
        // ptMissSum = ptMissSum + p[i].perp();
        momentumMissSum = momentumMissSum + p[i];
        nonu = false;
      }
    }
    else if ( invisibleParticles().size() != 0 ) {
      for ( vector<int>::const_iterator iID = invisibleParticles().begin(); iID != invisibleParticles().end(); ++iID ) {
        int iInt = *iID;
        if ( abs(ptype[i]->id())==iInt ) {
          // ptMissSum = ptMissSum + p[i].perp();
          momentumMissSum = momentumMissSum + p[i];
          nonu = false;
        }
      }
    }

  }

  if ( nonu ) return true;

  Energy ptMiss = momentumMissSum.perp();  

  double weight = 1.0;

  // if ( !parent->isInside<CutTypes::Momentum>(ptMissSum,ptMissMin(),ptMissMax(),weight) ) {
  if ( !parent->isInside<CutTypes::Momentum>(ptMiss,ptMissMin(),ptMissMax(),weight) ) {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

string MissingPtCut::doInvisibleParticles(string in) {
  istringstream ins(in);
  int first;
  ins >> first;
  theInvisibleParticles.push_back(first);
  return "";
}

void MissingPtCut::describe() const {

  CurrentGenerator::log()
    << "MissingPtCut '" << name() << "' matching "
    << "'" << matcher()->name() << "'";
  CurrentGenerator::log() << " within:\n";

  CurrentGenerator::log() 
    << "ptMiss  = " << ptMissMin()/GeV << " .. " << ptMissMax()/GeV << " GeV\n";

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MissingPtCut::persistentOutput(PersistentOStream & os) const {
//   os << ounit(thePtMissMin,GeV) << ounit(thePtMissMax,GeV);
  os << ounit(thePtMissMin,GeV) << ounit(thePtMissMax,GeV) 
     << theInvisibleParticles << theMatcher;
}

void MissingPtCut::persistentInput(PersistentIStream & is, int) {
//   is >> iunit(thePtMissMin,GeV) >> iunit(thePtMissMax,GeV);
  is >> iunit(thePtMissMin,GeV) >> iunit(thePtMissMax,GeV)
     >> theInvisibleParticles >> theMatcher;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MissingPtCut,MultiCutBase>
  describeHerwigMissingPtCut("Herwig::MissingPtCut", "HwMatchboxCuts.so");

void MissingPtCut::Init() {

  static ClassDocumentation<MissingPtCut> documentation
//    ("MissingPtCut implements a cut on the missing transverse momentum "
//     "of a set of outgoing particles, i.e. for now the total transverse momentum "
//     "of all outgoing neutrinos in an event.");
    ("MissingPtCut implements a cut on the transverse momentum of the four-momentum "
     "sum of a set of outgoing particles that cannot be detected. By default the three "
     "standard model neutrinos are considered. If at least one undetectable particle "
     "is specified through the InvisibleParticles interface, the default choice is "
     "nullified.");

  static Command<MissingPtCut> interfaceInvisibleParticles
    ("InvisibleParticles",
     "Insert the PDG code of a particle that cannot be detected. If no particle " 
     "is inserted at all, the three standard model neutrinos are considered by "
     "default. If at least one particle is inserted, the default choice is nullified.",
     &MissingPtCut::doInvisibleParticles, false);

  static Parameter<MissingPtCut,Energy> interfacePtMissMin
    ("PtMissMin",
     "The minimum missing pt required.",
     &MissingPtCut::thePtMissMin, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<MissingPtCut,Energy> interfacePtMissMax
    ("PtMissMax",
     "The maximum missing pt allowed.",
     &MissingPtCut::thePtMissMax, GeV, Constants::MaxEnergy, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Reference<MissingPtCut,MatcherBase> interfaceMatcher
    ("Matcher",
     "A matcher for particles to cut on.",
     &MissingPtCut::theMatcher, false, false, true, false, false);

}

#line 1 "./FrixionePhotonSeparationCut.cc"
// -*- C++ -*-
//
// FrixionePhotonSeparationCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FrixionePhotonSeparationCut class.
//

#include "FrixionePhotonSeparationCut.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace Herwig;

struct FrixionePartonInfo {
  double DeltaR;
  Energy pT;
  double f;
};

void FrixionePhotonSeparationCut::describe() const {
  CurrentGenerator::log() 
    << fullName()
    << " matching unresolved particles from '"
    << matcher()->name() << "':\n"
    << "DeltaZero = " << theDeltaZero << " \n"
    << "Exponent n = " << theExponentn << " \n"
    << "Efficiency = " << theEfficiency << " \n"
    << "Cut Type = " << theCutType << " \n\n";
}

IBPtr FrixionePhotonSeparationCut::clone() const {
  return new_ptr(*this);
}

IBPtr FrixionePhotonSeparationCut::fullclone() const {
  return new_ptr(*this);
}

bool FrixionePhotonSeparationCut::passCuts(tcCutsPtr parent, const tcPDVector & ptype,
					   const vector<LorentzMomentum> & p) const {
  if ( theDeltaZero <= 0 ) return true;

  double weight = 1.0;

  for ( int i = 0, N = ptype.size(); i < N; ++i ) 
    if ( ptype[i]->id() == ParticleID::gamma ) {
      vector<FrixionePartonInfo> partonvec;
      for ( int j = 0, M = ptype.size(); j < M; ++j ) {

	if ( !matcher()->check(*ptype[j]) ) continue;

        FrixionePartonInfo finfo;
        double dphi = abs(p[i].phi() - p[j].phi());
        if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
        finfo.DeltaR = sqrt(sqr(p[i].eta() - p[j].eta()) + sqr(dphi));
	if ( finfo.DeltaR < theDeltaZero ){
          finfo.pT = p[j].perp();
          finfo.f = pow((1-cos(finfo.DeltaR))/(1-cos(theDeltaZero)),theExponentn);
          partonvec.push_back(finfo);
	}

      }

      for ( unsigned int j = 0; j < partonvec.size(); ++j ) {
        Energy chidelta=ZERO;
	if (theCutType == 1) {
          for ( unsigned int k = 0; k < partonvec.size(); ++k ) 
	    if ( partonvec[k].DeltaR <= partonvec[j].DeltaR ) chidelta += partonvec[k].pT;
	}
	else if (theCutType == 2) {
          chidelta = partonvec[j].pT;
	}
        if ( !parent->isLessThan<CutTypes::Momentum>(chidelta,p[i].perp() * theEfficiency * partonvec[j].f,weight) ) {
	  parent->lastCutWeight(0.0);
	  return false;
	}
      }

    }


  parent->lastCutWeight(weight);
  return true;

}

void FrixionePhotonSeparationCut::persistentOutput(PersistentOStream & os) const {
  os << theDeltaZero << theExponentn << theEfficiency << theCutType << theMatcher;
}

void FrixionePhotonSeparationCut::persistentInput(PersistentIStream & is, int) {
  is >> theDeltaZero >> theExponentn >> theEfficiency >> theCutType >> theMatcher;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FrixionePhotonSeparationCut,MultiCutBase>
describeHerwigFrixionePhotonSeparationCut("Herwig::FrixionePhotonSeparationCut", "HwMatchboxCuts.so");

void FrixionePhotonSeparationCut::Init() {

  static ClassDocumentation<FrixionePhotonSeparationCut> documentation
    ("This class implements a separation criterium a la Frixione between "
     "final-state partons and photons.");

   static Parameter<FrixionePhotonSeparationCut,double> interfaceDeltaZero
   ("DeltaZero",
     "The maximal legoplot separation up to which partons are included in the criterium ",
     &FrixionePhotonSeparationCut::theDeltaZero, 0.7, 0.0, 10.0,
     false, false, Interface::limited);

   static Parameter<FrixionePhotonSeparationCut,double> interfaceExponentn
   ("Exponentn",
     "The exponent n of the algorithm ",
     &FrixionePhotonSeparationCut::theExponentn, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

   static Parameter<FrixionePhotonSeparationCut,double> interfaceEfficiency
   ("Efficiency",
     "The efficiency epsilon of the algorithm ",
     &FrixionePhotonSeparationCut::theEfficiency, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Switch<FrixionePhotonSeparationCut,int> interfaceCutType
    ("CutType",
     "Switch for controlling which definition of Frixione cut is used",
     &FrixionePhotonSeparationCut::theCutType, 1, false, false);
  static SwitchOption interfaceCutTypeVBFNLO
    (interfaceCutType,
     "VBFNLO",
     "Switch to Frixione cut a la VBFNLO",
     1);
  static SwitchOption interfaceCutTypeMCFM
    (interfaceCutType,
     "MCFM",
     "Switch to Frixione cut a la MCFM",
     2);

  static Reference<FrixionePhotonSeparationCut,MatcherBase> interfaceMatcher
    ("UnresolvedMatcher",
     "A matcher for particles to isolate on.",
     &FrixionePhotonSeparationCut::theMatcher, false, false, true, false, false);

  interfaceDeltaZero.setHasDefault(false);
  interfaceExponentn.setHasDefault(false);
  interfaceEfficiency.setHasDefault(false);
  interfaceCutType.setHasDefault(false);

}

#line 1 "./InvariantMassCut.cc"
// -*- C++ -*-
//
// InvariantMassCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InvariantMassCut class.
//

#include "InvariantMassCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace Herwig;

void InvariantMassCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching distances between: '"
    << theFirstMatcher->name() << "' and '"
    << theSecondMatcher->name() << "':\n"
    << "M = " << theMinMass/GeV << " .. " << theMaxMass/GeV << " GeV\n"
    << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
    << "opposite sign only = " << (theOppositeSignOnly?"Yes":"No") << " \n\n";
}

IBPtr InvariantMassCut::clone() const {
  return new_ptr(*this);
}

IBPtr InvariantMassCut::fullclone() const {
  return new_ptr(*this);
}

bool InvariantMassCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
  if ( !match ||
       ( theMinMass == ZERO && theMaxMass == Constants::MaxEnergy ) ) return true;
  if ( inci || incj ) return true;  

  if ( sameFlavourOnly() || oppositeSignOnly() ) {

    int fam1 = family(pitype->id());
    int fam2 = family(pjtype->id());

    if ( fam1 && fam2 ) {
      if ( sameFlavourOnly() && ( abs(fam1) != abs(fam2) ) ) return true;
      if ( oppositeSignOnly() && ( fam1*fam2 > 0 ) ) return true;
    }

  }


  double weight = 1.0;

  Energy minv = (pi+pj).m();

  if ( minv < ZERO ||
       !parent->isInside<CutTypes::Energy>(minv,minMass(),maxMass(),weight) ) {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

int InvariantMassCut::family(long id) const {

  int sign = (id>0)?-1:1;

  switch ( id ) {
    case ParticleID::u:
    case ParticleID::ubar:
    case ParticleID::d:
    case ParticleID::dbar:
      return 1*sign; break;
    case ParticleID::c:
    case ParticleID::cbar:
    case ParticleID::s:
    case ParticleID::sbar:
      return 2*sign; break;
    case ParticleID::t:
    case ParticleID::tbar:
    case ParticleID::b:
    case ParticleID::bbar:
      return 3*sign; break;
    case ParticleID::eminus:
    case ParticleID::eplus:
    case ParticleID::nu_e:
    case ParticleID::nu_ebar:
      return 11*sign; break;
    case ParticleID::muminus:
    case ParticleID::muplus:
    case ParticleID::nu_mu:
    case ParticleID::nu_mubar:
      return 12*sign; break;
    case ParticleID::tauminus:
    case ParticleID::tauplus:
    case ParticleID::nu_tau:
    case ParticleID::nu_taubar:
      return 13*sign; break;
  }

  return 0;

}

void InvariantMassCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinMass,GeV) << ounit(theMaxMass,GeV) 
     << theSameFlavourOnly << theOppositeSignOnly
     << theFirstMatcher << theSecondMatcher;
}

void InvariantMassCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinMass,GeV) >> iunit(theMaxMass,GeV) 
     >> theSameFlavourOnly >> theOppositeSignOnly
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<InvariantMassCut,TwoCutBase>
  describeHerwigInvariantMassCut("Herwig::InvariantMassCut", "HwMatchboxCuts.so");

void InvariantMassCut::Init() {

  static ClassDocumentation<InvariantMassCut> documentation
    ("This class implements an invariant mass cut between "
     "final-state particles.");

  static Parameter<InvariantMassCut,Energy> interfaceMinMass
   ("MinMass",
     "The minimal allowed invariant mass ",
     &InvariantMassCut::theMinMass, GeV, 0*GeV, 0*GeV, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Parameter<InvariantMassCut,Energy> interfaceMaxMass
   ("MaxMass",
     "The maximal allowed invariant mass ",
     &InvariantMassCut::theMaxMass, GeV, Constants::MaxEnergy, 0*GeV, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Switch<InvariantMassCut,bool> interfaceSameFlavourOnly
   ("SameFlavourOnly",
     "Whether cut works on fermion pairs of the same flavour only ",
     &InvariantMassCut::theSameFlavourOnly, true, false, false);
  static SwitchOption interfaceSameFlavourOnlyYes
    (interfaceSameFlavourOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceSameFlavourOnlyNo
    (interfaceSameFlavourOnly,
     "No",
     "No",
     false);

  static Switch<InvariantMassCut,bool> interfaceOppositeSignOnly
   ("OppositeSignOnly",
     "Whether cut works on fermion pairs of opposite sign only ",
     &InvariantMassCut::theOppositeSignOnly, true, false, false);
  static SwitchOption interfaceOppositeSignOnlyYes
    (interfaceOppositeSignOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceOppositeSignOnlyNo
    (interfaceOppositeSignOnly,
     "No",
     "No",
     false);

  static Reference<InvariantMassCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &InvariantMassCut::theFirstMatcher, true, false, true, true, false);

  static Reference<InvariantMassCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &InvariantMassCut::theSecondMatcher, true, false, true, true, false);

}

#line 1 "./PairPtCut.cc"
// -*- C++ -*-
//
// PairPtCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PairPtCut class.
//

#include "PairPtCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

using namespace Herwig;

void PairPtCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching distances between: '"
    << theFirstMatcher->name() << "' and '"
    << theSecondMatcher->name() << "':\n"
    << "pT = " << theMinPt/GeV << " .. " << theMaxPt/GeV << " GeV\n"
    << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
    << "opposite sign only = " << (theOppositeSignOnly?"Yes":"No") << " \n\n";
}

IBPtr PairPtCut::clone() const {
  return new_ptr(*this);
}

IBPtr PairPtCut::fullclone() const {
  return new_ptr(*this);
}

bool PairPtCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
  if ( !match ||
       ( theMinPt == ZERO && theMaxPt == Constants::MaxEnergy ) ) return true;
  if ( inci || incj ) return true;  

  if ( sameFlavourOnly() || oppositeSignOnly() ) {

    int fam1 = family(pitype->id());
    int fam2 = family(pjtype->id());

    if ( fam1 && fam2 ) {
      if ( sameFlavourOnly() && ( abs(fam1) != abs(fam2) ) ) return true;
      if ( oppositeSignOnly() && ( fam1*fam2 > 0 ) ) return true;
    }

  }


  double weight = 1.0;

  Energy minv = (pi+pj).perp();

  if ( !parent->isInside<CutTypes::Energy>(minv,minPt(),maxPt(),weight) ) 
  {
    parent->lastCutWeight(0.0);
    return false;
  }

  parent->lastCutWeight(weight);
  return true;

}

int PairPtCut::family(long id) const {

  int sign = (id>0)?-1:1;

  switch ( id ) {
    case ParticleID::u:
    case ParticleID::ubar:
    case ParticleID::d:
    case ParticleID::dbar:
      return 1*sign; break;
    case ParticleID::c:
    case ParticleID::cbar:
    case ParticleID::s:
    case ParticleID::sbar:
      return 2*sign; break;
    case ParticleID::t:
    case ParticleID::tbar:
    case ParticleID::b:
    case ParticleID::bbar:
      return 3*sign; break;
    case ParticleID::eminus:
    case ParticleID::eplus:
    case ParticleID::nu_e:
    case ParticleID::nu_ebar:
      return 11*sign; break;
    case ParticleID::muminus:
    case ParticleID::muplus:
    case ParticleID::nu_mu:
    case ParticleID::nu_mubar:
      return 12*sign; break;
    case ParticleID::tauminus:
    case ParticleID::tauplus:
    case ParticleID::nu_tau:
    case ParticleID::nu_taubar:
      return 13*sign; break;
  }

  return 0;

}

void PairPtCut::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinPt,GeV) << ounit(theMaxPt,GeV) 
     << theSameFlavourOnly << theOppositeSignOnly
     << theFirstMatcher << theSecondMatcher;
}

void PairPtCut::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinPt,GeV) >> iunit(theMaxPt,GeV) 
     >> theSameFlavourOnly >> theOppositeSignOnly
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PairPtCut,TwoCutBase>
  describeHerwigPairPtCut("Herwig::PairPtCut", "HwMatchboxCuts.so");

void PairPtCut::Init() {

  static ClassDocumentation<PairPtCut> documentation
    ("This class implements a transverse momentum cut on lepton pairs of "
     "final-state particles.");

  static Parameter<PairPtCut,Energy> interfaceMinPt
   ("MinPt",
     "The minimal allowed transverse momentum of the particle pair ",
     &PairPtCut::theMinPt, GeV, 0*GeV, 0*GeV, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Parameter<PairPtCut,Energy> interfaceMaxPt
   ("MaxPt",
     "The maximal allowed transverse momentum of the particle pair ",
     &PairPtCut::theMaxPt, GeV, Constants::MaxEnergy, 0*GeV, Constants::MaxEnergy,
     false, false, Interface::limited);

  static Switch<PairPtCut,bool> interfaceSameFlavourOnly
   ("SameFlavourOnly",
     "Whether cut works on fermion pairs of the same flavour only ",
     &PairPtCut::theSameFlavourOnly, true, false, false);
  static SwitchOption interfaceSameFlavourOnlyYes
    (interfaceSameFlavourOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceSameFlavourOnlyNo
    (interfaceSameFlavourOnly,
     "No",
     "No",
     false);

  static Switch<PairPtCut,bool> interfaceOppositeSignOnly
   ("OppositeSignOnly",
     "Whether cut works on fermion pairs of opposite sign only ",
     &PairPtCut::theOppositeSignOnly, true, false, false);
  static SwitchOption interfaceOppositeSignOnlyYes
    (interfaceOppositeSignOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceOppositeSignOnlyNo
    (interfaceOppositeSignOnly,
     "No",
     "No",
     false);

  static Reference<PairPtCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairPtCut::theFirstMatcher, true, false, true, true, false);

  static Reference<PairPtCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairPtCut::theSecondMatcher, true, false, true, true, false);

}


#line 1 "./PairRapidityCut.cc"
// -*- C++ -*-
//
// PairRapidityCut.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PairRapidityCut class.
//

#include "PairRapidityCut.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/StandardMatchers.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PairRapidityCut::describe() const {
  CurrentGenerator::log() 
    << fullName() << "\n"
    << "matching between: '"
    << theFirstMatcher->name() << "' and '"
//     << theSecondMatcher->name() << "':\n"
    << theSecondMatcher->name() << "':\n";

//     << "y = " << theMinRapidity << " .. " << theMaxRapidity << "\n"
    for ( vector<pair<double,double> >::const_iterator r = yRanges().begin();
	  r != yRanges().end(); ++r ) {
      CurrentGenerator::log() << "y = " << r->first << " .. " << r->second << "\n";
    }

//     << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
  CurrentGenerator::log()
    << "same flavour only = " << (theSameFlavourOnly?"Yes":"No") << " \n"
    << "opposite sign only = " << (theOppositeSignOnly?"Yes":"No") << " \n\n";
}

PairRapidityCut::PairRapidityCut() 
  : thePseudo(false), theSameFlavourOnly(false), theOppositeSignOnly(false) {}

IBPtr PairRapidityCut::clone() const {
  return new_ptr(*this);
}

IBPtr PairRapidityCut::fullclone() const {
  return new_ptr(*this);
}

bool PairRapidityCut::passCuts(tcCutsPtr parent, tcPDPtr pitype, tcPDPtr pjtype,
		                 LorentzMomentum pi, LorentzMomentum pj,
		                 bool inci, bool incj) const {

  bool match = false;
  if ( theFirstMatcher->check(*pitype) && theSecondMatcher->check(*pjtype) ) match = true;
  if ( theFirstMatcher->check(*pjtype) && theSecondMatcher->check(*pitype) ) match = true;
//   if ( !match ||
//        ( theMinRapidity == -Constants::MaxRapidity && theMaxRapidity == Constants::MaxRapidity ) ) return true;
  if ( !match || theYRanges.empty() ) return true;
  if ( inci || incj ) return true;  

  if ( sameFlavourOnly() || oppositeSignOnly() ) {

    int fam1 = family(pitype->id());
    int fam2 = family(pjtype->id());

    if ( fam1 && fam2 ) {
      if ( sameFlavourOnly() && ( abs(fam1) != abs(fam2) ) ) return true;
      if ( oppositeSignOnly() && ( fam1*fam2 > 0 ) ) return true;
    }

  }


  double weight = 1.0;

//   double minv = (pi+pj).rapidity();
// 
//   if ( !parent->isInside<CutTypes::Rapidity>(minv,minRapidity(),maxRapidity(),weight) ) 
//   {
//     parent->lastCutWeight(0.0);
//     return false;
//   }

//// From IdentifiedParticleCut.cc
//   double y = p.rapidity() + parent->currentYHat();
//   for ( vector<pair<double,double> >::const_iterator dy = yRanges().begin();
// 	 dy != yRanges().end(); ++dy ) {
//     if ( !parent->isInside<CutTypes::Rapidity>(y,dy->first,dy->second,weight) ) {
//       parent->lastCutWeight(0.0);
//       return false;
//     }
//   }

  double y = (pi+pj).rapidity() + parent->currentYHat();

  // Actually, why not 
  // double y = (pi+pj).rapidity() + parent->currentYHat() + parent->Y(); 
  // as in ThePEG /Cuts/Cuts.cc, /Cuts/OneCutBase.cc, etc. ???

  for ( vector<pair<double,double> >::const_iterator dy = yRanges().begin();
	dy != yRanges().end(); ++dy ) {
    if (!thePseudo) {
      if ( !parent->isInside<CutTypes::Rapidity>(y,dy->first,dy->second,weight) ) {
        parent->lastCutWeight(0.0);
        return false;
      }
    } else if (thePseudo) {
      //// From ThePEG/Cuts/OneCutBase or ThePEG/Cuts/SimpleKTCut
      // if ( p.mt()*sinh(y) <= p.perp()*sinh(theMinEta) ) return false;
      // if ( p.mt()*sinh(y) >= p.perp()*sinh(theMaxEta) ) return false;
      if ( !parent->isInside<CutTypes::Rapidity>( (pi+pj).mt()*sinh(y)/GeV ,
                                                  (pi+pj).perp()*sinh(dy->first)/GeV ,
                                                  (pi+pj).perp()*sinh(dy->second)/GeV ,
                                                  weight) ) {
        parent->lastCutWeight(0.0);
        return false;
      }
    }
  }

  parent->lastCutWeight(weight);
  return true;

}

int PairRapidityCut::family(long id) const {

  int sign = (id>0)?-1:1;

  switch ( id ) {
    case ParticleID::u:
    case ParticleID::ubar:
    case ParticleID::d:
    case ParticleID::dbar:
      return 1*sign; break;
    case ParticleID::c:
    case ParticleID::cbar:
    case ParticleID::s:
    case ParticleID::sbar:
      return 2*sign; break;
    case ParticleID::t:
    case ParticleID::tbar:
    case ParticleID::b:
    case ParticleID::bbar:
      return 3*sign; break;
    case ParticleID::eminus:
    case ParticleID::eplus:
    case ParticleID::nu_e:
    case ParticleID::nu_ebar:
      return 11*sign; break;
    case ParticleID::muminus:
    case ParticleID::muplus:
    case ParticleID::nu_mu:
    case ParticleID::nu_mubar:
      return 12*sign; break;
    case ParticleID::tauminus:
    case ParticleID::tauplus:
    case ParticleID::nu_tau:
    case ParticleID::nu_taubar:
      return 13*sign; break;
  }

  return 0;

}

string PairRapidityCut::doYRange(string in) {
  istringstream ins(in);
  double first, second;
  ins >> first >> second;
  if ( first > second )
    swap(first,second);
  theYRanges.push_back(make_pair(first,second));
  return "";
}

void PairRapidityCut::persistentOutput(PersistentOStream & os) const {
//   os << theMinRapidity << theMaxRapidity 
  os << theYRanges << thePseudo
     << theSameFlavourOnly << theOppositeSignOnly
     << theFirstMatcher << theSecondMatcher;
}

void PairRapidityCut::persistentInput(PersistentIStream & is, int) {
//   is >> theMinRapidity >> theMaxRapidity
  is >> theYRanges >> thePseudo
     >> theSameFlavourOnly >> theOppositeSignOnly
     >> theFirstMatcher >> theSecondMatcher;
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PairRapidityCut,TwoCutBase>
  describeHerwigPairRapidityCut("Herwig::PairRapidityCut", "HwMatchboxCuts.so");

void PairRapidityCut::Init() {

  static ClassDocumentation<PairRapidityCut> documentation
    ("This class implements a rapidity cut on lepton pairs of "
     "final-state particles.");

//   static Parameter<PairRapidityCut,double> interfaceMinRapidity
//    ("MinRapidity",
//      "The minimal allowed rapditiy of the particle pair ",
//      &PairRapidityCut::theMinRapidity, -Constants::MaxRapidity, -Constants::MaxRapidity, Constants::MaxRapidity,
//      false, false, Interface::limited);
// 
//   static Parameter<PairRapidityCut,double> interfaceMaxRapidity
//    ("MaxRapidity",
//      "The maximal allowed rapidity of the particle pair ",
//      &PairRapidityCut::theMaxRapidity, Constants::MaxRapidity, -Constants::MaxRapidity, Constants::MaxRapidity,
//      false, false, Interface::limited);

  static Command<PairRapidityCut> interfaceYRange
    ("YRange",
     "Insert a rapidity range.",
     &PairRapidityCut::doYRange, false);

  static Switch<PairRapidityCut,bool> interfacePseudo
   ("Pseudo",
     "Use pseudo rapidity instead of rapidity ",
     &PairRapidityCut::thePseudo, false, false, false);
  static SwitchOption interfacePseudoNo
    (interfacePseudo,
     "No",
     "No",
     false);
  static SwitchOption interfacePseudoYes
    (interfacePseudo,
     "Yes",
     "Yes",
     true);

  static Switch<PairRapidityCut,bool> interfaceSameFlavourOnly
   ("SameFlavourOnly",
     "Whether cut works on fermion pairs of the same flavour only ",
     &PairRapidityCut::theSameFlavourOnly, true, false, false);
  static SwitchOption interfaceSameFlavourOnlyYes
    (interfaceSameFlavourOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceSameFlavourOnlyNo
    (interfaceSameFlavourOnly,
     "No",
     "No",
     false);

  static Switch<PairRapidityCut,bool> interfaceOppositeSignOnly
   ("OppositeSignOnly",
     "Whether cut works on fermion pairs of opposite sign only ",
     &PairRapidityCut::theOppositeSignOnly, true, false, false);
  static SwitchOption interfaceOppositeSignOnlyYes
    (interfaceOppositeSignOnly,
     "Yes",
     "Yes",
     true);
  static SwitchOption interfaceOppositeSignOnlyNo
    (interfaceOppositeSignOnly,
     "No",
     "No",
     false);

  static Reference<PairRapidityCut,MatcherBase> interfaceFirstMatcher
    ("FirstMatcher",
     "Matcher for first particle of type pitype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairRapidityCut::theFirstMatcher, true, false, true, true, false);

  static Reference<PairRapidityCut,MatcherBase> interfaceSecondMatcher
    ("SecondMatcher",
     "Matcher for second particle of type pjtype in the pair (pitype,pjtype). "
     "Only particles matching this object will be affected by the cut. ",
     &PairRapidityCut::theSecondMatcher, true, false, true, true, false);

}


