#line 1 "./FFggxDipole.cc"
// -*- C++ -*-
//
// FFggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFggxDipole class.
//

#include "FFggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightInvertedTildeKinematics.h"

using namespace Herwig;

FFggxDipole::FFggxDipole() 
  : SubtractionDipole() {}

IBPtr FFggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FFggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double res = 
    1./(1.-z*(1.-y)) + 1./(1.-(1.-z)*(1.-y)) - 2. + z*(1.-z);

  res *= -ccme2;

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < y )
    return 0.0;
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double diag = 1./(1.-z*(1.-y))+1./(1.-(1.-z)*(1.-y))-2.;
  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-diag,pc,prop/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFggxDipole::persistentOutput(PersistentOStream &) const {
}

void FFggxDipole::persistentInput(PersistentIStream &, int) {
}

void FFggxDipole::Init() {

  static ClassDocumentation<FFggxDipole> documentation
    ("FFggxDipole");

  DipoleRepository::registerDipole<0,FFggxDipole,FFLightTildeKinematics,FFLightInvertedTildeKinematics>
    ("FFggxDipole","FFLightTildeKinematics","FFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFggxDipole,SubtractionDipole>
describeHerwigFFggxDipole("Herwig::FFggxDipole", "Herwig.so");
#line 1 "./FFqgxDipole.cc"
// -*- C++ -*-
//
// FFqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFqgxDipole class.
//

#include "FFqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightInvertedTildeKinematics.h"

using namespace Herwig;

FFqgxDipole::FFqgxDipole() 
  : SubtractionDipole() {}

IBPtr FFqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FFqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - (1.+z) );

  res *= -ccme2;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < y )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - (1.+z) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FFqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FFqgxDipole::Init() {

  static ClassDocumentation<FFqgxDipole> documentation
    ("FFqgxDipole");

  DipoleRepository::registerDipole<0,FFqgxDipole,FFLightTildeKinematics,FFLightInvertedTildeKinematics>
    ("FFqgxDipole","FFLightTildeKinematics","FFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFqgxDipole,SubtractionDipole>
describeHerwigFFqgxDipole("Herwig::FFqgxDipole", "Herwig.so");
#line 1 "./FFqqxDipole.cc"
// -*- C++ -*-
//
// FFqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFqqxDipole class.
//

#include "FFqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFLightInvertedTildeKinematics.h"

using namespace Herwig;

FFqqxDipole::FFqqxDipole() 
  : SubtractionDipole() {}

IBPtr FFqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FFqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double res = 1.-2.*z*(1.-z);

  res *= -ccme2;

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
 
  if ( alpha() < y )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-1.,pc,-prop/4.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FFqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FFqqxDipole::Init() {

  static ClassDocumentation<FFqqxDipole> documentation
    ("FFqqxDipole");

  DipoleRepository::registerDipole<0,FFqqxDipole,FFLightTildeKinematics,FFLightInvertedTildeKinematics>
    ("FFqqxDipole","FFLightTildeKinematics","FFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFqqxDipole,SubtractionDipole>
describeHerwigFFqqxDipole("Herwig::FFqqxDipole", "Herwig.so");
#line 1 "./FIggxDipole.cc"
// -*- C++ -*-
//
// FIggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIggxDipole class.
//

#include "FIggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"

using namespace Herwig;

FIggxDipole::FIggxDipole() 
  : SubtractionDipole() {}

IBPtr FIggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = 
    1./((1.-z)+(1.-x)) + 1./(z+(1.-x)) - 2.+z*(1.-z) 
    + (1.-x)*(1.+x*z*(1.-z))
    ;

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < (1.-x) )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double diag = 
    1./(1.-z+1.-x) + 1./(z+1.-x) - 2. 
    //+ (1.-x)*(1.+x*z*(1.-z))
    ;
  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-diag,pc,prop/(2.*x));

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIggxDipole::persistentOutput(PersistentOStream &) const {
}

void FIggxDipole::persistentInput(PersistentIStream &, int) {
}

void FIggxDipole::Init() {

  static ClassDocumentation<FIggxDipole> documentation
    ("FIggxDipole");

  DipoleRepository::registerDipole<0,FIggxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
    ("FIggxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIggxDipole,SubtractionDipole>
describeHerwigFIggxDipole("Herwig::FIggxDipole", "Herwig.so");
#line 1 "./FIqgxDipole.cc"
// -*- C++ -*-
//
// FIqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqgxDipole class.
//

#include "FIqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"

using namespace Herwig;

FIqgxDipole::FIqgxDipole() 
  : SubtractionDipole() {}

IBPtr FIqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 
    2./(1.-z+(1.-x)) -(1.+z) 
    + (1.-x)*(1.+3.*x*z) 
    );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < (1.-x) )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 
    2./(1.-z+(1.-x)) - (1.+z) 
    //+ (1.-x)*(1.+3.*x*z) 
    );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FIqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FIqgxDipole::Init() {

  static ClassDocumentation<FIqgxDipole> documentation
    ("FIqgxDipole");

  DipoleRepository::registerDipole<0,FIqgxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
    ("FIqgxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIqgxDipole,SubtractionDipole>
describeHerwigFIqgxDipole("Herwig::FIqgxDipole", "Herwig.so");
#line 1 "./FIqqxDipole.cc"
// -*- C++ -*-
//
// FIqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqqxDipole class.
//

#include "FIqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"

using namespace Herwig;

FIqqxDipole::FIqqxDipole() 
  : SubtractionDipole() {}

IBPtr FIqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = 1.-2.*z*(1.-z);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  if ( alpha() < (1.-x) )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-1.,pc,-prop/(4.*x));

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FIqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FIqqxDipole::Init() {

  static ClassDocumentation<FIqqxDipole> documentation
    ("FIqqxDipole");

  DipoleRepository::registerDipole<0,FIqqxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
    ("FIqqxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIqqxDipole,SubtractionDipole>
describeHerwigFIqqxDipole("Herwig::FIqqxDipole", "Herwig.so");
#line 1 "./IFggxDipole.cc"
// -*- C++ -*-
//
// IFggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFggxDipole class.
//

#include "IFggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFggxDipole::IFggxDipole() 
  : SubtractionDipole() {}

IBPtr IFggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IFggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = 1./(1.-x+u) + (1.-x)/x - 1. + x*(1.-x);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  if ( alpha() < u )
    return 0.0;
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double diag = 1./(1.-x+u)-1.+x*(1.-x);
  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]/u -
    realEmissionME()->lastXComb().meMomenta()[realSpectator()]/(1.-u);

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= u*(1.-u)*(1.-x)/x;

  SpinCorrelationTensor corr(-diag,pc,sc);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFggxDipole::persistentOutput(PersistentOStream &) const {
}

void IFggxDipole::persistentInput(PersistentIStream &, int) {
}

void IFggxDipole::Init() {

  static ClassDocumentation<IFggxDipole> documentation
    ("IFggxDipole");

  DipoleRepository::registerDipole<0,IFggxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFggxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFggxDipole,SubtractionDipole>
describeHerwigIFggxDipole("Herwig::IFggxDipole", "Herwig.so");
#line 1 "./IFgqxDipole.cc"
// -*- C++ -*-
//
// IFgqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgqxDipole class.
//

#include "IFgqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFgqxDipole::IFgqxDipole() 
  : SubtractionDipole() {}

IBPtr IFgqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFgqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFgqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emitter]->id() == ParticleID::g &&
    abs(partons[emission]->id()) < 6 &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IFgqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFgqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  if ( alpha() < u )
    return 0.0;
  
  Energy2 prop =
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFgqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFgqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFgqxDipole::Init() {

  static ClassDocumentation<IFgqxDipole> documentation
    ("IFgqxDipole");

  DipoleRepository::registerDipole<0,IFgqxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFgqxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFgqxDipole,SubtractionDipole>
describeHerwigIFgqxDipole("Herwig::IFgqxDipole", "Herwig.so");
#line 1 "./IFqgxDipole.cc"
// -*- C++ -*-
//
// IFqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqgxDipole class.
//

#include "IFqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFqgxDipole::IFqgxDipole() 
  : SubtractionDipole() {}

IBPtr IFqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IFqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 
    2./(1.-x+u) - (1.+x) 
    + u*(1.+3.*x*(1.-u)) 
    );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  if ( alpha() < u )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 
    2./(1.-x+u) - (1.+x) 
    //+ u*(1.+3.*x*(1.-u)) 
    );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFqgxDipole::persistentOutput(PersistentOStream &) const {
}

void IFqgxDipole::persistentInput(PersistentIStream &, int) {
}

void IFqgxDipole::Init() {

  static ClassDocumentation<IFqgxDipole> documentation
    ("IFqgxDipole");

  DipoleRepository::registerDipole<0,IFqgxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFqgxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFqgxDipole,SubtractionDipole>
describeHerwigIFqgxDipole("Herwig::IFqgxDipole", "Herwig.so");
#line 1 "./IFqqxDipole.cc"
// -*- C++ -*-
//
// IFqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqqxDipole class.
//

#include "IFqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"

using namespace Herwig;

IFqqxDipole::IFqqxDipole() 
  : SubtractionDipole() {}

IBPtr IFqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() - partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IFqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = x + 2.*(1.-x)/x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  
  if ( alpha() < u )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]/u -
    realEmissionME()->lastXComb().meMomenta()[realSpectator()]/(1.-u);

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= u*(1.-u)*(1.-x)/x;

  SpinCorrelationTensor corr(-x,pc,sc/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFqqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFqqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFqqxDipole::Init() {

  static ClassDocumentation<IFqqxDipole> documentation
    ("IFqqxDipole");

  DipoleRepository::registerDipole<0,IFqqxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
    ("IFqqxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFqqxDipole,SubtractionDipole>
describeHerwigIFqqxDipole("Herwig::IFqqxDipole", "Herwig.so");
#line 1 "./IIggxDipole.cc"
// -*- C++ -*-
//
// IIggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIggxDipole class.
//

#include "IIggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightInvertedTildeKinematics.h"

using namespace Herwig;

IIggxDipole::IIggxDipole() 
  : SubtractionDipole() {}

IBPtr IIggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = x/(1.-x) + (1.-x)/x + x*(1.-x);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IIggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double diag = x/(1.-x)+x*(1.-x);
  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()] -
    v*realEmissionME()->lastXComb().meMomenta()[realSpectator()];

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= (1.-x)/(x*v);

  SpinCorrelationTensor corr(-diag,pc,sc);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIggxDipole::persistentOutput(PersistentOStream &) const {
}

void IIggxDipole::persistentInput(PersistentIStream &, int) {
}

void IIggxDipole::Init() {

  static ClassDocumentation<IIggxDipole> documentation
    ("IIggxDipole");

  DipoleRepository::registerDipole<0,IIggxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIggxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIggxDipole,SubtractionDipole>
describeHerwigIIggxDipole("Herwig::IIggxDipole", "Herwig.so");
#line 1 "./IIgqxDipole.cc"
// -*- C++ -*-
//
// IIgqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIgqxDipole class.
//

#include "IIgqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightInvertedTildeKinematics.h"

using namespace Herwig;

IIgqxDipole::IIgqxDipole() 
  : SubtractionDipole() {}

IBPtr IIgqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIgqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIgqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    partons[emitter]->id() == ParticleID::g &&
    abs(partons[emission]->id()) < 6 &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIgqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IIgqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIgqxDipole::persistentOutput(PersistentOStream &) const {
}

void IIgqxDipole::persistentInput(PersistentIStream &, int) {
}

void IIgqxDipole::Init() {

  static ClassDocumentation<IIgqxDipole> documentation
    ("IIgqxDipole");

  DipoleRepository::registerDipole<0,IIgqxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIgqxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIgqxDipole,SubtractionDipole>
describeHerwigIIgqxDipole("Herwig::IIgqxDipole", "Herwig.so");
#line 1 "./IIqgxDipole.cc"
// -*- C++ -*-
//
// IIqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqgxDipole class.
//

#include "IIqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightInvertedTildeKinematics.h"

using namespace Herwig;

IIqgxDipole::IIqgxDipole() 
  : SubtractionDipole() {}

IBPtr IIqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 2./(1.-x) - (1.+x);

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IIqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 2./(1.-x) - (1.+x);

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIqgxDipole::persistentOutput(PersistentOStream &) const {
}

void IIqgxDipole::persistentInput(PersistentIStream &, int) {
}

void IIqgxDipole::Init() {

  static ClassDocumentation<IIqgxDipole> documentation
    ("IIqgxDipole");

  DipoleRepository::registerDipole<0,IIqgxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIqgxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIqgxDipole,SubtractionDipole>
describeHerwigIIqgxDipole("Herwig::IIqgxDipole", "Herwig.so");
#line 1 "./IIqqxDipole.cc"
// -*- C++ -*-
//
// IIqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IIqqxDipole class.
//

#include "IIqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IILightInvertedTildeKinematics.h"

using namespace Herwig;

IIqqxDipole::IIqqxDipole() 
  : SubtractionDipole() {}

IBPtr IIqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IIqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IIqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator < 2 &&
    abs(partons[emission]->id()) < 6 &&
    abs(partons[emitter]->id()) < 6 &&
    partons[emission]->id() - partons[emitter]->id() == 0 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double IIqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res = (1.+sqr(1.-x))/x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IIqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  
  if ( alpha() < v )
    return 0.0;

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()] -
    v*realEmissionME()->lastXComb().meMomenta()[realSpectator()];

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= (1.-x)/(x*v);

  SpinCorrelationTensor corr(-x,pc,sc/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IIqqxDipole::persistentOutput(PersistentOStream &) const {
}

void IIqqxDipole::persistentInput(PersistentIStream &, int) {
}

void IIqqxDipole::Init() {

  static ClassDocumentation<IIqqxDipole> documentation
    ("IIqqxDipole");

  DipoleRepository::registerDipole<0,IIqqxDipole,IILightTildeKinematics,IILightInvertedTildeKinematics>
    ("IIqqxDipole","IILightTildeKinematics","IILightInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IIqqxDipole,SubtractionDipole>
describeHerwigIIqqxDipole("Herwig::IIqqxDipole", "Herwig.so");
#line 1 "./FFMggxDipole.cc"
// -*- C++ -*-
//
// FFMggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMggxDipole class.
//

#include "FFMggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FFMggxDipole::FFMggxDipole() 
  : SubtractionDipole() {}

IBPtr FFMggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double FFMggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  // masses, g->gg all masses zero except spectator
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms, viji = 1
  double vijk = sqrt( sqr(2.*muj2+(1.-muj2)*(1.-y))-4.*muj2 ) / ((1.-muj2)*(1.-y));

  Energy2 prop =
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
  (realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double zp = 0.5*(1.+vijk);
  double zm = 0.5*(1.-vijk);

  double res = -ccme2;

  // extra mass terms all = 0.
  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= (1./(1-z*(1-y))+1./(1-(1.-z)*(1.-y))+(z*(1.-z)-zm*zp-2.)/vijk);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFMggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses, g->gg all masses zero except spectator
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muj2)*(1.-y))-4.*muj2 ) / ((1.-muj2)*(1.-y));
  
  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double diag = 1./(1-z*(1-y))+1./(1-(1.-z)*(1.-y))-2./vijk; // kappa=0
  double zim = z-0.5*(1.-vijk), zjm = (1.-z)-0.5*(1.-vijk);
  Lorentz5Momentum pc = 
    zim*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    zjm*realEmissionME()->lastXComb().meMomenta()[realEmission()];
  
  SpinCorrelationTensor corr(-diag,pc,prop/2.*vijk);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  // extra mass terms all = 0.
  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFMggxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMggxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMggxDipole::Init() {

  static ClassDocumentation<FFMggxDipole> documentation
    ("FFMggxDipole");

  DipoleRepository::registerDipole<0,FFMggxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMggxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMggxDipole,SubtractionDipole>
describeHerwigFFMggxDipole("Herwig::FFMggxDipole", "Herwig.so");
#line 1 "./FFMqgxDipole.cc"
// -*- C++ -*-
//
// FFMqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMqgxDipole class.
//

#include "FFMqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FFMqgxDipole::FFMqgxDipole() 
  : SubtractionDipole() {}

IBPtr FFMqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
//     abs(partons[emitter]->id()) < 7 &&
    ( abs(partons[emitter]->id()) < 7 || abs(partons[emitter]->id()) == 1000021 ) &&
    !(partons[emitter]->hardProcessMass() == ZERO &&
      partons[spectator]->hardProcessMass() == ZERO);
}

double FFMqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muQ2)+sqr(muj2)-2.*(muQ2+muj2+muQ2*muj2) ) / (1.-muQ2-muj2);

  Energy2 prop =
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
  (realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g;k} subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel: mi2+m2-Mi2 = mQ2+0-mQ2
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( (1.+z) + muQ2*sqr(lastDipoleScale())*2./prop ) );

  res *= -ccme2;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;
}

double FFMqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muQ2)+sqr(muj2)-2.*(muQ2+muj2+muQ2*muj2) ) / (1.-muQ2-muj2);

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g;k} subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel: mi2+m2-Mi2 = mQ2+0-mQ2
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( (1.+z) + muQ2*sqr(lastDipoleScale())*2./prop ) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFMqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMqgxDipole::Init() {

  static ClassDocumentation<FFMqgxDipole> documentation
    ("FFMqgxDipole");

  DipoleRepository::registerDipole<0,FFMqgxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMqgxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMqgxDipole,SubtractionDipole>
describeHerwigFFMqgxDipole("Herwig::FFMqgxDipole", "Herwig.so");
#line 1 "./FFMsqgxDipole.cc"
// -*- C++ -*-
//
// FFMsqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMsqgxDipole class.
//

#include "FFMsqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FFMsqgxDipole::FFMsqgxDipole() 
  : SubtractionDipole() {}

IBPtr FFMsqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMsqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMsqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    ((abs(partons[emitter]->id())> 1000000 && abs(partons[emitter]->id())< 1000007) ||
     (abs(partons[emitter]->id())> 2000000 && abs(partons[emitter]->id())< 2000007)) &&
    !(partons[emitter]->hardProcessMass() == ZERO &&
      partons[spectator]->hardProcessMass() == ZERO);
}


double FFMsqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  // masses
  double muSQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muSQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muSQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muSQ2)+sqr(muj2)-2.*(muSQ2+muj2+muSQ2*muj2) ) / (1.-muSQ2-muj2);

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( 2. + muSQ2*sqr(lastDipoleScale())*2./prop ));

  res *= -ccme2;

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;


}

double FFMsqgxDipole::me2() const {
  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double muSQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-muSQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-muSQ2-muj2)*(1.-y));
  double vbar = sqrt( 1.+sqr(muSQ2)+sqr(muj2)-2.*(muSQ2+muj2+muSQ2*muj2) ) / (1.-muSQ2-muj2);

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  // extra mass terms cancel: mi2+m2-Mi2 = mQ2+0-mQ2
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= ( 2./(1.-z*(1.-y)) - vbar/vijk * ( 2. + muSQ2*sqr(lastDipoleScale())*2./prop ) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  
  return res;

}

void FFMsqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMsqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMsqgxDipole::Init() {

  static ClassDocumentation<FFMsqgxDipole> documentation
    ("FFMsqgxDipole");

  DipoleRepository::registerDipole<0,FFMsqgxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMsqgxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMsqgxDipole,SubtractionDipole>
describeHerwigFFMsqgxDipole("Herwig::FFMsqgxDipole", "Herwig.so");
#line 1 "./FFMqqxDipole.cc"
// -*- C++ -*-
//
// FFMqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMqqxDipole class.
//

#include "FFMqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FFMqqxDipole::FFMqqxDipole() 
  : SubtractionDipole() {}

IBPtr FFMqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FFMqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FFMqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator > 1 &&
    abs(partons[emission]->id()) < 7 &&
    abs(partons[emitter]->id()) < 7 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
    !(partons[emission]->hardProcessMass() == ZERO &&
      partons[emitter]->hardProcessMass() == ZERO &&
      partons[spectator]->hardProcessMass() == ZERO);
}

double FFMqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  Energy2 mQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass() );
  // massive extra terms
  double t = 1.-2.*muQ2-muj2;
  double vijk = sqrt( sqr(2.*muj2+t*(1.-y))-4.*muj2 ) / (t*(1.-y));
  double viji = sqrt( sqr(t*y) - 4.*sqr(muQ2) ) / ( t*y + 2.*muQ2);

  Energy2 prop =
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
  (realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double zp = 0.5*(1.+viji*vijk);
  double zm = 0.5*(1.-viji*vijk);

  // kappa=0 -- otherwise: extra term

  double res = -ccme2;

  res *= (1.-2.*(z*(1-z)-zp*zm));

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/ ((prop+2.*mQ2)*vijk);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FFMqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  
  // masses
  double muQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass() / lastDipoleScale() );
  double muj2 = sqr( realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass() / lastDipoleScale() );
  Energy2 mQ2 = sqr( realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass() );
  // massive extra terms
  double vijk = sqrt( sqr(2.*muj2+(1.-2.*muQ2-muj2)*(1.-y))-4.*muj2 ) / ((1.-2.*muQ2-muj2)*(1.-y));

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]));

  double zim = z-0.5*(1.-vijk), zjm = (1.-z)-0.5*(1.-vijk);
  Lorentz5Momentum pc = 
    zim*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    zjm*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  // kappa=0 -- otherwise: extra diagonal term (instead of just -1.)
  SpinCorrelationTensor corr(-1.,pc,-(prop+2.*mQ2)/4.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/ ((prop+2.*mQ2)*vijk);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FFMqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FFMqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FFMqqxDipole::Init() {

  static ClassDocumentation<FFMqqxDipole> documentation
    ("FFMqqxDipole");

  DipoleRepository::registerDipole<0,FFMqqxDipole,FFMassiveTildeKinematics,FFMassiveInvertedTildeKinematics>
    ("FFMqqxDipole","FFMassiveTildeKinematics","FFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMqqxDipole,SubtractionDipole>
describeHerwigFFMqqxDipole("Herwig::FFMqqxDipole", "Herwig.so");
#line 1 "./FIMqgxDipole.cc"
// -*- C++ -*-
//
// FIqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMqgxDipole class.
//

#include "FIMqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightTildeKinematics.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/FILightInvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FIMqgxDipole::FIMqgxDipole() 
  : SubtractionDipole() {}

IBPtr FIMqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIMqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIMqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
//     abs(partons[emitter]->id()) < 7 &&
    ( abs(partons[emitter]->id()) < 7 || abs(partons[emitter]->id()) == 1000021 ) &&
    partons[emitter]->hardProcessMass() != ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIMqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g}^a subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term taken from FIqgxDipole implementation
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-z+(1.-x)) - (1.+z) 
    // + (1.-x)*(1.+3.*x*z) 
    - sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x
    );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIMqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  if ( realEmissionME()->lastXComb().mePartonData()[realEmitter()]->id() == 1000021 )
    CF = SM().Nc(); // For the SUSY D_{gluino,g}^a subtraction dipole we need to replace CF by CA=Nc

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term taken from FIqgxDipole implementation
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-z+(1.-x)) - (1.+z) 
    // + (1.-x)*(1.+3.*x*z) 
    - sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x
    );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIMqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FIMqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FIMqgxDipole::Init() {

  static ClassDocumentation<FIMqgxDipole> documentation
    ("FIMqgxDipole");

//  DipoleRepository::registerDipole<0,FIMqgxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
//    ("FIMqgxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,FIMqgxDipole,FIMassiveTildeKinematics,FIMassiveInvertedTildeKinematics>
    ("FIMqgxDipole","FIMassiveTildeKinematics","FIMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMqgxDipole,SubtractionDipole>
describeHerwigFIMqgxDipole("Herwig::FIMqgxDipole", "Herwig.so");
#line 1 "./FIMsqgxDipole.cc"
// -*- C++ -*-
//
// FIqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMsqgxDipole class.
//

#include "FIMsqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FIMsqgxDipole::FIMsqgxDipole() 
  : SubtractionDipole() {}

IBPtr FIMsqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIMsqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIMsqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    partons[emission]->id() == ParticleID::g &&
    ((abs(partons[emitter]->id())> 1000000 && abs(partons[emitter]->id())< 1000007) ||
     (abs(partons[emitter]->id())> 2000000 && abs(partons[emitter]->id())< 2000007)) &&
    partons[emitter]->hardProcessMass() != ZERO &&
    partons[spectator]->hardProcessMass() == ZERO;
}

double FIMsqgxDipole::me2Avg(double ccme2) const {
  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

//   // NOTE: extra term taken from FIqgxDipole implementation??
//   res *= ( 2./(1.-z+(1.-x)) - 2. +(1.-x)*(1.+3.*x*z) -
//   	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass()) / prop * 2.*x);
  // NOTE: CR: extra term switched off in massive implementation for the moment,
  //           mass of realEmission changed to mass of realEmitter
  res *= ( 2./(1.-z+(1.-x)) - 2. -
  	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x);

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIMsqgxDipole::me2() const {
  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  // extra mass terms cancel
  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

//   // NOTE: extra term taken from FIqgxDipole implementation
//   res *= ( 2./(1.-z+(1.-x)) -2. +(1.-x)*(1.+3.*x*z) -
//   	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmission()]->hardProcessMass()) / prop * 2.*x);
  // NOTE: CR: extra term switched off in massive implementation for the moment,
  //           mass of realEmission changed to mass of realEmitter
  res *= ( 2./(1.-z+(1.-x)) - 2. -
  	   sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass()) / prop * 2.*x);

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIMsqgxDipole::persistentOutput(PersistentOStream &) const {
}

void FIMsqgxDipole::persistentInput(PersistentIStream &, int) {
}

void FIMsqgxDipole::Init() {

  static ClassDocumentation<FIMsqgxDipole> documentation
    ("FIMsqgxDipole");

  DipoleRepository::registerDipole<0,FIMsqgxDipole,FIMassiveTildeKinematics,FIMassiveInvertedTildeKinematics>
    ("FIMsqgxDipole","FIMassiveTildeKinematics","FIMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMsqgxDipole,SubtractionDipole>
describeHerwigFIMsqgxDipole("Herwig::FIMsqgxDipole", "Herwig.so");
#line 1 "./FIMqqxDipole.cc"
// -*- C++ -*-
//
// FIMqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIqqxDipole class.
//

#include "FIMqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/FIMassiveInvertedTildeKinematics.h"

using namespace Herwig;

FIMqqxDipole::FIMqqxDipole() 
  : SubtractionDipole() {}

IBPtr FIMqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr FIMqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool FIMqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter > 1 && spectator < 2 &&
    abs(partons[emission]->id()) < 7 &&
    abs(partons[emitter]->id()) < 7 &&
    partons[emission]->id() + partons[emitter]->id() == 0 &&
//    !(partons[emitter]->hardProcessMass() == ZERO &&
//      partons[emission]->hardProcessMass() == ZERO &&
//      partons[spectator]->hardProcessMass() == ZERO);
    !(partons[emitter]->hardProcessMass() == ZERO &&
      partons[emission]->hardProcessMass() == ZERO) &&
      partons[spectator]->hardProcessMass() == ZERO;
}

double FIMqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Energy2 mQ2 = sqr(realEmissionME()->lastXComb().mePartonData()[realEmitter()]->hardProcessMass());
//  double muQ2 = x * mQ2 /
  double muQ2 = 0.5 * z * mQ2 /
    ((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
     (realEmissionME()->lastXComb().meMomenta()[realSpectator()]));

  // mu_ij=0, mu_i=mu_j=mu_Q.
//  double zm = ( 1.-x - sqrt( sqr(1.-x-2.*muQ2) - 4.*muQ2 ) ) / ( 2.*(1.-x) );
//  double zp = ( 1.-x + sqrt( sqr(1.-x-2.*muQ2) - 4.*muQ2 ) ) / ( 2.*(1.-x) );
  double zm = ( 1.-x - sqrt( sqr(1.-x-2.*muQ2) - 4.*sqr(muQ2) ) ) / ( 2.*(1.-x) );
  double zp = ( 1.-x + sqrt( sqr(1.-x-2.*muQ2) - 4.*sqr(muQ2) ) ) / ( 2.*(1.-x) );

  double res = 1.-2.*(z-zm)*(zp-z);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/( prop+2.*mQ2*x );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double FIMqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Energy2 mQ2 = sqr((realEmissionME()->lastXComb().mePartonData()[realEmitter()])->hardProcessMass());

  Lorentz5Momentum pc = 
    z*realEmissionME()->lastXComb().meMomenta()[realEmitter()] -
    (1.-z)*realEmissionME()->lastXComb().meMomenta()[realEmission()];

  SpinCorrelationTensor corr(-1.,pc,-(prop+2.*mQ2*x)/(4.*x));

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 4.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/( prop+2.*mQ2*x );

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void FIMqqxDipole::persistentOutput(PersistentOStream &) const {
}

void FIMqqxDipole::persistentInput(PersistentIStream &, int) {
}

void FIMqqxDipole::Init() {

  static ClassDocumentation<FIMqqxDipole> documentation
    ("FIMqqxDipole");

//  DipoleRepository::registerDipole<0,FIMqqxDipole,FILightTildeKinematics,FILightInvertedTildeKinematics>
//    ("FIMqqxDipole","FILightTildeKinematics","FILightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,FIMqqxDipole,FIMassiveTildeKinematics,FIMassiveInvertedTildeKinematics>
    ("FIMqqxDipole","FIMassiveTildeKinematics","FIMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMqqxDipole,SubtractionDipole>
describeHerwigFIMqqxDipole("Herwig::FIMqqxDipole", "Herwig.so");
#line 1 "./IFMggxDipole.cc"
// -*- C++ -*-
//
// IFMggxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFggxDipole class.
//

#include "IFMggxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

IFMggxDipole::IFMggxDipole() 
  : SubtractionDipole() {}

IBPtr IFMggxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMggxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMggxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    partons[emitter]->id() == ParticleID::g &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMggxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double muj2 = sqr( (realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass()) ) /
    (2.* (realEmissionME()->lastXComb().meMomenta()[bornSpectator()])*
     (realEmissionME()->lastXComb().meMomenta()[realEmitter()]) );

  double res = 1./(1.-x+u) + (1.-x)/x - 1. + x*(1.-x) -
    muj2/x*u/(1.-u);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFMggxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double diag = 1./(1.-x+u)-1.+x*(1.-x);
  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]/u -
    realEmissionME()->lastXComb().meMomenta()[realSpectator()]/(1.-u);

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= u*(1.-u)*(1.-x)/x;

  SpinCorrelationTensor corr(-diag,pc,sc);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  res *= 16.*Constants::pi*SM().Nc()*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFMggxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMggxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMggxDipole::Init() {

  static ClassDocumentation<IFMggxDipole> documentation
    ("IFMggxDipole");

  DipoleRepository::registerDipole<0,IFMggxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMggxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMggxDipole,SubtractionDipole>
describeHerwigIFMggxDipole("Herwig::IFMggxDipole", "Herwig.so");
#line 1 "./IFMgqxDipole.cc"
// -*- C++ -*-
//
// IFMgqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFgqxDipole class.
//

#include "IFMgqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

IFMgqxDipole::IFMgqxDipole() 
  : SubtractionDipole() {}

IBPtr IFMgqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMgqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMgqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emitter]->id() == ParticleID::g &&
    abs(partons[emission]->id()) < 7 &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMgqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFMgqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double res =
    8.*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= .5 * ( 1.-2.*x*(1.-x) );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFMgqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMgqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMgqxDipole::Init() {

  static ClassDocumentation<IFMgqxDipole> documentation
    ("IFMgqxDipole");

  DipoleRepository::registerDipole<0,IFMgqxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMgqxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMgqxDipole,SubtractionDipole>
describeHerwigIFMgqxDipole("Herwig::IFMgqxDipole", "Herwig.so");
#line 1 "./IFMqgxDipole.cc"
// -*- C++ -*-
//
// IFMqgxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqgxDipole class.
//

#include "IFMqgxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

IFMqgxDipole::IFMqgxDipole() 
  : SubtractionDipole() {}

IBPtr IFMqgxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqgxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMqgxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
    partons[emission]->id() == ParticleID::g &&
    abs(partons[emitter]->id()) < 7 &&
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMqgxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term same as in IFqgxDipole
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-x+u) - (1.+x) 
    // + u*(1.+3.*x*(1.-u)) 
    );

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFMqgxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  double res =
    8.*Constants::pi*CF*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  // NOTE: extra term same as in IFqgxDipole
  // NOTE: extra term switched off for the moment in the massive case
  res *= ( 
    2./(1.-x+u) - (1.+x) 
    // + u*(1.+3.*x*(1.-u)) 
    );

  res *= -underlyingBornME()->colourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()));

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFMqgxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMqgxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMqgxDipole::Init() {

  static ClassDocumentation<IFMqgxDipole> documentation
    ("IFMqgxDipole");

//  DipoleRepository::registerDipole<0,IFMqgxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
//    ("IFMqgxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,IFMqgxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMqgxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMqgxDipole,SubtractionDipole>
describeHerwigIFMqgxDipole("Herwig::IFMqgxDipole", "Herwig.so");
#line 1 "./IFMqqxDipole.cc"
// -*- C++ -*-
//
// IFMqqxDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFqqxDipole class.
//

#include "IFMqqxDipole.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "Herwig/MatrixElement/Matchbox/Base/DipoleRepository.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinCorrelationTensor.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightTildeKinematics.h"
//#include "Herwig/MatrixElement/Matchbox/Phasespace/IFLightInvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/IFMassiveInvertedTildeKinematics.h"

using namespace Herwig;

IFMqqxDipole::IFMqqxDipole() 
  : SubtractionDipole() {}

IBPtr IFMqqxDipole::clone() const {
  return new_ptr(*this);
}

IBPtr IFMqqxDipole::fullclone() const {
  return new_ptr(*this);
}

bool IFMqqxDipole::canHandle(const cPDVector& partons,
			    int emitter, int emission, int spectator) const {
  return
    emitter < 2 && spectator > 1 &&
//    abs(partons[emission]->id()) < 6 &&
//    abs(partons[emitter]->id()) < 6 &&
    abs(partons[emission]->id()) < 7 &&
    abs(partons[emitter]->id()) < 7 &&
    partons[emission]->id() - partons[emitter]->id() == 0 &&
//    !(partons[emitter]->hardProcessMass() == ZERO &&
//      partons[emission]->hardProcessMass() == ZERO &&
//      partons[spectator]->hardProcessMass() == ZERO);
    partons[emitter]->hardProcessMass() == ZERO &&
    partons[emission]->hardProcessMass() == ZERO &&
    partons[spectator]->hardProcessMass() != ZERO;
}

double IFMqqxDipole::me2Avg(double ccme2) const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  double muj2 = sqr( (realEmissionME()->lastXComb().mePartonData()[realSpectator()]->hardProcessMass()) ) /
    (2.* (realEmissionME()->lastXComb().meMomenta()[bornSpectator()])*
     (realEmissionME()->lastXComb().meMomenta()[realEmitter()]) );

  double res = x + 2.*(1.-x)/x -
    2.*muj2/x*u/(1.-u);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= -ccme2;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

double IFMqqxDipole::me2() const {

  if ( jacobian() == 0.0 )
    return 0.0;

  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

  Energy2 prop = 
    2.*((realEmissionME()->lastXComb().meMomenta()[realEmitter()])*
	(realEmissionME()->lastXComb().meMomenta()[realEmission()]))*x;

  Lorentz5Momentum pc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]/u -
    realEmissionME()->lastXComb().meMomenta()[realSpectator()]/(1.-u);

  Energy2 sc = 
    realEmissionME()->lastXComb().meMomenta()[realEmission()]*
    realEmissionME()->lastXComb().meMomenta()[realSpectator()];
  sc /= u*(1.-u)*(1.-x)/x;

  SpinCorrelationTensor corr(-x,pc,sc/2.);

  double res = -underlyingBornME()->spinColourCorrelatedME2(make_pair(bornEmitter(),bornSpectator()),
							    corr);

  double CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());

  res *= 8.*CF*Constants::pi*(realEmissionME()->lastXComb().lastSHat())*
    (underlyingBornME()->lastXComb().lastAlphaS())/prop;

  res *= 
    pow(realEmissionME()->lastXComb().lastSHat() / underlyingBornME()->lastXComb().lastSHat(),
	underlyingBornME()->lastXComb().mePartonData().size()-4.);

  res *=
    realEmissionME()->finalStateSymmetry() /
    underlyingBornME()->finalStateSymmetry();

  return res;

}

void IFMqqxDipole::persistentOutput(PersistentOStream &) const {
}

void IFMqqxDipole::persistentInput(PersistentIStream &, int) {
}

void IFMqqxDipole::Init() {

  static ClassDocumentation<IFMqqxDipole> documentation
    ("IFMqqxDipole");

//  DipoleRepository::registerDipole<0,IFMqqxDipole,IFLightTildeKinematics,IFLightInvertedTildeKinematics>
//    ("IFMqqxDipole","IFLightTildeKinematics","IFLightInvertedTildeKinematics");
  DipoleRepository::registerDipole<0,IFMqqxDipole,IFMassiveTildeKinematics,IFMassiveInvertedTildeKinematics>
    ("IFMqqxDipole","IFMassiveTildeKinematics","IFMassiveInvertedTildeKinematics");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMqqxDipole,SubtractionDipole>
describeHerwigIFMqqxDipole("Herwig::IFMqqxDipole", "Herwig.so");
#line 1 "./SubtractionDipole.cc"
// -*- C++ -*-
//
// SubtractionDipole.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SubtractionDipole class.
//

#include "SubtractionDipole.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDF/PartonBin.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/TildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/InvertedTildeKinematics.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"
#include "Herwig/MatrixElement/Matchbox/Utility/DiagramDrawer.h"

#include <iterator>
using std::ostream_iterator;

using namespace Herwig;

SubtractionDipole::SubtractionDipole() 
  : MEBase(), theSplitting(false), theApply(true), theSubtractionTest(false),
    theIgnoreCuts(false),
    theRealEmitter(-1), theRealEmission(-1), theRealSpectator(-1), 
    lastRealEmissionKey(realEmissionKey(cPDVector(),-1,-1,-1)),
    lastUnderlyingBornKey(underlyingBornKey(cPDVector(),-1,-1)),
    theBornEmitter(-1), theBornSpectator(-1),
    theLastSubtractionScale(ZERO), theLastSplittingScale(ZERO),
    theLastSubtractionPt(ZERO), theLastSplittingPt(ZERO),
    theLastSubtractionZ(0.0), theLastSplittingZ(0.0),
    theRealShowerSubtraction(false), theVirtualShowerSubtraction(false),
    theLoopSimSubtraction(false), theRealEmissionScales(false),
    theShowerHardScale(ZERO), theShowerScale(ZERO), 
    theIsInShowerPhasespace(false), theIsAboveCutoff(false) {}

SubtractionDipole::~SubtractionDipole() {}

double SubtractionDipole::alpha() const{
  return factory()->alphaParameter();
}

void SubtractionDipole::clearBookkeeping() {
  theRealEmitter = -1;
  theRealEmission = -1;
  theRealSpectator = -1;
  theBornEmitter = -1;
  theBornSpectator = -1;
  theMergingMap.clear();
  theSplittingMap.clear();
  theIndexMap.clear();
  theUnderlyingBornDiagrams.clear();
  theRealEmissionDiagrams.clear();
  theBornToRealDiagrams.clear();
  theRealToBornDiagrams.clear();
}

void SubtractionDipole::setupBookkeeping(const map<Ptr<DiagramBase>::ptr,SubtractionDipole::MergeInfo>& mergeInfo,bool slim) {

  theMergingMap.clear();
  theSplittingMap.clear();

  theUnderlyingBornDiagrams.clear();
  theRealEmissionDiagrams.clear();

  theBornToRealDiagrams.clear();
  theRealToBornDiagrams.clear();

  int xemitter = -1;
  int xspectator = -1;
  map<int,int> mergeLegs;
  map<int,int> remapLegs;
  map<int,int> realBornMap;
  map<int,int> bornRealMap;

  set<Ptr<DiagramBase>::cptr> usedDiagrams;

  for ( map<Ptr<DiagramBase>::ptr,MergeInfo>::const_iterator mit = mergeInfo.begin();
	mit != mergeInfo.end(); ++mit ) {

    DiagramVector::const_iterator bd = 
      theUnderlyingBornME->diagrams().end();

    // work out the most similar underlying Born diagram
    map<int,int> xRemapLegs;
    int nomapScore = 0;
    for ( DiagramVector::const_iterator b = 
	    theUnderlyingBornME->diagrams().begin();
	  b != theUnderlyingBornME->diagrams().end(); ++b ) {
      map<int,int> theRemapLegs;
      if ( mit->second.diagram->isSame(*b,theRemapLegs) &&
	   usedDiagrams.find(*b) == usedDiagrams.end() ) {
	int theNomapScore = 0;
	for ( map<int,int>::const_iterator m = theRemapLegs.begin();
	      m != theRemapLegs.end(); ++m )
	  if ( m->first == m->second )
	    theNomapScore += 1;
	if ( theNomapScore >= nomapScore ) {
	  nomapScore = theNomapScore;
	  xRemapLegs = theRemapLegs;
	  bd = b;
	}
      }
    }

    // no underlying Born
    if ( bd == theUnderlyingBornME->diagrams().end() )
      continue;

    // as we deal with one splitting only we now mark this diagram as used
    // since we fixed the overall remapping of the process from the first
    // occurence, see below. TODO: This confuses this code even more, and
    // clearly calls for a cleanup. This is just grown historically and got
    // messed up with experiencing different processes and setups.
    usedDiagrams.insert(*bd);

    if ( xemitter == -1 ) {

      xemitter = mit->second.emitter;
      mergeLegs = mit->second.mergeLegs;
      remapLegs = xRemapLegs;

      assert(remapLegs.find(xemitter) != remapLegs.end());
      xemitter = remapLegs[xemitter];

      // work out the leg remapping real -> born
      for ( map<int,int>::const_iterator k = mergeLegs.begin();
	    k != mergeLegs.end(); ++k ) {
	assert(remapLegs.find(k->second) != remapLegs.end());
	realBornMap[k->first] = remapLegs[k->second];
      }

      // work out the leg remapping born -> real
      for ( map<int,int>::const_iterator k = realBornMap.begin();
	    k != realBornMap.end(); ++k ) {
	bornRealMap[k->second] = k->first;
      }

      // work out the spectator
      assert(mergeLegs.find(realSpectator()) != mergeLegs.end());
      assert(remapLegs.find(mergeLegs[realSpectator()]) != remapLegs.end());
      xspectator = realBornMap[realSpectator()];

    }

    RealEmissionKey realKey = realEmissionKey((*mit->first).partons(),realEmitter(),realEmission(),realSpectator());
    UnderlyingBornKey bornKey = underlyingBornKey((**bd).partons(),xemitter,xspectator);
    if ( theMergingMap.find(realKey) == theMergingMap.end() )
      theMergingMap.insert(make_pair(realKey,make_pair(bornKey,realBornMap)));
    RealEmissionInfo realInfo = make_pair(realKey,bornRealMap);
    bool gotit = false;
    typedef multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator spIterator;
    pair<spIterator,spIterator> range = theSplittingMap.equal_range(bornKey);
    for ( ; range.first != range.second; ++range.first )
      if ( range.first->second == realInfo ) {
	gotit = true;
	break;
      }
    if ( !gotit )
      theSplittingMap.insert(make_pair(bornKey,realInfo));
    theUnderlyingBornDiagrams[process(realKey)].push_back(*bd);
    theRealEmissionDiagrams[process(bornKey)].push_back(mit->first);
    theBornToRealDiagrams[*bd] = mit->first;
    theRealToBornDiagrams[mit->first] = *bd;

  }
  
  
  if (slim) {
    theIndexMap.clear();
    theSplittingMap.clear();
    theBornToRealDiagrams.clear();
    theRealEmissionDiagrams.clear();
  }

  if ( theSplittingMap.empty() )
    return;

  theIndexMap.clear();

  for ( multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator s =
	  theSplittingMap.begin(); s != theSplittingMap.end(); ++s ) {
    theIndexMap[process(s->first)] = make_pair(emitter(s->first),spectator(s->first));
  }

}

void SubtractionDipole::subtractionBookkeeping() {
  /*
  if ( theMergingMap.empty() )
    setupBookkeeping();
  */
  assert(!theMergingMap.empty());
  lastRealEmissionKey = 
    realEmissionKey(lastHeadXComb().mePartonData(),realEmitter(),realEmission(),realSpectator());
  map<RealEmissionKey,UnderlyingBornInfo>::const_iterator k =
    theMergingMap.find(lastRealEmissionKey);
  if ( k == theMergingMap.end() ) {
    theApply = false;
    return;
  }
  theApply = true;
  lastUnderlyingBornKey = k->second.first;
  bornEmitter(emitter(lastUnderlyingBornKey));
  bornSpectator(spectator(lastUnderlyingBornKey));
}

void SubtractionDipole::splittingBookkeeping() {
  /*
  if ( theMergingMap.empty() )
    setupBookkeeping();
  */
  assert(!theMergingMap.empty());
  map<cPDVector,pair<int,int> >::const_iterator esit =
    theIndexMap.find(lastHeadXComb().mePartonData());
  if ( esit == theIndexMap.end() ) {
    theApply = false;
    return;
  }
  theApply = true;
  pair<int,int> es = esit->second;
  bornEmitter(es.first);
  bornSpectator(es.second);
  lastUnderlyingBornKey = underlyingBornKey(lastHeadXComb().mePartonData(),bornEmitter(),bornSpectator());
  typedef multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator spit;
  pair<spit,spit> kr = theSplittingMap.equal_range(lastUnderlyingBornKey);
  assert(kr.first != kr.second);
  lastRealEmissionInfo = kr.first;
  for ( ; lastRealEmissionInfo != kr.second; ++lastRealEmissionInfo )
    if ( process(lastRealEmissionInfo->second.first) == lastXComb().mePartonData() )
      break;
  assert(lastRealEmissionInfo != kr.second);
  lastRealEmissionKey = lastRealEmissionInfo->second.first;
  realEmitter(emitter(lastRealEmissionKey));
  realEmission(emission(lastRealEmissionKey));
  realSpectator(spectator(lastRealEmissionKey));
}

StdXCombPtr SubtractionDipole::makeXComb(Energy newMaxEnergy, const cPDPair & inc,
					 tEHPtr newEventHandler,tSubHdlPtr newSubProcessHandler,
					 tPExtrPtr newExtractor,	tCascHdlPtr newCKKW,
					 const PBPair & newPartonBins, tCutsPtr newCuts,
					 const DiagramVector & newDiagrams, bool mir,
					 const PartonPairVec& allBins,
					 tStdXCombPtr newHead,
					 tMEPtr newME) {

  if ( !newME )
    newME = this;

  if ( !splitting() ) {
    return
      underlyingBornME()->makeXComb(newMaxEnergy, inc,
				    newEventHandler, newSubProcessHandler,
				    newExtractor, newCKKW,
				    newPartonBins, newCuts,
				    newDiagrams, mir, allBins,
				    newHead, newME);
  }

  return
    realEmissionME()->makeXComb(newMaxEnergy, inc,
				newEventHandler, newSubProcessHandler,
				newExtractor, newCKKW,
				newPartonBins, newCuts,
				newDiagrams, mir, allBins,
				newHead, newME);

}

StdXCombPtr SubtractionDipole::makeXComb(tStdXCombPtr newHead,
					 const PBPair & newPartonBins,
					 const DiagramVector & newDiagrams,
					 tMEPtr newME) {

  if ( !newME )
    newME = this;

  if ( !splitting() ) {
    return
      underlyingBornME()->makeXComb(newHead, newPartonBins,
				    newDiagrams, newME);
  }

  return
    realEmissionME()->makeXComb(newHead, newPartonBins,
				newDiagrams, newME);

}

StdXCombPtr SubtractionDipole::makeBornXComb(tStdXCombPtr realXC) {

  const cPDVector& proc = const_cast<const StandardXComb&>(*realXC).mePartonData();

  lastRealEmissionKey = 
    realEmissionKey(proc,realEmitter(),realEmission(),realSpectator());
  map<RealEmissionKey,UnderlyingBornInfo>::const_iterator k =
    theMergingMap.find(lastRealEmissionKey);

  if ( k == theMergingMap.end() )
    return StdXCombPtr();

  PartonPairVec pbs = realXC->pExtractor()->getPartons(realXC->maxEnergy(), 
						       realXC->particles(),
						       *(realXC->cuts()));

  DiagramVector bornDiags = underlyingBornDiagrams(proc);
  assert(!bornDiags.empty());

  PartonPairVec::iterator ppit = pbs.begin();
  for ( ; ppit != pbs.end(); ++ppit ) {
    if ( ppit->first->parton() == bornDiags.front()->partons()[0] &&
	 ppit->second->parton() == bornDiags.front()->partons()[1] )
      break;
  }

  assert(ppit != pbs.end());

  return
    underlyingBornME()->makeXComb(realXC,*ppit,bornDiags,this);

}

vector<StdXCombPtr> SubtractionDipole::makeRealXCombs(tStdXCombPtr bornXC) {

  const cPDVector& proc = const_cast<const StandardXComb&>(*bornXC).mePartonData();

  map<cPDVector,pair<int,int> >::const_iterator esit = theIndexMap.find(proc);
  if ( esit == theIndexMap.end() ) 
    return vector<StdXCombPtr>();
  pair<int,int> es = esit->second;
  bornEmitter(es.first);
  bornSpectator(es.second);
  lastUnderlyingBornKey = underlyingBornKey(proc,bornEmitter(),bornSpectator());

  if ( theSplittingMap.find(lastUnderlyingBornKey) == theSplittingMap.end() )
    return vector<StdXCombPtr>();

  PartonPairVec pbs = bornXC->pExtractor()->getPartons(bornXC->maxEnergy(), 
						       bornXC->particles(),
						       *(bornXC->cuts()));

  DiagramVector realDiags = realEmissionDiagrams(proc);
  assert(!realDiags.empty());

  vector<StdXCombPtr> res;

  map<cPDVector,DiagramVector> realProcs;

  for ( MEBase::DiagramVector::const_iterator d = realDiags.begin();
	d != realDiags.end(); ++d ) {
    realProcs[(**d).partons()].push_back(*d);
  }

  for ( map<cPDVector,DiagramVector>::const_iterator pr =
	  realProcs.begin(); pr != realProcs.end(); ++pr ) {

    PartonPairVec::iterator ppit = pbs.begin();
    for ( ; ppit != pbs.end(); ++ppit ) {
      if ( ppit->first->parton() == pr->second.front()->partons()[0] &&
	   ppit->second->parton() == pr->second.front()->partons()[1] )
	break;
    }

    assert(ppit != pbs.end());

    StdXCombPtr rxc =
      realEmissionME()->makeXComb(bornXC,*ppit,pr->second,this);

    res.push_back(rxc);

  }

  return res;

}

const MEBase::DiagramVector& SubtractionDipole::underlyingBornDiagrams(const cPDVector& real) const {
  static DiagramVector empty;
  map<cPDVector,DiagramVector>::const_iterator k = theUnderlyingBornDiagrams.find(real);
  if (k == theUnderlyingBornDiagrams.end() )
    return empty;
  return k->second;
}

tcDiagPtr SubtractionDipole::underlyingBornDiagram(tcDiagPtr realDiag) const {
  map<tcDiagPtr,tcDiagPtr>::const_iterator it = theRealToBornDiagrams.find(realDiag);
  assert(it != theRealToBornDiagrams.end());
  return it->second;
}

const MEBase::DiagramVector& SubtractionDipole::realEmissionDiagrams(const cPDVector& born) const {
  static DiagramVector empty;
  map<cPDVector,DiagramVector>::const_iterator k = theRealEmissionDiagrams.find(born);
  if ( k == theRealEmissionDiagrams.end() )
    return empty;
  return k->second;
}

tcDiagPtr SubtractionDipole::realEmissionDiagram(tcDiagPtr bornDiag) const {
  map<tcDiagPtr,tcDiagPtr>::const_iterator it = theBornToRealDiagrams.find(bornDiag);
  assert(it != theBornToRealDiagrams.end());
  return it->second;
}

void SubtractionDipole::getDiagrams() const {
  if ( splitting() ) {
    realEmissionME()->diagrams();
    useDiagrams(realEmissionME());
  } else {
    underlyingBornME()->diagrams();
    useDiagrams(underlyingBornME());
  }
}

Selector<MEBase::DiagramIndex> SubtractionDipole::diagrams(const DiagramVector & dv) const {
  Ptr<MatchboxMEBase>::tcptr me = 
    splitting() ?
    realEmissionME() :
    underlyingBornME();
  if ( me->phasespace() ) {
    me->phasespace()->setXComb(lastXCombPtr());
    me->phasespace()->clearDiagramWeights();
    me->phasespace()->fillDiagramWeights();
  }
  return 
    me->diagrams(dv);
}

Selector<const ColourLines *>
SubtractionDipole::colourGeometries(tcDiagPtr diag) const {
  return 
    splitting() ?
    realEmissionME()->colourGeometries(diag) :
    underlyingBornME()->colourGeometries(diag);
}

const ColourLines &
SubtractionDipole::selectColourGeometry(tcDiagPtr diag) const {
  return 
    splitting() ?
    realEmissionME()->selectColourGeometry(diag) :
    underlyingBornME()->selectColourGeometry(diag);
}

void SubtractionDipole::flushCaches() {
  theUnderlyingBornME->flushCaches();
  theRealEmissionME->flushCaches();
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator r =
	  reweights().begin(); r != reweights().end(); ++r ) {
    (**r).flushCaches();
  }
}

void SubtractionDipole::setXComb(tStdXCombPtr xc) {
  if ( !xc ) {
    theApply = false;
    return;
  } else {
    theApply = true;
  }
  lastMatchboxXComb(xc);
  MEBase::setXComb(xc); 
  if ( splitting() ) {
    realEmissionME()->setXComb(xc);
    underlyingBornME()->setXComb(xc->head());
    splittingBookkeeping();
  } else {
    realEmissionME()->setXComb(xc->head());
    underlyingBornME()->setXComb(xc);
    subtractionBookkeeping();
  }
  if ( !apply() )
    return;
}

void SubtractionDipole::setKinematics() {
  MEBase::setKinematics(); 
  if ( splitting() )
    realEmissionME()->setKinematics();
  else
    underlyingBornME()->setKinematics();
}

bool SubtractionDipole::generateKinematics(const double * r) {
  if ( lastXCombPtr()->kinematicsGenerated() )
    return true;
  if ( splitting() ) {
    if ( !generateRadiationKinematics(r) )
      return false;
    if( ! realEmissionME()->lastXCombPtr()->setIncomingPartons())
      return false;
    realEmissionME()->setScale();
    double jac = jacobian();
    jac *= pow(underlyingBornME()->lastXComb().lastSHat() / realEmissionME()->lastXComb().lastSHat(),
	       realEmissionME()->lastXComb().mePartonData().size()-4.);
    jacobian(jac);
    assert(lastXCombPtr() == realEmissionME()->lastXCombPtr());
    lastXCombPtr()->didGenerateKinematics();
    return true;
  }
  if ( !generateTildeKinematics() ){ return false;}
  if( ! underlyingBornME()->lastXCombPtr()->setIncomingPartons() )
    return false;
  underlyingBornME()->setScale();
  assert(lastXCombPtr() == underlyingBornME()->lastXCombPtr());
  if( ! underlyingBornME()->lastXCombPtr()->setIncomingPartons() )
      return false;
  // need to have the scale and x's available for checking shower phase space
  if ( showerApproximation() &&
       lastXCombPtr()->willPassCuts() )
    showerApproximation()->getShowerVariables();
  lastXCombPtr()->didGenerateKinematics();
  return true;
}

int SubtractionDipole::nDim() const {
  if ( !splitting() )
    return underlyingBornME()->nDim();
  return underlyingBornME()->nDim() + nDimRadiation();
}

void SubtractionDipole::clearKinematics() {
  MEBase::clearKinematics(); 
  if ( splitting() )
    realEmissionME()->clearKinematics();
  else
    underlyingBornME()->clearKinematics();
}

void SubtractionDipole::tildeKinematics(Ptr<TildeKinematics>::tptr tk) { 
  theTildeKinematics = tk;
}

bool SubtractionDipole::generateTildeKinematics() {

  assert(!splitting());

  Ptr<TildeKinematics>::tptr kinematics = theTildeKinematics;
  if ( showerApproximation() ) {
    showerApproximation()->setBornXComb(lastXCombPtr());
    showerApproximation()->setRealXComb(realEmissionME()->lastXCombPtr());
    showerApproximation()->setDipole(this);
    showerApproximation()->checkCutoff();
    if ( showerApproximation()->showerTildeKinematics() &&
	 isAboveCutoff() &&
	 realShowerSubtraction() )
      kinematics = showerApproximation()->showerTildeKinematics();
  }

  if ( !kinematics ) {
    jacobian(0.0);
    return false;
  }

  kinematics->prepare(lastHeadXCombPtr(),lastXCombPtr());

  if ( !kinematics->doMap() ) {
    jacobian(0.0);
    return false;
  }

  theLastSubtractionScale = kinematics->lastScale();
  theLastSubtractionPt = kinematics->lastPt();
  theLastSubtractionZ = kinematics->lastZ();

  meMomenta().resize(lastHeadXComb().meMomenta().size() - 1);

  assert(mergingMap().find(lastRealEmissionKey) != mergingMap().end());
  map<int,int>& trans = theMergingMap[lastRealEmissionKey].second;

  int n = lastHeadXComb().meMomenta().size();
  for ( int k = 0; k < n; ++k ) {
    if ( k == realEmitter() || k == realEmission() || k == realSpectator() )
      continue;
    meMomenta()[trans[k]] = lastHeadXComb().meMomenta()[k];
    if ( kinematics->doesTransform() && k > 1 )
      meMomenta()[trans[k]] = kinematics->transform(meMomenta()[trans[k]]);
  }

  meMomenta()[bornEmitter()] = 
    const_cast<const TildeKinematics&>(*kinematics).bornEmitterMomentum();
  meMomenta()[bornSpectator()] = 
    const_cast<const TildeKinematics&>(*kinematics).bornSpectatorMomentum();

  cPDVector::const_iterator pd = mePartonData().begin();
  vector<Lorentz5Momentum>::iterator p = meMomenta().begin();
  for ( ; pd != mePartonData().end(); ++pd, ++p ) {
    p->setMass((**pd).hardProcessMass());
    p->rescaleRho();
  }

  jacobian(realEmissionME()->lastXComb().jacobian());

  logGenerateTildeKinematics();

  return true;

}

void SubtractionDipole::invertedTildeKinematics(Ptr<InvertedTildeKinematics>::tptr itk) { 
  theInvertedTildeKinematics = itk;
}

int SubtractionDipole::nDimRadiation() const {
  return invertedTildeKinematics() ? 
    invertedTildeKinematics()->nDimRadiation() :
    0;
}

bool SubtractionDipole::generateRadiationKinematics(const double * r) {

  assert(splitting());

  Ptr<InvertedTildeKinematics>::tptr kinematics = theInvertedTildeKinematics;
  if ( showerApproximation() ) {
    showerApproximation()->setBornXComb(lastHeadXCombPtr());
    showerApproximation()->setRealXComb(lastXCombPtr());
    showerApproximation()->setDipole(this);
    if ( showerApproximation()->showerInvertedTildeKinematics() ) {
      kinematics = showerApproximation()->showerInvertedTildeKinematics();
    }
  }

  if ( !kinematics ) {
    jacobian(0.0);
    return false;
  }

  kinematics->prepare(lastXCombPtr(),lastHeadXCombPtr());

  if ( !kinematics->doMap(r) ) {
    jacobian(0.0);
    return false;
  }

  theLastSplittingScale = kinematics->lastScale();
  theLastSplittingPt = kinematics->lastPt();
  theLastSplittingZ = kinematics->lastZ();

  meMomenta().resize(lastHeadXComb().meMomenta().size() + 1);

  assert(splittingMap().find(lastUnderlyingBornKey) != splittingMap().end());
  map<int,int>& trans = const_cast<map<int,int>&>(lastRealEmissionInfo->second.second);

  int n = lastHeadXComb().meMomenta().size();
  for ( int k = 0; k < n; ++k ) {
    if ( k == bornEmitter() || k == bornSpectator() )
      continue;
    meMomenta()[trans[k]] = lastHeadXComb().meMomenta()[k];
    if ( kinematics->doesTransform() && k > 1 )
      meMomenta()[trans[k]] = kinematics->transform(meMomenta()[trans[k]]);
  }

  meMomenta()[realEmitter()] = 
    const_cast<const InvertedTildeKinematics&>(*kinematics).realEmitterMomentum();
  meMomenta()[realEmission()] = 
    const_cast<const InvertedTildeKinematics&>(*kinematics).realEmissionMomentum();
  meMomenta()[realSpectator()] = 
    const_cast<const InvertedTildeKinematics&>(*kinematics).realSpectatorMomentum();

  cPDVector::const_iterator pd = mePartonData().begin();
  vector<Lorentz5Momentum>::iterator p = meMomenta().begin();
  for ( ; pd != mePartonData().end(); ++pd, ++p ) {
    p->setMass((**pd).hardProcessMass());
    p->rescaleRho();
  }

  jacobian(underlyingBornME()->lastXComb().jacobian() *
	   kinematics->jacobian());

  logGenerateRadiationKinematics(r);

  return true;

}

void SubtractionDipole::ptCut(Energy cut) {
  theInvertedTildeKinematics->ptCut(cut);
}

CrossSection SubtractionDipole::dSigHatDR(Energy2 factorizationScale) const {

  double pdfweight = 1.;

  double jac = jacobian();

  if ( splitting() && jac == 0.0 ) {
    lastMECrossSection(ZERO);
    return ZERO;
  }

  if ( factorizationScale == ZERO ) {
    factorizationScale = underlyingBornME()->lastScale();
  }

  if ( havePDFWeight1() ) {
    pdfweight *= realEmissionME()->pdf1(factorizationScale);
  }
 
  if ( havePDFWeight2() ) {
    pdfweight *= realEmissionME()->pdf2(factorizationScale);
  }

  lastMEPDFWeight(pdfweight);

  bool needTheDipole = true;
  CrossSection shower = ZERO;

  double lastThetaMu = 1.0;

  double showerFactor = 1.;

  if ( showerApproximation() ) {
    assert(!splitting());
    showerApproximation()->setBornXComb(lastXCombPtr());
    showerApproximation()->setRealXComb(realEmissionME()->lastXCombPtr());
    showerApproximation()->setDipole(const_cast<SubtractionDipole*>(this));
    if ( !isAboveCutoff() ) {
      showerApproximation()->wasBelowCutoff();
      lastThetaMu = 0.0;
    } else {
      lastThetaMu = 1.0;
    }
    if ( lastThetaMu > 0.0 && isInShowerPhasespace() ) {
      if ( realShowerSubtraction() )
	shower = showerApproximation()->dSigHatDR()*lastThetaMu;
      if ( virtualShowerSubtraction() || loopSimSubtraction() )
	shower = -showerApproximation()->dSigHatDR()*lastThetaMu;
      if ( virtualShowerSubtraction() &&
	   isAboveCutoff() &&
	   showerApproximation()->showerTildeKinematics() ) {
	// map shower to dipole kinematics; we are always above the
	// cutoff in this case
	showerFactor *= 
	  showerApproximation()->showerTildeKinematics()->jacobianRatio();
      }
      shower *= showerFactor;
    }
    if ( realShowerSubtraction() && lastThetaMu == 1.0 )
      needTheDipole = false;
    if ( virtualShowerSubtraction() && lastThetaMu == 0.0 )
      needTheDipole = false;
    if ( factory()->loopSimCorrections() ||
	 factory()->meCorrectionsOnly() )
      needTheDipole = false;
  }

  double xme2 = 0.0;

  if ( needTheDipole )
    xme2 = me2();

  if ( factory()->loopSimCorrections() ||
       factory()->meCorrectionsOnly() ) {

    assert(showerApproximation());
    xme2 = realEmissionME()->me2() * showerApproximation()->channelWeight();

    double rws =
      pow(underlyingBornME()->lastXComb().lastAlphaS()/
	  realEmissionME()->lastXComb().lastAlphaS(),
	  realEmissionME()->orderInAlphaS());

    xme2 *= rws;

    double rwe =
      pow(underlyingBornME()->lastXComb().lastAlphaEM()/
	  realEmissionME()->lastXComb().lastAlphaEM(),
	  underlyingBornME()->orderInAlphaEW());

    xme2 *= rwe;

  }

  if ( realShowerSubtraction() )
    xme2 *= 1. - lastThetaMu;
  if ( virtualShowerSubtraction() || loopSimSubtraction() )
    xme2 *= lastThetaMu;

  double coupl = lastMECouplings();
  coupl *= underlyingBornME()->lastXComb().lastAlphaS();
  lastMECouplings(coupl);

  CrossSection res = 
    sqr(hbarc) * jac * pdfweight * xme2 /
    (2. * realEmissionME()->lastXComb().lastSHat());

  if ( !showerApproximation() && xme2 != 0.0 ) {
    double weight = 0.0;
    bool applied = false;
    for ( vector<Ptr<MatchboxReweightBase>::ptr>::const_iterator rw =
	    theReweights.begin(); rw != theReweights.end(); ++rw ) {
      (**rw).setXComb(theRealEmissionME->lastXCombPtr());
      if ( !(**rw).apply() )
	continue;
      weight += (**rw).evaluate();
      applied = true;
    }
    if ( applied )
      res *= weight;
  }

  lastMECrossSection(-res-shower);

  logDSigHatDR(jac);

  return lastMECrossSection();
}


bool SubtractionDipole::aboveAlpha() const{return theTildeKinematics->aboveAlpha();}



CrossSection SubtractionDipole::prefactor(Energy2 factorizationScale)const{
  
  const double jac = jacobian();
  assert( factorizationScale != ZERO );
  assert (! splitting());
  double pdfweight = 1.;
  if ( havePDFWeight1() ) pdfweight *= realEmissionME()->pdf1(factorizationScale);
  if ( havePDFWeight2() ) pdfweight *= realEmissionME()->pdf2(factorizationScale);

  
  return sqr(hbarc) * jac * pdfweight /  (2. * realEmissionME()->lastXComb().lastSHat());
  
}




CrossSection SubtractionDipole::ps(Energy2 factorizationScale,Ptr<ColourBasis>::tptr largeNBasis) const {

  double ccme2 =underlyingBornME()->me2()*
                underlyingBornME()->
                largeNColourCorrelatedME2(
                  make_pair(bornEmitter(),bornSpectator()),largeNBasis)/
                underlyingBornME()->largeNME2(largeNBasis);
  
  return prefactor(factorizationScale) * me2Avg(ccme2);
}




pair<CrossSection,CrossSection> SubtractionDipole::dipandPs(Energy2 factorizationScale,Ptr<ColourBasis>::tptr largeNBasis) const {
  
  CrossSection  factor= prefactor(factorizationScale);

  double ccme2 =underlyingBornME()->me2()*
                underlyingBornME()->
                largeNColourCorrelatedME2(
                              make_pair(bornEmitter(),bornSpectator()),largeNBasis)/
                underlyingBornME()->largeNME2(largeNBasis);

  double ps = me2Avg(ccme2);
  double dip = me2();
  
  return make_pair(factor*dip,factor*ps);
}

CrossSection SubtractionDipole::dip(Energy2 factorizationScale) const {
  CrossSection  factor= prefactor(factorizationScale);
  double dip = me2();
  return factor*dip;
}


void SubtractionDipole::print(ostream& os) const {

  os << "--- SubtractionDipole setup ----------------------------------------------------\n";

  os << " subtraction '" << name() << "'\n for real emission '"
     << theRealEmissionME->name() << "'\n using underlying Born '"
     << theUnderlyingBornME->name() << "'\n";

  os << " tilde kinematics are '"
     << (theTildeKinematics ? theTildeKinematics->name() : "") 
     << " '\n inverted tilde kinematics are '"
     << (theInvertedTildeKinematics ? theInvertedTildeKinematics->name() : "") << "'\n";

  os << " the following subtraction mappings have been found:\n";

  for ( map<RealEmissionKey,UnderlyingBornInfo>::const_iterator m =
	  theMergingMap.begin(); m != theMergingMap.end(); ++m ) {
    os << " " << process(m->second.first)[0]->PDGName() << " "
       << process(m->second.first)[1]->PDGName() << " -> ";
    for ( cPDVector::const_iterator p = process(m->second.first).begin() + 2;
	  p != process(m->second.first).end(); ++p ) {
      os << (**p).PDGName() << " ";
    }
    os << "[" << emitter(m->second.first) << "," << spectator(m->second.first) << "] <=> ";
    os << process(m->first)[0]->PDGName() << " "
       << process(m->first)[1]->PDGName() << " -> ";
    for ( cPDVector::const_iterator p = process(m->first).begin() + 2;
	  p != process(m->first).end(); ++p ) {
      os << (**p).PDGName() << " ";
    }
    os << "[(" << emitter(m->first) << "," << emission(m->first) << ")," << spectator(m->first) << "]\n"
       << " non-dipole momenta ( ";
    for ( map<int,int>::const_iterator k = m->second.second.begin();
	  k != m->second.second.end(); ++k ) {
      if ( k->first == spectator(m->first) )
	continue;
      os << k->second << " ";
    }
    os << ") <=> ( ";
    for ( map<int,int>::const_iterator k = m->second.second.begin();
	  k != m->second.second.end(); ++k ) {
      if ( k->first == spectator(m->first) )
	continue;
      os << k->first << " ";
    }
    os << ")\n";
  }

  os << "--------------------------------------------------------------------------------\n";

  os << flush;

}

void SubtractionDipole::printLastEvent(ostream& os) const {

  os << "--- SubtractionDipole last event information -----------------------------------\n";

  os << " for dipole '" << name() << "' applying [" 
     << bornEmitter() << "," << bornSpectator() << "] <=> [("
     << realEmitter() << "," << realEmission() << ")," << realSpectator() << "]\n"
     << " evaluated the cross section/nb " << (lastMECrossSection()/nanobarn) << "\n"
     << " with subtraction parameters x[0] = " << subtractionParameters()[0]
     << " x[1] = " << subtractionParameters()[1] << "\n";

  os << " the last real emission event was:\n";
  realEmissionME()->printLastEvent(os);

  os << " the last underlying Born event was:\n";
  underlyingBornME()->printLastEvent(os);


  os << "--- end SubtractionDipole last event information -------------------------------\n";

  os << flush;

}

void SubtractionDipole::logME2() const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  tcStdXCombPtr bornxc = splitting() ? lastHeadXCombPtr() : lastXCombPtr();
  tcStdXCombPtr realxc = splitting() ? lastXCombPtr() : lastHeadXCombPtr();

  generator()->log() << "'" << name() << "' evaluated me2 using\n"
		     << "Born XComb " << bornxc << " real XComb " << realxc << "\n";

  generator()->log() << "subtraction parameters: ";
  copy(subtractionParameters().begin(),subtractionParameters().end(),
       ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n";

  generator()->log() << "Born phase space point (in GeV):\n";

  vector<Lorentz5Momentum>::const_iterator pit = bornxc->meMomenta().begin();
  cPDVector::const_iterator dit = bornxc->mePartonData().begin();

  for ( ; pit != bornxc->meMomenta().end() ; ++pit, ++dit )
    generator()->log() << (**dit).PDGName() << " : "
		       << (*pit/GeV) << "\n";

  generator()->log() << "with x1 = " << bornxc->lastX1() << " x2 = " << bornxc->lastX2() << "\n"
		     << "sHat/GeV2 = " << (bornxc->lastSHat()/GeV2) << "\n";

  generator()->log() << "Real emission phase space point (in GeV):\n";

  pit = realxc->meMomenta().begin();
  dit = realxc->mePartonData().begin();

  for ( ; pit != realxc->meMomenta().end() ; ++pit, ++dit )
    generator()->log() << (**dit).PDGName() << " : "
		       << (*pit/GeV) << "\n";

  generator()->log() << "with x1 = " << realxc->lastX1() << " x2 = " << realxc->lastX2() << "\n"
		     << "sHat/GeV2 = " << (realxc->lastSHat()/GeV2) << "\n";

}

void SubtractionDipole::logDSigHatDR(double effectiveJac) const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  tcStdXCombPtr bornxc = splitting() ? lastHeadXCombPtr() : lastXCombPtr();
  tcStdXCombPtr realxc = splitting() ? lastXCombPtr() : lastHeadXCombPtr();

  generator()->log() << "'" << name() << "' evaluated cross section using\n"
		     << "Born XComb " << bornxc << " real XComb " << realxc << "\n"
		     << "Jacobian = " << jacobian()
		     << " effective Jacobian = " << effectiveJac << "\n"
		     << "Born sHat/GeV2 = " << (bornxc->lastSHat()/GeV2)
		     << " real sHat/GeV2 = " << (realxc->lastSHat()/GeV2)
		     << " dsig/nb = "
		     << (lastMECrossSection()/nanobarn) << "\n" << flush;

}

void SubtractionDipole::logGenerateTildeKinematics() const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  generator()->log() << "'" << name() << "' generating tilde kinematics.\n"
		     << "configuration: [" << bornEmitter() << ","
		     << bornSpectator() << "] => "
		     << "[(" << realEmitter() << "," << realEmission() << "),"
		     << realSpectator() << "]\n"
		     << "with real xcomb " << lastHeadXCombPtr() << " born xcomb "
		     << lastXCombPtr() << "\n"
		     << "from real emission phase space point:\n";
  Lorentz5Momentum rSum;
  vector<Lorentz5Momentum>::const_iterator pr = lastHeadXComb().meMomenta().begin();
  cPDVector::const_iterator dr = lastHeadXComb().mePartonData().begin();
  size_t count = 0;
  for ( ; pr != lastHeadXComb().meMomenta().end(); ++pr,++dr ) {
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";
    if ( count < 2 ) {
      rSum -= *pr;
    } else {
      rSum += *pr;
    }
    ++count;

  }
  generator()->log() << "sum : " << (rSum/GeV) << "\n";

  generator()->log() << "subtraction parameters: ";
  copy(subtractionParameters().begin(),subtractionParameters().end(),
       ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n"
		     << "with scale/GeV = " << (theLastSubtractionScale/GeV)
		     << "and pt/GeV = " << (theLastSubtractionPt/GeV) << "\n";

  generator()->log() << "generated tilde kinematics:\n";
  pr = lastXComb().meMomenta().begin();
  dr = lastXComb().mePartonData().begin();
  count = 0;
  Lorentz5Momentum bSum;
  for ( ; pr != lastXComb().meMomenta().end(); ++pr,++dr ) {
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";
    if ( count < 2 ) {
      bSum -= *pr;
    } else {
      bSum += *pr;
    }
    ++count;
  }
  generator()->log() << "sum : " << (bSum/GeV) << "\n";

  generator()->log() << "Jacobian = " << jacobian() << "\n" << flush;

}


void SubtractionDipole::logGenerateRadiationKinematics(const double * r) const {

  if ( !realEmissionME()->verbose() &&
       !underlyingBornME()->verbose() )
    return;

  generator()->log() << "'" << name() << "' generating radiation kinematics.\n"
		     << "configuration: [" << bornEmitter() << ","
		     << bornSpectator() << "] => "
		     << "[(" << realEmitter() << "," << realEmission() << "),"
		     << realSpectator() << "]\n"
		     << "with born xcomb " << lastHeadXCombPtr() << " real xcomb "
		     << lastXCombPtr() << "\n"
		     << "from random numbers:\n";
  copy(r,r+nDimRadiation(),ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n";
  generator()->log() << "and born phase space point:\n";
  vector<Lorentz5Momentum>::const_iterator pr = lastHeadXComb().meMomenta().begin();
  cPDVector::const_iterator dr = lastHeadXComb().mePartonData().begin();
  for ( ; pr != lastHeadXComb().meMomenta().end(); ++pr,++dr )
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";

  generator()->log() << "subtraction parameters: ";
  copy(subtractionParameters().begin(),subtractionParameters().end(),
       ostream_iterator<double>(generator()->log()," "));
  generator()->log() << "\n" << flush;

  generator()->log() << "scales: scale/GeV = " << (theLastSplittingScale/GeV)
		     << " pt/GeV = " << (theLastSplittingPt/GeV) << "\n" << flush;

  generator()->log() << "generated real emission kinematics:\n";
  pr = lastXComb().meMomenta().begin();
  dr = lastXComb().mePartonData().begin();
  for ( ; pr != lastXComb().meMomenta().end(); ++pr,++dr )
    generator()->log() << (**dr).PDGName() << " : "
		       << (*pr/GeV) << "\n";

  generator()->log() << "Jacobian = "
		     << jacobian() << " = "
		     << underlyingBornME()->lastXComb().jacobian()
		     << "|Born * "
		     << invertedTildeKinematics()->jacobian()
		     << "|Radiation\n" << flush;

}

void SubtractionDipole::doinit() {
  MEBase::doinit();
  if ( underlyingBornME() ) {
    theUnderlyingBornME->init();
  }
  if ( realEmissionME() ) {
    theRealEmissionME->init();
  }
  if ( tildeKinematics() ) {
    theTildeKinematics->init();
  }
  if ( invertedTildeKinematics() ) {
    theInvertedTildeKinematics->init();
  }
  if ( showerApproximation() ) {
    theShowerApproximation->init();
  }
  for ( vector<Ptr<SubtractionDipole>::tptr>::iterator p = thePartners.begin();
	p != thePartners.end(); ++p ) {
    (**p).init();
  }
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    (**rw).init();
  }
}

void SubtractionDipole::doinitrun() {
  MEBase::doinitrun();
  if ( underlyingBornME() ) {
    theUnderlyingBornME->initrun();
  }
  if ( realEmissionME() ) {
    theRealEmissionME->initrun();
  }
  if ( tildeKinematics() ) {
    theTildeKinematics->initrun();
  }
  if ( invertedTildeKinematics() ) {
    theInvertedTildeKinematics->initrun();
  }
  if ( showerApproximation() ) {
    theShowerApproximation->initrun();
  }
  for ( vector<Ptr<SubtractionDipole>::tptr>::iterator p = thePartners.begin();
	p != thePartners.end(); ++p ) {
    (**p).initrun();
  }
  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    (**rw).initrun();
  }
}

void SubtractionDipole::cloneDependencies(const std::string& prefix,bool slim) {

  if ( underlyingBornME() ) {
    Ptr<MatchboxMEBase>::ptr myUnderlyingBornME = underlyingBornME()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myUnderlyingBornME->name();
    if ( ! (generator()->preinitRegister(myUnderlyingBornME,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Matrix element " << pname.str() << " already existing." << Exception::runerror;
    myUnderlyingBornME->cloneDependencies(pname.str(),slim);
    underlyingBornME(myUnderlyingBornME);
  }

  if ( realEmissionME()&& !slim ) {
    Ptr<MatchboxMEBase>::ptr myRealEmissionME = realEmissionME()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myRealEmissionME->name();
    if ( ! (generator()->preinitRegister(myRealEmissionME,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Matrix element " << pname.str() << " already existing." << Exception::runerror;
    myRealEmissionME->cloneDependencies(pname.str());
    realEmissionME(myRealEmissionME);
  }

  if ( tildeKinematics() ) {
    Ptr<TildeKinematics>::ptr myTildeKinematics = tildeKinematics()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myTildeKinematics->name();
    if ( ! (generator()->preinitRegister(myTildeKinematics,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Tilde kinematics " << pname.str() << " already existing." << Exception::runerror;
    myTildeKinematics->dipole(this);
    tildeKinematics(myTildeKinematics);
  }

  if ( invertedTildeKinematics()&& !slim ) {
    Ptr<InvertedTildeKinematics>::ptr myInvertedTildeKinematics = invertedTildeKinematics()->cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << myInvertedTildeKinematics->name();
    if ( ! (generator()->preinitRegister(myInvertedTildeKinematics,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Inverted tilde kinematics " << pname.str() << " already existing." << Exception::runerror;
    myInvertedTildeKinematics->dipole(this);
    invertedTildeKinematics(myInvertedTildeKinematics);
  }

  for ( vector<Ptr<MatchboxReweightBase>::ptr>::iterator rw =
	  theReweights.begin(); rw != theReweights.end(); ++rw ) {
    Ptr<MatchboxReweightBase>::ptr myReweight = (**rw).cloneMe();
    ostringstream pname;
    pname << (prefix == "" ? fullName() : prefix) << "/" << (**rw).name();
    if ( ! (generator()->preinitRegister(myReweight,pname.str()) ) )
      throw Exception() << "SubtractionDipole::cloneDependencies(): Reweight " << pname.str() << " already existing." << Exception::runerror;
    myReweight->cloneDependencies(pname.str());
    *rw = myReweight;
  }

}

void SubtractionDipole::constructVertex(tSubProPtr sub) {
  if ( splitting() )
    realEmissionME()->constructVertex(sub);
  else 
    underlyingBornME()->constructVertex(sub);
}

void SubtractionDipole::constructVertex(tSubProPtr sub, const ColourLines* cl) {
  if ( splitting() )
    realEmissionME()->constructVertex(sub,cl);
  else 
    underlyingBornME()->constructVertex(sub,cl);
}

void SubtractionDipole::generateSubCollision(SubProcess & sub) {
  if ( splitting() )
    realEmissionME()->generateSubCollision(sub);
  else
    underlyingBornME()->generateSubCollision(sub);
}

void SubtractionDipole::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theSplitting << theApply << theSubtractionTest 
     << theIgnoreCuts << theRealEmissionME << theUnderlyingBornME 
     << thePartners << theTildeKinematics << theInvertedTildeKinematics 
     << theReweights << theRealEmitter << theRealEmission << theRealSpectator 
     << theSubtractionParameters << theMergingMap << theSplittingMap 
     << theIndexMap << theUnderlyingBornDiagrams << theRealEmissionDiagrams 
     << theBornToRealDiagrams << theRealToBornDiagrams
     << lastRealEmissionKey << lastUnderlyingBornKey 
     << theBornEmitter << theBornSpectator << ounit(theLastSubtractionScale,GeV) 
     << ounit(theLastSplittingScale,GeV) << ounit(theLastSubtractionPt,GeV) 
     << ounit(theLastSplittingPt,GeV) << theLastSubtractionZ
     << theLastSplittingZ << theShowerApproximation 
     << theRealShowerSubtraction << theVirtualShowerSubtraction 
     << theLoopSimSubtraction << theRealEmissionScales
     << ounit(theShowerHardScale,GeV) << ounit(theShowerScale,GeV) 
     << theShowerParameters << theIsInShowerPhasespace << theIsAboveCutoff;
}

void SubtractionDipole::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theSplitting >> theApply >> theSubtractionTest 
     >> theIgnoreCuts >> theRealEmissionME >> theUnderlyingBornME 
     >> thePartners >> theTildeKinematics >> theInvertedTildeKinematics 
     >> theReweights >> theRealEmitter >> theRealEmission >> theRealSpectator 
     >> theSubtractionParameters >> theMergingMap >> theSplittingMap 
     >> theIndexMap >> theUnderlyingBornDiagrams >> theRealEmissionDiagrams 
     >> theBornToRealDiagrams >> theRealToBornDiagrams
     >> lastRealEmissionKey >> lastUnderlyingBornKey 
     >> theBornEmitter >> theBornSpectator >> iunit(theLastSubtractionScale,GeV) 
     >> iunit(theLastSplittingScale,GeV) >> iunit(theLastSubtractionPt,GeV) 
     >> iunit(theLastSplittingPt,GeV) >> theLastSubtractionZ
     >> theLastSplittingZ >> theShowerApproximation 
     >> theRealShowerSubtraction >> theVirtualShowerSubtraction 
     >> theLoopSimSubtraction >> theRealEmissionScales
     >> iunit(theShowerHardScale,GeV) >> iunit(theShowerScale,GeV) 
     >> theShowerParameters >> theIsInShowerPhasespace >> theIsAboveCutoff;
  lastMatchboxXComb(theLastXComb);
  typedef multimap<UnderlyingBornKey,RealEmissionInfo>::const_iterator spit;
  pair<spit,spit> kr = theSplittingMap.equal_range(lastUnderlyingBornKey);
  lastRealEmissionInfo = kr.first;
  for ( ; lastRealEmissionInfo != kr.second; ++lastRealEmissionInfo )
    if ( process(lastRealEmissionInfo->second.first) == lastXComb().mePartonData() )
      break;
}

Ptr<MatchboxFactory>::tptr SubtractionDipole::factory() const {
  return MatchboxFactory::currentFactory();
}

void SubtractionDipole::Init() {

  static ClassDocumentation<SubtractionDipole> documentation
    ("SubtractionDipole represents a dipole subtraction "
     "term in the formalism of Catani and Seymour.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<SubtractionDipole,MEBase>
describeSubtractionDipole("Herwig::SubtractionDipole", "Herwig.so");
