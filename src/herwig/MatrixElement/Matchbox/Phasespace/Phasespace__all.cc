#line 1 "./FFLightInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// FFLightInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFLightInvertedTildeKinematics class.
//

#include "FFLightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FFLightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFLightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFLightInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;

  double z = ptz.second;
  double y = sqr(pt/lastScale())/(z*(1.-z));
  if ( y < 0. || y > 1. ||
       z < 0. || z > 1. ) {
    jacobian(0.0);
    return false;
  }

  mapping *= (1.-y)/(z*(1.-z));
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt
    = getKt(emitter,spectator,pt,phi);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;

  realEmitterMomentum() = z*emitter + y*(1.-z)*spectator + kt;
  realEmissionMomentum() = (1.-z)*emitter + y*z*spectator - kt;
  realSpectatorMomentum() = (1.-y)*spectator;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FFLightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(y*z*(1.-z));

}

double FFLightInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

Energy FFLightInvertedTildeKinematics::ptMax() const {
  return lastScale()/2.;
}

pair<double,double> FFLightInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFLightInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFLightInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFLightInvertedTildeKinematics::Init() {

  static ClassDocumentation<FFLightInvertedTildeKinematics> documentation
    ("FFLightInvertedTildeKinematics inverts the final-final tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFLightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFFLightInvertedTildeKinematics("Herwig::FFLightInvertedTildeKinematics", "Herwig.so");
#line 1 "./FFLightTildeKinematics.cc"
// -*- C++ -*-
//
// FFLightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFLightTildeKinematics class.
//

#include "FFLightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FFLightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFLightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FFLightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double y = emission*emitter / (emission*emitter + emission*spectator + emitter*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = z;

  bornEmitterMomentum() = emitter+emission-(y/(1.-y))*spectator;
  bornSpectatorMomentum() = spectator/(1.-y);

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FFLightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(y*z*(1.-z));

}


Energy FFLightTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {
  Energy scale =  (emitter+emission+spectator).m();
  double y = emission*emitter/(emission*emitter + emission*spectator + emitter*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);
  Energy ret = scale * sqrt( y  * z*(1.-z) );
  return ret;
}

pair<double,double> FFLightTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}


double FFLightTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFLightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFLightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFLightTildeKinematics::Init() {

  static ClassDocumentation<FFLightTildeKinematics> documentation
    ("FFLightTildeKinematics implements the 'tilde' kinematics for "
     "a final-final subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFLightTildeKinematics,TildeKinematics>
describeHerwigFFLightTildeKinematics("Herwig::FFLightTildeKinematics", "Herwig.so");
#line 1 "./FFMassiveInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// FFMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMassiveInvertedTildeKinematics class.
//

#include "FFMassiveInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

IBPtr FFMassiveInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

// Matches Stephen Webster's thesis
bool FFMassiveInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  vector<double> values(6);
  pair<Energy,double> ptz = generatePtZ(mapping, r, &values);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  Energy2 pt2 = sqr(pt);
  double z = ptz.second;

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / lastScale() );
  double muj2 = sqr( realEmissionData()->hardProcessMass() / lastScale() );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / lastScale() );
  double Muij2 = sqr( bornEmitterData()->hardProcessMass() / lastScale() );
  double Muk2 = sqr( bornSpectatorData()->hardProcessMass() / lastScale() );

  // Define the scale
  Energy2 Qijk = sqr(lastScale());

  // Most of the required values have been calculated in ptzAllowed
  double y = values[0];
  double zi = values[1];
  double xk = values[2];
  double xij = values[3];
  double QijN2 = values[4];
  double sijkN = values[5];
  double sijkN2 = sqr(sijkN);

  // Construct reference momenta nk, nij
  Lorentz5Momentum nij = ( sijkN2 / (sijkN2-Muij2*Muk2) )
    * (emitter - (Muij2/sijkN)*spectator);
  Lorentz5Momentum nk = ( sijkN2 / (sijkN2-Muij2*Muk2) )
    * (spectator - (Muk2/sijkN)*emitter);

  // Construct qt
  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum qt = getKt(emitter,spectator,pt,phi);

  // Construct qij, qk, qi and qj
  Lorentz5Momentum qij = xij*nij + (Muij2/(xij*sijkN))*nk;
  Lorentz5Momentum spe = xk*nk + (Muk2/(xk*sijkN))*nij;

  Lorentz5Momentum em = z*qij
    + ((pt2/Qijk + mui2 - z*z*Muij2)/(xij*sijkN*z))*nk + qt;
  Lorentz5Momentum emm = (1.-z)*qij
    + ((pt2/Qijk + muj2 - sqr(1.-z)*Muij2)/(xij*sijkN*(1.-z)))*nk - qt;

  em.setMass(realEmitterData()->hardProcessMass());
  em.rescaleEnergy();
  emm.setMass(realEmissionData()->hardProcessMass());
  emm.rescaleEnergy();
  spe.setMass(realSpectatorData()->hardProcessMass());
  spe.rescaleEnergy();

  // book
  realEmitterMomentum() = em;
  realEmissionMomentum() = emm;
  realSpectatorMomentum() = spe;

  // Calculate the jacobian
  double bar = 1.-mui2-muj2-muk2;
  
  // mapFactor defined as dy dz = mapFactor * dpt2/sqr(lastScale()) dz
  double mapFactor = 0.0;
  mapFactor = y*(sqr(lastScale()) / (pt2 + sqr(1.-z)*mui2*Qijk + sqr(z)*muj2*Qijk))
      * abs(1. - 2.*Muk2*QijN2 / (bar*(1.-y)*xij*xk*sijkN));
  
  // Mapping includes only the variable changes/jacobians
  mapping *= mapFactor;
  jacobian( (Qijk*sqr(bar)/rootOfKallen(1.,Muij2,Muk2)) * mapping
            * (1.-y)/(16.*sqr(Constants::pi)) / sHat() );
  
  // Store the parameters
  subtractionParameters().resize(3);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = zi;
  subtractionParameters()[2] = z;
  
  return true;
}


Energy FFMassiveInvertedTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[2];
  
  Energy ret = scale * sqrt( y * (1.-mui2-muj2-muk2) * z*(1.-z) - sqr(1.-z)*mui2 - sqr(z)*muj2 );

  return ret;
}

double FFMassiveInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[2];
}    

Energy FFMassiveInvertedTildeKinematics::ptMax() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  Energy ptmax = rootOfKallen( mui2, muj2, sqr(1.-sqrt(muk2)) ) /
    ( 2.-2.*sqrt(muk2) ) * scale;
  
  return ptmax > 0.*GeV ? ptmax : 0.*GeV;
}

// NOTE: bounds calculated at this step may be too loose
pair<double,double> FFMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {

  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double zp = ( 1.+mui2-muj2+muk2-2.*sqrt(muk2) +
    rootOfKallen(mui2,muj2,sqr(1.-sqrt(muk2))) *
    sqrt( 1.-sqr(pt/hardPt) ) ) /
    ( 2.*sqr(1.-sqrt(muk2)) );
  double zm = ( 1.+mui2-muj2+muk2-2.*sqrt(muk2) -
    rootOfKallen(mui2,muj2,sqr(1.-sqrt(muk2))) *
    sqrt( 1.-sqr(pt/hardPt) ) ) /
    ( 2.*sqr(1.-sqrt(muk2)) );
    
  return make_pair(zm,zp);
}

// Matches Stephen Webster's thesis
bool FFMassiveInvertedTildeKinematics::ptzAllowed(pair<Energy,double> ptz, vector<double>* values) const {

  Energy pt = ptz.first;
  Energy2 pt2 = sqr(pt);
  double z = ptz.second;

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / lastScale() );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / lastScale() );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / lastScale() );
  double Muij2 = sqr( bornEmitterData()->hardProcessMass() / lastScale() );
  double Muk2 = sqr( bornSpectatorData()->hardProcessMass() / lastScale() );
  
  // Calculate the scales that we need
  Energy2 Qijk = sqr(lastScale());
  double QijN2 = (pt2/Qijk + (1.-z)*mui2 + z*muj2) / z / (1.-z);
  double sijkN = 0.5*( 1. - Muij2 - Muk2 + sqrt( sqr(1.-Muij2-Muk2) - 4.*Muij2*Muk2 ) );
  double bar = 1.-mui2-muj2-muk2;

  double y = ( pt2/Qijk + sqr(1.-z)*mui2 + z*z*muj2 ) /
      (z*(1.-z)*bar);

  // Calculate the scaling factors, xk and xij
  double lambdaK = 1. + (Muk2/sijkN);
  double lambdaIJ = 1. + (Muij2/sijkN);
  double fac1 = lambdaIJ*lambdaK + (muk2 - QijN2)/sijkN;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*muk2/sijkN ) )
    / 2. / lambdaK ;
  double xij = 1. - muk2*(1.-xk) / xk / sijkN;

  // Calculate zi
  double zi =
    ( z*xij*xk*sijkN + muk2*(pt2/Qijk + mui2) / (z*xij*xk*sijkN) )
    / (1.-y) / bar;
  
  // Limits on zi
  double facA = (2.*mui2+bar*y)/2./(mui2 + muj2 + bar*y);
  double facB =
    sqrt( (sqr(2.*muk2 + bar*(1.-y)) - 4.*muk2) *
          (sqr(bar)*sqr(y) - 4.*mui2*muj2))
    / bar / (1.-y) / (bar*y + 2.*mui2);
  double zim = facA * (1. - facB);
  double zip = facA * (1. + facB);
  
  // check (y,z) phase space boundary
  double ym = 2.*sqrt(mui2)*sqrt(muj2)/bar;
  double yp = 1. - 2.*sqrt(muk2)*(1.-sqrt(muk2))/bar;

  if ( y<ym || y>yp || zi<zim || zi>zip ) return false;

  assert( (*values).size() == 6);
  (*values)[0] = y;
  (*values)[1] = zi;
  (*values)[2] = xk;
  (*values)[3] = xij;
  (*values)[4] = QijN2;
  (*values)[5] = sijkN;

  return true;
}


// This is used to generate pt and z
pair<Energy,double> FFMassiveInvertedTildeKinematics::generatePtZ(double& jac, const double * r, vector<double> * values) const {

  double kappaMin = 
    ptCut() != ZERO ?
    sqr(ptCut()/ptMax()) :
    sqr(0.1*GeV/GeV);

  double kappa;

  using namespace RandomHelpers;

  if ( ptCut() > ZERO ) {
    pair<double,double> kw =
      generate(inverse(0.,kappaMin,1.),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  } else {
    pair<double,double> kw =
      generate((piecewise(),
		flat(1e-4,kappaMin),
		match(inverse(0.,kappaMin,1.))),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  }

  Energy pt = sqrt(kappa)*ptMax();

  pair<double,double> zLims = zBounds(pt);

  pair<double,double> zw =
    generate(inverse(0.,zLims.first,zLims.second)+
	     inverse(1.,zLims.first,zLims.second),r[1]);

  double z = zw.first;
  jac *= zw.second;

  jac *= sqr(ptMax()/lastScale());
  
  if( !ptzAllowed(make_pair(pt,z), values )) jac = 0.;

  return make_pair(pt,z);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMassiveInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFMassiveInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {  
}

void FFMassiveInvertedTildeKinematics::Init() {

  static ClassDocumentation<FFMassiveInvertedTildeKinematics> documentation
    ("FFMassiveInvertedTildeKinematics inverts the final-final tilde "
     "kinematics involving a massive particle.");
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFFMassiveInvertedTildeKinematics("Herwig::FFMassiveInvertedTildeKinematics", "Herwig.so");
#line 1 "./FFMassiveTildeKinematics.cc"
// -*- C++ -*-
//
// FFMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFMassiveTildeKinematics class.
//

#include "FFMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FFMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FFMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

// Matches Stephen Webster's thesis
bool FFMassiveTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  // Compute y
  double y = emission*emitter / (emission*emitter + emission*spectator + emitter*spectator);

  // Calculate the scale
  Lorentz5Momentum pTot = emitter+emission+spectator;
  Energy scale = pTot.m();

  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  double Muij2 = sqr( bornEmitterData()->hardProcessMass() / scale );
  double Muk2 = sqr( bornSpectatorData()->hardProcessMass() / scale );
 
  // Calculate the invariants
  Energy2 Qijk = sqr(scale);
  double sijkN = 0.5*( 1. - Muij2 - Muk2 + sqrt( sqr(1.-Muij2-Muk2) - 4.*Muij2*Muk2 ) );
  double bar = 1. - mui2 - muj2 - muk2;
  
  // Calculate Qij2
  double QijN2 = sqr(emitter + emission)/Qijk;
  
  // Calculate the scale factors, xk and xij
  double lambdaK = 1. + (Muk2/sijkN);
  double lambdaIJ = 1. + (Muij2/sijkN);
  double fac1 = lambdaIJ*lambdaK + (muk2 - QijN2)/sijkN;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*muk2/sijkN ) )
    / 2. / lambdaK ;
  double xij = 1. - muk2*(1.-xk) / xk / sijkN;
  
  // Calculate z = qi.nk / (qi+qj).nk from qi.qk and y
  double l = Muk2 / (2.*xk*xij*sijkN);
  double a = (xij*xk*sijkN/2.) - l*(bar*y + mui2 + muj2);
  double b = (-emitter*spectator)/Qijk + l*(bar*y + 2.*mui2);
  double z = -b/a;

  // Calculate zi
  double zi = emitter*spectator / ((emitter + emission)*spectator);

  // Store the variables
  subtractionParameters().resize(3);
  subtractionParameters()[0] = y;
  subtractionParameters()[1] = zi;
  subtractionParameters()[2] = z;
  
  // Calculate nij and nk from qi, q, qj
  double L = QijN2/xij/xk/sijkN;
  Lorentz5Momentum nij = (emitter + emission - L*spectator)/(xij - L*Muk2/xk/sijkN);
  Lorentz5Momentum nk = (1./xk)*(spectator - Muk2*nij/(xk*sijkN));
    
  // Calculate the born momenta from nij and nk
  bornSpectatorMomentum() = nk + (Muk2/sijkN)*nij;
  bornEmitterMomentum() = pTot - bornSpectatorMomentum();

  bornEmitterMomentum().setMass( bornEmitterData()->hardProcessMass() );
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass( bornSpectatorData()->hardProcessMass() );
  bornSpectatorMomentum().rescaleEnergy();
  
  return true;

}

Energy FFMassiveTildeKinematics::lastPt() const {
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
  
  double y = subtractionParameters()[0];
  double z = subtractionParameters()[2];
  
  // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );

  Energy ret = scale * sqrt( y * (1.-mui2-muj2-muk2) * z*(1.-z) - sqr(1.-z)*mui2 - sqr(z)*muj2 );
  
  return ret;
  
}

// Matches Stephen Webster's thesis
Energy FFMassiveTildeKinematics::lastPt(Lorentz5Momentum emitter,
					Lorentz5Momentum emission,
					Lorentz5Momentum spectator)const {


  // Compute y
  double y = emission*emitter / (emission*emitter 
			       + emission*spectator 
			       + emitter*spectator);

  // Calculate the scale
  Lorentz5Momentum pTot = emitter+emission+spectator;
  Energy scale = pTot.m();

  // masses
  double mui2 = sqr( emitter.mass() / scale );
  double muj2  = sqr( emission.mass() / scale );
  double muk2 = sqr( spectator.mass() / scale );
   // TODO: here we assume a gluon
  bool isgluon= emitter.mass()==emission.mass();
  double Muij2 = sqr(( isgluon?ZERO:max(emission.mass(),emitter.mass()) )/ scale );
  double Muk2 = sqr( spectator.mass() / scale );
 
  // Calculate the invariants
  Energy2 Qijk = sqr(scale);
  double sijkN = 0.5*( 1. - Muij2 - Muk2 + sqrt( sqr(1.-Muij2-Muk2) - 4.*Muij2*Muk2 ) );
  double bar = 1. - mui2 - muj2 - muk2;
  
  // Calculate Qij2
  double QijN2 = sqr(emitter + emission)/Qijk;
  
  // Calculate the scale factors, xk and xij
  double lambdaK = 1. + (Muk2/sijkN);
  double lambdaIJ = 1. + (Muij2/sijkN);
  double fac1 = lambdaIJ*lambdaK + (muk2 - QijN2)/sijkN;
  double xk =
    ( fac1 + sqrt( sqr(fac1) - 4.*lambdaIJ*lambdaK*muk2/sijkN ) )
    / 2. / lambdaK ;
  double xij = 1. - muk2*(1.-xk) / xk / sijkN;
  
  // Calculate z = qi.nk / (qi+qj).nk from qi.qk and y
  double l = Muk2 / (2.*xk*xij*sijkN);
  double a = (xij*xk*sijkN/2.) - l*(bar*y + mui2 + muj2);
  double b = (-emitter*spectator)/Qijk + l*(bar*y + 2.*mui2);
  double z = -b/a;

  Energy ret = scale * sqrt( z*(1.-z)*QijN2 - (1.-z)*mui2 - z*muj2 );
    
  return ret;
  
}


  // NOTE: bounds calculated at this step may be too loose
pair<double,double> FFMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  
  if(pt>hardPt) return make_pair(0.5,0.5);
  
  Energy scale = (bornEmitterMomentum()+bornSpectatorMomentum()).m();
    // masses
  double mui2 = sqr( realEmitterData()->hardProcessMass() / scale );
  double muj2  = sqr( realEmissionData()->hardProcessMass() / scale );
  double muk2 = sqr( realSpectatorData()->hardProcessMass() / scale );
  
  double zp = ( 1.+mui2-muj2+muk2-2.*sqrt(muk2) +
               rootOfKallen(mui2,muj2,sqr(1.-sqrt(muk2))) *
               sqrt( 1.-sqr(pt/hardPt) ) ) /
  ( 2.*sqr(1.-sqrt(muk2)) );
  double zm = ( 1.+mui2-muj2+muk2-2.*sqrt(muk2) -
               rootOfKallen(mui2,muj2,sqr(1.-sqrt(muk2))) *
               sqrt( 1.-sqr(pt/hardPt) ) ) /
  ( 2.*sqr(1.-sqrt(muk2)) );
  
  return make_pair(zm,zp);
}



double FFMassiveTildeKinematics::lastZ() const {
  return subtractionParameters()[2];
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FFMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FFMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FFMassiveTildeKinematics::Init() {

  static ClassDocumentation<FFMassiveTildeKinematics> documentation
    ("FFMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a final-final subtraction dipole involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FFMassiveTildeKinematics,TildeKinematics>
describeHerwigFFMassiveTildeKinematics("Herwig::FFMassiveTildeKinematics", "Herwig.so");
#line 1 "./FILightInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// FILightInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FILightInvertedTildeKinematics class.
//

#include "FILightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FILightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FILightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FILightInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;

  double z = ptz.second;
  double y = sqr(pt/lastScale())/(z*(1.-z));
  double x = 1./(1.+y);

  if ( x < spectatorX() || x > 1. ||
       z < 0. || z > 1. ) {
    jacobian(0.0);
    return false;
  }

  // This should (and does) have a factor of 1/x relative to
  // the dipole shower jacobian. 
  mapping /= z*(1.-z);
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt
    = getKt(spectator,emitter,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  realEmitterMomentum() = z*emitter + (1.-z)*((1.-x)/x)*spectator + kt;
  realEmissionMomentum() = (1.-z)*emitter + z*((1.-x)/x)*spectator - kt;
  realSpectatorMomentum() = (1./x)*spectator;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FILightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(z*(1.-z)*(1.-x)/x);

}

double FILightInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

Energy FILightInvertedTildeKinematics::ptMax() const {
  double x = spectatorX();
  return sqrt((1.-x)/x)*lastScale()/2.;
}

pair<double,double> FILightInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}

void FILightInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FILightInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FILightInvertedTildeKinematics::Init() {

  static ClassDocumentation<FILightInvertedTildeKinematics> documentation
    ("FILightInvertedTildeKinematics inverts the final-initial tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FILightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFILightInvertedTildeKinematics("Herwig::FILightInvertedTildeKinematics", "Herwig.so");
#line 1 "./FILightTildeKinematics.cc"
// -*- C++ -*-
//
// FILightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FILightTildeKinematics class.
//

#include "FILightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FILightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FILightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FILightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = 
    (- emission*emitter + emission*spectator + emitter*spectator) / 
    (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  bornEmitterMomentum() = emitter+emission-(1.-x)*spectator;
  bornSpectatorMomentum() = x*spectator;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FILightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return scale * sqrt(z*(1.-z)*(1.-x)/x);

}

Energy FILightTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {

  
  Energy2 scale =  - (emitter+emission-spectator).m2();
  
  
  double x =
  (- emission*emitter + emission*spectator + emitter*spectator) /
  (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);
  
  return sqrt( z*(1.-z)*(1.-x)/x*scale ) ;
}

pair<double,double> FILightTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*(1.-s),0.5*(1.+s));
}



double FILightTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FILightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FILightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FILightTildeKinematics::Init() {

  static ClassDocumentation<FILightTildeKinematics> documentation
    ("FILightTildeKinematics implements the 'tilde' kinematics for "
     "a final-initial subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FILightTildeKinematics,TildeKinematics>
describeHerwigFILightTildeKinematics("Herwig::FILightTildeKinematics", "Herwig.so");
#line 1 "./FIMassiveInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// FIMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMassiveInvertedTildeKinematics class.
//

#include "FIMassiveInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FIMassiveInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FIMassiveInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FIMassiveInvertedTildeKinematics::doMap(const double * r) { 
  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;

  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());

  Energy2 scale=2.*emitter*spectator;
  double y = (pt*pt+(1.-z)*mi2+z*m2-z*(1.-z)*Mi2) / (z*(1.-z)*scale);
  double x = 1./(1.+y);

  if ( x < spectatorX() ) {
    jacobian(0.0);
    return false;
  }

  // SW 05/12/2016: Checked this is correct
  // This should appear to have a factor of 1/x relative
  // to the dipole shower jacobian. It is cancelled by ratios of
  // real and born cross sections in the units.
  mapping /= z*(1.-z);
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(spectator,emitter,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  realEmitterMomentum() = z*emitter +
    (sqr(pt)+mi2-z*z*Mi2)/(z*scale)*spectator + kt;
  realEmissionMomentum() = (1.-z)*emitter +
    (pt*pt+m2-sqr(1.-z)*Mi2)/((1.-z)*scale)*spectator - kt;
  realSpectatorMomentum() = (1.+y)*spectator;

  double mui2 = x*mi2/scale;
  double mu2  = x*m2/scale;
  double Mui2 = x*Mi2/scale;
  double xp = 1. + Mui2 - sqr(sqrt(mui2)+sqrt(mu2));
  double zm = .5*( 1.-x+Mui2+mui2-mu2 -
  		   sqrt( sqr(1.-x+Mui2-mui2-mu2)-4.*mui2*mu2 ) ) / (1.-x+Mui2);
  double zp = .5*( 1.-x+Mui2+mui2-mu2 +
  		   sqrt( sqr(1.-x+Mui2-mui2-mu2)-4.*mui2*mu2 ) ) / (1.-x+Mui2);

  if ( x > xp || z < zm || z > zp ) {
    jacobian(0.0);
    return false;
  }
  realEmitterMomentum().setMass(sqrt(mi2));
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(sqrt(m2));
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();
  return true;

}

Energy FIMassiveInvertedTildeKinematics::lastPt() const {
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());

  Energy2 scale = Mi2 - (realEmitterMomentum()+realEmissionMomentum()-realSpectatorMomentum()).m2();
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];

  return sqrt( z*(1.-z)*(1.-x)/x*scale -
	       ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) );

}

double FIMassiveInvertedTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

Energy FIMassiveInvertedTildeKinematics::ptMax() const {
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());
  double x = spectatorX();
  // s^star/x
  Energy2 scale=2.*bornEmitterMomentum()*bornSpectatorMomentum();
  Energy2 s = scale * (1.-x)/x + Mi2;
  Energy ptmax = .5 * sqrt(s) * rootOfKallen( s/s, mi2/s, m2/s );
  return ptmax > 0.*GeV ? ptmax : 0.*GeV;
}

pair<double,double> FIMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());
  // s^star/x
  Energy2 scale=2.*bornEmitterMomentum()*bornSpectatorMomentum();
  Energy2 s = scale * (1.-spectatorX())/spectatorX() +  Mi2;

  double zm = .5*( 1.+(mi2-m2)/s - rootOfKallen(s/s,mi2/s,m2/s) *
		    sqrt( 1.-sqr(pt/hardPt) ) );
  double zp = .5*( 1.+(mi2-m2)/s + rootOfKallen(s/s,mi2/s,m2/s) *
		    sqrt( 1.-sqr(pt/hardPt) ) );
  return make_pair(zm, zp);

}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMassiveInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FIMassiveInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FIMassiveInvertedTildeKinematics::Init() {

  static ClassDocumentation<FIMassiveInvertedTildeKinematics> documentation
    ("FIMassiveInvertedTildeKinematics inverts the final-initial tilde "
     "kinematics involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigFIMassiveInvertedTildeKinematics("Herwig::FIMassiveInvertedTildeKinematics", "Herwig.so");
#line 1 "./FIMassiveTildeKinematics.cc"
// -*- C++ -*-
//
// FIMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FIMassiveTildeKinematics class.
//

#include "FIMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FIMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr FIMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool FIMassiveTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());

  double x = 
    (- emission*emitter + emission*spectator + emitter*spectator +
     0.5*(Mi2-mi2-m2)) / 
    (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = z;

  bornEmitterMomentum() = emitter+emission-(1.-x)*spectator;
  bornSpectatorMomentum() = x*spectator;

  bornEmitterMomentum().setMass(bornEmitterData()->hardProcessMass());
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(bornSpectatorData()->hardProcessMass());
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy FIMassiveTildeKinematics::lastPt() const {

  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());

  Energy2 scale = Mi2 - (realEmitterMomentum()+realEmissionMomentum()-realSpectatorMomentum()).m2();
  double x = subtractionParameters()[0];
  double z = subtractionParameters()[1];
  return sqrt( z*(1.-z)*(1.-x)/x*scale -
	       ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) );

}



Energy FIMassiveTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {
    // g->QQ or Q -> Qg
  Energy2 Mi2 = emitter.m()==emission.m()?0.*GeV2:max(emitter.m2(),emission.m2());
  Energy2 mi2 = emitter.m2();
  Energy2 m2  = emission.m2();
  
  Energy2 scale = Mi2 - (emitter+emission-spectator).m2();
  
  double x =
  (- emission*emitter + emission*spectator + emitter*spectator +
   0.5*(Mi2-mi2-m2)) /
  (emitter*spectator + emission*spectator);
  double z = emitter*spectator / (emitter*spectator + emission*spectator);
  
  return sqrt( z*(1.-z)*(1.-x)/x*scale -
              ((1.-z)*mi2+z*m2-z*(1.-z)*Mi2) );
}


pair<double,double> FIMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  Energy2 mi2 = sqr(realEmitterData()->hardProcessMass());
  Energy2 m2  = sqr(realEmissionData()->hardProcessMass());
  Energy2 Mi2 = sqr(bornEmitterData()->hardProcessMass());
    // s^star/x
  Energy2 scale=2.*bornEmitterMomentum()*bornSpectatorMomentum();
  Energy2 s = scale * (1.-spectatorX())/spectatorX() +  Mi2;
  
  double zm = .5*( 1.+(mi2-m2)/s - rootOfKallen(s/s,mi2/s,m2/s) *
                  sqrt( 1.-sqr(pt/hardPt) ) );
  double zp = .5*( 1.+(mi2-m2)/s + rootOfKallen(s/s,mi2/s,m2/s) *
                  sqrt( 1.-sqr(pt/hardPt) ) );
  return make_pair(zm, zp);
  
}



double FIMassiveTildeKinematics::lastZ() const {
  return subtractionParameters()[1];
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FIMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void FIMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void FIMassiveTildeKinematics::Init() {

  static ClassDocumentation<FIMassiveTildeKinematics> documentation
    ("FIMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a final-initial subtraction dipole involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FIMassiveTildeKinematics,TildeKinematics>
describeHerwigFIMassiveTildeKinematics("Herwig::FIMassiveTildeKinematics", "Herwig.so");
#line 1 "./IFLightInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// IFLightInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFLightInvertedTildeKinematics class.
//

#include "IFLightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr IFLightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFLightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFLightInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;

  double ratio = sqr(pt/lastScale());
  double rho = 1. - 4.*ratio*z*(1.-z) / sqr(1. - z + ratio);
  if ( rho < 0. ) {
    jacobian(0.0);
    return false;
  }

  double x = 0.5*(1./ratio)*(1.-z+ratio)*(1.-sqrt(rho));
  double u = 0.5*(1./(1.-z))*(1.-z+ratio)*(1.-sqrt(rho));

  if ( x < emitterX() || x > 1. || 
       u < 0. || u > 1. ) {
    jacobian(0.0);
    return false;
  }

  // This jacobian is (1/x^2)*dx*du
  mapping *= (1.-x)/((1.-z)*(z*(1.-z)+sqr(x-z)));
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi,true);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;

  realEmitterMomentum() = (1./x)*emitter;
  realEmissionMomentum() = ((1.-x)*(1.-u)/x)*emitter + u*spectator + kt;
  realSpectatorMomentum() = ((1.-x)*u/x)*emitter + (1.-u)*spectator - kt;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFLightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return scale * sqrt(u*(1.-u)*(1.-x)/x);

}

Energy IFLightInvertedTildeKinematics::ptMax() const {
  double x = emitterX();
  return sqrt((1.-x)/x)*lastScale()/2.;
}

pair<double,double> IFLightInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double x = emitterX();
  return make_pair(0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s));
}

double IFLightInvertedTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return 1. - (1.-x)*(1.-u);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFLightInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFLightInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFLightInvertedTildeKinematics::Init() {

  static ClassDocumentation<IFLightInvertedTildeKinematics> documentation
    ("IFLightInvertedTildeKinematics inverts the initial-final tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFLightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIFLightInvertedTildeKinematics("Herwig::IFLightInvertedTildeKinematics", "Herwig.so");
#line 1 "./IFLightTildeKinematics.cc"
// -*- C++ -*-
//
// IFLightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFLightTildeKinematics class.
//

#include "IFLightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr IFLightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFLightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFLightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = 
    (- emission*spectator + emitter*spectator + emitter*emission) / 
    (emitter*emission + emitter*spectator);
  double u = emitter*emission / (emitter*emission + emitter*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;

  bornEmitterMomentum() = x*emitter;
  bornSpectatorMomentum() = spectator + emission - (1.-x)*emitter;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFLightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return scale * sqrt(u*(1.-u)*(1.-x)/x);

}

Energy IFLightTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {

  
  double x =
  (- emission*spectator + emitter*spectator + emitter*emission) /
  (emitter*emission + emitter*spectator);
  double u = emitter*emission / (emitter*emission + emitter*spectator);
  
  Energy scale = sqrt(2.*(emission*emitter-emission*spectator+emitter*spectator));
  
  return scale * sqrt(u*(1.-u)*(1.-x)/x);
}

pair<double,double> IFLightTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double x = emitterX();
  return make_pair(0.5*(1.+x-(1.-x)*s),0.5*(1.+x+(1.-x)*s));
}



double IFLightTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
  return 1. - (1.-x)*(1.-u);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFLightTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFLightTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFLightTildeKinematics::Init() {

  static ClassDocumentation<IFLightTildeKinematics> documentation
    ("IFLightTildeKinematics implements the 'tilde' kinematics for "
     "a initial-final subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFLightTildeKinematics,TildeKinematics>
describeHerwigIFLightTildeKinematics("Herwig::IFLightTildeKinematics", "Herwig.so");
#line 1 "./IFMassiveInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// IFMassiveInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMassiveInvertedTildeKinematics class.
//

#include "IFMassiveInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr IFMassiveInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFMassiveInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFMassiveInvertedTildeKinematics::doMap(const double * r) { 
  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  // Compute dipole scale
  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();
  Energy2 scale = 2.*(spectator*emitter);

  // Generate pt and z
  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r);
  if ( mapping == 0.0 ){
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;
 
    // Compute x and u
    double ratio = sqr(pt)/scale;
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    double rho = 1. - 4.*ratio*(1.-muk2)*z*(1.-z)/sqr(1.-z+ratio);

    double x = 0.5*((1.-z+ratio)/(ratio*(1.-muk2))) * (1. - sqrt(rho));
    double u = x*ratio / (1.-z);
      
    // Following Catani-Seymour paper
    double muk2CS = x*muk2;
    double up = (1.-x) /
      ( 1.-x + muk2CS );
      
    if ( x < emitterX() || x > 1. ||
	 u < 0. || u > up ) {
      jacobian(0.0);
      return false;
    }
   

    // Store x and u
    subtractionParameters().resize(2);
    subtractionParameters()[0] = x;
    subtractionParameters()[1] = u;
    
    // jac = sajk*(1./x^2)*dx*du
    // Note - lastScale() is not equal to scale!!!!!!!
    double jac = u/x/(u + x - 2.*u*x*(1.-muk2))*scale/sqr(pt);
    mapping *= jac;
    jacobian( mapping*(sqr(lastScale())/sHat()) / (16.*sqr(Constants::pi)) );

    // Compute the new momenta
    double phi = 2.*Constants::pi*r[2];
    Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi,true);
    
    realEmitterMomentum() = (1./x)*emitter;
    realEmissionMomentum() = ((1.-x)*(1.-u)/x - 2.*u*muk2)*emitter + u*spectator + kt;
    realSpectatorMomentum() = ((1.-x)*u/x + 2.*u*muk2)*emitter + (1.-u)*spectator - kt;

    realEmitterMomentum().setMass(ZERO);
    realEmitterMomentum().rescaleEnergy();
    realEmissionMomentum().setMass(ZERO);
    realEmissionMomentum().rescaleEnergy();
    realSpectatorMomentum().setMass(bornSpectatorData()->hardProcessMass());
    realSpectatorMomentum().rescaleEnergy();
    return true;
    
  }

Energy IFMassiveInvertedTildeKinematics::lastPt() const {
    Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    double x = subtractionParameters()[0];
    double u = subtractionParameters()[1];
    return sqrt(scale * ( u*(1.-u)*(1.-x)/x - u*u*muk2 ));
}

double IFMassiveInvertedTildeKinematics::lastZ() const {
    Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    double x = subtractionParameters()[0];
    double u = subtractionParameters()[1];  
    return u + x + u*x*(muk2-1.);
}

Energy IFMassiveInvertedTildeKinematics::ptMax() const {
    double xe = emitterX();
    Energy2 scale = 2.*(bornEmitterMomentum()*bornSpectatorMomentum());
    Energy2 A = scale*(1.-xe)/xe;
    Energy2 mk2 = sqr(bornSpectatorData()->hardProcessMass());
    return 0.5*A/sqrt(mk2+A);
}

pair<double,double> IFMassiveInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double xe = emitterX();
  return make_pair(0.5*(1.+xe-(1.-xe)*s),0.5*(1.+xe+(1.-xe)*s));
}



// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMassiveInvertedTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFMassiveInvertedTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFMassiveInvertedTildeKinematics::Init() {

  static ClassDocumentation<IFMassiveInvertedTildeKinematics> documentation
    ("IFMassiveInvertedTildeKinematics inverts the initial-final tilde "
     "kinematics involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMassiveInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIFMassiveInvertedTildeKinematics("Herwig::IFMassiveInvertedTildeKinematics", "Herwig.so");
#line 1 "./IFMassiveTildeKinematics.cc"
// -*- C++ -*-
//
// IFMassiveTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IFMassiveTildeKinematics class.
//

#include "IFMassiveTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr IFMassiveTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IFMassiveTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IFMassiveTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = 
    (- emission*spectator + emitter*spectator + emitter*emission) / 
    (emitter*emission + emitter*spectator);
  double u = emitter*emission / (emitter*emission + emitter*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = u;

  bornEmitterMomentum() = x*emitter;
  bornSpectatorMomentum() = spectator + emission - (1.-x)*emitter;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(bornSpectatorData()->hardProcessMass());
  bornSpectatorMomentum().rescaleEnergy();

  return true;

}

Energy IFMassiveTildeKinematics::lastPt() const {
  Energy2 scale = 2.*(realEmissionMomentum()*realEmitterMomentum()
-realEmissionMomentum()*realSpectatorMomentum()
+realEmitterMomentum()*realSpectatorMomentum());
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];

    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    return sqrt(scale * ( u*(1.-u)*(1.-x)/x - u*u*muk2 ));
   }

Energy IFMassiveTildeKinematics::lastPt(Lorentz5Momentum emitter,Lorentz5Momentum emission,Lorentz5Momentum spectator)const {
  Energy2 scale = 2.*(emission*emitter-emission*spectator+emitter*spectator);
  double x = 0.5*scale / (emitter*emission + emitter*spectator);
  double u = emitter*emission / (emitter*emission + emitter*spectator);
  
    double muk2 = sqr(spectator.mass())/scale;
    return sqrt(scale * ( u*(1.-u)*(1.-x)/x - u*u*muk2 ));
  }
  
pair<double,double> IFMassiveTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double s = sqrt(1.-sqr(pt/hardPt));
  double xe = emitterX();
  return make_pair(0.5*(1.+xe-(1.-xe)*s),0.5*(1.+xe+(1.-xe)*s));
}

double IFMassiveTildeKinematics::lastZ() const {
Energy2 scale = 2.*(realEmissionMomentum()*realEmitterMomentum()
-realEmissionMomentum()*realSpectatorMomentum()
+realEmitterMomentum()*realSpectatorMomentum());
  double x = subtractionParameters()[0];
  double u = subtractionParameters()[1];
    double muk2 = sqr(bornSpectatorData()->hardProcessMass())/scale;
    return u + x - u*x*(1.-muk2);
  }


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IFMassiveTildeKinematics::persistentOutput(PersistentOStream &) const {
}

void IFMassiveTildeKinematics::persistentInput(PersistentIStream &, int) {
}

void IFMassiveTildeKinematics::Init() {

  static ClassDocumentation<IFMassiveTildeKinematics> documentation
    ("IFMassiveTildeKinematics implements the 'tilde' kinematics for "
     "a initial-final subtraction dipole involving a massive particle.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IFMassiveTildeKinematics,TildeKinematics>
describeHerwigIFMassiveTildeKinematics("Herwig::IFMassiveTildeKinematics", "Herwig.so");
#line 1 "./IILightInvertedTildeKinematics.cc"
// -*- C++ -*-
//
// IILightInvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IILightInvertedTildeKinematics class.
//

#include "IILightInvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr IILightInvertedTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IILightInvertedTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

bool IILightInvertedTildeKinematics::doMap(const double * r) {

  if ( ptMax() < ptCut() ) {
    jacobian(0.0);
    return false;
  }

  Lorentz5Momentum emitter = bornEmitterMomentum();
  Lorentz5Momentum spectator = bornSpectatorMomentum();

  double mapping = 1.0;
  pair<Energy,double> ptz = generatePtZ(mapping,r,2.0);
  if ( mapping == 0.0 ) {
    jacobian(0.0);
    return false;
  }

  Energy pt = ptz.first;
  double z = ptz.second;

  double ratio = sqr(pt/lastScale());
  double x = z*(1.-z) / ( 1. - z + ratio );
  double v = ratio * z / ( 1. - z + ratio );

  if ( x < emitterX() || x > 1. ||
       v < 0. || v > 1.-x ) {
    jacobian(0.0);
    return false;
  }

  mapping /= z*(1.-z);
  jacobian(mapping*(sqr(lastScale())/sHat())/(16.*sqr(Constants::pi)));

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = v;

  double phi = 2.*Constants::pi*r[2];
  Lorentz5Momentum kt = getKt(emitter,spectator,pt,phi);

  realEmitterMomentum() = (1./x)*emitter;
  realEmissionMomentum() = ((1.-x-v)/x)*emitter+v*spectator+kt;
  realSpectatorMomentum() = spectator;

  realEmitterMomentum().setMass(ZERO);
  realEmitterMomentum().rescaleEnergy();
  realEmissionMomentum().setMass(ZERO);
  realEmissionMomentum().rescaleEnergy();
  realSpectatorMomentum().setMass(ZERO);
  realSpectatorMomentum().rescaleEnergy();

  K = realEmitterMomentum() + realSpectatorMomentum() - realEmissionMomentum();
  K2 = K.m2();

  Ktilde = emitter + spectator;
  KplusKtilde = K + Ktilde;

  KplusKtilde2 = KplusKtilde.m2();

  return true;

}

Energy IILightInvertedTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return scale * sqrt(v*(1.-x-v)/x);

}

double IILightInvertedTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return x + v;
}

Energy IILightInvertedTildeKinematics::ptMax() const {
  return 0.5*(1.-emitterX())/sqrt(emitterX())*lastScale();
}

pair<double,double> IILightInvertedTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  hardPt = hardPt == ZERO ? ptMax() : min(hardPt,ptMax());
  if(pt>hardPt) return make_pair(0.5,0.5);
  double root = (1.-emitterX())*sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*( 1.+emitterX() - root),0.5*( 1.+emitterX() + root));
}

void IILightInvertedTildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << ounit(K,GeV) << ounit(K2,GeV2) << ounit(Ktilde,GeV)
     << ounit(KplusKtilde,GeV) << ounit(KplusKtilde2,GeV2);
}

void IILightInvertedTildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> iunit(K,GeV) >> iunit(K2,GeV2) >> iunit(Ktilde,GeV)
     >> iunit(KplusKtilde,GeV) >> iunit(KplusKtilde2,GeV2);
}

void IILightInvertedTildeKinematics::Init() {

  static ClassDocumentation<IILightInvertedTildeKinematics> documentation
    ("IILightInvertedTildeKinematics inverts the initial-initial tilde "
     "kinematics.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IILightInvertedTildeKinematics,InvertedTildeKinematics>
describeHerwigIILightInvertedTildeKinematics("Herwig::IILightInvertedTildeKinematics", "Herwig.so");
#line 1 "./IILightTildeKinematics.cc"
// -*- C++ -*-
//
// IILightTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the IILightTildeKinematics class.
//

#include "IILightTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr IILightTildeKinematics::clone() const {
  return new_ptr(*this);
}

IBPtr IILightTildeKinematics::fullclone() const {
  return new_ptr(*this);
}

Lorentz5Momentum IILightTildeKinematics::transform(const Lorentz5Momentum& k) const {

  LorentzMomentum res =
    k - 2.*((k*(K+Ktilde)/(K+Ktilde).m2())*(K+Ktilde)-((k*K)/(K.m2()))*Ktilde);

  return res;

}

bool IILightTildeKinematics::doMap() {

  Lorentz5Momentum emitter = realEmitterMomentum();
  Lorentz5Momentum emission = realEmissionMomentum();
  Lorentz5Momentum spectator = realSpectatorMomentum();

  double x = (emitter*spectator - emitter*emission - spectator*emission)/(emitter*spectator);
  double v = (emitter*emission)/(emitter*spectator);

  subtractionParameters().resize(2);
  subtractionParameters()[0] = x;
  subtractionParameters()[1] = v;

  bornEmitterMomentum() = x * emitter;
  bornSpectatorMomentum() = spectator;

  bornEmitterMomentum().setMass(ZERO);
  bornEmitterMomentum().rescaleEnergy();
  bornSpectatorMomentum().setMass(ZERO);
  bornSpectatorMomentum().rescaleEnergy();

  K = emitter + spectator - emission;
  Ktilde = x * emitter + spectator;

  return true;

}

Energy IILightTildeKinematics::lastPt() const {

  Energy scale = sqrt(2.*(bornEmitterMomentum()*bornSpectatorMomentum()));
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return scale * sqrt(v*(1.-x-v)/x);

}


Energy IILightTildeKinematics::lastPt(Lorentz5Momentum ,Lorentz5Momentum emission,Lorentz5Momentum )const {
  return emission.perp();
}

pair<double,double> IILightTildeKinematics::zBounds(Energy pt, Energy hardPt) const {
  if(pt>hardPt) return make_pair(0.5,0.5);
  double root = (1.-emitterX())*sqrt(1.-sqr(pt/hardPt));
  return make_pair(0.5*( 1.+emitterX() - root),0.5*( 1.+emitterX() + root));
}


double IILightTildeKinematics::lastZ() const {
  double x = subtractionParameters()[0];
  double v = subtractionParameters()[1];
  return x + v;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void IILightTildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << ounit(K,GeV) << ounit(Ktilde,GeV);
}

void IILightTildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> iunit(K,GeV) >> iunit(Ktilde,GeV);
}

void IILightTildeKinematics::Init() {

  static ClassDocumentation<IILightTildeKinematics> documentation
    ("IILightTildeKinematics implements the 'tilde' kinematics for "
     "a initial-initial subtraction dipole.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<IILightTildeKinematics,TildeKinematics>
describeHerwigIILightTildeKinematics("Herwig::IILightTildeKinematics", "Herwig.so");
#line 1 "./InvertedTildeKinematics.cc"
// -*- C++ -*-
//
// InvertedTildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the InvertedTildeKinematics class.
//

#include <limits>

#include "InvertedTildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;

InvertedTildeKinematics::InvertedTildeKinematics() 
  : HandlerBase(), theJacobian(0.0), thePtCut(0.0*GeV) {}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

Lorentz5Momentum InvertedTildeKinematics::getKt(const Lorentz5Momentum& p1,
						const Lorentz5Momentum& p2,
						Energy pt,
						double phi,
						bool spacelike) const {

  Lorentz5Momentum P;
  if ( !spacelike )
    P = p1 + p2;
  else
    P = p1 - p2;

  Energy2 Q2 = abs(P.m2());

  Lorentz5Momentum Q = 
    !spacelike ? 
    Lorentz5Momentum(ZERO,ZERO,ZERO,sqrt(Q2),sqrt(Q2)) :
    Lorentz5Momentum(ZERO,ZERO,sqrt(Q2),ZERO,-sqrt(Q2));

  if ( spacelike && Q.z() < P.z() )
    Q.setZ(-Q.z());

  bool boost =
    abs((P-Q).vect().mag2()/GeV2) > 1e-10 ||
    abs((P-Q).t()/GeV) > 1e-5;
  boost &= (P*Q-Q.mass2())/GeV2 > 1e-8;

  Lorentz5Momentum inFrame1;
  if ( boost )
    inFrame1 = p1 + ((P*p1-Q*p1)/(P*Q-Q.mass2()))*(P-Q);
  else
    inFrame1 = p1;

  Energy ptx = inFrame1.x();
  Energy pty = inFrame1.y();
  Energy q = 2.*inFrame1.z();

  Energy Qp = sqrt(4.*(sqr(ptx)+sqr(pty))+sqr(q));
  Energy Qy = sqrt(4.*sqr(pty)+sqr(q));

  double cPhi = cos(phi);
  double sPhi = sqrt(1.-sqr(cPhi));
  if ( phi > Constants::pi )
    sPhi = -sPhi;

  Lorentz5Momentum kt;

  if ( !spacelike ) {
    kt.setT(ZERO);
    kt.setX(pt*Qy*cPhi/Qp);
    kt.setY(-pt*(4*ptx*pty*cPhi/Qp+q*sPhi)/Qy);
    kt.setZ(2.*pt*(-ptx*q*cPhi/Qp + pty*sPhi)/Qy);
  } else {
    kt.setT(2.*pt*(ptx*q*cPhi+pty*Qp*sPhi)/(q*Qy));
    kt.setX(pt*(Qp*q*cPhi+4.*ptx*pty*sPhi)/(q*Qy));
    kt.setY(pt*Qy*sPhi/q);
    kt.setZ(ZERO);
  }

  if ( boost )
    kt = kt + ((P*kt-Q*kt)/(P*Q-Q.mass2()))*(P-Q);
  kt.setMass(-pt);
  kt.rescaleRho();

  return kt;

}

Energy InvertedTildeKinematics::lastScale() const {
  if ( ( theDipole->bornEmitter() < 2 && theDipole->bornSpectator() > 1 ) ||
       ( theDipole->bornEmitter() > 1 && theDipole->bornSpectator() < 2 ) ) {
    return -(bornEmitterMomentum()-bornSpectatorMomentum()).m();
  }
  return (bornEmitterMomentum()+bornSpectatorMomentum()).m();
}

pair<Energy,double> InvertedTildeKinematics::generatePtZ(double& jac, const double * r,
							 double pow, vector<double>* ) const {

  double kappaMin = 
    ptCut() != ZERO ?
    sqr(ptCut()/ptMax()) :
    sqr(0.1*GeV/GeV);

  double kappa;

  using namespace RandomHelpers;

  if ( ptCut() > ZERO ) {
    pair<double,double> kw = pow==1. ?
      generate(inverse(0.,kappaMin,1.),r[0]) :
      generate(power(0.,-pow,kappaMin,1.),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  } else {
    pair<double,double> kw =
      generate((piecewise(),
		flat(1e-4,kappaMin),
		match(inverse(0.,kappaMin,1.))),r[0]);
    kappa = kw.first;
    jac *= kw.second;
  }

  Energy pt = sqrt(kappa)*ptMax();

  pair<double,double> zLims = zBounds(pt);

  pair<double,double> zw(0,0);// =
  //  generate(inverse(0.,zLims.first,zLims.second)+
	//     inverse(1.,zLims.first,zLims.second),r[1]);

  // FlatZ = 1
  if ( theDipole->samplingZ() == 1 ) {
    zw = generate(flat(zLims.first,zLims.second),r[1]);
  }
  // OneOverZ = 2
  if ( theDipole->samplingZ() == 2 ) {
    zw = generate(inverse(0.0,zLims.first,zLims.second),r[1]);
  }
  // OneOverOneMinusZ = 3
  if ( theDipole->samplingZ() == 3 ) {
    zw = generate(inverse(1.0,zLims.first,zLims.second),r[1]);
  }
  // OneOverZOneMinusZ = 4
  if ( theDipole->samplingZ() == 4 ) {
    zw = generate(inverse(0.0,zLims.first,zLims.second) +
                                inverse(1.0,zLims.first,zLims.second),r[1]);
  }

  double z = zw.first;
  jac *= zw.second;

  jac *= sqr(ptMax()/lastScale());

  return make_pair(pt,z);

}

void InvertedTildeKinematics::rebind(const TranslationMap & trans) {
  theDipole = trans.translate(theDipole);
  HandlerBase::rebind(trans);
}

IVector InvertedTildeKinematics::getReferences() {
  IVector ret = HandlerBase::getReferences();
  ret.push_back(theDipole);
  return ret;
}

void InvertedTildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << theDipole << theRealXComb << theBornXComb
     << ounit(theRealEmitterMomentum,GeV) << ounit(theRealEmissionMomentum,GeV)
     << ounit(theRealSpectatorMomentum,GeV) << theJacobian
     << ounit(thePtCut,GeV);
}

void InvertedTildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theDipole >> theRealXComb >> theBornXComb
     >> iunit(theRealEmitterMomentum,GeV) >> iunit(theRealEmissionMomentum,GeV)
     >> iunit(theRealSpectatorMomentum,GeV) >> theJacobian
     >> iunit(thePtCut,GeV);
}

void InvertedTildeKinematics::Init() {

  static ClassDocumentation<InvertedTildeKinematics> documentation
    ("InvertedTildeKinematics is the base class for the inverted 'tilde' "
     "kinematics being used for subtraction terms in the "
     "formalism of Catani and Seymour.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<InvertedTildeKinematics,HandlerBase>
describeInvertedTildeKinematics("Herwig::InvertedTildeKinematics", "Herwig.so");
#line 1 "./MatchboxPhasespace.cc"
// -*- C++ -*-
//
// MatchboxPhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxPhasespace class.
//

#include "MatchboxPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/ProcessData.h"
#include "Herwig/MatrixElement/Matchbox/MatchboxFactory.h"

using namespace Herwig;

MatchboxPhasespace::MatchboxPhasespace() 
  : singularCutoff(10*GeV), theUseMassGenerators(false),
    theLoopParticleIdMin(200001), theLoopParticleIdMax(200100) {}

void MatchboxPhasespace::cloneDependencies(const std::string&) {}

Ptr<MatchboxFactory>::tcptr MatchboxPhasespace::factory() const {
  return MatchboxFactory::currentFactory();
}

Ptr<ProcessData>::tptr MatchboxPhasespace::processData() const {
  return factory()->processData();
}

double MatchboxPhasespace::generateKinematics(const double* r,
					      vector<Lorentz5Momentum>& momenta) {

  diagramWeights().clear();

  cPDVector::const_iterator pd = mePartonData().begin() + 2;
  vector<Lorentz5Momentum>::iterator p = momenta.begin() + 2;

  double massJacobian = 1.;
  Energy summ = ZERO;

  if ( useMassGenerators() ) {
    Energy gmass = ZERO;
    tGenericMassGeneratorPtr mgen;
    Energy maxMass = 
      (!haveX1X2() && momenta.size() > 3) ? 
      sqrt(lastSHat()) : sqrt(lastS());
    for ( ; pd != mePartonData().end(); ++pd, ++p ) {
      mgen = processData()->massGenerator(*pd);
      if ( mgen && !isInvertible() ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	gmass = mgen->mass(massJacobian,**pd,massMin,massMax,r[0]);
	++r;
      } else if ( (**pd).hardProcessWidth() != ZERO ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	// use a standard Breit Wigner here which we can invert
	// see invertKinematics as well
	double bwILow = 
	  atan((sqr(massMin)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	double bwIUp = 
	  atan((sqr(massMax)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	gmass = sqrt(sqr((**pd).hardProcessMass()) + 
		     (**pd).hardProcessMass()*(**pd).hardProcessWidth()*tan(bwILow+r[0]*(bwIUp-bwILow)));
	++r;
      } else {
	gmass = (**pd).hardProcessMass();
      }
      maxMass -= gmass;
      p->setMass(gmass);
      summ += gmass;
    }
  } else {
    for ( ; pd != mePartonData().end(); ++pd, ++p ) {
      summ += (**pd).hardProcessMass();
      p->setMass((**pd).hardProcessMass());
    }
  }

  if ( momenta.size() > 3 && !haveX1X2() ) {
    if ( summ > (momenta[0]+momenta[1]).m() )
      return 0.0;
  }

  double weight = momenta.size() > 3 ? 
    generateTwoToNKinematics(r,momenta) : 
    generateTwoToOneKinematics(r,momenta);

  fillDiagramWeights();

  return weight*massJacobian;

}

double MatchboxPhasespace::generateTwoToOneKinematics(const double* r,
						      vector<Lorentz5Momentum>& momenta) {

  double tau = momenta[2].mass2()/lastXCombPtr()->lastS();
  double ltau = log(tau)/2.;
  //old:  y = ltau - 2.*r[0]*ltau; x1 = sqrt(tau)*exp(y); x2 = sqrt(tau)*exp(-y);
  double x1=pow(tau,1.-r[0]);
  double x2=pow(tau,r[0]);

  // Due to the proton mass and P1.e() + P2.e() == lastS() we multiply here 
  // with the correction factor abs(P1.e()/P1.z()) to produce incoming 
  // p1/2 = (e1/2,0,0,+/- e1/2)  
  Lorentz5Momentum P1 = lastXCombPtr()->lastParticles().first->momentum();
  ThreeVector<Energy> p1 =  x1 * (P1.vect()) * abs(P1.e()/P1.z());

  Lorentz5Momentum P2 = lastXCombPtr()->lastParticles().second->momentum();
  ThreeVector<Energy> p2 =  x2 * (P2.vect()) * abs(P2.e()/P2.z());

  ThreeVector<Energy> q = p1 + p2;

  momenta[0] = Lorentz5Momentum(momenta[0].mass(),p1);
  momenta[1] = Lorentz5Momentum(momenta[1].mass(),p2);
  momenta[2] = Lorentz5Momentum(momenta[2].mass(),q);
  
  // check for energy conservation:
  if ((momenta[0]+momenta[1]-momenta[2]).e()>pow(10,-9)*GeV)
    generator()->log() 
    << "Warning: Momentum conservation in generateTwoToOneKinematics not precise.\n"
    << flush;
  

  lastXCombPtr()->lastX1X2({x1,x2});
  lastXCombPtr()->lastSHat((momenta[0]+momenta[1]).m2());

  return -4.*Constants::pi*ltau;

}

double MatchboxPhasespace::invertKinematics(const vector<Lorentz5Momentum>& momenta,
					    double* r) const {

  if ( useMassGenerators() ) {

    Energy gmass = ZERO;
    Energy maxMass = 
      (!haveX1X2() && momenta.size() > 3) ? 
      sqrt((momenta[0]+momenta[1]).m2()) : sqrt(lastS());

    cPDVector::const_iterator pd = mePartonData().begin() + 2;
    vector<Lorentz5Momentum>::const_iterator p = momenta.begin() + 2;

    for ( ; pd != mePartonData().end(); ++pd, ++p ) {

      if ( (**pd).hardProcessWidth() != ZERO ) {
	Energy massMax = min((**pd).massMax(),maxMass);
	Energy massMin = (**pd).massMin();
	if ( massMin > massMax )
	  return 0.0;
	double bwILow = 
	  atan((sqr(massMin)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	double bwIUp = 
	  atan((sqr(massMax)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	gmass = p->mass();
	double bw =
	  atan((sqr(gmass)-sqr((**pd).hardProcessMass()))/((**pd).hardProcessMass() * (**pd).hardProcessWidth()));
	r[0] = (bw-bwILow)/(bwIUp-bwILow);
	++r;
      } else {
	gmass = (**pd).hardProcessMass();
      }
      maxMass -= gmass;

    }

  }

  return momenta.size() > 3 ? 
    invertTwoToNKinematics(momenta,r) : 
    invertTwoToOneKinematics(momenta,r);
}

double MatchboxPhasespace::invertTwoToOneKinematics(const vector<Lorentz5Momentum>& momenta,
						    double* r) const {

  double tau = momenta[2].mass2()/lastXCombPtr()->lastS();
  double ltau = log(tau)/2.;

  r[0] = (ltau - (momenta[0]+momenta[1]).rapidity())/(2.*ltau);

  return -4.*Constants::pi*ltau;

}

void MatchboxPhasespace::setCoupling(long a, long b, long c, 
				     double coupling, bool includeCrossings) {
  cPDPtr A = getParticleData(a);
  cPDPtr B = getParticleData(b);
  cPDPtr C = getParticleData(c);
  if ( !A || !B || !C ) {
    generator()->log() << "Warning: could not determine particle data for ids "
		       << a << " " << b << " " << c << " when setting coupling in MatchboxPhasespace.\n"
		       << flush;
    return;
  }
  if ( !includeCrossings ) {
    theCouplings->couplings()[LTriple(a,b,c)] = coupling;
    return;
  }
  if ( A->CC() ) {
    theCouplings->couplings()[LTriple(-a,b,c)] = coupling;
    theCouplings->couplings()[LTriple(-a,c,b)] = coupling;
  } else {
    theCouplings->couplings()[LTriple(a,b,c)] = coupling;
    theCouplings->couplings()[LTriple(a,c,b)] = coupling;
  }
  if ( B->CC() ) {
    theCouplings->couplings()[LTriple(-b,a,c)] = coupling;
    theCouplings->couplings()[LTriple(-b,c,a)] = coupling;
  } else {
    theCouplings->couplings()[LTriple(b,a,c)] = coupling;
    theCouplings->couplings()[LTriple(b,c,a)] = coupling;
  }
  if ( C->CC() ) {
    theCouplings->couplings()[LTriple(-c,a,b)] = coupling;
    theCouplings->couplings()[LTriple(-c,b,a)] = coupling;
  } else {
    theCouplings->couplings()[LTriple(c,a,b)] = coupling;
    theCouplings->couplings()[LTriple(c,b,a)] = coupling;
  }
}

string MatchboxPhasespace::doSetCoupling(string in) {
  istringstream is(in);
  long a,b,c; double coupling;
  is >> a >> b >> c >> coupling;
  if ( !is )
    return "MatchboxPhasespace: error in setting coupling.";
  setCoupling(a,b,c,coupling,true);
  return "";
}

string MatchboxPhasespace::doSetPhysicalCoupling(string in) {
  istringstream is(in);
  long a,b,c; double coupling;
  is >> a >> b >> c >> coupling;
  if ( !is )
    return "MatchboxPhasespace: error in setting coupling.";
  setCoupling(a,b,c,coupling,false);
  return "";
}


pair<double,Lorentz5Momentum> 
MatchboxPhasespace::timeLikeWeight(const Tree2toNDiagram& diag,
				   int branch, double flatCut) const {

  pair<int,int> children = diag.children(branch);

  if ( children.first == -1 ) {
    return make_pair(1.,meMomenta()[diag.externalId(branch)]);
  }

  pair<double,Lorentz5Momentum> res
    = timeLikeWeight(diag,children.first,flatCut);

  pair<double,Lorentz5Momentum> other
    = timeLikeWeight(diag,children.second,flatCut);

  res.first *= other.first;
  res.second += other.second;

  LTriple vertexKey(diag.allPartons()[branch]->id(),
		    diag.allPartons()[children.first]->id(),
		    diag.allPartons()[children.second]->id());
  map<LTriple,double>::const_iterator cit = theCouplings->couplings().find(vertexKey);
  if ( cit != theCouplings->couplings().end() ){
    res.first *= cit->second;
  }

  Energy2 mass2 = sqr(diag.allPartons()[branch]->hardProcessMass());
  Energy2 width2 = sqr(diag.allPartons()[branch]->hardProcessWidth());

  if ( abs(diag.allPartons()[branch]->id()) >= theLoopParticleIdMin
       && abs(diag.allPartons()[branch]->id()) <= theLoopParticleIdMax ) { // "loop particle"

    if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut ) {
      res.first /=
	abs((res.second.m2()-mass2)/GeV2);
      res.first *=
	log(abs((res.second.m2()-mass2)/GeV2)); // normal. of the argument in the log?
    }

  } else {

    if ( width2 == ZERO ) {
      if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut )
	res.first /=
	  abs((res.second.m2()-mass2)/GeV2);
    } else {
      res.first /=
	(sqr((res.second.m2()-mass2)/GeV2) +
	 mass2*width2/sqr(GeV2))/(abs(res.second.m2()/GeV2));
    }

  }

  return res;

}

double MatchboxPhasespace::spaceLikeWeight(const Tree2toNDiagram& diag,
					   const Lorentz5Momentum& incoming,
					   int branch, double flatCut) const {

  if ( branch == -1 )
    return 1.;

  pair<int,int> children = diag.children(branch);

  pair<double,Lorentz5Momentum> res =
    timeLikeWeight(diag,children.second,flatCut);

  
  LTriple vertexKey(diag.allPartons()[branch]->id(),
		    diag.allPartons()[children.first]->id(),
		    diag.allPartons()[children.second]->id());
  if ( children.first == diag.nSpace() - 1 ) {
    if ( diag.allPartons()[children.first]->CC() )
      vertexKey = LTriple(diag.allPartons()[branch]->id(),
			  diag.allPartons()[children.second]->id(),
			  diag.allPartons()[children.first]->CC()->id());
    else
      vertexKey = LTriple(diag.allPartons()[branch]->id(),
			  diag.allPartons()[children.second]->id(),
			  diag.allPartons()[children.first]->id());
  }
  map<LTriple,double>::const_iterator cit = theCouplings->couplings().find(vertexKey);
  if ( cit != theCouplings->couplings().end() ){
    res.first *= cit->second;
  }
  if ( children.first == diag.nSpace() - 1 ) {
    return res.first;
  }

  res.second = incoming - res.second;

  Energy2 mass2 = sqr(diag.allPartons()[children.first]->hardProcessMass());
  Energy2 width2 = sqr(diag.allPartons()[children.first]->hardProcessWidth());

  if ( abs(diag.allPartons()[children.first]->id()) >= theLoopParticleIdMin
       && (diag.allPartons()[children.first]->id()) <= theLoopParticleIdMax ) { // "loop particle"

    if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut ) {
      res.first /=
	abs((res.second.m2()-mass2)/GeV2);
      res.first *=
	log(abs((res.second.m2()-mass2)/GeV2)); // normal. of the argument in the log?
    }

  } else {

    if ( width2 == ZERO ) {
      if ( abs((res.second.m2()-mass2)/lastSHat()) > flatCut )
	res.first /=
	  abs((res.second.m2()-mass2)/GeV2);
    } else {
      res.first /=
	(sqr((res.second.m2()-mass2)/GeV2) +
	 mass2*width2/sqr(GeV2))/(abs(res.second.m2()/GeV2));
    }

  }

  return
    res.first * spaceLikeWeight(diag,res.second,children.first,flatCut);

}

void MatchboxPhasespace::fillDiagramWeights(double flatCut) {

  if ( !diagramWeights().empty() )
    return;

  for ( auto & d : lastXComb().diagrams() ) {
    diagramWeights()[d->id()] =
      spaceLikeWeight(dynamic_cast<const Tree2toNDiagram&>(*d),meMomenta()[0],0,flatCut);
  }

}

Selector<MEBase::DiagramIndex> 
MatchboxPhasespace::selectDiagrams(const MEBase::DiagramVector& diags) const {
  Selector<MEBase::DiagramIndex> ret;
  for ( MEBase::DiagramIndex d = 0; d < diags.size(); ++d ) {
    ret.insert(diagramWeight(dynamic_cast<const Tree2toNDiagram&>(*diags[d])),d);
  }
  return ret;
}

bool MatchboxPhasespace::matchConstraints(const vector<Lorentz5Momentum>& momenta) {

  if ( singularLimits().empty() )
    return true;

  lastSingularLimit() = singularLimits().begin();

  for ( ; lastSingularLimit() != singularLimits().end(); ++lastSingularLimit() ) {
    if ( lastSingularLimit()->first == lastSingularLimit()->second &&
	 momenta[lastSingularLimit()->first].t() < singularCutoff )
      break;
    if ( lastSingularLimit()->first != lastSingularLimit()->second &&
	 sqrt(momenta[lastSingularLimit()->first]*
	      momenta[lastSingularLimit()->second]) < singularCutoff ) {
      bool match = true;
      for ( set<pair<size_t,size_t> >::const_iterator other =
	      singularLimits().begin(); other != singularLimits().end(); ++other ) {
	if ( other == lastSingularLimit() )
	  continue;
	if ( other->first == other->second &&
	     momenta[other->first].t() < singularCutoff ) {
	  match = false;
	  break;
	}
	if ( other->first != other->second &&
	     sqrt(momenta[other->first]*
		  momenta[other->second]) < singularCutoff ) {
	  match = false;
	  break;
	}
      }
      if ( match )
	break;
    }
  }

  return lastSingularLimit() != singularLimits().end();

}

int MatchboxPhasespace::nDim(const cPDVector& data) const {

  int ndimps = nDimPhasespace(data.size()-2);

  if ( useMassGenerators() ) {
    for ( cPDVector::const_iterator pd = data.begin();
	  pd != data.end(); ++pd ) {
      if ( (**pd).massGenerator() ||
	   (**pd).hardProcessWidth() != ZERO ) {
	++ndimps;
      }
    }
  }

  return ndimps;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void MatchboxPhasespace::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb
     << ounit(singularCutoff,GeV) << theUseMassGenerators 
     << theLoopParticleIdMin << theLoopParticleIdMax
     << theCouplings;
}

void MatchboxPhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb
     >> iunit(singularCutoff,GeV) >> theUseMassGenerators
     >> theLoopParticleIdMin >> theLoopParticleIdMax
     >> theCouplings;
  lastMatchboxXComb(theLastXComb);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxPhasespace,HandlerBase>
  describeMatchboxPhasespace("Herwig::MatchboxPhasespace", "Herwig.so");

void MatchboxPhasespace::Init() {

  static ClassDocumentation<MatchboxPhasespace> documentation
    ("MatchboxPhasespace defines an abstract interface to a phase "
     "space generator.");


  static Parameter<MatchboxPhasespace,Energy> interfaceSingularCutoff
    ("SingularCutoff",
     "[debug] Cutoff below which a region is considered singular.",
     &MatchboxPhasespace::singularCutoff, GeV, 10.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  interfaceSingularCutoff.rank(-1);

  /*
  static Switch<MatchboxPhasespace,bool> interfaceUseMassGenerators
    ("UseMassGenerators",
     "Use mass generators instead of fixed masses.",
     &MatchboxPhasespace::theUseMassGenerators, false, false, false);
  static SwitchOption interfaceUseMassGeneratorsYes
    (interfaceUseMassGenerators,
     "Yes",
     "Use mass generators.",
     true);
  static SwitchOption interfaceUseMassGeneratorsNo
    (interfaceUseMassGenerators,
     "No",
     "Do not use mass generators.",
     false);
  */

  static Command<MatchboxPhasespace> interfaceSetCoupling
    ("SetCoupling",
     "",
     &MatchboxPhasespace::doSetCoupling, false);

  static Command<MatchboxPhasespace> interfaceSetPhysicalCoupling
    ("SetPhysicalCoupling",
     "",
     &MatchboxPhasespace::doSetPhysicalCoupling, false);

  static Parameter<MatchboxPhasespace,int> interfaceLoopParticleIdMin
    ("LoopParticleIdMin",
     "First id in a range of id's meant to denote fictitious "
     "'ghost' particles to be used by the diagram generator "
     "in loop induced processes.",
     &MatchboxPhasespace::theLoopParticleIdMin, 200001, 0, 0,
     false, false, Interface::lowerlim);
  interfaceLoopParticleIdMin.rank(-1);

  static Parameter<MatchboxPhasespace,int> interfaceLoopParticleIdMax
    ("LoopParticleIdMax",
     "Last id in a range of id's meant to denote fictitious "
     "'ghost' particles to be used by the diagram generator "
     "in loop induced processes.",
     &MatchboxPhasespace::theLoopParticleIdMax, 200100, 0, 0,
     false, false, Interface::lowerlim);
  interfaceLoopParticleIdMax.rank(-1);


  static Reference<MatchboxPhasespace,PhasespaceCouplings> interfaceCouplingData
    ("CouplingData",
     "Set the storage for the couplings.",
     &MatchboxPhasespace::theCouplings, false, false, true, false, false);
  interfaceCouplingData.rank(-1);

}

#line 1 "./MatchboxRambo.cc"
// -*- C++ -*-
//
// MatchboxRambo.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxRambo class.
//

#include "MatchboxRambo.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxRambo::MatchboxRambo() 
  : needToReshuffle(false), theMakeReferenceSample(false),
    referenceSample(0) {}

IBPtr MatchboxRambo::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxRambo::fullclone() const {
  return new_ptr(*this);
}

static double weights[7] = {

  -1.,-1.,
  0.039788735772973833942,
  0.00012598255637968550463,
  1.3296564302788840628E-7,
  7.0167897579949011130E-11,
  2.2217170114046130768E-14

};

void MatchboxRambo::setXComb(tStdXCombPtr xc) {
  MatchboxPhasespace::setXComb(xc);
  needToReshuffle = false;
  if ( xc ) {
    for ( cPDVector::const_iterator d = mePartonData().begin();
	  d != mePartonData().end(); ++d ) {
      if ( (**d).hardProcessMass() != ZERO ) {
	needToReshuffle = true;
	break;
      }
    }
  }
}

void MatchboxRambo::dumpReference(const vector<Lorentz5Momentum>& momenta, double weight) const {
  *referenceSample << lastX1() << " " << lastX2() << " ";
  Boost toLab = (lastPartons().first->momentum() + 
		 lastPartons().second->momentum()).boostVector();
  for ( vector<Lorentz5Momentum>::const_iterator p = momenta.begin();
	p != momenta.end(); ++p ) {
    Lorentz5Momentum pl = *p;
    if ( toLab.mag2() > Constants::epsilon )
      pl.boost(toLab);
    *referenceSample 
      << (pl.x()/GeV) << " "
      << (pl.y()/GeV) << " "
      << (pl.z()/GeV) << " "
      << (pl.t()/GeV) << " "
      << (pl.mass()/GeV) << " ";
  }
  double ymax = lastCuts().yHatMax();
  double ymin = lastCuts().yHatMin();
  double km = log(lastCuts().sHatMax()/lastCuts().sHatMin());
  ymax = min(ymax, log(lastCuts().x1Max()*sqrt(lastS()/lastSHat())));
  ymin = max(ymin, -log(lastCuts().x2Max()*sqrt(lastS()/lastSHat())));
  *referenceSample << weight*km*(ymax-ymin)/(lastX1()*lastX2()) << "\n" << flush;
}

double MatchboxRambo::generateTwoToNKinematics(const double* r,
					       vector<Lorentz5Momentum>& momenta) {

  if ( theMakeReferenceSample ) {
    map<cPDVector,ofstream*>::iterator ref =
      referenceSamples.find(mePartonData());
    if ( ref == referenceSamples.end() ) {
      ostringstream refname;
      for ( cPDVector::const_iterator p = mePartonData().begin();
	    p != mePartonData().end(); ++p ) {
	refname << (**p).PDGName();
      }
      refname << ".rambo";
      referenceSamples[mePartonData()] = new ofstream(refname.str().c_str(),std::ios_base::app);
      ref = referenceSamples.find(mePartonData());
      *(ref->second) << setprecision(26);
    }
    assert(ref != referenceSamples.end());
    referenceSample = ref->second;
  }

  size_t offset = 2;
  if ( lastXCombPtr() )
    offset = dynamic_cast<const Tree2toNDiagram&>(*lastXComb().diagrams().front()).nSpace() > 0 ? 2 : 1;

  Energy w = sqrt(lastSHat());
  size_t count = 0;
  Lorentz5Momentum Q;
  for ( vector<Lorentz5Momentum>::iterator k = momenta.begin() + offset;
	k != momenta.end(); ++k ) {
    Energy q = -w*log(r[count]*r[count+1]);
    double ct = 2.*r[count+2]-1.;
    double st = sqrt(1.-sqr(ct));
    double phi = 2.*Constants::pi*r[count+3];
    double cphi = cos(phi);
    double sphi = sqrt(1.-sqr(cphi));
    if ( phi > Constants::pi )
      sphi = -sphi;
    (*k).setMass(ZERO);
    (*k).setT(q);
    (*k).setX(q*cphi*st);
    (*k).setY(q*sphi*st);
    (*k).setZ(q*ct);
    count += 4;
    Q += *k;
  }

  Energy M = sqrt(Q.m2());
  double x = w/M;
  Boost beta = -(Q.vect() * (1./M));
  double gamma = Q.t()/M;
  double a = 1./(1.+gamma);

  for ( vector<Lorentz5Momentum>::iterator k = momenta.begin() + offset;
	k != momenta.end(); ++k ) {
    Energy q = (*k).t();
    Energy bq = beta*(*k).vect();
    (*k).setT(x*(gamma*q+bq));
    (*k).setVect(x*((*k).vect()+(q+a*bq)*beta));
  }

  size_t n = momenta.size()-offset;
  double weight = weights[n];

  if ( !needToReshuffle ) {
    if ( !matchConstraints(momenta) )
      return 0.;
    if ( theMakeReferenceSample )
      dumpReference(momenta, weight);
    return weight;
  }

  double xi;

  ReshuffleEquation solve(w,mePartonData().begin()+offset,mePartonData().end(),
			  momenta.begin()+2,momenta.end());

  GSLBisection solver(1e-10,1e-8,10000);

  try {
    xi = solver.value(solve,0.0,1.1);
  } catch (GSLBisection::GSLerror) {
    return 0.;
  } catch (GSLBisection::IntervalError) {
    return 0.;
  }

  weight *= pow(xi,3.*(n-1.));

  Energy num = ZERO;
  Energy den = ZERO;

  cPDVector::const_iterator d = mePartonData().begin()+offset;
  for ( vector<Lorentz5Momentum>::iterator k = momenta.begin()+offset;
	k != momenta.end(); ++k, ++d ) {
    num += (*k).vect().mag2()/(*k).t();
    Energy q = (*k).t();
    (*k).setT(sqrt(sqr((**d).hardProcessMass())+xi*xi*sqr((*k).t())));
    (*k).setVect(xi*(*k).vect());
    weight *= q/(*k).t();
    den += (*k).vect().mag2()/(*k).t();
    (*k).setMass((**d).hardProcessMass());
  }

  if ( !matchConstraints(momenta) )
    return 0.;

  weight *= num/den;

  if ( theMakeReferenceSample )
    dumpReference(momenta, weight);
  
  return weight;

}

Energy MatchboxRambo::ReshuffleEquation::operator() (double xi) const {
  cPDVector::const_iterator d = dataBegin;
  vector<Lorentz5Momentum>::const_iterator p = momentaBegin;
  Energy res = -w;
  for ( ; d != dataEnd; ++d, ++p ) {
    res += sqrt(sqr((**d).hardProcessMass()) +
		xi*xi*sqr(p->t()));
  }
  return res;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxRambo::persistentOutput(PersistentOStream & os) const {
  os << needToReshuffle << theMakeReferenceSample;
}

void MatchboxRambo::persistentInput(PersistentIStream & is, int) {
  is >> needToReshuffle >> theMakeReferenceSample;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxRambo,MatchboxPhasespace>
  describeHerwigMatchboxRambo("Herwig::MatchboxRambo", "Herwig.so");

void MatchboxRambo::Init() {

  static ClassDocumentation<MatchboxRambo> documentation
    ("MatchboxRambo implements RAMBO phase space generation.");


  static Switch<MatchboxRambo,bool> interfaceMakeReferenceSample
    ("MakeReferenceSample",
     "Switch on generation of a reference sample of phase space points.",
     &MatchboxRambo::theMakeReferenceSample, false, false, false);
  static SwitchOption interfaceMakeReferenceSampleYes
    (interfaceMakeReferenceSample,
     "Yes",
     "Generate a reference sample.",
     true);
  static SwitchOption interfaceMakeReferenceSampleNo
    (interfaceMakeReferenceSample,
     "No",
     "Do not generate a reference sample.",
     false);
  interfaceMakeReferenceSample.rank(-1);

}

#line 1 "./PhasespaceHelpers.cc"
// -*- C++ -*-
//
// PhasespaceHelpers.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "PhasespaceHelpers.h"
#include "Herwig/MatrixElement/Matchbox/Phasespace/RandomHelpers.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;
using namespace Herwig::RandomHelpers;

Energy PhasespaceInfo::generateMass(tcPDPtr data, 
                                    const pair<Energy,Energy>& range) {

  double xlow = sqr(range.first)/sHat;
  if ( range.first < ZERO )
    xlow = -xlow;

  double xup = sqr(range.second)/sHat;
  if ( range.second < ZERO )
    xup = -xup;

  double mu = sqr(data->hardProcessMass())/sHat;
  double gamma = sqr(data->hardProcessWidth())/sHat;

  if ( gamma < 1e-14 )
    gamma = 0.0;

  if ( M0 != ZERO )
    x0 = M0/sqrtSHat;
  if ( Mc != ZERO )
    xc = Mc/sqrtSHat;

  double r = rnd();

  pair<double,double> event;
  
  if ( gamma == 0. ) {

    if ( mu < xlow || mu > xup ) {
      if ( abs(xlow-mu) < xc )
        xlow = mu;
      if ( abs(xup-mu) < xc )
        xup = mu;
    }

    if ( mu < xlow || mu > xup ) {
      event =
      generate(inverse(mu,xlow,xup),r);
    } else {
      pair<double,double> pLeft(xlow,xlow < mu-x0 ? mu-x0 : xlow);
      pair<double,double> pRight(xup > mu+x0 ? mu+x0 : xup,xup);
      pair<double,double> fLeft(pLeft.second < mu-x0 ? mu-x0 : pLeft.second,mu-xc);
      if ( fLeft.first >= fLeft.second ) fLeft.first = fLeft.second;
      pair<double,double> fRight(mu+xc,pRight.first > mu+x0 ? mu+x0 : pRight.first);
      if ( fRight.first >= fRight.second ) fRight.second = fRight.first;

      
      const bool pL= abs( pLeft.first - pLeft.second ) < 1e-14;
      const bool fL= abs( fLeft.first - fLeft.second ) < 1e-14;
      const bool fR= abs( fRight.first - fRight.second ) < 1e-14;
      const bool pR= abs( pRight.first - pRight.second ) < 1e-14;


      if ( !pL && !fL && !fR && !pR ) {
        event =
        generate((piecewise(),
                  inverse(mu,pLeft.first,pLeft.second),
                  match(flat(fLeft.first,fLeft.second))) +
                 match((piecewise(),
                        flat(fRight.first,fRight.second),
                        match(inverse(mu,pRight.first,pRight.second)))),
                 r);
      } else if ( pL && !fL && !fR && !pR ) {
        event =
        generate(flat(fLeft.first,fLeft.second) +
                 match((piecewise(),
                        flat(fRight.first,fRight.second),
                        match(inverse(mu,pRight.first,pRight.second)))), r);
      } else if ( !pL && !fL && !fR && pR ) {
        event =
        generate((piecewise(),
                  inverse(mu,pLeft.first,pLeft.second),
                  match(flat(fLeft.first,fLeft.second))) +
                 match(flat(fRight.first,fRight.second)),
                 r);
      } else if ( pL && fL && !fR && !pR ) {
        event =
        generate((piecewise(),flat(fRight.first,fRight.second),
                  match(inverse(mu,pRight.first,pRight.second))), r);
      } else if ( !pL && !fL && fR && pR ) {
        event =
        generate((piecewise(),
                  inverse(mu,pLeft.first,pLeft.second),
                  match(flat(fLeft.first,fLeft.second))),r);
      } else if ( pL && !fL && !fR && pR ) {
        event = generate(flat(fLeft.first,fLeft.second) +
                 match(flat(fRight.first,fRight.second)),r);
      } else if ( pL && fL && !fR && pR ) {
        event = generate(flat(fRight.first,fRight.second),r);
      } else if ( pL && !fL && fR && pR ) {
        event = generate(flat(fLeft.first,fLeft.second),r);
      } else if ( pL && fL && fR && pR ) {
        throw Veto();
      } else assert(false);
    }
  } else {
    event = generate(breitWigner(mu,gamma,xlow,xup),r);
  }

  if ( abs(event.first) < xc )
    throw Veto();

  weight *= event.second;

  Energy res = sqrt(abs(event.first)*sHat);
  if ( event.first < 0. )
    res = -res;

  return res;
}

Lorentz5Momentum PhasespaceInfo::generateKt(const Lorentz5Momentum& p1,
                                            const Lorentz5Momentum& p2,
                                            Energy pt) {

  double phi = 2.*Constants::pi*rnd();
  weight *= 2.*Constants::pi;

  Lorentz5Momentum P = p1 + p2;

  Energy2 Q2 = abs(P.m2());

  Lorentz5Momentum Q =  Lorentz5Momentum(ZERO,ZERO,ZERO,sqrt(Q2),sqrt(Q2));

  bool boost =
    abs((P-Q).vect().mag2()/GeV2) > 1e-8 ||
    abs((P-Q).t()/GeV) > 1e-4;
  boost &= (P*Q-Q.mass2())/GeV2 > 1e-8;

  Lorentz5Momentum inFrame1;
  if ( boost )
    inFrame1 = p1 + ((P*p1-Q*p1)/(P*Q-Q.mass2()))*(P-Q);
  else
    inFrame1 = p1;

  Energy ptx = inFrame1.x();
  Energy pty = inFrame1.y();
  Energy q = 2.*inFrame1.z();

  Energy Qp = sqrt(4.*(sqr(ptx)+sqr(pty))+sqr(q));
  Energy Qy = sqrt(4.*sqr(pty)+sqr(q));

  double cPhi = cos(phi);
  double sPhi = sqrt(1.-sqr(cPhi));
  if ( phi > Constants::pi )
    sPhi = -sPhi;

  Lorentz5Momentum kt;

  kt.setT(ZERO);
  kt.setX(pt*Qy*cPhi/Qp);
  kt.setY(-pt*(4*ptx*pty*cPhi/Qp+q*sPhi)/Qy);
  kt.setZ(2.*pt*(-ptx*q*cPhi/Qp + pty*sPhi)/Qy);

  if ( boost )
    kt = kt + ((P*kt-Q*kt)/(P*Q-Q.mass2()))*(P-Q);
  kt.setMass(-pt);
  kt.rescaleRho();

  return kt;

}

void PhasespaceTree::setup(const Tree2toNDiagram& diag, 
                           int pos) {
  doMirror = false;

  pair<int,int> dchildren =  diag.children(pos);

  data = diag.allPartons()[pos];

  spacelike = pos < diag.nSpace();

  if ( pos == 0 )
    externalId = 0;

  if ( dchildren.first == -1 ) {
    externalId = diag.externalId(pos);
    leafs.insert(externalId);
    return;
  }

  children.push_back(PhasespaceTree());
  children.back().setup(diag,dchildren.first);
  children.push_back(PhasespaceTree());
  children.back().setup(diag,dchildren.second);

  if ( !children[0].children.empty() &&
      children[1].children.empty() &&
      !spacelike )
    swap(children[0],children[1]);
  if ( spacelike &&
      !children[0].spacelike )
    swap(children[0],children[1]);

  copy(children[0].leafs.begin(),children[0].leafs.end(),
       inserter(leafs,leafs.begin()));
  copy(children[1].leafs.begin(),children[1].leafs.end(),
       inserter(leafs,leafs.begin()));

}

void PhasespaceTree::setupMirrored(const Tree2toNDiagram& diag, 
                                   int pos) {

  doMirror = true;

  spacelike = pos < diag.nSpace();

  pair<int,int> dchildren;
  if (pos != 0 && spacelike)
    dchildren = {diag.parent(pos), diag.children(diag.parent(pos)).second};
  else if ( !spacelike ) dchildren = diag.children(pos);
  else                   dchildren = {-1,-1};

  data = diag.allPartons()[pos];

  if ( pos == diag.nSpace() - 1 )
    externalId = 1;

  if ( dchildren.first == -1 ) {
    externalId = diag.externalId(pos);
    leafs.insert(externalId);
    return;
  }
  


  children.push_back(PhasespaceTree());
  children.back().setupMirrored(diag,dchildren.first);
  children.push_back(PhasespaceTree());
  children.back().setupMirrored(diag,dchildren.second);

  if ( !children[0].children.empty() &&
      children[1].children.empty() &&
      !spacelike )
    swap(children[0],children[1]);
  if ( spacelike &&
      !children[0].spacelike ) {
    assert (false);
  }

  copy(children[0].leafs.begin(),children[0].leafs.end(),
       inserter(leafs,leafs.begin()));
  copy(children[1].leafs.begin(),children[1].leafs.end(),
       inserter(leafs,leafs.begin()));

}

void PhasespaceTree::print(int in) {
  for (int i = 0; i != in; i++)
    cerr << "   ";
  cerr << " |- "  << data->PDGName() << " " << externalId << "\n" << flush;
  if ( !children.empty() ) {
    children[1].print(in+1);
    children[0].print(in+int(!spacelike));
  }
  else {
    cerr << "\n";
  }
}

void PhasespaceTree::init(const vector<Lorentz5Momentum>& meMomenta) {

  if ( children.empty() ) {
    massRange.first = meMomenta[externalId].mass();
    massRange.second = massRange.first;
    if ( !doMirror && externalId == 1 )
      momentum = meMomenta[1];
    if ( doMirror && externalId == 0 )
      momentum = meMomenta[0];
    momentum.setMass(meMomenta[externalId].mass());
    return;
  }

  children[0].init(meMomenta);
  children[1].init(meMomenta);

  if ( !children[0].spacelike &&
      !children[1].spacelike ) {
    massRange.first =
    children[0].massRange.first +
    children[1].massRange.first;
  }

}

void PhasespaceTree::generateKinematics(PhasespaceInfo& info,
                                        vector<Lorentz5Momentum>& meMomenta) {
  
  if ( !doMirror && externalId == 0 ) {
    init(meMomenta);
    Energy2 s = (meMomenta[0]+meMomenta[1]).m2();
    double sign = meMomenta[0].z() >= ZERO ? 1. : -1;
    momentum = Lorentz5Momentum(ZERO,ZERO,sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
    backwardMomentum = Lorentz5Momentum(ZERO,ZERO,-sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
  }
  else if ( doMirror && externalId == 1) {
    init(meMomenta);
    Energy2 s = (meMomenta[0]+meMomenta[1]).m2();
    double sign = meMomenta[0].z() >= ZERO ? 1. : -1;
    momentum = Lorentz5Momentum(ZERO,ZERO,-sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
    backwardMomentum = Lorentz5Momentum(ZERO,ZERO,sign*sqrt(s)/2.,sqrt(s)/2.,ZERO);
  }

  if ( children.empty() ) {
    if ( ( !doMirror && externalId != 1 )
        || ( doMirror && externalId !=0 ) )
      meMomenta[externalId] = momentum;
    return;
  }

    // s-channel
  if ( ( !doMirror && externalId == 0 &&
          children[0].externalId == 1 )
     || ( doMirror && externalId == 1 &&
          children[0].externalId == 0 ) ) {
        children[1].momentum = meMomenta[0] + meMomenta[1];
        children[1].momentum.setMass(info.sqrtSHat);
        children[1].momentum.rescaleEnergy();
        children[1].generateKinematics(info,meMomenta);
        return;
      }

  if ( !spacelike ) {

    Energy mij = momentum.mass();
    Energy mi,mj;
    
      // work out the mass for the first child
    if ( !children[0].children.empty() ) {
      Energy sumOthers = ZERO;
      for ( size_t k = 2; k < meMomenta.size(); ++k )
        if ( children[1].leafs.find(k) != children[1].leafs.end() )
          sumOthers += meMomenta[k].mass();
      children[0].massRange.second = momentum.mass() - sumOthers;
      if ( children[0].massRange.second < children[0].massRange.first )
        throw Veto();
      if ( children[0].massRange.second > momentum.mass() )
        throw Veto();
      mi = info.generateMass(children[0].data,children[0].massRange);
      children[0].momentum.setMass(mi);
    } else {
      mi = children[0].momentum.mass();
    }

      // work out the mass for the second child
    if ( !children[1].children.empty() ) {
      children[1].massRange.second = momentum.mass()-children[0].momentum.mass();
      if ( children[1].massRange.second < children[1].massRange.first )
        throw Veto();
      mj = info.generateMass(children[1].data,children[1].massRange);
      children[1].momentum.setMass(mj);
    } else {
      mj = children[1].momentum.mass();
    }

    Energy2 mij2 = sqr(mij);
    Energy2 mi2 = sqr(mi);
    Energy2 mj2 = sqr(mj);

      // perform the decay
    Energy4 lambda2 = sqr(mij2-mi2-mj2)-4.*mi2*mj2;
    if ( lambda2 <= ZERO )
      throw Veto();
    Energy2 lambda = sqrt(lambda2);
    double phi = 2.*Constants::pi*info.rnd();
    double cosPhi = cos(phi);
    double sinPhi = sqrt(1.-sqr(cosPhi));
    if ( phi > Constants::pi )
      sinPhi = -sinPhi;
    info.weight *= Constants::pi*lambda/(2.*mij2);
    double cosTheta = 2.*info.rnd() - 1.;
    double sinTheta = sqrt(1.-sqr(cosTheta));
    Energy p = lambda/(2.*mij);
    children[0].momentum.setX(p*cosPhi*sinTheta);
    children[0].momentum.setY(p*sinPhi*sinTheta);
    children[0].momentum.setZ(p*cosTheta);
    children[0].momentum.rescaleEnergy();
    if ( momentum.m2() <= ZERO )
      throw Veto();
    Boost out = momentum.boostVector();
    if ( out.mag2() > Constants::epsilon ) {
      children[0].momentum.boost(out);
    }
    children[1].momentum = momentum - children[0].momentum;
    children[1].momentum.setMass(mj);
    children[1].momentum.rescaleEnergy();

      // go on with next branchings
    children[0].generateKinematics(info,meMomenta);
    children[1].generateKinematics(info,meMomenta);

    return;

  }

    // get the minimum mass of the `W' system
  Energy Wmin = ZERO;
  PhasespaceTree* current = &children[0];
  while ( !(current->children.empty()) ) {
    Wmin += current->children[1].massRange.first;
    current = &(current->children[0]);
  }

    // get the CM energy avaialble
  Energy2 s = (momentum+backwardMomentum).m2();
  if ( s <= ZERO )
    throw Veto();

    // generate a mass for the timelike child
  Energy mi;
  if ( !children[1].children.empty() ) {
    children[1].massRange.second = sqrt(s)-Wmin;
    if ( children[1].massRange.second < children[1].massRange.first )
      throw Veto();
    mi = info.generateMass(children[1].data,children[1].massRange);
    children[1].momentum.setMass(mi);
  } else {
    mi = children[1].momentum.mass();
  }
  Energy2 mi2 = sqr(mi);
 
    // wether or not this is the last 2->2 scatter
  bool lastScatter = children[0].children[0].children.empty();

    // `W' mass relevant for the other boundaries
  Energy MW = Wmin;

    // generate a mass for second outgoing leg, if needed
  if ( lastScatter )
    if ( !children[0].children[1].children.empty() ) {
        // get the maximum `W' mass
      Energy Wmax = sqrt(s)-children[1].momentum.mass();
      children[0].children[1].massRange.second = Wmax;
      if ( children[0].children[1].massRange.second <
          children[0].children[1].massRange.first )
        throw Veto();
      MW = info.generateMass(children[0].children[1].data,
                             children[0].children[1].massRange);
      children[0].children[1].momentum.setMass(MW);
    }
  Energy2 MW2 = sqr(MW);

  Energy ma = momentum.mass();
  Energy2 ma2 = sqr(ma);
  if ( ma < ZERO )
    ma2 = -ma2;
  Energy mb = backwardMomentum.mass();
  Energy2 mb2 = sqr(mb);
  if ( mb < ZERO )
    mb2 = -mb2;

    // pick the ys variable
  Energy2 ys = ZERO;
  if ( !lastScatter ) {
    ys = info.rnd()*(sqr(sqrt(s)-mi)-MW2);
    info.weight *= (sqr(sqrt(s)-mi)-MW2)/info.sHat;
  }

  Energy4 lambda2 = sqr(s-ma2-mb2)-4.*ma2*mb2;
  if ( lambda2 <= ZERO ) {
    throw Veto();
  }
  Energy2 lambda = sqrt(lambda2);
  info.weight *= info.sHat/(4.*lambda);

    // get the boundaries on the momentum transfer
  Energy4 rho2 = sqr(s-ys-MW2-mi2)-4.*mi2*(ys+MW2);
  if ( rho2 < ZERO )
    throw Veto();
  Energy2 rho = sqrt(rho2);
  Energy4 tau2 =
  ys*(ma2-mb2+s)
  - sqr(s)+s*(ma2+mb2+mi2+MW2)-(mi2-MW2)*(ma2-mb2);
  pair<Energy2,Energy2> tBounds
  ((tau2-rho*lambda)/(2.*s),(tau2+rho*lambda)/(2.*s));
  children[0].massRange.first = sqrt(abs(tBounds.first));
  if ( tBounds.first < ZERO )
    children[0].massRange.first = -children[0].massRange.first;
  children[0].massRange.second = sqrt(abs(tBounds.second));
  if ( tBounds.second < ZERO )
    children[0].massRange.second = -children[0].massRange.second;

    // generate a momentum transfer
  Energy mai = info.generateMass(children[0].data,children[0].massRange);
  children[0].momentum.setMass(mai);
  Energy2 t = sqr(mai);
  if ( mai < ZERO )
    t = -t;

  Energy2 u = -s -t + ys + ma2 + mb2 + mi2 + MW2;
  Energy2 st = s - ma2 - mb2;
  Energy2 tt = t - mi2 - ma2;
  Energy2 ut = u - mi2 - mb2;

    // get the timelike momentum
  double xa = (-st*ut+2.*mb2*tt)/lambda2;
  double xb = (-st*tt+2.*ma2*ut)/lambda2;
  Energy2 pt2 = (st*tt*ut-ma2*sqr(ut)-mb2*sqr(tt)-mi2*sqr(st)+4.*ma2*mb2*mi2)/lambda2;
  if ( pt2 < ZERO )
    throw Veto();
  Energy pt = sqrt(pt2);

  children[1].momentum =
  xa*momentum + xb*backwardMomentum 
  + info.generateKt(momentum,backwardMomentum,pt);
  children[1].momentum.setMass(mi);
  children[1].momentum.rescaleEnergy();

  children[0].momentum =
  momentum - children[1].momentum;
  children[0].momentum.setMass(mai);
  bool changeSign = false;
  if ( children[0].momentum.t() < ZERO && mai > ZERO ) changeSign = true;
  if ( mai < ZERO )
    children[0].momentum.rescaleRho();
  else
    children[0].momentum.rescaleEnergy();
  if ( changeSign ) children[0].momentum.setT(-children[0].momentum.t());
 
  children[1].generateKinematics(info,meMomenta);

  if ( !lastScatter ) {
    children[0].backwardMomentum = backwardMomentum;
    children[0].generateKinematics(info,meMomenta);
  } else {
    children[0].children[1].momentum =
    backwardMomentum + children[0].momentum;
    children[0].children[1].momentum.setMass(MW);
    children[0].children[1].momentum.rescaleEnergy();
    children[0].children[1].generateKinematics(info,meMomenta);
  }

}


void PhasespaceTree::put(PersistentOStream& os) const {
  os << children.size();
  if ( !children.empty() ) {
    children[0].put(os);
    children[1].put(os);
  }
  os << data << externalId << leafs << spacelike << doMirror;
}

void PhasespaceTree::get(PersistentIStream& is) {
  size_t nc; is >> nc;
  assert(nc == 0 || nc == 2);
  if ( nc == 2 ) {
    children.resize(2,PhasespaceTree());
    children[0].get(is);
    children[1].get(is);
  }
  is >> data >> externalId >> leafs >> spacelike >> doMirror;
}
#line 1 "./TildeKinematics.cc"
// -*- C++ -*-
//
// TildeKinematics.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TildeKinematics class.
//

#include "TildeKinematics.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

Energy TildeKinematics::lastScale() const {
  if ( ( theDipole->bornEmitter() < 2 && theDipole->bornSpectator() > 1 ) ||
       ( theDipole->bornEmitter() > 1 && theDipole->bornSpectator() < 2 ) ) {
    return -(bornEmitterMomentum()-bornSpectatorMomentum()).m();
  }
  return (bornEmitterMomentum()+bornSpectatorMomentum()).m();
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void TildeKinematics::rebind(const TranslationMap & trans) {
  theDipole = trans.translate(theDipole);
  HandlerBase::rebind(trans);
}

IVector TildeKinematics::getReferences() {
  IVector ret = HandlerBase::getReferences();
  ret.push_back(theDipole);
  return ret;
}



void TildeKinematics::persistentOutput(PersistentOStream & os) const {
  os << theDipole << theRealXComb << theBornXComb
     << ounit(theBornEmitterMomentum,GeV) << ounit(theBornSpectatorMomentum,GeV);
}

void TildeKinematics::persistentInput(PersistentIStream & is, int) {
  is >> theDipole >> theRealXComb >> theBornXComb
     >> iunit(theBornEmitterMomentum,GeV) >> iunit(theBornSpectatorMomentum,GeV);
}

void TildeKinematics::Init() {

  static ClassDocumentation<TildeKinematics> documentation
    ("TildeKinematics is the base class for the 'tilde' "
     "kinematics being used for subtraction terms in the "
     "formalism of Catani and Seymour.");

}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<TildeKinematics,HandlerBase>
describeTildeKinematics("Herwig::TildeKinematics", "Herwig.so");
#line 1 "./TreePhasespace.cc"
// -*- C++ -*-
//
// TreePhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TreePhasespace class.
//

#include <sstream> 
#include <string> 
#include "TreePhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/DiagramDrawer.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;

TreePhasespace::TreePhasespace() 
  : x0(0.01), xc(1e-4), M0(ZERO), Mc(ZERO) {
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
  theIncludeMirrored = true;
}

IBPtr TreePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr TreePhasespace::fullclone() const {
  return new_ptr(*this);
}

void TreePhasespace::setXComb(tStdXCombPtr xco) {

  MatchboxPhasespace::setXComb(xco);

  lastChannelsIterator = channelMap().find(lastXCombPtr());

  if ( lastChannelsIterator == channelMap().end() ) {
    map<Ptr<Tree2toNDiagram>::ptr,pair<PhasespaceTree, PhasespaceTree> > channels;
    for ( auto const & d : lastXComb().diagrams()) {
      PhasespaceTree tree;
      Ptr<Tree2toNDiagram>::ptr diag =
	dynamic_ptr_cast<Ptr<Tree2toNDiagram>::ptr>(d);
      tree.setup(*diag);
      PhasespaceTree treeMirror;
      treeMirror.setupMirrored(*diag, diag->nSpace() - 1);
      channels[diag] = make_pair(tree,treeMirror);
    }
    channelMap()[lastXCombPtr()] = channels;
    lastChannelsIterator = channelMap().find(lastXCombPtr());
  }

}

double TreePhasespace::generateTwoToNKinematics(const double* random,
						vector<Lorentz5Momentum>& momenta) {

  lastPhasespaceInfo.sHat = lastXComb().lastSHat();
  lastPhasespaceInfo.sqrtSHat = sqrt(lastXComb().lastSHat());
  lastPhasespaceInfo.weight = 1.;

  size_t nchannels = lastXComb().diagrams().size();
  bool doMirror = (UseRandom::rnd() < 0.5) && theIncludeMirrored;
  map<Ptr<Tree2toNDiagram>::ptr,
      pair <PhasespaceHelpers::PhasespaceTree, PhasespaceHelpers::PhasespaceTree> >::iterator ds =
    lastChannels().begin();

  size_t i = (size_t)(random[0]*nchannels);
  advance(ds,i);

  Ptr<Tree2toNDiagram>::ptr channel = ds->first;
  ++random;
  
  lastPhasespaceInfo.rnd.numbers = random;
  lastPhasespaceInfo.rnd.nRnd = 3*momenta.size() - 10;
    
  try {
    if ( !doMirror )
      lastChannels()[channel].first.generateKinematics(lastPhasespaceInfo,momenta);
    else 
      lastChannels()[channel].second.generateKinematics(lastPhasespaceInfo,momenta);
  } catch (Veto) {
    return 0.;
  }
    
  if ( !matchConstraints(momenta) )
    return 0.;

  double flatCut = x0;
  if ( M0 != ZERO )
    flatCut = M0/sqrt(lastSHat());

  fillDiagramWeights(flatCut);

  double sum = 0.;
  for ( auto const & d : lastChannels())
    sum += diagramWeight(*(d.first));

  double piWeight = pow(2.*Constants::pi,(double)(3*(momenta.size()-2)-4));

  for ( auto & k : momenta )
    k.rescaleRho();

  return nchannels*lastPhasespaceInfo.weight*diagramWeight(*channel)/(sum*piWeight);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

void TreePhasespace::doinit() {
  MatchboxPhasespace::doinit();
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
}

void TreePhasespace::doinitrun() {
  MatchboxPhasespace::doinitrun();
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
}

void TreePhasespace::persistentOutput(PersistentOStream & os) const {
  os << theChannelMap << x0 << xc 
     << ounit(M0,GeV) << ounit(Mc,GeV)
     << theIncludeMirrored
     << theLastXComb;
}

void TreePhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theChannelMap >> x0 >> xc 
     >> iunit(M0,GeV) >> iunit(Mc,GeV)
     >> theIncludeMirrored
     >> theLastXComb;
  lastPhasespaceInfo.x0 = x0;
  lastPhasespaceInfo.xc = xc;
  lastPhasespaceInfo.M0 = M0;
  lastPhasespaceInfo.Mc = Mc;
  lastChannelsIterator = channelMap().find(lastXCombPtr());
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TreePhasespace,MatchboxPhasespace>
  describeHerwigTreePhasespace("Herwig::TreePhasespace", "Herwig.so");

void TreePhasespace::Init() {

  static ClassDocumentation<TreePhasespace> documentation
    ("TreePhasespace is a multi-channel phase space generator "
     "adapting to singularity structures as determined from the matrix "
     "elements diagrams.");


  static Reference<TreePhasespace,TreePhasespaceChannels> interfaceChannelMap
    ("ChannelMap",
     "Set the object storing the channels.",
     &TreePhasespace::theChannelMap, false, false, true, false, false);
  interfaceChannelMap.rank(-1);


  static Parameter<TreePhasespace,double> interfaceX0
    ("X0",
     "Set the cut below which flat virtuality sampling is imposed.",
     &TreePhasespace::x0, 0.01, 0.0, 0,
     false, false, Interface::lowerlim);


  static Parameter<TreePhasespace,double> interfaceXC
    ("XC",
     "Set the cut below which no virtualities are generated.",
     &TreePhasespace::xc, 1e-4, 0.0, 0,
     false, false, Interface::lowerlim);


  static Parameter<TreePhasespace,Energy> interfaceM0
    ("M0",
     "Set the cut below which flat virtuality sammpling is imposed.",
     &TreePhasespace::M0, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Parameter<TreePhasespace,Energy> interfaceMC
    ("MC",
     "Set the cut below which no virtualities are generated.",
     &TreePhasespace::Mc, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);

  static Switch<TreePhasespace,bool> interfaceIncludeMirrored
    ("IncludeMirrored",
     "Choose whether to include mirrored diagrams for PS generation",
     &TreePhasespace::theIncludeMirrored, true, true, false);
  static SwitchOption interfaceIncludeMirroredYes
    (interfaceIncludeMirrored,
     "Yes",
     "Use unmirrored and mirrored diagrams",
     true);
  static SwitchOption interfaceIncludeMirroredNo
    (interfaceIncludeMirrored,
     "No",
     "Use only unmirrored diagrams",
     false);
  interfaceIncludeMirrored.rank(-1);

}

#line 1 "./TreePhasespaceChannels.cc"
// -*- C++ -*-
//
// TreePhasespaceChannels.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TreePhasespaceChannels class.
//

#include "TreePhasespaceChannels.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace Herwig::PhasespaceHelpers;

IBPtr TreePhasespaceChannels::clone() const {
  return new_ptr(*this);
}

IBPtr TreePhasespaceChannels::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TreePhasespaceChannels::persistentOutput(PersistentOStream & os) const {
  os << theChannelMap.size();
  for ( map<tStdXCombPtr,map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceTree,PhasespaceTree> > >::const_iterator k =
          theChannelMap.begin(); k != theChannelMap.end(); ++k ) {
    os << k->first << k->second.size();
    for ( map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceTree,PhasespaceTree> >::const_iterator l = 
            k->second.begin(); l != k->second.end(); ++l ) {
      os << l->first;
      l->second.first.put(os);
      l->second.second.put(os);
    }
  }
}

void TreePhasespaceChannels::persistentInput(PersistentIStream & is, int) {
  size_t nk; is >> nk;
  for ( size_t k = 0; k < nk; ++k ) {
    tStdXCombPtr xc; is >> xc;
    size_t nl; is >> nl;
    map<Ptr<Tree2toNDiagram>::ptr,pair <PhasespaceTree,PhasespaceTree> > cm;
    for ( size_t l = 0; l < nl; ++l ) {
      Ptr<Tree2toNDiagram>::ptr ci; is >> ci;
      pair<PhasespaceTree,PhasespaceTree> cp; cp.first.get(is); cp.second.get(is);
      cm[ci] = cp;
    }
    theChannelMap[xc] = cm;
  }
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<TreePhasespaceChannels,HandlerBase>
  describeHerwigTreePhasespaceChannels("Herwig::TreePhasespaceChannels", "Herwig.so");

void TreePhasespaceChannels::Init() {

  static ClassDocumentation<TreePhasespaceChannels> documentation
    ("Store channels for the tree phase space.");

}

#line 1 "./MatchboxReference.cc"
// -*- C++ -*-
//
// MatchboxReference.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxReference class.
//

#include "MatchboxReference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr MatchboxReference::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxReference::fullclone() const {
  return new_ptr(*this);
}

double MatchboxReference::generateTwoToNKinematics(const double*,
						   vector<Lorentz5Momentum>& momenta) {


  map<cPDVector,ifstream*>::iterator ref =
    referenceSamples.find(mePartonData());
  if ( ref == referenceSamples.end() ) {
    ostringstream refname;
    for ( cPDVector::const_iterator p = mePartonData().begin();
	  p != mePartonData().end(); ++p ) {
      refname << (**p).PDGName();
    }
    refname << ".rambo";
    referenceSamples[mePartonData()] = new ifstream(refname.str().c_str());
    ref = referenceSamples.find(mePartonData());
  }
  assert(ref != referenceSamples.end());

  ifstream& in = *(ref->second);
  assert(in);

  double x1,x2;
  double x,y,z,t,m;
  double weight;

  in >> x1 >> x2;
  for ( vector<Lorentz5Momentum>::iterator p = momenta.begin();
	p != momenta.end(); ++p ) {
    in >> x >> y >> z >> t >> m;
    *p = Lorentz5Momentum(x*GeV,y*GeV,z*GeV,t*GeV,m*GeV);
  }

  in >> weight;

  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat((momenta[0]+momenta[1]).m2());

  return weight;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxReference::persistentOutput(PersistentOStream &) const {}

void MatchboxReference::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxReference,MatchboxPhasespace>
  describeHerwigMatchboxReference("Herwig::MatchboxReference", "Herwig.so");

void MatchboxReference::Init() {

  static ClassDocumentation<MatchboxReference> documentation
    ("MatchboxReference implements reference sample phase space generation.");

}

#line 1 "./FlatInvertiblePhasespace.cc"
// -*- C++ -*-
//
// FlatInvertiblePhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatInvertiblePhasespace class.
//

#include "FlatInvertiblePhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr FlatInvertiblePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr FlatInvertiblePhasespace::fullclone() const {
  return new_ptr(*this);
}

double FlatInvertiblePhasespace::bisect(double v, double n, 
					double target, double maxLevel) const {

  if ( v != 0.0 && v != 1.0 ) {

    double level = 0;
    double left = 0;
    double right = 1;

    double checkV = -1.;
    double u = -1;

    while ( level < maxLevel ) {

      u = (left+right)*pow(0.5,level+1.);
      checkV = 
	pow(u,n+1.)*(n+2.-(n+1.)*u);

      if ( log10(abs(1.-checkV/v)) <= target )
	break;

      left *= 2.;
      right *= 2.;

      if ( v <= checkV ) {
	right -= 1.;
	++level;
      }

      if ( v > checkV ) {
	left += 1.;
	++level;
      }

    }

    return u;

  }

  return v;

}

double FlatInvertiblePhasespace::generateIntermediates(vector<Energy>& K,
						       const double* r) const {

  size_t n = K.size() + 1;

  
  for ( size_t i = 2; i <= n-1; ++i ) {
    double u = bisect(r[i-2],n-1-i);
    K[i-1] = sqrt(u*sqr(K[i-2]));
  } 
  
  int kap = K.size() + 1;
  
  return flatWeights(kap);

}

double FlatInvertiblePhasespace::invertIntermediates(const vector<Energy>& K,
						     double* r) const {

  size_t n = K.size() + 1;
  for ( size_t i = 2; i <= n-1; ++i ) {
    double u = sqr(K[i-1]/K[i-2]);
    r[i-2] = (n+1-i)*pow(u,(double)(n-i)) - (n-i)*pow(u,(double)(n+1-i));
  } 

  int kap = K.size() + 1;
  return flatWeights(kap);

}

double FlatInvertiblePhasespace::generateIntermediates(vector<Energy>& M,
						       const vector<Energy>& m,
						       const double* r) const {

  size_t n = M.size() + 1;

  vector<Energy> K = M;
  for ( size_t i = 1; i <= n; ++i )
    K[0] -= m[i-1];

  double w0 = generateIntermediates(K,r);

  M = K;
  for ( size_t i = 1; i <= n-1; ++i ) {
    for ( size_t k = i; k <= n; ++k )
      M[i-1] += m[k-1];
  }

  double weight = 8.*w0*rho(M[n-2],m[n-1],m[n-2]);

  for ( size_t i = 2; i <= n-1; ++i ) {
    weight *= 
      (rho(M[i-2],M[i-1],m[i-2])/rho(K[i-2],K[i-1],ZERO)) * (M[i-1]/K[i-1]);
  }

  weight *= pow(K[0]/M[0],2.*n-4.);
  
  return weight;

}

double FlatInvertiblePhasespace::invertIntermediates(const vector<Energy>& M,
						     const vector<Energy>& m,
						     double* r) const {

  size_t n = M.size() + 1;

  vector<Energy> K = M;
  for ( size_t i = 1; i <= n-1; ++i ) {
    for ( size_t k = i; k <= n; ++k )
      K[i-1] -= m[k-1];
  }

  double w0 = invertIntermediates(K,r);

  double weight = 8.*w0*rho(M[n-2],m[n-1],m[n-2]);

  for ( size_t i = 2; i <= n-1; ++i ) {
    weight *= 
      (rho(M[i-2],M[i-1],m[i-2])/rho(K[i-2],K[i-1],ZERO)) * (M[i-1]/K[i-1]);
  }

  weight *= pow(K[0]/M[0],2.*n-4.);

  return weight;

}


double FlatInvertiblePhasespace::generateKinematics(vector<Lorentz5Momentum>& P,
						    Energy Ecm,
						    const double* r) const {

  vector<Energy> m;
  for ( vector<Lorentz5Momentum>::const_iterator p =
	  P.begin() + 2; p != P.end(); ++p )
    m.push_back(p->mass());

  size_t n = P.size() - 2;
  vector<Energy> M(n-1);
  M[0] = Ecm;

  double weight = generateIntermediates(M,m,r);

  M.push_back(m.back());

  Lorentz5Momentum Q(M[0]);
  Lorentz5Momentum nextQ;

  for ( size_t i = 2; i <= n; ++i ) {

    Energy q = 4.*M[i-2]*rho(M[i-2],M[i-1],m[i-2]);

    double c = 2.*r[n-6+2*i]-1.;
    double s = sqrt(1.-sqr(c));
    double phi = 2.*Constants::pi*r[n-5+2*i];
    double cphi = cos(phi);
    double sphi = sqrt(1.-sqr(cphi));
    if ( phi > Constants::pi )
      sphi = -sphi;

    P[i].setX(q*cphi*s);
    P[i].setY(q*sphi*s);
    P[i].setZ(q*c);
    P[i].rescaleEnergy();
    P[i].boost(Q.boostVector());
    P[i].rescaleEnergy();

    nextQ = Q - P[i];
    nextQ.setMass(M[i-1]);
    nextQ.rescaleEnergy();

    Q = nextQ;

  }

  P.back() = Q;

  return weight;

}

double FlatInvertiblePhasespace::invertKinematics(const vector<Lorentz5Momentum>& P,
						  Energy Ecm,
						  double* r) const {

  vector<Energy> m;
  for ( vector<Lorentz5Momentum>::const_iterator p =
	  P.begin() + 2; p != P.end(); ++p )
    m.push_back(p->mass());

  size_t n = P.size() - 2;
  vector<Energy> M(n-1);
  M[0] = Ecm;

  vector<Lorentz5Momentum> Q(n-1);
  Q[0] = Lorentz5Momentum(M[0]);

  for ( size_t i = 2; i <= n-1; ++i ) {
    for ( size_t k = i; k <= n; ++k )
      Q[i-1] += P[k+1];
    M[i-1] = Q[i-1].m();
  }

  double weight = invertIntermediates(M,m,r);

  for ( size_t i = 2; i <= n; ++i ) {
    Lorentz5Momentum p = P[i];
    p.boost(-Q[i-2].boostVector());
    r[n-6+2*i] = (p.cosTheta()+1.)/2.;
    double phi = p.phi();
    if ( phi < 0. )
      phi = 2.*Constants::pi + phi;
    r[n-5+2*i] = phi/(2.*Constants::pi);
  }

  return weight;

}


double FlatInvertiblePhasespace::generateTwoToNKinematics(const double* r,
							  vector<Lorentz5Momentum>& momenta) {

  double weight = generateKinematics(momenta,sqrt(lastXCombPtr()->lastSHat()),r);
  return weight;

}

long double FlatInvertiblePhasespace::flatWeights(int k) const{
  using Constants::pi;
  if(k<2) { return -1;
  } else return pow((pi/2),(k-1)) * pow((2*pi),(4-3*k))/factorial(k-1)/factorial(k-2);
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FlatInvertiblePhasespace::persistentOutput(PersistentOStream &) const {}

void FlatInvertiblePhasespace::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FlatInvertiblePhasespace,MatchboxPhasespace>
  describeHerwigFlatInvertiblePhasespace("Herwig::FlatInvertiblePhasespace", "Herwig.so");

void FlatInvertiblePhasespace::Init() {

  static ClassDocumentation<FlatInvertiblePhasespace> documentation
    ("FlatInvertiblePhasespace implements flat, invertible phase space generation.");

}

#line 1 "./FlatInvertibleLabframePhasespace.cc"
// -*- C++ -*-
//
// FlatInvertiblePhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatInvertibleLabframePhasespace class.
//

#include "FlatInvertibleLabframePhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

FlatInvertibleLabframePhasespace::FlatInvertibleLabframePhasespace()
  : theLogSHat(false) {}

IBPtr FlatInvertibleLabframePhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr FlatInvertibleLabframePhasespace::fullclone() const {
  return new_ptr(*this);
}

double FlatInvertibleLabframePhasespace::invertTwoToNKinematics(const vector<Lorentz5Momentum>& momenta,
								double* r) const {

  double weight = 1.;

  Energy finalstatemass = 0*GeV;
  for ( vector<Lorentz5Momentum>::const_iterator p =
        momenta.begin()+2; p != momenta.end(); ++p )
    finalstatemass += p->mass();

  Lorentz5Momentum pinitial = momenta[0]+momenta[1];
  Energy2 sh = pinitial.m2();
  double tau = sh/lastS();
  Energy2 shmax = lastCuts().sHatMax();
  Energy2 shmin = max(lastCuts().sHatMin(),sqr(finalstatemass));
  if (theLogSHat) {
    r[0] = log(sh/shmin)/log(shmax/shmin);
    weight *= tau*log(shmax/shmin);
  } else {
    r[0] = (sh-shmin)/(shmax-shmin);
    weight *= (shmax-shmin)/lastS();
  }
  double ltau = log(tau);
  r[1] = 0.5 - pinitial.rapidity()/ltau;
  weight *= -ltau;

  vector<Lorentz5Momentum> Pcms = momenta;
  Boost toCMS = pinitial.findBoostToCM();
  for ( vector<Lorentz5Momentum>::iterator pit =
        Pcms.begin(); pit != Pcms.end(); ++pit )
    pit->boost(toCMS);
  
  weight *= FlatInvertiblePhasespace::invertTwoToNKinematics(Pcms, r+2);

  return weight;

}


double FlatInvertibleLabframePhasespace::generateTwoToNKinematics(const double* r,
							          vector<Lorentz5Momentum>& momenta) {

  double weight = 1.;

  Energy finalstatemass = 0*GeV;
  for ( vector<Lorentz5Momentum>::const_iterator p =
        momenta.begin()+2; p != momenta.end(); ++p )
    finalstatemass += p->mass();

  Energy beamenergy = sqrt(lastS())/2.;
  Energy2 shmax = lastCuts().sHatMax();
  Energy2 shmin = max(lastCuts().sHatMin(),sqr(finalstatemass));
  Energy2 sh;
  double tau; 
  if (theLogSHat) {
    sh = shmin*pow(shmax/shmin, r[0]);
    tau = sh/lastS(); 
    weight *= tau*log(shmax/shmin);
  } else {
    sh = r[0]*(shmax-shmin)+shmin;
    tau = sh/lastS(); 
    weight *= (shmax-shmin)/lastS();
  }
  double ltau = log(tau);
  double y = ltau*(0.5 - r[1]);
  weight *= -ltau;

  double x1 = sqrt(tau)*exp(y);
  double x2 = sqrt(tau)*exp(-y);
  momenta[0] = Lorentz5Momentum(0*GeV,0*GeV,+x1*beamenergy,x1*beamenergy);
  momenta[1] = Lorentz5Momentum(0*GeV,0*GeV,-x2*beamenergy,x2*beamenergy);
  lastXCombPtr()->lastX1X2(make_pair(x1,x2));
  lastXCombPtr()->lastSHat(sh);

  weight *= FlatInvertiblePhasespace::generateTwoToNKinematics(r+2, momenta);

  // find boost to the relevant partonic frame note final state kinematics are
  // always generated in the CMS for this phase space algorithm
  Boost boostinitial = (momenta[0]+momenta[1]).findBoostToCM();
  for ( vector<Lorentz5Momentum>::iterator pit =
          momenta.begin()+2; pit != momenta.end(); ++pit )
    pit->boost(-boostinitial);

  return weight;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FlatInvertibleLabframePhasespace::persistentOutput(PersistentOStream & os) const {
  os << theLogSHat;
}

void FlatInvertibleLabframePhasespace::persistentInput(PersistentIStream & is, int) {
  is >> theLogSHat;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FlatInvertibleLabframePhasespace,MatchboxPhasespace>
  describeHerwigFlatInvertibleLabframePhasespace("Herwig::FlatInvertibleLabframePhasespace", "Herwig.so");

void FlatInvertibleLabframePhasespace::Init() {

  static ClassDocumentation<FlatInvertibleLabframePhasespace> documentation
    ("FlatInvertibleLabframePhasespace implements flat, invertible phase space generation in the lab frame.");

  static Switch<FlatInvertibleLabframePhasespace,bool> interfaceLogSHat
    ("LogSHat",
     "Generate a flat distribution in \\f$\\log(\\hat{s})\\f$.",
     &FlatInvertibleLabframePhasespace::theLogSHat, false, false, false);

  static SwitchOption interfaceLogSHatYes
    (interfaceLogSHat,
     "Yes", "Generate flat in \\f$\\log(\\hat{s})\\f$", true);

  static SwitchOption interfaceLogSHatNo
    (interfaceLogSHat,
     "No", "Generate flat in \\f$\\hat{s}\\f$", false);

  interfaceLogSHat.rank(-1);

}

#line 1 "./PhasespaceCouplings.cc"
// -*- C++ -*-
//
// PhasespaceCouplings.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PhasespaceCouplings class.
//

#include "PhasespaceCouplings.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr PhasespaceCouplings::clone() const {
  return new_ptr(*this);
}

IBPtr PhasespaceCouplings::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void PhasespaceCouplings::persistentOutput(PersistentOStream & os) const {
  os << theCouplings.size();
  for ( map<LTriple,double>::const_iterator cit = 
	  theCouplings.begin(); cit != theCouplings.end(); ++cit )
    os << cit->first << cit->second;
}

void PhasespaceCouplings::persistentInput(PersistentIStream & is, int) {
  theCouplings.clear();
  size_t size;
  LTriple k;
  is >> size;
  while ( size-- && is ) {
    is >> k;
    is >> theCouplings[k];
  }
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<PhasespaceCouplings,HandlerBase>
  describeHerwigPhasespaceCouplings("Herwig::PhasespaceCouplings", "Herwig.so");

void PhasespaceCouplings::Init() {

  static ClassDocumentation<PhasespaceCouplings> documentation
    ("Store couplings for the phase space generator.");

}

