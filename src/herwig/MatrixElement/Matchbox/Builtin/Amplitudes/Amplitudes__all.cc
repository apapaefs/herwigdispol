#line 1 "./MatchboxCurrents.cc"
// - * - C++ - * -
//
// MatchboxCurrents.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "MatchboxCurrents.h"
#include "Herwig/Utilities/Maths.h"

using namespace Herwig;
using namespace Herwig::Math;

using Constants::pi;

namespace {

  static const LorentzVector<Complex> czero(0.,0.,0.,0.);

  inline Complex csqr(const Complex & a) {
    return a * a;
  }

  inline double theta(const double x) {
    if ( x >= 0. )
      return 1.;
    return 0.;
  }
  inline double sign(const double x) {
    if ( x >= 0. )
      return 1.;
    return -1.;
  }

}


void MatchboxCurrents::setupLeptons(const int l,    const Lorentz5Momentum& pl,
				    const int lbar, const Lorentz5Momentum& plbar) {

  const Energy4 Delta = (sqr(pl*plbar) - (pl*pl)*(plbar*plbar));
  const Energy2 prod = pl*plbar;

  // Variable to contain the sign of pl*plbar
  double sgn;
  if (prod < ZERO ) {sgn = -1;}
  else if (prod > ZERO) {sgn = 1;}
  else {sgn = 0;}

  InvEnergy2 fact = 0.5/(sgn*sqrt(Delta));

  Lorentz5Momentum lmassless    = ( double(fact*(sgn*sqrt(Delta) + prod))*pl    - double(fact*(   pl*pl))*plbar );
  Lorentz5Momentum lbarmassless = ( double(fact*(sgn*sqrt(Delta) + prod))*plbar - double(fact*(plbar*plbar))*pl );

  lmassless.setMass(ZERO); lmassless.rescaleEnergy();
  lbarmassless.setMass(ZERO); lbarmassless.rescaleEnergy();

  if ( pl.t() < ZERO )
    lmassless.setT(-lmassless.t());

  if ( plbar.t() < ZERO )
    lbarmassless.setT(-lbarmassless.t());

  momentum(l,lmassless,true,pl.mass());
  momentum(lbar,lbarmassless,true,plbar.mass());
}


void MatchboxCurrents::setupQuarks(const int q,    const Lorentz5Momentum& pq,
				   const int qbar, const Lorentz5Momentum& pqbar) {

  const Energy4 Delta = (sqr(pq*pqbar) - (pq*pq)*(pqbar*pqbar));
  const Energy2 prod = pq*pqbar;

  // Variable to contain the sign of pq*pqbar
  double sgn;
  if (prod < ZERO) {sgn = -1;}
  else if (prod > ZERO) {sgn = 1;}
  else {sgn = 0;}
  
  InvEnergy2 fact =  0.5/(sgn*sqrt(Delta));

  Lorentz5Momentum qmassless    = ( double(fact*(sgn*sqrt(Delta) + prod))*pq    - double(fact*(pq*pq))*pqbar );
  Lorentz5Momentum qbarmassless = ( double(fact*(sgn*sqrt(Delta) + prod))*pqbar - double(fact*(pqbar*pqbar))*pq );

  qmassless.setMass(ZERO); qmassless.rescaleEnergy();
  qbarmassless.setMass(ZERO); qbarmassless.rescaleEnergy();

  if ( pq.t() < ZERO )
    qmassless.setT(-qmassless.t());

  if ( pqbar.t() < ZERO )
    qbarmassless.setT(-qbarmassless.t());

  momentum(q,qmassless,true,pq.mass());
  momentum(qbar,qbarmassless,true,pqbar.mass());
}


const LorentzVector<Complex>& MatchboxCurrents::llbarLeftCurrent(const int l,    const int lHel,
								 const int lbar, const int lbarHel) {

  if ( getCurrent(hash<0>(1,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,2.) * mass(lbar)/plusProduct(l,lbar)) * momentum(l));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,-2.) * mass(l)/minusProduct(l,lbar)) * momentum(lbar));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,1.) * mass(l) * mass(lbar)/invariant(l,lbar)) * minusCurrent(lbar,l));
  }

  return cachedCurrent();
}

const LorentzVector<Complex>& MatchboxCurrents::llbarRightCurrent(const int l,    const int lHel,
								  const int lbar, const int lbarHel) {
    
  if ( getCurrent(hash<0>(2,1,l,lHel,lbar,lbarHel)) ) {
    if ( lHel == 1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,1.) * mass(l) * mass(lbar)/invariant(l,lbar)) * minusCurrent(l,lbar));
    if ( lHel == 1 && lbarHel == -1 )
      cacheCurrent((Complex(0.,-2.) * mass(l)/plusProduct(l,lbar)) * momentum(lbar));
    if ( lHel == -1 && lbarHel == 1 )
      cacheCurrent((Complex(0.,2.) * mass(lbar)/minusProduct(l,lbar)) * momentum(l));
    if ( lHel == -1 && lbarHel == -1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(lbar,l));
  }

  return cachedCurrent();
}


const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftCurrent(const int q,    const int qHel,
								 const int qbar, const int qbarHel) {

  if ( getCurrent(hash<1>(1,1,q,qHel,qbar,qbarHel)) ) {
    if ( qHel == 1 && qbarHel == 1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(q,qbar));
    if ( qHel == 1 && qbarHel == -1 )
      cacheCurrent((Complex(0.,2.) * mass(qbar)/plusProduct(q,qbar)) * momentum(q));
    if ( qHel == -1 && qbarHel == 1 )
      cacheCurrent((Complex(0.,-2.) * mass(q)/minusProduct(q,qbar)) * momentum(qbar));
    if ( qHel == -1 && qbarHel == -1 )
      cacheCurrent((Complex(0.,1.) * mass(q) * mass(qbar)/invariant(q,qbar)) * minusCurrent(qbar,q));
  }

#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  
  return cachedCurrent();
}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightCurrent(const int q,    const int qHel,
								  const int qbar, const int qbarHel) {
    
  if ( getCurrent(hash<1>(2,1,q,qHel,qbar,qbarHel)) ) {
    if ( qHel == 1 && qbarHel == 1 )
      cacheCurrent((Complex(0.,1.) * mass(q) * mass(qbar)/invariant(q,qbar)) * minusCurrent(q,qbar));
    if ( qHel == 1 && qbarHel == -1 )
      cacheCurrent((Complex(0.,-2.) * mass(q)/plusProduct(q,qbar)) * momentum(qbar));
    if ( qHel == -1 && qbarHel == 1 )
      cacheCurrent((Complex(0.,2.) * mass(qbar)/minusProduct(q,qbar)) * momentum(q));
    if ( qHel == -1 && qbarHel == -1 )
      cacheCurrent(Complex(0.,1.) * minusCurrent(qbar,q));
  }
 
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif
  
  return cachedCurrent();
}


const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftCurrent(const int q,    const int qHel,
								  const int qbar, const int qbarHel,
								  const int g,    const int gHel) {

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) {

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);

      // 2*factor from the spinor definition of the negative helicity gluon 
      // Note that the gluon is outgoing so the polarisation vector of the hel=+1 gluon is conjugated to give the hel=-1 vector
      const Complex cminus = sqrt(2.0) / minusProduct(g,q);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (minusProduct(qbar,q)*plusProduct(g,qbar)/den_j))*minusCurrent(q, qbar) - (minusProduct(g,q)*plusProduct(g,qbar)/den_j)*minusCurrent(q,g) ) ); 
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*(-mass(qbar)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (plusProduct(qbar,g)*minusProduct(q,qbar)/den_j))*2.*momentum(q) + (invariant(q,g)/den_j)*minusCurrent(q,g) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(q)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (plusProduct(g,qbar)*minusProduct(qbar,q)/den_j))*2.*momentum(qbar) - (minusProduct(g,q)*plusProduct(g,qbar)/den_j)*minusCurrent(qbar,g) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (minusProduct(q,qbar)*plusProduct(qbar,g)/den_j))*minusCurrent(qbar,q) + (invariant(q,g)/den_j)*minusCurrent(qbar,g) ) );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(1,1,q,qHel,qbar,qbarHel,g,gHel)) ) { 

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);

      // 2*factor from the spinor definition of the positive helicity gluon
      const Complex cplus = sqrt(2.0) / plusProduct(q,g);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (minusProduct(qbar,g)*plusProduct(q,qbar)/den_j))*minusCurrent(q, qbar) - (invariant(q,g)/den_i)*minusCurrent(g,qbar) ) ); 
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*(-mass(qbar)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (plusProduct(qbar,q)*minusProduct(g,qbar)/den_j))*2.*momentum(q) - (invariant(q,g)/den_i)*minusCurrent(g,q) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(q)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (minusProduct(qbar,g)*plusProduct(q,qbar)/den_j))*2.*momentum(qbar) + (minusProduct(qbar,g)*plusProduct(q,g)/den_i)*minusCurrent(g,qbar) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (plusProduct(qbar,q)*minusProduct(g,qbar)/den_j))*minusCurrent(qbar, q) + (minusProduct(qbar,g)*plusProduct(q,g)/den_i)*minusCurrent(g,q) ) );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  return czero;
}

const LorentzVector<Complex>& MatchboxCurrents::qqbargRightCurrent(const int q,    const int qHel,
								   const int qbar, const int qbarHel,
								   const int g,    const int gHel) {

  if ( gHel == 1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);
     
      // 2*factor from the spinor definition of the positive helicity gluon
      const Complex cminus = sqrt(2.0) / minusProduct(g,q);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (minusProduct(qbar,q)*plusProduct(g,qbar)/den_j))*plusCurrent(qbar, q) + (plusProduct(qbar,g)*minusProduct(q,g)/den_i)*plusCurrent(g,q) ) );
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*(mass(q)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(qbar,g)/(plusProduct(qbar,q)*den_i)) - (plusProduct(qbar,g)*minusProduct(q,qbar)/den_j))*2.*momentum(qbar) + (plusProduct(qbar,g)*minusProduct(q,g)/den_i)*plusCurrent(g,qbar) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cminus*(-mass(qbar)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (minusProduct(qbar,q)*plusProduct(g,qbar)/den_j))*2.*momentum(q) - (invariant(q,g)/den_i)*plusCurrent(g,q) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cminus*( ((sqr(mass(q))*plusProduct(g,qbar)/(plusProduct(q,qbar)*den_i)) - (plusProduct(qbar,g)*minusProduct(q,qbar)/den_j))*plusCurrent(q, qbar) - (invariant(q,g)/den_i)*plusCurrent(g,qbar) ) ); 
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  if ( gHel == -1 ) {
    if ( getCurrent(hash<2>(2,1,q,qHel,qbar,qbarHel,g,gHel)) ) {

      // Invariant products from propagator denominators
      const Complex den_i = invariant(q,g) + (sqr(mass(q))/invariant(q,qbar))*invariant(qbar,g);
      const Complex den_j = invariant(qbar,g) + (sqr(mass(qbar))/invariant(q,qbar))*invariant(q,g);

      // 2*factor from the spinor definition of the positive helicity gluon
      const Complex cplus = sqrt(2.0) / plusProduct(q,g);

      if ( qHel == 1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(qbar)*mass(q)/(invariant(q,qbar))) * ( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (plusProduct(q,qbar)*minusProduct(qbar,g)/den_j))*plusCurrent(qbar, q) + (invariant(q,g)/den_j)*plusCurrent(qbar,g) ) );
      if ( qHel == 1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*(mass(q)/plusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(g,qbar)/(minusProduct(q,qbar)*den_i)) - (minusProduct(g,qbar)*plusProduct(qbar,q)/den_j))*2.*momentum(qbar) - (plusProduct(g,q)*minusProduct(g,qbar)/den_j)*plusCurrent(qbar,g) ) );
      if ( qHel == -1 && qbarHel == 1 )
	cacheCurrent( Complex(0.,1.)*cplus*(-mass(qbar)/minusProduct(qbar,q)) * ( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (minusProduct(qbar,g)*plusProduct(q,qbar)/den_j))*2.*momentum(q) + (invariant(q,g)/den_j)*plusCurrent(q,g) ) );
      if ( qHel == -1 && qbarHel == -1 )
	cacheCurrent( Complex(0.,1.)*cplus*( ((sqr(mass(q))*minusProduct(qbar,g)/(minusProduct(qbar,q)*den_i)) - (plusProduct(qbar,q)*minusProduct(g,qbar)/den_j))*plusCurrent(q, qbar) - (plusProduct(g,q)*minusProduct(g,qbar)/den_j)*plusCurrent(q,g) ) ); 
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbargRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g));
#endif

    return cachedCurrent();
  }

  return czero;
}


LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralLeftCurrent(const int i, const int,
								   const int j, const int,
								   const int k, const int g1Hel,
								   const int l, const int g2Hel,
								   const int n) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kl = plusProduct(k,l);
  const Complex plusP_kn = plusProduct(k,n);
  const Complex plusP_ln = plusProduct(l,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_kn = minusProduct(k,n);
  const Complex minusP_ln = minusProduct(l,n);

  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);
  const LorentzVector<Complex> & minusC_il = minusCurrent(i,l);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);
  const LorentzVector<Complex> & minusC_kl = minusCurrent(k,l);

  const LorentzVector<Complex> & minusC_lj = minusCurrent(l,j);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2) * plusP_jl * plusP_kl * minusC_ik)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ik)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_il)/
      (kl * (jk + jl + kl)) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_ij * minusP_in)/
      (kl * (ik + il + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_il * minusP_in)/
      (ik * jl * minusP_kn) + 
      (Complex(0,2) * sqr(plusP_kl) * minusC_kj * minusP_in)/
      (kl * (ik + il + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_il * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_ij * minusP_in)/
      (kl * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * sqr(plusP_kl) * minusC_lj * minusP_in)/
      (kl * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ij * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ik * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_ij * sqr(minusP_in))/
      (ik * (ik + il + kl) * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_kj * sqr(minusP_in))/
      (ik * (ik + il + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_ij * minusP_in * minusP_jn)/
      (ik * jl * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ij * sqr(minusP_jn))/
      (jl * (jk + jl + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ik * minusP_kn)/
      (kl * (jk + jl + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_ln)/
      (jl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_ln)/
      (kl * (jk + jl + kl) * minusP_kn);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * plusP_jk * plusP_jn * minusC_ik * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_kn * minusC_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_in * minusC_ij * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_kj * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * minusC_lj * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_jn * minusC_ij * minusP_in * minusP_jl)/
      (ik * jl * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_jk * plusP_jn * minusC_ij * minusP_jl * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_ij * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_kl * plusP_kn * minusC_lj * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jk * plusP_kn * minusC_ij * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_ik * plusP_kn * minusC_ij * minusP_ik * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_lj * minusP_ik * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ij * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_il * plusP_kn * minusC_ij * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_kj * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_kl * minusC_lj * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * minusP_kn) + 
      (Complex(0,1) * plusP_jk * plusP_kn * minusC_ij * minusP_jk * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ij * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_jl * plusP_kn * minusC_ij * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_jk * plusP_jn * minusC_il * minusP_jl * minusP_ln)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ik * minusP_kl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_jl * plusP_kn * minusC_ik * minusP_kl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_jk * plusP_kn * minusC_il * minusP_kl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,2) * plusP_in * plusP_jl * minusC_il * minusP_ik)/
      (ik * jl * plusP_kn) + 
      (Complex(0,2) * plusP_jl * minusC_kl * minusP_ik)/(ik * jl) - 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_il * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_il * plusP_in * minusC_ij * minusP_ik * minusP_in)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_kj * minusP_ik * minusP_in)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_in * plusP_jl * minusC_ij * minusP_ik * minusP_jn)/
      (ik * jl * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * minusC_kj * minusP_ik * minusP_jn)/
      (ik * jl * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_jn * minusC_ij * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * plusP_ln * minusC_ij * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_kl * plusP_ln * minusC_kj * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_ij * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * plusP_jn * minusC_il * minusP_jn * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * minusC_ij * minusP_ik * minusP_kn)/
      (ik * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ij * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_ik * plusP_ln * minusC_ij * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_kl * minusC_kj * minusP_ik * minusP_kn)/
      (ik * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_kl * minusC_kj * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * minusP_ln) - 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_lj * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_il * plusP_ln * minusC_ij * minusP_il * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_kj * minusP_il * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ij * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_jk * plusP_ln * minusC_ij * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_jl * plusP_ln * minusC_ij * minusP_jl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_jl * plusP_ln * minusC_ik * minusP_kl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_il * minusP_kl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_jk * plusP_ln * minusC_il * minusP_kl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2) * sqr(plusP_in) * minusC_ij * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_kj * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_lj * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_jn * minusC_ij * minusP_ik * minusP_jl)/
      (ik * jl * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_jn * minusC_kj * minusP_ik * minusP_jl)/
      (ik * jl * plusP_ln) + 
      (Complex(0,2) * sqr(plusP_jn) * minusC_ij * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_ij * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_ij * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_kn * minusC_kj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_kn * minusC_kj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_in * minusC_ij * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) + 
      (Complex(0,2) * minusC_kj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ln * minusC_lj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_jn * minusC_ij * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_jn * minusC_ij * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * sqr(plusP_jn) * minusC_il * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_jn * minusC_ik * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_jn * minusC_il * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_ln);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedLeftCurrent(const int i, const int,
								 const int j, const int,
								 const int k, const int g1Hel,
								 const int l, const int g2Hel) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);

  const Complex plusP_kl = plusProduct(k,l);

  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_kl = minusProduct(k,l);

  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);
  const LorentzVector<Complex> & minusC_il = minusCurrent(i,l);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);
  const LorentzVector<Complex> & minusC_kl = minusCurrent(k,l);

  const LorentzVector<Complex> & minusC_lj = minusCurrent(l,j);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,-2) * plusP_jl * plusP_kl * minusC_ik)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ik)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_il)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_ij)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ij * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_ik) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_il * minusP_ij)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ij * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ik * minusP_ij)/
      (jl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ij * sqr(minusP_ij))/
      (jl * (jk + jl + kl) * minusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ik * minusP_ik)/
      (kl * (jk + jl + kl) * minusP_il) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_il)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_il * minusP_il)/
      (kl * (jk + jl + kl) * minusP_ik);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-1) * sqr(plusP_ik) * minusC_ij * minusP_il)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_lj * minusP_il)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,1) * plusP_ik * minusC_ij * sqr(minusP_il))/
      (kl * (ik + il + kl) * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_kj * sqr(minusP_il))/
      (kl * (ik + il + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_kl * minusC_lj * sqr(minusP_il))/
      (kl * (ik + il + kl) * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_jk * minusC_ij * minusP_il * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_ik * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_ij * minusP_ij * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_jl * minusC_ij * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ij * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_il * minusP_il * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ik * plusP_jk * minusC_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_jk * minusC_ij * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_jl * minusC_ik * minusP_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ik * minusP_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_jk * minusC_il * minusP_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,1) * sqr(plusP_il) * minusC_ij * minusP_ik)/
      (kl * (ik + il + kl) * plusP_ik) + 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_kj * minusP_ik)/
      (kl * (ik + il + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_jl * minusC_kl * minusP_ik)/(ik * jl) + 
      (Complex(0,2) * plusP_jl * minusC_kj * minusP_ij * minusP_ik)/
      (ik * jl * minusP_il) - 
      (Complex(0,2) * plusP_il * minusC_ij * sqr(minusP_ik))/
      (ik * (ik + il + kl) * minusP_il) - 
      (Complex(0,1) * plusP_il * minusC_ij * sqr(minusP_ik))/
      (kl * (ik + il + kl) * minusP_il) - 
      (Complex(0,2) * plusP_kl * minusC_kj * sqr(minusP_ik))/
      (ik * (ik + il + kl) * minusP_il) - 
      (Complex(0,2) * plusP_kl * minusC_kj * sqr(minusP_ik))/
      (kl * (ik + il + kl) * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_lj * sqr(minusP_ik))/
      (kl * (ik + il + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * plusP_jl * minusC_ij * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_jk * minusC_ij * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ij * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_jl * minusC_ij * minusP_ik * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_il * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_il * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_ij * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ij * plusP_jl * minusC_il * minusP_ij * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_jl * minusC_ik * minusP_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_jk * minusC_il * minusP_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_il * minusP_ik * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * plusP_ij * minusC_kj * minusP_ik * minusP_jl)/
      (ik * jl * plusP_il) + 
      (Complex(0,2) * sqr(plusP_ij) * minusC_ij * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ik * plusP_il) + 
      (Complex(0,2) * plusP_ik * minusC_kj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ik * minusC_kj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * minusC_lj * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * minusC_kj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_il * minusC_lj * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * minusC_ij * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * plusP_ij * minusC_ij * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * sqr(plusP_ij) * minusC_il * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik * plusP_il) - 
      (Complex(0,2) * plusP_ij * minusC_ik * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_ij * minusC_il * sqr(minusP_kl))/
      (kl * (jk + jl + kl) * plusP_il);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggGeneralRightCurrent(const int i, const int,
								    const int j, const int,
								    const int k, const int g1Hel,
								    const int l, const int g2Hel,
								    const int n) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kl = plusProduct(k,l);
  const Complex plusP_kn = plusProduct(k,n);
  const Complex plusP_ln = plusProduct(l,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_kn = minusProduct(k,n);
  const Complex minusP_ln = minusProduct(l,n);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);
  const LorentzVector<Complex> & minusC_jl = minusCurrent(j,l);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  const LorentzVector<Complex> & minusC_li = minusCurrent(l,i);
  const LorentzVector<Complex> & minusC_lk = minusCurrent(l,k);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jk)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_ji * minusP_in)/
      (kl * (ik + il + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_jl * minusP_in)/
      (ik * (ik + il + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ji * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * sqr(plusP_kl) * minusC_ki * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_ji * minusP_in)/
      (ik * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_ji * minusP_in)/
      (kl * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_jk * minusP_in)/
      (ik * (ik + il + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ji * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_jk * minusP_jn)/
      (ik * jl * minusP_ln) + 
      (Complex(0,2) * sqr(plusP_kl) * minusC_li * minusP_jn)/
      (kl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_il * minusC_ji * sqr(minusP_in))/
      (ik * (ik + il + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_ji * minusP_in * minusP_jn)/
      (ik * jl * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ji * sqr(minusP_jn))/
      (jl * (jk + jl + kl) * minusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_li * sqr(minusP_jn))/
      (jl * (jk + jl + kl) * minusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_kn)/
      (ik * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_kn)/
      (kl * (ik + il + kl) * minusP_ln) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jl * minusP_ln)/
      (kl * (ik + il + kl) * minusP_kn);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * plusP_ik * plusP_kn * minusC_ji * minusP_il)/
      (ik * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_ik * plusP_jn * minusC_jk * minusP_jl)/
      (ik * jl * plusP_ln) + 
      (Complex(0,2) * plusP_ik * minusC_lk * minusP_jl)/(ik * jl) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_jk * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_jk * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_ik * plusP_in * minusC_ji * minusP_il * minusP_in)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_ik * plusP_jn * minusC_ji * minusP_in * minusP_jl)/
      (ik * jl * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_ik * minusC_li * minusP_in * minusP_jl)/
      (ik * jl * minusP_kn) - 
      (Complex(0,2) * plusP_jk * plusP_jn * minusC_ji * minusP_jl * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_li * minusP_jl * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_kn * minusC_ji * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_ik * plusP_in * minusC_jk * minusP_in * minusP_kl)/
      (ik * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_jk * plusP_kn * minusC_ji * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_kl * plusP_kn * minusC_li * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_ik * plusP_kn * minusC_ji * minusP_ik * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ji * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_il * plusP_kn * minusC_ji * minusP_il * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_jk * plusP_kn * minusC_ji * minusP_jk * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_li * minusP_jk * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,2) * plusP_jk * minusC_ji * minusP_jl * minusP_ln)/
      (jl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ji * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_jl * plusP_kn * minusC_ji * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_kl * plusP_kn * minusC_ki * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_jl * minusP_ln)/
      (jl * (jk + jl + kl) * minusP_kn) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_jl * minusP_ln)/
      (kl * (jk + jl + kl) * minusP_kn) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_jk * minusP_kl * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) + 
      (Complex(0,1) * plusP_il * plusP_kn * minusC_jk * minusP_kl * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn) - 
      (Complex(0,1) * plusP_ik * plusP_kn * minusC_jl * minusP_kl * minusP_ln)/
      (kl * (ik + il + kl) * plusP_ln * minusP_kn);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,-2) * plusP_il * plusP_in * minusC_jl * minusP_ik)/
      (ik * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_il * plusP_ln * minusC_jl * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_il * plusP_in * minusC_ji * minusP_ik * minusP_in)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_in * plusP_jl * minusC_ji * minusP_ik * minusP_jn)/
      (ik * jl * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_jn * minusC_ji * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_jl * minusC_ki * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * minusP_ln) - 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_li * minusP_jk * minusP_jn)/
      (jl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * plusP_ln * minusC_ji * minusP_in * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jl * plusP_ln * minusC_ji * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_kl * plusP_ln * minusC_ki * minusP_jn * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_ji * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_ik * plusP_ln * minusC_ji * minusP_ik * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,2) * plusP_il * plusP_in * minusC_jk * minusP_ik * minusP_kn)/
      (ik * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_il * plusP_ln * minusC_ji * minusP_il * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_jn * plusP_kl * minusC_ji * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_jk * plusP_ln * minusC_ji * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_kl * minusC_ki * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * minusP_ln) + 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_li * minusP_jk * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_jl * plusP_ln * minusC_ji * minusP_jl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_kl * plusP_ln * minusC_ki * minusP_jl * minusP_kn)/
      (kl * (jk + jl + kl) * plusP_kn * minusP_ln) - 
      (Complex(0,1) * plusP_il * plusP_ln * minusC_jk * minusP_kl * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,2) * plusP_in * plusP_kl * minusC_jl * minusP_kl * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln) + 
      (Complex(0,1) * plusP_ik * plusP_ln * minusC_jl * minusP_kl * minusP_kn)/
      (kl * (ik + il + kl) * plusP_kn * minusP_ln);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2) * sqr(plusP_in) * minusC_ji * minusP_ik * minusP_il)/
      (ik * (ik + il + kl) * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_in * plusP_jn * minusC_ji * minusP_ik * minusP_jl)/
      (ik * jl * plusP_kn * plusP_ln) - 
      (Complex(0,2) * plusP_in * minusC_li * minusP_ik * minusP_jl)/
      (ik * jl * plusP_kn) + 
      (Complex(0,2) * sqr(plusP_jn) * minusC_ji * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_jn * minusC_ki * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ln) + 
      (Complex(0,2) * plusP_jn * minusC_li * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_in * minusC_ji * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ln) + 
      (Complex(0,2) * sqr(plusP_in) * minusC_jk * minusP_ik * minusP_kl)/
      (ik * (ik + il + kl) * plusP_kn * plusP_ln) + 
      (Complex(0,2) * plusP_in * minusC_ji * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_jn * minusC_ji * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * plusP_kn * minusC_ki * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ln) - 
      (Complex(0,2) * minusC_li * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_jn * minusC_ji * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_jn * minusC_ji * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_ln * minusC_li * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_ln * minusC_li * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_kn) + 
      (Complex(0,2) * plusP_in * minusC_jk * sqr(minusP_kl))/
      (kl * (ik + il + kl) * plusP_kn) - 
      (Complex(0,2) * plusP_in * minusC_jl * sqr(minusP_kl))/
      (kl * (ik + il + kl) * plusP_ln);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbarggFixedRightCurrent(const int i, const int,
								  const int j, const int,
								  const int k, const int g1Hel,
								  const int l, const int g2Hel) {
  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jl = plusProduct(j,l);

  const Complex plusP_kl = plusProduct(k,l);

  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_kl = minusProduct(k,l);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);
  const LorentzVector<Complex> & minusC_jl = minusCurrent(j,l);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  const LorentzVector<Complex> & minusC_li = minusCurrent(l,i);
  const LorentzVector<Complex> & minusC_lk = minusCurrent(l,k);

  if ( g1Hel == 1 && g2Hel == 1 ) {
    return
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jk)/
      (kl * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (ik * (ik + il + kl)) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jl)/
      (kl * (ik + il + kl)) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_ji * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * sqr(plusP_kl) * minusC_ki * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * plusP_jk * plusP_kl * minusC_ji * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_il) - 
      (Complex(0,2) * plusP_ik * plusP_jl * minusC_jk * minusP_ij)/
      (ik * jl * minusP_il) + 
      (Complex(0,2) * sqr(plusP_kl) * minusC_li * minusP_ij)/
      (kl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,2) * plusP_jk * plusP_jl * minusC_ji * sqr(minusP_ij))/
      (jl * (jk + jl + kl) * minusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_jl * plusP_kl * minusC_li * sqr(minusP_ij))/
      (jl * (jk + jl + kl) * minusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_ik)/
      (ik * (ik + il + kl) * minusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_jk * minusP_ik)/
      (kl * (ik + il + kl) * minusP_il) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_jl * minusP_il)/
      (kl * (ik + il + kl) * minusP_ik);
  }

  if ( g1Hel == 1 && g2Hel == -1 ) {
    return
      (Complex(0,-2) * sqr(plusP_ik) * minusC_ji * minusP_il)/
      (ik * (ik + il + kl) * plusP_il) - 
      (Complex(0,1) * sqr(plusP_ik) * minusC_ji * minusP_il)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,1) * plusP_ik * minusC_ji * sqr(minusP_il))/
      (kl * (ik + il + kl) * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_jk * minusC_ji * minusP_il * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_li * minusP_il * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_ik * minusC_jk * minusP_jl)/
      (ik * jl * plusP_il) + 
      (Complex(0,2) * plusP_ik * minusC_lk * minusP_jl)/(ik * jl) - 
      (Complex(0,2) * plusP_ij * plusP_jk * minusC_ji * minusP_ij * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_li * minusP_ij * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_jk * minusC_ji * minusP_il * minusP_jl)/
      (jl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,1) * plusP_ik * plusP_jl * minusC_ji * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ji * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,1) * plusP_ik * plusP_kl * minusC_ki * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_il * minusP_jl)/
      (jl * (jk + jl + kl) * minusP_ik) + 
      (Complex(0,2) * plusP_kl * minusC_li * minusP_il * minusP_jl)/
      (kl * (jk + jl + kl) * minusP_ik) - 
      (Complex(0,2) * sqr(plusP_ik) * minusC_jk * minusP_kl)/
      (ik * (ik + il + kl) * plusP_il) - 
      (Complex(0,2) * sqr(plusP_ik) * minusC_jk * minusP_kl)/
      (kl * (ik + il + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ik * plusP_jk * minusC_ji * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) - 
      (Complex(0,2) * plusP_ik * plusP_kl * minusC_li * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il * minusP_ik) + 
      (Complex(0,1) * plusP_ik * minusC_jk * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * minusP_ik) - 
      (Complex(0,1) * sqr(plusP_ik) * minusC_jl * minusP_il * minusP_kl)/
      (kl * (ik + il + kl) * plusP_il * minusP_ik);
  }

  if ( g1Hel == -1 && g2Hel == 1 ) {
    return
      (Complex(0,1) * sqr(plusP_il) * minusC_ji * minusP_ik)/
      (kl * (ik + il + kl) * plusP_ik) - 
      (Complex(0,1) * plusP_il * minusC_ji * sqr(minusP_ik))/
      (kl * (ik + il + kl) * minusP_il) - 
      (Complex(0,2) * plusP_ij * plusP_jl * minusC_ji * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * plusP_jl * minusC_ki * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * minusP_il) - 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_li * minusP_ij * minusP_jk)/
      (jl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_jk * minusC_ji * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_ij * plusP_kl * minusC_ji * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_kl * minusC_ki * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * minusP_il) + 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_li * minusP_ik * minusP_jk)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_jl * minusC_ji * minusP_ik * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * plusP_il * plusP_kl * minusC_ki * minusP_ik * minusP_jl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,2) * sqr(plusP_il) * minusC_jl * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ik) + 
      (Complex(0,2) * plusP_il * plusP_jl * minusC_ji * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) + 
      (Complex(0,2) * plusP_il * plusP_kl * minusC_ki * minusP_ij * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik * minusP_il) - 
      (Complex(0,1) * sqr(plusP_il) * minusC_jk * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * plusP_ik * minusP_il) + 
      (Complex(0,1) * plusP_il * minusC_jl * minusP_ik * minusP_kl)/
      (kl * (ik + il + kl) * minusP_il);
  }

  if ( g1Hel == -1 && g2Hel == -1 ) {
    return
      (Complex(0,2) * sqr(plusP_ij) * minusC_ji * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ik * plusP_il) + 
      (Complex(0,2) * plusP_ij * minusC_ki * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_il) + 
      (Complex(0,2) * plusP_ij * minusC_li * minusP_jk * minusP_jl)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * minusC_ji * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * plusP_ik * minusC_ki * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_il) - 
      (Complex(0,2) * minusC_li * minusP_jk * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_ij * minusC_ji * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_ij * minusC_ji * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl)) - 
      (Complex(0,2) * minusC_ki * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl)) - 
      (Complex(0,2) * plusP_il * minusC_li * minusP_jl * minusP_kl)/
      (jl * (jk + jl + kl) * plusP_ik) - 
      (Complex(0,2) * plusP_il * minusC_li * minusP_jl * minusP_kl)/
      (kl * (jk + jl + kl) * plusP_ik);
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggLeftCurrent(const int q,    const int qHel,
								   const int qbar, const int qbarHel,
								   const int g1,   const int g1Hel,
								   const int g2,   const int g2Hel) {
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<3>(1,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,q);
    LorentzVector<Complex> nj = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,qbar);
    LorentzVector<Complex> nl = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,0);
    LorentzVector<Complex> nlbar = qqbarggGeneralLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,1);
    LorentzVector<Complex> fixed = qqbarggFixedLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
    ostream& ncheck = checkStream("qqbarggLeftCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(qqbarggFixedLeftCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarggLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1)+momentum(g2));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarggRightCurrent(const int q,    const int qHel,
								    const int qbar, const int qbarHel,
								    const int g1,   const int g1Hel,
								    const int g2,   const int g2Hel) {

  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<3>(2,1,q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,q);
    LorentzVector<Complex> nj = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,qbar);
    LorentzVector<Complex> nl = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,0);
    LorentzVector<Complex> nlbar = qqbarggGeneralRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel,1);
    LorentzVector<Complex> fixed = qqbarggFixedRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
    ostream& ncheck = checkStream("qqbarggRightCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(qqbarggFixedRightCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,g2,g2Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbarggRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1)+momentum(g2));
#endif
  return cachedCurrent();


}

const LorentzVector<Complex>& MatchboxCurrents::qqbarqqbarLeftCurrent(const int q,    const int qHel,
								      const int qbar, const int qbarHel,
								      const int k,    const int kHel,
								      const int kbar, const int kbarHel) {

  if ( qHel != 1 || qbarHel != 1 ||
       abs(kHel+kbarHel) != 2 )
    return czero;

  const int i = q; const int j = qbar; const int l = kbar;

  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_kj = plusProduct(k,j);
  const Complex plusP_kl = plusProduct(k,l);
  const Complex plusP_lj = plusProduct(l,j);
  const Complex plusP_lk = plusProduct(l,k);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_il = minusProduct(i,l);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_ki = minusProduct(k,i);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_li = minusProduct(l,i);
  const Complex minusP_lk = minusProduct(l,k);

  
  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);
  const LorentzVector<Complex> & minusC_il = minusCurrent(i,l);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);

  const LorentzVector<Complex> & minusC_lj = minusCurrent(l,j);

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusP_ki * plusP_il * minusC_ij+
		     minusP_ik * plusP_lk * minusC_kj)/
		    (kl+il+ik)-
		    (minusP_jk * plusP_lj * minusC_ij+
		     minusP_lk * plusP_lj * minusC_il)/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarqqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }


  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<4>(1,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusP_li * plusP_ik * minusC_ij+
		     minusP_il * plusP_kl * minusC_lj)/
		    (kl+il+ik)-
		    (minusP_jl * plusP_kj * minusC_ij+
		     minusP_kl * plusP_kj * minusC_ik)/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarqqbarLeftCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbarqqbarRightCurrent(const int q,    const int qHel,
								       const int qbar, const int qbarHel,
								       const int k,    const int kHel,
								       const int kbar, const int kbarHel) {

  if ( qHel != -1 || qbarHel != -1 ||
       abs(kHel+kbarHel) != 2 )
    return czero;

  const int i = q; const int j = qbar; const int l = kbar;

  const double ik = invariant(i,k);
  const double il = invariant(i,l);
  const double jk = invariant(j,k);
  const double jl = invariant(j,l);
  const double kl = invariant(k,l);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_il = plusProduct(i,l);

  const Complex plusP_ki = plusProduct(k,i);
  const Complex plusP_kj = plusProduct(k,j);
  const Complex plusP_kl = plusProduct(k,l);

  const Complex plusP_li = plusProduct(l,i);
  const Complex plusP_lj = plusProduct(l,j);
  const Complex plusP_lk = plusProduct(l,k);

  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jl = minusProduct(j,l);
  const Complex minusP_ki = minusProduct(k,i);
  const Complex minusP_kl = minusProduct(k,l);
  const Complex minusP_li = minusProduct(l,i);
  const Complex minusP_lk = minusProduct(l,k);

  
  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);
  const LorentzVector<Complex> & minusC_jl = minusCurrent(j,l);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  const LorentzVector<Complex> & minusC_li = minusCurrent(l,i);

  if ( kHel == 1 && kbarHel == 1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusP_ki * plusP_il * minusC_ji+
		     minusP_lk * plusP_li * minusC_jl)/
		    (kl+il+ik)-
		    (minusP_jk * plusP_lj * minusC_ji+
		     minusP_jk * plusP_lk * minusC_ki)/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarqqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  if ( kHel == -1 && kbarHel == -1 ) {
    if ( getCurrent(hash<4>(2,1,q,qHel,qbar,qbarHel,k,kHel,kbar,kbarHel)) ) {
      cacheCurrent((Complex(0.,-2.)/kl)*
		   ((minusP_li * plusP_ik * minusC_ji+
		     minusP_kl * plusP_ki * minusC_jk)/
		    (kl+il+ik)-
		    (minusP_jl * plusP_kj * minusC_ji+
		     minusP_jl * plusP_kl * minusC_li)/
		    (kl+jl+jk)));
    }
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarqqbarRightCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(k)+momentum(kbar));
#endif
    return cachedCurrent();
  }

  return czero;

}


// Definition of sqrt to enable calculation of the sqrt of a negative double
inline Complex sqrt1 (double a) {

  if (a > 0.) { return Complex(sqrt(a), 0.) ;}
  else if (a < 0.) { return Complex(0., sqrt(abs(a))) ;}
  else { return Complex(0., 0.); }
}

// Definition of sqrt to enable calculation of the sqrt of Complex arguments
inline Complex sqrt1 (Complex a) {
  
  const double real_part = sqrt(abs(a))*cos(0.5*arg(a));
  const double imag_part = sqrt(abs(a))*sin(0.5*arg(a)); 
  return Complex(real_part, imag_part) ;
}

// Definition of log to enable continuation of the log of a negative double
inline Complex log1 (double a) {

  if (a < 0.) { return Complex(log(abs(a)), Constants::pi) ;}
  else { return Complex(log(a), 0.) ;}
}

// Definition of log to enable continuation of the log of a Complex argument with a negative real part
inline Complex log1 (Complex a) {

  return Complex(log(abs(a)), arg(a)) ;
}


const LorentzVector<Complex>& MatchboxCurrents::qqbarLeftOneLoopCurrent(const int q,    const int qHel,
									const int qbar, const int qbarHel) {

  // Note this cannot currently handle the case of one massive quark and one massless quark
  assert( (mass(q) == 0 && mass(qbar) == 0) || (mass(q) != 0 && mass(qbar) != 0) );

  // Massless quarks
  if ( mass(q) == 0 && mass(qbar) == 0 ) {

    if ( qHel != 1 || qbarHel != 1 )
      return czero;

    const LorentzVector<Complex>& tree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(1,2,q,qHel,qbar,qbarHel)) ) {
      cacheCurrent( 0.5*CF*( -8. - 3.*log1(-1./invariant(q,qbar)) - sqr(log1(-1./invariant(q,qbar))) ) * tree );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarLeftOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }

  // Massive quarks
  else {
    const LorentzVector<Complex>& momQ = momentum(q) + (sqr(mass(q))/invariant(q,qbar))*momentum(qbar);
    const LorentzVector<Complex>& momQbar = momentum(qbar) + (sqr(mass(qbar))/invariant(q,qbar))*momentum(q);

    const Complex s = (momQ+momQbar).dot(momQ+momQbar);
    const Complex inv12 = s - sqr(mass(q)) - sqr(mass(qbar));

    // Coefficient of the left-handed born-level current
    const Complex coeffLeftTree = -1.0*log1(1./sqr(mass(q)))-1.0*log1(1./sqr(mass(qbar)))-4.0 + 0.5*((2.*log1(sqr(mass(q))/sqr(mass(qbar)))*(0.5*inv12+sqr(mass(q))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))-(2.*inv12*Li2(0.5-(0.25*inv12)/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*sqr(mass(qbar)))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*Li2((0.5*(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(1.*inv12*log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*log1((2.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))*log1((-0.5*inv12-1.*sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(1.*inv12*log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(4.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(-0.25*(3.+2.*log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*sqr(inv12)+sqr(mass(q))*sqr(mass(qbar))-0.5*inv12*(1.+log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*(sqr(mass(q))+sqr(mass(qbar)))))/((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))));

    // Coefficient of the right-handed born-level current
    const Complex coeffRightTree = (2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*mass(q)*mass(qbar))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)));      
       
    const LorentzVector<Complex>& leftTree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);
    const LorentzVector<Complex>& rightTree = qqbarRightCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(1,2,q,qHel,qbar,qbarHel)) ) {

      if ( qHel == 1 && qbarHel == 1 ) {
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree) );
      }

      if ( qHel == 1 && qbarHel == -1 ) {
	// Coefficients of the left and right handed products of massive spinors
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffRightProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * mass(q)*mass(qbar)/plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == 1 ){
	// Coefficients of the left and right handed products of massive spinors
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffRightProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * mass(q)*mass(qbar)/minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == -1 ){
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree ) );
      }

    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarLeftOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }
}

const LorentzVector<Complex>& MatchboxCurrents::qqbarRightOneLoopCurrent(const int q,    const int qHel,
									 const int qbar, const int qbarHel) {

  // Note this cannot currently handle the case of one massive quark and one massless quark
  assert( (mass(q) == 0 && mass(qbar) == 0) || (mass(q) != 0 && mass(qbar) != 0) );

  // Massless quarks
  if ( mass(q) == 0 && mass(qbar) ==0 ) {

    if ( qHel != -1 || qbarHel != -1 )
      return czero;

    const LorentzVector<Complex>& tree = qqbarRightCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(2,2,q,qHel,qbar,qbarHel)) ) {
      cacheCurrent( 0.5*CF*( -8. - 3.*log1(-1./invariant(q,qbar)) - sqr(log1(-1./invariant(q,qbar))) ) * tree );
    }

#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarRightOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }

  // Massive quarks
  else {
    const LorentzVector<Complex>& momQ = momentum(q) + (sqr(mass(q))/invariant(q,qbar))*momentum(qbar);
    const LorentzVector<Complex>& momQbar = momentum(qbar) + (sqr(mass(qbar))/invariant(q,qbar))*momentum(q);

    const Complex s = (momQ+momQbar).dot(momQ+momQbar);
    const Complex inv12 = s - sqr(mass(q)) - sqr(mass(qbar));

    // Coefficient of the right-handed born-level current
    const Complex coeffRightTree = -1.0*log1(1./sqr(mass(q)))-1.0*log1(1./sqr(mass(qbar)))-4.0 + 0.5*((2.*log1(sqr(mass(q))/sqr(mass(qbar)))*(0.5*inv12+sqr(mass(q))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))-(2.*inv12*Li2(0.5-(0.25*inv12)/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*sqr(mass(qbar)))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*Li2((0.5*(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(1.*inv12*log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(2.*inv12*log1((2.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))*log1((-0.5*inv12-1.*sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(1.*inv12*log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1(-((0.5*inv12+sqr(mass(q))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(0.5*inv12*sqr(log1((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(qbar))-1.*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1(-((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))/(0.5*inv12+sqr(mass(q))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))-(0.5*inv12*sqr(log1((0.5*inv12+sqr(mass(qbar))+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar))))/(1.*inv12+sqr(mass(q))+sqr(mass(qbar))))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))+(4.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(-0.25*(3.+2.*log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*sqr(inv12)+sqr(mass(q))*sqr(mass(qbar))-0.5*inv12*(1.+log1(1./(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))))*(sqr(mass(q))+sqr(mass(qbar)))))/((1.*inv12+sqr(mass(q))+sqr(mass(qbar)))*sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))));

    // Coefficient of the left-handed born-level current
    const Complex coeffLeftTree = (2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*mass(q)*mass(qbar))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)));      
       
    const LorentzVector<Complex>& leftTree = qqbarLeftCurrent(q,qHel,qbar,qbarHel);
    const LorentzVector<Complex>& rightTree = qqbarRightCurrent(q,qHel,qbar,qbarHel);

    if ( getCurrent(hash<1>(2,2,q,qHel,qbar,qbarHel)) ) {

      if ( qHel == 1 && qbarHel == 1 ) {
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree ) );
      }

      if ( qHel == 1 && qbarHel == -1 ) {
	// Coefficients of the right and left handed products of massive spinors
	const LorentzVector<Complex>& coeffRightProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * mass(q)*mass(qbar)/plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == 1 ){
	// Coefficients of the right and left handed products of massive spinors
	const LorentzVector<Complex>& coeffRightProd = ( (mass(qbar)*(-2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(2.*momQ+momQbar)+(3.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))-(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(2.*momQ+momQbar)*sqr(inv12)-momQbar*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*((5.*momQ+2.*momQbar)*sqr(mass(q))+momQ*sqr(mass(qbar)))+(2.*momQ+momQbar)*sqr(sqr(mass(q)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );
	const LorentzVector<Complex>& coeffLeftProd = ( (mass(q)*(2.*(momQ+momQbar)*(1.*inv12+sqr(mass(q))+sqr(mass(qbar)))+log1(sqr(mass(q))/sqr(mass(qbar)))*(1.*inv12*(momQ+2.*momQbar)+momQbar*sqr(mass(q))+(2.*momQ+3.*momQbar)*sqr(mass(qbar)))+(2.*log1((-1.*mass(q)*mass(qbar))/(0.5*inv12+sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))*(0.5*(momQ+2.*momQbar)*sqr(inv12)-momQ*sqr(mass(q))*sqr(mass(qbar))+0.5*inv12*(momQbar*sqr(mass(q))+(2.*momQ+5.*momQbar)*sqr(mass(qbar)))+(momQ+2.*momQbar)*sqr(sqr(mass(qbar)))))/sqrt1(0.25*sqr(inv12)-sqr(mass(q))*sqr(mass(qbar)))))/sqr(1.*inv12+sqr(mass(q))+sqr(mass(qbar))) );

	const Complex leftProd = Complex(0.,1.) * mass(q)*mass(qbar)/minusProduct(q,qbar);
	const Complex rightProd = Complex(0.,1.) * plusProduct(q,qbar);

	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree + coeffLeftProd*leftProd + coeffRightProd*rightProd ) );
      }

      if ( qHel == -1 && qbarHel == -1 ){
	cacheCurrent( 0.5*CF*( coeffLeftTree*leftTree + coeffRightTree*rightTree ) );
      }

    }
  
#ifdef CHECK_MatchboxCurrents
    checkCurrent("qqbarRightOneLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar));
#endif

    return cachedCurrent();
  }
}


// ln(s(a+i0))
inline Complex log(double s, double a) {
  return
    s < 0. ?
	Complex(log(abs(a)),-pi * theta(a)) :
    Complex(log(abs(a)),pi * theta(-a));
}

// ln(s(a+i0)/(b+i0))
inline Complex log(double s, double a, double b) {
  return
    s < 0. ?
	Complex(log(abs(a/b)),-pi * theta(a/b) * sign(b-a)) :
    Complex(log(abs(a/b)),pi * theta(-a/b) * sign(b-a));
}


// Li2(-(a+i0)/(b+i0))
inline Complex Li2(double a, double b) {
  if ( -a/b < 1. )
    return Complex(Herwig::Math::ReLi2(-a/b),0.0);
  return Complex(Herwig::Math::ReLi2(-a/b),-pi * log(-a/b) * sign(b-a));
}

Complex MatchboxCurrents::box6(const int i, const int j, const int k) {

  const double sij = invariant(i,j);
  const double sik = invariant(i,k);
  const double sjk = invariant(j,k);

  return
    -( Li2(sik+sjk,sij) + Li2(sik+sij,sjk) + 0.5 * csqr(log(1.,sij,sjk)) + sqr(pi)/6. )/8.;
      
}

void MatchboxCurrents::qqbargLoopCoefficients(const int i, const int j, const int k) {

  // use a dummy cache entry to check if we need to get some work done
  static Complex dummy;

  if ( getAmplitude(hash<5>(1,2,i,0,j,0,k,0)) ) {
    dummy = 0.;
    cacheAmplitude(dummy);
    cachedAmplitude();
  } else {
    cachedAmplitude();
    return;
  }

  qqbargLoops.resize(13);

  // get the transcendentals
  const double ij = invariant(i,j);  
  const double ij2 = sqr(ij);
  const double ij3 = ij2 * ij;

  const double ik = invariant(i,k);
  const double ik2 = sqr(ik);
  //const double ik3 = ik2 * ik;

  const double jk = invariant(j,k);
  const double jk2 = sqr(jk);
  const double jk3 = jk2 * jk;

  const double ij_ik = ij + ik;
  const double ij_ik_2 = sqr(ij_ik);

  const double ij_jk = ij + jk;
  const double ij_jk_2 = sqr(ij_jk);

  const double ik_jk = ik + jk;
  const double ik_jk_2 = sqr(ik_jk);


  const double Q2 = ij + ik + jk;
  // checked for LEP that virtuals + I operator are mu2 independent
  //double xmu2 = 10 * GeV2/sqr(amplitudeScale());
  const double xmu2 = 1.;

  const Complex Lijk = log(1.,-xmu2/Q2);

  const Complex Lij = log(1.,Q2,ij);
  const Complex Lik = log(1.,Q2,ik);
  const Complex Ljk = log(1.,Q2,jk);

  const Complex Box6ijk = box6(i,j,k);
  const Complex Box6ikj = box6(i,k,j);
  const Complex Box6jik = box6(j,i,k);

  // get the coefficients
  qqbargLoops[0] = (
		    (2 * CF * ij2) 			- 
		    (32 * CA * Box6ijk * ij2)		+ 
		    (64 * CF * Box6ijk * ij2)		- 
		    (8 * CA * Box6jik * ij2)		+ 
		    (16 * CF * Box6jik * ij2)		+ 
		    (2 * CA * Lij * ij2)		- 
		    (4 * CF * Lij * ij2)		- 
		    (CA * Lik * ij2)			- 
		    (2 * CF * Lik * ij2)		- 
		    (4 * CF * Ljk * ij2)		- 
		    (16 * CA * Box6ijk * ij3)		/ ik + 
		    (32 * CF * Box6ijk * ij3)		/ ik + 
		    (CA * Lij * ij3)			/ ik - 
		    (2 * CF * Lij * ij3)		/ ik + 
		    (2 * CF * ij * ik)			- 
		    (16 * CA * Box6ijk * ij * ik)	+ 
		    (32 * CF * Box6ijk * ij * ik)	- 
		    (16 * CA * Box6jik * ij * ik)	+ 
		    (32 * CF * Box6jik * ij * ik)	+ 
		    (CA * Lij * ij * ik)		- 
		    (2 * CF * Lij * ij * ik)		- 
		    (2 * CA * Lik * ij * ik)		- 
		    (4 * CF * Lik * ij * ik)		- 
		    (4 * CF * Ljk * ij * ik)		- 
		    (8 * CA * Box6jik * ik2)		+ 
		    (16 * CF * Box6jik * ik2)		- 
		    (CA * Lik * ik2)			- 
		    (2 * CF * Lik * ik2)		- 
		    (8 * CA * Box6jik * ij3)		/ jk + 
		    (16 * CF * Box6jik * ij3)		/ jk - 
		    (16 * CA * Box6jik * ij2 * ik)	/ jk + 
		    (32 * CF * Box6jik * ij2 * ik)	/ jk - 
		    (8 * CA * Box6jik * ij * ik2)	/ jk + 
		    (16 * CF * Box6jik * ij * ik2)	/ jk + 
		    (2 * CF * ij * jk)			- 
		    (40 * CA * Box6ijk * ij * jk)	+ 
		    (80 * CF * Box6ijk * ij * jk)	+ 
		    (24 * CA * Box6ikj * ij * jk)	+ 
		    (2 * CA * Lij * ij * jk)		- 
		    (4 * CF * Lij * ij * jk)		- 
		    (CA * Lik * ij * jk)		- 
		    (4 * CF * Lik * ij * jk)		- 
		    (12 * CF * Ljk * ij * jk)		- 
		    (8 * CA * Box6ijk * ij3 * jk)	/ ik2 + 
		    (16 * CF * Box6ijk * ij3 * jk)	/ ik2 - 
		    (32 * CA * Box6ijk * ij2 * jk)	/ ik + 
		    (64 * CF * Box6ijk * ij2 * jk)	/ ik + 
		    (CA * Lij * ij2 * jk)		/ ik - 
		    (2 * CF * Lij * ij2 * jk)		/ ik + 
		    (CA * Ljk * ij2 * jk)		/ ik - 
		    (2 * CF * Ljk * ij2 * jk)		/ ik + 
		    (2 * CF * ik * jk)			- 
		    (16 * CA * Box6ijk * ik * jk)	+ 
		    (32 * CF * Box6ijk * ik * jk)	+ 
		    (48 * CA * Box6ikj * ik * jk)	+ 
		    (CA * Lij * ik * jk)		- 
		    (2 * CF * Lij * ik * jk)		- 
		    (2 * CA * Lik * ik * jk)		- 
		    (8 * CF * Lik * ik * jk)		- 
		    (CA * Ljk * ik * jk)		- 
		    (8 * CF * Ljk * ik * jk)		+ 
		    (24 * CA * Box6ikj * ik2 * jk)	/ ij - 
		    (CA * Lik * ik2 * jk)		/ ij - 
		    (4 * CF * Lik * ik2 * jk)		/ ij - 
		    (8 * CA * Box6ijk * jk2)		+ 
		    (16 * CF * Box6ijk * jk2)		+ 
		    (24 * CA * Box6ikj * jk2)		- 
		    (8 * CF * Ljk * jk2)		- 
		    (8 * CA * Box6ijk * ij2 * jk2)	/ ik2 + 
		    (16 * CF * Box6ijk * ij2 * jk2)	/ ik2 - 
		    (16 * CA * Box6ijk * ij * jk2)	/ ik + 
		    (32 * CF * Box6ijk * ij * jk2)	/ ik + 
		    (CA * Ljk * ij * jk2)		/ ik - 
		    (2 * CF * Ljk * ij * jk2)		/ ik + 
		    (48 * CA * Box6ikj * ik * jk2)	/ ij - 
		    (CA * Ljk * ik * jk2)		/ ij - 
		    (4 * CF * Ljk * ik * jk2)		/ ij + 
		    (24 * CA * Box6ikj * ik2 * jk2)	/ ij2
		    ) / (ij_ik_2 * ij_jk);

  qqbargLoops[1] = (
		    (-2 * CF * ij2)			+ 
		    (8 * CA * Box6ijk * ij2)		- 
		    (16 * CF * Box6ijk * ij2)		+ 
		    (32 * CA * Box6jik * ij2)		- 
		    (64 * CF * Box6jik * ij2)		- 
		    (2 * CA * Lij * ij2)		+ 
		    (4 * CF * Lij * ij2)		+ 
		    (4 * CF * Lik * ij2)		+ 
		    (CA * Ljk * ij2)			+ 
		    (2 * CF * Ljk * ij2)		+ 
		    (8 * CA * Box6ijk * ij3)		/ ik - 
		    (16 * CF * Box6ijk * ij3)		/ ik - 
		    (2 * CF * ij * ik)			- 
		    (24 * CA * Box6ikj * ij * ik)	+ 
		    (40 * CA * Box6jik * ij * ik)	- 
		    (80 * CF * Box6jik * ij * ik)	- 
		    (2 * CA * Lij * ij * ik)		+ 
		    (4 * CF * Lij * ij * ik)		+ 
		    (12 * CF * Lik * ij * ik)		+ 
		    (CA * Ljk * ij * ik)		+ 
		    (4 * CF * Ljk * ij * ik)		- 
		    (24 * CA * Box6ikj * ik2)		+ 
		    (8 * CA * Box6jik * ik2)		- 
		    (16 * CF * Box6jik * ik2)		+ 
		    (8 * CF * Lik * ik2)		+ 
		    (8 * CA * Box6jik * ij3 * ik)	/ jk2 - 
		    (16 * CF * Box6jik * ij3 * ik)	/ jk2 + 
		    (8 * CA * Box6jik * ij2 * ik2)	/ jk2 - 
		    (16 * CF * Box6jik * ij2 * ik2)	/ jk2 + 
		    (16 * CA * Box6jik * ij3)		/ jk - 
		    (32 * CF * Box6jik * ij3)		/ jk - 
		    (CA * Lij * ij3)			/ jk + 
		    (2 * CF * Lij * ij3)		/ jk + 
		    (32 * CA * Box6jik * ij2 * ik)	/ jk - 
		    (64 * CF * Box6jik * ij2 * ik)	/ jk - 
		    (CA * Lij * ij2 * ik)		/ jk + 
		    (2 * CF * Lij * ij2 * ik)		/ jk - 
		    (CA * Lik * ij2 * ik)		/ jk + 
		    (2 * CF * Lik * ij2 * ik)		/ jk + 
		    (16 * CA * Box6jik * ij * ik2)	/ jk - 
		    (32 * CF * Box6jik * ij * ik2)	/ jk - 
		    (CA * Lik * ij * ik2)		/ jk + 
		    (2 * CF * Lik * ij * ik2)		/ jk - 
		    (2 * CF * ij * jk)			+ 
		    (16 * CA * Box6ijk * ij * jk)	- 
		    (32 * CF * Box6ijk * ij * jk)	+ 
		    (16 * CA * Box6jik * ij * jk)	- 
		    (32 * CF * Box6jik * ij * jk)	- 
		    (CA * Lij * ij * jk)		+ 
		    (2 * CF * Lij * ij * jk)		+ 
		    (4 * CF * Lik * ij * jk)		+ 
		    (2 * CA * Ljk * ij * jk)		+ 
		    (4 * CF * Ljk * ij * jk)		+ 
		    (16 * CA * Box6ijk * ij2 * jk)	/ ik - 
		    (32 * CF * Box6ijk * ij2 * jk)	/ ik - 
		    (2 * CF * ik * jk)			- 
		    (48 * CA * Box6ikj * ik * jk)	+ 
		    (16 * CA * Box6jik * ik * jk)	- 
		    (32 * CF * Box6jik * ik * jk)	- 
		    (CA * Lij * ik * jk)		+ 
		    (2 * CF * Lij * ik * jk)		+ 
		    (CA * Lik * ik * jk)		+ 
		    (8 * CF * Lik * ik * jk)		+ 
		    (2 * CA * Ljk * ik * jk)		+ 
		    (8 * CF * Ljk * ik * jk)		- 
		    (48 * CA * Box6ikj * ik2 * jk)	/ ij + 
		    (CA * Lik * ik2 * jk)		/ ij + 
		    (4 * CF * Lik * ik2 * jk)		/ ij + 
		    (8 * CA * Box6ijk * jk2)		- 
		    (16 * CF * Box6ijk * jk2)		+ 
		    (CA * Ljk * jk2)			+ 
		    (2 * CF * Ljk * jk2)		+ 
		    (8 * CA * Box6ijk * ij * jk2)	/ ik - 
		    (16 * CF * Box6ijk * ij * jk2)	/ ik - 
		    (24 * CA * Box6ikj * ik * jk2)	/ ij + 
		    (CA * Ljk * ik * jk2)		/ ij + 
		    (4 * CF * Ljk * ik * jk2)		/ ij - 
		    (24 * CA * Box6ikj * ik2 * jk2)	/ ij2
		    ) / (ij_ik * ij_jk_2);

  qqbargLoops[2] = -3 * CF * Lijk + (
				     (-4 * CA * Box6jik * ij3)		+ 
				     (8 * CF * Box6jik * ij3)		+ 
				     (CA * Lij * ij3)			/ 2. -
				     (CF * Lij * ij3)			+ 
				     (CA * ij2 * ik)			- 
				     (9 * CF * ij2 * ik)			+ 
				     (8 * CA * Box6ijk * ij2 * ik)	- 
				     (16 * CF * Box6ijk * ij2 * ik)	- 
				     (8 * CA * Box6ikj * ij2 * ik)	- 
				     (8 * CA * Box6jik * ij2 * ik)	+ 
				     (16 * CF * Box6jik * ij2 * ik)	+ 
				     (CA * Lij * ij2 * ik)		/ 2. -
				     (CF * Lij * ij2 * ik)		+ 
				     (CA * Lik * ij2 * ik)		/ 2. -
				     (CF * Lik * ij2 * ik)		+ 
				     (CA * ij * ik2)			- 
				     (9 * CF * ij * ik2)			+ 
				     (8 * CA * Box6ijk * ij * ik2)	- 
				     (16 * CF * Box6ijk * ij * ik2)	- 
				     (8 * CA * Box6ikj * ij * ik2)	- 
				     (4 * CA * Box6jik * ij * ik2)	+ 
				     (8 * CF * Box6jik * ij * ik2)	+ 
				     (CA * Lik * ij * ik2)		/ 2. -
				     (CF * Lik * ij * ik2)		- 
				     (4 * CA * Box6jik * ij3 * ik)	/ jk + 
				     (8 * CF * Box6jik * ij3 * ik)	/ jk - 
				     (4 * CA * Box6jik * ij2 * ik2)	/ jk + 
				     (8 * CF * Box6jik * ij2 * ik2)	/ jk + 
				     (CA * ij2 * jk)			- 
				     (9 * CF * ij2 * jk)			+ 
				     (12 * CA * Box6ijk * ij2 * jk)	- 
				     (24 * CF * Box6ijk * ij2 * jk)	- 
				     (8 * CA * Box6ikj * ij2 * jk)	- 
				     (4 * CA * Box6jik * ij2 * jk)	+ 
				     (8 * CF * Box6jik * ij2 * jk)	+ 
				     (CA * Lik * ij2 * jk)		/ 2. -
				     (CF * Lik * ij2 * jk)		- 
				     (CA * Ljk * ij2 * jk)		/ 2. +
				     (CF * Ljk * ij2 * jk)		+ 
				     (4 * CA * Box6ijk * ij3 * jk)	/ ik - 
				     (8 * CF * Box6ijk * ij3 * jk)	/ ik - 
				     (CA * Lij * ij3 * jk)		/ (2. * ik) + 
				     (CF * Lij * ij3 * jk)		/ ik + 
				     (2 * CA * ij * ik * jk)		- 
				     (18 * CF * ij * ik * jk)		+ 
				     (16 * CA * Box6ijk * ij * ik * jk)	- 
				     (32 * CF * Box6ijk * ij * ik * jk)	- 
				     (28 * CA * Box6ikj * ij * ik * jk)	- 
				     (4 * CA * Box6jik * ij * ik * jk)	+ 
				     (8 * CF * Box6jik * ij * ik * jk)	+ 
				     (CA * Lij * ij * ik * jk)		/ 2. -
				     (CF * Lij * ij * ik * jk)		+ 
				     (CA * Lik * ij * ik * jk)		- 
				     (CF * Lik * ij * ik * jk)		- 
				     (CA * Ljk * ij * ik * jk)		/ 2. +
				     (3 * CF * Ljk * ij * ik * jk)	+ 
				     (CA * ik2 * jk)			- 
				     (9 * CF * ik2 * jk)			+ 
				     (8 * CA * Box6ijk * ik2 * jk)	- 
				     (16 * CF * Box6ijk * ik2 * jk)	- 
				     (20 * CA * Box6ikj * ik2 * jk)	+ 
				     (CA * Lik * ik2 * jk)		/ 2. +
				     (CA * ij * jk2)			- 
				     (9 * CF * ij * jk2)			+ 
				     (12 * CA * Box6ijk * ij * jk2)	- 
				     (24 * CF * Box6ijk * ij * jk2)	- 
				     (20 * CA * Box6ikj * ij * jk2)	- 
				     (CA * Lij * ij * jk2)		/ 2. + 
				     (CF * Lij * ij * jk2)		+ 
				     (CA * Lik * ij * jk2)		/ 2. - 
				     (CA * Ljk * ij * jk2)		+ 
				     (4 * CF * Ljk * ij * jk2)		+ 
				     (4 * CA * Box6ijk * ij3 * jk2)	/ ik2 - 
				     (8 * CF * Box6ijk * ij3 * jk2)	/ ik2 + 
				     (8 * CA * Box6ijk * ij2 * jk2)	/ ik - 
				     (16 * CF * Box6ijk * ij2 * jk2)	/ ik - 
				     (CA * Lij * ij2 * jk2)		/ (2. * ik) + 
				     (CF * Lij * ij2 * jk2)		/ ik - 
				     (CA * Ljk * ij2 * jk2)		/ (2. * ik) + 
				     (CF * Ljk * ij2 * jk2)		/ ik + 
				     (CA * ik * jk2)			- 
				     (9 * CF * ik * jk2)			+ 
				     (8 * CA * Box6ijk * ik * jk2)	- 
				     (16 * CF * Box6ijk * ik * jk2)	- 
				     (32 * CA * Box6ikj * ik * jk2)	+ 
				     (CA * Lik * ik * jk2)		/ 2. - 
				     (CA * Ljk * ik * jk2)		/ 2. + 
				     (3 * CF * Ljk * ik * jk2)		- 
				     (12 * CA * Box6ikj * ik2 * jk2)	/ ij - 
				     (12 * CA * Box6ikj * jk3)		- 
				     (CA * Ljk * jk3)			/ 2. +
				     (3 * CF * Ljk * jk3)		+ 
				     (4 * CA * Box6ijk * ij2 * jk3)	/ ik2 - 
				     (8 * CF * Box6ijk * ij2 * jk3)	/ ik2 + 
				     (4 * CA * Box6ijk * ij * jk3)	/ ik - 
				     (8 * CF * Box6ijk * ij * jk3)	/ ik - 
				     (CA * Ljk * ij * jk3)		/ (2. * ik) + 
				     (CF * Ljk * ij * jk3)		/ ik - 
				     (12 * CA * Box6ikj * ik * jk3)	/ ij
				     ) / (ij_ik * ij_jk * ik_jk);

  qqbargLoops[3] = 3 * CF * Lijk + (
				    (8 * CF * ij2)			- 
				    (8 * CA * Box6ijk * ij2)		+ 
				    (16 * CF * Box6ijk * ij2)		+ 
				    (8 * CA * Box6ikj * ij2)		- 
				    (8 * CA * Box6jik * ij2)		+ 
				    (16 * CF * Box6jik * ij2)		+ 
				    (CA * Lij * ij2)			/ 2. - 
				    (CF * Lij * ij2)			+ 
				    (8 * CF * ij * ik)			- 
				    (8 * CA * Box6ijk * ij * ik)	+ 
				    (16 * CF * Box6ijk * ij * ik)	+ 
				    (8 * CA * Box6ikj * ij * ik)	- 
				    (12 * CA * Box6jik * ij * ik)	+ 
				    (24 * CF * Box6jik * ij * ik)	+ 
				    (CA * Lij * ij * ik)		/ 2. - 
				    (CF * Lij * ij * ik)		+ 
				    (CA * Lik * ij * ik)		/ 2. - 
				    (CF * Lik * ij * ik)		- 
				    (4 * CA * Box6jik * ik2)		+ 
				    (8 * CF * Box6jik * ik2)		+ 
				    (CA * Lik * ik2)			/ 2. - 
				    (CF * Lik * ik2)			- 
				    (4 * CA * Box6jik * ij2 * ik)	/ jk + 
				    (8 * CF * Box6jik * ij2 * ik)	/ jk - 
				    (4 * CA * Box6jik * ij * ik2)	/ jk + 
				    (8 * CF * Box6jik * ij * ik2)	/ jk + 
				    (8 * CF * ij * jk)			- 
				    (12 * CA * Box6ijk * ij * jk)	+ 
				    (24 * CF * Box6ijk * ij * jk)	+ 
				    (8 * CA * Box6ikj * ij * jk)	- 
				    (8 * CA * Box6jik * ij * jk)	+ 
				    (16 * CF * Box6jik * ij * jk)	+ 
				    (CA * Lij * ij * jk)		/ 2. - 
				    (CF * Lij * ij * jk)		+ 
				    (CA * Ljk * ij * jk)		/ 2. - 
				    (CF * Ljk * ij * jk)		- 
				    (4 * CA * Box6ijk * ij2 * jk)	/ ik + 
				    (8 * CF * Box6ijk * ij2 * jk)	/ ik + 
				    (8 * CF * ik * jk)			- 
				    (8 * CA * Box6ijk * ik * jk)	+ 
				    (16 * CF * Box6ijk * ik * jk)	- 
				    (4 * CA * Box6ikj * ik * jk)	- 
				    (8 * CA * Box6jik * ik * jk)	+ 
				    (16 * CF * Box6jik * ik * jk)	+ 
				    (CA * Lij * ik * jk)		/ 2. - 
				    (CF * Lij * ik * jk)		+ 
				    (CA * Lik * ik * jk)		/ 2. + 
				    (2 * CF * Lik * ik * jk)		+ 
				    (CA * Ljk * ik * jk)		/ 2. + 
				    (2 * CF * Ljk * ik * jk)		- 
				    (12 * CA * Box6ikj * ik2 * jk)	/ ij + 
				    (CA * Lik * ik2 * jk)		/ (2. * ij) + 
				    (2 * CF * Lik * ik2 * jk)		/ ij - 
				    (4 * CA * Box6ijk * jk2)		+ 
				    (8 * CF * Box6ijk * jk2)		+ 
				    (CA * Ljk * jk2)			/ 2. - 
				    (CF * Ljk * jk2)			- 
				    (4 * CA * Box6ijk * ij * jk2)	/ ik + 
				    (8 * CF * Box6ijk * ij * jk2)	/ ik - 
				    (12 * CA * Box6ikj * ik * jk2)	/ ij + 
				    (CA * Ljk * ik * jk2)		/ (2. * ij) + 
				    (2 * CF * Ljk * ik * jk2)		/ ij - 
				    (12 * CA * Box6ikj * ik2 * jk2)	/ ij2
				    ) / (ij_ik * ij_jk);

  qqbargLoops[4] = -3 * CF * Lijk + (
				     (-8 * CF * ij2)			+ 
				     (8 * CA * Box6ijk * ij2)		- 
				     (16 * CF * Box6ijk * ij2)		- 
				     (8 * CA * Box6ikj * ij2)		+ 
				     (8 * CA * Box6jik * ij2)		- 
				     (16 * CF * Box6jik * ij2)		- 
				     (CA * Lij * ij2)			/ 2. + 
				     (CF * Lij * ij2)			- 
				     (8 * CF * ij * ik)			+ 
				     (8 * CA * Box6ijk * ij * ik)	- 
				     (16 * CF * Box6ijk * ij * ik)	- 
				     (8 * CA * Box6ikj * ij * ik)	+ 
				     (12 * CA * Box6jik * ij * ik)	- 
				     (24 * CF * Box6jik * ij * ik)	- 
				     (CA * Lij * ij * ik)		/ 2. + 
				     (CF * Lij * ij * ik)		- 
				     (CA * Lik * ij * ik)		/ 2. + 
				     (CF * Lik * ij * ik)		+ 
				     (4 * CA * Box6jik * ik2)		- 
				     (8 * CF * Box6jik * ik2)		- 
				     (CA * Lik * ik2)			/ 2. + 
				     (CF * Lik * ik2)			+ 
				     (4 * CA * Box6jik * ij2 * ik)	/ jk - 
				     (8 * CF * Box6jik * ij2 * ik)	/ jk + 
				     (4 * CA * Box6jik * ij * ik2)	/ jk - 
				     (8 * CF * Box6jik * ij * ik2)	/ jk - 
				     (8 * CF * ij * jk)			+ 
				     (12 * CA * Box6ijk * ij * jk)	- 
				     (24 * CF * Box6ijk * ij * jk)	- 
				     (8 * CA * Box6ikj * ij * jk)	+ 
				     (8 * CA * Box6jik * ij * jk)	- 
				     (16 * CF * Box6jik * ij * jk)	- 
				     (CA * Lij * ij * jk)		/ 2. + 
				     (CF * Lij * ij * jk)		- 
				     (CA * Ljk * ij * jk)		/ 2. + 
				     (CF * Ljk * ij * jk)		+ 
				     (4 * CA * Box6ijk * ij2 * jk)	/ ik - 
				     (8 * CF * Box6ijk * ij2 * jk)	/ ik - 
				     (8 * CF * ik * jk)			+ 
				     (8 * CA * Box6ijk * ik * jk)	- 
				     (16 * CF * Box6ijk * ik * jk)	+ 
				     (4 * CA * Box6ikj * ik * jk)	+ 
				     (8 * CA * Box6jik * ik * jk)	- 
				     (16 * CF * Box6jik * ik * jk)	- 
				     (CA * Lij * ik * jk)		/ 2. + 
				     (CF * Lij * ik * jk)		- 
				     (CA * Lik * ik * jk)		/ 2. - 
				     (2 * CF * Lik * ik * jk)		- 
				     (CA * Ljk * ik * jk)		/ 2. - 
				     (2 * CF * Ljk * ik * jk)		+ 
				     (12 * CA * Box6ikj * ik2 * jk)	/ ij - 
				     (CA * Lik * ik2 * jk)		/ (2. * ij) - 
				     (2 * CF * Lik * ik2 * jk)		/ ij + 
				     (4 * CA * Box6ijk * jk2)		- 
				     (8 * CF * Box6ijk * jk2)		- 
				     (CA * Ljk * jk2)			/ 2. + 
				     (CF * Ljk * jk2)			+ 
				     (4 * CA * Box6ijk * ij * jk2)	/ ik - 
				     (8 * CF * Box6ijk * ij * jk2)	/ ik + 
				     (12 * CA * Box6ikj * ik * jk2)	/ ij - 
				     (CA * Ljk * ik * jk2)		/ (2. * ij) - 
				     (2 * CF * Ljk * ik * jk2)		/ ij + 
				     (12 * CA * Box6ikj * ik2 * jk2)	/ ij2
				     ) / (ij_ik * ij_jk);

  qqbargLoops[5] = 3 * CF * Lijk + (
				    (-4 * CA * Box6jik * ij2)		+ 
				    (8 * CF * Box6jik * ij2)		+ 
				    (CA * Lij * ij2)			/ 2. - 
				    (CF * Lij * ij2)			- 
				    (CA * ij * ik)			+ 
				    (9 * CF * ij * ik)			- 
				    (8 * CA * Box6ijk * ij * ik)	+ 
				    (16 * CF * Box6ijk * ij * ik)	+ 
				    (8 * CA * Box6ikj * ij * ik)	- 
				    (4 * CA * Box6jik * ij * ik)	+ 
				    (8 * CF * Box6jik * ij * ik)	+ 
				    (CA * Lij * ij * ik)		/ 2. - 
				    (CF * Lij * ij * ik)		+ 
				    (CA * Lik * ij * ik)		/ 2. - 
				    (CF * Lik * ij * ik)		- 
				    (CA * ik2)				+ 
				    (9 * CF * ik2)			- 
				    (8 * CA * Box6ijk * ik2)		+ 
				    (16 * CF * Box6ijk * ik2)		+ 
				    (8 * CA * Box6ikj * ik2)		+ 
				    (CA * Lik * ik2)			/ 2. - 
				    (CF * Lik * ik2)			- 
				    (4 * CA * Box6jik * ij2 * ik)	/ jk + 
				    (8 * CF * Box6jik * ij2 * ik)	/ jk - 
				    (4 * CA * Box6jik * ij * ik2)	/ jk + 
				    (8 * CF * Box6jik * ij * ik2)	/ jk - 
				    (CA * ij * jk)			+ 
				    (9 * CF * ij * jk)			- 
				    (4 * CA * Box6ijk * ij * jk)	+ 
				    (8 * CF * Box6ijk * ij * jk)	+ 
				    (8 * CA * Box6ikj * ij * jk)	- 
				    (CA * Lij * ij * jk)		/ 2. + 
				    (CF * Lij * ij * jk)		+ 
				    (CA * Lik * ij * jk)		/ 2. - 
				    (CF * Lik * ij * jk)		- 
				    (CA * Ljk * ij * jk)		/ 2. + 
				    (CF * Ljk * ij * jk)		+ 
				    (4 * CA * Box6ijk * ij2 * jk)	/ ik - 
				    (8 * CF * Box6ijk * ij2 * jk)	/ ik - 
				    (CA * Lij * ij2 * jk)		/ (2. * ik) + 
				    (CF * Lij * ij2 * jk)		/ ik - 
				    (CA * ik * jk)			+ 
				    (9 * CF * ik * jk)			- 
				    (8 * CA * Box6ijk * ik * jk)	+ 
				    (16 * CF * Box6ijk * ik * jk)	+ 
				    (20 * CA * Box6ikj * ik * jk)	+ 
				    (CA * Lik * ik * jk)		/ 2. - 
				    (CF * Lik * ik * jk)		- 
				    (CA * Ljk * ik * jk)		/ 2. - 
				    (2 * CF * Ljk * ik * jk)		+ 
				    (12 * CA * Box6ikj * ik2 * jk)	/ ij + 
				    (12 * CA * Box6ikj * jk2)		- 
				    (CA * Ljk * jk2)			/ 2. - 
				    (2 * CF * Ljk * jk2)		+ 
				    (4 * CA * Box6ijk * ij2 * jk2)	/ ik2 - 
				    (8 * CF * Box6ijk * ij2 * jk2)	/ ik2 + 
				    (4 * CA * Box6ijk * ij * jk2)	/ ik - 
				    (8 * CF * Box6ijk * ij * jk2)	/ ik - 
				    (CA * Ljk * ij * jk2)		/ (2. * ik) + 
				    (CF * Ljk * ij * jk2)		/ ik + 
				    (12 * CA * Box6ikj * ik * jk2)	/ ij
				    ) / (ij_ik * ik_jk);

  qqbargLoops[6] = (
		    (-2 * CF * ij)			+ 
		    (32 * CA * Box6ijk * ij)		- 
		    (64 * CF * Box6ijk * ij)		- 
		    (4 * CA * Lij * ij)			+ 
		    (8 * CF * Lij * ij)			+ 
		    (4 * CF * Ljk * ij)			+ 
		    (16 * CA * Box6ijk * ij2)		/ ik - 
		    (32 * CF * Box6ijk * ij2)		/ ik - 
		    (2 * CA * Lij * ij2)		/ ik + 
		    (4 * CF * Lij * ij2)		/ ik - 
		    (2 * CF * ik)			+ 
		    (16 * CA * Box6ijk * ik)		- 
		    (32 * CF * Box6ijk * ik)		- 
		    (2 * CA * Lij * ik)			+ 
		    (4 * CF * Lij * ik)			+ 
		    (4 * CF * Ljk * ik)			+ 
		    (16 * CA * Box6ijk * jk)		- 
		    (32 * CF * Box6ijk * jk)		- 
		    (2 * CA * Ljk * jk)			+ 
		    (6 * CF * Ljk * jk)			+ 
		    (16 * CA * Box6ijk * ij2 * jk)	/ ik2 - 
		    (32 * CF * Box6ijk * ij2 * jk)	/ ik2 + 
		    (32 * CA * Box6ijk * ij * jk)	/ ik - 
		    (64 * CF * Box6ijk * ij * jk)	/ ik - 
		    (2 * CA * Ljk * ij * jk)		/ ik + 
		    (4 * CF * Ljk * ij * jk)		/ ik
		    ) / ij_ik_2;

  qqbargLoops[7] = (
		    (8 * CA * Box6jik * ij)		- 
		    (16 * CF * Box6jik * ij)		+ 
		    (CA * Lij * ij)			- 
		    (2 * CF * Lij * ij)			+ 
		    (CA * Lik * ij)			+ 
		    (2 * CF * Lik * ij)			+ 
		    (CA * Lij * ij2)			/ ik - 
		    (2 * CF * Lij * ij2)		/ ik + 
		    (8 * CA * Box6jik * ik)		- 
		    (16 * CF * Box6jik * ik)		+ 
		    (CA * Lik * ik)			+ 
		    (2 * CF * Lik * ik)			+ 
		    (8 * CA * Box6jik * ij2)		/ jk - 
		    (16 * CF * Box6jik * ij2)		/ jk + 
		    (8 * CA * Box6jik * ij * ik)	/ jk - 
		    (16 * CF * Box6jik * ij * ik)	/ jk - 
		    (24 * CA * Box6ikj * jk)		+ 
		    (CA * Lij * jk)			- 
		    (2 * CF * Lij * jk)			+ 
		    (CA * Lik * jk)			+ 
		    (4 * CF * Lik * jk)			+ 
		    (CA * Ljk * jk)			+ 
		    (4 * CF * Ljk * jk)			- 
		    (8 * CA * Box6ijk * ij2 * jk)	/ ik2 + 
		    (16 * CF * Box6ijk * ij2 * jk)	/ ik2 - 
		    (8 * CA * Box6ijk * ij * jk)	/ ik + 
		    (16 * CF * Box6ijk * ij * jk)	/ ik + 
		    (CA * Lij * ij * jk)		/ ik - 
		    (2 * CF * Lij * ij * jk)		/ ik + 
		    (CA * Ljk * ij * jk)		/ ik - 
		    (2 * CF * Ljk * ij * jk)		/ ik - 
		    (24 * CA * Box6ikj * ik * jk)	/ ij + 
		    (CA * Lik * ik * jk)		/ ij + 
		    (4 * CF * Lik * ik * jk)		/ ij - 
		    (24 * CA * Box6ikj * jk2)		/ ij + 
		    (CA * Ljk * jk2)			/ ij + 
		    (4 * CF * Ljk * jk2)		/ ij - 
		    (8 * CA * Box6ijk * ij * jk2)	/ ik2 + 
		    (16 * CF * Box6ijk * ij * jk2)	/ ik2 - 
		    (8 * CA * Box6ijk * jk2)		/ ik + 
		    (16 * CF * Box6ijk * jk2)		/ ik + 
		    (CA * Ljk * jk2)			/ ik - 
		    (2 * CF * Ljk * jk2)		/ ik - 
		    (24 * CA * Box6ikj * ik * jk2)	/ ij2
		    ) / (ij_ik * ij_jk);

  qqbargLoops[8] = (
		    (-8 * CA * Box6ijk * ij)		+ 
		    (16 * CF * Box6ijk * ij)		- 
		    (CA * Lij * ij)			+ 
		    (2 * CF * Lij * ij)			- 
		    (CA * Ljk * ij)			- 
		    (2 * CF * Ljk * ij)			- 
		    (8 * CA * Box6ijk * ij2)		/ ik + 
		    (16 * CF * Box6ijk * ij2)		/ ik + 
		    (24 * CA * Box6ikj * ik)		- 
		    (CA * Lij * ik)			+ 
		    (2 * CF * Lij * ik)			- 
		    (CA * Lik * ik)			- 
		    (4 * CF * Lik * ik)			- 
		    (CA * Ljk * ik)			- 
		    (4 * CF * Ljk * ik)			+ 
		    (24 * CA * Box6ikj * ik2)		/ ij - 
		    (CA * Lik * ik2)			/ ij - 
		    (4 * CF * Lik * ik2)		/ ij + 
		    (8 * CA * Box6jik * ij2 * ik)	/ jk2 - 
		    (16 * CF * Box6jik * ij2 * ik)	/ jk2 + 
		    (8 * CA * Box6jik * ij * ik2)	/ jk2 - 
		    (16 * CF * Box6jik * ij * ik2)	/ jk2 - 
		    (CA * Lij * ij2)			/ jk + 
		    (2 * CF * Lij * ij2)		/ jk + 
		    (8 * CA * Box6jik * ij * ik)	/ jk - 
		    (16 * CF * Box6jik * ij * ik)	/ jk - 
		    (CA * Lij * ij * ik)		/ jk + 
		    (2 * CF * Lij * ij * ik)		/ jk - 
		    (CA * Lik * ij * ik)		/ jk + 
		    (2 * CF * Lik * ij * ik)		/ jk + 
		    (8 * CA * Box6jik * ik2)		/ jk - 
		    (16 * CF * Box6jik * ik2)		/ jk - 
		    (CA * Lik * ik2)			/ jk + 
		    (2 * CF * Lik * ik2)		/ jk - 
		    (8 * CA * Box6ijk * jk)		+ 
		    (16 * CF * Box6ijk * jk)		- 
		    (CA * Ljk * jk)			- 
		    (2 * CF * Ljk * jk)			- 
		    (8 * CA * Box6ijk * ij * jk)	/ ik + 
		    (16 * CF * Box6ijk * ij * jk)	/ ik + 
		    (24 * CA * Box6ikj * ik * jk)	/ ij - 
		    (CA * Ljk * ik * jk)		/ ij - 
		    (4 * CF * Ljk * ik * jk)		/ ij + 
		    (24 * CA * Box6ikj * ik2 * jk)	/ ij2
		    ) / (ij_ik * ij_jk);

  qqbargLoops[9] = (
		    (2 * CF * ij) 			- 
		    (32 * CA * Box6jik * ij) 		+ 
		    (64 * CF * Box6jik * ij) 		+ 
		    (4 * CA * Lij * ij) 		- 
		    (8 * CF * Lij * ij) 		- 
		    (4 * CF * Lik * ij) 		- 
		    (16 * CA * Box6jik * ik)		+ 
		    (32 * CF * Box6jik * ik)		+ 
		    (2 * CA * Lik * ik) 		- 
		    (6 * CF * Lik * ik) 		- 
		    (16 * CA * Box6jik * ij2 * ik)	/ jk2 + 
		    (32 * CF * Box6jik * ij2 * ik)	/ jk2 - 
		    (16 * CA * Box6jik * ij2)		/ jk + 
		    (32 * CF * Box6jik * ij2)		/ jk + 
		    (2 * CA * Lij * ij2)		/ jk - 
		    (4 * CF * Lij * ij2)		/ jk - 
		    (32 * CA * Box6jik * ij * ik)	/ jk + 
		    (64 * CF * Box6jik * ij * ik)	/ jk + 
		    (2 * CA * Lik * ij * ik)		/ jk - 
		    (4 * CF * Lik * ij * ik)		/ jk + 
		    (2 * CF * jk) - 
		    (16 * CA * Box6jik * jk) + 
		    (32 * CF * Box6jik * jk) + 
		    (2 * CA * Lij * jk) - 
		    (4 * CF * Lij * jk) - 
		    (4 * CF * Lik * jk)
		    ) / ij_jk_2;

  qqbargLoops[10] = (
		     (-8 * CA * Box6ijk * ij2 * jk)	+ 
		     (16 * CF * Box6ijk * ij2 * jk)	+ 
		     (2 * CA * Lij * ij2 * jk)		- 
		     (4 * CF * Lij * ij2 * jk)		- 
		     (CA * ij * ik * jk)			+ 
		     (2 * CF * ij * ik * jk)		- 
		     (8 * CA * Box6ijk * ij * ik * jk)	+ 
		     (16 * CF * Box6ijk * ij * ik * jk)	+ 
		     (3 * CA * Lij * ij * ik * jk)	- 
		     (6 * CF * Lij * ij * ik * jk)	+ 
		     (CA * Ljk * ij * ik * jk)		- 
		     (2 * CF * Ljk * ij * ik * jk)	- 
		     (CA * ik2 * jk)			+ 
		     (2 * CF * ik2 * jk)			+ 
		     (CA * Lij * ik2 * jk)		- 
		     (2 * CF * Lij * ik2 * jk)		+ 
		     (CA * Ljk * ik2 * jk)		- 
		     (CF * Ljk * ik2 * jk)		- 
		     (CA * ij * jk2)			+ 
		     (2 * CF * ij * jk2)			- 
		     (16 * CA * Box6ijk * ij * jk2)	+ 
		     (32 * CF * Box6ijk * ij * jk2)	+ 
		     (2 * CA * Lij * ij * jk2)		- 
		     (4 * CF * Lij * ij * jk2)		+ 
		     (2 * CA * Ljk * ij * jk2)		- 
		     (4 * CF * Ljk * ij * jk2)		- 
		     (16 * CA * Box6ijk * ij2 * jk2)	/ ik + 
		     (32 * CF * Box6ijk * ij2 * jk2)	/ ik + 
		     (CA * Lij * ij2 * jk2)		/ ik - 
		     (2 * CF * Lij * ij2 * jk2)		/ ik - 
		     (CA * ik * jk2)			+ 
		     (2 * CF * ik * jk2)			+ 
		     (CA * Lij * ik * jk2)		- 
		     (2 * CF * Lij * ik * jk2)		+ 
		     (2 * CA * Ljk * ik * jk2)		- 
		     (2 * CF * Ljk * ik * jk2)		+ 
		     (CA * Ljk * jk3)			- 
		     (CF * Ljk * jk3)			- 
		     (8 * CA * Box6ijk * ij2 * jk3)	/ ik2 + 
		     (16 * CF * Box6ijk * ij2 * jk3)	/ ik2 - 
		     (8 * CA * Box6ijk * ij * jk3)	/ ik + 
		     (16 * CF * Box6ijk * ij * jk3)	/ ik + 
		     (CA * Ljk * ij * jk3)		/ ik - 
		     (2 * CF * Ljk * ij * jk3)		/ ik
		     ) / (ij_ik * ik_jk_2);

  qqbargLoops[11] = (
		     (16 * CA * Box6jik * ij2 * ik)	- 
		     (32 * CF * Box6jik * ij2 * ik)	- 
		     (CA * Lij * ij2 * ik)		+ 
		     (2 * CF * Lij * ij2 * ik)		+ 
		     (8 * CA * Box6jik * ij * ik2)	- 
		     (16 * CF * Box6jik * ij * ik2)	- 
		     (CA * Lik * ij * ik2)		+ 
		     (2 * CF * Lik * ij * ik2)		+ 
		     (8 * CA * Box6jik * ij2 * ik2)	/ jk - 
		     (16 * CF * Box6jik * ij2 * ik2)	/ jk + 
		     (8 * CA * Box6jik * ij2 * jk)	- 
		     (16 * CF * Box6jik * ij2 * jk)	- 
		     (2 * CA * Lij * ij2 * jk)		+ 
		     (4 * CF * Lij * ij2 * jk)		+ 
		     (CA * ij * ik * jk)			- 
		     (2 * CF * ij * ik * jk)		+ 
		     (16 * CA * Box6jik * ij * ik * jk)	- 
		     (32 * CF * Box6jik * ij * ik * jk)	- 
		     (2 * CA * Lij * ij * ik * jk)	+ 
		     (4 * CF * Lij * ij * ik * jk)	- 
		     (2 * CA * Lik * ij * ik * jk)	+ 
		     (4 * CF * Lik * ij * ik * jk)	- 
		     (CA * Lik * ik2 * jk)		+ 
		     (CF * Lik * ik2 * jk)		+ 
		     (CA * ij * jk2)			- 
		     (2 * CF * ij * jk2)			+ 
		     (8 * CA * Box6jik * ij * jk2)	- 
		     (16 * CF * Box6jik * ij * jk2)	- 
		     (3 * CA * Lij * ij * jk2)		+ 
		     (6 * CF * Lij * ij * jk2)		- 
		     (CA * Lik * ij * jk2)		+ 
		     (2 * CF * Lik * ij * jk2)		+ 
		     (CA * ik * jk2)			- 
		     (2 * CF * ik * jk2)			- 
		     (CA * Lij * ik * jk2)		+ 
		     (2 * CF * Lij * ik * jk2)		- 
		     (2 * CA * Lik * ik * jk2)		+ 
		     (2 * CF * Lik * ik * jk2)		+ 
		     (CA * jk3)				- 
		     (2 * CF * jk3)			- 
		     (CA * Lij * jk3)			+ 
		     (2 * CF * Lij * jk3)		- 
		     (CA * Lik * jk3)		 	+ 
		     (CF * Lik * jk3)
		     ) / (ij_jk * ik_jk_2);

  qqbargLoops[12] = -3 * CF * Lijk + (
				      (CA * ij2 * ik)			- 
				      (9 * CF * ij2 * ik)			+ 
				      (8 * CA * Box6ijk * ij2 * ik)	- 
				      (16 * CF * Box6ijk * ij2 * ik)	- 
				      (8 * CA * Box6ikj * ij2 * ik)	+ 
				      (CA * ij * ik2)			- 
				      (9 * CF * ij * ik2)			+ 
				      (8 * CA * Box6ijk * ij * ik2)	- 
				      (16 * CF * Box6ijk * ij * ik2)	- 
				      (8 * CA * Box6ikj * ij * ik2)	+ 
				      (CA * ij2 * jk)			- 
				      (9 * CF * ij2 * jk)			- 
				      (8 * CA * Box6ikj * ij2 * jk)	+ 
				      (8 * CA * Box6jik * ij2 * jk)	- 
				      (16 * CF * Box6jik * ij2 * jk)	+ 
				      (2 * CA * ij * ik * jk)		- 
				      (18 * CF * ij * ik * jk)		+ 
				      (8 * CA * Box6ijk * ij * ik * jk)	- 
				      (16 * CF * Box6ijk * ij * ik * jk)	- 
				      (40 * CA * Box6ikj * ij * ik * jk)	+ 
				      (8 * CA * Box6jik * ij * ik * jk)	- 
				      (16 * CF * Box6jik * ij * ik * jk)	+ 
				      (3 * CF * Lik * ij * ik * jk)	+ 
				      (3 * CF * Ljk * ij * ik * jk)	+ 
				      (CA * ik2 * jk)			- 
				      (9 * CF * ik2 * jk)			+ 
				      (8 * CA * Box6ijk * ik2 * jk)	- 
				      (16 * CF * Box6ijk * ik2 * jk)	- 
				      (32 * CA * Box6ikj * ik2 * jk)	+ 
				      (3 * CF * Lik * ik2 * jk)		+ 
				      (CA * ij * jk2)			- 
				      (9 * CF * ij * jk2)			- 
				      (8 * CA * Box6ikj * ij * jk2)	+ 
				      (8 * CA * Box6jik * ij * jk2)	- 
				      (16 * CF * Box6jik * ij * jk2)	+ 
				      (CA * ik * jk2)			- 
				      (9 * CF * ik * jk2)			- 
				      (32 * CA * Box6ikj * ik * jk2)	+ 
				      (8 * CA * Box6jik * ik * jk2)	- 
				      (16 * CF * Box6jik * ik * jk2)	+ 
				      (3 * CF * Ljk * ik * jk2)		- 
				      (24 * CA * Box6ikj * ik2 * jk2)	/ ij
				      ) / (ij_ik * ij_jk * ik_jk);

  /* // idendities implied by gauge invariance and current conservation; checked analytically and numerically
     Complex c1 = qqbargLoops[0] + qqbargLoops[6] + qqbargLoops[7];
     Complex c2 = qqbargLoops[1] + qqbargLoops[8] + qqbargLoops[9];
     Complex c3 = qqbargLoops[3] + qqbargLoops[4];
     Complex c4 = qqbargLoops[2] + qqbargLoops[5] + qqbargLoops[10] + qqbargLoops[11];
     Complex c5 = 
     2. * qqbargLoops[3]/ik +
     2. * qqbargLoops[5]/jk +
     qqbargLoops[6] * (1.+ij/ik) +
     qqbargLoops[8] * (jk+ij)/ik +
     2. * qqbargLoops[10] * (1./ik+1./jk) +
     2. * qqbargLoops[12] * (1./ik+1./jk);
     Complex c6 = 
     2. * qqbargLoops[4]/jk +
     2. * qqbargLoops[5]/jk +
     qqbargLoops[7] * (ik+ij)/jk +
     qqbargLoops[9] * (1.+ij/jk) +
     2. * qqbargLoops[11] * (ik/jk2+1./jk);
     Complex c7 =
     0.5 * qqbargLoops[0] * (ij+ik) +
     0.5 * qqbargLoops[1] * (ij+jk) +
     qqbargLoops[2] * (1.+ik/jk) -
     qqbargLoops[12] * (1.+ik/jk);

     double x1 = c1 != 0. ? log(abs(real(c1 * conj(c1)))) : 0.;
     double x2 = c2 != 0. ? log(abs(real(c2 * conj(c2)))) : 0.;
     double x3 = c3 != 0. ? log(abs(real(c3 * conj(c3)))) : 0.;
     double x4 = c4 != 0. ? log(abs(real(c4 * conj(c4)))) : 0.;
     double x5 = c5 != 0. ? log(abs(real(c5 * conj(c5)))) : 0.;
     double x6 = c6 != 0. ? log(abs(real(c6 * conj(c6)))) : 0.;
     double x7 = c7 != 0. ? log(abs(real(c7 * conj(c7)))) : 0.;

     cerr << x1 << " " << x2 << " " << x3 << " " << x4 << " "
     << x5 << " " << x6 << " " << x7 << "\n";
  */

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralLeftLoopCurrent(const int i, const int,
								      const int j, const int,
								      const int k, const int gHel,
								      const int n) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kn = plusProduct(k,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kn = minusProduct(k,n);
  
  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);

  const LorentzVector<Complex> & minusC_nk = minusCurrent(n,k);
  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);
  const LorentzVector<Complex> & minusC_kn = minusCurrent(k,n);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      (sqrt(2) * c6 * plusP_jk * minusC_nk * minusP_ik)/(jk * minusP_kn) + 
      (sqrt(2) * c1 * plusP_jk * momentum(i) * minusP_in)/minusP_kn + 
      (sqrt(2) * c2 * plusP_jk * momentum(j) * minusP_in)/minusP_kn + 
      (2 * sqrt(2) * c3 * plusP_jk * momentum(k) * minusP_in)/(jk * minusP_kn) + 
      (sqrt(2) * c4 * plusP_ik * minusC_ij * minusP_in)/(ik * minusP_kn) - 
      (sqrt(2) * c7 * plusP_ik * plusP_jk * momentum(i) * minusP_ik * minusP_in)/(ik * minusP_kn) - 
      (sqrt(2) * c9 * plusP_ik * plusP_jk * momentum(j) * minusP_ik * minusP_in)/(ik * minusP_kn) - 
      (2 * sqrt(2) * c11 * plusP_ik * plusP_jk * momentum(k) * minusP_ik * minusP_in)/(ik * jk * minusP_kn) + 
      (sqrt(2) * c5 * plusP_jk * minusC_ij * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c8 * sqr(plusP_jk) * momentum(i) * minusP_ik * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c10 * sqr(plusP_jk) * momentum(j) * minusP_ik * minusP_jn)/(jk * minusP_kn) - 
      (2 * sqrt(2) * c12 * sqr(plusP_jk) * momentum(k) * minusP_ik * minusP_jn)/(sqr(jk) * minusP_kn);
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2) * c1 * plusP_jn * momentum(i) * minusP_ik)/plusP_kn) - 
      (sqrt(2) * c2 * plusP_jn * momentum(j) * minusP_ik)/plusP_kn - 
      (2 * sqrt(2) * c3 * plusP_jn * momentum(k) * minusP_ik)/(jk * plusP_kn) - 
      (sqrt(2) * c4 * plusP_in * minusC_ij * minusP_ik)/(ik * plusP_kn) + 
      (sqrt(2) * c13 * minusC_kj * minusP_ik)/ik + (sqrt(2) * c13 * minusC_kj * minusP_ik)/jk - 
      (sqrt(2) * c6 * plusP_jk * minusC_kn * minusP_ik)/(jk * plusP_kn) + 
      (sqrt(2) * c7 * plusP_in * plusP_jk * momentum(i) * sqr(minusP_ik))/(ik * plusP_kn) + 
      (sqrt(2) * c9 * plusP_in * plusP_jk * momentum(j) * sqr(minusP_ik))/(ik * plusP_kn) + 
      (2 * sqrt(2) * c11 * plusP_in * plusP_jk * momentum(k) * sqr(minusP_ik))/(ik * jk * plusP_kn) - 
      (sqrt(2) * c5 * plusP_jn * minusC_ij * minusP_jk)/(jk * plusP_kn) + 
      (sqrt(2) * c8 * plusP_jk * plusP_jn * momentum(i) * minusP_ik * minusP_jk)/(jk * plusP_kn) + 
      (sqrt(2) * c10 * plusP_jk * plusP_jn * momentum(j) * minusP_ik * minusP_jk)/(jk * plusP_kn) + 
      (2 * sqrt(2) * c12 * plusP_jk * plusP_jn * momentum(k) * minusP_ik * minusP_jk)/(sqr(jk) * plusP_kn);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedLeftLoopCurrent(const int i, const int,
								    const int j, const int,
								    const int k, const int gHel) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_jk = plusProduct(j,k);

  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_ik = minusProduct(i,k);

  const LorentzVector<Complex> & minusC_ij = minusCurrent(i,j);
  const LorentzVector<Complex> & minusC_ik = minusCurrent(i,k);

  const LorentzVector<Complex> & minusC_kj = minusCurrent(k,j);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2) * c6 * plusP_jk * minusC_ik)/jk) 
      - (sqrt(2) * c8 * sqr(plusP_jk) * momentum(i) * minusP_ij)/jk - 
      (sqrt(2) * c10 * sqr(plusP_jk) * momentum(j) * minusP_ij)/jk - 
      (2 * sqrt(2) * c12 * sqr(plusP_jk) * momentum(k) * minusP_ij)/sqr(jk) + 
      (sqrt(2) * c5 * plusP_jk * minusC_ij * minusP_ij)/(jk * minusP_ik);
  }

  if ( gHel == -1 ) {
    return
      (sqrt(2) * c4 * plusP_ij * minusC_ij * minusP_ik)/(ik * plusP_jk) + 
      (sqrt(2) * c13 * minusC_kj * minusP_ik)/ik + (sqrt(2) * c13 * minusC_kj * minusP_ik)/jk + 
      (sqrt(2) * c6 * minusC_kj * minusP_ik)/jk - (sqrt(2) * c7 * plusP_ij * momentum(i)*
						   sqr(minusP_ik))/ik - 
      (sqrt(2) * c9 * plusP_ij * momentum(j) * sqr(minusP_ik))/ik - 
      (2 * sqrt(2) * c11 * plusP_ij * momentum(k) * sqr(minusP_ik))/(ik * jk);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargGeneralRightLoopCurrent(const int i,    const int,
								       const int j,    const int,
								       const int k,    const int gHel,
								       const int n) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ik = plusProduct(i,k);
  const Complex plusP_in = plusProduct(i,n);

  const Complex plusP_jk = plusProduct(j,k);
  const Complex plusP_jn = plusProduct(j,n);

  const Complex plusP_kn = plusProduct(k,n);

  const Complex minusP_ik = minusProduct(i,k);
  const Complex minusP_in = minusProduct(i,n);
  const Complex minusP_jk = minusProduct(j,k);
  const Complex minusP_jn = minusProduct(j,n);
  const Complex minusP_kn = minusProduct(k,n);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);

  const LorentzVector<Complex> & minusC_nk = minusCurrent(n,k);
  const LorentzVector<Complex> & minusC_kn = minusCurrent(k,n);

  Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2) * c13 * plusP_ik * minusC_jk)/ik) - (sqrt(2) * c13 * plusP_ik * minusC_jk)/jk + 
      (sqrt(2) * c4 * plusP_ik * minusC_ji * minusP_in)/(ik * minusP_kn) + 
      (sqrt(2) * c6 * plusP_ik * minusC_nk * minusP_jk)/(jk * minusP_kn) - 
      (sqrt(2) * c7 * sqr(plusP_ik) * momentum(i) * minusP_in * minusP_jk)/(ik * minusP_kn) - 
      (sqrt(2) * c9 * sqr(plusP_ik) * momentum(j) * minusP_in * minusP_jk)/(ik * minusP_kn) - 
      (2 * sqrt(2) * c11 * sqr(plusP_ik) * momentum(k) * minusP_in * minusP_jk)/(ik * jk * minusP_kn) + 
      (sqrt(2) * c1 * plusP_ik * momentum(i) * minusP_jn)/minusP_kn + 
      (sqrt(2) * c2 * plusP_ik * momentum(j) * minusP_jn)/minusP_kn + 
      (2 * sqrt(2) * c3 * plusP_ik * momentum(k) * minusP_jn)/(jk * minusP_kn) + 
      (sqrt(2) * c5 * plusP_jk * minusC_ji * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c8 * plusP_ik * plusP_jk * momentum(i) * minusP_jk * minusP_jn)/(jk * minusP_kn) - 
      (sqrt(2) * c10 * plusP_ik * plusP_jk * momentum(j) * minusP_jk * minusP_jn)/(jk * minusP_kn) - 
      (2 * sqrt(2) * c12 * plusP_ik * plusP_jk * momentum(k) * minusP_jk * minusP_jn)/(sqr(jk) * minusP_kn);
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2) * c4 * plusP_in * minusC_ji * minusP_ik)/(ik * plusP_kn)) - 
      (sqrt(2) * c1 * plusP_in * momentum(i) * minusP_jk)/plusP_kn - 
      (sqrt(2) * c2 * plusP_in * momentum(j) * minusP_jk)/plusP_kn - 
      (2 * sqrt(2) * c3 * plusP_in * momentum(k) * minusP_jk)/(jk * plusP_kn) - 
      (sqrt(2) * c5 * plusP_jn * minusC_ji * minusP_jk)/(jk * plusP_kn) - 
      (sqrt(2) * c6 * plusP_ik * minusC_kn * minusP_jk)/(jk * plusP_kn) + 
      (sqrt(2) * c7 * plusP_ik * plusP_in * momentum(i) * minusP_ik * minusP_jk)/(ik * plusP_kn) + 
      (sqrt(2) * c9 * plusP_ik * plusP_in * momentum(j) * minusP_ik * minusP_jk)/(ik * plusP_kn) + 
      (2 * sqrt(2) * c11 * plusP_ik * plusP_in * momentum(k) * minusP_ik * minusP_jk)/
      (ik * jk * plusP_kn) + 
      (sqrt(2) * c8 * plusP_ik * plusP_jn * momentum(i) * sqr(minusP_jk))/(jk * plusP_kn) + 
      (sqrt(2) * c10 * plusP_ik * plusP_jn * momentum(j) * sqr(minusP_jk))/(jk * plusP_kn) + 
      (2 * sqrt(2) * c12 * plusP_ik * plusP_jn * momentum(k) * sqr(minusP_jk))/(sqr(jk) * plusP_kn);
  }

  return czero;

}

LorentzVector<Complex> MatchboxCurrents::qqbargFixedRightLoopCurrent(const int i, const int,
								     const int j, const int,
								     const int k, const int gHel) {

  qqbargLoopCoefficients(i,j,k);

  const double ik = invariant(i,k);
  const double jk = invariant(j,k);

  const Complex plusP_ij = plusProduct(i,j);
  const Complex plusP_ik = plusProduct(i,k);
  
  const Complex minusP_ij = minusProduct(i,j);
  const Complex minusP_jk = minusProduct(j,k);

  const LorentzVector<Complex> & minusC_ji = minusCurrent(j,i);
  const LorentzVector<Complex> & minusC_jk = minusCurrent(j,k);

  const LorentzVector<Complex> & minusC_ki = minusCurrent(k,i);

  //Complex c1  = qqbargLoops[0]; Complex c2  = qqbargLoops[1];  Complex c3  = qqbargLoops[2];
  Complex c4  = qqbargLoops[3]; Complex c5  = qqbargLoops[4];  Complex c6  = qqbargLoops[5];
  Complex c7  = qqbargLoops[6]; Complex c8  = qqbargLoops[7];  Complex c9  = qqbargLoops[8];
  Complex c10 = qqbargLoops[9]; Complex c11 = qqbargLoops[10]; Complex c12 = qqbargLoops[11];
  Complex c13 = qqbargLoops[12];

  if ( gHel == 1 ) {
    return
      -((sqrt(2) * c13 * plusP_ik * minusC_jk)/ik) - 
      (sqrt(2) * c13 * plusP_ik * minusC_jk)/jk - 
      (sqrt(2) * c6 * plusP_ik * minusC_jk)/jk + 
      (sqrt(2) * c7 * sqr(plusP_ik) * momentum(i) * minusP_ij)/ik + 
      (sqrt(2) * c9 * sqr(plusP_ik) * momentum(j) * minusP_ij)/ik + 
      (2 * sqrt(2) * c11 * sqr(plusP_ik) * momentum(k) * minusP_ij)/(ik * jk) - 
      (sqrt(2) * c4 * plusP_ik * minusC_ji * minusP_ij)/(ik * minusP_jk);
  }

  if ( gHel == -1 ) {
    return
      -((sqrt(2) * c5 * plusP_ij * minusC_ji * minusP_jk)/(jk * plusP_ik)) + 
      (sqrt(2) * c6 * minusC_ki * minusP_jk)/jk + 
      (sqrt(2) * c8 * plusP_ij * momentum(i) * sqr(minusP_jk))/jk + 
      (sqrt(2) * c10 * plusP_ij * momentum(j) * sqr(minusP_jk))/jk + 
      (2 * sqrt(2) * c12 * plusP_ij * momentum(k) * sqr(minusP_jk))/sqr(jk);
  }

  return czero;

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargLeftOneLoopCurrent(const int q,    const int qHel,
									 const int qbar, const int qbarHel,
									 const int g1,   const int g1Hel) {
  if ( qHel != 1 || qbarHel != 1 )
    return czero;

  if ( getCurrent(hash<2>(1,2,q,qHel,qbar,qbarHel,g1,g1Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,q);
    LorentzVector<Complex> nj = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,qbar);
    LorentzVector<Complex> nl = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,0);
    LorentzVector<Complex> nlbar = Complex(0.,0.5) * qqbargGeneralLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,1);
    LorentzVector<Complex> fixed = Complex(0.,0.5) * qqbargFixedLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
    ostream& ncheck = checkStream("qqbargLeftLoopCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(Complex(0.,0.5) * qqbargFixedLeftLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargLeftLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1));
#endif
  return cachedCurrent();

}

const LorentzVector<Complex>& MatchboxCurrents::qqbargRightOneLoopCurrent(const int q,    const int qHel,
									  const int qbar, const int qbarHel,
									  const int g1,   const int g1Hel) {

  if ( qHel != -1 || qbarHel != -1 )
    return czero;

  if ( getCurrent(hash<2>(2,2,q,qHel,qbar,qbarHel,g1,g1Hel)) ) {
#ifdef CHECK_MatchboxCurrents
    LorentzVector<Complex> ni = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,q);
    LorentzVector<Complex> nj = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,qbar);
    LorentzVector<Complex> nl = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,0);
    LorentzVector<Complex> nlbar = Complex(0.,0.5) * qqbargGeneralRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel,1);
    LorentzVector<Complex> fixed = Complex(0.,0.5) * qqbargFixedRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel);
    LorentzVector<Complex> x1 = fixed - ni;
    LorentzVector<Complex> x2 = fixed - nj;
    LorentzVector<Complex> x3 = fixed - nl;
    LorentzVector<Complex> x4 = fixed - nlbar;
    double c1 = 
      real(x1.t() * conj(x1.t())) + real(x1.x() * conj(x1.x())) + real(x1.y() * conj(x1.y())) + real(x1.z() * conj(x1.z()));
    double c2 = 
      real(x2.t() * conj(x2.t())) + real(x2.x() * conj(x2.x())) + real(x2.y() * conj(x2.y())) + real(x2.z() * conj(x2.z()));
    double c3 = 
      real(x3.t() * conj(x3.t())) + real(x3.x() * conj(x3.x())) + real(x3.y() * conj(x3.y())) + real(x3.z() * conj(x3.z()));
    double c4 = 
      real(x4.t() * conj(x4.t())) + real(x4.x() * conj(x4.x())) + real(x4.y() * conj(x4.y())) + real(x4.z() * conj(x4.z()));
    ostream& ncheck = checkStream("qqbargRightLoopCurrentNChoice");
    ncheck << (c1 != 0. ? log10(abs(c1)) : 0.) << " "
	   << (c2 != 0. ? log10(abs(c2)) : 0.) << " "
	   << (c3 != 0. ? log10(abs(c3)) : 0.) << " "
	   << (c4 != 0. ? log10(abs(c4)) : 0.) << " "
	   << "\n" << flush;
#endif
    cacheCurrent(Complex(0.,0.5) * qqbargFixedRightLoopCurrent(q,qHel,qbar,qbarHel,g1,g1Hel));
  }
#ifdef CHECK_MatchboxCurrents
  checkCurrent("qqbargRightLoopCurrent",cachedCurrent(),momentum(q)+momentum(qbar)+momentum(g1));
#endif
  return cachedCurrent();

}

#ifdef CHECK_MatchboxCurrents

map<string,ofstream * >& MatchboxCurrents::checkStreams() {
  static map<string,ofstream * > theMap;
  return theMap;
}

ostream& MatchboxCurrents::checkStream(const string& id) {
  map<string,ofstream * >::iterator ret = checkStreams().find(id);
  if ( ret == checkStreams().end() ) {
    checkStreams()[id] = new ofstream(id.c_str());
    ret = checkStreams().find(id);
  }
  return *(ret->second);
}

void MatchboxCurrents::checkCurrent(const string& id,
				    const LorentzVector<Complex>& current,
				    const LorentzVector<double>& q) {
  
  Complex c = current.dot(q);
  double ac = abs(real(conj(c) * c));

  if ( ! isfinite(ac) ) {
    cerr << "ooops ... nan encountered in current conservation\n" << flush;
    return;
  }

  checkStream(id) << (ac > 0. ? log10(ac) : 0.) << "\n";

}

#endif // CHECK_MatchboxCurrents
#line 1 "./MatchboxZGammaAmplitude.cc"
// -*- C++ -*-
//
// MatchboxZGammaAmplitude.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxZGammaAmplitude class.
//

#include "MatchboxZGammaAmplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MatchboxZGammaAmplitude::MatchboxZGammaAmplitude() 
  : MatchboxAmplitude(), theIncludeZ(true), theIncludeGamma(true) {}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxZGammaAmplitude::persistentOutput(PersistentOStream & os) const {
  os << theIncludeZ << theIncludeGamma;
}

void MatchboxZGammaAmplitude::persistentInput(PersistentIStream & is, int) {
  is >> theIncludeZ >> theIncludeGamma;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeAbstractClass<MatchboxZGammaAmplitude,MatchboxAmplitude>
  describeHerwigMatchboxZGammaAmplitude("Herwig::MatchboxZGammaAmplitude", "HwMatchboxBuiltin.so");

void MatchboxZGammaAmplitude::Init() {

  static ClassDocumentation<MatchboxZGammaAmplitude> documentation
    ("There is no documentation for the MatchboxZGammaAmplitude class");

  static Switch<MatchboxZGammaAmplitude,bool> interfaceIncludeZ
    ("IncludeZ",
     "Include the Z contribution.",
     &MatchboxZGammaAmplitude::theIncludeZ, true, false, false);
  static SwitchOption interfaceIncludeZYes
    (interfaceIncludeZ,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIncludeZNo
    (interfaceIncludeZ,
     "No",
     "",
     false);

  static Switch<MatchboxZGammaAmplitude,bool> interfaceIncludeGamma
    ("IncludeGamma",
     "Include the photon contribution.",
     &MatchboxZGammaAmplitude::theIncludeGamma, true, false, false);
  static SwitchOption interfaceIncludeGammaYes
    (interfaceIncludeGamma,
     "Yes",
     "",
     true);
  static SwitchOption interfaceIncludeGammaNo
    (interfaceIncludeGamma,
     "No",
     "",
     false);

}

#line 1 "./MatchboxAmplitudellbarqqbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudellbarqqbar class.
//

#include "MatchboxAmplitudellbarqqbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudellbarqqbar::MatchboxAmplitudellbarqqbar() {}

IBPtr MatchboxAmplitudellbarqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudellbarqqbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudellbarqqbar::doinit() {
  MatchboxZGammaAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(4);
}

void MatchboxAmplitudellbarqqbar::doinitrun() {
  MatchboxZGammaAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(4);
}

bool MatchboxAmplitudellbarqqbar::canHandle(const PDVector& proc) const {
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator lepton = xproc.begin();
  long leptonId = 0;
  for ( ; lepton != xproc.end(); ++lepton )
    if ( (**lepton).id() == 11 ||
	 (**lepton).id() == 13 ||
	 (**lepton).id() == 15 ) {
      break;
    }
  if ( lepton == xproc.end() )
    return false;
  leptonId = (**lepton).id();
  xproc.erase(lepton);
  PDVector::iterator antiLepton = xproc.begin();
  for ( ; antiLepton != xproc.end(); ++antiLepton )
    if ( (**antiLepton).id() == -leptonId ) {
      break;
    }
  if ( antiLepton == xproc.end() )
    return false;
  xproc.erase(antiLepton);

  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 7 &&
	 (**quark).id() > 0 ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);

  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  return xproc.empty();
}

void MatchboxAmplitudellbarqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxZGammaAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));

  setupQuarks(2,amplitudeMomentum(2),
 	      3,amplitudeMomentum(3));

  MatchboxZGammaAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudellbarqqbar::evaluate(size_t, const vector<int>& hel, Complex& largeN) {

  if ( abs(hel[2])+abs(hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]);
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]);

  const LorentzVector<Complex>& quarkLeft
    = qqbarLeftCurrent(2,hel[2],3,hel[3]);
  const LorentzVector<Complex>& quarkRight
    = qqbarRightCurrent(2,hel[2],3,hel[3]);

  Complex LL = leptonLeft.dot( quarkLeft );
  Complex RL = leptonRight.dot( quarkLeft );
  Complex LR = leptonLeft.dot( quarkRight );
  Complex RR = leptonRight.dot( quarkRight );

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma = 0.0;
  if ( includeGamma() )
    gamma = Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
      (LL + RL + LR + RR)/bProp;

  bool up = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  Complex Z = 0.0;
  if ( includeZ() )
    Z = Complex(0.,-1.)*
      (standardModel()->le()*(up ? standardModel()->lu() : standardModel()->ld())*LL + 
       standardModel()->re()*(up ? standardModel()->lu() : standardModel()->ld())*RL + 
       standardModel()->le()*(up ? standardModel()->ru() : standardModel()->rd())*LR + 
       standardModel()->re()*(up ? standardModel()->ru() : standardModel()->rd())*RR)/
      Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex res = 4.*Constants::pi*SM().alphaEMMZ()*(gamma+Z);
  largeN = res;
  return res;

}

Complex MatchboxAmplitudellbarqqbar::evaluateOneLoop(size_t, const vector<int>& hel) {

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]);
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]);

  const LorentzVector<Complex>& quarkLeft
    = qqbarLeftOneLoopCurrent(2,hel[2],3,hel[3]);
  const LorentzVector<Complex>& quarkRight
    = qqbarRightOneLoopCurrent(2,hel[2],3,hel[3]);

  Complex LL = leptonLeft.dot( quarkLeft );
  Complex RL = leptonRight.dot( quarkLeft );
  Complex LR = leptonLeft.dot( quarkRight );
  Complex RR = leptonRight.dot( quarkRight );

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma = 0.0;
  if ( includeGamma() )
    gamma = Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
      (LL + RL + LR + RR)/bProp;

  bool up = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  Complex Z = 0.0;
  if ( includeZ() )
    Z = Complex(0.,-1.)*
      (standardModel()->le()*(up ? standardModel()->lu() : standardModel()->ld())*LL +
       standardModel()->re()*(up ? standardModel()->lu() : standardModel()->ld())*RL +
       standardModel()->le()*(up ? standardModel()->ru() : standardModel()->rd())*LR +
       standardModel()->re()*(up ? standardModel()->ru() : standardModel()->rd())*RR)/
      Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex res = (SM().alphaS()/(2.*Constants::pi))*4.*Constants::pi*SM().alphaEMMZ()*(gamma+Z);
  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudellbarqqbar::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudellbarqqbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudellbarqqbar,MatchboxZGammaAmplitude>
describeHerwigMatchboxAmplitudellbarqqbar("Herwig::MatchboxAmplitudellbarqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudellbarqqbar::Init() {

  static ClassDocumentation<MatchboxAmplitudellbarqqbar> documentation
    ("MatchboxAmplitudellbarqqbar");

}

#line 1 "./MatchboxAmplitudellbarqqbarg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudellbarqqbarg class.
//

#include "MatchboxAmplitudellbarqqbarg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudellbarqqbarg::MatchboxAmplitudellbarqqbarg() {}

IBPtr MatchboxAmplitudellbarqqbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudellbarqqbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudellbarqqbarg::doinit() {
  MatchboxZGammaAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(5);
}

void MatchboxAmplitudellbarqqbarg::doinitrun() {
  MatchboxZGammaAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(5);
}

bool MatchboxAmplitudellbarqqbarg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 5 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator lepton = xproc.begin();
  long leptonId = 0;
  for ( ; lepton != xproc.end(); ++lepton )
    if ( (**lepton).id() == 11 ||
	 (**lepton).id() == 13 ||
	 (**lepton).id() == 15 ) {
      break;
    }
  if ( lepton == xproc.end() )
    return false;
  leptonId = (**lepton).id();
  xproc.erase(lepton);
  PDVector::iterator antiLepton = xproc.begin();
  for ( ; antiLepton != xproc.end(); ++antiLepton )
    if ( (**antiLepton).id() == -leptonId ) {
      break;
    }
  if ( antiLepton == xproc.end() )
    return false;
  xproc.erase(antiLepton);

  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 7 &&
	 (**quark).id() > 0 ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);

  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  if ( xproc.size() != 1 )
    return false;
  return xproc[0]->id() == 21;
}


void MatchboxAmplitudellbarqqbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxZGammaAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));

  setupQuarks(2,amplitudeMomentum(2),
	      3,amplitudeMomentum(3));

  momentum(4,amplitudeMomentum(4));

  MatchboxZGammaAmplitude::prepareAmplitudes(me);

}


Complex MatchboxAmplitudellbarqqbarg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {

  assert(nPoints() == 5);

  if ( abs(hel[2])+abs(hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]); 

  const LorentzVector<Complex>& quarkLeft
    = qqbargLeftCurrent(2,hel[2],3,hel[3],4,hel[4]);
  const LorentzVector<Complex>& quarkRight
    = qqbargRightCurrent(2,hel[2],3,hel[3],4,hel[4]);

  Complex LL = leptonLeft.dot( quarkLeft );
  Complex RL = leptonRight.dot( quarkLeft );
  Complex LR = leptonLeft.dot( quarkRight );
  Complex RR = leptonRight.dot( quarkRight );

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma = 0.0;
  if ( includeGamma() )
    gamma = Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
      (LL + RL + LR + RR)/bProp;

  bool up = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  Complex Z = 0.0;
  if ( includeZ() )
    Z = Complex(0.,-1.)*
      (standardModel()->le()*(up ? standardModel()->lu() : standardModel()->ld())*LL +
       standardModel()->re()*(up ? standardModel()->lu() : standardModel()->ld())*RL +
       standardModel()->le()*(up ? standardModel()->ru() : standardModel()->rd())*LR +
       standardModel()->re()*(up ? standardModel()->ru() : standardModel()->rd())*RR)/
      Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat()); 

  Complex res = 4.*Constants::pi*SM().alphaEMMZ()*sqrt(4.*Constants::pi*SM().alphaS())*(gamma+Z);
  largeN = res;
  return res;

}


Complex MatchboxAmplitudellbarqqbarg::evaluateOneLoop(size_t, const vector<int>& hel) {

  if ( abs(hel[2]+hel[3]) != 2 ) {
    return 0.;
  }

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]); 

  Complex LL =
    hel[2] ==  1 ? leptonLeft.dot( qqbargLeftOneLoopCurrent (2,hel[2],3,hel[3],4,hel[4]))  : 0.;
  Complex RL =
    hel[2] ==  1 ? leptonRight.dot(qqbargLeftOneLoopCurrent (2,hel[2],3,hel[3],4,hel[4]))  : 0.;
  Complex LR =
    hel[2] == -1 ? leptonLeft.dot( qqbargRightOneLoopCurrent(2,hel[2],3,hel[3],4,hel[4]))  : 0.;
  Complex RR =
    hel[2] == -1 ? leptonRight.dot(qqbargRightOneLoopCurrent(2,hel[2],3,hel[3],4,hel[4]))  : 0.;

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma = 0.0;
  if ( includeGamma() )
    gamma = Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
      (LL + RL + LR + RR)/bProp;

  bool up = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  Complex Z = 0.0;
  if ( includeZ() )
    Z = Complex(0.,-1.)*
      (standardModel()->le()*(up ? standardModel()->lu() : standardModel()->ld())*LL +
       standardModel()->re()*(up ? standardModel()->lu() : standardModel()->ld())*RL +
       standardModel()->le()*(up ? standardModel()->ru() : standardModel()->rd())*LR +
       standardModel()->re()*(up ? standardModel()->ru() : standardModel()->rd())*RR)/
      Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex res = 
    (SM().alphaS()/(2.*Constants::pi))*
    4.*Constants::pi*SM().alphaEMMZ()*sqrt(4.*Constants::pi*SM().alphaS())*(gamma+Z);
  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudellbarqqbarg::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudellbarqqbarg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudellbarqqbarg,MatchboxZGammaAmplitude>
  describeHerwigMatchboxAmplitudellbarqqbarg("Herwig::MatchboxAmplitudellbarqqbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudellbarqqbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudellbarqqbarg> documentation
    ("MatchboxAmplitudellbarqqbarg");

}

#line 1 "./MatchboxAmplitudellbarqqbargg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbargg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudellbarqqbargg class.
//

#include "MatchboxAmplitudellbarqqbargg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudellbarqqbargg::MatchboxAmplitudellbarqqbargg() {}

IBPtr MatchboxAmplitudellbarqqbargg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudellbarqqbargg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudellbarqqbargg::doinit() {
  MatchboxZGammaAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

void MatchboxAmplitudellbarqqbargg::doinitrun() {
  MatchboxZGammaAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

bool MatchboxAmplitudellbarqqbargg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 6 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator lepton = xproc.begin();
  long leptonId = 0;
  for ( ; lepton != xproc.end(); ++lepton )
    if ( (**lepton).id() == 11 ||
	 (**lepton).id() == 13 ||
	 (**lepton).id() == 15 ) {
      break;
    }
  if ( lepton == xproc.end() )
    return false;
  leptonId = (**lepton).id();
  xproc.erase(lepton);
  PDVector::iterator antiLepton = xproc.begin();
  for ( ; antiLepton != xproc.end(); ++antiLepton )
    if ( (**antiLepton).id() == -leptonId ) {
      break;
    }
  if ( antiLepton == xproc.end() )
    return false;
  xproc.erase(antiLepton);
  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0 &&
	 (**quark).hardProcessMass() == ZERO ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);
  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  if ( xproc.size() != 2 )
    return false;
  return xproc[0]->id() == 21 && xproc[1]->id() == 21;
}

void MatchboxAmplitudellbarqqbargg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxZGammaAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));

  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  momentum(5,amplitudeMomentum(5));

  MatchboxZGammaAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudellbarqqbargg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]); 

  assert(amplitudeToColourMap()[2] == 0 &&
	 amplitudeToColourMap()[3] == 1);

  int g1,hg1,g2,hg2;
  if ( amplitudeToColourMap()[4] == 2 &&
       amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else if ( a == 1 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else assert(false);
  } else if ( amplitudeToColourMap()[4] == 3 &&
	      amplitudeToColourMap()[5] == 2 ) {
    if ( a == 0 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else if ( a == 1 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else assert(false);
  } else assert(false);

  Complex LL =
    hel[2] ==  1 ? leptonLeft.dot( qqbarggLeftCurrent (2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;
  Complex RL =
    hel[2] ==  1 ? leptonRight.dot(qqbarggLeftCurrent (2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;
  Complex LR =
    hel[2] == -1 ? leptonLeft.dot( qqbarggRightCurrent(2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;
  Complex RR =
    hel[2] == -1 ? leptonRight.dot(qqbarggRightCurrent(2,hel[2],3,hel[3],g1,hg1,g2,hg2))  : 0.;

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma = 0.0;
  if ( includeGamma() )
    gamma = Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
      (LL + RL + LR + RR)/bProp;

  bool up = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  Complex Z = 0.0;
  if ( includeZ() )
    Z = Complex(0.,-1.)*
      (standardModel()->le()*(up ? standardModel()->lu() : standardModel()->ld())*LL +
       standardModel()->re()*(up ? standardModel()->lu() : standardModel()->ld())*RL +
       standardModel()->le()*(up ? standardModel()->ru() : standardModel()->rd())*LR +
       standardModel()->re()*(up ? standardModel()->ru() : standardModel()->rd())*RR)/
      Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex res = sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*(gamma+Z);
  largeN = res;
  return res;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudellbarqqbargg::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudellbarqqbargg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudellbarqqbargg,MatchboxZGammaAmplitude>
  describeHerwigMatchboxAmplitudellbarqqbargg("Herwig::MatchboxAmplitudellbarqqbargg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudellbarqqbargg::Init() {

  static ClassDocumentation<MatchboxAmplitudellbarqqbargg> documentation
    ("MatchboxAmplitudellbarqqbargg");

}

#line 1 "./MatchboxAmplitudellbarqqbarqqbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudellbarqqbarqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudellbarqqbarqqbar class.
//

#include "MatchboxAmplitudellbarqqbarqqbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudellbarqqbarqqbar::MatchboxAmplitudellbarqqbarqqbar() {}

IBPtr MatchboxAmplitudellbarqqbarqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudellbarqqbarqqbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudellbarqqbarqqbar::doinit() {
  MatchboxZGammaAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

void MatchboxAmplitudellbarqqbarqqbar::doinitrun() {
  MatchboxZGammaAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  nPoints(6);
}

bool MatchboxAmplitudellbarqqbarqqbar::canHandle(const PDVector& proc) const {
  if ( proc.size() != 6 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator lepton = xproc.begin();
  long leptonId = 0;
  for ( ; lepton != xproc.end(); ++lepton )
    if ( (**lepton).id() == 11 ||
	 (**lepton).id() == 13 ||
	 (**lepton).id() == 15 ) {
      break;
    }
  if ( lepton == xproc.end() )
    return false;
  leptonId = (**lepton).id();
  xproc.erase(lepton);
  PDVector::iterator antiLepton = xproc.begin();
  for ( ; antiLepton != xproc.end(); ++antiLepton )
    if ( (**antiLepton).id() == -leptonId ) {
      break;
    }
  if ( antiLepton == xproc.end() )
    return false;
  xproc.erase(antiLepton);
  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0 &&
	 (**quark).hardProcessMass() == ZERO ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);
  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  quark = xproc.begin();
  quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0 &&
	 (**quark).hardProcessMass() == ZERO ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);
  antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  return xproc.empty();
}

inline bool leftNonZero(int heli, int helj, int helk, int hell) {
  return 
    heli == 1 && helj == 1 &&
    abs(helk+hell) == 2;
}

inline bool rightNonZero(int heli, int helj, int helk, int hell) {
  return 
    heli == -1 && helj == -1 &&
    abs(helk+hell) == 2;
}

void MatchboxAmplitudellbarqqbarqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxZGammaAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));

  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  momentum(5,amplitudeMomentum(5));

  MatchboxZGammaAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudellbarqqbarqqbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  const LorentzVector<Complex>& leptonLeft
    = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& leptonRight
    = llbarRightCurrent(0,hel[0],1,hel[1]); 

  Complex LL2345 =
    leftNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonLeft.dot(qqbarqqbarLeftCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex LL4523 =
    leftNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonLeft.dot(qqbarqqbarLeftCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex LL2543 =
    leftNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonLeft.dot(qqbarqqbarLeftCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex LL4325 =
    leftNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonLeft.dot(qqbarqqbarLeftCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  Complex LR2345 =
    rightNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonLeft.dot(qqbarqqbarRightCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex LR4523 =
    rightNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonLeft.dot(qqbarqqbarRightCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex LR2543 =
    rightNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonLeft.dot(qqbarqqbarRightCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex LR4325 =
    rightNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonLeft.dot(qqbarqqbarRightCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  Complex RL2345 =
    leftNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonRight.dot(qqbarqqbarLeftCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex RL4523 =
    leftNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonRight.dot(qqbarqqbarLeftCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex RL2543 =
    leftNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonRight.dot(qqbarqqbarLeftCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex RL4325 =
    leftNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonRight.dot(qqbarqqbarLeftCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  Complex RR2345 =
    rightNonZero(hel[2],hel[3],hel[4],hel[5]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonRight.dot(qqbarqqbarRightCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex RR4523 =
    rightNonZero(hel[4],hel[5],hel[2],hel[3]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonRight.dot(qqbarqqbarRightCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex RR2543 =
    rightNonZero(hel[2],hel[5],hel[4],hel[3]) &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonRight.dot(qqbarqqbarRightCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex RR4325 =
    rightNonZero(hel[4],hel[3],hel[2],hel[5]) &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonRight.dot(qqbarqqbarRightCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;

  double bProp = (amplitudeMomentum(0)+amplitudeMomentum(1)).m2()/lastSHat();

  Complex gamma2345 =
    Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
    (LL2345 + RL2345 + LR2345 + RR2345)/bProp;
  Complex gamma2543 =
    Complex(0.,-1.)*(-amplitudePartonData()[2]->iCharge()/3.)*
    (LL2543 + RL2543 + LR2543 + RR2543)/bProp;
  Complex gamma4523 =
    Complex(0.,-1.)*(-amplitudePartonData()[4]->iCharge()/3.)*
    (LL4523 + RL4523 + LR4523 + RR4523)/bProp;
  Complex gamma4325 =
    Complex(0.,-1.)*(-amplitudePartonData()[4]->iCharge()/3.)*
    (LL4325 + RL4325 + LR4325 + RR4325)/bProp;

  bool up2 = abs(amplitudePartonData()[2]->id()) % 2 == 0;
  bool up4 = abs(amplitudePartonData()[4]->id()) % 2 == 0;

  Complex Z2345 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up2 ? standardModel()->lu() : standardModel()->ld())*LL2345 +
     standardModel()->re()*(up2 ? standardModel()->lu() : standardModel()->ld())*RL2345 +
     standardModel()->le()*(up2 ? standardModel()->ru() : standardModel()->rd())*LR2345 +
     standardModel()->re()*(up2 ? standardModel()->ru() : standardModel()->rd())*RR2345)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());
  Complex Z2543 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up2 ? standardModel()->lu() : standardModel()->ld())*LL2543 +
     standardModel()->re()*(up2 ? standardModel()->lu() : standardModel()->ld())*RL2543 +
     standardModel()->le()*(up2 ? standardModel()->ru() : standardModel()->rd())*LR2543 +
     standardModel()->re()*(up2 ? standardModel()->ru() : standardModel()->rd())*RR2543)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());
  Complex Z4523 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up4 ? standardModel()->lu() : standardModel()->ld())*LL4523 +
     standardModel()->re()*(up4 ? standardModel()->lu() : standardModel()->ld())*RL4523 +
     standardModel()->le()*(up4 ? standardModel()->ru() : standardModel()->rd())*LR4523 +
     standardModel()->re()*(up4 ? standardModel()->ru() : standardModel()->rd())*RR4523)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());
  Complex Z4325 =
    Complex(0.,-1.)*
    (standardModel()->le()*(up4 ? standardModel()->lu() : standardModel()->ld())*LL4325 +
     standardModel()->re()*(up4 ? standardModel()->lu() : standardModel()->ld())*RL4325 +
     standardModel()->le()*(up4 ? standardModel()->ru() : standardModel()->rd())*LR4325 +
     standardModel()->re()*(up4 ? standardModel()->ru() : standardModel()->rd())*RR4325)/
    Complex(bProp-sqr(MZ)/lastSHat(),MZ*GZ/lastSHat());

  Complex sum2345 = 0.0;
  Complex sum2543 = 0.0;
  Complex sum4523 = 0.0;
  Complex sum4325 = 0.0;

  if ( includeGamma() ) {
    sum2345 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma2345;
    sum2543 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma2543;
    sum4523 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma4523;
    sum4325 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*gamma4325;
  }

  if ( includeZ() ) {
    sum2345 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z2345;
    sum2543 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z2543;
    sum4523 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z4523;
    sum4325 += sqr(4.*Constants::pi)*SM().alphaEMMZ()*SM().alphaS()*Z4325;
  }

  double Nc = SM().Nc();

  Complex resLeading = 0.;
  Complex resSubLeading = 0.;

  if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 1 &&
       amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else if ( a == 1 ) {  //(25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { // (25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else if ( a == 1 ) { // (23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else if ( a == 1 ) { //(25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 1 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(25)(43)
      resLeading = sum2345 + sum4523;
      resSubLeading = sum2543 + sum4325;
    } else if ( a == 1 ) { //(23)(45)
      resLeading = sum2543 + sum4325;
      resSubLeading = sum2345 + sum4523;
    } else assert(false);
  } else assert(false);

  resSubLeading *= -1./Nc;
  largeN = resLeading/2.;

  return (resLeading + resSubLeading)/2.;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudellbarqqbarqqbar::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudellbarqqbarqqbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudellbarqqbarqqbar,MatchboxZGammaAmplitude>
  describeHerwigMatchboxAmplitudellbarqqbarqqbar("Herwig::MatchboxAmplitudellbarqqbarqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudellbarqqbarqqbar::Init() {

  static ClassDocumentation<MatchboxAmplitudellbarqqbarqqbar> documentation
    ("MatchboxAmplitudellbarqqbarqqbar");

}

#line 1 "./MatchboxAmplitudelnuqqbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbar class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbar.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbar::MatchboxAmplitudelnuqqbar() 
  : theDiagonal(false) {}

void MatchboxAmplitudelnuqqbar::doinit() {
  MatchboxAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(4);  
}

void MatchboxAmplitudelnuqqbar::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(4);
}

IBPtr MatchboxAmplitudelnuqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbar::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbar::canHandle(const PDVector& proc) const {
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator elektron = xproc.begin();
  for ( ; elektron != xproc.end(); ++elektron ) 
    if ( abs((*elektron)->id()) >= 11 && abs((*elektron)->id()) <= 16 && abs((*elektron)->id()) % 2 == 1 ) break;
  if ( elektron == xproc.end() ) return false;
  PDPtr e = *elektron;
  xproc.erase(elektron);
  PDVector::iterator neutrino = xproc.begin();
  for ( ; neutrino != xproc.end(); ++neutrino ) 
    if ( abs((*neutrino)->id()) >= 11 && abs((*neutrino)->id()) <= 16 && abs((*neutrino)->id()) % 2 == 0 ) break;
  if ( neutrino == xproc.end() ) return false;
  PDPtr n = *neutrino;
  xproc.erase(neutrino);
  PDVector::iterator quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 1 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr d = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 0 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr u = *quark;
  xproc.erase(quark);
  if ( u->iCharge() + d->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
  if ( theDiagonal && SU2Helper::family(u) != SU2Helper::family(d) ) return false;
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  setupLeptons(0,amplitudeMomentum(0),1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudelnuqqbar::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbarLeftCurrent(2,hel[2],3,hel[3]); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = current*wVertices*wPropergator;
  largeN = res;
  return res;
}

Complex MatchboxAmplitudelnuqqbar::evaluateOneLoop(size_t, const vector<int>& hel) {
  if ( abs(hel[2]+hel[3]) != 2 ) return 0.;
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbarLeftOneLoopCurrent(2,hel[2],3,hel[3]); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = (SM().alphaS()/(2.*Constants::pi))*current*wVertices*wPropergator;
  return res;
}

void MatchboxAmplitudelnuqqbar::persistentOutput(PersistentOStream & os) const {
  os <<  theDiagonal << theCKM ;
}

void MatchboxAmplitudelnuqqbar::persistentInput(PersistentIStream & is, int) {
  is >>  theDiagonal >> theCKM ;
}

DescribeClass<MatchboxAmplitudelnuqqbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbar("Herwig::MatchboxAmplitudelnuqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbar::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbar> documentation
    ("MatchboxAmplitudelnuqqbar");
  static Switch<MatchboxAmplitudelnuqqbar,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbar::theDiagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);
}
#line 1 "./MatchboxAmplitudelnuqqbarg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbarg class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbarg.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbarg::MatchboxAmplitudelnuqqbarg() 
  : theDiagonal(false) {}

void MatchboxAmplitudelnuqqbarg::doinit() {
  MatchboxAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(5);
}

void MatchboxAmplitudelnuqqbarg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(5);
}

IBPtr MatchboxAmplitudelnuqqbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbarg::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbarg::canHandle(const PDVector& proc) const {
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator elektron = xproc.begin();
  for ( ; elektron != xproc.end(); ++elektron ) 
    if ( abs((*elektron)->id()) >= 11 && abs((*elektron)->id()) <= 16 && abs((*elektron)->id()) % 2 == 1 ) break;
  if ( elektron == xproc.end() ) return false;
  PDPtr e = *elektron;
  xproc.erase(elektron);
  PDVector::iterator neutrino = xproc.begin();
  for ( ; neutrino != xproc.end(); ++neutrino ) 
    if ( abs((*neutrino)->id()) >= 11 && abs((*neutrino)->id()) <= 16 && abs((*neutrino)->id()) % 2 == 0 ) break;
  if ( neutrino == xproc.end() ) return false;
  PDPtr n = *neutrino;
  xproc.erase(neutrino);
  PDVector::iterator quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 1 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr d = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 0 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr u = *quark;
  xproc.erase(quark);
  PDVector::iterator gluon = xproc.begin();
  for ( ; gluon != xproc.end(); ++gluon ) 
    if ( (*gluon)->id() == 21 ) break; 
  if ( gluon == xproc.end() ) return false;
  xproc.erase(gluon);
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  if ( u->iCharge() + d->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
  if ( theDiagonal && SU2Helper::family(u) != SU2Helper::family(d) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  setupLeptons(0,amplitudeMomentum(0),1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudelnuqqbarg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  Complex sVertex =
          sqrt(4.*Constants::pi*SM().alphaS());
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbargLeftCurrent(2,hel[2],3,hel[3],4,hel[4]); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = current*wVertices*wPropergator*sVertex;
  largeN = res;
  return res;
}

Complex MatchboxAmplitudelnuqqbarg::evaluateOneLoop(size_t, const vector<int>& hel) {
  if ( abs(hel[2]+hel[3]) != 2 ) return 0.;
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  Complex sVertex =
          sqrt(4.*Constants::pi*SM().alphaS());
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbargLeftOneLoopCurrent(2,hel[2],3,hel[3],4,hel[4]); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = (SM().alphaS()/(2.*Constants::pi))*current*wVertices*wPropergator*sVertex;
  return res;
}

void MatchboxAmplitudelnuqqbarg::persistentOutput(PersistentOStream & os) const {
  os << theDiagonal << theCKM ;
}

void MatchboxAmplitudelnuqqbarg::persistentInput(PersistentIStream & is, int) {
  is >> theDiagonal >> theCKM ;
}

DescribeClass<MatchboxAmplitudelnuqqbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbarg("Herwig::MatchboxAmplitudelnuqqbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbarg::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbarg> documentation
    ("MatchboxAmplitudelnuqqbarg");
  static Switch<MatchboxAmplitudelnuqqbarg,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbarg::theDiagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);
}

#line 1 "./MatchboxAmplitudelnuqqbargg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbargg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbargg class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbargg.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbargg::MatchboxAmplitudelnuqqbargg() 
  : theDiagonal(false) {}

void MatchboxAmplitudelnuqqbargg::doinit() {
  MatchboxAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(6);
}

void MatchboxAmplitudelnuqqbargg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(6);
}

IBPtr MatchboxAmplitudelnuqqbargg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbargg::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbargg::canHandle(const PDVector& proc) const {
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator elektron = xproc.begin();
  for ( ; elektron != xproc.end(); ++elektron ) 
    if ( abs((*elektron)->id()) >= 11 && abs((*elektron)->id()) <= 16 && abs((*elektron)->id()) % 2 == 1 ) break;
  if ( elektron == xproc.end() ) return false;
  PDPtr e = *elektron;
  xproc.erase(elektron);
  PDVector::iterator neutrino = xproc.begin();
  for ( ; neutrino != xproc.end(); ++neutrino ) 
    if ( abs((*neutrino)->id()) >= 11 && abs((*neutrino)->id()) <= 16 && abs((*neutrino)->id()) % 2 == 0 ) break;
  if ( neutrino == xproc.end() ) return false;
  PDPtr n = *neutrino;
  xproc.erase(neutrino);
  PDVector::iterator quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 1 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr d = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 && abs((*quark)->id()) % 2 == 0 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr u = *quark;
  xproc.erase(quark);
  PDVector::iterator gluon = xproc.begin();
  for ( ; gluon != xproc.end(); ++gluon ) 
    if ( (*gluon)->id() == 21 ) break; 
  if ( gluon == xproc.end() ) return false;
  xproc.erase(gluon);
  gluon = xproc.begin();
  for ( ; gluon != xproc.end(); ++gluon ) 
    if ( (*gluon)->id() == 21 ) break; 
  if ( gluon == xproc.end() ) return false;
  xproc.erase(gluon);
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  if ( u->iCharge() + d->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
  if ( theDiagonal && SU2Helper::family(u) != SU2Helper::family(d) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbargg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  setupLeptons(0,amplitudeMomentum(0),1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  momentum(5,amplitudeMomentum(5));
  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudelnuqqbargg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {
  if ( abs(hel[2]+hel[3]) != 2 ) {
    largeN = 0.;
    return 0.;
  }
  assert ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 1 );
  int g1,hg1,g2,hg2;
  if ( amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else if ( a == 1 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else assert ( false );
  } else if ( amplitudeToColourMap()[4] == 3 && amplitudeToColourMap()[5] == 2 ) {
    if ( a == 0 ) {
      g1 = 5; hg1 = hel[5];
      g2 = 4; hg2 = hel[4];
    } else if ( a == 1 ) {
      g1 = 4; hg1 = hel[4];
      g2 = 5; hg2 = hel[5];
    } else assert ( false );
  } else assert ( false );
  Complex ckmelement = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp(
      SU2Helper::family(amplitudePartonData()[2]),
      SU2Helper::family(amplitudePartonData()[3]));
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp.first,tmp.second);
    ckmelement = theCKM[tmp.first][tmp.second];
    if ( !wPlus ) ckmelement = conj(ckmelement);
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement;
  Complex sVertices =
          4.*Constants::pi*SM().alphaS();
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  const LorentzVector<Complex>& quarkCurrent = qqbarggLeftCurrent(2,hel[2],3,hel[3],g1,hg1,g2,hg2); 
  Complex current = hel[2] == 1 ? Complex(0.,-1)*leptonCurrent.dot(quarkCurrent): 0.;
  Complex res = current*wVertices*wPropergator*sVertices;
  largeN = res;
  return res;  
}

void MatchboxAmplitudelnuqqbargg::persistentOutput(PersistentOStream & os) const {
  os << theDiagonal << theCKM ;
}

void MatchboxAmplitudelnuqqbargg::persistentInput(PersistentIStream & is, int) {
  is >> theDiagonal >> theCKM ;
}

DescribeClass<MatchboxAmplitudelnuqqbargg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbargg("Herwig::MatchboxAmplitudelnuqqbargg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbargg::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbargg> documentation
    ("MatchboxAmplitudelnuqqbargg");
  static Switch<MatchboxAmplitudelnuqqbargg,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbargg::theDiagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);
}

#line 1 "./MatchboxAmplitudelnuqqbarqqbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudelnuqqbarqqbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudelnuqqbarqqbar class.
//

#include "Herwig/MatrixElement/Matchbox/Builtin/Amplitudes/MatchboxAmplitudelnuqqbarqqbar.h"

using namespace Herwig;

MatchboxAmplitudelnuqqbarqqbar::MatchboxAmplitudelnuqqbarqqbar() 
  : theDiagonal(false) {}

void MatchboxAmplitudelnuqqbarqqbar::doinit() {
  MatchboxAmplitude::doinit();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  theCKM = standardCKM(SM())->getUnsquaredMatrix(6);
  nPoints(6);
}

void MatchboxAmplitudelnuqqbarqqbar::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MZ = getParticleData(ParticleID::Z0)->hardProcessMass();
  GZ = getParticleData(ParticleID::Z0)->hardProcessWidth();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  GW = getParticleData(ParticleID::Wplus)->hardProcessWidth();
  CA = SM().Nc();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(6);
}

IBPtr MatchboxAmplitudelnuqqbarqqbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudelnuqqbarqqbar::fullclone() const {
  return new_ptr(*this);
}

bool MatchboxAmplitudelnuqqbarqqbar::canHandle(const PDVector& proc) const {
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  //  Charge charge(ZERO);
  PDVector::iterator elektron = xproc.begin();
  for ( ; elektron != xproc.end(); ++elektron ) 
    if ( abs((*elektron)->id()) >= 11 && abs((*elektron)->id()) <= 16 && abs((*elektron)->id()) % 2 == 1 ) break;
  if ( elektron == xproc.end() ) return false;
  PDPtr e = *elektron;
  xproc.erase(elektron);
  PDVector::iterator neutrino = xproc.begin();
  for ( ; neutrino != xproc.end(); ++neutrino ) 
    if ( abs((*neutrino)->id()) >= 11 && abs((*neutrino)->id()) <= 16 && abs((*neutrino)->id()) % 2 == 0 ) break;
  if ( neutrino == xproc.end() ) return false;
  PDPtr n = *neutrino;
  xproc.erase(neutrino);
  PDVector::iterator quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q1 = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q2 = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q3 = *quark;
  xproc.erase(quark);
  quark = xproc.begin();
  for ( ; quark != xproc.end(); ++quark ) 
    if ( abs((*quark)->id()) >= 1 && abs((*quark)->id()) <= 6 ) {
      assert( (*quark)->hardProcessMass() == ZERO );
      break;
    } 
  if ( quark == xproc.end() ) return false;
  PDPtr q4 = *quark;
  xproc.erase(quark);
  if ( q1->id() == -q2->id() ) {
    if ( q3->iCharge() + q4->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q3) != SU2Helper::family(q4) ) return false;
  } else if ( q1->id() == -q3->id() ) {
    if ( q2->iCharge() + q4->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q2) != SU2Helper::family(q4) ) return false;
  } else  if ( q1->id() == -q4->id() ) {
    if ( q2->iCharge() + q3->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q2) != SU2Helper::family(q3) ) return false;
  } else  if ( q2->id() == -q3->id() ) {
    if ( q1->iCharge() + q4->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q1) != SU2Helper::family(q4) ) return false;
  } else  if ( q2->id() == -q4->id() ) {
    if ( q1->iCharge() + q3->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q1) != SU2Helper::family(q3) ) return false;
  } else  if ( q3->id() == -q4->id() ) {
    if ( q1->iCharge() + q2->iCharge() + e->iCharge() != PDT::ChargeNeutral ) return false;
    if ( theDiagonal && SU2Helper::family(q1) != SU2Helper::family(q2) ) return false;
  } else return false;
  if ( SU2Helper::family(e) != SU2Helper::family(n) ) return false;
  return xproc.empty();
}

void MatchboxAmplitudelnuqqbarqqbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  setupLeptons(0,amplitudeMomentum(0),1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));
  momentum(5,amplitudeMomentum(5));
  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudelnuqqbarqqbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {
  const LorentzVector<Complex>& leptonCurrent = llbarLeftCurrent(0,hel[0],1,hel[1]); 
  Complex Current2345 =
    hel[2] == 1 && hel[3] == 1 && abs(hel[4]+hel[5]) == 2 &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    leptonCurrent.dot(qqbarqqbarLeftCurrent(2,hel[2],3,hel[3],4,hel[4],5,hel[5])) : 0.;
  Complex Current4523 =
    hel[4] == 1 && hel[5] == 1 && abs(hel[2]+hel[3]) == 2 &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    leptonCurrent.dot(qqbarqqbarLeftCurrent(4,hel[4],5,hel[5],2,hel[2],3,hel[3])) : 0.;
  Complex Current2543 =
    hel[2] == 1 && hel[5] == 1 && abs(hel[4]+hel[3]) == 2 &&
    abs(amplitudePartonData()[4]->id()) == abs(amplitudePartonData()[3]->id()) ? 
    -leptonCurrent.dot(qqbarqqbarLeftCurrent(2,hel[2],5,hel[5],4,hel[4],3,hel[3])) : 0.;
  Complex Current4325 =
    hel[4] == 1 && hel[3] == 1 && abs(hel[2]+hel[5]) == 2 &&
    abs(amplitudePartonData()[2]->id()) == abs(amplitudePartonData()[5]->id()) ? 
    -leptonCurrent.dot(qqbarqqbarLeftCurrent(4,hel[4],3,hel[3],2,hel[2],5,hel[5])) : 0.;
  Complex ckmelement23 = 1.;
  Complex ckmelement45 = 1.;
  if ( !theDiagonal ) {
    bool wPlus = ( abs(amplitudePartonData()[0]->id()) % 2 == 0 ) ?
      amplitudePartonData()[1]->id() < 0:
      amplitudePartonData()[0]->id() < 0;
    pair<int,int> tmp23(
      SU2Helper::family(amplitudePartonData()[2])-1,
      SU2Helper::family(amplitudePartonData()[3])-1);
    pair<int,int> tmp45(
      SU2Helper::family(amplitudePartonData()[4])-1,
      SU2Helper::family(amplitudePartonData()[5])-1);
    if ( amplitudePartonData()[3]->id() < 0 ) swap(tmp23.first,tmp23.second);
    if ( amplitudePartonData()[5]->id() < 0 ) swap(tmp45.first,tmp45.second);
    ckmelement23 = theCKM[tmp23.first][tmp23.second];
    ckmelement45 = theCKM[tmp45.first][tmp45.second];
    if ( !wPlus ) {
      ckmelement23 = conj(ckmelement23);
      ckmelement45 = conj(ckmelement45);
    }
  }
  Complex wPropergator =
          1./Complex(((amplitudeMomentum(0)+amplitudeMomentum(1)).m2()-sqr(MW))/lastSHat(),MW*GW/lastSHat());
  Complex wVertices23 = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement23;
  Complex wVertices45 = 
          2.*SM().alphaEMMZ()*Constants::pi/SM().sin2ThetaW()*ckmelement45;
  Complex sVertices =
          4.*Constants::pi*SM().alphaS();
  Complex res2345 =
    Complex(0.,-1.)*wPropergator*sVertices*Current2345*wVertices23;
  Complex res2543 =
    Complex(0.,-1.)*wPropergator*sVertices*Current2543*wVertices23;
  Complex res4523 =
    Complex(0.,-1.)*wPropergator*sVertices*Current4523*wVertices45;
  Complex res4325 =
    Complex(0.,-1.)*wPropergator*sVertices*Current4325*wVertices45;
  double Nc = SM().Nc();
  Complex resLeading = 0.;
  Complex resSubLeading = 0.;
  if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 1 &&
       amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else if ( a == 1 ) {  //(25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 0 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 2 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { // (25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else if ( a == 1 ) { // (23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 3 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 1 ) {
    if ( a == 0 ) { //(23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else if ( a == 1 ) { //(25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else assert(false);
  } else if ( amplitudeToColourMap()[2] == 2 && amplitudeToColourMap()[3] == 1 &&
	      amplitudeToColourMap()[4] == 0 && amplitudeToColourMap()[5] == 3 ) {
    if ( a == 0 ) { //(25)(43)
      resLeading = res2345 + res4523;
      resSubLeading = res2543 + res4325;
    } else if ( a == 1 ) { //(23)(45)
      resLeading = res2543 + res4325;
      resSubLeading = res2345 + res4523;
    } else assert(false);
  } else assert(false);
  resSubLeading *= -1./Nc;
  largeN = resLeading/2.;
  return (resLeading + resSubLeading)/2.;
}

void MatchboxAmplitudelnuqqbarqqbar::persistentOutput(PersistentOStream & os) const {
  os << theDiagonal << theCKM ;
}

void MatchboxAmplitudelnuqqbarqqbar::persistentInput(PersistentIStream & is, int) {
  is >> theDiagonal >> theCKM ;
}

DescribeClass<MatchboxAmplitudelnuqqbarqqbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudelnuqqbarqqbar("Herwig::MatchboxAmplitudelnuqqbarqqbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudelnuqqbarqqbar::Init() {
  static ClassDocumentation<MatchboxAmplitudelnuqqbarqqbar> documentation
    ("MatchboxAmplitudelnuqqbarqqbar");
  static Switch<MatchboxAmplitudelnuqqbarqqbar,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &MatchboxAmplitudelnuqqbarqqbar::theDiagonal, false, false, false);
  static SwitchOption interfaceDiagonalYes
    (interfaceDiagonal,
     "Yes",
     "Use a diagonal CKM matrix.",
     true);
  static SwitchOption interfaceDiagonalNo
    (interfaceDiagonal,
     "No",
     "Use the CKM object as used by the StandardModel.",
     false);
}

#line 1 "./MatchboxAmplitudehbbbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudehbbbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehbbbar class.
//

#include "MatchboxAmplitudehbbbar.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudehbbbar::MatchboxAmplitudehbbbar() {}

IBPtr MatchboxAmplitudehbbbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehbbbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehbbbar::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(3);;
}

void MatchboxAmplitudehbbbar::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  CF = (SM().Nc()*SM().Nc()-1.)/(2.*SM().Nc());
  nPoints(3);
}

bool MatchboxAmplitudehbbbar::canHandle(const PDVector& proc) const {
  if ( proc.size() != 3 )
    return false;
  PDVector xproc = proc;
  PDVector::iterator q=xproc.begin();
  for (; q!=xproc.end(); ++q){
      if ((**q).id()==1 || 
          (**q).id()==2 ||
          (**q).id()==3 ||
          (**q).id()==4 ||
          (**q).id()==5) 
          break;
  }   
  if(q==xproc.end()) return false;
  int qid = (**q).id();
  xproc.erase(q);
  PDVector::iterator qb=xproc.begin();
  for (; qb!=xproc.end(); ++qb){
      if ((**qb).id()==-qid) break;
  }
  if(qb==xproc.end()) return false;
  xproc.erase(qb);
  PDVector::iterator h0=xproc.begin();
  if (xproc.size()==1 && (**h0).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehbbbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  //cout<<"prepare erreicht"<<endl;
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));
  //cout<<"momenta werden gesetzt"<<endl;
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehbbbar::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  unsigned int q=0;
  unsigned int qbar=0;
  cPDVector x=amplitudePartonData();   
  /*for (int i=0;i<3;++i){
    cout<<"x["<<i<<"]: "<<x[i]->id()<<endl;
  }*/
  //cout<<"teilchenidentifizierung";
  for (;q<amplitudePartonData().size();++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} //cout<<"x[q]"<<x[q]->id()<<" hel: "<<hel[q]<<endl;
  for (;qbar<amplitudePartonData().size();++qbar){if (x[qbar]->id() ==-x[q]->id()) break;} //cout<<"x[qbar]"<<x[qbar]->id()<<" hel: "<<hel[qbar]<<endl;
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  
  long id=x[q]->id();
  Energy Mf = 0*GeV;
  switch(id){
    case 1 : Mf = interfaceDMass; /*cout<<"d"<<ounit(Mf,GeV)<<endl;*/ break;
    case 2 : Mf = interfaceUMass; /*cout<<"u"<<ounit(Mf,GeV)<<endl;*/ break;
    case 3 : Mf = interfaceSMass; /*cout<<"s"<<ounit(Mf,GeV)<<endl;*/ break;
    case 4 : Mf = interfaceCMass; /*cout<<"c"<<ounit(Mf,GeV)<<endl;*/ break;
    case 5 : Mf = interfaceBMass; /*cout<<"b"<<ounit(Mf,GeV)<<endl;*/ break;
  }
  if (Mf==0*GeV) 
    throw Exception() << "Invalid settings in MatchboxAmplitudehbbbar -- zero fermion mass."
		      << Exception::runerror; 
  double im=gw*Mf/2/MW;
  Complex c = Complex(0.,im);
  //cout<<"c: "<<c<<endl; 
  //cout<<"largeN wird berechnet werden";
  //cout<<"plusproduct"<<plusProduct(qbar,q)<<endl;
  //cout<<"minusproduct"<<minusProduct(qbar,q)<<endl;
  if (hel[qbar]==-hel[q]&& hel[q]==1){
    largeN = 0;
    //cout<<"largeN plus hat geklappt"<<largeN;
    return(largeN);
  }
  if (hel[qbar]==-hel[q]&& hel[q]==-1){
    largeN = 0;
    //cout<<"largeN plus hat geklappt"<<largeN;
    return(largeN);
  }
  if (hel[qbar]==hel[q] && hel[q]==1){
    largeN = c*(plusProduct(qbar,q));
    //cout<<"largeN plus hat geklappt"<<largeN;
    return(largeN);
  }
  if (hel[qbar]==hel[q] && hel[q]==-1){
    largeN = c*(minusProduct(qbar,q));
    //cout<<"largeN minus hat geklappt"<<largeN;
    return(largeN);
  }
  assert(false);
  return 0;
}

Complex MatchboxAmplitudehbbbar::evaluateOneLoop(size_t, const vector<int>& hel) {
  Complex res = 0;
  int q=0;
  int qbar=0;
  cPDVector x=amplitudePartonData();   
  
  
  for (;q<3;++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} 
  for (;qbar<3;++qbar){if (x[qbar]->id() ==-x[q]->id()) break;} 
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  
  long id=x[q]->id();
  Energy Mf = 0*GeV;
  switch(id){
    case 1 : Mf = interfaceDMass; break;
    case 2 : Mf = interfaceUMass; break;
    case 3 : Mf = interfaceSMass; break;
    case 4 : Mf = interfaceCMass; break;
    case 5 : Mf = interfaceBMass; break;
  }
  if (Mf==0*GeV) 
    throw Exception() << "Invalid settings in MatchboxAmplitudehbbbar -- zero fermion mass."
		      << Exception::runerror; 
  
  double loop = SM().alphaS()*CF/2/Constants::pi ; //one-loop-Factor
  double bornim = gw*Mf/2/MW; //constant factor from born
  double real = loop*bornim;
  Complex c = Complex(real,0.);
  
  if (hel[qbar]==-hel[q]&& hel[q]==1){
    res = 0;
    return(res);
  }
  if (hel[qbar]==-hel[q]&& hel[q]==-1){
    res = 0;
    return(res);
  }
  if (hel[qbar]==hel[q] && hel[q]==1){
    res = c*(plusProduct(qbar,q));
    return(res);
  }
  if (hel[qbar]==hel[q] && hel[q]==-1){
    res = c*(minusProduct(qbar,q));
    return(res);
  }
  assert(false);
  return 0;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehbbbar::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceUMass,GeV)<< ounit(interfaceDMass,GeV)<< ounit(interfaceSMass,GeV)<< ounit(interfaceCMass,GeV)<< ounit(interfaceBMass,GeV);
}
void MatchboxAmplitudehbbbar::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceUMass,GeV)>> iunit(interfaceDMass,GeV)>> iunit(interfaceSMass,GeV)>> iunit(interfaceCMass,GeV)>> iunit(interfaceBMass,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehbbbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehbbbar("Herwig::MatchboxAmplitudehbbbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehbbbar::Init() {

  static ClassDocumentation<MatchboxAmplitudehbbbar> documentation
    ("MatchboxAmplitudehbbbar");
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceUMass
    ("interfaceUMass",
     "The up quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceUMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceDMass
    ("interfaceDMass",
     "The down quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceDMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceSMass
    ("interfaceSMass",
     "The strange quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceSMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceCMass
    ("interfaceCMass",
     "The charm quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceCMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbar,Energy> interfaceBMass
    ("interfaceBMass",
     "The bottom quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbar::interfaceBMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  

}

#line 1 "./MatchboxAmplitudehbbbarg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudehbbbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehbbbarg class.
//

#include "MatchboxAmplitudehbbbarg.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudehbbbarg::MatchboxAmplitudehbbbarg() {}

IBPtr MatchboxAmplitudehbbbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehbbbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehbbbarg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);;
}

void MatchboxAmplitudehbbbarg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

bool MatchboxAmplitudehbbbarg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  /*for (PDVector::iterator x=xproc.begin(); x!=xproc.end();  ++x){
    cout<<"xproc: "<<(**x).id()<<endl;
  }
  cout<<endl;*/
  PDVector::iterator q=xproc.begin();
  for (; q!=xproc.end(); ++q){
      if ((**q).id()==1 || 
          (**q).id()==2 ||
          (**q).id()==3 ||
          (**q).id()==4 ||
          (**q).id()==5) 
          break;
      //cout<<"kein Quark gefunden";
  }
  if(q==xproc.end()) return false;
  int qid = (**q).id(); //cout<<"qid ="<<qid<<endl;
  xproc.erase(q);
  PDVector::iterator qb=xproc.begin();
  for (; qb!=xproc.end(); ++qb){
      if ((**qb).id()==-qid /*|| (**qb).id()==qid*/ ){/*cout<<"qbid = "<<(**qb).id()<<endl; cout<<"kein antiquark gefunden!"<<endl;*/ break;}
  }
  if(qb==xproc.end()){/*cout<<"kein antiquark gefunden"<<endl;*/ return false;}
  xproc.erase(qb);
  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21){xproc.erase(g); break; }
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehbbbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  /*setupLeptons(0,amplitudeMomentum(0),
	       1,amplitudeMomentum(1));*/

  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehbbbarg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  unsigned int q=0;
  unsigned int qbar=0;
  unsigned int g=0;
  cPDVector x=amplitudePartonData(); 
  for (;q<amplitudePartonData().size();++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} 
  for (;qbar<amplitudePartonData().size();++qbar){if (x[qbar]->id() ==-x[q]->id()) break;}
  for (;g<amplitudePartonData().size();++g){if (x[g]->id() == 21) break;}
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double gs = sqrt(4*Constants::pi*SM().alphaS());
  //double alphaS = SM().alphaS();
  //double v= 2*MW/gw/sqrt(mu2()) ; //Auf mu normierter VakuumErwartungswert
  //double c = gs*alphaS/3/Constants::pi/v;
  long id=x[q]->id();
  Energy Mf = 0*GeV;
  switch(id){
    case 1 : Mf = interfaceDMass; /*cout<<"d"<<ounit(Mf,GeV)<<endl;*/ break;
    case 2 : Mf = interfaceUMass; /*cout<<"u"<<ounit(Mf,GeV)<<endl;*/ break;
    case 3 : Mf = interfaceSMass; /*cout<<"s"<<ounit(Mf,GeV)<<endl;*/ break;
    case 4 : Mf = interfaceCMass; /*cout<<"c"<<ounit(Mf,GeV)<<endl;*/ break;
    case 5 : Mf = interfaceBMass; /*cout<<"b"<<ounit(Mf,GeV)<<endl;*/ break;
  }
  if (Mf==0*GeV) 
    throw Exception() << "Invalid settings in MatchboxAmplitudehbbbarg -- zero fermion mass."
		      << Exception::runerror; 
  Complex c2= Complex(0.,gw*gs*Mf/MW/sqrt(2));
  
  //qghq Prozesse mit Quark als "Zwischenteilchen" aus qghq.nb
  if(hel[qbar]==+1 && hel[q]==+1 && hel[g]==-1){
    largeN = -c2*minusProduct(q,qbar)*minusProduct(q,qbar)/minusProduct(q,g)/minusProduct(qbar,g);
    return(largeN);
  }
  if(hel[qbar]==+1 && hel[q]==+1 && hel[g]==+1){
    largeN = c2*(minusProduct(qbar,g)/plusProduct(q,g)+
                 minusProduct(q,g)/plusProduct(qbar,g)-
                 invariant(qbar,q)/plusProduct(q,g)/plusProduct(qbar,g));
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==-1 && hel[g]==+1){
    largeN = c2*plusProduct(q,qbar)*plusProduct(q,qbar)/plusProduct(q,g)/plusProduct(qbar,g);
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==-1 && hel[g]==-1){
    largeN = -c2*(plusProduct(qbar,g)/minusProduct(q,g)+
                  plusProduct(q,g)/minusProduct(qbar,g)-
                  invariant(qbar,q)/minusProduct(q,g)/minusProduct(qbar,g));
    return(largeN);
  }
  if(hel[qbar]==-hel[q]){
    largeN = 0;
    return(largeN);
  }
  
//  Viel zu hohe Ordnnung in alphaS!!! Vielleicht spaeter per Interface dazuschaltbar??  
//  //mpm verdrehte vorzeichen wegen simons spinordef q=p qbar=q; g=r gluonhelizitaeten bleiben bleiben
//  if(hel[qbar]==+1 && hel[q]==-1 && hel[g]==-1){
//      largeN = c*minusProduct(q,qbar)*plusProduct(q,g)*plusProduct(q,g)/(invariant(q,qbar));
//      return(largeN);
//  }
//  //mpp
//  if(hel[qbar]==1 && hel[q]==-1 && hel[g]==1){
//      largeN = c*minusProduct(qbar,g)*minusProduct(qbar,g)*plusProduct(q,qbar)/(invariant(q,qbar));
//      return(largeN);
//  }
//  //pmm
//  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==-1){
//      largeN = -c*minusProduct(q,qbar)*plusProduct(qbar,g)*plusProduct(qbar,g)/(invariant(q,qbar));
//      return(largeN);
//  }
//  //pmp
//  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==1){
//      largeN = -c*minusProduct(q,g)*minusProduct(q,g)*plusProduct(q,qbar)/(invariant(q,qbar));
//      return(largeN);
//  }
  assert(false);

  return 0;
}

/*Complex MatchboxAmplitudehbbbarg::evaluateOneLoop(size_t, const vector<int>& hel) {

}*/

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehbbbarg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceUMass,GeV)<< ounit(interfaceDMass,GeV)<< ounit(interfaceSMass,GeV)<< ounit(interfaceCMass,GeV)<< ounit(interfaceBMass,GeV);
}
void MatchboxAmplitudehbbbarg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceUMass,GeV)>> iunit(interfaceDMass,GeV)>> iunit(interfaceSMass,GeV)>> iunit(interfaceCMass,GeV)>> iunit(interfaceBMass,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehbbbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehbbbarg("Herwig::MatchboxAmplitudehbbbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehbbbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudehbbbarg> documentation
    ("MatchboxAmplitudehbbbarg");
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceUMass
    ("interfaceUMass",
     "The up quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceUMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceDMass
    ("interfaceDMass",
     "The down quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceDMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceSMass
    ("interfaceSMass",
     "The strange quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceSMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceCMass
    ("interfaceCMass",
     "The charm quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceCMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  static Parameter<MatchboxAmplitudehbbbarg,Energy> interfaceBMass
    ("interfaceBMass",
     "The bottom quark mass to be used in the amplitude.",
     &MatchboxAmplitudehbbbarg::interfaceBMass, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  

}

#line 1 "./MatchboxAmplitudehgg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudehgg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehgg class.
//



#include "MatchboxAmplitudehgg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinorHelicity.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"


using namespace Herwig;

MatchboxAmplitudehgg::MatchboxAmplitudehgg()
  : interfaceTHooft(126*GeV) {}

IBPtr MatchboxAmplitudehgg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehgg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehgg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(3);
}

void MatchboxAmplitudehgg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(3);
}

bool MatchboxAmplitudehgg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 3 ) return false;
  PDVector xproc = proc;
  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
    if ((**g).id()==21) {xproc.erase(g); --g;}
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehgg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {
  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehgg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {

  unsigned int g1=0;
  unsigned int g2=0;
  cPDVector x=amplitudePartonData();
  for (;g1<amplitudePartonData().size();++g1){
    if (x[g1]->id()==21) {
      for (g2=g1+1;g2<amplitudePartonData().size();++g2){if (x[g2]->id()==21) break;}
      break;
    }
  }
  // wrong assignement of the particles. g1 and g2 have to be different gluons.
  assert(g1!=g2);
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double v= 2*MW/gw/sqrt(lastSHat()); 
  Complex c = Complex (0.,-1.*SM().alphaS()/3/Constants::pi/v); 
 
  if (hel[g1]==-hel[g2]){
    largeN=0;
    return largeN;
  }
  if (hel[g1]==hel[g2] && hel[g1]==-1){
    largeN=c*plusProduct(g1,g2)*plusProduct(g1,g2);
    return largeN;
  }
  if (hel[g1]==hel[g2] && hel[g1]==1){
    largeN=c*minusProduct(g1,g2)*minusProduct(g1,g2);
    return largeN;
  }
  // unknown helicity configuration
  assert(false);
  return(0.);
}

Complex MatchboxAmplitudehgg::evaluateOneLoop(size_t a , const vector<int>& Hel) {
  Complex E = SM().alphaS()/Constants::pi*11/4;
  Complex largeN = 0.;
  return (E*evaluate(a,Hel,largeN));   
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehgg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceTHooft,GeV);
}

void MatchboxAmplitudehgg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceTHooft,GeV);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehgg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehgg("Herwig::MatchboxAmplitudehgg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehgg::Init() {

  static ClassDocumentation<MatchboxAmplitudehgg> documentation
    ("MatchboxAmplitudehgg");

  /* // not used guess leftover from validation (mu2() variation)
  static Parameter<MatchboxAmplitudehgg,Energy> interfaceTHooft
    ("interfaceTHooft",
     "The THooft Mass.",
     &MatchboxAmplitudehgg::interfaceTHooft, GeV, 115.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  */

}

#line 1 "./MatchboxAmplitudehggg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudehggg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehggg class.
//

#include "MatchboxAmplitudehggg.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/MatrixElement/Matchbox/Utility/SpinorHelicity.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudehggg::MatchboxAmplitudehggg()
  : interfaceTHooft(126*GeV) {}

IBPtr MatchboxAmplitudehggg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehggg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehggg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

void MatchboxAmplitudehggg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

bool MatchboxAmplitudehggg::canHandle(const PDVector& proc) const {

  if ( proc.size() != 4 ) return false;
  PDVector xproc = proc;

  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21) {xproc.erase(g); --g;}
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehggg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));

  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehggg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {
  
  unsigned int p=0;
  unsigned int q=0;
  unsigned int r=0;
  cPDVector x=amplitudePartonData();
  for (;p<amplitudePartonData().size();++p){
    if (x[p]->id()==21) {
      for (q=p+1;q<amplitudePartonData().size();++q){
        if (x[q]->id()==21){
         for (r=q+1;r<amplitudePartonData().size();++r){if (x[r]->id()==21) break;}
          } 
        break;}
      break;
    }
  }
  // Wrong particle assignment. There have to be three distinct Gluons p, q and r.
  assert(!((p==q) || (p==r) || (q==r)));
 
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double gs = sqrt(4*Constants::pi*SM().alphaS());
  double v= 2*MW/gw/sqrt(lastSHat()) ;
  Complex c = Complex (0.,0.);                                                                                                                    
  
  // Assertion makes sure that there is no crossing, which causes a relative minus sign. If Assert -> new, more general code is needed                                                                                                                                    
  assert(amplitudeToColourMap()[1]==0 && amplitudeToColourMap()[2]==1 && amplitudeToColourMap()[3]==2 );

  if (a==0){
    c = Complex (0.,-1.)*sqrt(2.)*SM().alphaS()/3./Constants::pi/v*gs ; 
  } 
  else {
    if (a==1){
      c = Complex (0.,+1.)*sqrt(2.)*SM().alphaS()/3./Constants::pi/v*gs ; 
    }
    else{ 
      //The Colourbasis a is not appropriate for this process. 
      //  hggg ~ f^{abc} -> ~ tr(t^a,t^b,t^c) - tr(t^a,t^c,t^b) -> a in {0,1}
      assert(true);
    }
  }
  if(hel[p]==+1 && hel[q]==+1 && hel[r]==+1){
    largeN = c*(invariant(p,q)*invariant(p,q)+2*invariant(p,q)*invariant(p,r)
            +invariant(p,r)*invariant(p,r)+2*invariant(p,q)*invariant(q,r)
            +2*invariant(p,r)*invariant(q,r)+invariant(q,r)*invariant(q,r))
            /plusProduct(p,q)/plusProduct(p,r)/plusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==+1 && hel[q]==+1 && hel[r]==-1){
    largeN = -c*minusProduct(p,q)*minusProduct(p,q)*minusProduct(p,q)/minusProduct(p,r)/minusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==+1 && hel[q]==-1 && hel[r]==+1){
    largeN = -c*minusProduct(p,r)*minusProduct(p,r)*minusProduct(p,r)/minusProduct(p,q)/minusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==+1 && hel[q]==-1 && hel[r]==-1){
    largeN = c*plusProduct(q,r)*plusProduct(q,r)*plusProduct(q,r)/plusProduct(p,q)/plusProduct(p,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==+1 && hel[r]==+1){
    largeN = -c*minusProduct(q,r)*minusProduct(q,r)*minusProduct(q,r)/minusProduct(p,q)/minusProduct(p,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==+1 && hel[r]==-1){
    largeN = c*plusProduct(p,r)*plusProduct(p,r)*plusProduct(p,r)/plusProduct(p,q)/plusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==-1 && hel[r]==+1){
    largeN = c*plusProduct(p,q)*plusProduct(p,q)*plusProduct(p,q)/plusProduct(p,r)/plusProduct(q,r);
    return(largeN);
  }
  if(hel[p]==-1 && hel[q]==-1 && hel[r]==-1){
    largeN = -c*(invariant(p,q)*invariant(p,q)+2*invariant(p,q)*invariant(p,r)
            +invariant(p,r)*invariant(p,r)+2*invariant(p,q)*invariant(q,r)
            +2*invariant(p,r)*invariant(q,r)+invariant(q,r)*invariant(q,r))
            /plusProduct(p,q)/plusProduct(p,r)/plusProduct(q,r);
    return(largeN);
  }
  // Unknown helicity configuration
  assert(false);
  return(0.);
}

/*Complex MatchboxAmplitudehggg::evaluateOneLoop(size_t, const vector<int>& hel) {

}*/

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehggg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceTHooft,GeV);
}

void MatchboxAmplitudehggg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceTHooft,GeV);
}

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehggg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehggg("Herwig::MatchboxAmplitudehggg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehggg::Init() {

  static ClassDocumentation<MatchboxAmplitudehggg> documentation
    ("MatchboxAmplitudehggg");

  /*  // not used guess leftover from validation (mu2() variation)
  static Parameter<MatchboxAmplitudehggg,Energy> interfaceTHooft
    ("interfaceTHooft",
     "The THooft Mass.",
     &MatchboxAmplitudehggg::interfaceTHooft, GeV, 115.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  */
}

 
#line 1 "./MatchboxAmplitudehqqbarg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudehqqbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudehqqbarg class.
//

#include "MatchboxAmplitudehqqbarg.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudehqqbarg::MatchboxAmplitudehqqbarg() 
  : interfaceTHooft(126*GeV) {}

IBPtr MatchboxAmplitudehqqbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudehqqbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudehqqbarg::doinit() {
  MatchboxAmplitude::doinit();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);;
}

void MatchboxAmplitudehqqbarg::doinitrun() {
  MatchboxAmplitude::doinitrun();
  MW = getParticleData(ParticleID::Wplus)->hardProcessMass();
  nPoints(4);
}

bool MatchboxAmplitudehqqbarg::canHandle(const PDVector& proc) const {
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();
  PDVector::iterator q=xproc.begin();
  for (; q!=xproc.end(); ++q){
    if ((**q).id()==1 || 
        (**q).id()==2 ||
        (**q).id()==3 ||
        (**q).id()==4 ||
        (**q).id()==5) 
    break;
  }
  if(q==xproc.end()) return false;
  int qid = (**q).id(); 
  xproc.erase(q);
  PDVector::iterator qb=xproc.begin();
  for (; qb!=xproc.end(); ++qb){
    if ((**qb).id()==-qid ){break;}
  }
  if(qb==xproc.end()){ return false;}
  xproc.erase(qb);
  for (PDVector::iterator g=xproc.begin(); g!=xproc.end(); ++g){
      if ((**g).id()==21){xproc.erase(g); break; }
  }   
  if (xproc.size()==1 && (**xproc.begin()).id()==25) {return true;}
  return false;
}

void MatchboxAmplitudehqqbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }

  amplitudeScale(sqrt(lastSHat()));
 
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudehqqbarg::evaluate(size_t, const vector<int>& hel, Complex& largeN) {
  unsigned int q=0;
  unsigned int qbar=0;
  unsigned int g=0;
  cPDVector x=amplitudePartonData();
  for (;q<amplitudePartonData().size();++q){if (x[q]->id()!= 25 && x[q]->id()>0 ) break;} 
  for (;qbar<amplitudePartonData().size();++qbar){if (x[qbar]->id() ==-x[q]->id()) break;}
  for (;g<amplitudePartonData().size();++g){if (x[g]->id() == 21) break;}
  double gw = sqrt(4*Constants::pi*SM().alphaEMMZ()) / sqrt(SM().sin2ThetaW());
  double gs = sqrt(4*Constants::pi*SM().alphaS());
  double alphaS = SM().alphaS();
  double v= 2*MW/gw/sqrt(lastSHat()) ;
  double c = gs*alphaS/3/sqrt(2.)/Constants::pi/v;
  
  if(hel[qbar]==hel[q]){
    largeN = 0;
    return(largeN);
  }
  
  if(hel[qbar]==+1 && hel[q]==-1 && hel[g]==-1){
    largeN = -c*plusProduct(qbar,g)*plusProduct(qbar,g)/(plusProduct(q,qbar));
    return(largeN);
  }
  if(hel[qbar]==1 && hel[q]==-1 && hel[g]==1){
    largeN = -c*minusProduct(q,g)*minusProduct(q,g)/(minusProduct(q,qbar));
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==-1){
    largeN = c*plusProduct(q,g)*plusProduct(q,g)/(plusProduct(q,qbar));
    return(largeN);
  }
  if(hel[qbar]==-1 && hel[q]==1 && hel[g]==1){
    largeN = c*minusProduct(qbar,g)*minusProduct(qbar,g)/(minusProduct(q,qbar));
    return(largeN);
  }
  // Unknown helicity configuration
  assert(false);
  return(0.);
}

/*Complex MatchboxAmplitudehqqbarg::evaluateOneLoop(size_t, const vector<int>& hel) {

}*/

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudehqqbarg::persistentOutput(PersistentOStream &os) const {
  os << ounit(interfaceTHooft,GeV);
}

void MatchboxAmplitudehqqbarg::persistentInput(PersistentIStream &is, int) {
  is >> iunit(interfaceTHooft,GeV);
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudehqqbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudehqqbarg("Herwig::MatchboxAmplitudehqqbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudehqqbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudehqqbarg> documentation
    ("MatchboxAmplitudehqqbarg");

  /*  // not used guess leftover from validation (mu2() variation)
  static Parameter<MatchboxAmplitudehqqbarg,Energy> interfaceTHooft
    ("interfaceTHooft",
     "The THooft Mass.",
     &MatchboxAmplitudehqqbarg::interfaceTHooft, GeV, 115.0*GeV, 0.0*GeV, 0*GeV,
     false, false, Interface::lowerlim);
  */
}

#line 1 "./MatchboxAmplitudeqqbarttbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudeqqbarttbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudeqqbarttbar class.
//

#include "MatchboxAmplitudeqqbarttbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudeqqbarttbar::MatchboxAmplitudeqqbarttbar() {}

IBPtr MatchboxAmplitudeqqbarttbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudeqqbarttbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudeqqbarttbar::setupParams(map<string, double> &MGParams){
  // create parameter map for adapted Madgraph amplitude to use
  MGParams["aS"]     = SM().alphaS();
  MGParams["MZ"]     = getParticleData(ParticleID::Z0)      -> hardProcessMass()/GeV;
  MGParams["MW"]     = getParticleData(ParticleID::Wplus)   -> hardProcessMass()/GeV;
  MGParams["MH"]     = getParticleData(ParticleID::h0)      -> hardProcessMass()/GeV;
  MGParams["MT"]     = getParticleData(ParticleID::t)       -> hardProcessMass()/GeV;
  MGParams["MB"]     = getParticleData(ParticleID::b)       -> hardProcessMass()/GeV;
  MGParams["MTA"]    = getParticleData(ParticleID::tauminus)-> hardProcessMass()/GeV;
  MGParams["WW"]     = getParticleData(ParticleID::Wplus)   ->hardProcessWidth()/GeV;
  MGParams["WZ"]     = getParticleData(ParticleID::Z0)      ->hardProcessWidth()/GeV;
  MGParams["WH"]     = getParticleData(ParticleID::h0)      ->hardProcessWidth()/GeV;
  MGParams["WT"]     = getParticleData(ParticleID::t)       ->hardProcessWidth()/GeV;
  MGParams["WB"]     = getParticleData(ParticleID::b)       ->hardProcessWidth()/GeV;
  MGParams["GF"]     = SM().fermiConstant()*GeV2;
  MGParams["aEWM1"]  = 1./SM().alphaEMMZ();
  return;
}

void MatchboxAmplitudeqqbarttbar::doinit() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinit();
  nPoints(4);
}

void MatchboxAmplitudeqqbarttbar::doinitrun() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinitrun();
  nPoints(4);
}

bool MatchboxAmplitudeqqbarttbar::canHandle(const PDVector& proc) const {
  // check process is qqbar > ttbar
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  PDVector::iterator top = xproc.begin();
  long topId = 0;
  for ( ; top != xproc.end(); ++top )
    if ( (**top).id() == 6 ) {
      break;
    } 
  if ( top == xproc.end() )
    return false;
  topId = (**top).id();
  xproc.erase(top);
  PDVector::iterator antiTop = xproc.begin();
  for ( ; antiTop != xproc.end(); ++antiTop )
    if ( (**antiTop).id() == -topId ) {
      break;
    }
  if ( antiTop == xproc.end() )
    return false;
  xproc.erase(antiTop);
  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0  ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);
  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    }
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  return xproc.empty();
}


void MatchboxAmplitudeqqbarttbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  
  amplitudeScale(sqrt(lastSHat()));

  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudeqqbarttbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  // set up momenta to pass into Madgraph amplitude
  vector<double*> p;
  p.push_back(mom0);  p.push_back(mom1);  p.push_back(mom2);  p.push_back(mom3);
  for (int ip=0; ip<4; ++ip){
    p[ip][0] = abs(momentum(ip).e())<1.e-13 ? 0.0:double(momentum(ip).e()*amplitudeScale()/GeV);
    p[ip][1] = abs(momentum(ip).x())<1.e-13 ? 0.0:double(momentum(ip).x()*amplitudeScale()/GeV);
    p[ip][2] = abs(momentum(ip).y())<1.e-13 ? 0.0:double(momentum(ip).y()*amplitudeScale()/GeV);
    p[ip][3] = abs(momentum(ip).z())<1.e-13 ? 0.0:double(momentum(ip).z()*amplitudeScale()/GeV);
  } 

  // calculate amplitudes
  vector<complex<double> > amplitudes;
  MG_qqx2ttx process;
  process.initProc(MGParams_);  
  process.setMomenta(p);
  process.sigmaKin(amplitudes, hel);

  // calculate colour flows
  complex<double> matrixElement;
  if      (a==0) matrixElement = amplitudes[0]*( 1./6.);
  else if (a==1) matrixElement = amplitudes[0]*(-1./2.);
  else assert(false);

  largeN = matrixElement;
  return matrixElement;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudeqqbarttbar::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudeqqbarttbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudeqqbarttbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudeqqbarttbar("Herwig::MatchboxAmplitudeqqbarttbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudeqqbarttbar::Init() {

  static ClassDocumentation<MatchboxAmplitudeqqbarttbar> documentation
    ("MatchboxAmplitudeqqbarttbar");

}

#line 1 "./MatchboxAmplitudeqqbarttbarg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudeqqbarttbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudeqqbarttbarg class.
//

#include "MatchboxAmplitudeqqbarttbarg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudeqqbarttbarg::MatchboxAmplitudeqqbarttbarg() {}

IBPtr MatchboxAmplitudeqqbarttbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudeqqbarttbarg::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudeqqbarttbarg::setupParams(map<string, double> &MGParams){
  // create parameter map for adapted Madgraph amplitude to use
  MGParams["aS"]     = SM().alphaS();
  MGParams["MZ"]     = getParticleData(ParticleID::Z0)    -> hardProcessMass()/GeV;
  MGParams["MW"]     = getParticleData(ParticleID::Wplus) -> hardProcessMass()/GeV;
  MGParams["MH"]     = getParticleData(ParticleID::h0)    -> hardProcessMass()/GeV;
  MGParams["WW"]     = getParticleData(ParticleID::Wplus) ->hardProcessWidth()/GeV;
  MGParams["WZ"]     = getParticleData(ParticleID::Z0)    ->hardProcessWidth()/GeV;
  MGParams["WH"]     = getParticleData(ParticleID::h0)    ->hardProcessWidth()/GeV;
  MGParams["MT"]     = getParticleData(ParticleID::t)     -> hardProcessMass()/GeV;
  MGParams["MB"]     = getParticleData(ParticleID::b)     -> hardProcessMass()/GeV;
  MGParams["WT"]     = getParticleData(ParticleID::t)     ->hardProcessWidth()/GeV;
  MGParams["WB"]     = getParticleData(ParticleID::b)     ->hardProcessWidth()/GeV;
  MGParams["MTA"]    = getParticleData(ParticleID::tauminus)-> hardProcessMass()/GeV;
  MGParams["GF"]     = SM().fermiConstant()*GeV2;
  MGParams["aEWM1"]  = 1./SM().alphaEMMZ();
  return;
}


void MatchboxAmplitudeqqbarttbarg::doinit() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinit();
  nPoints(5);
}

void MatchboxAmplitudeqqbarttbarg::doinitrun() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinitrun();
  nPoints(5);
}

bool MatchboxAmplitudeqqbarttbarg::canHandle(const PDVector& proc) const {
  // check process is qqbar > ttbarg or some crossing of this
  if ( proc.size() != 5 )
    return false;
  PDVector xproc = proc;
  if ( xproc[0]->CC() )
    xproc[0] = xproc[0]->CC();
  if ( xproc[1]->CC() )
    xproc[1] = xproc[1]->CC();  
  PDVector::iterator top = xproc.begin();
  long topId = 0;
  for ( ; top != xproc.end(); ++top )
    if ( (**top).id() == 6 ) {
      break;
    }
  if ( top == xproc.end() )
    return false;
  topId = (**top).id();
  xproc.erase(top);
  PDVector::iterator antiTop = xproc.begin();
  for ( ; antiTop != xproc.end(); ++antiTop )
    if ( (**antiTop).id() == -topId ) {
      break;
    }
  if ( antiTop == xproc.end() )
    return false;
  xproc.erase(antiTop);
  PDVector::iterator quark = xproc.begin();
  long quarkId = 0;
  for ( ; quark != xproc.end(); ++quark )
    if ( abs((**quark).id()) < 6 &&
	 (**quark).id() > 0  ) {
      break;
    }
  if ( quark == xproc.end() )
    return false;
  quarkId = (**quark).id();
  xproc.erase(quark);
  PDVector::iterator antiQuark = xproc.begin();
  for ( ; antiQuark != xproc.end(); ++antiQuark )
    if ( (**antiQuark).id() == -quarkId ) {
      break;
    } 
  if ( antiQuark == xproc.end() )
    return false;
  xproc.erase(antiQuark);
  return xproc[0]->id() == 21;
}


void MatchboxAmplitudeqqbarttbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));

  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudeqqbarttbarg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  // check if doing qqbar, qg or qbarg initiated process
  int crossed;
  cPDVector partons = mePartonData();
  if (abs(partons[0]->id())< 6 && 
      partons[1]->id()==-partons[0]->id()) crossed=1;
  else if (((partons[0]->id()>0 && partons[0]->id()<6) && 
  	    partons[1]->id()==ParticleID::g) ||
  	   ((partons[1]->id()>0 && partons[1]->id()<6) && 
  	    partons[0]->id()==ParticleID::g)) crossed=2;
  else if (((partons[0]->id()<0 && partons[0]->id()>-6) && 
  	    partons[1]->id()==ParticleID::g) ||
  	   ((partons[1]->id()<0 && partons[1]->id()>-6) && 
  	    partons[0]->id()==ParticleID::g)) crossed=3;
  else 
    throw Exception() << "MatchboxAmplitudeqqbarststbarg::evaluate(): Unrecognised process\n"
  		      << Exception::runerror;


  // set up momenta to pass into Madgraph amplitude
  vector<double*> p;
  p.push_back(mom0);  p.push_back(mom1);  p.push_back(mom2);  p.push_back(mom3); p.push_back(mom4);
  for (int ip=0; ip<5; ++ip){
    p[ip][0] = abs(momentum(ip).e())<1.e-13 ? 0.0:double(momentum(ip).e()*amplitudeScale()/GeV);
    p[ip][1] = abs(momentum(ip).x())<1.e-13 ? 0.0:double(momentum(ip).x()*amplitudeScale()/GeV);
    p[ip][2] = abs(momentum(ip).y())<1.e-13 ? 0.0:double(momentum(ip).y()*amplitudeScale()/GeV);
    p[ip][3] = abs(momentum(ip).z())<1.e-13 ? 0.0:double(momentum(ip).z()*amplitudeScale()/GeV);
  } 
 
  // calculate amplitudes
  vector<complex<double> > amplitudes;
  MG_qqx2ttxg process;
  double factor=lastSHat()/GeV2;
  process.initProc(MGParams_);  
  process.setMomenta(p);
  process.sigmaKin(amplitudes, hel, crossed);

  for (int iamp=0; iamp<int(amplitudes.size()); ++iamp)
    amplitudes[iamp]*=sqrt(factor);

  complex<double> matrixElement;
  complex<double> i = complex<double>(0,1);

  // calculate colour flows
  if (a==0)
    matrixElement = (1./6.)*(amplitudes[1] + amplitudes[2]);
  else if (a==1)
    matrixElement =  1./2.*( i*amplitudes[0] - amplitudes[1] - amplitudes[3]);  
  else if (a==2) 
    matrixElement = (1./6.)*( amplitudes[3] + amplitudes[4]);
  else if (a==3)
    matrixElement =  1./2.*(-i*amplitudes[0] - amplitudes[2] - amplitudes[4]);
  else assert(false);

  largeN = matrixElement;
  return matrixElement;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudeqqbarttbarg::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudeqqbarttbarg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudeqqbarttbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudeqqbarttbarg("Herwig::MatchboxAmplitudeqqbarttbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudeqqbarttbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudeqqbarttbarg> documentation
    ("MatchboxAmplitudeqqbarttbarg");

}

#line 1 "./MatchboxAmplitudeggttbar.cc"
// -*- C++ -*-
//
// MatchboxAmplitudeggttbar.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudeggttbar class.
//

#include "MatchboxAmplitudeggttbar.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudeggttbar::MatchboxAmplitudeggttbar() {}

IBPtr MatchboxAmplitudeggttbar::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudeggttbar::fullclone() const {
  return new_ptr(*this);
}

void MatchboxAmplitudeggttbar::setupParams(map<string, double> &MGParams){
  // create parameter map for adapted Madgraph amplitude to use
  MGParams["aS"]     = SM().alphaS();
  MGParams["MZ"]     = getParticleData(ParticleID::Z0)    -> hardProcessMass()/GeV;
  MGParams["MW"]     = getParticleData(ParticleID::Wplus) -> hardProcessMass()/GeV;
  MGParams["MH"]     = getParticleData(ParticleID::h0)    -> hardProcessMass()/GeV;
  MGParams["WW"]     = getParticleData(ParticleID::Wplus) ->hardProcessWidth()/GeV;
  MGParams["WZ"]     = getParticleData(ParticleID::Z0)    ->hardProcessWidth()/GeV;
  MGParams["WH"]     = getParticleData(ParticleID::h0)    ->hardProcessWidth()/GeV;
  MGParams["MT"]     = getParticleData(ParticleID::t)     -> hardProcessMass()/GeV;
  MGParams["MB"]     = getParticleData(ParticleID::b)     -> hardProcessMass()/GeV;
  MGParams["WT"]     = getParticleData(ParticleID::t)     ->hardProcessWidth()/GeV;
  MGParams["WB"]     = getParticleData(ParticleID::b)     ->hardProcessWidth()/GeV;
  MGParams["MTA"]    = getParticleData(ParticleID::tauminus)-> hardProcessMass()/GeV;
  MGParams["GF"]     = SM().fermiConstant()*GeV2;
  MGParams["aEWM1"]  = 1./SM().alphaEMMZ();
  return;
}


void MatchboxAmplitudeggttbar::doinit() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinit();
  nPoints(4);
}

void MatchboxAmplitudeggttbar::doinitrun() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinitrun();
  nPoints(4);
}

bool MatchboxAmplitudeggttbar::canHandle(const PDVector& proc) const {
  // check process is gg > ttbar
  if ( proc.size() != 4 )
    return false;
  PDVector xproc = proc;
  PDVector::iterator top = xproc.begin();
  long topId = 0;
  for ( ; top != xproc.end(); ++top )
    if ( (**top).id() == 6 ) {
      break;
    }
  if ( top == xproc.end() )
    return false;
  topId = (**top).id();
  xproc.erase(top);
  PDVector::iterator antiTop = xproc.begin();
  for ( ; antiTop != xproc.end(); ++antiTop )
    if ( (**antiTop).id() == -topId ) {
      break;
    }
  if ( antiTop == xproc.end() )
    return false;
  xproc.erase(antiTop);
  return (xproc[0]->id()==21 && xproc[1]->id()==21);
}


void MatchboxAmplitudeggttbar::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  
  amplitudeScale(sqrt(lastSHat()));
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));

  MatchboxAmplitude::prepareAmplitudes(me);
}

Complex MatchboxAmplitudeggttbar::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  // set up momenta to pass into Madgraph amplitude
  vector<double*> p;
  p.push_back(mom0);  p.push_back(mom1);  p.push_back(mom2);  p.push_back(mom3);
  for (int ip=0; ip<4; ++ip){
    p[ip][0] = abs(momentum(ip).e())<1.e-13 ? 0.0:double(momentum(ip).e()*amplitudeScale()/GeV);
    p[ip][1] = abs(momentum(ip).x())<1.e-13 ? 0.0:double(momentum(ip).x()*amplitudeScale()/GeV);
    p[ip][2] = abs(momentum(ip).y())<1.e-13 ? 0.0:double(momentum(ip).y()*amplitudeScale()/GeV);
    p[ip][3] = abs(momentum(ip).z())<1.e-13 ? 0.0:double(momentum(ip).z()*amplitudeScale()/GeV);
  } 

  // calculate amplitudes
  vector<complex<double> > amplitudes;
  MG_gg2ttx process;
  process.initProc(MGParams_);
  process.setMomenta(p);
  process.sigmaKin(amplitudes, hel);

  complex<double> matrixElement;
  complex<double> i = complex<double> (0, 1);

  // calculate colour flows
  if      (a==0) matrixElement = 0.;
  else if (a==1) matrixElement = i*amplitudes[0] - amplitudes[1];
  else if (a==2) matrixElement =-i*amplitudes[0] - amplitudes[2];
  else assert(false); 

  largeN = matrixElement;
  return matrixElement;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudeggttbar::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudeggttbar::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudeggttbar,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudeggttbar("Herwig::MatchboxAmplitudeggttbar", "HwMatchboxBuiltin.so");

void MatchboxAmplitudeggttbar::Init() {

  static ClassDocumentation<MatchboxAmplitudeggttbar> documentation
    ("MatchboxAmplitudeggttbar");

}

#line 1 "./MatchboxAmplitudeggttbarg.cc"
// -*- C++ -*-
//
// MatchboxAmplitudeggttbarg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MatchboxAmplitudeggttbarg class.
//

#include "MatchboxAmplitudeggttbarg.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace Herwig;

MatchboxAmplitudeggttbarg::MatchboxAmplitudeggttbarg() {}

IBPtr MatchboxAmplitudeggttbarg::clone() const {
  return new_ptr(*this);
}

IBPtr MatchboxAmplitudeggttbarg::fullclone() const {
  return new_ptr(*this);
}


void MatchboxAmplitudeggttbarg::setupParams(map<string, double> &MGParams){
  // create parameter map for adapted Madgraph amplitude to use
  MGParams["aS"]     = SM().alphaS();
  MGParams["MZ"]     = getParticleData(ParticleID::Z0)    -> hardProcessMass()/GeV;
  MGParams["MW"]     = getParticleData(ParticleID::Wplus) -> hardProcessMass()/GeV;
  MGParams["MH"]     = getParticleData(ParticleID::h0)    -> hardProcessMass()/GeV;
  MGParams["WW"]     = getParticleData(ParticleID::Wplus) ->hardProcessWidth()/GeV;
  MGParams["WZ"]     = getParticleData(ParticleID::Z0)    ->hardProcessWidth()/GeV;
  MGParams["WH"]     = getParticleData(ParticleID::h0)    ->hardProcessWidth()/GeV;
  MGParams["MT"]     = getParticleData(ParticleID::t)     -> hardProcessMass()/GeV;
  MGParams["MB"]     = getParticleData(ParticleID::b)     -> hardProcessMass()/GeV;
  MGParams["WT"]     = getParticleData(ParticleID::t)     ->hardProcessWidth()/GeV;
  MGParams["WB"]     = getParticleData(ParticleID::b)     ->hardProcessWidth()/GeV;
  MGParams["MTA"]    = getParticleData(ParticleID::tauminus)-> hardProcessMass()/GeV;
  MGParams["GF"]     = SM().fermiConstant()*GeV2;
  MGParams["aEWM1"]  = 1./SM().alphaEMMZ();
  return;
}


void MatchboxAmplitudeggttbarg::doinit() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinit();
  nPoints(5);
}

void MatchboxAmplitudeggttbarg::doinitrun() {
  setupParams(MGParams_);
  MatchboxAmplitude::doinitrun();
  nPoints(5);
}

bool MatchboxAmplitudeggttbarg::canHandle(const PDVector& proc) const {
  // check process is gg > ttbarg
  if ( proc.size() != 5 )
    return false;
  PDVector xproc = proc;
  PDVector::iterator top = xproc.begin();
  long topId = 0;
  for ( ; top != xproc.end(); ++top )
    if ( (**top).id() == 6 ) {
      break;
    }
  if ( top == xproc.end() )
    return false;
  topId = (**top).id();
  xproc.erase(top);
  PDVector::iterator antiTop = xproc.begin();
  for ( ; antiTop != xproc.end(); ++antiTop )
    if ( (**antiTop).id() == -topId ) {
      break;
    }
  if ( antiTop == xproc.end() )
    return false;
  xproc.erase(antiTop);
  if ( xproc.size() != 3 )
    return false;
  return (xproc[0]->id() == 21 && xproc[1]->id() == 21 && xproc[2]->id() == 21);
}


void MatchboxAmplitudeggttbarg::prepareAmplitudes(Ptr<MatchboxMEBase>::tcptr me) {

  if ( !calculateTreeAmplitudes() ) {
    MatchboxAmplitude::prepareAmplitudes(me);
    return;
  }
  amplitudeScale(sqrt(lastSHat()));
  momentum(0,amplitudeMomentum(0));
  momentum(1,amplitudeMomentum(1));
  momentum(2,amplitudeMomentum(2));
  momentum(3,amplitudeMomentum(3));
  momentum(4,amplitudeMomentum(4));

  MatchboxAmplitude::prepareAmplitudes(me);

}

Complex MatchboxAmplitudeggttbarg::evaluate(size_t a, const vector<int>& hel, Complex& largeN) {

  // set up momenta to pass into Madgraph amplitude
  vector<double*> p;
  p.push_back(mom0);  p.push_back(mom1);  p.push_back(mom2);  p.push_back(mom3); p.push_back(mom4);
  for (int ip=0; ip<5; ++ip){
    p[ip][0] = abs(momentum(ip).e())<1.e-13 ? 0.0:double(momentum(ip).e()*amplitudeScale()/GeV);
    p[ip][1] = abs(momentum(ip).x())<1.e-13 ? 0.0:double(momentum(ip).x()*amplitudeScale()/GeV);
    p[ip][2] = abs(momentum(ip).y())<1.e-13 ? 0.0:double(momentum(ip).y()*amplitudeScale()/GeV);
    p[ip][3] = abs(momentum(ip).z())<1.e-13 ? 0.0:double(momentum(ip).z()*amplitudeScale()/GeV);
  } 

  // calculate amplitudes
  vector<complex<double> > amplitudes;
  MG_gg2ttxg process;
  double factor=lastSHat()/GeV2;
  process.initProc(MGParams_);
  process.setMomenta(p);
  process.sigmaKin(amplitudes, hel);

  for (int iamp=0; iamp<int(amplitudes.size()); ++iamp)
    amplitudes[iamp]*=sqrt(factor);
  complex<double> matrixElement;
  complex<double> i = complex<double>(0,1);
  
  // calculate colour flows
  if (a<5)          matrixElement = 0.;
  else if (a==5)    
    matrixElement = (-amplitudes[0]  -  amplitudes[5]  + amplitudes[14] - 
  		      amplitudes[17] +  amplitudes[15] + 
		   i*(amplitudes[2]  +  amplitudes[4]));
  else if (a==6)
    matrixElement = (-amplitudes[3]  +  amplitudes[11] - amplitudes[14] - 
		      amplitudes[16] -  amplitudes[15] + 
		  i*(-amplitudes[4]  +  amplitudes[10])); 
  else if (a==7)
    matrixElement = ( amplitudes[0]  -   amplitudes[11] - amplitudes[12] + 
  		      amplitudes[17] +   amplitudes[16] + 
		  i*(-amplitudes[2]  +   amplitudes[9]));
  else if (a==8)
    matrixElement = (-amplitudes[6]  +   amplitudes[11] - amplitudes[14] - 
		      amplitudes[16] -   amplitudes[15] 
	         + i*(amplitudes[7]  -   amplitudes[9]));
  else if (a==9)
    matrixElement = (amplitudes[0]  -   amplitudes[11] - amplitudes[13] + 
		     amplitudes[17] +   amplitudes[16] 
		+ i*(amplitudes[1]  -   amplitudes[10]));
  else if (a==10)
    matrixElement = (-amplitudes[0]  -   amplitudes[8] + amplitudes[14] - 
		      amplitudes[17] +   amplitudes[15] 
		+ i*(-amplitudes[1]  -   amplitudes[7]));
  else assert(false);
  
  largeN = matrixElement;
  return matrixElement;
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MatchboxAmplitudeggttbarg::persistentOutput(PersistentOStream &) const {}

void MatchboxAmplitudeggttbarg::persistentInput(PersistentIStream &, int) {}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MatchboxAmplitudeggttbarg,MatchboxAmplitude>
  describeHerwigMatchboxAmplitudeggttbarg("Herwig::MatchboxAmplitudeggttbarg", "HwMatchboxBuiltin.so");

void MatchboxAmplitudeggttbarg::Init() {

  static ClassDocumentation<MatchboxAmplitudeggttbarg> documentation
    ("MatchboxAmplitudeggttbarg");

}

#line 1 "./HelAmps_sm.cc"
//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include "HelAmps_sm.h"
#include <iostream>
namespace MG5_sm 
{


double Sgn(double a, double b)
{
  return (b < 0)? - abs(a):abs(a); 
}


void ixxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  if (nsf==1 && p[0]>=0.)       i2xxxx(p, fmass, nhel, nsf, fi);
  else if (nsf==-1 && p[0]>=0.) i2xxxx(p, fmass, nhel, nsf, fi);
  else if (nsf==1 && p[0]<0.){
    p[0]*=-1; p[1]*=-1; p[2]*=-1; p[3]*=-1;
    o2xxxx(p, fmass, nhel, -1.*nsf, fi);
    p[0]*=-1; p[1]*=-1; p[2]*=-1; p[3]*=-1;
  }
  else if (nsf==-1 && p[0]<0.){
    p[0]*=-1; p[1]*=-1; p[2]*=-1; p[3]*=-1;
    o2xxxx(p, fmass, nhel, -1.*nsf, fi);
    p[0]*=-1; p[1]*=-1; p[2]*=-1; p[3]*=-1;
  }
  return;
}

  void oxxxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  if (nsf==1 && p[0]>=0.)       o2xxxx(p, fmass, nhel, nsf, fo);
  else if (nsf==-1 && p[0]>=0.) o2xxxx(p, fmass, nhel, nsf, fo);
  else if (nsf==1 && p[0]<0.){
    p[0]*=-1.; p[1]*=-1.; p[2]*=-1.; p[3]*=-1.;
    i2xxxx(p, fmass, nhel, -1.*nsf, fo);
    p[0]*=-1.; p[1]*=-1.; p[2]*=-1.; p[3]*=-1.;
  }
  else if (nsf==-1 && p[0]<0.){
    p[0]*=-1.; p[1]*=-1.; p[2]*=-1.; p[3]*=-1.;
    i2xxxx(p, fmass, nhel, -1.*nsf, fo);
    p[0]*=-1.; p[1]*=-1.; p[2]*=-1.; p[3]*=-1.;
  }
  return;
}

void sxxxxx(double p[4], int nss, complex<double> sc[3])
{
  s2xxxx(p, nss, sc);
  return; 
}

void vxxxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  v2xxxx(p, vmass, nhel, nsv, vc);
  return;
}

  void i2xxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fi[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomega[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int ip, im, nh; 
  fi[0] = complex<double> (-p[0] * nsf, -p[3] * nsf); 
  fi[1] = complex<double> (-p[1] * nsf, -p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.0)
  {
    pp = min(p[0], pow((pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2)), 0.5)); 
    if (pp == 0.0)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      fi[2] = ip * sqm[ip]; 
      fi[3] = im * nsf * sqm[ip]; 
      fi[4] = ip * nsf * sqm[im]; 
      fi[5] = im * sqm[im]; 
    }
    else
    {
      sf[0] = (1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = (1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomega[0] = sf[0] * omega[ip]; 
      sfomega[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.0); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0); 
      if (pp3 == 0.0)
      {
        chi[1] = complex<double> (-nh, 0); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fi[2] = sfomega[0] * chi[im]; 
      fi[3] = sfomega[0] * chi[ip]; 
      fi[4] = sfomega[1] * chi[im]; 
      fi[5] = sfomega[1] * chi[ip]; 
    }
  }
  else
  {
    if (p[1] == 0.0 and p[2] == 0.0 and p[3] < 0.0)
    {
      sqp0p3 = 0.0; 
    }
    else
    {     
      sqp0p3 = pow(max(p[0] + p[3], 0.0), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.0); 
    if (sqp0p3 == 0.0)
    {
      chi[1] = complex<double> (-nhel * pow(2.0 * p[0], 0.5), 0.0); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], p[2])/sqp0p3; 
    }
    if (nh == 1)
    {
      fi[2] = complex<double> (0.0, 0.0); 
      fi[3] = complex<double> (0.0, 0.0); 
      fi[4] = chi[0]; 
      fi[5] = chi[1]; 
    }
    else
    {
      fi[2] = chi[1]; 
      fi[3] = chi[0]; 
      fi[4] = complex<double> (0.0, 0.0); 
      fi[5] = complex<double> (0.0, 0.0); 
    }
  }
  return; 
}

void o2xxxx(double p[4], double fmass, int nhel, int nsf, complex<double> fo[6])
{
  complex<double> chi[2]; 
  double sf[2], sfomeg[2], omega[2], pp, pp3, sqp0p3, sqm[2]; 
  int nh, ip, im; 
  fo[0] = complex<double> (p[0] * nsf, p[3] * nsf); 
  fo[1] = complex<double> (p[1] * nsf, p[2] * nsf); 
  nh = nhel * nsf; 
  if (fmass != 0.000)
  {
    pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
    if (pp == 0.000)
    {
      sqm[0] = pow(abs(fmass), 0.5); 
      sqm[1] = Sgn(sqm[0], fmass); 
      ip = -((1 - nh)/2) * nhel; 
      im = (1 + nh)/2 * nhel; 
      fo[2] = im * sqm[im]; 
      fo[3] = ip * nsf * sqm[im]; 
      fo[4] = im * nsf * sqm[ - ip]; 
      fo[5] = ip * sqm[ - ip]; 
    }
    else
    {
      pp = min(p[0], pow(pow(p[1], 2) + pow(p[2], 2) + pow(p[3], 2), 0.5)); 
      sf[0] = double(1 + nsf + (1 - nsf) * nh) * 0.5; 
      sf[1] = double(1 + nsf - (1 - nsf) * nh) * 0.5; 
      omega[0] = pow(p[0] + pp, 0.5); 
      omega[1] = fmass/omega[0]; 
      ip = (1 + nh)/2; 
      im = (1 - nh)/2; 
      sfomeg[0] = sf[0] * omega[ip]; 
      sfomeg[1] = sf[1] * omega[im]; 
      pp3 = max(pp + p[3], 0.00); 
      chi[0] = complex<double> (pow(pp3 * 0.5/pp, 0.5), 0.00); 
      if (pp3 == 0.00)
      {
        chi[1] = complex<double> (-nh, 0.00); 
      }
      else
      {
        chi[1] = complex<double> (nh * p[1], -p[2])/pow(2.0 * pp * pp3, 0.5); 
      }
      fo[2] = sfomeg[1] * chi[im]; 
      fo[3] = sfomeg[1] * chi[ip]; 
      fo[4] = sfomeg[0] * chi[im]; 
      fo[5] = sfomeg[0] * chi[ip]; 
    }
  }
  else
  {
    if((p[1] == 0.00) and (p[2] == 0.00) and (p[3] < 0.00))
    {
      sqp0p3 = 0.00; 
    }
    else
    {
      sqp0p3 = pow(max(p[0] + p[3], 0.00), 0.5) * nsf; 
    }
    chi[0] = complex<double> (sqp0p3, 0.00); 
    if(sqp0p3 == 0.000)
    {
      chi[1] = complex<double> (-nhel, 0.00) * pow(2.0 * p[0], 0.5); 
    }
    else
    {
      chi[1] = complex<double> (nh * p[1], -p[2])/sqp0p3; 
    }
    if(nh == 1)
    {
      fo[2] = chi[0]; 
      fo[3] = chi[1]; 
      fo[4] = complex<double> (0.00, 0.00); 
      fo[5] = complex<double> (0.00, 0.00); 
    }
    else
    {
      fo[2] = complex<double> (0.00, 0.00); 
      fo[3] = complex<double> (0.00, 0.00); 
      fo[4] = chi[1]; 
      fo[5] = chi[0]; 
    }
  }
  return; 
}

void s2xxxx(double p[4], int nss, complex<double> sc[3])
{
  sc[2] = complex<double> (1.00, 0.00); 
  sc[0] = complex<double> (p[0] * nss, p[3] * nss); 
  sc[1] = complex<double> (p[1] * nss, p[2] * nss); 
  return; 
}

void v2xxxx(double p[4], double vmass, int nhel, int nsv, complex<double> vc[6])
{
  double hel, hel0, pt, pt2, pp, pzpt, emp, sqh; 
  int nsvahl; 
  sqh = pow(0.5, 0.5); 
  hel = double(nhel); 
  nsvahl = nsv * abs(hel); 
  pt2 = pow(p[1], 2) + pow(p[2], 2); 
  pp = min(p[0], pow(pt2 + pow(p[3], 2), 0.5)); 
  pt = min(pp, pow(pt2, 0.5)); 
  vc[0] = complex<double> (p[0] * nsv, p[3] * nsv); 
  vc[1] = complex<double> (p[1] * nsv, p[2] * nsv); 
  if (vmass != 0.0)
  {
    hel0 = 1.0 - abs(hel); 
    if(pp == 0.0)
    {
      vc[2] = complex<double> (0.0, 0.0); 
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsvahl * sqh); 
      vc[5] = complex<double> (hel0, 0.0); 
    }
    else
    {
      emp = p[0]/(vmass * pp); 
      vc[2] = complex<double> (hel0 * pp/vmass, 0.0); 
      vc[5] = complex<double> (hel0 * p[3] * emp + hel * pt/pp * sqh, 0.0); 
      if (pt != 0.0)
      {
        pzpt = p[3]/(pp * pt) * sqh * hel; 
        vc[3] = complex<double> (hel0 * p[1] * emp - p[1] * pzpt, -nsvahl *
            p[2]/pt * sqh);
        vc[4] = complex<double> (hel0 * p[2] * emp - p[2] * pzpt, nsvahl *
            p[1]/pt * sqh);
      }
      else
      {
        vc[3] = complex<double> (-hel * sqh, 0.0); 
        vc[4] = complex<double> (0.0, nsvahl * Sgn(sqh, p[3])); 
      }
    }
  }
  else
  {
    pp = p[0]; 
    pt = pow(pow(p[1], 2) + pow(p[2], 2), 0.5); 
    vc[2] = complex<double> (0.0, 0.0); 
    vc[5] = complex<double> (hel * pt/pp * sqh, 0.0); 
    if (pt != 0.0)
    {
      pzpt = p[3]/(pp * pt) * sqh * hel; 
      vc[3] = complex<double> (-p[1] * pzpt, -nsv * p[2]/pt * sqh); 
      vc[4] = complex<double> (-p[2] * pzpt, nsv * p[1]/pt * sqh); 
    }
    else
    {
      vc[3] = complex<double> (-hel * sqh, 0.0); 
      vc[4] = complex<double> (0.0, nsv * Sgn(sqh, p[3])); 
    }
  }
  return; 
}

void FFV1_0(complex<double> F1[], complex<double> F2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP3; 
  TMP3 = (F1[2] * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI * (V3[4]))) +
      (F1[3] * (F2[4] * (V3[3] - cI * (V3[4])) + F2[5] * (V3[2] - V3[5])) +
      (F1[4] * (F2[2] * (V3[2] - V3[5]) - F2[3] * (V3[3] + cI * (V3[4]))) +
      F1[5] * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] + V3[5])))));
  vertex = COUP * - cI * TMP3; 
}

void FFV1P0_3(complex<double> F1[], complex<double> F2[], complex<double> COUP,
    double M3, double W3, complex<double> V3[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P3[4]; 
  complex<double> denom; 
  V3[0] = +F1[0] + F2[0]; 
  V3[1] = +F1[1] + F2[1]; 
  P3[0] = -V3[0].real(); 
  P3[1] = -V3[1].real(); 
  P3[2] = -V3[1].imag(); 
  P3[3] = -V3[0].imag(); 
  denom = COUP/(pow(P3[0], 2) - pow(P3[1], 2) - pow(P3[2], 2) - pow(P3[3], 2) -
      M3 * (M3 - cI * W3));
  V3[2] = denom * - cI * (F1[2] * F2[4] + F1[3] * F2[5] + F1[4] * F2[2] + F1[5]
      * F2[3]);
  V3[3] = denom * - cI * (F1[4] * F2[3] + F1[5] * F2[2] - F1[2] * F2[5] - F1[3]
      * F2[4]);
  V3[4] = denom * - cI * (-cI * (F1[2] * F2[5] + F1[5] * F2[2]) + cI * (F1[3] *
      F2[4] + F1[4] * F2[3]));
  V3[5] = denom * - cI * (F1[3] * F2[5] + F1[4] * F2[2] - F1[2] * F2[4] - F1[5]
      * F2[3]);
}
void VVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  double P3[4]; 
  complex<double> TMP5; 
  complex<double> TMP4; 
  complex<double> denom; 
  complex<double> TMP3; 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  V1[0] = +V2[0] + V3[0]; 
  V1[1] = +V2[1] + V3[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP5 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP4 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP1 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  TMP3 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP2 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (TMP5 * (-cI * (P2[0]) + cI * (P3[0])) + (V2[2] * (-cI *
      (TMP1) + cI * (TMP2)) + V3[2] * (-cI * (TMP4) + cI * (TMP3))));
  V1[3] = denom * (TMP5 * (-cI * (P2[1]) + cI * (P3[1])) + (V2[3] * (-cI *
      (TMP1) + cI * (TMP2)) + V3[3] * (-cI * (TMP4) + cI * (TMP3))));
  V1[4] = denom * (TMP5 * (-cI * (P2[2]) + cI * (P3[2])) + (V2[4] * (-cI *
      (TMP1) + cI * (TMP2)) + V3[4] * (-cI * (TMP4) + cI * (TMP3))));
  V1[5] = denom * (TMP5 * (-cI * (P2[3]) + cI * (P3[3])) + (V2[5] * (-cI *
      (TMP1) + cI * (TMP2)) + V3[5] * (-cI * (TMP4) + cI * (TMP3))));
}

void FFV1_2(complex<double> F1[], complex<double> V3[], complex<double> COUP,
    double M2, double W2, complex<double> F2[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P2[4]; 
  complex<double> denom; 
  F2[0] = +F1[0] + V3[0]; 
  F2[1] = +F1[1] + V3[1]; 
  P2[0] = -F2[0].real(); 
  P2[1] = -F2[1].real(); 
  P2[2] = -F2[1].imag(); 
  P2[3] = -F2[0].imag(); 
  denom = COUP/(pow(P2[0], 2) - pow(P2[1], 2) - pow(P2[2], 2) - pow(P2[3], 2) -
      M2 * (M2 - cI * W2));
  F2[2] = denom * cI * (F1[2] * (P2[0] * (V3[2] + V3[5]) + (P2[1] * - 1. *
      (V3[3] + cI * (V3[4])) + (P2[2] * (+cI * (V3[3]) - V3[4]) - P2[3] *
      (V3[2] + V3[5])))) + (F1[3] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] *
      (V3[5] - V3[2]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] * (+cI
      * (V3[4]) - V3[3])))) + M2 * (F1[4] * (V3[2] - V3[5]) + F1[5] * (+cI *
      (V3[4]) - V3[3]))));
  F2[3] = denom * - cI * (F1[2] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) +
      (P2[1] * (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) - P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[3] * (P2[0] * (V3[5] - V3[2]) + (P2[1] *
      (V3[3] - cI * (V3[4])) + (P2[2] * (V3[4] + cI * (V3[3])) + P2[3] * (V3[5]
      - V3[2])))) + M2 * (F1[4] * (V3[3] + cI * (V3[4])) - F1[5] * (V3[2] +
      V3[5]))));
  F2[4] = denom * - cI * (F1[4] * (P2[0] * (V3[5] - V3[2]) + (P2[1] * (V3[3] +
      cI * (V3[4])) + (P2[2] * (V3[4] - cI * (V3[3])) + P2[3] * (V3[5] -
      V3[2])))) + (F1[5] * (P2[0] * (V3[3] - cI * (V3[4])) + (P2[1] * - 1. *
      (V3[2] + V3[5]) + (P2[2] * (+cI * (V3[2] + V3[5])) + P2[3] * (V3[3] - cI
      * (V3[4]))))) + M2 * (F1[2] * - 1. * (V3[2] + V3[5]) + F1[3] * (+cI *
      (V3[4]) - V3[3]))));
  F2[5] = denom * cI * (F1[4] * (P2[0] * - 1. * (V3[3] + cI * (V3[4])) + (P2[1]
      * (V3[2] - V3[5]) + (P2[2] * (-cI * (V3[5]) + cI * (V3[2])) + P2[3] *
      (V3[3] + cI * (V3[4]))))) + (F1[5] * (P2[0] * (V3[2] + V3[5]) + (P2[1] *
      (+cI * (V3[4]) - V3[3]) + (P2[2] * - 1. * (V3[4] + cI * (V3[3])) - P2[3]
      * (V3[2] + V3[5])))) + M2 * (F1[2] * (V3[3] + cI * (V3[4])) + F1[3] *
      (V3[2] - V3[5]))));
}

void FFV1_1(complex<double> F2[], complex<double> V3[], complex<double> COUP,
    double M1, double W1, complex<double> F1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> denom; 
  F1[0] = +F2[0] + V3[0]; 
  F1[1] = +F2[1] + V3[1]; 
  P1[0] = -F1[0].real(); 
  P1[1] = -F1[1].real(); 
  P1[2] = -F1[1].imag(); 
  P1[3] = -F1[0].imag(); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  F1[2] = denom * cI * (F2[2] * (P1[0] * (V3[5] - V3[2]) + (P1[1] * (V3[3] - cI
      * (V3[4])) + (P1[2] * (V3[4] + cI * (V3[3])) + P1[3] * (V3[5] - V3[2]))))
      + (F2[3] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] * - 1. * (V3[2] +
      V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (V3[3] + cI *
      (V3[4]))))) + M1 * (F2[4] * (V3[2] + V3[5]) + F2[5] * (V3[3] + cI *
      (V3[4])))));
  F1[3] = denom * - cI * (F2[2] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] - V3[5]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) + P1[3] *
      (V3[3] - cI * (V3[4]))))) + (F2[3] * (P1[0] * (V3[2] + V3[5]) + (P1[1] *
      - 1. * (V3[3] + cI * (V3[4])) + (P1[2] * (+cI * (V3[3]) - V3[4]) - P1[3]
      * (V3[2] + V3[5])))) + M1 * (F2[4] * (+cI * (V3[4]) - V3[3]) + F2[5] *
      (V3[5] - V3[2]))));
  F1[4] = denom * - cI * (F2[4] * (P1[0] * (V3[2] + V3[5]) + (P1[1] * (+cI *
      (V3[4]) - V3[3]) + (P1[2] * - 1. * (V3[4] + cI * (V3[3])) - P1[3] *
      (V3[2] + V3[5])))) + (F2[5] * (P1[0] * (V3[3] + cI * (V3[4])) + (P1[1] *
      (V3[5] - V3[2]) + (P1[2] * (-cI * (V3[2]) + cI * (V3[5])) - P1[3] *
      (V3[3] + cI * (V3[4]))))) + M1 * (F2[2] * (V3[5] - V3[2]) + F2[3] *
      (V3[3] + cI * (V3[4])))));
  F1[5] = denom * cI * (F2[4] * (P1[0] * (+cI * (V3[4]) - V3[3]) + (P1[1] *
      (V3[2] + V3[5]) + (P1[2] * - 1. * (+cI * (V3[2] + V3[5])) + P1[3] * (+cI
      * (V3[4]) - V3[3])))) + (F2[5] * (P1[0] * (V3[5] - V3[2]) + (P1[1] *
      (V3[3] + cI * (V3[4])) + (P1[2] * (V3[4] - cI * (V3[3])) + P1[3] * (V3[5]
      - V3[2])))) + M1 * (F2[2] * (+cI * (V3[4]) - V3[3]) + F2[3] * (V3[2] +
      V3[5]))));
}

void VVV1_0(complex<double> V1[], complex<double> V2[], complex<double> V3[],
    complex<double> COUP, complex<double> & vertex)
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP2; 
  complex<double> TMP1; 
  double P1[4]; 
  double P2[4]; 
  complex<double> TMP7; 
  double P3[4]; 
  complex<double> TMP6; 
  complex<double> TMP5; 
  complex<double> TMP4; 
  complex<double> TMP9; 
  complex<double> TMP3; 
  complex<double> TMP8; 
  P1[0] = V1[0].real(); 
  P1[1] = V1[1].real(); 
  P1[2] = V1[1].imag(); 
  P1[3] = V1[0].imag(); 
  P2[0] = V2[0].real(); 
  P2[1] = V2[1].real(); 
  P2[2] = V2[1].imag(); 
  P2[3] = V2[0].imag(); 
  P3[0] = V3[0].real(); 
  P3[1] = V3[1].real(); 
  P3[2] = V3[1].imag(); 
  P3[3] = V3[0].imag(); 
  TMP9 = (V1[2] * P3[0] - V1[3] * P3[1] - V1[4] * P3[2] - V1[5] * P3[3]); 
  TMP8 = (P2[0] * V1[2] - P2[1] * V1[3] - P2[2] * V1[4] - P2[3] * V1[5]); 
  TMP5 = (V3[2] * V1[2] - V3[3] * V1[3] - V3[4] * V1[4] - V3[5] * V1[5]); 
  TMP4 = (P1[0] * V2[2] - P1[1] * V2[3] - P1[2] * V2[4] - P1[3] * V2[5]); 
  TMP7 = (V3[2] * V2[2] - V3[3] * V2[3] - V3[4] * V2[4] - V3[5] * V2[5]); 
  TMP6 = (V2[2] * P3[0] - V2[3] * P3[1] - V2[4] * P3[2] - V2[5] * P3[3]); 
  TMP1 = (V2[2] * V1[2] - V2[3] * V1[3] - V2[4] * V1[4] - V2[5] * V1[5]); 
  TMP3 = (V3[2] * P2[0] - V3[3] * P2[1] - V3[4] * P2[2] - V3[5] * P2[3]); 
  TMP2 = (V3[2] * P1[0] - V3[3] * P1[1] - V3[4] * P1[2] - V3[5] * P1[3]); 
  vertex = COUP * (TMP1 * (-cI * (TMP2) + cI * (TMP3)) + (TMP5 * (-cI * (TMP6)
      + cI * (TMP4)) + TMP7 * (-cI * (TMP8) + cI * (TMP9))));
}

void VVVV1P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  double P1[4]; 
  complex<double> TMP0; 
  complex<double> denom; 
  complex<double> TMP3; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP0 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  TMP3 = (V2[2] * V3[2] - V2[3] * V3[3] - V2[4] * V3[4] - V2[5] * V3[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V4[2] * TMP3) + cI * (V3[2] * TMP0)); 
  V1[3] = denom * (-cI * (V4[3] * TMP3) + cI * (V3[3] * TMP0)); 
  V1[4] = denom * (-cI * (V4[4] * TMP3) + cI * (V3[4] * TMP0)); 
  V1[5] = denom * (-cI * (V4[5] * TMP3) + cI * (V3[5] * TMP0)); 
}


void VVVV4P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> TMP0; 
  complex<double> denom; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP1 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP0 = (V4[2] * V2[2] - V4[3] * V2[3] - V4[4] * V2[4] - V4[5] * V2[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V3[2] * TMP0) + cI * (V2[2] * TMP1)); 
  V1[3] = denom * (-cI * (V3[3] * TMP0) + cI * (V2[3] * TMP1)); 
  V1[4] = denom * (-cI * (V3[4] * TMP0) + cI * (V2[4] * TMP1)); 
  V1[5] = denom * (-cI * (V3[5] * TMP0) + cI * (V2[5] * TMP1)); 
}

void VVVV3P0_1(complex<double> V2[], complex<double> V3[], complex<double>
    V4[], complex<double> COUP, double M1, double W1, complex<double> V1[])
{
  complex<double> cI = complex<double> (0., 1.); 
  complex<double> TMP1; 
  double P1[4]; 
  complex<double> denom; 
  complex<double> TMP3; 
  V1[0] = +V2[0] + V3[0] + V4[0]; 
  V1[1] = +V2[1] + V3[1] + V4[1]; 
  P1[0] = -V1[0].real(); 
  P1[1] = -V1[1].real(); 
  P1[2] = -V1[1].imag(); 
  P1[3] = -V1[0].imag(); 
  TMP1 = (V4[2] * V3[2] - V4[3] * V3[3] - V4[4] * V3[4] - V4[5] * V3[5]); 
  TMP3 = (V2[2] * V3[2] - V2[3] * V3[3] - V2[4] * V3[4] - V2[5] * V3[5]); 
  denom = COUP/(pow(P1[0], 2) - pow(P1[1], 2) - pow(P1[2], 2) - pow(P1[3], 2) -
      M1 * (M1 - cI * W1));
  V1[2] = denom * (-cI * (V4[2] * TMP3) + cI * (V2[2] * TMP1)); 
  V1[3] = denom * (-cI * (V4[3] * TMP3) + cI * (V2[3] * TMP1)); 
  V1[4] = denom * (-cI * (V4[4] * TMP3) + cI * (V2[4] * TMP1)); 
  V1[5] = denom * (-cI * (V4[5] * TMP3) + cI * (V2[5] * TMP1)); 
}

 // end namespace $(namespace)s_s
}
#line 1 "./Parameters_sm.cc"
//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph 5 v. 1.5.7, 2013-01-15
// By the MadGraph Development Team
// Please visit us at https://launchpad.net/madgraph5
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_sm.h"

// Initialize static instance
Parameters_sm * Parameters_sm::instance = 0; 

// Function to get static instance - only one instance per program
Parameters_sm * Parameters_sm::getInstance()
{
  if (instance == 0)
    instance = new Parameters_sm(); 

  return instance; 
}

void Parameters_sm::setIndependentCouplings()
{
  GC_1 =  -(ee * complexi)/3.; 
  GC_2 =   (2. * ee * complexi)/3.; 
  GC_3 =  -(ee * complexi); 
  GC_50 = -(cw * ee * complexi)/(2. * sw); 
  GC_51 =  (cw * ee * complexi)/(2. * sw); 
  GC_58 = -(ee * complexi * sw)/(6. * cw); 
  GC_59 =  (ee * complexi * sw)/(2. * cw); 
}

void Parameters_sm::setDependentParameters()
{
  sqrt__aS  = sqrt(aS); 
  G         = 2. * sqrt__aS * sqrt(M_PI); 
  G__exp__2 = pow(G, 2.); 
}

void Parameters_sm::setDependentCouplings()
{
  GC_10 = -G; 
  GC_11 = complexi * G;
  GC_12 = complexi * G__exp__2; 
}

void Parameters_sm::setIndependentParameters(map<string, double> &MGParams)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  WH = MGParams.find("WH") ->second;
  WW = MGParams.find("WW") ->second;
  WZ = MGParams.find("WZ") ->second;
  WT = MGParams.find("WT") ->second;
  Gf = MGParams.find("GF") ->second;
  MH = MGParams.find("MH") ->second;
  MZ = MGParams.find("MZ") ->second;
  MW = MGParams.find("MW") ->second;
  MTA =MGParams.find("MTA") ->second;
  MT = MGParams.find("MT") ->second;
  MB = MGParams.find("MB") ->second;
  aS = MGParams.find("aS") ->second;
  aEWM1 = MGParams.find("aEWM1") ->second;
  ymtau = MTA;
  ymt   = MT;
  ymb   = MT;
  cw    = MW/MZ;
  cw__exp__2 = pow(cw, 2.); 
  sw = sqrt(1. - cw__exp__2);
  sw__exp_2 = pow(sw, 2.);
  conjg__CKM1x1 = 1.; 
  conjg__CKM3x3 = 1.; 
  CKM3x3 = 1.; 
  complexi = std::complex<double> (0., 1.); 
  MZ__exp__2 = pow(MZ, 2.); 
  MZ__exp__4 = pow(MZ, 4.); 
  sqrt__2 = sqrt(2.); 
  MH__exp__2 = pow(MH, 2.); 
  aEW = 1./aEWM1; 
  sqrt__aEW = sqrt(aEW); 
  ee = 2. * sqrt__aEW * sqrt(M_PI); 
  MW__exp__2 = pow(MW, 2.);
  g1 = ee/cw; 
  gw = ee/sw; 
  vev = (2. * MW * sw)/ee; 
  vev__exp__2 = pow(vev, 2.); 
  lam = MH__exp__2/(2. * vev__exp__2); 
  yb = (ymb * sqrt__2)/vev; 
  yt = (ymt * sqrt__2)/vev; 
  ytau = (ymtau * sqrt__2)/vev; 
  muH = sqrt(lam * vev__exp__2); 
  I1x33 = yb * conjg__CKM3x3; 
  I2x33 = yt * conjg__CKM3x3; 
  I3x33 = CKM3x3 * yt; 
  I4x33 = CKM3x3 * yb; 
  ee__exp__2 = pow(ee, 2.); 
  sw__exp__2 = pow(sw, 2.); 
  cw__exp__2 = pow(cw, 2.); 
}
#line 1 "./MG_qqx2ttx.cc"
// -*- C++ -*-
//
// MG_qqx2ttx.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_qqx2ttx class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15
 

#include "MG_qqx2ttx.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 

void MG_qqx2ttx::initProc(map<string, double> & MGParams) 
{
  // Instantiate the model class and set parameters using Herwig values
  pars = Parameters_sm::getInstance(); 
  pars->setIndependentParameters(MGParams); 
  pars->setIndependentCouplings(); 

  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->MT); 
}

void MG_qqx2ttx::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel) 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 

  // Reset the matrix elements
  for(int i = 0; i < namplitudes; i++ )
    amp[i] = 0.; 
 
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
    perm[i] = i; 

  // Set vector of helicities
  int helicities[1][nexternal] = {{hel[0], hel[1], hel[2], hel[3]}};

  // Calculate amplitudes
  calculate_wavefunctions(perm, helicities[0]);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;
}

void MG_qqx2ttx::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate all wavefunctions
  ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  swap(w[0],w[1]);
  oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 

  FFV1P0_3(w[0], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[4]); 
  // Calculate all amplitudes
  FFV1_0(w[3], w[2], w[4], pars->GC_11, amp[0]); 
}
#line 1 "./MG_qqx2ttxg.cc"
// -*- C++ -*-
//
// MG_qqx2ttxg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_qqx2ttxg class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15

#include "MG_qqx2ttxg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 

void MG_qqx2ttxg::initProc(map<string, double> & MGParams) 
{
  // Instantiate the model class and set parameters using Herwig values
  pars = Parameters_sm::getInstance(); 
  pars->setIndependentParameters(MGParams); 
  pars->setIndependentCouplings(); 

  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->ZERO); 
}

void MG_qqx2ttxg::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel, int crossed) 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 

  // Reset the matrix elements
  for(int i = 0; i < namplitudes; i++ )
    amp[i] = 0.; 
 
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
    perm[i] = i;

  // Set vector of helicities
  int helicities[1][nexternal] = {{hel[0], hel[1], hel[2], hel[3], hel[4]}};

  // Calculate amplitudes
  calculate_wavefunctions(perm, helicities[0], crossed);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;
}

void MG_qqx2ttxg::calculate_wavefunctions(const int perm[], const int hel[], int crossed)
{
  // Calculate all wavefunctions 
  // q qbar initiated
  if (crossed==1){   
    ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
    oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
    swap(w[0],w[1]);
    oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
    ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
    vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  }
  // q g initiated
  else if (crossed==2){
    oxxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]); 
    oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
    swap(w[0],w[1]);
    oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
    ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
    vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 
  }
  // qbar g initiated
  else if (crossed==3){
    ixxxxx(p[perm[0]], mME[0], hel[0], +1, w[0]);
    ixxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
    swap(w[0],w[1]);
    oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
    ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]); 
    vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  }

  FFV1P0_3(w[0], w[1], pars->GC_11, pars->ZERO, pars->ZERO, w[5]); 
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[6]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->MT, pars->ZERO, w[7]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->MT, pars->ZERO, w[8]); 
  FFV1_2(w[0], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[9]); 
  FFV1_1(w[1], w[4], pars->GC_11, pars->ZERO, pars->ZERO, w[10]); 
  // Calculate all amplitudes
  VVV1_0(w[5], w[6], w[4], pars->GC_10,  amp[0]); 
  FFV1_0(w[3], w[7], w[5], pars->GC_11,  amp[1]); 
  FFV1_0(w[8], w[2], w[5], pars->GC_11,  amp[2]); 
  FFV1_0(w[9], w[1], w[6], pars->GC_11,  amp[3]); 
  FFV1_0(w[0], w[10], w[6], pars->GC_11, amp[4]); 
}
#line 1 "./MG_gg2ttx.cc"
// -*- C++ -*-
//
// MG_gg2ttx.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_gg2ttx class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15


#include "MG_gg2ttx.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 


void MG_gg2ttx::initProc(map<string, double> & MGParams) 
{
  // Instantiate the model class and set parameters using Herwig values
  pars = Parameters_sm::getInstance(); 
  pars->setIndependentParameters(MGParams); 
  pars->setIndependentCouplings(); 

  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->MT); 
}


void MG_gg2ttx::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel) 
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 

  // Reset the matrix elements
  for(int i = 0; i < namplitudes; i++ )
    amp[i] = 0.; 
 
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
    perm[i] = i; 

  // Set vector of helicities
  int helicities[1][nexternal] = {{hel[0], hel[1], hel[2], hel[3]}};

  // Calculate amplitudes
  calculate_wavefunctions(perm, helicities[0]);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;
}


void MG_gg2ttx::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate all wavefunctions
  vxxxxx(p[perm[2]], mME[0], hel[2], +1, w[0]); 
  vxxxxx(p[perm[3]], mME[1], hel[3], +1, w[1]); 
  swap(w[0],w[1]);
  oxxxxx(p[perm[0]], mME[2], hel[0], +1, w[2]); 
  ixxxxx(p[perm[1]], mME[3], hel[1], -1, w[3]); 

  VVV1P0_1(w[0], w[1], pars->GC_10, pars->ZERO, pars->ZERO, w[4]); 
  FFV1_1(w[2], w[0], pars->GC_11, pars->MT, pars->ZERO, w[5]); 
  FFV1_2(w[3], w[0], pars->GC_11, pars->MT, pars->ZERO, w[6]); 

  // Calculate all amplitudes
  FFV1_0(w[3], w[2], w[4], pars->GC_11, amp[0]); 
  FFV1_0(w[3], w[5], w[1], pars->GC_11, amp[1]); 
  FFV1_0(w[6], w[2], w[1], pars->GC_11, amp[2]); 

}
#line 1 "./MG_gg2ttxg.cc"
// -*- C++ -*-
//
// MG_gg2ttxg.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MG_gg2ttxg class.
//
//

// Adapted from file automatically generated for C++ Standalone by
// MadGraph 5 v. 1.5.15


#include "MG_gg2ttxg.h"
#include "HelAmps_sm.h"
#include <iostream>

using namespace MG5_sm; 


void MG_gg2ttxg::initProc(map<string, double> & MGParams) 
{
  // Instantiate the model class and set parameters using Herwig values
  pars = Parameters_sm::getInstance(); 
  pars->setIndependentParameters(MGParams); 
  pars->setIndependentCouplings(); 

  // Set external particle masses for this matrix element
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->ZERO); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->MT); 
  mME.push_back(pars->ZERO); 
  
}


void MG_gg2ttxg::sigmaKin(vector<complex<double> >& amps, const vector<int>& hel)
{
  // Set the parameters which change event by event
  pars->setDependentParameters(); 
  pars->setDependentCouplings(); 

  // Reset the matrix elements
  for(int i = 0; i < namplitudes; i++ )
    amp[i] = 0.; 
 
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
    perm[i] = i;

  // Set vector of helicities
  int helicities[1][nexternal] = {{hel[0], hel[1], hel[2], hel[3], hel[4]}};

  // Calculate amplitudes
  calculate_wavefunctions(perm, helicities[0]);   
  for (int ir=0; ir<namplitudes; ++ir)
    amps.push_back(amp[ir]);
    
  return;    
}


void MG_gg2ttxg::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate all wavefunctions
  vxxxxx(p[perm[2]], mME[0], hel[0], +1, w[0]); 
  vxxxxx(p[perm[3]], mME[1], hel[1], +1, w[1]); 
  oxxxxx(p[perm[0]], mME[2], hel[2], +1, w[2]); 
  ixxxxx(p[perm[1]], mME[3], hel[3], -1, w[3]); 
  vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]); 

  VVV1P0_1(w[0], w[1], pars->GC_10, pars->ZERO, pars->ZERO, w[5]); 
  FFV1P0_3(w[3], w[2], pars->GC_11, pars->ZERO, pars->ZERO, w[6]); 
  FFV1_1(w[2], w[4], pars->GC_11, pars->MT, pars->ZERO, w[7]); 
  FFV1_2(w[3], w[4], pars->GC_11, pars->MT, pars->ZERO, w[8]); 
  FFV1_1(w[2], w[0], pars->GC_11, pars->MT, pars->ZERO, w[9]); 
  FFV1_2(w[3], w[1], pars->GC_11, pars->MT, pars->ZERO, w[10]); 
  VVV1P0_1(w[1], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[11]); 
  FFV1_2(w[3], w[0], pars->GC_11, pars->MT, pars->ZERO, w[12]); 
  FFV1_1(w[2], w[1], pars->GC_11, pars->MT, pars->ZERO, w[13]); 
  VVV1P0_1(w[0], w[4], pars->GC_10, pars->ZERO, pars->ZERO, w[14]); 
  VVVV1P0_1(w[0], w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[15]); 
  VVVV3P0_1(w[0], w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[16]); 
  VVVV4P0_1(w[0], w[1], w[4], pars->GC_12, pars->ZERO, pars->ZERO, w[17]); 

  // Calculate all amplitudes
  VVV1_0(w[5], w[6], w[4], pars->GC_10, amp[0]); 
  FFV1_0(w[3], w[7], w[5], pars->GC_11, amp[1]); 
  FFV1_0(w[8], w[2], w[5], pars->GC_11, amp[2]); 
  FFV1_0(w[10], w[9], w[4], pars->GC_11, amp[3]); 
  FFV1_0(w[3], w[9], w[11], pars->GC_11, amp[4]); 
  FFV1_0(w[8], w[9], w[1], pars->GC_11, amp[5]); 
  FFV1_0(w[12], w[13], w[4], pars->GC_11, amp[6]); 
  FFV1_0(w[12], w[2], w[11], pars->GC_11, amp[7]); 
  FFV1_0(w[12], w[7], w[1], pars->GC_11, amp[8]); 
  FFV1_0(w[3], w[13], w[14], pars->GC_11, amp[9]); 
  FFV1_0(w[10], w[2], w[14], pars->GC_11, amp[10]); 
  VVV1_0(w[14], w[1], w[6], pars->GC_10, amp[11]); 
  FFV1_0(w[8], w[13], w[0], pars->GC_11, amp[12]); 
  FFV1_0(w[10], w[7], w[0], pars->GC_11, amp[13]); 
  VVV1_0(w[0], w[11], w[6], pars->GC_10, amp[14]); 
  FFV1_0(w[3], w[2], w[15], pars->GC_11, amp[15]); 
  FFV1_0(w[3], w[2], w[16], pars->GC_11, amp[16]); 
  FFV1_0(w[3], w[2], w[17], pars->GC_11, amp[17]); 

}
