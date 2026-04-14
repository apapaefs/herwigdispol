// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HwMEBase class.
//

#include "HwMEBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/PDT/GenericMassGenerator.h"
#include "ThePEG/Cuts/Cuts.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "ThePEG/PDF/PolarizedPartonExtractor.h"

using namespace Herwig;

void HwMEBase::persistentOutput(PersistentOStream & os) const {
  os << massOption_ << rescaleOption_;
}

void HwMEBase::persistentInput(PersistentIStream & is, int) {
  is >> massOption_ >> rescaleOption_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<HwMEBase,MEBase>
describeHerwigHwMEBase("Herwig::HwMEBase", "Herwig.so");

void HwMEBase::Init() {

  static ClassDocumentation<HwMEBase> documentation
    ("The HwMEBase class is the base class for matrix elements in Herwig"
     " and provides the virtual members for hard radiation corrections in the"
     " shower.");

}

int HwMEBase::nDim() const {
  unsigned ndim = 1;
  for(unsigned int ix=0;ix<massOption_.size();++ix)
    if(massOption_[ix]==2) ++ndim;
  return ndim;
}

CrossSection HwMEBase::dSigHatDR() const {
  return me2()*jacobian()/(16.0*sqr(Constants::pi)*sHat())*sqr(hbarc);
}

void HwMEBase::setKinematics() {
  MEBase::setKinematics();
  lastTHat_ = (meMomenta()[0] - meMomenta()[2]).m2();
  lastUHat_ = (meMomenta()[1] - meMomenta()[2]).m2();
  lastPhi_ = meMomenta()[2].phi();
}

bool HwMEBase::generateMasses(vector<Energy> & masses, double & mjac,
			      const double *r) {
  assert(massOption_.size()+2==mePartonData().size());
  mjac = 1.;
  masses.clear();
  masses.resize(massOption_.size(),ZERO);
  Energy ecm = sqrt(sHat());
  Energy emin(ZERO);
  int noff(0);
  for(unsigned int ix=0;ix<massOption_.size();++ix) {
    if(massOption_[ix]==1) {
      masses[ix] = mePartonData()[ix+2]->hardProcessMass();
      emin += masses[ix];
    }
    else if (massOption_[ix]==2) {
      emin += mePartonData()[ix+2]->massMin();
      ++noff;
    }
  }
  // check allowed
  if(emin>ecm) return false;
  // if nothing off-shell return
  if(noff==0) return true;
  int iloc = nDim()-noff;
  emin = ecm - emin;
  // generate the masses
  for(unsigned int ix=0;ix<massOption_.size();++ix) {
    if(massOption_[ix]!=2) continue;
    Energy mmin = mePartonData()[ix+2]->massMin();
    emin += mmin;
    Energy mmax = min(mePartonData()[ix+2]->massMax(),emin);
    if(mmin>mmax) return false;
    tGenericMassGeneratorPtr gen = mePartonData()[ix+2]->massGenerator() ?
      dynamic_ptr_cast<tGenericMassGeneratorPtr>(mePartonData()[ix+2]->massGenerator()) :
      tGenericMassGeneratorPtr();
    if(gen) {
      double jtemp(0.);
      masses[ix] = gen->mass(jtemp,*mePartonData()[ix+2],mmin,mmax,r[iloc]);
      mjac *= jtemp;
    }
    else {
      Energy mon(mePartonData()[ix+2]->hardProcessMass());
      Energy width(mePartonData()[ix+2]->width());
      double rhomin = atan2((sqr(mmin)-sqr(mon)), mon*width);
      double rhomax = atan2((sqr(mmax)-sqr(mon)), mon*width);
      masses[ix] = sqrt(mon*width*tan(rhomin+r[iloc]*(rhomax-rhomin))+sqr(mon));
      mjac *= (rhomax-rhomin)/Constants::pi;
    }
    emin -= masses[ix];
    if(emin<ZERO) return false;
    ++iloc;
  }
  return true;
}

bool HwMEBase::setupTwoToTwoKinematics(const double * r,
                                       TwoToTwoKinematicsSetup & setup) {
  jacobian(1.);
  if(!generateMasses(setup.masses, setup.mjac, r)) return false;
  // set up the momenta
  for ( int i = 2, N = meMomenta().size(); i < N; ++i ) {
    meMomenta()[i] = Lorentz5Momentum(setup.masses[i-2]);
  }
  setup.ctmin = -1.0;
  setup.ctmax = 1.0;
  try {
    setup.q = SimplePhaseSpace::
      getMagnitude(sHat(), meMomenta()[2].mass(), meMomenta()[3].mass());
  } 
  catch ( ImpossibleKinematics & e) {
    return false;
  }

  setup.e = sqrt(sHat())/2.0;
     	    
  setup.m22 = meMomenta()[2].mass2();
  setup.m32 = meMomenta()[3].mass2();
  setup.e0e2 = 2.0*setup.e*sqrt(sqr(setup.q) + setup.m22);
  setup.e1e2 = 2.0*setup.e*sqrt(sqr(setup.q) + setup.m22);
  setup.e0e3 = 2.0*setup.e*sqrt(sqr(setup.q) + setup.m32);
  setup.e1e3 = 2.0*setup.e*sqrt(sqr(setup.q) + setup.m32);
  setup.pq = 2.0*setup.e*setup.q;

  Energy2 thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[2]);
  if ( thmin > ZERO )
    setup.ctmax = min(setup.ctmax, (setup.e0e2 - setup.m22 - thmin)/setup.pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[2]);
  if ( thmin > ZERO )
    setup.ctmin = max(setup.ctmin, (thmin + setup.m22 - setup.e1e2)/setup.pq);

  thmin = lastCuts().minTij(mePartonData()[1], mePartonData()[3]);
  if ( thmin > ZERO )
    setup.ctmax = min(setup.ctmax, (setup.e1e3 - setup.m32 - thmin)/setup.pq);

  thmin = lastCuts().minTij(mePartonData()[0], mePartonData()[3]);
  if ( thmin > ZERO )
    setup.ctmin = max(setup.ctmin, (thmin + setup.m32 - setup.e0e3)/setup.pq);

  Energy ptmin = max(lastCuts().minKT(mePartonData()[2]),
   		     lastCuts().minKT(mePartonData()[3]));
  if ( ptmin > ZERO ) {
    double ctm = 1.0 - sqr(ptmin/setup.q);
    if ( ctm <= 0.0 ) return false;
    setup.ctmin = max(setup.ctmin, -sqrt(ctm));
    setup.ctmax = min(setup.ctmax, sqrt(ctm));
  }

  double ymin2 = lastCuts().minYStar(mePartonData()[2]);
  double ymax2 = lastCuts().maxYStar(mePartonData()[2]);
  double ymin3 = lastCuts().minYStar(mePartonData()[3]);
  double ymax3 = lastCuts().maxYStar(mePartonData()[3]);
  double ytot = lastCuts().Y() + lastCuts().currentYHat();
  if ( ymin2 + ytot > -0.9*Constants::MaxRapidity )
    setup.ctmin = max(setup.ctmin, sqrt(sqr(setup.q) + setup.m22)*tanh(ymin2)/setup.q);
  if ( ymax2 + ytot < 0.9*Constants::MaxRapidity )
    setup.ctmax = min(setup.ctmax, sqrt(sqr(setup.q) + setup.m22)*tanh(ymax2)/setup.q);
  if ( ymin3 + ytot > -0.9*Constants::MaxRapidity )
    setup.ctmax = min(setup.ctmax, sqrt(sqr(setup.q) + setup.m32)*tanh(-ymin3)/setup.q);
  if ( ymax3 + ytot < 0.9*Constants::MaxRapidity )
    setup.ctmin = max(setup.ctmin, sqrt(sqr(setup.q) + setup.m32)*tanh(-ymax3)/setup.q);
  
  return setup.ctmin < setup.ctmax;
}

bool HwMEBase::finishTwoToTwoKinematics(const TwoToTwoKinematicsSetup & setup,
                                        double cth,
                                        bool * cutsRejected) {
  if ( cutsRejected ) *cutsRejected = false;
    
  Energy pt = setup.q*sqrt(1.0-sqr(cth));
  phi(rnd(2.0*Constants::pi));
  meMomenta()[2].setVect(Momentum3( pt*sin(phi()),  pt*cos(phi()),  setup.q*cth));
  meMomenta()[3].setVect(Momentum3(-pt*sin(phi()), -pt*cos(phi()), -setup.q*cth));

  meMomenta()[2].rescaleEnergy();
  meMomenta()[3].rescaleEnergy();

  vector<LorentzMomentum> out(2);
  out[0] = meMomenta()[2];
  out[1] = meMomenta()[3];
  tcPDVector tout(2);
  tout[0] = mePartonData()[2];
  tout[1] = mePartonData()[3];
  if ( !lastCuts().passCuts(tout, out, mePartonData()[0], mePartonData()[1]) ) {
    if ( cutsRejected ) *cutsRejected = true;
    return false;
  }

  tHat(setup.pq*cth + setup.m22 - setup.e0e2);
  uHat(setup.m22 + setup.m32 - sHat() - tHat());
  jacobian((setup.pq/sHat())*Constants::pi*jacobian()*setup.mjac);
  // compute the rescaled momenta
  return rescaleMomenta(meMomenta(),mePartonData());
}

bool HwMEBase::generateKinematics(const double * r) {
  TwoToTwoKinematicsSetup setup;
  if(!setupTwoToTwoKinematics(r, setup)) return false;
  return finishTwoToTwoKinematics(setup,
                                  getCosTheta(setup.ctmin, setup.ctmax, r[0]));
}

bool HwMEBase::rescaleMomenta(const vector<Lorentz5Momentum> & momenta,
			      const cPDVector & data) {
  assert(momenta.size()==4&&data.size()==4);
  // default just use the ones we generated
  rescaledMomenta_=momenta;
  if(rescaleOption_==1) return true;
  Energy mnew[2] = {0*MeV, ZERO};
  if(rescaleOption_==0) {
    mnew[0] = ZERO;
    mnew[1] = ZERO;
  }
  else if(rescaleOption_==2) {
    mnew[0] = data[2]->hardProcessMass();
    mnew[1] = data[3]->hardProcessMass();
  }
  else if(rescaleOption_==3) {
    if(abs(data[2]->id())!=abs(data[3]->id())) return true;
    mnew[0] = 0.5*(momenta[2].mass()+momenta[3].mass());
    mnew[1] = mnew[0];
  } 
  else {
    assert(false);
  }
  Lorentz5Momentum pcm(momenta[2]+momenta[3]);
  Energy m0=pcm.m();
  if(m0<mnew[0]+mnew[1]) return false;
  Boost bv = pcm.boostVector();
  rescaledMomenta_[2].boost(bv);
  rescaledMomenta_[2].setMass(mnew[0]);
  rescaledMomenta_[2].setE(0.5*(sqr(m0)+sqr(mnew[0])-sqr(mnew[1]))/m0);
  if(rescaledMomenta_[2].t()-rescaledMomenta_[2].mass()>1e-10*(rescaledMomenta_[2].t()+rescaledMomenta_[2].mass()))
    rescaledMomenta_[2].rescaleRho();
  else {
    rescaledMomenta_[2].setX(ZERO);
    rescaledMomenta_[2].setY(ZERO);
    rescaledMomenta_[2].setZ(ZERO);
  }
  rescaledMomenta_[2].boost(-bv);
  rescaledMomenta_[3].boost(bv);
  rescaledMomenta_[3].setMass(mnew[1]);
  rescaledMomenta_[3].setE(0.5*(sqr(m0)-sqr(mnew[0])+sqr(mnew[1]))/m0);
  if(rescaledMomenta_[3].t()-rescaledMomenta_[3].mass()>1e-10*(rescaledMomenta_[3].t()+rescaledMomenta_[3].mass()))
    rescaledMomenta_[3].rescaleRho();
  else {
    rescaledMomenta_[3].setX(ZERO);
    rescaledMomenta_[3].setY(ZERO);
    rescaledMomenta_[3].setZ(ZERO);
  }
  rescaledMomenta_[3].boost(-bv);
  return true;
}

double HwMEBase::getCosTheta(double ctmin, double ctmax, const double r) {
  double cth = 0.0;
  static const double eps = 1.0e-6;
  if ( 1.0 + ctmin <= eps && 1.0 - ctmax <= eps ) {
    jacobian(jacobian()*(ctmax - ctmin));
    cth = ctmin + r*(ctmax - ctmin);
  } else if (  1.0 + ctmin <= eps ) {
    cth = 1.0 - (1.0 - ctmax)*pow((1.0 - ctmin)/(1.0 - ctmax), r);
    jacobian(jacobian()*log((1.0 - ctmin)/(1.0 - ctmax))*(1.0 - cth));
  } else if (  1.0 - ctmax <= eps ) {
    cth = -1.0 + (1.0 + ctmin)*pow((1.0 + ctmax)/(1.0 + ctmin), r);
    jacobian(jacobian()*log((1.0 + ctmax)/(1.0 + ctmin))*(1.0 + cth));
  } else {
    double zmin = 0.5*(1.0 - ctmax);
    double zmax = 0.5*(1.0 - ctmin);
    double A1 = -ctmin/(zmax*(1.0-zmax));
    double A0 = -ctmax/(zmin*(1.0-zmin));
    double A = r*(A1 - A0) + A0;
    double z = A < 2.0? 2.0/(sqrt(sqr(A) + 4.0) + 2 - A):
      0.5*(A - 2.0 + sqrt(sqr(A) + 4.0))/A;
    cth = 1.0 - 2.0*z;
    jacobian(jacobian()*2.0*(A1 - A0)*sqr(z)*sqr(1.0 - z)/(sqr(z) + sqr(1.0 - z)));
  }
  return cth;
}

bool HwMEBase::softMatrixElementVeto(PPtr, PPtr,
				     const bool &,
				     const Energy & ,
				     const vector<tcPDPtr> & ,
				     const double & ,
				     const Energy & ,
				     const Energy & ) {
  assert(false);
  return false;
}

RealEmissionProcessPtr HwMEBase::generateHardest(RealEmissionProcessPtr,ShowerInteraction) {
  assert(false);
  return RealEmissionProcessPtr();
}

RealEmissionProcessPtr HwMEBase::applyHardMatrixElementCorrection(RealEmissionProcessPtr) {
  assert(false);
  return RealEmissionProcessPtr();
}

void HwMEBase::initializeMECorrection(RealEmissionProcessPtr , double & ,
				      double & ) {
  assert(false);
}

pair<RhoDMatrix,RhoDMatrix> HwMEBase::getRhoMatrices() const {
  // default unpolarized output
  pair<RhoDMatrix,RhoDMatrix> output = make_pair(RhoDMatrix(mePartonData()[0]->iSpin()), 
                                                 RhoDMatrix(mePartonData()[1]->iSpin()));
  // if already determined not polarized return
  if(!_polarized) return output;
  // check if polarized extractor and if not return
  Ptr<PolarizedPartonExtractor>::tptr pExtractor = dynamic_ptr_cast<Ptr<PolarizedPartonExtractor>::tptr>(lastExtractor());
  if(! pExtractor) {
    _polarized = false;
    return output;
  }
  output = pExtractor->getRhoMatrices(lastXComb().partonBinInstances());
  if(lastXComb().mirror()) swap(output.first,output.second);
  return output;
}

ThreeVector<double> HwMEBase::getBeamPolarization(bool firstOrSecond) const {
  // default unpolarized output
  ThreeVector<double> output = {0,0,0};
  // if already determined not polarized return
  if(!_polarized) return output;
  // check if polarized extractor and if not return
  Ptr<PolarizedPartonExtractor>::tptr pExtractor = dynamic_ptr_cast<Ptr<PolarizedPartonExtractor>::tptr>(lastExtractor());
  if(! pExtractor) {
    _polarized = false;
    return output;
  }
  output = pExtractor->polarization(firstOrSecond);
  return output;
}
