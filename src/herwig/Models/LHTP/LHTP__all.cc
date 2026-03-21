#line 1 "./LHTPModel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPModel class.
//

#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "gsl/gsl_multiroots.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include <algorithm>

using namespace Herwig;
using namespace ThePEG;

// equations for top parameters for GSL
namespace {

// struct to provide the model parameters to the functions
struct tparams {
  Energy v;
  Energy f;
  Energy mt;
  double tan2a;
};

// equations defining tan 2alpha and mt expressed in form f(lambda1,lambda2)=0
// to be solved to give lambda1 and lambda2 by gsl
int top_equation(const gsl_vector * x, void *params, gsl_vector *f ) {
  // yukawa and check top mass
  const double lam1 = gsl_vector_get(x,0);
  const double lam2 = gsl_vector_get(x,1);
  Energy fs = ((struct tparams *) params)->f;
  Energy v = ((struct tparams *) params)->v;
  double sv = sin(sqrt(2.)*v/fs);
  double cv = cos(sqrt(2.)*v/fs);
  Energy mt    = ((struct tparams *) params)->mt;
  double tan2a = ((struct tparams *) params)->tan2a;
  double f1 = 4.*lam1*lam2*(1.+cv)/(4.*sqr(lam2)-sqr(lam1)*(2.*sqr(sv)+sqr(1.+cv)))
    -tan2a;
  double delta = 0.5*(sqr(lam2)+0.5*sqr(lam1)*(sqr(sv)+0.5*sqr(1.+cv)));
  double f2 = sqr(fs/mt)*delta*(1.-sqrt(1.-0.5*sqr(lam1*lam2*sv/delta)))-1.;
  if(lam1*lam2<0.) f1+=1e10;
  if(lam1*lam2<0.) f2+=1e10;
  gsl_vector_set(f,0,f1);
  gsl_vector_set(f,1,f2);
  return GSL_SUCCESS;
}

}

LHTPModel::LHTPModel()
  : f_(0.5*TeV), salpha_(sqrt(0.5)), calpha_(sqrt(0.5)), sbeta_(0.), cbeta_(0.), 
    sL_(0.), cL_(1.), sR_(0.), cR_(0.),
    kappaQuark_(1.), kappaLepton_(1.), mh_(125.*GeV), v_(246.*GeV),
    g_(sqrt(0.43)), gp_(sqrt(0.12)), approximate_(false)
{}

IBPtr LHTPModel::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPModel::fullclone() const {
  return new_ptr(*this);
}

void LHTPModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(f_,TeV) << salpha_ << calpha_ << sbeta_ << cbeta_
     << kappaQuark_ << kappaLepton_ << ounit(v_,GeV) 
     << g_ << gp_ << sthetaH_ << cthetaH_ << approximate_
     << sL_ << cL_ << sR_ << cR_ << WHHVertex_;
}

void LHTPModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(f_,TeV) >> salpha_ >> calpha_ >> sbeta_ >> cbeta_
     >> kappaQuark_ >> kappaLepton_>> iunit(v_,GeV)
     >> g_ >> gp_ >> sthetaH_ >> cthetaH_ >> approximate_
     >> sL_ >> cL_ >> sR_ >> cR_ >> WHHVertex_;
}


// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPModel,StandardModel>
describeHerwigLHTPModel("Herwig::LHTPModel", "HwLHTPModel.so");

void LHTPModel::Init() {

  static ClassDocumentation<LHTPModel> documentation
    ("The LHTPModel class implements the Little Higgs model"
     " with T-parity");

  static Parameter<LHTPModel,Energy> interfacef
    ("f",
     "The scale of the non-linear sigma-model",
     &LHTPModel::f_, TeV, 1.*TeV, 0.0*TeV, 10.0*TeV,
     true, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceSinAlpha
    ("SinAlpha",
     "The parameter controlling the mixing in the top quark sector of the model",
     &LHTPModel::salpha_, sqrt(0.5), 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceKappaQuark
    ("KappaQuark",
     "The parameter controlling the masses of the T-odd quarks",
     &LHTPModel::kappaQuark_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,double> interfaceKappaLepton
    ("KappaLepton",
     "The parameter controlling the masses of the T-odd leptons",
     &LHTPModel::kappaLepton_, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHTPModel,Energy> interfaceHiggsMass
    ("HiggsMass",
     "The mass of the lightest Higgs boson",
     &LHTPModel::mh_, GeV, 120.0*GeV, 100.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Switch<LHTPModel,bool> interfaceApproximate
    ("Approximate",
     "Whether to use the full expression for the mases of the top quark"
     " and its partners or the second-order expansion in v/f.",
     &LHTPModel::approximate_, false, false, false);
  static SwitchOption interfaceApproximateYes
    (interfaceApproximate,
     "Yes",
     "Approximate",
     true);
  static SwitchOption interfaceApproximateNo
    (interfaceApproximate,
     "No",
     "Don't approximate",
     false);

  static Reference<LHTPModel,AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/WHH",
     "Vertex for the interactions of the electroweak gauge"
     " bosons and two Higgs bosons.",
     &LHTPModel::WHHVertex_, false, false, true, false, false);

}

void LHTPModel::doinit() {
  addVertex(WHHVertex_);
  StandardModel::doinit();
  using Constants::pi;
  // compute the parameters of the model
  // W and Z masses
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  // Energy mz(getParticleData(ParticleID::Z0)->mass());
  // couplings g and g'
  // double ee = sqrt(4.*pi*alphaEM(sqr(mz)));
  double ee = sqrt(4.*pi*alphaEMMZ());
  double sw(sqrt(sin2ThetaW())),cw(sqrt(1.-sin2ThetaW()));
  g_  = ee/sw;
  gp_ = ee/cw;
  // vev
  v_ = 2.*mw/g_;
  double vf(sqr(v_/f_));
  // calculate masses of the new particles from input
  // and SM parameters
  // masses of the new gauge bosons (MWH = MZH)
  Energy MAH = gp_*f_*sqrt(0.2)*(1.-0.625*vf);
  Energy MZH = g_ *f_*          (1.-0.125*vf);
  // mixings
  sthetaH_ = 1.25*g_*gp_/(5.*sqr(g_)-sqr(gp_))*vf;
  cthetaH_ = sqrt(1.-sqr(sthetaH_));
  // masses of the new top quarks
  Energy MTp,MTm;
  topMixing(MTp,MTm);
  // mixings in the top sector
  sL_ = sqr(salpha_)*v_/f_;
  cL_ = sqrt(1.-sqr(sL_));
  sR_ = salpha_*(1.-0.5*sqr(calpha_)*(sqr(calpha_)-sqr(salpha_))*vf);
  cR_ = sqrt(1.-sqr(sR_));
  // masses of the T-odd fermions
  Energy Mdm = sqrt(2.)*kappaQuark_ *f_;
  Energy Mum = sqrt(2.)*kappaQuark_ *f_*(1.-0.125*vf);
  Energy Mlm = sqrt(2.)*kappaLepton_*f_;
  Energy Mnm = sqrt(2.)*kappaLepton_*f_*(1.-0.125*vf);
  // masses of the triplet higgs
  Energy MPhi = sqrt(2.)*mh_*f_/v_;
  // set the masses of the new particles
  // new gauge bosons
  resetMass( 32      , MAH  );
  resetMass( 33      , MZH  );
  resetMass( 34      , MZH  );
  resetMass(-34      , MZH  );
  // masses of the new top quarks
  resetMass(  8      , MTp  );
  resetMass( -8      , MTp  );
  resetMass(  4000008, MTm  );
  resetMass( -4000008, MTm  );
  //  masses of the Higgs bosons
  resetMass( 25      , mh_  );
  resetMass( 35      , MPhi );
  resetMass( 36      , MPhi );
  resetMass( 37      , MPhi );
  resetMass(-37      , MPhi );
  resetMass( 38      , MPhi );
  resetMass(-38      , MPhi );
  // masses of the T-odd quarks
  resetMass( 4000001, Mdm  );
  resetMass(-4000001, Mdm  );
  resetMass( 4000002, Mum  );
  resetMass(-4000002, Mum  );
  resetMass( 4000003, Mdm  );
  resetMass(-4000003, Mdm  );
  resetMass( 4000004, Mum  );
  resetMass(-4000004, Mum  );
  resetMass( 4000005, Mdm  );
  resetMass(-4000005, Mdm  );
  resetMass( 4000006, Mum  );
  resetMass(-4000006, Mum  );
  // masses of the T-odd leptons
  resetMass( 4000011, Mlm  );
  resetMass(-4000011, Mlm  );
  resetMass( 4000012, Mnm  );
  resetMass(-4000012, Mnm  );
  resetMass( 4000013, Mlm  );
  resetMass(-4000013, Mlm  );
  resetMass( 4000014, Mnm  );
  resetMass(-4000014, Mnm  );
  resetMass( 4000015, Mlm  );
  resetMass(-4000015, Mlm  );
  resetMass( 4000016, Mnm  );
  resetMass(-4000016, Mnm  );
}

void LHTPModel::topMixing(Energy & MTp, Energy & MTm) {
  double vf(sqr(v_/f_));
  Energy mt = getParticleData(ParticleID::t)->mass();
  calpha_ = sqrt(1.-sqr(salpha_));
  double sv(sin(sqrt(2.)*v_/f_)),cv(cos(sqrt(2.)*v_/f_));
  // first guess for Yukawa's based on second-order in v/f expansion
  double lambda1 = mt/v_/calpha_*(1.+(2.-3.*pow(salpha_,4))*vf/6.);
  double lambda2 = mt/v_/salpha_*(1.+(2.-3.*pow(calpha_,4))*vf/6.);
  // first guess for masses
  MTp = sqrt(sqr(lambda1)+sqr(lambda2))*f_*(1-0.5*vf*sqr(calpha_*salpha_));
  MTm = lambda2*f_;
  if(!approximate_) {
    // special case where denominator of tan 2 alpha eqn is zero
    if(abs(salpha_-sqrt(0.5))<1e-4) {
      double a = 0.25*(2.*sqr(sv)+sqr(1.+cv));
      double b = 0.5*(a+0.5*(sqr(sv)+0.5*sqr(1.+cv)));
      lambda1 = mt/f_*sqrt(1./b/(1.-sqrt(1.-0.5*a*sqr(sv/b))));
      lambda2 = sqrt(a)*lambda1;
    }
    // general case using GSL
    else {
      double ca = sqrt(1.-sqr(salpha_));
      double ta = salpha_/ca;
      double tan2a = 2.*ta/(1.-sqr(ta));
      const gsl_multiroot_fsolver_type *T;
      gsl_multiroot_fsolver *s;
      int status;
      size_t iter=0;
      const size_t n=2;
      struct tparams p = {v_,f_,mt,tan2a};
      gsl_multiroot_function f = {&top_equation, n, &p};
      gsl_vector *x = gsl_vector_alloc(n);
      gsl_vector_set(x,0,lambda1);
      gsl_vector_set(x,1,lambda2);
      T = gsl_multiroot_fsolver_hybrids;
      s = gsl_multiroot_fsolver_alloc(T,2);
      gsl_multiroot_fsolver_set(s, &f,x);
      do {
	iter++;
	status = gsl_multiroot_fsolver_iterate(s);
	if(status) break;
	status = gsl_multiroot_test_residual(s->f,1e-7);
      }
      while (status==GSL_CONTINUE && iter < 1000);
      gsl_multiroot_fsolver_free(s);
      lambda1 = gsl_vector_get(s->x,0);
      lambda2 = gsl_vector_get(s->x,1);
      gsl_vector_free(x);
    }
    // calculate the heavy top masses using full result
    double delta = 0.5*(sqr(lambda2)+0.5*sqr(lambda1)*(sqr(sv)+0.5*sqr(1.+cv)));
    double det = sqrt(1.-0.5*sqr(lambda1*lambda2*sv/delta));
    MTp = sqrt(sqr(f_)*delta*(1.+det));
    MTm = lambda2*f_;
  }
  // beta mixing angle
  double beta = 0.5*atan(2.*sqrt(2.)*sqr(lambda1)*sv*(1.+cv)/
			 (4.*sqr(lambda2)+sqr(1.+cv)*sqr(lambda1)-2.*sqr(lambda1)*sv));
  sbeta_ = sin(beta);
  cbeta_ = cos(beta);
}
#line 1 "./LHTPWWHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWWHVertex class.
//

#include "LHTPWWHVertex.h"
#include "LHTPModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPWWHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPWWHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(coup_,GeV);
}

void LHTPWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coup_,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPWWHVertex,VVSVertex>
describeHerwigLHTPWWHVertex("Herwig::LHTPWWHVertex", "HwLHTPModel.so");

void LHTPWWHVertex::Init() {

  static ClassDocumentation<LHTPWWHVertex> documentation
    ("The LHTPWWHVertex class implements the coupling of two electroweak"
     " gauge bosons to a Higgs boson in the Little Higgs Model with T-Parity"
     "including the additional heavy photon, Z and W bosons and the "
     "triplet Higgs bosons.");

}

LHTPWWHVertex::LHTPWWHVertex() : coupLast_(0.), q2Last_(0.*GeV2) {
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHTPWWHVertex::doinit() {
  // W_L W_L H
  addToList(  24,  -24,    25);
  // Z_L Z_L H
  addToList(  23,   23,    25);
  // W_H W_H H
  addToList(  34,  -34,    25);
  // Z_H Z_H H
  addToList(  33,   33,    25);
  // A_H A_H H
  addToList(  32,   32,    25);
  // Z_H A_H H
  addToList(  33,   32,    25);
  // Z_L Z_H Phi0
  addToList(  23,   33,    35);
  // A_H Z_L Phi0
  addToList(  32,   23,    35);
  // W_L W_H PhiP
  addToList(  24,  -34,    36);
  addToList(  34,  -24,    36);
  // W_H Z_L Phi+/- 
  addToList(  34,   23,   -37);
  addToList( -34,   23,    37);
  // W_L A_H Phi+/-
  addToList(  24,   32,   -37);
  addToList( -24,   32,    37);
  // W_L Z_H Phi+/- 
  addToList(  24,   33,   -37);
  addToList( -24,   33,    37);
  // W_H A_L Phi+/- 
  addToList(  34,   22,   -37);
  addToList( -34,   22,    37);
  // W_L W_H Phi --/++
  addToList(  24,   34,   -38);
  addToList( -24,  -34,    38);
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHTPModel "
			  << " in LHTPWWHVertex::doinit()"
			  << Exception::runerror;
  // base class
  VVSVertex::doinit();
  // calculate the couplings for the different combinations of particles
  Energy fact = 0.5*model->vev()/model->sin2ThetaW();
  double sw(sqrt(model->sin2ThetaW())),cw(sqrt(1.-model->sin2ThetaW()));
  double vf(model->vev()/model->f());
  double r2(sqrt(2.));
  coup_.resize(14);
  // H 
  coup_[ 0] = fact        *(1.-sqr(vf)/3.);
  coup_[ 1] = fact/sqr(cw)*(1.-sqr(vf)/3.);
  coup_[ 2] =-fact;
  coup_[ 3] =-fact;
  coup_[ 4] =-fact*sqr(sw/cw);
  coup_[ 5] =-fact/cw*sw;
  // PhiP 
  coup_[ 6] = r2*fact*vf/3.;
  // Phi0
  coup_[ 7] =-fact*vf/r2/cw;
  coup_[ 8] = fact*vf/r2*sw/sqr(cw);
  // Phi+
  coup_[ 9] = fact*vf/6./cw*(1.+2.*sqr(sw));
  coup_[10] = fact*vf*sw/cw*0.5;
  coup_[11] = fact*vf*5./6.;
  coup_[12] =-fact*vf*sw/3.;
  // Phi++
  coup_[13] =-fact*vf;
}

void LHTPWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,
					      tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = sqr(electroMagneticCoupling(q2));
    q2Last_=q2;
  }
  long ih = abs(c->id());
  long ibos[2]={abs(a->id()),abs(b->id())};
  
  if(ih == 25) {
    if( ibos[0] == 24 && ibos[1] == 24) 
      norm(UnitRemoval::InvE *coupLast_*coup_[0]);
    else if( ibos[0] == 23 && ibos[1] == 23     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[1]);
    else if( ibos[0] == 34 && ibos[1] == 34     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[2]);
    else if( ibos[0] == 33 && ibos[1] == 33     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[3]);
    else if( ibos[0] == 32 && ibos[1] == 32     ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[4]);
    else if((ibos[0] == 33 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 33)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[5]);
    else 
      assert(false);
  }
  else if(ih == 36) {
    if( a->id() == 34 || b->id() == 34) 
      norm( Complex(0.,1.)*UnitRemoval::InvE *coupLast_*coup_[6]);
    else
      norm(-Complex(0.,1.)*UnitRemoval::InvE *coupLast_*coup_[6]);
  }
  else if(ih == 35) {
    if((ibos[0] == 23 && ibos[1] == 33) ||
       (ibos[0] == 33 && ibos[1] == 23)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[7]);
    else if((ibos[0] == 23 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 23)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[8]);
    else 
      assert(false);
  }
  else if(ih == 37) {
    if((ibos[0] == 34 && ibos[1] == 23) ||
       (ibos[0] == 23 && ibos[1] == 34)    ) {
      norm(UnitRemoval::InvE *coupLast_*coup_[ 9]);
    }
    else if((ibos[0] == 24 && ibos[1] == 32) ||
	    (ibos[0] == 32 && ibos[1] == 24)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[10]);
    else if((ibos[0] == 24 && ibos[1] == 33) ||
	    (ibos[0] == 33 && ibos[1] == 24)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[11]);
    else if((ibos[0] == 34 && ibos[1] == 22) ||
	    (ibos[0] == 22 && ibos[1] == 34)    ) 
      norm(UnitRemoval::InvE *coupLast_*coup_[12]);
    else
      assert(false);
  }
  else if(ih == 38) {
    norm(UnitRemoval::InvE *coupLast_*coup_[13]);
  }
  else
    assert(false);
}
#line 1 "./LHTPFFGVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFGVertex class.
//

#include "LHTPFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPFFGVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFGVertex::fullclone() const {
  return new_ptr(*this);
}

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<LHTPFFGVertex,FFVVertex>
describeHerwigLHTPFFGVertex("Herwig::LHTPFFGVertex", "HwLHTPModel.so");

void LHTPFFGVertex::Init() {

  static ClassDocumentation<LHTPFFGVertex> documentation
    ("The LHTPFFGVertex class implements the couples of the fermions "
     "to the gluons in the Little Higgs model with T-parity.");

}

LHTPFFGVertex::LHTPFFGVertex() 
  : coupLast_(0.), q2Last_(0.*GeV2) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

void LHTPFFGVertex::doinit() {
  // SM quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix, 21);
  }
  // additional top quark
  addToList(-8,  8, 21);
  // T odd quarks
  for(long ix = 4000001; ix <= 4000006; ++ix) {
    addToList(-ix,     ix, 21);
  }
  addToList(-4000008,   4000008, 21);
  FFVVertex::doinit();
}

// coupling for FFG vertex
void LHTPFFGVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -strongCoupling(q2);
    q2Last_=q2;
  }
  norm(coupLast_);
  // the left and right couplings
  int iferm=abs(a->id());
  if( iferm > 8 ) iferm -= 4000000;
  if((iferm>=1 && iferm<=8)) {
    left (1.);
    right(1.);
  }
  else
    assert(false);
}
#line 1 "./LHTPFFPVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFPVertex class.
//

#include "LHTPFFPVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "LHTPModel.h"

using namespace Herwig;

IBPtr LHTPFFPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFPVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << charge_ << coupd_ << coupu_ << coupe_ << coupnu_ 
     << TPreFactor_ << sL_ << cL_ << sR_ << cR_;
}

void LHTPFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> charge_ >> coupd_ >> coupu_ >> coupe_ >> coupnu_ 
     >> TPreFactor_ >> sL_ >> cL_ >> sR_ >> cR_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFPVertex,FFVVertex>
describeHerwigLHTPFFPVertex("Herwig::LHTPFFPVertex", "HwLHTPModel.so");

void LHTPFFPVertex::Init() {

  static ClassDocumentation<LHTPFFPVertex> documentation
    ("The LHTPFFPVertex class implements the coupling"
     " of the charged fermions to the photon in the Little Higgs"
     " model with T-parity.");

}

LHTPFFPVertex::LHTPFFPVertex() : 
  charge_(37,0.0), coupLast_(0.), q2Last_(-1.*GeV2),
  coupd_(0.), coupu_(0.), coupe_(0.), coupnu_(0.),
  TPreFactor_(0.), sL_(0.), cL_(1.), sR_(0.), cR_(1.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHTPFFPVertex::doinit() {
  // interactions with the photon
  // the quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix, 22);
  }
  // the leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix,    ix, 22);
  }
  // extra top quark
  addToList(-8,   8, 22);
  // the T-odd quarks
  for(long ix = 4000001;ix < 4000007; ++ix) {
    addToList(-ix,     ix, 22);
  }
  // the T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix,    ix, 22);
  }
  // extra top quark
  addToList(-4000008,  4000008, 22);
  // interactions with A_H
  // quarks and T-odd quark
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix - 4000000, ix          , 32);
    addToList(-ix          , ix + 4000000, 32);
  }
  // leptons and T-odd leptons (both charged leptons and neutrinos)
  for(int ix = 11; ix < 17; ++ix ) {
    addToList(-ix - 4000000, ix          , 32);
    addToList(-ix          , ix + 4000000, 32);
  }
  // T+T-A_H
  addToList(-4000008,       8, 32);
  addToList(      -8, 4000008, 32);
  // T-tA_H
  addToList(-4000008,       6, 32);
  addToList(      -6, 4000008, 32);
  // T-tA_H
  addToList(-4000006,       8, 32);
  addToList(      -8, 4000006, 32);
  // charges
  for(int ix = 1; ix < 16; ++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) charge_[ix] = double(ptemp->iCharge())/3.;
  }
  for(int ix = 4000001; ix < 4000016; ++ix) {
    tcPDPtr ptemp = getParticleData(ix);
    if(ptemp) charge_[ix-3999980] = double(ptemp->iCharge())/3.;
  }
  // couplings to A_H
  double sw = generator()->standardModel()->sin2ThetaW();
  double cw = sqrt(1.-sw);
  sw = sqrt(sw);
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFPVertex::doinit()"
				   << Exception::runerror;
  double cH = model->cosThetaH();
  double sH = model->sinThetaH();
  sL_ = model->sinThetaL();
  cL_ = model->cosThetaL();
  sR_ = model->sinThetaR();
  cR_ = model->cosThetaR();
  // couplings of fermion T-odd fermion A_H
  coupd_  = -0.1*(cH/cw-5.*sH/sw);
  coupu_  = -0.1*(cH/cw+5.*sH/sw);
  coupe_  = -0.1*(cH/cw-5.*sH/sw);
  coupnu_ = -0.1*(cH/cw+5.*sH/sw);
  // couplings of T+T- A_H
  TPreFactor_ = 0.4*cH/cw;
  FFVVertex::doinit();
}

// coupling for FFP vertex
void LHTPFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -electroMagneticCoupling(q2);
    q2Last_=q2;
  }
  norm(coupLast_);
  // the left and right couplings
  long iferm=abs(a->id()),ibos(c->id());
  if(ibos == ParticleID::gamma) {
    if(iferm < 20) {
      left (charge_[iferm]);
      right(charge_[iferm]);
    }
    else {
      iferm-=3999980;
      left (charge_[iferm]);
      right(charge_[iferm]);
    }
  }
  else if(ibos == 32) {
    long ianti = abs(b->id());
    if(iferm>4000000) swap(iferm,ianti);
    assert(iferm<4000000&&ianti>4000000);
    if( iferm == 6 || iferm == 8 ) {
      if     (iferm==6&&ianti==4000006) {
	left (cL_*coupu_);
	right(0.);
      }
      else if(iferm==6&&ianti==4000008) {
	left (-TPreFactor_*sL_);
	right(-TPreFactor_*sR_);
      }
      else if(iferm==8&&ianti==4000006) {
	left (sL_*coupu_);
	right(0.);
      }
      else if(iferm==8&&ianti==4000008) {
	left ( TPreFactor_*cL_);
	right( TPreFactor_*cR_);
      }
      else
	assert(false);
    }
    // quarks (inclding top)
    else if(iferm <= 5) {
      if(iferm % 2 == 0) left(coupu_);
      else               left(coupd_);
      right(0.);
    }
    // leptons
    else {
      if(iferm %2 == 0) left(coupnu_);
      else              left(coupe_ );
      right(0.);
    }
  }
  else
    assert(false);
}
#line 1 "./LHTPFFWVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFWVertex class.
//

#include "LHTPFFWVertex.h"
#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/CKMBase.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

IBPtr LHTPFFWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFWVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << ckm_ << sL_ << cL_;
}

void LHTPFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> ckm_ >> sL_ >> cL_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFWVertex,FFVVertex>
describeHerwigLHTPFFWVertex("Herwig::LHTPFFWVertex", "HwLHTPModel.so");

void LHTPFFWVertex::Init() {

  static ClassDocumentation<LHTPFFWVertex> documentation
    ("The LHTPFFWVertex class implements the couplings of the W"
     " and W_H bosons of the Little Higgss model with T-parity to the fermions.");

}

void LHTPFFWVertex::doinit() {
  // particles for outgoing W-
  // quarks
  for(int ix = 1; ix < 6; ix += 2) {
    for(int iy = 2; iy < 7; iy += 2) {
      addToList(-ix,      iy,      -24);
    }
  }
  // additional T quark
  addToList(-5,   8,  -24);
  // leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix,    ix + 1,    -24);
  }
  // T-odd quarks
  for(long ix = 4000002; ix < 4000007; ix += 2) {
    addToList(-ix + 1,    ix,    -24);
  }
  // T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix,    ix + 1,    -24);
  }
  // particles for outgoing W+
  // quarks
  for(int ix = 2; ix < 7; ix += 2) {
    for(int iy = 1; iy < 6; iy += 2) {
      addToList(-ix,       iy,       24);
    }
  }
  // additional T quark
  addToList(-8,   5,  24);
  // leptons
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix - 1,    ix,    24);
  }
  // T-odd quarks
  for(long ix = 4000002; ix < 4000009; ix += 2) {
    addToList(-ix,    ix-1,    24);
  }
  // T-odd leptons
  for(long ix = 4000011; ix < 4000017; ix += 2) {
    addToList(-ix-1,    ix,    24);
  }
  // particles for W_H-
  // quark and T-odd quark
  for(int ix = 1; ix < 6; ix += 2) {
    addToList(-ix-4000000,    ix+1,    -34);
    addToList(-ix,    ix+1+4000000,    -34);
  }
  addToList(-4000005,   8,    -34);
  // lepton and T-odd lepton
  for(int ix = 11;ix < 17; ix += 2) {
    addToList(-ix-4000000,    ix+1,    -34);
    addToList(-ix,    ix+1+4000000,    -34);
  }
  // particles for w_h+
  // quark and T-odd quark
  for(int ix = 1;ix < 6;ix += 2) {
    addToList(ix + 4000000,    -ix - 1,    34);
    addToList(ix,    -ix - 4000001,     34);
  }
  addToList(4000005,    -8,    34);
  // leptons and T-odd lepton
  for(int ix = 11; ix < 17; ix += 2) {
    addToList(-ix - 4000001,    ix,    34);
    addToList(-ix - 1,    ix + 4000000,    34);
  }
  ThePEG::Helicity::FFVVertex::doinit();
  Ptr<CKMBase>::transient_pointer CKM = generator()->standardModel()->CKM();
  // cast the CKM object to the HERWIG one
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    hwCKM = ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(CKM);
  if(hwCKM) {
    vector< vector<Complex > > CKM;
    CKM = hwCKM->getUnsquaredMatrix(generator()->standardModel()->families());
    for(unsigned int ix=0;ix<3;++ix) {
      for(unsigned int iy=0;iy<3;++iy) {
	ckm_[ix][iy]=CKM[ix][iy];
      }
    }
  }
  else {
    throw InitException() << "Must have access to the Herwig::StandardCKM object"
			  << "for the CKM matrix in LHTPFFWVertex::doinit()"
			  << Exception::runerror;
  }
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFWVertex::doinit()"
				   << Exception::runerror;
  sL_ = model->sinThetaL();
  cL_ = model->cosThetaL();
}

LHTPFFWVertex::LHTPFFWVertex()
  : sL_(0.), cL_(1.), ckm_(3,vector<Complex>(3,0.0)),
    coupLast_(0.),q2Last_(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

// coupling for FFW vertex
void LHTPFFWVertex::setCoupling(Energy2 q2, tcPDPtr a,
				tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -sqrt(0.5)*weakCoupling(q2);
    q2Last_   = q2;
  }
  norm(coupLast_);
  long ia(abs(a->id())),ib(abs(b->id()));
  // SM W boson
  if(abs(c->id())==ParticleID::Wplus) {
    // quarks
    if(ia >= 1 && ia <= 8 && ib >= 1 && ib <= 8 ) {
      int iu,id;
      // up type first
      if(ia % 2 == 0) {
	iu = ia/2;
	id = (ib+1)/2;
      }
      // down type first
      else {
	iu = ib/2;
	id = (ia+1)/2;
      }
      if(iu==4) iu=3;
      assert( iu>=1 && iu<=3 && id>=1 && id<=3);
      if      ( ia==6 || ib==6 ) left(ckm_[iu-1][id-1]*cL_);
      else if ( ia==8 || ib==8 ) left(ckm_[iu-1][id-1]*sL_);
      else                       left(ckm_[iu-1][id-1]    );
      right(0.);
    }
    // leptons
    else if( ia >= 11 && ia <= 16) {
      left(1.);
      right(0.);
    }
    else {
      left (1.);
      right(1.);
    }
  }
  else {
    if(ia==6||ib==6)      left(-cL_);
    else if(ia==8||ib==8) left(-sL_);
    else                  left(-1. );
    right(0.);
  }
}
#line 1 "./LHTPFFZVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFZVertex class.
//

#include "LHTPFFZVertex.h"
#include "LHTPModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr LHTPFFZVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPFFZVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << gl_ << gr_ << tl_ << tr_ << coupd_ << coupu_ << coupe_ << coupnu_ 
     << sL_ << cL_ << sR_ << cR_;
}

void LHTPFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> gl_ >> gr_ >> tl_ >> tr_ >> coupd_ >> coupu_ >> coupe_ >> coupnu_
     >> sL_ >> cL_ >> sR_ >> cR_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFZVertex,FFVVertex>
describeHerwigLHTPFFZVertex("Herwig::LHTPFFZVertex", "HwLHTPModel.so");

void LHTPFFZVertex::Init() {

  static ClassDocumentation<LHTPFFZVertex> documentation
    ("The LHTPFFZVertex class implements the couplings of "
     "the fermions to the Z boson and its heavy partner in the"
     " Little Higgs model with T-parity.");

}

LHTPFFZVertex::LHTPFFZVertex() 
  : gl_(37,0.0), gr_(37,0.0), tl_( 6,0.0), tr_( 6,0.0),
    coupd_(0.), coupu_(0.), coupe_(0.), coupnu_(0.),
    sL_(0.), cL_(1.), sR_(0.), cR_(1.),
    coupLast_(0.0), q2Last_(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHTPFFZVertex::doinit() {
  // Z
  // the quarks
  for(int ix = 1; ix < 7; ++ix) {
    addToList(-ix,    ix,    23);
  }
  // T+T+
  addToList(-8,  +8,  23);
  //T+t
  addToList(-6,  +8,  23);
  addToList(-8,  +6,  23);
  // the leptons
  for(int ix = 11; ix < 17; ++ix) {
    addToList(-ix,    ix,    23);
  }
  // the T-odd quarks
  for(long ix = 4000001; ix < 4000007; ++ix) {
    addToList(-ix,    ix,    23);
  }
  addToList(-4000008,  +4000008,  23);
  // the T-odd leptons
  for(long ix = 11;ix<17;++ix) {
    addToList(-ix-4000000,    ix+4000000,    23);
  }
  // Z_H
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix-4000000,    ix,    33);
    addToList(-ix,    ix+4000000,    33);
  }  
  addToList(      -8,  +4000008,  33);
  addToList(       8,  -4000008,  33);
  addToList(      -6,  +4000008,  33);
  addToList(       6,  -4000008,  33);
  // the leptons
  for(int ix=11;ix<17;++ix) {
    addToList(-ix-4000000,    ix,    33);
    addToList(-ix,    ix+4000000,    33);
  }
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHTPModel "
				   << " in LHTPFFPVertex::doinit()"
				   << Exception::runerror;
  double sw = model->sin2ThetaW();
  double cw = sqrt(1.-sw);
  sw = sqrt(sw);
  double fact = 0.25/sw/cw;
  for(int ix=1;ix<4;++ix) {
    // SM fermions
    gl_[2*ix-1]  = fact*(model->vd()  + model->ad() );
    gl_[2*ix ]   = fact*(model->vu()  + model->au() );
    gl_[2*ix+9 ] = fact*(model->ve()  + model->ae() );
    gl_[2*ix+10] = fact*(model->vnu() + model->anu());
    gr_[2*ix-1]  = fact*(model->vd()  - model->ad() );
    gr_[2*ix ]   = fact*(model->vu()  - model->au() );
    gr_[2*ix+9 ] = fact*(model->ve()  - model->ae() );
    gr_[2*ix+10] = fact*(model->vnu() - model->anu());
    // T-odd fermions
    gl_[2*ix-1 +20] = fact*(model->vd()  + model->ad() );
    gl_[2*ix   +20] = fact*(model->vu()  + model->au() );
    gl_[2*ix+9 +20] = fact*(model->ve()  + model->ae() );
    gl_[2*ix+10+20] = fact*(model->vnu() + model->anu());
    gr_[2*ix-1 +20] = gl_[2*ix-1 +20];
    gr_[2*ix   +20] = gl_[2*ix   +20];
    gr_[2*ix+9 +20] = gl_[2*ix+9 +20];
    gr_[2*ix+10+20] = gl_[2*ix+10+20];
  }
  // couplngis to Z for extended top sector
  tl_[0] = (0.5*sqr(model->cosThetaL())-2./3.*sqr(sw))/cw/sw;
  tr_[0] = -2./3.*sw/cw;
  tl_[1] = (0.5*sqr(model->sinThetaL())-2./3.*sqr(sw))/cw/sw;
  tr_[1] = -2./3.*sw/cw;
  tl_[2] =  0.5/sw/cw*model->sinThetaL()*model->cosThetaL();
  tr_[2] = 0.;
  // couplings to the Z_H of T-odd fermions
  double cH = model->cosThetaH();
  double sH = model->sinThetaH();
  sL_ = model->sinThetaL();
  cL_ = model->cosThetaL();
  sR_ = model->sinThetaR();
  cR_ = model->cosThetaR();
  coupd_  = 0.1*(sH/cw+5.*cH/sw);
  coupu_  = 0.1*(sH/cw-5.*cH/sw);
  coupe_  = 0.1*(sH/cw+5.*cH/sw);
  coupnu_ = 0.1*(sH/cw-5.*cH/sw);
  tl_[5] = -2./3.*sw/cw;
  tr_[5] = -2./3.*sw/cw;
  // couplings of T-odd top to the Z_H
  tl_[3] = 0.4*sH*sL_/cw;
  tr_[3] = 0.4*sH*sR_/cw;
  tl_[4] =-0.4*sH*cL_/cw;
  tr_[4] =-0.4*sH*cR_/cw;
  // base class initialisation
  FFVVertex::doinit();
}

void LHTPFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=q2Last_) {
    coupLast_ = -electroMagneticCoupling(q2);
    q2Last_=q2;
  }
  norm(coupLast_);
  // the left and right couplings
  long iferm = abs(a->id());
  long ianti = abs(b->id());
  long ibos  = c->id();
  if(ibos == ParticleID::Z0) {
    if(iferm == 8 || iferm == 6) {
      if(iferm == 6 && ianti == 6) {
	left (tl_[0]);
	right(tr_[0]);
      }
      else if(iferm == 8 && ianti == 8) {
	left (tl_[1]);
	right(tr_[1]);
      }
      else {
	left (tl_[2]);
	right(tr_[2]);
      }
    }
    else if(iferm == 4000008) {
      left (tl_[5]);
      right(tr_[5]);
    }
    else if((iferm >= 1 && iferm <= 6)|| (iferm >= 11 && iferm <= 16)) {
      left (gl_[iferm]);
      right(gr_[iferm]);
    }
    else {
      iferm = (iferm % 4000000) + 20;
      left (gl_[iferm]);
      right(gr_[iferm]);
    }
  }
  else if(ibos == 33) {
    if(iferm>4000000) swap(iferm,ianti);
    assert(iferm<4000000&&ianti>4000000);
    if( iferm == 6 || iferm == 8 ) {
      if     (iferm==6&&ianti==4000006) {
	left (cL_*coupu_);
	right(0.);
      }
      else if(iferm==6&&ianti==4000008) {
	left (tl_[3]);
	right(tr_[3]);
      }
      else if(iferm==8&&ianti==4000006) {
	left (sL_*coupu_);
	right(0.);
      }
      else if(iferm==8&&ianti==4000008) {
	left ( tl_[4]);
	right( tr_[4]);
      }
      else
	assert(false);
    }
    else {
      right(0.);
      if(iferm <= 6) {
	if(iferm % 2 == 0) left( coupu_ );
	else               left( coupd_ );
      }
      else {
	if(iferm % 2 == 0) left( coupnu_ );
	else               left( coupe_  );
      }
    }
  }
  else
    assert(false);
}
#line 1 "./LHTPWWWVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWWWVertex class.
//

#include "LHTPWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "LHTPModel.h"

using namespace Herwig;

LHTPWWWVertex::LHTPWWWVertex() : coupLast_(0.), q2Last_(ZERO),
				 couplings_(3 ,0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHTPWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << couplings_;
}

void LHTPWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> couplings_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPWWWVertex,VVVVertex>
describeHerwigLHTPWWWVertex("Herwig::LHTPWWWVertex", "HwLHTPModel.so");

void LHTPWWWVertex::Init() {

  static ClassDocumentation<LHTPWWWVertex> documentation
    ("The LHTPWWWVertex class implements the coupling of three "
     "electroweak gauge bosons and their heavy partners in the "
     "Little Higgs model with T-parity.");

}

void LHTPWWWVertex::doinit() {
  //SM interactions
  addToList( 24, -24,  22);
  addToList( 24, -24,  23);
  //LHTP
  //W_H W_H A_L
  addToList( 34, -34,  22);
  //W_H W_H Z_L
  addToList( 34, -34,  23);
  //W_H W_L A_H
  addToList( 34, -24,  32);
  addToList( 24, -34,  32); 
  //W_H W_L Z_H
  addToList( 34, -24,  33);
  addToList( 24, -34,  33);
  VVVVertex::doinit();
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if( !model )
    throw InitException() 
      << "LHTPWWWVertex::doinit() - Model pointer must be of LHTPModel"
      << "type, cannot continue without this."
      << Exception::abortnow;
  double sw(sqrt(model->sin2ThetaW()));
  double cw(sqrt(1. - model->sin2ThetaW()));
  //W W Z
  couplings_[0] = cw/sw;
  //W_L W_H A_H
  couplings_[1] = model->sinThetaH()/sw;
  //W_L W_H Z_H
  couplings_[2] = 1./sw;

}

void LHTPWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  if(q2 != q2Last_) {
    coupLast_ = electroMagneticCoupling(q2);
    q2Last_ = q2;
  }
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // get the PDG code for the neutral boson
  long boson(0);
  if(!a->charged())      {
    boson = ida;
    ida = 22;
  }
  else if(abs(ida) !=ParticleID::Wplus) {
    ida = ida > 0 ? 24 : -24;
  }
  if(!b->charged()) {
    boson = idb;
    idb = 22;
  }
  else if(abs(idb) !=ParticleID::Wplus) {
    idb = idb > 0 ? 24 : -24;
  }
  if(!c->charged()) {
    boson = idc;
    idc = 22;
  }
  else if(abs(idc) !=ParticleID::Wplus) {
    idc = idc > 0 ? 24 : -24;
  }
  assert( boson ==22 || boson==23 || boson==32 || boson==33);
  // get the prefactor
  double pre(0.);
  switch (boson) {
  case 22:
    pre = 1.;
    break;
  case 23:
    pre = couplings_[0];
    break;
  case 32:
    pre = couplings_[1];
    break;
  case 33:
    pre = couplings_[2];
    break;
  default:
    assert(false);
  };
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || 
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )      norm( coupLast_*pre);
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) ) norm(-coupLast_*pre);
  else 
    throw Helicity::HelicityConsistencyError() 
      << "LHTPWWWVertex::setCoupling - Incorrect particles in LHTPWWWVertex. " 
      << a->id() << " " << b->id() << " " << c->id() << '\n'
      << Exception::runerror;
}
#line 1 "./LHTPHHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPHHHVertex class.
//

#include "LHTPHHHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHTPHHHVertex::LHTPHHHVertex() : ratio_(ZERO), coupLast_(0.), q2Last_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHTPHHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPHHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(ratio_,GeV);
}

void LHTPHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ratio_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<LHTPHHHVertex,SSSVertex>
describeHerwigLHTPHHHVertex("Herwig::LHTPHHHVertex", "HwLHTPModel.so");

void LHTPHHHVertex::Init() {

  static ClassDocumentation<LHTPHHHVertex> documentation
    ("The LHTPHHHVertex class implements the trilinear Higgs boson"
     " self couplings in the Little Higgs model with T-parity");

}

void LHTPHHHVertex::doinit() {
  addToList(25,25, 25);
  addToList(25,35, 35);
  addToList(25,36, 36);
  addToList(25,37,-37);
  SSSVertex::doinit();
  ratio_ = sqr(getParticleData(ParticleID::h0)->mass())/
    getParticleData(ParticleID::Wplus)->mass();
}

void LHTPHHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr) {
  if(q2!=q2Last_||coupLast_==0.) {
    coupLast_ = weakCoupling(q2)*ratio_*UnitRemoval::InvE;
    q2Last_=q2;
  }
  long id = part2->id()!=ParticleID::h0 ? abs(part2->id()) : abs(part1->id());
  if(id==ParticleID::h0) norm(    -coupLast_);
  else if(id==35)        norm(  3.*coupLast_);
  else if(id==36)        norm(  3.*coupLast_);
  else if(id==37)        norm(1.5*coupLast_);
  else assert(false);
}
#line 1 "./LHTPWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPWHHVertex class.
//

#include "LHTPWHHVertex.h"
#include "LHTPModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHTPWHHVertex::LHTPWHHVertex() : 
  coupLast_(0.), q2Last_(ZERO), coup_(11) {
  orderInGs(0);
  orderInGem(1);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHTPWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHTPWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHTPWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void LHTPWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPWHHVertex,VSSVertex>
describeHerwigLHTPWHHVertex("Herwig::LHTPWHHVertex", "HwLHTPModel.so");

void LHTPWHHVertex::Init() {

  static ClassDocumentation<LHTPWHHVertex> documentation
    ("The LHTPWHHVertex class implements the coupling of a pair of Higgs"
     " bosons to an electroweak gauge boson in the Little"
     " Higgs model with T-parity.");

}

void LHTPWHHVertex::doinit() {
  // photon
  addToList( 22, 37,-37);
  addToList( 22, 38,-38);
  // Z0
  addToList( 23, 37,-37);
  addToList( 23, 38,-38);
  addToList( 23, 35, 36);
  // W+
  addToList( 24, 35,-37);
  addToList( 24, 36,-37);
  addToList( 24, 37,-38);
  // W-
  addToList(-24, 35, 37);
  addToList(-24, 36, 37);
  addToList(-24,-37, 38);
  // A_H
  addToList( 32, 25, 36);
  // Z_H
  addToList( 33, 25, 36);
  // W_H
  addToList( 34, 25,-37);
  addToList(-34, 25, 37);
  VSSVertex::doinit();
  // model
  cLHTPModelPtr model = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double vf(model->vev()/model->f());
  coup_[ 0] = 1.;
  coup_[ 1] = 2.;
  coup_[ 2] =-sw/cw;
  coup_[ 3] = (1.-2.*sw2)/cw/sw;
  coup_[ 4] =-Complex(0.,1.)/cw/sw;
  coup_[ 5] = sqrt(0.5)/sw;
  coup_[ 6] = Complex(0.,1.)/sw*sqrt(0.5);
  coup_[ 7] = 1./sw;
  coup_[ 8] = Complex(0.,1)*sqrt(0.5)*vf/3./cw;
  coup_[ 9] = Complex(0.,1)*sqrt(0.5)*vf/3./sw;
  coup_[10] =-vf/6./sw;
}

void LHTPWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
				tcPDPtr particle2, tcPDPtr particle3) {
  if( q2 != q2Last_ || coupLast_==0.) {
    q2Last_ = q2;
    coupLast_ = electroMagneticCoupling(q2);
  }
  int ibos = particle1->id();
  int isc1 = particle2->id();
  int isc2 = particle3->id();
  if(ibos==ParticleID::gamma) {
    if(isc1==37) 
      norm(coup_[0]*coupLast_);
    else if(isc1==38)
      norm(coup_[1]*coupLast_);
    else if(isc1==-37) 
      norm(-coup_[0]*coupLast_);
    else if(isc1==-38)
      norm(-coup_[1]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Z0) {
    if(isc1==37) 
      norm(coup_[2]*coupLast_);
    else if(isc1==38)
      norm(coup_[3]*coupLast_);
    else if(isc1==-37) 
      norm(-coup_[2]*coupLast_);
    else if(isc1==-38)
      norm(-coup_[3]*coupLast_);
    else if(isc1==35)
      norm(coup_[4]*coupLast_);
    else if(isc2==35)
      norm(-coup_[4]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wplus) {
    if(isc1==35)
      norm(coup_[5]*coupLast_);
    else if(isc1==36)
      norm(coup_[6]*coupLast_);
    else if(isc1==-38)
      norm(-coup_[7]*coupLast_);
    else if(isc2==35)
      norm(-coup_[5]*coupLast_);
    else if(isc2==36)
      norm(-coup_[6]*coupLast_);
    else if(isc2==-38)
      norm( coup_[7]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wminus) {
    if(isc1==35)
      norm(conj(coup_[5])*coupLast_);
    else if(isc1==36)
      norm(conj(coup_[6])*coupLast_);
    else if(isc1==38)
      norm(-conj(coup_[7])*coupLast_);
    else if(isc2==35)
      norm(-conj(coup_[5])*coupLast_);
    else if(isc2==36)
      norm(-conj(coup_[6])*coupLast_);
    else if(isc2==38)
      norm( conj(coup_[7])*coupLast_);
    else
      assert(false);
  }
  else if(ibos==32) {
    if(isc1==25)
      norm( coup_[8]*coupLast_);
    else if(isc2==25)
      norm(-coup_[8]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==33) {
    if(isc1==25)
      norm( coup_[9]*coupLast_);
    else if(isc2==25)
      norm(-coup_[9]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==34) {
    if(isc1==25)
      norm( coup_[10]*coupLast_);
    else if(isc2==25)
      norm(-coup_[10]*coupLast_);
    else
      assert(false);
  }
  else if(ibos==-34) {
    if(isc1==25)
      norm( conj(coup_[10])*coupLast_);
    else if(isc2==25)
      norm(-conj(coup_[10])*coupLast_);
    else
      assert(false);
  }
  else
    assert(false);
}
#line 1 "./LHTPFFHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHTPFFHVertex class.
//

#include "LHTPFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHTPFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(cL_,1./GeV) << ounit(cR_,1./GeV) << model_;
}

void LHTPFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(cL_,1./GeV) >> iunit(cR_,1./GeV) >> model_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHTPFFHVertex,FFSVertex>
describeHerwigLHTPFFHVertex("Herwig::LHTPFFHVertex", "HwLHTPModel.so");

void LHTPFFHVertex::Init() {

  static ClassDocumentation<LHTPFFHVertex> documentation
    ("The LHTPFFHVertex class implements the interaction of the fermions"
     " and the Higgs bosons in the Little Higgs model with T-parity");

}

LHTPFFHVertex::LHTPFFHVertex() 
  : q2Last_(ZERO) {
  orderInGem(1);
  orderInGs(0);
  massLast_[0] = 0.*GeV; 
  massLast_[1] = 0.*GeV;
  idLast_[0] = 0;
  idLast_[1] = 0;
  colourStructure(ColourStructure::DELTA);
}

void LHTPFFHVertex::doinit() {
  // SM like higgs
  addToList(  -3,   3,  25);
  addToList(  -4,   4,  25);
  addToList(  -5,   5,  25);
  addToList(  -6,   6,  25);
  addToList(  -6,   8,  25);
  addToList(  -8,   6,  25);
  addToList(  -8,   8,  25);
  addToList( -13,  13,  25);
  addToList( -15,  15,  25);
  addToList( -4000002, 4000002,  25);
  addToList( -4000004, 4000004,  25);
  addToList( -4000006, 4000006,  25);
  addToList( -4000012, 4000012,  25);
  addToList( -4000014, 4000014,  25);
  addToList( -4000016, 4000016,  25);
  // phi0
  addToList( -3      , 4000003,  35);
  addToList( -4      , 4000004,  35);
  addToList( -5      , 4000005,  35);
  addToList( -4000003,       3,  35);
  addToList( -4000004,       4,  35);
  addToList( -4000005,       5,  35);
  addToList( -6      , 4000006,  35);
  addToList( -8      , 4000006,  35);
  addToList( -4000006,       6,  35);
  addToList( -4000006,       8,  35);
  // phiP
  addToList( -2      , 4000002,  36);
  addToList( -3      , 4000003,  36);
  addToList( -4      , 4000004,  36);
  addToList( -5      , 4000005,  36);
  addToList( -4000002,       2,  36);
  addToList( -4000003,       3,  36);
  addToList( -4000004,       4,  36);
  addToList( -4000005,       5,  36);
  addToList( -12     , 4000012,  36);
  addToList( -14     , 4000014,  36);
  addToList( -16     , 4000016,  36);
  addToList( -4000012,      12,  36);
  addToList( -4000014,      14,  36);
  addToList( -4000016,      16,  36);
  addToList( -6      , 4000006,  36);
  addToList( -6      , 4000008,  36);
  addToList( -8      , 4000006,  36);
  addToList( -8      , 4000008,  36);
  addToList( -4000006,       6,  36);
  addToList( -4000008,       6,  36);
  addToList( -4000006,       8,  36);
  addToList( -4000008,       8,  36);
  // phi +/-
  addToList( -1      , 4000002, -37);
  addToList( -3      , 4000004, -37);
  addToList( -5      , 4000006, -37);
  addToList( -4000001,       2, -37);
  addToList( -4000003,       4, -37);
  addToList( -4000005,       6, -37);
  addToList( -4000005,       8, -37);
  addToList( -4000002,       1,  37);
  addToList( -4000004,       3,  37);
  addToList( -4000006,       5,  37);
  addToList( -2      , 4000001,  37);
  addToList( -4      , 4000003,  37);
  addToList( -6      , 4000005,  37);
  addToList( -8      , 4000005,  37);
  addToList( -11     , 4000012, -37);
  addToList( -13     , 4000014, -37);
  addToList( -15     , 4000016, -37);
  addToList( -4000011,      12, -37);
  addToList( -4000013,      14, -37);
  addToList( -4000015,      16, -37);
  addToList( -4000012, 11     ,  37);
  addToList( -4000014, 13     ,  37);
  addToList( -4000016, 15     ,  37);
  addToList(      -12, 4000011,  37);
  addToList(      -14, 4000013,  37);
  addToList(      -16, 4000015,  37);
  model_ = 
    dynamic_ptr_cast<cLHTPModelPtr>(generator()->standardModel());
  if(!model_)   throw InitException() << "Must be using the LHModel "
				      << " in LHFFPVertex::doinit()"
				      << Exception::runerror;
  cL_  .resize(18); 
  cR_  .resize(18);
  Energy v  = model_->vev();
  Energy f  = model_->f();
  double vf = model_->vev()/model_->f();
  double sa = model_->sinAlpha();
  double ca = model_->cosAlpha();
  // lightest higgs couplings
  // coupling of light SM fermions
  cL_[0]   =  cR_[0] = 1./v;
  // couplings to top quarks
  cL_[1]   =  cR_[1] = sa*ca/f;
  cL_[2]   = -sa/ca/v;
  cR_[2]   =  sqr(ca)*vf/v;
  // couplings to T-odd quarks
  cL_[3]   =  cR_[3] = 0.5*sqrt(0.5)/f*model_->kappaQuark();
  // couplings to T-odd leptons
  cL_[4]   =  cR_[4] = 0.5*sqrt(0.5)/f*model_->kappaLepton();
  // Phi0
  // quark, T-odd quark
  cL_[5]   = sqrt(0.5)/f;
  cR_[5]   = ZERO;
  // and top quarks
  cL_[6]   = sqrt(0.5)*model_->cosThetaR()/f/ca;
  cR_[6]   = ZERO;
  cL_[7]   = sqrt(0.5)*model_->sinThetaR()/f/ca;
  cR_[7]   = ZERO;
  // PhiP
  // quark, T-odd quark
  cL_[8] = vf/f*model_->kappaQuark()/12;
  cR_[8] = sqrt(0.5)/f;
  // lepton, T-odd lepton
  cL_[9] = vf/f*model_->kappaLepton()/12;
  cR_[9] = ZERO;
  // top, T_-
  cL_[10] = model_->cosThetaR()*sqrt(2.)*vf/f/ca/3;
  cR_[10] = ZERO;
  cL_[11] = model_->sinThetaR()*sqrt(2.)*vf/f/ca/3;
  cR_[11] = ZERO;
  // top, t_-
  cL_[12] = model_->cosThetaL()*vf/f/12.*model_->kappaQuark();
  cR_[12] = model_->cosThetaR()*sqrt(0.5)/f/ca;
  cL_[13] = model_->sinThetaL()*vf/f/12.*model_->kappaQuark();
  cR_[13] = model_->sinThetaR()*sqrt(0.5)/f/ca;
  // Phi +/-
  cL_[14] = vf/f*model_->kappaLepton()/24.;
  cR_[14] = ZERO;
  // quark T-odd quark
  cL_[15] = vf/f*model_->kappaQuark() /24.;
  cR_[15] =-vf*sqrt(0.5)/v;
  cL_[16] = model_->cosThetaL()*vf/f*model_->kappaQuark() /24.;
  cR_[16] =-model_->cosThetaR()*vf*sqrt(0.5)/v/ca;
  cL_[17] = model_->sinThetaL()*vf/f*model_->kappaQuark() /24.;
  cR_[17] =-model_->sinThetaR()*vf*sqrt(0.5)/v/ca;
  FFSVertex::doinit();
}

void LHTPFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  norm(1.);
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  //  int ihigg=abs(c->id());
  // left and right couplings set to one
  // SM like higgs
  if(c->id()==ParticleID::h0) {
    // to SM fermions and T
    if(iferm<=16&&ianti<=16) {
      // running masses
      if(q2!=q2Last_||idLast_[0]!=iferm||idLast_[1]!=ianti) {
	q2Last_ = q2;
	idLast_[0] = iferm;
	assert((idLast_[0]>=1  && idLast_[0]<=8 ) || 
	       (idLast_[0]>=11 && idLast_[0]<=16));
	if(idLast_[0]!=8)
	  massLast_[0] = model_->mass(q2,a);
	else
	  massLast_[0] = model_->mass(q2,getParticleData(ParticleID::t));
	idLast_[1] = ianti;
	assert((idLast_[1]>=1  && idLast_[1]<=8 ) || 
	       (idLast_[1]>=11 && idLast_[1]<=16));
	if(idLast_[0]!=idLast_[1]) {
	  if(idLast_[1]!=8)
	    massLast_[1] = model_->mass(q2,a);
	  else
	    massLast_[1] = model_->mass(q2,getParticleData(ParticleID::t));
	}
	else {
	  massLast_[1] = massLast_[0];
	}
      }
      if(iferm<6||iferm>8) {
	left (-Complex(cL_[0]*massLast_[0]));
	right(-Complex(cR_[0]*massLast_[0]));
      }
      else {
	if(iferm==8&&ianti==8) {
	  left ( Complex(cL_[1]*massLast_[0]));
	  right( Complex(cR_[1]*massLast_[0]));
	}
	else {
	  if(a->id()==ParticleID::tbar||b->id()==ParticleID::tbar) {
	    left (-Complex(cL_[2]*massLast_[0]));
	    right(-Complex(cR_[2]*massLast_[0]));
	  }
	  else {
	    left (-Complex(cR_[2]*massLast_[0]));
	    right(-Complex(cL_[2]*massLast_[0]));
	  } 
	}
      }
    }
    else {
      if(iferm<=4000006) {
	left ( Complex(cL_[3]*model_->vev()));
	right( Complex(cR_[3]*model_->vev()));
      }
      else {
	left ( Complex(cL_[4]*model_->vev()));
	right( Complex(cR_[4]*model_->vev()));
      }
    }
  }
  // Phi0
  else if(c->id()==ParticleID::H0 ||
	  c->id()==ParticleID::A0) {
    tcPDPtr ferm = a;
    if(iferm>4000000) {
      swap(iferm,ianti);
      ferm = b;
    }
    if(q2!=q2Last_||idLast_[0]!=iferm) {
      q2Last_ = q2;
      idLast_[0] = iferm;
      assert((idLast_[0]>=1  && idLast_[0]<=8 ) || 
	     (idLast_[0]>=11 && idLast_[0]<=16));
      if(idLast_[0]!=8)
	massLast_[0] = model_->mass(q2,ferm);
      else
	massLast_[0] = model_->mass(q2,getParticleData(ParticleID::t));
    }
    if(c->id()==ParticleID::H0 ) {
      unsigned int      iloc = 5;
      if(iferm==6)      iloc = 6;
      else if(iferm==8) iloc = 7;
      if( (a->id()>=1&&a->id()<=8) || (b->id()>=1&&b->id()<=8) ) {
	left ( Complex(cR_[iloc]*massLast_[0]));
	right( Complex(cL_[iloc]*massLast_[0]));
      }
      else {
	left ( Complex(cL_[iloc]*massLast_[0]));
	right( Complex(cR_[iloc]*massLast_[0]));
      }
    }
    // PhiP
    else if(c->id()==ParticleID::A0) {
      if(iferm<=5) {
	if( (a->id()>=1&&a->id()<=5) || (b->id()>=1&&b->id()<=5) ) {
	  if(iferm%2==0) {
	    right(Complex(0., 1.)*model_->vev()*cL_[8]);
	    left (Complex(0.,-1.)*massLast_[0] *cR_[8]);
	  }
	  else {
	    right(Complex(ZERO));
	    left (Complex(0., 1.)*massLast_[0] *cR_[8]);
	  }
	}
	else {
	  if(iferm%2==0) {
	    left (Complex(0.,-1.)*model_->vev()*cL_[8]);
	    right(Complex(0., 1.)*massLast_[0] *cR_[8]);
	  }
	  else {
	    left (Complex(ZERO));
	    right(Complex(0.,-1.)*massLast_[0] *cR_[8]);
	  }
	}
      }
      else if(iferm>=12) {
	if( (a->id()>=11&&a->id()<=16) || (b->id()>=11&&b->id()<=16) ) {
	  right(Complex(0., 1.)*model_->vev()*cL_[9]);
	  left (Complex(0.,-1.)*massLast_[0] *cR_[9]);
	}
	else {
	  right(Complex(0.,-1.)*massLast_[0] *cR_[9]);
	  left (Complex(0., 1.)*model_->vev()*cL_[9]);
	}
      }
      else {
	if(ianti==4000008) {
	  unsigned int iloc = (iferm+14)/2;
	  if( (a->id()==6||a->id()==8) || (b->id()==6||b->id()==8) ) {
	    left (Complex(0., 1.)*massLast_[0]*cR_[iloc]);
	    right(Complex(0., 1.)*massLast_[0]*cL_[iloc]);
	  }
	  else {
	    left (Complex(0.,-1.)*massLast_[0]*cL_[iloc]);
	    right(Complex(0.,-1.)*massLast_[0]*cR_[iloc]);
	  }
	}
	else {
	  unsigned int iloc = (iferm+18)/2;
	  if( (a->id()==6||a->id()==8) || (b->id()==6||b->id()==8) ) {
	    left (Complex(0., 1.)*model_->vev()*cL_[iloc]);
	    right(Complex(0.,-1.)*massLast_[0] *cR_[iloc]);
	  }
	  else {
	    left (Complex(0., 1.)*massLast_[0] *cR_[iloc]);
	    right(Complex(0.,-1.)*model_->vev()*cL_[iloc]);
	  }
	}
      }
    }
  }
  else if(abs(c->id())==ParticleID::Hplus) {
    tcPDPtr ferm = a;
    if(iferm>4000000) {
      swap(iferm,ianti);
      ferm = b;
    }
    if(q2!=q2Last_||idLast_[0]!=iferm) {
      q2Last_ = q2;
      idLast_[0] = iferm;
      assert((idLast_[0]>=1  && idLast_[0]<=8 ) || 
	     (idLast_[0]>=11 && idLast_[0]<=16));
      if(idLast_[0]!=8)
	massLast_[0] = model_->mass(q2,ferm);
      else
	massLast_[0] = model_->mass(q2,getParticleData(ParticleID::t));
    }
    Complex cleft(0.),cright(0.);
    // lepton and T-odd lepton
    if(iferm>=11&&iferm<=16) {
      cright = Complex(cR_[14]*massLast_[0]);
      cleft  = Complex(cL_[14]*model_->vev()); 
    }
    else if(iferm>=1&&iferm<=6) {
      cright = Complex(cR_[15]*massLast_[0]);
      cleft  = Complex(cL_[15]*model_->vev());
    }
    else if(iferm==6) {
      cright = Complex(cR_[16]*massLast_[0]);
      cleft  = Complex(cL_[16]*model_->vev()); 
    }
    else if(iferm==8) {
      cright = Complex(cR_[17]*massLast_[0]);
      cleft  = Complex(cL_[17]*model_->vev()); 
    }
    if((a->id()>=1&&a->id()<=16) ||(b->id()>=1&&b->id()<=16) ) {
      swap(cleft,cright);
      cleft  *= -1.;
      cright *= -1.;
    }
    if(c->id()==ParticleID::Hminus) {
      cleft  *= -1.;
      cright *= -1.;
    }
    left (cleft );
    right(cright);
  }
}
