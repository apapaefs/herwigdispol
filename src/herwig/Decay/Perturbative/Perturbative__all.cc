#line 1 "./SMWDecayer.cc"
// -*- C++ -*-
//
// SMWDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWDecayer class.
//

#include "SMWDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMWDecayer::EPS_=0.00000001;

SMWDecayer::SMWDecayer()
  : quarkWeight_(6,0.), leptonWeight_(3,0.), CF_(4./3.),
    NLO_(false) {
  quarkWeight_[0]  = 1.01596;
  quarkWeight_[1]  = 0.0537308;
  quarkWeight_[2]  = 0.0538085;
  quarkWeight_[3]  = 1.01377;
  quarkWeight_[4]  = 1.45763e-05;
  quarkWeight_[5]  = 0.0018143;
  leptonWeight_[0] = 0.356594;
  leptonWeight_[1] = 0.356593;
  leptonWeight_[2] = 0.356333;
  // intermediates
  generateIntermediates(false);
}

void SMWDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig StandardModel object in"
				  << "SMWDecayer::doinit()"
				  << Exception::runerror;
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  WWWVertex_ = hwsm->vertexWWW();
  FFPVertex_ = hwsm->vertexFFP();
  // make sure they are initialized
  FFGVertex_->init();
  FFWVertex_->init();
  WWWVertex_->init();
  FFPVertex_->init();
  // now set up the decay modes
  // W modes
  tPDPtr Wp = getParticleData(ParticleID::Wplus);
  // loop for the quarks
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(!FFWVertex_->allowed(-ix,iy,ParticleID::Wminus))
	throw InitException() << "SMWDecayer::doinit() the W vertex" 
			      << "cannot handle all the quark modes" 
			      << Exception::abortnow;
      tPDVector out = {getParticleData(-ix),
		       getParticleData( iy)};
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(Wp,out,quarkWeight_[iz]));
      addMode(mode);
      ++iz;
    }
  }
  // loop for the leptons
  for(int ix=11;ix<17;ix+=2) {
    tPDVector out = {getParticleData(-ix  ),
		     getParticleData( ix+1)};
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(Wp,out,leptonWeight_[(ix-11)/2]));
    addMode(mode);
  }
  gluon_ = getParticleData(ParticleID::g);
}

int SMWDecayer::modeNumber(bool & cc,tcPDPtr parent, 
			    const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  int id0=parent->id();
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if(abs(id0)!=ParticleID::Wplus) return imode;
  int idd(0),idu(0);
  if(abs(id1)%2==1&&abs(id2)%2==0) {
    idd=abs(id1);
    idu=abs(id2);
  }
  else if(abs(id1)%2==0&&abs(id2)%2==1) {
    idd=abs(id2);
    idu=abs(id1);
  }
  if(idd==0&&idu==0) {
    return imode;
  }
  else if(idd<=5) {
    imode=idd+idu/2-2;
  }
  else {
    imode=(idd-1)/2+1;
  }
  cc= (id0==ParticleID::Wminus);
  return imode;
}

ParticleVector SMWDecayer::decay(const Particle & parent,
				 const tPDVector & children) const {
  // generate the decay
  bool cc;
  unsigned int imode = modeNumber(cc,parent.dataPtr(),children);
  ParticleVector output = generate(false,cc,imode,parent);
  if(output[0]->hasColour())      output[0]->antiColourNeighbour(output[1]);
  else if(output[1]->hasColour()) output[1]->antiColourNeighbour(output[0]);
  return output;
}

void SMWDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFWVertex_ << quarkWeight_ << leptonWeight_
     << FFGVertex_ << gluon_ << NLO_
     << WWWVertex_ << FFPVertex_;  
}

void SMWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> quarkWeight_ >> leptonWeight_
     >> FFGVertex_ >> gluon_ >> NLO_
     >> WWWVertex_ >> FFPVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWDecayer,PerturbativeDecayer>
describeHerwigSMWDecayer("Herwig::SMWDecayer", "HwPerturbativeDecay.so");

void SMWDecayer::Init() {

  static ClassDocumentation<SMWDecayer> documentation
    ("The SMWDecayer class is the implementation of the decay"
     " of the W boson to the Standard Model fermions.");

  static ParVector<SMWDecayer,double> interfaceWquarkMax
    ("QuarkMax",
     "The maximum weight for the decay of the W to quarks",
     &SMWDecayer::quarkWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMWDecayer,double> interfaceWleptonMax
    ("LeptonMax",
     "The maximum weight for the decay of the W to leptons",
     &SMWDecayer::leptonWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static Switch<SMWDecayer,bool> interfaceNLO
    ("NLO",
     "Whether to return the LO or NLO result",
     &SMWDecayer::NLO_, false, false, false);
  static SwitchOption interfaceNLOLO
    (interfaceNLO,
     "No",
     "Leading-order result",
     false);
  static SwitchOption interfaceNLONLO
    (interfaceNLO,
     "Yes",
     "NLO result",
     true);

}
void SMWDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

// return the matrix element squared
double SMWDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME()) 
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(outgoing[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
  					       const_ptr_cast<tPPtr>(&part),
  					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  SpinorBarWaveFunction wbar(momenta[iferm],outgoing[iferm],Helicity::outgoing);
  SpinorWaveFunction    w   (momenta[ianti],outgoing[ianti],Helicity::outgoing);
  wavebar_.resize(2);
  wave_   .resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wbar.reset(ihel);
    wavebar_[ihel] = wbar;
    w.reset(ihel);
    wave_   [ihel] = w;
  }
  // compute the matrix element
  Energy2 scale(sqr(part.mass()));
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
  	if(iferm>ianti) (*ME())(vhel,ia,ifm)=
  	  FFWVertex_->evaluate(scale,wave_[ia],wavebar_[ifm],vectors_[vhel]);
  	else            (*ME())(vhel,ifm,ia)=
  	  FFWVertex_->evaluate(scale,wave_[ia],wavebar_[ifm],vectors_[vhel]);
      }
    }
  }
  double output=(ME()->contract(rho_)).real()*UnitRemoval::E2/scale;
  if(abs(outgoing[0]->id())<=6) output*=3.;
  // // leading-order result
  if(!NLO_) return output;
  // check decay products coloured, otherwise return
  if(!outgoing[0]->coloured()) return output;
  // inital masses, couplings  etc
  // W mass
  mW_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mW_));
  // reduced mass
  double mu1  = outgoing[0]->mass()/mW_;
  double mu2  = outgoing[1]->mass()/mW_;
  // scale
  scale_ = sqr(mW_);
  // now for the nlo loop correction
  double virt = CF_*aS_/Constants::pi;
  // now for the real correction
  double realFact=0.;
  for(int iemit=0;iemit<2;++iemit) {
    double phi  = UseRandom::rnd()*Constants::twopi;
    // set the emitter and the spectator
    double muj  = iemit==0 ? mu1 : mu2;
    double muk  = iemit==0 ? mu2 : mu1;
    double muj2 = sqr(muj);
    double muk2 = sqr(muk);
    // calculate y
    double yminus = 0.; 
    double yplus  = 1.-2.*muk*(1.-muk)/(1.-muj2-muk2);
    double y = yminus + UseRandom::rnd()*(yplus-yminus);
    double v = sqrt(sqr(2.*muk2 + (1.-muj2-muk2)*(1.-y))-4.*muk2)
      /(1.-muj2-muk2)/(1.-y);
    double zplus  = (1.+v)*(1.-muj2-muk2)*y/2./(muj2+(1.-muj2-muk2)*y);
    double zminus = (1.-v)*(1.-muj2-muk2)*y/2./(muj2+(1.-muj2-muk2)*y);
    double z = zminus + UseRandom::rnd()*(zplus-zminus);
    double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
    // calculate x1,x2,x3,xT
    double x2 = 1.-y*(1.-muj2-muk2)-muj2+muk2;
    double x1 = 1.+muj2-muk2-z*(x2-2.*muk2);
    // copy the particle objects over for calculateRealEmission
    tcPDVector oTemp = {part.dataPtr(),outgoing[0],outgoing[1],gluon_};
    vector<Lorentz5Momentum> mom = {part.momentum(),momenta[0],momenta[1]};
    realFact += 0.25*jac*sqr(1.-muj2-muk2)/
      sqrt((1.-sqr(muj-muk))*(1.-sqr(muj+muk)))/Constants::twopi
      *2.*CF_*aS_*calculateRealEmission(x1, x2, oTemp, mom, phi, 
  					muj, muk, iemit, true);
  }
  // the born + virtual + real
  output *= (1. + virt + realFact);
  return output;
}

void SMWDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<6) quarkWeight_ [ix]=mode(ix)->maxWeight();
      else     leptonWeight_[ix-6]=mode(ix)->maxWeight();
    }
  }
}

void SMWDecayer::dataBaseOutput(ofstream & output,
				 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<quarkWeight_.size();++ix) {
    output << "newdef " << name() << ":QuarkMax " << ix << " "
	   << quarkWeight_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<leptonWeight_.size();++ix) {
    output << "newdef " << name() << ":LeptonMax " << ix << " "
	   << leptonWeight_[ix] << "\n";
  }
  // parameters for the PerturbativeDecayer base class
  PerturbativeDecayer::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}


void SMWDecayer::
initializeMECorrection(RealEmissionProcessPtr born, double & initial,
		       double & final) {
  // get the quark and antiquark
  ParticleVector qq; 
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix)
    qq.push_back(born->bornOutgoing()[ix]);
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // centre of mass energy
  d_Q_ = (qq[0]->momentum() + qq[1]->momentum()).m();
  // quark mass
  d_m_ = 0.5*(qq[0]->momentum().m()+qq[1]->momentum().m());
  // set the other parameters
  setRho(sqr(d_m_/d_Q_));
  setKtildeSymm();
  // otherwise can do it
  initial=1.;
  final  =1.;
}

bool SMWDecayer::softMatrixElementVeto(PPtr parent,
				       PPtr progenitor,
				       const bool & ,
				       const Energy & highestpT,
				       const vector<tcPDPtr> & ids,
				       const double & d_z,
				       const Energy & d_qt,
				       const Energy & ) {
  // check we should be applying the veto
  if(parent->id()!=progenitor->id()||
     ids[0]!=ids[1]||
     ids[2]->id()!=ParticleID::g) return false;
  // calculate pt
  Energy2 d_m2 = parent->momentum().m2();
  Energy2 pPerp2 = sqr(d_z*d_qt) - d_m2;
  if(pPerp2<ZERO) return true;
  Energy pPerp = (1.-d_z)*sqrt(pPerp2);
  // if not hardest so far don't apply veto
  if(pPerp<highestpT) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight and return
  return !UseRandom::rndbool(weight);
}


void SMWDecayer::setRho(double r) 
{ 
  d_rho_ = r;
  d_v_ = sqrt(1.-4.*d_rho_);
}

void SMWDecayer::setKtildeSymm() { 
  d_kt1_ = (1. + sqrt(1. - 4.*d_rho_))/2.;
  setKtilde2();
}

void SMWDecayer::setKtilde2() { 
   double num = d_rho_ * d_kt1_ + 0.25 * d_v_ *(1.+d_v_)*(1.+d_v_);
   double den = d_kt1_ - d_rho_;
   d_kt2_ = num/den;
}

double SMWDecayer::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho_);
  return uval + num/den;
}

double SMWDecayer::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double SMWDecayer::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    - 8.*d_rho_*(1.+2.*d_rho_);
  double den = (1.+2.*d_rho_)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMWDecayer::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    + 2.*d_rho_*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho_);
  double den = d_v_*d_v_*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMWDecayer::u(double x2) {
  return 0.5*(1. + d_rho_/(1.-x2+d_rho_));
}

void SMWDecayer::
getXXbar(double kti, double z, double &x, double &xbar) {
  double w = sqr(d_v_) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z);
  if (w < 0) {
    x = -1.; 
    xbar = -1;
  } else {
    x = (1. + sqr(d_v_)*(-1. + z) + sqr(kti*(-1. + z))*z*z*z 
	 + z*sqrt(w)
	 - kti*(-1. + z)*z*(2. + z*(-2 + sqrt(w))))/
      (1. - kti*(-1. + z)*z + sqrt(w));
    xbar = 1. + kti*(-1. + z)*z;
  }
}

double SMWDecayer::qWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the quark emission zone?
  if(k1 < d_kt1_) {
    rval = MEV(x, xbar)/PS(x, xbar);
    // is it also in the anti-quark emission zone?
    if(k2 < d_kt2_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMWDecayer::qbarWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the antiquark emission zone?
  if(k2 < d_kt2_) {
    rval = MEV(x, xbar)/PS(xbar, x);
    // is it also in the quark emission zone?
    if(k1 < d_kt1_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMWDecayer::qWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, x, xb);
  // if exceptionally out of phase space, leave this emission, as there 
  // is no good interpretation for the soft ME correction. 
  if (x < 0 || xb < 0) return 1.0; 
  return qWeight(x, xb); 
}

double SMWDecayer::qbarWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, xb, x);
  // see above in qWeightX. 
  if (x < 0 || xb < 0) return 1.0; 
  return qbarWeight(x, xb); 
}

double SMWDecayer::PS(double x, double xbar) {
  double u = 0.5*(1. + d_rho_ / (1.-xbar+d_rho_));
  double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho_);
  double brack = (1.+z*z)/(1.-z)- 2.*d_rho_/(1-xbar);
  // interesting: the splitting function without the subtraction
  // term. Actually gives a much worse approximation in the collinear
  // limit.  double brack = (1.+z*z)/(1.-z);
  double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho_);
  return brack/den;
}

double SMWDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
				      const ParticleVector & decay3, MEOption,
				      ShowerInteraction inter) {
  // extract partons and LO momentas
  vector<tcPDPtr> partons(1,inpart.dataPtr());
  vector<Lorentz5Momentum> lomom(1,inpart.momentum());
  for(unsigned int ix=0;ix<2;++ix) {
    partons.push_back(decay2[ix]->dataPtr());
    lomom.push_back(decay2[ix]->momentum());
  }
  vector<Lorentz5Momentum> realmom(1,inpart.momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==2) partons.push_back(decay3[ix]->dataPtr());
    realmom.push_back(decay3[ix]->momentum());
  }
  if(partons[0]->id()<0) {
    swap(partons[1],partons[2]);
    swap(lomom[1],lomom[2]);
    swap(realmom[1],realmom[2]);
  }
  scale_ = sqr(inpart.mass());
  double     lome = loME(partons,lomom);
  InvEnergy2 reme = realME(partons,realmom,inter);
  double ratio = reme/lome*sqr(inpart.mass())*4.*Constants::pi;
  if(inter==ShowerInteraction::QCD) ratio *= CF_;
  return ratio;
}

double SMWDecayer::meRatio(tcPDVector partons, 
			   vector<Lorentz5Momentum> momenta,
			   unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  double lome(0.);
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;    
    Energy2 pipj = momenta[3      ] * momenta[1+iemit ];
    Energy2 pipk = momenta[3      ] * momenta[1+ispect];
    Energy2 pjpk = momenta[1+iemit] * momenta[1+ispect];
    double y = pipj/(pipj+pipk+pjpk); 
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[1+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[1+ispect].mass()))*
			(Q2-sqr(mij-momenta[1+ispect].mass())));
    Energy2 Qpk = q*momenta[1+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[1+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[1+ispect].mass())-sqr(momenta[1+ispect].mass()))*q;
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[1+iemit ].mass()/sqrt(Q2);
    double muk = momenta[1+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[1+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(3);
    lomom[0] = momenta[0];
    if(iemit==0) {
      lomom[1] = pijt;
      lomom[2] = pkt ;
    }
    else {
      lomom[2] = pijt;
      lomom[1] = pkt ;
    }
    if(iemit==0) lome  = loME(partons,lomom);
  }
  InvEnergy2 ratio = realME(partons,momenta,ShowerInteraction::QCD)/lome*abs(D[iemitter])
    /(abs(D[0])+abs(D[1]));
  if(subtract)
    return Q2*(ratio-2.*D[iemitter]);
  else
    return Q2*ratio;
}

double SMWDecayer::loME(const vector<tcPDPtr> & partons, 
			const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>    vin;
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  VectorWaveFunction    win  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    win.reset(ix);
    vin.push_back(win);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total(0.);
  for(unsigned int inhel=0;inhel<3;++inhel) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	Complex diag1 = FFWVertex_->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMWDecayer::realME(const vector<tcPDPtr> & partons, 
			      const vector<Lorentz5Momentum> & momenta,
			      ShowerInteraction inter) const {
  // compute the spinors
  vector<VectorWaveFunction>     vin;
  vector<SpinorWaveFunction>     aout;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction>     gout;
  VectorWaveFunction    win  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  VectorWaveFunction    gluon(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  for(unsigned int ix=0;ix<3;++ix){
    win.reset(ix);
    vin.push_back(win);
  }
  vector<Complex> diag(3,0.);

  double total(0.);

  AbstractFFVVertexPtr vertex = inter==ShowerInteraction::QCD ? FFGVertex_ : FFPVertex_;
  
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    vertex->evaluate(scale_,3,partons[1]->CC(),fout[outhel1],gout[outhel3]);
	  diag[0] = FFWVertex_->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);
	  
	  SpinorWaveFunction off2 = 
	    vertex->evaluate(scale_,3,partons[2]->CC(),aout[outhel2],gout[outhel3]);
	  diag[1] = FFWVertex_->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

	  if(inter==ShowerInteraction::QED) {
	    VectorWaveFunction off3 =
	      WWWVertex_->evaluate(scale_,3,partons[0],vin[inhel1],gout[outhel3]);
	    diag[2] = FFWVertex_->evaluate(scale_,aout[outhel2],fout[outhel1],off3);
	  }
	  
	  // sum of diagrams
	  Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  // me2
	  total += norm(sum);
	}
      }
    }
  }
  // divide out the coupling
  total /= norm(vertex->norm());
  // double g = sqrt(2.)*abs(FFWVertex_->norm());
  // double xg = 2.*momenta[3].t()/momenta[0].mass();
  // double xe,mue2;
  // if(abs(partons[1]->id())==ParticleID::eminus) {
  //   xe = 2.*momenta[1].t()/momenta[0].mass();
  //   mue2 = sqr(momenta[1].mass()/momenta[0].mass());
  // }
  // else {
  //   xe = 2.*momenta[2].t()/momenta[0].mass();
  //   mue2 = sqr(momenta[2].mass()/momenta[0].mass());
  // }
  // double cg = -4. * g * g * (-pow(mue2, 3.) / 2. + (xg * xg / 4. + (xe / 2. + 1.) * xg + 5. / 2. * xe - 2.) * mue2 * mue2
  // 			     + (pow(xg, 3.) / 4. + (xe / 4. - 5. / 4.) * xg * xg + (-7. / 2. * xe + 3.) * xg - 3. * xe * xe
  // 				+ 11. / 2. * xe - 7. / 2.) * mue2 + (xg * xg / 2. + (xe - 2.) * xg + xe * xe - 2. * xe + 2.) * (-1. + xg + xe)) * (xe - mue2 - 1.) *
  //   pow(xg, -2.) * pow(-1. + xg + xe - mue2, -2.);
  
  // cerr << "real " << cg/total << "\n";
  // return the total
  return total*UnitRemoval::InvE2;
}

double SMWDecayer::calculateRealEmission(double x1, double x2, 
					 tcPDVector partons,
					 vector<Lorentz5Momentum> pin,
					 double phi, double muj,
					 double muk, int iemit, 
					 bool subtract) const {
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3)-0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1)-4.*sqr(muk)+4.*sqr(muj))
		       /(sqr(x2)-4.*sqr(muk))));
  // calculate the momenta
  Energy M = mW_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*sqr(muk),0.)),
			  0.5*M*x2,M*muk); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*sqr(muj),0.)),
			  0.5*M*x1,M*muj);
  Lorentz5Momentum pgluon(0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // boost and rotate momenta
  LorentzRotation eventFrame( ( pin[1] +
				pin[2] ).findBoostToCM() );
  unsigned int ispect = iemit==0 ? 2 : 1;
  Lorentz5Momentum spectator = eventFrame*pin[ispect];
  eventFrame.rotateZ( -spectator.phi()    );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  vector<Lorentz5Momentum> momenta(3);
  momenta[0]   = pin[0];
  momenta[ispect ] = eventFrame*pspect;
  momenta[iemit+1] = eventFrame*pemit ;
  momenta.push_back(eventFrame*pgluon);
  // calculate the weight
  double realwgt(0.);
  if(1.-x1>1e-5 && 1.-x2>1e-5) 
    realwgt = meRatio(partons,momenta,iemit,subtract);
  return realwgt;
}
#line 1 "./SMZDecayer.cc"
// -*- C++ -*-
//
// SMZDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMZDecayer class.
//

#include "SMZDecayer.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG::Helicity;

const double SMZDecayer::EPS_=0.00000001;

SMZDecayer::SMZDecayer() 
  : quarkWeight_(5,0.), leptonWeight_(6,0.), CF_(4./3.),
    NLO_(false) {
   quarkWeight_[0]  = 0.488029;
   quarkWeight_[1]  = 0.378461;
   quarkWeight_[2]  = 0.488019;
   quarkWeight_[3]  = 0.378027;
   quarkWeight_[4]  = 0.483207;
   leptonWeight_[0] = 0.110709;
   leptonWeight_[1] = 0.220276;
   leptonWeight_[2] = 0.110708;
   leptonWeight_[3] = 0.220276;
   leptonWeight_[4] = 0.110458;
   leptonWeight_[5] = 0.220276;
   // intermediates
   generateIntermediates(false);
   // QED corrections
  hasRealEmissionME(true);
  hasOneLoopME(true);
}

void SMZDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig StandardModel object in"
  				  << "SMZDecayer::doinit()"
  				  << Exception::runerror;
  // cast the vertices
  FFZVertex_ = dynamic_ptr_cast<FFVVertexPtr>(hwsm->vertexFFZ());
  FFZVertex_->init();
  FFGVertex_ = hwsm->vertexFFG();
  FFGVertex_->init();
  FFPVertex_ = hwsm->vertexFFP();
  FFPVertex_->init();
  gluon_ = getParticleData(ParticleID::g);
  // now set up the decay modes
  tPDPtr Z0 = getParticleData(ParticleID::Z0);
  // loop over the  quarks and the leptons
  for(int istep=0;istep<11;istep+=10) {
    for(int ix=1;ix<7;++ix) {
      int iy=istep+ix;
      if(iy==6) continue;
      double maxWeight = iy<=6 ? quarkWeight_.at(ix-1) : leptonWeight_.at(iy-11);
      tPDVector out = {getParticleData(-iy),getParticleData( iy)};
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(Z0,out,maxWeight));
      addMode(mode);
    }
  }
}

int SMZDecayer::modeNumber(bool & cc,tcPDPtr parent, 
			    const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2) return imode;
  int id0=parent->id();
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  // Z to quarks or leptons
  cc =false;
  if(id0!=ParticleID::Z0) return imode;
  if(abs(id1)<6&&id1==-id2) {
    imode=abs(id1)-1;
  }
  else if(abs(id1)>=11&&abs(id1)<=16&&id1==-id2) {
    imode=abs(id1)-6;
  }
  cc = false;
  return imode;
}

void SMZDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ << FFPVertex_ << FFGVertex_
     << quarkWeight_ << leptonWeight_ << NLO_
     << gluon_;
}

void SMZDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFGVertex_
     >> quarkWeight_ >> leptonWeight_ >> NLO_
     >> gluon_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMZDecayer,PerturbativeDecayer>
describeHerwigSMZDecayer("Herwig::SMZDecayer", "HwPerturbativeDecay.so");

void SMZDecayer::Init() {

  static ClassDocumentation<SMZDecayer> documentation
    ("The SMZDecayer class is the implementation of the decay"
     " Z boson to the Standard Model fermions.");

  static ParVector<SMZDecayer,double> interfaceZquarkMax
    ("QuarkMax",
     "The maximum weight for the decay of the Z to quarks",
     &SMZDecayer::quarkWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static ParVector<SMZDecayer,double> interfaceZleptonMax
    ("LeptonMax",
     "The maximum weight for the decay of the Z to leptons",
     &SMZDecayer::leptonWeight_,
     0, 0, 0, -10000, 10000, false, false, true);

  static Switch<SMZDecayer,bool> interfaceNLO
    ("NLO",
     "Whether to return the LO or NLO result",
     &SMZDecayer::NLO_, false, false, false);
  static SwitchOption interfaceNLOLO
    (interfaceNLO,
     "No",
     "Leading-order result",
     false);
  static SwitchOption interfaceNLONLO
    (interfaceNLO,
     "Yes",
     "NLO result",
     true);

}

ParticleVector SMZDecayer::decay(const Particle & parent,
				 const tPDVector & children) const {
  // generate the decay
  bool cc;
  unsigned int imode = modeNumber(cc,parent.dataPtr(),children);
  ParticleVector output = generate(false,false,imode,parent);
  if(output[0]->hasColour())      output[0]->antiColourNeighbour(output[1]);
  else if(output[1]->hasColour()) output[1]->antiColourNeighbour(output[0]);
  return output;
}

void SMZDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  VectorWaveFunction::constructSpinInfo(_vectors,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  SpinorBarWaveFunction::
    constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
}

// return the matrix element squared
double SMZDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(outgoing[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(_vectors,_rho,
  					       const_ptr_cast<tPPtr>(&part),
  					       incoming,false);
    // fix rho if no correlations
    fixRho(_rho);
  }
  SpinorBarWaveFunction wbar(momenta[iferm],outgoing[iferm],Helicity::outgoing);
  SpinorWaveFunction    w   (momenta[ianti],outgoing[ianti],Helicity::outgoing);
  _wavebar.resize(2);
  _wave   .resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wbar.reset(ihel);
    _wavebar[ihel] = wbar;
    w.reset(ihel);
    _wave   [ihel] = w;
  }
  // compute the matrix element
  Energy2 scale(sqr(part.mass()));
  unsigned int ifm,ia,vhel;
  for(ifm=0;ifm<2;++ifm) {
    for(ia=0;ia<2;++ia) {
      for(vhel=0;vhel<3;++vhel) {
  	if(iferm>ianti) (*ME())(vhel,ia,ifm)=
  	  FFZVertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
  	else            (*ME())(vhel,ifm,ia)=
  	  FFZVertex_->evaluate(scale,_wave[ia],_wavebar[ifm],_vectors[vhel]);
      }
    }
  }
  double output=(ME()->contract(_rho)).real()*UnitRemoval::E2/scale;
  if(abs(outgoing[0]->id())<=6) output*=3.;
  // if LO return
  if(!NLO_) return output;  // check decay products coloured, otherwise return
  if(!outgoing[0]->coloured()) return output;
  // inital masses, couplings  etc
  // fermion mass
  Energy particleMass = outgoing[0]->mass();
  // Z mass
  mZ_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mZ_));
  // reduced mass
  mu_  = particleMass/mZ_;
  mu2_ = sqr(mu_);
  // scale
  scale_ = sqr(mZ_);
  // compute the spinors
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  vector<VectorWaveFunction>    vin;
  SpinorBarWaveFunction qkout(momenta[0],outgoing[0],Helicity::outgoing);
  SpinorWaveFunction    qbout(momenta[1],outgoing[1],Helicity::outgoing);
  VectorWaveFunction    zin  (part.momentum()     ,part.dataPtr()     ,incoming);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total=0.;
  if(mu_!=0.) {
    LorentzPolarizationVector momDiff = 
      (momenta[0]-momenta[1])/2./(momenta[0].mass()+momenta[1].mass());
    // scalars
    Complex scalar1 = zin.wave().dot(momDiff);
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {		
  	for(unsigned int inhel=0;inhel<3;++inhel) {
  	  // first the LO bit
  	  Complex diag1 = FFZVertex_->evaluate(scale_,aout[outhel2],
  					       fout[outhel1],vin[inhel]);
  	  // extra stuff for NLO 
  	  LorentzPolarizationVector left  = 
  	    aout[outhel2].wave().leftCurrent(fout[outhel1].wave());
  	  LorentzPolarizationVector right = 
  	    aout[outhel2].wave().rightCurrent(fout[outhel1].wave());
  	  Complex scalar = 
  	    aout[outhel2].wave().scalar(fout[outhel1].wave());
  	  // nlo specific pieces
  	  Complex diag3 =
  	    Complex(0.,1.)*FFZVertex_->norm()*
  	    (FFZVertex_->right()*( left.dot(zin.wave())) +
  	     FFZVertex_-> left()*(right.dot(zin.wave())) -
  	     ( FFZVertex_-> left()+FFZVertex_->right())*scalar1*scalar);
  	  // nlo piece
  	  total += real(diag1*conj(diag3) + diag3*conj(diag1));
  	}
      }
    }
    // rescale
    total *= UnitRemoval::E2/scale_;
  }
  else {
    total = ZERO;
  }
  // now for the NLO bit
  double mu4 = sqr(mu2_);
  double lmu = mu_!=0. ? log(mu_) : 0.;
  double v = sqrt(1.-4.*mu2_),v2(sqr(v));
  double omv = 4.*mu2_/(1.+v);
  double f1,f2,fNS,VNS;
  double r = omv/(1.+v);
  double lr = mu_!=0. ? log(r) : 0.;
  // normal form
  if(mu_>1e-4) {
    f1 = CF_*aS_/Constants::pi*
      ( +1. + 3.*log(0.5*(1.+v)) - 1.5*log(0.5*(1.+v2)) + sqr(Constants::pi)/6.
  	- 0.5*sqr(lr) - (1.+v2)/v*(lr*log(1.+v2) + sqr(Constants::pi)/12. 
  				       -0.5*log(4.*mu2_)*lr + 0.25*sqr(lr)));
    fNS = -0.5*(1.+2.*v2)*lr/v + 1.5*lr - 2./3.*sqr(Constants::pi) + 0.5*sqr(lr)
      + (1.+v2)/v*(Herwig::Math::ReLi2(r) + sqr(Constants::pi)/3. - 0.25*sqr(lr) 
  		   + lr*log((2.*v/ (1.+v))));
    VNS = 1.5*log(0.5*(1.+v2)) 
      + 0.5*(1.+v2)/v*( 2.*lr*log(2.*(1.+v2)/sqr(1.+v))  
  			+ 2.*Herwig::Math::ReLi2(sqr(r)) 
  			- 2.*Herwig::Math::ReLi2(2.*v/(1.+v)) - sqr(Constants::pi)/6.)
      + log(1.-mu_) - 2.*log(1.-2.*mu_) - 4.*mu2_/(1.+v2)*log(mu_/(1.-mu_)) 
      - mu_/(1.-mu_)
      + 4.*(2.*mu2_-mu_)/(1.+v2) + 0.5*sqr(Constants::pi); 
    f2 = CF_*aS_/Constants::pi*mu2_*lr/v;
  }
  // small mass limit
  else {
    f1 = -CF_*aS_/Constants::pi/6.*
      ( - 6. - 24.*lmu*mu2_ - 15.*mu4 - 12.*mu4*lmu - 24.*mu4*sqr(lmu) 
  	+ 2.*mu4*sqr(Constants::pi) - 12.*mu2_*mu4 - 96.*mu2_*mu4*sqr(lmu) 
  	+ 8.*mu2_*mu4*sqr(Constants::pi) - 80.*mu2_*mu4*lmu);
    fNS = - mu2_/18.*( + 36.*lmu - 36. - 45.*mu2_ + 216.*lmu*mu2_ - 24.*mu2_*sqr(Constants::pi) 
  		      + 72.*mu2_*sqr(lmu) - 22.*mu4 + 1032.*mu4 * lmu
  		      - 96.*mu4*sqr(Constants::pi) + 288.*mu4*sqr(lmu));
    VNS = - mu2_/1260.*(-6930. + 7560.*lmu + 2520.*mu_ - 16695.*mu2_ 
  		       + 1260.*mu2_*sqr(Constants::pi) 
  		       + 12600.*lmu*mu2_ + 1344.*mu_*mu2_ - 52780.*mu4 + 36960.*mu4*lmu 
  		       + 5040.*mu4*sqr(Constants::pi) - 12216.*mu_*mu4);
    f2 = CF_*aS_*mu2_/Constants::pi*( 2.*lmu + 4.*mu2_*lmu + 2.*mu2_ + 12.*mu4*lmu + 7.*mu4);
  }
  // add up bits for f1
  f1 += CF_*aS_/Constants::pi*(fNS+VNS);
  double realFact(0.); 
  for(int iemit=0;iemit<2;++iemit) {
    // now for the real correction
    double phi  = UseRandom::rnd()*Constants::twopi;
    // calculate y
    double yminus = 0.; 
    double yplus  = 1.-2.*mu_*(1.-mu_)/(1.-2*mu2_);
    double y = yminus + UseRandom::rnd()*(yplus-yminus);
    // calculate z
    double v1  = sqrt(sqr(2.*mu2_+(1.-2.*mu2_)*(1.-y))-4.*mu2_)/(1.-2.*mu2_)/(1.-y);
    double zplus  = (1.+v1)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
    double zminus = (1.-v1)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
    double z = zminus + UseRandom::rnd()*(zplus-zminus);
    double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
    // calculate x1,x2
    double x2 = 1. - y*(1.-2.*mu2_);
    double x1 = 1. - z*(x2-2.*mu2_);
    // copy the particle objects over for calculateRealEmission
    tcPDVector oTemp = {part.dataPtr(),outgoing[0],outgoing[1],gluon_};
    vector<Lorentz5Momentum> mom = {part.momentum(),momenta[0],momenta[1]};
    // total real emission contribution
    realFact += 0.25*jac*sqr(1.-2.*mu2_)/
      sqrt(1.-4.*mu2_)/Constants::twopi
      *2.*CF_*aS_*calculateRealEmission(x1, x2, oTemp, mom,  phi,
  					iemit, true);
  }
  // the born + virtual + real
  output = output*(1. + f1 + realFact) + f2*total;
  return output;
}

void SMZDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<5)       quarkWeight_ [ix   ]=mode(ix)->maxWeight();
      else if(ix<11) leptonWeight_[ix-5 ]=mode(ix)->maxWeight();
    }
  }
}

void SMZDecayer::dataBaseOutput(ofstream & output,
				 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<quarkWeight_.size();++ix) {
    output << "newdef " << name() << ":QuarkMax " << ix << " "
	   << quarkWeight_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<leptonWeight_.size();++ix) {
    output << "newdef " << name() << ":LeptonMax " << ix << " "
	   << leptonWeight_[ix] << "\n";
  }
  // parameters for the PerturbativeDecayer base class
  PerturbativeDecayer::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

InvEnergy2 SMZDecayer::
realEmissionME(unsigned int,const Particle &parent, 
	       ParticleVector &children,
	       unsigned int iemitter,
	       double ctheta, double stheta,
	       const LorentzRotation & rot1,
	       const LorentzRotation & rot2) {
  // check the number of products and parent
  assert(children.size()==3 && parent.id()==ParticleID::Z0);
  // the electric charge
  double e = sqrt(SM().alphaEM()*4.*Constants::pi);
  // azimuth of the photon
  double phi = children[2]->momentum().phi();
  // wavefunctions for the decaying particle in the rotated dipole frame
  vector<VectorWaveFunction> vec1 = _vectors;
  for(unsigned int ix=0;ix<vec1.size();++ix) {
    vec1[ix].transform(rot1);
    vec1[ix].transform(rot2);
  }
  // wavefunctions for the decaying particle in the rotated rest frame
  vector<VectorWaveFunction> vec2 = _vectors;
  for(unsigned int ix=0;ix<vec1.size();++ix) {
    vec2[ix].transform(rot2);
  }
  // find the outgoing particle and antiparticle
  unsigned int iferm(0),ianti(1);
  if(children[iferm]->id()<0) swap(iferm,ianti);
  // wavefunctions for the particles before the radiation
  // wavefunctions for the outgoing fermion
  SpinorBarWaveFunction wavebartemp;
  Lorentz5Momentum ptemp =  - _wavebar[0].momentum();
  ptemp *= rot2;
  if(ptemp.perp()/ptemp.e()<1e-10) {
    ptemp.setX(ZERO);
    ptemp.setY(ZERO);
  }
  wavebartemp = SpinorBarWaveFunction(ptemp,_wavebar[0].particle(),outgoing);
  // wavefunctions for the outgoing antifermion
  SpinorWaveFunction wavetemp;
  ptemp =  - _wave[0].momentum();
  ptemp *= rot2;
  if(ptemp.perp()/ptemp.e()<1e-10) {
    ptemp.setX(ZERO);
    ptemp.setY(ZERO);
  }
  wavetemp = SpinorWaveFunction(ptemp,_wave[0].particle(),outgoing);
  // loop over helicities
  vector<SpinorWaveFunction> wf_old;
  vector<SpinorBarWaveFunction> wfb_old;
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wavetemp.reset(ihel);
    wf_old.push_back(wavetemp);
    wavebartemp.reset(ihel);
    wfb_old.push_back(wavebartemp);
  }
  // calculate the wave functions for the new fermions
  // ensure the momenta have pT=0
  for(unsigned int ix=0;ix<2;++ix) {
    Lorentz5Momentum ptemp = children[ix]->momentum();
    if(ptemp.perp()/ptemp.e()<1e-10) {
      ptemp.setX(ZERO);
      ptemp.setY(ZERO);
      children[ix]->set5Momentum(ptemp);
    }
  }
  // calculate the wavefunctions
  vector<SpinorBarWaveFunction> wfb;
  SpinorBarWaveFunction::calculateWaveFunctions(wfb,children[iferm],outgoing);
  vector<SpinorWaveFunction> wf;
  SpinorWaveFunction::calculateWaveFunctions   (wf ,children[ianti],outgoing);
  // wave functions for the photons
  vector<VectorWaveFunction> photon;
  VectorWaveFunction::calculateWaveFunctions(photon,children[2],outgoing,true);
  // loop to calculate the matrix elements
  Complex lome[3][2][2],diffme[3][2][2][2],summe[3][2][2][2];
  Energy2 scale(sqr(parent.mass()));
  Complex diff[2]={0.,0.};
  Complex sum [2]={0.,0.};
  for(unsigned int ifm=0;ifm<2;++ifm) {
    for(unsigned int ia=0;ia<2;++ia) {
      for(unsigned int vhel=0;vhel<3;++vhel) {
	// calculation of the leading-order matrix element
	Complex loamp  = FFZVertex_->evaluate(scale,wf_old[ia],
					      wfb_old[ifm],vec2[vhel]);
	Complex lotemp = FFZVertex_->evaluate(scale,wf[ia],
					      wfb[ifm],vec1[vhel]);
	if(iferm>ianti) lome[vhel][ia][ifm] = loamp;
	else            lome[vhel][ifm][ia] = loamp;
	// photon loop for the real emmision terms
	for(unsigned int phel=0;phel<2;++phel) {
	  // radiation from the antifermion
	  // normal case with small angle treatment
	  if(children[2    ]->momentum().z()/
	     children[iferm]->momentum().z()>=ZERO && iemitter == iferm ) {
	    Complex dipole = e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    // sum and difference
	    SpinorBarWaveFunction foff =
	      FFPVertex_->evaluateSmall(ZERO,3,children[iferm]->dataPtr()->CC(),
					wfb[ifm],photon[2*phel],
					ifm,2*phel,ctheta,phi,stheta,false);
	    diff[0] = FFZVertex_->evaluate(scale,wf[ia],foff,vec1[vhel]) +
	      e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*(lotemp-loamp)*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    sum [0] = diff[0]+2.*dipole;
	  }
	  // special if fermion backwards
	  else {
	    SpinorBarWaveFunction foff = 
	      FFPVertex_->evaluate(ZERO,3,children[iferm]->dataPtr()->CC(),
				   wfb[ifm],photon[2*phel]);
	    Complex diag = 
	      FFZVertex_->evaluate(scale,wf[ia],foff,vec1[vhel]);
	    Complex dipole = e*double(children[iferm]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[iferm]->momentum()*photon[2*phel].wave())/
	      (children[iferm]->momentum()*children[2]->momentum());
	    diff[0] = diag-dipole;
	    sum [0] = diag+dipole;
	  }
	  // radiation from the anti fermion 
	  // small angle case in general
	  if(children[2    ]->momentum().z()/
	     children[ianti]->momentum().z()>=ZERO && iemitter == ianti ) {
	    Complex dipole = e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    // sum and difference
	    SpinorWaveFunction foff =
	      FFPVertex_->evaluateSmall(ZERO,3,children[ianti]->dataPtr()->CC(),
					wf[ia],photon[2*phel],
					ia,2*phel,ctheta,phi,stheta,false);
	    diff[1] = FFZVertex_->evaluate(scale,foff ,wfb[ifm],vec1[vhel]) +
	      e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*(lotemp-loamp)*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    sum [1] = diff[1]+2.*dipole;
	  }	    
	  // special if fermion backwards after radiation
	  else {
	    SpinorWaveFunction foff = 
	      FFPVertex_->evaluate(ZERO,3,children[ianti]->dataPtr()->CC(),
				   wf[ia],photon[2*phel]);
	    Complex diag = 
	      FFZVertex_->evaluate(scale,foff ,wfb[ifm],vec1[vhel]);
	    Complex dipole = e*double(children[ianti]->dataPtr()->iCharge())/3.*
	      UnitRemoval::E*loamp*
	      (children[ianti]->momentum()*photon[2*phel].wave())/
	      (children[ianti]->momentum()*children[2]->momentum());
	    // sum and difference
	    diff[1] = diag - dipole;
	    sum [1] = diag + dipole;
	  }
	  // add to me
	  if(iferm>ianti) {
	    diffme[vhel][ia][ifm][phel] = diff[0] + diff[1];
	    summe [vhel][ia][ifm][phel] = sum[0]  + sum[1] ;
	  }
	  else {
	    diffme  [vhel][ifm][ia][phel] = diff[0] + diff[1];
	    summe   [vhel][ifm][ia][phel] = sum[0]  + sum[1] ;
	  }
	}
      }
    }
  }
//   cerr << parent << "\n";
//   for(unsigned int ix=0;ix<children.size();++ix) {
//     cerr << *children[ix] << "\n";
//   }
//   _rho = RhoDMatrix(PDT::Spin1);
  Complex lo(0.),difference(0.);
  for(unsigned int vhel1=0;vhel1<3;++vhel1) {
    for(unsigned int vhel2=0;vhel2<3;++vhel2) {
      for(unsigned int ifm=0;ifm<2;++ifm) {
	for(unsigned int ia=0;ia<2;++ia) {
	  lo += _rho(vhel1,vhel2)*lome[vhel1][ifm][ia]*conj(lome[vhel2][ifm][ia]);
	  for(unsigned int phel=0;phel<2;++phel) {
	    difference += 
	      _rho(vhel1,vhel2)*diffme[vhel1][ifm][ia][phel]*conj(summe[vhel2][ifm][ia][phel]);
	  }
	}
      }
    }
  }
//   // analytic result
//   double iCharge = children[0]->dataPtr()->iCharge()*
//     children[1]->dataPtr()->iCharge()/9.;
//   Energy2 ubar = 2.*children[0]->momentum()*children[2]->momentum();
//   Energy2 tbar = 2.*children[1]->momentum()*children[2]->momentum();
//   double mu2 = sqr(children[1]->mass()/parent.mass());
//   double gL = (FFZVertex_->left() *FFZVertex_->norm()).real();
//   double gR = (FFZVertex_->right()*FFZVertex_->norm()).real();
//   Energy2 den = sqr(parent.mass())*(((sqr(gL)+sqr(gR))*(1-mu2)+6.*mu2*gL*gR));

//   InvEnergy2 anal =  -iCharge*( 2.*(ubar/tbar+tbar/ubar)/sqr(parent.mass())+
// 				4.*mu2/den*((sqr(gL)+sqr(gR))*(1+ubar/tbar+tbar/ubar)
// 					    -2.*gL*gR*(1.+2.*(ubar/tbar+tbar/ubar))));
//   cerr << "testing ratio " << parent.PDGName() 
//        << " " << difference.real()/sqr(e)/lo.real()*UnitRemoval::InvE2/(anal) << "\n"
//        << stheta << " " << ctheta << "\n";
  return difference.real()/sqr(e)/lo.real()*UnitRemoval::InvE2;
}

double SMZDecayer::oneLoopVirtualME(unsigned int,
				    const Particle & parent, 
				    const ParticleVector & children) {
  assert(children.size()==2);
  // velocities of the particles
  double beta = sqrt(1.-4.*sqr(children[0]->mass()/parent.mass()));
  double opb = 1.+beta;
  double omb = 4.*sqr(children[0]->mass()/parent.mass())/opb;
  // couplings
  double gL = (FFZVertex_->left() *FFZVertex_->norm()).real();
  double gR = (FFZVertex_->right()*FFZVertex_->norm()).real();
  double gA = 0.5*(gL-gR);
  double gV = 0.5*(gL+gR);
  // correction terms
  double ln = log(omb/opb);
  double f1 = 1. + ln*beta;
  double fA = 1. + ln/beta;
  InvEnergy f2 = 0.5*sqrt(omb*opb)/parent.mass()/beta*ln;
  // momentum difference for the loop
  Lorentz5Momentum q = children[0]->momentum()-children[1]->momentum();
  if(children[0]->id()<0) q *= -1.;
  // spinors
  vector<LorentzSpinor   <SqrtEnergy> > sp;
  vector<LorentzSpinorBar<SqrtEnergy> > sbar;
  for(unsigned int ix=0;ix<2;++ix) {
    sp  .push_back(   _wave[ix].dimensionedWave());
    sbar.push_back(_wavebar[ix].dimensionedWave());
  }
  // polarization vectors
  vector<LorentzPolarizationVector> pol;
  for(unsigned int ix=0;ix<3;++ix)
    pol.push_back(_vectors[ix].wave());
  // matrix elements
  complex<Energy> lome[3][2][2],loopme[3][2][2];
  for(unsigned int vhel=0;vhel<3;++vhel) {
    for(unsigned int ihel1=0;ihel1<2;++ihel1) {
      for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	complex<Energy> vector = 
	  sp[ihel1].generalCurrent(sbar[ihel2], 1.,1.).dot(pol[vhel]);
	complex<Energy>  axial = 
	  sp[ihel1].generalCurrent(sbar[ihel2],-1.,1.).dot(pol[vhel]);
	complex<Energy2> scalar =
	  sp[ihel1].scalar(sbar[ihel2])*(q*pol[vhel]);
	lome  [vhel][ihel1][ihel2] = gV*   vector-gA*   axial;
	loopme[vhel][ihel1][ihel2] = gV*f1*vector-gA*fA*axial+scalar*f2*gV;
      }
    }
  }
  // sum sums
  complex<Energy2> den(ZERO),num(ZERO);
  for(unsigned int vhel1=0;vhel1<3;++vhel1) {
    for(unsigned int vhel2=0;vhel2<3;++vhel2) {
      for(unsigned int ihel1=0;ihel1<2;++ihel1) {
	for(unsigned int ihel2=0;ihel2<2;++ihel2) {
	  num += _rho(vhel1,vhel2)*
	    (  lome[vhel1][ihel1][ihel2]*conj(loopme[vhel2][ihel1][ihel2])+
	     loopme[vhel1][ihel1][ihel2]*conj(  lome[vhel2][ihel1][ihel2]));
	  den += _rho(vhel1,vhel2)*
	    lome[vhel1][ihel1][ihel2]*conj(lome[vhel2][ihel1][ihel2]);
	}
      }
    }
  }
  // prefactor
  double iCharge = children[0]->dataPtr()->iCharge()*
                   children[1]->dataPtr()->iCharge()/9.;
  double pre = 0.5*SM().alphaEM()*iCharge/Constants::pi;
  // output
  return pre*num.real()/den.real();
}


void SMZDecayer::
initializeMECorrection(RealEmissionProcessPtr born, double & initial,
		       double & final) {
  // get the quark and antiquark
  ParticleVector qq; 
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix)
    qq.push_back(born->bornOutgoing()[ix]);
  // ensure quark first
  if(qq[0]->id()<0) swap(qq[0],qq[1]);
  // centre of mass energy
  d_Q_ = (qq[0]->momentum() + qq[1]->momentum()).m();
  // quark mass
  d_m_ = 0.5*(qq[0]->momentum().m()+qq[1]->momentum().m());
  // set the other parameters
  setRho(sqr(d_m_/d_Q_));
  setKtildeSymm();
  // otherwise can do it
  initial=1.;
  final  =1.;
}

bool SMZDecayer::softMatrixElementVeto(PPtr parent,
				       PPtr progenitor,
				       const bool & ,
				       const Energy & highestpT,
				       const vector<tcPDPtr> & ids,
				       const double & d_z,
				       const Energy & d_qt,
				       const Energy &) {
  // check we should be applying the veto
  if(parent->id()!=progenitor->id()||
     ids[0]->id()!=ids[1]->id()||
     ids[2]->id()!=ParticleID::g) return false;
  // calculate pt
  Energy2 d_m2 = parent->momentum().m2();
  Energy pPerp = (1.-d_z)*sqrt( sqr(d_z*d_qt) - d_m2);
  // if not hardest so far don't apply veto
  if(pPerp<highestpT) return false;
  // calculate the weight
  double weight = 0.;
  if(parent->id()>0) weight = qWeightX(d_qt, d_z);
  else weight = qbarWeightX(d_qt, d_z);
  // compute veto from weight and return 
  return !UseRandom::rndbool(weight);
}


void SMZDecayer::setRho(double r) 
{ 
  d_rho_ = r;
  d_v_ = sqrt(1.-4.*d_rho_);
}

void SMZDecayer::setKtildeSymm() { 
  d_kt1_ = (1. + sqrt(1. - 4.*d_rho_))/2.;
  setKtilde2();
}

void SMZDecayer::setKtilde2() { 
   double num = d_rho_ * d_kt1_ + 0.25 * d_v_ *(1.+d_v_)*(1.+d_v_);
   double den = d_kt1_ - d_rho_;
   d_kt2_ = num/den;
}

double SMZDecayer::getZfromX(double x1, double x2) {
  double uval = u(x2);
  double num = x1 - (2. - x2)*uval;
  double den = sqrt(x2*x2 - 4.*d_rho_);
  return uval + num/den;
}

double SMZDecayer::getKfromX(double x1, double x2) {
   double zval = getZfromX(x1, x2);
   return (1.-x2)/(zval*(1.-zval));
}

double SMZDecayer::MEV(double x1, double x2) {
  // Vector part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    - 8.*d_rho_*(1.+2.*d_rho_);
  double den = (1.+2.*d_rho_)*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMZDecayer::MEA(double x1, double x2) {
  // Axial part
  double num = (x1+2.*d_rho_)*(x1+2.*d_rho_) + (x2+2.*d_rho_)*(x2+2.*d_rho_) 
    + 2.*d_rho_*((5.-x1-x2)*(5.-x1-x2) - 19.0 + 4*d_rho_);
  double den = d_v_*d_v_*(1.-x1)*(1.-x2);
  return (num/den - 2.*d_rho_/((1.-x1)*(1.-x1)) 
	  - 2*d_rho_/((1.-x2)*(1.-x2)))/d_v_;
}

double SMZDecayer::u(double x2) {
  return 0.5*(1. + d_rho_/(1.-x2+d_rho_));
}

void SMZDecayer::
getXXbar(double kti, double z, double &x, double &xbar) {
  double w = sqr(d_v_) + kti*(-1. + z)*z*(2. + kti*(-1. + z)*z);
  if (w < 0) {
    x = -1.; 
    xbar = -1;
  } else {
    x = (1. + sqr(d_v_)*(-1. + z) + sqr(kti*(-1. + z))*z*z*z 
	 + z*sqrt(w)
	 - kti*(-1. + z)*z*(2. + z*(-2 + sqrt(w))))/
      (1. - kti*(-1. + z)*z + sqrt(w));
    xbar = 1. + kti*(-1. + z)*z;
  }
}

double SMZDecayer::qWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the quark emission zone?
  if(k1 < d_kt1_) {
    rval = MEV(x, xbar)/PS(x, xbar);
    // is it also in the anti-quark emission zone?
    if(k2 < d_kt2_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMZDecayer::qbarWeight(double x, double xbar) {
  double rval; 
  double xg = 2. - xbar - x;
  // always return one in the soft gluon region
  if(xg < EPS_) return 1.0;
  // check it is in the phase space
  if((1.-x)*(1.-xbar)*(1.-xg) < d_rho_*xg*xg) return 0.0;
  double k1 = getKfromX(x, xbar);
  double k2 = getKfromX(xbar, x);
  // Is it in the antiquark emission zone?
  if(k2 < d_kt2_) {
    rval = MEV(x, xbar)/PS(xbar, x);
    // is it also in the quark emission zone?
    if(k1 < d_kt1_) rval *= 0.5;
    return rval;
  }
  return 1.0;
}

double SMZDecayer::qWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, x, xb);
  // if exceptionally out of phase space, leave this emission, as there 
  // is no good interpretation for the soft ME correction. 
  if (x < 0 || xb < 0) return 1.0; 
  return qWeight(x, xb); 
}

double SMZDecayer::qbarWeightX(Energy qtilde, double z) {
  double x, xb;
  getXXbar(sqr(qtilde/d_Q_), z, xb, x);
  // see above in qWeightX. 
  if (x < 0 || xb < 0) return 1.0; 
  return qbarWeight(x, xb); 
}

double SMZDecayer::PS(double x, double xbar) {
  double u = 0.5*(1. + d_rho_ / (1.-xbar+d_rho_));
  double z = u + (x - (2.-xbar)*u)/sqrt(xbar*xbar - 4.*d_rho_);
  double brack = (1.+z*z)/(1.-z)- 2.*d_rho_/(1-xbar);
  // interesting: the splitting function without the subtraction
  // term. Actually gives a much worse approximation in the collinear
  // limit.  double brack = (1.+z*z)/(1.-z);
  double den = (1.-xbar)*sqrt(xbar*xbar - 4.*d_rho_);
  return brack/den;
}

double SMZDecayer::matrixElementRatio(const Particle & part, const ParticleVector & decay2,
				      const ParticleVector & decay3, MEOption,
				      ShowerInteraction inter) {
  // extract partons and LO momentas
  tcPDVector partons(1,part.dataPtr());
  vector<Lorentz5Momentum> lomom(1,part.momentum());
  for(unsigned int ix=0;ix<2;++ix) {
    partons.push_back(decay2[ix]->dataPtr());
    lomom.push_back(decay2[ix]->momentum());
  }
  vector<Lorentz5Momentum> realmom(1,part.momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==2) partons.push_back(decay3[ix]->dataPtr());
    realmom.push_back(decay3[ix]->momentum());
  }
  if(partons[0]->id()<0) {
    swap(partons[1],partons[2]);
    swap(lomom[1],lomom[2]);
    swap(realmom[1],realmom[2]);
  }
  scale_ = sqr(part.mass());
  double     lome = loME(partons,lomom);
  InvEnergy2 reme = realME(partons,realmom,inter);
  double ratio = reme/lome*sqr(part.mass())*4.*Constants::pi;
  if(inter==ShowerInteraction::QCD) ratio *= CF_;
  return ratio;
}

double SMZDecayer::meRatio(tcPDVector partons, 
			   vector<Lorentz5Momentum> momenta,
			   unsigned int iemitter, bool subtract) const {
  Lorentz5Momentum q = momenta[1]+momenta[2]+momenta[3];
  Energy2 Q2=q.m2();
  Energy2 lambda = sqrt((Q2-sqr(momenta[1].mass()+momenta[2].mass()))*
			(Q2-sqr(momenta[1].mass()-momenta[2].mass())));
  InvEnergy2 D[2];
  double lome[2];
  for(unsigned int iemit=0;iemit<2;++iemit) {
    unsigned int ispect = iemit==0 ? 1 : 0;    
    Energy2 pipj = momenta[3      ] * momenta[1+iemit ];
    Energy2 pipk = momenta[3      ] * momenta[1+ispect];
    Energy2 pjpk = momenta[1+iemit] * momenta[1+ispect];
    double y = pipj/(pipj+pipk+pjpk); 
    double z = pipk/(     pipk+pjpk);
    Energy mij = sqrt(2.*pipj+sqr(momenta[1+iemit].mass()));
    Energy2 lamB = sqrt((Q2-sqr(mij+momenta[1+ispect].mass()))*
			(Q2-sqr(mij-momenta[1+ispect].mass())));
    Energy2 Qpk = q*momenta[1+ispect];
    Lorentz5Momentum pkt = 
      lambda/lamB*(momenta[1+ispect]-Qpk/Q2*q)
      +0.5/Q2*(Q2+sqr(momenta[1+ispect].mass())-sqr(momenta[1+ispect].mass()))*q;
    Lorentz5Momentum pijt = 
      q-pkt;
    double muj = momenta[1+iemit ].mass()/sqrt(Q2);
    double muk = momenta[1+ispect].mass()/sqrt(Q2);
    double vt = sqrt((1.-sqr(muj+muk))*(1.-sqr(muj-muk)))/(1.-sqr(muj)-sqr(muk));
    double v  = sqrt(sqr(2.*sqr(muk)+(1.-sqr(muj)-sqr(muk))*(1.-y))-4.*sqr(muk))
      /(1.-y)/(1.-sqr(muj)-sqr(muk));
    // dipole term
    D[iemit] = 0.5/pipj*(2./(1.-(1.-z)*(1.-y))
			 -vt/v*(2.-z+sqr(momenta[1+iemit].mass())/pipj));
    // matrix element
    vector<Lorentz5Momentum> lomom(3);
    lomom[0] = momenta[0];
    if(iemit==0) {
      lomom[1] = pijt;
      lomom[2] = pkt ;
    }
    else {
      lomom[2] = pijt;
      lomom[1] = pkt ;
    }
    lome[iemit]  = loME(partons,lomom);
  }
  InvEnergy2 ratio = realME(partons,momenta,ShowerInteraction::QCD)*abs(D[iemitter])
    /(abs(D[0]*lome[0])+abs(D[1]*lome[1]));
  if(subtract)
    return Q2*(ratio-2.*D[iemitter]);
  else
    return Q2*ratio;
}

double SMZDecayer::loME(const tcPDVector & partons, 
			const vector<Lorentz5Momentum> & momenta) const {
  // compute the spinors
  vector<VectorWaveFunction>    vin;
  vector<SpinorWaveFunction>    aout;
  vector<SpinorBarWaveFunction> fout;
  VectorWaveFunction    zin  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  // temporary storage of the different diagrams
  // sum over helicities to get the matrix element
  double total(0.);
  for(unsigned int inhel=0;inhel<3;++inhel) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	Complex diag1 = FFZVertex_->evaluate(scale_,aout[outhel2],fout[outhel1],vin[inhel]);
	total += norm(diag1);
      }
    }
  }
  // return the answer
  return total;
}
 
InvEnergy2 SMZDecayer::realME(const tcPDVector & partons, 
			      const vector<Lorentz5Momentum> & momenta,
			      ShowerInteraction inter) const {
  AbstractFFVVertexPtr vertex = inter==ShowerInteraction::QCD ?
    FFGVertex_ : FFPVertex_;
  // compute the spinors
  vector<VectorWaveFunction>     vin;
  vector<SpinorWaveFunction>     aout;
  vector<SpinorBarWaveFunction>  fout;
  vector<VectorWaveFunction>     gout;
  VectorWaveFunction    zin  (momenta[0],partons[0],incoming);
  SpinorBarWaveFunction qkout(momenta[1],partons[1],outgoing);
  SpinorWaveFunction    qbout(momenta[2],partons[2],outgoing);
  VectorWaveFunction    gluon(momenta[3],partons[3],outgoing);
  for(unsigned int ix=0;ix<2;++ix){
    qkout.reset(ix);
    fout.push_back(qkout);
    qbout.reset(ix);
    aout.push_back(qbout);
    gluon.reset(2*ix);
    gout.push_back(gluon);
  }
  for(unsigned int ix=0;ix<3;++ix){
    zin.reset(ix);
    vin.push_back(zin);
  }
  vector<Complex> diag(2,0.);

  double total(0.);
  for(unsigned int inhel1=0;inhel1<3;++inhel1) {
    for(unsigned int outhel1=0;outhel1<2;++outhel1) {
      for(unsigned int outhel2=0;outhel2<2;++outhel2) {
	for(unsigned int outhel3=0;outhel3<2;++outhel3) {
	  SpinorBarWaveFunction off1 =
	    vertex->evaluate(scale_,3,partons[1]->CC(),fout[outhel1],gout[outhel3]);
	  diag[0] = FFZVertex_->evaluate(scale_,aout[outhel2],off1,vin[inhel1]);

	  SpinorWaveFunction off2 = 
	    vertex->evaluate(scale_,3,partons[2]->CC(),aout[outhel2],gout[outhel3]);
	  diag[1] = FFZVertex_->evaluate(scale_,off2,fout[outhel1],vin[inhel1]);

	  // sum of diagrams
	  Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	  // me2
	  total += norm(sum);
	}
      }
    }
  }

  // divide out the coupling
  total /= norm(vertex->norm());
  // return the total
  return total*UnitRemoval::InvE2;
}

double SMZDecayer::calculateRealEmission(double x1, double x2, 
					 vector<PPtr> hardProcess,
					 double phi,
					 bool subtract) const {
  // make partons data object for meRatio
  tcPDVector partons (3);
  for(int ix=0; ix<3; ++ix)
    partons[ix] = hardProcess[ix]->dataPtr();
  partons.push_back(gluon_);
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2_)));
  // calculate the momenta
  Energy M = mZ_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2_,0.)),0.5*M*x2,M*mu_); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
  			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2_,0.)),0.5*M*x1,M*mu_);
  Lorentz5Momentum pgluon(0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
			  0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // loop over the possible emitting partons
  double realwgt(0.);
  for(unsigned int iemit=0;iemit<2;++iemit) {
    // boost and rotate momenta
    LorentzRotation eventFrame( ( hardProcess[1]->momentum() +
				  hardProcess[2]->momentum() ).findBoostToCM() );
    Lorentz5Momentum spectator = eventFrame*hardProcess[iemit+1]->momentum();
    eventFrame.rotateZ( -spectator.phi()    );
    eventFrame.rotateY( -spectator.theta()  );
    eventFrame.invert();
    vector<Lorentz5Momentum> momenta(3);
    momenta[0]   = hardProcess[0]->momentum();
    if(iemit==0) {
      momenta[2] = eventFrame*pspect;
      momenta[1] = eventFrame*pemit ;
    }
    else {
      momenta[1] = eventFrame*pspect;
      momenta[2] = eventFrame*pemit ;
    }
    momenta.push_back(eventFrame*pgluon);
    // calculate the weight
    if(1.-x1>1e-5 && 1.-x2>1e-5) 
      realwgt += meRatio(partons,momenta,iemit,subtract);
  }
  
  // total real emission contribution
  return realwgt;
}

double SMZDecayer::calculateRealEmission(double x1, double x2,
					 tcPDVector partons,
					 vector<Lorentz5Momentum> pin,
					 double phi,
					 bool subtract,
					 int emitter) const {
  // calculate x3
  double x3 = 2.-x1-x2;
  double xT = sqrt(max(0.,sqr(x3) -0.25*sqr(sqr(x2)+sqr(x3)-sqr(x1))/(sqr(x2)-4.*mu2_)));
  // calculate the momenta
  Energy M = mZ_;
  Lorentz5Momentum pspect(ZERO,ZERO,-0.5*M*sqrt(max(sqr(x2)-4.*mu2_,0.)),0.5*M*x2,M*mu_); 
  Lorentz5Momentum pemit (-0.5*M*xT*cos(phi),-0.5*M*xT*sin(phi),
  			  0.5*M*sqrt(max(sqr(x1)-sqr(xT)-4.*mu2_,0.)),0.5*M*x1,M*mu_);
  Lorentz5Momentum pgluon( 0.5*M*xT*cos(phi), 0.5*M*xT*sin(phi),
  			   0.5*M*sqrt(max(sqr(x3)-sqr(xT),0.)),0.5*M*x3,ZERO);
  if(abs(pspect.z()+pemit.z()-pgluon.z())/M<1e-6) 
    pgluon.setZ(-pgluon.z());
  else if(abs(pspect.z()-pemit.z()+pgluon.z())/M<1e-6) 
    pemit .setZ(- pemit.z());
  // boost and rotate momenta
  LorentzRotation eventFrame( ( pin[1]+pin[2] ).findBoostToCM() );
  unsigned int ispect =  emitter==0 ? 2 : 1;
  Lorentz5Momentum spectator = eventFrame*pin[ispect];
  eventFrame.rotateZ( -spectator.phi()    );
  eventFrame.rotateY( -spectator.theta()  );
  eventFrame.invert();
  vector<Lorentz5Momentum> momenta(3);
  momenta[0]   = pin[0];
  momenta[ispect   ] = eventFrame*pspect;
  momenta[emitter+1] = eventFrame*pemit ;
  momenta.push_back(eventFrame*pgluon);
  // calculate the weight
  double realwgt(0.);
  if(1.-x1>1e-5 && 1.-x2>1e-5) 
    realwgt = meRatio(partons,momenta,emitter,subtract);  
  // total real emission contribution
  return realwgt;
}
#line 1 "./SMTopDecayer.cc"
// -*- C++ -*-
//
// SMTopDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMTopDecayer class.
//

#include "SMTopDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOn1IntegralCalculator.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMTopDecayer::SMTopDecayer() 
  : _wquarkwgt(6,0.),_wleptonwgt(3,0.), _xg_sampling(1.5), 
    _initialenhance(1.),  _finalenhance(2.3) {
  _wleptonwgt[0] = 0.302583;
  _wleptonwgt[1] = 0.301024;
  _wleptonwgt[2] = 0.299548;
  _wquarkwgt[0]  = 0.851719;
  _wquarkwgt[1]  = 0.0450162;
  _wquarkwgt[2]  = 0.0456962;
  _wquarkwgt[3]  = 0.859839;
  _wquarkwgt[4]  = 3.9704e-06;
  _wquarkwgt[5]  = 0.000489657; 
  generateIntermediates(true);
}
  
bool SMTopDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  if(abs(parent->id()) != ParticleID::t) return false;
  int id0(0),id1(0),id2(0);
  for(tPDVector::const_iterator it = children.begin();
      it != children.end();++it) {
    int id=(**it).id(),absid(abs(id));
    if(absid==ParticleID::b&&double(id)/double(parent->id())>0) {
      id0=id;
    }
    else {
      switch (absid) {
      case ParticleID::nu_e: 
      case ParticleID::nu_mu:
      case ParticleID::nu_tau:
	id1 = id;
	break;
      case ParticleID::eminus:
      case ParticleID::muminus:
      case ParticleID::tauminus:
	id2 = id;
	break;
      case ParticleID::b:
      case ParticleID::d:
      case ParticleID::s:
	id1 = id;
	break;
      case ParticleID::u:
      case ParticleID::c:
	id2=id;
	break;
      default :
	break;
      }
    }
  }
  if(id0==0||id1==0||id2==0) return false;
  if(double(id1)/double(id2)>0) return false;
  return true;
}
  
ParticleVector SMTopDecayer::decay(const Particle & parent,
				   const tPDVector & children) const {
  int id1(0),id2(0);
  for(tPDVector::const_iterator it = children.begin();
      it != children.end();++it) {
    int id=(**it).id(),absid=abs(id);
    if(absid == ParticleID::b && double(id)/double(parent.id())>0) continue;
    //leptons
    if(absid > 10 && absid%2==0) id1=absid;
    if(absid > 10 && absid%2==1) id2=absid;
    //quarks
    if(absid < 10 && absid%2==0) id2=absid;
    if(absid < 10 && absid%2==1) id1=absid;
  }
  unsigned int imode(0);
  if(id2 >=11 && id2<=16) imode = (id1-12)/2;
  else imode = id1+1+id2/2;
  bool cc = parent.id() == ParticleID::tbar;
  ParticleVector out(generate(true,cc,imode,parent));
  //arrange colour flow
  PPtr pparent=const_ptr_cast<PPtr>(&parent);
  out[1]->incomingColour(pparent,out[1]->id()<0);
  ParticleVector products = out[0]->children();
  if(products[0]->hasColour())
    products[0]->colourNeighbour(products[1],true);
  else if(products[0]->hasAntiColour())
    products[0]->colourNeighbour(products[1],false);
  return out;
}   
 
void SMTopDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFWVertex_ << FFGVertex_ << FFPVertex_ << WWWVertex_
     << _wquarkwgt << _wleptonwgt << _wplus
     << _initialenhance << _finalenhance << _xg_sampling;
}
  
void SMTopDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFWVertex_ >> FFGVertex_ >> FFPVertex_ >> WWWVertex_
     >> _wquarkwgt >> _wleptonwgt >> _wplus
     >> _initialenhance >> _finalenhance >> _xg_sampling;
}
  
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMTopDecayer,PerturbativeDecayer>
describeHerwigSMTopDecayer("Herwig::SMTopDecayer", "HwPerturbativeDecay.so");
  
void SMTopDecayer::Init() {
    
  static ClassDocumentation<SMTopDecayer> documentation
    ("This is the implementation of the SMTopDecayer which "
     "decays top quarks into bottom quarks and either leptons  "
     "or quark-antiquark pairs including the matrix element for top decay",
     "The matrix element correction for top decay \\cite{Hamilton:2006ms}.",
     "%\\cite{Hamilton:2006ms}\n"
     "\\bibitem{Hamilton:2006ms}\n"
     "  K.~Hamilton and P.~Richardson,\n"
     "  ``A simulation of QCD radiation in top quark decays,''\n"
     "  JHEP {\\bf 0702}, 069 (2007)\n"
     "  [arXiv:hep-ph/0612236].\n"
     "  %%CITATION = JHEPA,0702,069;%%\n");
  
  static ParVector<SMTopDecayer,double> interfaceQuarkWeights
    ("QuarkWeights",
     "Maximum weights for the hadronic decays",
     &SMTopDecayer::_wquarkwgt, 6, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  
  static ParVector<SMTopDecayer,double> interfaceLeptonWeights
    ("LeptonWeights",
     "Maximum weights for the semi-leptonic decays",
     &SMTopDecayer::_wleptonwgt, 3, 1.0, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<SMTopDecayer,double> interfaceEnhancementFactor
    ("InitialEnhancementFactor",
     "The enhancement factor for initial-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one.",
     &SMTopDecayer::_initialenhance, 1.0, 1.0, 10000.0,
     false, false, Interface::limited);

  static Parameter<SMTopDecayer,double> interfaceFinalEnhancementFactor
    ("FinalEnhancementFactor",
     "The enhancement factor for final-state radiation in the shower to ensure"
     " the weight for the matrix element correction is less than one",
     &SMTopDecayer::_finalenhance, 1.6, 1.0, 1000.0,
     false, false, Interface::limited);

  static Parameter<SMTopDecayer,double> interfaceSamplingTopHardMEC
    ("SamplingTopHardMEC",
     "The importance sampling power for choosing an initial xg, "
     "to sample xg according to xg^-_xg_sampling",
     &SMTopDecayer::_xg_sampling, 1.5, 1.2, 2.0,
     false, false, Interface::limited);
}

void SMTopDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()>0) {
    SpinorWaveFunction::
      constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    SpinorWaveFunction   ::constructSpinInfo(_outHalf   ,decay[1],outgoing,true);
    SpinorBarWaveFunction::constructSpinInfo(_outHalfBar,decay[2],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    SpinorBarWaveFunction::constructSpinInfo(_outHalfBar,decay[1],outgoing,true);
    SpinorWaveFunction   ::constructSpinInfo(_outHalf   ,decay[2],outgoing,true);
  }
}

double SMTopDecayer::me2(const int,const Particle & part,
			 const tPDVector & outgoing,
			 const vector<Lorentz5Momentum> & momenta,
			 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
  					 PDT::Spin1Half,PDT::Spin1Half)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }
  if ( ( momenta[1] + momenta[2] ).m()
       < outgoing[1]->constituentMass() + outgoing[2]->constituentMass() )
    return 0.0;

  // spinors for the decay product
  _outHalf.resize(2);
  _outHalfBar.resize(2);
  if(part.id()>0) {
    SpinorBarWaveFunction w0(momenta[0],outgoing[0],Helicity::outgoing);
    SpinorWaveFunction    w1(momenta[1],outgoing[1],Helicity::outgoing);
    SpinorBarWaveFunction w2(momenta[2],outgoing[2],Helicity::outgoing);
    _inHalfBar.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel) {
      w0.reset(ihel);      _inHalfBar [ihel] = w0;
      w1.reset(ihel);      _outHalf   [ihel] = w1;
      w2.reset(ihel);      _outHalfBar[ihel] = w2;
    }
  }
  else {
    SpinorWaveFunction    w0(momenta[0],outgoing[0],Helicity::outgoing);
    SpinorBarWaveFunction w1(momenta[1],outgoing[1],Helicity::outgoing);
    SpinorWaveFunction    w2(momenta[2],outgoing[2],Helicity::outgoing);
    _inHalf.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel) {
      w0.reset(ihel);      _inHalf    [ihel] = w0;
      w1.reset(ihel);      _outHalfBar[ihel] = w1;
      w2.reset(ihel);      _outHalf   [ihel] = w2;
    }
  }
  Energy2 scale(sqr(part.mass()));
  if(part.id() == ParticleID::t) {
    //Define intermediate vector wave-function for Wplus 
    tcPDPtr Wplus(getParticleData(ParticleID::Wplus));
    VectorWaveFunction inter;
    unsigned int thel,bhel,fhel,afhel;
    for(thel = 0;thel<2;++thel){
      for(bhel = 0;bhel<2;++bhel){	  
  	inter = FFWVertex_->evaluate(scale,1,Wplus,_inHalf[thel],
  				   _inHalfBar[bhel]);
  	for(afhel=0;afhel<2;++afhel){
  	  for(fhel=0;fhel<2;++fhel){
  	    (*ME())(thel,bhel,afhel,fhel) = 
  	      FFWVertex_->evaluate(scale,_outHalf[afhel],
  				 _outHalfBar[fhel],inter);
  	  }
  	}
      }
    }
  }
  else if(part.id() == ParticleID::tbar) {
    VectorWaveFunction inter;
    tcPDPtr Wminus(getParticleData(ParticleID::Wminus));
    unsigned int tbhel,bbhel,afhel,fhel;
    for(tbhel = 0;tbhel<2;++tbhel){
      for(bbhel = 0;bbhel<2;++bbhel){
  	inter = FFWVertex_->
  	  evaluate(scale,1,Wminus,_inHalf[bbhel],_inHalfBar[tbhel]);
  	for(afhel=0;afhel<2;++afhel){
  	  for(fhel=0;fhel<2;++fhel){
  	    (*ME())(tbhel,bbhel,fhel,afhel) = 
  	      FFWVertex_->evaluate(scale,_outHalf[afhel],
  				 _outHalfBar[fhel],inter);
  	  }
  	}
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  if(abs(outgoing[1]->id())<=6) output *=3.;
  return output;
}

void SMTopDecayer::doinit() {
  PerturbativeDecayer::doinit();
  //get vertices from SM object
  tcHwSMPtr hwsm = dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) throw InitException() << "Must have Herwig::StandardModel in "
  				  << "SMTopDecayer::doinit()";
  FFWVertex_ = hwsm->vertexFFW();
  FFGVertex_ = hwsm->vertexFFG();
  FFPVertex_ = hwsm->vertexFFP();
  WWWVertex_ = hwsm->vertexWWW();
  //initialise
  FFWVertex_->init();
  FFGVertex_->init();
  FFPVertex_->init();
  WWWVertex_->init();
  //set up decay modes
  _wplus = getParticleData(ParticleID::Wplus);
  tPDPtr top = getParticleData(ParticleID::t);
  tPDPtr bot = getParticleData(ParticleID::b);
  //lepton modes
  for(int i=11; i<17;i+=2) {
    tPDVector out = {bot,getParticleData(-i),getParticleData(i+1)};
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(top,out,_wleptonwgt[(i-11)/2]));
    mode->addChannel((PhaseSpaceChannel(mode),0,_wplus,0,1,1,2,1,3));
    addMode(mode);
  }
  //quark modes
  unsigned int iz=0;
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<6;iy+=2) {
      // check that the combination of particles is allowed
      if(FFWVertex_->allowed(-ix,iy,ParticleID::Wminus)) {
	tPDVector out = {bot,getParticleData(-ix),getParticleData( iy)};
	PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(top,out,_wquarkwgt[iz]));
	mode->addChannel((PhaseSpaceChannel(mode),0,_wplus,0,1,1,2,1,3));
	addMode(mode);
 	++iz;
      }
      else {
  	throw InitException() << "SMTopDecayer::doinit() the W vertex" 
  			      << "cannot handle all the quark modes" 
  			      << Exception::abortnow;
      }
    }
  }
}

void SMTopDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the PerturbativeDecayer base class
  for(unsigned int ix=0;ix<_wquarkwgt.size();++ix) {
    os << "newdef " << name() << ":QuarkWeights " << ix << " "
	   << _wquarkwgt[ix] << "\n";
  }
  for(unsigned int ix=0;ix<_wleptonwgt.size();++ix) {
    os << "newdef " << name() << ":LeptonWeights " << ix << " "
	   << _wleptonwgt[ix] << "\n";
  }
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMTopDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(ix<3) _wleptonwgt[ix  ] = mode(ix)->maxWeight();
      else     _wquarkwgt [ix-3] = mode(ix)->maxWeight();
    }
  }
}

WidthCalculatorBasePtr SMTopDecayer::threeBodyMEIntegrator(const DecayMode & dm) const {
  // identify W decay products
  int sign = dm.parent()->id() > 0 ? 1 : -1;
  int iferm(0),ianti(0);
  for(ParticleMSet::const_iterator pit=dm.products().begin();
      pit!=dm.products().end();++pit) {
    int id = (**pit).id();
    if(id*sign != ParticleID::b) {
      if   (id*sign > 0 ) iferm = id*sign;
      else                ianti = id*sign;
    }
  }
  assert(iferm!=0&&ianti!=0);
  // work out which mode we are doing
  int imode(-1);
  for(unsigned int ix=0;ix<numberModes();++ix) {
    if(mode(ix)->outgoing()[1]->id() == ianti &&
       mode(ix)->outgoing()[2]->id() == iferm ) {
      imode = ix;
      break;
    }
  }
  assert(imode>=0);
  // get the masses we need
  Energy m[3] = {mode(imode)->outgoing()[0]->mass(),
		 mode(imode)->outgoing()[2]->mass(),
		 mode(imode)->outgoing()[1]->mass()};
  return 
    new_ptr(ThreeBodyAllOn1IntegralCalculator<SMTopDecayer>
	    (3,_wplus->mass(),_wplus->width(),0.0,*this,imode,m[0],m[1],m[2]));
}

InvEnergy SMTopDecayer::threeBodydGammads(const int imode, const Energy2 mt2,
					  const Energy2 mffb2, const Energy mb,
					  const Energy mf, const Energy mfb) const {
  Energy mffb(sqrt(mffb2));
  Energy mw(_wplus->mass());
  Energy2 mw2(sqr(mw)),gw2(sqr(_wplus->width()));
  Energy mt(sqrt(mt2));
  Energy Eb  = 0.5*(mt2-mffb2-sqr(mb))/mffb;
  Energy Ef  = 0.5*(mffb2-sqr(mfb)+sqr(mf))/mffb;
  Energy Ebm = sqrt(sqr(Eb)-sqr(mb));
  Energy Efm = sqrt(sqr(Ef)-sqr(mf));
  Energy2 upp = sqr(Eb+Ef)-sqr(Ebm-Efm);
  Energy2 low = sqr(Eb+Ef)-sqr(Ebm+Efm);
  InvEnergy width=(dGammaIntegrand(mffb2,upp,mt,mb,mf,mfb,mw)-
		   dGammaIntegrand(mffb2,low,mt,mb,mf,mfb,mw))
    /32./mt2/mt/8/pow(Constants::pi,3)/(sqr(mffb2-mw2)+mw2*gw2);
  // couplings
  width *= 0.25*sqr(4.*Constants::pi*generator()->standardModel()->alphaEM(mt2)/
		    generator()->standardModel()->sin2ThetaW());
  width *= generator()->standardModel()->CKM(*mode(imode)->incoming().first,
					     *mode(imode)->outgoing()[0]);
  if(abs(mode(imode)->outgoing()[1]->id())<=6) {
    width *=3.;
    if(abs(mode(imode)->outgoing()[1]->id())%2==0)
      width *=generator()->standardModel()->CKM(*mode(imode)->outgoing()[1],
						*mode(imode)->outgoing()[2]);
    else
      width *=generator()->standardModel()->CKM(*mode(imode)->outgoing()[2],
						*mode(imode)->outgoing()[1]);
  }
  // final spin average
  assert(!std::isnan(double(width*MeV)));
  return 0.5*width;
}

Energy6 SMTopDecayer::dGammaIntegrand(Energy2 mffb2, Energy2 mbf2, Energy mt,
				      Energy mb, Energy mf, Energy mfb, Energy mw) const {
  Energy2 mt2(sqr(mt)) ,mb2(sqr(mb)) ,mf2(sqr(mf )),mfb2(sqr(mfb )),mw2(sqr(mw ));
  Energy4 mt4(sqr(mt2)),mb4(sqr(mb2)),mf4(sqr(mf2)),mfb4(sqr(mfb2)),mw4(sqr(mw2));
  return -mbf2 * ( + 6 * mb2 * mf2 * mfb2 * mffb2    +   6 * mb2 * mt2 * mfb2 * mffb2 
		   + 6 * mb2 * mt2 * mf2  * mffb2    +  12 * mb2 * mt2 * mf2 * mfb2 
		   - 3  * mb2 * mfb4  * mffb2        +   3 * mb2 * mf2 * mffb2 * mffb2 
		   - 3  * mb2 * mf4   * mffb2        -   6 * mb2 * mt2 * mfb4 
		   - 6  * mb2 * mt2 * mf4            -   3 * mb4 * mfb2 * mffb2 
		   - 3  * mb4 * mf2 * mffb2          -   6 * mb4 * mf2 * mfb2
		   + 3  * mt4 * mf4                  +   3 * mb4 * mfb4 
		   + 3  * mb4 * mf4                  +   3 * mt4 * mfb4
		   + 3  * mb2 * mfb2 * mffb2 * mffb2 +   3 * mt2 * mfb2 * mffb2 * mffb2 
		   - 3  * mt2 * mfb4 * mffb2         +   3 * mt2 * mf2 * mffb2 * mffb2 
		   - 3  * mt2 * mf4 * mffb2          -   3 * mt4 * mfb2 * mffb2 
		   - 3  * mt4 * mf2 * mffb2          -   6 * mt4 * mf2 * mfb2 
		   + 6  * mt2 * mf2 * mfb2 * mffb2   +  12 * mt2 * mf2 * mw4 
		   + 12 * mb2 * mfb2 * mw4           +  12 * mb2 * mt2 * mw4 
		   + 6  * mw2 * mt2 * mfb2 * mbf2    -  12 * mw2 * mt2 * mf2 * mffb2 
		   - 6  * mw2 * mt2 * mf2 * mbf2     -  12 * mw2 * mt2 * mf2 * mfb2 
		   - 12 * mw2 * mb2  * mfb2 * mffb2  -   6 * mw2 * mb2 * mfb2 * mbf2 
		   + 6  * mw2 * mb2  * mf2 * mbf2    -  12 * mw2 * mb2 * mf2 * mfb2 
		   - 12 * mw2 * mb2 * mt2 * mfb2     -  12 * mw2 * mb2 * mt2 * mf2 
		   + 12 * mf2 * mfb2 * mw4           +   4 * mbf2 * mbf2 * mw4 
		   -  6 * mfb2 * mbf2 * mw4          -   6 * mf2 * mbf2 * mw4 
		   -  6 * mt2 * mbf2 * mw4           -   6 * mb2 * mbf2 * mw4 
		   + 12 * mw2 * mt2 * mf4            +  12 * mw2 * mt4 * mf2 
		   + 12 * mw2 * mb2 * mfb4           +  12 * mw2 * mb4 * mfb2) /mw4 / 3.;
}

void SMTopDecayer::initializeMECorrection(RealEmissionProcessPtr born, double & initial,
					  double & final) {
  if(born->bornOutgoing().size()!=2) return;
  // check the outgoing particles
  PPtr part[2];
  for(unsigned int ix=0;ix<born->bornOutgoing().size();++ix) {
    part[ix]= born->bornOutgoing()[ix];
  }
  // check the final-state particles and get the masses
  if(abs(part[0]->id())==ParticleID::Wplus&&abs(part[1]->id())==ParticleID::b) {
    _ma=part[0]->mass();
    _mc=part[1]->mass();
  }
  else if(abs(part[1]->id())==ParticleID::Wplus&&abs(part[0]->id())==ParticleID::b) {
    _ma=part[1]->mass();
    _mc=part[0]->mass();
  }
  else {
    return;
  }
  // set the top mass
  _mt=born->bornIncoming()[0]->mass();
  // set the gluon mass
  _mg=getParticleData(ParticleID::g)->constituentMass();
  // set the radiation enhancement factors
  initial = _initialenhance;
  final   = _finalenhance;
  // reduced mass parameters
  _a=sqr(_ma/_mt);
  _g=sqr(_mg/_mt);
  _c=sqr(_mc/_mt);
  double lambda = sqrt(1.+sqr(_a)+sqr(_c)-2.*_a-2.*_c-2.*_a*_c);
  _ktb = 0.5*(3.-_a+_c+lambda);
  _ktc = 0.5*(1.-_a+3.*_c+lambda);
  useMe();
}


bool SMTopDecayer::softMatrixElementVeto(PPtr parent,
					 PPtr progenitor,
					 const bool & ,
					 const Energy & highestpT,
					 const vector<tcPDPtr> &,
					 const double & z,
					 const Energy & scale,
					 const Energy & pt) {
  // check if we need to apply the full correction
  // the initial-state correction
  if(abs(progenitor->id())==ParticleID::t&&abs(parent->id())==ParticleID::t) {
    // check if hardest so far
    // if not just need to remove effect of enhancement
    bool veto(false);
    // if not hardest so far
    if(pt<highestpT)
      veto=!UseRandom::rndbool(1./_initialenhance);
    // if hardest so far do calculation
    else {
      // values of kappa and z
      double kappa(sqr(scale/_mt));
      // parameters for the translation
      double w(1.-(1.-z)*(kappa-1.)),u(1.+_a-_c-(1.-z)*kappa),v(sqr(u)-4.*_a*w*z);
      // veto if outside phase space
      if(v<0.) 
	veto=true;
      // otherwise calculate the weight
      else {
	v = sqrt(v);
	double xa((0.5*(u+v)/w+0.5*(u-v)/z)),xg((1.-z)*kappa);
	double f(me(xa,xg)),
	  J(0.5*(u+v)/sqr(w)-0.5*(u-v)/sqr(z)+_a*sqr(w-z)/(v*w*z));
	double wgt(f*J*2./kappa/(1.+sqr(z)-2.*z/kappa)/_initialenhance);
	// This next `if' prevents the hardest emission from the 
	// top shower ever entering the so-called T2 region of the
	// phase space if that region is to be populated by the hard MEC.
	if(useMEforT2()&&xg>xgbcut(_ktb)) wgt = 0.;
	if(wgt>1.) {
	  generator()->log() << "Violation of maximum for initial-state "
			     << " soft veto in "
			     << "SMTopDecayer::softMatrixElementVeto"
			     << "xg = " << xg << " xa = " << xa 
			     << "weight =  " << wgt << "\n";
	  wgt=1.;
	}
	// compute veto from weight
	veto = !UseRandom::rndbool(wgt);
      }
    }
    // return the veto
    return veto;
  }
  // final-state correction
  else if(abs(progenitor->id())==ParticleID::b&&abs(parent->id())==ParticleID::b) {
    // check if hardest so far
    // if not just need to remove effect of enhancement
    // if not hardest so far
    if(pt<highestpT) return !UseRandom::rndbool(1./_finalenhance);
    // if hardest so far do calculation
    // values of kappa and z
    double kappa(sqr(scale/_mt));
    // momentum fractions
    double xa(1.+_a-_c-z*(1.-z)*kappa),r(0.5*(1.+_c/(1.+_a-xa))),root(sqr(xa)-4.*_a);
    if(root<0.) {
      generator()->log() << "Imaginary root for final-state veto in "
			 << "SMTopDecayer::softMatrixElementVeto"
			 << "\nz =  " << z  << "\nkappa = " << kappa
			 << "\nxa = " << xa 
			 << "\nroot^2= " << root;
      return true;
    } 
    root=sqrt(root);
    double xg((2.-xa)*(1.-r)-(z-r)*root);
    // xfact (below) is supposed to equal xg/(1-z). 
    double xfact(z*kappa/2./(z*(1.-z)*kappa+_c)*(2.-xa-root)+root);
    // calculate the full result
    double f(me(xa,xg));
    // jacobian
    double J(z*root);
    double wgt(f*J*2.*kappa/(1.+sqr(z)-2.*_c/kappa/z)/sqr(xfact)/_finalenhance);
    if(wgt>1.) {
      generator()->log() << "Violation of maximum for final-state  soft veto in "
			 << "SMTopDecayer::softMatrixElementVeto"
			 << "xg = " << xg << " xa = " << xa 
			 << "weight =  " << wgt << "\n";
      wgt=1.;
    }
    // compute veto from weight and return
    return !UseRandom::rndbool(wgt);
  }
  // otherwise don't veto
  else return !UseRandom::rndbool(1./_finalenhance);
}

double SMTopDecayer::me(double xw,double xg) {
  double prop(1.+_a-_c-xw),xg2(sqr(xg));
  double lambda=sqrt(1.+_a*_a+_c*_c-2.*_a-2.*_c-2.*_a*_c);
  double denom=(1.-2*_a*_a+_a+_c*_a+_c*_c-2.*_c);
  double wgt=-_c*xg2/prop+(1.-_a+_c)*xg-(prop*(1 - xg)+xg2)
    +(0.5*(1.+2.*_a+_c)*sqr(prop-xg)*xg+2.*_a*prop*xg2)/denom;
  return wgt/(lambda*prop);
}

// xgbcut is the point along the xg axis where the upper bound on the 
// top quark (i.e. b) emission phase space goes back on itself in the 
// xa vs xg plane i.e. roughly mid-way along the xg axis in
// the xa vs xg Dalitz plot.
double SMTopDecayer::xgbcut(double kt) { 
  double lambda2 = 1.+_a*_a+_c*_c-2.*_a-2.*_c-2.*_a*_c; 
  double num1    = kt*kt*(1.-_a-_c);
  double num2    = 2.*kt*sqrt(_a*(kt*kt*_c+lambda2*(kt-1.)));
  return (num1-num2)/(kt*kt-4.*_a*(kt-1.));
}

double SMTopDecayer::loME(const Particle & inpart, const ParticleVector & decay) {
  // spinors
  vector<SpinorWaveFunction   > swave;
  vector<SpinorBarWaveFunction> awave;
  vector<VectorWaveFunction> vwave;
  tPPtr Wboson = abs(decay[0]->id())==ParticleID::Wplus ? decay[0] : decay[1];
  tPPtr bquark = abs(decay[0]->id())==ParticleID::Wplus ? decay[1] : decay[0];
  // spinors 
  if(inpart.id()>0) {
    SpinorWaveFunction   ::calculateWaveFunctions(swave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorBarWaveFunction::calculateWaveFunctions(awave,bquark,outgoing);
  }
  else {
    SpinorBarWaveFunction::calculateWaveFunctions(awave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorWaveFunction   ::calculateWaveFunctions(swave,bquark,outgoing);
  }
  // polarization vectors
  VectorWaveFunction::calculateWaveFunctions(vwave,Wboson,outgoing,false);
  Energy2 scale(sqr(inpart.mass()));
  double me=0.;
  if(inpart.id() == ParticleID::t) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel) {
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  Complex diag = FFWVertex_->evaluate(scale,swave[thel],awave[bhel],vwave[whel]);
	  me += norm(diag);
	}
      }
    }
  }
  else if(inpart.id() == ParticleID::tbar) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel){ 
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  Complex diag = FFWVertex_->evaluate(scale,swave[bhel],awave[thel],vwave[whel]);
	  me += norm(diag);
 	}
      }
    }
  }
  return me;
}


double SMTopDecayer::realME(const Particle & inpart, const ParticleVector & decay,
			    ShowerInteraction inter) {
  // vertex for emission from fermions
  AbstractFFVVertexPtr vertex = inter==ShowerInteraction::QCD ? FFGVertex_ : FFPVertex_;
  // spinors
  vector<SpinorWaveFunction   > swave;
  vector<SpinorBarWaveFunction> awave;
  vector<VectorWaveFunction> vwave,gwave;
  tPPtr Wboson = abs(decay[0]->id())==ParticleID::Wplus ? decay[0] : decay[1];
  tPPtr bquark = abs(decay[0]->id())==ParticleID::Wplus ? decay[1] : decay[0];
  // spinors 
  if(inpart.id()>0) {
    SpinorWaveFunction   ::calculateWaveFunctions(swave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorBarWaveFunction::calculateWaveFunctions(awave,bquark,outgoing);
  }
  else {
    SpinorBarWaveFunction::calculateWaveFunctions(awave,const_ptr_cast<tPPtr>(&inpart),
						  incoming);
    SpinorWaveFunction   ::calculateWaveFunctions(swave,bquark,outgoing);
  }
  // polarization vectors
  VectorWaveFunction::calculateWaveFunctions(vwave,Wboson,outgoing,false);
  VectorWaveFunction::calculateWaveFunctions(gwave,decay[2],outgoing,true );
  Energy2 scale(sqr(inpart.mass()));
  double me=0.;
  vector<Complex> diag(3,0.);
  if(inpart.id() == ParticleID::t) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel) {
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  for(unsigned int ghel =0; ghel <3; ghel+=2) {
	    // emission from top
	    SpinorWaveFunction interF = vertex->evaluate(scale,3,inpart.dataPtr(),swave[thel],gwave[ghel]);
	    diag[0] = FFWVertex_->evaluate(scale,interF,awave[bhel],vwave[whel]);
	    // emission from bottom
	    SpinorBarWaveFunction  interB = vertex->evaluate(scale,3,bquark->dataPtr()->CC(),awave[bhel],gwave[ghel]);
	    diag[1] = FFWVertex_->evaluate(scale,swave[thel],interB,vwave[whel]);
	    // emission from W
	    if(inter==ShowerInteraction::QED) {
	      VectorWaveFunction interV = WWWVertex_->evaluate(scale,3,Wboson->dataPtr()->CC(),vwave[whel],gwave[ghel]);
	      diag[1] = FFWVertex_->evaluate(scale,swave[thel],awave[bhel],interV);
	    }
	    Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    me += norm(sum);
	  }
	}
      }
    }
  }
  else if(inpart.id() == ParticleID::tbar) {
    for(unsigned int thel = 0; thel < 2; ++thel) {
      for(unsigned int bhel = 0; bhel < 2; ++bhel){ 
	for(unsigned int whel = 0; whel < 3; ++whel) {
	  for(unsigned int ghel =0; ghel <3; ghel+=2) {
	    // emission from top
	    SpinorBarWaveFunction  interB = vertex->evaluate(scale,3,inpart.dataPtr(),awave[thel],gwave[ghel]);
	    diag[1] = FFWVertex_->evaluate(scale,swave[bhel],interB,vwave[whel]);
	    // emission from bottom
	    SpinorWaveFunction interF = vertex->evaluate(scale,3,bquark->dataPtr()->CC(),swave[bhel],gwave[ghel]);
	    diag[0] = FFWVertex_->evaluate(scale,interF,awave[thel],vwave[whel]);
	    // emission from W
	    if(inter==ShowerInteraction::QED) {
	      VectorWaveFunction interV = WWWVertex_->evaluate(scale,3,Wboson->dataPtr()->CC(),vwave[whel],gwave[ghel]);
	      diag[1] = FFWVertex_->evaluate(scale,swave[bhel],awave[thel],interV);
	    }
	    Complex sum = std::accumulate(diag.begin(),diag.end(),Complex(0.));
	    me += norm(sum);
	  }
	}
      }
    }
  }
  // divide out the coupling
  me /= norm(vertex->norm());
  // return the total
  return me;
}

double SMTopDecayer::matrixElementRatio(const Particle & inpart,
					const ParticleVector & decay2,
					const ParticleVector & decay3,
					MEOption ,
					ShowerInteraction inter) {
  double Nc = standardModel()->Nc();
  double Cf = (sqr(Nc) - 1.) / (2.*Nc);  
  // if(inter==ShowerInteraction::QED) return 0.;
  // double f  = (1. + sqr(e2()) - 2.*sqr(s2()) + s2() + s2()*e2() - 2.*e2());
  // 
  // 
  // double B  = f/s2();
  
  // Energy2 PbPg = decay3[0]->momentum()*decay3[2]->momentum();
  // Energy2 PtPg = inpart.momentum()*decay3[2]->momentum();
  // Energy2 PtPb = inpart.momentum()*decay3[0]->momentum();

  // double R = Cf *((-4.*sqr(mb())*f/s2()) * ((sqr(mb())*e2()/sqr(PbPg)) + 
  // 		  (sqr(mb())/sqr(PtPg)) - 2.*(PtPb/(PtPg*PbPg))) +
  // 		  (16. + 8./s2() + 8.*e2()/s2()) * ((PtPg/PbPg) + (PbPg/PtPg)) -
  // 		  (16./s2()) * (1. + e2()));
  // return R/B*Constants::pi;
  double Bnew = loME(inpart,decay2);
  double Rnew = realME(inpart,decay3,inter);
  double output = Rnew/Bnew*4.*Constants::pi*sqr(inpart.mass())*UnitRemoval::InvE2;
  if(inter==ShowerInteraction::QCD) output *= Cf;
  return output;
}
#line 1 "./SMHiggsFermionsDecayer.cc"
// -*- C++ -*-
//
// SMHiggsFermionsDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsFermionsDecayer class.
//

#include "SMHiggsFermionsDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/ScalarSpinInfo.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "Herwig/Utilities/Maths.h"
#include "Herwig/Shower/RealEmissionProcess.h"
#include "Herwig/Shower/ShowerAlpha.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMHiggsFermionsDecayer::SMHiggsFermionsDecayer() :
  CF_(4./3.), NLO_(false) {
  _maxwgt.resize(9);
  _maxwgt[0]=0.;
  _maxwgt[1]=0;		
  _maxwgt[2]=0;		
  _maxwgt[3]=0.0194397;	
  _maxwgt[4]=0.463542;	
  _maxwgt[5]=0.;		
  _maxwgt[6]=6.7048e-09; 
  _maxwgt[7]=0.00028665; 
  _maxwgt[8]=0.0809643;  
}

void SMHiggsFermionsDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm)
    throw InitException() << "SMHiggsFermionsDecayer needs the StandardModel class"
			  << " to be either the Herwig one or a class inheriting"
			  << " from it";
  _hvertex = hwsm->vertexFFH();
  // make sure they are initialized
  _hvertex->init();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  // set up the decay modes
  unsigned int imode=0;
  for(unsigned int istep=0;istep<11;istep+=10) {
    for(unsigned ix=1;ix<7;++ix) {
      if(istep<10||ix%2!=0) {
	int iy = ix+istep;
	tPDVector out = {getParticleData( iy),
			 getParticleData(-iy)};
	addMode(new_ptr(PhaseSpaceMode(higgs,out,_maxwgt[imode])));
	++imode;
      }
    }
  }
//   Energy quarkMass = getParticleData(ParticleID::b )->mass();
//   Energy higgsMass = getParticleData(ParticleID::h0)->mass();
//   double mu = quarkMass/higgsMass;
//   double beta = sqrt(1.-4.*sqr(mu));
//   double beta2 = sqr(beta);
//   double aS = SM().alphaS(sqr(higgsMass));
//   double L = log((1.+beta)/(1.-beta));
//   cerr << "testing " << beta << " " << mu << "\n";
//   cerr << "testing " << aS << " " << L << "\n";
//   double fact = 
//     6.-0.75*(1.+beta2)/beta2+12.*log(mu)-8.*log(beta)
//     +(5./beta-2.*beta+0.375*sqr(1.-beta2)/beta2/beta)*L
//     +(1.+beta2)/beta*(4.*L*log(0.5*(1.+beta)/beta)
// 		      -2.*log(0.5*(1.+beta))*log(0.5*(1.-beta))
// 		      +8.*Herwig::Math::ReLi2((1.-beta)/(1.+beta))
// 		      -4.*Herwig::Math::ReLi2(0.5*(1.-beta)));
//   cerr << "testing correction " 
//        << 1.+4./3.*aS/Constants::twopi*fact
//        << "\n"; 
//   double real = 4./3.*aS/Constants::twopi*
//     (8.-0.75*(1.+beta2)/beta2+8.*log(mu)-8.*log(beta)
//      +(3./beta+0.375*sqr(1.-beta2)/pow(beta,3))*L
//      +(1.+beta2)/beta*(-0.5*sqr(L)+4.*L*log(0.5*(1.+beta))
// 		       -2.*L*log(beta)-2.*log(0.5*(1.+beta))*log(0.5*(1.-beta))
// 		       +6.*Herwig::Math::ReLi2((1.-beta)/(1.+beta))
// 		       -4.*Herwig::Math::ReLi2(0.5*(1.-beta))
// 		       -2./3.*sqr(Constants::pi)));
//   double virt = 4./3.*aS/Constants::twopi*
//     (-2.+4.*log(mu)+(2./beta-2.*beta)*L
//      +(1.+beta2)/beta*(0.5*sqr(L)-2.*L*log(beta)+2.*sqr(Constants::pi)/3.
// 		       +2.*Herwig::Math::ReLi2((1.-beta)/(1.+beta))));
//   cerr << "testing real " << real << "\n";
//   cerr << "testing virtual " << virt << "\n";
//   cerr << "testing total no mb corr " << 1.+real+virt << "\n";
//   cerr << "testing total    mb corr " << 1.+real+virt +(8./3. - 2.*log(sqr(mu)))*aS/Constants::pi << "\n";
//   InvEnergy2 Gf = 1.166371e-5/GeV2;
//   Gf = sqrt(2.)*4*Constants::pi*SM().alphaEM(sqr(higgsMass))/8./SM().sin2ThetaW()/
//     sqr(getParticleData(ParticleID::Wplus)->mass());
//   cerr << "testing GF " << Gf*GeV2 << "\n";
//   Energy LO = (3./8./Constants::pi)*sqrt(2)*sqr(quarkMass)*Gf*higgsMass*beta*beta*beta;
//   cerr << "testing LO " << LO/GeV << "\n";
//   cerr << "testing quark mass " << quarkMass/GeV << "\n";
//   cerr << "testing gamma " << (1.+real+virt)*LO/MeV << "\n";
}
  
bool SMHiggsFermionsDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  if(parent->id()!=ParticleID::h0||children.size()!=2) return false;
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if(id1==-id2&&(abs(id1)<=6||(abs(id1)>=11&&abs(id1)<=16)))
    return true;
  else
    return false;
}

ParticleVector SMHiggsFermionsDecayer::decay(const Particle & parent,
					     const tPDVector & children) const {
  // id's of the decaying particles
  tPDVector::const_iterator pit(children.begin());
  int id1((**pit).id());
  int imode=-1;
  if(abs(id1)<=6)                     imode = abs(id1)-1;
  else if(abs(id1)>=11&&abs(id1)<=16) imode = (abs(id1)-11)/2+6;
  ParticleVector output(generate(false,false,imode,parent));
  // set up the colour flow
  if(output[0]->hasColour())      output[0]->antiColourNeighbour(output[1]);
  else if(output[1]->hasColour()) output[1]->antiColourNeighbour(output[0]);
  return output;
}


void SMHiggsFermionsDecayer::persistentOutput(PersistentOStream & os) const {
  os << _maxwgt << _hvertex << NLO_;
}

void SMHiggsFermionsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _maxwgt >> _hvertex >> NLO_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHiggsFermionsDecayer,PerturbativeDecayer>
describeHerwigSMHiggsFermionsDecayer("Herwig::SMHiggsFermionsDecayer", "HwPerturbativeHiggsDecay.so");

void SMHiggsFermionsDecayer::Init() {

  static ClassDocumentation<SMHiggsFermionsDecayer> documentation
    ("The SMHiggsFermionsDecayer class implements the decat of the Standard Model"
     " Higgs boson to the Standard Model fermions");

  static ParVector<SMHiggsFermionsDecayer,double> interfaceMaxWeights
    ("MaxWeights",
     "Maximum weights for the various decays",
     &SMHiggsFermionsDecayer::_maxwgt, 9, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  
  static Switch<SMHiggsFermionsDecayer,bool> interfaceNLO
    ("NLO",
     "Whether to return the LO or NLO result",
     &SMHiggsFermionsDecayer::NLO_, false, false, false);
  static SwitchOption interfaceNLOLO
    (interfaceNLO,
     "No",
     "Leading-order result",
     false);
  static SwitchOption interfaceNLONLO
    (interfaceNLO,
     "Yes",
     "NLO result",
     true);

}

void SMHiggsFermionsDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  int iferm(1),ianti(0);
  if(decay[0]->id()>0) swap(iferm,ianti);
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(_wavebar,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(_wave   ,decay[ianti],outgoing,true);
}

// return the matrix element squared
double SMHiggsFermionsDecayer::me2(const int,const Particle & part,
				   const tPDVector & outgoing,
				   const vector<Lorentz5Momentum> & momenta,
				   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  int iferm(1),ianti(0);
  if(outgoing[0]->id()>0) swap(iferm,ianti);
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }
  SpinorBarWaveFunction wbar(momenta[iferm],outgoing[iferm],Helicity::outgoing);
  SpinorWaveFunction    w   (momenta[ianti],outgoing[ianti],Helicity::outgoing);
  _wavebar.resize(2);
  _wave.resize(2);
  for(unsigned int ihel=0;ihel<2;++ihel) {
    wbar.reset(ihel);
    _wavebar[ihel] = wbar;
    w.reset(ihel);
    _wave[ihel] = w;
  }
  Energy2 scale(sqr(part.mass()));
  unsigned int ifm,ia;
  for(ifm=0;ifm<2;++ifm) {
    for(ia=0;ia<2;++ia) {
      if(iferm>ianti)
	(*ME())(0,ia,ifm)=_hvertex->evaluate(scale,_wave[ia],
					  _wavebar[ifm],_swave);
      else
	(*ME())(0,ifm,ia)=_hvertex->evaluate(scale,_wave[ia],
					  _wavebar[ifm],_swave);
    }
  }
  int id = abs(outgoing[0]->id());
  double output=(ME()->contract(_rho)).real()*UnitRemoval::E2/scale;
  if(id <=6) output*=3.;
  // test of the partial width
  //   Ptr<Herwig::StandardModel>::transient_const_pointer 
  //     hwsm=dynamic_ptr_cast<Ptr<Herwig::StandardModel>::transient_const_pointer>(standardModel());
  //   double g2(hwsm->alphaEM(scale)*4.*Constants::pi/hwsm->sin2ThetaW());
  //   Energy mass(hwsm->mass(scale,outgoing[0])),
  //     mw(getParticleData(ParticleID::Wplus)->mass());
  //   double beta(sqrt(1.-4.*decay[0]->mass()*decay[0]->mass()/scale));
  //   cerr << "testing alpha " << hwsm->alphaEM(scale) << "\n";
  //   Energy test(g2*mass*mass*beta*beta*beta*part.mass()/32./Constants::pi/mw/mw);
  //   if(abs(decay[0]->id())<=6){test *=3.;}
  //   cout << "testing the answer " << output << "     " 
  //        << test/GeV
  //        << endl;
  // leading-order result
  if(!NLO_) return output;
  // fermion mass
  Energy particleMass = outgoing[0]->mass();
  // check decay products coloured, otherwise return
  if(!outgoing[0]->coloured()||
     particleMass==ZERO) return output;
  // inital masses, couplings  etc
  // higgs mass
  mHiggs_ = part.mass();
  // strong coupling
  aS_ = SM().alphaS(sqr(mHiggs_));
  // reduced mass
  mu_  = particleMass/mHiggs_;
  mu2_ = sqr(mu_);
  // generate y
  double yminus = 0.; 
  double yplus  = 1.-2.*mu_*(1.-mu_)/(1.-2*mu2_);
  double y = yminus + UseRandom::rnd()*(yplus-yminus);
  //generate z for D31,2
  double v  = sqrt(sqr(2.*mu2_+(1.-2.*mu2_)*(1.-y))-4.*mu2_)/(1.-2.*mu2_)/(1.-y);
  double zplus  = (1.+v)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
  double zminus = (1.-v)*(1.-2.*mu2_)*y/2./(mu2_ +(1.-2.*mu2_)*y);
  double z = zminus + UseRandom::rnd()*(zplus-zminus);
  // map y,z to x1,x2 for both possible emissions
  double x2 = 1. - y*(1.-2.*mu2_);
  double x1 = 1. - z*(x2-2.*mu2_);
  //get the dipoles
  InvEnergy2 D1 = dipoleSubtractionTerm( x1, x2); 
  InvEnergy2 D2 = dipoleSubtractionTerm( x2, x1); 
  InvEnergy2 dipoleSum = abs(D1) + abs(D2);
  //jacobian
  double jac = (1.-y)*(yplus-yminus)*(zplus-zminus);
  //calculate real
  Energy2 realPrefactor = 0.25*sqr(mHiggs_)*sqr(1.-2.*mu2_)
    /sqrt(calculateLambda(1,mu2_,mu2_))/sqr(Constants::twopi);
  InvEnergy2 realEmission = 4.*Constants::pi*aS_*CF_*calculateRealEmission( x1, x2);
  // calculate the virtual
  double virtualTerm = calculateVirtualTerm();
  // running mass correction
  virtualTerm += (8./3. - 2.*log(mu2_))*aS_/Constants::pi;
  //answer = (born + virtual + real)/born * LO
  output *= 1. + virtualTerm + 2.*jac*realPrefactor*(realEmission*abs(D1)/dipoleSum  - D1);
  // return the answer
  return output;
}

void SMHiggsFermionsDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the PerturbativeDecayer base class
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    os << "newdef " << name() << ":MaxWeights " << ix << " "
	   << _maxwgt[ix] << "\n";
  }
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}

void SMHiggsFermionsDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _maxwgt[ix] = mode(ix)->maxWeight();
    }
  }
}


//calculate lambda
double SMHiggsFermionsDecayer::calculateLambda(double x, double y, double z) const{
  return sqr(x)+sqr(y)+sqr(z)-2.*x*y-2.*x*z-2.*y*z;
}

//calculates the dipole subtraction term for x1, D31,2 (Dij,k),
// 2 is the spectator anti-fermion and 3 is the gluon
InvEnergy2 SMHiggsFermionsDecayer::
dipoleSubtractionTerm(double x1, double x2) const{
  InvEnergy2 commonPrefactor = CF_*8.*Constants::pi*aS_/sqr(mHiggs_);
  return commonPrefactor/(1.-x2)*
    (2.*(1.-2.*mu2_)/(2.-x1-x2)- 
     sqrt((1.-4.*mu2_)/(sqr(x2)-4.*mu2_))*
     (x2-2.*mu2_)*(2.+(x1-1.)/(x2-2.*mu2_)+2.*mu2_/(1.-x2))/(1.-2.*mu2_));
}

//return ME for real emission
InvEnergy2 SMHiggsFermionsDecayer::
calculateRealEmission(double x1, double x2) const {
  InvEnergy2 prefactor = 2./sqr(mHiggs_)/(1.-4.*mu2_);
  return prefactor*(2. + (1.-x1)/(1.-x2) + (1.-x2)/(1.-x1) 
                    + 2.*(1.-2.*mu2_)*(1.-4.*mu2_)/(1.-x1)/(1.-x2)
                    - 2.*(1.-4.*mu2_)*(1./(1.-x2)+1./(1.-x1)) 
                    - 2.*mu2_*(1.-4.*mu2_)*(1./sqr(1.-x2)+1./sqr(1.-x1)));
}

double SMHiggsFermionsDecayer::
calculateVirtualTerm() const {
  // logs and prefactors
  double beta = sqrt(1.-4.*mu2_);
  double L = log((1.+beta)/(1.-beta));
  double prefactor = CF_*aS_/Constants::twopi;
  // non-singlet piece
  double nonSingletTerm = calculateNonSingletTerm(beta, L);
  double virtualTerm = 
    -2.+4.*log(mu_)+(2./beta - 2.*beta)*L 
    + (2.-4.*mu2_)/beta*(0.5*sqr(L) - 2.*L*log(beta)
			 + 2.*Herwig::Math::ReLi2((1.-beta)/(1.+beta)) 
			 + 2.*sqr(Constants::pi)/3.);
  double iEpsilonTerm = 
    2.*(3.-sqr(Constants::pi)/2. + 0.5*log(mu2_) - 1.5*log(1.-2.*mu2_)
	-(1.-2.*mu2_)/beta*(0.5*sqr(L)+sqr(Constants::pi)/6.
			    -2.*L*log(1.-2.*mu2_))
	+ nonSingletTerm);
  return prefactor*(virtualTerm+iEpsilonTerm);
}

//non-singlet piece of I(epsilon) insertion operator
double SMHiggsFermionsDecayer::
calculateNonSingletTerm(double beta, double L) const {
  return  1.5*log(1.-2.*mu2_)  
    + (1.-2.*mu2_)/beta*(- 2.*L*log(4.*(1.-2.*mu2_)/sqr(1.+beta))+
			 + 2.*Herwig::Math::ReLi2(sqr((1.-beta)/(1.+beta)))
			 - 2.*Herwig::Math::ReLi2(2.*beta/(1.+beta)) 
			 - sqr(Constants::pi)/6.) 
    + log(1.-mu_) 
    - 2.*log(1.-2.*mu_) 
    - 2.*mu2_/(1.-2.*mu2_)*log(mu_/(1.-mu_))
    - mu_/(1.-mu_)
    + 2.*mu_*(2*mu_-1.)/(1.-2.*mu2_)
    + 0.5*sqr(Constants::pi);
}

double SMHiggsFermionsDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
						  const ParticleVector & decay3, MEOption,
						  ShowerInteraction inter) {
  mHiggs_ = inpart.mass();
  mu_ = decay2[0]->mass()/mHiggs_;
  mu2_ = sqr(mu_);
  double x1 = 2.*decay3[0]->momentum().t()/mHiggs_;
  double x2 = 2.*decay3[1]->momentum().t()/mHiggs_;
  double pre = inter==ShowerInteraction::QCD ? CF_ : sqr(double(decay2[0]->dataPtr()->iCharge())/3.);
  return pre*calculateRealEmission(x1,x2)*4.*Constants::pi*sqr(mHiggs_);
}
#line 1 "./SMHiggsGGHiggsPPDecayer.cc"
// -*- C++ -*-
//
// SMHiggsGGHiggsPPDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsGGHiggsPPDecayer class.
//

#include "SMHiggsGGHiggsPPDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "Herwig/Utilities/HiggsLoopFunctions.h"

using namespace Herwig;
using namespace Herwig::HiggsLoopFunctions;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<SMHiggsGGHiggsPPDecayer,PerturbativeDecayer>
describeHerwigSMHiggsGGHiggsPPDecayer("Herwig::SMHiggsGGHiggsPPDecayer",
				      "HwPerturbativeHiggsDecay.so");

bool SMHiggsGGHiggsPPDecayer::accept(tcPDPtr parent,
				       const tPDVector & children) const {
  int idp = parent->id();
  int id0 = children[0]->id();
  int id1 = children[1]->id();
  if((idp == ParticleID::h0 && 
      id0 == ParticleID::g     && id1 == ParticleID::g) ||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::gamma)||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::Z0 && id1 == ParticleID::gamma)||
     (idp == ParticleID::h0 && 
      id0 == ParticleID::gamma && id1 == ParticleID::Z0)) 
    return true;
  else
    return false;
}

ParticleVector SMHiggsGGHiggsPPDecayer::decay(const Particle & parent,
					      const tPDVector & children) const {
  int imode(2);
  if(children[0]->id() == ParticleID::gamma && 
     children[1]->id() == ParticleID::gamma)
    imode = 1;
  else if(children[0]->id() ==ParticleID::g)
    imode = 0;
  ParticleVector out(generate(true,false,imode,parent));
  //colour flow
  if(children[0]->id() == ParticleID::g &&
     children[1]->id() == ParticleID::g) {
    out[0]->colourNeighbour(out[1]);
    out[0]->antiColourNeighbour(out[1]);
  }
  return out;
}

void SMHiggsGGHiggsPPDecayer::persistentOutput(PersistentOStream & os) const {
  os << _hggvertex << _hppvertex << _hzpvertex  << _h0wgt
     << _minloop << _maxloop << _massopt;
}

void SMHiggsGGHiggsPPDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _hggvertex >> _hppvertex >> _hzpvertex >> _h0wgt
     >> _minloop >> _maxloop >> _massopt;
}

void SMHiggsGGHiggsPPDecayer::Init() {

  static ClassDocumentation<SMHiggsGGHiggsPPDecayer> documentation
    ("This is an implentation of h0->gg or h0->gamma,gamma "
     "decayer using the SMHGGVertex.");
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHGGVertex
    ("SMHGGVertex",
     "Pointer to SMHGGVertex",
     &SMHiggsGGHiggsPPDecayer::_hggvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHPPVertex
    ("SMHPPVertex",
     "Pointer to SMHPPVertex",
     &SMHiggsGGHiggsPPDecayer::_hppvertex, false, false, true, 
     false, false);
  
  static Reference<SMHiggsGGHiggsPPDecayer,AbstractVVSVertex> 
    interfaceSMHZPVertex
    ("SMHZPVertex",
     "Pointer to SMHZPVertex",
     &SMHiggsGGHiggsPPDecayer::_hzpvertex, false, false, true, 
     false, false);
  
  static ParVector<SMHiggsGGHiggsPPDecayer,double> interfaceMaxWeights
    ("MaxWeights",
     "Maximum weights for the various decays",
     &SMHiggsGGHiggsPPDecayer::_h0wgt, 3, 1.0, 0.0, 10.0,
     false, false, Interface::limited);
  static Parameter<SMHiggsGGHiggsPPDecayer,int> interfaceMinimumInLoop
    ("MinimumInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHiggsGGHiggsPPDecayer::_minloop, 6, 4, 6,
     false, false, Interface::limited);

  static Parameter<SMHiggsGGHiggsPPDecayer,int> interfaceMaximumInLoop
    ("MaximumInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHiggsGGHiggsPPDecayer::_maxloop, 6, 4, 6,
     false, false, Interface::limited);

  static Switch<SMHiggsGGHiggsPPDecayer,unsigned int> interfaceMassOption
    ("MassOption",
     "Option for the treatment of the masses in the loop diagrams",
     &SMHiggsGGHiggsPPDecayer::_massopt, 0, false, false);
  static SwitchOption interfaceMassOptionFull
    (interfaceMassOption,
     "Full",
     "Include the full mass dependence",
     0);
  static SwitchOption interfaceMassOptionLarge
    (interfaceMassOption,
     "Large",
     "Use the heavy mass limit",
     1);

}

void SMHiggsGGHiggsPPDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::constructSpinInfo(_vwave[ix],decay[ix],
					  outgoing,true,
					  decay[ix]->id()!=ParticleID::Z0);
}

double SMHiggsGGHiggsPPDecayer::me2(const int,const Particle & part,
				    const tPDVector & outgoing,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }
  for(unsigned int ix=0;ix<2;++ix) {
    VectorWaveFunction::
      calculateWaveFunctions(_vwave[ix],momenta[ix],outgoing[ix],
			     Helicity::outgoing,outgoing[ix]->id()!=ParticleID::Z0);


  }
  //Set up decay matrix
  Energy2 scale(sqr(part.mass()));
  unsigned int v1hel,v2hel;
  AbstractVVSVertexPtr vertex;
  unsigned int vstep1(2),vstep2(2);
  double sym(1.);
  if(outgoing[0]->id() == ParticleID::g &&
     outgoing[1]->id() == ParticleID::g) {
    vertex = _hggvertex;
    sym = 2.;
  }
  else if(outgoing[0]->id() == ParticleID::gamma &&
  	  outgoing[1]->id() == ParticleID::gamma) {
    vertex = _hppvertex;
    sym = 2.;
  }
  else if(outgoing[0]->id() == ParticleID::Z0 &&
  	  outgoing[1]->id() == ParticleID::gamma) {
    vertex = _hzpvertex;
    vstep1 = 1;
  }
  else if(outgoing[1]->id() == ParticleID::Z0 &&
  	  outgoing[0]->id() == ParticleID::gamma) {
    vertex = _hzpvertex;
    vstep2 = 1;
  }
  else
    assert(false);
  // loop over the helicities of the outgoing bosons
  for(v1hel = 0;v1hel < 3;v1hel+=vstep1) {
    for(v2hel = 0;v2hel < 3;v2hel+=vstep2) {
      (*ME())(0,v1hel,v2hel) = vertex->evaluate(scale,_vwave[0][v1hel],
  						_vwave[1][v2hel],_swave);
    }
  }
  //store matrix element
  double output = ME()->contract(_rho).real()*UnitRemoval::E2/scale;
  //colour factor (N^2 - 1)/4
  if(outgoing[0]->id() == ParticleID::g) output *= 8.;
  //symmetric final states
  output /= sym;
  // return the answer
  return output;
}

void SMHiggsGGHiggsPPDecayer::doinit() {
  PerturbativeDecayer::doinit();
  if(_hggvertex) _hggvertex->init();
  else {
    throw InitException() << "SMHiggsGGHiggsPPDecayer::doinit() - " 
  			  << "_hggvertex is null";
  }
  if(_hppvertex) _hppvertex->init();
  else {
    throw InitException() << "SMHiggsGGHiggsPPDecayer::doinit() - " 
  			  << "_hppvertex is null";
  }
  if(_hzpvertex) _hzpvertex->init();
  //set up decay modes
  tPDPtr higgs  = getParticleData(ParticleID::h0);
  tPDPtr gluon  = getParticleData(ParticleID::g);
  tPDPtr photon = getParticleData(ParticleID::gamma);
  tPDPtr Z0     = getParticleData(ParticleID::Z0);
  // glu,glu mode
  addMode(new_ptr(PhaseSpaceMode(higgs,{gluon,gluon},_h0wgt[0])));
  // gamma,gamma mode
  addMode(new_ptr(PhaseSpaceMode(higgs,{photon,photon},_h0wgt[1])));
  // Z0,gamma mode
  addMode(new_ptr(PhaseSpaceMode(higgs,{Z0,photon},_h0wgt[2])));
}

void SMHiggsGGHiggsPPDecayer::doinitrun() {
  _hggvertex->initrun();
  _hppvertex->initrun();
  _hzpvertex->initrun();
  PerturbativeDecayer::doinitrun();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      _h0wgt[ix] = mode(ix)->maxWeight();
    }
  }
}

void SMHiggsGGHiggsPPDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  // parameters for the PerturbativeDecayer base class
  for(unsigned int ix=0;ix<_h0wgt.size();++ix) {
    os << "newdef " << name() << ":MaxWeights " << ix << " "
	   << _h0wgt[ix] << "\n";
  }
  os << "newdef " << name() << ":SMHGGVertex " << _hggvertex->fullName() << "\n";
  os << "newdef " << name() << ":SMHPPVertex " << _hppvertex->fullName() << "\n";
  os << "newdef " << name() << ":SMHZPVertex " << _hzpvertex->fullName() << "\n";
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" 
		<< fullName() << "\";" << endl;
}

double SMHiggsGGHiggsPPDecayer::matrixElementRatio(const Particle & inpart, const ParticleVector & decay2,
						   const ParticleVector & decay3, MEOption,
#ifndef NDEBUG
						   ShowerInteraction inter) {
#else
						   ShowerInteraction ) {
#endif
  assert(inter==ShowerInteraction::QCD);
  // extract partons and LO momentas
  vector<cPDPtr> partons(1,inpart.dataPtr());
  vector<Lorentz5Momentum> lomom(1,inpart.momentum());
  for(unsigned int ix=0;ix<2;++ix) {
    partons.push_back(decay2[ix]->dataPtr());
    lomom.push_back(decay2[ix]->momentum());
  }
  vector<Lorentz5Momentum> realmom(1,inpart.momentum());
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==2) partons.push_back(decay3[ix]->dataPtr());
    realmom.push_back(decay3[ix]->momentum());
  }
  Energy2 scale = sqr(inpart.mass());
  Energy2 lome = loME(inpart.mass());
  double reme = realME(partons,realmom);
  double ratio = reme/lome*scale;
  // // analytic value for mt -> infinity
  // double x1 = 2.*decay3[0]->momentum().t()/inpart.mass();
  // double x2 = 2.*decay3[1]->momentum().t()/inpart.mass();
  // double x3 = 2.*decay3[2]->momentum().t()/inpart.mass();
  // double test = 8.*Constants::pi*3.*(1.+pow(1-x1,4)+pow(1-x2,4)+pow(1-x3,4))
  //   /(1.-x1)/(1.-x2)/(1.-x3);
  // generator()->log() << "TESTING RATIO " << test << " " << ratio << " " << ratio/test << "\n";
  // remember the symmetry factor
  return ratio/3.;
}

double SMHiggsGGHiggsPPDecayer::realME(//const vector<cPDPtr> & partons, 
				       const vector<cPDPtr> &, 
				       const vector<Lorentz5Momentum> & momenta) const {
  // using std::norm;
  // ScalarWaveFunction  hout(momenta[0],partons[0],outgoing);
  // LorentzPolarizationVector g[3][2];
  // // calculate the polarization vectors for the gluons
  // for(unsigned int iw=0;iw<3;++iw) {
  //   VectorWaveFunction gwave(momenta[iw+1],partons[iw+1],outgoing);
  //   for(unsigned int ix=0;ix<2;++ix) {
  //     //if(iw==2) gwave.reset(10);
  //     //else
  //     gwave.reset(2*ix);
  //     g[iw][ix] = gwave.wave();
  //   }
  // }
  Energy2 mh2 = momenta[0].mass2();
  Energy2 s = (momenta[1]+momenta[2]).m2();
  Energy2 t = (momenta[1]+momenta[3]).m2();
  Energy2 u = (momenta[2]+momenta[3]).m2();
  // calculate the loop functions
  Complex A4stu(0.),A2stu(0.),A2tsu(0.),A2ust(0.);
  for(int ix=_minloop;ix<=_maxloop;++ix) {
    // loop functions
    if(_massopt==0) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      A4stu+=A4(s,t,u,mf2);
      A2stu+=A2(s,t,u,mf2);
      A2tsu+=A2(t,s,u,mf2);
      A2ust+=A2(u,s,t,mf2);
    }
    else {
      A4stu=-1./3.;
      A2stu=-sqr(s/mh2)/3.;
      A2tsu=-sqr(t/mh2)/3.;
      A2ust=-sqr(u/mh2)/3.;
    }
  }
  // Complex A3stu=0.5*(A2stu+A2ust+A2tsu-A4stu);
  // // compute the dot products for the matrix element
  // // and polarization vector * momenta
  // Energy2 pdot[3][3];
  // complex<InvEnergy> eps[3][3][2]; 
  // for(unsigned int ig=0;ig<3;++ig) {
  //   for(unsigned int ip=0;ip<3;++ip) {
  //     pdot[ig][ip]=momenta[ig+1]*momenta[ip+1];
  //     for(unsigned int ih=0;ih<2;++ih) {
  // 	if(ig!=ip)
  // 	  eps[ig][ip][ih]=g[ig][ih].dot(momenta[ip+1])/pdot[ig][ip];
  // 	else
  // 	  eps[ig][ip][ih]=ZERO;
  //     }
  //   }
  // }
  // prefactors
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  // Energy3 pre=sqr(mh2)/mw;
  // // compute the matrix element
  // double output(0.);
  // complex<InvEnergy2> wdot[3][3];
  //  for(unsigned int ghel1=0;ghel1<2;++ghel1) {
  //    for(unsigned int ghel2=0;ghel2<2;++ghel2) {
  //      for(unsigned int ghel3=0;ghel3<2;++ghel3) {
  // 	 wdot[0][1]=g[0][ghel1].dot(g[1][ghel2])/pdot[0][1];
  // 	 wdot[0][2]=g[0][ghel1].dot(g[2][ghel3])/pdot[0][2];
  // 	 wdot[1][0]=wdot[0][1];
  // 	 wdot[1][2]=g[1][ghel2].dot(g[2][ghel3])/pdot[1][2];
  // 	 wdot[2][0]=wdot[0][2];
  // 	 wdot[2][1]=wdot[1][2];
  // 	 // last piece
  // 	 Complex diag=pre*A3stu*(eps[0][2][ghel1]*eps[1][0][ghel2]*eps[2][1][ghel3]-
  // 				 eps[0][1][ghel1]*eps[1][2][ghel2]*eps[2][0][ghel3]+
  // 				 (eps[2][0][ghel3]-eps[2][1][ghel3])*wdot[0][1]+
  // 				 (eps[1][2][ghel2]-eps[1][0][ghel2])*wdot[0][2]+
  // 				 (eps[0][1][ghel1]-eps[0][2][ghel1])*wdot[1][2]);
  // 	 // first piece
  // 	 diag+=pre*(+A2stu*(eps[0][1][ghel1]*eps[1][0][ghel2]-wdot[0][1])*
  // 		    (eps[2][0][ghel3]-eps[2][1][ghel3])
  // 		    +A2tsu*(eps[0][2][ghel1]*eps[2][0][ghel3]-wdot[0][2])*
  // 		    (eps[1][2][ghel2]-eps[1][0][ghel2])
  // 		    +A2ust*(eps[1][2][ghel2]*eps[2][1][ghel3]-wdot[1][2])*
  // 		    (eps[0][1][ghel1]-eps[0][2][ghel1]));
  // 	 output+=norm(diag);
  //      }
  //    }
  //  }
  //  // colour factor and what's left of the prefactor
  //  output *= 6.;
   double me=4.*24./s/t/u*pow<4,1>(mh2)/sqr(mw)*
     (norm(A2stu)+norm(A2ust)+norm(A2tsu)+norm(A4stu));
   return me;
}

Energy2 SMHiggsGGHiggsPPDecayer::loME(Energy mh) const {
  Complex loop(0.);
  Energy2 mh2(sqr(mh));
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  for(int ix=_minloop;ix<=_maxloop;++ix) {
    // loop functions
    if(_massopt==0) {
      Energy2 mf2=sqr(getParticleData(ix)->mass());
      loop += A1(mh2,mf2);
    }
    else {
      loop += 2./3.;
    }
  }
  return 1./Constants::pi*sqr(mh2)/sqr(mw)*norm(loop);
}
#line 1 "./SMHiggsWWDecayer.cc"
// -*- C++ -*-
//
// SMHiggsWWDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHiggsWWDecayer class.
//

#include "SMHiggsWWDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/ParticleData.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
typedef Selector<tDMPtr> DecaySelector;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHiggsWWDecayer,PerturbativeDecayer>
describeHerwigSMHiggsWWDecayer("Herwig::SMHiggsWWDecayer", "HwPerturbativeHiggsDecay.so");

void SMHiggsWWDecayer::Init() {

  static ClassDocumentation<SMHiggsWWDecayer> documentation
    ("The SMHiggsWWDecayer class performs the decay of the Standard Model Higgs"
     " boson to W+W- and Z0Z0");

  static ParVector<SMHiggsWWDecayer,double> interfaceWMaximum
    ("WMaximum",
     "The maximum weight for H-> W+W- decays",
     &SMHiggsWWDecayer::_wmax, 2, 1.0, 0.0, 10000.0,
     false, false, Interface::limited);

  static ParVector<SMHiggsWWDecayer,double> interfaceZMaximum
    ("ZMaximum",
     "The maximum weight for H-> Z0Z0 decays",
     &SMHiggsWWDecayer::_zmax, 2, 1.0, 0.0, 10000.0,
     false, false, Interface::limited);
}

SMHiggsWWDecayer::SMHiggsWWDecayer() : _wmax(2,1.00), _zmax(2,1.00)
{}

void SMHiggsWWDecayer::doinit() {
  PerturbativeDecayer::doinit();
  // get the vertices from the Standard Model object
  tcHwSMPtr hwsm=dynamic_ptr_cast<tcHwSMPtr>(standardModel());
  if(!hwsm) 
    throw InitException() << "SMHiggsWWDecayer needs the StandardModel class"
  			  << " to be either the Herwig one or a class inheriting"
  			  << " from it";
  _theFFWVertex = hwsm->vertexFFW();
  _theFFZVertex = hwsm->vertexFFZ();
  _theHVVVertex = hwsm->vertexWWH();
  // get the width generator for the higgs
  tPDPtr higgs = getParticleData(ParticleID::h0);
  // the W+W- decays
  for(unsigned int ix=0;ix<2;++ix) {
    tPDPtr wplus  = getParticleData(ParticleID::Wplus);
    tPDPtr wminus = getParticleData(ParticleID::Wminus);
    DecaySelector wpDecay =  wplus->decaySelector();
    DecaySelector wmDecay = wminus->decaySelector();
    unsigned int imode=0;
    for(DecaySelector::const_iterator wp=wpDecay.begin();wp!=wpDecay.end();++wp) {
      // extract the decay products of W+
      tPDVector prod=(*wp).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      tPDVector out={prod[0],prod[1],tPDPtr(),tPDPtr()};
      for(DecaySelector::const_iterator wm=wmDecay.begin();wm!=wmDecay.end();++wm) {
	// extract the decay products of W-
	tPDVector prod=(*wm).second->orderedProducts();
	if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
	out[2] = prod[0];
	out[3] = prod[1];
	// create the new mode
	PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(higgs,out,_wmax[ix]));
	// create the phase space channel
	PhaseSpaceChannel phase((PhaseSpaceChannel(mode),0,wplus,0,wminus,1,1,1,2,2,3,2,4));
	mode->addChannel(phase);
	if(ix==0) {
	  phase.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	  phase.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	}
	addMode(mode);
	// insert mode into selector
	_ratio.push_back(wp->second->brat()*wm->second->brat());
	if(ix==0) _wdecays.insert (_ratio.back(),imode);
  	++imode;
      }
    }
    // the Z0Z0 decays
    tPDPtr Z0=getParticleData(ParticleID::Z0);
    DecaySelector Z0Decay = Z0->decaySelector();
    for(DecaySelector::const_iterator z1=Z0Decay.begin();z1!=Z0Decay.end();++z1) {
      // extract the decay products of Z0
      tPDVector prod=(*z1).second->orderedProducts();
      if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
      tPDVector out = {prod[0],prod[1],tPDPtr(),tPDPtr()};
      for(DecaySelector::const_iterator z2=Z0Decay.begin();z2!=Z0Decay.end();++z2) {
	// extract the decay products of Z0
	tPDVector prod=(*z2).second->orderedProducts();
	if(prod[0]->id()<prod[1]->id()) swap(prod[0],prod[1]);
	out[2] = prod[0];
	out[3] = prod[1];
	// create the new mode
	PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(higgs,out,_zmax[ix]));
	// create the phase space channel
	PhaseSpaceChannel phase((PhaseSpaceChannel(mode),0,Z0,0,Z0,1,1,1,2,2,3,2,4));
	mode->addChannel(phase);
	if(ix==0) {
	  phase.setJacobian(1,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	  phase.setJacobian(2,PhaseSpaceChannel::PhaseSpaceResonance::Power,0.);
	}
	addMode(mode);
	// insert mode into selector
	_ratio.push_back(z1->second->brat()*z2->second->brat());
	if(ix==0) _zdecays.insert (_ratio.back(),imode);
	++imode;
      }
    }
  }
}

bool SMHiggsWWDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  // if not two decay products return false
  if(children.size()!=2) return false;
  // if not decaying higgs return false
  if(parent->id()!=ParticleID::h0) return false;
  tPDVector::const_iterator pit = children.begin();
  int id1=(**pit).id();
  ++pit;
  int id2=(**pit).id();
  if((id1==-id2&&abs(id1)==ParticleID::Wplus)||
     (id1== id2&&    id1 ==ParticleID::Z0))
    return true;
  else
    return false;
}

void SMHiggsWWDecayer::persistentOutput(PersistentOStream & os) const {
  os << _theFFWVertex << _theFFZVertex << _theHVVVertex 
     << _wdecays << _zdecays << _ratio << _wmax << _zmax;
}

void SMHiggsWWDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _theFFWVertex >> _theFFZVertex >> _theHVVVertex 
     >> _wdecays >> _zdecays >> _ratio >> _wmax >> _zmax;
}

ParticleVector SMHiggsWWDecayer::decay(const Particle & parent,
				       const tPDVector & children) const {
  // select the decay modes of the gauge bosons
  unsigned int imode;
  if(abs(children[0]->id())==ParticleID::Wplus)
    imode=_wdecays.select(UseRandom::rnd());
  else
    imode=_zdecays.select(UseRandom::rnd());
  // use different phase space for low/high mass higgs
  if(parent.mass()>1.8*children[0]->mass()) 
    imode+=_wdecays.size()+_zdecays.size();
  // generate the kinematics
  ParticleVector decay = generate(true,false,imode,parent);
  // set up the colour flows
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->children()[0]->coloured()) {
      decay[ix]->children()[0]->antiColourNeighbour(decay[ix]->children()[1]);
    }
  }
  return decay;
}

void SMHiggsWWDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
  					  incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(_fwave1,decay[0],outgoing,true);
  SpinorWaveFunction   ::
    constructSpinInfo(_awave1,decay[1],outgoing,true);
  SpinorBarWaveFunction::
    constructSpinInfo(_fwave2,decay[2],outgoing,true);
  SpinorWaveFunction   ::
    constructSpinInfo(_awave2,decay[3],outgoing,true);
}

double SMHiggsWWDecayer::me2(const int,const Particle & part,
				   const tPDVector & outgoing,
				   const vector<Lorentz5Momentum> & momenta,
				   MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half,
  					 PDT::Spin1Half,PDT::Spin1Half)));
  // check if Z or W decay
  bool Z0=outgoing[0]->id()==-outgoing[1]->id();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(_rho,const_ptr_cast<tPPtr>(&part),incoming);
    _swave = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(_rho);
  }

  SpinorBarWaveFunction fw1(momenta[0],outgoing[0],Helicity::outgoing);
  SpinorBarWaveFunction fw2(momenta[2],outgoing[2],Helicity::outgoing);
  SpinorWaveFunction    aw1(momenta[1],outgoing[1],Helicity::outgoing);
  SpinorWaveFunction    aw2(momenta[3],outgoing[3],Helicity::outgoing);
  _fwave1.resize(2);
  _fwave2.resize(2);
  _awave1.resize(2);
  _awave2.resize(2);  
  for(unsigned int ihel=0;ihel<2;++ihel) {
    fw1.reset(ihel);   _fwave1[ihel] = fw1;
    fw2.reset(ihel);   _fwave2[ihel] = fw2;
    aw1.reset(ihel);   _awave1[ihel] = aw1;
    aw2.reset(ihel);   _awave2[ihel] = aw2;
  }  
  // get the intermediates and vertex
  tcPDPtr inter[2];
  AbstractFFVVertexPtr vert;
  if(Z0) {
    inter[0]=getParticleData(ParticleID::Z0);
    inter[1]=inter[0];
    vert=_theFFZVertex;
  }
  else {
    inter[0]=getParticleData(ParticleID::Wplus);
    inter[1]=getParticleData(ParticleID::Wminus);
    vert=_theFFWVertex;
  }
  // construct the spinors for the outgoing particles
  Energy2 scale0(sqr(part.mass()));
  Energy2 scale1((momenta[0]+momenta[1]).m2());
  Energy2 scale2((momenta[2]+momenta[3]).m2());
  // for decays to quarks ensure boson is massive enough to
  // put quarks on constituent mass-shell
  if(scale1<sqr(outgoing[0]->constituentMass()+
  		outgoing[1]->constituentMass())) return 0.;
  if(scale2<sqr(outgoing[2]->constituentMass()+
  		outgoing[3]->constituentMass())) return 0.;
  // compute the boson currents
  VectorWaveFunction curr1[2][2],curr2[2][2];
  unsigned int ohel1,ohel2,ohel3,ohel4;
  for(ohel1=0;ohel1<2;++ohel1) {
    for(ohel2=0;ohel2<2;++ohel2) {
      curr1[ohel1][ohel2]=vert->evaluate(scale1,1,inter[0],
  					 _awave1[ohel2],_fwave1[ohel1]);
      curr2[ohel1][ohel2]=vert->evaluate(scale2,1,inter[1],
  					 _awave2[ohel2],_fwave2[ohel1]);
    }
  }
  // compute the matrix element
  for(ohel1=0;ohel1<2;++ohel1) {
    for(ohel2=0;ohel2<2;++ohel2) {
      for(ohel3=0;ohel3<2;++ohel3) {
  	for(ohel4=0;ohel4<2;++ohel4) {
  	  (*ME())(0,ohel1,ohel2,ohel3,ohel4)=
  	    _theHVVVertex->evaluate(scale0,curr1[ohel1][ohel2],
  				    curr2[ohel3][ohel4],_swave);
  	}
      }
    }
  }
  double output=(ME()->contract(_rho)).real()*scale0*UnitRemoval::InvE2;
  // colour factors
  if(outgoing[0]->coloured()) output *= 3.;
  if(outgoing[2]->coloured()) output *= 3.;
  // divide out the gauge boson branching ratios
  output/=_ratio[imode()];
  // if Z0 decays identical particle factor
  if(Z0) output*=0.5;
  // return the answer
  return output;
}

void SMHiggsWWDecayer::dataBaseOutput(ofstream & os,bool header) const {
  if(header) os << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<2;++ix) {
    os << "newdef " << name() << ":WMaximum "    << ix << " " << _wmax[ix]  << "\n";
    os << "newdef " << name() << ":ZMaximum "    << ix << " " << _zmax[ix]  << "\n";
  }
  PerturbativeDecayer::dataBaseOutput(os,false);
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

void SMHiggsWWDecayer::doinitrun() {
  PerturbativeDecayer::doinitrun();
  for(unsigned int ix=0;ix<2;++ix) {
    _zmax[ix]=0.;
    _wmax[ix]=0.;
  }
  unsigned int ntest=_wdecays.size()+_zdecays.size();
  if(initialize()) {
    for(unsigned int ix=0;ix<numberModes();++ix) {
      unsigned int iloc = ix<ntest ? 0 : 1;
      if(mode(ix)->outgoing()[0]->iCharge()==
	 -mode(ix)->outgoing()[1]->iCharge()) {
	_zmax[iloc]=max(mode(ix)->maxWeight(),_zmax[iloc]);
      }
      else {
	_wmax[iloc]=max(mode(ix)->maxWeight(),_wmax[iloc]);
      }
    }
  }
}
