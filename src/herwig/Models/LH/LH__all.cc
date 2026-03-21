#line 1 "./LHModel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHModel class.
//

#include "LHModel.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHModel::LHModel() 
  : _cott(1.), _tantp(1.),
    _v(246.*GeV), _lamratio(1.), _mH(120.*GeV), _vacratio(0.05),
    _f(3.*TeV), _lambda1(0.), _lambda2(0.),
    _s(0.), _c(0.), _sp(0.), _cp(0.)
{}

void LHModel::persistentOutput(PersistentOStream & os) const {
  os << _cott << _tantp << ounit(_v,GeV) << _lamratio
     << ounit(_mH,GeV) << _vacratio << ounit(_f,GeV) 
     << _s0 << _c0 << _sP << _cP << _sPlus << _cPlus
     << _lambda1 << _lambda2 << _s << _c << _sp << _cp
     << WHHVertex_;
}

void LHModel::persistentInput(PersistentIStream & is, int) {
  is >> _cott >> _tantp >> iunit(_v,GeV)  >> _lamratio
     >> iunit(_mH,GeV) >> _vacratio >> iunit(_f,GeV)  
     >> _s0 >> _c0 >> _sP >> _cP >> _sPlus >> _cPlus
     >> _lambda1 >> _lambda2 >> _s >> _c >> _sp >> _cp
     >> WHHVertex_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHModel,StandardModel>
describeThePEGLHModel("Herwig::LHModel", "HwLHModel.so");

void LHModel::Init() {

  static ClassDocumentation<LHModel> documentation
    ("The LHModel class");

  static Parameter<LHModel,double> interfaceCotTheta
    ("CotTheta",
     "The cotangent of the theta mixing angle",
     &LHModel::_cott, 1.0, 0.1, 10.,
     false, false, Interface::limited);

  static Parameter<LHModel,double> interfaceTanThetaPrime
    ("TanThetaPrime",
     "The tangent of the theta' mixing angle",
     &LHModel::_tantp, 1.0, 0.1, 10.0,
     false, false, Interface::limited);

  static Parameter<LHModel,Energy> interfacef
    ("f",
     "The scale of the non-linear sigma-model",
     &LHModel::_f, TeV, 3.*TeV, 0.0*TeV, 100.0*TeV,
     true, false, Interface::limited);

  static Parameter<LHModel,double> interfaceLambdaRatio
    ("LambdaRatio",
     "The ratio lambda_1/lambda_2 of the top Yukawa couplings.",
     &LHModel::_lamratio, 1.0, 0.01, 100.,
     false, false, Interface::limited);

  static Parameter<LHModel,double> interfaceVEVRatio
    ("VEVRatio",
     "The ratio of the vacuum expection values v'/v",
     &LHModel::_vacratio, 0.05, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<LHModel,Energy> interfacemH
    ("mH",
     "The mass of the lightest Higgs",
     &LHModel::_mH, GeV, 120.0*GeV, 100.0*GeV, 1000.0*GeV,
     false, false, Interface::limited);

  static Reference<LHModel,AbstractVSSVertex> interfaceVertexWHH
    ("Vertex/WHH",
     "Pointer to the WHH vertex",
     &LHModel::WHHVertex_, false, false, true, false, false);

}

void LHModel::doinit() {
  // additional vertices
  if(WHHVertex_) addVertex(WHHVertex_);
  // stuff from the base class
  BSMModel::doinit();
  // compute the parameters of the model
  // W and Z masses
  Energy mw(getParticleData(ParticleID::Wplus)->mass());
  Energy mz(getParticleData(ParticleID::Z0)   ->mass());
  // SM couplings
  double e  = sqrt(4.*Constants::pi*alphaEM(sqr(mz)));
  double sw2(sin2ThetaW()),cw2(1.-sin2ThetaW());
  double g  = e/sqrt(sw2);
  double gp = e/sqrt(cw2);
  // vev
  _v = 2.*mw/g;
  // cos and sin of mixing angles
  double theta (atan(1./_cott));
  _c = cos(theta );
  _s = sin(theta );
  double thetap(atan(_tantp  ));
  _cp = cos(thetap);
  _sp = sin(thetap);
  // xH (Eqn A35)
  double xH = 2.5*g*gp*_s*_c*_sp*_cp*(sqr(_c*_sp)+sqr(_s*_cp))/
    (5.*sqr(g*_sp*_cp)-sqr(gp*_s*_c));
  double vf(sqr(_v/_f));
  // masses of the neutral gauge bosons (Eqn 21,22,A37)
  Energy2 MAH2 = sqr(mz)*sw2*(0.2/sqr(_sp*_cp)/vf-1.+0.25*xH*cw2/sqr(_s*_c)/sw2);
  Energy2 MZH2 = sqr(mw)*(1./sqr(_s*_c)/vf-1.-xH*sw2/sqr(_sp*_cp)/cw2);
  // mass of the heavy charged gauge boson (Eqn. 19/A33) 
  Energy2 MWH2 = sqr(mw)*(1./sqr(_s*_c)/vf-1.);
  // top and heavy top yukawas (from Eqns 26,27)
  Energy mt = getParticleData(ParticleID::t)->mass();
  Energy MT = _f/_v*(1.+sqr(_lamratio))/_lamratio*mt;
  _lambda2 = MT/sqrt(1.+sqr(_lamratio))/_f;
  _lambda1 = _lamratio*_lambda2;
  // masses of the Higgs bosons (Eqns 12 and 13)
  double r = 8.*_f/_v*_vacratio;
  double lamh = 2.*sqr(_mH/_v)/(1.-0.25*sqr(r));
  if(lamh<0.) {
    throw Exception() << "Higgs trilinear coupling negative, reduce f or v'\n"
		      << Exception::runerror;
  }
  Energy2 MPhi2 = lamh*sqr(_f);
  // from Eqn A27
  _sP    = 2.*sqrt(2.)*_vacratio;
  _cP    = 1.-4.*sqr(_vacratio);
  _sPlus = 2.*_vacratio;
  _cPlus = 1.-2.*sqr(_vacratio);
  // from Eqn A28
  _s0 = 2.*sqrt(2.)*_vacratio;
  _c0 = 1.-4.*sqr(_vacratio);
  // set the masses of the new particles
  resetMass( 32,sqrt(MAH2));
  resetMass( 33,sqrt(MZH2));
  resetMass( 34,sqrt(MWH2));
  resetMass(-34,sqrt(MWH2));
  resetMass(  8,MT);
  resetMass( -8,MT);
  resetMass( 35,sqrt(MPhi2));
  resetMass( 36,sqrt(MPhi2));
  resetMass( 37,sqrt(MPhi2));
  resetMass(-37,sqrt(MPhi2));
  resetMass( 38,sqrt(MPhi2));
  resetMass(-38,sqrt(MPhi2));
}
#line 1 "./LHFFZVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFZVertex class.
//

#include "LHFFZVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr << _glH << _grH;
}

void LHFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr >> _glH >> _grH;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFZVertex,FFVVertex>
describeHerwigLHFFZVertex("Herwig::LHFFZVertex", "HwLHModel.so");

void LHFFZVertex::Init() {

  static ClassDocumentation<LHFFZVertex> documentation
    ("The LHFFZVertex class implements the couplings of the Z and Z_H in"
     " the Little Higgs model to the fermions, both of the Standard Model"
     " and the additional heavy top.");

}

LHFFZVertex::LHFFZVertex() : _couplast(0.0), _q2last(0.*GeV2) {
  // set order in the couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHFFZVertex::doinit() {
  for(int ib=23;ib<34;ib+=10) {
    // the quarks
    for(int ix=1;ix<7;++ix) {
      addToList(-ix,    ix,    ib);
    }
    addToList( -8,   8,  ib);
    addToList( -6,   8,  ib);
    addToList( -8,   6,  ib);
    // the leptons
    for(int ix=11;ix<17;++ix) {
      addToList(-ix,    ix,   ib);
    }
  }
  FFVVertex::doinit();
  cLHModelPtr model = dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHModel "
				   << " in LHFFZVertex::doinit()"
				   << Exception::runerror;
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double pre =-0.5/sw/cw;
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double sp2(sqr(sp)),cp2(sqr(cp));
  // from Eqn A35
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  double xB(-2.5/sw*sp*cp*(cp2-sp2));
  double yu  = -0.4, ye  =  0.6;
  double vf(model->vev()/model->f());
  double xL(sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2())));
  double vu  = pre*( 0.5-4./3.*sw2-sqr(vf)*(+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+7./15.  -cp2/6.)));
  double vd  = pre*(-0.5+2./3.*sw2-sqr(vf)*(-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+11./15. +cp2/6.))); 
  double ve  = pre*(-0.5+2.*   sw2-sqr(vf)*(-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*ye-9./5. +1.5*cp2))); 
  double vv  = pre*(+0.5          -sqr(vf)*(+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(   ye-4./5. +0.5*cp2)));
  double au  = pre*(-0.5-sqr(vf)*(-0.5*cw*xW*c/s+sw*xB/sp/cp*(0.2-0.5*cp2)));
  double ad  = pre*( 0.5-sqr(vf)*(+0.5*cw*xW*c/s-sw*xB/sp/cp*(0.2-0.5*cp2)));
  double ae  = pre*( 0.5-sqr(vf)*(+0.5*cw*xW*c/s-sw*xB/sp/cp*(0.2-0.5*cp2)));
  double av  = pre*(-0.5-sqr(vf)*(-0.5*cw*xW*c/s+sw*xB/sp/cp*(0.2-0.5*cp2)));
  double vtl = pre*( 0.5-4./3.*sw2-sqr(vf)*(-0.5*sqr(xL)+0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(2.*yu+9./5.-1.5*cp2
							  +(7./15.-2.*cp2/3.)*xL)));
  double atl = pre*(-0.5          -sqr(vf)*(+0.5*sqr(xL)-0.5*cw*xW*c/s
					    +sw*xB/sp/cp*(+0.2-0.5*cp2-0.2*xL)));
  double vth = 2./3.*sw/cw;
  double ath = 0.;
  double vtm = 0.25*xL*vf/cw/sw;
  double atm = -vtm;
  _gl.resize(17);
  _gr.resize(17);
  for(unsigned ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = vd - ad;
    _gl[2*ix ]   = vu - au;
    _gl[2*ix+9 ] = ve - ae;
    _gl[2*ix+10] = vv - av;
    _gr[2*ix-1]  = vd + ad;
    _gr[2*ix ]   = vu + au;
    _gr[2*ix+9 ] = ve + ae;
    _gr[2*ix+10] = vv + av;
  }
  _gl[6] = vtl - atl;
  _gr[6] = vtl + atl;
  _gl[7] = vtm - atm;
  _gr[7] = vtm + atm;
  _gl[8] = vth - ath;
  _gr[8] = vth + ath;
  // heavy Z
  double fact = 0.25*c/s/sw;
  vu  =  fact;
  vd  = -fact;
  ve  = -fact;
  vv  =  fact;
  au  = -fact;
  ad  =  fact;
  ae  =  fact;
  av  = -fact;
  vtl =  fact;
  atl = -fact;
  vth =  0.;
  ath =  0.;
  vtm =  -0.25*xL*vf*c/s/sw;
  atm = -vtm;
  _glH.resize(17);
  _grH.resize(17);
  for(unsigned ix=1;ix<4;++ix) {
    _glH[2*ix-1]  = vd - ad;
    _glH[2*ix ]   = vu - au;
    _glH[2*ix+9 ] = ve - ae;
    _glH[2*ix+10] = vv - av;
    _grH[2*ix-1]  = vd + ad;
    _grH[2*ix ]   = vu + au;
    _grH[2*ix+9 ] = ve + ae;
    _grH[2*ix+10] = vv + av;
  }
  _glH[6] = vtl - atl;
  _grH[6] = vtl + atl;
  _glH[7] = vtm - atm;
  _grH[7] = vtm + atm;
  _glH[8] = vth - ath;
  _grH[8] = vth + ath;
}

void LHFFZVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm = abs(a->id());
  int ianti = abs(b->id());
  assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)|| iferm == 8);
  // Z0
  if(c->id()==ParticleID::Z0) {
    if(ianti==iferm) {
      left (_gl[iferm]);
      right(_gr[iferm]);
    }
    else {
      left (_gl[7]);
      right(_gr[7]);
    }
  }
  else {
    if(ianti==iferm) {
      left (_glH[iferm]);
      right(_grH[iferm]);
    }
    else {
      left (_glH[7]);
      right(_grH[7]);
    }
  }
}
#line 1 "./LHFFPVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFPVertex class.
//

#include "LHFFPVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHFFPVertex::LHFFPVertex() 
  : _couplast(0.), _q2last(-1.*GeV2) {
  // order in strong and em coupling
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge <<  _gl << _gr;
}

void LHFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge >>  _gl >> _gr;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFPVertex,FFVVertex>
describeHerwigLHFFPVertex("Herwig::LHFFPVertex", "HwLHModel.so");

void LHFFPVertex::Init() {

  static ClassDocumentation<LHFFPVertex> documentation
    ("The LHFFPVertex class implements the couplings of the fermions to"
     " the photon and A_H in the Little Higgs model");

}

void LHFFPVertex::doinit() {
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,    ix,    22);
    addToList(-ix,    ix,    32);
  }
  addToList( -8,   8,  22);
  addToList( -8,   8,  32);
  addToList( -6,   8,  32);
  addToList( -8,   6,  32);
  // the leptons
  for(int ix=11;ix<17;++ix) {
    if(ix%2!=0) addToList(-ix,    ix,    22);
    addToList(-ix,    ix,    32);
  }
  FFVVertex::doinit();
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHFFPVertex::doinit()"
			  << Exception::runerror;
  // charges
  _charge.resize(17);  
  for(int ix=1;ix<4;++ix) {
    _charge[2*ix-1]  = model->ed();
    _charge[2*ix ]   = model->eu();
    _charge[2*ix+9 ] = model->ee();
    _charge[2*ix+10] = model->enu();
  }
  _charge[8] =  model->eu();
  // couplings for the heavy photon taken from table IX
  double cw  = sqrt(1.-sin2ThetaW());
  double xL  = sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2()));
  double cp2 = sqr(model->cosThetaPrime());
  double yu  = -0.4;
  double ye  =  0.6;
  // prefactor after removal of -e
  double pre = -0.5/cw/model->cosThetaPrime()/model->sinThetaPrime();
  // down type quarks
  double gvd   = pre*(2.*yu+11./15.+cp2/6);
  double gad   = pre*(-0.2+0.5*cp2);
  // up type quarks
  double gvu   = pre*(2.*yu+17./15.-5./6.*cp2);
  double gau   = pre*( 0.2-0.5*cp2);
  // charged leptons
  double gve   = pre*(2.*ye-9./5.+1.5*cp2);
  double gae   = pre*(-0.2+0.5*cp2);
  // neutrinos
  double gvv   = pre*(-0.2+0.5*cp2);
  double gav   = pre*( 0.2-0.5*cp2);
  // light top
  double gvtll = pre*(2.*yu+17./15.-5./6.*cp2-0.2*xL);
  double gatll = pre*(0.2-0.5*cp2-0.2*xL);
  // mixed top
  double gvtlh = pre*0.2*model->lambda1()*model->lambda2()/
    (sqr(model->lambda1())+sqr(model->lambda2()));
  double gatlh = gvtlh;
  // heavy top
  double gvthh = pre*(2.*yu+14./15.-4./3.*cp2+0.2*xL);
  double gathh = pre*0.2*xL;
  _gl.resize(17);
  _gr.resize(17);
  for(unsigned int ix=1;ix<4;++ix) {
    _gr[2*ix-1]  = gvd+gad;
    _gl[2*ix-1]  = gvd-gad;
    _gr[2*ix ]   = gvu+gau;
    _gl[2*ix ]   = gvu-gau;
    _gr[2*ix+9 ] = gve+gae;
    _gl[2*ix+9 ] = gve-gae;
    _gr[2*ix+10] = gvv+gav;
    _gl[2*ix+10] = gvv-gav;
  }
  // light top
  _gr[6] = gvtll+gatll;
  _gl[6] = gvtll-gatll;
  // mixed top
  _gr[7] = gvtlh+gatlh;
  _gl[7] = gvtlh-gatlh;
  // heavy top
  _gr[8] = gvthh+gathh;
  _gl[8] = gvthh-gathh;
}

// coupling for FFP vertex
void LHFFPVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int iferm=abs(a->id());
  assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)||iferm==8);
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  // photon
  if(c->id()==ParticleID::gamma) {
    left (_charge[iferm]);
    right(_charge[iferm]);
  }
  // heavy photon
  else {
    assert(c->id()==32);
    int ianti = abs(b->id());
    if(ianti==iferm) {
      left (_gl[iferm]);
      right(_gr[iferm]);
    }
    else {
      left (_gl[7]);
      right(_gr[7]);
    }
  }
}
#line 1 "./LHFFGVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFGVertex class.
//

#include "LHFFGVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

// Static variable needed for the type description system in ThePEG.
DescribeNoPIOClass<LHFFGVertex,FFVVertex>
describeHerwigLHFFGVertex("Herwig::LHFFGVertex", "HwLHModel.so");

void LHFFGVertex::Init() {

  static ClassDocumentation<LHFFGVertex> documentation
    ("The LHFFGVertex class implements the coupling of the quarks"
     " to the gluon in the Little Higgs model");

}

// coupling for FFG vertex
void LHFFGVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
			      tcPDPtr a,
#else
			      tcPDPtr ,
#endif
			      tcPDPtr,tcPDPtr) {
  // check allowed
#ifndef NDEBUG
  int iferm=abs(a->id());
#endif
  assert((iferm>=1 && iferm<=6) || iferm==8);
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = -strongCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  left(1.);
  right(1.);
}

void LHFFGVertex::doinit() {
  // PDG codes for the particles
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 21);
  }
  addToList(-8, 8, 21);
  FFVVertex::doinit();
}

LHFFGVertex::LHFFGVertex() : _couplast(0.), _q2last(0.*GeV2) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}
#line 1 "./LHFFWVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFWVertex class.
//

#include "LHFFWVertex.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;

LHFFWVertex::LHFFWVertex() 
  : _ckm(3,vector<Complex>(3,0.0)), _couplast(0.), _q2last(0.*GeV2),
    _corrL(0.),_corrH(0.),_tcorrL(0.),_tcorrH(0.),_tHcorrL(0.), _tHcorrH(0.) {
  // order of vertex in couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LHFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _ckm << _corrL << _corrH << _tcorrL << _tcorrH << _tHcorrL << _tHcorrH;
}

void LHFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _ckm >> _corrL >> _corrH >> _tcorrL >> _tcorrH >> _tHcorrL >> _tHcorrH;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFWVertex,FFVVertex>
describeHerwigLHFFWVertex("Herwig::LHFFWVertex", "HwLHModel.so");

void LHFFWVertex::Init() {

  static ClassDocumentation<LHFFWVertex> documentation
    ("The LHFFWVertex class implements the vertices for"
     " the coupling of the W and heavy W to the Standard Model "
     "fermions and the heavy top quark in the Little Higgs model");

}
  
void LHFFWVertex::doinit() {
  // particles for outgoing W-
  // quarks
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<7;iy+=2) {
      addToList(-ix,      iy,      -24);
      addToList(-ix,      iy,      -34);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix,    ix+1,    -24);
    addToList(-ix,    ix+1,    -34);
  }
  // particles for outgoing W+
  // quarks
  for(int ix=2;ix<7;ix+=2) {
    for(int iy=1;iy<6;iy+=2) {
      addToList(-ix,      iy,      24);
      addToList(-ix,      iy,      34);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix-1,    ix,    24);
    addToList(-ix-1,    ix,    34);
  }
  // couplings to new heavy quark
  addToList(-5,  8,  -24);
  addToList(-5,  8,  -34);
  addToList(-8,  5,  24);
  addToList(-8,  5,  34);
  ThePEG::Helicity::FFVVertex::doinit();
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHFFWVertex::doinit()"
			  << Exception::runerror;
  // cast the CKM object to the HERWIG one
  ThePEG::Ptr<Herwig::StandardCKM>::transient_const_pointer 
    hwCKM=ThePEG::dynamic_ptr_cast< ThePEG::Ptr<Herwig::StandardCKM>::
    transient_const_pointer>(generator()->standardModel()->CKM());
  if(!hwCKM) 
    throw InitException() << "Must have access to the Herwig::StandardCKM object"
			  << "for the CKM matrix in LHFFWVertex::doinit()"
			  << Exception::runerror;
  _ckm = hwCKM->getUnsquaredMatrix(model->families());
  // compute the correction factors
  double s2(sqr(model->sinTheta())),c2(sqr(model->cosTheta()));
  double vf(model->vev()/model->f());
  double xL = sqr(model->lambda1())/(sqr(model->lambda1())+sqr(model->lambda2()));
  // from Table VIII with -sign to agree with our SM conventions
  _corrL   =  1.-0.5*sqr(vf)*c2*(c2-s2);
  _corrH   = -model->cosTheta()/model->sinTheta();
  _tcorrL  =  1.-0.5*sqr(vf)*(c2*(c2-s2)+sqr(xL));
  _tcorrH  = -model->cosTheta()/model->sinTheta();
  _tHcorrL = -vf*xL;
  _tHcorrH =  vf*xL*model->cosTheta()/model->sinTheta();
}

void LHFFWVertex::setCoupling(Energy2 q2, tcPDPtr a, 
				       tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast    = -sqrt(0.5)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  right(0.);
  // the left and right couplings
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  bool heavy(false);
  if(iferm==8) {
    iferm = 6;
    heavy = true;
  }
  if(ianti==8) {
    ianti = 6;
    heavy = true;
  }
  // quarks
  if(iferm>=1 && iferm <=6) {
    int iu,id;
    // up type first
    if(iferm%2==0) {
      iu = iferm/2;
      id = (ianti+1)/2;
    }
    // down type first
    else {
      iu = ianti/2;
      id = (iferm+1)/2;
    }
    assert( iu>=1 && iu<=3 && id>=1 && id<=3);
    left(_ckm[iu-1][id-1]);
  }
  // leptons
  else if(iferm>=11 && iferm <=16) {
    left(1.);
  }
  else assert(false);
  // correction factors
  // light W
  if(abs(c->id())==ParticleID::Wplus) {
    // light quarks or leptons
    if(iferm<6&&ianti<6) {
      left(_corrL*left());
    }
    // light top quark
    else if(!heavy) {
      left(_tcorrL*left());
    }
    // heavy top quark
    else {
      left(_tHcorrL*left());
    }
  }
  // heavy W
  else {
    // light quarks or leptons
    if(iferm<6&&ianti<6) {
      left(_corrH*left());
    }
    // light top quark
    else if(!heavy) {
      left(_tcorrH*left());
    }
    // heavy top quark
    else {
      left(_tHcorrH*left());
    }
  }
}
#line 1 "./LHWWWVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWWVertex class.
//

#include "LHWWWVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _corr; 
}

void LHWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _corr; 
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWWVertex,VVVVertex>
describeHerwigLHWWWVertex("Herwig::LHWWWVertex", "HwLHModel.so");

void LHWWWVertex::Init() {

  static ClassDocumentation<LHWWWVertex> documentation
    ("The LHWWWVertex class implements the triple electroweak"
     " gauge boson couplings in the Little Higgs model.");

}

LHWWWVertex::LHWWWVertex() : _couplast(0.),_q2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHWWWVertex::doinit() {
  // particles
  addToList(24,  -24,  22);
  addToList(24,  -24,  23);
  addToList(24,  -24,  32);
  addToList(24,  -24,  33);

  addToList(34,  -24,  23);
  addToList(34,  -24,  32);
  addToList(34,  -24,  33);

  addToList(24,  -34,  23);
  addToList(24,  -34,  32);
  addToList(24,  -34,  33);

  addToList(34,  -34,  22);
  addToList(34,  -34,  23);
  addToList(34,  -34,  32);
  addToList(34,  -34,  33);
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWVertex::doinit()"
			  << Exception::runerror;
  // correction factors for the different interactions
  double sw(sqrt(model->sin2ThetaW())),cw(sqrt(1.-model->sin2ThetaW()));
  double vf(sqr(model->vev()/model->f()));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double xB(-2.5/sw*sp*cp*(sqr(cp)-sqr(sp)));
  double xH(2.5/sw/cw*s*c*sp*cp*(sqr(c*sp)+sqr(s*cp))/
	    (5.*sqr(sp*cp/sw)-sqr(s*c/cw)));
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  _corr.resize(12);
  // W_L W_L A_L
  _corr[ 0] = -1.;
  // W_L W_L A_H
  _corr[ 1] = cw/sw*vf*xB;
  // W_L W_H A_L
  _corr[ 2] = 0.;
  // W_L W_H A_H
  _corr[ 3] = -vf/sw*xH;
  // W_H W_H A_L
  _corr[ 4] = -1.;
  // W_H W_H A_H
  _corr[ 5] = vf/sw*(xH*(sqr(c)-sqr(s))/s/c+cw*xB);
  // W_L W_L Z_L
  _corr[ 6] = -cw/sw;
  // W_L W_L Z_H
  _corr[ 7] = vf/sw*(cw*xW+s*c*(sqr(c)-sqr(s)));
  // W_L W_H Z_L
  _corr[ 8] = -vf/sw*xW;
  // W_L W_H Z_H
  _corr[ 9] = -1./sw;
  // W_H W_H Z_L
  _corr[10] = -cw/sw;
  // W_H W_H Z_H
  _corr[11] = (sqr(c)-sqr(s))/s/c/sw;
  VVVVertex::doinit();
}

// couplings for the WWW vertex
void LHWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  int ia(a->iCharge()/3),ib(b->iCharge()/3),ic(c->iCharge()/3);
  int ida(a->id()), idb(b->id()), idc(c->id());
  // get the particles in the interaction
  int ineut,nh(0);
  if(ia==0) {
    ineut=ida;
    if(abs(idb)==34) ++nh;
    if(abs(idc)==34) ++nh;
  }
  else if(ib==0) {
    ineut=idb;
    if(abs(ida)==34) ++nh;
    if(abs(idc)==34) ++nh;
  }
  else {
    ineut=idc;
    if(abs(ida)==34) ++nh;
    if(abs(idb)==34) ++nh;
  }
  if(nh==0) {
    if     (ineut==22) norm(_corr[ 0]*_couplast);
    else if(ineut==23) norm(_corr[ 6]*_couplast);
    else if(ineut==32) norm(_corr[ 1]*_couplast);
    else if(ineut==33) norm(_corr[ 7]*_couplast);
  }
  else if(nh==1) {
    if     (ineut==22) norm(_corr[ 2]*_couplast);
    else if(ineut==23) norm(_corr[ 8]*_couplast);
    else if(ineut==32) norm(_corr[ 3]*_couplast);
    else if(ineut==33) norm(_corr[ 9]*_couplast);
  }
  else if(nh==2) {
    if     (ineut==22) norm(_corr[ 4]*_couplast);
    else if(ineut==23) norm(_corr[10]*_couplast);
    else if(ineut==32) norm(_corr[ 5]*_couplast);
    else if(ineut==33) norm(_corr[11]*_couplast);
  }
  // check the order for the overall sign 
  if((ia<0  && ib>0  && ic==0) || 
     (ia==0 && ib<0  && ic>0 ) || 
     (ia>0  && ib==0 && ic<0 ) ) norm(-norm());
}
#line 1 "./LHWWWWVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWWWVertex class.
//

#include "LHWWWWVertex.h"
#include "LHModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHWWWWVertex::LHWWWWVertex() : 
  _couplast(0.0), _q2last(sqr(Constants::MaxEnergy)), _coup(36,0.) {
  // order in the couplings
  orderInGem(2);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHWWWWVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHWWWWVertex::fullclone() const {
  return new_ptr(*this);
}

void LHWWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _coup;
}

void LHWWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _coup;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWWWVertex,VVVVVertex>
describeHerwigLHWWWWVertex("Herwig::LHWWWWVertex", "HwLHModel.so");

void LHWWWWVertex::Init() {

  static ClassDocumentation<LHWWWWVertex> documentation
    ("The LHWWWWVertex class implements the quartic electroweak"
     " boson couplings in the Little Higgs Model");

}

void LHWWWWVertex::doinit() {
  // all charge W's
  addToList(24, -24, 24, -24);
  addToList(34, -34, 34, -34);
  addToList(24, -24, 34, -34);
  addToList(24, -24, 24, -34);
  addToList(24, -24, 34, -24);
  addToList(34, -24, 34, -24);
  addToList(24, -34, 24, -34);
  addToList(34, -34, 24, -34);
  addToList(34, -34, 34, -24);
  // two neutral and 2 W_L
  addToList(22,  24, 22, -24);
  addToList(23,  24, 23, -24);
  addToList(22,  24, 23, -24);
  addToList(22,  24, 32, -24);
  addToList(22,  24, 33, -24);
  addToList(23,  24, 33, -24);
  addToList(23,  24, 32, -24);
  addToList(33,  24, 33, -24);
  addToList(33,  24, 32, -24);
  // two neutral and 2 W_H
  addToList(22,  34, 22, -34);
  addToList(23,  34, 23, -34);
  addToList(22,  34, 23, -34);
  addToList(22,  34, 32, -34);
  addToList(22,  34, 33, -34);
  addToList(23,  34, 33, -34);
  addToList(23,  34, 32, -34);
  addToList(33,  34, 33, -34);
  addToList(33,  34, 32, -34);
  // two neutral W_L W_H
  addToList(23,  24, 23, -34);
  addToList(23,  24, 22, -34);
  addToList(22,  24, 32, -34);
  addToList(23,  24, 32, -34);
  addToList(33,  24, 33, -34);
  addToList(33,  24, 32, -34);
  addToList(22,  24, 33, -34);
  addToList(23,  24, 33, -34);
  addToList(23,  34, 23, -24);
  addToList(23,  34, 22, -24);
  addToList(22,  34, 32, -24);
  addToList(23,  34, 32, -24);
  addToList(33,  34, 33, -24);
  addToList(33,  34, 32, -24);
  addToList(22,  34, 33, -24);
  addToList(23,  34, 33, -24);
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  // correction factors for the different interactions
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double vf(sqr(model->vev()/model->f()));
  double xB(-2.5/sw*sp*cp*(sqr(cp)-sqr(sp)));
  double xW(-0.5/cw*s*c*(sqr(c)-sqr(s)));
  double xH(2.5/sw/cw*s*c*sp*cp*(sqr(c*sp)+sqr(s*cp))/
	    (5.*sqr(sp*cp/sw)-sqr(s*c/cw)));
  // 4 W's
  _coup[ 0] =-1./sw2;
  _coup[ 1] =-1./sw2;
  _coup[ 2] = 0.5/sw2*(sqr(c)-sqr(s))/c/s;
  _coup[ 3] =-0.25/sw2*vf*s*c*(sqr(c)-sqr(s));
  _coup[ 4] =-0.25/sw2;
  _coup[ 5] =-1./sw2*(pow(c,6)+pow(s,c))/sqr(s*c);
  // 2 W_L
  _coup[ 6] = 1.;
  _coup[ 7] = sqr(cw/sw);
  _coup[ 8] = cw/sw;
  _coup[ 9] =-cw/sw*xB*vf;
  _coup[10] =-cw/sw*xW*vf+0.5/sw*s*c*(sqr(c)-sqr(s))*vf;
  _coup[11] =-(sqr(cw)-sw2)/sw2*xW*vf;
  _coup[12] =-sqr(cw/sw)*xB*vf;
  _coup[13] = 0.;
  _coup[14] = 1./sw2;
  _coup[15] = xH*vf/sw2;
  _coup[16] = 1.;
  _coup[17] = sqr(cw/sw);
  _coup[18] = cw/sw;
  _coup[19] =-cw/sw*xB*vf-xH/sw*vf*(sqr(c)-sqr(s))/s/c;
  _coup[20] =-1./sw*(sqr(c)-sqr(s))/s/c;
  _coup[21] =-cw/sw2*(sqr(c)-sqr(s))/s/c;
  _coup[22] =-sqr(cw/sw)*xB*vf-cw/sw2*xH*vf*(sqr(c)-sqr(s))/c/s;
  _coup[23] = 0.;
  _coup[24] = (pow(c,6)+pow(s,6))/sqr(s*c)/sw2;
  _coup[25] = xH/sw2*vf*(pow(c,6)+pow(s,6))/sqr(s*c)
             +cw/sw2*xB*vf*(sqr(c)-sqr(s))/s/c;
  _coup[26] = 0.;
  _coup[27] = 2.*cw/sw2*xW*vf;
  _coup[28] = 0.;
  _coup[29] = xH*vf/sw;
  _coup[30] = xH*vf*cw/sw2;
  _coup[31] = xW*vf/sw;
  _coup[32] =-(sqr(c)-sqr(s))/s/c/sw2;
  _coup[33] =-xH*vf*(sqr(c)-sqr(s))/s/c-cw/sw2*xB*vf;
  _coup[34] = 1./sw;
  _coup[35] = cw/sw2;
  VVVVVertex::doinit();
}

void LHWWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
			       tcPDPtr c,tcPDPtr d) {
  // id's of the particles
  long int id[4]={a->id(),b->id(),c->id(),d->id()};
  // order the particles
  int ngamma(0),nz(0);
  int iorder[4];
  for(int ix=0;ix<4;++ix) {
    if      (id[ix]==22||id[ix]==32) ++ngamma;
    else if (id[ix]==23||id[ix]==33) ++nz;
  }
  // if photons or Z's
  if(ngamma!=0 || nz!=0) {
    int iy=0;
    // put the photons first
    for(int ix=0;iy<ngamma&&ix<4;++ix) {
      if(id[ix]==22||id[ix]==32) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the Z bosons
    for(int ix=0;iy<ngamma+nz&&ix<4;++ix) {
      if(id[ix]==23||id[ix]==33) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24||id[ix]==34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==3);
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24||id[ix]==-34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==4);
  }
  else {
    int iy=0;
    // first the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24||id[ix]==34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==2);
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24||id[ix]==-34) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==4);
  }
  setOrder(iorder[0],iorder[1],iorder[2],iorder[3]);
  setType(2);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = sqr(electroMagneticCoupling(q2));
    _q2last=q2;
  }
  // ids of the particles
  for(unsigned int ix=0;ix<4;++ix) {
    if     (iorder[ix]==0) id[ix] = abs(a->id());
    else if(iorder[ix]==1) id[ix] = abs(b->id());
    else if(iorder[ix]==2) id[ix] = abs(c->id());
    else if(iorder[ix]==3) id[ix] = abs(d->id());
  }
  if( ngamma == 0 && nz == 0 ) {
    if(id[0]==id[1]) {
      if(id[2]==id[3]) {
	if(id[0]==24&&id[2]==24)
	  norm(_couplast*_coup[0]);
	else if(id[0]==34&&id[2]==34)
	  norm(_couplast*_coup[5]);
	else
	  norm(_couplast*_coup[1]);
      }
      else {
	if(id[0]==24)
	  norm(_couplast*_coup[3]);
	else
	  norm(_couplast*_coup[2]);
      }
    }
    else {
      if(id[2]==id[3]) {
	if(id[2]==24)
	  norm(_couplast*_coup[3]);
	else
	  norm(_couplast*_coup[2]);
      }
      else
	norm(_couplast*_coup[4]);
    } 
  }
  else {
    if(id[2]==id[3]) {
      unsigned int ioff = id[2]==24 ? 0 : 10;
      if(id[0]==22&&id[1]==22)
	norm(_couplast*_coup[6+ioff]);
      else if(id[0]==23&&id[1]==23)
	norm(_couplast*_coup[7+ioff]);
      else if((id[0]==22&&id[1]==23) || (id[0]==23&&id[1]==22))
	norm(_couplast*_coup[8+ioff]);
      else if((id[0]==22&&id[1]==32) || (id[0]==32&&id[1]==22))
	norm(_couplast*_coup[9+ioff]);
      else if((id[0]==22&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[10+ioff]);
      else if((id[0]==23&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[11+ioff]);
      else if((id[0]==23&&id[1]==32) || (id[0]==32&&id[1]==23))
	norm(_couplast*_coup[12+ioff]);
      else if( id[0]==33&&id[1]==33)
	norm(_couplast*_coup[14+ioff]);
      else if((id[0]==32&&id[1]==33) || (id[0]==33&&id[1]==32))
	norm(_couplast*_coup[15+ioff]);
      else
	assert(false);
    }
    else { 
      if(id[0]==23&&id[1]==23)
	norm(_couplast*_coup[27]);
      else if((id[0]==22&&id[1]==23) || (id[0]==23&&id[1]==22))
	norm(_couplast*_coup[28]);
      else if((id[0]==22&&id[1]==32) || (id[0]==32&&id[1]==22))
	norm(_couplast*_coup[29]);
      else if((id[0]==23&&id[1]==32) || (id[0]==32&&id[1]==23))
	norm(_couplast*_coup[30]);
      else if( id[0]==33&&id[1]==33)
	norm(_couplast*_coup[31]);
      else if((id[0]==32&&id[1]==33) || (id[0]==33&&id[1]==32))
	norm(_couplast*_coup[32]);
      else if((id[0]==22&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[33]);
      else if((id[0]==23&&id[1]==33) || (id[0]==33&&id[1]==22))
	norm(_couplast*_coup[34]);
      else
	assert(false);
    }
  }
}
#line 1 "./LHFFHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHFFHVertex class.
//

#include "LHFFHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coup,1./GeV) << _model;
}

void LHFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coup,1./GeV) >> _model;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHFFHVertex,FFSVertex>
describeHerwigLHFFHVertex("Herwig::LHFFHVertex", "HwLHModel.so");

void LHFFHVertex::Init() {

  static ClassDocumentation<LHFFHVertex> documentation
    ("The LHFFHVertex class implements the interaction of the fermions"
     " and the Higgs bosons in the Little Higgs model");

}

LHFFHVertex::LHFFHVertex() 
  : _q2last(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
  _masslast[0] = 0.*GeV; 
  _masslast[1] = 0.*GeV;
  _idlast[0] = 0;
  _idlast[1] = 0;
  colourStructure(ColourStructure::DELTA);
}

void LHFFHVertex::doinit() {
  // SM like higgs
  for (int ix=1;ix<=6;++ix) {
    addToList(  -ix,  ix,  25);
  }
  addToList(  -6,   8,  25);
  addToList(  -8,   6,  25);
  addToList(  -8,   8,  25);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix,  ix,  25);
  }
  // phi0
  for (int ix=1;ix<=6;++ix) {
    addToList(  -ix,  ix,  35);
  }
  addToList(  -6,   8,  35);
  addToList(  -8,   6,  35);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix,  ix,  35);
  }
  // phiP
  for (int ix=1;ix<=6;++ix) {
    addToList(  -ix,  ix,  36);
  }
  addToList(  -6,   8,  36);
  addToList(  -8,   6,  36);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix,  ix,  36);
  }
  // phi +/-
  for(int ix=1;ix<6;ix+=2) {
    addToList( -ix-1,   ix,  37);
    addToList( -ix  , ix+1, -37);
  }
  addToList( -8 ,   5,  37);
  addToList( -5 ,   8, -37);
  for(int ix=11;ix<16;ix+=2) {
    addToList( -ix-1,   ix,  37);
    addToList( -ix  , ix+1, -37);
  }
  _model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!_model)   throw InitException() << "Must be using the LHModel "
				      << " in LHFFPVertex::doinit()"
				      << Exception::runerror;
  _coup.resize(11);
  Energy v     = _model->vev();
  double s0    = _model->sinTheta0();
  double sP    = _model->sinThetaP();
  double sPlus = _model->sinThetaPlus();
  double s02   = sqr(s0);
  double vf    = _model->vev()/_model->f();
  double xL    = sqr(_model->lambda1())/(sqr(_model->lambda1())+sqr(_model->lambda2()));
  double xR    = sqr(_model->lambda1())/sqrt(sqr(_model->lambda1())+sqr(_model->lambda2()));
  Energy mT    = getParticleData(8)->mass();
  // lightest higgs couplings
  // coupling of light SM fermions
  _coup[0] = (1.-0.5*s02+vf*s0/sqrt(2.)-2./3.*sqr(vf))/v;
  // couplings to top quark
  _coup[1] = (1.-0.5*s02+vf*s0/sqrt(2.)-2./3.*sqr(vf)+sqr(vf)*xL*(1.+xL))/v;
  // couplings to the T quark
  _coup[2] =-xR*(1.+xL)*vf/mT;
  // couplings to tT
  _coup[3] = xR/mT;
  _coup[4] = vf/v*(1.+xL);
  // phi 0
  // light particles
  _coup[5] = sqrt(0.5)/v*(vf-sqrt(2.)*s0);
  // mixed
  _coup[6] = sqrt(0.5)/v*(vf-sqrt(2.)*s0)*_model->lambda1()/_model->lambda2();
  // phi P
  _coup[7] = Complex(0.,1.)*sqrt(0.5)/v*(vf-sqrt(2.)*sP);
  _coup[8] = Complex(0.,1.)*sqrt(0.5)/v*(vf-sqrt(2.)*sP)*_model->lambda1()/_model->lambda2();
  // phi +/-
  _coup[9] = -sqrt(0.5)/v*(vf-2.*sPlus);
  _coup[9] = -sqrt(0.5)/v*(vf-2.*sPlus)*_model->lambda1()/_model->lambda2();
  FFSVertex::doinit();
}

void LHFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int iferm=abs(a->id());
  int ianti=abs(b->id());
  // left and right couplings set to one
  left (1.);
  right(1.);
  // first the overall normalisation
  if(q2!=_q2last||_idlast[0]!=iferm||_idlast[1]!=ianti) {
    _q2last = q2;
    _idlast[0] = iferm;
    if(_idlast[0]==8) _idlast[0]=6;
    assert((_idlast[0]>=1  && _idlast[0]<=6 ) || 
	   (_idlast[0]>=11 && _idlast[0]<=16));
    if(iferm==_idlast[0])
      _masslast[0] = _model->mass(q2,a);
    else
      _masslast[0] = _model->mass(q2,getParticleData(ParticleID::t));
    _idlast[1] = ianti;
    if(_idlast[1]==8) _idlast[1]=6;
    assert((_idlast[1]>=1  && _idlast[1]<=6 ) || 
	   (_idlast[1]>=11 && _idlast[1]<=16));
    if(_idlast[0]!=_idlast[1]) {
      if(ianti==_idlast[1])
	_masslast[1] = _model->mass(q2,a);
      else
	_masslast[1] = _model->mass(q2,getParticleData(ParticleID::t));
    }
    else {
      _masslast[1] = _masslast[0];
    }
  }
  // SM like higgs
  if(c->id()==ParticleID::h0) {
    if(iferm==ianti) {
      if((iferm>=1  && iferm<=5 ) || 
	 (iferm>=11 && iferm<=16)) {
	norm(-Complex(_coup[0]*_masslast[0]));
      }
      else if(iferm==6) {
	norm(-Complex(_coup[1]*_masslast[0]));
      }
      else if(iferm==8) {
	norm(-Complex(_coup[2]*a->mass()));
      }
      else assert(false);
    }
    else {
      assert( (iferm == 6 && ianti == 8 ) ||
	      (ianti == 6 && iferm == 8 ));
      Complex cleft,cright;
      if(iferm==6) {
	cleft  = Complex(-_coup[3]*b->mass());
	cright = Complex(-_coup[4]*_masslast[0]);
      }
      else {
	cleft  = Complex(-_coup[3]*a->mass());
	cright = Complex(-_coup[4]*_masslast[0]);
      }
      if(b->id()==ParticleID::tbar || c->id()==ParticleID::tbar) {
	cright = conj(cleft);
	cleft = 0.;
      }
      left (cleft );
      right(cright);
      norm(1.);
    }
  }
  else if(c->id()==ParticleID::H0) {
    if(iferm==ianti) {
      if((iferm>=1  && iferm<=6 ) || 
	 (iferm>=11 && iferm<=16)) {
	norm(-Complex(_coup[5]*_masslast[0]));
      }
      else assert(false);
    }
    else {
      assert( (iferm == 6 && ianti == 8 ) ||
	      (iferm == 8 && ianti == 6 ) );
      Complex cleft  = Complex(_coup[6]*_masslast[0]);
      Complex cright = 0.;
      if(b->id()==ParticleID::tbar || c->id()==ParticleID::tbar) {
	cright = conj(cleft);
	cleft = 0.;
      }
      left (cleft );
      right(cright);
      norm(1.);
    }
  }
  else if(c->id()==ParticleID::A0) {
    left(-1.);
    right(1.);
    if(iferm==ianti) {
      if((iferm>=1  && iferm<=6 ) || 
	 (iferm>=11 && iferm<=16)) {
	if(iferm%2==0)
	  norm(-Complex( _coup[7]*_masslast[0]));
	else
	  norm(-Complex(-_coup[7]*_masslast[0]));
      }
      else assert(false);
    }
    else {
      assert( (iferm == 6 && ianti == 8 ) ||
	      (iferm == 8 && ianti == 6 ));
      Complex cleft  = Complex(_coup[8]*_masslast[0]);
      Complex cright = 0.;
      if(b->id()==ParticleID::tbar || c->id()==ParticleID::tbar) {
	cright = conj(cleft);
	cleft = 0.;
      }
      left (cleft );
      right(cright);
      norm(1.);
    }
  }
  else if(c->id()==ParticleID::Hplus) {
    norm(1.);
    Complex cleft(0.),cright(0.);
    if(iferm%2==0) {
      if(iferm==ParticleID::t) {
	cleft  = Complex(_masslast[0]*_coup[ 9]);
      }
      else {
	cleft  = Complex(_masslast[0]*_coup[10]);
	cright = Complex(_masslast[1]*_coup[10]);
      }
    }
    else {
      if(ianti==ParticleID::t) {
	cleft  = Complex(_masslast[1]*_coup[ 9]);
      }
      else {
	cleft  = Complex(_masslast[1]*_coup[10]);
	cright = Complex(_masslast[0]*_coup[10]);
      }
    }
    left ( cleft);
    right(cright);
  }
  else if(c->id()==ParticleID::Hminus) {
    norm(1.);
    Complex cleft(0.),cright(0.);
    if(iferm%2==0) {
      if(iferm==ParticleID::t) {
	cright = Complex(_masslast[0]*_coup[ 9]);
      }
      else {
	cright = Complex(_masslast[0]*_coup[10]);
	cleft  = Complex(_masslast[1]*_coup[10]);
      }
    }
    else {
      if(ianti==ParticleID::t) {
	cright = Complex(_masslast[1]*_coup[ 9]);
      }
      else {
	cright = Complex(_masslast[1]*_coup[10]);
	cleft  = Complex(_masslast[0]*_coup[10]);
      }
    }
    left ( cleft);
    right(cright);
  }
}
#line 1 "./LHWWHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWHVertex class.
//

#include "LHWWHVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LHWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_coup,GeV);
}

void LHWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_coup,GeV);
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWHVertex,VVSVertex>
describeHerwigLHWWHVertex("Herwig::LHWWHVertex", "HwLHModel.so");

void LHWWHVertex::Init() {

  static ClassDocumentation<LHWWHVertex> documentation
    ("The LHWWHVertex class implements the coupling of two electroweak"
     " gauge bosons to a Higgs boson in the Little Higgs Model including the "
     "additional heavy photon, Z and W bosons and the triplet Higgs bosons.");

}

LHWWHVertex::LHWWHVertex() 
  : _couplast(0.), _q2last(0.*GeV2) {
  // order in the couplings
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void LHWWHVertex::doinit() {
  // W_L W_L H
  addToList(  24,  -24,    25);
  // Z_L Z_L H
  addToList(  23,   23,    25);
  // W_L W_H H
  addToList(  24,  -34,    25);
  addToList(  34,  -24,    25);
  // Z_L A_H H
  addToList(  23,   32,    25);
  // W_H W_H H
  addToList(  34,  -34,    25);
  // Z_H Z_H H
  addToList(  33,   33,    25);
  // A_H A_H H
  addToList(  32,   32,    25);
  // Z_H Z_L H
  addToList(  23,   33,    25);
  // Z_H A_H H
  addToList(  33,   32,    25);
  // W_L W_L Phi0
  addToList(  24,  -24,    35);
  // W_L W_H Phi0
  addToList(  24,  -34,    35);
  addToList(  34,  -24,    35);
  // Z_L Z_L Phi0
  addToList(  23,   23,    35);
  // Z_L Z_H Phi0
  addToList(  23,   33,    35);
  // W_H W_H Phi0
  addToList(  34,  -34,    35);
  // Z_H Z_H Phi0
  addToList(  33,   33,    35);
  // A_H Z_H Phi0
  addToList(  32,   33,    35);
  // A_H Z_L Phi0
  addToList(  32,   23,    35);
  // A_H A_H Phi0
  addToList(  32,   32,    35);
  // W_L Z_L Phi-
  addToList(  24,   23,   -37);
  addToList( -24,   23,    37);
  // W_L A_H Phi-
  addToList(  24,   32,   -37);
  addToList( -24,   32,    37);
  // W_L Z_H Phi-
  addToList(  24,   33,   -37);
  addToList( -24,   33,    37);
  // W_H Z_L Phi-
  addToList(  34,   23,   -37);
  addToList( -34,   23,    37);
  // W_H A_H Phi-
  addToList(  34,   32,   -37);
  addToList( -34,   32,    37);
  // W_H Z_H Phi-
  addToList(  34,   33,   -37);
  addToList( -34,   33,    37);
  // W_L W_L Phi--
  addToList(  24,   24,   -38);
  addToList( -24,  -24,    38);
  // W_H W_H Phi--
  addToList(  34,   34,   -38);
  addToList( -34,  -34,    38);
  // W_L W_H Phi--
  addToList(  24,   34,   -38);
  addToList( -24,  -34,    38);
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) throw InitException() << "Must be using the LHModel "
				   << " in LHWWHVertex::doinit()"
				   << Exception::runerror;
  // base class
  VVSVertex::doinit();
  // calculate the couplings for the different combinations of particles
  double sw(sqrt(sin2ThetaW())),cw(sqrt(1.-sin2ThetaW()));
  Energy fact = getParticleData(ParticleID::Wplus)->mass()/sw;
  double vf(sqr(model->vev()/model->f()));
  double vr(model->vevPrime()/model->vev());
  double r2(sqrt(2.));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double sPlus(model->sinTheta0());
  double s0(model->sinTheta0());
  _coup.resize(27);
  // couplings to SM higgs
  _coup[ 0] = fact        *(1.-vf/3.+0.5*vf* sqr(sqr(c)-sqr(s))
			    -0.5*sqr(s0)-2.*r2*s0*vr);
  _coup[ 1] = fact/sqr(cw)*(1.-vf/3.-0.5*vf*(sqr(sqr(c)-sqr(s))+5.*sqr(sqr(cp)-sqr(sp)))
			    -0.5*sqr(s0)+4.*r2*s0*vr);
  _coup[ 2] =-fact;
  _coup[ 3] =-fact;
  _coup[ 4] =-fact*sqr(sw/cw);
  _coup[ 5] =-fact*   0.5*(sqr(c)-sqr(s))/s/c;
  _coup[ 6] =-fact/cw*0.5*(sqr(c)-sqr(s))/s/c;
  _coup[ 7] =-fact/sqr(cw)*sw*0.5*(sqr(cp)-sqr(sp))/sp/cp;
  _coup[ 8] =-fact/cw*sw*0.5/s/c/sp/cp*(sqr(c*sp)+sqr(s*cp));
  _coup[ 9] =-fact*(s0-2.*r2*vr);
  _coup[10] = fact*(s0-2.*r2*vr);
  _coup[11] = fact*(s0-2.*r2*vr)*0.5*(sqr(c)-sqr(s))/s/c;
  _coup[12] =-fact/sqr(cw)*(s0-4.*r2*vr);
  _coup[13] = fact*(s0+sqr(sqr(c)-sqr(s))/sqr(s*c)*r2*vr);
  _coup[14] = fact/cw*0.5*(sqr(c)-sqr(s))/s/c*(s0-4.*r2*vr);
  _coup[15] = fact*sw/sqr(cw)*0.5*(sqr(cp)-sqr(sp))/sp/cp*(s0-4.*r2*vr);
  _coup[16] = fact*sw/cw*0.5/s/c/sp/cp*(s0*(sqr(c*sp)+sqr(s*cp))
					+2.*r2*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))*vr);
  _coup[17] = fact*sqr(sw/cw)*(s0+r2*vr*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp));
  _coup[18] =-2.*fact/cw*vr;
  _coup[19] = fact/cw*(sqr(c)-sqr(s))/s/c*vr;
  _coup[20] =-fact*sw/cw*0.5*(sqr(cp)-sqr(sp))/sp/cp*(sPlus-4.*vr);
  _coup[21] =-fact*sw/cw*(sqr(c*cp)+sqr(s*sp))/s/c/sp/cp*vr;
  _coup[22] = fact*(sqr(c)-sqr(s))/s/c*vr;
  _coup[23] =-fact*(pow(c,4)+pow(s,4))/sqr(s*c)*vr;
  _coup[24] = fact*4.*vr;
  _coup[25] = fact*2.*(pow(c,4)+pow(s,4))/sqr(s*c)*vr;
  _coup[26] =-fact*2.*vr*(sqr(c)-sqr(s))/s/c;
}

void LHWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  int ih = abs(c->id());
  long int ibos[2]={abs(a->id()),abs(b->id())};
  if(ih==25) {
    if     ( ibos[0]==24&&ibos[1]==24     ) norm(UnitRemoval::InvE *_couplast*_coup[0]);
    else if( ibos[0]==23&&ibos[1]==23     ) norm(UnitRemoval::InvE *_couplast*_coup[1]);
    else if( ibos[0]==34&&ibos[1]==34     ) norm(UnitRemoval::InvE *_couplast*_coup[2]);
    else if( ibos[0]==33&&ibos[1]==33     ) norm(UnitRemoval::InvE *_couplast*_coup[3]);
    else if( ibos[0]==32&&ibos[1]==32     ) norm(UnitRemoval::InvE *_couplast*_coup[4]);
    else if((ibos[0]==24&&ibos[1]==34) ||
	    (ibos[0]==34&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[5]);
    else if((ibos[0]==23&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[6]);
    else if((ibos[0]==23&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[7]);
    else if((ibos[0]==33&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==33)    ) norm(UnitRemoval::InvE *_couplast*_coup[8]);
    else assert(false);
  }
  else if(ih==35) {
    if     ( ibos[0]==24&&ibos[1]==24     ) norm(UnitRemoval::InvE *_couplast*_coup[ 9]);
    else if( ibos[0]==34&&ibos[1]==34     ) norm(UnitRemoval::InvE *_couplast*_coup[10]);
    else if((ibos[0]==24&&ibos[1]==34) ||
	    (ibos[0]==34&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[11]);
    else if( ibos[0]==23&&ibos[1]==23     ) norm(UnitRemoval::InvE *_couplast*_coup[12]);
    else if( ibos[0]==33&&ibos[1]==33     ) norm(UnitRemoval::InvE *_couplast*_coup[13]);
    else if((ibos[0]==23&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[14]);
    else if((ibos[0]==23&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==23)    ) norm(UnitRemoval::InvE *_couplast*_coup[15]);
    else if((ibos[0]==33&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==33)    ) norm(UnitRemoval::InvE *_couplast*_coup[16]);
    else if((ibos[0]==32&&ibos[1]==32)    ) norm(UnitRemoval::InvE *_couplast*_coup[17]);
    else assert(false);
  }
  else if(ih==37) {
    if     ((ibos[0]==24&&ibos[1]==23) ||
	    (ibos[0]==23&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[18]);
    else if((ibos[0]==34&&ibos[1]==23) ||
	    (ibos[0]==23&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[19]);
    else if((ibos[0]==24&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[20]);
    else if((ibos[0]==34&&ibos[1]==32) ||
	    (ibos[0]==32&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[21]);
    else if((ibos[0]==24&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[22]);
    else if((ibos[0]==34&&ibos[1]==33) ||
	    (ibos[0]==33&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[23]);
    else assert(false);
  }
  else if(ih==38) {
    if     ((ibos[0]==24&&ibos[1]==24)    ) norm(UnitRemoval::InvE *_couplast*_coup[24]);
    else if((ibos[0]==34&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[24]);
    else if((ibos[0]==34&&ibos[1]==24) ||
	    (ibos[0]==24&&ibos[1]==34)    ) norm(UnitRemoval::InvE *_couplast*_coup[24]);
    else assert(false);
  }
  else assert(false);
}
#line 1 "./LHWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWHHVertex class.
//

#include "LHWHHVertex.h"
#include "LHModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHWHHVertex::LHWHHVertex() : 
  couplast_(0.), q2last_(ZERO), coup_(24) {
  orderInGs(0);
  orderInGem(1);
  // neutral
  addToList( 22, 37,-37);
  addToList( 22, 38,-38);
  addToList( 32, 37,-37);
  addToList( 32, 38,-38);
  addToList( 23, 37,-37);
  addToList( 23, 38,-38);
  addToList( 33, 37,-37);
  addToList( 33, 38,-38);
  addToList( 32, 25, 36);
  addToList( 32, 35, 36);
  addToList( 23, 25, 36);
  addToList( 23, 35, 36);
  addToList( 33, 25, 36);
  addToList( 33, 35, 36);
  // W+
  addToList( 24, 25,-37);
  addToList( 24, 35,-37);
  addToList( 24, 36,-37);
  addToList( 24, 37,-38);
  addToList( 34, 25,-37);
  addToList( 34, 35,-37);
  addToList( 34, 36,-37);
  addToList( 34, 37,-38);
  // W-
  addToList(-24, 25, 37);
  addToList(-24, 35, 37);
  addToList(-24, 36, 37);
  addToList(-24,-37, 38);
  addToList(-34, 25, 37);
  addToList(-34, 35, 37);
  addToList(-34, 36, 37);
  addToList(-34,-37, 38);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void LHWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWHHVertex,VSSVertex>
describeHerwigLHWHHVertex("Herwig::LHWHHVertex", "HwLHModel.so");

void LHWHHVertex::Init() {

  static ClassDocumentation<LHWHHVertex> documentation
    ("The LHWHHVertex class implements the coupling of a pair of Higgs"
     " bosons to an electroweak gauge boson in the Little Higgs model.");

}

void LHWHHVertex::doinit() {
  VSSVertex::doinit();
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWHHVertex::doinit()"
			  << Exception::runerror;
  double sw2(sin2ThetaW());
  double sw(sqrt(sw2)),cw(sqrt(1.-sw2));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double s0   (model->sinTheta0());
  double sP   (model->sinThetaP());
  double sPlus(model->sinThetaPlus());
  coup_[ 0] = 0.5/sw*(sqrt(2.)*s0-sPlus);
  coup_[ 1] = sqrt(0.5)/sw;
  coup_[ 2] = Complex(0.,1.)/sw*sqrt(0.5);
  coup_[ 3] = 1./sw;
  coup_[ 4] = 0.;
  coup_[ 5] = 0.;
  coup_[ 6] = 1.;
  coup_[ 7] = 2.;
  coup_[ 8] = Complex(0.,0.5)/cw/sw*(sP-2.*s0);	       
  coup_[ 9] =-Complex(0.,1.)/cw/sw;
  coup_[10] =-sw/cw;
  coup_[11] = (1.-2.*sw2)/cw/sw;
  coup_[12] =-0.25/sw*(sqr(c)-sqr(s))/s/c*(sqrt(2.)*s0-sPlus);
  coup_[13] =-sqrt(0.5)/sw*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[14] =-Complex(0.,1.)*sqrt(0.5)*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[15] =-0.5*(sqr(c)-sqr(s))/s/c/sw;
  coup_[16] =-Complex(0.,0.5)/cw*0.5*(sqr(cp)-sqr(sp))/sp/cp*(sP-2.*s0);
  coup_[17] = Complex(0.,1.)/cw*0.5*(sqr(cp)-sqr(sp))/sp/cp;
  coup_[18] =-0.5*(sqr(cp)-sqr(sp))/sp/cp/cw;
  coup_[19] =-0.5*(sqr(cp)-sqr(sp))/sp/cp/cw;
  coup_[20] =-Complex(0.,0.5)/sw*0.5*(sqr(c)-sqr(s))/s/c*(sP-2.*s0);
  coup_[21] = Complex(0.,1.)/sw*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[22] = 0.;
  coup_[23] =-0.5/sw*(sqr(c)-sqr(s))/s/c;
}

void LHWHHVertex::setCoupling(Energy2 q2, tcPDPtr particle1,
			      tcPDPtr particle2, tcPDPtr particle3) {
  if( q2 != q2last_ || couplast_==0.) {
    q2last_ = q2;
    couplast_ = -electroMagneticCoupling(q2);
  }
  int ibos = particle1->id();
  int isc1 = particle2->id();
  int isc2 = particle3->id();
  if(ibos==ParticleID::gamma) {
    if(isc1==37) 
      norm(coup_[6]*couplast_);
    else if(isc1==38)
      norm(coup_[7]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[6]*couplast_);
    else if(isc1==-38)
      norm(-coup_[7]*couplast_);
    else
      assert(false);
  }
  else if(ibos==32) {
    if(isc1==37) 
      norm(coup_[18]*couplast_);
    else if(isc1==38)
      norm(coup_[19]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[18]*couplast_);
    else if(isc1==-38)
      norm(-coup_[19]*couplast_);
    else if(isc1==25)
      norm(coup_[16]*couplast_);
    else if(isc1==35)
      norm(coup_[17]*couplast_);
    else if(isc2==25)
      norm(-coup_[16]*couplast_);
    else if(isc2==35)
      norm(-coup_[17]*couplast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Z0) {
    if(isc1==37) 
      norm(coup_[10]*couplast_);
    else if(isc1==38)
      norm(coup_[11]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[10]*couplast_);
    else if(isc1==-38)
      norm(-coup_[11]*couplast_);
    else if(isc1==25)
      norm(coup_[8]*couplast_);
    else if(isc1==35)
      norm(coup_[9]*couplast_);
    else if(isc2==25)
      norm(-coup_[8]*couplast_);
    else if(isc2==35)
      norm(-coup_[9]*couplast_);
    else
      assert(false);
  }
  else if(ibos==33) {
    if(isc1==37) 
      norm(coup_[22]*couplast_);
    else if(isc1==38)
      norm(coup_[23]*couplast_);
    else if(isc1==-37) 
      norm(-coup_[22]*couplast_);
    else if(isc1==-38)
      norm(-coup_[23]*couplast_);
    else if(isc1==25)
      norm(coup_[20]*couplast_);
    else if(isc1==35)
      norm(coup_[21]*couplast_);
    else if(isc2==25)
      norm(-coup_[20]*couplast_);
    else if(isc2==35)
      norm(-coup_[21]*couplast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wplus) {
    if(isc1==25)
      norm(coup_[0]*couplast_);
    else if(isc1==35)
      norm(coup_[1]*couplast_);
    else if(isc1==36)
      norm(coup_[2]*couplast_);
    else if(isc1==37)
      norm(coup_[3]*couplast_);
    else if(isc2==25)
      norm(-coup_[0]*couplast_);
    else if(isc2==35)
      norm(-coup_[1]*couplast_);
    else if(isc2==36)
      norm(-coup_[2]*couplast_);
    else if(isc2==37)
      norm(-coup_[3]*couplast_);
    else
      assert(false);
  }
  else if(ibos==34) {
    if(isc1==25)
      norm(coup_[12]*couplast_);
    else if(isc1==35)
      norm(coup_[13]*couplast_);
    else if(isc1==36)
      norm(coup_[14]*couplast_);
    else if(isc1==37)
      norm(coup_[15]*couplast_);
    else if(isc2==25)
      norm(-coup_[12]*couplast_);
    else if(isc2==35)
      norm(-coup_[13]*couplast_);
    else if(isc2==36)
      norm(-coup_[14]*couplast_);
    else if(isc2==37)
      norm(-coup_[15]*couplast_);
    else
      assert(false);
  }
  else if(ibos==ParticleID::Wminus) {
    if(isc1==25)
      norm(conj(coup_[0])*couplast_);
    else if(isc1==35)
      norm(conj(coup_[1])*couplast_);
    else if(isc1==36)
      norm(conj(coup_[2])*couplast_);
    else if(isc1==37)
      norm(conj(coup_[3])*couplast_);
    else if(isc2==25)
      norm(-conj(coup_[0])*couplast_);
    else if(isc2==35)
      norm(-conj(coup_[1])*couplast_);
    else if(isc2==36)
      norm(-conj(coup_[2])*couplast_);
    else if(isc2==37)
      norm(-conj(coup_[3])*couplast_);
    else
      assert(false);
  }
  else if(ibos==-34) {
    if(isc1==25)
      norm(conj(coup_[12])*couplast_);
    else if(isc1==35)
      norm(conj(coup_[13])*couplast_);
    else if(isc1==36)
      norm(conj(coup_[14])*couplast_);
    else if(isc1==-37)
      norm(conj(coup_[15])*couplast_);
    else if(isc2==25)
      norm(-conj(coup_[12])*couplast_);
    else if(isc2==35)
      norm(-conj(coup_[13])*couplast_);
    else if(isc2==36)
      norm(-conj(coup_[14])*couplast_);
    else if(isc2==-37)
      norm(-conj(coup_[15])*couplast_);
    else
      assert(false);
  }
}
#line 1 "./LHWWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHWWHHVertex class.
//

#include "LHWWHHVertex.h"
#include "LHModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

LHWWHHVertex::LHWWHHVertex() : 
  couplast_(0.), q2last_(ZERO), coup_(107) {
  orderInGs(0);
  orderInGem(2);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr LHWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LHWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void LHWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << coup_;
}

void LHWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> coup_;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<LHWWHHVertex,VVSSVertex>
describeHerwigLHWWHHVertex("Herwig::LHWWHHVertex", "HwLHModel.so");

void LHWWHHVertex::Init() {

  static ClassDocumentation<LHWWHHVertex> documentation
    ("The LHWWHHVertex class implements the couplings of a pair"
     " of electroweak gauge bosons and a pair of Higgs bosons in"
     " the Little Higgs model.");

}

void LHWWHHVertex::doinit() {
  // VVHH
  addToList(  24,  -24,  25,  25);
  addToList(  23,   23,  25,  25);
  addToList(  24,  -34,  25,  25);
  addToList(  34,  -24,  25,  25);
  addToList(  23,   32,  25,  25);
  addToList(  34,  -34,  25,  25);
  addToList(  33,   33,  25,  25);
  addToList(  32,   32,  25,  25);
  addToList(  23,   33,  25,  25);
  addToList(  32,   33,  25,  25);
  // VVH Phi_0
  addToList(  24,  -24,  25,  35);
  addToList(  23,   23,  25,  35);
  addToList(  24,  -34,  25,  35);
  addToList(  34,  -24,  25,  35);
  addToList(  23,   32,  25,  35);
  addToList(  34,  -34,  25,  35);
  addToList(  33,   33,  25,  35);
  addToList(  32,   32,  25,  35);
  addToList(  23,   33,  25,  35);
  addToList(  32,   33,  25,  35);
  // VV Phi_0 Phi_0
  addToList(  24,  -24,  35,  35);
  addToList(  23,   23,  35,  35);
  addToList(  24,  -34,  35,  35);
  addToList(  34,  -24,  35,  35);
  addToList(  23,   32,  35,  35);
  addToList(  34,  -34,  35,  35);
  addToList(  33,   33,  35,  35);
  addToList(  32,   32,  35,  35);
  addToList(  23,   33,  35,  35);
  addToList(  32,   33,  35,  35);
  // VV Phi_P Phi_P
  addToList(  24,  -24,  36,  36);
  addToList(  23,   23,  36,  36);
  addToList(  24,  -34,  36,  36);
  addToList(  34,  -24,  36,  36);
  addToList(  23,   32,  36,  36);
  addToList(  34,  -34,  36,  36);
  addToList(  33,   33,  36,  36);
  addToList(  32,   32,  36,  36);
  addToList(  23,   33,  36,  36);
  addToList(  32,   33,  36,  36);
  // VV Phi+ Phi-
  addToList(  24,  -24,  37, -37);
  addToList(  23,   23,  37, -37);
  addToList(  22,   22,  37, -37);
  addToList(  22,   23,  37, -37);
  addToList(  24,  -34,  37, -37);
  addToList(  34,  -24,  37, -37);
  addToList(  34,  -34,  37, -37);
  addToList(  33,   33,  37, -37);
  addToList(  32,   32,  37, -37);
  addToList(  32,   33,  37, -37);
  addToList(  22,   32,  37, -37);
  addToList(  23,   33,  37, -37);
  addToList(  23,   32,  37, -37);
  // VV Phi++ Phi--
  addToList(  24,  -24,  38, -38);
  addToList(  23,   23,  38, -38);
  addToList(  22,   22,  38, -38);
  addToList(  22,   23,  38, -38);
  addToList(  24,  -34,  38, -38);
  addToList(  34,  -24,  38, -38);
  addToList(  34,  -34,  38, -38);
  addToList(  33,   33,  38, -38);
  addToList(  32,   32,  38, -38);
  addToList(  32,   33,  38, -38);
  addToList(  22,   32,  38, -38);
  addToList(  23,   33,  38, -38);
  addToList(  23,   32,  38, -38);
  // VV H phi-  + cc
  addToList(  24,   22,  25, -37);
  addToList(  24,   23,  25, -37);
  addToList(  24,   32,  25, -37);
  addToList(  24,   33,  25, -37);
  addToList(  34,   22,  25, -37);
  addToList(  34,   23,  25, -37);
  addToList(  34,   32,  25, -37);
  addToList(  34,   33,  25, -37);
  addToList( -24,   22,  25,  37);
  addToList( -24,   23,  25,  37);
  addToList( -24,   32,  25,  37);
  addToList( -24,   33,  25,  37);
  addToList( -34,   22,  25,  37);
  addToList( -34,   23,  25,  37);
  addToList( -34,   32,  25,  37);
  addToList( -34,   33,  25,  37);
  // VV phi0  phi-  + cc
  addToList(  24,   22,  35, -37);
  addToList(  24,   23,  35, -37);
  addToList(  24,   32,  35, -37);
  addToList(  24,   33,  35, -37);
  addToList(  34,   22,  35, -37);
  addToList(  34,   23,  35, -37);
  addToList(  34,   32,  35, -37);
  addToList(  34,   33,  35, -37);
  addToList( -24,   22,  35,  37);
  addToList( -24,   23,  35,  37);
  addToList( -24,   32,  35,  37);
  addToList( -24,   33,  35,  37);
  addToList( -34,   22,  35,  37);
  addToList( -34,   23,  35,  37);
  addToList( -34,   32,  35,  37);
  addToList( -34,   33,  35,  37);
  // VV phiP  phi-  + cc
  addToList(  24,   22,  36, -37);
  addToList(  24,   23,  36, -37);
  addToList(  24,   32,  36, -37);
  addToList(  24,   33,  36, -37);
  addToList(  34,   22,  36, -37);
  addToList(  34,   23,  36, -37);
  addToList(  34,   32,  36, -37);
  addToList(  34,   33,  36, -37);
  addToList( -24,   22,  36,  37);
  addToList( -24,   23,  36,  37);
  addToList( -24,   32,  36,  37);
  addToList( -24,   33,  36,  37);
  addToList( -34,   22,  36,  37);
  addToList( -34,   23,  36,  37);
  addToList( -34,   32,  36,  37);
  addToList( -34,   33,  36,  37);
  // VV phi+ phi -- + cc
  addToList(  24,   22,  37, -38);
  addToList(  24,   23,  37, -38);
  addToList(  24,   32,  37, -38);
  addToList(  24,   33,  37, -38);
  addToList(  34,   22,  37, -38);
  addToList(  34,   23,  37, -38);
  addToList(  34,   32,  37, -38);
  addToList(  34,   33,  37, -38);
  addToList( -24,   22, -37,  38);
  addToList( -24,   23, -37,  38);
  addToList( -24,   32, -37,  38);
  addToList( -24,   33, -37,  38);
  addToList( -34,   22, -37,  38);
  addToList( -34,   23, -37,  38);
  addToList( -34,   32, -37,  38);
  addToList( -34,   33, -37,  38);
  // VV H phi-- + cc
  addToList(  24,   24,  25, -38);
  addToList( -24,  -24,  25,  38);
  addToList(  24,   34,  25, -38);
  addToList( -24,  -34,  25,  38);
  addToList(  34,   34,  25, -38);
  addToList( -34,  -34,  25,  38);
  // VV phi0 phi-- + cc
  addToList(  24,   24,  35, -38);
  addToList( -24,  -24,  35,  38);
  addToList(  24,   34,  35, -38);
  addToList( -24,  -34,  35,  38);
  addToList(  34,   34,  35, -38);
  addToList( -34,  -34,  35,  38);
  // VV phiP  phi-- + cc
  addToList(  24,   24,  36, -38);
  addToList( -24,  -24,  36,  38);
  addToList(  24,   34,  36, -38);
  addToList( -24,  -34,  36,  38);
  addToList(  34,   34,  36, -38);
  addToList( -34,  -34,  36,  38);
  VVSSVertex::doinit();
  // model
  cLHModelPtr model = 
    dynamic_ptr_cast<cLHModelPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must be using the LHModel "
			  << " in LHWWWWVertex::doinit()"
			  << Exception::runerror;
  double sw2(sin2ThetaW()),cw2(1.-sw2);
  double sw(sqrt(sw2)),cw(sqrt(cw2));
  double s (model->sinTheta()     ),c (model->cosTheta()     );
  double sp(model->sinThetaPrime()),cp(model->cosThetaPrime());
  double s0   (model->sinTheta0());
  double sPlus(model->sinThetaPlus());
  // VV HH
  coup_[  0] = 0.5/sw2;
  coup_[  1] = 0.5/sw2/cw2;
  coup_[  2] = 0.;
  coup_[  3] =-0.25/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[  4] =-0.25/sw/cw2*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[  5] =-0.5/sw2;
  coup_[  6] =-0.5/sw2;
  coup_[  7] =-0.5/cw2;
  coup_[  8] =-0.25/sw2/cw*(sqr(c)-sqr(s))/c/s;
  coup_[  9] =-0.25/sw/cw*(sqr(c*sp)+sqr(s*cp))/c/s/sp/cp;
  // VV H Phi_0
  coup_[ 10] = 0.5*s0/sw2;
  coup_[ 11] = 1.5*s0/sw2/cw2;
  coup_[ 12] = 0.;
  coup_[ 13] =-0.25/sw2*(sqr(c)-sqr(s))/c/s*s0;
  coup_[ 14] =-0.75/sw/cw2*(sqr(cp)-sqr(sp))/cp/sp*s0;
  coup_[ 15] = -0.5/sw2*s0;
  coup_[ 16] = 0.5/sw2*(1.+sqr(sqr(c)-sqr(s))/sqr(s*c))*s0;
  coup_[ 17] = 0.5/cw2*(1.+sqr(sqr(cp)-sqr(sp))/sqr(sp*cp))*s0;
  coup_[ 18] =-0.75/cw/sw2*(sqr(c)-sqr(s))/c/s*s0;
  coup_[ 19] = 0.25/sw/cw/c/s/sp/cp*((sqr(c*sp)+sqr(s*cp))
				    +2.*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp)))*s0;
  // VV phi0 phi0
  coup_[ 20] = 1./sw2;
  coup_[ 21] = 2./cw2/sw2;
  coup_[ 22] = 0.;
  coup_[ 23] =-0.5/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 24] =-1./sw/cw2*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 25] =-1./sw2;
  coup_[ 26] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 26] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 28] =-1./cw/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 29] = 0.5/cw/sw*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))/s/c/sp/cp;
  // VV phi_P phi_P
  coup_[ 30] = 1./sw2;
  coup_[ 31] = 2./sw2/cw2;
  coup_[ 32] = 0.;
  coup_[ 33] =-0.5/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 34] =-1./sw/cw2*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 35] =-1./sw2;
  coup_[ 36] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 37] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 38] =-1./cw/sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 39] = 0.5/cw/sw*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))/s/c/sp/cp;
  // VV phi+ phi-
  coup_[ 40] = 2./sw2;
  coup_[ 41] = 2.*sw2/cw2;
  coup_[ 42] = 2.;
  coup_[ 43] =-2.*sw/cw;
  coup_[ 44] =-1./sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 45] = 0.;
  coup_[ 46] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 47] =-0.5/sw2/sqr(s*c);
  coup_[ 48] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 49] = 0.;
  coup_[ 50] =-1./cw*(sqr(cp)-sqr(sp))/sp/cp;
  coup_[ 51] = 0.;
  coup_[ 52] = sw/cw2*(sqr(cp)-sqr(sp))/sp/cp;
  // VV phi++ phi--
  coup_[ 53] = 1./sw2;
  coup_[ 54] = 2./cw2/sw2*sqr(1.-2.*sw2);
  coup_[ 55] = 8.;
  coup_[ 56] = 4./sw/cw*(1.-2.*sw2);
  coup_[ 57] =-0.5/sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 58] = 2./sw*(sqr(c)-sqr(s))/s/c;
  coup_[ 59] =-1./sw2;
  coup_[ 60] = 0.5/sw2*sqr(sqr(c)-sqr(s))/sqr(s*c);
  coup_[ 61] = 0.5/cw2*sqr(sqr(cp)-sqr(sp))/sqr(sp*cp);
  coup_[ 62] =-0.5/cw/sw*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp))/s/c/sp/cp;
  coup_[ 63] =-2./cw*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 64] = 1./sw2/cw*(sqr(c )-sqr(s ))/s /c *(1.-2.*sw2);
  coup_[ 65] =-1./cw2/sw*(sqr(cp)-sqr(sp))/sp/cp*(1.-2.*sw2);
  // VV h phi-
  coup_[ 66] =-0.5/sw*(sPlus-sqrt(2.)*s0);
  coup_[ 67] = 0.5/cw/sw2*(sPlus*sw2-sqrt(2.)*s0*(1.+sw2));
  coup_[ 68] =-0.25/sw/cw*(sqr(cp)-sqr(sp))/cp/sp*(sPlus-2.*sqrt(2.)*s0);
  coup_[ 69] = 0.25/sw2*(sqr(c)-sqr(s))/s/c*s0;
  coup_[ 70] = 0.25/sw*(sqr(c)-sqr(s))/s/c*(sPlus-sqrt(2.)*s0);
  coup_[ 71] =-0.25/sw2/cw*(sqr(c)-sqr(s))/s/c*(sPlus*sw2-sqrt(2.)*s0*(1.+sw2));
  coup_[ 72] =-0.25/sw/cw/s/c/sp/cp*(sPlus*(sqr(c*sp)+sqr(s*cp))
				    +sqrt(2.)*(sqr(c)-sqr(s))*(sqr(cp)-sqr(sp)));
  coup_[ 73] =-0.25/sw2*s0*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phi0 phi-
  coup_[ 74] =-sqrt(0.5)/sw;
  coup_[ 75] =-sqrt(0.5)/sw2/cw*(1.+sw2);
  coup_[ 76] = sqrt(0.5)/sw/cw*(sqr(cp)-sqr(sp))/sp/cp;
  coup_[ 77] = 0.5*sqrt(0.5)/sw2*(sqr(c)-sqr(s))/c/s;
  coup_[ 78] = 0.5*sqrt(0.5)/sw*(sqr(c)-sqr(s))/c/s;
  coup_[ 79] = 0.5*sqrt(0.5)/sw2/cw*(sqr(c)-sqr(s))/c/s*(1.+sw2);
  coup_[ 80] =-0.5*sqrt(0.5)/sw/cw*(sqr(cp)-sqr(sp))/cp/sp*(sqr(c)-sqr(s))/c/s;
  coup_[ 81] =-0.5*sqrt(0.5)/sw2*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phi+ phi--
  coup_[ 82] = 3./sw;
  coup_[ 83] = (1.-3.*sw2)/cw/sw2;
  coup_[ 84] = 1./sw/cw*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 85] = Complex(0.,1.)*0.5*sqrt(0.5)/sw2*(sqr(c)-sqr(s))/s/c;
  coup_[ 86] =-3./sw*(sqr(c)-sqr(s))/s/c;
  coup_[ 87] =-0.5/sw2/cw*(sqr(c)-sqr(s))/s/c*(1.-3.*sw2);
  coup_[ 88] =-0.5/sw/cw*(sqr(c)-sqr(s))/s/c*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 89] =-Complex(0.,1.)*0.5*sqrt(0.5)/sw2*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phip phi-
  coup_[ 90] =-Complex(0.,1.)/sw*sqrt(0.5);
  coup_[ 91] =-Complex(0.,1.)/sw2/cw*(1.+sw2);
  coup_[ 92] = Complex(0.,1.)/sw/cw*sqrt(0.5)*(sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 93] = Complex(0.,1.)/sw2*sqrt(0.5)*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[ 94] =-Complex(0.,1.)*sqrt(0.5)*0.5/sw*(sqr(c)-sqr(s))/s/c;
  coup_[ 95] = Complex(0.,1.)*sqrt(0.5)*0.5/sw2/cw*(sqr(c)-sqr(s))/s/c*(1.+sw2);
  coup_[ 96] =-Complex(0.,1.)*sqrt(0.5)*0.5/sw/cw*(sqr(c )-sqr(s ))/s/c*
                                                 (sqr(cp)-sqr(sp))/cp/sp;
  coup_[ 97] =-Complex(0.,1.)/sw2*sqrt(0.5)*0.5*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV H phi--
  coup_[ 98] = sqrt(2.)/sw2*s0;
  coup_[ 99] =-sqrt(2.)/sw2*0.5*(sqr(c)-sqr(s))/s/c*s0;
  coup_[100] = sqrt(2.)/sw2*0.5*(pow(c,4)+pow(s,4))/sqr(s*c)*s0;
  // VV phi0 phi--
  coup_[101] = sqrt(2.)/sw2;
  coup_[102] =-sqrt(2.)/sw2*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[103] = sqrt(2.)/sw2*0.5*(pow(c,4)+pow(s,4))/sqr(s*c);
  // VV phip phi--
  coup_[104] = Complex(0.,1.)*sqrt(2.)/sw2;
  coup_[105] =-Complex(0.,1.)*sqrt(2.)/sw2*0.5*(sqr(c)-sqr(s))/s/c;
  coup_[106] = Complex(0.,1.)*sqrt(2.)/sw2*0.5*(pow(c,4)+pow(s,4))/sqr(s*c);
}

void LHWWHHVertex::setCoupling(Energy2 q2,
			       tcPDPtr part1,tcPDPtr part2,
			       tcPDPtr part3,tcPDPtr part4) {
  if( q2 != q2last_ || couplast_==0.) {
    q2last_ = q2;
    couplast_ = sqr(electroMagneticCoupling(q2));
  }
  int ibos1 = part1->id();
  int ibos2 = part2->id();
  int isca1 = part3->id();
  int isca2 = part4->id();
  if( isca1 == isca2 || 
      (isca1==25&&isca2==35) || (isca1==35&&isca2==25)) {
    unsigned int ioff = 0;
    if     (isca1!=isca2) ioff = 10;
    else if(isca1==35   ) ioff = 20;
    else if(isca1==36   ) ioff = 30;
    if(ibos1==23&&ibos2==23)
      norm(coup_[1+ioff]*couplast_);
    else if(ibos1==33&&ibos2==33)
      norm(coup_[6+ioff]*couplast_);
    else if(ibos1==33&&ibos2==33)
      norm(coup_[7+ioff]*couplast_);
    else if(abs(ibos1)==24&&abs(ibos2)==24)
      norm(coup_[0+ioff]*couplast_);
    else if(abs(ibos1)==34&&abs(ibos2)==34)
      norm(coup_[5+ioff]*couplast_);
    else if(( abs(ibos1) == 24 && abs(ibos2) == 34) ||
	    ( abs(ibos1) == 34 && abs(ibos2) == 24))
      norm(coup_[3+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 32) ||
	    ( ibos1 == 32 && ibos2 == 23))
      norm(coup_[4+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 23))
      norm(coup_[8+ioff]*couplast_);
    else if(( ibos1 == 32 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 32))
      norm(coup_[9+ioff]*couplast_);
    else
      assert(false);
  }
  else if(isca1==-isca2) {
    unsigned int ioff = abs(isca1) == 37 ? 40 : 53;
    if(abs(ibos1)==24&&abs(ibos2)==24)
      norm(coup_[0+ioff]*couplast_);
    else if(ibos1==23&&ibos2==23)
      norm(coup_[1+ioff]*couplast_);
    else if(ibos1==22&&ibos2==22)
      norm(coup_[2+ioff]*couplast_);
    else if(( ibos1 == 22 && ibos2 == 23) ||
	    ( ibos1 == 23 && ibos2 == 22))
      norm(coup_[3+ioff]*couplast_);
    else if(( abs(ibos1) == 24 && abs(ibos2) == 34) ||
	    ( abs(ibos1) == 34 && abs(ibos2) == 24))
      norm(coup_[4+ioff]*couplast_);
    else if(( ibos1 == 22 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 22))
      norm(coup_[5+ioff]*couplast_);
    else if(abs(ibos1)==34&&abs(ibos2)==34)
      norm(coup_[6+ioff]*couplast_);
    else if(ibos1==33&&ibos2==33)
      norm(coup_[7+ioff]*couplast_);
    else if(ibos1==32&&ibos2==32)
      norm(coup_[8+ioff]*couplast_);
    else if(( ibos1 == 32 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 32))
      norm(coup_[9+ioff]*couplast_);
    else if(( ibos1 == 22 && ibos2 == 32) ||
	    ( ibos1 == 32 && ibos2 == 22))
      norm(coup_[10+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 33) ||
	    ( ibos1 == 33 && ibos2 == 23))
      norm(coup_[11+ioff]*couplast_);
    else if(( ibos1 == 23 && ibos2 == 32) ||
	    ( ibos1 == 32 && ibos2 == 23))
      norm(coup_[12+ioff]*couplast_);
    else
      assert(false);
  }
  else if(((abs(ibos1) == 24 || abs(ibos1) == 34) && 
	   (abs(ibos2) != 24 && abs(ibos2) != 34)) ||
	  ((abs(ibos2) == 24 || abs(ibos2) == 34) && 
	   (abs(ibos1) != 24 && abs(ibos1) != 34))) {
    int iw,ineut;
    if(abs(ibos1) == 24 || abs(ibos1) == 34) {
      iw    = abs(ibos1);
      ineut = ibos2;
    }
    else {
      iw    = abs(ibos2);
      ineut = ibos1;
    }
    unsigned int ioff = 66;
    if((isca1 == 35 && abs(isca2) == 37) ||
       (isca2 == 35 && abs(isca1) == 37)) {
      ioff += 8;
    }
    else if ((abs(isca1) == 37 && abs(isca2) == 38) ||
	     (abs(isca2) == 37 && abs(isca1) == 38)) {
      ioff += 16;
    } 
    else if ((isca1 == 35 && abs(isca2) == 37) ||
	     (isca2 == 35 && abs(isca1) == 37)) {
      ioff += 24;
    }
    else
      assert(false);
    if(iw==34) ioff += 4;
    if(ineut==22)
      norm(coup_[0+ioff]*couplast_);
    else if(ineut==23)
      norm(coup_[1+ioff]*couplast_);
    else if(ineut==32)
      norm(coup_[2+ioff]*couplast_);
    else if(ineut==33)
      norm(coup_[3+ioff]*couplast_);
    else
      assert(false);
  }
  else {
    unsigned int ioff = 98;
    if(isca1==25||isca2==25)
      ioff += 0;
    else if(isca1==35||isca2==35)
      ioff += 3;
    else if(isca1==36||isca2==36)
      ioff += 6;
    else
      assert(false);
    if(ibos1==ibos2) {
      if(abs(ibos1)==24) {
	norm(coup_[0+ioff]*couplast_);
      }
      else if(abs(ibos1)==34) {
	norm(coup_[1+ioff]*couplast_);
      }
      else
	assert(false);
    }
    else {
      norm(coup_[2+ioff]*couplast_);
    }
  }
}
