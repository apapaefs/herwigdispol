#line 1 "./BallZwickyScalarFormFactor.cc"
// -*- C++ -*-
//
// BallZwickyScalarFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyScalarFormFactor class.
//

#include "BallZwickyScalarFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

BallZwickyScalarFormFactor::BallZwickyScalarFormFactor()  
  : _r10(7), _r20(7), _r1plus(7), _r2plus(7), _r1T(7), _r2T(7), _m120(7),
    _mfit20(7), _m12plus(7), _mfit2plus(7), _m12T(7), _mfit2T(7) {
  // parameters for the B to pi  form-factors
  for(unsigned int ix=0;ix<4;++ix) {
    _r10[ix]     = 0.            ; _r20[ix]       = 0.258; 
    _m120[ix]    = -1.*GeV2      ; _mfit20[ix]    = 33.81*GeV2; 
    _r1plus[ix]  = 0.744         ; _r2plus[ix]    = -0.486; 
    _m12plus[ix] = sqr(5.32)*GeV2; _mfit2plus[ix] = 40.73*GeV2; 
    _r1T[ix]     = 1.387         ; _r2T[ix]       = -1.134; 
    _m12T[ix]    = sqr(5.32)*GeV2; _mfit2T[ix]    = 32.22*GeV2; 
  }
  addFormFactor(-521, 111,0,-2,5,2);
  addFormFactor(-511, 111,0,-2,5,1);
  addFormFactor(-511, 211,0,-2,5,2);  
  addFormFactor(-521, 211,0,-2,5,1);
  for(unsigned int ix=0;ix<2;++ix) {
    double fact(sqrt(0.5));
    if(ix==1) fact *= -1.;
    _r20[ix] *= fact; _r1plus[ix] *= fact; _r2plus[ix] *= fact; 
    _r1T[ix] *= fact; _r2T[ix]    *= fact; 
  }
  // parameters for the B to K   form-factors
  addFormFactor(-521,-321,0,-2,5,3);
  addFormFactor(-511,-311,0,-2,5,3);
  for(unsigned int ix=4;ix<6;++ix) {
    _r10[ix]     = 0.            ; _r20[ix]       = 0.330; 
    _m120[ix]    = -1.*GeV2      ; _mfit20[ix]    = 37.46*GeV2; 
    _r1plus[ix]  = 0.162         ; _r2plus[ix]    = 0.173; 
    _m12plus[ix] = sqr(5.41)*GeV2; _mfit2plus[ix] = -1.*GeV2; 
    _r1T[ix]     = 0.161         ; _r2T[ix]       = 0.198; 
    _m12T[ix]    = sqr(5.41)*GeV2; _mfit2T[ix]    = -1.*GeV2; 
  }
  // parameters for the B to eta form-factors
  addFormFactor(521,221,0,2,-5,-2); 
  _r10[6]     = 0.            ; _r20[6]       = 0.273; 
  _m120[6]    = -1.*GeV2      ; _mfit20[6]    = 31.03*GeV2; 
  _r1plus[6]  = 0.122         ; _r2plus[6]    = 0.155; 
  _m12plus[6] = sqr(5.32)*GeV2; _mfit2plus[6] = -1.*GeV2; 
  _r1T[6]     = 0.111         ; _r2T[6]       = 0.175; 
  _m12T[6]    = sqr(5.32)*GeV2; _mfit2T[6]    = -1.*GeV2; 
  // initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta = -Constants::pi/9.;
}

void BallZwickyScalarFormFactor::doinit() {
  ScalarFormFactor::doinit();
  // check all the vectors have the same size
  unsigned int isize=numberOfFactors();
  if(isize!=_r10.size()||isize!=_r20.size()||isize!=_r1plus.size()||
     isize!=_r2plus.size()||isize!=_r1T.size()||
     isize!=_r2T.size()||isize!=_m120.size()||isize!=_mfit20.size()||
     isize!=_m12plus.size()||isize!=_mfit2plus.size()||
     isize!=_m12T.size()||isize!=_mfit2T.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "BallZwickyScalarFormFactor::doinit()" 
			  << Exception::abortnow;
  // output some graphs to check the answers
//   int id0,id1;
//   unsigned int iz;
//   Energy m0,m1; 
//   Energy2 q2,step(14./100.*GeV2);
//   tcPDPtr in,out;
//   Complex f0,fp,ft;
//   ofstream output("Ball.top");
//   for(unsigned int ix=0;ix<numberOfFactors();++ix) {
//     particleID(ix,id0,id1);
//     in = getParticleData(id0);
//     m0=in->mass();
//     out= getParticleData(id1);
//     m1=out->mass();
//     output << "new frame " << endl;
//     output << "newdef font duplex" << endl;
//     output << "title top \"" << in->PDGName() << " to " << out->PDGName() 
// 	   << " scalar form factors \"" << endl;
//     output << "newdef limits x 0 14. y 0 1" << endl;
//     double rt(sqrt(2.));
//     for(iz=0;iz<3;++iz) {
//       q2=ZERO;
//       for( ;q2<14.*GeV2+step;q2+=step) {
// 	ScalarScalarFormFactor(q2,ix,id0,id1,m0,m1,f0,fp);
// 	ScalarScalarSigmaFormFactor(q2,ix,id0,id1,m0,m1,ft);
// 	if(id1==111) {
// 	  if((abs(id0)%100)/10==1) {
// 	    f0*=-rt;
// 	    fp*=-rt;
// 	    ft*=-rt;
// 	  }
// 	  else {
// 	    f0*=rt;
// 	    fp*=rt;
// 	    ft*=rt;
// 	  }
// 	}
// 	else if(id1==221) {
// 	  double fact(cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.));
// 	  f0/=fact;
// 	  fp/=fact;
// 	  ft/=fact;
// 	}
// 	if(iz==0)      output << q2/GeV2 << "   " << f0.real() << "\n";
// 	else if(iz==1) output << q2/GeV2 << "   " << fp.real() << "\n";
// 	else if(iz==2) output << q2/GeV2 << "   " << ft.real() << "\n";
//       }
//       if(iz==0)      output << "join red  \n";
//       else if(iz==1) output << "join blue \n";
//       else if(iz==2) output << "join green\n";
//     }
//   }
}

void BallZwickyScalarFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _r10 << _r20 << _r1plus << _r2plus << _r1T << _r2T 
     << ounit(_m120,GeV2) << ounit(_mfit20,GeV2) 
     << ounit(_m12plus,GeV2) << ounit(_mfit2plus,GeV2) 
     << ounit(_m12T,GeV2) << ounit(_mfit2T,GeV2) << _thetaeta;
}
  
void BallZwickyScalarFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _r10 >> _r20 >> _r1plus >> _r2plus >> _r1T >> _r2T 
     >> iunit(_m120,GeV2) >> iunit(_mfit20,GeV2) 
     >> iunit(_m12plus,GeV2) >> iunit(_mfit2plus,GeV2) 
     >> iunit(_m12T,GeV2) >> iunit(_mfit2T,GeV2) >> _thetaeta;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BallZwickyScalarFormFactor,ScalarFormFactor>
describeHerwigBallZwickyScalarFormFactor("Herwig::BallZwickyScalarFormFactor", "HwFormFactors.so");

void BallZwickyScalarFormFactor::Init() {

  static ClassDocumentation<BallZwickyScalarFormFactor> documentation
    ("The BallZwickyScalarFormFactor class implements the form-factors"
     " of PRD71 014015 (2005) for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson",
     "The form factors of \\cite{Ball:2004ye} for $B\\to\\pi, K, \\eta$ were used.",
     "\\bibitem{Ball:2004ye} P.~Ball and R.~Zwicky,\n "
     "Phys.\\ Rev.\\  D {\\bf 71} (2005) 014015 [arXiv:hep-ph/0406232].\n"
     "%%CITATION = PHRVA,D71,014015;%%\n");

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_10
    ("r_10",
     "The r_1 coefficient for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_r10,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_20
    ("r_20",
     "The r_2 coefficient for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_r20,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_1plus
    ("r_1plus",
     "The r_1 coefficient for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_r1plus,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_2plus
    ("r_2plus",
     "The r_2 coefficient for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_r2plus,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_1T
    ("r_1T",
     "The r_1 coefficient for the f_T form-factor",
     &BallZwickyScalarFormFactor::_r1T,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,double> interfacer_2T
    ("r_2T",
     "The r_2 coefficient for the f_T form-factor",
     &BallZwickyScalarFormFactor::_r2T,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem120
    ("m_120",
     "The value of m_1^2 for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_m120,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit20
    ("mfit20",
     "The value of m_fit^2 for the f_0 form-factor",
     &BallZwickyScalarFormFactor::_mfit20,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem12plus
    ("m_12plus",
     "The value of m_1^2 for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_m12plus,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit2plus
    ("mfit2plus",
     "The value of m_fit^2 for the f_+ form-factor",
     &BallZwickyScalarFormFactor::_mfit2plus,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacem12T
    ("m_12T",
     "The value of m_1^2 for the f_T form-factor",
     &BallZwickyScalarFormFactor::_m12T,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyScalarFormFactor,Energy2> interfacemfit2T
    ("mfit2T",
     "The value of m_fit^2 for the f_T form-factor",
     &BallZwickyScalarFormFactor::_mfit2T,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static Parameter<BallZwickyScalarFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &BallZwickyScalarFormFactor::_thetaeta, -Constants::pi/9.,
     -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void BallZwickyScalarFormFactor::
ScalarScalarFormFactor(Energy2 q2,unsigned  int mode,
		       int, int id1, Energy, Energy,
		       Complex & f0, Complex & fp) const {
  useMe();
  // the F_0 form-factor
  if(_m120[mode]<ZERO) {
    f0=_r20[mode]/(1.-q2/_mfit20[mode]);
  }
  else if(_mfit20[mode]<ZERO) {
    f0=(_r10[mode]+_r20[mode]/(1.-q2/_m120[mode]))/(1.-q2/_m120[mode]);
  }
  else {
    f0=_r10[mode]/(1.-q2/_m120[mode])+_r20[mode]/(1.-q2/_mfit20[mode]);
  }
  // the F_1 form-factor
  if(_m12plus[mode]<ZERO) {
    fp = _r2plus[mode]/(1.-q2/_mfit2plus[mode]);
  }
  else if(_mfit2plus[mode]<ZERO) {
    fp = (_r1plus[mode]+_r2plus[mode]/(1.-q2/_m12plus[mode]))/(1.-q2/_m12plus[mode]);
  }
  else {
    fp =_r1plus[mode]/(1.-q2/_m12plus[mode])+_r2plus[mode]/(1.-q2/_mfit2plus[mode]);
  }
  if(id1==ParticleID::eta) {
    double fact(cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.));
    fp *= fact;
    f0 *= fact;
  }
}

void BallZwickyScalarFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,
							     unsigned int mode,int,
							     int id1,Energy,
							     Energy,
							     Complex & fT) const {
  useMe();
  // the F_T form-factor
  if(_m12T[mode]<ZERO) {
    fT = _r2T[mode]/(1.-q2/_mfit2T[mode]);
  }
  else if(_mfit2T[mode]<ZERO) {
    fT = (_r1T[mode]+_r2T[mode]/(1.-q2/_m12T[mode]))/(1.-q2/_m12T[mode]);
  }
  else {
    fT =_r1T[mode]/(1.-q2/_m12T[mode])+_r2T[mode]/(1.-q2/_mfit2T[mode]);
  }
  if(id1==ParticleID::eta) {
    fT *=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
  }
}

void BallZwickyScalarFormFactor::dataBaseOutput(ofstream & output,bool header,
						bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BallZwickyScalarFormFactor "
		    << name() << " \n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":r_10 " << ix << " " << _r10[ix] << "\n";
      output << "newdef " << name() << ":r_20 " << ix << " " << _r20[ix] << "\n";
      output << "newdef " << name() << ":r_1plus " << ix << " " << _r1plus[ix] << "\n";
      output << "newdef " << name() << ":r_2plus " << ix << " " << _r2plus[ix] << "\n";
      output << "newdef " << name() << ":r_1T " << ix << " " << _r1T[ix] << "\n";
      output << "newdef " << name() << ":r_2T " << ix << " " << _r2T[ix] << "\n";
      output << "newdef " << name() << ":m_120 " 
	     << ix << " " << _m120[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":mfit20 " 
	     << ix << " " << _mfit20[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":m_12plus " 
	     << ix << " " << _m12plus[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":mfit2plus " 
	     << ix << " " << _mfit2plus[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":m_12T " 
	     << ix << " " << _m12T[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":mfit2T " 
	     << ix << " " << _mfit2T[ix]/GeV2 << "\n";
    }
    else {
      output << "insert " << name() << ":r_10 " 
	     << ix << " " << _r10[ix] << "\n";
      output << "insert " << name() << ":r_20 " 
	     << ix << " " << _r20[ix] << "\n";
      output << "insert " << name() << ":r_1plus " 
	     << ix << " " << _r1plus[ix] << "\n";
      output << "insert " << name() << ":r_2plus " 
	     << ix << " " << _r2plus[ix] << "\n";
      output << "insert " << name() << ":r_1T " << ix << " " << _r1T[ix] << "\n";
      output << "insert " << name() << ":r_2T " << ix << " " << _r2T[ix] << "\n";
      output << "insert " << name() << ":m_120 " 
	     << ix << " " << _m120[ix]/GeV2 << "\n";
      output << "insert " << name() << ":mfit20 " 
	     << ix << " " << _mfit20[ix]/GeV2 << "\n";
      output << "insert " << name() << ":m_12plus " 
	     << ix << " " << _m12plus[ix]/GeV2 << "\n";
      output << "insert " << name() << ":mfit2plus " 
	     << ix << " " << _mfit2plus[ix]/GeV2 << "\n";
      output << "insert " << name() << ":m_12T " 
	     << ix << " " << _m12T[ix]/GeV2 << "\n";
      output << "insert " << name() << ":mfit2T " 
	     << ix << " " << _mfit2T[ix]/GeV2 << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./BallZwickyVectorFormFactor.cc"
// -*- C++ -*-
//
// BallZwickyVectorFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BallZwickyVectorFormFactor class.
//

#include "BallZwickyVectorFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

BallZwickyVectorFormFactor::BallZwickyVectorFormFactor() 
  : _Vr1(9), _Vr2(9), _A0r1(9), _A0r2(9), _A1r1(9), _A1r2(9), 
    _A2r1(9), _A2r2(9), _T1r1(9), _T1r2(9), _T2r1(9), _T2r2(9), 
    _T3r1(9), _T3r2(9),  _VmR2(9), _Vmfit2(9), _A0mR2(9), _A0mfit2(9), 
    _A1mR2(9), _A1mfit2(9), _A2mR2(9), _A2mfit2(9), _T1mR2(9), _T1mfit2(9), 
    _T2mR2(9), _T2mfit2(9), _T3mR2(9), _T3mfit2(9) {
  double ort(1./sqrt(2.));
  // parameters for the different form-factors
  // B to rho
  addFormFactor(-521, 113,1,-2,5,2);
  addFormFactor(-511, 113,1,-2,5,1);
  addFormFactor(-511, 213,1,-2,5,2);  
  addFormFactor(-521, 213,1,-2,5,1);
  for(unsigned int ix=0;ix<4;++ix) {
    _Vr1[ix]   = 1.045         ; _Vr2[ix]     = -0.721; 
    _VmR2[ix]  = sqr(5.32)*GeV2; _Vmfit2[ix]  = 38.34*GeV2; 
    _A0r1[ix]  = 1.527         ; _A0r2[ix]    = -1.220; 
    _A0mR2[ix] = sqr(5.28)*GeV2; _A0mfit2[ix] = 33.36*GeV2; 
    _A1r1[ix]  = 0.240         ; _A1r2[ix]    = 0.; 
    _A1mR2[ix] = -1.0*GeV2     ; _A1mfit2[ix] = 37.51*GeV2; 
    _A2r1[ix]  = 0.009         ; _A2r2[ix]    = 0.212; 
    _A2mR2[ix] = -1.0*GeV2     ; _A2mfit2[ix] = 40.82*GeV2; 
    _T1r1[ix]  = 0.897         ; _T1r2[ix]    = -0.629; 
    _T1mR2[ix] = sqr(5.32)*GeV2; _T1mfit2[ix] = 38.04*GeV2; 
    _T2r1[ix]  = 0.267         ; _T2r2[ix]    = 0.; 
    _T2mR2[ix] = -1.0*GeV2     ; _T2mfit2[ix] = 38.59*GeV2; 
    _T3r1[ix]  = 0.022         ; _T3r2[ix]    = 0.245; 
    _T3mR2[ix] = -1.0*GeV2     ; _T3mfit2[ix] = 40.88*GeV2;
  }
  for(unsigned int ix=0;ix<2;++ix) {
    double fact = ix==0 ? ort : -ort;
    _Vr1[ix]  *= fact; _Vr2[ix]  *= fact;
    _A0r1[ix] *= fact; _A0r2[ix] *= fact;
    _A1r1[ix] *= fact; _A1r2[ix] *= fact;
    _A2r1[ix] *= fact; _A2r2[ix] *= fact;
    _T1r1[ix] *= fact; _T1r2[ix] *= fact;
    _T2r1[ix] *= fact; _T2r2[ix] *= fact;
    _T3r1[ix] *= fact; _T3r2[ix] *= fact;
  }
  // parameters for the B to K   form-factors
  addFormFactor(-521,-323,1,-2,5,3);  
  addFormFactor(-511,-313,1,-2,5,3); 
  for(unsigned int ix=4;ix<6;++ix) {
    _Vr1[ix]   = 0.923         ; _Vr2[ix] = -0.511; 
    _VmR2[ix]  = sqr(5.32)*GeV2; _Vmfit2[ix] = 49.40*GeV2; 
    _A0r1[ix]  = 1.364         ; _A0r2[ix] = -0.990; 
    _A0mR2[ix] = sqr(5.28)*GeV2; _A0mfit2[ix] = 36.78*GeV2; 
    _A1r1[ix]  = 0.290         ; _A1r2[ix] = 0.0; 
    _A1mR2[ix] = -1.0*GeV2     ; _A1mfit2[ix] = 40.38*GeV2; 
    _A2r1[ix]  = -0.084        ; _A2r2[ix] = 0.342; 
    _A2mR2[ix] = -1.0*GeV2     ; _A2mfit2[ix] = 52.00*GeV2; 
    _T1r1[ix]  = 0.823         ; _T1r2[ix] = -0.491; 
    _T1mR2[ix] = sqr(5.32)*GeV2; _T1mfit2[ix] = 46.31*GeV2; 
    _T2r1[ix]  = 0.333         ; _T2r2[ix] = 0.; 
    _T2mR2[ix] = -1.0*GeV2     ; _T2mfit2[ix] = 41.41*GeV2; 
    _T3r1[ix]  = -0.036        ; _T3r2[ix] = 0.369; 
    _T3mR2[ix] = -1.0*GeV2     ; _T3mfit2[ix] = 48.10*GeV2;
  }
  // B to omega
  addFormFactor(-521,223,1,-2,5,2);
  _Vr1[6]   = 1.006*ort     ; _Vr2[6]     = -0.713*ort; 
  _VmR2[6]  = 5.32*5.32*GeV2; _Vmfit2[6]  = 37.45*GeV2; 
  _A0r1[6]  = 1.321*ort     ; _A0r2[6]    = -1.040*ort; 
  _A0mR2[6] = 5.28*5.28*GeV2; _A0mfit2[6] = 34.47*GeV2; 
  _A1r1[6]  = 0.217*ort     ; _A1r2[6]    = 0.; 
  _A1mR2[6] = -1.0*GeV2     ; _A1mfit2[6] = 37.01*GeV2; 
  _A2r1[6]  = 0.006*ort     ; _A2r2[6]    = 0.192*ort; 
  _A2mR2[6] = -1.0*GeV2     ; _A2mfit2[6] = 41.24*GeV2; 
  _T1r1[6]  = 0.865*ort     ; _T1r2[6]    = -0.622*ort; 
  _T1mR2[6] = 5.32*5.32*GeV2; _T1mfit2[6] = 37.19*GeV2; 
  _T2r1[6]  = 0.242*ort     ; _T2r2[6]    = 0.; 
  _T2mR2[6] = -1.0*GeV2     ; _T2mfit2[6] = 37.95*GeV2; 
  _T3r1[6]  = 0.023*ort     ; _T3r2[6]    = 0.219*ort; 
  _T3mR2[6] = -1.0*GeV2     ; _T3mfit2[6] = 40.87*GeV2; 
  // B_s to K*
  addFormFactor(-531,323,1,-3,5,2); 
  _Vr1[7]   = 2.351         ; _Vr2[7]     = -2.039; 
  _VmR2[7]  = sqr(5.42)*GeV2; _Vmfit2[7]  = 33.10*GeV2; 
  _A0r1[7]  = 2.813         ; _A0r2[7]    = -2.450; 
  _A0mR2[7] = sqr(5.37)*GeV2; _A0mfit2[7] = 31.58*GeV2; 
  _A1r1[7]  = 0.231         ; _A1r2[7]    = 0.; 
  _A1mR2[7] = -1.0*GeV2     ; _A1mfit2[7] = 32.94*GeV2; 
  _A2r1[7]  = -0.011        ; _A2r2[7]    = 0.192; 
  _A2mR2[7] = -1.0*GeV2     ; _A2mfit2[7] = 40.14*GeV2; 
  _T1r1[7]  = 2.047         ; _T1r2[7]    = -1.787; 
  _T1mR2[7] = sqr(5.42)*GeV2; _T1mfit2[7] = 32.83*GeV2; 
  _T2r1[7]  = 0.260         ; _T2r2[7]    = 0.; 
  _T2mR2[7] = -1.0*GeV2     ; _T2mfit2[7] = 33.01*GeV2; 
  _T3r1[7]  = 0.043         ; _T3r2[7]    = 0.217; 
  _T3mR2[7] = -1.0*GeV2     ; _T3mfit2[7] = 39.38*GeV2; 
  // B_s to phi
  addFormFactor(-531,333,1,-3,5,3);
  _Vr1[8]   = 1.484         ; _Vr2[8]     = -1.049; 
  _VmR2[8]  = sqr(5.42)*GeV2; _Vmfit2[8]  = 39.52*GeV2; 
  _A0r1[8]  = 3.310         ; _A0r2[8]    = -2.835; 
  _A0mR2[8] = sqr(5.37)*GeV2; _A0mfit2[8] = 31.57*GeV2; 
  _A1r1[8]  = 0.308         ; _A1r2[8]    = 0.; 
  _A1mR2[8] = -1.0*GeV2     ; _A1mfit2[8] = 36.54*GeV2; 
  _A2r1[8]  = -0.054        ; _A2r2[8]    = 0.288; 
  _A2mR2[8] = -1.0*GeV2     ; _A2mfit2[8] = 48.94*GeV2; 
  _T1r1[8]  = 1.303         ; _T1r2[8]    = -0.954; 
  _T1mR2[8] = sqr(5.42)*GeV2; _T1mfit2[8] = 38.28*GeV2; 
  _T2r1[8]  = 0.349         ; _T2r2[8]    = 0.; 
  _T2mR2[8] = -1.0*GeV2     ; _T2mfit2[8] = 37.21*GeV2; 
  _T3r1[8]  = 0.027         ; _T3r2[8]    = 0.322; 
  _T3mR2[8] = -1.0*GeV2     ; _T3mfit2[8] = 45.56*GeV2; 
  initialModes(numberOfFactors());
  // cut-off parameter
  _cutoff=0.01*GeV2;
}

void BallZwickyVectorFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_Vr1.size()||isize!=_Vr2.size()||isize!=_A0r1.size()||isize!=_A0r2.size()||
     isize!=_A1r1.size()||isize!=_A1r2.size()||isize!=_A2r1.size()||isize!=_A2r2.size()||
     isize!=_T1r1.size()||isize!=_T1r2.size()||isize!=_T2r1.size()||isize!=_T2r2.size()||
     isize!=_T3r1.size()||isize!=_T3r2.size()||isize!=_VmR2.size()||
     isize!=_Vmfit2.size()||isize!=_A0mR2.size()||isize!=_A0mfit2.size()||
     isize!=_A1mR2.size()||isize!=_A1mfit2.size()||isize!=_A2mR2.size()||
     isize!=_A2mfit2.size()||isize!=_T1mR2.size()||isize!=_T1mfit2.size()||
     isize!=_T2mR2.size()||isize!=_T2mfit2.size()||isize!=_T3mR2.size()||
     isize!=_T3mfit2.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "BallZwickyScalarFormFactor::doinit()" 
			  << Exception::abortnow;
  // output some graphs to check the answers
//   int id0,id1;
//   unsigned int iz;
//   Energy m0,m1;
//   Energy2 q2,step(14./100.*GeV2);
//   tcPDPtr in,out;
//   Complex A0,A1,A2,V;
//   ofstream output("Ball.top");
//   for(unsigned int ix=0;ix<numberOfFactors();++ix) {
//     particleID(ix,id0,id1);
//     in = getParticleData(id0);
//     m0=in->mass();
//     out= getParticleData(id1);
//     m1=out->mass();
//     output << "new frame " << endl;
//     output << "newdef font duplex" << endl;
//     output << "title top \"" << in->PDGName() << " to " << out->PDGName() 
// 	   << " vector form factors \"" << endl;
//     output << "newdef limits x 0 14. y 0 1" << endl;
//     double rt(sqrt(2.));
//     for(iz=0;iz<4;++iz) {
//       q2=ZERO;
//       for( ;q2<14.*GeV2+step;q2+=step) {
// 	ScalarVectorFormFactor(q2,ix,id0,id1,m0,m1,A0,A1,A2,V);
// 	if(id1==113||id1==223) {
//  	  if((abs(id0)%100)/10==1) {
// 	    A0*=-rt;
// 	    A1*=-rt;
// 	    A2*=-rt;
// 	    V *=-rt;
// 	  }
// 	  else {
// 	    A0*=rt;
// 	    A1*=rt;
// 	    A2*=rt;
// 	    V*=rt;
// 	  }
// 	}
// 	if(iz==0)      output << q2/GeV2 << "   " << A0.real() << endl;
// 	else if(iz==1) output << q2/GeV2 << "   " << A1.real() << endl;
// 	else if(iz==2) output << q2/GeV2 << "   " << A2.real() << endl;
// 	else if(iz==3) output << q2/GeV2 << "   " << V.real()  << endl;
//       }
//       if(iz==0)      output << "join red"    << endl;
//       else if(iz==1) output << "join blue"   << endl;
//       else if(iz==2) output << "join green"  << endl;
//       else if(iz==3) output << "join yellow" << endl;
//     }
//     output << "new frame " << endl;
//     output << "newdef font duplex" << endl;
//     output << "title top \"" << in->PDGName() << " to " << out->PDGName() 
// 	   << " penguin form factors\" " << endl;
//     output << "newdef limits x 0 14. y 0 1" << endl;
//     for(iz=0;iz<3;++iz) {
//       q2=ZERO;
//       for( ;q2<14.*GeV2+step;q2+=step) {
// 	ScalarVectorSigmaFormFactor(q2,ix,id0,id1,m0,m1,A0,A1,A2);
// 	if(id1==113||id1==223) {
//  	  if((abs(id0)%100)/10==1) {
// 	    A0*=-rt;
// 	    A1*=-rt;
// 	    A2*=-rt;
// 	  }
// 	  else {
// 	    A0*=rt;
// 	    A1*=rt;
// 	    A2*=rt;
// 	  }
// 	}
// 	if(iz==0)      output << q2/GeV2 << "   " << A0.real() << endl;
// 	else if(iz==1) output << q2/GeV2 << "   " << A1.real() << endl;
// 	else if(iz==2) output << q2/GeV2 << "   " << A2.real() << endl;
//       }
//       if(iz==0){output      << "join red"    << endl;}
//       else if(iz==1){output << "join blue"   << endl;}
//       else if(iz==2){output << "join green"  << endl;}
//     }
//   }
}

void BallZwickyVectorFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _Vr1 << _Vr2 << _A0r1 << _A0r2 << _A1r1 << _A1r2 << _A2r1 << _A2r2 << _T1r1
     << _T1r2 << _T2r1 << _T2r2 << _T3r1 << _T3r2 
     << ounit(_VmR2,GeV2) << ounit(_Vmfit2,GeV2) << ounit(_A0mR2,GeV2) 
     << ounit(_A0mfit2,GeV2) << ounit(_A1mR2,GeV2) << ounit(_A1mfit2,GeV2) 
     << ounit(_A2mR2,GeV2) << ounit(_A2mfit2,GeV2) << ounit(_T1mR2,GeV2) 
     << ounit(_T1mfit2,GeV2) << ounit(_T2mR2,GeV2) << ounit(_T2mfit2,GeV2) 
     << ounit(_T3mR2,GeV2) << ounit(_T3mfit2,GeV2) << ounit(_cutoff,GeV2);
}

void BallZwickyVectorFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _Vr1 >> _Vr2 >> _A0r1 >> _A0r2 >> _A1r1 >> _A1r2 >> _A2r1 >> _A2r2 >> _T1r1
     >> _T1r2 >> _T2r1 >> _T2r2 >> _T3r1 >> _T3r2 
     >> iunit(_VmR2,GeV2) >> iunit(_Vmfit2,GeV2) >> iunit(_A0mR2,GeV2)
     >> iunit(_A0mfit2,GeV2) >> iunit(_A1mR2,GeV2) >> iunit(_A1mfit2,GeV2) 
     >> iunit(_A2mR2,GeV2) >> iunit(_A2mfit2,GeV2) >> iunit(_T1mR2,GeV2) 
     >> iunit(_T1mfit2,GeV2) >> iunit(_T2mR2,GeV2) >> iunit(_T2mfit2,GeV2) 
     >> iunit(_T3mR2,GeV2) >> iunit(_T3mfit2,GeV2) >> iunit(_cutoff,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BallZwickyVectorFormFactor,ScalarFormFactor>
describeHerwigBallZwickyVectorFormFactor("Herwig::BallZwickyVectorFormFactor", "HwFormFactors.so");

void BallZwickyVectorFormFactor::Init() {

  static ClassDocumentation<BallZwickyVectorFormFactor> documentation
    ("The BallZwickyVectorFormFactor class implements the vector form"
     " factors of hep-ph/0412079 for the form-factor for the decay of a B-meson to a"
     " light pseudoscalar meson",
     "The form factors of \\cite{Ball:2004rg} for $B_{d,s}\\to\\rho,\\omega,K^*,\\phi$"
     " were used.",
     "\\bibitem{Ball:2004rg} P.~Ball and R.~Zwicky, \n"
     "Phys.\\ Rev.\\  D {\\bf 71} (2005) 014029 [arXiv:hep-ph/0412079].\n"
     "%%CITATION = PHRVA,D71,014029;%%\n");

  static ParVector<BallZwickyVectorFormFactor,double> interfaceVr_1
    ("Vr_1",
     "The r_1 coefficient for the V form-factor",
     &BallZwickyVectorFormFactor::_Vr1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceVr_2
    ("Vr_2",
     "The r_2 coefficient for the V form-factor",
     &BallZwickyVectorFormFactor::_Vr2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA0r_1
    ("A0r_1",
     "The r_1 coefficient for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA0r_2
    ("A0r_2",
     "The r_2 coefficient for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA1r_1
    ("A1r_1",
     "The r_1 coefficient for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA1r_2
    ("A1r_2",
     "The r_2 coefficient for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA2r_1
    ("A2r_1",
     "The r_1 coefficient for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceA2r_2
    ("A2r_2",
     "The r_2 coefficient for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT1r_1
    ("T1r_1",
     "The r_1 coefficient for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT1r_2
    ("T1r_2",
     "The r_2 coefficient for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT2r_1
    ("T2r_1",
     "The r_1 coefficient for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT2r_2
    ("T2r_2",
     "The r_2 coefficient for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT3r_1
    ("T3r_1",
     "The r_1 coefficient for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3r1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,double> interfaceT3r_2
    ("T3r_2",
     "The r_2 coefficient for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3r2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceVmR2
    ("VmR2",
     "The value of m_R^2 for the V form-factor",
     &BallZwickyVectorFormFactor::_VmR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceVmfit2 
    ("Vmfit2",
     "The value of m_fit^2 for the V form-factor",
     &BallZwickyVectorFormFactor::_Vmfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA0mR2
    ("A0mR2",
     "The value of m_R^2 for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA0mfit2 
    ("A0mfit2",
     "The value of m_fit^2 for the A_0 form-factor",
     &BallZwickyVectorFormFactor::_A0mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA1mR2
    ("A1mR2",
     "The value of m_R^2 for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA1mfit2 
    ("A1mfit2",
     "The value of m_fit^2 for the A_1 form-factor",
     &BallZwickyVectorFormFactor::_A1mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA2mR2
    ("A2mR2",
     "The value of m_R^2 for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceA2mfit2 
    ("A2mfit2",
     "The value of m_fit^2 for the A_2 form-factor",
     &BallZwickyVectorFormFactor::_A2mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT1mR2
    ("T1mR2",
     "The value of m_R^2 for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT1mfit2 
    ("T1mfit2",
     "The value of m_fit^2 for the T_1 form-factor",
     &BallZwickyVectorFormFactor::_T1mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT2mR2
    ("T2mR2",
     "The value of m_R^2 for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT2mfit2 
    ("T2mfit2",
     "The value of m_fit^2 for the T_2 form-factor",
     &BallZwickyVectorFormFactor::_T2mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT3mR2
    ("T3mR2",
     "The value of m_R^2 for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3mR2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static ParVector<BallZwickyVectorFormFactor,Energy2> interfaceT3mfit2 
    ("T3mfit2",
     "The value of m_fit^2 for the T_3 form-factor",
     &BallZwickyVectorFormFactor::_T3mfit2,
     GeV2, -1, 25.*GeV2, -2.*GeV2, 100.*GeV2, false, false, true);

  static Parameter<BallZwickyVectorFormFactor,Energy2> interfaceCutOff
    ("CutOff",
     "Parameter controlling the value of q^2 where we switch from the fit "
     "to a small q^2 expansion for numerical stability.",
     &BallZwickyVectorFormFactor::_cutoff, GeV2, 2.0*GeV2, ZERO, 10.0*GeV2,
     false, false, true);

}

// form-factor for scalar to vector
void BallZwickyVectorFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
							int,int,Energy, 
							Energy,
							Complex & A0,Complex & A1,
							Complex & A2,Complex & V) const {
  useMe();
  // the form-factors
  // A_0
  if(_A0mR2[mode]<ZERO) {
    A0 = (_A0r1[mode]+_A0r2[mode]/(1.-q2/_A0mfit2[mode]))/(1.-q2/_A0mfit2[mode]);
  }
  else {
    A0 = _A0r1[mode]/(1.-q2/_A0mR2[mode])+_A0r2[mode]/(1.-q2/_A0mfit2[mode]);
  }
  // A_1
  if(_A1mR2[mode]<ZERO) {
    A1 = (_A1r1[mode]+_A1r2[mode]/(1.-q2/_A1mfit2[mode]))/(1.-q2/_A1mfit2[mode]);
  }
  else {
    A1 = _A1r1[mode]/(1.-q2/_A1mR2[mode])+_A1r2[mode]/(1.-q2/_A1mfit2[mode]);
  }
  // A_2
  if(_A2mR2[mode]<ZERO) {
    A2 = (_A2r1[mode]+_A2r2[mode]/(1.-q2/_A2mfit2[mode]))/(1.-q2/_A2mfit2[mode]);
  }
  else {
    A2 = _A2r1[mode]/(1.-q2/_A2mR2[mode])+_A2r2[mode]/(1.-q2/_A2mfit2[mode]);
  }
  // V
  if(_VmR2[mode]<ZERO) {
    V = (_Vr1[mode]+_Vr2[mode]/(1.-q2/_Vmfit2[mode]))/(1.-q2/_Vmfit2[mode]);
  }
  else {
    V = _Vr1[mode]/(1.-q2/_VmR2[mode])+_Vr2[mode]/(1.-q2/_Vmfit2[mode]);
  }
}

void BallZwickyVectorFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,
							     unsigned int mode,int,
							     int,Energy m0,Energy m1,
							     Complex & T1,Complex & T2,
							     Complex & T3) const {
  useMe();
  // T_1
  if(_T1mR2[mode]<ZERO) {
    T1 = (_T1r1[mode]+_T1r2[mode]/(1.-q2/_T1mfit2[mode]))/(1.-q2/_T1mfit2[mode]);
  }
  else {
    T1 = _T1r1[mode]/(1.-q2/_T1mR2[mode])+_T1r2[mode]/(1.-q2/_T1mfit2[mode]);
  }
  // T_2
  if(_T2mR2[mode]<ZERO) {
    T2 = (_T2r1[mode]+_T2r2[mode]/(1.-q2/_T2mfit2[mode]))/(1.-q2/_T2mfit2[mode]);
  }
  else {
    T2 = _T2r1[mode]/(1.-q2/_T2mR2[mode])+_T2r2[mode]/(1.-q2/_T2mfit2[mode]);
  }
  // T_3
  if(q2>_cutoff) {
    if(_T3mR2[mode]<ZERO) {
      T3 = (_T3r1[mode]+_T3r2[mode]/(1.-q2/_T3mfit2[mode]))/(1.-q2/_T3mfit2[mode]);
    }
    else {
      T3 = _T3r1[mode]/(1.-q2/_T3mR2[mode])+_T3r2[mode]/(1.-q2/_T3mfit2[mode]);
    }
    // convert for T_3tilde to T_3
    T3 = Complex((m0*m0-m1*m1)/q2*(T3-T2));
  }
  else {
    InvEnergy2 smallT2,smallT3;
    if(_T2mR2[mode]<ZERO) {
      double a(q2/_T2mfit2[mode]);
      smallT2=1./_T2mfit2[mode]*
	(_T2r1[mode]+2.*_T2r2[mode]+a*(_T2r1[mode]+3.*_T2r2[mode]+
				       a*(_T2r1[mode]+4.*_T2r2[mode]+
					  a*(_T2r1[mode]+5.*_T2r2[mode]))));
    }
    else {
      smallT2=(_T2r1[mode]/_T2mR2[mode]+_T2r2[mode]/_T2mfit2[mode])
	+q2*(+_T2r1[mode]/_T2mR2[mode]/_T2mR2[mode]
	     +_T2r2[mode]/_T2mfit2[mode]/_T2mfit2[mode])
	+q2*q2*(+_T2r1[mode]/_T2mR2[mode]/_T2mR2[mode]/_T2mR2[mode]
		+_T2r2[mode]/_T2mfit2[mode]/_T2mfit2[mode]/_T2mfit2[mode]);
    }
    if(_T3mR2[mode]<ZERO) {
      double a(q2/_T3mfit2[mode]);
      smallT3=1./_T3mfit2[mode]*
	(_T3r1[mode]+2.*_T3r2[mode]+a*(_T3r1[mode]+3.*_T3r2[mode]+
				       a*(_T3r1[mode]+4.*_T3r2[mode]+
					  a*(_T3r1[mode]+5.*_T3r2[mode]))));
    }
    else {
      smallT3=(_T3r1[mode]/_T3mR2[mode]+_T3r2[mode]/_T3mfit2[mode])
	+q2*(+_T3r1[mode]/_T3mR2[mode]/_T3mR2[mode]
	     +_T3r2[mode]/_T3mfit2[mode]/_T3mfit2[mode])
	+q2*q2*(+_T3r1[mode]/_T3mR2[mode]/_T3mR2[mode]/_T3mR2[mode]
		+_T3r2[mode]/_T3mfit2[mode]/_T3mfit2[mode]/_T3mfit2[mode]);
    }
    T3 = (m0+m1)*(m0-m1)*(smallT3-smallT2);
  }
}

void BallZwickyVectorFormFactor::dataBaseOutput(ofstream & output,bool header,
						bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BallZwickyVectorFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":CutOff " << _cutoff/GeV2 << "\n";
  for(unsigned int ix=0;ix<_Vr1.size();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":Vr_1 "  << ix << "  " 
	     << _Vr1[ix]  << "\n";
      output << "newdef " << name() << ":Vr_2 "  << ix << "  " 
	     << _Vr2[ix]  << "\n";
      output << "newdef " << name() << ":A0r_1 "  << ix << "  " 
	     << _A0r1[ix] << "\n";
      output << "newdef " << name() << ":A0r_2 "  << ix << "  " 
	     << _A0r2[ix] << "\n";
      output << "newdef " << name() << ":A1r_1 "  << ix << "  " 
	     << _A1r1[ix] << "\n";
      output << "newdef " << name() << ":A1r_2 "  << ix << "  " 
	     << _A1r2[ix] << "\n";
      output << "newdef " << name() << ":A2r_1 "  << ix << "  " 
	     << _A2r1[ix] << "\n";
      output << "newdef " << name() << ":A2r_2 "  << ix << "  " 
	     << _A2r2[ix] << "\n";
      output << "newdef " << name() << ":T1r_1 "  << ix << "  " 
	     << _T1r1[ix] << "\n";
      output << "newdef " << name() << ":T1r_2 "  << ix << "  " 
	     << _T1r2[ix] << "\n";
      output << "newdef " << name() << ":T2r_1 "  << ix << "  " 
	     << _T2r1[ix] << "\n";
      output << "newdef " << name() << ":T2r_2 "  << ix << "  " 
	     << _T2r2[ix] << "\n";
      output << "newdef " << name() << ":T3r_1 "  << ix << "  " 
	     << _T3r1[ix] << "\n";
      output << "newdef " << name() << ":T3r_2 "  << ix << "  " 
	     << _T3r2[ix] << "\n";
      output << "newdef " << name() << ":VmR2 "  << ix 
	     << "  " << _VmR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":Vmfit2 "  << ix 
	     << "  " << _Vmfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mR2 "  << ix 
	     << "  " << _A0mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mfit2 "  << ix 
	     << "  " << _A0mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mR2 "  << ix 
	     << "  " << _A1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mfit2 "  << ix 
	     << "  " << _A1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mR2 "  << ix 
	     << "  " << _A2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mfit2 "  << ix 
	     << "  " << _A2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mR2 "  << ix 
	     << "  " << _T1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mfit2 "  << ix 
	     << "  " << _T1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mR2 "  << ix 
	     << "  " << _T2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mfit2 "  << ix 
	     << "  " << _T2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mR2 "  << ix 
	     << "  " << _T3mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mfit2 "  << ix 
	     << "  " << _T3mfit2[ix]/GeV2 << "\n";
    }
    else {
      output << "newdef " << name() << ":Vr_1 "  << ix << "  " 
	     << _Vr1[ix] << "\n";
      output << "newdef " << name() << ":Vr_2 "  << ix << "  " 
	     << _Vr2[ix] << "\n";
      output << "newdef " << name() << ":A0r_1 "  << ix << "  " 
	     << _A0r1[ix] << "\n";
      output << "newdef " << name() << ":A0r_2 "  << ix << "  " 
	     << _A0r2[ix] << "\n";
      output << "newdef " << name() << ":A1r_1 "  << ix << "  " 
	     << _A1r1[ix] << "\n";
      output << "newdef " << name() << ":A1r_2 "  << ix << "  " 
	     << _A1r2[ix] << "\n";
      output << "newdef " << name() << ":A2r_1 "  << ix << "  " 
	     << _A2r1[ix] << "\n";
      output << "newdef " << name() << ":A2r_2 "  << ix << "  " 
	     << _A2r2[ix] << "\n";
      output << "newdef " << name() << ":T1r_1 "  << ix << "  " 
	     << _T1r1[ix] << "\n";
      output << "newdef " << name() << ":T1r_2 "  << ix << "  " 
	     << _T1r2[ix] << "\n";
      output << "newdef " << name() << ":T2r_1 "  << ix << "  " 
	     << _T2r1[ix] << "\n";
      output << "newdef " << name() << ":T2r_2 "  << ix << "  " 
	     << _T2r2[ix] << "\n";
      output << "newdef " << name() << ":T3r_1 "  << ix << "  " 
	     << _T3r1[ix] << "\n";
      output << "newdef " << name() << ":T3r_2 "  << ix << "  " 
	     << _T3r2[ix] << "\n";
      output << "newdef " << name() << ":VmR2 "  << ix 
	     << "  " << _VmR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":Vmfit2 "  << ix 
	     << "  " << _Vmfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mR2 "  << ix 
	     << "  " << _A0mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A0mfit2 "  << ix 
	     << "  " << _A0mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mR2 "  << ix 
	     << "  " << _A1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A1mfit2 "  << ix 
	     << "  " << _A1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mR2 "  << ix 
	     << "  " << _A2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":A2mfit2 "  << ix 
	     << "  " << _A2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mR2 "  << ix 
	     << "  " << _T1mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T1mfit2 "  << ix 
	     << "  " << _T1mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mR2 "  << ix 
	     << "  " << _T2mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T2mfit2 "  << ix 
	     << "  " << _T2mfit2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mR2 "  << ix 
	     << "  " << _T3mR2[ix]/GeV2 << "\n";
      output << "newdef " << name() << ":T3mfit2 "  << ix 
	     << "  " << _T3mfit2[ix]/GeV2 << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./BaryonSimpleFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonSimpleFormFactor class.
//

#include "BaryonSimpleFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

BaryonSimpleFormFactor::BaryonSimpleFormFactor() {
  // axial vector coupling in beta decay
  _gA=1.25;
  // SU(3) breaking paramters
  // D/(D+F) ratio
  _alphaD=0.6;
  // eta_V parameter
  _etaV=0.970;
  // eta_A parameter
  _etaA=1.080;
  // SU(3) breaking factors for the electric dipole moment
  _rhoE=0.094;
  // SU(3) breaking factors for the magnetic dipole moment
  _rhoM=0.860;
  // the various decay modes handled by the model
  addFormFactor(2112,2212,2,2,2,1,1,2);
  addFormFactor(3222,3122,2,2,3,2,1,2);
  addFormFactor(3112,3122,2,2,3,1,1,2);
  addFormFactor(3112,3212,2,2,3,1,1,2);
  addFormFactor(3212,3222,2,2,3,2,1,2);
  addFormFactor(3312,3322,2,2,3,3,1,2);
  addFormFactor(3122,2212,2,2,2,1,3,2);
  addFormFactor(3212,2212,2,2,2,1,3,2);
  addFormFactor(3112,2112,2,2,1,1,3,2);
  addFormFactor(3312,3122,2,2,3,1,3,2);
  addFormFactor(3312,3212,2,2,3,1,3,2);
  addFormFactor(3322,3222,2,2,3,2,3,2);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void BaryonSimpleFormFactor::doinit() {
  BaryonFormFactor::doinit();
  _f1.clear();
  _f2.clear();
  _g1.clear();
  _g2.clear();
  _f1.resize(numberOfFactors());
  _f2.resize(numberOfFactors());
  _g1.resize(numberOfFactors());
  _g2.resize(numberOfFactors());
  // calculate the couplings for the different modes
  int id0,id1;
  double root23(sqrt(2./3.)),root2(sqrt(2.)),root32(sqrt(3./2.));
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    // get the particle ids for the mode
    particleID(ix,id0,id1);
    // work out the couplings
    // neutron beta decay
    if(id0==2112&&id1==2212) {
      _f1[ix] = 1. ;
      _g1[ix] = _gA;
      _f2[ix] = 3.7*_gA-1.0;
      _g2[ix] = 0.;
    }
    // sigma+ to Lambda
    else if(id0==3222&&id1==3122) {
      _f1[ix] = 0.;
      _g1[ix] = _gA*root23*_alphaD;
      _f2[ix] = 4.55*_g1[ix]-1.0*_f1[ix];
      _g2[ix] = -0.03*_g1[ix];
    }
    // sigma
    else if(id0==3112&&id1==3122) {
      _f1[ix] = 0.;
      _g1[ix] = _gA*root23*_alphaD;
      _f2[ix] = 4.55*_g1[ix]-1.0*_f1[ix];
      _g2[ix] = -0.03*_g1[ix];
    }
    else if(id0==3112&&id1==3212) {
      _f1[ix] = root2;
      _g1[ix] = _gA*root2*(1.-_alphaD);
      _f2[ix] = 4.69*_g1[ix]-_f1[ix];
      _g2[ix] = 0.;
    }
    else if(id0==3212&&id1==3222) {
      _f1[ix] = -root2;
      _g1[ix] = -_gA*root2*(1.-_alphaD);
      _f2[ix] = 4.69*_g1[ix]-_f1[ix];
      _g2[ix] = 0.;
    }
    else if(id0==3312&&id1==3322) {
      _f1[ix] = -1.;
      _g1[ix] = -_gA*(1.-2.*_alphaD);
      _f2[ix] = 5.21*_g1[ix]-_f1[ix];
      _g2[ix] = 0.;
    }
    else if(id0==3122&&id1==2212) {
      _f1[ix] = -root32*_etaV;
      _g1[ix] = -_gA*root32*_etaA*(1.-2.*_alphaD/3.);
      _f2[ix] = 4.05*_rhoM*_g1[ix]-1.01*_f1[ix];
      _g2[ix] = (4.05*_rhoE-0.09)*_g1[ix];
    }
    else if(id0==3212&&id1==2212) {
      _f1[ix] = -_etaV/root2;
      _g1[ix] = -_gA/root2*_etaA*(1.-2.*_alphaD);
      _f2[ix] = 4.2*_rhoM*_g1[ix]-1.03*_f1[ix];
      _g2[ix] = (4.2*_rhoE-0.12)*_g1[ix];
    }
    else if(id0==3112&&id1==2112) {
      _f1[ix] = -_etaV;
      _g1[ix] = -_gA*_etaA*(1.-2.*_alphaD);
      _f2[ix] = 4.2*_rhoM*_g1[ix]-1.03*_f1[ix];
      _g2[ix] = (4.2*_rhoE-0.12)*_g1[ix];
    }
    else if(id0==3312&&id1==3122) {
      _f1[ix] = root32*_etaV;
      _g1[ix] = _gA*root32*_etaA*(1.-4./3.*_alphaD);
      _f2[ix] = 4.8*_rhoM*_g1[ix]-1.01*_f1[ix];
      _g2[ix] = (4.8*_rhoE-0.08)*_g1[ix];
    }
    else if(id0==3312&&id1==3212) {
      _f1[ix] = _etaV/root2;
      _g1[ix] = _gA*_etaA/root2;
      _f2[ix] = 4.95*_rhoM*_g1[ix]-1.*_f1[ix];
      _g2[ix] = (4.95*_rhoE-0.05)*_g1[ix];
    }
    else if(id0==3322&&id1==3222) {
      _f1[ix] = _etaV;
      _g1[ix] = _gA*_etaA;
      _f2[ix] = 4.95*_rhoM*_g1[ix]-1.*_f1[ix];
      _g2[ix] = (4.95*_rhoE-0.05)*_g1[ix];
    }
    else {
      throw InitException() << "Mode not recognised in BaryonSimpleFormFactor "
			    << id0 << "  " << id1 
			    << Exception::abortnow;
    }
  }
}

void BaryonSimpleFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _gA << _alphaD << _etaV << _etaA << _rhoE << _rhoM 
     << _f1 << _f2 << _g1 << _g2;
}

void BaryonSimpleFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _gA >> _alphaD >> _etaV >> _etaA >> _rhoE >> _rhoM 
     >> _f1 >> _f2 >> _g1 >> _g2;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BaryonSimpleFormFactor,BaryonFormFactor>
describeHerwigBaryonSimpleFormFactor("Herwig::BaryonSimpleFormFactor",
				     "HwFormFactors.so");

void BaryonSimpleFormFactor::Init() {

  static ClassDocumentation<BaryonSimpleFormFactor> documentation
    ("The BaryonSimpleFormFactor class implements the"
     " quark model calculation of the form-factors from PRD25, 206",
     "The BaryonSimpleFormFactor class which implements the results"
     " of \\cite{Donoghue:1981uk}was used for the weak decays of the light baryons.",
     "\\bibitem{Donoghue:1981uk}\n"
     "J.~F.~Donoghue and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 25} (1982) 206.\n"
     "%%CITATION = PHRVA,D25,206;%%\n");

  static Parameter<BaryonSimpleFormFactor,double> interfaceg_A
    ("g_A",
     "The axial-vector coupling in neutron beta decay.",
     &BaryonSimpleFormFactor::_gA, 1.25, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacealpha_D
    ("alpha_D",
     "SU(3) breaking parameter which is the ratio D/(D+F). ",
     &BaryonSimpleFormFactor::_alphaD, 0.6, 0.0, 1.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfaceeta_V
    ("eta_V",
     "The eta_V SU(3) breaking parameter",
     &BaryonSimpleFormFactor::_etaV, .97, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfaceeta_A
    ("eta_A",
     "The eta_A SU(3) breaking parameter",
     &BaryonSimpleFormFactor::_etaA, 1.08, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacerho_E
    ("rho_E",
     "The SU(3) breaking parameter for the electric dipole moment.",
     &BaryonSimpleFormFactor::_rhoE, 0.094, 0.0, 10.0,
     false, false, true);

  static Parameter<BaryonSimpleFormFactor,double> interfacerho_M
    ("rho_M",
     "The SU(3) breaking parameter for the magentic dipole moment.",
     &BaryonSimpleFormFactor::_rhoM, 0.860, 0.0, 10.0,
     false, false, true);
}

void BaryonSimpleFormFactor::
SpinHalfSpinHalfFormFactor(Energy2,int iloc, int,int,Energy,Energy,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo,
			   Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  f1v =  _f1[iloc];
  f1a = -_g1[iloc];
  f2v = -_f2[iloc];
  f2a =  _g2[iloc];
  f3v = 0.;
  f3a = 0.;
}

void BaryonSimpleFormFactor::dataBaseOutput(ofstream& output,bool header,
					    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BaryonSimpleFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":g_A " <<  _gA << " \n";
  output << "newdef " << name() << ":alpha_D " << _alphaD  << " \n";
  output << "newdef " << name() << ":eta_V " << _etaV  << " \n";
  output << "newdef " << name() << ":eta_A " <<  _etaA << " \n";
  output << "newdef " << name() << ":rho_E " << _rhoE  << " \n";
  output << "newdef " << name() << ":rho_M " << _rhoM  << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./BaryonThreeQuarkModelFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonThreeQuarkModelFormFactor class.
//

#include "BaryonThreeQuarkModelFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/GaussianIntegrator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

// function for the integral of the partial width
namespace {
using namespace Herwig;

struct BaryonCMatrixElement {

typedef ThePEG::Ptr<BaryonThreeQuarkModelFormFactor>::pointer 
BaryonThreeQuarkModelFormFactorPtr;

BaryonCMatrixElement(BaryonThreeQuarkModelFormFactorPtr in,
		     Energy m0, Energy m1, int type, int mass,
		     int id0,int id1) {
  _formFactor=in;
  _m0=m0;
  _m1=m1;
  _type=type;
  _mass=mass;
  _id0=id0;
  _id1=id1;
}
  
// calculate the integrand  
Energy operator() (double x) const {
  return _formFactor->widthIntegrand(x,_m0,_m1,_type,_mass,_id0,_id1);
} 
  /** Argument type for GaussianIntegrator */
  typedef double ArgType;
  /** Return type for GaussianIntegrator */
  typedef Energy ValType;

private:
  // private variables
  BaryonThreeQuarkModelFormFactorPtr _formFactor;
  Energy _m0,_m1;
  int _type,_mass;
  int _id0,_id1;
};

}

BaryonThreeQuarkModelFormFactor::BaryonThreeQuarkModelFormFactor() 
  : _initialize(false), _order(50),_mlight(420*MeV),_mstrange(570*MeV),
    _LambdaQ(2.5*GeV),_Lambdaqq(0.71*GeV),_Lambdasq(850*MeV),_Lambdass(1.0*GeV) {
  // modes handled by this form factor
  // lambda_b
  addFormFactor(5122,4122,2,2,1,2,5,4);
  // xi_b
  addFormFactor(5232,4232,2,2,2,3,5,4);
  addFormFactor(5132,4132,2,2,1,3,5,4);
  // sigma_b
  addFormFactor(5112,4112,2,2,1,1,5,4);
  addFormFactor(5212,4212,2,2,2,1,5,4);
  addFormFactor(5222,4222,2,2,2,2,5,4);
  // omega_b
  addFormFactor(5332,4332,2,2,3,3,5,4);
  // sigma_b-> sigma_c*
  addFormFactor(5112,4114,2,4,1,1,5,4);
  addFormFactor(5212,4214,2,4,2,1,5,4);
  addFormFactor(5222,4224,2,4,2,2,5,4);
  // omega_b -> omega_c*
  addFormFactor(5332,4334,2,4,3,3,5,4);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void BaryonThreeQuarkModelFormFactor::doinit() {
  BaryonFormFactor::doinit();
  // initialization in needed
  if(_initialize) {
    GaussianIntegrator integrator;
    _C0.clear();_C1.clear();_C2.clear();
    double pre(0.),root(2.*sqrt(6.));
    double gamma1(1),gamma2(1),gamma3(sqrt(acos(-1.)));
    unsigned int ix,iy;
    for(iy=0;iy<2;++iy) {
      _mu2 = iy==0 ? sqr(_mlight  /_LambdaQ) : sqr(_mstrange/_LambdaQ);
      for(ix=0;ix<=_order;++ix) {
	if(ix>0)    gamma1*=ix;
	if(ix%2==1) { gamma2*=(ix+1)/2.0; gamma3*=ix/2.0; }
	if(ix%2==0) pre=pow(root,double(ix))/12.*gamma2/gamma1;
	else        pre=pow(root,double(ix))/12.*gamma3/gamma1;
	// for the xi_0 function
	_a=_mu2;_b=2.;_N=ix;
	_C0.push_back(pre*integrator.value(*this,0.,1.));
	// for the xi_1 function
	_a=_mu2;_b=1.;
	_C1.push_back(pre*integrator.value(*this,0.,1.));
	// for the xi_2 function
	_a=0.;
	_b=0.;
	_C2.push_back(pre*integrator.value(*this,0.,1.));
      }
    }
    generator()->log() << "Checking results of BaryonThreeQuarkModelFormFactor"
		       << "vs results from the original paper\n";
    // first matrix element
    Energy m0=getParticleData(5122)->mass();
    Energy m1=getParticleData(4122)->mass();
    double omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    BaryonCMatrixElement int1(this,m0,m1,1,0,5122,4122);
    Energy width = integrator.value(int1,1.,omegamax); 
    generator()->log() << "Lambda_b0->Lambda_c+ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // second matrix element
    m0=getParticleData(5222)->mass();
    m1=getParticleData(4222)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,2,0,5222,4222);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Sigma_b+->Sigma_c++ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // third matrix element
    m0=getParticleData(5232)->mass();
    m1=getParticleData(4232)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,1,1,5232,4232);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Xi_b0->Xi_c+ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // fourth matrix element
    m0=getParticleData(5332)->mass();
    m1=getParticleData(4332)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,2,2,5332,4332);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Omega_b-->Omega_c0 decay" 
		       << width/6.582119E-22/MeV << "\n";
    // fifth matrix element
    m0=getParticleData(5222)->mass();
    m1=getParticleData(4224)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,3,0,5222,4224);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Sigma_vb+->Sigma_c*++ decay" 
		       << width/6.582119E-22/MeV << "\n";
    // fourth matrix element
    m0=getParticleData(5332)->mass();
    m1=getParticleData(4334)->mass();
    omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
    int1 = BaryonCMatrixElement(this,m0,m1,3,2,5332,4334);
    width = integrator.value(int1,1.,omegamax);
    generator()->log() << "Omega_b-->Omega_c*0 decay" 
		       << width/6.582119E-22/MeV << "\n";
    // output some plots for testing
    double lambdabar = -999.999;
    ofstream output("ThreeQuark.top");
    output << "newdef font duplex" << endl;
    output << "title top \"Figure 3 from paper \"" << endl;
    output << "newdef limits x 1 1.44 y 0.5 1" << endl;
    for(unsigned int ix=0;ix<5;++ix) {
      double omegamin(1.),omegamax(1.44),
	step((omegamax-omegamin)/100.),omega(1.),xi;
      unsigned int ioff(0);
      if(ix==0){lambdabar=600*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=650*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=710*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=750*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=800*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step) {
	vector<double> phi(phiFunction(omega));
	double power(1.),numer(0.),denom(0.);
	for(unsigned int iy=0;iy<=_order;++iy) {
	  numer+=phi[iy]*power*_C0[iy+ioff];
	  denom+=power*_C0[iy+ioff];
	  power*=lambdabar;
	}
	xi=numer/denom;
	output << omega << "   " << xi << endl; 
      }
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
    output << "new frame " << endl;
    output << "newdef font duplex" << endl;
    output << "title top \"Figure 6 from paper \"" << endl;
    output << "newdef limits x 1 1.4 y 0.5 1" << endl;
    for(unsigned int ix=0;ix<5;++ix) {
      double omegamin(1.),omegamax(1.45),step((omegamax-omegamin)/100.),omega(1.);
      unsigned int ioff(0);
      if(ix==0){lambdabar=600*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=650*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=710*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=750*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=800*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step) {
	vector<double> phi(phiFunction(omega));
	double power(1.),numer[2]={0.,0.},denom(0.);
	for(unsigned int iy=0;iy<=_order;++iy) {
	  numer[0]+=phi[iy]*power*_C1[iy+ioff];
	  denom+=power*_C1[iy+ioff];
	  numer[1]+=power*_C2[iy+ioff]*(phi[iy]-phi[iy+2]);
	  power*=lambdabar;
	}
	numer[1]/=(omega-1.);
	double xi1(numer[0]/denom);
	output << omega << "   " << xi1 << endl; 
      }
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
    output << "new frame " << endl;
    output << "newdef font duplex" << endl;
    output << "title top \"Figure 7 from paper \"" << endl;
    output << "newdef limits x 1 1.33 y 0.4 1" << endl;
    for(unsigned int ix=0;ix<5;++ix) {
      double omegamin(1.),omegamax(1.45),step((omegamax-omegamin)/100.),omega(1.);
      unsigned int ioff(_order+1);
      if(ix==0){lambdabar=800*MeV/_LambdaQ;}
      else if(ix==1){lambdabar=900*MeV/_LambdaQ;}
      else if(ix==2){lambdabar=1000*MeV/_LambdaQ;}
      else if(ix==3){lambdabar=1050*MeV/_LambdaQ;}
      else if(ix==4){lambdabar=1100*MeV/_LambdaQ;}
      for(;omega<omegamax;omega+=step) {
	vector<double> phi(phiFunction(omega));
	double power(1.),numer[2]={0.,0.},denom(0.);
	for(unsigned int iy=0;iy<=_order;++iy) {
	  numer[0]+=phi[iy]*power*_C1[iy+ioff];
	  denom+=power*_C1[iy+ioff];
	  numer[1]+=power*_C2[iy+ioff]*(phi[iy]-phi[iy+2]);
	  power*=lambdabar;
	}
	numer[1]/=(omega-1.);
	double xi1(numer[0]/denom);
	output << omega << "   " << xi1 << endl; 
      }
      if(ix==0){output << "join red" << endl;}
      else if(ix==1){output << "join green" << endl;}
      else if(ix==2){output << "join blue" << endl;}
      else if(ix==3){output << "join cyan" << endl;}
      else if(ix==4){output << "join magenta" << endl;}
    }
  }
}

void BaryonThreeQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _initialize << _order << ounit(_mlight,MeV) << ounit(_mstrange,MeV) 
     << ounit(_LambdaQ,MeV) << ounit(_Lambdaqq,MeV) 
     << ounit(_Lambdasq,MeV) << ounit(_Lambdass,MeV) << _C0 << _C1 << _C2;
}

void BaryonThreeQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _initialize >> _order >> iunit(_mlight,MeV) >> iunit(_mstrange,MeV) 
     >> iunit(_LambdaQ,MeV) >> iunit(_Lambdaqq,MeV)
     >> iunit(_Lambdasq,MeV) >> iunit(_Lambdass,MeV) >> _C0 >> _C1 >> _C2;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BaryonThreeQuarkModelFormFactor,BaryonFormFactor>
describeHerwigBaryonThreeQuarkModelFormFactor("Herwig::BaryonThreeQuarkModelFormFactor", "HwFormFactors.so");

void BaryonThreeQuarkModelFormFactor::Init() {

  static ClassDocumentation<BaryonThreeQuarkModelFormFactor> documentation
    ("The BaryonThreeQuarkModelFormFactor class implements"
     " the form-factors for semi-leptonic decay of baryon containing a"
     " heavy quark from PRD56, 348.",
     "The form factors from \\cite{Ivanov:1996fj} were used.",
     "%\\cite{Ivanov:1996fj}\n"
     "\\bibitem{Ivanov:1996fj}\n"
     "  M.~A.~Ivanov, V.~E.~Lyubovitskij, J.~G.~Korner and P.~Kroll,\n"
     "  ``Heavy baryon transitions in a relativistic three-quark model,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 56} (1997) 348\n"
     "  [arXiv:hep-ph/9612463].\n"
     "  %%CITATION = PHRVA,D56,348;%%\n"
     );

  static Parameter<BaryonThreeQuarkModelFormFactor,unsigned int> interfaceOrder
    ("Order",
     "The order of terms to include in the series expansion of the form-factor.",
     &BaryonThreeQuarkModelFormFactor::_order, 10, 0, 1000,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLightMass
    ("LightMass",
     "The mass of the light quark",
     &BaryonThreeQuarkModelFormFactor::_mlight, GeV, .42*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark",
     &BaryonThreeQuarkModelFormFactor::_mstrange, GeV, .57*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdaQ
    ("LambdaQ",
     "Heavy Baryon Size Parameter",
     &BaryonThreeQuarkModelFormFactor::_LambdaQ, GeV, 2.5*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdaqq
    ("Lambdaqq",
     "The size parameter for light quarks",
     &BaryonThreeQuarkModelFormFactor::_Lambdaqq, GeV, 0.71*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdasq
    ("Lambdasq",
     "The size parameter for one strange quark",
     &BaryonThreeQuarkModelFormFactor::_Lambdasq, GeV, 0.85*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static Parameter<BaryonThreeQuarkModelFormFactor,Energy> interfaceLambdass
    ("Lambdass",
     "The size parameter with two strange quarks.",
     &BaryonThreeQuarkModelFormFactor::_Lambdass, GeV, 1.0*GeV, ZERO, 2.0*GeV,
     false, false, true);

  static ParVector<BaryonThreeQuarkModelFormFactor,double> interfaceC0
    ("C0",
     "The coefficient of the expansion for xi_0.",
     &BaryonThreeQuarkModelFormFactor::_C0,
     0, 0, 0, -1E20, 1E20, false, false, true);

  static ParVector<BaryonThreeQuarkModelFormFactor,double> interfaceC1
    ("C1",
     "The coefficient of the expansion for xi_1.",
     &BaryonThreeQuarkModelFormFactor::_C1,
     0, 0, 0, -1E20, 1E20, false, false, true);

  static ParVector<BaryonThreeQuarkModelFormFactor,double> interfaceC2
    ("C2",
     "The coefficient of the expansion for xi_2.",
     &BaryonThreeQuarkModelFormFactor::_C2,
     0, 0, 0, -1E20, 1E20, false, false, true);


  static Switch<BaryonThreeQuarkModelFormFactor,bool> interfaceInitialize
    ("Initialize",
     "Initialize the coefficient for the expansion of the form-factor",
     &BaryonThreeQuarkModelFormFactor::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialize
    (interfaceInitialize,
     "Yes",
     "Perform the initialize",
     true);
  static SwitchOption interfaceInitializeNoInitialize
    (interfaceInitialize,
     "No",
     "No initialization",
     false);
}

void BaryonThreeQuarkModelFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int,int id0,int id1,Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo ,
			   Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  // this model is based on heavy quark theory
  // therefore most of the factors are zero
  Complex g1v(0.),g1a(0.),g2v(0.),g2a(0.),g3a(0.),g3v(0.);
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4332) {
    lambdabar=_Lambdass/_LambdaQ;
    ioff=_order+1;
  }
  else if(abs(id1)==4232||abs(id1)==4132) {
    lambdabar=_Lambdasq/_LambdaQ;
    ioff=_order+1;
  }
  else {
    lambdabar=_Lambdaqq/_LambdaQ;
  }
  // the omega value
  double omega = 0.5/m0/m1*(m0*m0+m1*m1-q2);
  // calculate the form factor
  // the phi functions
  vector<double> phi(phiFunction(omega));
  // use the xi0 functions
  if(abs(id0)==5122||abs(id0)==5232||abs(id0)==5132) {
    double power(1.),numer(0.),denom(0.);
    for(unsigned int ix=0;ix<=_order;++ix) {
      numer+=phi[ix]*power*_C0[ix+ioff];
      denom+=power*_C0[ix+ioff];
      power*=lambdabar;
    }
    g1v=numer/denom;
    g1a=g1v;
  }
  else {
    double power(1.),numer[2]={0.,0.},denom(0.);
    for(unsigned int ix=0;ix<=_order;++ix) {
      numer[0]+=phi[ix]*power*_C1[ix+ioff];
      denom+=power*_C1[ix+ioff];
      numer[1]+=power*_C2[ix+ioff]*(phi[ix]-phi[ix+2]);
      power*=lambdabar;
    }
    numer[1]/=(omega-1.);
    double xi1(numer[0]/denom),xi2(numer[1]/denom);
    // the couplings in the velocity form
    g1v = -(omega*xi1-(sqr(omega)-1.)*xi2)/3.;
    g1a = g1v;
    g2v = 2./3.*(xi1-(omega-1.)*xi2);
    g3v = g2v;
    g2a = 2./3.*(xi1-(omega+1.)*xi2);
    g3a =-g2a;
  }
  // convert to our form
  f1v = g1v+Complex(0.5*(m0+m1)*(g2v/m0+g3v/m1));
  f1a =-g1a+Complex(0.5*(m0-m1)*(g2a/m0+g3a/m1));
  f2v = Complex(0.5*(m0+m1)*( g2v/m0+g3v/m1));
  f3v = Complex(0.5*(m0+m1)*( g2v/m0-g3v/m1));
  f2a =-Complex(0.5*(m0+m1)*( g2a/m0+g3a/m1));
  f3a = Complex(0.5*(m0+m1)*(-g2a/m0+g3a/m1));
}

void  BaryonThreeQuarkModelFormFactor::
SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int,int,int id1,Energy m0,
				Energy m1, Complex & f1v,Complex & f2v,
				Complex & f3v,Complex & f4v,Complex & f1a,
				Complex & f2a,Complex & f3a,Complex & f4a,
				FlavourInfo ,
				Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  // work out which light quark constant to use
  double lambdabar;unsigned int ioff(0);
  if(abs(id1)==4334) {
    lambdabar=_Lambdass/_LambdaQ;
    ioff=_order+1;
  }
  else if(abs(id1)==4234||abs(id1)==4134||abs(id1)==3324) {
    lambdabar=_Lambdasq/_LambdaQ;
    ioff=_order+1;
  }
  else {
    lambdabar=_Lambdaqq/_LambdaQ;
  }
  // the omega value
  double omega=0.5/m0/m1*(m0*m0+m1*m1-q2);
  // calculate the form factor
  // the phi functions
  vector<double> phi=phiFunction(omega);
  double power(1.),numer[2]={0.,0.},denom(0.);
  for(unsigned int ix=0;ix<=_order;++ix) {
    numer[0]+=phi[ix]*power*_C1[ix+ioff];
    denom+=power*_C1[ix+ioff];
    numer[1]+=power*_C2[ix+ioff]*(phi[ix]-phi[ix+2]);
    power*=lambdabar;
  }
  numer[1]/=(omega-1.);
  double xi1=numer[0]/denom;
  double xi2=numer[1]/denom;
  // couplings in the velocity form
  Complex g1v,g2v,g3v,g4v,g1a,g2a,g3a,g4a;
  double orr(1./sqrt(3.));
  Energy msum(m0+m1);
  Energy2 msum2(msum*msum);
  g1v = 2.*orr*xi1;
  g1a = -g1v;
  g2v = orr*(xi1-(omega-1)*xi2);
  g2a = orr*(xi1-(omega+1)*xi2);
  g3v = 0.;
  g3a = 0.;
  g4v = -2.*orr*xi2;
  g4a = -g4v;
  // convert to our form
  f1v = g1v;
  f1a =-g1a;
  f2v = g2v*msum/m0;
  f2a =-g2a*msum/m0;
  f3v = Complex(msum2/m0*(g3v/m0+g4v/m1));
  f3a =-Complex(msum2/m0*(g3a/m0+g4a/m1));
  f4v = Complex(msum2/m0/m0*g3v);
  f4a =-Complex(msum2/m0/m0*g3a);
}

void BaryonThreeQuarkModelFormFactor::dataBaseOutput(ofstream & output,bool header,
						     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BaryonThreeQuarkModelFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":Order       " << _order        << " \n";
  output << "newdef " << name() << ":LightMass   " << _mlight/GeV   << " \n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "newdef " << name() << ":LambdaQ     " << _LambdaQ/GeV  << " \n";
  output << "newdef " << name() << ":Lambdaqq    " << _Lambdaqq/GeV << " \n";
  output << "newdef " << name() << ":Lambdasq    " << _Lambdasq/GeV << " \n";
  output << "newdef " << name() << ":Lambdass    " << _Lambdass/GeV << " \n";
  // the number of terms to include in the sum for the form-factors
  for(unsigned int ix=0;ix<_C0.size();++ix)
    output << "insert " << name() << ":C0 " << ix << "   " << _C0[ix] << " \n";
  for(unsigned int ix=0;ix<_C1.size();++ix)
    output << "insert " << name() << ":C1 " << ix << "   " << _C1[ix] << " \n";
  for(unsigned int ix=0;ix<_C2.size();++ix)
    output << "insert " << name() << ":C2 " << ix << "   " << _C2[ix] << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}

double BaryonThreeQuarkModelFormFactor::operator ()(double x) const {
  // convert the integration variable
  double y=(1.-x)/x;
  double output =exp(-24.*_mu2*y)*y;
  // the integrals
  double I,Im2;
  SN(y,_N,Im2,I);
  double Nfact=0.5*_N+1.0;
  output *= (_a+Nfact/24./(1.+y))*I+1./12./(1.+y)*(_b-Nfact)*Im2;
  return output;
}

void BaryonThreeQuarkModelFormFactor::SN(double y, int N, double & SNm2,
					 double & SN) const {
  // special cases for the low lying values
  if(N==0) {
    double root=sqrt((1.+y)*(3.+y));
    SN   = 0.5/y*sqrt((1.+y)/(3.+y))*log((root+y)/(root-y));
    SNm2 = 0.5/(y+3.)*(SN+(1.+y)/(3.+4.*y)); 
  }
  else if(N==1) {
    SN   = sqrt(1.+y)/y*asin(y/sqrt(1.+y)/sqrt(3.+y));
    SNm2 = 1./(3.+y)*sqrt((1.+y)/(3.+4.*y)); 
  }
  else if(N==2) {
    double root=sqrt((1.+y)*(3.+y));
    SN   = 1.;
    SNm2 = 0.5/y*sqrt((1.+y)/(3.+y))*log((root+y)/(root-y));
  }
  // the general case
  else {
    int ix; double root;
    if(N%2==0) {
      SN=1.;
      ix=2;
      root=1.;
    }
    else {
      SN = sqrt(1.+y)/y*asin(y/sqrt(1.+y)/sqrt(3.+y));
      ix=1;
      root=sqrt((1.+y)/(3.+4.*y));
    }
    do {
      SNm2=SN;
      ix+=2;
      root*=(3.+4.*y)/(1.+y);
      SN=1.0/(ix-1.0)*(root+(ix-2.0)*(y+3.)*SNm2);
    }
    while(ix<N);
  }
}

// return the phi_N functions calculated using recursion
vector<double> BaryonThreeQuarkModelFormFactor::phiFunction(double omega) {
  vector<double> output;
  double root(sqrt(omega*omega-1.));
  output.push_back(1./root*log(omega+root));
  if(omega<1.00001) output.back()=1.;
  if(_order>0) output.push_back(2./(omega+1.));
  if(_order<2) return output;
  for(unsigned int ix=2;ix<=_order+2;++ix) {
    output.push_back(2./ix/(omega+1.)*(1.+(ix-1)*output[ix-2]));
  }
  return output;
}

Energy BaryonThreeQuarkModelFormFactor::widthIntegrand(double omega,Energy m0,
						       Energy m1, int type,
						       int ,int id0,
						       int id1) {
  // prefactors
  double omegamax=0.5*(m0*m0+m1*m1)/m0/m1;
  double pi=acos(-1.);
  InvEnergy kw=sqr(generator()->standardModel()->fermiConstant())
    /8./pi/pi/pi*m1*m1*m1/6.*(omegamax-omega)*sqrt(omega*omega-1.);
  Energy2 q2 = sqr(m0)+sqr(m1)-2.*m0*m1*omega;
  if(type<=2) {
    Complex f1v,f2v,f3v,f1a,f2a,f3a;
    SpinHalfSpinHalfFormFactor(q2,0,id0,id1,m0,m1,f1v,f2v,f3v,f1a,f2a,f3a,FlavourInfo());
    Complex left  =f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
    Complex right =f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
    double g1v = 0.5*( left+right).real();
    double g2v = m0*(f2v+f3v).real()/(m0+m1);
    double g3v = m1*(f2v-f3v).real()/(m0+m1);
    double g1a = -0.5*(+right-left).real();
    double g2a = -m0*(f2a+f3a).real()/(m0+m1);
    double g3a = -m1*(f2a-f3a).real()/(m0+m1);
    Energy Hpp = -2.*sqrt(m0*m1*(omega-1.))*g1v+2.*sqrt(m0*m1*(omega+1))*g1a;
    Energy Hmm = -2.*sqrt(m0*m1*(omega-1.))*g1v-2.*sqrt(m0*m1*(omega+1))*g1a;
    Energy Hp0 = 
      (sqrt(2.*m0*m1*(omega-1))*((m0+m1)*g1v+m1*(omega+1)*g2v+m0*(omega+1)*g3v)-
       sqrt(2.*m0*m1*(omega+1))*((m0-m1)*g1a-m1*(omega-1)*g2a-m0*(omega-1)*g3a))/sqrt(q2);
    Energy Hm0 = 
      (sqrt(2.*m0*m1*(omega-1))*((m0+m1)*g1v+m1*(omega+1)*g2v+m0*(omega+1)*g3v)+
       sqrt(2.*m0*m1*(omega+1))*((m0-m1)*g1a-m1*(omega-1)*g2a-m0*(omega-1)*g3a))/sqrt(q2);
    return kw*sqr(0.04)*(sqr(Hpp)+sqr(Hmm)+sqr(Hp0)+sqr(Hm0));
  }
  else {
    Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
    double  g1v,g2v,g3v,g4v,g1a,g2a,g3a,g4a;
    SpinHalfSpinThreeHalfFormFactor(q2,0,id0,id1,m0,m1,
				    f1v,f2v,f3v,f4v,
				    f1a,f2a,f3a,f4a,FlavourInfo());
    g1v = f1v.real();
    g1a = -f1a.real();
    g2v = m0/(m0+m1)*f2v.real();
    g2a =-m0/(m0+m1)*f2a.real();
    g3v = sqr(m0/(m0+m1))*f4v.real();
    g3a =-sqr(m0/(m0+m1))*f4a.real();
    g4v = m0*m1/sqr(m0+m1)*(f3v.real()-f4v.real());
    g4a =-m0*m1/sqr(m0+m1)*(f3a.real()-f4a.real());
    Energy HppC = 
      +sqrt(2./3.)*sqrt(m0*m1*(omega-1))*(g1v-2.*(omega+1)*g2v)
      -sqrt(2./3.)*sqrt(m0*m1*(omega+1))*(g1a-2.*(omega-1)*g2a);
    Energy HmmC = 
      -sqrt(2./3.)*sqrt(m0*m1*(omega-1))*(g1v-2.*(omega+1)*g2v)
      -sqrt(2./3.)*sqrt(m0*m1*(omega+1))*(g1a-2.*(omega-1)*g2a);
    Energy HppbC = 
      -sqrt(2.*m0*m1*(omega-1))*g1v
      -sqrt(2.*m0*m1*(omega+1))*g1a;
    Energy HmmbC = 
      sqrt(2.*m0*m1*(omega-1))*g1v
      -sqrt(2.*m0*m1*(omega+1))*g1a;
    Energy Hp0C = (
		   -2./sqrt(3.)*sqrt(m0*m1*(omega-1))*
		   ((m0*omega-m1)*g1v-(m0-m1)*(omega+1)*g2v+m1*(sqr(omega)-1)*g3v+m0*(sqr(omega)-1)*g4v)
		   -2./sqrt(3.)*sqrt(m0*m1*(omega+1))*
		   ((m0*omega-m1)*g1a+(m0+m1)*(omega-1)*g2a+m1*(sqr(omega)-1)*g3a+m0*(sqr(omega)-1)*g4a))/sqrt(q2);
    Energy Hm0C = (
		   +2./sqrt(3.)*sqrt(m0*m1*(omega-1))*
		   ((m0*omega-m1)*g1v-(m0-m1)*(omega+1)*g2v+m1*(sqr(omega)-1)*g3v+m0*(sqr(omega)-1)*g4v)
		   -2./sqrt(3.)*sqrt(m0*m1*(omega+1))*
		   ((m0*omega-m1)*g1a+(m0+m1)*(omega-1)*g2a+m1*(sqr(omega)-1)*g3a+m0*(sqr(omega)-1)*g4a))/sqrt(q2);
    return kw*sqr(0.04)*(sqr(HppC)+sqr(HmmC)+sqr(HppbC)+sqr(HmmbC)+sqr(Hp0C)+sqr(Hm0C));
  }
}
#line 1 "./ChengHeavyBaryonFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ChengHeavyBaryonFormFactor class.
//

#include "ChengHeavyBaryonFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h" 
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ChengHeavyBaryonFormFactor::ChengHeavyBaryonFormFactor() {
  // consituent quark masses
  _mu = 338*MeV;
  _md = 322*MeV;
  _ms = 510*MeV;
  _mc = 1.6*GeV;
  _mb = 5.0*GeV;
  // masses for the q^2 dependence
  _mVbc = 6.34*GeV;
  _mVbs = 5.42*GeV;
  _mVcs = 2.11*GeV;
  _mVbd = 5.32*GeV;
  _mVcu = 2.01*GeV;
  _mAbc = 6.73*GeV;
  _mAbs = 5.86*GeV;
  _mAcs = 2.54*GeV;
  _mAbd = 5.71*GeV;
  _mAcu = 2.42*GeV;
  double one3(1./sqrt(3.)),one2(1./sqrt(2.));
  // lambda_b to lambda_c
  addFormFactor(5122,4122,2,2,1,2,5,4);_Nfi.push_back(1.     );_eta.push_back(1.);
  // lambda_b to lambda
  addFormFactor(5122,3122,2,2,1,2,5,3);_Nfi.push_back(one3   );_eta.push_back(1.);
  // lambda_b to n
  addFormFactor(5122,2112,2,2,1,2,5,1);_Nfi.push_back(one2   );_eta.push_back(1.);
  // xi_b to xi_c
  addFormFactor(5232,4232,2,2,2,3,5,4);_Nfi.push_back(1.     );_eta.push_back(1.);
  addFormFactor(5132,4132,2,2,1,3,5,4);_Nfi.push_back(1.     );_eta.push_back(1.);
  // xi_b to xi
  addFormFactor(5232,3322,2,2,2,3,5,3);_Nfi.push_back(one2   );_eta.push_back(1.);
  addFormFactor(5132,3312,2,2,1,3,5,3);_Nfi.push_back(one2   );_eta.push_back(1.);
  // xi_b to sigma
  addFormFactor(5232,3212,2,2,2,3,5,1);_Nfi.push_back(0.5    );_eta.push_back(1.);
  addFormFactor(5132,3112,2,2,1,3,5,1);_Nfi.push_back(0.5    );_eta.push_back(1.);
  // xi_b to lambda
  addFormFactor(5232,3122,2,2,2,3,5,1);_Nfi.push_back(one3/2.);_eta.push_back(1.);
  // omega_b to omega_c
  addFormFactor(5332,4332,2,2,3,3,5,4);_Nfi.push_back(1.     );_eta.push_back(-1./3.);
  // omega_b to xi
  addFormFactor(5332,3312,2,2,3,3,5,1);_Nfi.push_back(one3   );_eta.push_back(-1./3.);
  // omega_b to omega_c*
  addFormFactor(5332,4334,2,4,3,3,5,4);_Nfi.push_back(1.     );_eta.push_back(0.);
  // omega_b to omega
  addFormFactor(5332,3334,2,4,3,3,5,3);_Nfi.push_back(1.     );_eta.push_back(0.);
  // omega_b to xi*
  addFormFactor(5332,3314,2,4,3,3,5,1);_Nfi.push_back(one3   );_eta.push_back(0.);
  // omega_c to omega
  addFormFactor(4332,3334,2,4,3,3,4,3);_Nfi.push_back(1.     );_eta.push_back(0.);
  // omega_c to xi*
  addFormFactor(4332,3324,2,4,3,3,4,2);_Nfi.push_back(one3   );_eta.push_back(0.);
  // lambda_c to lambda_0
  addFormFactor(4122,3122,2,2,1,2,4,3);_Nfi.push_back(1./sqrt(3.));_eta.push_back(1.);
  // xi_c to xi
  addFormFactor(4232,3322,2,2,2,3,4,3);_Nfi.push_back(1./sqrt(3.));_eta.push_back(1.);
  addFormFactor(4132,3312,2,2,1,3,4,3);_Nfi.push_back(1./sqrt(3.));_eta.push_back(1.);
  // initial number of form factors
  initialModes(numberOfFactors());
}

void ChengHeavyBaryonFormFactor::doinit() {
  BaryonFormFactor::doinit();
  // check the parameters are consistent
  unsigned int isize(numberOfFactors());
  if(isize!=_eta.size()||isize!=_Nfi.size())
    {throw InitException() << "Inconsistent paramters in ChengHeavyBaryon" 
			   << "FormFactor::doinit() " << Exception::abortnow;}
  Energy mi,mf,mq,mQ,lambda,delta,msum;
  int id0,id1,inspin,outspin,isp1,isp2,inq,outq;
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      // ids of the external particles
      particleID(ix,id0,id1);
      formFactorInfo(ix,inspin,outspin,isp1,isp2,inq,outq);
      id0=abs(id0);id1=abs(id1);
      mi=getParticleData(id0)->mass();
      mf=getParticleData(id1)->mass();
      msum=mi+mf;
      // masses of the incoming and outgoing quarks
      if((id0==4122&&id1==3122)||(id0==4232&&id1==3322)||(id0==4132&&id1==3312)||
	 (id0==4332&&id1==3334))
	{mq=_ms;mQ=_mc;}
      else if((id0==4332&&id1==3322)||(id0==4332&&id1==3324))
	{mq=_mu;mQ=_mc;}
      else if((id0==5122&&id1==4122)||(id0==5232&&id1==4232)||(id0==5132&&id1==4132)||
	      (id0==5332&&id1==4332)||(id0==5332&&id1==4334))
	{mq=_mc;mQ=_mb;}
      else if((id0==5122&&id1==3122)||(id0==5132&&id1==3312)||(id0==5232&&id1==3322)||
	      (id0==5332&&id1==3334))
	{mq=_ms;mQ=_mb;}
      else if((id0==5122&&id1==2112)||(id0==5132&&id1==3112)||(id0==5232&&id1==3212)||
	      (id0==5332&&id1==3312)||(id0==5232&&id1==3122)||(id0==5332&&id1==3314))
	{mq=_md;mQ=_mb;}
      else
	{throw InitException() << "Unknown decay in ChengHeavyBaryon" 
			       << "FormFactor::doinit() " << Exception::abortnow;}
      // parameters
      lambda = mf-mq;
      delta  = mi-mf;
      // compute the form-factors
      if(inspin==2&&outspin==2)
	{
	  _f1.push_back(_Nfi[ix]*(1.-0.5*delta/mi
				  +0.25*delta/mi/mq*(1.-0.5*lambda/mf)*
				  (mi+mf-_eta[ix]*delta)
				  -0.125*delta/mi/mf/mQ*lambda*(mi+mf+_eta[ix]*delta)));
	  _f2.push_back(_Nfi[ix]*msum*(0.5/mi+0.25/mi/mq*(1.-0.5*lambda/mf)*
				       (delta-(mi+mf)*_eta[ix])
				       -0.125*lambda/mi/mf/mQ*(delta+(mi+mf)*_eta[ix])));
	  _f3.push_back(_Nfi[ix]*msum*(0.5/mi-0.25/mi/mq*(1.-0.5*lambda/mf)*
				       (mi+mf-_eta[ix]*delta)
				       +0.125*lambda/mi/mf/mQ*(mi+mf+_eta[ix]*delta)));
	  _g1.push_back(_Nfi[ix]*_eta[ix]*(1.+0.25*delta*lambda*(1./mi/mq-1./mf/mQ)));
	  _g2.push_back(-0.25*msum*_Nfi[ix]*_eta[ix]*lambda*(1./mi/mq-1./mf/mQ));
	  _g3.push_back(-0.25*msum*_Nfi[ix]*_eta[ix]*lambda*(1./mi/mq+1./mf/mQ));
	}
      else if(inspin==2&&outspin==4)
	{
	  _f1.push_back(2.*_Nfi[ix]/sqrt(3.)*(1.+0.5*lambda*(1./mq+1./mQ)));
	  _f2.push_back(_Nfi[ix]*msum/sqrt(3.)/mi*(1.+0.5*lambda*(1./mq+1./mQ)));
	  _f3.push_back(-_Nfi[ix]*msum*msum/mi/mf/sqrt(3.)*
			(1.+0.5*lambda*(1./mq+1./mQ)));
	  _g1.push_back(-2./sqrt(3.)*_Nfi[ix]);
	  _g2.push_back(-_Nfi[ix]*msum/sqrt(3.)*lambda/mq/mi);
	  _g3.push_back(-_f3.back());
	}
      else
	{throw InitException() << "Unknown spin combination in ChengHeavyBaryon" 
			       << "FormFactor::doinit() "  << Exception::abortnow;}
    }
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      int id0,id1;
      particleID(ix,id0,id1);
      tcPDPtr part0=getParticleData(id0); Energy m0=part0->mass();
      tcPDPtr part1=getParticleData(id1); Energy m1=part1->mass();
      if ( part1->iSpin() == 2 ) {
	Complex f1v,f2v,f3v,f1a,f2a,f3a; // dummy variables
	SpinHalfSpinHalfFormFactor(ZERO,ix,id0,id1,m0,m1,f1v,f2v,f3v,f1a,f2a,f3a,FlavourInfo());
      }
      else {
	Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a; // dummy variables
	SpinHalfSpinThreeHalfFormFactor(ZERO,ix,id0,id1,m0,m1,f1v,f2v,f3v,
					f4v,f1a,f2a,f3a,f4a,FlavourInfo());
      }
    }
}
  
void ChengHeavyBaryonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mu,MeV) << ounit(_md,MeV) << ounit(_ms,MeV) << ounit(_mc,MeV) << ounit(_mb,MeV) 
     << _Nfi << _eta << _f1 << _f2 << _f3 
     << _g1 << _g2 << _g3 << ounit(_mVbc,MeV) << ounit(_mVbs,MeV) << ounit(_mVcs,MeV) 
     << ounit(_mVbd,MeV) << ounit(_mVcu,MeV) << ounit(_mAbc,MeV) 
     << ounit(_mAbs,MeV) << ounit(_mAcs,MeV) << ounit(_mAbd,MeV) << ounit(_mAcu,MeV);
 }
  

void ChengHeavyBaryonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mu,MeV) >> iunit(_md,MeV) >> iunit(_ms,MeV) >> iunit(_mc,MeV) >> iunit(_mb,MeV) 
     >> _Nfi >> _eta >> _f1 >> _f2 >> _f3 
     >> _g1 >> _g2 >> _g3 >> iunit(_mVbc,MeV) >> iunit(_mVbs,MeV) >> iunit(_mVcs,MeV) 
     >> iunit(_mVbd,MeV) >> iunit(_mVcu,MeV) >> iunit(_mAbc,MeV) 
     >> iunit(_mAbs,MeV) >> iunit(_mAcs,MeV) >> iunit(_mAbd,MeV) >> iunit(_mAcu,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ChengHeavyBaryonFormFactor,BaryonFormFactor>
describeHerwigChengHeavyBaryonFormFactor("Herwig::ChengHeavyBaryonFormFactor", "HwFormFactors.so");

void ChengHeavyBaryonFormFactor::Init() {

  static ClassDocumentation<ChengHeavyBaryonFormFactor> documentation
    ("The ChengHeavyBaryonFormFactor class is the implementation"
     " of the form-factors of PRD53, 1457 and PRD56, 2799 for the weak decay of"
     "baryons containing a heavy quark. This model can be used for either"
     "semi-leptonic decays, or with the factorization approximation for"
     " non-leptonic weak decays",
     "The weak decay of baryons containing a heavy quark used form factors from "
     "\\cite{Cheng:1995fe,Cheng:1996cs}.",
     "%\\cite{Cheng:1995fe}\n"
     "\\bibitem{Cheng:1995fe}\n"
     "  H.~Y.~Cheng and B.~Tseng,\n"
     "  %``1/M corrections to baryonic form-factors in the quark model,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 53} (1996) 1457\n"
     "  [Erratum-ibid.\\  D {\\bf 55} (1997) 1697]\n"
     "  [arXiv:hep-ph/9502391].\n"
     "  %%CITATION = PHRVA,D53,1457;%%\n"
     "%\\cite{Cheng:1996cs}\n"
     "\\bibitem{Cheng:1996cs}\n"
     "  H.~Y.~Cheng,\n"
     "  %``Nonleptonic weak decays of bottom baryons,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 56} (1997) 2799\n"
     "  [arXiv:hep-ph/9612223].\n"
     "  %%CITATION = PHRVA,D56,2799;%%\n"
     );

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceUpMass
    ("DownMass",
     "The consituent mass of the down quark",
     &ChengHeavyBaryonFormFactor::_md, GeV, 0.322*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceDownMass
    ("UpMass",
     "The consituent mass of the up quark",
     &ChengHeavyBaryonFormFactor::_mu, GeV, 0.338*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfacStrangeMass
    ("StrangeMass",
     "The consituent mass of the strange quark",
     &ChengHeavyBaryonFormFactor::_ms, GeV, 0.510*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The consituent mass of the charm quark",
     &ChengHeavyBaryonFormFactor::_mc, GeV, 1.6*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The consituent mass of the bottom quark",
     &ChengHeavyBaryonFormFactor::_mb, GeV, 5.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMassbc
    ("VectorMassbc",
     "The vector mass for the b->c transitions.",
     &ChengHeavyBaryonFormFactor::_mVbc, GeV, 6.34*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMassbc
    ("AxialMassbc",
     "The axial-vector mass for the b->c transitions.",
     &ChengHeavyBaryonFormFactor::_mAbc, GeV, 6.73*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMassbs
    ("VectorMassbs",
     "The vector mass for the b->s transitions.",
     &ChengHeavyBaryonFormFactor::_mVbs, GeV, 5.42*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMassbs
    ("AxialMassbs",
     "The axial-vector mass for the b->s transitions.",
     &ChengHeavyBaryonFormFactor::_mAbs, GeV, 5.86*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMassbd
    ("VectorMassbd",
     "The vector mass for the b->d transitions.",
     &ChengHeavyBaryonFormFactor::_mVbd, GeV, 5.32*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMassbd
    ("AxialMassbd",
     "The axial-vector mass for the b->d transitions.",
     &ChengHeavyBaryonFormFactor::_mAbd, GeV, 5.71*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMasscs
    ("VectorMasscs",
     "The vector mass for the c->s transitions.",
     &ChengHeavyBaryonFormFactor::_mVcs, GeV, 2.11*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMasscs
    ("AxialMasscs",
     "The axial-vector mass for the c->s transitions.",
     &ChengHeavyBaryonFormFactor::_mAcs, GeV, 2.54*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceVectorMasscu
    ("VectorMasscu",
     "The vector mass for the c->u transitions.",
     &ChengHeavyBaryonFormFactor::_mVcu, GeV, 2.01*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ChengHeavyBaryonFormFactor,Energy> interfaceAxialMasscu
    ("AxialMasscu",
     "The axial-vector mass for the c->u transitions.",
     &ChengHeavyBaryonFormFactor::_mAcu, GeV, 2.42*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static ParVector<ChengHeavyBaryonFormFactor,double> interfaceNfi
    ("Nfi",
     "The prefactor for a given form factor",
     &ChengHeavyBaryonFormFactor::_Nfi, -1, 1.0, -10.0, 10.0,
     false, false, true);

  static ParVector<ChengHeavyBaryonFormFactor,double> interfaceEta
    ("Eta",
     "The eta parameter for the form factor",
     &ChengHeavyBaryonFormFactor::_eta, -1, 0.0, -10.0, 10.0,
     false, false, true);

}

// form factor for spin-1/2 to spin-1/2
void ChengHeavyBaryonFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc,int id0,int id1, Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo ,
			   Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  id0=abs(id0);
  id1=abs(id1);
  // masses for the energy dependence of the form-factors
  Energy mV(ZERO),mA(ZERO);
   if((id0==4122&&id1==3122)||(id0==4232&&id1==3322)||(id0==4132&&id1==3312)||
      (id0==4332&&id1==3334))
     {mA=_mAcs;mV=_mVcs;}
   else if((id0==4332&&id1==3322)||(id0==4332&&id1==3324))
     {mA=_mAcu;mV=_mVcu;}
   else if((id0==5122&&id1==4122)||(id0==5232&&id1==4232)||(id0==5132&&id1==4132)||
	   (id0==5332&&id1==4332)||(id0==5332&&id1==4334))
     {mA=_mAbc;mV=_mVbc;}
   else if((id0==5122&&id1==3122)||(id0==5132&&id1==3312)||(id0==5232&&id1==3322)||
	   (id0==5332&&id1==3334))
     {mA=_mAbs;mV=_mVbs;}
   else if((id0==5122&&id1==2112)||(id0==5132&&id1==3112)||(id0==5232&&id1==3212)||
	   (id0==5332&&id1==3312)||(id0==5232&&id1==3122)||(id0==5332&&id1==3314))
     {mA=_mAbd;mV=_mVbd;}
   Energy delta=m0-m1;
   double Vfact = (1.-delta*delta/mV/mV)/(1.-q2/mV/mV);
   Vfact *=Vfact;
   double Afact = (1.-delta*delta/mA/mA)/(1.-q2/mA/mA);
   Afact *=Afact;
   f1v = _f1[iloc]*Vfact;
   f2v = _f2[iloc]*Vfact;
   f3v = _f3[iloc]*Vfact;
   f1a =-_g1[iloc]*Afact;
   f2a =-_g2[iloc]*Afact;
   f3a =-_g3[iloc]*Afact;
}

// form factor for spin-1/2 to spin-3/2
void ChengHeavyBaryonFormFactor::
SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int iloc, int id0, int id1,
				Energy m0, Energy m1,
				Complex & g1v,Complex & g2v,Complex & g3v,
				Complex & g4v,Complex & g1a,Complex & g2a,
				Complex & g3a,Complex & g4a,
				FlavourInfo , Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  id0=abs(id0);
  id1=abs(id1);
  // masses for the energy dependence of the form-factors
  Energy mV(ZERO),mA(ZERO);
   if((id0==4122&&id1==3122)||(id0==4232&&id1==3322)||(id0==4132&&id1==3312)||
      (id0==4332&&id1==3334))
     {mA=_mAcs;mV=_mVcs;}
   else if((id0==4332&&id1==3322)||(id0==4332&&id1==3324))
     {mA=_mAcu;mV=_mVcu;}
   else if((id0==5122&&id1==4122)||(id0==5232&&id1==4232)||(id0==5132&&id1==4132)||
	   (id0==5332&&id1==4332)||(id0==5332&&id1==4334))
     {mA=_mAbc;mV=_mVbc;}
   else if((id0==5122&&id1==3122)||(id0==5132&&id1==3312)||(id0==5232&&id1==3322)||
	   (id0==5332&&id1==3334))
     {mA=_mAbs;mV=_mVbs;}
   else if((id0==5122&&id1==2112)||(id0==5132&&id1==3112)||(id0==5232&&id1==3212)||
	   (id0==5332&&id1==3312)||(id0==5232&&id1==3122)||(id0==5332&&id1==3314))
     {mA=_mAbd;mV=_mVbd;}
   Energy delta=m0-m1;
   double Vfact = (1.-delta*delta/mV/mV)/(1.-q2/mV/mV);
   Vfact *=Vfact;
   double Afact = (1.-delta*delta/mA/mA)/(1.-q2/mA/mA);
   Afact *=Afact;
   g1v = _f1[iloc]*Vfact;
   g2v = _f2[iloc]*Vfact;
   g3v = _f3[iloc]*Vfact;
   g4v = 0.;
   g1a =-_g1[iloc]*Afact;
   g2a =-_g2[iloc]*Afact;
   g3a =-_g3[iloc]*Afact;
   g4a = 0.;
}

// output the information for the database
void ChengHeavyBaryonFormFactor::dataBaseOutput(ofstream& output,bool header,
						bool create) const
{
  if(header){output << "update decayers set parameters=\"";}
  if(create)
    {output << "create Herwig::ChengHeavyBaryonFormFactor " << name() << " \n";}
  output << "newdef " << name() << ":DownMass     " << _md/GeV << " \n";
  output << "newdef " << name() << ":UpMass       " << _mu/GeV << " \n";
  output << "newdef " << name() << ":StrangeMass  " << _ms/GeV << " \n";
  output << "newdef " << name() << ":CharmMass    " << _mc/GeV << " \n";
  output << "newdef " << name() << ":BottomMass   " << _mb/GeV << " \n";
  output << "newdef " << name() << ":VectorMassbc " << _mVbc/GeV << " \n";
  output << "newdef " << name() << ":AxialMassbc  " << _mAbc/GeV << " \n";
  output << "newdef " << name() << ":VectorMassbs " << _mVbs/GeV << " \n";
  output << "newdef " << name() << ":AxialMassbs  " << _mAbs/GeV << " \n";
  output << "newdef " << name() << ":VectorMassbd " << _mVbd/GeV << " \n";
  output << "newdef " << name() << ":AxialMassbd  " << _mAbd/GeV << " \n";
  output << "newdef " << name() << ":VectorMasscs " << _mVcs/GeV << " \n";
  output << "newdef " << name() << ":AxialMasscs  " << _mAcs/GeV << " \n";
  output << "newdef " << name() << ":VectorMasscu " << _mVcu/GeV << " \n";
  output << "newdef " << name() << ":AxialMasscu  " << _mAcu/GeV << " \n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix)
    {
      if(ix<initialModes())
	{
	  output << "newdef " << name() << ":Nfi " << ix << "  " 
		<< _Nfi[ix] << endl;
	  output << "newdef " << name() << ":Eta " << ix << "  " 
		<< _eta[ix] << endl;
	}
      else
	{
	  output << "insert " << name() << ":Nfi " << ix << "  " 
		<< _Nfi[ix] << endl;
	  output << "insert " << name() << ":Eta " << ix << "  " 
		<< _eta[ix] << endl;
	}
    }
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;}
}
#line 1 "./ISGW2FormFactor.cc"
// -*- C++ -*-
//
// ISGW2FormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGW2FormFactor class.
//

#include "ISGW2FormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ISGW2FormFactor::ISGW2FormFactor() { 
  // include the a(omega) piece (off by default)
  _includeaW = false;
  // values of a_S at matching scale
  _alphamuQM=0.6;
  // the quark masses
  _mdown    = 0.33*GeV;
  _mup      = 0.33*GeV;
  _mstrange = 0.55*GeV;
  _mcharm   = 1.82*GeV;
  _mbottom  = 5.20*GeV;
  // 1 1S0 parameters
  // the wavefunction parameters
  _beta1S0ud = 0.41*GeV;
  _beta1S0us = 0.44*GeV;
  _beta1S0ss = 0.53*GeV;
  _beta1S0cu = 0.45*GeV;
  _beta1S0cs = 0.56*GeV;
  _beta1S0ub = 0.43*GeV;
  _beta1S0sb = 0.54*GeV;
  _beta1S0cc = 0.88*GeV;
  _beta1S0bc = 0.92*GeV;
  // 1 3S1 parameters
  // the wavefunction parameters
  _beta3S1ud = 0.30*GeV;
  _beta3S1us = 0.33*GeV;
  _beta3S1ss = 0.37*GeV;
  _beta3S1cu = 0.38*GeV;
  _beta3S1cs = 0.44*GeV;
  _beta3S1ub = 0.40*GeV;
  _beta3S1sb = 0.49*GeV;
  _beta3S1cc = 0.62*GeV;
  _beta3S1bc = 0.75*GeV;
  // 1P parameters
  // the wavefunction parameters
  _beta1Pud = 0.28*GeV;
  _beta1Pus = 0.30*GeV;
  _beta1Pss = 0.33*GeV;
  _beta1Pcu = 0.33*GeV;
  _beta1Pcs = 0.38*GeV;
  _beta1Pub = 0.35*GeV;
  _beta1Psb = 0.41*GeV;
  _beta1Pcc = 0.52*GeV;
  _beta1Pbc = 0.60*GeV;
  // relativistic correction factors
  _CfDrho     = 0.889;
  _CfDKstar   = 0.928;
  _CfDsKstar  = 0.873;
  _CfDsphi    = 0.911;
  _CfBrho     = 0.905;
  _CfBDstar   = 0.989;
  _CfBsKstar  = 0.892;
  _CfBsDstar  = 0.984;
  _CfBcDstar  = 0.868;
  _CfBcpsi    = 0.967;
  _CfBcBsstar = 1.0;
  _CfBcBstar  = 1.0;
  // eta-eta' mixing angle
  _thetaeta = -Constants::pi/9.;
  // B_c to d cbar
  addFormFactor(-541,-411  ,0,-4, 5, 1);
  addFormFactor(-541,-413  ,1,-4, 5, 1);
  addFormFactor(-541,-415  ,2,-4, 5, 1);
  addFormFactor(-541, 10413,1,-2, 5, 1);
  addFormFactor(-541,-20413,1,-4, 5, 1);
  addFormFactor(-541, 10411,0, 4, 5, 1);
  // B_c to u cbar
  addFormFactor(-541,-421  ,0,-4, 5, 2);
  addFormFactor(-541,-423  ,1,-4, 5, 2);
  addFormFactor(-541,-425  ,2,-4, 5, 2);
  addFormFactor(-541,-10423,1,-4, 5, 2);
  addFormFactor(-541,-20423,1,-4, 5, 2);
  addFormFactor(-541,-10421,0,-4, 5, 2);
  // B_c to s cbar
  addFormFactor(-541,-431  ,0,-4, 5, 3);
  addFormFactor(-541,-433  ,1,-4, 5, 3);
  addFormFactor(-541,-435  ,2,-4, 5, 3);
  addFormFactor(-541,-10433,1,-4, 5, 3);
  addFormFactor(-541,-20433,1,-4, 5, 3);
  addFormFactor(-541,-10431,0, 4, 5, 3);
  // B_c decays to c cbar
  addFormFactor(-541, 441  ,0,-4, 5, 4);
  addFormFactor(-541, 443  ,1,-4, 5, 4);
  addFormFactor(-541, 445  ,2,-4, 5, 4);
  addFormFactor(-541, 10443,1,-4, 5, 4);
  addFormFactor(-541, 20443,1,-4, 5, 4);
  addFormFactor(-541, 10441,0, 4, 5, 4);
  // B_c to b dbar
  addFormFactor( 541, 511  ,0, 5,-4,-1);
  addFormFactor( 541, 513  ,1, 5,-4,-1);
  addFormFactor( 541, 515  ,2, 5,-4,-1);
  addFormFactor( 541, 10513,1, 5,-4,-1); 
  addFormFactor( 541, 20513,1, 5,-4,-1);
  addFormFactor( 541, 10511,0, 5,-4,-1);
  // B_c to b ubar
  addFormFactor( 541, 521  ,0, 5,-4,-2);
  addFormFactor( 541, 523  ,1, 5,-4,-2);
  addFormFactor( 541, 525  ,2, 5,-4,-2);
  addFormFactor( 541, 10523,1, 5,-4,-2); 
  addFormFactor( 541, 20523,1, 5,-4,-2);
  addFormFactor( 541, 10521,0, 5,-4,-2);
  // B_c decays to s bbar
  addFormFactor( 541, 531  ,0, 5,-4,-3);
  addFormFactor( 541, 533  ,1, 5,-4,-3);
  addFormFactor( 541, 535  ,2, 5,-4,-3);
  addFormFactor( 541, 10533,1, 5,-4,-3);
  addFormFactor( 541, 20533,1, 5,-4,-3);
  addFormFactor( 541, 10531,0, 5,-4,-3);
  // B_s to d sbar
  addFormFactor(-531, 311  ,0, 3,-5,-1);
  addFormFactor(-531, 313  ,1, 3,-5,-1);
  addFormFactor(-531, 315  ,2, 3,-5,-1);
  addFormFactor(-531, 10313,1, 3,-5,-1);
  addFormFactor(-531, 20313,1, 3,-5,-1);
  addFormFactor(-531, 10311,0, 3,-5,-1);
  // B_s to u sbar
  addFormFactor(-531, 321  ,0, 3,-5,-2);
  addFormFactor(-531, 323  ,1, 3,-5,-2);
  addFormFactor(-531, 325  ,2, 3,-5,-2);
  addFormFactor(-531, 10323,1, 3,-5,-2);
  addFormFactor(-531, 20323,1, 3,-5,-2);
  addFormFactor(-531, 10321,0, 3,-5,-2);
  // B_s to s sbar
  addFormFactor(-531, 221  ,0, 3,-5,-3);
  addFormFactor(-531, 331  ,0, 3,-5,-3);
  addFormFactor(-531, 333  ,1, 3,-5,-3);
  addFormFactor(-531, 335  ,2, 3,-5,-3);
  addFormFactor(-531, 10333,1, 3,-5,-3);
  addFormFactor(-531, 20333,1, 3,-5,-3);
  addFormFactor(-531, 10331,0, 3,-5,-3);
  // B_s decays to c sbar
  addFormFactor(-531, 431  ,0, 3,-5,-4);
  addFormFactor(-531, 433  ,1, 3,-5,-4);
  addFormFactor(-531, 435  ,2, 3,-5,-4);
  addFormFactor(-531, 10433,1, 3,-5,-4);
  addFormFactor(-531, 20433,1, 3,-5,-4);
  addFormFactor(-531, 10431,0, 3,-5,-4);
  // B_u decays to d ubar
  addFormFactor(-521,-211  ,0,-2, 5, 1);
  addFormFactor(-521,-213  ,1,-2, 5, 1);
  addFormFactor(-521,-215  ,2,-2, 5, 1);
  addFormFactor(-521,-10213,1,-2, 5, 1);
  addFormFactor(-521,-20213,1,-2, 5, 1);
  addFormFactor(-521,-10211,0,-2, 5, 1);
  // B_u to uu (I=0)
  addFormFactor(-521, 221  ,0,-2, 5, 2);
  addFormFactor(-521, 331  ,0,-2, 5, 2);
  addFormFactor(-521, 223  ,1,-2, 5, 2);
  addFormFactor(-521, 225  ,2,-2, 5, 2);
  addFormFactor(-521, 10223,1,-2, 5, 2);
  addFormFactor(-521, 20223,1,-2, 5, 2);
  addFormFactor(-521, 10221,0,-2, 5, 2);
  // B_u to uu (I=1)
  addFormFactor(-521, 111  ,0,-2, 5, 2);
  addFormFactor(-521, 113  ,1,-2, 5, 2);
  addFormFactor(-521, 115  ,2,-2, 5, 2);
  addFormFactor(-521, 10113,1,-2, 5, 2);
  addFormFactor(-521, 20113,1,-2, 5, 2);
  addFormFactor(-521, 10111,0,-2, 5, 2);
  // B_u decays to s ubar
  addFormFactor(-521,-321  ,0,-2, 5, 3);
  addFormFactor(-521,-323  ,1,-2, 5, 3);
  addFormFactor(-521,-325  ,2,-2, 5, 3);
  addFormFactor(-521,-10323,1,-2, 5, 3);
  addFormFactor(-521,-20323,1,-2, 5, 3);
  addFormFactor(-521,-10321,0,-2, 5, 3);
  // B_u decays to c ubar
  addFormFactor(-521, 421  ,0,-2, 5, 4);
  addFormFactor(-521, 423  ,1,-2, 5, 4);
  addFormFactor(-521, 425  ,2,-2, 5, 4);
  addFormFactor(-521, 10423,1,-2, 5, 4);
  addFormFactor(-521, 20423,1,-2, 5, 4);
  addFormFactor(-521, 10421,0,-2, 5, 4);
  // B_d decays to d dbar (I=0)
  addFormFactor(-511, 221  ,0, 1,-5,-1); 
  addFormFactor(-511, 331  ,0, 1,-5,-1); 
  addFormFactor(-511, 223  ,1, 1,-5,-1); 
  addFormFactor(-511, 225  ,2, 1,-5,-1); 
  addFormFactor(-511, 10223,1, 1,-5,-1); 
  addFormFactor(-511, 20223,1, 1,-5,-1);
  addFormFactor(-511, 10221,0, 1,-5,-1);
  // B_d decays to d dbar (I=1)
  addFormFactor(-511, 111  ,0, 1,-5,-1); 
  addFormFactor(-511, 113  ,1, 1,-5,-1); 
  addFormFactor(-511, 115  ,2, 1,-5,-1); 
  addFormFactor(-511, 10113,1, 1,-5,-1); 
  addFormFactor(-511, 20113,1, 1,-5,-1);
  addFormFactor(-511, 10111,0, 1,-5,-1);
  // B_d decays to u dbar
  addFormFactor(-511, 211  ,0, 1,-5,-2); 
  addFormFactor(-511, 213  ,1, 1,-5,-2); 
  addFormFactor(-511, 215  ,2, 1,-5,-2); 
  addFormFactor(-511, 10213,1, 1,-5,-2); 
  addFormFactor(-511, 20213,1, 1,-5,-2);
  addFormFactor(-511, 10211,0, 1,-5,-2);
  // B_d decays to  s dbar
  addFormFactor(-511, 311  ,0, 1,-5,-3);
  addFormFactor(-511, 313  ,1, 1,-5,-3);
  addFormFactor(-511, 315  ,2, 1,-5,-3);
  addFormFactor(-511, 10313,1, 1,-5,-3);
  addFormFactor(-511, 20313,1, 1,-5,-3);
  addFormFactor(-511, 10311,0, 1,-5,-3);
  // B_d decays to  c dbar
  addFormFactor(-511, 411  ,0, 1,-5,-4);
  addFormFactor(-511, 413  ,1, 1,-5,-4);
  addFormFactor(-511, 415  ,2, 1,-5,-4);
  addFormFactor(-511, 10413,1, 1,-5,-4);
  addFormFactor(-511, 20413,1, 1,-5,-4);
  addFormFactor(-511, 10411,0, 1,-5,-4);
  // D_s to d sbar
  addFormFactor( 431,   311,0,-3, 4, 1);
  addFormFactor( 431,   313,1,-3, 4, 1);
  addFormFactor( 431,   315,2,-3, 4, 1);
  addFormFactor( 431, 10313,1,-3, 4, 1);
  addFormFactor( 431, 20313,1,-3, 4, 1);
  addFormFactor( 431, 10311,0,-3, 4, 1);
  // D_s to u sbar
  addFormFactor( 431,   321,0,-3, 4, 2);
  addFormFactor( 431,   323,1,-3, 4, 2);
  addFormFactor( 431,   325,2,-3, 4, 2);
  addFormFactor( 431, 10323,1,-3, 4, 2);
  addFormFactor( 431, 20323,1,-3, 4, 2);
  addFormFactor( 431, 10321,0,-3, 4, 2);
  // D_s to s sbar
  addFormFactor( 431, 221  ,0,-3, 4, 3);
  addFormFactor( 431, 331  ,0,-3, 4, 3);
  addFormFactor( 431, 333  ,1,-3, 4, 3);
  addFormFactor( 431, 335  ,2,-3, 4, 3);
  addFormFactor( 431, 10333,1,-3, 4, 3);
  addFormFactor( 431, 20333,1,-3, 4, 3);
  addFormFactor( 431, 10331,0,-3, 4, 3);
  // D0 to d ubar
  addFormFactor( 421,-211  ,0,-2, 4, 1);
  addFormFactor( 421,-213  ,1,-2, 4, 1);
  addFormFactor( 421,-215  ,2,-2, 4, 1);
  addFormFactor( 421,-10213,1,-2, 4, 1);
  addFormFactor( 421,-20213,1,-2, 4, 1);
  addFormFactor( 421,-10211,0,-2, 4, 1);
  // D0 to u ubar (I=1)
  addFormFactor( 421, 111  ,0,-2, 4, 2);
  addFormFactor( 421, 113  ,1,-2, 4, 2);
  addFormFactor( 421, 115  ,2,-2, 4, 2);
  addFormFactor( 421, 10113,1,-2, 4, 2);
  addFormFactor( 421, 20113,1,-2, 4, 2);
  addFormFactor( 421, 10111,0,-2, 4, 2);
  // D0 to u ubar (I=0)
  addFormFactor( 421, 221  ,0,-2, 4, 2);
  addFormFactor( 421, 331  ,0,-2, 4, 2);
  addFormFactor( 421, 223  ,1,-2, 4, 2);
  addFormFactor( 421, 225  ,2,-2, 4, 2);
  addFormFactor( 421, 10223,1,-2, 4, 2);
  addFormFactor( 421, 20223,1,-2, 4, 2);
  addFormFactor( 421, 10221,0,-2, 4, 2);
  // D0 to s ubar
  addFormFactor( 421,-321  ,0,-2, 4, 3);
  addFormFactor( 421,-323  ,1,-2, 4, 3);
  addFormFactor( 421,-325  ,2,-2, 4, 3);
  addFormFactor( 421,-10323,1,-2, 4, 3);
  addFormFactor( 421,-20323,1,-2, 4, 3);
  addFormFactor( 421,-10321,0,-2, 4, 3);
  // D+ to d dbar I=0
  addFormFactor( 411, 221  ,0,-1, 4, 1);
  addFormFactor( 411, 331  ,0,-1, 4, 1);
  addFormFactor( 411, 223  ,1,-1, 4, 1);
  addFormFactor( 411, 225  ,2,-1, 4, 1);
  addFormFactor( 411, 10223,1,-1, 4, 1);
  addFormFactor( 411, 20223,1,-1, 4, 1);
  addFormFactor( 411, 10221,0,-1, 4, 1);
  // D+ to d dbar I=1
  addFormFactor( 411, 111  ,0,-1, 4, 1);
  addFormFactor( 411, 113  ,1,-1, 4, 1);
  addFormFactor( 411, 115  ,2,-1, 4, 1);
  addFormFactor( 411, 10113,1,-1, 4, 1);
  addFormFactor( 411, 20113,1,-1, 4, 1);
  addFormFactor( 411, 10111,0,-1, 4, 1);
  // D+ to u dbar
  addFormFactor( 411, 211  ,0,-1, 4, 2);
  addFormFactor( 411, 213  ,1,-1, 4, 2);
  addFormFactor( 411, 215  ,2,-1, 4, 2);
  addFormFactor( 411, 10213,1,-1, 4, 2);
  addFormFactor( 411, 20213,1,-1, 4, 2);
  addFormFactor( 411, 10211,0,-1, 4, 2);
  // D+ to s dbar
  addFormFactor( 411,-311  ,0,-1, 4, 3);
  addFormFactor( 411,-313  ,1,-1, 4, 3);
  addFormFactor( 411,-315  ,2,-1, 4, 3);
  addFormFactor( 411,-10313,1,-1, 4, 3);
  addFormFactor( 411,-20313,1,-1, 4, 3);
  addFormFactor( 411,-10311,0,-1, 4, 3);
  // set the initial number of modes
  initialModes(numberOfFactors());
}			     

void ISGW2FormFactor::doinit() {
  ScalarFormFactor::doinit();
  // set up the quark masses
  _mquark.resize(5);
  _mquark[0]=_mdown;
  _mquark[1]=_mup;
  _mquark[2]=_mstrange;
  _mquark[3]=_mcharm;
  _mquark[4]=_mbottom;
  // value of alpha_S at the quark masses
  _alphaQ.resize(5);
  for(unsigned int ix=0;ix<5;++ix) {
    _alphaQ[ix]=alphaS(_mquark[ix],sqr(_mquark[ix]));
  }
  _beta1S0.resize(5,vector<Energy>(5));
  _mass1S0.resize(5,vector<Energy>(5));
  _beta3S1.resize(5,vector<Energy>(5));
  _beta1P .resize(5,vector<Energy>(5));
  _massPoh.resize(5,vector<Energy>(5));
  _massPth.resize(5,vector<Energy>(5));
  // set up the beta values
  _beta1S0[0][0] = _beta1S0ud;_beta3S1[0][0] = _beta3S1ud;_beta1P[0][0] = _beta1Pud;
  _beta1S0[1][0] = _beta1S0ud;_beta3S1[1][0] = _beta3S1ud;_beta1P[1][0] = _beta1Pud;
  _beta1S0[2][0] = _beta1S0us;_beta3S1[2][0] = _beta3S1us;_beta1P[2][0] = _beta1Pus;
  _beta1S0[3][0] = _beta1S0cu;_beta3S1[3][0] = _beta3S1cu;_beta1P[3][0] = _beta1Pcu;
  _beta1S0[4][0] = _beta1S0ub;_beta3S1[4][0] = _beta3S1ub;_beta1P[4][0] = _beta1Pub;
  _beta1S0[0][1] = _beta1S0ud;_beta3S1[0][1] = _beta3S1ud;_beta1P[0][1] = _beta1Pud;
  _beta1S0[1][1] = _beta1S0ud;_beta3S1[1][1] = _beta3S1ud;_beta1P[1][1] = _beta1Pud;
  _beta1S0[2][1] = _beta1S0us;_beta3S1[2][1] = _beta3S1us;_beta1P[2][1] = _beta1Pus;
  _beta1S0[3][1] = _beta1S0cu;_beta3S1[3][1] = _beta3S1cu;_beta1P[3][1] = _beta1Pcu;
  _beta1S0[4][1] = _beta1S0ub;_beta3S1[4][1] = _beta3S1ub;_beta1P[4][1] = _beta1Pub;
  _beta1S0[0][2] = _beta1S0us;_beta3S1[0][2] = _beta3S1us;_beta1P[0][2] = _beta1Pus;
  _beta1S0[1][2] = _beta1S0us;_beta3S1[1][2] = _beta3S1us;_beta1P[1][2] = _beta1Pus;
  _beta1S0[2][2] = _beta1S0ss;_beta3S1[2][2] = _beta3S1ss;_beta1P[2][2] = _beta1Pss;
  _beta1S0[3][2] = _beta1S0cs;_beta3S1[3][2] = _beta3S1cs;_beta1P[3][2] = _beta1Pcs;
  _beta1S0[4][2] = _beta1S0sb;_beta3S1[4][2] = _beta3S1sb;_beta1P[4][2] = _beta1Psb;
  _beta1S0[0][3] = _beta1S0cu;_beta3S1[0][3] = _beta3S1cu;_beta1P[0][3] = _beta1Pcu;
  _beta1S0[1][3] = _beta1S0cu;_beta3S1[1][3] = _beta3S1cu;_beta1P[1][3] = _beta1Pcu;
  _beta1S0[2][3] = _beta1S0cs;_beta3S1[2][3] = _beta3S1cs;_beta1P[2][3] = _beta1Pcs;
  _beta1S0[3][3] = _beta1S0cc;_beta3S1[3][3] = _beta3S1cc;_beta1P[3][3] = _beta1Pcc;
  _beta1S0[4][3] = _beta1S0bc;_beta3S1[4][3] = _beta3S1bc;_beta1P[4][3] = _beta1Pbc;
  _beta1S0[0][4] = _beta1S0ub;_beta3S1[0][4] = _beta3S1ub;_beta1P[0][4] = _beta1Pub;
  _beta1S0[1][4] = _beta1S0ub;_beta3S1[1][4] = _beta3S1ub;_beta1P[1][4] = _beta1Pub;
  _beta1S0[2][4] = _beta1S0sb;_beta3S1[2][4] = _beta3S1sb;_beta1P[2][4] = _beta1Psb;
  _beta1S0[3][4] = _beta1S0bc;_beta3S1[3][4] = _beta3S1bc;_beta1P[3][4] = _beta1Pbc;
  _beta1S0[4][4] = ZERO    ;_beta3S1[4][4] = ZERO    ;_beta1P[4][4] = ZERO   ;
  // set up the values of mbar
  // get the particle data objects
  tcPDPtr p1S0[5][5],p3S1[5][5];
  p1S0[0][0] = getParticleData(111); p3S1[0][0] = getParticleData(113);
  p1S0[1][0] = getParticleData(211); p3S1[1][0] = getParticleData(213);
  p1S0[2][0] = getParticleData(311); p3S1[2][0] = getParticleData(313);
  p1S0[3][0] = getParticleData(411); p3S1[3][0] = getParticleData(413);
  p1S0[4][0] = getParticleData(511); p3S1[4][0] = getParticleData(513);
  p1S0[0][1] = getParticleData(211); p3S1[0][1] = getParticleData(213);
  p1S0[1][1] = getParticleData(111); p3S1[1][1] = getParticleData(113);
  p1S0[2][1] = getParticleData(321); p3S1[2][1] = getParticleData(323);
  p1S0[3][1] = getParticleData(421); p3S1[3][1] = getParticleData(423);
  p1S0[4][1] = getParticleData(521); p3S1[4][1] = getParticleData(523);
  p1S0[0][2] = getParticleData(311); p3S1[0][2] = getParticleData(313);
  p1S0[1][2] = getParticleData(321); p3S1[1][2] = getParticleData(323);
  p1S0[2][2] = getParticleData(331); p3S1[2][2] = getParticleData(333);
  p1S0[3][2] = getParticleData(431); p3S1[3][2] = getParticleData(433);
  p1S0[4][2] = getParticleData(531); p3S1[4][2] = getParticleData(533);
  p1S0[0][3] = getParticleData(411); p3S1[0][3] = getParticleData(413);
  p1S0[1][3] = getParticleData(421); p3S1[1][3] = getParticleData(423);
  p1S0[2][3] = getParticleData(431); p3S1[2][3] = getParticleData(433);
  p1S0[3][3] = getParticleData(441); p3S1[3][3] = getParticleData(443);
  p1S0[4][3] = getParticleData(541); p3S1[4][3] = getParticleData(543);
  p1S0[0][4] = getParticleData(511); p3S1[0][4] = getParticleData(513);
  p1S0[1][4] = getParticleData(521); p3S1[1][4] = getParticleData(523);
  p1S0[2][4] = getParticleData(531); p3S1[2][4] = getParticleData(533);
  p1S0[3][4] = getParticleData(541); p3S1[3][4] = getParticleData(543);
  p1S0[4][4] = getParticleData(551); p3S1[4][4] = getParticleData(553);
  tcPDPtr p3P0[5][5],p3P1[5][5],p3P2[5][5],p1P1[5][5];
  p3P0[0][0] = getParticleData(10111); p3P1[0][0] = getParticleData(20113);
  p3P0[1][0] = getParticleData(10211); p3P1[1][0] = getParticleData(20213);
  p3P0[2][0] = getParticleData(10311); p3P1[2][0] = getParticleData(20313);
  p3P0[3][0] = getParticleData(10411); p3P1[3][0] = getParticleData(20413);
  p3P0[4][0] = getParticleData(10511); p3P1[4][0] = getParticleData(20513);
  p3P0[0][1] = getParticleData(10211); p3P1[0][1] = getParticleData(20213);
  p3P0[1][1] = getParticleData(10111); p3P1[1][1] = getParticleData(20113);
  p3P0[2][1] = getParticleData(10321); p3P1[2][1] = getParticleData(20323);
  p3P0[3][1] = getParticleData(10421); p3P1[3][1] = getParticleData(20423);
  p3P0[4][1] = getParticleData(10521); p3P1[4][1] = getParticleData(20523);
  p3P0[0][2] = getParticleData(10311); p3P1[0][2] = getParticleData(20313);
  p3P0[1][2] = getParticleData(10321); p3P1[1][2] = getParticleData(20323);
  p3P0[2][2] = getParticleData(10331); p3P1[2][2] = getParticleData(20333);
  p3P0[3][2] = getParticleData(10431); p3P1[3][2] = getParticleData(20433);
  p3P0[4][2] = getParticleData(10531); p3P1[4][2] = getParticleData(20533);
  p3P0[0][3] = getParticleData(10411); p3P1[0][3] = getParticleData(20413);
  p3P0[1][3] = getParticleData(10421); p3P1[1][3] = getParticleData(20423);
  p3P0[2][3] = getParticleData(10431); p3P1[2][3] = getParticleData(20433);
  p3P0[3][3] = getParticleData(10441); p3P1[3][3] = getParticleData(20443);
  p3P0[4][3] = getParticleData(10541); p3P1[4][3] = getParticleData(20543);
  p3P0[0][4] = getParticleData(10511); p3P1[0][4] = getParticleData(20513);
  p3P0[1][4] = getParticleData(10521); p3P1[1][4] = getParticleData(20523);
  p3P0[2][4] = getParticleData(10531); p3P1[2][4] = getParticleData(20533);
  p3P0[3][4] = getParticleData(10541); p3P1[3][4] = getParticleData(20543);
  p3P0[4][4] = getParticleData(10551); p3P1[4][4] = getParticleData(20553);
  p1P1[0][0]=getParticleData(10113); p3P2[0][0]=getParticleData(115);
  p1P1[1][0]=getParticleData(10213); p3P2[1][0]=getParticleData(215);
  p1P1[2][0]=getParticleData(10313); p3P2[2][0]=getParticleData(315);
  p1P1[3][0]=getParticleData(10413); p3P2[3][0]=getParticleData(415);
  p1P1[4][0]=getParticleData(10513); p3P2[4][0]=getParticleData(515);
  p1P1[0][1]=getParticleData(10213); p3P2[0][1]=getParticleData(215);
  p1P1[1][1]=getParticleData(10113); p3P2[1][1]=getParticleData(115);
  p1P1[2][1]=getParticleData(10323); p3P2[2][1]=getParticleData(325);
  p1P1[3][1]=getParticleData(10423); p3P2[3][1]=getParticleData(425);
  p1P1[4][1]=getParticleData(10523); p3P2[4][1]=getParticleData(525);
  p1P1[0][2]=getParticleData(10313); p3P2[0][2]=getParticleData(315);
  p1P1[1][2]=getParticleData(10323); p3P2[1][2]=getParticleData(325);
  p1P1[2][2]=getParticleData(10333); p3P2[2][2]=getParticleData(335);
  p1P1[3][2]=getParticleData(10433); p3P2[3][2]=getParticleData(435);
  p1P1[4][2]=getParticleData(10533); p3P2[4][2]=getParticleData(535);
  p1P1[0][3]=getParticleData(10413); p3P2[0][3]=getParticleData(415);
  p1P1[1][3]=getParticleData(10423); p3P2[1][3]=getParticleData(425);
  p1P1[2][3]=getParticleData(10433); p3P2[2][3]=getParticleData(435);
  p1P1[3][3]=getParticleData(10443); p3P2[3][3]=getParticleData(445);
  p1P1[4][3]=getParticleData(10543); p3P2[4][3]=getParticleData(545);
  p1P1[0][4]=getParticleData(10513); p3P2[0][4]=getParticleData(515);
  p1P1[1][4]=getParticleData(10523); p3P2[1][4]=getParticleData(525);
  p1P1[2][4]=getParticleData(10533); p3P2[2][4]=getParticleData(535);
  p1P1[3][4]=getParticleData(10543); p3P2[3][4]=getParticleData(545);
  p1P1[4][4]=getParticleData(10553); p3P2[4][4]=getParticleData(555);
  // calculate the masses
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      Energy m1S0,m3S1,m3P0,m3P1,m3P2,m1P1;
      if(!p1S0[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 1S0 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m1S0 = ZERO;
      }
      else {
	m1S0 = p1S0[ix][iy]->mass();
      }
      if(!p3S1[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3S1 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3S1 = ZERO;
      }
      else {
	m3S1 = p3S1[ix][iy]->mass();
      }
      if(!p3P0[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3P0 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3P0 = ZERO;
      }
      else {
	m3P0 = p3P0[ix][iy]->mass();
      }
      if(!p3P1[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3P1 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3P1 = ZERO;
      }
      else {
	m3P1 = p3P1[ix][iy]->mass();
      }
      if(!p3P2[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 3P2 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m3P2 = ZERO;
      }
      else {
	m3P2 = p3P2[ix][iy]->mass();
      }
      if(!p1P1[ix][iy]) {
	generator()->log() << "Error in ISGW2FormFactor::doinit don't have "
			   << "ParticleData object for 1P1 " << ix << " " << iy 
			   << " setting mass to zero\n";
	m1P1 = ZERO;
      }
      else {
	m1P1 = p1P1[ix][iy]->mass();
      }
      // 1S0
      _mass1S0[ix][iy] = 0.75 *m3S1+0.25 *m1S0;
      //  1p 1/2
      _massPoh[ix][iy] = 0.75 *m3P1+0.25 *m3P0;
      //  1p 3/2
      _massPth[ix][iy] = 0.625*m3P2+0.375*m1P1;
    }
  }
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      Energy mtemp = (4.*_massPoh[ix][iy]+8.*_massPth[ix][iy])/12.;
      _massPoh[ix][iy]=mtemp;
      _massPth[ix][iy]=mtemp;
    }
  }
}

void ISGW2FormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mdown,GeV)  << ounit(_mup,GeV) << ounit(_mstrange,GeV) 
     << ounit(_mcharm,GeV) << ounit(_mbottom,GeV) << ounit(_beta1S0ud,GeV)
     << ounit(_beta1S0us,GeV) << ounit(_beta1S0ss,GeV) << ounit(_beta1S0cu,GeV) 
     << ounit(_beta1S0cs,GeV) << ounit(_beta1S0ub,GeV) << ounit(_beta1S0sb,GeV) 
     << ounit(_beta1S0cc,GeV) << ounit(_beta1S0bc,GeV) << ounit(_beta3S1ud,GeV) 
     << ounit(_beta3S1us,GeV) << ounit(_beta3S1ss,GeV) << ounit(_beta3S1cu,GeV) 
     << ounit(_beta3S1cs,GeV) << ounit(_beta3S1ub,GeV) << ounit(_beta3S1sb,GeV) 
     << ounit(_beta3S1cc,GeV) << ounit(_beta3S1bc,GeV) << ounit(_beta1Pud ,GeV) 
     << ounit(_beta1Pus ,GeV) << ounit(_beta1Pss ,GeV) << ounit(_beta1Pcu ,GeV) 
     << ounit(_beta1Pcs ,GeV) << ounit(_beta1Pub ,GeV) << ounit(_beta1Psb ,GeV) 
     << ounit(_beta1Pcc ,GeV) << ounit(_beta1Pbc ,GeV)
     << _alphamuQM  << _CfDrho << _CfDKstar << _CfDsKstar << _CfDsphi 
     << _CfBrho << _CfBDstar << _CfBsKstar << _CfBsDstar << _CfBcDstar << _CfBcpsi
     << _CfBcBsstar << _CfBcBstar << _thetaeta << ounit(_mquark,GeV) << _alphaQ 
     << ounit(_beta1S0,GeV) << ounit(_mass1S0,GeV) << ounit(_beta3S1,GeV) 
     << ounit(_beta1P,GeV) << ounit(_massPoh,GeV) << ounit(_massPth,GeV) << _includeaW;
}

void ISGW2FormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mdown,GeV) >> iunit(_mup,GeV) >> iunit(_mstrange,GeV) 
     >> iunit(_mcharm,GeV) >> iunit(_mbottom,GeV) >> iunit(_beta1S0ud,GeV) 
     >> iunit(_beta1S0us,GeV) >> iunit(_beta1S0ss,GeV) >> iunit(_beta1S0cu,GeV) 
     >> iunit(_beta1S0cs,GeV) >> iunit(_beta1S0ub,GeV) >> iunit(_beta1S0sb,GeV) 
     >> iunit(_beta1S0cc,GeV) >> iunit(_beta1S0bc,GeV) >> iunit(_beta3S1ud,GeV) 
     >> iunit(_beta3S1us,GeV) >> iunit(_beta3S1ss,GeV) >> iunit(_beta3S1cu,GeV) 
     >> iunit(_beta3S1cs,GeV) >> iunit(_beta3S1ub,GeV) >> iunit(_beta3S1sb,GeV) 
     >> iunit(_beta3S1cc,GeV) >> iunit(_beta3S1bc,GeV) >> iunit(_beta1Pud ,GeV) 
     >> iunit(_beta1Pus ,GeV) >> iunit(_beta1Pss ,GeV) >> iunit(_beta1Pcu ,GeV) 
     >> iunit(_beta1Pcs ,GeV) >> iunit(_beta1Pub ,GeV) >> iunit(_beta1Psb ,GeV) 
     >> iunit(_beta1Pcc ,GeV) >> iunit(_beta1Pbc ,GeV) 
     >> _alphamuQM >> _CfDrho >> _CfDKstar >> _CfDsKstar >> _CfDsphi 
     >> _CfBrho >> _CfBDstar >> _CfBsKstar >> _CfBsDstar >> _CfBcDstar >> _CfBcpsi
     >> _CfBcBsstar >> _CfBcBstar >> _thetaeta >> iunit(_mquark,GeV) >> _alphaQ 
     >> iunit(_beta1S0,GeV) >> iunit(_mass1S0,GeV) >> iunit(_beta3S1,GeV) 
     >> iunit(_beta1P,GeV) >> iunit(_massPoh,GeV) >> iunit(_massPth,GeV) >> _includeaW;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ISGW2FormFactor,ScalarFormFactor>
describeHerwigISGW2FormFactor("Herwig::ISGW2FormFactor", "HwFormFactors.so");

void ISGW2FormFactor::Init() {

  static ClassDocumentation<ISGW2FormFactor> documentation
    ("The ISGW2FormFactor class implements the ISGW2 model of "
     "PRD52, 2783.",
     "The ISGW2 form factor model \\cite{Scora:1995ty} was used.",
     "\\bibitem{Scora:1995ty} D.~Scora and N.~Isgur,"
     "Phys.\\ Rev.\\  D {\\bf 52} (1995) 2783 [arXiv:hep-ph/9503486].\n"
     "%%CITATION = PHRVA,D52,2783;%%\n");

  static Parameter<ISGW2FormFactor,Energy> interfaceDownMass
    ("DownMass",
     "The mass of the down quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mdown, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceUpMass
    ("UpMass",
     "The mass of the up quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mup, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mstrange, GeV, 0.55*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mcharm, GeV, 1.82*GeV, ZERO, 3.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark in the ISGW model (this is a consituent mass)",
     &ISGW2FormFactor::_mbottom, GeV, 5.20*GeV, 3.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ud
    ("Beta1S0ud",
     "The beta wavefunction parameter for the ud meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ud, GeV, 0.41*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0us
    ("Beta1S0us",
     "The beta wavefunction parameter for the us meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0us, GeV, 0.44*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ss
    ("Beta1S0ss",
     "The beta wavefunction parameter for the ss meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ss, GeV, 0.53*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cu
    ("Beta1S0cu",
     "The beta wavefunction parameter for the cu meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cu, GeV, 0.45*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cs
    ("Beta1S0cs",
     "The beta wavefunction parameter for the cs meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cs, GeV, 0.56*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0ub
    ("Beta1S0ub",
     "The beta wavefunction parameter for the ub meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0ub, GeV, 0.43*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0sb
    ("Beta1S0sb",
     "The beta wavefunction parameter for the sb meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0sb, GeV, 0.54*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0cc
    ("Beta1S0cc",
     "The beta wavefunction parameter for the cc meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0cc, GeV, 0.88*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1S0bc
    ("Beta1S0bc",
     "The beta wavefunction parameter for the bc meson in the 1 1S0 level",
     &ISGW2FormFactor::_beta1S0bc, GeV, 0.92*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pud
    ("Beta1Pud",
     "The beta wavefunction parameter for the ud meson in the 1P level",
     &ISGW2FormFactor::_beta1Pud, GeV, 0.28*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pus
    ("Beta1Pus",
     "The beta wavefunction parameter for the us meson in the 1P level",
     &ISGW2FormFactor::_beta1Pus, GeV, 0.30*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pss
    ("Beta1Pss",
     "The beta wavefunction parameter for the ss meson in the 1P level",
     &ISGW2FormFactor::_beta1Pss, GeV, 0.33*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcu
    ("Beta1Pcu",
     "The beta wavefunction parameter for the cu meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcu, GeV, 0.33*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcs
    ("Beta1Pcs",
     "The beta wavefunction parameter for the cs meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcs, GeV, 0.38*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pub
    ("Beta1Pub",
     "The beta wavefunction parameter for the ub meson in the 1P level",
     &ISGW2FormFactor::_beta1Pub, GeV, 0.35*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Psb
    ("Beta1Psb",
     "The beta wavefunction parameter for the sb meson in the 1P level",
     &ISGW2FormFactor::_beta1Psb, GeV, 0.41*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pcc
    ("Beta1Pcc",
     "The beta wavefunction parameter for the cc meson in the 1P level",
     &ISGW2FormFactor::_beta1Pcc, GeV, 0.52*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta1Pbc
    ("Beta1Pbc",
     "The beta wavefunction parameter for the bc meson in the 1P level",
     &ISGW2FormFactor::_beta1Pbc, GeV, 0.60*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ud
    ("Beta3S1ud",
     "The beta wavefunction parameter for the ud meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ud, GeV, 0.30*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1us
    ("Beta3S1us",
     "The beta wavefunction parameter for the us meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1us, GeV, 0.33*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ss
    ("Beta3S1ss",
     "The beta wavefunction parameter for the ss meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ss, GeV, 0.37*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cu
    ("Beta3S1cu",
     "The beta wavefunction parameter for the cu meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cu, GeV, 0.38*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cs
    ("Beta3S1cs",
     "The beta wavefunction parameter for the cs meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cs, GeV, 0.44*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1ub
    ("Beta3S1ub",
     "The beta wavefunction parameter for the ub meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1ub, GeV, 0.40*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1sb
    ("Beta3S1sb",
     "The beta wavefunction parameter for the sb meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1sb, GeV, 0.49*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1cc
    ("Beta3S1cc",
     "The beta wavefunction parameter for the cc meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1cc, GeV, 0.62*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,Energy> interfaceBeta3S1bc
    ("Beta3S1bc",
     "The beta wavefunction parameter for the bc meson in the 3S1 level",
     &ISGW2FormFactor::_beta3S1bc, GeV, 0.75*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceAlphaCutOff
    ("AlphaCutOff",
     "The value of the strong coupling constnat at the cut-off",
     &ISGW2FormFactor::_alphamuQM, 0.6, 0.0, 10.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDrho
    ("CfDrho",
     "The relativistic correction factor for D -> rho",
     &ISGW2FormFactor::_CfDrho, 0.889, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDKstar
    ("CfDKstar",
     "The relativistic correction factor for D -> Kstar",
     &ISGW2FormFactor::_CfDKstar, 0.928, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDsKstar
    ("CfDsKstar",
     "The relativistic correction factor for Ds -> Kstar",
     &ISGW2FormFactor::_CfDsKstar, 0.873, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfDsphi
    ("CfDsphi",
     "The relativistic correction factor for Ds -> phi",
     &ISGW2FormFactor::_CfDsphi, 0.911, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBrho
    ("CfBrho",
     "The relativistic correction factor for B -> rho",
     &ISGW2FormFactor::_CfBrho, 0.905, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBDstar
    ("CfBDstar",
     "The relativistic correction factor for B -> Dstar",
     &ISGW2FormFactor::_CfBDstar, 0.989, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBsKstar
    ("CfBsKstar",
     "The relativistic correction factor for Bs -> Kstar",
     &ISGW2FormFactor::_CfBsKstar, 0.892, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBsDstar
    ("CfBsDstar",
     "The relativistic correction factor for Bs -> Dstar",
     &ISGW2FormFactor::_CfBsDstar, 0.984, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcDstar
    ("CfBcDstar",
     "The relativistic correction factor for Bc -> Dstar",
     &ISGW2FormFactor::_CfBcDstar, 0.868, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcpsi
    ("CfBcpsi",
     "The relativistic correction factor for Bc -> psi",
     &ISGW2FormFactor::_CfBcpsi, 0.967, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcBsstar
    ("CfBcBsstar",
     "The relativistic correction factor for Bc -> Bsstar",
     &ISGW2FormFactor::_CfBcBsstar, 1.000, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceCfBcBstar
    ("CfBcBstar",
     "The relativistic correction factor for Bc -> Bstar",
     &ISGW2FormFactor::_CfBcBstar, 1.000, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGW2FormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &ISGW2FormFactor::_thetaeta, -Constants::pi/9., -Constants::pi, Constants::pi,
     false, false, true);

  static Switch<ISGW2FormFactor,bool> interfaceIncludeaW
    ("IncludeaW",
     "Include the a(omega) piece of the Cji factor",
     &ISGW2FormFactor::_includeaW, true, false, false);
  static SwitchOption interfaceIncludeaWInclude
    (interfaceIncludeaW,
     "Yes",
     "Include the factor",
     true);
  static SwitchOption interfaceIncludeaWDoNot
    (interfaceIncludeaW,
     "No",
     "Do not include the factor",
     false);

}

// member which does the work
void ISGW2FormFactor::formFactor(Energy2 q2, unsigned int iloc, int, int id1,
				 Energy mY,
				 Energy mX, Complex & f1,Complex & f2,Complex & f3,
				 Complex & f4) const {
  useMe();
  // get the flavours of the quarks etc
  int jspin,spect,inquark,outquark;
  formFactorInfo(iloc,jspin,spect,inquark,outquark);
  int ifl0(abs(inquark)),ifl1(abs(outquark)),ifls(abs(spect));
  // determine the multiplet
  int ispin(abs(id1)/1000);
  // masses of the quarks
  Energy mQ(_mquark[ifl0-1]),mq(_mquark[ifl1-1]),ms(_mquark[ifls-1]);
  // of the mesons
  Energy mtildeX(mq+ms),mtildeY(mQ+ms),mup(mq*mQ/(mQ+mq)),mum(mq*mQ/(mQ-mq));
  // wavefunction parameters for the mesons
  Energy betaX(ZERO),mbarX(ZERO),
    betaY(_beta1S0[ifl0-1][ifls-1]),mbarY(_mass1S0[ifl0-1][ifls-1]);
  double Cf(1.);
  // the wavefunction parameter for the outgoing meson
  // 1S0
  if(ispin==0&&jspin==0) {
    betaX=_beta1S0[ifl1-1][ifls-1];
    mbarX=_mass1S0[ifl1-1][ifls-1];
  }
  // 3S1
  else if(ispin==0&&jspin==1) {
    betaX = _beta3S1[ifl1-1][ifls-1];
    mbarX = _mass1S0[ifl1-1][ifls-1];
    // set the relativistic correction parameter
    // decaying b
    if(ifl0==5) {
      if(ifls<3)       Cf = ifl1<3  ? _CfBrho    : _CfBDstar;
      else if(ifls==3) Cf = ifl1==4 ? _CfBsDstar : _CfBsKstar;
      else if(ifls==4) Cf = ifl1==4 ? _CfBcpsi   : _CfBcDstar;
    }
    // decaying D
    else if(ifl0==4) {
      if(ifls<3)       Cf = ifl1<3  ? _CfDrho    : _CfDKstar;
      else if(ifls==3) Cf = ifl1<3  ? _CfDsKstar : _CfDsphi;
      else if(ifls==5) Cf = ifl1<3  ? _CfBcBstar : _CfBcBsstar;
    } 
  }
  else if(ispin==10&&jspin==0) {
    betaX=_beta1P[ifl1-1][ifls-1];
    mbarX=_massPoh[ifl1-1][ifls-1];
  }
  // 1 3/2 P 1 (1 P1)
  else if((ispin==0&&jspin==2)||(ispin==10&&jspin==1)) {
    betaX = _beta1P[ifl1-1][ifls-1];
    mbarX=_massPth[ifl1-1][ifls-1];
  }
  // 1 1/2 P1 ( 3 P1) 
  else if(ispin==20&&jspin==1) {
    betaX = _beta1P[ifl1-1][ifls-1];
    mbarX=_massPoh[ifl1-1][ifls-1];
  }
  else {
    throw Exception() << "ISGWS2FormFactor::formFactor" 
		      << " unknown multiplet" << Exception::abortnow;
  }
  Energy2 beta2XY(0.5*(betaX*betaX+betaY*betaY));
  // number of active flavours
  int Nf  = ifl0-1;
  int Nfp = ifl1==2 ? 0 : ifl1-1;
  // first piece of the f_n function
  double betar(betaX*betaY/beta2XY),fn(sqrt(mtildeX/mtildeY)*betar*sqrt(betar));
  // q dependent piece
  Energy2 tm((mY-mX)*(mY-mX)),tmmt(tm-q2);
  // radius parameter
  InvEnergy2 r2(0.75/mQ/mq+1.5*ms*ms/mbarX/mbarY/beta2XY
		+16./mbarX/mbarY/(33.-2.*Nfp)*log(_alphamuQM/_alphaQ[ifl1-1]));
  // the parameters for the form-factor depenedent piece
  double rmbmtY(sqrt(mbarY/mtildeY)),rmbmtX(sqrt(mbarX/mtildeX));
  // work out wtilde
  double wt(1.+0.5*tmmt/mbarX/mbarY);
  // storage of the form factors
  Energy f(ZERO);
  InvEnergy g(ZERO),appam(ZERO),apmam(ZERO);
  InvEnergy2 h(ZERO),bp(ZERO),bm(ZERO);
  double fpmfm(0.),fppfm(0.),k(0.);
  // scalar and vector from 1S levels
  if(ispin==0) {
    // parameters for the beta functions
    double asopi(alphaS(mq,mq*mQ)/Constants::pi),w(1.+0.5*tmmt/mX/mY);
    double aI(-6./(33.-2.*Nf)),rw(1./sqrt(w*w-1)*log(w+sqrt(w*w-1.)));
    double aLw(8./(33.-2.*Nfp)*(w*rw-1.)); 
    double cji(pow(_alphaQ[ifl0-1]/_alphaQ[ifl1-1],aI));
    if(_includeaW) cji*=pow(_alphaQ[ifl1-1]/_alphamuQM,aLw);
    double zji(mq/mQ); 
    double gamji(-2.*zji/(1.-zji)*log(zji)-2.),chiji(-1.-gamji/(1.-zji));
    // scalar
    if(jspin==0) {
      double fact((1.+1./12.*r2*tmmt));
      fn/=(fact*fact);
      fact = (1.-0.5*ms*mq/mup/mtildeX*betaY*betaY/beta2XY);
      fppfm = fn*rmbmtX/rmbmtY*cji*(1.+asopi*(gamji-2./3.*chiji))*
	(2.-mtildeX/mq*fact);
      fpmfm = fn*rmbmtY/rmbmtX*cji*(1.+asopi*(gamji+2./3.*chiji))*mtildeY/mq*fact;
    }
    else if(jspin==1) {
      // factors for the F and R functions
      double fact((1.+1./12.*r2*tmmt));
      fn/=(fact*fact);
      double betaapmam=1./3.-4./3./(1-zji)-chiji
	+gamji*(1.-2./3.*(1.+zji)/(1.-zji)/(1.-zji));
      double ftemp  = Cf*fn*rmbmtX*rmbmtY*cji*(1.+asopi*(-2./3.+gamji));
      double gtemp  = fn/rmbmtX/rmbmtY*cji*(1.+asopi*( 2./3.+gamji));
      double aptemp = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY)*cji;
      double amtemp = fn/rmbmtX/rmbmtY*cji*(1.+asopi*betaapmam);
      // rest of the calculation
      f     =    ftemp*mtildeY*(1.+wt+0.5*ms*(wt-1.)/mup);
      g     =0.5*gtemp*(1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY);
      appam = aptemp*(ms*betaX*betaX/(1.+wt)/mq/mQ/beta2XY*
		      (1.-0.5*ms*betaX*betaX/mtildeY/beta2XY)
		      +asopi/mtildeY*(-1.-chiji+4./3./(1.-zji)
				      +2./3.*(1.+zji)/sqr(1.-zji)*gamji));
      apmam =-amtemp/mtildeX*(mtildeY/mQ
			      -0.5*ms*betaX*betaX/mup/beta2XY
			      +wt*ms*mtildeY*betaX*betaX/(wt+1.)/mq/mQ/beta2XY*
			      (1.-0.5*ms/mtildeY*betaX*betaX/beta2XY)); 
    }
    else if(jspin==2) {
      // factors for the F function
      double fact((1.+1./18.*r2*tmmt));
      fn*=betar/(fact*fact*fact);
      double htemp = fn/rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
      double ktemp = fn*rmbmtX/rmbmtY;
      double bptemp(fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY*rmbmtY*rmbmtY));
      double bmtemp(fn/rmbmtX/(rmbmtY*rmbmtY*rmbmtY));
      // functions themselves
      double or2(sqrt(0.5));
      h = 0.5*htemp*ms*or2/mtildeY/betaY*(1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY);
      k = or2*ktemp*ms/betaY*(1.+wt);
      InvEnergy2 bppbm = 0.25*bptemp*or2*ms*ms/mq/mQ/mtildeY/betaY*betaX*betaX/beta2XY*
	(1.-0.5*ms/mtildeY*betaX*betaX/beta2XY);
      InvEnergy2 bpmbm = -or2*bmtemp*ms/mQ/mtildeX/betaY*
	(1.-0.5*ms*mQ/mup/mtildeY*betaX*betaX/beta2XY
	 +0.25*ms/mq*betaX*betaX/beta2XY*(1.-0.5*ms/mtildeY*betaX*betaX/beta2XY));
      // conversion
      bp = 0.5*(bppbm+bpmbm);
      bm = 0.5*(bppbm-bpmbm);
    }
  }
  // 1 3P0
  else if(ispin==10&&jspin==0) {
    fn*=betar;
    double fact=(1.+1./18.*r2*tmmt);
    fn/=(fact*fact*fact);
    fn *= sqrt(2./3.)*ms/betaY;
    fppfm =-fn*rmbmtX/rmbmtY;
    fpmfm = fn*rmbmtY/rmbmtX*mtildeY/mtildeX;
  }
  // 1 3/2 P1 ( 1 P1) 
  else if(ispin==10&&jspin==1) {
    // factors for the F and R functions
    double fact=(1.+1./18.*r2*tmmt);
    fn*=betar/(fact*fact*fact);
    double ftemp  = fn*rmbmtX*rmbmtY;
    double gtemp  = fn/rmbmtX/rmbmtY;
    double aptemp = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
    double amtemp = fn/rmbmtX/rmbmtY;
    // light meson or onium
    if((ifls<3&&ifl1<3)||(ifls==ifl1)) {
      double oor2(sqrt(0.5));
      f     = oor2*ftemp*mtildeY*betaY*(1./mup
				       +ms*mtildeX/3./mq/betaY/betaY*(wt-1.)*(wt-1.));
      g     = oor2*gtemp *(0.25*mtildeY*betaY/mQ/mq/mtildeX+(wt-1.)*ms/6./mtildeX/betaY);
      appam = oor2*aptemp*ms/mtildeY/betaY*(1.-ms/mq+0.5*ms/mup*betaY*betaY/beta2XY);
      apmam = oor2*amtemp*ms/mq/betaY*((4.-wt)/3.
				       -0.5*ms*mq/mtildeX/mup*betaY*betaY/beta2XY);
    }
    // heavy meson
    else {
      double oor3(1./sqrt(3.));
      f     =-2.*ftemp*oor3*mtildeY*betaY*
	(1./mq+0.5*mtildeX*ms*(wt-1.)/betaY/betaY*
	 (0.5*(wt+1.)/mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY)); 
      g     =-0.5*gtemp*oor3*(0.5*(1.+wt)+0.5*betaY*betaY*mtildeY/ms/mq/mQ)*ms
	/betaY/mtildeX;
      appam =-0.5*aptemp/oor3*ms/betaY/mtildeY*
	(1.-ms/3./mq-ms/3.*betaY*betaY/beta2XY*
	 (0.5/mum-1./mup));
      apmam =-0.5*amtemp*oor3*ms/betaY/mtildeX*
	((2.-wt)*mtildeX/mq+ms*betaY*betaY/beta2XY*(0.5/mum-1./mup));
    }
  }
  // 1 1/2 P 1 (3 P1)
  else if(ispin==20&&jspin==1) {
    // factors for the F and R functions
    double fact=(1.+1./18.*r2*tmmt);
    fn*=betar/(fact*fact*fact);
    double ftemp  = fn*rmbmtX*rmbmtY;
    double gtemp  = fn/rmbmtX/rmbmtY;
    double aptemp = fn*rmbmtX/(rmbmtY*rmbmtY*rmbmtY);
    double amtemp = fn/rmbmtX/rmbmtY;
    // light meson
    if( ( ifls<3 && ifl1<3 ) || ( ifl1==ifls ) ) {
      f     =-ftemp *mtildeY*betaY*(1./mum
				   +ms*mtildeX*(wt-1.)/betaY/betaY*
				   ((5.+wt)/6./mq-0.5/mum*ms/mtildeX*
				    betaY*betaY/beta2XY)); 
      g     =-gtemp *0.5*ms/mtildeX/betaY*(5.+wt)/6.;
      appam =-aptemp*0.5*ms*mtildeX/mq/mtildeY/betaY*
	(1.-0.5*ms*mq/mtildeX/mum*betaY*betaY/beta2XY); 
      apmam = amtemp*0.5*ms/mq/betaY*((wt+2.)/3.
				      -0.5*ms*mq/mtildeX/mum*betaY*betaY/beta2XY);
    }
    // heavy meson
    else {
      double r2o3(sqrt(2./3.));
      f     =  ftemp*r2o3*mtildeY*betaY*(0.5/mq-1.5/mQ+ms*mtildeX*(wt-1.)/betaY/betaY*
					 (1./mq-0.5*ms*betaY*betaY/mum/mtildeX/beta2XY));
      g     =  gtemp *0.5*r2o3*ms/betaY/mtildeX*(1.-0.25*betaY*betaY*mtildeY/ms/mq/mQ);
      appam =  aptemp*0.5*r2o3*ms*ms*betaX*betaX/mtildeY/mq/betaY/beta2XY; 
      apmam = -amtemp*r2o3*ms/mtildeX/betaY*(1.+0.5*ms*betaX*betaX/mq/beta2XY);
    }
  }
  else {
    throw Exception() << "ISGWS2FormFactor::formFactor" 
		      << " unknown multiplet" << Exception::abortnow;
  }
  // the final manipulations
  if(jspin==0) {
    double fp,fm;
    fp = 0.5*(fppfm+fpmfm);
    fm = 0.5*(fppfm-fpmfm);
    // convert to the standard form
    f1 = q2/(mY*mY-mX*mX)*fm+fp;
    f2 = fp;
  }
  else if(jspin==1) {
    InvEnergy ap(0.5*(appam+apmam)),am(0.5*(appam-apmam));
    // convert to the standard notation
    Energy msum(mX+mY),mdiff(mY-mX);
    Complex ii(0.,1.);
    f2 = -ii*f/msum;
    f3 = +Complex(ii*ap*msum);
    f1 = -Complex(ii*0.5/mX*(am*q2+ii*msum*f2-ii*mdiff*f3));
    f4 =  Complex(ii*g*msum);
  }
  else if(jspin==2) {
    Energy msum(mX+mY);
    f1 = h*sqr(msum);
    f2 = k;
    f3 = bp*sqr(msum);
    f4 = bm*sqr(msum);
  }
  // special for mixing
  double fact;
  if(id1==ParticleID::eta) {
    if(ifl1==3&&ifls==3){fact=-2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
    else{fact=cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);}
    f1*=fact;f2*=fact;f3*=fact;f4*=fact;
  }
  else if(id1==ParticleID::etaprime) {
    if(ifl1==3&&ifls==3){fact=-2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
    else{fact=sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);}
    f1*=fact;f2*=fact;f3*=fact;f4*=fact;
  }
  else if(ifl1==ifls&&ifl1<3) {
    if(abs(ifl1)==1&&int(id1/10)%100==1){fact=-sqrt(0.5);}
    else{fact=sqrt(0.5);}
    f1*=fact;f2*=fact;f3*=fact;f4*=fact;
  }
}

// form-factor for scalar to scalar
void ISGW2FormFactor::ScalarScalarFormFactor(Energy2 q2, unsigned int iloc,int id0,
					     int id1, Energy mY, Energy mX,
					     Complex & f0,Complex & fp) const {
  Complex d1(0.),d2(0.);
  formFactor(q2,iloc,id0,id1,mY,mX,f0,fp,d1,d2);
}

// form-factor for scalar to vector
void ISGW2FormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0, 
					     int id1, Energy mY, Energy mX,
					     Complex & A0,Complex & A1,
					     Complex & A2,Complex & V) const {
  formFactor(q2,iloc,id0,id1,mY,mX,A0,A1,A2,V);
}


// form-factor for scalar to tensor
void ISGW2FormFactor::
ScalarTensorFormFactor(Energy2 q2, unsigned int iloc, int id0, int id1, 
		       Energy mY, Energy mX, complex<InvEnergy2> & h,
		       Complex & k,complex<InvEnergy2> & bp,
		       complex<InvEnergy2> & bm) const {
  Complex f1,f2,f3,f4;
  formFactor(q2,iloc,id0,id1,mY,mX,f1,f2,f3,f4);
  Energy msum(mX+mY);
  h = f1/sqr(msum);
  k = f2;
  bp = f3/sqr(msum);
  bm = f4/sqr(msum);
}

void ISGW2FormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ISGW2FormFactor " << name() << "\n";
  output << "newdef " << name() << ":DownMass "    << _mdown/GeV    << "\n";
  output << "newdef " << name() << ":UpMass "      << _mup/GeV      << "\n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << "\n";
  output << "newdef " << name() << ":CharmMass "   << _mcharm/GeV   << "\n";
  output << "newdef " << name() << ":BottomMass "  << _mbottom/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0ud " << _beta1S0ud/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0us " << _beta1S0us/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0ss " << _beta1S0ss/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0cu " << _beta1S0cu/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0cs " << _beta1S0cs/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0ub " << _beta1S0ub/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0sb " << _beta1S0sb/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0cc " << _beta1S0cc/GeV  << "\n";
  output << "newdef " << name() << ":Beta1S0bc " << _beta1S0bc/GeV  << "\n";
  output << "newdef " << name() << ":Beta1Pud  " << _beta1Pud/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pus  " << _beta1Pus/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pss  " << _beta1Pss/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pcu  " << _beta1Pcu/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pcs  " << _beta1Pcs/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pub  " << _beta1Pub/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Psb  " << _beta1Psb/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pcc  " << _beta1Pcc/GeV   << "\n";
  output << "newdef " << name() << ":Beta1Pbc  " << _beta1Pbc/GeV   << "\n";
  output << "newdef " << name() << ":Beta3S1ud " << _beta3S1ud/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1us " << _beta3S1us/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1ss " << _beta3S1ss/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1cu " << _beta3S1cu/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1cs " << _beta3S1cs/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1ub " << _beta3S1ub/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1sb " << _beta3S1sb/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1cc " << _beta3S1cc/GeV  << "\n";
  output << "newdef " << name() << ":Beta3S1bc " << _beta3S1bc/GeV  << "\n";
  output << "newdef " << name() << ":AlphaCutOff " << _alphamuQM    << "\n";
  output << "newdef " << name() << ":CfDrho "      << _CfDrho       << "\n";
  output << "newdef " << name() << ":CfDKstar "    << _CfDKstar     << "\n";
  output << "newdef " << name() << ":CfDsKstar "   << _CfDsKstar    << "\n";
  output << "newdef " << name() << ":CfDsphi "     << _CfDsphi      << "\n";
  output << "newdef " << name() << ":CfBrho "      << _CfBrho       << "\n";
  output << "newdef " << name() << ":CfBDstar "    << _CfBDstar     << "\n";
  output << "newdef " << name() << ":CfBsKstar "   << _CfBsKstar    << "\n";
  output << "newdef " << name() << ":CfBsDstar "   << _CfBsDstar    << "\n";
  output << "newdef " << name() << ":CfBcDstar "   << _CfBcDstar    << "\n";
  output << "newdef " << name() << ":CfBcpsi "     << _CfBcpsi      << "\n";
  output << "newdef " << name() << ":CfBcBsstar "  << _CfBcBsstar   << "\n";
  output << "newdef " << name() << ":CfBcBstar "   << _CfBcBstar    << "\n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./ISGWFormFactor.cc"
// -*- C++ -*-
//
// ISGWFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISGWFormFactor class.
//

#include "ISGWFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ISGWFormFactor::ISGWFormFactor() {
  // default values of the parameters
  // fudge factor
  _kappa=0.7;
  // quark masses
  _mdown   =0.33*GeV;
  _mup     =0.33*GeV;
  _mstrange=0.55*GeV;
  _mcharm  =1.82*GeV;
  _mbottom =5.12*GeV;
  // beta values
  _betaSud = 0.31*GeV;
  _betaSus = 0.34*GeV;
  _betaSuc = 0.39*GeV;
  _betaSub = 0.41*GeV;
  _betaPud = 0.27*GeV;
  _betaPus = 0.30*GeV;
  _betaPuc = 0.34*GeV;
  // the mixing for eta/eta'
  _thetaeta=-Constants::pi/18.;
  // B_u decays to d ubar
  addFormFactor(-521,-211  ,0,-2, 5, 1);
  addFormFactor(-521,-213  ,1,-2, 5, 1);
  addFormFactor(-521,-215  ,2,-2, 5, 1);
  addFormFactor(-521,-10213,1,-2, 5, 1);
  addFormFactor(-521,-20213,1,-2, 5, 1);
  addFormFactor(-521,-10211,0,-2, 5, 1);
  // B_u to uu (I=0)
  addFormFactor(-521, 221  ,0,-2, 5, 2);
  addFormFactor(-521, 331  ,0,-2, 5, 2);
  addFormFactor(-521, 223  ,1,-2, 5, 2);
  addFormFactor(-521, 225  ,2,-2, 5, 2);
  addFormFactor(-521, 10223,1,-2, 5, 2);
  addFormFactor(-521, 20223,1,-2, 5, 2);
  addFormFactor(-521, 10221,0,-2, 5, 2);
  // B_u to uu (I=1)
  addFormFactor(-521, 111  ,0,-2, 5, 2);
  addFormFactor(-521, 113  ,1,-2, 5, 2);
  addFormFactor(-521, 115  ,2,-2, 5, 2);
  addFormFactor(-521, 10113,1,-2, 5, 2);
  addFormFactor(-521, 20113,1,-2, 5, 2);
  addFormFactor(-521, 10111,0,-2, 5, 2);
  // B_u decays to d ubar
  addFormFactor(-521,-321  ,0,-2, 5, 3);
  addFormFactor(-521,-323  ,1,-2, 5, 3);
  addFormFactor(-521,-325  ,2,-2, 5, 3);
  addFormFactor(-521,-10323,1,-2, 5, 3);
  addFormFactor(-521,-20323,1,-2, 5, 3);
  addFormFactor(-521,-10321,0,-2, 5, 3);
  // B_u decays to c ubar
  addFormFactor(-521, 421  ,0,-2, 5, 4);
  addFormFactor(-521, 423  ,1,-2, 5, 4);
  addFormFactor(-521, 425  ,2,-2, 5, 4);
  addFormFactor(-521, 10423,1,-2, 5, 4);
  addFormFactor(-521, 20423,1,-2, 5, 4);
  addFormFactor(-521, 10421,0,-2, 5, 4);
  // B_d to d dbar (I=0)
  addFormFactor(-511, 221  ,0, 1,-5,-1);
  addFormFactor(-511, 331  ,0, 1,-5,-1);
  addFormFactor(-511, 223  ,1, 1,-5,-1);
  addFormFactor(-511, 225  ,2, 1,-5,-1);
  addFormFactor(-511, 10223,1, 1,-5,-1);
  addFormFactor(-511, 20223,1, 1,-5,-1);
  addFormFactor(-511, 10221,0, 1,-5,-1);
  // B_d to d dbar (I=1)
  addFormFactor(-511, 111  ,0, 1,-5,-1);
  addFormFactor(-511, 113  ,1, 1,-5,-1);
  addFormFactor(-511, 115  ,2, 1,-5,-1);
  addFormFactor(-511, 10113,1, 1,-5,-1);
  addFormFactor(-511, 20113,1, 1,-5,-1);
  addFormFactor(-511, 10111,0, 1,-5,-1);
  // B_d to u dbar
  addFormFactor(-511, 211  ,0, 1,-5,-2);
  addFormFactor(-511, 213  ,1, 1,-5,-2);
  addFormFactor(-511, 215  ,2, 1,-5,-2);
  addFormFactor(-511, 10213,1, 1,-5,-2);
  addFormFactor(-511, 20213,1, 1,-5,-2);
  addFormFactor(-511, 10211,0, 1,-5,-2);
  // B_d to s dbar 
  addFormFactor(-511, 311  ,0, 1,-5,-3);
  addFormFactor(-511, 313  ,1, 1,-5,-3);
  addFormFactor(-511, 315  ,2, 1,-5,-3);
  addFormFactor(-511, 10313,1, 1,-5,-3);
  addFormFactor(-511, 20313,1, 1,-5,-3);
  addFormFactor(-511, 10311,0, 1,-5,-3);
  // B_d decays to  c dbar
  addFormFactor(-511, 411  ,0, 1,-5,-4);
  addFormFactor(-511, 413  ,1, 1,-5,-4);
  addFormFactor(-511, 415  ,2, 1,-5,-4);
  addFormFactor(-511, 10413,1, 1,-5,-4);
  addFormFactor(-511, 20413,1, 1,-5,-4);
  addFormFactor(-511, 10411,0, 1,-5,-4);
  // D0 to d ubar
  addFormFactor( 421,-211  ,0,-2, 4, 1);
  addFormFactor( 421,-213  ,1,-2, 4, 1);
  addFormFactor( 421,-215  ,2,-2, 4, 1);
  addFormFactor( 421,-10213,1,-2, 4, 1);
  addFormFactor( 421,-20213,1,-2, 4, 1);
  addFormFactor( 421,-10211,0,-2, 4, 1);
  // D0 to d ubar (I=1)
  addFormFactor( 421, 111  ,0,-2, 4, 2);
  addFormFactor( 421, 113  ,1,-2, 4, 2);
  addFormFactor( 421, 115  ,2,-2, 4, 2);
  addFormFactor( 421, 10113,1,-2, 4, 2);
  addFormFactor( 421, 20113,1,-2, 4, 2);
  addFormFactor( 421, 10111,0,-2, 4, 2);
  // D0 to d ubar (I=0)
  addFormFactor( 421, 221  ,0,-2, 4, 2);
  addFormFactor( 421, 331  ,0,-2, 4, 2);
  addFormFactor( 421, 223  ,1,-2, 4, 2);
  addFormFactor( 421, 225  ,2,-2, 4, 2);
  addFormFactor( 421, 10223,1,-2, 4, 2);
  addFormFactor( 421, 20223,1,-2, 4, 2);
  addFormFactor( 421, 10221,0,-2, 4, 2);
  // D0 to s ubar
  addFormFactor( 421,-321  ,0,-2, 4, 3);
  addFormFactor( 421,-323  ,1,-2, 4, 3);
  addFormFactor( 421,-325  ,2,-2, 4, 3);
  addFormFactor( 421,-10323,1,-2, 4, 3);
  addFormFactor( 421,-20323,1,-2, 4, 3);
  addFormFactor( 421,-10321,0,-2, 4, 3);
  // D+ to d dbar I=0
  addFormFactor( 411, 221  ,0,-1, 4, 1);
  addFormFactor( 411, 331  ,0,-1, 4, 1);
  addFormFactor( 411, 223  ,1,-1, 4, 1);
  addFormFactor( 411, 225  ,2,-1, 4, 1); 
  addFormFactor( 411, 10223,1,-1, 4, 1); 
  addFormFactor( 411, 20223,1,-1, 4, 1); 
  addFormFactor( 411, 10221,0,-1, 4, 1);
  // D+ to d dbar I=1
  addFormFactor( 411, 111  ,0,-1, 4, 1);
  addFormFactor( 411, 113  ,1,-1, 4, 1); 
  addFormFactor( 411, 115  ,2,-1, 4, 1); 
  addFormFactor( 411, 10113,1,-1, 4, 1); 
  addFormFactor( 411, 20113,1,-1, 4, 1); 
  addFormFactor( 411, 10111,0,-1, 4, 1);
  // D+ to u dbar
  addFormFactor( 411, 211  ,0,-1, 4, 2);
  addFormFactor( 411, 213  ,1,-1, 4, 2);
  addFormFactor( 411, 215  ,2,-1, 4, 2);
  addFormFactor( 411, 10213,1,-1, 4, 2);
  addFormFactor( 411, 20213,1,-1, 4, 2);
  addFormFactor( 411, 10211,0,-1, 4, 2);
  // D+ to s dbar
  addFormFactor( 411,-311  ,0,-1, 4, 3);
  addFormFactor( 411,-313  ,1,-1, 4, 3);
  addFormFactor( 411,-315  ,2,-1, 4, 3);
  addFormFactor( 411,-10313,1,-1, 4, 3);
  addFormFactor( 411,-20313,1,-1, 4, 3);
  addFormFactor( 411,-10311,0,-1, 4, 3);
  // set the initial number of modes
  initialModes(numberOfFactors());
}

void ISGWFormFactor::doinit() {
  ScalarFormFactor::doinit();
  // set up the quark masses
  _mquark.resize(5);
  _mquark[0]=_mdown;
  _mquark[1]=_mup;
  _mquark[2]=_mstrange;
  _mquark[3]=_mcharm;
  _mquark[4]=_mbottom;
  // and the beta values
  _betaS.resize(5,vector<Energy>(5));
  _betaP.resize(5,vector<Energy>(5));
  _betaS[0][0] = _betaSud;_betaP[0][0] = _betaPud;
  _betaS[1][0] = _betaSud;_betaP[1][0] = _betaPud;
  _betaS[2][0] = _betaSus;_betaP[2][0] = _betaPus;
  _betaS[3][0] = _betaSuc;_betaP[3][0] = _betaPuc;
  _betaS[4][0] = _betaSub;_betaP[4][0] = ZERO  ;
  _betaS[0][1] = _betaSud;_betaP[0][1] = _betaPud;
  _betaS[1][1] = _betaSud;_betaP[1][1] = _betaPud;
  _betaS[2][1] = _betaSus;_betaP[2][1] = _betaPus;
  _betaS[3][1] = _betaSuc;_betaP[3][1] = _betaPuc;
  _betaS[4][1] = _betaSub;_betaP[4][1] = ZERO  ;
  _betaS[0][2] = _betaSus;_betaP[0][2] = _betaPus;
  _betaS[1][2] = _betaSus;_betaP[1][2] = _betaPus;
  _betaS[2][2] = ZERO  ;_betaP[2][2] = ZERO  ;
  _betaS[3][2] = ZERO  ;_betaP[3][2] = ZERO  ;
  _betaS[4][2] = ZERO  ;_betaP[4][2] = ZERO  ;
  _betaS[0][3] = _betaSuc;_betaP[0][3] = _betaPuc;
  _betaS[1][3] = _betaSuc;_betaP[1][3] = _betaPuc;
  _betaS[2][3] = ZERO  ;_betaP[2][3] = ZERO  ;
  _betaS[3][3] = ZERO  ;_betaP[3][3] = ZERO  ;
  _betaS[4][3] = ZERO  ;_betaP[4][3] = ZERO  ;
  _betaS[0][4] = ZERO  ;_betaP[0][4] = ZERO  ;
  _betaS[1][4] = ZERO  ;_betaP[1][4] = ZERO  ;
  _betaS[2][4] = ZERO  ;_betaP[2][4] = ZERO  ;
  _betaS[3][4] = ZERO  ;_betaP[3][4] = ZERO  ;
  _betaS[4][4] = ZERO  ;_betaP[4][4] = ZERO  ;
}

void ISGWFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _kappa << ounit(_mdown,GeV) << ounit(_mup,GeV) << ounit(_mstrange,GeV) 
     << ounit(_mcharm,GeV) << ounit(_mbottom,GeV) << ounit(_betaSud,GeV) 
     << ounit(_betaSus,GeV) << ounit(_betaSuc,GeV) << ounit(_betaSub,GeV) 
     << ounit(_betaPud,GeV) << ounit(_betaPus,GeV) << ounit(_betaPuc,GeV)
     << _thetaeta << ounit(_mquark,GeV) << ounit(_betaS,GeV) << ounit(_betaP,GeV);
}

void ISGWFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _kappa >> iunit(_mdown,GeV) >> iunit(_mup,GeV) >> iunit(_mstrange,GeV) 
     >> iunit(_mcharm,GeV) >> iunit(_mbottom,GeV) >> iunit(_betaSud,GeV) 
     >> iunit(_betaSus,GeV) >> iunit(_betaSuc,GeV) >> iunit(_betaSub,GeV) 
     >> iunit(_betaPud,GeV) >> iunit(_betaPus,GeV) >> iunit(_betaPuc,GeV)
     >> _thetaeta >> iunit(_mquark,GeV) >> iunit(_betaS,GeV) >> iunit(_betaP,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ISGWFormFactor,ScalarFormFactor>
describeHerwigISGWFormFactor("Herwig::ISGWFormFactor", "HwFormFactors.so");

void ISGWFormFactor::Init() {

  static ClassDocumentation<ISGWFormFactor> documentation
    ("The ISGWFormFactor class implements the ISGW model of"
     "Phys. Rev. D39, 799 (1989) for the scalar meson form-factors.",
     "The form factor model of ISGW \\cite{Isgur:1988gb} together with the"
     "form factors for the term which are supressed by the lepton mass from"
     "\\cite{Scora:1989ys,Isgur:1990jf}",
     "\\bibitem{Isgur:1988gb} N.~Isgur, D.~Scora, B.~Grinstein and M.~B.~Wise,\n"
     "Phys.\\ Rev.\\  D {\\bf 39} (1989) 799.\n"
     "%%CITATION = PHRVA,D39,799;%%\n"
     "\\bibitem{Scora:1989ys} D.~Scora and N.~Isgur, \n"
     "Phys.\\ Rev.\\  D {\\bf 40} (1989) 1491.\n"
     "%%CITATION = PHRVA,D40,1491;%%\n"
     "\\bibitem{Isgur:1990jf} N.~Isgur and M.~B.~Wise,\n"
     "Phys.\\ Rev.\\  D {\\bf 43} (1991) 819.\n"
     "%%CITATION = PHRVA,D43,819;%%\n");

  static Parameter<ISGWFormFactor,double> interfaceKappa
    ("Kappa",
     "The relavistic compensation factor of the ISGW model",
     &ISGWFormFactor::_kappa, 0.7, 0.0, 1.0,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceDownMass
    ("DownMass",
     "The mass of the down quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mdown, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceUpMass
    ("UpMass",
     "The mass of the up quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mup, GeV, 0.33*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mstrange, GeV, 0.55*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mcharm, GeV, 1.82*GeV, ZERO, 3.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark in the ISGW model (this is a consituent mass)",
     &ISGWFormFactor::_mbottom, GeV, 5.12*GeV, 3.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSud
    ("BetaSud",
     "The variational parameter for s-wave ud mesons",
     &ISGWFormFactor::_betaSud, GeV, 0.31*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSus
    ("BetaSus",
     "The variational parameter for s-wave us mesons",
     &ISGWFormFactor::_betaSus, GeV, 0.34*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSuc
    ("BetaSuc",
     "The variational parameter for s-wave uc mesons",
     &ISGWFormFactor::_betaSuc, GeV, 0.39*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaSub
    ("BetaSub",
     "The variational parameter for s-wave ub mesons",
     &ISGWFormFactor::_betaSub, GeV, 0.41*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPud
    ("BetaPud",
     "The variational parameter for p-wave ud mesons",
     &ISGWFormFactor::_betaPud, GeV, 0.27*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPus
    ("BetaPus",
     "The variational parameter for p-wave us mesons",
     &ISGWFormFactor::_betaPus, GeV, 0.30*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,Energy> interfaceBetaPuc
    ("BetaPuc",
     "The variational parameter for p-wave uc mesons",
     &ISGWFormFactor::_betaPuc, GeV, 0.34*GeV, ZERO, 1.0*GeV,
     false, false, true);

  static Parameter<ISGWFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &ISGWFormFactor::_thetaeta, -Constants::pi/18., -Constants::pi, Constants::pi,
     false, false, true);

}


// form-factor for scalar to scalar
void ISGWFormFactor::ScalarScalarFormFactor(Energy2 q2, unsigned int iloc,int id0,
					    int id1,Energy mY,Energy mX,
					    Complex & f0,Complex & fp) const {
  Complex d1(0.),d2(0.);
  formFactor(q2,iloc,id0,id1,mY,mX,f0,fp,d1,d2);
}

// form-factor for scalar to vector
void ISGWFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int id0,
					    int id1,Energy mY, Energy mX,
					    Complex & A0,Complex & A1,
					    Complex & A2,Complex & V) const {
  formFactor(q2,iloc,id0,id1,mY,mX,A0,A1,A2,V);
}


// form-factor for scalar to tensor
void ISGWFormFactor::ScalarTensorFormFactor(Energy2 q2, unsigned int iloc, int id0,
					    int id1, Energy mY, Energy mX,
					    complex<InvEnergy2> & h,Complex & k,
					    complex<InvEnergy2> & bp,
					    complex<InvEnergy2> & bm) const {
  Complex f1,f2,f3,f4;
  formFactor(q2,iloc,id0,id1,mY,mX,f1,f2,f3,f4);
  Energy msum(mX+mY);
  h = f1/sqr(msum);
  k = f2;
  bp = f3/sqr(msum);
  bm = f4/sqr(msum);
}

// member which does the work
void ISGWFormFactor::formFactor(Energy2 q2, unsigned int iloc, int, int id1,
				Energy mY,Energy mX, Complex & f1,Complex & f2,
				Complex & f3, Complex & f4) const {
  useMe();
  // work out the flavour of the heavy quarks etc
  int jspin,spect,inquark,outquark;
  formFactorInfo(iloc,jspin,spect,inquark,outquark);
  int ifl0(abs(inquark)),ifl1(abs(outquark)),ifls(abs(spect));
  // determine the multiplet
  int ispin(abs(id1)/1000);
  // masses of the quarks
  Energy mQ(_mquark[ifl0-1]),mq(_mquark[ifl1-1]),ms(_mquark[ifls-1]);
  Energy mtildeX(mq+ms),mtildeY(mQ+ms);
  // wavefunction parameters for the mesons
  Energy betaX(ZERO),betaY(_betaS[ifl0-1][ifls-1]);
  // spin-0 outgoing mesons
  if(ispin==0&&jspin<2) {
    betaX=_betaS[ifl1-1][ifls-1];
  }
  else {
    betaX=_betaP[ifl1-1][ifls-1];
  }
  // compute the F_n function we will need
  Energy2 beta2XY(0.5*(betaX*betaX+betaY*betaY)),tm((mY-mX)*(mY-mX));
  double betar(betaX*betaY/beta2XY),kappa2(_kappa*_kappa),
    slope((tm-q2)/(kappa2*beta2XY));
  Energy mup(mq*mQ/(mQ+mq)),mum(mq*mQ/(mQ-mq));
  double fn(sqrt(mtildeX/mtildeY)*betar*sqrt(betar)*
	    exp(-0.25*ms*ms/(mtildeX*mtildeY)*slope));
  // now we can compute the form-factors
  // for scalars
  if(jspin==0) {
    Complex fp,fm;
    // 1 1S0
    if(ispin==0) {
      double yratio(ms/mtildeX*betaY*betaY/beta2XY);
      fp = fn*(1.+0.5*mQ/mum-0.25*mQ*mq/mup/mum*yratio);
      fm = fn*(1.-(mtildeX+mtildeY)*(0.5/mq-0.25/mup*yratio));
    }
    // 1 3P0
    else if(ispin<100) {
      // extra power of beta factors
      fn*=betar;
      fp = fn*ms*mQ*mq/sqrt(6.)/betaY/mtildeX/mum;
      fm =-fn*ms/betaY/sqrt(6.)*(mtildeX+mtildeY)/mtildeX;
    }
    // 2 1S0
    else {
      throw Exception() << "ISGWFormFactor::formFactor" 
			<< " 2S not implemented" << Exception::abortnow;
    }
    // convert to the standard form
    f1 = Complex(q2/(mY*mY-mX*mX)*fm)+fp;
    f2 = fp;
  }
  // for vectors
  else if(jspin==1) {
    complex<Energy> f;
    complex<InvEnergy> g,ap,am;
    Energy2 betaX2(betaX*betaX),betaY2(betaY*betaY);
    //  1 3S1
    if(ispin==0) {
      f  = 2.*mtildeY*fn;
      g  = 0.5*fn*(1./mq-0.5/mum*ms/mtildeX*betaY2/beta2XY);
      ap =-0.5*fn/mtildeX*(1.+ms/mQ*(betaY2-betaX2)/(betaX2+betaY2)
			   -0.25*ms*ms/mum/mtildeY*betaX2*betaX2/beta2XY/beta2XY);
      am = 0.5*fn/mtildeX*(1.+ms/mQ+ms*ms/mq/mQ*betaX2/beta2XY*
 			   (1.-0.25*(mtildeX+mtildeY)/mtildeY*betaX2/beta2XY));
    }
    // 1 3P1
    else if(ispin==20) {
      fn*=betar;
      f  =-fn*mtildeY*betaY*(1./mum+0.5*ms/mtildeY*beta2XY*slope/betaY2*
			     (1./mq-0.5/mum*ms/mtildeX*betaY2/beta2XY));
      g  = 0.5*fn*ms/mtildeX/betaY;
      ap = 0.25*fn*ms*mQ/mtildeY/betaY/mum*(1.-0.5*ms*mq/mtildeX/mum*betaY2/beta2XY);
      am = -0.25*fn*ms*(mtildeX+mtildeY)/mq/betaY/mtildeY*
 	(1.-0.5*ms*mq/mtildeX/mum*betaY2/beta2XY);
    }
    //  1 1P1
    else if(ispin==10) {
      fn*=betar;
      double ort(1./sqrt(2.));
      f  = fn*ort*mtildeY*betaY/mup;
      g  = 0.25*fn*mtildeY*betaY*ort/mq/mQ/mtildeX;
      Energy msum(mtildeX+mtildeY);
      ap = fn*ms*ort/betaY/mtildeY*(1.+0.5*mQ/mum
				    -0.25*mq*mQ*ms/mum/mup/mtildeX*betaY2/beta2XY);
      am = fn*ms*ort/betaY/mtildeY*(1.-0.5/mq*msum
				    +0.25*ms*betaY2/mup/beta2XY*msum/mtildeX);
    }
    // 2 1S0
    else {
      throw Exception() << "ISGWFormFactor::formFactor" 
			<< " 2S not implemented" << Exception::abortnow;
    }
    // convert to the standard notation
    Energy msum(mX+mY),mdiff(mY-mX);
    Complex ii(0.,1.);
    f2 = -ii*f/msum;
    f3 = +Complex(ii*ap*msum);
    f1 = -Complex(ii*0.5/mX*(am*q2+ii*msum*f2-ii*mdiff*f3));
    f4 =  Complex(ii*g*msum);
  }
  // for tensors
  else if(jspin==2) {
    Energy msum(mX+mY);
    fn *=betar/sqrt(2.);
    double betaXb2(betaX*betaX/beta2XY);
    // 1 3P2
    if(ispin==0) {
      f1 = 0.5*fn*ms/mtildeY/betaY*(1./mq
				    -0.5*ms/mtildeX/mum*betaY*betaY/beta2XY)*sqr(msum);
      f2 = 2.*fn*ms/betaY;
      f3 =-0.5*fn*ms/mtildeX/mQ/betaY*
	(1.-0.5*ms*mQ/mup/mtildeY*betaXb2
	 +0.25*ms*mQ/mtildeY/mum*betaXb2*(1.-0.5*ms*betaXb2/mtildeY))* 
	sqr(msum);
      f4 = 0.5*fn*ms/mtildeX/mQ/betaY*
 	(1.-0.5*ms*mQ/mup/mtildeY*betaXb2+
 	 0.25*ms*betaXb2/mq*(mtildeX+mtildeY)/mtildeY*(1.-0.5*ms*betaXb2/mtildeY))*
	sqr(msum);
    }
  }
  else {
    throw Exception() << "ISGWFormFactor::FormFactor spin = " << jspin 
		      << " but spins higher than 2 not implemented"
		      << Exception::runerror;
  }
  // special for mixing
  double fact(1.);
  if(id1==ParticleID::eta) {
    if(ifl1==3&&ifls==3) fact = -2.*cos(_thetaeta)/sqrt(6.) - sin(_thetaeta)/sqrt(3.);
    else                 fact =     cos(_thetaeta)/sqrt(6.) - sin(_thetaeta)/sqrt(3.);
   
  }
  else if(id1==ParticleID::etaprime) {
    if(ifl1==3&&ifls==3) fact = -2.*sin(_thetaeta)/sqrt(6.) + cos(_thetaeta)/sqrt(3.);
    else                 fact =     sin(_thetaeta)/sqrt(6.) + cos(_thetaeta)/sqrt(3.);
  }
  else if(ifl1==ifls&&ifl1<3) {
    if(abs(ifl1)==1&&int(id1/10)%100==1) fact = -sqrt(0.5);
    else                                 fact =  sqrt(0.5);
  }
  f1*=fact;
  f2*=fact;
  f3*=fact;
  f4*=fact;
}

void ISGWFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ISGWFormFactor " << name() << "\n";
  output << "newdef " << name() << ":Kappa    "    << _kappa        << "\n";
  output << "newdef " << name() << ":DownMass "    << _mdown/GeV    << "\n";
  output << "newdef " << name() << ":UpMass "      << _mup/GeV      << "\n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << "\n";
  output << "newdef " << name() << ":CharmMass "   << _mcharm/GeV   << "\n";
  output << "newdef " << name() << ":BottomMass "  << _mbottom/GeV  << "\n";
  output << "newdef " << name() << ":BetaSud "     << _betaSud/GeV  << "\n";
  output << "newdef " << name() << ":BetaSus "     << _betaSus/GeV  << "\n";
  output << "newdef " << name() << ":BetaSuc "     << _betaSuc/GeV  << "\n";
  output << "newdef " << name() << ":BetaSub "     << _betaSub/GeV  << "\n";
  output << "newdef " << name() << ":BetaPud "     << _betaPud/GeV  << "\n";
  output << "newdef " << name() << ":BetaPus "     << _betaPus/GeV  << "\n";
  output << "newdef " << name() << ":BetaPuc "     << _betaPuc/GeV  << "\n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./LambdabExcitedLambdacSumRuleFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LambdabExcitedLambdacSumRuleFormFactor class.
//

#include "LambdabExcitedLambdacSumRuleFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

LambdabExcitedLambdacSumRuleFormFactor::LambdabExcitedLambdacSumRuleFormFactor() {
  _xi1=0.29;
  _rho2=2.01;
  // modes handled by this form-factor
  // lambda_b to lambda_c1
  addFormFactor(5122,14122,2,2,1,2,5,4);
  // lambda_b to lambda_c1*
  addFormFactor(5122,4124 ,2,4,1,2,5,4);
}

void LambdabExcitedLambdacSumRuleFormFactor::
persistentOutput(PersistentOStream & os) const {
  os << _xi1 << _rho2;
}

void LambdabExcitedLambdacSumRuleFormFactor::
persistentInput(PersistentIStream & is, int) {
  is >> _xi1 >> _rho2;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LambdabExcitedLambdacSumRuleFormFactor,BaryonFormFactor>
describeHerwigLambdabExcitedLambdacSumRuleFormFactor("Herwig::LambdabExcitedLambdacSumRuleFormFactor", "HwFormFactors.so");

void LambdabExcitedLambdacSumRuleFormFactor::Init() {

  static ClassDocumentation<LambdabExcitedLambdacSumRuleFormFactor> documentation
    ("The LambdabExcitedLambdacSumRuleFormFactor class implements the"
     " form-factors for Lambda_b to Lambda_c1(*) from hep-ph/0012114.",
     "Lambda_b to Lambda_c1(*) used the formfactors from \\cite{Huang:2000xw}.",
     "%\\cite{Huang:2000xw}\n"
     "\\bibitem{Huang:2000xw}\n"
     "  M.~Q.~Huang, J.~P.~Lee, C.~Liu and H.~S.~Song,\n"
     "  %``Leading Isgur-Wise form factor of Lambda/b to Lambda/c1 transition  using\n"
     "  %QCD sum rules,''\n"
     "  Phys.\\ Lett.\\  B {\\bf 502}, 133 (2001)\n"
     "  [arXiv:hep-ph/0012114].\n"
     "  %%CITATION = PHLTA,B502,133;%%\n"
     );

  static Parameter<LambdabExcitedLambdacSumRuleFormFactor,double> interfaceXi
    ("Xi",
     "The intercept for the Isgur-Wise form-factor",
     &LambdabExcitedLambdacSumRuleFormFactor::_xi1, 0.29, 0.0, 10.0,
     false, false, true);

  static Parameter<LambdabExcitedLambdacSumRuleFormFactor,double> interfaceRho2
    ("Rho2",
     "The slope parameter for the form-factor.",
     &LambdabExcitedLambdacSumRuleFormFactor::_rho2, 2.01, -10.0, 10.0,
     false, false, true);

}

void LambdabExcitedLambdacSumRuleFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int,int,int,Energy m0,Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo ,
			   Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  double omega(.5/m0/m1*(m0*m0+m1*m1-q2)),orr(1./sqrt(3.));
  // the universal form-factor
  double xi=_xi1*(1.-_rho2*(omega-1.));
  // the couplings in the velocity form
  Complex g1v,g1a,g2v,g2a,g3a,g3v;
  g1v = orr*(omega-1.)*xi;
  g1a = orr*(omega+1.)*xi;
  g2v =-2.*orr*xi;
  g3v = 0.;
  g2a =-2.*orr*xi;
  g3a = 0.;
  // convert to our form
  f1a = g1v-0.5*(m0-m1)*(g2v/m0+g3v/m1);
  f1v =-g1a-0.5*(m0+m1)*(g2a/m0+g3a/m1);
  f2a = Complex(0.5*(m0+m1)*( g2v/m0+g3v/m1));
  f2v =-Complex(0.5*(m0+m1)*( g2a/m0+g3a/m1));
  f3a = Complex(0.5*(m0+m1)*( g2v/m0-g3v/m1));
  f3v = Complex(0.5*(m0+m1)*(-g2a/m0+g3a/m1));
}

void  LambdabExcitedLambdacSumRuleFormFactor::
SpinHalfSpinThreeHalfFormFactor(Energy2 q2,int,int,int,Energy m0,Energy m1,
				Complex & f1v,Complex & f2v,
				Complex & f3v,Complex & f4v,
				Complex & f1a,Complex & f2a,
				Complex & f3a,Complex & f4a,
				FlavourInfo ,
				Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  // the omega value
  double omega(.5/m0/m1*(m0*m0+m1*m1-q2));
  // the universal form-factor
  double xi(_xi1*(1.-_rho2*(omega-1.)));
  // calculate the form factor
  // the couplings in the velocity form
  Complex N1,N2,N3,N4,K1,K2,K3,K4; 
  Energy msum(m0+m1);Energy2 msum2(msum*msum);
  // in the form of the heavy quark papers
  N1 = xi;
  K1 = xi;
  N2 = 0.;
  K2 = 0.;
  N3 = 0.;
  K3 = 0.;
  N4 = 0.;
  K4 = 0.;
  // convert to our form
  f1v =-N4;
  f1a = K4;
  f2v =-N1*msum/m0;
  f2a = K1*msum/m0;
  f3v =-Complex(msum2/m0*(N2/m0+N3/m1));
  f3a = Complex(msum2/m0*(K2/m0+K3/m1));
  f4v =-Complex(msum2/m0/m0*N2);
  f4a = Complex(msum2/m0/m0*K2);
}

void LambdabExcitedLambdacSumRuleFormFactor::dataBaseOutput(ofstream & output,
							    bool header,
							    bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::LambdabExcitedLambdacSumRuleFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":Xi          " << _xi1          << " \n";
  output << "newdef " << name() << ":Rho2        " << _rho2         << " \n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./LightBaryonQuarkModelFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LightBaryonQuarkModelFormFactor class.
//

#include "LightBaryonQuarkModelFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void LightBaryonQuarkModelFormFactor::doinit() {
  BaryonFormFactor::doinit();
  // check that the parameters are consistent
  unsigned int isize=numberOfFactors();
  if(isize!=_f1.size()||isize!=_f2.size()||isize!=_g1.size()||isize!=_g2.size()||
     isize!=_Lambdaf1.size()||isize!=_Lambdaf2.size()||
     isize!=_Lambdag1.size()||isize!=_Lambdag2.size())
    throw InitException() << "Inconsistent parameters in "
			  << "LightBaryonQuarkModelFormFactor::doinit()" 
			  << Exception::abortnow;
}

LightBaryonQuarkModelFormFactor::LightBaryonQuarkModelFormFactor() {
  // the various decay modes handled by the model and the parameters
  // neutron to proton
  addFormFactor(2112,2212,2,2,2,1,1,2);
  _f1.push_back(1.00);_f2.push_back(1.81/GeV);
  _g1.push_back(1.25);_g2.push_back(ZERO);
  _Lambdaf1.push_back(0.69*GeV);_Lambdaf2.push_back(0.96*GeV);
  _Lambdag1.push_back(0.76*GeV);_Lambdag2.push_back(1.04*GeV);
  // sigma+  to lambda
  addFormFactor(3222,3122,2,2,3,2,2,1);
  _f1.push_back(0.00);_f2.push_back(1.04/GeV);
  _g1.push_back(0.60);_g2.push_back(ZERO);
  _Lambdaf1.push_back(-0.32*GeV);_Lambdaf2.push_back(-1.72*GeV);
  _Lambdag1.push_back( 0.77*GeV);_Lambdag2.push_back( 1.05*GeV);
  // sigma-  to lambda
  addFormFactor(3112,3122,2,2,3,1,1,2);
  _f1.push_back(0.00);_f2.push_back(1.04/GeV);
  _g1.push_back(0.60);_g2.push_back(ZERO);
  _Lambdaf1.push_back(-0.32*GeV);_Lambdaf2.push_back(-1.72*GeV);
  _Lambdag1.push_back( 0.77*GeV);_Lambdag2.push_back( 1.05*GeV);
  // sigma-  to sigma0
  addFormFactor(3112,3212,2,2,3,1,1,2);
  _f1.push_back(1.41);_f2.push_back(0.76/GeV);
  _g1.push_back(0.69);_g2.push_back(ZERO);
  _Lambdaf1.push_back(0.60*GeV);_Lambdaf2.push_back(0.81*GeV);
  _Lambdag1.push_back(0.77*GeV);_Lambdag2.push_back(1.04*GeV);
  // sigma0  to sigma+
  addFormFactor(3212,3222,2,2,3,2,1,2);
  _f1.push_back(-1.41);_f2.push_back(-0.76/GeV);
  _g1.push_back(-0.69);_g2.push_back( ZERO);
  _Lambdaf1.push_back(0.60*GeV);_Lambdaf2.push_back(0.81*GeV);
  _Lambdag1.push_back(0.77*GeV);_Lambdag2.push_back(1.04*GeV);
  // Xi- to Xi0
  addFormFactor(3312,3322,2,2,3,3,1,2);
  _f1.push_back(-1.00);_f2.push_back(0.73/GeV);
  _g1.push_back( 0.24);_g2.push_back(ZERO);
  _Lambdaf1.push_back(0.56*GeV);_Lambdaf2.push_back(0.71*GeV);
  _Lambdag1.push_back(0.76*GeV);_Lambdag2.push_back(1.04*GeV);
  // lambda to proton
  addFormFactor(3122,2212,2,2,1,2,3,2);
  _f1.push_back(-1.19);_f2.push_back(-0.850/GeV);
  _g1.push_back(-0.99);_g2.push_back(-0.025/GeV);
  _Lambdaf1.push_back(0.71*GeV);_Lambdaf2.push_back(0.98*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.12*GeV);
  // sigma0 to proton
  addFormFactor(3212,2212,2,2,2,1,3,2);
  _f1.push_back(-0.69);_f2.push_back(0.44/GeV);
  _g1.push_back( 0.19);_g2.push_back(0.0043/GeV);
  _Lambdaf1.push_back(0.64*GeV);_Lambdaf2.push_back(0.84*GeV);
  _Lambdag1.push_back(0.83*GeV);_Lambdag2.push_back(1.16*GeV);
  // sigma- to neutron
  addFormFactor(3112,2112,2,2,1,1,3,2);
  _f1.push_back(-0.97);_f2.push_back(0.62/GeV);
  _g1.push_back(0.27);_g2.push_back(0.0061/GeV);
  _Lambdaf1.push_back(0.64*GeV);_Lambdaf2.push_back(0.90*GeV);
  _Lambdag1.push_back(0.83*GeV);_Lambdag2.push_back(1.16*GeV);
  // xi- to lambda
  addFormFactor(3312,3122,2,2,3,1,3,2);
  _f1.push_back(1.19);_f2.push_back(0.07/GeV);
  _g1.push_back(0.33);_g2.push_back(0.0076/GeV);
  _Lambdaf1.push_back(0.68*GeV);_Lambdaf2.push_back(0.89*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.10*GeV);
  // xi- to sigma0
  addFormFactor(3312,3212,2,2,3,1,3,2);
  _f1.push_back(0.69);_f2.push_back(0.98/GeV);
  _g1.push_back(0.94);_g2.push_back(0.022/GeV);
  _Lambdaf1.push_back(0.75*GeV);_Lambdaf2.push_back(1.05*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.12*GeV);
  // xi0 to sigma+
  addFormFactor(3322,3222,2,2,3,2,3,2);
  _f1.push_back(0.98);_f2.push_back(1.38/GeV);
  _g1.push_back(1.33);_g2.push_back(0.031/GeV);
  _Lambdaf1.push_back(0.75*GeV);_Lambdaf2.push_back(1.05*GeV);
  _Lambdag1.push_back(0.81*GeV);_Lambdag2.push_back(1.12*GeV);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void LightBaryonQuarkModelFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _f1 << ounit(_f2,1/GeV) << _g1 << ounit(_g2,1/GeV) 
     << ounit(_Lambdaf1,GeV) << ounit(_Lambdaf2,GeV) 
     << ounit(_Lambdag1,GeV) << ounit(_Lambdag2,GeV);
}

void LightBaryonQuarkModelFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _f1 >> iunit(_f2,1/GeV) >> _g1 >> iunit(_g2,1/GeV) 
     >> iunit(_Lambdaf1,GeV) >> iunit(_Lambdaf2,GeV) 
     >> iunit(_Lambdag1,GeV) >> iunit(_Lambdag2,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LightBaryonQuarkModelFormFactor,BaryonFormFactor>
describeHerwigLightBaryonQuarkModelFormFactor("Herwig::LightBaryonQuarkModelFormFactor", "HwFormFactors.so");

void LightBaryonQuarkModelFormFactor::Init() {

  static ClassDocumentation<LightBaryonQuarkModelFormFactor> documentation
    ("The LightBaryonQuarkModelFormFactor class implements"
     " the quark model calculation of hep-ph/9409272 for the form-factors"
     " for the light quarks",
     "The quark model calculation of \\cite{Schlumpf:1994fb} was used"
     "for the weak decay of the light baryons",
     "\\bibitem{Schlumpf:1994fb}\n"
     "F.~Schlumpf,\n"
     "Phys.\\ Rev.\\  D {\\bf 51} (1995) 2262 [arXiv:hep-ph/9409272].\n"
     "%%CITATION = PHRVA,D51,2262;%%\n");

  static ParVector<LightBaryonQuarkModelFormFactor,double> interfacef1
    ("f1",
     "The form-factor f1 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_f1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,double> interfaceg1
    ("g1",
     "The form-factor g1 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_g1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,InvEnergy> interfacef2
    ("f2",
     "The form-factor f2 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_f2,
     1./GeV, 0, ZERO, -10./GeV, 10./GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,InvEnergy> interfaceg2
    ("g2",
     "The form-factor g2 at zero q^2",
     &LightBaryonQuarkModelFormFactor::_g2,
     1./GeV, 0, ZERO, -10./GeV, 10./GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdaf1
    ("Lambdaf1",
     "The first mass for the energy dependence of the f1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdaf1,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdaf2
    ("Lambdaf2",
     "The second mass for the energy dependence of the f1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdaf2,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdag1
    ("Lambdag1",
     "The first mass for the energy dependence of the g1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdag1,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<LightBaryonQuarkModelFormFactor,Energy> interfaceLambdag2
    ("Lambdag2",
     "The second mass for the energy dependence of the g1 form-factor.",
     &LightBaryonQuarkModelFormFactor::_Lambdag2,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void LightBaryonQuarkModelFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int mode,int, int, Energy m0, Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo ,
			   Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  // f_3 is zero
  f3v = 0.;
  f3a = 0.;
  // energy dependence of the f1 and g1 factors
  InvEnergy2 lam1,lam2;
  if(_Lambdaf1[mode]>ZERO) {
    lam1 = 1./sqr(_Lambdaf1[mode]);
    lam2 = 1./sqr(_Lambdaf2[mode]);
    f1v= _f1[mode]/(1.-q2*(lam1-q2*sqr(lam2)));
  }
  else {
    f1v = _f1[mode]*(1.+_Lambdaf1[mode]/GeV*q2/GeV2+_Lambdaf2[mode]/GeV*sqr(q2/GeV2));
  }
  lam1 = 1./sqr(_Lambdag1[mode]);
  lam2 = 1./sqr(_Lambdag2[mode]);
  f1a=-_g1[mode]/(1.-q2*(lam1-q2*sqr(lam2)));
  // the f2 and g2 factors
  f2v =-(m0+m1)*_f2[mode];
  f2a = (m0+m1)*_g2[mode]; 
}

void LightBaryonQuarkModelFormFactor::dataBaseOutput(ofstream& output,bool header,
						     bool create) const  {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::LightBaryonQuarkModelFormFactor " 
		    << name() << " \n";
  for(unsigned int ix=0;ix<_f1.size();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":f1 " << ix << "  " << _f1[ix] << "\n";
      output << "newdef " << name() << ":g1 " << ix << "  " << _g1[ix] << "\n";
      output << "newdef " << name() << ":f2 " << ix << "  " << _f2[ix]*GeV << "\n";
      output << "newdef " << name() << ":g2 " << ix << "  " << _g2[ix]*GeV << "\n";
      output << "newdef " << name() << ":Lambdaf1 " << ix << "  " 
	     << _Lambdaf1[ix]/GeV << "\n";
      output << "newdef " << name() << ":Lambdaf2 " << ix << "  " 
	     << _Lambdaf2[ix]/GeV << "\n";
      output << "newdef " << name() << ":Lambdag1 " << ix << "  " 
	     << _Lambdag1[ix]/GeV << "\n";
      output << "newdef " << name() << ":Lambdag2 " << ix << "  " 
	     << _Lambdag2[ix]/GeV << "\n";
    }
    else {
      output << "insert " << name() << ":f1 " << ix << "  " << _f1[ix] << "\n";
      output << "insert " << name() << ":g1 " << ix << "  " << _g1[ix] << "\n";
      output << "insert " << name() << ":f2 " << ix << "  " << _f2[ix]*GeV << "\n";
      output << "insert " << name() << ":g2 " << ix << "  " << _g2[ix]*GeV << "\n";
      output << "insert " << name() << ":Lambdaf1 " << ix << "  " 
	     << _Lambdaf1[ix]/GeV << "\n";
      output << "insert " << name() << ":Lambdaf2 " << ix << "  " 
	     << _Lambdaf2[ix]/GeV << "\n";
      output << "insert " << name() << ":Lambdag1 " << ix << "  " 
	     << _Lambdag1[ix]/GeV << "\n";
      output << "insert " << name() << ":Lambdag2 " << ix << "  " 
	     << _Lambdag2[ix]/GeV << "\n";
    }
  }
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SingletonFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SingletonFormFactor class.
//

#include "SingletonFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

SingletonFormFactor::SingletonFormFactor() {
  // the charm and strange quark masses
  _mcharm   = 1.80*GeV;
  _mstrange = 0.51*GeV;
  // mixing angles
  using Constants::pi;
  _thetalambda = 0.25*pi;
  _thetasigma  = 0.25*pi;
  _thetaxi     = 0.25*pi;
  _thetaxip    = 0.25*pi;
  // the particles handled and the masses for them 
  // lambda_b
  addFormFactor(5122,4122,2,2,1,2,5,4);_polemass.push_back(6.0*GeV);
  // sigma_b
  addFormFactor(5112,4112,2,2,1,1,5,4);_polemass.push_back(6.0*GeV);
  addFormFactor(5212,4212,2,2,2,1,5,4);_polemass.push_back(6.0*GeV);
  addFormFactor(5222,4222,2,2,2,2,5,4);_polemass.push_back(6.0*GeV);
  // omega_b
  addFormFactor(5332,4332,2,2,3,3,5,4);_polemass.push_back(6.4*GeV);
  // xi_b
  addFormFactor(5232,4232,2,2,2,3,5,4);_polemass.push_back(6.0*GeV);
  addFormFactor(5132,4132,2,2,1,3,5,4);_polemass.push_back(6.0*GeV);
  // lambda_c
  addFormFactor(4122,3122,2,2,1,2,4,3);_polemass.push_back(2.5*GeV);
  // sigma_c
  addFormFactor(4112,3112,2,2,1,1,4,3);_polemass.push_back(2.8*GeV);
  addFormFactor(4212,3212,2,2,2,1,4,3);_polemass.push_back(2.8*GeV);
  addFormFactor(4222,3222,2,2,2,2,4,3);_polemass.push_back(2.8*GeV);
  // xi_c
  addFormFactor(4232,3322,2,2,2,3,4,3);_polemass.push_back(2.8*GeV);
  addFormFactor(4132,3312,2,2,1,3,4,3);_polemass.push_back(2.8*GeV);
  // set the inital number of form factors
  initialModes(numberOfFactors());
}

void SingletonFormFactor::doinit() {
  BaryonFormFactor::doinit();
  if(numberOfFactors()!=_polemass.size())
    throw InitException() << "Inconsistent parameters in SingletonFormFactor::doinit()"
			  << Exception::abortnow;
  // calculate the constants for the form-factors
  int id0,id1;
  _xi.clear();
  _nmM.clear();
  _mquark.clear();
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    // id codes for the particles
    particleID(ix,id0,id1);
    if((abs(id0)==5122&&abs(id1)==4122)||(abs(id0)==5132&&abs(id1)==4132)||
       (abs(id0)==5232&&abs(id1)==4232)) {
      _mquark.push_back(_mcharm);
      _xi.push_back(1.);
      _nmM.push_back(1.);
    }
    else if((abs(id0)==5222&&abs(id1)==4222)||(abs(id0)==5212&&abs(id1)==4212)||
	    (abs(id0)==5112&&abs(id1)==4112)||(abs(id0)==5332&&abs(id1)==4332)||
	    (abs(id0)==5312&&abs(id1)==4312)||(abs(id0)==5322&&abs(id1)==4322)) {
      _mquark.push_back(_mcharm);
      _xi.push_back(-1./3.);
      _nmM.push_back(1.);
    }
    else if(abs(id0)==4122&&abs(id1)==3122) {
      _mquark.push_back(_mstrange);
      _xi.push_back(1.);
      _nmM.push_back(sqrt(2./3.)*sin(_thetalambda));
    }
    else if((abs(id0)==4222&&abs(id1)==3222)||(abs(id0)==4212&&abs(id1)==3212)||
	    (abs(id0)==4112&&abs(id1)==3112)) {
      _mquark.push_back(_mstrange);
      _xi.push_back(-1./3.);
      _nmM.push_back(sqrt(2./3.)*cos(_thetasigma));
    }
    else if((abs(id0)==4132&&abs(id1)==3312)||(abs(id0)==4232&&abs(id1)==3322)) {
      _mquark.push_back(_mstrange);
      _xi.push_back(1.);
      _nmM.push_back(1./sqrt(2.)*sin(_thetaxi));
    }
    else if((abs(id0)==4312&&abs(id1)==3322)||(abs(id0)==4322&&abs(id1)==3312)) {
      _mquark.push_back(_mstrange);
      _xi.push_back(-1./3.);
      _nmM.push_back(1./sqrt(6.)*cos(_thetaxip));
    }
    else
      throw InitException() << "Unknown mode in SingletonFormFactor::doinit()"
			    << Exception::abortnow;
  }
}

void SingletonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mcharm,GeV) << ounit(_mstrange,GeV) <<  _thetalambda 
     << _thetasigma << _thetaxi 
     << _thetaxip << ounit(_polemass,GeV) << _xi << _nmM << ounit(_mquark,GeV);
}

void SingletonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mcharm,GeV) >> iunit(_mstrange,GeV) >>  _thetalambda 
     >> _thetasigma >> _thetaxi 
     >> _thetaxip >> iunit(_polemass,GeV) >> _xi >> _nmM >> iunit(_mquark,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SingletonFormFactor,BaryonFormFactor>
describeHerwigSingletonFormFactor("Herwig::SingletonFormFactor", "HwFormFactors.so");

void SingletonFormFactor::Init() {

  static ClassDocumentation<SingletonFormFactor> documentation
    ("The SingletonFormFactor class implements the"
     " form-factors of PRD43, 2939 for the decay of spin-1/2 baryons"
     " containing one heavy quark.",
     "Spin-1/2 baryons with one heavy quark were decayed using "
     "the form factors in \\cite{Singleton:1990ye}.",
     "%\\cite{Singleton:1990ye}\n"
     "\\bibitem{Singleton:1990ye}\n"
     "  R.~L.~Singleton,\n"
     "  %``Semileptonic baryon decays with a heavy quark,''\n"
     "  Phys.\\ Rev.\\  D {\\bf 43} (1991) 2939.\n"
     "  %%CITATION = PHRVA,D43,2939;%%\n"
     );

  static Parameter<SingletonFormFactor,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark",
     &SingletonFormFactor::_mcharm, GeV, 1.8*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<SingletonFormFactor,Energy> interfaceStrangeMass
    ("StrangeMass",
     "The mass of the strange quark",
     &SingletonFormFactor::_mstrange, GeV, 0.51*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaLambda
    ("ThetaLambda",
     "The mixing angle for the Lambda",
     &SingletonFormFactor::_thetalambda, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaSigma
    ("ThetaSigma",
     "The mixing angle for the Sigma",
     &SingletonFormFactor::_thetasigma, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaXi
    ("ThetaXi",
     "The mixing angle for the Xi",
     &SingletonFormFactor::_thetaxi, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static Parameter<SingletonFormFactor,double> interfaceThetaXiPrime
    ("ThetaXiPrime",
     "The mixing angle for the Xi'",
     &SingletonFormFactor::_thetaxip, 0.785398163, 0.0, 6.28318507,
     false, false, true);

  static ParVector<SingletonFormFactor,Energy> interfacePoleMass
    ("PoleMass",
     "The mass for the energy dependence of the form-factors.",
     &SingletonFormFactor::_polemass,
     1.*GeV, 0, ZERO, -10.*GeV, 10.*GeV, false, false, true);
}

// form factor for spin-1/2 to spin-1/2
void SingletonFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc,int, int, Energy m0, Energy m1,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo, Virtuality virt) {
  assert(virt==SpaceLike);
  useMe();
  InvEnergy ratio(0.5/m0);
  // all factors divided by sqrt(4.*m0*m1) normalisation  
  double gbar(m1*_xi[iloc]*_nmM[iloc]/_mquark[iloc]);
  double abar(_xi[iloc]*_nmM[iloc]);
  InvEnergy apbar(-ratio*(m1/_mquark[iloc]-1.)*_xi[iloc]*_nmM[iloc]);
  InvEnergy ambar(apbar);
  InvEnergy gpbar(-ratio*_nmM[iloc]*((_xi[iloc]*m1/_mquark[iloc]-1.)
				  +0.5*(m1-m0)/_mquark[iloc]*(1.-_xi[iloc])));
  InvEnergy gmbar(-ratio*_nmM[iloc]*((_xi[iloc]*m1/_mquark[iloc]-1.)
				  +0.5*(m1+m0)/_mquark[iloc]*(1.-_xi[iloc])));
  // energy dependence
  double ymax(sqr(1.-m1/m0));
  double yres(sqr(_polemass[iloc]/m0));
  double y(q2/sqr(m0)),efact((ymax-yres)/(y-yres));
  f1v = efact*(gbar+(m0+m1)*gpbar);
  f2v = efact*gpbar*(m0+m1);
  f3v = efact*gmbar*(m0+m1);
  f1a = efact*(-abar+(m0-m1)*apbar);
  f2a =-efact*apbar*(m0+m1);
  f3a =-efact*ambar*(m0+m1);
}

void SingletonFormFactor::dataBaseOutput(ofstream & output,bool header,
					 bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::SingletonFormFactor " << name() << " \n";
  output << "newdef " << name() << ":CharmMass " << _mcharm/GeV << " \n";
  output << "newdef " << name() << ":StrangeMass " << _mstrange/GeV << " \n";
  output << "newdef " << name() << ":ThetaLambda " << _thetalambda << " \n";
  output << "newdef " << name() << ":ThetaSigma " << _thetasigma << " \n";
  output << "newdef " << name() << ":ThetaXi " << _thetaxi << " \n";
  output << "newdef " << name() << ":ThetaXiPrime " << _thetaxip << " \n";
  for(unsigned int ix=0;ix<_polemass.size();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":PoleMass " << ix << "  " 
	     << _polemass[ix]/GeV << endl;
    }
    else {
      output << "insert " << name() << ":PoleMass "<< ix << "  " 
	     << _polemass[ix]/GeV << endl;
    }
  }
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./WSBFormFactor.cc"
// -*- C++ -*-
//
// WSBFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WSBFormFactor class.
//

#include "WSBFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Repository/EventGenerator.h"

using namespace Herwig;

void WSBFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_F0.size() ||isize!=_V.size()  ||isize!=_A0.size() ||
     isize!=_A1.size() ||isize!=_A2.size() ||isize!=_mS0.size()||
     isize!=_mS1.size()||isize!=_mV0.size()||isize!=_mV1.size())
    throw InitException() << "Inconsistent parameters in WSBFormFactor::doinit()" 
			  << Exception::abortnow;
}


WSBFormFactor::WSBFormFactor() 
  : _F0(51), _V(51), _A0(51), _A1(51), _A2(51), 
    _mS0(51), _mS1(51), _mV0(51), _mV1(51) {
  // modes handled by this and the parameters model
  // K to pi 
  addFormFactor(-321, 111,0,-2,3,2);
  addFormFactor(-311, 211,0,-1,3,2);
  addFormFactor(-321,-211,0,-2,3,1);
  addFormFactor(-311, 111,0,-1,3,1);
  for(unsigned int ix=0;ix<4;++ix) {
    _F0[ix] = 0.992; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 0.494*GeV; _mV0[ix] = 0.892*GeV; 
    _mS1[ix] = 1.430*GeV; _mV1[ix] = 1.273*GeV; 
  }
  // D to K 
  addFormFactor(421,-321,0,-2,4,3);
  addFormFactor(411,-311,0,-1,4,3);
  for(unsigned int ix=4;ix<6;++ix) {
    _F0[ix] = 0.762; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.97*GeV; _mV0[ix] = 2.11*GeV; 
    _mS1[ix] = 2.60*GeV; _mV1[ix] = 2.53*GeV;
  }
  // D to pi
  addFormFactor(421,-211,0,-2,4,1);
  addFormFactor(421, 111,0,-2,4,2);
  addFormFactor(411, 111,0,-1,4,1);
  addFormFactor(411, 211,0,-1,4,2);
  for(unsigned int ix=6;ix<10;++ix) {
    _F0[ix] = 0.692; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to eta
  addFormFactor(421,221,0,-2,4,2);
  addFormFactor(411,221,0,-1,4,1);
  for(unsigned int ix=10;ix<12;++ix) {
    _F0[ix] = 0.681; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to eta'
  addFormFactor(421,331,0,-2,4,2);
  addFormFactor(411,331,0,-1,4,1);
  for(unsigned int ix=12;ix<14;++ix) {
    _F0[ix] = 0.655; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to K*
  addFormFactor(421,-323,1,-2,4,3);
  addFormFactor(411,-313,1,-1,4,3);
  for(unsigned int ix=14;ix<16;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.226; _A0[ix] = 0.733; 
    _A1[ix] = 0.880; _A2[ix] = 1.147; 
    _mS0[ix] = 1.97*GeV; _mV0[ix] = 2.11*GeV; 
    _mS1[ix] = 2.60*GeV; _mV1[ix] = 2.53*GeV;
  } 
  // D to rho
  addFormFactor(421,-213,1,-2,4,1);
  addFormFactor(421, 113,1,-2,4,2);
  addFormFactor(411, 113,1,-1,4,1);
  addFormFactor(411, 213,1,-1,4,2);
  for(unsigned int ix=16;ix<20;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.225; _A0[ix] = 0.669; 
    _A1[ix] = 0.775; _A2[ix] = 0.923; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D to omega
  addFormFactor(411,223,1,-1,4,1); 
  addFormFactor(421,223,1,-2,4,2); 
  for(unsigned int ix=20;ix<22;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.236; _A0[ix] = 0.669; 
    _A1[ix] = 0.772; _A2[ix] = 0.920; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D_s to eta
  addFormFactor(431,221,0,-3,4,3);
  _F0[22] = 0.723; _V[22]  = 0.000; _A0[22] = 0.000; 
  _A1[22] = 0.000; _A2[22] = 0.000; 
  _mS0[22] = 1.97*GeV; _mV0[22] = 2.11*GeV; 
  _mS1[22] = 2.60*GeV; _mV1[22] = 2.53*GeV; 
  // D_s to eta'
  addFormFactor(431,331,0,-3,4,3);
  _F0[23] = 0.704; _V[23]  = 0.000; _A0[23] = 0.000; 
  _A1[23] = 0.000; _A2[23] = 0.000; 
  _mS0[23] = 1.97*GeV; _mV0[23] = 2.11*GeV; 
  _mS1[23] = 2.60*GeV; _mV1[23] = 2.53*GeV; 
  // D_s to K
  addFormFactor(431,311,0,-3,4,1);
  addFormFactor(431,321,0,-3,4,2);
  for(unsigned int ix=24;ix<26;++ix) {
    _F0[ix] = 0.643; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D_s to K*
  addFormFactor(431,313,1,-3,4,1);
  addFormFactor(431,323,1,-3,4,2);
  for(unsigned int ix=26;ix<28;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 1.250; _A0[ix] = 0.634; 
    _A1[ix] = 0.717; _A2[ix] = 0.853; 
    _mS0[ix] = 1.87*GeV; _mV0[ix] = 2.01*GeV; 
    _mS1[ix] = 2.47*GeV; _mV1[ix] = 2.42*GeV;
  }
  // D_s to phi
  addFormFactor(431,333,1,-3,4,3);
  _F0[28] = 0.000; _V[28]  = 1.319; _A0[28] = 0.700; 
  _A1[28] = 0.820; _A2[28] = 1.076; 
  _mS0[28] = 1.97*GeV; _mV0[28] = 2.11*GeV; 
  _mS1[28] = 2.60*GeV; _mV1[28] = 2.53*GeV; 
  // B to D
  addFormFactor(-521,421,0,-2,5,4);
  addFormFactor(-511,411,0,-2,5,4);
  for(unsigned int ix=29;ix<31;++ix) {
    _F0[ix] = 0.690; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 6.30*GeV; _mV0[ix] = 6.34*GeV; 
    _mS1[ix] = 6.80*GeV; _mV1[ix] = 6.73*GeV;
  }
  // B to K 
  addFormFactor(-521,-321,0,-2,5,3);
  addFormFactor(-511,-311,0,-1,5,3);
  for(unsigned int ix=31;ix<33;++ix) {
    _F0[ix] = 0.379; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.38*GeV; _mV0[ix] = 5.43*GeV; 
    _mS1[ix] = 5.89*GeV; _mV1[ix] = 5.82*GeV;
  }
  // B to pi
  addFormFactor(-521, 111,0,-2,5,2);
  addFormFactor(-511, 211,0,-1,5,2);
  addFormFactor(-521,-211,0,-2,5,1);
  addFormFactor(-511, 111,0,-1,5,1);
  for(unsigned int ix=33;ix<37;++ix) {
    _F0[ix] = 0.333; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to eta
  addFormFactor(-521,221,0,-2,5,2);
  addFormFactor(-511,221,0,-1,5,1);
  for(unsigned int ix=37;ix<39;++ix) {
    _F0[ix] = 0.307; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to eta'
  addFormFactor(-521,331,0,-2,5,2);
  addFormFactor(-511,331,0,-1,5,1);
  for(unsigned int ix=39;ix<41;++ix) {
    _F0[ix] = 0.254; _V[ix]  = 0.000; _A0[ix] = 0.000; 
    _A1[ix] = 0.000; _A2[ix] = 0.000; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to D*
  addFormFactor(-521,423,1,-2,5,4);
  addFormFactor(-511,413,1,-1,5,4);
  for(unsigned int ix=41;ix<43;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 0.705; _A0[ix] = 0.623; 
    _A1[ix] = 0.651; _A2[ix] = 0.686; 
    _mS0[ix] = 6.30*GeV; _mV0[ix] = 6.34*GeV; 
    _mS1[ix] = 6.80*GeV; _mV1[ix] = 6.73*GeV;
  }
  // B to K* 
  addFormFactor(-521,-323,1,-2,5,3);
  addFormFactor(-511,-313,1,-1,5,3);
  for(unsigned int ix=43;ix<45;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 0.369; _A0[ix] = 0.321; 
    _A1[ix] = 0.328; _A2[ix] = 0.331; 
    _mS0[ix] = 5.38*GeV; _mV0[ix] = 5.43*GeV; 
    _mS1[ix] = 5.89*GeV; _mV1[ix] = 5.82*GeV;
  }
  // B to rho
  addFormFactor(-521, 113,1,-2,5,2);
  addFormFactor(-511, 213,1,-1,5,2);
  addFormFactor(-521,-213,1,-2,5,1);
  addFormFactor(-511, 113,1,-1,5,1); 
  for(unsigned int ix=45;ix<49;++ix) {
    _F0[ix] = 0.000; _V[ix]  = 0.329; _A0[ix] = 0.281; 
    _A1[ix] = 0.283; _A2[ix] = 0.283; 
 _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  }
  // B to omega
  addFormFactor(-521,223,1,-2,5,2); 
  addFormFactor(-511,223,1,-1,5,1);
  for(unsigned int ix=49;ix<51;++ix) {
    _F0[ix] = 0.000; _V[ix] = 0.328; _A0[ix] = 0.280; 
    _A1[ix] = 0.281; _A2[ix] = 0.281; 
    _mS0[ix] = 5.27*GeV; _mV0[ix] = 5.32*GeV; 
    _mS1[ix] = 5.78*GeV; _mV1[ix] = 5.71*GeV;
  } 
  // set the initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta=-0.194;
}

void WSBFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _F0 << _V << _A0 << _A1 << _A2 << ounit(_mS0,GeV) 
     << ounit(_mS1,GeV) << ounit(_mV0,GeV) << ounit(_mV1,GeV) << _thetaeta;
}

void WSBFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _F0 >> _V >> _A0 >> _A1 >> _A2 >> iunit(_mS0,GeV) 
     >> iunit(_mS1,GeV) >> iunit(_mV0,GeV) >> iunit(_mV1,GeV) >> _thetaeta;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<WSBFormFactor,ScalarFormFactor>
describeHerwigWSBFormFactor("Herwig::WSBFormFactor", "HwFormFactors.so");

void WSBFormFactor::Init() {

  static ClassDocumentation<WSBFormFactor> documentation
    ("The WSBFormFactor class is the implementation of the form-factors of "
     "Z.Phys.C29,637.",
     "The form factor model of \\cite{Bauer:1986bm,Wirbel:1985ji} was used"
     " for either semi-leptonic or hadronic weak decays",
     "\\bibitem{Bauer:1986bm} M.~Bauer, B.~Stech and M.~Wirbel,\n"
     "Z.\\ Phys.\\  C {\\bf 34} (1987) 103.\n"
     "%%CITATION = ZEPYA,C34,103;%%\n"
     "\\bibitem{Wirbel:1985ji} M.~Wirbel, B.~Stech and M.~Bauer,"
     "Z.\\ Phys.\\  C {\\bf 29} (1985) 637.\n"
     "%%CITATION = ZEPYA,C29,637;%%\n");

  static ParVector<WSBFormFactor,double> interfaceF0
    ("F0",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_F0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceV
    ("V",
     "The form-factor V at zero q^2",
     &WSBFormFactor::_V,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA0
    ("A0",
     "The form-factor A0 at zero q^2",
     &WSBFormFactor::_A0,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA1
    ("A1",
     "The form-factor A1 at zero q^2",
     &WSBFormFactor::_A1,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,double> interfaceA2
    ("A2",
     "The form-factor F0 at zero q^2",
     &WSBFormFactor::_A2,
     0, 0, 0, -10., 10., false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceScalarMass
    ("ScalarMass",
     "The scalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS0,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoScalarMass
    ("PseudoScalarMass",
     "The pseudoscalar mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mS1,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfaceVectorMass
    ("VectorMass",
     "The vector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV0,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static ParVector<WSBFormFactor,Energy> interfacePseudoVectorMass
    ("PseudoVectorMass",
     "The pseudovector mass for the energy dependence of the form-factors.",
     &WSBFormFactor::_mV1,
     1.*GeV, -1, 5.*GeV, -10.*GeV, 10.*GeV, false, false, true);

  static Parameter<WSBFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     &WSBFormFactor::_thetaeta, -0.194, -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void WSBFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
					   int,int id1,
					   Energy, Energy,Complex & f0,
					   Complex & fp) const {
  useMe();
  f0 = _F0[mode]/(1.-q2/sqr(_mS1[mode]));
  fp = _F0[mode]/(1.-q2/sqr(_mV0[mode]));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)) {
    double fact;
    if(id1==ParticleID::eta) {
      if(abs(outquark)==3) fact = -2.*cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
      else                 fact =     cos(_thetaeta)/sqrt(6.)-sin(_thetaeta)/sqrt(3.);
    }
    else if(id1==ParticleID::etaprime) {
      if(abs(outquark)==3) fact = -2.*sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
      else                 fact =     sin(_thetaeta)/sqrt(6.)+cos(_thetaeta)/sqrt(3.);
    }
    else if(id1==ParticleID::pi0&&abs(outquark)==1) fact=-sqrt(0.5);
    else                                            fact= sqrt(0.5);
    f0*=fact;
    fp*=fact;
  }
}

void WSBFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
					   int, int id1, 
					   Energy, Energy,Complex & A0,
					   Complex & A1,Complex & A2,Complex & V) const {
  A0 = -_A0[mode]/(1.-q2/_mS0[mode]/_mS0[mode]);
  A1 = -_A1[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  A2 = -_A2[mode]/(1.-q2/_mV1[mode]/_mV1[mode]);
  V  =   _V[mode]/(1.-q2/_mV0[mode]/_mV0[mode]);
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3) {
    double fact = id1==ParticleID::rho0&&abs(outquark)==1 ? -sqrt(0.5) : sqrt(0.5);
    A0 *= fact;
    A1 *= fact;
    A2 *= fact;
    V  *= fact;
  }
}

void WSBFormFactor::dataBaseOutput(ofstream & output,bool header,bool create) const {
  useMe();
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::WSBFormFactor " << name() << " \n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":F0 " 
	     << ix << "  " << _F0[ix] << "\n";
      output << "newdef " << name() << ":V  " 
	     << ix << "  " << _V[ix]  << "\n";
      output << "newdef " << name() << ":A0 " 
	     << ix << "  " << _A0[ix] << "\n";
      output << "newdef " << name() << ":A1 " 
	     << ix << "  " << _A1[ix] << "\n";
      output << "newdef " << name() << ":A2 " 
	     << ix << "  " << _A2[ix] << "\n";
      output << "newdef " << name() << ":ScalarMass " 
	     << ix << "  " << _mS0[ix]/GeV << "\n";
      output << "newdef " << name() << ":PseudoScalarMass " 
	     << ix << "  " << _mS1[ix]/GeV << "\n";
      output << "newdef " << name() << ":VectorMass " 
	     << ix << "  " << _mV0[ix]/GeV << "\n";
      output << "newdef " << name() << ":PseudoVectorMass " 
	     << ix << "  " << _mV1[ix]/GeV << "\n";
    }
    else {
      output << "insert " << name() << ":F0 " 
	     << ix << "  " << _F0[ix] << "\n";
      output << "insert " << name() << ":V  " 
	     << ix << "  " << _V[ix]  << "\n";
      output << "insert " << name() << ":A0 " 
	     << ix << "  " << _A0[ix] << "\n";
      output << "insert " << name() << ":A1 " 
	     << ix << "  " << _A1[ix] << "\n";
      output << "insert " << name() << ":A2 " 
	     << ix << "  " << _A2[ix] << "\n";
      output << "insert " << name() << ":ScalarMass " 
	     << ix << "  " << _mS0[ix]/GeV << "\n";
      output << "insert " << name() << ":PseudoScalarMass " 
	     << ix << "  " << _mS1[ix]/GeV << "\n";
      output << "insert " << name() << ":VectorMass " 
	     << ix << "  " << _mV0[ix]/GeV << "\n";
      output << "insert " << name() << ":PseudoVectorMass " 
	     << ix << "  " << _mV1[ix]/GeV << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./KiselevBcFormFactor.cc"
// -*- C++ -*-
//
// KiselevBcFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KiselevBcFormFactor class.
//
#include "KiselevBcFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KiselevBcFormFactor::KiselevBcFormFactor() :
  _fp  (16,0.    ), _fm  (16,0.    ), _FV( 16,ZERO),
  _F0A (16,ZERO), _FpA (16,ZERO), _FmA (16,ZERO),
  _Mfp (16,ZERO), _Mfm (16,ZERO), _MFV (16,ZERO),
  _MF0A(16,ZERO), _MFpA(16,ZERO), _MFmA(16,ZERO) {
  // B_c to B_s
  addFormFactor(541,531,0,-5,4,3);
  _fp[0] =   1.30    ;_Mfp[0] =   1.8*GeV;
  _fm[0] =  -5.80    ;_Mfm[0] =   1.8*GeV;
  // B_c to B_s*
  addFormFactor(541,533,1,-5,4,3);
  _FV[1] =   1.10/GeV;_MFV[1] =   1.8*GeV;
  _F0A[1] =  8.10*GeV;_MF0A[1] =  1.8*GeV;
  _FpA[1] =  0.20/GeV;_MFpA[1] =  1.8*GeV;
  _FmA[1] =  1.80/GeV;_MFmA[1] =  1.8*GeV;
  // B_c to B
  addFormFactor(541,511,0,-5,4,1);
  addFormFactor(541,521,0,-5,4,2);
  for(unsigned int ix=2;ix<4;++ix) {
    _fp[ix] =   1.27    ;_Mfp[ix] =   1.7*GeV;
    _fm[ix] =  -7.30    ;_Mfm[ix] =   1.7*GeV;
  }
  // B_c to B*
  addFormFactor(541,513,1,-5,4,1);
  addFormFactor(541,523,1,-5,4,2);
  for(unsigned int ix=4;ix<6;++ix) {
    _FV[ix] =   1.35/GeV; _MFV[ix] =   2.2*GeV;
    _F0A[ix] =  9.80*GeV; _MF0A[ix] =  3.2*GeV;
    _FpA[ix] =  0.35/GeV; _MFpA[ix] =  2.2*GeV;
    _FmA[ix] =  2.50/GeV; _MFmA[ix] =  3.2*GeV;
  }
  // B_c to D
  addFormFactor(-541,-411,0,-4,5,1);
  addFormFactor(-541,-421,0,-4,5,2);
  for(unsigned int ix=6;ix<8;++ix) {
    _fp[ix] =   0.32    ; _Mfp[ix] =   5.0*GeV;
    _fm[ix] =  -0.34    ; _Mfm[ix] =   5.0*GeV;
  }
  // B_c to D*
  addFormFactor(-541,-413,1,-4,5,1);
  addFormFactor(-541,-423,1,-4,5,2);
  for(unsigned int ix=8;ix<10;++ix) {
    _FV[ix] =   0.20 /GeV; _MFV[ix] =   6.2*GeV;
    _F0A[ix] =  3.60 *GeV; _MF0A[ix] = -1.0*GeV;
    _FpA[ix] = -0.062/GeV; _MFpA[ix] =  6.2*GeV;
    _FmA[ix] =  0.10 /GeV; _MFmA[ix] =  6.2*GeV;
  }
  // B_c to D_s
  addFormFactor(-541,-431,0,-4,5,3);
  _fp[10] =   0.45    ; _Mfp[10] =   5.0*GeV;
  _fm[10] =  -0.43    ; _Mfm[10] =   5.0*GeV;
  // B_c to D_s
  addFormFactor(-541,-433,1,-4,5,3);
  _FV[11] =   0.24 /GeV; _MFV[11] =   6.2*GeV;
  _F0A[11] =  4.70 *GeV; _MF0A[11] = -1.0*GeV;
  _FpA[11] = -0.077/GeV; _MFpA[11] =  6.2*GeV;
  _FmA[11] =  0.13 /GeV; _MFmA[11] =  6.2*GeV;
  // B_c to eta_c
  addFormFactor(-541,441,0,-4,5,4);
  _fp[12] =   0.66    ; _Mfp[12] =   4.5*GeV;
  _fm[12] =  -0.36    ; _Mfm[12] =   4.5*GeV;
  // B_c to J/psi
  addFormFactor(-541,443,1,-4,5,4);
  _FV[13] =   0.11 /GeV; _MFV[13] =   5.5*GeV;
  _F0A[13] =  5.90 *GeV; _MF0A[13] =  5.5*GeV;
  _FpA[13] = -0.074/GeV; _MFpA[13] =  5.5*GeV;
  _FmA[13] =  0.12 /GeV; _MFmA[13] =  5.5*GeV;
  // B_c to eta_c(2S)
  addFormFactor(-541,100441,0,-4,5,4);
  _fp[14] =   0.17    ; _Mfp[14] =   4.5*GeV;
  _fm[14] =  -0.16    ; _Mfm[14] =   4.5*GeV;
  // B_c to J/psi(2S)
  addFormFactor(-541,100443,1,-4,5,4);
  _FV[15] =   0.035/GeV; _MFV[15] =   4.2*GeV;
  _F0A[15] =  1.686*GeV; _MF0A[15] =  4.2*GeV;
  _FpA[15] = -0.015/GeV; _MFpA[15] =  4.2*GeV;
  _FmA[15] =  0.052/GeV; _MFmA[15] =  4.2*GeV;
  // set the inital number of modes
  initialModes(numberOfFactors());
}

void KiselevBcFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_fp.size() ||isize!=_fm.size()  ||isize!=_FV.size()  ||isize!=_F0A.size()||
     isize!=_FpA.size()||isize!=_FmA.size() ||isize!=_Mfp.size() ||isize!=_Mfm.size()||
     isize!=_MFV.size()||isize!=_MF0A.size()||isize!=_MFpA.size()||isize!=_MFmA.size())
    throw InitException() << "Inconsistent parameters in KiselevBcFormFactor::doinit()" 
			  << Exception::abortnow;

}

void KiselevBcFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _fp << _fm << ounit(_FV,1/GeV) << ounit(_F0A,GeV) << ounit(_FpA,1/GeV) 
     << ounit(_FmA,1/GeV) << ounit(_Mfp,GeV) << ounit(_Mfm,GeV) << ounit(_MFV,GeV) 
     << ounit(_MF0A,GeV) << ounit(_MFpA,GeV) << ounit(_MFmA,GeV);
}

void KiselevBcFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _fp >> _fm >> iunit(_FV,1/GeV) >> iunit(_F0A,GeV) >> iunit(_FpA,1/GeV) 
     >> iunit(_FmA,1/GeV) >> iunit(_Mfp,GeV) >> iunit(_Mfm,GeV) >> iunit(_MFV,GeV) 
     >> iunit(_MF0A,GeV) >> iunit(_MFpA,GeV) >> iunit(_MFmA,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KiselevBcFormFactor,ScalarFormFactor>
describeHerwigKiselevBcFormFactor("Herwig::KiselevBcFormFactor", "HwFormFactors.so");

void KiselevBcFormFactor::Init() {

  static ClassDocumentation<KiselevBcFormFactor> documentation
    ("The KiselevBcFormFactor class implements the form factors from hep-ph/0211021"
     " for the decay of the B_c",
     "The form factors of \\cite{Kiselev:2002vz} for the decay of the $B_c$ meson"
     " were used.",
     "\\bibitem{Kiselev:2002vz} V.~V.~Kiselev, arXiv:hep-ph/0211021.\n"
     "%%CITATION = HEP-PH/0211021;%%");

  static ParVector<KiselevBcFormFactor,double> interfaceFplus
    ("Fplus",
     "The value of the f_+ form factor at q^2=0",
     &KiselevBcFormFactor::_fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<KiselevBcFormFactor,double> interfaceFminus
    ("Fminus",
     "The value of the f_- form factor at q^2=0",
     &KiselevBcFormFactor::_fm, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFV
    ("FV",
     "The value of the F_V form factor at q^2=0",
     &KiselevBcFormFactor::_FV, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceF0A
    ("F0A",
     "The value of the F_0^A form factor at q^2=0",
     &KiselevBcFormFactor::_F0A, 1.*GeV, -1, ZERO, -10.*GeV, 10.*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFpA
    ("FplusA",
     "The value of the F_+^A form factor at q^2=0",
     &KiselevBcFormFactor::_FpA, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,InvEnergy> interfaceFmA
    ("FminusA",
     "The value of the F_-^A form factor at q^2=0",
     &KiselevBcFormFactor::_FmA, 1./GeV, -1, ZERO, -10./GeV, 10./GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFplus
    ("MpoleFplus",
     "The pole mass for the f_+ form factor",
     &KiselevBcFormFactor::_Mfp, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFminus
    ("MpoleFminus",
     "The pole mass for the f_- form factor",
     &KiselevBcFormFactor::_Mfm, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFV
    ("MpoleFV",
     "The pole mass for the f_V form factor",
     &KiselevBcFormFactor::_MFV, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleF0A
    ("MpoleF0A",
     "The pole mass for the f_0^A form factor",
     &KiselevBcFormFactor::_MF0A, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFpA
    ("MpoleFplusA",
     "The pole mass for the f_+^A form factor",
     &KiselevBcFormFactor::_MFpA, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KiselevBcFormFactor,Energy> interfaceMpoleFmA
    ("MpoleFminusA",
     "The pole mass for the f_-^A form factor",
     &KiselevBcFormFactor::_MFmA, GeV, -1, ZERO, -2.0*GeV, 10.0*GeV,
     false, false, true);
}

void KiselevBcFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int iloc,
						 int,int,Energy m0,Energy m1,
						 Complex & f0,Complex & fp) const {
  useMe();
  fp = _fp[iloc]/(1.-q2/_Mfp[iloc]/_Mfp[iloc]);
  f0 = Complex(q2/(m0+m1)/(m0-m1)*_fm[iloc]/(1.-q2/_Mfm[iloc]/_Mfm[iloc]))+fp;
}

void KiselevBcFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int iloc, int,
						 int,Energy m0, Energy m1,
						 Complex & A0,Complex & A1,Complex & A2,
						 Complex & V) const {
  useMe();
  InvEnergy fv,fp,fm;
  Energy f0;
  if(_MFV[iloc]>ZERO)  fv = _FV[iloc]/(1.-q2/_MFV[iloc]/_MFV[iloc]);
  else                  fv = _FV[iloc];
  if(_MFmA[iloc]>ZERO) fm = _FmA[iloc]/(1.-q2/_MFmA[iloc]/_MFmA[iloc]);
  else                  fm = _FmA[iloc];
  if(_MFpA[iloc]>ZERO) fp = _FpA[iloc]/(1.-q2/_MFpA[iloc]/_MFpA[iloc]);
  else                  fp = _FpA[iloc];
  if(_MF0A[iloc]>ZERO) f0 = _F0A[iloc]/(1.-q2/_MF0A[iloc]/_MF0A[iloc]);
  else                  f0 = _F0A[iloc];
  Energy msum(m0+m1);
  V  = fv*msum;
  A1 =-f0/msum;
  A2 = fp*msum;
  A0 =-0.5/m1*(f0+msum*(m0-m1)*fp+q2*fm);
}

void KiselevBcFormFactor::dataBaseOutput(ofstream & output,
					 bool header,bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KiselevBcFormFactor " << name() << " \n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":Fplus "  << ix << "  " 
	     << _fp[ix]  << "\n";
      output << "newdef " << name() << ":Fminus "  << ix << "  " 
	     << _fm[ix]  << "\n";
      output << "newdef " << name() << ":FV "  << ix << "  " 
	     << _FV[ix]*GeV  << "\n";
      output << "newdef " << name() << ":F0A "  << ix << "  " 
	     << _F0A[ix]/GeV  << "\n";
      output << "newdef " << name() << ":FplusA "  << ix << "  " 
	     << _FpA[ix]*GeV  << "\n";
      output << "newdef " << name() << ":FminusA "  << ix << "  " 
	     << _FmA[ix]*GeV  << "\n";
      output << "newdef " << name() << ":MpoleFplus "  << ix << "  " 
	     << _Mfp[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFminus "  << ix << "  " 
	     << _Mfm[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFV "  << ix << "  " 
	     << _MFV[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleF0A "  << ix << "  " 
	     << _MF0A[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFplusA "  << ix << "  " 
	     << _MFpA[ix]/GeV  << "\n";
      output << "newdef " << name() << ":MpoleFminusA "  << ix << "  " 
	     << _MFmA[ix]/GeV  << "\n";
    }
    else {
      output << "insert " << name() << ":Fplus "  << ix << "  " 
	     << _fp[ix]  << "\n";
      output << "insert " << name() << ":Fminus "  << ix << "  " 
	     << _fm[ix]  << "\n";
      output << "insert " << name() << ":FV "  << ix << "  " 
	     << _FV[ix]*GeV  << "\n";
      output << "insert " << name() << ":F0A "  << ix << "  " 
	     << _F0A[ix]/GeV  << "\n";
      output << "insert " << name() << ":FplusA "  << ix << "  " 
	     << _FpA[ix]*GeV  << "\n";
      output << "insert " << name() << ":FminusA "  << ix << "  " 
	     << _FmA[ix]*GeV  << "\n";
      output << "insert " << name() << ":MpoleFplus "  << ix << "  " 
	     << _Mfp[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFminus "  << ix << "  " 
	     << _Mfm[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFV "  << ix << "  " 
	     << _MFV[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleF0A "  << ix << "  " 
	     << _MF0A[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFplusA "  << ix << "  " 
	     << _MFpA[ix]/GeV  << "\n";
      output << "insert " << name() << ":MpoleFminusA "  << ix << "  " 
	     << _MFmA[ix]/GeV  << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./MelikhovFormFactor.cc"
// -*- C++ -*-
//
// MelikhovFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MelikhovFormFactor class.
//

#include "MelikhovFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

MelikhovFormFactor::MelikhovFormFactor() : 
  _ifit(1), _Rplus0(0), _Mplus(ZERO), _nplus(0), _RV0(0), _MV(ZERO), _nV(0),
  _R10(0), _M1(ZERO), _n1(0), _R20(0), _M2(ZERO), _n2(0) {
  // the possible modes
  // B to rho
  addFormFactor(-521, 113,1,2,5,2);
  addFormFactor(-511, 213,1,1,5,2);
  addFormFactor(-521,-213,1,2,5,1);
  addFormFactor(-511, 113,1,1,5,1);
  // B to pi
  addFormFactor(-521, 111,0,2,5,2);
  addFormFactor(-511, 211,0,1,5,2);
  addFormFactor(-521,-211,0,2,5,1);
  addFormFactor(-511, 111,0,1,5,1);
  // set the initial number of modes
  initialModes(numberOfFactors());
}

void MelikhovFormFactor::doinit() {
  ScalarFormFactor::doinit();
  // the parameters for the different fits
  double Rplus0[4]={0.29    ,0.20    ,0.21    ,0.26    };
  Energy Mplus[4] ={6.29*GeV,6.22*GeV,5.90*GeV,5.44*GeV};
  double nplus[4] ={2.35    ,2.45    ,2.33    ,1.72    };
  double RV[4]    ={0.30    ,0.20    ,0.21    ,0.29    };
  Energy MV[4]    ={6.28*GeV,6.22*GeV,5.90*GeV,5.46*GeV};
  double nv[4]    ={2.36    ,2.46    ,2.35    ,1.73    };
  double R10[4]   ={0.27    ,0.20    ,0.21    ,0.29    };
  Energy M1[4]    ={7.07*GeV,6.78*GeV,6.50*GeV,5.68*GeV};
  double n1[4]    ={2.65    ,2.65    ,2.70    ,1.67    };
  double R20[4]   ={0.25    ,0.19    ,0.20    ,0.28    };
  Energy M2[4]    ={6.13*GeV,6.00*GeV,5.90*GeV,5.36*GeV};
  double n2[4]    ={2.17    ,2.34    ,2.45    ,1.67    };
  // set the values
  _Rplus0=Rplus0[_ifit-1];
  _Mplus=Mplus[_ifit-1];
  _nplus=nplus[_ifit-1];
  _RV0=RV[_ifit-1];
  _MV=MV[_ifit-1];
  _nV=nv[_ifit-1];
  _R10=R10[_ifit-1];
  _M1=M1[_ifit-1];
  _n1=n1[_ifit-1];
  _R20=R20[_ifit-1];
  _M2=M2[_ifit-1];
  _n2=n2[_ifit-1];
}

void MelikhovFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _ifit << _Rplus0 << ounit(_Mplus,GeV) << _nplus << _RV0 << ounit(_MV,GeV) 
     << _nV << _R10 << ounit(_M1,GeV) << _n1 << _R20 << ounit(_M2,GeV) << _n2; 
}

void MelikhovFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _ifit >> _Rplus0 >> iunit(_Mplus,GeV) >> _nplus >> _RV0 >> iunit(_MV,GeV) 
     >> _nV >> _R10 >> iunit(_M1,GeV) >> _n1 >> _R20 >> iunit(_M2,GeV) >> _n2; 
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MelikhovFormFactor,ScalarFormFactor>
describeHerwigMelikhovFormFactor("Herwig::MelikhovFormFactor", "HwFormFactors.so");

void MelikhovFormFactor::Init() {

  static ClassDocumentation<MelikhovFormFactor> documentation
    ("The MelikhovFormFactor class implements the form factors from"
     " Phys. Lett. B 380 (1996) 363 for B to pi,rho",
     "The form factors of \\cite{Melikhov:1996ge} were used for $B\\to \\pi,\\rho$.",
     "\\bibitem{Melikhov:1996ge} D.~Melikhov, Phys.\\ Lett.\\  B {\\bf 380} (1996) 363\n"
     "[arXiv:hep-ph/9603340]. %%CITATION = PHLTA,B380,363;%%\n");

  static Parameter<MelikhovFormFactor,unsigned int> interfaceFit
    ("Fit",
     "Which of the fits from hep-ph/9603340 to use",
     &MelikhovFormFactor::_ifit, 1, 1, 4,
     false, false, true);

}

// form-factor for scalar to scalar
void MelikhovFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
						int,int id1,
						Energy, Energy,Complex & f0,
						Complex & fp) const {
  useMe();
  fp = _Rplus0/pow((1.-q2/_Mplus/_Mplus),_nplus);
  f0 = fp;
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)) {
    double fact = (id1==ParticleID::pi0&&abs(outquark)==1) ? -sqrt(0.5) : sqrt(0.5);
    f0 *= fact;
    fp *= fact;
  }
}

void MelikhovFormFactor::ScalarVectorFormFactor(Energy2 q2,unsigned int mode,
						int, int id1, Energy mY, Energy mX,
						Complex & A0,Complex & A1,Complex & A2,
						Complex & V) const {
  useMe();
  // constants
  double r(mX/mY),y(q2/mY/mY);
  // V form-factor
  V  =-double((1.+r)*_RV0/pow((1.-q2/_MV/_MV),_nV));
  // A_1 form-factor
  A1 = (1.+r*r-y)/(1.+r)*_R10/pow((1.-q2/_M1/_M1),_n1);
  // A_2 form-factor
  A2 = (1.+r)*(1.-r*r-y)/((1.+r)*(1.+r)-y)*_R20/pow((1.-q2/_M2/_M2),_n2);
  // set the a_- factor to zero
//   A0 = 0.5/mX*((mY+mX)*A1-(mY-mX)*A2);
  A0 = 0.;
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)==abs(spect)&&abs(spect)<3) {
    double fact = (id1==ParticleID::rho0&&abs(outquark)==1) ? -sqrt(0.5) : sqrt(0.5);
    A0 *= fact;
    A1 *= fact;
    A2 *= fact;
    V  *= fact;
  }
}

void MelikhovFormFactor::dataBaseOutput(ofstream & output,bool header,
					bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::MelikhovFormFactor " << name() << " \n";
  output << "newdef " << name() << ":Fit " << _ifit << " \n";
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./MelikhovStechFormFactor.cc"
// -*- C++ -*-
//
// MelikhovStechFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MelikhovStechFormFactor class.
//

#include "MelikhovStechFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Herwig;

MelikhovStechFormFactor::MelikhovStechFormFactor() 
: _fplus0(42,0.),_sigma1fp(42,0.),_sigma2fp(42,0.),
  _f00   (42,0.),_sigma1f0(42,0.),_sigma2f0(42,0.),
  _fT0   (42,0.),_sigma1fT(42,0.),_sigma2fT(42,0.),
  _V0    (42,0.),_sigma1V0(42,0.),_sigma2V0(42,0.),
  _A00   (42,0.),_sigma1A0(42,0.),_sigma2A0(42,0.),
  _A10   (42,0.),_sigma1A1(42,0.),_sigma2A1(42,0.),
  _A20   (42,0.),_sigma1A2(42,0.),_sigma2A2(42,0.),
  _T10   (42,0.),_sigma1T1(42,0.),_sigma2T1(42,0.),
  _T20   (42,0.),_sigma1T2(42,0.),_sigma2T2(42,0.),
  _T30   (42,0.),_sigma1T3(42,0.),_sigma2T3(42,0.),
  _massP(42,ZERO), _massV(42,ZERO) {
  // form factors for D to K
  addFormFactor(421,-321,0,-2,4,3);
  addFormFactor(411,-311,0,-1,4,3);
  for(unsigned int ix=0;ix<2;++ix) {
    _fplus0[ix] = 0.78    ; _sigma1fp[ix] = 0.24    ;_sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.78    ; _sigma1f0[ix] = 0.38    ;_sigma2f0[ix] = 0.46;
    _fT0[ix]    = 0.75    ; _sigma1fT[ix] = 0.27    ;_sigma2fT[ix] = 0.00;
    _massP[ix]  = 1.97*GeV; _massV[ix]    = 2.11*GeV;
  }
  // form factors for D to K*
  addFormFactor(421,-323,1,-2,4,3);
  addFormFactor(411,-313,1,-1,4,3);
  for(unsigned int ix=2;ix<4;++ix) {
    _V0[ix]    = 1.03    ; _sigma1V0[ix] = 0.27    ; _sigma2V0[ix] = 0.00;
    _A00[ix]   = 0.76    ; _sigma1A0[ix] = 0.17    ; _sigma2A0[ix] = 0.00;
    _A10[ix]   = 0.66    ; _sigma1A1[ix] = 0.30    ; _sigma2A1[ix] = 0.20;
    _A20[ix]   = 0.49    ; _sigma1A2[ix] = 0.67    ; _sigma2A2[ix] = 0.16;
    _T10[ix]   = 0.78    ; _sigma1T1[ix] = 0.25    ; _sigma2T1[ix] = 0.00;
    _T20[ix]   = 0.78    ; _sigma1T2[ix] = 0.02    ; _sigma2T2[ix] = 1.80;
    _T30[ix]   = 0.45    ; _sigma1T3[ix] = 1.23    ; _sigma2T3[ix] = 0.34;
    _massP[ix] = 1.97*GeV; _massV[ix]    = 2.11*GeV;
  }
  // form factors for D to pi
  addFormFactor(421,-211,0,-2,4,1);
  addFormFactor(421, 111,0,-2,4,2);
  addFormFactor(411, 111,0,-1,4,1);
  addFormFactor(411, 211,0,-1,4,2);
  for(unsigned int ix=4;ix<8;++ix) {
    _fplus0[ix] = 0.69    ; _sigma1fp[ix] = 0.30; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.69    ; _sigma1f0[ix] = 0.54; _sigma2f0[ix] = 0.32;
    _fT0[ix]    = 0.60    ; _sigma1fT[ix] = 0.34; _sigma2fT[ix] = 0.00;
    _massP[ix]  = 1.87*GeV; _massV[ix]    = 2.01*GeV;
  }
  // form factors for D to rho
  addFormFactor(421,-213,1,-2,4,1);
  addFormFactor(421, 113,1,-2,4,2);
  addFormFactor(411, 113,1,-1,4,1);
  addFormFactor(411, 213,1,-1,4,2);
  for(unsigned int ix=8;ix<12;++ix) {
    _V0[ix]  = 0.90; _sigma1V0[ix] = 0.46; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.66; _sigma1A0[ix] = 0.36; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.59; _sigma1A1[ix] = 0.50; _sigma2A1[ix] = 0.00;
    _A20[ix] = 0.49; _sigma1A2[ix] = 0.89; _sigma2A2[ix] = 0.00;
    _T10[ix] = 0.66; _sigma1T1[ix] = 0.44; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.66; _sigma1T2[ix] = 0.38; _sigma2T2[ix] = 0.50;
    _T30[ix] = 0.31; _sigma1T3[ix] = 1.10; _sigma2T3[ix] = 0.17;
    _massP[ix] = 1.87*GeV; _massV[ix] = 2.01*GeV;
  }
  // B to D
  addFormFactor(-521,421,0,-2,5,4);
  addFormFactor(-511,411,0,-2,5,4);
  for(unsigned int ix=12;ix<14;++ix) {
    _fplus0[ix] = 0.67; _sigma1fp[ix] = 0.57; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.67; _sigma1f0[ix] = 0.78; _sigma2f0[ix] = 0.00;
    _fT0[ix]    = 0.69; _sigma1fT[ix] = 0.56; _sigma2fT[ix] = 0.00;
    _massP[ix]  = 6.4*GeV; _massV[ix] = 6.4*GeV;
  }
  // B to D*
  addFormFactor(-521,423,1,-2,5,4);
  addFormFactor(-511,413,1,-1,5,4);
  for(unsigned int ix=14;ix<16;++ix) {
    _V0[ix]  = 0.76; _sigma1V0[ix] = 0.57; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.69; _sigma1A0[ix] = 0.58; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.66; _sigma1A1[ix] = 0.78; _sigma2A1[ix] = 0.00;
    _A20[ix] = 0.62; _sigma1A2[ix] = 1.40; _sigma2A2[ix] = 0.41;
    _T10[ix] = 0.68; _sigma1T1[ix] = 0.57; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.68; _sigma1T2[ix] = 0.64; _sigma2T2[ix] = 0.00;
    _T30[ix] = 0.33; _sigma1T3[ix] = 1.46; _sigma2T3[ix] = 0.50;
    _massP[ix] = 6.4*GeV; _massV[ix] = 6.4*GeV;
  }
  // B to K
  addFormFactor(-521,-321,0,-2,5,3);
  addFormFactor(-511,-311,0,-1,5,3);
  for(unsigned int ix=16;ix<18;++ix) {
    _fplus0[ix] = 0.36; _sigma1fp[ix] = 0.43; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.36; _sigma1f0[ix] = 0.70; _sigma2f0[ix] = 0.27;
    _fT0[ix]    = 0.35; _sigma1fT[ix] = 0.43; _sigma2fT[ix] = 0.00;
    _massP[ix] = 5.37*GeV; _massV[ix] = 5.42*GeV;
  }
  // B to K*
  addFormFactor(-521,-323,1,-2,5,3);
  addFormFactor(-511,-313,1,-1,5,3);
  for(unsigned int ix=18;ix<20;++ix) {
    _V0[ix]  = 0.44; _sigma1V0[ix] = 0.45; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.45; _sigma1A0[ix] = 0.46; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.36; _sigma1A1[ix] = 0.64; _sigma2A1[ix] = 0.36;
    _A20[ix] = 0.32; _sigma1A2[ix] = 1.23; _sigma2A2[ix] = 0.38;
    _T10[ix] = 0.39; _sigma1T1[ix] = 0.45; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.39; _sigma1T2[ix] = 0.72; _sigma2T2[ix] = 0.62;
    _T30[ix] = 0.27; _sigma1T3[ix] = 1.31; _sigma2T3[ix] = 0.41;
    _massP[ix] = 5.37*GeV; _massV[ix] = 5.42*GeV;
  }
  // B to pi
  addFormFactor(-521, 111,0,-2,5,2);
  addFormFactor(-511, 211,0,-1,5,2);
  addFormFactor(-521,-211,0,-2,5,1);
  addFormFactor(-511, 111,0,-1,5,1);
  for(unsigned int ix=20;ix<24;++ix) {
    _fplus0[ix] = 0.29; _sigma1fp[ix] = 0.48; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.29; _sigma1f0[ix] = 0.76; _sigma2f0[ix] = 0.28;
    _fT0[ix]    = 0.28; _sigma1fT[ix] = 0.48; _sigma2fT[ix] = 0.00;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // B to rho
  addFormFactor(-521, 113,1,-2,5,2);
  addFormFactor(-511, 213,1,-1,5,2);
  addFormFactor(-521,-213,1,-2,5,1);
  addFormFactor(-511, 113,1,-1,5,1);
  for(unsigned int ix=24;ix<28;++ix) {
    _V0[ix]  = 0.31; _sigma1V0[ix] = 0.59; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.30; _sigma1A0[ix] = 0.54; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.26; _sigma1A1[ix] = 0.73; _sigma2A1[ix] = 0.10;
    _A20[ix] = 0.24; _sigma1A2[ix] = 1.40; _sigma2A2[ix] = 0.50;
    _T10[ix] = 0.27; _sigma1T1[ix] = 0.60; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.27; _sigma1T2[ix] = 0.74; _sigma2T2[ix] = 0.19;
    _T30[ix] = 0.19; _sigma1T3[ix] = 1.42; _sigma2T3[ix] = 0.51;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // D_s to K
  addFormFactor(431,311,0,-3,4,1);
  addFormFactor(431,321,0,-3,4,2);
  for(unsigned int ix=28;ix<30;++ix) {
    _fplus0[ix] = 0.72; _sigma1fp[ix] = 0.20; _sigma2fp[ix] = 0.00;
    _f00[ix]    = 0.72; _sigma1f0[ix] = 0.41; _sigma2f0[ix] = 0.70;
    _fT0[ix]    = 0.77; _sigma1fT[ix] = 0.24; _sigma2fT[ix] = 0.00;
    _massP[ix] = 1.87*GeV; _massV[ix] = 2.01*GeV;
  }
  // D_s to K*
  addFormFactor(431,313,1,-3,4,1);
  addFormFactor(431,323,1,-3,4,2);
  for(unsigned int ix=30;ix<32;++ix) {
    _V0[ix]  = 1.04; _sigma1V0[ix] =  0.24; _sigma2V0[ix] = 0.00;
    _A00[ix] = 0.67; _sigma1A0[ix] =  0.20; _sigma2A0[ix] = 0.00;
    _A10[ix] = 0.57; _sigma1A1[ix] =  0.29; _sigma2A1[ix] = 0.42;
    _A20[ix] = 0.42; _sigma1A2[ix] =  0.58; _sigma2A2[ix] = 0.00;
    _T10[ix] = 0.71; _sigma1T1[ix] =  0.22; _sigma2T1[ix] = 0.00;
    _T20[ix] = 0.71; _sigma1T2[ix] = -0.06; _sigma2T2[ix] = 0.44;
    _T30[ix] = 0.45; _sigma1T3[ix] =  1.08; _sigma2T3[ix] = 0.68;
    _massP[ix] = 1.87*GeV; _massV[ix] = 2.01*GeV;
  }
  // D_s to eta
  addFormFactor(431,221,0,-3,4,3);
  _fplus0[32] = 0.78; _sigma1fp[32] = 0.23; _sigma2fp[32] = 0.00;
  _f00[32]    = 0.78; _sigma1f0[32] = 0.33; _sigma2f0[32] = 0.38;
  _fT0[32]    = 0.80; _sigma1fT[32] = 0.24; _sigma2fT[32] = 0.00;
  _massP[32] = 1.97*GeV; _massV[32] = 2.11*GeV;
  // D_s to eta'
  addFormFactor(431,331,0,-3,4,3);
  _fplus0[33] = 0.78; _sigma1fp[33] = 0.23; _sigma2fp[33] = 0.00;
  _f00[33]    = 0.78; _sigma1f0[33] = 0.21; _sigma2f0[33] = 0.76;
  _fT0[33]    = 0.94; _sigma1fT[33] = 0.24; _sigma2fT[33] = 0.00;
  _massP[33] = 1.97*GeV; _massV[33] = 2.11*GeV;
  // D_s to phi
  addFormFactor(431,333,1,-3,4,3);
  _V0[34]  = 1.10; _sigma1V0[34] = 0.26; _sigma2V0[34] = 0.00;
  _A00[34] = 0.73; _sigma1A0[34] = 0.10; _sigma2A0[34] = 0.00;
  _A10[34] = 0.64; _sigma1A1[34] = 0.29; _sigma2A1[34] = 0.00;
  _A20[34] = 0.47; _sigma1A2[34] = 0.63; _sigma2A2[34] = 0.00;
  _T10[34] = 0.77; _sigma1T1[34] = 0.25; _sigma2T1[34] = 0.00;
  _T20[34] = 0.77; _sigma1T2[34] = 0.02; _sigma2T2[34] = 2.01;
  _T30[34] = 0.46; _sigma1T3[34] = 1.34; _sigma2T3[34] = 0.45;
  _massP[34] = 1.97*GeV; _massV[34] = 2.11*GeV;
  // B_s to K
  addFormFactor(-531, 311  ,0,-3,5,1);
  addFormFactor(-531, 321  ,0,-3,5,2);
  for(unsigned int ix=35;ix<37;++ix) {
    _fplus0[ix] = 0.31; _sigma1fp[ix] = 0.63; _sigma2fp[ix] = 0.33;
    _f00[ix]    = 0.31; _sigma1f0[ix] = 0.93; _sigma2f0[ix] = 0.70;
    _fT0[ix]    = 0.31; _sigma1fT[ix] = 0.61; _sigma2fT[ix] = 0.30;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // B_s to K*
  addFormFactor(-531, 313  ,1,-3,5,1);
  addFormFactor(-531, 323  ,1,-3,5,2);
  for(unsigned int ix=37;ix<39;++ix) {
    _V0[ix]  = 0.38; _sigma1V0[ix] = 0.66; _sigma2V0[ix] = 0.30;
    _A00[ix] = 0.37; _sigma1A0[ix] = 0.60; _sigma2A0[ix] = 0.16;
    _A10[ix] = 0.29; _sigma1A1[ix] = 0.86; _sigma2A1[ix] = 0.60;
    _A20[ix] = 0.26; _sigma1A2[ix] = 1.32; _sigma2A2[ix] = 0.54;
    _T10[ix] = 0.32; _sigma1T1[ix] = 0.66; _sigma2T1[ix] = 0.31;
    _T20[ix] = 0.32; _sigma1T2[ix] = 0.98; _sigma2T2[ix] = 0.90;
    _T30[ix] = 0.23; _sigma1T3[ix] = 1.42; _sigma2T3[ix] = 0.62;
    _massP[ix] = 5.27*GeV; _massV[ix] = 5.32*GeV;
  }
  // B_s to eta
  addFormFactor(-531, 221  ,0,-3,5,3);
  _fplus0[39] = 0.36; _sigma1fp[39] = 0.60; _sigma2fp[39] = 0.20;
  _f00[39]    = 0.36; _sigma1f0[39] = 0.80; _sigma2f0[39] = 0.40;
  _fT0[39]    = 0.36; _sigma1fT[39] = 0.58; _sigma2fT[39] = 0.18;
  _massP[39] = 5.37*GeV; _massV[39] = 5.42*GeV;
  // B_s to eta'
  addFormFactor(-531, 331  ,0,-3,5,3);
  _fplus0[40] = 0.36; _sigma1fp[40] = 0.60; _sigma2fp[40] = 0.20;
  _f00[40]    = 0.36; _sigma1f0[40] = 0.80; _sigma2f0[40] = 0.45;
  _fT0[40]    = 0.39; _sigma1fT[40] = 0.58; _sigma2fT[40] = 0.18;
  _massP[40] = 5.37*GeV; _massV[40] = 5.42*GeV;
  // B_s to phi
  addFormFactor(-531, 333  ,1,-3,5,3);
  _V0[41]  = 0.44; _sigma1V0[41] = 0.62; _sigma2V0[41] = 0.20;
  _A00[41] = 0.42; _sigma1A0[41] = 0.55; _sigma2A0[41] = 0.12;
  _A10[41] = 0.34; _sigma1A1[41] = 0.73; _sigma2A1[41] = 0.42;
  _A20[41] = 0.31; _sigma1A2[41] = 1.30; _sigma2A2[41] = 0.52;
  _T10[41] = 0.38; _sigma1T1[41] = 0.62; _sigma2T1[41] = 0.20;
  _T20[41] = 0.38; _sigma1T2[41] = 0.83; _sigma2T2[41] = 0.71;
  _T30[41] = 0.26; _sigma1T3[41] = 1.41; _sigma2T3[41] = 0.57;
  _massP[41] = 5.37*GeV; _massV[41] = 5.42*GeV;
  // set the initial number of modes
  initialModes(numberOfFactors());
  // eta-eta' mixing angle
  _thetaeta = 2./9.*Constants::pi;
}

void MelikhovStechFormFactor::doinit() {
  ScalarFormFactor::doinit();
  unsigned int isize(numberOfFactors());
  if(isize!=_fplus0.size()||isize!=_sigma1fp.size()||isize!=_sigma2fp.size()||
     isize!=_f00.size()   ||isize!=_sigma1f0.size()||isize!=_sigma2f0.size()||
     isize!=_fT0.size()   ||isize!=_sigma1fT.size()||isize!=_sigma2fT.size()||
     isize!=_V0.size()    ||isize!=_sigma1V0.size()||isize!=_sigma2V0.size()||
     isize!=_A00.size()   ||isize!=_sigma1A0.size()||isize!=_sigma2A0.size()||
     isize!=_A10.size()   ||isize!=_sigma1A1.size()||isize!=_sigma2A1.size()||
     isize!=_A20.size()   ||isize!=_sigma1A2.size()||isize!=_sigma2A2.size()||
     isize!=_T10.size()   ||isize!=_sigma1T1.size()||isize!=_sigma2T1.size()||
     isize!=_T20.size()   ||isize!=_sigma1T2.size()||isize!=_sigma2T2.size()||
     isize!=_T30.size()   ||isize!=_sigma1T3.size()||isize!=_sigma2T3.size()||
     isize!=_massP.size() ||isize!=_massV.size())
    throw InitException() << "Inconsistent parameters in " 
			  << "MelikhovStechFormFactor::doinit()" 
			  << Exception::abortnow;
}

void MelikhovStechFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _fplus0 << _sigma1fp << _sigma2fp << _f00 << _sigma1f0 << _sigma2f0 << _fT0 
     << _sigma1fT << _sigma2fT << _V0 << _sigma1V0 << _sigma2V0 << _A00 << _sigma1A0 
     << _sigma2A0 << _A10 << _sigma1A1 << _sigma2A1 << _A20 << _sigma1A2 << _sigma2A2 
     << _T10 << _sigma1T1 << _sigma2T1 << _T20 << _sigma1T2 << _sigma2T2 << _T30 
     << _sigma1T3 << _sigma2T3 << ounit(_massP,GeV) << ounit(_massV,GeV) << _thetaeta;
}

void MelikhovStechFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _fplus0 >> _sigma1fp >> _sigma2fp >> _f00 >> _sigma1f0 >> _sigma2f0 >> _fT0 
     >> _sigma1fT >> _sigma2fT >> _V0 >> _sigma1V0 >> _sigma2V0 >> _A00 >> _sigma1A0 
     >> _sigma2A0 >> _A10 >> _sigma1A1 >> _sigma2A1 >> _A20 >> _sigma1A2 >> _sigma2A2 
     >> _T10 >> _sigma1T1 >> _sigma2T1 >> _T20 >> _sigma1T2 >> _sigma2T2 >> _T30 
     >> _sigma1T3 >> _sigma2T3 >> iunit(_massP,GeV) >> iunit(_massV,GeV) >> _thetaeta;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<MelikhovStechFormFactor,ScalarFormFactor>
describeHerwigMelikhovStechFormFactor("Herwig::MelikhovStechFormFactor", "HwFormFactors.so");

void MelikhovStechFormFactor::Init() {

  static ClassDocumentation<MelikhovStechFormFactor> documentation
    ("The MelikhovStechFormFactor class is the implementation of the"
     " form factors from Phys. Rev. D62  014006 (2000).",
     "The form factors of \\cite{Melikhov:2000yu} were used.",
     "\\bibitem{Melikhov:2000yu} D.~Melikhov and B.~Stech,\n"
     "Phys.\\ Rev.\\  D {\\bf 62} (2000) 014006 [arXiv:hep-ph/0001113].\n"
     "%%CITATION = PHRVA,D62,014006;%%\n");

  static ParVector<MelikhovStechFormFactor,double> interfaceFPlus0
    ("FPlus0",
     "The value of the F_+ form factor at q^2=0",
     &MelikhovStechFormFactor::_fplus0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFpsigma_1
    ("F+sigma_1",
     "The value of sigma_1 for the F_+ form factor",
     &MelikhovStechFormFactor::_sigma1fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFpsigma_2
    ("F+sigma_2",
     "The value of sigma_2 for the F_+ form factor",
     &MelikhovStechFormFactor::_sigma2fp, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF00
    ("F00",
     "The value of the F_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_f00, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF0sigma_1
    ("F0sigma_1",
     "The value of sigma_1 for the F_0 form factor",
     &MelikhovStechFormFactor::_sigma1f0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceF0sigma_2
    ("F0sigma_2",
     "The value of sigma_2 for the F_0 form factor",
     &MelikhovStechFormFactor::_sigma2f0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFT0
    ("FT0",
     "The value of the F_T form factor at q^2=0",
     &MelikhovStechFormFactor::_fT0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFTsigma_1
    ("FTsigma_1",
     "The value of sigma_1 for the F_T form factor",
     &MelikhovStechFormFactor::_sigma1fT, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceFTsigma_2
    ("FTsigma_2",
     "The value of sigma_2 for the F_T form factor",
     &MelikhovStechFormFactor::_sigma2fT, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV00
    ("V00",
     "The value of the V_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV0sigma_1
    ("V0sigma_1",
     "The value of sigma_1 for the V_0 form factor",
     &MelikhovStechFormFactor::_sigma1V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceV0sigma_2
    ("V0sigma_2",
     "The value of sigma_2 for the V_0 form factor",
     &MelikhovStechFormFactor::_sigma2V0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA00
    ("A00",
     "The value of the A_0 form factor at q^2=0",
     &MelikhovStechFormFactor::_A00, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA0sigma_1
    ("A0sigma_1",
     "The value of sigma_1 for the A_0 form factor",
     &MelikhovStechFormFactor::_sigma1A0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA0sigma_2
    ("A0sigma_2",
     "The value of sigma_2 for the A_0 form factor",
     &MelikhovStechFormFactor::_sigma2A0, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA10
    ("A10",
     "The value of the A_1 form factor at q^2=0",
     &MelikhovStechFormFactor::_A10, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA1sigma_1
    ("A1sigma_1",
     "The value of sigma_1 for the A_1 form factor",
     &MelikhovStechFormFactor::_sigma1A1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA1sigma_2
    ("A1sigma_2",
     "The value of sigma_2 for the A_1 form factor",
     &MelikhovStechFormFactor::_sigma2A1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA20
    ("A20",
     "The value of the A_2 form factor at q^2=0",
     &MelikhovStechFormFactor::_A20, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA2sigma_1
    ("A2sigma_1",
     "The value of sigma_1 for the A_2 form factor",
     &MelikhovStechFormFactor::_sigma1A2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceA2sigma_2
    ("A2sigma_2",
     "The value of sigma_2 for the A_2 form factor",
     &MelikhovStechFormFactor::_sigma2A2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT10
    ("T10",
     "The value of the T_1 form factor at q^2=0",
     &MelikhovStechFormFactor::_T10, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT1sigma_1
    ("T1sigma_1",
     "The value of sigma_1 for the T_1 form factor",
     &MelikhovStechFormFactor::_sigma1T1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT1sigma_2
    ("T1sigma_2",
     "The value of sigma_2 for the T_1 form factor",
     &MelikhovStechFormFactor::_sigma2T1, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT20
    ("T20",
     "The value of the T_2 form factor at q^2=0",
     &MelikhovStechFormFactor::_T20, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT2sigma_1
    ("T2sigma_1",
     "The value of sigma_1 for the T_2 form factor",
     &MelikhovStechFormFactor::_sigma1T2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT2sigma_2
    ("T2sigma_2",
     "The value of sigma_2 for the T_2 form factor",
     &MelikhovStechFormFactor::_sigma2T2, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT30
    ("T30",
     "The value of the T_3 form factor at q^2=0",
     &MelikhovStechFormFactor::_T30, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT3sigma_1
    ("T3sigma_1",
     "The value of sigma_1 for the T_3 form factor",
     &MelikhovStechFormFactor::_sigma1T3, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,double> interfaceT3sigma_2
    ("T3sigma_2",
     "The value of sigma_2 for the T_3 form factor",
     &MelikhovStechFormFactor::_sigma2T3, -1, 0.0, -10.0, 10.0,
     false, false, true);

  static ParVector<MelikhovStechFormFactor,Energy> interfaceMassP
    ("MassP",
     "The mass of the pseudoscalar for the q^2 dependence of the form factors.",
     &MelikhovStechFormFactor::_massP, GeV, -1, ZERO, ZERO, 10.*GeV,
     false, false, false);

  static ParVector<MelikhovStechFormFactor,Energy> interfaceMassV
    ("MassV",
     "The mass of the vector for the q^2 dependence of the form factors.",
     &MelikhovStechFormFactor::_massV, GeV, -1, ZERO, ZERO, 10.*GeV,
     false, false, false);

  static Parameter< MelikhovStechFormFactor,double> interfaceThetaEtaEtaPrime
    ("ThetaEtaEtaPrime",
     "The eta-eta' mixing angle",
     & MelikhovStechFormFactor::_thetaeta, 2.*Constants::pi/9.,
     -Constants::pi, Constants::pi,
     false, false, true);
}

// form-factor for scalar to scalar
void MelikhovStechFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int mode,
						     int,int id1,
						     Energy, Energy,Complex & f0,
						     Complex & fp) const {
  useMe();
  double ratio(q2/sqr(_massV[mode]));
  fp = _fplus0[mode]/(1.-ratio)/(1.-ratio*(_sigma1fp[mode]-_sigma2fp[mode]*ratio));
  f0 = _f00[mode]              /(1.-ratio*(_sigma1f0[mode]-_sigma2f0[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>3) return;
  double fact;
  if(id1==ParticleID::eta) {
    if(abs(outquark)==3) fact =-sin(_thetaeta);
    else                 fact = cos(_thetaeta)*sqrt(0.5);
  }
  else if(id1==ParticleID::etaprime) {
    if(abs(outquark)==3)  fact=cos(_thetaeta);
    else                  fact=sin(_thetaeta);
  }
  else if(id1==ParticleID::pi0&&abs(outquark)==1) fact=-sqrt(0.5);
  else                                            fact= sqrt(0.5);
  f0 *= fact;
  fp *= fact;
}
  
void MelikhovStechFormFactor::ScalarScalarSigmaFormFactor(Energy2 q2,unsigned int mode,
							  int,int id1,
							  Energy, Energy,
							  Complex & fT) const {
  useMe();
  double ratio(q2/sqr(_massV[mode]));
  fT = _fT0[mode]/(1.-ratio)/(1.-ratio*(_sigma1fT[mode]-_sigma2fT[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>3) return;
  double fact;
  if(id1==ParticleID::eta) {
    if(abs(outquark)==3) fact = -sin(_thetaeta);
    else                 fact =  cos(_thetaeta)*sqrt(0.5);
  }
  else if(id1==ParticleID::etaprime) {
    if(abs(outquark)==3) fact = cos(_thetaeta);
    else                 fact = sin(_thetaeta);
  }
  else if(id1==ParticleID::pi0&&abs(outquark)==1) fact=-sqrt(0.5);
  else                                            fact= sqrt(0.5);
  fT*=fact;
}

void MelikhovStechFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int mode,
						     int, int id1,
						     Energy, Energy,Complex & A0,
						     Complex & A1,Complex & A2,
						     Complex & V) const {
  useMe();
  double ratio(q2/sqr(_massV[mode])),ratioP(q2/sqr(_massP[mode]));
  A0= _A00[mode]/(1.-ratioP)/(1.-ratioP*(_sigma1A0[mode]-_sigma2A0[mode]*ratioP));
  A1= _A10[mode]            /(1.-ratio *(_sigma1A1[mode]-_sigma2A1[mode]*ratio));
  A2= _A20[mode]            /(1.-ratio *(_sigma1A2[mode]-_sigma2A2[mode]*ratio));
  V =-_V0[mode] /(1.-ratio )/(1.-ratio *(_sigma1V0[mode]-_sigma2V0[mode]*ratio ));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>2) return;
  double fact(sqrt(0.5));
  if(id1==ParticleID::rho0&&abs(outquark)==1) fact=-fact;
  A0 *= fact;
  A1 *= fact;
  A2 *= fact;
  V  *= fact;
}

void MelikhovStechFormFactor::ScalarVectorSigmaFormFactor(Energy2 q2,unsigned int mode,
							  int,int id1,
							  Energy, Energy,
							  Complex & T1,Complex & T2,
							  Complex & T3) const {
  useMe();
  double ratio(q2/sqr(_massV[mode]));
  T1= _T10[mode]/(1.-ratio)/(1.-ratio*(_sigma1T1[mode]-_sigma2T1[mode]*ratio));
  T2= _T20[mode]           /(1.-ratio*(_sigma1T2[mode]-_sigma2T2[mode]*ratio));
  T3= _T30[mode]           /(1.-ratio*(_sigma1T3[mode]-_sigma2T3[mode]*ratio));
  int jspin,spect,inquark,outquark;
  formFactorInfo(mode,jspin,spect,inquark,outquark);
  if(abs(outquark)!=abs(spect)) return;
  if(abs(spect)>2) return;
  double fact(sqrt(0.5));
  if(id1==ParticleID::rho0&&abs(outquark)==1){fact=-fact;}
  T1 *= fact;
  T2 *= fact;
  T3 *= fact;
}

void MelikhovStechFormFactor::dataBaseOutput(ofstream & output,bool header,
					     bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::MelikhovStechFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":ThetaEtaEtaPrime " << _thetaeta  << "\n";
  for(unsigned int ix=0;ix<numberOfFactors();++ix) {
    if(ix<initialModes()) {
      output << "newdef " << name() << ":FPlus0 " 
	     << ix << "  " << _fplus0[ix] << "\n";
      output << "newdef " << name() << ":F+sigma_1 " 
	     << ix << "  " << _sigma1fp[ix] << "\n";
      output << "newdef " << name() << ":F+sigma_2 " 
	     << ix << "  " << _sigma2fp[ix] << "\n";
      output << "newdef " << name() << ":F00 " 
	     << ix << "  " << _f00[ix] << "\n";
      output << "newdef " << name() << ":F0sigma_1 " 
	     << ix << "  " << _sigma1f0[ix] << "\n";
      output << "newdef " << name() << ":F0sigma_2 " 
	     << ix << "  " << _sigma2f0[ix] << "\n";
      output << "newdef " << name() << ":FT0 " 
	     << ix << "  " << _fT0[ix] << "\n";
      output << "newdef " << name() << ":FTsigma_1 " 
	     << ix << "  " << _sigma1fT[ix] << "\n";
      output << "newdef " << name() << ":FTsigma_2 " 
	     << ix << "  " << _sigma2fT[ix] << "\n";
      output << "newdef " << name() << ":V00 " 
	     << ix << "  " << _V0[ix] << "\n";
      output << "newdef " << name() << ":V0sigma_1 " 
	     << ix << "  " << _sigma1V0[ix] << "\n";
      output << "newdef " << name() << ":V0sigma_2 " 
	     << ix << "  " << _sigma2V0[ix] << "\n";
      output << "newdef " << name() << ":A00 " 
	     << ix << "  " << _A00[ix] << "\n";
      output << "newdef " << name() << ":A0sigma_1 " 
	     << ix << "  " << _sigma1A0[ix] << "\n";
      output << "newdef " << name() << ":A0sigma_2 " 
	     << ix << "  " << _sigma2A0[ix] << "\n";
      output << "newdef " << name() << ":A10 " 
	     << ix << "  " << _A10[ix] << "\n";
      output << "newdef " << name() << ":A1sigma_1 " 
	     << ix << "  " << _sigma1A1[ix] << "\n";
      output << "newdef " << name() << ":A1sigma_2 " 
	     << ix << "  " << _sigma2A1[ix] << "\n";
      output << "newdef " << name() << ":A20 " 
	     << ix << "  " << _A20[ix] << "\n";
      output << "newdef " << name() << ":A2sigma_1 " 
	     << ix << "  " << _sigma1A2[ix] << "\n";
      output << "newdef " << name() << ":A2sigma_2 " 
	     << ix << "  " << _sigma2A2[ix] << "\n";
      output << "newdef " << name() << ":T10 " 
	     << ix << "  " << _T10[ix] << "\n";
      output << "newdef " << name() << ":T1sigma_1 " 
	     << ix << "  " << _sigma1T1[ix] << "\n";
      output << "newdef " << name() << ":T1sigma_2 " 
	     << ix << "  " << _sigma2T1[ix] << "\n";
      output << "newdef " << name() << ":T20 " 
	     << ix << "  " << _T20[ix] << "\n";
      output << "newdef " << name() << ":T2sigma_1 " 
	     << ix << "  " << _sigma1T2[ix] << "\n";
      output << "newdef " << name() << ":T2sigma_2 " 
	     << ix << "  " << _sigma2T2[ix] << "\n";
      output << "newdef " << name() << ":T30 " 
	     << ix << "  " << _T30[ix] << "\n";
      output << "newdef " << name() << ":T3sigma_1 " 
	     << ix << "  " << _sigma1T3[ix] << "\n";
      output << "newdef " << name() << ":T3sigma_2 " 
	     << ix << "  " << _sigma2T3[ix] << "\n";
      output << "newdef " << name() << ":MassP " 
	     << ix << "  " << _massP[ix]/GeV << "\n";
      output << "newdef " << name() << ":MassV " 
	     << ix << "  " << _massV[ix]/GeV << "\n";
    }
    else {
      output << "insert " << name() << ":FPlus0 " 
	     << ix << "  " << _fplus0[ix] << "\n";
      output << "insert " << name() << ":F+sigma_1 " 
	     << ix << "  " << _sigma1fp[ix] << "\n";
      output << "insert " << name() << ":F+sigma_2 " 
	     << ix << "  " << _sigma2fp[ix] << "\n";
      output << "insert " << name() << ":F00 " 
	     << ix << "  " << _f00[ix] << "\n";
      output << "insert " << name() << ":F0sigma_1 " 
	     << ix << "  " << _sigma1f0[ix] << "\n";
      output << "insert " << name() << ":F0sigma_2 " 
	     << ix << "  " << _sigma2f0[ix] << "\n";
      output << "insert " << name() << ":FT0 " 
	     << ix << "  " << _fT0[ix] << "\n";
      output << "insert " << name() << ":FTsigma_1 " 
	     << ix << "  " << _sigma1fT[ix] << "\n";
      output << "insert " << name() << ":FTsigma_2 " 
	     << ix << "  " << _sigma2fT[ix] << "\n";
      output << "insert " << name() << ":V00 " 
	     << ix << "  " << _V0[ix] << "\n";
      output << "insert " << name() << ":V0sigma_1 " 
	     << ix << "  " << _sigma1V0[ix] << "\n";
      output << "insert " << name() << ":V0sigma_2 " 
	     << ix << "  " << _sigma2V0[ix] << "\n";
      output << "insert " << name() << ":A00 " 
	     << ix << "  " << _A00[ix] << "\n";
      output << "insert " << name() << ":A0sigma_1 " 
	     << ix << "  " << _sigma1A0[ix] << "\n";
      output << "insert " << name() << ":A0sigma_2 " 
	     << ix << "  " << _sigma2A0[ix] << "\n";
      output << "insert " << name() << ":A10 " 
	     << ix << "  " << _A10[ix] << "\n";
      output << "insert " << name() << ":A1sigma_1 " 
	     << ix << "  " << _sigma1A1[ix] << "\n";
      output << "insert " << name() << ":A1sigma_2 " 
	     << ix << "  " << _sigma2A1[ix] << "\n";
      output << "insert " << name() << ":A20 " 
	     << ix << "  " << _A20[ix] << "\n";
      output << "insert " << name() << ":A2sigma_1 " 
	     << ix << "  " << _sigma1A2[ix] << "\n";
      output << "insert " << name() << ":A2sigma_2 " 
	     << ix << "  " << _sigma2A2[ix] << "\n";
      output << "insert " << name() << ":T10 " 
	     << ix << "  " << _T10[ix] << "\n";
      output << "insert " << name() << ":T1sigma_1 " 
	     << ix << "  " << _sigma1T1[ix] << "\n";
      output << "insert " << name() << ":T1sigma_2 " 
	     << ix << "  " << _sigma2T1[ix] << "\n";
      output << "insert " << name() << ":T20 " 
	     << ix << "  " << _T20[ix] << "\n";
      output << "insert " << name() << ":T2sigma_1 " 
	     << ix << "  " << _sigma1T2[ix] << "\n";
      output << "insert " << name() << ":T2sigma_2 " 
	     << ix << "  " << _sigma2T2[ix] << "\n";
      output << "insert " << name() << ":T30 " 
	     << ix << "  " << _T30[ix] << "\n";
      output << "insert " << name() << ":T3sigma_1 " 
	     << ix << "  " << _sigma1T3[ix] << "\n";
      output << "insert " << name() << ":T3sigma_2 " 
	     << ix << "  " << _sigma2T3[ix] << "\n";
      output << "insert " << name() << ":MassP " 
	     << ix << "  " << _massP[ix]/GeV << "\n";
      output << "insert " << name() << ":MassV " 
	     << ix << "  " << _massV[ix]/GeV << "\n";
    }
  }
  ScalarFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./BtoSGammaFlatEnergy.cc"
// -*- C++ -*-
//
// BtoSGammaFlatEnergy.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaFlatEnergy class.
//

#include "BtoSGammaFlatEnergy.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<BtoSGammaFlatEnergy,BtoSGammaHadronicMass>
describeHerwigBtoSGammaFlatEnergy("Herwig::BtoSGammaFlatEnergy", "HwFormFactors.so");

void BtoSGammaFlatEnergy::Init() {

  static ClassDocumentation<BtoSGammaFlatEnergy> documentation
    ("The BtoSGammaFlatEnergy class implements a model of the hadronic mass"
     " in B to s gamma decays designed to give a flat energy spectrum for"
     " testing purposes.");

}

Energy BtoSGammaFlatEnergy::hadronicMass(Energy mb,Energy mquark) {
  Energy upper(min(mb,maxMass())),lower(max(minMass(),mquark));
  double rand(UseRandom::rnd());
  return sqrt(upper*upper*rand+(1.-rand)*lower*lower);
}

void BtoSGammaFlatEnergy::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BtoSGammaFlatEnergy " 
		    << name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() 
		    << "\";" << endl;
}
#line 1 "./BtoSGammaKagan.cc"
// -*- C++ -*-
//
// BtoSGammaKagan.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BtoSGammaKagan class.
//

#include "BtoSGammaKagan.h"
#include "Herwig/Utilities/Maths.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/GaussianIntegrator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;
using Herwig::Math::Li2;

DescribeClass<BtoSGammaKagan,BtoSGammaHadronicMass>
describeHerwigBtoSGammaKagan("Herwig::BtoSGammaKagan",
			     "HwFormFactors.so");

HERWIG_INTERPOLATOR_CLASSDESC(BtoSGammaKagan1,double,double)

HERWIG_INTERPOLATOR_CLASSDESC(BtoSGammaKagan2,InvEnergy,Energy)


BtoSGammaKagan::BtoSGammaKagan() 
  : _initialize(false),_mt(175.*GeV),_mb(4.8*GeV),
    _mc(1.392*GeV),_ms(ZERO),_msovermb(1./50.),_zratio(0.),
    _lambda2(0.12*GeV2),_mw(80.425*GeV),_mz(91.1876*GeV),
    _MB(5279.4*MeV),_c20(0.),_c70(0.),_c80(0.),
    _beta0(23./3.),_beta1(116./3.),_alpha(1./137.036),
    _alphaSZ(0.118),_mub(4.8*GeV),_alphaSM(0.),
    _ckm(0.976),_delta(0.),_spectmax(0.00025/GeV),_maxtry(100),
    _fermilambda(ZERO),_fermia(0.),_ferminorm(1./GeV),
    _fermilambda1(-0.3*GeV2),_ycut(0.9999999999),
    _y(0.),_deltacut(0.9),_nsfunct(100),_nspect(100),_iopt(9999) {
  Energy mHin[100]={0*GeV,0.0505907*GeV,0.101181*GeV,0.151772*GeV,0.202363*GeV,
		    0.252953*GeV,0.303544*GeV,0.354135*GeV,0.404726*GeV,0.455316*GeV,
		    0.505907*GeV,0.556498*GeV,0.607088*GeV,0.657679*GeV,0.70827*GeV,
		    0.75886*GeV,0.809451*GeV,0.860042*GeV,0.910632*GeV,0.961223*GeV,
		    1.01181*GeV,1.0624*GeV,1.113*GeV,1.16359*GeV,1.21418*GeV,
		    1.26477*GeV,1.31536*GeV,1.36595*GeV,1.41654*GeV,1.46713*GeV,
		    1.51772*GeV,1.56831*GeV,1.6189*GeV,1.66949*GeV,1.72008*GeV,
		    1.77067*GeV,1.82126*GeV,1.87186*GeV,1.92245*GeV,1.97304*GeV,
		    2.02363*GeV,2.07422*GeV,2.12481*GeV,2.1754*GeV,2.22599*GeV,
		    2.27658*GeV,2.32717*GeV,2.37776*GeV,2.42835*GeV,2.47894*GeV,
		    2.52953*GeV,2.58013*GeV,2.63072*GeV,2.68131*GeV,2.7319*GeV,
		    2.78249*GeV,2.83308*GeV,2.88367*GeV,2.93426*GeV,2.98485*GeV,
		    3.03544*GeV,3.08603*GeV,3.13662*GeV,3.18721*GeV,3.2378*GeV,
		    3.2884*GeV,3.33899*GeV,3.38958*GeV,3.44017*GeV,3.49076*GeV,
		    3.54135*GeV,3.59194*GeV,3.64253*GeV,3.69312*GeV,3.74371*GeV,
		    3.7943*GeV,3.84489*GeV,3.89548*GeV,3.94607*GeV,3.99666*GeV,
		    4.04726*GeV,4.09785*GeV,4.14844*GeV,4.19903*GeV,4.24962*GeV,
		    4.30021*GeV,4.3508*GeV,4.40139*GeV,4.45198*GeV,4.50257*GeV,
		    4.55316*GeV,4.60375*GeV,4.65434*GeV,4.70493*GeV,4.75553*GeV,
		    4.80612*GeV,4.85671*GeV,4.9073*GeV,4.95789*GeV,5.00848*GeV};
  InvEnergy spin[100]={0./GeV,3.40885e-10/GeV,1.06258e-08/GeV,7.30539e-08/GeV,
		       2.75462e-07/GeV,7.49796e-07/GeV,1.66662e-06/GeV,3.22116e-06/GeV,
		       5.62241e-06/GeV,9.06604e-06/GeV,1.37419e-05/GeV,1.98035e-05/GeV,
		       2.73352e-05/GeV,3.6404e-05/GeV,4.69896e-05/GeV,5.90282e-05/GeV,
		       7.23334e-05/GeV,8.67468e-05/GeV,0.000101999/GeV,0.000117792/GeV,
		       0.000133824/GeV,0.0001497/GeV,0.000165089/GeV,0.000179618/GeV,
		       0.000193034/GeV,0.000205017/GeV,0.000215324/GeV,0.000223761/GeV,
		       0.000230224/GeV,0.00023456/GeV,0.000236774/GeV,0.000236858/GeV,
		       0.000234965/GeV,0.000231179/GeV,0.00022566/GeV,0.000218597/GeV,
		       0.000210199/GeV,0.000200691/GeV,0.000190323/GeV,0.000179277/GeV,
		       0.000167797/GeV,0.000156088/GeV,0.000144322/GeV,0.000132705/GeV,
		       0.000121364/GeV,0.00011042/GeV,9.99745e-05/GeV,9.01017e-05/GeV,
		       8.08564e-05/GeV,7.22729e-05/GeV,6.43679e-05/GeV,5.71471e-05/GeV,
		       5.05892e-05/GeV,4.46771e-05/GeV,3.93802e-05/GeV,3.46595e-05/GeV,
		       3.04797e-05/GeV,2.67948e-05/GeV,2.3561e-05/GeV,2.07344e-05/GeV,
		       1.82726e-05/GeV,1.61351e-05/GeV,1.42839e-05/GeV,1.26838e-05/GeV,
		       1.13031e-05/GeV,1.01119e-05/GeV,9.08476e-06/GeV,8.19854e-06/GeV,
		       7.43307e-06/GeV,6.77062e-06/GeV,6.19614e-06/GeV,5.69634e-06/GeV,
		       5.26005e-06/GeV,4.87775e-06/GeV,4.54138e-06/GeV,4.24416e-06/GeV,
		       3.9805e-06/GeV,3.74572e-06/GeV,3.53602e-06/GeV,3.34835e-06/GeV,
		       3.18023e-06/GeV,3.02972e-06/GeV,2.89566e-06/GeV,2.77756e-06/GeV,
		       2.67492e-06/GeV,2.58796e-06/GeV,2.52423e-06/GeV,2.51738e-06/GeV,
		       2.58615e-06/GeV,2.72927e-06/GeV,2.94112e-06/GeV,3.22002e-06/GeV,
		       3.57015e-06/GeV,4.00189e-06/GeV,4.53288e-06/GeV,5.19051e-06/GeV,
		       6.0167e-06/GeV,7.0773e-06/GeV,8.4804e-06/GeV,1.04156e-05/GeV};
  _mHinter=vector<Energy>(mHin,mHin+100);
  _spectrum=vector<InvEnergy>(spin,spin+100);
}

void BtoSGammaKagan::doinitrun() {
  BtoSGammaHadronicMass::doinitrun();
  _pmHinter = new_ptr(Interpolator<InvEnergy,Energy>(_spectrum,_mHinter,3));
}

void BtoSGammaKagan::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mt,GeV) << ounit(_mb,GeV) << ounit(_mc,GeV) << ounit(_ms,GeV) 
     << _msovermb << _zratio << ounit(_lambda2,GeV2) << ounit(_mw,GeV) << ounit(_mz,GeV) 
     << ounit(_MB,GeV) << _c20 << _c70 << _c80 << _beta0 << _beta1 << _alpha 
     << _alphaSZ << ounit(_mub,GeV) << _alphaSM << _ckm << _delta 
     << ounit(_mHinter,GeV) << ounit(_spectrum,1./GeV) 
     << ounit(_spectmax,1./GeV) << ounit(_fermilambda,GeV)
     << _fermia << ounit(_ferminorm,1./GeV) << ounit(_fermilambda1,GeV2)
     <<_ycut << _deltacut << _nsfunct 
     << _nspect << _maxtry << _initialize
     << _s22inter << _s27inter << _pmHinter;
}

void BtoSGammaKagan::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mt,GeV) >> iunit(_mb,GeV) >> iunit(_mc,GeV) >> iunit(_ms,GeV) 
     >> _msovermb >> _zratio >> iunit(_lambda2,GeV2) >> iunit(_mw,GeV) >> iunit(_mz,GeV) 
     >> iunit(_MB,GeV) >> _c20 >> _c70 >> _c80 >> _beta0 >> _beta1 >> _alpha 
     >> _alphaSZ >> iunit(_mub,GeV) >> _alphaSM >> _ckm >> _delta 
     >> iunit(_mHinter,GeV) >> iunit(_spectrum,1./GeV)
     >> iunit(_spectmax,1./GeV) >> iunit(_fermilambda,GeV) 
     >> _fermia >> iunit(_ferminorm,1./GeV) >> iunit(_fermilambda1,GeV2)
     >>_ycut >> _deltacut >> _nsfunct 
     >> _nspect >> _maxtry >> _initialize
     >> _s22inter >> _s27inter >> _pmHinter;
}

void BtoSGammaKagan::Init() {

  static ClassDocumentation<BtoSGammaKagan> documentation
    ("The BtoSGammaKagan class implements the calculation of hep-ph/9805303 for the"
     " hadronic mass spectrum in b to s gamma decays.",
     "The decay $b\\to s\\gamma$ was simulated using the hadronic mass spectrum from"
     "\\cite{Kagan:1998ym}.\n",
     "\\bibitem{Kagan:1998ym} A.~L.~Kagan and M.~Neubert,\n"
     "Eur.\\ Phys.\\ J.\\  C {\\bf 7} (1999) 5 [arXiv:hep-ph/9805303].\n"
     "%%CITATION = EPHJA,C7,5;%%\n");

  static Switch<BtoSGammaKagan,bool> interfaceInitialize
    ("Initialize",
     "Initialize the interpolation tables for the hadronic mass",
     &BtoSGammaKagan::_initialize, false, false, false);
  static SwitchOption interfaceInitializeInitialize
    (interfaceInitialize,
     "Yes",
     "Perform initialization",
     true);
  static SwitchOption interfaceInitializeNoInitialization
    (interfaceInitialize,
     "No",
     "No initialization is performed",
     false);

  static Parameter<BtoSGammaKagan,Energy> interfaceTopMass
    ("TopMass",
     "The mass of the top quark",
     &BtoSGammaKagan::_mt, GeV, 175.*GeV, 100.0*GeV, 200.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceBottomMass
    ("BottomMass",
     "The mass of the bottom quark",
     &BtoSGammaKagan::_mb, GeV, 4.8*GeV, 1.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceCharmMass
    ("CharmMass",
     "The mass of the charm quark",
     &BtoSGammaKagan::_mc, GeV, 1.392*GeV, 1.0*GeV, 2.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceStrangeMassRatio
    ("StrangeMassRatio",
     "The ratio of the strange quark mass to that of the bottom",
     &BtoSGammaKagan::_msovermb, 1./50., 0.001, 0.1,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceWMass
    ("WMass",
     "The mass of the W boson",
     &BtoSGammaKagan::_mw, GeV, 80.425*GeV, 75.*GeV, 85.*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceZMass
    ("ZMass",
     "The mass of the Z boson",
     &BtoSGammaKagan::_mz, GeV, 91.1876*GeV, 85.*GeV, 95.*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy2> interfaceLambda2
    ("Lambda2",
     "Hadronic parameter from hep-ph/9805303",
     &BtoSGammaKagan::_lambda2, GeV2, 0.12*GeV2, ZERO, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceBMesonMass
    ("BMesonMass",
     "The mass of the decaying B meson",
     &BtoSGammaKagan::_MB, GeV, 5.2794*GeV, 5.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy> interfaceMu
    ("Mu",
     "The renormalisation scale",
     &BtoSGammaKagan::_mub, GeV, 4.8*GeV, 3.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceDelta
    ("Delta",
     "The cut off on the photon energy",
     &BtoSGammaKagan::_deltacut, 0.9, 0.0, 1.0,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,Energy2> interfaceLambda1
    ("Lambda1",
     "Hadronic scale for the fermi motion function",
     &BtoSGammaKagan::_fermilambda1, GeV2, -0.3*GeV2, -10.0*GeV2, 10.0*GeV2,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfacealpha
    ("alpha",
     "The fine structure constant",
     &BtoSGammaKagan::_alpha, 1./137.036, 0.005, 0.01,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfacealphaSZ
    ("alphaSZ",
     "The strong coupling at the Z mass",
     &BtoSGammaKagan::_alphaSZ, 0.118, 0.1, 0.3,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceCKM
    ("CKM",
     "The CKM prefactor for the decay",
     &BtoSGammaKagan::_ckm, 0.976, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<BtoSGammaKagan,Energy> interfacemHValues
    ("mHValues",
     "The mH values for the interpolation of the spectrum",
     &BtoSGammaKagan::_mHinter, GeV, -1, 1.*GeV, ZERO, 10.*GeV,
     false, false, Interface::limited);

  static ParVector<BtoSGammaKagan,InvEnergy> interfaceSpectrum
    ("Spectrum",
     "Values of the spectrum for interpolation",
     &BtoSGammaKagan::_spectrum, 1./GeV, -1, ZERO, ZERO, 1./GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,InvEnergy> interfaceFermiNormalisation
    ("FermiNormalisation",
     "The normalisation for the fermi motion function",
     &BtoSGammaKagan::_ferminorm, 1./GeV, 1./GeV, 0*1./GeV, 0*1./GeV,
     false, false, Interface::nolimits);

  static Parameter<BtoSGammaKagan,unsigned int> interfaceMaximumTries
    ("MaximumTries",
     "Maximum number of attempts to unweight the spectrum",
     &BtoSGammaKagan::_maxtry, 100, 10, 10000,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,InvEnergy> interfaceSpectrumMaximum
    ("SpectrumMaximum",
     "The maximum value of the spectrum for unweighting",
     &BtoSGammaKagan::_spectmax, 1./GeV, 1./GeV, ZERO, 10000.0/GeV,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,double> interfaceycut
    ("ycut",
     "Limit of the value of y to avoid singualarities in integrals",
     &BtoSGammaKagan::_ycut, 0.9999999999, 0.999, 1.,
     false, false, Interface::limited);

  static Parameter<BtoSGammaKagan,unsigned int> interfaceNSpectrum
    ("NSpectrum",
     "Number of spectrum points to calculate for interpolation",
     &BtoSGammaKagan::_nspect, 100, 10, 1000,
     false, false, Interface::limited);
}

void BtoSGammaKagan::calculateWilsonCoefficients() {
  // strong coupling at various scales
  double alphaSMW(alphaS(_mw)),alphaSMT(alphaS(_mt));
  _alphaSM=alphaS(_mub);
  // coupling running top mass and ratio to W mass
  Energy mtbar=_mt*pow(alphaSMW/alphaSMT,12./23.)*
    (1.+12./23.*(253./72.-58./56.)/pi*(alphaSMW-alphaSMT)-4./3.*alphaSMT/pi);
  double xt(mtbar/_mw);xt*=xt;
  // leading order coefficients
  double eta(alphaSMW/_alphaSM);
  _c20=0.5*(pow(eta,-12./23.)+pow(eta,6./23.));
  double c7w=1./pow(xt-1.,3)*(0.25*xt*xt*(3.*xt-2.)/(xt-1)*log(xt)
			      +xt/24.*(7.-5.*xt-8.*xt*xt));
  double c8w=0.25*xt/pow(xt-1.,3)*(-3.*xt/(xt-1.)*log(xt)+0.5*(2.+5.*xt-xt*xt));
  // corrections
  double c7c=626126./272277.*pow(eta,14./23.)-56281./51730.*pow(eta,16./23.)
    -3./7.*pow(eta,6./23.)-1./14.*pow(eta,-12./23.)-0.6494*pow(eta,0.4086)
    -0.0380*pow(eta,-.4230)-0.0186*pow(eta,-0.8994)-0.0057*pow(eta,0.1456);
  double c8c=313063./363036.*pow(eta,14./23.)-0.9135*pow(eta,0.4086)
    +0.0873*pow(eta,-0.4230)-0.0571*pow(eta,-0.8994)+0.0209*pow(eta,0.1456);
  // total
  _c70=pow(eta,16./23.)*c7w+8./3.*(pow(eta,14./23.)-pow(eta,16./23.))*c8w+c7c;
  _c80=pow(eta,14./23.)*c8w+c8c;
  // corrections
  double sp(real(Li2(Complex(1.-1./xt))));
  double xt2(xt*xt),xt3(xt2*xt),xt4(xt3*xt),xt5(xt4*xt);
  double c71w=1./pow(xt-1.,4)*( sp/9.*xt*(-8.+80.*xt-122.*xt2-16.*xt3)
				+pow(log(xt),2)/3./(xt-1.)*xt2*(-28.+46.*xt+6.*xt2)
				+log(xt)/81./(xt-1.)*(208.-1364.*xt+3244.*xt2
						      -2262.*xt3-588.*xt4-102.*xt5)
				+1./486.*(-436.+2509.*xt-10740.*xt2+12205.*xt3+1646.*xt4)
				);
  double c81w=1./pow(xt-1.,4)*(sp/6.*(xt+41.*xt2+40.*xt3-4.*xt4)
			       +0.5*pow(log(xt),2)/(xt-1.)*(-17.*xt3-31.*xt2)
			       +log(xt)/216./(xt-1.)*(280.-1994.*xt+2857.*xt2
						      +4893.*xt3+1086.*xt4-210.*xt5)
			       +1./1296.*(-508.+610.*xt-28209.*xt2-14102.*xt3+737.*xt4));
  // coefficients for the series in eta
  double ai[8]={14./23.,16./23.,6./23.,-12./23.,0.4086,-0.4230,-0.8994,0.1456};
  double ei[8]={4661194./816831.,-8516./2217.,0.,0.,-1.9043,-0.1008,0.1216,0.0183};
  double fi[8]={-17.3023,8.5027,4.5508,0.7519,2.0040,0.7476,-0.5385,0.0914};
  double gi[8]={14.8088,-10.8090,-0.8740,0.4218,-2.9347,0.3971,0.1600,0.0225};
  double Ex=-2./3.*log(xt)+log(xt)/6./pow(1.-xt,4)*xt2*(15.-16.*xt+4.*xt2)
    +xt/12./pow(1.-xt,3)*(18.-11.*xt-xt2);
  double c71c=(297664./14283.*pow(eta,16./23.)-7164416./357075.*pow(eta,14./23.)
	       +256868./14283.*pow(eta,37./23.)-6698884./357075.*pow(eta,39./23.))*c8w
    +37208./4761.*(pow(eta,39./23.)-pow(eta,16./23.))*c7w;
  for(unsigned int ix=0;ix<8;++ix) {
    c71c+=(ei[ix]*eta*Ex+fi[ix]+gi[ix]*eta)*pow(eta,ai[ix]);
  }
  // total correction
  double c71(pow(eta,39./23.)*c71w+8./3.*(pow(eta,37./23.)-pow(eta,39./23.))*c81w+c71c);
  // electromagnetic corrections
  double c7em=
    (32./75.*pow(eta,-9./23.)-40./69.*pow(eta,-7./23.)+88./575.*pow(eta,16./23.))*c7w
    +(-32./575.*pow(eta,-9./23.)+32./1449.*pow(eta,-7./23.)+640./1449.*pow(eta,14./23.)
      -704./1725.*pow(eta,16./23.))*c8w
    -190./8073.*pow(eta,-35./23.)-359./3105.*pow(eta,-17./23.)
    +4276./121095.*pow(eta,-12./23.)+350531./1009125.*pow(eta,-9./23.)
    +2./4347.*pow(eta,-7./23.)-5956./15525.*pow(eta,6./23.)
    +38380./169533.*pow(eta,14./23.)-748./8625.*pow(eta,16./23.);
  // nlo correction to semi-leptonic decay
  double kSL(12./23.*(1./eta -1.));
  // corrections 
  double r7(-10./3.-8.*pi*pi/9.),gam27(416./81.),gam77(32./3.),gam87(-32./9.);
  double fz(semiLeptonicf());
  double kappa(3.382-4.14*(sqrt(_zratio)-0.29));
  double realr2(-4.092+12.78*(sqrt(_zratio)-0.29));
  double realr8(44./9.-8.*pi*pi/27.);
  // correction to the various terms (get multiplied by Delta(y))
  double delta77=1.+0.5*_alphaSM/pi*(r7+gam77*log(_mb/_mub)-16./3.)
    +(pow(1.-_zratio,4)/fz-1.)*6.*_lambda2/_mb/_mb
    +0.5*_alphaSM/pi*kappa;
  double delta27=0.5*_alphaSM/pi*(realr2+gam27*log(_mb/_mub))-_lambda2/9./_mc/_mc;
  double delta78=0.5*_alphaSM/pi*(realr8+gam87*log(_mb/_mub));
  // combination for the genuine NLO terms
  _delta = delta77*_c70*_c70+delta27*_c20*_c70+delta78*_c70*_c80
    +0.5*_alphaSM/pi*c71*_c70
    +_alpha/_alphaSM*(2.*c7em*_c70-kSL*_c70*_c70);
}

void BtoSGammaKagan::doinit() {
  BtoSGammaHadronicMass::doinit();
  // quark masses
  _ms=_msovermb*_mb;
  _zratio=sqr(_mc/_mb);
  // parameters for the fermi motion function
  _fermilambda=_MB-_mb;
  _fermia=-3.*sqr(_fermilambda)/_fermilambda1-1.;
  if(_initialize) {
    // calculate the wilson coefficients etc.
    calculateWilsonCoefficients();
    // calculate the interpolation tables for the s functions
    vector<double> sfunct,yvalues;
    // s_22 function
    sfunct.push_back(0.),yvalues.push_back(0.);
    double step(1./_nsfunct);
    _y=-0.5*step;
    // perform the integrals
    GaussianIntegrator integrator;
    _iopt=0;
    for(unsigned int ix=0;ix<_nsfunct;++ix) {
      _y+=step;
      yvalues.push_back(_y);
      sfunct.push_back(integrator.value(*this,0.,_y));
    }
    _s22inter = new_ptr(Interpolator<double,double>(sfunct,yvalues,3));
    // s_27 function;
    _iopt=1;
    sfunct.clear();yvalues.clear();
    sfunct.push_back(0.),yvalues.push_back(0.);
    _y=-0.5*step;
    // perform integrals
    for(unsigned int ix=0;ix<_nsfunct;++ix) {
      _y+=step;
      yvalues.push_back(_y);
      sfunct.push_back(integrator.value(*this,0.,_y));
    }
    _s27inter = new_ptr(Interpolator<double,double>(sfunct,yvalues,3));
    // compute the normalisation constant
    KaganIntegrand integrand(this);
    _iopt=0;
    _ferminorm*=1./integrator.value(integrand,_MB*(1.-_deltacut)-_mb,_MB-_mb);
    // now for the spectrum
    _mHinter.clear();
    _spectrum.clear();
    _spectmax=ZERO;
    // limits on the mass
    Energy minegamma(0.5*_MB*(1. - _deltacut)),maxegamma(0.5*_MB);
    Energy minhadronmass(max(minMass(),sqrt(_MB*_MB-2.*_MB*maxegamma)));
    Energy maxhadronmass(min(maxMass(),sqrt(_MB*_MB-2.*_MB*minegamma)));
    Energy hstep=(maxhadronmass-minhadronmass)/(_nspect-1);
    Energy mhadron(minhadronmass);
    // function to be integrated
    _iopt=1;
    // prefactor
    InvEnergy2 pre(6.*0.105*2./_MB/_MB*_alpha/pi/semiLeptonicf()*_ckm*_ckm);
    // compute the table
    for(unsigned int ix=0;ix<_nspect;++ix) {
      // calculate y
      _y=1.-mhadron*mhadron/_MB/_MB;
      // perform the integral
      _spectrum.push_back(pre*mhadron*integrator.value(integrand,_MB*_y-_mb,_MB-_mb));
      _spectmax=max(_spectmax,_spectrum.back());
      _mHinter.push_back(mhadron);
      // increment the loop
      mhadron+=hstep;
    }
  }
}

Energy BtoSGammaKagan::hadronicMass(Energy mb,Energy mquark) {
  useMe();
  Energy minmass(max(minMass(),mquark)),maxmass(min(maxMass(),mb)),mass;
  Energy minegamma(0.5*_MB*(1. - _deltacut));//,maxegamma(0.5*_MB);
  minmass=max(minmass,ZERO); // ZERO==sqrt(_MB*_MB-2.*_MB*maxegamma));
  maxmass=min(maxmass,sqrt(_MB*_MB-2.*_MB*minegamma));
  unsigned int ntry(0);
  double rand;
  do {
    rand=UseRandom::rnd();
    mass = minmass*(1.-rand)+rand*maxmass;
    ++ntry;
  }
  while(ntry<_maxtry&&(*_pmHinter)(mass)<UseRandom::rnd()*_spectmax);
  if(ntry>=_maxtry) throw Exception() 
    << "Unweighting failed in BtoSGammaKagan::hadronicMass()" 
    << Exception::eventerror;
  return mass;
}

void BtoSGammaKagan::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::BtoSGammaKagan " << name() << " \n";
  output << "newdef " << name() << ":TopMass "    << _mt/GeV << " \n";
  output << "newdef " << name() << ":BottomMass " << _mb/GeV << " \n";
  output << "newdef " << name() << ":CharmMass "  << _mc/GeV << " \n";
  output << "newdef " << name() << ":StrangeMassRatio " << _msovermb << " \n";
  output << "newdef " << name() << ":WMass " << _mw/GeV << " \n";
  output << "newdef " << name() << ":ZMass " << _mz/GeV << " \n";
  output << "newdef " << name() << ":Lambda2 " << _lambda2/GeV2 << " \n";
  output << "newdef " << name() << ":BMesonMass " << _MB/GeV << " \n";
  output << "newdef " << name() << ":Mu " << _mub/GeV << " \n";
  output << "newdef " << name() << ":Delta " << _deltacut << " \n";
  output << "newdef " << name() << ":Lambda1 " << _fermilambda1/GeV2 << " \n";
  output << "newdef " << name() << ":alpha " << _alpha << " \n";
  output << "newdef " << name() << ":CKM " << _ckm << " \n";
  output << "newdef " << name() << ":FermiNormalisation " << _ferminorm*GeV << " \n";
  output << "newdef " << name() << ":MaximumTries " << _maxtry << " \n";
  output << "newdef " << name() << ":ycut " << _ycut << " \n";
  output << "newdef " << name() << ":NSpectrum " <<  _nspect << " \n";
  for(unsigned int ix=0;ix<_mHinter.size();++ix) {
    if(ix<100) output << "newdef " << name() << ":mHValues " << ix << " " 
		      << _mHinter[ix]/GeV << " \n";
    else       output << "insert " << name() << ":mHValues " << ix << " " 
		      << _mHinter[ix]/GeV << " \n";
  }
  for(unsigned int ix=0;ix<_spectrum.size();++ix) {
    if(ix<100) output << "newdef " << name() << ":Spectrum " << ix << " " 
		      << _spectrum[ix]*GeV << " \n";
    else       output << "insert " << name() << ":Spectrum " << ix << " " 
		      << _spectrum[ix]*GeV << " \n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./HQETFormFactor.cc"
// -*- C++ -*-
//
// HQETFormFactor.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HQETFormFactor class.
//

#include "HQETFormFactor.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

HQETFormFactor::HQETFormFactor() {
  // parameters from arXiv:0705.4008
  // F(1) from F(1)V_cb and our value of V_cb
  _f1scalar = 1.0269328;
  _f1vector = 0.84;
  // slope
  _rho2scalar = 1.17;
  _rho2vector = 1.179;
  // R_1(1)
  _r1 = 1.417;
  // R_2(1)
  _r2 = 0.836;
  // allowed form factors
  addFormFactor(-521, 421  ,0,-2, 5, 4);
  addFormFactor(-521, 423  ,1,-2, 5, 4);
  addFormFactor(-511, 411  ,0, 1, 5, 4);
  addFormFactor(-511, 413  ,1, 1, 5, 4);
  // set the initial number of modes
  initialModes(numberOfFactors());
}

void HQETFormFactor::persistentOutput(PersistentOStream & os) const {
  os << _f1scalar << _f1vector << _r1 << _r2 << _rho2scalar << _rho2vector;
}

void HQETFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> _f1scalar >> _f1vector >> _r1 >> _r2 >> _rho2scalar >> _rho2vector;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<HQETFormFactor,ScalarFormFactor>
describeHerwigHQETFormFactor("Herwig::HQETFormFactor", "HwFormFactors.so");

void HQETFormFactor::Init() {

  static ClassDocumentation<HQETFormFactor> documentation
    ("The HQETFormFactor class uses the parameterisation of hep-ph/9712417"
     " of the form factor in the heavy quark limit.",
     "The parameterisation of \\cite{Caprini:1997mu} was used for the "
     "$B\\to D^{(*)}$ form factors together with the parameters from"
     "\\cite{Snyder:2007qn} for the $D$ and \\cite{Aubert:2007rs}"
     " for the $D^*$",
     "\\bibitem{Caprini:1997mu} I.~Caprini, L.~Lellouch and M.~Neubert,"
     "Nucl.\\ Phys.\\  B {\\bf 530} (1998) 153 [arXiv:hep-ph/9712417].\n"
     "%%CITATION = NUPHA,B530,153;%%\n"
     "\\bibitem{Aubert:2007rs} B.~Aubert {\\it et al.}  [BABAR Collaboration],"
     "arXiv:0705.4008 [hep-ex]. %%CITATION = ARXIV:0705.4008;%%\n"
     "\\bibitem{Snyder:2007qn} A.~E.~Snyder, [arXiv:hep-ex/0703035].\n"
     "%%CITATION = ECONF,C0610161,015;%%\n");

  static Parameter<HQETFormFactor,double> interfaceF1Scalar
    ("F1Scalar",
     "The normalisation factor for the scalar form factor",
     &HQETFormFactor::_f1scalar, 1.0269328, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceF1Vector
    ("F1Vector",
     "The normalisation factor for the vector form factor",
     &HQETFormFactor::_f1vector, 0.84,  0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceRho2Scalar
    ("Rho2Scalar",
     "The slope parameter for the scalar form factor",
     &HQETFormFactor::_rho2scalar, 1.17, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceRho2Vector
    ("Rho2Vector",
     "The slope parameter for the vector form factor",
     &HQETFormFactor::_rho2vector, 1.179, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceR1
    ("R1",
     "The ratio R_1 at omega=1",
     &HQETFormFactor::_r1, 1.417, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<HQETFormFactor,double> interfaceR2
    ("R2",
     "The ratio R_2 at omega=1",
     &HQETFormFactor::_r2, 0.836, 0.0, 10.0,
     false, false, Interface::limited);

}

void HQETFormFactor::ScalarScalarFormFactor(Energy2 q2,unsigned int ,int,int,
					    Energy m0,Energy m1,Complex & f0,
					    Complex & fp) const {
  useMe();
  double omega = 0.5*(sqr(m0)+sqr(m1)-q2)/m0/m1;
  double root = sqrt(1.+omega),rt2=sqrt(2.);
  double z = (root-rt2)/(root+rt2);
  double Rs = 2.*sqrt(m0*m1)/(m0+m1);
  fp = 1.-8.*_rho2scalar*z+((51.*_rho2scalar-10.)-(252*_rho2scalar-84.)*z)*sqr(z);
  fp *=_f1scalar/Rs;
  f0 = fp*(1.-q2/sqr(m0+m1));
}

void HQETFormFactor::ScalarVectorFormFactor(Energy2 q2, unsigned int,
					    int, int, Energy m0, Energy m1,
					    Complex & A0, Complex & A1,Complex & A2,
					    Complex & V) const {
  useMe();
  double omega = 0.5*(sqr(m0)+sqr(m1)-q2)/m0/m1;
  double root = sqrt(1.+omega),rt2=sqrt(2.);
  double z = (root-rt2)/(root+rt2);
  double hA1 = _f1vector*(1.-8.*_rho2vector*z+((53.*_rho2vector-15.)
					       -(231.*_rho2vector-91.)*z)*sqr(z));
  double wmo=omega-1.;
  double R1 = _r1-0.12*wmo+0.05*sqr(wmo);
  double R2 = _r2+0.11*wmo-0.06*sqr(wmo);
  double Rs = 2.*sqrt(m0*m1)/(m0+m1);
  A1 = 0.5*(omega+1.)*Rs*hA1;
  A2 = R2*hA1/Rs;
  V  =-R1*hA1/Rs;
  Complex A3 = 0.5/m1*((m0+m1)*A1-(m0-m1)*A2);
  A0 = A3+0.5*A2*q2/m1/(m0+m1);
}

void HQETFormFactor::dataBaseOutput(ofstream & os,bool header,bool create) const {
  if(header) os << "update decayers set parameters=\"";
  if(create) os << "create Herwig::HQETFormFactor " << name() << "\n";
  ScalarFormFactor::dataBaseOutput(os,false,false);
  os << "newdef " << name() << ":F1Scalar   " << _f1scalar   << "\n";
  os << "newdef " << name() << ":F1Vector   " << _f1vector   << "\n";
  os << "newdef " << name() << ":Rho2Scalar " << _rho2scalar << "\n";
  os << "newdef " << name() << ":Rho2Vector " << _rho2vector << "\n";
  os << "newdef " << name() << ":R1         " << _r1         << "\n";
  os << "newdef " << name() << ":R2         " << _r2         << "\n";
  if(header) os << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./CzyzNucleonFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CzyzNucleonFormFactor class.
//

#include "CzyzNucleonFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/ResonanceHelpers.h"

using namespace Herwig;

CzyzNucleonFormFactor::CzyzNucleonFormFactor() {
  // masses and widths of rho resonances
  rhoMasses_   = {775.49*MeV, 1465*MeV, 1720*MeV, 2.12*GeV,2.32647*GeV};
  rhoWidths_   = {149.10*MeV,  400*MeV,  250*MeV, 0.3 *GeV,0.4473 *GeV};
  // masses and width of omega resonances
  omegaMasses_ = {782.65*MeV, 1425*MeV, 1670*MeV, 2.0707 *GeV, 2.34795 *GeV};
  omegaWidths_ = {8.49  *MeV,  215*MeV,  315*MeV, 1.03535*GeV, 1.173975*GeV};
  // c_1 couplings
  c1Re_ = {1.,-0.467,-0.177, 0.301};
  c1Im_ = {0.,-0.385, 0.149, 0.264};
  // c_2 couplings
  c2Re_ = {1., 0.0521, -0.00308,-0.348};
  c2Im_ = {0.,-3.04, 2.38,-0.104};
  // c_3 couplings
  c3Re_ = {1.,-7.88,10.2};
  c3Im_ = {0., 5.67, -1.94};
  // c_4 couplings
  c4Re_ = {1.,-0.832, 0.405};
  c4Im_ = {0., 0.308,-0.25};
  // Magnetic moments
  mup_ =  2.792847356;
  mun_ = -1.9130427;
  // initialise a_ and b_
  a_ = mup_-mun_-1.;
  b_ = -mup_-mun_+1.;
  // set up the form factors
  addFormFactor(2212,2212,2,2,2,1,1,1);
  addFormFactor(2112,2112,2,2,2,1,1,1);
  initialModes(numberOfFactors());
}

IBPtr CzyzNucleonFormFactor::clone() const {
  return new_ptr(*this);
}

IBPtr CzyzNucleonFormFactor::fullclone() const {
  return new_ptr(*this);
}

void CzyzNucleonFormFactor::doinit() {
  BaryonFormFactor::doinit();
  static const Complex ii(0.,1.);
  // calculate c_1
  c1_.clear();
  assert(c1Re_.size()==4 && c1Im_.size()==4);
  // c_1 1 -> 4
  complex<Energy2> fact(ZERO);
  for(unsigned int ix=0;ix<4;++ix) {
    c1_.push_back(c1Re_[ix]+ii*c1Im_[ix]);
    fact += c1_[ix]*sqr(omegaMasses_[ix]);
  }
  c1_.push_back(-fact/sqr(omegaMasses_[4]));
  // calculate c_2
  c2_.clear();
  assert(c2Re_.size()==4 && c2Im_.size()==4);
  // c_2 1 -> 4
  fact = ZERO;
  for(unsigned int ix=0;ix<4;++ix) {
    c2_.push_back(c2Re_[ix]+ii*c2Im_[ix]);
    fact += c2_[ix]*sqr(rhoMasses_[ix]);
  }
  c2_.push_back(-fact/sqr(rhoMasses_[4]));
  // calculate c_3
  c3_.clear();
  assert(c3Re_.size()==3 && c3Im_.size()==3);
  // c_3 1 -> 4
  fact = ZERO;
  complex<Energy4> fact2(ZERO);
  for(unsigned int ix=0;ix<3;++ix) {
    c3_.push_back(c3Re_[ix]+ii*c3Im_[ix]);
    fact += c3_[ix]*sqr(omegaMasses_[ix]);
    fact2 += c3_[ix]*sqr(omegaMasses_[ix])*
      (sqr(omegaMasses_[ix])-sqr(omegaMasses_[4])
       +ii*(omegaMasses_[4]*omegaWidths_[4]-omegaMasses_[ix]*omegaWidths_[ix]));
  }
  c3_.push_back(fact2/sqr(omegaMasses_[3])/
		(sqr(omegaMasses_[4])-sqr(omegaMasses_[3])
		 +ii*(omegaMasses_[3]*omegaWidths_[3]-omegaMasses_[4]*omegaWidths_[4])) );
  fact += c3_[3]*sqr(omegaMasses_[3]);
  c3_.push_back(-fact/sqr(omegaMasses_[4]));
  // c_4 1 -> 4
  fact = ZERO;
  fact2 = ZERO;
  for(unsigned int ix=0;ix<3;++ix) {
    c4_.push_back(c4Re_[ix]+ii*c4Im_[ix]);
    fact += c4_[ix]*sqr(rhoMasses_[ix]);
    fact2 += c4_[ix]*sqr(rhoMasses_[ix])*
      (sqr(rhoMasses_[ix])-sqr(rhoMasses_[4])
       +ii*(rhoMasses_[4]*rhoWidths_[4]-rhoMasses_[ix]*rhoWidths_[ix]));
  }
  c4_.push_back(fact2/sqr(rhoMasses_[3])/
		(sqr(rhoMasses_[4])-sqr(rhoMasses_[3])
		 +ii*(rhoMasses_[3]*rhoWidths_[3]-rhoMasses_[4]*rhoWidths_[4])) );
  fact += c4_[3]*sqr(rhoMasses_[3]);
  c4_.push_back(-fact/sqr(rhoMasses_[4]));
  // a and b parameters
  a_ = mup_-mun_-1.;
  b_ = -mup_-mun_+1.;
}

void CzyzNucleonFormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(rhoMasses_,GeV) << ounit(rhoWidths_,GeV)
     << ounit(omegaMasses_,GeV) << ounit(omegaWidths_,GeV)
     << c1Re_ << c1Im_ << c2Re_ << c2Im_
     << c3Re_ << c3Im_ << c4Re_ << c4Im_
     << c1_ << c2_ << c3_ << c4_
     <<  mup_ << mun_ << a_ << b_;
}

void CzyzNucleonFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rhoMasses_,GeV) >> iunit(rhoWidths_,GeV)
     >> iunit(omegaMasses_,GeV) >> iunit(omegaWidths_,GeV)
     >> c1Re_ >> c1Im_ >> c2Re_ >> c2Im_
     >> c3Re_ >> c3Im_ >> c4Re_ >> c4Im_
     >> c1_ >> c2_ >> c3_ >> c4_
     >> mup_ >> mun_ >> a_ >> b_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<CzyzNucleonFormFactor,BaryonFormFactor>
describeHerwigCzyzNucleonFormFactor("Herwig::CzyzNucleonFormFactor",
				    "HwFormFactors.so");

void CzyzNucleonFormFactor::Init() {

  static ClassDocumentation<CzyzNucleonFormFactor> documentation
    ("The CzyzNucleonFormFactor class implements the model of "
     "Phys.Rev. D90 (2014) no.11, 114021 for the nucleon form factor",
     "The nucleon form factor model of \\cite{Czyz:2014sha} was used",
     "\\bibitem{Czyz:2014sha}\n"
     "H.~Czyż, J.~H.~Kühn and S.~Tracz,\n"
     "%``Nucleon form factors and final state radiative corrections to $e^+e^-  \\to p\\bar{p}γ$,''\n"
     "Phys.\\ Rev.\\ D {\\bf 90} (2014) no.11,  114021\n"
     "doi:10.1103/PhysRevD.90.114021\n"
     "[arXiv:1407.7995 [hep-ph]].\n"
     "%%CITATION = doi:10.1103/PhysRevD.90.114021;%%\n"
     "%5 citations counted in INSPIRE as of 25 Aug 2018\n");

  static ParVector<CzyzNucleonFormFactor,Energy> interfaceRhoMasses
    ("RhoMasses",
     "The masses of the rho mesons for the form factor",
     &CzyzNucleonFormFactor::rhoMasses_, GeV, 5, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,Energy> interfaceRhoWidths
    ("RhoWidths",
     "The widths of the rho mesons for the form factor",
     &CzyzNucleonFormFactor::rhoWidths_, GeV, 5, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,Energy> interfaceOmegaMasses
    ("OmegaMasses",
     "The masses of the omega mesons for the form factor",
     &CzyzNucleonFormFactor::omegaMasses_, GeV, 5, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,Energy> interfaceOmegaWidths
    ("OmegaWidths",
     "The widths of the omega mesons for the form factor",
     &CzyzNucleonFormFactor::omegaWidths_, GeV, 5, 1.0*GeV, 0.0*GeV, 10.0*GeV,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec1Real
    ("c1Real",
     "The real part of the c_1 coupling",
     &CzyzNucleonFormFactor::c1Re_, 4, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec1Imag
    ("c1Imag",
     "The imaginary part of the c_1 coupling",
     &CzyzNucleonFormFactor::c1Im_, 4, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec2Real
    ("c2Real",
     "The real part of the c_2 coupling",
     &CzyzNucleonFormFactor::c2Re_, 4, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec2Imag
    ("c2Imag",
     "The imaginary part of the c_2 coupling",
     &CzyzNucleonFormFactor::c2Im_, 4, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec3Real
    ("c3Real",
     "The real part of the c_3 coupling",
     &CzyzNucleonFormFactor::c3Re_, 3, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec3Imag
    ("c3Imag",
     "The imaginary part of the c_3 coupling",
     &CzyzNucleonFormFactor::c3Im_, 3, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec4Real
    ("c4Real",
     "The real part of the c_4 coupling",
     &CzyzNucleonFormFactor::c4Re_, 3, 1., -100., 100.0,
     false, false, Interface::limited);

  static ParVector<CzyzNucleonFormFactor,double> interfacec4Imag
    ("c4Imag",
     "The imaginary part of the c_4 coupling",
     &CzyzNucleonFormFactor::c4Im_, 3, 1., -100., 100.0,
     false, false, Interface::limited);

  static Parameter<CzyzNucleonFormFactor,double> interfaceMuProton
    ("MuProton",
     "The proton magnetic moment",
     &CzyzNucleonFormFactor::mup_, 2.792, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<CzyzNucleonFormFactor,double> interfaceMuNeutron
    ("MuNeutron",
     "The proton magnetic moment",
     &CzyzNucleonFormFactor::mun_, -1.913, 0.0, 10.0,
     false, false, Interface::limited);
}


void CzyzNucleonFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int id0,int id1,Energy,Energy,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo flavour,
			   Virtuality virt) {
  f1a = f2a = f3a = f1v = f2v = f3v = 0.;
  assert(abs(id0)==abs(id1));
  if(iloc==0) assert(abs(id0)==2212);
  else        assert(abs(id0)==2112);
  assert(virt==TimeLike);
  if(flavour.strange != Strangeness::Unknown and flavour.strange != Strangeness::Zero) return;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return;
  bool I0 = true, I1 = true;
  if(flavour.I!=IsoSpin::IUnknown) {
    if(flavour.I==IsoSpin::IZero) {
      I1=false;
      if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
    }
    else if(flavour.I==IsoSpin::IOne) {
      I0=false;
      if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
    }
  }
  if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
  // calculate the form factors
  Complex F1S(0.),F1V(0.),F2S(0.),F2V(0.);
  Complex n1(0.),n2(0.),n3(0.),n4(0.);
  for(unsigned int ix=0;ix<5;++ix) {
    if(I0) {
      F1S += c1_[ix]*Resonance::BreitWignerFW(q2,omegaMasses_[ix],omegaWidths_[ix]);
      F2S += c3_[ix]*Resonance::BreitWignerFW(q2,omegaMasses_[ix],omegaWidths_[ix]);
    }
    if(I1) {
      F1V += c2_[ix]*Resonance::BreitWignerFW(q2,  rhoMasses_[ix],  rhoWidths_[ix]);
      F2V += c4_[ix]*Resonance::BreitWignerFW(q2,  rhoMasses_[ix],  rhoWidths_[ix]);
    }
    n1 += c1_[ix];
    n2 += c2_[ix];
    n3 += c3_[ix];
    n4 += c4_[ix];
  }
  F1S *=  0.5   /n1;
  F1V *=  0.5   /n2;
  F2S *= -0.5*b_/n3;
  F2V *=  0.5*a_/n4;
  if(iloc==0) {
    f1v =  F1S + F1V;
    f2v =  F2S + F2V;
  }
  else {
    f1v = F1S - F1V;
    f2v = F2S - F2V;
  }
}

void CzyzNucleonFormFactor::
dataBaseOutput(ofstream& output,bool header,
	       bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::CzyzNucleonFormFactor "
		    << name() << " \n";
  for(unsigned int ix=0;ix<5;++ix) {
    output << "newdef " << name() << ":RhoMasses " << ix << " "
	   << rhoMasses_[ix]/GeV << "\n";
    output << "newdef " << name() << ":RhoWidths " << ix << " "
	   << rhoWidths_[ix]/GeV << "\n";
    output << "newdef " << name() << ":OmegaMasses " << ix << " "
	   << omegaMasses_[ix]/GeV << "\n";
    output << "newdef " << name() << ":OmegaWidths " << ix << " "
	   << omegaWidths_[ix]/GeV << "\n";
  }
  for(unsigned int ix=0;ix<4;++ix) {
    output << "newdef " << name() << ":c1Real " << ix << " "
	   << c1Re_[ix] << "\n";
    output << "newdef " << name() << ":c1Imag " << ix << " "
	   << c1Im_[ix] << "\n";
    output << "newdef " << name() << ":c2Real " << ix << " "
	   << c2Re_[ix] << "\n";
    output << "newdef " << name() << ":c2Imag " << ix << " "
	   << c2Im_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<3;++ix) {
    output << "newdef " << name() << ":c3Real " << ix << " "
	   << c3Re_[ix] << "\n";
    output << "newdef " << name() << ":c3Imag " << ix << " "
	   << c3Im_[ix] << "\n";
    output << "newdef " << name() << ":c4Real " << ix << " "
	   << c4Re_[ix] << "\n";
    output << "newdef " << name() << ":c4Imag " << ix << " "
	   << c4Im_[ix] << "\n";
  }
  output << "newdef " << name() << ":MuProton  " << mup_ << "\n";
  output << "newdef " << name() << ":MuNeutron " << mun_ << "\n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}
#line 1 "./KornerKurodaFormFactor.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KornerKurodaFormFactor class.
//

#include "KornerKurodaFormFactor.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

KornerKurodaFormFactor::KornerKurodaFormFactor() : includeNucleon_(true),
						   mRho_(775.26*MeV), mOmega_(782.65*MeV), mPhi_(1019.461*MeV),
						   aPrime_(1./GeV2) {
  const double o3=1./3., o6=1./6., o9=1./6.;
  c1Rho_    = {  0.5,  -0.5,    1.,   0.,    -1.,     0.5,   -0.5,  0., 0.};
  c1Omega_  = {  0.5,   0.5,    o3,   o3,     o3,      o6,     o6,  o3, 0.};
  c1Phi_    = {   0.,   0. ,   -o3,  -o3,    -o3,  -2.*o3,- 2.*o3, -o3, 0.};
  c12Rho_   = {5.*o6,  1.25, 2.*o3,   0.,     2.,    0.25,   -0.5,  0., 1.};
  c12Omega_ = {   o6, -0.25, 2.*o9,2.*o3, -2.*o3, 0.25*o3,     o6,  0., 0.};
  c12Phi_   = {   0.,    0.,    o9,   o3,    -o3,   2.*o3,  4.*o3,  1., 0.};
  mu_       = {2.792847,-1.192304,2.458,0.930949,-1.160,-1.250,-0.6507,-0.613,1.61};
  // set up the form factors
  addFormFactor(2212,2212,2,2,2,2,1,1);
  addFormFactor(2112,2112,2,2,2,1,1,1);
  addFormFactor(3222,3222,2,2,2,2,3,3);
  addFormFactor(3212,3212,2,2,2,1,3,3);
  addFormFactor(3112,3112,2,2,1,1,3,3);
  addFormFactor(3322,3322,2,2,2,3,3,3);
  addFormFactor(3312,3312,2,2,1,3,3,3);
  addFormFactor(3122,3122,2,2,2,1,3,3);
  addFormFactor(3122,3212,2,2,2,1,3,3);
  initialModes(numberOfFactors());
}

IBPtr KornerKurodaFormFactor::clone() const {
  return new_ptr(*this);
}

IBPtr KornerKurodaFormFactor::fullclone() const {
  return new_ptr(*this);
}

void KornerKurodaFormFactor::persistentOutput(PersistentOStream & os) const {
  os << ounit(mRho_,GeV) << ounit(mOmega_,GeV) <<  ounit(mPhi_, GeV)
     << ounit(aPrime_,1./GeV2) << includeNucleon_
     << c1Rho_ << c1Omega_ << c1Phi_ << c12Rho_ << c12Omega_ << c12Phi_
     << c2Rho_ << c2Omega_ << c2Phi_ << mu_;
}

void KornerKurodaFormFactor::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mRho_,GeV) >> iunit(mOmega_,GeV) >>  iunit(mPhi_, GeV)
     >> iunit(aPrime_,1./GeV2) >> includeNucleon_
     >> c1Rho_ >> c1Omega_ >> c1Phi_ >> c12Rho_ >> c12Omega_ >> c12Phi_
     >> c2Rho_ >> c2Omega_ >> c2Phi_ >> mu_;
}

void KornerKurodaFormFactor::doinit() {
  // adjust units from nuclear magnetrons
  Energy mp = getParticleData(ParticleID::pplus)->mass();
  mu_[2] *= mp/getParticleData(ParticleID::Sigmaplus)->mass();
  mu_[3] *= mp/getParticleData(ParticleID::Sigma0)->mass();
  mu_[4] *= mp/getParticleData(ParticleID::Sigmaminus)->mass();
  mu_[5] *= mp/getParticleData(ParticleID::Xi0)->mass();
  mu_[6] *= mp/getParticleData(ParticleID::Ximinus)->mass();
  mu_[7] *= mp/getParticleData(ParticleID::Lambda0)->mass();
  mu_[8] *= mp/getParticleData(ParticleID::Lambda0)->mass();
  for(unsigned int ix=0;ix<c1Rho_.size();++ix) {
    c2Rho_  .push_back(mu_[ix]*c12Rho_  [ix]-c1Rho_  [ix]);
    c2Omega_.push_back(mu_[ix]*c12Omega_[ix]-c1Omega_[ix]);
    c2Phi_  .push_back(mu_[ix]*c12Phi_  [ix]-c1Phi_  [ix]);
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KornerKurodaFormFactor,BaryonFormFactor>
describeThePEGKornerKurodaFormFactor("Herwig::KornerKurodaFormFactor",
				     "HwFormFactors.so");

void KornerKurodaFormFactor::Init() {

  static ClassDocumentation<KornerKurodaFormFactor> documentation
    ("Simple mode of the nucelon form-factors based on  Phys. Rev. D 16 (1977) 2165",
     "Simple model of the nucleon form factors based on \\cite{Korner:1976hv}.",
     "\\bibitem{Korner:1976hv}"
     "J.~G.~Korner and M.~Kuroda,"
     "%``E+ e- Annihilation Into Baryon-anti-Baryon Pairs,''"
     "Phys.\\ Rev.\\ D {\\bf 16} (1977) 2165."
     "doi:10.1103/PhysRevD.16.2165"
     "%%CITATION = doi:10.1103/PhysRevD.16.2165;%%"
     "%108 citations counted in INSPIRE as of 31 Oct 2019");

  static Parameter<KornerKurodaFormFactor,Energy> interfacemRho
    ("mRho",
     "Mass of the rho meson for the form factor",
     &KornerKurodaFormFactor::mRho_, GeV, 0.77526*GeV, 0.6*GeV, 1.0*GeV,
     false, false, Interface::limited);

  static Parameter<KornerKurodaFormFactor,Energy> interfacemOmega
    ("mOmega",
     "Mass of the omega meson for the form factor",
     &KornerKurodaFormFactor::mOmega_, GeV, 0.78265*GeV, 0.6*GeV, 1.0*GeV,
     false, false, Interface::limited);

  static Parameter<KornerKurodaFormFactor,Energy> interfacemPhi
    ("mPhi",
     "Mass of the phi meson for the form factor",
     &KornerKurodaFormFactor::mPhi_, GeV, 1.019461*GeV, 0.9*GeV, 1.5*GeV,
     false, false, Interface::limited);
  
  static Parameter<KornerKurodaFormFactor,InvEnergy2> interfacealphaPrime
    ("alphaPrime",
     "The regge slope",
     &KornerKurodaFormFactor::aPrime_, 1./GeV2, 1./GeV2, 0.5/GeV2, 1.5/GeV2,
     false, false, Interface::limited);

  
  static Switch<KornerKurodaFormFactor,bool> interfaceIncludeNucleon
    ("IncludeNucleon",
     "Whether or not to include the nucleon (i.e. proton and neutron) modes or not."
     " Often omit these as other currents do a better job.",
     &KornerKurodaFormFactor::includeNucleon_, true, false, false);
  static SwitchOption interfaceIncludeNucleonYes
    (interfaceIncludeNucleon,
     "Yes",
     "Include proton and neutron modes",
     true);
  static SwitchOption interfaceIncludeNucleonNo
    (interfaceIncludeNucleon,
     "No",
     "Don't include neutron and proton modes",
     false);

}

void KornerKurodaFormFactor::
SpinHalfSpinHalfFormFactor(Energy2 q2,int iloc, int ,int ,Energy,Energy,
			   Complex & f1v,Complex & f2v,Complex & f3v,
			   Complex & f1a,Complex & f2a,Complex & f3a,
			   FlavourInfo flavour,
			   Virtuality virt) {
  useMe();
  assert(virt==TimeLike);
  f1a = f2a = f3a = f1v = f2v = f3v = 0.;
  if(flavour.charm   != Charm::Unknown       and flavour.charm   != Charm::Zero      ) return;
  if(flavour.bottom  != Beauty::Unknown      and flavour.bottom  !=Beauty::Zero      ) return;
  if(!includeNucleon_&&iloc<=1) return;
  bool rho = true, omega = true, phi = true;
  // strange content
  if(flavour.strange != Strangeness::Unknown) {
    if(flavour.strange == Strangeness::Zero) phi = false;
    else if(flavour.strange == Strangeness::ssbar) phi = true;
    else return;
  }
  if(flavour.I!=IsoSpin::IUnknown) {
    // I=0, i.e omega content
    if(flavour.I==IsoSpin::IZero) {
      rho=false;
      if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
    }
    // I=1, i.e rho content
    else if(flavour.I==IsoSpin::IOne) {
      omega=false;
      if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
    }
  }
  if(flavour.I3!=IsoSpin::I3Unknown and flavour.I3!=IsoSpin::I3Zero) return;
  // form factors
  // rho component
  if(rho) {
    Energy2 mRho2[3] = {sqr(mRho_), mRho2[0]+1./aPrime_, mRho2[0]+2./aPrime_}; 
    f1v += c1Rho_[iloc]/(1.-q2/mRho2[0])/(1.-q2/mRho2[1]);
    f2v += c2Rho_[iloc]/(1.-q2/mRho2[0])/(1.-q2/mRho2[1])/(1.-q2/mRho2[2]);
  }
  // omega component
  if(omega) {
    Energy2 mOmega2[3] = {sqr(mOmega_), mOmega2[0]+1./aPrime_, mOmega2[0]+2./aPrime_};
    f1v += c1Omega_[iloc]/(1.-q2/mOmega2[0])/(1.-q2/mOmega2[1]);
    f2v += c2Omega_[iloc]/(1.-q2/mOmega2[0])/(1.-q2/mOmega2[1])/(1.-q2/mOmega2[2]);
  }
  // phi component
  if(phi) {
    Energy2 mPhi2[3] = {sqr(mPhi_), mPhi2[0]+1./aPrime_, mPhi2[0]+2./aPrime_}; 
    f1v += c1Phi_[iloc]/(1.-q2/mPhi2[0])/(1.-q2/mPhi2[1]);
    f2v += c2Phi_[iloc]/(1.-q2/mPhi2[0])/(1.-q2/mPhi2[1])/(1.-q2/mPhi2[2]);
  }
}

void KornerKurodaFormFactor::
dataBaseOutput(ofstream& output,bool header,
	       bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::KornerKurodaFormFactor " 
		    << name() << " \n";
  output << "newdef " << name() << ":mRho "   << mRho_/GeV   << "\n";
  output << "newdef " << name() << ":mOmega " << mOmega_/GeV << "\n";
  output << "newdef " << name() << ":mPhi "   << mPhi_/GeV   << "\n";
  output << "newdef " << name() << ":alphaPrime " << aPrime_*GeV2 << "\n";
  BaryonFormFactor::dataBaseOutput(output,false,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./KPiIThreeHalfFOCUSKMatrix.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiIThreeHalfFOCUSKMatrix class.
//

#include "KPiIThreeHalfFOCUSKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KPiIThreeHalfFOCUSKMatrix::KPiIThreeHalfFOCUSKMatrix()
  : KMatrix(FlavourInfo(IsoSpin::IThreeHalf,IsoSpin::I3Unknown,
			Strangeness::PlusOne,Charm::Zero,
			Beauty::Zero),
	    vector<Channels>({KMatrix::KPi})),
    D_({-0.22147,0.026637,-0.00092057}),
    sThreeHalf_(0.27*GeV2)
{}

IBPtr KPiIThreeHalfFOCUSKMatrix::clone() const {
  return new_ptr(*this);
}

IBPtr KPiIThreeHalfFOCUSKMatrix::fullclone() const {
  return new_ptr(*this);
}

void KPiIThreeHalfFOCUSKMatrix::persistentOutput(PersistentOStream & os) const {
  os << D_ << ounit(sThreeHalf_,GeV2) << ounit(sNorm_,GeV2);
}

void KPiIThreeHalfFOCUSKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> D_ >> iunit(sThreeHalf_,GeV2) >> iunit(sNorm_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiIThreeHalfFOCUSKMatrix,KMatrix>
describeHerwigKPiIThreeHalfFOCUSKMatrix("Herwig::KPiIThreeHalfFOCUSKMatrix", "HwFormFactors.so");

void KPiIThreeHalfFOCUSKMatrix::Init() {

  static ClassDocumentation<KPiIThreeHalfFOCUSKMatrix> documentation
    ("The KPiIThreeHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration (Phys.Lett. B653 (2007) 1-11) for the I=3/2 "
     "component of the Kpi K-matrix.",
     "The KPiIThreeHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration \\cite{Pennington:2007se} for the $I=3/2$ "
     "component of the $K\\pi$ K-matrix.",
     "\\bibitem{Pennington:2007se}"
     "J.~M.~Link {\\it et al.} [FOCUS Collaboration],"
     "%``Dalitz plot analysis of the $D^{+} \\to K^{-} \\pi^{+} \\pi^{+}$ decay in the FOCUS experiment,''"
     "Phys.\\ Lett.\\ B {\\bf 653} (2007) 1"
     "doi:10.1016/j.physletb.2007.06.070"
     "[arXiv:0705.2248 [hep-ex]]."
     "%%CITATION = doi:10.1016/j.physletb.2007.06.070;%%"
     "%79 citations counted in INSPIRE as of 14 Jan 2020");

}

void KPiIThreeHalfFOCUSKMatrix::doinit() {
  KMatrix::doinit();
  Energy mK = getParticleData(ParticleID::Kplus)->mass();
  Energy mpi= getParticleData(ParticleID::piplus)->mass();
  sNorm_ = sqr(mK)+sqr(mpi);
}

boost::numeric::ublas::matrix<double> KPiIThreeHalfFOCUSKMatrix::K(Energy2 s,bool) const {
  double st = s/sNorm_-1.;
  double param=1.;
  boost::numeric::ublas::matrix<double> output =
    boost::numeric::ublas::zero_matrix<double>(1,1);
  for(unsigned int ix=0;ix<D_.size();++ix) {
    output(0,0) += D_[ix]*param;
    param *= st;
  }
  output *=(s-sThreeHalf_)/sNorm_;
  return output;
}
#line 1 "./KPiIHalfFOCUSKMatrix.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KPiIHalfFOCUSKMatrix class.
//

#include "KPiIHalfFOCUSKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

KPiIHalfFOCUSKMatrix::KPiIHalfFOCUSKMatrix()
  : KMatrix(FlavourInfo(IsoSpin::IHalf,IsoSpin::I3Unknown,
			Strangeness::PlusOne,Charm::Zero,
			Beauty::Zero),
	    vector<Channels>({KMatrix::KPi,KMatrix::KEtaPrime}),
	    vector<Energy2>({1.7919*GeV2}),
	    vector<vector<Energy> >(1,vector<Energy>({0.31072*GeV,-0.02323*GeV}))),
    C11_({0.79299,-0.15099,0.00811}),
    C22_({0.17054,-0.0219,0.00085655}),
    C12_({0.15040,-0.038266,0.0022596}),
    sHalf_(0.23*GeV2)
{}

IBPtr KPiIHalfFOCUSKMatrix::clone() const {
  return new_ptr(*this);
}

IBPtr KPiIHalfFOCUSKMatrix::fullclone() const {
  return new_ptr(*this);
}

void KPiIHalfFOCUSKMatrix::persistentOutput(PersistentOStream & os) const {
  os << C11_ << C22_ << C12_ << ounit(sHalf_,GeV2) << ounit(sNorm_,GeV2);
}

void KPiIHalfFOCUSKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> C11_ >> C22_ >> C12_ >> iunit(sHalf_,GeV2) >> iunit(sNorm_,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KPiIHalfFOCUSKMatrix,KMatrix>
describeHerwigKPiIHalfFOCUSKMatrix("Herwig::KPiIHalfFOCUSKMatrix", "HwFormFactors.so");

void KPiIHalfFOCUSKMatrix::Init() {

  static ClassDocumentation<KPiIHalfFOCUSKMatrix> documentation
    ("The KPiIHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration (Phys.Lett. B653 (2007) 1-11) for the I=1/2 "
     "component of the Kpi K-matrix.",
     "The KPiIHalfFOCUSKMatrix class implements the K-matrix fit of "
     "the FOCUS collaboration \\cite{Pennington:2007se} for the $I=1/2$ "
     "component of the $K\\pi$ K-matrix.",
     "\\bibitem{Pennington:2007se}"
     "J.~M.~Link {\\it et al.} [FOCUS Collaboration],"
     "%``Dalitz plot analysis of the $D^{+} \\to K^{-} \\pi^{+} \\pi^{+}$ decay in the FOCUS experiment,''"
     "Phys.\\ Lett.\\ B {\\bf 653} (2007) 1"
     "doi:10.1016/j.physletb.2007.06.070"
     "[arXiv:0705.2248 [hep-ex]]."
     "%%CITATION = doi:10.1016/j.physletb.2007.06.070;%%"
     "%79 citations counted in INSPIRE as of 14 Jan 2020");

}

void KPiIHalfFOCUSKMatrix::doinit() {
  KMatrix::doinit();
  Energy mK = getParticleData(ParticleID::Kplus)->mass();
  Energy mpi= getParticleData(ParticleID::piplus)->mass();
  sNorm_ = sqr(mK)+sqr(mpi);
}

boost::numeric::ublas::matrix<double> KPiIHalfFOCUSKMatrix::K(Energy2 s, bool multiplyByPoles) const {
  double st = s/sNorm_-1.;
  double pre = (s-sHalf_)/sNorm_;
  Energy2 denom = !multiplyByPoles ? poles()[0]-s : poles()[0];
  boost::numeric::ublas::matrix<double> output =
    boost::numeric::ublas::zero_matrix<double>(2,2);
  output(0,0) = poleCouplings()[0][0]*poleCouplings()[0][0]/denom;
  output(0,1) = poleCouplings()[0][0]*poleCouplings()[0][1]/denom;
  output(1,1) = poleCouplings()[0][1]*poleCouplings()[0][1]/denom;
  double param = !multiplyByPoles ? 1. : (1.-s/poles()[0]);
  for(unsigned int ix=0;ix<C11_.size();++ix) {
    output(0,0) += C11_[ix]*param;
    output(1,1) += C22_[ix]*param;
    output(0,1) += C12_[ix]*param;
    param *= st;
  }
  output(1,0) = output(0,1);
  output *= pre;
  return output;
}
#line 1 "./PiPiAnisovichKMatrix.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PiPiAnisovichKMatrix class.
//

#include "PiPiAnisovichKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

PiPiAnisovichKMatrix::PiPiAnisovichKMatrix()
  : KMatrix(FlavourInfo(IsoSpin::IZero,IsoSpin::I3Unknown,
			Strangeness::Zero,Charm::Zero,
			Beauty::Zero),
	    vector<Channels>({KMatrix::PiPi,KMatrix::KK,KMatrix::FourPi,KMatrix::EtaEta,KMatrix::EtaEtaPrime}),
	    vector<Energy2>({sqr(0.65100*GeV), sqr(1.20720*GeV), sqr(1.56122*GeV), sqr(1.21257*GeV), sqr(1.81746*GeV) }),
	    vector<vector<Energy> >({{0.24844*GeV,-0.52523*GeV,0.00000*GeV,-0.38878*GeV,-0.36397*GeV},
				     {0.91779*GeV,0.55427*GeV,0.00000*GeV,0.38705*GeV,0.29448*GeV},
				     {0.37024*GeV,0.23591*GeV,0.62605*GeV,0.18409*GeV,0.18923*GeV},
				     {0.34501*GeV,0.39642*GeV,0.97644*GeV,0.19746*GeV,0.00357*GeV},
				     {0.15770*GeV,-0.17915*GeV,-0.90100*GeV,-0.00931*GeV,0.20689*GeV}})),
    s0Scatt_(-3.30564*GeV2), f1a_({0.26681,0.16583,-0.19840,0.32808,0.31193}),sA_(1.),sA0_(-0.2)
{}

IBPtr PiPiAnisovichKMatrix::clone() const {
  return new_ptr(*this);
}

IBPtr PiPiAnisovichKMatrix::fullclone() const {
  return new_ptr(*this);
}

void PiPiAnisovichKMatrix::persistentOutput(PersistentOStream & os) const {
  os << ounit(s0Scatt_,GeV2) << f1a_ << sA_ << sA0_ << ounit(mPi_,GeV);
}

void PiPiAnisovichKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> iunit(s0Scatt_,GeV2) >> f1a_ >> sA_ >> sA0_ >> iunit(mPi_,GeV);
}

void PiPiAnisovichKMatrix::doinit() {
  KMatrix::doinit();
  mPi_ = getParticleData(211)->mass();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PiPiAnisovichKMatrix,KMatrix>
describeHerwigPiPiAnisovichKMatrix("Herwig::PiPiAnisovichKMatrix", "HwFormFactors.so");

void PiPiAnisovichKMatrix::Init() {

  static ClassDocumentation<PiPiAnisovichKMatrix> documentation
    ("The PiPiAnisovichKMatrix class implements the K matrix pararmeterization from Eur.Phys.J.A 16 (2003) 229-258",
     "The K-matrix parameterization of Anisovich and Sarantsev \\cite{Anisovich:2002ij} was used.",
     "\bibitem{Anisovich:2002ij}\n"
     "V.~V.~Anisovich and A.~V.~Sarantsev,"
     "%``K matrix analysis of the (I J**(PC) = 00++)-wave in the mass region below 1900 MeV,''\n"
     "Eur. Phys. J. A \\textbf{16} (2003), 229-258\n"
     "doi:10.1140/epja/i2002-10068-x\n"
     "[arXiv:hep-ph/0204328 [hep-ph]].\n"
     "%170 citations counted in INSPIRE as of 20 Jun 2022\n");
  
  static Parameter<PiPiAnisovichKMatrix,Energy2> interfaces0Scatt
    ("s0Scatt",
     "The s0 parameters for scattering",
     &PiPiAnisovichKMatrix::s0Scatt_, GeV2, -3.30564*GeV2, -100.0*GeV2, 100.0*GeV2,
     false, false, Interface::limited);

  static ParVector<PiPiAnisovichKMatrix,double> interfaceF1A
    ("F1A",
     "The non-pole coefficients",
     &PiPiAnisovichKMatrix::f1a_, 5, 1.0, -100, 100.,
     false, false, Interface::limited);

  static Parameter<PiPiAnisovichKMatrix,double> interfacesA
    ("sA",
     "The sA parameter",
     &PiPiAnisovichKMatrix::sA_, 1., -10., 10.0,
     false, false, Interface::limited);

  static Parameter<PiPiAnisovichKMatrix,double> interfacesA0
    ("sA0",
     "The sA0 parameter",
     &PiPiAnisovichKMatrix::sA0_, -0.2, -10., 10.0,
     false, false, Interface::limited);
}

boost::numeric::ublas::matrix<double> PiPiAnisovichKMatrix::K(Energy2 s, bool multiplyByPoles) const {
  double pre = (s-0.5*sA_*sqr(mPi_))*(1.-sA0_)/(s-sA0_*GeV2);
  double coeff = (GeV2-s0Scatt_)/(s-s0Scatt_);
  if(multiplyByPoles)
    for (Energy2 pole : poles() ) coeff *= 1.-s/pole;
  boost::numeric::ublas::matrix<double> output =
    boost::numeric::ublas::zero_matrix<double>(5,5);
  for(unsigned int im=0;im<poles().size();++im) {
    InvEnergy2 term;
    if(multiplyByPoles) {
      term = 1./poles()[im];
      for(unsigned int iz=0;iz<poles().size();++iz) {
	if(iz==im) continue;
	term *= 1. - s/poles()[iz];
      }
    }
    else
      term = 1./(poles()[im]-s);
    for(unsigned int ix=0;ix<5;++ix)
      for(unsigned int iy=ix;iy<5;++iy)
	output(ix,iy) += term*poleCouplings()[im][ix]*poleCouplings()[im][iy];
  }
  for(unsigned int iy=0;iy<5;++iy) output(0,iy) += coeff*f1a_[iy];
  for(unsigned int ix=0;ix<5;++ix)
    for(unsigned int iy=ix+1;iy<5;++iy)
      output(iy,ix)=output(ix,iy);
  output *= pre;
  return output;
}
#line 1 "./AnalyticOmnesFunction.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnalyticOmnesFunction class.
//

#include "AnalyticOmnesFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"

using namespace Herwig;

IBPtr AnalyticOmnesFunction::clone() const {
  return new_ptr(*this);
}

IBPtr AnalyticOmnesFunction::fullclone() const {
  return new_ptr(*this);
}

void AnalyticOmnesFunction::doinit() {
  OmnesFunction::doinit();
  // set the parameters
  tPDPtr rho(getParticleData(ParticleID::rho0));
  if(!localParameters_) {
    mRho_=rho->mass();
    rhoWidth_=rho->width();
  }
  mPi_=getParticleData(ParticleID::piplus)->mass();
  Energy pcm(Kinematics::pstarTwoBodyDecay(mRho_,mPi_,mPi_));
  rhoConst_=sqr(mRho_)*rhoWidth_/pow<3,1>(pcm);
}

void AnalyticOmnesFunction::persistentOutput(PersistentOStream & os) const {
  os << ounit(fPi_,MeV) << ounit(mRho_,MeV) << ounit(rhoWidth_,MeV)
     << rhoConst_ << ounit(mPi_,MeV) << localParameters_;
}
void AnalyticOmnesFunction::persistentInput(PersistentIStream & is, int) {
  is >> iunit(fPi_,MeV) >> iunit(mRho_,MeV) >> iunit(rhoWidth_,MeV)
     >> rhoConst_ >> iunit(mPi_,MeV) >> localParameters_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<AnalyticOmnesFunction,OmnesFunction>
describeHerwigAnalyticOmnesFunction("Herwig::AnalyticOmnesFunction", "HwFormFactors.so");

void AnalyticOmnesFunction::Init() {

  static ClassDocumentation<AnalyticOmnesFunction> documentation
    ("The AnalyticOmnesFunction class implements the analytic version of the Omnes function from hep-ph/0112150",
     "The AnalyticOmnesFunction class implementing the analytic version of the Omnes function "
     "from \\cite{Holstein:2001bt} was used.",
     "\bibitem{Holstein:2001bt}\n"
     "B.~R.~Holstein,\n"
     "%``Allowed eta decay modes and chiral symmetry,''\n"
     "Phys. Scripta T \textbf{99} (2002), 55-67\n"
     "doi:10.1238/Physica.Topical.099a00055\n"
     "[arXiv:hep-ph/0112150 [hep-ph]].\n");

  static Parameter<AnalyticOmnesFunction,Energy> interfacefpi
    ("fpi",
     "The pion decay constant",
     &AnalyticOmnesFunction::fPi_, MeV, 130.7*MeV, ZERO, 200.*MeV,
     false, false, false); 

  static Parameter<AnalyticOmnesFunction,Energy> interfaceRhoMass
    ("RhoMass",
     "The mass of the rho",
     &AnalyticOmnesFunction::mRho_, MeV, 771.1*MeV, 400.*MeV, 1000.*MeV,
     false, false, false);

  static Parameter<AnalyticOmnesFunction,Energy> interfaceRhoWidth
    ("RhoWidth",
     "The width of the rho",
     &AnalyticOmnesFunction::rhoWidth_, MeV, 149.2*MeV, 100.*MeV, 300.*MeV,
     false, false, false);

  static Switch<AnalyticOmnesFunction,bool> interfaceLocalParameters
    ("LocalParameters",
     "Use local values of the rho mass and width",
     &AnalyticOmnesFunction::localParameters_, true, false, false);
  static SwitchOption interfaceLocalParametersLocal
    (interfaceLocalParameters,
     "Local",
     "Use local parameters",
     true);
  static SwitchOption interfaceLocalParametersParticleData
    (interfaceLocalParameters,
     "ParticleData",
     "Use values from the particle data objects",
     false);
}

Complex AnalyticOmnesFunction::D(Energy2 s) const {
    Energy2 mpi2(mPi_*mPi_),mrho2(mRho_*mRho_);
    double root, pi2 = sqr(Constants::pi);
    Complex f,ii(0.,1.);
    double pre(mpi2/12./pi2/fPi_/fPi_);
    if(s>4.*mpi2) {
      // real piece
      root=sqrt(1.-4.*mpi2/s);
      f=(1.-0.25*s/mpi2)*root*log((root+1.)/(-root+1.))-2.;
      f *=pre;
      // imaginary piece
      f += ii*s/mrho2*rhoConst_/8.*pow(root,3);
    }
    else {
      root=sqrt(4.*mpi2/s-1.);
      f=2.*(1.-0.25*s/mpi2)*root*atan2(1.,root)-2.;
      f *=pre;
    }
    return 1.-s/mrho2-s/48./pi2/fPi_/fPi_*log(mrho2/mpi2)-f;
}

void AnalyticOmnesFunction::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::AnalyticOmnesFunction " << name() << " \n";
  output << "newdef " << name() << ":fpi             " << fPi_/MeV         << "\n";
  output << "newdef " << name() << ":RhoMass         " << mRho_/MeV        << "\n";
  output << "newdef " << name() << ":RhoWidth        " << rhoWidth_/MeV    << "\n";
  output << "newdef " << name() << ":LocalParameters " << localParameters_ << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./ExperimentalOmnesFunction.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExperimentalOmnesFunction class.
//

#include "ExperimentalOmnesFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

ExperimentalOmnesFunction::ExperimentalOmnesFunction() {
  // initialization of the experimental function
  initialize_ =false;
  nPoints_=100;
  energy_ = {300*MeV, 320*MeV, 340*MeV, 360*MeV, 380*MeV,
	     400*MeV, 420*MeV, 440*MeV, 460*MeV, 480*MeV,
	     500*MeV, 520*MeV, 540*MeV, 560*MeV, 580*MeV,
	     600*MeV, 620*MeV, 640*MeV, 660*MeV, 680*MeV,
	     700*MeV, 720*MeV, 740*MeV, 760*MeV, 780*MeV,
	     800*MeV, 820*MeV, 840*MeV, 860*MeV, 880*MeV,
	     900*MeV, 920*MeV, 940*MeV, 960*MeV, 980*MeV};
  phase_={  00.1,     0.4,     0.7,     1.0,     1.5,
	    02.0,     2.5,     3.2,     4.0,     4.9,
	    05.9,     7.1,     8.5,    10.1,    12.1,
	    14.4,    17.3,    20.9,    25.4,    31.2,
	    38.7,    48.4,    60.6,    74.9,    90.0, 
	    103.8,   115.3,   124.3,   131.3,   136.7,
	    141.0,   144.5,   147.3,   149.7,   151.8};
  // experimental omnes function 
  omnesEnergy_ ={ 282.534*MeV, 289.32 *MeV, 296.106*MeV, 302.893*MeV, 309.679*MeV,
		  316.466*MeV, 323.252*MeV, 330.038*MeV, 336.825*MeV, 343.611*MeV,
		  350.398*MeV, 357.184*MeV, 363.97 *MeV, 370.757*MeV, 377.543*MeV,
		  384.33 *MeV, 391.116*MeV, 397.902*MeV, 404.689*MeV, 411.475*MeV,
		  418.261*MeV, 425.048*MeV, 431.834*MeV, 438.621*MeV, 445.407*MeV, 
		  452.193*MeV, 458.98 *MeV, 465.766*MeV, 472.553*MeV, 479.339*MeV,
		  486.125*MeV, 492.912*MeV, 499.698*MeV, 506.485*MeV, 513.271*MeV,
		  520.057*MeV, 526.844*MeV, 533.63 *MeV, 540.417*MeV, 547.203*MeV,
		  553.989*MeV, 560.776*MeV, 567.562*MeV, 574.349*MeV, 581.135*MeV,
		  587.921*MeV, 594.708*MeV, 601.494*MeV, 608.281*MeV, 615.067*MeV,
		  621.853*MeV, 628.64 *MeV, 635.426*MeV, 642.213*MeV, 648.999*MeV,
		  655.785*MeV, 662.572*MeV, 669.358*MeV, 676.145*MeV, 682.931*MeV,
		  689.717*MeV, 696.504*MeV, 703.29 *MeV, 710.077*MeV, 716.863*MeV,
		  723.649*MeV, 730.436*MeV, 737.222*MeV, 744.009*MeV, 750.795*MeV,
		  757.581*MeV, 764.368*MeV, 771.154*MeV, 777.94 *MeV, 784.727*MeV, 
		  791.513*MeV, 798.3  *MeV, 805.086*MeV, 811.872*MeV, 818.659*MeV, 
		  825.445*MeV, 832.232*MeV, 839.018*MeV, 845.804*MeV, 852.591*MeV, 
		  859.377*MeV, 866.164*MeV, 872.95 *MeV, 879.736*MeV, 886.523*MeV,
		  893.309*MeV, 900.096*MeV, 906.882*MeV, 913.668*MeV, 920.455*MeV, 
		  927.241*MeV, 934.028*MeV, 940.814*MeV, 947.6  *MeV, 954.387*MeV};
  omnesFunctionRe_ = { 0.860676  , 0.851786  , 0.843688  , 0.835827  , 0.828031  ,
		       0.820229  , 0.812370  , 0.804424  , 0.796354  , 0.788143  ,
		       0.779698  , 0.770939  , 0.761692  , 0.752707  , 0.743823  ,
		       0.735004  , 0.726091  , 0.717047  , 0.707862  , 0.698439  ,
		       0.688685  , 0.678510  , 0.668518  , 0.658481  , 0.648344  ,
		       0.638219  , 0.627989  , 0.617603  , 0.607222  , 0.596711  ,
		       0.586026  , 0.575280  , 0.564282  , 0.553067  , 0.541923  ,
		       0.530574  , 0.519112  , 0.507690  , 0.496055  , 0.484313  ,
		       0.472496  , 0.460245  , 0.447943  , 0.435766  , 0.423390  ,
		       0.410997  , 0.398510  , 0.385479  , 0.372458  , 0.359520  ,
		       0.346129  , 0.332837  , 0.319623  , 0.305858  , 0.292238  ,
		       0.278690  , 0.264391  , 0.250316  , 0.236400  , 0.221655  ,
		       0.207196  , 0.192956  , 0.177745  , 0.162833  , 0.148209  ,
		       0.132603  , 0.117202  , 0.102090  , 0.0862283 , 0.0703392 ,     
		       0.0545317 , 0.0383762 , 0.0219486 , 0.00518648,-0.0113217 ,    
		       -0.0280201 ,-0.045445  ,-0.0625479 ,-0.079748  ,-0.0978819 ,    
		       -0.11569   ,-0.133447  ,-0.152117  ,-0.170608  ,-0.189137  ,     
		       -0.208597  ,-0.227864  ,-0.247185  ,-0.267306  ,-0.287382  ,
		       -0.307707  ,-0.328882  ,-0.350103  ,-0.37178   ,-0.394464  ,
		       -0.417228  ,-0.440561  ,-0.464976  ,-0.490278  ,-0.517527};
  omnesFunctionIm_ = { 0.00243346, 0.000894972,-0.000612496,-0.00209178,-0.00354344,
		       -0.00496737,-0.00636316 ,-0.00773022 ,-0.00906769,-0.0103569 ,
		       -0.0116108 ,-0.0128658  ,-0.0145424  ,-0.0165746 ,-0.0186438 ,
		       -0.0206363 ,-0.0225379  ,-0.0243827  ,-0.0261488 ,-0.0278572 ,
		       -0.0295317 ,-0.0316349  ,-0.0339321  ,-0.0362345 ,-0.0386555 ,
		       -0.0410799 ,-0.0434534  ,-0.0459509  ,-0.0484302 ,-0.0508376 ,
		       -0.0533398 ,-0.0557937  ,-0.0581587  ,-0.0608612 ,-0.0635382 ,
		       -0.0661231 ,-0.068983   ,-0.0717604  ,-0.0744215 ,-0.0772635 ,
		       -0.0799845 ,-0.0825991  ,-0.0857537  ,-0.0888139 ,-0.0917441 ,
		       -0.0948263 ,-0.0977055  ,-0.100462   ,-0.103773  ,-0.106912  ,
		       -0.109931  ,-0.113413   ,-0.116647   ,-0.119722  ,-0.123282  ,
		       -0.126521  ,-0.129593   ,-0.133324   ,-0.136691  ,-0.139854  ,
		       -0.143729  ,-0.14718    ,-0.150356   ,-0.154353  ,-0.157926  ,
		       -0.161133  ,-0.165174   ,-0.168899   ,-0.172212  ,-0.176116  ,
		       -0.179892  ,-0.183445   ,-0.187134   ,-0.190947  ,-0.195144  ,
		       -0.198771  ,-0.202443   ,-0.206906   ,-0.210561  ,-0.214207  ,
		       -0.218943  ,-0.222806   ,-0.226551   ,-0.231273  ,-0.235267  ,
		       -0.239178  ,-0.244082   ,-0.24836    ,-0.252492  ,-0.257394  ,
		       -0.261812  ,-0.266156   ,-0.271161   ,-0.275849  ,-0.280675  ,
		       -0.286275  ,-0.291716   ,-0.297353   ,-0.303621  ,-0.310452  };
  // integration cut parameter
  epsCut_=0.4*MeV;
  // size of the arrays
  nsizea_ = energy_.size(); nsizeb_ = omnesEnergy_.size();
}


IBPtr ExperimentalOmnesFunction::clone() const {
  return new_ptr(*this);
}

IBPtr ExperimentalOmnesFunction::fullclone() const {
  return new_ptr(*this);
}

void ExperimentalOmnesFunction::doinit() {
  using Constants::pi;
  OmnesFunction::doinit();
  if(initialize_) {
    // convert the phase shift into radians
    vector<double> radphase;
    for(unsigned int ix=0;ix<phase_.size();++ix) {
      radphase.push_back(phase_[ix]/180.*Constants::pi);
    }
    // set up an interpolator for this
    interpolator_ = make_InterpolatorPtr(radphase,energy_,3);
    // limits and step sizes
    Energy mPi = getParticleData(ParticleID::piplus)->mass();
    Energy mEta(getParticleData(ParticleID::etaprime)->mass());
    Energy moff(2.*mPi),upp(mEta),step((mEta-moff)/nPoints_);
    // intergrator
    GaussianIntegrator integrator;
    // integrand
    // loop for integrals
    double D1real,D1imag;
    Complex ii(0.,1.),answer;
    moff+=0.5*step;
    omnesFunctionRe_.clear();
    omnesFunctionIm_.clear();
    omnesEnergy_.clear();
    for( ;moff<upp;moff+=step) {
      s_ = sqr(moff);
      // piece between 0 and 1 GeV
      Energy2 moff2(sqr(moff)),eps2(sqr(epsCut_));
      D1real=-moff2*(integrator.value(*this,4.*sqr(mPi),moff2-eps2)+
		     integrator.value(*this,moff2+eps2,upp*upp))/pi;
      D1imag=-(*interpolator_)(moff);
      // piece above 1 GeV
      D1real+=-(*interpolator_)(upp)/pi*log(upp*upp/(upp*upp-moff*moff));
      // calculate the answer
      answer = exp(D1real+ii*D1imag);
      // put into the arrays
      omnesFunctionRe_.push_back(answer.real());
      omnesFunctionIm_.push_back(answer.imag());
      omnesEnergy_.push_back(moff);
    }
  }
}

void ExperimentalOmnesFunction::persistentOutput(PersistentOStream & os) const {
  os << ounit(energy_,MeV) << ounit(omnesEnergy_,MeV) 
     << phase_ << omnesFunctionRe_ << omnesFunctionIm_ << initialize_
     << nPoints_ << ounit(epsCut_,MeV);
}

void ExperimentalOmnesFunction::persistentInput(PersistentIStream & is, int) {
  is >> iunit(energy_,MeV) >> iunit(omnesEnergy_,MeV) 
     >> phase_ >> omnesFunctionRe_ >> omnesFunctionIm_ >> initialize_
     >> nPoints_ >> iunit(epsCut_,MeV);
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<ExperimentalOmnesFunction,OmnesFunction>
describeHerwigExperimentalOmnesFunction("Herwig::ExperimentalOmnesFunction", "HwFormFactors.so");
HERWIG_INTERPOLATOR_CLASSDESC(ExperimentalOmnesFunction,double,Energy)

void ExperimentalOmnesFunction::Init() {

  static ClassDocumentation<ExperimentalOmnesFunction> documentation
    ("The ExperimentalOmnesFunction class implements the Omnes function"
     " integrating the experimental measurement of the phase");

  static Switch<ExperimentalOmnesFunction,bool> interfaceInitializeOmnes
    ("Initialize",
     "Initialize the experimental version of the Omnes function.",
     &ExperimentalOmnesFunction::initialize_, false, false, false);
  static SwitchOption interfaceInitializeOmnesInitialize
    (interfaceInitializeOmnes,
     "Yes",
     "Perform the initialization",
     true);
  static SwitchOption interfaceInitializeOmnesNoInitialization
    (interfaceInitializeOmnes,
     "No",
     "No initialization",
     false);

  static ParVector<ExperimentalOmnesFunction,Energy> interfacePhase_Energy
    ("Phase_Energy",
     "The energy values for the phase shift for the experimental version of the"
     " Omnes function",
     &ExperimentalOmnesFunction::energy_, MeV, -1, 1.0*MeV, 300.0*MeV, 2000.0*MeV,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,double> interfacePhase_Shift
    ("Phase_Shift",
     "The experimental values of the phase shift for the experimental version"
     " of the Omnes function",
     &ExperimentalOmnesFunction::phase_, 1.0, -1, 0.0, 0.0, 1000.0,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,Energy> interfaceOmnesEnergy
    ("OmnesEnergy",
     "The energy values for the interpolation of the experimental Omnes function",
     &ExperimentalOmnesFunction::omnesEnergy_, MeV, -1, 1.*MeV, 250.0*MeV, 2000.*MeV,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,double> interfaceOmnesReal
    ("OmnesReal",
     "The real part of the experimental Omnes function for the interpolation.",
     &ExperimentalOmnesFunction::omnesFunctionRe_, 1., -1, 1., -100.,
     100.,
     false, false, true);

  static ParVector<ExperimentalOmnesFunction,double> interfaceOmnesImag
    ("OmnesImag",
     "The imaginary part of the experimental Omnes function for the interpolation.",
     &ExperimentalOmnesFunction::omnesFunctionIm_, 1., -1, 1., -100.,
     100.,
     false, false, true);

  static Parameter<ExperimentalOmnesFunction,unsigned int> interfaceOmnesPoints
    ("OmnesPoints",
     "The number of points for the interpolation table for the experimental"
     " Omnes function.",
     &ExperimentalOmnesFunction::nPoints_, 100, 50, 200,
     false, false, true);

  static Parameter<ExperimentalOmnesFunction,Energy> interfaceOmnesCut
    ("OmnesCut",
     "The cut parameter for the integral in the experimental Omnes function.",
     &ExperimentalOmnesFunction::epsCut_, MeV, 0.1*MeV, 0.001*MeV, 1.0*MeV,
     false, false, true);

}

Complex ExperimentalOmnesFunction::D(Energy2 s) const {
  if(!oRe_) {
    oRe_ = make_InterpolatorPtr(omnesFunctionRe_,omnesEnergy_,3);
    oIm_ = make_InterpolatorPtr(omnesFunctionIm_,omnesEnergy_,3);
  }
  Energy q(sqrt(s)); Complex ii(0.,1.);
  return (*oRe_)(q)+ii*(*oIm_)(q);
}

void ExperimentalOmnesFunction::dataBaseOutput(ofstream & output,bool header,
					   bool create) const {
  if(header) output << "update decayers set parameters=\"";
  if(create) output << "create Herwig::ExperimentalOmnesFunction " << name() << " \n";
  output << "newdef " << name() << ":Initialize " << initialize_      << "\n";
  output << "newdef " << name() << ":OmnesPoints     " << nPoints_         << "\n";
  output << "newdef " << name() << ":OmnesCut        " << epsCut_/MeV      << "\n";
  for(unsigned int ix=0;ix<energy_.size();++ix) {
    if(ix<nsizea_) {
      output << "newdef " << name() << ":Phase_Energy " << ix << "  " 
  	     << energy_[ix]/MeV << "\n";
      output << "newdef " << name() << ":Phase_Shift  " << ix << "  " 
  	     << phase_[ix]  << "\n";
    }
    else {
      output << "insert " << name() << ":Phase_Energy " << ix << "  " 
  	     << energy_[ix]/MeV << "\n";
      output << "insert " << name() << ":Phase_Shift  " << ix << "  " 
  	     << phase_[ix]  << "\n";
    }
  }
  for(unsigned int ix=0;ix<omnesEnergy_.size();++ix) {
    if(ix<nsizeb_) {
      output << "newdef " << name() << ":OmnesEnergy " << ix << "  " 
	     << omnesEnergy_[ix]/MeV << "\n";
      output << "newdef " << name() << ":OmnesReal " << ix << "  " 
	     << omnesFunctionRe_[ix] << "\n";
      output << "newdef " << name() << ":OmnesImag " << ix << "  " 
	     << omnesFunctionIm_[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":OmnesEnergy " << ix << "  " 
	     << omnesEnergy_[ix]/MeV << "\n";
      output << "insert " << name() << ":OmnesReal " << ix << "  " 
	     << omnesFunctionRe_[ix] << "\n";
      output << "insert " << name() << ":OmnesImag " << ix << "  " 
	     << omnesFunctionIm_[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
