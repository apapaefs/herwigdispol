#line 1 "./LeptoquarkModel.cc"

// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModel class.
//

#include "LeptoquarkModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void LeptoquarkModel::doinit()  {
  addVertex(_theSLQSLQGVertex);
  addVertex(_theSLQSLQGGVertex);
  addVertex(_theSLQFFVertex);
  
  BSMModel::doinit();
}

LeptoquarkModel::LeptoquarkModel() :  _CouplFF(0.312), 
				      _leftcoup(1.0), 
				      _rightcoup(1.0), 
				      _rightcouptilde(1.0), 
				      _leftcoup1(1.0) , 
				      _leftcoup12(1.0), 
				      _rightcoup12(1.0), 
				      _leftcoup12t(1.0), 
				      _dleftcoup(1.0), 
				      _drightcoup(1.0), 
				      _drightcouptilde(1.0), 
				      _dleftcoup1(1.0) , 
				      _dleftcoup12(1.0), 
				      _drightcoup12(1.0), 
				      _dleftcoup12t(1.0), 
				      _derivscalef(500.0*GeV) 
{}


IBPtr LeptoquarkModel::clone() const {
  return new_ptr(*this);
}
IBPtr LeptoquarkModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void LeptoquarkModel::persistentOutput(PersistentOStream & os) const {
  os <<  _theSLQSLQGGVertex
     << _theSLQSLQGVertex
     << _theSLQFFVertex
     << _CouplFF
     << _leftcoup
     << _rightcoup
     << _leftcoup1
     << _rightcouptilde
     << _leftcoup12
     << _rightcoup12
     << _leftcoup12t
     << _dleftcoup
     << _drightcoup
     << _dleftcoup1
     << _drightcouptilde
     << _dleftcoup12
     << _drightcoup12
     << _dleftcoup12t
     << ounit(_derivscalef,GeV);

    
  
}

void LeptoquarkModel::persistentInput(PersistentIStream & is, int) {
  is >> _theSLQSLQGGVertex
     >> _theSLQSLQGVertex
     >> _theSLQFFVertex
     >> _CouplFF
     >> _leftcoup
     >> _rightcoup
     >> _leftcoup1
     >> _rightcouptilde
     >> _leftcoup12
     >> _rightcoup12
     >> _leftcoup12t
     >> _dleftcoup
     >> _drightcoup
     >> _dleftcoup1
     >> _drightcouptilde
     >> _dleftcoup12
     >> _drightcoup12
     >> _dleftcoup12t
     >> iunit(_derivscalef,GeV);
    
  
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LeptoquarkModel,BSMModel>
describeHerwigLeptoquarkModel("Herwig::LeptoquarkModel", "HwLeptoquarkModel.so");

void LeptoquarkModel::Init() {
  
  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractVSSVertex> interfaceVertexSLQSLQG
  ("Vertex/SLQSLQG",
   "Reference to the scalar leptoquark-scalar leptoquark-gluon vertex",
   &LeptoquarkModel::_theSLQSLQGVertex, false, false, true, false, false);

  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractVVSSVertex> interfaceVertexSLQSLQGG
  ("Vertex/SLQSLQGG",
   "Reference to the scalar leptoquark-scalar leptoquark-gluon-gluon vertex",
   &LeptoquarkModel::_theSLQSLQGGVertex, false, false, true, false, false);

  static Reference<LeptoquarkModel,ThePEG::Helicity::AbstractFFSVertex> interfaceVertexSLQFF
  ("Vertex/SLQFF",
   "Reference to the scalar leptoquark-scalar-quark-lepton",
   &LeptoquarkModel::_theSLQFFVertex, false, false, true, false, false);

  static Parameter<LeptoquarkModel, double> interfaceLQCoupling
    ("LQCoupling",
     "The overall Leptoquark Coupling",
     &LeptoquarkModel::_CouplFF, 0.312, 0., 10.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_L
    ("g_S0_L",
     "The leptoquark S0 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_leftcoup, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_R
    ("g_S0_R",
     "The leptoquark S0 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_rightcoup, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_Rt
    ("g_S0t_R",
     "The leptoquark ~S0 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_rightcouptilde, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ_L1
    ("g_S1_L",
     "The leptoquark S1 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_leftcoup1, 1.0, 0., 1.0,
     false, false, Interface::limited);
  
    static Parameter<LeptoquarkModel, double> interfacegLQ12_L
    ("g_S12_L",
     "The leptoquark S1/2 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_leftcoup12, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegLQ12_R
    ("g_S12_R",
     "The leptoquark S1/2 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_rightcoup12, 1.0, 0., 1.0,
     false, false, Interface::limited);
  
  static Parameter<LeptoquarkModel, double> interfacegLQ12t_L
    ("g_S12t_L",
     "The leptoquark ~S1/2 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_leftcoup12t, 1.0, 0., 1.0,
     false, false, Interface::limited);


  static Parameter<LeptoquarkModel, double> interfacegdLQ_L
    ("g_dS0_L",
     "The leptoquark dS0 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_dleftcoup, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegdLQ_R
    ("g_dS0_R",
     "The leptoquark dS0 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_drightcoup, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegdLQ_Rt
    ("g_dS0t_R",
     "The leptoquark ~dS0 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_drightcouptilde, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegdLQ_L1
    ("g_dS1_L",
     "The leptoquark dS1 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_dleftcoup1, 1.0, 0., 1.0,
     false, false, Interface::limited);
  
    static Parameter<LeptoquarkModel, double> interfacegdLQ12_L
    ("g_dS12_L",
     "The leptoquark dS1/2 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_dleftcoup12, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, double> interfacegdLQ12_R
    ("g_dS12_R",
     "The leptoquark dS1/2 coupling LQ-lepton_right-quark_left",
     &LeptoquarkModel::_drightcoup12, 1.0, 0., 1.0,
     false, false, Interface::limited);
  
  static Parameter<LeptoquarkModel, double> interfacegdLQ12t_L
    ("g_dS12t_L",
     "The leptoquark ~dS1/2 coupling LQ-lepton_left-quark_right",
     &LeptoquarkModel::_dleftcoup12t, 1.0, 0., 1.0,
     false, false, Interface::limited);

  static Parameter<LeptoquarkModel, Energy> interfaceDerivativeScale
    ("derivscale",
     "The suppression scale for the derivatively coupled leptoquarks",
     &LeptoquarkModel::_derivscalef, GeV, 500.0*GeV, ZERO, 10000.0*GeV,
     false, false, Interface::limited);


  static ClassDocumentation<LeptoquarkModel> documentation
    ("There is no documentation for the LeptoquarkModel class");

}

#line 1 "./LeptoquarkModelSLQSLQGGVertex.cc"
// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQSLQGGVertex class.
//

#include "LeptoquarkModelSLQSLQGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"


using namespace Herwig;
using namespace ThePEG;

LeptoquarkModelSLQSLQGGVertex::LeptoquarkModelSLQSLQGGVertex() : q2last_(ZERO),
								 couplast_(0.) {
  orderInGs(2);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TTFUNDS);
}

void LeptoquarkModelSLQSLQGGVertex::doinit() {
  addToList(21,21,9941551,-9941551);
  addToList(21,21,9911561,-9911561);
  addToList(21,21,9921551,-9921551);
  addToList(21,21,9931561,-9931561);
  addToList(21,21,9931551,-9931551);
  addToList(21,21,9931661,-9931661);
  addToList(21,21,9941561,-9941561);

  addToList(21,21,9951551,-9951551);
  addToList(21,21,9951651,-9951651);
  addToList(21,21,9961551,-9961551);

  addToList(21,21,9971561,-9971561);
  addToList(21,21,9981561,-9981561);
  addToList(21,21,9981551,-9981551);
  addToList(21,21,9981651,-9981651);

  addToList(21,21,9991551,-9991551);
  addToList(21,21,9991561,-9991561);
  addToList(21,21,9901561,-9901561);
  addToList(21,21,9901661,-9901661);
  VVSSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<LeptoquarkModelSLQSLQGGVertex,VVSSVertex>
describeHerwigLeptoquarkModelSLQSLQGGVertex("Herwig::LeptoquarkModelSLQSLQGGVertex", "HwLeptoquarkModel.so");

void LeptoquarkModelSLQSLQGGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGGVertex> documentation
    ("The LeptoquarkModelSLQSLQGGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr ,tcPDPtr, tcPDPtr ) {
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = sqr(strongCoupling(q2));  
    q2last_ = q2;
  }
  norm(couplast_);
}
#line 1 "./LeptoquarkModelSLQFFVertex.cc"
// -*- C++ -*-
//
// LeptoquarkModelSLQFFVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQFFVertex class.
//

#include "LeptoquarkModelSLQFFVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr LeptoquarkModelSLQFFVertex::clone() const {
  return new_ptr(*this);
}

IBPtr LeptoquarkModelSLQFFVertex::fullclone() const {
  return new_ptr(*this);
}

LeptoquarkModelSLQFFVertex::LeptoquarkModelSLQFFVertex() {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void LeptoquarkModelSLQFFVertex::doinit() {
  //S0
  addToList( 15, 6,-9911561);
  addToList(-15,-6, 9911561);
  
  addToList(-16,-5, 9911561);
  addToList( 16, 5,-9911561);

  //~S0
  addToList(-15,-5, 9921551);
  addToList( 15, 5,-9921551);

  //S1 triplet
  //S1p
  addToList(-15,-5, 9931551);
  addToList( 15, 5,-9931551);
  //S1z
  addToList(-15,-6, 9931561);
  addToList( 15, 6,-9931561);
  addToList(-16,-5, 9931561);
  addToList( 16, 5,-9931561);
  //S1m
  addToList(-16,-6, 9931661);
  addToList( 16, 6,-9931661);

  //S1/2 doublet
  addToList( 15,-6, 9941561);
  addToList(-15, 6,-9941561);
  
  addToList(-15, 5,-9941551);
  addToList(-16, 6,-9941551);
  addToList( 15,-5, 9941551);
  addToList( 16,-6, 9941551);


  //S1/2 tilde doublet
  addToList( 5,-16,-9951651);
  addToList(-5, 16, 9951651);

  addToList(-5, 15, 9951551);
  addToList( 5,-15,-9951551);


  //dS0
  addToList( 15,-5, 9961551);
  addToList(-15, 5,-9961551);

  addToList( 16,-6, 9961551);
  addToList(-16, 6,-9961551);

  //~dS0
  addToList( 15,-6, 9971561);
  addToList(-15, 6,-9971561);


  //dS1 triplet

  //dS1p
  addToList( 15,-6, 9981561);
  addToList(-15, 6,-9981561);

  //dS1z
  addToList( 16,-6, 9981551);
  addToList(-16, 6,-9981551);

  addToList( 15,-5, 9981551);
  addToList(-15, 5,-9981551);

  //dS1m
  addToList( 16,-5, 9981651);
  addToList(-16, 5,-9981651);

  //dS1/2 doublet
  addToList(-15,-5, 9991551);
  addToList( 15, 5,-9991551);

  addToList(-15,-6, 9991561);
  addToList( 15, 6,-9991561);

  addToList(-16,-5, 9991561);
  addToList( 16, 5,-9991561);

  //dS1/2 tilde doublet
  addToList(-15,-6, 9901561);
  addToList( 15, 6,-9901561);

  addToList(-16,-6, 9901661);
  addToList( 16, 6,-9901661);


  _theModel = generator()->standardModel();
  tcHwLeptoquarkPtr hwLeptoquark=dynamic_ptr_cast<tcHwLeptoquarkPtr>(_theModel);
  if(hwLeptoquark){
    _CFF=hwLeptoquark->cfermion();
    _cL0 =hwLeptoquark->cleft();
    _cR0 =hwLeptoquark->cright();
    _cR0t = hwLeptoquark->crighttilde();
    _cL1 =hwLeptoquark->cleft1(); 
    _cL12 =hwLeptoquark->cleft12(); 
    _cR12 =hwLeptoquark->cright12(); 
    _cL12t =hwLeptoquark->cleft12tilde(); 
    
    
    _derivscale = hwLeptoquark->fscale();
    _dcL0 =hwLeptoquark->dcleft();
    _dcR0 =hwLeptoquark->dcright();
    _dcR0t = hwLeptoquark->dcrighttilde();
    _dcL1 =hwLeptoquark->dcleft1(); 
    _dcL12 =hwLeptoquark->dcleft12(); 
    _dcR12 =hwLeptoquark->dcright12(); 
    _dcL12t =hwLeptoquark->dcleft12tilde(); 
    
  }
  FFSVertex::doinit();
}

void LeptoquarkModelSLQFFVertex::persistentOutput(PersistentOStream & os) const {
  os << _CFF << _cL0 << _cR0 << _cR0t 
     << _cL1 << _cL12 << _cR12 << _cL12t 
     << _dcL0 << _dcR0 << _dcR0t 
     << _dcL1 << _dcL12 << _dcR12 << _dcL12t 
     << ounit(_derivscale,GeV);
}

void LeptoquarkModelSLQFFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _CFF >> _cL0 >> _cR0 >> _cR0t 
     >> _cL1 >> _cL12 >> _cR12 >> _cL12t 
     >>_dcL0 >> _dcR0 >> _dcR0t 
     >> _dcL1 >> _dcL12 >> _dcR12 >> _dcL12t 
     >> iunit(_derivscale,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<LeptoquarkModelSLQFFVertex,FFSVertex>
describeHerwigLeptoquarkModelSLQFFVertex("Herwig::LeptoquarkModelSLQFFVertex", "Herwig.so");


void LeptoquarkModelSLQFFVertex::Init() {
  
  static ClassDocumentation<LeptoquarkModelSLQFFVertex> documentation
    ("The LeptoquarkModelSLQFFVertex class is the implementation"
     " of the helicity amplitude calculation of the Leptoquark"
     " quark-lepton vertex.");
}


void LeptoquarkModelSLQFFVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  long isc(cc->id()), ism(aa->id()), 
    ichg(bb->id());
  long lqid = isc;
  long smid_1 = ism;
  long smid_2 = ichg;
  if(abs(lqid) < 9900000) { 
    lqid = ism; 
    smid_1 = ichg;
    smid_2 = isc; 
  }
  if(abs(lqid) < 9900000) {
    smid_1 = ism;
    smid_2 = isc;
  }
  if( abs(smid_1) > abs(smid_2) ) { swap(smid_1, smid_2); }

  const Energy denom = sqrt(2.) * _derivscale;
  const double mtop = getParticleData(ParticleID::t)->mass()        / denom;
  const double mbot = getParticleData(ParticleID::b)->mass()        / denom;
  const double mtau = getParticleData(ParticleID::tauminus)->mass() / denom;

  //set the couplings to left and right 
  //S0
  if( abs(isc) == 9911561 || abs(ism) == 9911561 || abs(ichg) == 9911561 ) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cL = -_cL0; _cR = Complex(0.);
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cL = _cL0;
      _cR = _cR0;
    }
  }
  //~S0
  if( abs(isc) == 9921551 || abs(ism) == 9921551 || abs(ichg) == 9921551 ) {
    _cL = Complex(0.); _cR = _cR0t;
  }
  
  //S1 triplet
  //Q = + 4/3
  if( abs(isc) == 9931551 || abs(ism) == 9931551 || abs(ichg) == 9931551 ) {
    _cL = sqrt(2.)* _cL1; _cR = Complex(0.);
  }
  //Q = + 1/3
  if( abs(isc) == 9931561 || abs(ism) == 9931561 || abs(ichg) == 9931561 ) {
    _cL = - _cL1; _cR = Complex(0.);
  }
  //Q = - 2/3
  if( abs(isc) == 9931661 || abs(ism) == 9931661 || abs(ichg) == 9931661 ) {
    _cL = sqrt(2.) * _cL1; _cR = Complex(0.);
  }
  
  //S1/2 doublet

  //Q = + 5/3 
  if( abs(isc) == 9941561 || abs(ism) == 9941561 || abs(ichg) == 9941561 ) {
    _cR = _cL12; _cL = _cR12;
  }
  
  
  //Q = + 2/3 
  if( abs(isc) == 9941551 || abs(ism) == 9941551 || abs(ichg) == 9941551 ) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cR = Complex(0.); _cL = - _cR12;
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cL = Complex(0.); _cR = _cL12;
    }
  }

  //S1/2 tilde doublet

  //Q = + 2/3 
  if( abs(isc) == 9951551 || abs(ism) == 9951551 || abs(ichg) == 9951551 ) {
    _cR = _cL12t; _cL = Complex(0.);
  }
  
  
  //Q = - 1/3 
  if( abs(isc) == 9951651 || abs(ism) == 9951651 || abs(ichg) == 9951651 ) {
    _cR = _cL12t; _cL = Complex(0.);
  }



  //dS0
  if( abs(isc) == 9961551 || abs(ism) == 9961551 || abs(ichg) == 9961551) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cR = _dcL0 * mbot +_dcR0 * mtau; 
      _cL = _dcR0 * mbot + _dcL0 * mtau;
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cR = _dcL0 * mtop; 
      _cL = Complex(0.);
    }
  }

  //d~S0
  if( abs(isc) == 9971561 || abs(ism) == 9971561 || abs(ichg) ==  9971561) {
    _cR = _dcR0t * mtau;
    _cL = _dcR0t * mtop;
  }

  //dS1 triplet
  if( abs(isc) == 9981561 || abs(ism) == 9981561 || abs(ichg) ==  9981561) {
    _cR = sqrt(2.) * _dcL1 * mtop;
    _cL = sqrt(2.) * _dcL1 * mtau;
  }
  if( abs(isc) == 9981551 || abs(ism) == 9981551 || abs(ichg) ==  9981551) {
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cR = -_dcL1 * mbot; 
      _cL = -_dcL1 * mtau;
    }
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cR = _dcL1 * mtop;
      _cL = Complex(0.);
    }
  }

  if( abs(isc) == 9981651 || abs(ism) == 9981651 || abs(ichg) ==  9981651) {
    _cL = sqrt(2.) * _dcL1 * mbot;
    _cR = Complex(0.);
  }
  
  
  //dS1/2 doublet
  if( abs(isc) == 9991551 || abs(ism) == 9991551 || abs(ichg) == 9991551 ) {
    _cL = _dcL12 * mbot + _dcR12 * mtau;
    _cR = _dcR12 * mbot + _dcL12 * mtau;
  }

  if( abs(isc) == 9991561 || abs(ism) == 9991561 || abs(ichg) == 9991561 ) {
    if(abs(isc) == 6 || abs(ism) == 6 || abs(ichg) == 6) { 
      _cL = _dcR12 * mtau; 
      _cR = _dcR12 * mtop;
    }
    if(abs(isc) == 5 || abs(ism) == 5 || abs(ichg) == 5) { 
      _cL = _dcL12 * mbot;
    }    
  }

  //dS1/2 tilde doublet
  if( abs(isc) == 9901561 || abs(ism) == 9901561  || abs(ichg) == 9901561 ) {
    _cL = _dcL12t * mtop; 
    _cR = _dcL12t * mtau;
  }
  
  if( abs(isc) == 9901661 || abs(ism) == 9901661  || abs(ichg) == 9901661 ) {
    _cL = _dcL12t * mtop;
    _cR = Complex(0.);
  }


  if(smid_1 > 0) { 
    left(conj(_cR)); 
    right(conj(_cL));
  } 
  else { 
    left(_cL); 
    right(_cR); 
  } 

  norm(_CFF);
}
#line 1 "./LeptoquarkModelSLQSLQGVertex.cc"
// -*- C++ -*-
//
// LeptoquarkModelSLQSLQGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LeptoquarkModelSLQSLQGVertex class.
//

#include "LeptoquarkModelSLQSLQGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

LeptoquarkModelSLQSLQGVertex::LeptoquarkModelSLQSLQGVertex() : q2last_(ZERO),
							       couplast_(0.) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}

void LeptoquarkModelSLQSLQGVertex::doinit() {
  addToList(21,9941551,-9941551);
  addToList(21,9911561,-9911561);
  addToList(21,9921551,-9921551);
  addToList(21,9931561,-9931561);
  addToList(21,9931551,-9931551);
  addToList(21,9931661,-9931661);
  addToList(21,9941561,-9941561);
  addToList(21,9951551,-9951551);
  addToList(21,9951651,-9951651);
  addToList(21,9961551,-9961551);
  addToList(21,9971561,-9971561);
  addToList(21,9981561,-9981561);
  addToList(21,9981551,-9981551);
  addToList(21,9981651,-9981651);
  addToList(21,9991551,-9991551);
  addToList(21,9991561,-9991561);
  addToList(21,9901561,-9901561);
  addToList(21,9901661,-9901661);
  VSSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<LeptoquarkModelSLQSLQGVertex,VSSVertex>
describeHerwigLeptoquarkModelSLQSLQGVertex("Herwig::LeptoquarkModelSLQSLQGVertex", "HwLeptoquarkModel.so");

void LeptoquarkModelSLQSLQGVertex::Init() {
  static ClassDocumentation<LeptoquarkModelSLQSLQGVertex> documentation
    ("The LeptoquarkModelSLQSLQGVertex class is the implementation of"
     " the LeptoquarkModel scalar LQ-scalar LQ-gluon vertex");
  
}

void LeptoquarkModelSLQSLQGVertex::setCoupling(Energy2 q2,tcPDPtr ,tcPDPtr p2,tcPDPtr ) { 
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = strongCoupling(q2);
    q2last_ = q2;
  }
  if(p2->id()<0)
    norm( couplast_);
  else
    norm(-couplast_);
}

