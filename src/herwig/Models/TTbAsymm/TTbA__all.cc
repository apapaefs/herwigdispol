#line 1 "./TTbAModel.cc"

// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModel class.
//

#include "TTbAModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;


void TTbAModel::doinit()  {
  addVertex(_theWPTDVertex);
  addVertex(_theZPQQVertex);
  addVertex(_theAGQQVertex);
  addVertex(_theSU2XVertex);
  StandardModel::doinit();
}

TTbAModel::TTbAModel(): _gWPTD_L(1.0), _gWPTD_R(1.0),_gZPTU_L(1.0), _gZPTU_R(1.0),_gZPUU_L(1.0), _gZPUU_R(1.0),_gZPCC_L(1.0), _gZPCC_R(1.0),_gAGQQ_L(1.0), _gAGQQ_R(1.0),_gAGTT_L(1.0), _gAGTT_R(1.0), _alphaXparam(0.060), _costhetaXparam(0.95), _modelselect(1) {}

IBPtr TTbAModel::clone() const {
  return new_ptr(*this);
}
IBPtr TTbAModel::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void TTbAModel::persistentOutput(PersistentOStream & os) const {
  os << _theWPTDVertex
     << _theZPQQVertex
     << _theAGQQVertex
     << _theSU2XVertex
     << _gWPTD_L
     << _gWPTD_R
     << _gZPTU_L
     << _gZPTU_R
     << _gZPUU_L
     << _gZPUU_R
     << _gZPCC_L
     << _gZPCC_R
     << _gAGQQ_L
     << _gAGQQ_R
     << _gAGTT_L
     << _gAGTT_R
     << _modelselect;
}

void TTbAModel::persistentInput(PersistentIStream & is, int) {
  is >> _theWPTDVertex
     >> _theZPQQVertex
     >> _theAGQQVertex
     >> _theSU2XVertex
     >> _gWPTD_L
     >> _gWPTD_R
     >> _gZPTU_L
     >> _gZPTU_R
     >> _gZPUU_L
     >> _gZPUU_R
     >> _gZPCC_L
     >> _gZPCC_R
     >> _gAGQQ_L
     >> _gAGQQ_R
     >> _gAGTT_L
     >> _gAGTT_R
     >> _modelselect;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TTbAModel,BSMModel>
describeHerwigTTbAModel("Herwig::TTbAModel", "HwTTbAModel.so");

void TTbAModel::Init() {
  
  static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexWPTD
  ("Vertex/WPTD",
   "Reference to the W prime Top Down vertex",
   &TTbAModel::_theWPTDVertex, false, false, true, false, false);

 static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexZPQQ
  ("Vertex/ZPQQ",
   "Reference to the Z prime Quark-Antiquark vertex",
   &TTbAModel::_theZPQQVertex, false, false, true, false, false);

 static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexAGQQ
  ("Vertex/AGQQ",
   "Reference to the Axigluon Quark-Antiquark vertex",
   &TTbAModel::_theAGQQVertex, false, false, true, false, false);

 static Reference<TTbAModel,ThePEG::Helicity::AbstractFFVVertex> interfaceVertexSU2X
  ("Vertex/SU2X",
   "Reference to the non-Abelian SU(2)_X vertex",
   &TTbAModel::_theSU2XVertex, false, false, true, false, false);




  static Parameter<TTbAModel, double> interfaceWPTDLCoupling
    ("WPTDLCoupling",
     "The left-handed W prime coupling to top down",
     &TTbAModel::_gWPTD_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceWPTDRCoupling
    ("WPTDRCoupling",
     "The right-handed W prime coupling to top down",
     &TTbAModel::_gWPTD_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceAGQQLCoupling
    ("AGQQLCoupling",
     "The left-handed axigluon coupling to q-qbar",
     &TTbAModel::_gAGQQ_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

static Parameter<TTbAModel, double> interfaceAGQQRCoupling
    ("AGQQRCoupling",
     "The right-handed axigluon coupling to q-qbar",
     &TTbAModel::_gAGQQ_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);


  static Parameter<TTbAModel, double> interfaceAGTTLCoupling
    ("AGTTLCoupling",
     "The left-handed axigluon coupling to t-tbar",
     &TTbAModel::_gAGTT_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

static Parameter<TTbAModel, double> interfaceAGTTRCoupling
    ("AGTTRCoupling",
     "The right-handed axigluon coupling to t-tbar",
     &TTbAModel::_gAGTT_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPTULCoupling
    ("ZPTULCoupling",
     "The left-handed Z prime coupling to top up",
     &TTbAModel::_gZPTU_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPTURCoupling
    ("ZPTURCoupling",
     "The right-handed Z prime coupling to top up",
     &TTbAModel::_gZPTU_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPUULCoupling
    ("ZPUULCoupling",
     "The left-handed Z prime coupling to up upbar",
     &TTbAModel::_gZPUU_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPUURCoupling
    ("ZPUURCoupling",
     "The right-handed Z prime coupling to up upbar",
     &TTbAModel::_gZPUU_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPCCLCoupling
    ("ZPCCLCoupling",
     "The left-handed Z prime coupling to char charmbar",
     &TTbAModel::_gZPCC_L, 1.0, -10.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceZPCCRCoupling
    ("ZPCCRCoupling",
     "The right-handed Z prime coupling to charm charmbar",
     &TTbAModel::_gZPCC_R, 1.0, -10.0, 10.0,
     false, false, Interface::limited);


  static Parameter<TTbAModel, double> interfaceSU2Xcostheta
    ("SU2Xcostheta",
     "Misalignment parameter of SU(2)_X model",
     &TTbAModel::_costhetaXparam, 0.95, -1.0, 1.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, double> interfaceSU2Xalpha
    ("SU2Xalpha",
     "alphaX coupling constant",
     &TTbAModel::_alphaXparam, 0.060, 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<TTbAModel, int> interfacemodelselect
    ("modelselect",
     "Selet which model to run",
     &TTbAModel::_modelselect, 0, 0, 4,
     false, false, Interface::limited);


  static ClassDocumentation<TTbAModel> documentation
    ("There is no documentation for the TTbAModel class");

}

#line 1 "./TTbAModelWPTDVertex.cc"
// -*- C++ -*-
//
// TTbAModelWPTDVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelWPTDVertex class.
//

#include "TTbAModelWPTDVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr TTbAModelWPTDVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelWPTDVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelWPTDVertex::TTbAModelWPTDVertex()  {
  addToList(-1,6,-34);
  addToList(-6,1,34);
  orderInGem(1);
  orderInGs(1);
  colourStructure(ColourStructure::DELTA);
}

void TTbAModelWPTDVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _cWPTD_R =hwTTbA->_cWPTD_right();
    _cWPTD_L =hwTTbA->_cWPTD_left();
    _models =hwTTbA->_model();
  }
  FFVVertex::doinit();
}

void TTbAModelWPTDVertex::persistentOutput(PersistentOStream & os) const {
  os << _cWPTD_R << _cWPTD_L << _models;
}

void TTbAModelWPTDVertex::persistentInput(PersistentIStream & is, int) {
  is >> _cWPTD_R >> _cWPTD_L >> _models; 
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TTbAModelWPTDVertex,FFVVertex>
describeHerwigTTbAModelWPTDVertex("Herwig::TTbAModelWPTDVertex", "Herwig.so");


void TTbAModelWPTDVertex::Init() {
  
  static ClassDocumentation<TTbAModelWPTDVertex> documentation
    ("The TTbAModelWPTDVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " quark-lepton vertex.");
}


void TTbAModelWPTDVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  
  double _cL = 0, _cR = 0;
  
  if(abs(aa->id()) == 34 || abs(bb->id()) == 34 || abs(cc->id()) == 34) {
    _cR = _cWPTD_R; 
    _cL = _cWPTD_L; 
  }
  if(_models!=0) { _cL = 1E-10; _cR = 1E-10; }
  left(_cL);
  right(_cR);
  
  norm(1.0);
}
#line 1 "./TTbAModelZPQQVertex.cc"
// -*- C++ -*-
//
// TTbAModelZPQQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelZPQQVertex class.
//

#include "TTbAModelZPQQVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;



IBPtr TTbAModelZPQQVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelZPQQVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelZPQQVertex::TTbAModelZPQQVertex()  {
  addToList(-2,6,32);
  addToList(-6,2,32);
  addToList(-2,2,32);
  addToList(-4,4,32);
  orderInGem(1);
  orderInGs(1);
  colourStructure(ColourStructure::DELTA);
}

void TTbAModelZPQQVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _cZPTU_R =hwTTbA->_cZPTU_right();
    _cZPTU_L =hwTTbA->_cZPTU_left();
    _cZPUU_R =hwTTbA->_cZPUU_right();
    _cZPUU_L =hwTTbA->_cZPUU_left();  
    _cZPCC_R =hwTTbA->_cZPCC_right();
    _cZPCC_L =hwTTbA->_cZPCC_left();
    _models =hwTTbA->_model();

  }
  FFVVertex::doinit();
}

void TTbAModelZPQQVertex::persistentOutput(PersistentOStream & os) const {
  os << _cZPTU_R << _cZPTU_L << _cZPUU_R << _cZPUU_L << _cZPCC_R << _cZPCC_L << _models;
}

void TTbAModelZPQQVertex::persistentInput(PersistentIStream & is, int) {
  is >> _cZPTU_R >> _cZPTU_L >> _cZPUU_R >> _cZPUU_L >> _cZPCC_R >> _cZPCC_L >> _models;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TTbAModelZPQQVertex,FFVVertex>
describeHerwigTTbAModelZPQQVertex("Herwig::TTbAModelZPQQVertex", "Herwig.so");


void TTbAModelZPQQVertex::Init() {
  
  static ClassDocumentation<TTbAModelZPQQVertex> documentation
    ("The TTbAModelZPQQVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " Z prime Quark-antiQuark vertex.");
}

void TTbAModelZPQQVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  double _cR = 0, _cL = 0;
  if( abs(aa->id()) == 6 || abs(bb->id()) == 6 || abs(cc->id()) == 6) { 
    _cR = _cZPTU_R; 
    _cL = _cZPTU_L; 
  } else {
    if( abs(aa->id()) != 4 && abs(bb->id()) != 4 && abs(cc->id()) != 4) { 
      _cR = _cZPUU_R; 
      _cL = _cZPUU_L;
    }
    if( abs(aa->id()) == 4 || abs(bb->id()) == 4 || abs(cc->id()) == 4) { 
      _cR = _cZPCC_R; 
      _cL = _cZPCC_L;
    }
  }
  
  if(_models!=1) { _cL = 1E-10; _cR = 1E-10; }
  right(_cR);
  left(_cL);

  norm(1.0);
}
#line 1 "./TTbAModelAGQQVertex.cc"
// -*- C++ -*-
//
// TTbAModelAGQQVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelAGQQVertex class.
//

#include "TTbAModelAGQQVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr TTbAModelAGQQVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelAGQQVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelAGQQVertex::TTbAModelAGQQVertex()  {
  orderInGem(1);
  orderInGs(1);
  addToList(-1,1,63);
  addToList(-2,2,63);
  addToList(-3,3,63);
  addToList(-4,4,63);
  addToList(-5,5,63);
  addToList(-6,6,63);
  colourStructure(ColourStructure::DELTA);
}

void TTbAModelAGQQVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _cAGQQ_R =hwTTbA->_cAGQQ_right();
    _cAGQQ_L =hwTTbA->_cAGQQ_left();
    _cAGTT_R =hwTTbA->_cAGTT_right();
    _cAGTT_L =hwTTbA->_cAGTT_left();
    _models = hwTTbA->_model();
  }
  FFVVertex::doinit();
}

void TTbAModelAGQQVertex::persistentOutput(PersistentOStream & os) const {
  os << _cAGQQ_R << _cAGQQ_L << _cAGTT_R << _cAGTT_L << _models;
}

void TTbAModelAGQQVertex::persistentInput(PersistentIStream & is, int) {
  is >> _cAGQQ_R >> _cAGQQ_L >>_cAGTT_R >> _cAGTT_L >> _models;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TTbAModelAGQQVertex,FFVVertex>
describeHerwigTTbAModelAGQQVertex("Herwig::TTbAModelAGQQVertex", "Herwig.so");


void TTbAModelAGQQVertex::Init() {
  
  static ClassDocumentation<TTbAModelAGQQVertex> documentation
    ("The TTbAModelAGQQVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " quark-lepton vertex.");
}


void TTbAModelAGQQVertex::setCoupling(Energy2 q2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {
  
  double _cL = 0, _cR = 0;
  double gstrong = 1.0;
  gstrong = strongCoupling(q2);

  if(abs(aa->id()) == 63 || abs(bb->id()) == 63 || abs(cc->id()) == 63) {
    if(abs(aa->id()) !=6 && abs(bb->id()) !=6 && abs(cc->id()) != 6) {
      _cR = _cAGQQ_R; 
      _cL = _cAGQQ_L; 
    } else { 
      _cR = _cAGTT_R; 
      _cL = _cAGTT_L; 
    }
    
  }

  if(_models!=2) { _cL = 1E-10; _cR = 1E-10; }
  left(_cL);
  right(_cR);
  
  norm(gstrong);
}
#line 1 "./TTbAModelSU2XVertex.cc"
// -*- C++ -*-
//
// TTbAModelSU2XVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TTbAModelSU2XVertex class.
//

#include "TTbAModelSU2XVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;



IBPtr TTbAModelSU2XVertex::clone() const {
  return new_ptr(*this);
}

IBPtr TTbAModelSU2XVertex::fullclone() const {
  return new_ptr(*this);
}

TTbAModelSU2XVertex::TTbAModelSU2XVertex()  {
  orderInGem(1);
  orderInGs(1);
  
  addToList(-6,6,70);
  addToList(-2,2,70);
  addToList(-2,6,70);
  addToList(2,-6,70);

  addToList(-6,6,71);
  addToList(-2,2,71);
  addToList(-2,6,71);
  addToList(2,-6,71);

   
  addToList(-6,6,-71);
  addToList(-2,2,-71);
  addToList(-2,6,-71);
  addToList(2,-6,-71);

  colourStructure(ColourStructure::DELTA);
}

void TTbAModelSU2XVertex::doinit() {
  _theModel = generator()->standardModel();
  tcHwTTbAPtr hwTTbA=dynamic_ptr_cast<tcHwTTbAPtr>(_theModel);
  if(hwTTbA) {
    _alphaX =hwTTbA->_alphaX_value();
    _costhetaX =hwTTbA->_costhetaX_value();
    _models =hwTTbA->_model();

  }
  FFVVertex::doinit();
}

void TTbAModelSU2XVertex::persistentOutput(PersistentOStream & os) const {
  os << _alphaX << _costhetaX << _models;
}

void TTbAModelSU2XVertex::persistentInput(PersistentIStream & is, int) {
  is >> _alphaX >> _costhetaX >> _models;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TTbAModelSU2XVertex,FFVVertex>
describeHerwigTTbAModelSU2XVertex("Herwig::TTbAModelSU2XVertex", "Herwig.so");


void TTbAModelSU2XVertex::Init() {
  
  static ClassDocumentation<TTbAModelSU2XVertex> documentation
    ("The TTbAModelSU2XVertex class is the implementation"
     " of the helicity amplitude calculation of the TTbA"
     " SU(2)_X vertex.");
}

void TTbAModelSU2XVertex::setCoupling(Energy2,tcPDPtr aa ,tcPDPtr bb, tcPDPtr cc) {

  double _cR = 0, _fac = 0;
  _gX = sqrt( 4 * Constants::pi * _alphaX ); 
  double ct = _costhetaX;
  double st = sqrt(1 - pow(ct,2)); 

  //Vz
  if( abs(aa->id()) == 70 || abs(bb->id()) == 70 || abs(cc->id()) == 70) { 
    _fac = _gX / 2.0; 
    if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = pow(ct,2) - pow(st,2);
      }
    }
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
       if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
	_cR = pow(st,2) - pow(ct,2);
      }
    }
    if( abs(aa->id()) == 2 || abs(bb->id()) == 2 || abs(cc->id()) == 2) {
      if( abs(aa->id()) == 6 || abs(bb->id()) == 6 || abs(cc->id()) == 6) {
	_cR = 2 * ct * st;
      }
    }
  }
  
  //Ym
  if( aa->id() == -71 || bb->id() == -71 || cc->id() == -71 ) { 
    
    _fac = _gX / sqrt(2.0); 
    
    if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = - ct * st;
      }
    }
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
      if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
	_cR = ct * st;
      }
    }    
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = - pow(st,2); 
      }
    }
    if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
      if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
	_cR = pow(ct,2);
      }
    }   
  }


  //Yp 
  if( aa->id() == 71 || bb->id() == 71 || cc->id() == 71 ) { 

    _fac = _gX / sqrt(2.0); 

    if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR =  - ct * st;
      }
    }
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
       if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
	_cR = ct * st;
      }
    }
    
    if( aa->id() == 2 || bb->id() == 2 || cc->id() == 2) {
      if( aa->id() == -6 || bb->id() == -6 || cc->id() == -6) {
	_cR = pow(ct,2);
      }
    }
    if( aa->id() == -2 || bb->id() == -2 || cc->id() == -2) {
      if( aa->id() == 6 || bb->id() == 6 || cc->id() == 6) {
	_cR = - pow(st,2);
      }
    } 
  }

  //normalise according to Lagrangian factor
  _cR *= _fac;


  //If this model is not selected set coupling to zero.
  if(_models!=3) { _cR = 1E-10; }
 

  right(_cR);
  left(0.);
 
  norm(1.0);

}
