#line 1 "./NMSSM.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSM class.
//

#include "NMSSM.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Reference.h"

using namespace Herwig;

void NMSSM::persistentOutput(PersistentOStream & os) const {
  os << _lambda << _kappa << ounit(_theAlambda,GeV) 
     << ounit(_theAkappa, GeV) << ounit(_lambdaVEV, GeV)
     << ounit(_MQ3, GeV) << ounit(_MU2, GeV);
}

void NMSSM::persistentInput(PersistentIStream & is, int) {
  is >> _lambda >> _kappa >> iunit(_theAlambda,GeV) 
     >> iunit(_theAkappa, GeV) >> iunit(_lambdaVEV, GeV)
     >> iunit(_MQ3, GeV) >> iunit(_MU2, GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSM,MSSM>
describeHerwigNMSSM("Herwig::NMSSM", "HwSusy.so HwNMSSM.so");

void NMSSM::Init() {

  static ClassDocumentation<NMSSM> documentation
    ("The NMSSM class is the base class for the NMSSM model");

}

void NMSSM::extractParameters(bool checkmodel) {
  MSSM::extractParameters(false);
  if(checkmodel) {
    map<string,ParamMap>::const_iterator pit;
    pit = parameters().find("modsel");
    if(pit == parameters().end()) return;
    ParamMap::const_iterator it;
    // nmssm or mssm
    it = pit->second.find(3);
    int inmssm = (it != pit->second.end()) ? int(it->second) : 0;
    if(inmssm == 0) 
      throw Exception() << "R-parity, CP and flavour conserving NMSSM model"
			<< " used but MSSM read in." << Exception::runerror; 
    // RPV
    it = pit->second.find(4);
    int irpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(irpv != 0) throw Exception() << "NMSSM model does not support RPV"
				  << Exception::runerror; 
    // CPV
    it = pit->second.find(5);
    int icpv = (it != pit->second.end()) ? int(it->second) : 0;
    if(icpv != 0) throw Exception() << "NMSSM model does not support CPV" 
				  << Exception::runerror; 
    // flavour violation
    it = pit->second.find(6);
    int ifv = (it != pit->second.end()) ? int(it->second) : 0;
    if(ifv != 0) throw Exception() << "NMSSM model does not support "
				 << "flavour violation"
				 << Exception::runerror;
  }
  // get the NMSSM parameters
  map<string,ParamMap>::const_iterator pit;
  pit=parameters().find("msoft");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it;
    it = pit->second.find(43);
    if(it != pit->second.end()) _MQ3 = it->second*GeV;
    it = pit->second.find(46);
    if(it != pit->second.end()) _MU2 = it->second*GeV;
  }
  pit=parameters().find("nmssmrun");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it = pit->second.find(1);
    if(it != pit->second.end()) _lambda = it->second;
    it = pit->second.find(2);
    if(it != pit->second.end()) _kappa = it->second;
    it = pit->second.find(3);
    if(it != pit->second.end()) _theAlambda = it->second*GeV;
    it = pit->second.find(4);
    if(it != pit->second.end()) _theAkappa = it->second*GeV;
    it = pit->second.find(5);
    if(it != pit->second.end()) _lambdaVEV = it->second*GeV;
  }
  pit=parameters().find("extpar");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it = pit->second.find(61);
    if(_lambda==ZERO     && it != pit->second.end()) _lambda = it->second;
    it = pit->second.find(62);
    if(_kappa==ZERO      && it != pit->second.end()) _kappa = it->second;
    it = pit->second.find(63);
    if(_theAlambda==ZERO && it != pit->second.end()) _theAlambda = it->second*GeV;
    it = pit->second.find(64);
    if(_theAkappa==ZERO  && it != pit->second.end()) _theAkappa = it->second*GeV;
    it = pit->second.find(65);
    if(_lambdaVEV==ZERO  && it != pit->second.end()) _lambdaVEV = it->second*GeV;
    it = pit->second.find(43);
    if(_MQ3==ZERO        && it != pit->second.end()) _MQ3 = it->second*GeV;
    it = pit->second.find(46);
    if(_MU2==ZERO        && it != pit->second.end()) _MU2 = it->second*GeV;
  }
  else {
    throw Exception() << "NMSSM::extractParameters - There was no EXTPAR block "
		      << "in the extracted parameters list. The model cannot "
		      << "be used without these." << Exception::runerror;
  }
  pit=parameters().find("msoft");
  if( pit != parameters().end() ) {
    ParamMap::const_iterator it;
    if(_MQ3==ZERO) {
      it = pit->second.find(43);
      if(it != pit->second.end()) _MQ3 = it->second*GeV;
    }
    if(_MU2==ZERO) {
      it = pit->second.find(46);
      if(it != pit->second.end()) _MU2 = it->second*GeV;
    }
  }
}

void NMSSM::createMixingMatrices() {
  map<string,pair<MatrixSize, MixingVector> >::const_iterator it;
  for(it=mixings().begin();it!=mixings().end();++it) {
    string name=it->first;
    // pseudo-scalar higgs mixing
    if (name == "nmamix") {
      MixingMatrixPtr temp;
      createMixingMatrix(temp,name,it->second.second,it->second.first);
      CPoddHiggsMix(temp);
    }
  }
  // base class for neutralinos and charginos
  MSSM::createMixingMatrices();
}
#line 1 "./NMSSMFFHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMFFHVertex class.
//

#include "NMSSMFFHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMFFHVertex::NMSSMFFHVertex() : _mw(0.*MeV), _sinb(0.), _cosb(0.), 
				   _tanb(0.), _idlast(make_pair(0,0)),
				   _q2last(0.*MeV2), 
				   _masslast(make_pair(0.*MeV,0*MeV)),
				   _couplast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void NMSSMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _mixS << _mixP << ounit(_mw,GeV)
     << _sinb << _cosb << _tanb << _sw << _theSM;
}

void NMSSMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixS >> _mixP >> iunit(_mw,GeV)
     >> _sinb >> _cosb >> _tanb >> _sw >> _theSM;
}

void NMSSMFFHVertex::doinit() {
  // the quarks and neutral higgs
  int in[5]={25,35,45,36,46};
  for(unsigned int iy=0;iy<5;++iy)
    for(int ix=1;ix<7;++ix)
      addToList( -ix, ix, in[iy] );

  // leptons and neutral higgs
  for(unsigned int iy=0;iy<5;++iy)
    for(int ix=11;ix<17;ix+=2)
      addToList( -ix, ix, in[iy] );

  // the quarks  and the charged higgs
  //H-
  for(int ix=0;ix<3;++ix) 
    addToList(2*ix+2, -2*ix-1, -37);

  //H+
  for(int ix=0;ix<3;++ix)
    addToList(-(2*ix+2), 2*ix+1, 37);

  // the leptons and the charged higgs
  //H-
  for(int ix=0;ix<3;++ix)
    addToList( 2*ix+12, -2*ix-11, -37 );

  //H+
  for(int ix=0;ix<3;++ix)
    addToList( -(2*ix+12), 2*ix+11, 37 );
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMFFHVertex::doinit()"
			  << Exception::runerror;
  _theSM = model;
  // sin theta_W
  double sw2=_theSM->sin2ThetaW();
  _sw = sqrt(sw2);
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMFFHVertex::doinit()" 
				   << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMFFHVertex::doinit()" 
				   << Exception::runerror;
  // Mass of the W boson
  _mw=getParticleData(ParticleID::Wplus)->mass();
  // sin and cos beta
  _tanb = model->tanBeta();
  double beta = atan(_tanb);
  _sinb=sin(beta);
  _cosb=cos(beta);
  // base class
  FFSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMFFHVertex,FFSVertex>
describeHerwigNMSSMFFHVertex("Herwig::NMSSMFFHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMFFHVertex::Init() {

  static ClassDocumentation<NMSSMFFHVertex> documentation
    ("The NMSSMFFHVertex class implements the vertex for the couplings"
     " of the Higgs bosons of the NMSSM to Standard Model fermions");

}
//calulate the couplings
void NMSSMFFHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int ihiggs=c->id();
  int id(abs(a->id()));
  Complex output(1.);
  // neutral Higgs
  if(ihiggs==25||ihiggs==35||ihiggs==45||ihiggs==36||ihiggs==46) {
    if(_idlast.first!=id||q2!=_q2last) {
      _idlast.first=id;
      _masslast.first = _theSM->mass(q2,a);
    }
    output = _masslast.first/_mw;
    // CP-even 
    if(ihiggs==25||ihiggs==35||ihiggs==45) {
      int iloc = (ihiggs-25)/10;
      output *= (id%2==0) ? (*_mixS)(iloc,1)/_sinb : (*_mixS)(iloc,0)/_cosb;
      left(1.); right(1.);
    } 
    // CP-odd
    else {
      int iloc = (ihiggs-36)/10;
      output *= (id%2==0) ? (*_mixP)(iloc,1)/_sinb : (*_mixP)(iloc,0)/_cosb;
      left(1.); right(-1.);
      output *= Complex(0., 1.);
    }
  }
  // Charged higgs
  else if(abs(ihiggs)==37) {
    output *= -sqrt(2.);
    int id2=abs(b->id());
    if(id2<id) {
      swap(id,id2);
      swap(a,b);
    }
    if(_idlast.first!=id||_idlast.second!=id2||q2!=_q2last) {
      _idlast.first =id ;
      _idlast.second=id2;
      _masslast.first  = _theSM->mass(q2,a);
      _masslast.second = _theSM->mass(q2,b);
    }
    double rgt = _masslast.first *_tanb/_mw;
    double lft = _masslast.second/_tanb/_mw;
    if(ihiggs>0) swap(lft,rgt);
    right(rgt);
    left (lft);
  }
  else {
    throw Exception() << "Unknown Higgs boson, PDG code = " << ihiggs 
		      << "in NMSSMFFHVertex::setCoupling()"
		      << Exception::runerror;
  }
  // prefactor
  if(q2!=_q2last) {
    _couplast = 0.5*weakCoupling(q2);
    _q2last=q2;
  }
  norm(-_couplast*output);
}
#line 1 "./NMSSMWWHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWWHVertex class.
//

#include "NMSSMWWHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMWWHVertex::NMSSMWWHVertex() 
  : _couplast(0.), _q2last(), _mw(), _zfact(0.), _sinb(0.),_cosb(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMWWHVertex::doinit() {
  int id[3]={25,35,45};
  // PDG codes for the particles in the vertex
  for(unsigned int ix=0;ix<3;++ix) {
    // Higgs WW
    addToList( 24, -24, id[ix] );
    //Higgs ZZ
    addToList( 23, 23, id[ix] );
  }
  // SM parameters
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  double sw2 = sin2ThetaW();
  _zfact = 1./(1.-sw2);
  // NMSSM parameters
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMWWHVertex::doinit()"
			  << Exception::runerror;
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMWWHVertex::doinit()" 
				   << Exception::runerror;
  // sin and cos beta
  double beta = atan(model->tanBeta());
  _sinb=sin(beta);
  _cosb=cos(beta);
  // base class
  VVSVertex::doinit();
}

void NMSSMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw,GeV) << _zfact << _sinb << _cosb << _mixS;
}

void NMSSMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw,GeV) >> _zfact >> _sinb >> _cosb >> _mixS;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMWWHVertex,VVSVertex>
describeHerwigNMSSMWWHVertex("Herwig::NMSSMWWHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMWWHVertex::Init() {

  static ClassDocumentation<NMSSMWWHVertex> documentation
    ("The NMSSMWWHVertex class implements the coupling of two electroweak gauge"
     " bosons with the Higgs bosons of the NMSSM");

}
//calulate couplings
void NMSSMWWHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr, tcPDPtr c) {
  // ID of gauge bosons
  int ibos=abs(a->id());
  // ID of Higgs
  int ihiggs = (c->id()-25)/10;
  // first the overall normalisation
  if(q2!=_q2last) {
    _couplast = weakCoupling(q2)*_mw*UnitRemoval::InvE;
    _q2last=q2;
  }
  // higgs mixing factor
  Complex hmix = _cosb*(*_mixS)(ihiggs,0)+_sinb*(*_mixS)(ihiggs,1);
  // couplings
  if(ibos==24)      norm(_couplast*hmix);
  else if(ibos==23) norm(_couplast*hmix*_zfact);
  else assert(false);
}
#line 1 "./NMSSMWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWHHVertex class.
//

#include "NMSSMWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMWHHVertex::NMSSMWHHVertex() : _sinb(0.), _cosb(0.), _sw(0.), _cw(0.),
				   _q2last(0.*MeV2), _couplast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMWHHVertex::doinit() {
  // codes for the neutral higgs
  //CP even
  int ieven[3]={25,35,45};
  //CP odd
  int iodd [2]={36,46};
  // Z CP even CP odd
  for(unsigned int ix=0;ix<3;++ix)
    for(unsigned int iy=0;iy<2;++iy)
      addToList( 23, ieven[ix], iodd[iy] );

  // W H+ CP even
  for(unsigned int ix=0;ix<3;++ix)
    addToList( -24, 37, ieven[ix] );

   // W+ H- CP even
  for(unsigned int ix=0;ix<3;++ix)
    addToList( 24, -37, ieven[ix] );

  // W H+ CP odd
  for(unsigned int ix=0;ix<2;++ix)
    addToList( -24, 37, iodd[ix] );

  //W+ H- CP odd
  for(unsigned int ix=0;ix<2;++ix)
    addToList( 24, -37, iodd[ix] );

  // Charged higgs Z/gamma
  addToList( 22, 37, -37 );
  addToList( 23, 37, -37 );
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMFFHVertex::doinit()"
			  << Exception::runerror;
  // sin theta_W
  double sw2 = sin2ThetaW();
  _sw = sqrt(sw2);
  _cw = sqrt(1.-sw2);
  // get the mixing matrices
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMWHHVertex::doinit()" 
				   << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMWHHVertex::doinit()" 
				   << Exception::runerror;
  // sin and cos beta
  double beta = atan(model->tanBeta());
  _sinb = sin(beta);
  _cosb = cos(beta);
  // base class
  VSSVertex::doinit();
}

void NMSSMWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << _sinb << _cosb << _sw << _cw << _mixS << _mixP;
}

void NMSSMWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _sinb >> _cosb >> _sw >> _cw >> _mixS >> _mixP;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMWHHVertex,VSSVertex>
describeHerwigNMSSMWHHVertex("Herwig::NMSSMWHHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMWHHVertex::Init() {

  static ClassDocumentation<NMSSMWHHVertex> documentation
    ("The NMSSMWHHVertex class implements the coupling of an electroweak"
     " gauge boson with two Higgs bosons in the NMSSM.");

}

//calulate the couplings
void NMSSMWHHVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,tcPDPtr c) {
  // weak coupling
  if(q2!=_q2last) {
    _couplast = weakCoupling(q2);
    _q2last=q2;
  }
  // gauge bosons
  int ibos= a->id();
  int ih1 = b->id();
  int ih2 = c->id();
  Complex fact;
  if(ibos==ParticleID::Z0) {
    fact = 0.5/_cw;
    // Z H+ H-
    if(abs(ih1)==37) {
      fact *= (sqr(_cw)-sqr(_sw));
      if(ih1<0) fact *=-1.;  
    }
    // Z CP even CP odd
    else {
      if(ih1%10==6) {
	fact *= -1.; 
	swap(ih1,ih2);
      }
      int is = (ih1-25)/10;
      int ip = (ih2-36)/10;
      fact *= Complex(0.,1.)*((*_mixS)(is,1)*(*_mixP)(ip,1)-
			      (*_mixS)(is,0)*(*_mixP)(ip,0));
    }
  }
  // gamma CP even CP odd
  else if(ibos==ParticleID::gamma) {
    fact = ih1>0 ? _sw : -_sw;  
  }
  // W boson
  else {
    fact = 0.5; 
    if(abs(ih2)==37) {
      swap(ih1,ih2);
      fact*=-1; 
    }
    // H+ CP even
    if(ih2%5==0) {
      int is = (ih2-25)/10;
      fact *= (_cosb*(*_mixS)(is,1)-_sinb*(*_mixS)(is,0));
      if(ibos<0) fact*=-1; 
    }
    // H+ CP odd
    else {
      int ip = (ih2-36)/10;
      fact *=-Complex(0.,1.)*(_cosb*(*_mixP)(ip,1)+_sinb*(*_mixP)(ip,0));
    }
  }
  //output the coupling
  norm(_couplast*fact);
}
#line 1 "./NMSSMHSFSFVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMHSFSFVertex class.
//

#include "NMSSMHSFSFVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMHSFSFVertex::NMSSMHSFSFVertex() : 
  _triTp(0.*MeV), _triBt(0.*MeV), _triTa(0.*MeV), _lambda(0.),
  _lambdaVEV(0.*MeV), _v1(0.*MeV), _v2(0.*MeV), _sw(0.), _cw(0.), 
  _mw(0.*MeV), _mz(0.*MeV), _sb(0.), _cb(0.), _tb(0.), _q2last(0.*MeV2), 
  _couplast(0.), _masslast(make_pair(0.*MeV,0.*MeV)), _idlast(make_pair(0,0)) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void NMSSMHSFSFVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _mixS << _mixP << _mixTp << _mixBt << _mixTa
     << ounit(_triTp,GeV) << ounit(_triBt,GeV) << ounit(_triTa,GeV) 
     << _lambda << ounit(_lambdaVEV,GeV) << ounit(_v1,GeV) << ounit(_v2,GeV)
     << _sw << _cw << ounit(_mw,GeV) << ounit(_mz,GeV) << _sb << _cb
     << _tb;
}


void NMSSMHSFSFVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _mixS >> _mixP >> _mixTp >> _mixBt >> _mixTa
     >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) >> iunit(_triTa,GeV) 
     >> _lambda >> iunit(_lambdaVEV,GeV) >> iunit(_v1,GeV) >> iunit(_v2,GeV) 
     >> _sw >> _cw >> iunit(_mw,GeV) >> iunit(_mz,GeV) >> _sb >> _cb >> _tb;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMHSFSFVertex,Helicity::SSSVertex>
describeHerwigNMSSMHSFSFVertex("Herwig::NMSSMHSFSFVertex", "HwSusy.so HwNMSSM.so");

void NMSSMHSFSFVertex::Init() {

  static ClassDocumentation<NMSSMHSFSFVertex> documentation
    ("The coupling of Higgs bosons to sfermions in the MSSM.");

}

void NMSSMHSFSFVertex::doinit() {
  //CP even
  int even[3] = {25, 35, 45};
  for(size_t h = 0; h < 3; ++h ) {
    //squarks
    for(long q = 1; q < 7; ++q) {
      //11
      addToList(even[h], -1000000 - q, 1000000 + q);
      //22
      addToList(even[h], -2000000 - q, 2000000 + q);
      //12
      addToList(even[h], -1000000 - q, 2000000 + q);
      //21
      addToList(even[h], -2000000 - q, 1000000 + q);
    }
    //sleptons
    for(long l = 11; l < 17; ++l) {
      //11
      addToList(even[h], -1000000 - l, 1000000 + l);
      //no right handed sneutrinos
      if( l % 2 != 0 ) {
	//22
	addToList(even[h], -2000000 - l, 2000000 + l);
	//12
	addToList(even[h], -1000000 - l, 2000000 + l);
	//21
	addToList(even[h], -2000000 - l, 1000000 + l);
      }
    }
  }
  //CP odd
  int odd[2] = {36, 46};
  for(size_t h = 0; h < 2; ++h ) {
    //squarks
    for(long q = 1; q < 7; ++q) {
      //12
      addToList(odd[h], -1000000 - q, 2000000 + q);
      //21
      addToList(odd[h], -2000000 - q, 1000000 + q);
    }
    //sleptons
    for(long l = 11; l < 16; l += 2) {
      //12
      addToList(odd[h], -1000000 - l, 2000000 + l);
      //21
      addToList(odd[h], -2000000 - l, 1000000 + l);
    }
  }
  //charged higgs
  //squarks
  for(long q = 1; q < 4; ++q ) {
    //H-
    //LL
    addToList(-37, -2*q - 999999, 2*q + 1000000);
    //RR
    addToList(-37, -2*q - 1999999, 2*q + 2000000);
    //LR
    addToList(-37, -2*q - 999999, 2*q + 2000000);
    //RL
    addToList(-37, -2*q - 1999999, 2*q + 1000000);
    //H+
    //LL
    addToList(37, -2*q - 1000000, 2*q + 999999);
    //RR
    addToList(37, -2*q - 2000000, 2*q + 1999999);
    //LR
    addToList(37, -2*q - 1000000, 2*q + 1999999);
    //RL
    addToList(37, -2*q - 2000000, 2*q + 999999);
  }
  //sleptons
  //easier as there are no right handed sneutrinos
  for(long l = 11; l <= 15; l +=2 ) {
    //H-
    //LL
    addToList(-37, -l - 1000000, l + 1000001);
    //RL
    addToList(-37, -l - 2000000, l + 1000001);
    //H+
    //LL
    addToList(+37, -l - 1000001, l + 1000000);
    //RL
    addToList(+37, -l - 1000001, l + 2000000);
  }
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if( !nmssm )
    throw InitException() << "NMSSMHSFSFVertex::doinit() - The model pointer "
			  << "in this vertex is not an NMSSM one as it "
			  << "should be." << Exception::runerror;
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixTp = nmssm->stopMix();
  _mixBt = nmssm->sbottomMix();
  _mixTa = nmssm->stauMix();
  if( !_mixS || !_mixP || !_mixTp || !_mixBt || !_mixTa ) 
    throw InitException() 
      << "NMSSMHSFSFVertex::doinit() - One of the mixing matrix pointers is "
      << "null, cannot continue. CP-even: " << _mixS << "  CP-odd: " << _mixP
      << "  ~t: " << _mixTp << "  ~b: " << _mixBt << "  ~tau: " << _mixTa
      << Exception::runerror;

  _triTp = nmssm->topTrilinear();
  _triBt = nmssm->bottomTrilinear();
  _triTa = nmssm->tauTrilinear();

  _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();

  _sw = sin2ThetaW();
  _cw = sqrt( 1. - _sw);
  _sw = sqrt(_sw);
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();

  _tb = nmssm->tanBeta();
  double beta = atan(_tb);
  _sb = sin(beta);
  _cb = cos(beta);
  
  _v1 = sqrt(2.)*_mw*_cb;
  _v2 = sqrt(2.)*_mw*_sb;
  SSSVertex::doinit();
}

void NMSSMHSFSFVertex::setCoupling(Energy2 q2,tcPDPtr part1,
				   tcPDPtr part2, tcPDPtr part3) {
  // extract particle ids
  long higgs(part1->id()), isf1(part2->id()), isf2(part3->id());
  // higgs first
  if(abs(isf1)<100) swap(higgs,isf1);
  if(abs(isf2)<100) swap(higgs,isf2);
  // squark second
  if(isf1<0) swap(isf1,isf2);
  // check higgs
  assert( higgs == 25 || higgs == 35 || higgs == 45 || 
	  higgs == 36 || higgs == 46 || abs(higgs) == 37 );
  // abs of antisquark and check
  isf2 *=-1;
  assert(isf1>0&&isf2>0);
  // running coupling
  if( q2 != _q2last ) {
    _q2last = q2;
    _couplast = weakCoupling(q2);
  }
 
  //charged higgs
  if( abs(higgs) == 37 ) {
    norm(_couplast*chargedHiggs(q2, isf1, isf2));
    return;
  }
  // neutral higgs
  // L/R states of the sfermions
  unsigned int alpha = ( isf1 > 2000000 ) ? 1 : 0;
  unsigned int beta  = ( isf2 > 2000000 ) ? 1 : 0;
  // id nad mass of corresponding SM fermion
  long smid = ( alpha == 0 ) ? isf1 - 1000000 : isf1 - 2000000;
  if( q2 != _q2last || smid != _idlast.first) {
    _idlast.first = smid;
    _masslast.first = _theSM->mass(q2, getParticleData(smid));
  }
  double f1 = _masslast.first/_mw;
  complex<Energy> af(ZERO);
  Complex m1a(0.), m1b(0.), m2a(0.), m2b(0.);
  // mixing for down type squarks and charged sleptons
  if( smid % 2 != 0 ) {
    f1 /= _cb;
    // sbottom
    if( smid == 5 ) {
      m1a = (*_mixBt)(alpha, 0);
      m1b = (*_mixBt)(alpha, 1);
      m2a = (*_mixBt)(beta , 0) ;
      m2b = (*_mixBt)(beta , 1);
      af = _triBt;
    }
    // stau
    else if( smid == 15 ) {
      m1a = (*_mixTa)(alpha, 0);
      m1b = (*_mixTa)(alpha, 1);
      m2a = (*_mixTa)(beta , 0) ;
      m2b = (*_mixTa)(beta , 1);
      af = _triTa;
    }
    // 1st 2 generations
    else {
      m1a = (alpha == 0) ? 1. : 0.;
      m1b = (alpha == 0) ? 0. : 1.;
      m2a = (beta  == 0) ? 1. : 0.;
      m2b = (beta  == 0) ? 0. : 1.;
      af = ZERO;
    }
  }
  // mixing for up type squarks and sneutrions
  else {
    f1 /= _sb;
    // stop
    if( smid == 6 ) {	  
      m1a = (*_mixTp)(alpha, 0);
      m1b = (*_mixTp)(alpha, 1);
      m2a = (*_mixTp)(beta , 0);
      m2b = (*_mixTp)(beta , 1);
      af = _triTp;
    }
    // everything else
    else {
      m1a = (alpha == 0) ? 1. : 0.;
      m1b = (alpha == 0) ? 0. : 1.;
      m2a = (beta == 0) ? 1. : 0.;
      m2b = (beta == 0) ? 0. : 1.;
      af = 0.*MeV;
    }
  }
  // scalar higgs bosons
  complex<Energy> fact(ZERO);
  if( higgs == 25 || higgs == 35 || higgs == 45 ) {
    int iloc = (higgs - 25)/10;
    complex<Energy> f2 = 0.5*_mz*( - _cb*(*_mixS)(iloc,0) 
				   + _sb*(*_mixS)(iloc,1))/_cw;

    // down type squarks and charged sleptons
    if( smid % 2 != 0 ) {
      double ef = (smid < 7) ? -1./3. : -1.;
      fact = - f2*( (1. + 2.*ef*sqr(_sw))*m1a*m2a - 2.*ef*sqr(_sw)*m1b*m2b)
	- f1*_masslast.first*(*_mixS)(iloc,0)*(m1a*m2a + m1b*m2b) 
	- 0.5*f1*(( - _lambdaVEV*(*_mixS)(iloc,1) 
		    - _lambda*_v2*(*_mixS)(iloc,2)/_couplast 
		    +  af*(*_mixS)(iloc,0)) * 
		  (m2a*m1b + m1a*m2b) );
    }
    // up type squarks and sneutrinos
    else {
      double ef = (smid < 7) ? 2./3. : 0.;
      fact =  +f2*( (1. - 2.*ef*sqr(_sw))*m1a*m2a + 2.*ef*sqr(_sw)*m1b*m2b )
	- f1*_masslast.first*(*_mixS)(iloc,1)*(m1a*m2a + m1b*m2b)  
	-  0.5*f1*(( - _lambdaVEV*(*_mixS)(iloc,0) 
		     - _lambda*_v1*(*_mixS)(iloc,2)/_couplast
		     +  af*(*_mixS)(iloc,1) ) *
		   (m2a*m1b + m1a*m2b));
    }
  }
  // pseudo scalar
  else if( higgs == 36 || higgs == 46 ) {
    int iloc = (higgs - 36)/10;
    // down type squarks and charged sleptons
    if( smid % 2 != 0 ) {
      fact = -0.5*f1*Complex(0.0,1.0)*
	( _lambdaVEV*(*_mixP)(iloc,1) +
	  _lambda*_v2*(*_mixP)(iloc,2)/_couplast +
	  af*(*_mixP)(iloc,0) );
    }
    // up-type squarks and sneutrinos
    else {
      fact =-0.5*f1*Complex(0.0,1.0)*
	( _lambdaVEV *(*_mixP)(iloc,0) +
	  _lambda*_v1*(*_mixP)(iloc,2)/_couplast +
	  af*(*_mixP)(iloc,1));
    }
    if(alpha<beta) fact *= -1.;
  }
  norm(_couplast*fact*UnitRemoval::InvE);
}

Complex NMSSMHSFSFVertex::chargedHiggs(Energy2 q2, long id1, long id2) {
  //have id1 as up-type
  if( id1 % 2 != 0) swap(id1, id2);
  // sfermion L/R states
  unsigned int alpha = ( id1/1000000 == 2 ) ? 1 : 0;
  unsigned int beta  = ( id2/1000000 == 2 ) ? 1 : 0;
  // type of quarks
  long utype = (alpha == 0) ? id1 - 1000000 : id1 - 2000000;
  long dtype = ( beta == 0) ? id2 - 1000000 : id2 - 2000000;
  // compute the running masses
  if( q2 != _q2last || id1 != _idlast.first || id2 != _idlast.second) {
    _idlast.first  = id1;
    _idlast.second = id2;
    _masslast.first  = _theSM->mass(q2, getParticleData(utype) );
    _masslast.second = _theSM->mass(q2, getParticleData(dtype) );
  }
  Energy2 facta = 2.*sqr(_mw)*_sb*_cb;
  complex<Energy2> coupling(ZERO);
  // sleptons
  if( dtype == 11 || dtype == 13 || dtype == 15) {
    Complex l1b = 0., l2b = 0.;
    complex<Energy> tri(ZERO);
    // 1st 2 generations
    if (dtype == 11 || dtype == 13) {
      l1b = (beta == 0) ? 1.0 : 0.0;
      l2b = (beta == 0) ? 0.0 : 1.0;
    }
    // stau
    else {
      l1b = (*_mixTa)(beta, 0) ;
      l2b = (*_mixTa)(beta, 1);
      tri = _triTa;
    }
    coupling = ( l1b*(sqr(_masslast.second)*_tb - facta) +
		 l2b*_masslast.second*(tri*_tb + _lambdaVEV) );
  }
  // squarks
  else {
    Complex q1a(0.0), q1b(0.0), q2a(0.0), q2b(0.0);
    complex<Energy> triD(ZERO), triU(ZERO);
    // up-type bit
    // stop
    if(utype == 6){
      q1a =  (*_mixTp)(alpha, 0) ;
      q2a =  (*_mixTp)(alpha, 1);
      triU = _triTp;
    }
    // light
    else{
      q1a = (alpha == 0) ? 1.0 : 0.0;
      q2a = (alpha == 0) ? 0.0 : 1.0;
    }
    // down-type bit
    // sbottom
    if(dtype == 5){
      q1b =  (*_mixBt)(beta, 0) ;
      q2b =  (*_mixBt)(beta, 1);
      triD = _triBt;
    }
    // light
    else{
      q1b = (beta == 0) ? 1.0 : 0.0;
      q2b = (beta == 0) ? 0.0 : 1.0;
    }
    Energy mfu = _masslast.first;
    Energy mfd = _masslast.second;
    coupling = ( q1a*q1b*((sqr(mfd)*_tb + sqr(mfu)/_tb) - facta)
		 + q2a*q2b*mfu*mfd*(_tb + (1./_tb))
		 + q1a*q2b*mfd*(triD*_tb + _lambdaVEV)
		 + q2a*q1b*mfu*(triU/_tb + _lambdaVEV));
  }
  return coupling * UnitRemoval::InvE/_mw/sqrt(2.);
}
#line 1 "./NMSSMGOGOHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGOGOHVertex class.
//

#include "NMSSMGOGOHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMGOGOHVertex::NMSSMGOGOHVertex() : _lambda(0.), _kappa(0.), _sinb(0.),
				       _cosb(0.), _sw(0.), _cw(0.),
				       _q2last(0.*MeV2), _couplast(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMGOGOHVertex::persistentOutput(PersistentOStream & os) const {
   os << _mixV << _mixU << _mixN << _mixS << _mixP << _lambda << _kappa << _sinb
      << _cosb << _sw << _cw;
}

void NMSSMGOGOHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _mixV >> _mixU >> _mixN >> _mixS >> _mixP >> _lambda >> _kappa >> _sinb
     >> _cosb >> _sw >> _cw;
}

void NMSSMGOGOHVertex::doinit() {
  int ieven[3]={25,35,45};
  int iodd [2]={36,46};
  long ichar[2]={1000024,1000037};
  long ineut[5]={1000022,1000023,1000025,1000035,1000045};
  // CP-even charginos
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      for(unsigned int iz=0;iz<3;++iz) {
	addToList(-ichar[ix], ichar[iy], ieven[iz]);
      }
    }
  }
  // CP-odd charginos
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      for(unsigned int iz=0;iz<2;++iz) {
	addToList(-ichar[ix], ichar[iy], iodd [iz]);

      }
    }
  }
  // CP-even neutralinos
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      for(unsigned int iz=0;iz<3;++iz) {
	addToList( ineut[ix], ineut[iy], ieven[iz]);
        }
      }
    }
  // CP-odd  neutralinos
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<5;++iy) {
      for(unsigned int iz=0;iz<2;++iz) {
	addToList( ineut[ix], ineut[iy], iodd[iz]);
        }
      }
    }

  // charged higgs
  for(unsigned int ix=0;ix<5;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      addToList(ineut[ix], -ichar[iy], 37);

      addToList(ineut[ix], ichar[iy], -37);

    }
  }

   tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
     // SM parameters
  // sin theta_W

  double sw2=sin2ThetaW();
  _cw=sqrt(1.0 - sw2);
  _sw=sqrt(sw2);
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in "
			  << "NMSSMGOGOHVertex::doinit()"
			  << Exception::runerror;
  // get the mixing matrices
  // higgs
  _mixS=model->CPevenHiggsMix();
  if(!_mixS) 
    throw InitException() << "Mixing matrix for CP-even neutral Higgs"
			  << " bosons is not set in NMSSMGOGOHVertex::doinit()" 
			  << Exception::runerror;
  _mixP=model->CPoddHiggsMix();
  if(!_mixP) 
    throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
			  << " bosons is not set in NMSSMGOGOHVertex::doinit()" 
			  << Exception::runerror;
  // charginos
  _mixU = model->charginoUMix();
  _mixV = model->charginoVMix();
  if(!_mixU || !_mixV)
    throw InitException() << "NMSSMGOGOHVertex::doinit - "
			  << "A mixing matrix pointer is null.  U: " 
			  << _mixU << "  V: " << _mixV
			  << Exception::abortnow;
  // neutralinos
  _mixN  = model->neutralinoMix();
  if(!_mixN)
    throw InitException() << "NMSSMGOGOHVertex::doinit - The neutralino "
			  << "mixing matrix pointer is null." 
			  << Exception::abortnow;
  // kappa and lambda couplings
  _lambda = model->lambda();
  _kappa  = model->kappa();
  // sin and cos beta
  double beta = atan(model->tanBeta());
  _sinb=sin(beta);
  _cosb=cos(beta);
  FFSVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMGOGOHVertex,FFSVertex>
describeHerwigNMSSMGOGOHVertex("Herwig::NMSSMGOGOHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMGOGOHVertex::Init() {

  static ClassDocumentation<NMSSMGOGOHVertex> documentation
    ("The NMSSMGOGOHVertex class implements the couplings of the Higgs bosons"
     " of the NMSSM and the electroweak gauginos");

}

void NMSSMGOGOHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
				   tcPDPtr part3) {
  long id1(part1->id()), id2(part2->id()), 
    id3(part3->id()), ihigg(0), ig1(0), ig2(0);
  if( abs(id1) == 25 || abs(id1) == 35 || abs(id1) == 45 || 
      abs(id1) == 36 || abs(id1) == 46 || abs(id1) == 37 ) {
    ihigg = id1;
    ig1 = id2;
    ig2 = id3;
  }
  else if( abs(id2) == 25 || abs(id2) == 35 || 
	   abs(id2) == 45 ||abs(id2) == 36 ||abs(id2) == 46 || abs(id2) == 37  ) {
    ihigg = id2;
    ig1 = id1;
    ig2 = id3;
  }
  else if( abs(id3) ==25 || abs(id3) == 35 || 
	   abs(id3) == 45 ||abs(id3) == 36 ||abs(id3) == 46 || abs(id3) == 37 ) {
    ihigg = id3;
    ig1 = id1;
    ig2 = id2;
  }
  else {
    throw HelicityConsistencyError() 
      << "NMSSMGOGOHVertex::setCoupling - There is no higgs particle in "
      << "this vertex. Particles: " << id1 << " " << id2 << " " << id3
      << Exception::runerror;
    return;
  }
  
  // weak coupling
  if(q2!=_q2last) {
    _couplast = weakCoupling(q2);
    _q2last = q2;
  }
  double rt = sqrt(0.5);
  // CP-even neutral higgs
  if(ihigg == 25 || ihigg == 35 || ihigg == 45) {
    int iloc = (ihigg - 25)/10;
    // chargino
    if(abs(ig1) == 1000024 || abs(ig1) == 1000037) {
      if( ig1 < 0 ) swap(ig1, ig2);
      int ic1 = abs(ig1)==1000024 ? 0 : 1;
      int ic2 = abs(ig2)==1000024 ? 0 : 1;
      Complex coupL = -_lambda*rt*conj((*_mixS)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1))
	-_couplast*rt*(conj((*_mixS)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0) +
			    (*_mixS)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1)));
      Complex coupR = -_lambda*rt*(*_mixS)(iloc,2)*(*_mixU)(ic2,1)*(*_mixV)(ic1,1)
	-_couplast*rt*((*_mixS)(iloc,0)*(*_mixU)(ic2,1)*(*_mixV)(ic1,0)+
		       (*_mixS)(iloc,1)*(*_mixU)(ic2,0)*(*_mixV)(ic1,1));
      left(coupL);
      right(coupR);
      norm(1.0);
    }
    // neutralino
    else  {
      int in1 = (ig1 < 1000024) ? (ig1 - 1000022) : (ig1 - 1000005)/10; 
      int in2 = (ig2 < 1000024) ? (ig2 - 1000022) : (ig2 - 1000005)/10;
      Complex us1 = (*_mixS)(iloc, 0), us2 = (*_mixS)(iloc, 1);
      Complex us3 = (*_mixS)(iloc, 2);
      Complex ni1 = (*_mixN)(in1,0), nj1 = (*_mixN)(in2,0);	  
      Complex ni2 = (*_mixN)(in1,1), nj2 = (*_mixN)(in2,1);
      Complex ni3 = (*_mixN)(in1,3), nj3 = (*_mixN)(in2,3); 
      Complex ni4 = (*_mixN)(in1,2), nj4 = (*_mixN)(in2,2); 
      Complex ni5 = (*_mixN)(in1,4), nj5 = (*_mixN)(in2,4);
      Complex YL =  
	- _lambda*rt*(us2*(ni4*nj5 + ni5*nj4) + 
		      us1*(ni3*nj5 + ni5*nj3) +
		      us3*(ni3*nj4 + ni4*nj3))
	+ sqrt(2.)*_kappa*us3*ni5*nj5 
	- _couplast*0.5*(us2*(ni2*nj3 + ni3*nj2) -
			 us1*(ni2*nj4 + ni4*nj2))
	+ _couplast*0.5*_sw*(us2*(ni1*nj3 + ni3*nj1) - 
			     us1*(ni1*nj4 + ni4*nj1) )/_cw;
      left(-conj(YL));
      right(-YL);
      norm(1.0);
    }
  }
  // CP-odd  neutral higgs
  else if(ihigg==36||ihigg==46) {
    int iloc = (ihigg-36)/10;
    // chargino
    if(abs(ig1)==1000024||abs(ig1)==1000037) {
	if( ig1 < 0 ) swap(ig1, ig2);
	int ic1 = abs(ig1)==1000024 ? 0 : 1;
	int ic2 = abs(ig2)==1000024 ? 0 : 1;
	Complex QL = Complex(0,-1.0)*
	  (_lambda*rt*conj((*_mixP)(iloc,2)*(*_mixU)(ic1,1)*(*_mixV)(ic2,1))
	   -_couplast*rt*conj(((*_mixP)(iloc,0)*(*_mixU)(ic1,1)*(*_mixV)(ic2,0) +
			       (*_mixP)(iloc,1)*(*_mixU)(ic1,0)*(*_mixV)(ic2,1))));
	Complex QR = Complex(0,-1.0)*
	  (_lambda*rt*(*_mixP)(iloc,2)*(*_mixU)(ic2,1)*(*_mixV)(ic1,1)
	   -_couplast*rt*((*_mixP)(iloc,0)*(*_mixU)(ic2,1)*(*_mixV)(ic1,0) +
			  (*_mixP)(iloc,1)*(*_mixU)(ic2,0)*(*_mixV)(ic1,1)));
	left(QL);
	right(-QR);
	norm(1.);
    }
    // neutralino
    else {
      int in1 = (ig1 < 1000024) ? (ig1 - 1000022) : (ig1 - 1000005)/10; 
      int in2 = (ig2 < 1000024) ? (ig2 - 1000022) : (ig2 - 1000005)/10;
      Complex up1 = (*_mixP)(iloc, 0), up2 = (*_mixP)(iloc, 1);
      Complex up3 = (*_mixP)(iloc, 2);
      Complex ni1 = (*_mixN)(in1,0), nj1 = (*_mixN)(in2,0);	  
      Complex ni2 = (*_mixN)(in1,1), nj2 = (*_mixN)(in2,1); 
      Complex ni3 = (*_mixN)(in1,2), nj3 = (*_mixN)(in2,2); 
      Complex ni4 = (*_mixN)(in1,3), nj4 = (*_mixN)(in2,3); 
      Complex ni5 = (*_mixN)(in1,4), nj5 = (*_mixN)(in2,4);
      Complex AL = 
	_lambda*rt*(up2*(ni3*nj5 + ni5*nj3) +
		    up1*(ni4*nj5 + ni5*nj4) + 
		    up3*(ni3*nj4 + ni4*nj3))
	- sqrt(2.)*_kappa*up3*ni5*nj5
	- _couplast*0.5*(up2*(ni2*nj4 + ni4*nj2) - 
			 up1*(ni2*nj3 + ni3*nj2))
	+ _couplast*0.5*_sw*(up2*(ni1*nj4 + ni4*nj1) - 
			     up1*(ni1*nj3 + ni3*nj1))/_cw;
      AL *= Complex(0.0, -1.0);
      left(conj(AL));
      right(AL);
      norm(1.);
    }
  }
  // charged higgs
  else {
    if (abs(ig1) == 1000024 || abs(ig1) == 1000037) swap (ig1,ig2);
    int in = (abs(ig1) < 1000024) ? (ig1-1000022) : (ig1-1000005)/10; 
    int ic = (abs(ig2) == 1000024) ? 0 : 1;
    Complex QpR = _lambda*_cosb*(*_mixU)(ic,1)*(*_mixN)(in,4)
      -_sinb*_couplast*(rt*(*_mixU)(ic,1)*(_sw*(*_mixN)(in,0)/_cw + (*_mixN)(in,1))
		 	 - (*_mixU)(ic,0)*(*_mixN)(in,2));
    Complex QpL = _lambda*_sinb*(*_mixV)(ic,1)*(*_mixN)(in,4)
      + _couplast*_cosb*(rt*(*_mixV)(ic,1)
			 *(_sw*(*_mixN)(in,0)/_cw + (*_mixN)(in,1))
			 + (*_mixV)(ic,0)*(*_mixN)(in,3));
    QpL = conj(QpL);
    if(ihigg > 0) {
      left (QpL);
      right(QpR);
      norm(-1.);
    }
    else {
      left (conj(QpR));
      right(conj(QpL));
      norm(-1.);
    }
  }
}
#line 1 "./NMSSMHHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMHHHVertex class.
//

#include "NMSSMHHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "NMSSM.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

NMSSMHHHVertex::NMSSMHHHVertex() : _mw(0.*MeV), _mz(0.*MeV), _sw2(0.),
				   _cw(0.), _lambda(0.), _kappa(0.) ,
				   _lambdaVEV(0.*MeV), _theAl(0.*MeV),
				   _theAk(0.*MeV), _sb(0.), _cb(0.),
				   _s2b(0.), _c2b(0.), _vu(0.*MeV),
				   _vd(0.*MeV), _s(0.*MeV), _q2last(0.*MeV2),
				   _glast(0.), _MQ3(0.*MeV), _MU2(0.*MeV),
				   _includeRadiative(false) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMHHHVertex::doinit() {
  // PDG codes for the particles in vertex _vd
  //CP-even Higgs
  addToList(25, 35, 45);
  for( unsigned int i = 25; i <= 45; i += 10 ) {
    addToList(i, i, 25);
    addToList(i, i, 35);
    addToList(i, i, 45);
    //Charged Higgs
    addToList(i, 37, -37);
    //CP-odd Higgs
    addToList(i, 36, 36);
    addToList(i, 36, 46);
    addToList(i, 46, 36);
    addToList(i, 46, 46);
  }
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  if( !nmssm ) 
    throw InitException() << "NMSSMHHHVertex::doinit - The model object is"
			  << "not the NMSSM object."
			  << Exception::runerror;  
  //SM parameters
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _sw2 = sin2ThetaW(); 
  _cw = sqrt(1. - _sw2);
  //NMSSM parameters
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  if( !_mixS || !_mixP ) 
    throw InitException() << "NMSSMHHHVertex::doinit - One of the mixing matrix "
			  << "pointers is null, cannot continue. S: "
			  << _mixS << "  P: " << _mixP << Exception::runerror;
  _lambda = nmssm->lambda();
  _kappa = nmssm->kappa();
  _lambdaVEV = nmssm->lambdaVEV();
  _theAl = nmssm->trilinearLambda();
  _theAk = nmssm->trilinearKappa();
  _MQ3 = nmssm->MQ3();  
  _MU2 = nmssm->MU2(); 
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);
  _vd = sqrt(2)*_mw*_cb;
  _vu = sqrt(2)*_mw*_sb;
  _s  = _lambdaVEV/_lambda;
  SSSVertex::doinit();
}

void NMSSMHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw, GeV) << ounit(_mz,GeV) 
     << _sw2 << _cw <<  _lambda << _includeRadiative
     << _kappa <<  ounit(_lambdaVEV,GeV) <<  ounit(_theAl, GeV) 
     << ounit(_theAk,GeV) <<  _sb <<  _cb << _s2b <<  _c2b
     << ounit(_vu,GeV) << ounit(_vd,GeV) << ounit(_s,GeV) << _mixS << _mixP
     << ounit(_MQ3,GeV) << ounit(_MU2,GeV) << _theSM;
}

void NMSSMHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw, GeV) >> iunit(_mz,GeV) 
     >> _sw2 >> _cw >>  _lambda >> _includeRadiative
     >> _kappa >>  iunit(_lambdaVEV,GeV) >>  iunit(_theAl, GeV) 
     >> iunit(_theAk,GeV) >>  _sb >>  _cb >> _s2b >>  _c2b
     >> iunit(_vu,GeV) >> iunit(_vd,GeV) >> iunit(_s,GeV)>> _mixS >> _mixP
     >> iunit(_MQ3,GeV) >> iunit(_MU2,GeV)  >> _theSM;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMHHHVertex,SSSVertex>
describeHerwigNMSSMHHHVertex("Herwig::NMSSMHHHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMHHHVertex::Init() {

  static ClassDocumentation<NMSSMHHHVertex> documentation
    ("This is the triple Higgs coupling in the NMSSM.");

  static Switch<NMSSMHHHVertex,bool> interfaceIncludeRadiativeCorrections
    ("IncludeRadiativeCorrections",
     "Include radiative corrections in the vertex",
     &NMSSMHHHVertex::_includeRadiative, false, false, false);
  static SwitchOption interfaceIncludeRadiativeCorrectionsYes
    (interfaceIncludeRadiativeCorrections,
     "Yes",
     "Include the radiative terms",
     true);
  static SwitchOption interfaceIncludeRadiativeCorrectionsNo
    (interfaceIncludeRadiativeCorrections,
     "No",
     "Don't include them",
     false);

}

//calulate couplings
void NMSSMHHHVertex::setCoupling(Energy2 q2,tcPDPtr p1,tcPDPtr p2, 
				 tcPDPtr p3) {
  using Constants::pi;
  long higgs[3] = {p1->id(), p2->id(), p3->id()};
  unsigned int ns(0), np(0), nc(0);
  for( int i = 0; i < 3; ++i ) {
    if( higgs[i] == 25 || higgs[i] == 35 || higgs[i] == 45 )
      ++ns;
    else if( higgs[i] == 36 || higgs[i] == 46 )
      ++np;
    else if( abs(higgs[i]) ==  37 )
      ++nc;
  }
  //check three Higgs in vertex
  assert( ns + np + nc == 3 );
  if( q2 != _q2last ) {
    _q2last = q2;
    _glast = weakCoupling(q2);
    _mb = _theSM->mass(q2,getParticleData(5));
    _mt = _theSM->mass(q2,getParticleData(6));
  }
  //define VEV's
  double rt = sqrt(0.5);
  Energy _mtpole = getParticleData(6)->mass();
  Energy2 Qstsb = _MQ3*_MU2;
  double radlog(0.);
  if(_includeRadiative) {
    radlog = Qstsb/sqr(_mtpole);
    assert(radlog!=0.);
    radlog = log(radlog);
  }
  complex<Energy> coupling;
  //CP even Higgs
  if( ns == 3 ) {
    unsigned int a = (higgs[0] - 25)/10;
    unsigned int b = (higgs[1] - 25)/10;
    unsigned int c = (higgs[2] - 25)/10;
    coupling = 
      sqr(_lambda)*rt*(_vu*(usMix(a,b,c,1,0,0) + usMix(a,b,c,1,2,2))/_glast +
		       _vd*(usMix(a,b,c,0,1,1) + usMix(a,b,c,0,2,2))/_glast +
		       _s *(usMix(a,b,c,2,1,1) + usMix(a,b,c,2,0,0)))
      - _lambda*_kappa*rt*(_vu*usMix(a,b,c,0,2,2)/_glast +
			   _vd*usMix(a,b,c,2,1,2)/_glast + 2.*_s*usMix(a,b,c,1,0,2))
      + sqr(_kappa)/rt*_s*usMix(a,b,c,2,2,2)
      - _lambda*_theAl*rt*usMix(a,b,c,1,0,2)
      + _kappa*_theAk*rt/3.*usMix(a,b,c,2,2,2)
      + sqr(_glast)*0.25*rt/sqr(_cw)*(_vu*(usMix(a,b,c,1,1,1) -
					   usMix(a,b,c,1,0,0))/_glast -
				      _vd*(usMix(a,b,c,0,1,1) -
					   usMix(a,b,c,0,0,0))/_glast);
    // additional radiative terms
    if(_includeRadiative) {
      complex <Energy> radtop = usMix(a,b,c,1,1,1)*3.0*sqrt(2.0)*radlog
	*sqr(_mt)*sqr(_mt)*sqr(_glast)*_glast/
	(16.0*sqr(pi)*_vu*_vu*_vu);
      complex <Energy> radbot= usMix(a,b,c,0,0,0)*3.0*sqrt(2.0)*radlog
	*sqr(_mb)*sqr(_mb)*sqr(_glast)*_glast
	/(16.0*sqr(pi)*_vd*_vd*_vd);
      coupling += radbot + radtop;
    }					
  }
  //CP even, CP odd Vertex
  else if(ns == 1 && np == 2) {
    unsigned int a(0), b(0), c(0);
    if( higgs[0] == 25 || higgs[0] == 35 || higgs[0] == 45 ) {
      a = (higgs[0] - 25)/10;
      b = (higgs[1] - 36)/10;
      c = (higgs[2] - 36)/10;
    }
    else if(higgs[1] == 25 || higgs[1] == 35 || higgs[1] == 45 ) {
      a = (higgs[1] - 25)/10;
      b = (higgs[0] - 36)/10;
      c = (higgs[2] - 36)/10;
    }
    else {
      a = (higgs[2] - 25)/10;
      b = (higgs[0] - 36)/10;
      c = (higgs[1] - 36)/10;
    }
    coupling =	
      sqr(_lambda)*rt*(_vu*(upMix(a,b,c,1,0,0) + upMix(a,b,c,1,2,2))/_glast +
		       _vd*(upMix(a,b,c,0,1,1) + upMix(a,b,c,0,2,2))/_glast +
		       _s *(upMix(a,b,c,2,1,1) + upMix(a,b,c,2,0,0)))
      + _lambda*_kappa*rt*(_vu*(upMix(a,b,c,0,2,2) 
				- 2.*upMix(a,b,c,2,0,2))/_glast +
			   _vd*(upMix(a,b,c,1,2,2)
				- 2.*upMix(a,b,c,2,1,2))/_glast 
			   + 2.*_s*(upMix(a,b,c,2,1,0)
				    - upMix(a,b,c,1,0,2) - upMix(a,b,c,0,1,2)))
      + sqr(_kappa)/rt*_s*upMix(a,b,c,2,2,2)
      +_lambda*_theAl*rt*(upMix(a,b,c,1,0,2)
			  + upMix(a,b,c,0,1,2) + upMix(a,b,c,2,1,0))
      - _kappa*_theAk*rt*upMix(a,b,c,2,2,2)
      + sqr(_glast)*0.25*rt/sqr(_cw)*(_vu*(upMix(a,b,c,1,1,1) -
					   upMix(a,b,c,1,0,0))/_glast - 
				      _vd*(upMix(a,b,c,0,1,1) -
					   upMix(a,b,c,0,0,0))/_glast);
    if(_includeRadiative) {
      complex <Energy> radtop = upMix(a,b,c,1,1,1)*3.0*sqrt(2.0)*radlog*
	sqr(_mt)*sqr(_mt)*sqr(_glast)*_glast/
	(16.0*sqr(pi)*_vu*_vu*_vu);
      complex <Energy> radbot= upMix(a,b,c,0,0,0)*3.0*sqrt(2.0)*radlog*
	sqr(_mb)*sqr(_mb)*sqr(_glast)*_glast
	/(16.0*sqr(pi)*_vd*_vd*_vd);
      coupling += radbot + radtop;
    }
  }
  //Charged Higgs
  else {
    unsigned int a(0);
    if( higgs[0] == 25 || higgs[0] == 35 || higgs[0] == 45 )
      a = (higgs[0] - 25)/10;
    else if(higgs[1] == 25 || higgs[1] == 35 || higgs[1] == 45 )
      a = (higgs[1] - 25)/10;
    else
      a = (higgs[2] - 25)/10;
    coupling = 	
      sqr(_lambda)*rt*2.*(_s*((*_mixS)(a,2)*sqr(_cb) + (*_mixS)(a,2)*sqr(_sb))
			  - (_vu*(*_mixS)(a,0)/_glast + 
			     _vd*(*_mixS)(a,1)/_glast)*_sb*_cb)
      +_lambda*_sb*_cb*2.*(*_mixS)(a,2)*(_kappa*_s/rt + rt*_theAl)
      + sqr(_glast)*0.5*rt*_sw2/sqr(_cw)*((_vu*(*_mixS)(a,1)/_glast - 
					   _vd*(*_mixS)(a,0)/_glast)*sqr(_cb) + 
					  (_vd*(*_mixS)(a,0)/_glast -
					   _vu*(*_mixS)(a,1)/_glast)*sqr(_sb))
      + sqr(_glast)*0.5*rt*(_vu*((*_mixS)(a,1)*sqr(_cb) + 
				 (*_mixS)(a,1)*sqr(_sb) + 
				 2.*(*_mixS)(a,0)*_cb*_sb)/_glast 
			    + _vd*((*_mixS)(a,0)*sqr(_cb) + 
				   (*_mixS)(a,0)*sqr(_sb) 
				   + 2.*(*_mixS)(a,1)*_sb*_cb)/_glast);
    if(_includeRadiative) {
      complex <Energy> radtop =(*_mixS)(a,1)*sqr(_sb)*6.0*sqrt(2.0)*radlog*
	sqr(_mt)*sqr(_mt)*sqr(_glast)*_glast/
	(16.0*sqr(pi)*_vu*_vu*_vu);
      complex <Energy> radbot=(*_mixS)(a,0)*sqr(_cb)*6.0*sqrt(2.0)*radlog*
	sqr(_mb)*sqr(_mb)*sqr(_glast)*_glast
	/(16.0*sqr(pi)*_vd*_vd*_vd);
      complex <Energy> temp2 = _vu*((*_mixS)(a,1)*sqr(_cb) + 
				    (*_mixS)(a,0)*_sb*_cb)/_glast+ 
	_vd*((*_mixS)(a,1)*_sb*_cb + 
	     (*_mixS)(a,0)*sqr(_sb))/_glast;				
      
      complex <Energy> radtopbot= temp2*6.0*sqrt(2.0)*radlog*
	sqr(_mt)*sqr(_mb)*sqr(_glast)*sqr(_glast)
	/(16.0*sqr(pi)*sqr(_vu)*sqr(_vd));
      coupling += radbot + radtop + radtopbot;					
    }
  }
  norm(-coupling * UnitRemoval::InvE);
}
#line 1 "./NMSSMGGHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMGGHVertex class.
//

#include "NMSSMGGHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/Susy/NMSSM/NMSSM.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

NMSSMGGHVertex::NMSSMGGHVertex() : _sw(0.), _cw(0.), _mw(0.*MeV),
	_mz(0.*MeV),_lambdaVEV(0.*MeV), _lambda(0.), _v1(0.*MeV),
	_v2(0.*MeV), _triTp(0.*MeV), _triBt(0.*MeV),
	_sb(0.), _cb(0.), _masslast(make_pair(0.*MeV,0.*MeV)),
	_q2last(0.*MeV2), _couplast(0.), _coup(0.),
    _hlast(0), _recalc(true) {
  orderInGem(1);
  orderInGs(2);
  colourStructure(ColourStructure::DELTA);
}

void NMSSMGGHVertex::doinit()  {
  addToList(21,21,25);
  addToList(21,21,35);
  addToList(21,21,36);
  addToList(21,21,45);
  addToList(21,21,46);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) {
    throw InitException() << "NMSSMGGHVertex::doinit - The SM pointer is null!"
			  << Exception::abortnow;
  }
  // SM parameters
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1. - sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _top = getParticleData(6);
  _bt = getParticleData(5);
  
  //NMSSM parameters
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixQt = nmssm->stopMix();
  _mixQb = nmssm->sbottomMix();
  
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);
  
  _v1 = sqrt(2.)*_mw*_cb;
  _v2 = sqrt(2.)*_mw*_sb;
 
  _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();

  _triTp = nmssm->topTrilinear();
  _triBt = nmssm->bottomTrilinear();

  // resize vectors here and use setNParticles method
  // to the set the actual number in the loop.
  // Also only the top mass hass to be calculated at runtime
  masses.resize(6, Energy());
  masses[0] = getParticleData(6)->mass();
  masses[1] = getParticleData(5)->mass();

  masses[2] = getParticleData(1000005)->mass();
  masses[3] = getParticleData(2000005)->mass();

  masses[4] = getParticleData(1000006)->mass();
  masses[5] = getParticleData(2000006)->mass();

  type.resize(6, PDT::Spin0);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1Half;
  couplings.resize(6);

  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void NMSSMGGHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _sw << _cw << ounit(_mw, GeV) << ounit(_mz, GeV)
     << ounit(_lambdaVEV,GeV) << _lambda << ounit(_v1,GeV) << ounit(_v2,GeV)
     << ounit(_triTp,GeV) << ounit(_triBt,GeV) 
     << _top << _bt << _mixS << _mixP << _mixQt << _mixQb << _sb << _cb; 
}


void NMSSMGGHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _sw >> _cw >> iunit(_mw, GeV) >> iunit(_mz, GeV)
     >> iunit(_lambdaVEV,GeV) >> _lambda >> iunit(_v1,GeV) >> iunit(_v2,GeV)
     >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) 
     >> _top >> _bt >> _mixS >> _mixP >> _mixQt >> _mixQb >> _sb >> _cb; 
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMGGHVertex,VVSLoopVertex>
describeHerwigNMSSMGGHVertex("Herwig::NMSSMGGHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMGGHVertex::Init() {

  static ClassDocumentation<NMSSMGGHVertex> documentation
    ("The effective coupling of a higgs to a pair of gluons in the "
     "NMSSM.");

}

void NMSSMGGHVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2,
				 tcPDPtr p3) {			 
  long hid(p3->id());
  if( q2 != _q2last ) {
    Looptools::clearcache();
    _couplast = sqr(strongCoupling(q2));
    _coup = weakCoupling(q2);
    _q2last = q2;
    _recalc = true;
  }
  norm(_couplast*_coup);
  // scalar higgs bosons  
  if( hid != _hlast ) {
    _hlast = hid;
    _recalc = true;
    if( hid % 5 == 0 ) {
      // location of the higgs
      int iloc = (hid - 25)/10;
      // 6 particles in the loop
      setNParticles(6);
      // top and bottom quark masses
      Energy mt = _theSM->mass(q2, _top);
      Energy mb = _theSM->mass(q2,  _bt);
      Complex c(0.);
      // couplings for the top quark loop
      c = -0.25*mt*(*_mixS)(iloc, 1)/_sb/_mw;
	
      couplings[0] = make_pair(c,c);
      masses[0] = mt;
      // couplings for the bottom quark loop
      c = -0.25*mb*(*_mixS)(iloc, 0)/_cb/_mw;
	
      couplings[1] = make_pair(c,c);	
      masses[1] = mb;
      // sbottoms
      double f1 = mb/_mw/_cb;
      complex<Energy>  f2 = 0.5*_mz/_cw*
	( - _cb*(*_mixS)(iloc,0) + _sb*(*_mixS)(iloc,1));
      complex<Energy> cpl;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl = -f2*( (1. - 2.*sqr(_sw)/3.)*(*_mixQb)(ix, 0)*(*_mixQb)(ix, 0)
		    + 2.*sqr(_sw)*(*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)/3.)
	  - f1*mb*(*_mixS)(iloc,0)
	  *((*_mixQb)(ix, 0)*(*_mixQb)(ix, 0) + (*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)) 
	  - 0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,1) - _lambda*_v2*(*_mixS)(iloc,2)/_coup 
		    +  _triBt*(*_mixS)(iloc,0))*((*_mixQb)(ix, 1)*(*_mixQb)(ix, 0)
						 + (*_mixQb)(ix, 0)*(*_mixQb)(ix, 1));

	couplings[2+ix] = make_pair(Complex(0.5*cpl*UnitRemoval::InvE),
				    Complex(0.5*cpl*UnitRemoval::InvE)); 
      }
      // stop
      f1 = mt/_mw/_sb;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl   =+f2*( (1. - 4.*sqr(_sw)/3.)*(*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
		     + 4.*sqr(_sw)*(*_mixQt)(ix, 1)*(*_mixQt)(ix, 1)/3.)
	  - f1*mt*(*_mixS)(iloc,1)
	  *((*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
	    + (*_mixQt)(ix, 1)*(*_mixQt)(ix, 1))  
	  -  0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,0) - _lambda*_v1*(*_mixS)(iloc,2)/_coup
		     + _triTp*(*_mixS)(iloc,1))*((*_mixQt)(ix, 1)*(*_mixQt)(ix, 0)
						 + (*_mixQt)(ix, 0)*(*_mixQt)(ix, 1));

	couplings[4+ix] = make_pair(Complex(0.5*cpl*UnitRemoval::InvE),
				    Complex(0.5*cpl*UnitRemoval::InvE));
      }
    }
    // pseudoscalar higgs bosons	
    else {
      // location of the higgs
      int iloc = (hid - 36)/10;
      // 2 particles in the loop
      setNParticles(2);
      // top and bottom quark masses
      Energy mt = _theSM->mass(q2, _top);
      Energy mb = _theSM->mass(q2,  _bt);
      Complex c(0.);
      // top quark couplings
      c = Complex(0.,-1.)*0.25*mt*(*_mixP)(iloc, 1)/_sb/_mw;
	  	  
      couplings[0] = make_pair(-c,c);
      masses[0] = mt;
      // bottom quark couplings
      c = Complex(0., -1.)*0.25*mb*(*_mixP)(iloc, 0)/_cb/_mw;
	  
      couplings[1] = make_pair(-c,c);	
      masses[1] = mb;
    }
  }

  if( _recalc ) {
    VVSLoopVertex::setCoupling(q2, p1, p2, p3);
    _recalc = false;
  } 
}
#line 1 "./NMSSMPPHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMPPHVertex class.
//

#include "NMSSMPPHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Models/Susy/NMSSM/NMSSM.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;

NMSSMPPHVertex::NMSSMPPHVertex() 
  : _sw(0.), _cw(0.), _mw(0.*MeV),
    _mz(0.*MeV),_lambdaVEV(0.*MeV), _lambda(0.),
    _triTp(0.*MeV), _triBt(0.*MeV),
    _sb(0.), _cb(0.), 
    _kappa(0.),_vu(ZERO),_vd(ZERO),_s(ZERO),_theAl(ZERO),
    _masslast(make_pair(0.*MeV,0.*MeV)),_q2last(0.*MeV2),
    _couplast(0.), _coup(0.), _hlast(0), _recalc(true) {
  orderInGem(3);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void NMSSMPPHVertex::doinit()  {
  addToList(22,22,25);
  addToList(22,22,35);
  addToList(22,22,36);
  addToList(22,22,45);
  addToList(22,22,46);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) {
    throw InitException() << "NMSSMPPHVertex::doinit - The SM pointer is null!"
			  << Exception::abortnow;
  }
  // SM parameters
  _sw = sqrt(sin2ThetaW());
  _cw = sqrt(1. - sin2ThetaW());
  _mw = getParticleData(24)->mass();
  _mz = getParticleData(23)->mass();
  _top = getParticleData(6);
  _bt  = getParticleData(5);
  _tau = getParticleData(15);
  
  //NMSSM parameters
  tcNMSSMPtr nmssm = dynamic_ptr_cast<tcNMSSMPtr>(_theSM);
  _mixS = nmssm->CPevenHiggsMix();
  _mixP = nmssm->CPoddHiggsMix();
  _mixQt = nmssm->stopMix();
  _mixQb = nmssm->sbottomMix();
  _mixLt = nmssm->stauMix();
  
  double beta = atan(nmssm->tanBeta());
  _sb = sin(beta);
  _cb = cos(beta);

  _lambda = nmssm->lambda();
  _lambdaVEV = nmssm->lambdaVEV();
  
  _triTp = nmssm->topTrilinear();
  _triBt = nmssm->bottomTrilinear();
  _triTa = nmssm->tauTrilinear();

  _vd = sqrt(2)*_mw*_cb;
  _vu = sqrt(2)*_mw*_sb;
  _s  = _lambdaVEV/_lambda;
  _theAl = nmssm->trilinearLambda();
  _kappa = nmssm->kappa();

  _mixU = nmssm->charginoUMix();
  _mixV = nmssm->charginoVMix();

  // resize vectors here and use setNParticles method
  // to the set the actual number in the loop.
  // Also only the top mass hass to be calculated at runtime
  masses.resize(13, Energy());
  masses[ 0] = getParticleData( 6)->mass();
  masses[ 1] = getParticleData( 5)->mass();
  masses[ 2] = getParticleData(15)->mass();
  masses[ 3] = getParticleData(ParticleID::SUSY_chi_1plus)->mass();
  masses[ 4] = getParticleData(ParticleID::SUSY_chi_2plus)->mass();
  masses[ 5] = _mw;
  masses[ 6] = getParticleData(ParticleID::Hplus)->mass();
  masses[ 7] = getParticleData(1000005)->mass();
  masses[ 8] = getParticleData(2000005)->mass();
  masses[ 9] = getParticleData(1000006)->mass();
  masses[10] = getParticleData(2000006)->mass();
  masses[11] = getParticleData(1000015)->mass();
  masses[12] = getParticleData(2000015)->mass();
  type.resize(13, PDT::Spin0);
  type[0] = PDT::Spin1Half;
  type[1] = PDT::Spin1Half;
  type[2] = PDT::Spin1Half;
  type[3] = PDT::Spin1Half;
  type[4] = PDT::Spin1Half;
  type[5] = PDT::Spin1;
  couplings.resize(13);
  VVSLoopVertex::doinit();
  if(loopToolsInitialized()) Looptools::ltexi();
}

void NMSSMPPHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << _sw << _cw << ounit(_mw, GeV) << ounit(_mz, GeV)
     << ounit(_lambdaVEV,GeV) << _lambda
     << ounit(_triTp,GeV) << ounit(_triBt,GeV) << ounit(_triTa,GeV)
     << _top << _bt << _tau << _mixS << _mixP << _mixU << _mixV
     << _mixQt << _mixQb << _mixLt << _sb << _cb <<  _kappa
     << ounit(_vu,GeV) << ounit(_vd,GeV) << ounit(_s,GeV) << ounit(_theAl,GeV);
}

void NMSSMPPHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> _sw >> _cw >> iunit(_mw, GeV) >> iunit(_mz, GeV)
     >> iunit(_lambdaVEV,GeV) >> _lambda
     >> iunit(_triTp,GeV) >> iunit(_triBt,GeV) >> iunit(_triTa,GeV)
     >> _top >> _bt >> _tau >> _mixS >> _mixP >> _mixU >> _mixV
     >> _mixQt >> _mixQb >> _mixLt >> _sb >> _cb >> _kappa 
     >> iunit(_vu,GeV) >> iunit(_vd,GeV) >> iunit(_s,GeV) >> iunit(_theAl,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMPPHVertex,VVSLoopVertex>
describeHerwigNMSSMPPHVertex("Herwig::NMSSMPPHVertex", "HwSusy.so HwNMSSM.so");

void NMSSMPPHVertex::Init() {

  static ClassDocumentation<NMSSMPPHVertex> documentation
    ("The effective coupling of a higgs to a pair of gluons in the "
     "NMSSM.");

}

void NMSSMPPHVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2,
				 tcPDPtr p3) {			 
  long hid(p3->id());
  double rt = sqrt(0.5);
  if( q2 != _q2last ) {
    Looptools::clearcache();
    _couplast = sqr(electroMagneticCoupling(q2));
    _coup = weakCoupling(q2);
    _q2last = q2;
    _recalc = true;
  }
  norm(_couplast*_coup);
  // scalar higgs bosons  
  if( hid != _hlast ) {
    _hlast = hid;
    _recalc = true;
    // top and bottom quark masses
    Energy mt   = _theSM->mass(q2, _top);
    Energy mb   = _theSM->mass(q2,  _bt);
    Energy mtau = _theSM->mass(q2, _tau);
    // scalar
    if( hid % 5 == 0 ) {
      // location of the higgs
      int iloc = (hid - 25)/10;
      // 6 particles in the loop
      setNParticles(13);
      Complex c(0.);
      // couplings for the top quark loop
      c = -1.5*sqr(_theSM->eu())*  mt*(*_mixS)(iloc, 1)/_sb/_mw;
      couplings[0] = make_pair(c,c);
      masses[0] = mt;
      // couplings for the bottom quark loop
      c = -1.5*sqr(_theSM->ed())*  mb*(*_mixS)(iloc, 0)/_cb/_mw;
      couplings[1] = make_pair(c,c);	
      masses[1] = mb;
      // couplings for the tau lepton loop
      c = -0.5*sqr(_theSM->ee())*mtau*(*_mixS)(iloc, 0)/_cb/_mw;
      couplings[2] = make_pair(c,c);	
      masses[2] = mtau;
      // charginos
      for(unsigned int ic=0;ic<2;++ic) {
	c = -_lambda/_coup*rt*(*_mixS)(iloc,2)*(*_mixU)(ic,1)*(*_mixV)(ic,1)
	  -rt*((*_mixS)(iloc,0)*(*_mixU)(ic,1)*(*_mixV)(ic,0) +
	       (*_mixS)(iloc,1)*(*_mixU)(ic,0)*(*_mixV)(ic,1));
	couplings[3+ic] = make_pair(c,c);
      }
      // W boson
      c = Complex(UnitRemoval::InvE*_mw*
		  (_cb*(*_mixS)(iloc,0)+_sb*(*_mixS)(iloc,1)));
      couplings[5] = make_pair(c,c);
      // charged Higgs
      complex<Energy> cpl;
      cpl = sqr(_lambda)*rt*2.*(_s*((*_mixS)(iloc,2)*sqr(_cb) + (*_mixS)(iloc,2)*sqr(_sb))
				- (_vu*(*_mixS)(iloc,0)/_coup + 
				   _vd*(*_mixS)(iloc,1)/_coup)*_sb*_cb)
	+_lambda*_sb*_cb*2*(*_mixS)(iloc,2)*(_kappa*_s/rt + rt*_theAl)
	+ sqr(_coup)*0.5*rt*sqr(_sw)/sqr(_cw)*((_vu*(*_mixS)(iloc,1)/_coup - 
						_vd*(*_mixS)(iloc,0)/_coup)*sqr(_cb) + 
					       (_vd*(*_mixS)(iloc,0)/_coup -
						_vu*(*_mixS)(iloc,1)/_coup)*sqr(_sb))
	+ sqr(_coup)*0.5*rt*(_vu*((*_mixS)(iloc,1)*sqr(_cb) + 
				  (*_mixS)(iloc,1)*sqr(_sb) + 
				  2.*(*_mixS)(iloc,0)*_cb*_sb)/_coup 
			     + _vd*((*_mixS)(iloc,0)*sqr(_cb) + 
				    (*_mixS)(iloc,0)*sqr(_sb) 
				    + 2.*(*_mixS)(iloc,1)*_sb*_cb)/_coup);
      cpl /= -_coup;
      couplings[6] = make_pair(Complex(cpl*UnitRemoval::InvE),
			       Complex(cpl*UnitRemoval::InvE));
      // sbottoms
      double f1 = mb/_mw/_cb;
      complex<Energy>  f2 = 0.5*_mz/_cw*
	( - _cb*(*_mixS)(iloc,0) + _sb*(*_mixS)(iloc,1));
      for(unsigned int ix=0;ix<2;++ix) {
	cpl = -f2*( (1. - 2.*sqr(_sw)/3.)*(*_mixQb)(ix, 0)*(*_mixQb)(ix, 0)
		    + 2.*sqr(_sw)*(*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)/3.)
	  - f1*mb*(*_mixS)(iloc,0)
	  *((*_mixQb)(ix, 0)*(*_mixQb)(ix, 0) + (*_mixQb)(ix, 1)*(*_mixQb)(ix, 1)) 
	  - 0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,1) - _lambda*_vu*(*_mixS)(iloc,2)/_coup 
		    +  _triBt*(*_mixS)(iloc,0))*((*_mixQb)(ix, 1)*(*_mixQb)(ix, 0)
						 + (*_mixQb)(ix, 0)*(*_mixQb)(ix, 1));
	cpl *= 3.*sqr(_theSM->ed());
	couplings[7+ix] = make_pair(Complex(cpl*UnitRemoval::InvE),Complex(cpl*UnitRemoval::InvE)); 
      }
      // stop
      f1 = mt/_mw/_sb;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl   =+f2*( (1. - 4.*sqr(_sw)/3.)*(*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
		     + 4.*sqr(_sw)*(*_mixQt)(ix, 1)*(*_mixQt)(ix, 1)/3.)
	  - f1*mt*(*_mixS)(iloc,1)
	  *((*_mixQt)(ix, 0)*(*_mixQt)(ix, 0)
	    + (*_mixQt)(ix, 1)*(*_mixQt)(ix, 1))  
	  -  0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,0) - _lambda*_vd*(*_mixS)(iloc,2)/_coup
		     + _triTp*(*_mixS)(iloc,1))*((*_mixQt)(ix, 1)*(*_mixQt)(ix, 0)
						 + (*_mixQt)(ix, 0)*(*_mixQt)(ix, 1));
	cpl *= 3.*sqr(_theSM->eu());
	couplings[9+ix] = make_pair(Complex(cpl*UnitRemoval::InvE),
				    Complex(cpl*UnitRemoval::InvE));
      } // sbottoms
      f1 = mtau/_mw/_cb;
      for(unsigned int ix=0;ix<2;++ix) {
	cpl = -f2*( (1. - 2.*sqr(_sw))*(*_mixLt)(ix, 0)*(*_mixLt)(ix, 0)
		    + 2.*sqr(_sw)*(*_mixLt)(ix, 1)*(*_mixLt)(ix, 1))
	  - f1*mtau*(*_mixS)(iloc,0)
	  *((*_mixLt)(ix, 0)*(*_mixLt)(ix, 0) + (*_mixLt)(ix, 1)*(*_mixLt)(ix, 1)) 
	  - 0.5*f1*(-_lambdaVEV*(*_mixS)(iloc,1) - _lambda*_vu*(*_mixS)(iloc,2)/_coup 
		    +  _triTa*(*_mixS)(iloc,0))*((*_mixLt)(ix, 1)*(*_mixLt)(ix, 0)
						 + (*_mixLt)(ix, 0)*(*_mixLt)(ix, 1));
	cpl *= sqr(_theSM->ee());
	couplings[11+ix] = make_pair(Complex(cpl*UnitRemoval::InvE),
				     Complex(cpl*UnitRemoval::InvE)); 
      }
    }
    // pseudoscalar higgs bosons	
    else {
      // location of the higgs
      int iloc = (hid - 36)/10;
      // 2 particles in the loop
      setNParticles(5);
      Complex c(0.);
      // top quark couplings
      c = Complex(0., 1.)*1.5*sqr(_theSM->eu())*  mt*(*_mixP)(iloc, 1)/_sb/_mw;
      couplings[0] = make_pair(c,-c);
      masses[0] = mt;
      // bottom quark couplings
      c = Complex(0., 1.)*1.5*sqr(_theSM->ed())*  mb*(*_mixP)(iloc, 0)/_cb/_mw;
      couplings[1] = make_pair(c,-c);	
      masses[1] = mb;
      // tau lepton couplings
      c = Complex(0., 1.)*0.5*sqr(_theSM->ee())*mtau*(*_mixP)(iloc, 0)/_cb/_mw;
      couplings[2] = make_pair(c,-c);	
      masses[2] = mtau;
      // charginos
      for(unsigned int ic=0;ic<2;++ic) {
	c = Complex(0,-1.0)*
	  (_lambda/_coup*rt*(*_mixP)(iloc,2)*(*_mixU)(ic,1)*(*_mixV)(ic,1)
	   -rt*((*_mixP)(iloc,0)*(*_mixU)(ic,1)*(*_mixV)(ic,0)
		     + (*_mixP)(iloc,1)*(*_mixU)(ic,0)*(*_mixV)(ic,1)));
			  couplings[3+ic] = make_pair(-c,c); 

      }
    }
  }

  if( _recalc ) {
    VVSLoopVertex::setCoupling(q2, p1, p2, p3);
    _recalc = false;
  } 
}
#line 1 "./NMSSMWWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NMSSMWWHHVertex class.
//

#include "NMSSMWWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "NMSSM.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NMSSMWWHHVertex::NMSSMWWHHVertex() : couplast_(0.),q2last_(ZERO),
				     sw_(0.), cw_(0.), sb_(0.), cb_(0.) {
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr NMSSMWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr NMSSMWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void NMSSMWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << sw_ << cw_ << sb_ << cb_ << mixS_ << mixP_;
}

void NMSSMWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> sw_ >> cw_ >> sb_ >> cb_ >> mixS_ >> mixP_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NMSSMWWHHVertex,Helicity::VVSSVertex>
describeHerwigNMSSMWWHHVertex("Herwig::NMSSMWWHHVertex", "NMSSMWWHHVertex.so");

void NMSSMWWHHVertex::Init() {

  static ClassDocumentation<NMSSMWWHHVertex> documentation
    ("The NMSSMWWHHVertex class implements the coupling of"
     " two electroweak and two Higgs bosons in the NMSSM.");

}

void NMSSMWWHHVertex::doinit() {
  int scalar[3]={25,35,45};
  int pseudo[2]={36,46};
  // scalar higgs bosons
  for(unsigned int i=0;i<3;++i) {
    // pair of scalars
    for(unsigned int j=0;j<3;++j) {
      addToList( 24,-24,scalar[i],scalar[j]);
      addToList( 23, 23,scalar[i],scalar[j]);
    }
    // scalar charged
    addToList( 22, 24,scalar[i],-37);
    addToList( 22,-24,scalar[i], 37);
    addToList( 23, 24,scalar[i],-37);
    addToList( 23,-24,scalar[i], 37);
  }
  // pair of pseudoscalars
  for(unsigned int i=0;i<2;++i) {
    for(unsigned int j=0;j<2;++j) {
      addToList( 24,-24,pseudo[i],pseudo[j]);
      addToList( 23, 23,pseudo[i],pseudo[j]);
    }
    // pseudo charged
    addToList( 22, 24,pseudo[i],-37);
    addToList( 22,-24,pseudo[i], 37);
    addToList( 23, 24,pseudo[i],-37);
    addToList( 23,-24,pseudo[i], 37);
  }
  addToList( 24,-24, 37,-37);
  addToList( 23, 23, 37,-37);
  addToList( 22, 23, 37,-37);
  addToList( 22, 22, 37,-37);
  // cast to NMSSM model
  tcNMSSMPtr model=dynamic_ptr_cast<tcNMSSMPtr>(generator()->standardModel());
  if(!model) 
    throw InitException() << "Must have the NMSSM Model in NMSSMWWHHVertex::doinit()"
			  << Exception::runerror;
  // sin and cos theta_W
  sw_ = sqrt(   sin2ThetaW());
  cw_ = sqrt(1.-sin2ThetaW());
  // get the mixing matrices
  mixS_ = model->CPevenHiggsMix();
  if(!mixS_) throw InitException() << "Mixing matrix for CP-even neutral Higgs"
				   << " bosons is not set in NMSSMWWHHVertex::doinit()" 
				   << Exception::runerror;
  mixP_ = model->CPoddHiggsMix();
  if(!mixP_) throw InitException() << "Mixing matrix for CP-odd neutral Higgs"
				   << " bosons is not set in NMSSMWWHHVertex::doinit()" 
				   << Exception::runerror;
  // sin and cos beta
  double beta = atan(model->tanBeta());
  sb_ = sin(beta);
  cb_ = cos(beta);
  // base class
  VVSSVertex::doinit();
}
 
void NMSSMWWHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,
				  tcPDPtr part3, tcPDPtr part4) {
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = sqr(electroMagneticCoupling(q2));
    q2last_=q2;
  }
  int ibos1 = part1->id(), ibos2 = part2->id();
  int isca1 = part3->id(), isca2 = part4->id();
  Complex fact(1.);
  if(abs(ibos1)==abs(ibos2)) {
    fact *= 0.5/sqr(sw_);
    if(ibos1==ParticleID::Z0) fact /= sqr(cw_);
    // charged
    if(abs(isca1)==37) {
      if(ibos1==ParticleID::Z0) fact *= sqr(sqr(cw_)-sqr(sw_));
      else if(ibos1==ParticleID::gamma) fact = 2.;
    }
    // pair of scalars
    else if((isca1-5)%10==0) {
      unsigned int i = (isca1-25)/10;
      unsigned int j = (isca2-25)/10;
      fact *= (*mixS_)(i,0)*(*mixS_)(j,0)+(*mixS_)(i,1)*(*mixS_)(j,1);
    }
    // pair of pseudoscalars
    else if((isca1-6)%10==0) {
      unsigned int i = (isca1-36)/10;
      unsigned int j = (isca2-36)/10;
      fact *= (*mixP_)(i,0)*(*mixP_)(j,0)+(*mixP_)(i,1)*(*mixP_)(j,1);
    }
    else
      assert(false);
  }
  else if(abs(ibos1)==ParticleID::Wplus ||
	  abs(ibos2)==ParticleID::Wplus) {
    if(abs(ibos1)==ParticleID::gamma ||
       abs(ibos2)==ParticleID::gamma) {
      fact *= -0.5/sw_;
    }
    else {
      fact *=  0.5/cw_;
    }
    if((isca1-5)%10==0) {
      unsigned int i = (isca1-25)/10;
      fact *= sb_*(*mixS_)(i,0) - cb_*(*mixS_)(i,1);
    }
    else if((isca2-5)%10==0) {
      unsigned int i = (isca2-25)/10;
      fact *= sb_*(*mixS_)(i,0) - cb_*(*mixS_)(i,1);
    }
    else if((isca1-6)%10==0) {
      unsigned int i = (isca1-36)/10;
      fact *= sb_*(*mixP_)(i,0) + cb_*(*mixP_)(i,1);
      fact *= isca2==ParticleID::Hplus ? -Complex(0.,1.) : Complex(0.,1.);
    }
    else if((isca2-6)%10==0) {
      unsigned int i = (isca2-36)/10;
      fact *= sb_*(*mixP_)(i,0) + cb_*(*mixP_)(i,1);
      fact *= isca1==ParticleID::Hplus ? -Complex(0.,1.) : Complex(0.,1.);
    }
    else
      assert(false);
  }
  else {
    fact = (sqr(cw_)-sqr(sw_))/cw_/sw_;
  }
  // set the coupling
  norm(couplast_*fact);
}
