#line 1 "./StandardModel.cc"
// -*- C++ -*-
//
// StandardModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModel class.
//

#include "StandardModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/General/ModelGenerator.h"
#include "ThePEG/Repository/BaseRepository.h"

using namespace Herwig;

StandardModel::StandardModel() {}

StandardModel::~StandardModel() {}

StandardModel::StandardModel(const StandardModel & x)
  : StandardModelBase(x), 
    FFZVertex_ (x.FFZVertex_),
    FFPVertex_ (x.FFPVertex_) , FFGVertex_ (x.FFGVertex_) ,
    FFWVertex_ (x.FFWVertex_) , FFHVertex_ (x.FFHVertex_) ,
    WWHVertex_ (x.WWHVertex_) ,
    GGGVertex_ (x.GGGVertex_) ,
    WWWVertex_ (x.WWWVertex_) , GGGGVertex_(x.GGGGVertex_),
    WWWWVertex_(x.WWWWVertex_), HGGVertex_ (x.HGGVertex_) ,
    HPPVertex_ (x.HPPVertex_) , HHHVertex_ (x.HHHVertex_) ,
    WWHHVertex_ (x.WWHHVertex_) ,
    vertexList_(x.vertexList_), extraVertices_(x.extraVertices_),
    runningMass_(x.runningMass_),modelGenerator_(x.modelGenerator_),
    couplings_(x.couplings_)
{}

IBPtr StandardModel::clone() const {
  return new_ptr(*this);
}

IBPtr StandardModel::fullclone() const {
  return new_ptr(*this);
}

void StandardModel::doinit() {
  if(runningMass_) {
    runningMass_->init();
  }
  //add Standard Model vertices
  if ( registerDefaultVertices() ) {
    addVertex(FFZVertex_);
    addVertex(FFPVertex_);
    addVertex(FFGVertex_);
    addVertex(FFWVertex_);
    addVertex(vertexFFH());
    addVertex(vertexWWH());
    addVertex(GGGVertex_);
    addVertex(WWWVertex_);
    addVertex(GGGGVertex_);
    addVertex(WWWWVertex_);
    addVertex(vertexHGG());
    addVertex(HPPVertex_);
    if(HHHVertex_ ) addVertex(HHHVertex_);
    if(WWHHVertex_) addVertex(WWHHVertex_);
  }
  if(couplings_.find("QED")==couplings_.end()) {
    couplings_["QED"] = make_pair(1,99);
  }
  if(couplings_.find("QCD")==couplings_.end()) {
    couplings_["QCD"] = make_pair(2,99);
  }
  StandardModelBase::doinit();
}

void StandardModel::persistentOutput(PersistentOStream & os) const {
  os << FFZVertex_ <<FFPVertex_ << FFGVertex_ << FFWVertex_ 
     << FFHVertex_ << WWHVertex_ << GGGGVertex_ << WWWWVertex_
     << GGGVertex_ << WWWVertex_  << HGGVertex_  << HPPVertex_ 
     << HHHVertex_ << WWHHVertex_ 
     << runningMass_ << vertexList_ << extraVertices_ << modelGenerator_;
}

void StandardModel::persistentInput(PersistentIStream & is, int) {
  is >> FFZVertex_ >> FFPVertex_ >> FFGVertex_ >> FFWVertex_
     >> FFHVertex_ >> WWHVertex_ >> GGGGVertex_ >> WWWWVertex_
     >> GGGVertex_ >> WWWVertex_ >> HGGVertex_  >> HPPVertex_ 
     >> HHHVertex_ >> WWHHVertex_ 
     >> runningMass_ >> vertexList_ >> extraVertices_ >> modelGenerator_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StandardModel,StandardModelBase>
describeHerwigStandardModel("Herwig::StandardModel", "Herwig.so");

void StandardModel::Init() {

  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFZ
    ("Vertex/FFZ",
     "Reference to the Standard Model FFZ Vertex",
     &StandardModel::FFZVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFP
    ("Vertex/FFP",
     "Reference to the Standard Model FFP Vertex",
     &StandardModel::FFPVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFG
    ("Vertex/FFG",
     "Reference to the Standard Model FFG Vertex",
     &StandardModel::FFGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFVVertex> interfaceVertexFFW
    ("Vertex/FFW",
     "Reference to the Standard Model FFW Vertex",
     &StandardModel::FFWVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractFFSVertex> interfaceVertexFFH
    ("Vertex/FFH",
     "Reference to the Standard Model FFH Vertex.",
     &StandardModel::FFHVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVertex> interfaceVertexGGG
    ("Vertex/GGG",
     "Reference to the Standard Model GGG Vertex",
     &StandardModel::GGGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVertex> interfaceVertexWWW
    ("Vertex/WWW",
     "Reference to the Standard Model WWW Vertex",
     &StandardModel::WWWVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVSVertex> interfaceVertexWWH
    ("Vertex/WWH",
     "Reference to the Standard Model WWH Vertex",
     &StandardModel::WWHVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVVertex> interfaceVertexWWWW
    ("Vertex/WWWW",
     "Reference to the Standard Model WWWW Vertex",
     &StandardModel::WWWWVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVVVVertex> interfaceVertexGGGG
    ("Vertex/GGGG",
     "Reference to the Standard Model GGGG Vertex",
     &StandardModel::GGGGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVSVertex> interfaceVertexHGG
    ("Vertex/HGG",
     "Reference to the StandardModel HGG Vertex",
     &StandardModel::HGGVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractVVSVertex> interfaceVertexHPP
    ("Vertex/HPP",
     "Reference to StandardModel HPPVertex",
     &StandardModel::HPPVertex_, false, false, true, false);
  
  static Reference<StandardModel,AbstractSSSVertex> interfaceVertexHHH
    ("Vertex/HHH",
     "Reference to the Standard Model HHHVertex",
     &StandardModel::HHHVertex_, false, false, true, true);

  static Reference<StandardModel,AbstractVVSSVertex> interfaceVertexWWHH
    ("Vertex/WWHH",
     "Reference to the Standard Model WWHHVertex",
     &StandardModel::WWHHVertex_, false, false, true, true);

  static RefVector<StandardModel,VertexBase> interfaceExtraVertices
    ("ExtraVertices",
     "Additional vertices to be considered in automatic ME construction.",
     &StandardModel::extraVertices_, -1, true, false, true, false, false);

  static Reference<StandardModel,RunningMassBase> interfaceRunningMass
    ("RunningMass",
     "Reference to the running mass object",
     &StandardModel::runningMass_, false, false, true, false);
  
  static Reference<StandardModel,Herwig::ModelGenerator> interfaceModelGenerator
    ("ModelGenerator",
     "Pointer to ModelGenerator class",
     &StandardModel::modelGenerator_, false, false, true, true);

  static ClassDocumentation<StandardModel> documentation
    ("The StandardModel class inherits from StandardModelBase"
     "and supplies additional couplings and access to the StandardModel"
     "vertices for helicity amplitude calculations" );

}

void StandardModel::resetMass(long id, Energy mass,tPDPtr part) {
  if(!part) part = getParticleData(id);
  if(!part) return;
  const InterfaceBase * ifb = BaseRepository::FindInterface(part, "NominalMass");
  ostringstream os;
  os << setprecision(12) << abs(mass/GeV);
  ifb->exec(*part, "set", os.str());
}
#line 1 "./RunningMassBase.cc"
// -*- C++ -*-
//
// RunningMassBase.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RunningMassBase class.
//

#include "RunningMassBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void RunningMassBase::persistentOutput(PersistentOStream & os) const {
  os << ounit(_theMass, GeV);
}

void RunningMassBase::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_theMass, GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<RunningMassBase,Interfaced>
describeHerwigRunningMassBase("Herwig::RunningMassBase", "Herwig.so");

void RunningMassBase::Init() {
 
  static ClassDocumentation<RunningMassBase> documentation
    ("The RunningMassBase class is the base class for running mass"
     "calculations");
  
}

void RunningMassBase::doinit() {
  _theMass = mass();
  Interfaced::doinit();
}
#line 1 "./RunningMass.cc"
// -*- C++ -*-
//
// RunningMass.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RunningMass class.
//

#include "RunningMass.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"

namespace Herwig {
using namespace ThePEG;

void RunningMass::persistentOutput(PersistentOStream & os) const {
  os << _theQCDOrder << _thePower << _theCoefficient << _theMaxFlav
     << _theStandardModel << _lightOption << _heavyOption;
}

void RunningMass::persistentInput(PersistentIStream & is, int) {
  is >> _theQCDOrder >> _thePower >> _theCoefficient >> _theMaxFlav 
     >> _theStandardModel >> _lightOption >> _heavyOption;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RunningMass,RunningMassBase>
describeHerwigRunningMass("Herwig::RunningMass", "Herwig.so");

void RunningMass::Init() {

  static Parameter<RunningMass,unsigned int> interfaceQCDOrder
    ("QCDOrder",
     "The order in alpha_S",
     &RunningMass::_theQCDOrder, 1, 1, 2,
     false, false, true);
    
  static Parameter<RunningMass,unsigned int> interfaceMaxFlav
    ("MaxFlav",
     "maximum number of flavours",
     &RunningMass::_theMaxFlav, 6, 3, 6,
     false, false, true);
  
  static ClassDocumentation<RunningMass> documentation
    ("The RunningMass class is the implementation of the"
     " QCD running mass to one or two loop in alpha_S");

  static Switch<RunningMass,unsigned int> interfaceLightQuarkMass
    ("LightQuarkMass",
     "Option for the treatment of light mass masses",
     &RunningMass::_lightOption, 1, false, false);
  static SwitchOption interfaceLightQuarkMassRunning
    (interfaceLightQuarkMass,
     "Running",
     "Use a running, probably zero mass",
     0);
  static SwitchOption interfaceLightQuarkMassPole
    (interfaceLightQuarkMass,
     "Pole",
     "Use the pole mass",
     1);

  static Switch<RunningMass,unsigned int> interfaceBottomCharmMass
    ("TopBottomCharmMass",
     "Option for using a running or pole mass for the top, bottom and charm quarks",
     &RunningMass::_heavyOption, 0, false, false);
  static SwitchOption interfaceBottomCharmMassRunning
    (interfaceBottomCharmMass,
     "Running",
     "Use the running mass",
     0);
  static SwitchOption interfaceBottomCharmMassPole
    (interfaceBottomCharmMass,
     "Pole",
     "Use the pole mass",
     1);

}
// Return the masses used.
vector<Energy> RunningMass::mass() const {
  using Constants::pi;
  vector<Energy> masses;
  for ( long f = 1; f <= long(_theMaxFlav); ++f ) {
    PDPtr p = getParticleData(f);
    Energy massf = p ? p->mass() : ZERO;
    if ( !_theStandardModel->alphaSPtr()->quarkMasses().empty() ) {
      if ( f < 3 ) {
	massf = 
	  f == 1 ?
	  _theStandardModel->alphaSPtr()->quarkMasses()[1] :
	  _theStandardModel->alphaSPtr()->quarkMasses()[0];
      } else {
	massf = _theStandardModel->alphaSPtr()->quarkMasses()[f-1];
      }
    }
    if((f<=3&&_lightOption==0) ||
       (f>3 &&_heavyOption==0)) {
      double coeff = _theQCDOrder==2 ? _theCoefficient[f-1]+4./3./pi : 0.;
      double as = _theStandardModel->alphaS(massf*massf);
      massf = as>0 ? massf/(1.+coeff*as)/pow(as,_thePower[f-1]) : ZERO;
    }
    masses.push_back(massf);
  }
  return masses;
}

// Return the running mass for a given scale and particle type.
Energy RunningMass::value(Energy2 scale, tcPDPtr part) const {
  Energy output;
  unsigned int id=abs(part->id());
  // calculate the running mass
  if(id<=_theMaxFlav) {
    if( (id <= 3 && _lightOption == 1 ) ||
	(id >= 4 && _heavyOption == 1 ) ) {
      output= part->mass();
    }
    else {
      // calculate the value of alphaS and number of flavours
      //unsigned int nf=  _theStandardModel->Nf(scale);
      unsigned int nf=id;
      double as = _theStandardModel->alphaS(scale);
      id=id-1;
      output = massElement(id)*pow(as,_thePower[nf-1]);
      if(_theQCDOrder==2){output*=(1.+as*_theCoefficient[nf-1]);}
    }
  }
  // 
  else {
    output= part->mass();
  }
  return output;
}

void RunningMass::doinit() {
  using Constants::pi;
  _theStandardModel = generator()->standardModel();
  _theStandardModel->alphaSPtr()->init();
  // coefficients for the calculation
  double c = 1./pi,cprime,b,bprime,power,coeff;
  for(unsigned int f=1;f<=_theMaxFlav;++f) {
    // the basic parameters for the running mass
    cprime =     c*(303.-10.*f)/72.;
    b      =     c*(33. -2. *f)/12.;
    bprime = 0.5*c*(153.-19.*f)/(33.-2.*f);
    power = c/b;
    coeff = c*(cprime-bprime)/b;
    _thePower.push_back(power);
    _theCoefficient.push_back(coeff);
  }
  RunningMassBase::doinit();
}
}
#line 1 "./StandardCKM.cc"
// -*- C++ -*-
//
// StandardCKM.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardCKM class.
//

#include "StandardCKM.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"


using namespace Herwig;
using namespace ThePEG;

IBPtr StandardCKM::clone() const {
  return new_ptr(*this);
}

IBPtr StandardCKM::fullclone() const {
  return new_ptr(*this);
}

vector< vector<double> > StandardCKM::getMatrix(unsigned int nFamilies) const {
  vector< vector<double> > ckm(nFamilies, vector<double>(nFamilies, 0.0));
  for ( unsigned int i = 0; i < nFamilies; ++i ) ckm[i][i] = 1.0;
  if ( nFamilies <= 1 ) return ckm;
  double s12 = sin(theta12);
  double c12 = cos(theta12);
  if ( nFamilies == 2 ) {
    ckm[0][0] = sqr(c12);
    ckm[0][1] = sqr(s12);
    ckm[1][0] = sqr(s12);
    ckm[1][1] = sqr(c12);
    return ckm;
  }
  double s13 = sin(theta13);
  double c13 = cos(theta13);
  double s23 = sin(theta23);
  double c23 = cos(theta23);
  double cd = cos(delta);
  ckm[0][0] = sqr(c12*c13);
  ckm[0][1] = sqr(s12*c13);
  ckm[0][2] = sqr(s13);
  ckm[1][0] = sqr(s12*c23)+sqr(c12*s23*s13)+2.0*s12*c23*c12*s23*s13*cd;
  ckm[1][1] = sqr(c12*c23)+sqr(s12*s23*s13)-2.0*c12*c23*s12*s23*s13*cd;
  ckm[1][2] = sqr(s23*c13);
  ckm[2][0] = sqr(s12*s23)+sqr(c12*c23*s13)-2.0*s12*s23*c12*c23*s13*cd;
  ckm[2][1] = sqr(c12*s23)+sqr(s12*c23*s13)+2.0*c12*s23*s12*c23*s13*cd;
  ckm[2][2] = sqr(c23*c13);
  return ckm;
}
vector< vector<complex<double> > > 
StandardCKM::getUnsquaredMatrix(unsigned int nFamilies) const {
  vector< vector<complex<double> > > ckm(nFamilies, vector<complex<double> >(nFamilies, 0.0));
  for ( unsigned int i = 0; i < nFamilies; ++i ) ckm[i][i] = 1.0;
  if ( nFamilies <= 1 ) return ckm;
  double s12 = sin(theta12);
  double c12 = cos(theta12);
  if ( nFamilies == 2 ) {
    ckm[0][0] = sqr(c12);
    ckm[0][1] = sqr(s12);
    ckm[1][0] = sqr(s12);
    ckm[1][1] = sqr(c12);
    return ckm;
  }
  double s13 = sin(theta13);
  double c13 = cos(theta13);
  double s23 = sin(theta23);
  double c23 = cos(theta23);
  double cd = cos(delta);
  double sd = sin(delta);
  complex<double> ii(0.,1.);
  complex<double> expid  = cd+ii*sd;
  complex<double> expmid = cd-ii*sd;
  ckm[0][0] =  c12*c13;
  ckm[0][1] =  s12*c13;
  ckm[0][2] =  s13*expmid;
  ckm[1][0] = -s12*c23-c12*s23*s13*expid;
  ckm[1][1] =  c12*c23-s12*s23*s13*expid;
  ckm[1][2] =  s23*c13;
  ckm[2][0] =  s12*s23-c12*c23*s13*expid;
  ckm[2][1] = -c12*s23-s12*c23*s13*expid;
  ckm[2][2] =  c23*c13;
  return ckm;
}

void StandardCKM::persistentOutput(PersistentOStream & os) const {
  os << theta12 << theta13 << theta23 << delta;
}

void StandardCKM::persistentInput(PersistentIStream & is, int) {
  is >> theta12 >> theta13 >> theta23 >> delta;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StandardCKM,CKMBase>
describeHerwigStandardCKM("Herwig::StandardCKM", "Herwig.so");

void StandardCKM::Init() {
  
  static Parameter<StandardCKM,double> interfaceTheta12
    ("theta_12",
     "The mixing angle between the first and second generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta12, 0.2262, 0.0, Constants::twopi, false, false, true);
  static Parameter<StandardCKM,double> interfaceTheta13
    ("theta_13",
     "The mixing angle between the first and third generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta13, 0.0037, 0.0, Constants::twopi, false, false, true);
  static Parameter<StandardCKM,double> interfaceTheta23
    ("theta_23",
     "The mixing angle between the second and third generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta23, 0.0413, 0.0, Constants::twopi, false, false, true);
  static Parameter<StandardCKM,double> interfaceDelta
    ("delta",
     "The phase angle in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::delta, 1.05, 0.0, Constants::twopi, false, false, true);
}
#line 1 "./O2AlphaS.cc"
// -*- C++ -*-
//
// O2AlphaS.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the O2AlphaS class.
//

#include "O2AlphaS.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void O2AlphaS::persistentOutput(PersistentOStream & os) const {
  os << ounit(_lambdaQCD,GeV) << _bcoeff << _ccoeff << ounit(_lambdas,GeV) 
     << ounit(_threshold,GeV) << _match << _copt;
}

void O2AlphaS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_lambdaQCD,GeV) >> _bcoeff >> _ccoeff >> iunit(_lambdas,GeV) 
     >> iunit(_threshold,GeV) >> _match >> _copt;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<O2AlphaS,AlphaSBase>
describeHerwigO2AlphaS("Herwig::O2AlphaS", "Herwig.so");

void O2AlphaS::Init() {

  static ClassDocumentation<O2AlphaS> documentation
    ("The O2AlphaS class implements the next-to-leading order alphaS in the same"
     " way as in FORTRAN HERWIG");

  static Parameter<O2AlphaS,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "The value of Lambda QCD",
     &O2AlphaS::_lambdaQCD, MeV, 180.*MeV, 50.*MeV, 500.0*MeV,
     false, false, Interface::limited);


  static Switch<O2AlphaS,unsigned int> interfaceLambdaType
    ("LambdaType",
     "Which type of Lambda to use",
     &O2AlphaS::_copt, 0, false, false);
  static SwitchOption interfaceLambdaTypeMonteCarlo
    (interfaceLambdaType,
     "MonteCarlo",
     "Use a Monte Carlo scheme as in the FORTRAN program",
     0);
  static SwitchOption interfaceLambdaTypeMSBar
    (interfaceLambdaType,
     "MSBar",
     "Use the MSBar scheme",
     1);
}

vector<Energy2> O2AlphaS::flavourThresholds() const {
  vector<Energy2> thresholds(_threshold.size());
  transform(_threshold.begin(), _threshold.end(),
	    thresholds.begin(),
	    sqr<Energy>);
  return thresholds;
}

void O2AlphaS::doinit() {
  // thresholds
  for ( int ix=1; ix<7; ++ix ) {
    if ( quarkMasses().empty() ) {
      tPDPtr p = getParticleData(ix);
      _threshold[ix-1] = p->mass();
    } else
      _threshold[ix-1] = quarkMasses()[ix-1];
  }
  // d is heavier than u, need to swap
  swap(_threshold[0],_threshold[1]);

  // beta function coefficients
  const double ca = generator()->standardModel()->Nc();
  const double cf = (sqr(ca)-1.)/2./ca;
  for(unsigned int ix=3;ix<7;++ix)
    {
      _bcoeff[ix-1]=(11.*ca-2.*ix)/(12.*Constants::pi);
      _ccoeff[ix-1]=(17.*sqr(ca)-ix*(5.*ca+3.*cf))/(24.*sqr(Constants::pi))/sqr(_bcoeff[ix-1]);
    }
  if(_copt==0)
    {
      double kfac(ca*(67./18.-sqr(Constants::pi)/6.)-25./9.);
      _lambdas[5]=_lambdaQCD*exp(kfac/(4.*Constants::pi*_bcoeff[4]))/sqrt(2.);
    }
  else{_lambdas[5]=_lambdaQCD;}
  // calculate the threshold matching
  double rho=2.*log(_threshold[5]/_lambdas[5]);
  double rat=log(rho)/rho;
  _match[5]=(_bcoeff[4]/(1.-_ccoeff[4]*rat)-_bcoeff[5]/(1.-_ccoeff[5]*rat))*rho;
  rho=2.*log(_threshold[4]/_lambdas[5]);
  rat=log(rho)/rho;
  _match[4]=(_bcoeff[4]/(1.-_ccoeff[4]*rat)-_bcoeff[3]/(1.-_ccoeff[3]*rat))*rho;
  rho=2.*log(_threshold[3]/_lambdas[5]);
  rat=log(rho)/rho;
  _match[3]=(_bcoeff[3]/(1.-_ccoeff[3]*rat)-_bcoeff[2]/(1.-_ccoeff[2]*rat))*rho
    +_match[4];
  // calculate the 4-flavour lambda
  _lambdas[4]=_lambdas[5]*pow(_threshold[4]/_lambdas[5],2./25.)*
    pow(2.*log(_threshold[4]/_lambdas[5]),963./14375.);
  // calculate the 3-flavour lambda
  double eps(1.e-6),d35(-1./(_bcoeff[2]*_match[3])),rlf,drh;
  unsigned int ix=0;
  do
    {
      rat=log(d35)/d35;
      rlf=_bcoeff[2]*d35/(1.-_ccoeff[2]*rat);
      drh=_bcoeff[2]*(rlf+_match[3])*sqr(d35)/
	((1.-2.*_ccoeff[2]*rat+_ccoeff[2]/d35)*sqr(rlf));
      d35=d35-drh;
      ++ix;
    }
  while(ix<100&&abs(drh)>eps*d35);
  _lambdas[3]=_lambdas[5]*exp(0.5*d35);
  AlphaSBase::doinit();
}

vector<Energy> O2AlphaS::LambdaQCDs() const
{
  vector<Energy> output(4,_lambdas[3]);
  output.push_back(_lambdas[4]);
  output.push_back(_lambdas[5]);
  output.push_back(_lambdas[5]);
  return output;
}



double O2AlphaS::value(Energy2 scale, const StandardModelBase &) const
{
  Energy rs=sqrt(scale);
  if(scale<sqr(_lambdas[5])) {
    generator()->logWarning(Exception()
			    << "O2AlphaS called with scale less than Lambda QCD "
			    << "scale = " << rs/MeV << " MeV and "
			    << "Lambda = " << _lambdas[5]/MeV << " MeV"
			    << Exception::warning);
    return 0.;
  }
  double rho=2.*log(rs/_lambdas[5]),rat(log(rho)/rho);
  double rlf;
  if(rs>_threshold[5])      rlf=_bcoeff[5]*rho/(1.-_ccoeff[5]*rat)+_match[5];
  else if(rs>_threshold[4]) rlf=_bcoeff[4]*rho/(1.-_ccoeff[4]*rat);
  else if(rs>_threshold[3]) rlf=_bcoeff[3]*rho/(1.-_ccoeff[3]*rat)+_match[4];
  else                      rlf=_bcoeff[2]*rho/(1.-_ccoeff[2]*rat)+_match[3];
  // must be possible
  if(rlf<=0.) {
    generator()->logWarning(Exception() << "O2AlphaS coupling less than zero"
			    << Exception::warning) ;
    return 0.;
  }
  return 1./rlf;
}


#line 1 "./AlphaEM.cc"
// -*- C++ -*-
//
// AlphaEM.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AlphaEM class.
//

#include "AlphaEM.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;

IBPtr AlphaEM::clone() const {
  return new_ptr(*this);
}

IBPtr AlphaEM::fullclone() const {
  return new_ptr(*this);
}

void AlphaEM::persistentOutput(PersistentOStream & os) const {
  os << ounit(_me,GeV2) << ounit(_mmu,GeV2) 
     << ounit(_mtau,GeV2) << ounit(_mtop,GeV2);
}

void AlphaEM::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_me,GeV2) >> iunit(_mmu,GeV2) 
     >> iunit(_mtau,GeV2) >> iunit(_mtop,GeV2);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<AlphaEM,AlphaEMBase>
describeHerwigAlphaEM("Herwig::AlphaEM", "Herwig.so");

void AlphaEM::Init() {

  static ClassDocumentation<AlphaEM> documentation
    ("This class implements a running \\f$\\alpha_{\\mbox{EM}}\\f$ according "
     "to Buckhardt et al.",
     "In the running of $\\alpha_{EM}$, the parametrization of "
     "H.~Buckhardt et al. was used. See \\cite{KLEISSCERN9808v3pp129}.",
     "\\bibitem{KLEISSCERN9808v3pp129} R.~Kleiss et al, "
     "CERN 89-08, vol.~3, pp 129-131.");

}

void AlphaEM::doinit() {
  AlphaEMBase::doinit();
  _me   = sqr(getParticleData(ParticleID::eminus)->mass());
  _mmu  = sqr(getParticleData(ParticleID::muminus)->mass());
  _mtau = sqr(getParticleData(ParticleID::tauminus)->mass());
  _mtop = sqr(getParticleData(ParticleID::t)->mass());
}

double AlphaEM::realPi(double r) const {
  static double fvthr=1.666666666666667e0,rmax=1.e6;
  // use assymptotic formula
  if(abs(r)<1e-3) {
    return -fvthr-log(r);
  }
  // return zero for large values
  else if(abs(r)>rmax) {
    return 0.;
  }
  else if(4.*r>1.) {
    double beta=sqrt(4.*r-1.);
    return 1./3.
      -(1.+2.*r)*(2.-beta*acos(1.-1./(2.*r)));
  }
  else {
    double beta=sqrt(1.-4.*r);
    return 1./3.
      -(1.+2.*r)*(2.+beta*log(abs((beta-1.)/(beta+1.))));
  }
}

double AlphaEM::value(Energy2 scale, const StandardModelBase & sm) const {
  useMe();
  static double eps=1e-6;
  static double a1=0.0    ,b1=0.00835,c1=1.000;
  static double a2=0.0    ,b2=0.00238,c2=3.927;
  static double a3=0.00165,b3=0.00299,c3=1.000;
  static double a4=0.00221,b4=0.00293,c4=1.000;
  // alpha_EM at Q^2=0
  double alem=sm.alphaEM();
  double aempi=alem/(3.*Constants::pi);
  // convert scale to GeV2
  double Q2=scale/GeV2;
  double x=abs(Q2);
  // return q^2=0 value for small scales
  if(x<eps) return alem;
  // leptonic component
  double repigg=aempi*(realPi(_me/scale)+realPi(_mmu/scale)+realPi(_mtau/scale));
  // Hadronic component from light quarks
  if(x<9e-2)      repigg+=a1+b1*log(1.+c1*x);
  else if(x<9.)   repigg+=a2+b2*log(1.+c2*x);
  else if(x<1.e4) repigg+=a3+b3*log(1.+c3*x);
  else            repigg+=a4+b4*log(1.+c4*x);
  // Top Contribution
  repigg+=aempi*realPi(_mtop/scale);
  // reutrn the answer
  return alem/(1.-repigg);
}
#line 1 "./SMFFGVertex.cc"
// -*- C++ -*-
//
// SMFFGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelFFGVertex class.
//

#include "SMFFGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Herwig;
using namespace ThePEG;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SMFFGVertex,FFVVertex>
describeHerwigSMFFGVertex("Herwig::SMFFGVertex", "Herwig.so");

void SMFFGVertex::Init() {

  static ClassDocumentation<SMFFGVertex> documentation
    ("The SMFFGVertex class is the implementation of"
     "the coupling of the gluon to the Standard Model fermions");
  
}

// coupling for FFG vertex
#ifndef NDEBUG
void SMFFGVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr bb,tcPDPtr) {
#else
void SMFFGVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,tcPDPtr) {
#endif
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -strongCoupling(q2);
    _q2last=q2;
  }
  // the left and right couplings
  assert( abs(aa->id()) >= 1 && abs(aa->id()) <= 6 );
  assert( aa->id()==-bb->id());
  if(aa->id()<0)
    norm( _couplast);
  else
    norm(-_couplast);
  left(1.);
  right(1.);
}

SMFFGVertex::SMFFGVertex() : _couplast(0.), _q2last(ZERO) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3TFUND);
}
  
void SMFFGVertex::doinit() {
  // PDG codes for the particles
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21);
  }
  FFVVertex::doinit();
}
#line 1 "./SMFFHVertex.cc"
// -*- C++ -*-
//
// SMFFHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFHVertex class.
//

#include "SMFFHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/Constants.h"

using namespace Herwig;

IBPtr SMFFHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMFFHVertex::fullclone() const {
  return new_ptr(*this);
}


SMFFHVertex::SMFFHVertex()  {
  // set up for the couplings
  _couplast=InvEnergy();
  _idlast=0;
  _q2last=ZERO;
  _masslast=ZERO;
  _mw=ZERO;
  _fermion = 0;
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SMFFHVertex::doinit() {
  if ( !_fermion ) {
    // PDG codes for the particles
    // the quarks
    for(int ix=1;ix<7;++ix) {
      addToList(-ix, ix, 25);
    }
    // the leptons
    for(int ix=11;ix<17;ix+=2) {
      addToList(-ix, ix, 25);
    }
  } else {
    addToList(-_fermion,_fermion,25);
  }
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if (!_theSM) 
    throw InitException();
  _mw= getParticleData(ThePEG::ParticleID::Wplus)->mass();
  FFSVertex::doinit();
}

void SMFFHVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _fermion;
}

void SMFFHVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _fermion;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMFFHVertex,FFSVertex>
describeHerwigSMFFHVertex("Herwig::SMFFHVertex", "Herwig.so");

void SMFFHVertex::Init() {

  static ClassDocumentation<SMFFHVertex> documentation
    ("The SMFFHVertex class is the implementation"
     " of the helicity amplitude calculation of the Standard Model Higgs"
     " fermion-antiferiom vertex.");

  static Parameter<SMFFHVertex,int> interfaceFermion
    ("Fermion",
     "The fermion to couple to. If not set all fermions are considered.",
     &SMFFHVertex::_fermion, 0, 0, 16,
     false, false, Interface::limited);
 
}

void SMFFHVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr, tcPDPtr) {
  int iferm=abs(aa->id());
  // left and right couplings set to one
  left (1.);
  right(1.);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==complex<InvEnergy>()) {
    _couplast = -0.5*weakCoupling(q2)/_mw;
    _q2last=q2;
    _idlast=iferm;
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    _masslast=_theSM->mass(q2,aa);
  }
  else if(iferm!=_idlast) {
    _idlast=iferm;
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    _masslast=_theSM->mass(q2,aa);
  }
  norm(_couplast*_masslast);
}
#line 1 "./SMFFPVertex.cc"
// -*- C++ -*-
//
// SMFFPVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardModelFFPVertex class.
//

#include "SMFFPVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMFFPVertex::SMFFPVertex()  : _charge(17,0.0), _couplast(0.), _q2last(0.*GeV2) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SMFFPVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 22);
  }
  // the leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix, 22);
  }
  for(int ix=1;ix<4;++ix) {
    _charge[2*ix-1]  = generator()->standardModel()->ed();
    _charge[2*ix ]   = generator()->standardModel()->eu();
    _charge[2*ix+9 ] = generator()->standardModel()->ee();
    _charge[2*ix+10] = generator()->standardModel()->enu();
  }
  FFVVertex::doinit();
}

void SMFFPVertex::persistentOutput(PersistentOStream & os) const {
  os << _charge;
}

void SMFFPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _charge;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SMFFPVertex,FFVVertex>
describeHerwigSMFFPVertex("Herwig::SMFFPVertex", "Herwig.so");

void SMFFPVertex::Init() {
 static ClassDocumentation<SMFFPVertex> documentation
    ("The SMFFPVertex class is the implementation of"
     "the coupling of the photon to the Standard Model fermions");
}

// coupling for FFP vertex
#ifndef NDEBUG
void SMFFPVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr bb,tcPDPtr) {
#else
void SMFFPVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,tcPDPtr) {
#endif
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
  assert(aa->id()==-bb->id());
  if(aa->id()<0) {
    left(_charge[iferm]);
    right(_charge[iferm]);
  }
  else {
    left(-_charge[iferm]);
    right(-_charge[iferm]);
  }
}
#line 1 "./SMFFWVertex.cc"
// -*- C++ -*-
//
// SMFFWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFWVertex class.
//

#include "SMFFWVertex.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
    
SMFFWVertex::SMFFWVertex() : 
  _diagonal(false), _ckm(3,vector<Complex>(3,0.0)), 
  _couplast(0.), _q2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SMFFWVertex::persistentOutput(PersistentOStream & os) const {
  os << _diagonal << _ckm;
}
  
void SMFFWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _diagonal >> _ckm;
}
  
void SMFFWVertex::doinit() {
  // particles for outgoing W-
  // quarks
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<7;iy+=2) {
      bool isOff = iy/2 != (ix+1)/2;
      if ( isOff && _diagonal )
	continue;
      addToList(-ix, iy, -24);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix+1, -24);
  }
  // particles for outgoing W+
  // quarks
  for(int ix=2;ix<7;ix+=2) {
    for(int iy=1;iy<6;iy+=2) {
      bool isOff = ix/2 != (iy+1)/2;
      if ( isOff && _diagonal )
	continue;
      addToList(-ix, iy, 24);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix-1, ix, 24);
  }
  ThePEG::Helicity::FFVVertex::doinit();
  if ( !_diagonal ) {
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
	  _ckm[ix][iy]=CKM[ix][iy];
	}
      }
    }
    else {
      throw Exception() << "Must have access to the Herwig::StandardCKM object"
			<< "for the CKM matrix in SMFFWVertex::doinit()"
			<< Exception::runerror;
    }
  } else {
    _ckm = vector< vector<Complex > >(3,vector<Complex >(3,1.0));
  }
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SMFFWVertex,FFVVertex>
describeHerwigSMFFWVertex("Herwig::SMFFWVertex", "Herwig.so");
  
void SMFFWVertex::Init() {
  static ClassDocumentation<SMFFWVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the W boson to the Standard Model fermions");


  static Switch<SMFFWVertex,bool> interfaceDiagonal
    ("Diagonal",
     "Use a diagonal CKM matrix (ignoring the CKM object of the StandardModel).",
     &SMFFWVertex::_diagonal, false, false, false);
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
  
// coupling for FFW vertex
void SMFFWVertex::setCoupling(Energy2 q2, tcPDPtr aa, tcPDPtr bb, tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -sqrt(0.5)*weakCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  int ianti=abs(bb->id());
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
    right(0.);
  }
  // leptons
  else if(iferm>=11 && iferm <=16) {
    left(1.);
    right(0.);
  }
  else 
    assert(false);
}







#line 1 "./SMFFZVertex.cc"
// -*- C++ -*-
//
// SMFFZVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMFFZVertex class.
//

#include "SMFFZVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void SMFFZVertex::persistentOutput(PersistentOStream & os) const {
  os << _gl << _gr;
}

void SMFFZVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gl >> _gr;
}

// Static variable needed for the type description system in ThePEG.
DescribeClass<SMFFZVertex,FFVVertex>
describeHerwigSMFFZVertex("Herwig::SMFFZVertex", "Herwig.so");

void SMFFZVertex::Init() {
  static ClassDocumentation<SMFFZVertex> documentation
    ("The SMFFZVertex class is the implementation of"
     "the coupling of the Z boson to the Standard Model fermions");
}

void SMFFZVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = -electroMagneticCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
  // the left and right couplings
  int iferm=abs(aa->id());
  if((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16)) {
    left(_gl[iferm]);
    right(_gr[iferm]);
  }
  else
    throw HelicityConsistencyError() << "SMFFZVertex::setCoupling "
				     << "Unknown particle in Z vertex" 
				     << Exception::runerror;
}

SMFFZVertex::SMFFZVertex() : _gl(17,0.0), _gr(17,0.0),
			     _couplast(0.0), _q2last(ZERO) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::DELTA);
}

void SMFFZVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for(int ix=1;ix<7;++ix) {
    addToList(-ix, ix, 23);
  }
  // the leptons
  for(int ix=11;ix<17;++ix) {
    addToList(-ix, ix, 23);
  }
  tcSMPtr sm = generator()->standardModel();
  double sw2 = sin2ThetaW();
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    _gl[2*ix-1]  = fact*(sm->vd()  + sm->ad() );
    _gl[2*ix ]   = fact*(sm->vu()  + sm->au() );
    _gl[2*ix+9 ] = fact*(sm->ve()  + sm->ae() );
    _gl[2*ix+10] = fact*(sm->vnu() + sm->anu());
    _gr[2*ix-1]  = fact*(sm->vd()  - sm->ad() );
    _gr[2*ix ]   = fact*(sm->vu()  - sm->au() );
    _gr[2*ix+9 ] = fact*(sm->ve()  - sm->ae() );
    _gr[2*ix+10] = fact*(sm->vnu() - sm->anu());
  }
  FFVVertex::doinit();
}
#line 1 "./SMGGGGVertex.cc"
// -*- C++ -*-
//
// SMGGGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGGVertex class.
//

#include "SMGGGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMGGGGVertex::SMGGGGVertex() : _couplast(0.),_q2last() {
  orderInGs(2);
  orderInGem(0);
  colourStructure(ColourStructure::SU3FF);
}

void SMGGGGVertex::doinit() {
  // particles
  addToList(21,21,21,21);
  VVVVVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SMGGGGVertex,VVVVVertex>
describeHerwigSMGGGGVertex("Herwig::SMGGGGVertex", "Herwig.so");

void SMGGGGVertex::Init() {
  static ClassDocumentation<SMGGGGVertex> documentation
    ("The SMGGGGVertex class is the implementation of the"
     " Standard Model quartic gluon coupling");
  
}


// couplings for the GGGG vertex
void SMGGGGVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
			       tcPDPtr,tcPDPtr) {
  // set the order and type
  setType(1);
  setOrder(0,1,2,3);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = sqr(strongCoupling(q2));
    _q2last=q2;
  }
  norm(_couplast);
}

#line 1 "./SMGGGVertex.cc"
// -*- C++ -*-
//
// SMGGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMGGGVertex class.
//

#include "SMGGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMGGGVertex::SMGGGVertex() : _couplast(0.), _q2last(0.*GeV2) {
  orderInGs(1);
  orderInGem(0);
  colourStructure(ColourStructure::SU3F);
}

void SMGGGVertex::doinit() {
  // the particles
  addToList(21,21,21);
  VVVVertex::doinit();
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SMGGGVertex,Helicity::VVVVertex>
describeHerwigSMGGGVertex("Herwig::SMGGGVertex", "Herwig.so");

void SMGGGVertex::Init() {
 static ClassDocumentation<SMGGGVertex> documentation
    ("The SMGGGVertex class is the implementation"
     " of the Standard Model triple gluon vertex.");
  
}

// couplings for the GGG vertex
void SMGGGVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr, tcPDPtr) {
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = strongCoupling(q2);
    _q2last=q2;
  }
  norm(_couplast);
}
#line 1 "./SMWWHVertex.cc"
// -*- C++ -*-
//
// SMWWHVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the WWHVertex class.
//
#include "SMWWHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SMWWHVertex::SMWWHVertex() 
  : _couplast(0.), _q2last(ZERO), _mw(ZERO), _zfact(0.) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SMWWHVertex::doinit() {
  addToList(24,-24, 25);
  addToList(23, 23, 25);
  // parameters
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _zfact = 1./(1.-sin2ThetaW());
  // base class
  VVSVertex::doinit();
}
    
void SMWWHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(_mw,GeV) << _zfact;
}

void SMWWHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_mw,GeV) >> _zfact;
}
    
// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWWHVertex,VVSVertex>
describeHerwigSMWWHVertex("Herwig::SMWWHVertex", "Herwig.so");

void SMWWHVertex::Init() {
  static ClassDocumentation<SMWWHVertex> documentation
    ("The SMWWHVertex class is the implementation"
     " of the helicity amplitude calculation for the coupling of the Standard"
     " Model electroweak gauge bosons to the Higgs.");
}

void SMWWHVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr, tcPDPtr) {
  int ibos=abs(aa->id());
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = weakCoupling(q2) * UnitRemoval::InvE * _mw;
    _q2last=q2;
  }
  if(ibos==24)      norm(_couplast);
  else if(ibos==23) norm(_couplast*_zfact);
  else
    throw HelicityConsistencyError() << "SMWWHVertex::setCoupling "
				     << "Invalid particles in WWH Vertex" 
				     << Exception::runerror;
}
#line 1 "./SMWWWVertex.cc"
// -*- C++ -*-
//
// SMWWWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWWVertex class.
//

#include "SMWWWVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SMWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _zfact; 
}

void SMWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _zfact;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWWWVertex,Helicity::VVVVertex>
describeHerwigSMWWWVertex("Herwig::SMWWWVertex", "Herwig.so");

void SMWWWVertex::Init() {
  static ClassDocumentation<SMWWWVertex> documentation
    ("The SMWWWVertex class is the implementation of the "
     "Standard Model triple electroweak boson coupling.");
  
}
    
// couplings for the WWW vertex
void SMWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b, tcPDPtr c) {
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = electroMagneticCoupling(q2);
    _q2last=q2;
  }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) || 
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )          norm(_couplast);
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) || 
          (ida== 22 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 22 && idc== 24) )     norm(-_couplast);
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) || 
          (ida== 23 && idb==-24 && idc== 24) || 
          (ida== 24 && idb== 23 && idc==-24) )     norm(_couplast*_zfact);
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) || 
          (ida== 23 && idb== 24 && idc==-24) || 
          (ida==-24 && idb== 23 && idc== 24) )     norm(-_couplast*_zfact);
  else
    throw Helicity::HelicityConsistencyError() 
      << "SMWWWVertex::setCoupling "
      << "Invalid particles in WWW Vertex"
      << a->PDGName() << " " << b->PDGName() << " " << c->PDGName() 
      << Exception::runerror;
}

SMWWWVertex::SMWWWVertex() : _zfact(0.),_couplast(0.), 
			     _q2last(sqr(Constants::MaxEnergy)) {
  orderInGem(1);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SMWWWVertex::doinit() {
  addToList(24, -24, 22);
  addToList(24, -24, 23);
  VVVVertex::doinit();
  // factor for the Z vertex
  double sw2=sin2ThetaW();
  _zfact = sqrt((1.-sw2)/sw2);
}
#line 1 "./SMWWWWVertex.cc"
// -*- C++ -*-
//
// SMWWWWVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWWWVertex class.
//

#include "SMWWWWVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SMWWWWVertex::SMWWWWVertex() 
  : _couplast(0.0), _q2last(sqr(Constants::MaxEnergy)), 
    _vfact(4,0.0), _sw2(0.), _cw2(0.) {
  orderInGem(2);
  orderInGs(0);
  colourStructure(ColourStructure::SINGLET);
}

void SMWWWWVertex::doinit() {
  // particles
  addToList(24, -24, 24, -24);
  addToList(23,  24, 23, -24);
  addToList(22,  24, 22, -24);
  addToList(22,  24, 23, -24);
  VVVVVertex::doinit();
  // couplings
  _sw2 = sin2ThetaW();
  _cw2 = 1.-_sw2;
  double sw = sqrt(_sw2);
  double cw = sqrt(_cw2);
  _vfact[0] = -1./_sw2;
  _vfact[1] = _cw2/_sw2;
  _vfact[2] = 1.;
  _vfact[3] = cw/sw;
  // pointer for intermediate particles
  _gamma  = getParticleData(ThePEG::ParticleID::gamma);
  _Z0     = getParticleData(ThePEG::ParticleID::Z0);
  _wplus  = getParticleData(ThePEG::ParticleID::Wplus);
  _wminus = getParticleData(ThePEG::ParticleID::Wminus);
}

void SMWWWWVertex::persistentOutput(PersistentOStream & os) const {
  os << _gamma << _Z0 << _wplus << _wminus
     << _vfact  << _sw2 << _cw2;
}

void SMWWWWVertex::persistentInput(PersistentIStream & is, int) {
  is >> _gamma >> _Z0 >> _wplus >> _wminus
     >> _vfact >> _sw2 >> _cw2;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWWWWVertex,VVVVVertex>
describeHerwigSMWWWWVertex("Herwig::SMWWWWVertex", "Herwig.so");

void SMWWWWVertex::Init() {
  static ClassDocumentation<SMWWWWVertex> documentation
    ("The SMWWWWVertex class is the implementation of the"
     " Standard Model quartic electroweka gauge boson coupling.");
  
}

// couplings for the WWWW vertex
void SMWWWWVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
			       tcPDPtr c,tcPDPtr d) {
  // id's of the particles
  long id[4]={a->id(),b->id(),c->id(),d->id()};
  // order the particles
  int ngamma(0),nz(0);
  int iorder[4];
  for(int ix=0;ix<4;++ix) {
    if      (id[ix]==22) ++ngamma;
    else if (id[ix]==23) ++nz;
  }
  // if photons or Z's
  if(ngamma!=0 || nz!=0) {
    int iy=0;
    // put the photons first
    for(int ix=0;iy<ngamma&&ix<4;++ix) {
      if(id[ix]==22) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the Z bosons
    for(int ix=0;iy<ngamma+nz&&ix<4;++ix) {
      if(id[ix]==23) {
	iorder[iy]=ix;
	++iy;
      }
    }
    // then the W+
    for(int ix=0;iy<3&&ix<4;++ix) {
      if(id[ix]==24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==3);
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24) {
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
      if(id[ix]==24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==2);
    // finally the W-
    for(int ix=0;iy<4&&ix<4;++ix) {
      if(id[ix]==-24) {
	iorder[iy]=ix;
	++iy;
      }
    }
    assert(iy==4);
    setIntermediate(_gamma,_Z0,_sw2,_cw2);
  }
  setOrder(iorder[0],iorder[1],iorder[2],iorder[3]);
  setType(2);
  // first the overall normalisation
  if(q2!=_q2last||_couplast==0.) {
    _couplast = sqr(electroMagneticCoupling(q2));
    _q2last=q2;
  }
  // id's of the first two particles
  int ida(0),idb(0);
  if(iorder[0]==0)      ida = abs(a->id());
  else if(iorder[0]==1) ida = abs(b->id());
  else if(iorder[0]==2) ida = abs(c->id());
  else if(iorder[0]==3) ida = abs(d->id());
  if(iorder[1]==0)      idb = abs(a->id());
  else if(iorder[1]==1) idb = abs(b->id());
  else if(iorder[1]==2) idb = abs(c->id());
  else if(iorder[1]==3) idb = abs(d->id());
  // WWWW coupling
  if(ida==24)               norm(_vfact[0]*_couplast);
  // ZZWW coupling
  else if(ida==23&&idb==23) norm(_vfact[1]*_couplast);
  // gamma gamma WW coupling
  else if(ida==22&&idb==22) norm(_couplast);
  // gamma  Z WW coupling
  else if(ida==22&&idb==23) norm(_vfact[3]*_couplast);
  else assert(false);
}
#line 1 "./SMHGGVertex.cc"
// -*- C++ -*-
//
// SMHGGVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHGGVertex class.
//

#include "SMHGGVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;
using namespace ThePEG;  

SMHGGVertex::SMHGGVertex()
  :_couplast(0.),_q2last(ZERO),_mw(),massopt(1),_minloop(6),
   _maxloop(6),_CoefRepresentation(1) {
  orderInGs(2);
  orderInGem(1);
  colourStructure(ColourStructure::DELTA);
}

void SMHGGVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(21,21,25);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if(!_theSM) 
    throw InitException();
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  VVSLoopVertex::doinit();
  // code to test the partial width
//   Energy mh = getParticleData(25)->mass();
//   Complex I(0.);
//   for(long ix=int(_minloop);ix<=int(_maxloop);++ix) {
//     tcPDPtr qrk = getParticleData(ix);
//     Energy mt = (2 == massopt) ? _theSM->mass(sqr(mh),qrk) : qrk->mass();
//     double lambda = sqr(mt/mh);
//     Complex fl;
//     if(lambda>=0.25) {
//       fl = -2.*sqr(asin(0.5/sqrt(lambda)));
//     }
//     else {
//       double etap = 0.5+sqrt(0.25-lambda);
//       double etam = 0.5-sqrt(0.25-lambda);
//       fl = 0.5*sqr(log(etap/etam))-0.5*sqr(Constants::pi)
// 	-Complex(0.,1.)*Constants::pi*log(etap/etam);
//     }
//     I += 3.*(2.*lambda+lambda*(4.*lambda-1)*fl);
//   }
//   Energy width = sqr(weakCoupling(sqr(mh))*sqr(strongCoupling(sqr(mh))))/36./8.*sqr(mh/_mw)*mh
//     /sqr(4.*sqr(Constants::pi))*std::norm(I)/Constants::pi;
//   cerr << "testing anal " << width/GeV << "\n";
  if(loopToolsInitialized()) Looptools::ltexi();
}

void SMHGGVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << massopt 
     << _minloop << _maxloop << _CoefRepresentation;
}

void SMHGGVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> massopt 
     >> _minloop >> _maxloop >> _CoefRepresentation;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHGGVertex,VVSLoopVertex>
describeHerwigSMHGGVertex("Herwig::SMHGGVertex", "Herwig.so");

void SMHGGVertex::Init() {
  
  static ClassDocumentation<SMHGGVertex> documentation
    ("This class implements the h->g,g vertex");

  static Parameter<SMHGGVertex,int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHGGVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHGGVertex,int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHGGVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);

  static Switch<SMHGGVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHGGVertex::massopt, 1, false, false);
  static SwitchOption interfaceHeavyMass
    (interfaceMassOption,
     "PoleMasses",
     "The loop is calculcated with the pole quark masses",
     1);
  static SwitchOption interfaceNormalMass
    (interfaceMassOption,
     "RunningMasses",
     "running quark masses are taken in the loop",
     2);
  static SwitchOption interfaceInfiniteTopMass
    (interfaceMassOption,
     "InfiniteTopMass",
     "the loop consists of an infinitely massive top quark",
     3);

  static Switch<SMHGGVertex,unsigned int> interfaceScheme
    ("CoefficientScheme",
     "Which scheme for the tensor coefficients is applied",
     &SMHGGVertex::_CoefRepresentation, 1, false, false);
  static SwitchOption interfaceSchemeSimplified
    (interfaceScheme,
     "Simplified",
     "Represection suitable for the simplified the H-g-g and H-gamma-gamma vertices",
     1);
  static SwitchOption interfaceSchemeGeneral
    (interfaceScheme,
     "General",
     "Represection suitable for the Passarino-Veltman tensor reduction scheme",
     2);
}

void SMHGGVertex::setCoupling(Energy2 q2, tcPDPtr part2, tcPDPtr part3, tcPDPtr part1) {
  assert(part1 && part2 && part3);
  assert(part1->id() == ParticleID::h0 &&
	 part2->id() == ParticleID::g  && part3->id() == ParticleID::g );
  int Qminloop = _minloop;
  int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) {
    Qmaxloop=_minloop;
    Qminloop=_maxloop;
  }
  if(massopt==3) {
    if(q2 != _q2last) {
      double g   = weakCoupling(q2);
      double gs2 = sqr(strongCoupling(q2));
      _couplast = UnitRemoval::E * gs2 * g / 16. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    norm(_couplast);
    Complex loop(2./3.);
    a00( loop);    a11(0.0);   a12(0.0);
    a21(-loop);    a22(0.0);   aEp(0.0);
    return;
  }
  switch (_CoefRepresentation) {
  case 1: {
    if(q2 != _q2last||_couplast==0.) {
      double g   = weakCoupling(q2);
      double gs2 = sqr(strongCoupling(q2));
      _couplast = UnitRemoval::E * gs2 * g / 16. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    norm(_couplast);
    Complex loop(0.);
    for ( int i = Qminloop; i <= Qmaxloop; ++i ) {
      tcPDPtr qrk = getParticleData(i);
      Energy mass = (2 == massopt) ? _theSM->mass(q2,qrk) : qrk->mass();
      loop += Af(sqr(mass)/invariant(0,0));
    }
    a00(loop);
    a11(0.0);
    a12(0.0);
    a21(-loop);
    a22(0.0);
    aEp(0.0);
    break;
  }
  case 2: {
    if (q2 != _q2last) {
      Looptools::clearcache();
      _couplast = 0.25*sqr(strongCoupling(q2))*weakCoupling(q2);
      _q2last = q2;
    }
    norm(_couplast);
    int delta = Qmaxloop - Qminloop + 1;
    type.resize(delta,PDT::SpinUnknown);
    masses.resize(delta,ZERO);
    couplings.clear();
    for (int i = 0; i < delta; ++i) {
      tcPDPtr q = getParticleData(_minloop+i);
      type[i] = PDT::Spin1Half;
      masses[i] = (2 == massopt) ? _theSM->mass(q2,q) : q->mass();
      const double ratio = masses[i]/_mw;
      couplings.push_back(make_pair(ratio, ratio));
    }
    setNParticles(delta);
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
    break;
  }
  }
}

Complex SMHGGVertex::Af(double tau) const {
  return tau*(4.- W2(tau)*(1.-4.*tau));
}

Complex SMHGGVertex::W2(double lambda) const {
  double pi = Constants::pi;
  if (0.0 == lambda)     return 0.0;
  else if (lambda < 0.0) return 4.*sqr(asinh(0.5*sqrt(-1./lambda)));

  double root(0.5*sqrt(1./lambda));
  Complex ac(0.);
    // formulae from NPB297,221
  if(root < 1.) {
    ac = -sqr(asin(root));
  } 
  else {
    double ex = acosh(root);
    ac = sqr(ex) - 0.25*sqr(pi) - pi*ex*Complex(0.,1.);
  }
  return 4.*ac;
}
#line 1 "./SMHZPVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHZPVertex class.
//

#include "SMHZPVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMHZPVertex::SMHZPVertex()
  :_couplast(0.),_q2last(),_mw(),_mz(),_massopt(1),
   _minloop(6),_maxloop(6) {
  orderInGs(0);
  orderInGem(3);
  kinematics(true);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SMHZPVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMHZPVertex::fullclone() const {
  return new_ptr(*this);
}

void SMHZPVertex::doinit() {
  GeneralVVSVertex::doinit();
  //PDG codes for particles at vertices
  addToList(23,22,25);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) 
    throw InitException() 
      << "SMHGGVertex::doinit() - The pointer to the SM object is null."
      << Exception::abortnow;
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  _mz = getParticleData(ThePEG::ParticleID::Z0)->mass();
  GeneralVVSVertex::doinit();
}

void SMHZPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << ounit(_mz,GeV) << _massopt
     << _minloop << _maxloop;
}

void SMHZPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> iunit(_mz,GeV) >> _massopt
     >> _minloop >> _maxloop;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<SMHZPVertex,GeneralVVSVertex>
describeHerwigSMHZPVertex("Herwig::SMHZPVertex", "libHerwig.so");

void SMHZPVertex::Init() {

  static ClassDocumentation<SMHZPVertex> documentation
    ("The SMHZPVertex class provides a simple implementation of the "
     "Higgs-Z-Photon loop looping to allow the calculation of the "
     "associated Higgs decay mode H -> Z gamma.");


  static Parameter<SMHZPVertex,int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHZPVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHZPVertex,int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHZPVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);

  static Switch<SMHZPVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHZPVertex::_massopt, 2, false, false);
  static SwitchOption interfaceHeavyMass
    (interfaceMassOption,
     "PoleMasses",
     "The loop is calculcated with the pole quark masses",
     1);
  static SwitchOption interfaceNormalMass
    (interfaceMassOption,
     "RunningMasses",
     "running quark masses are taken in the loop",
     2);
}

void SMHZPVertex::setCoupling(Energy2 q2, tcPDPtr part2,
			      tcPDPtr part3,
#ifndef NDEBUG
			      tcPDPtr part1) {
#else
			      tcPDPtr) {
#endif
  if(part3->id()==ParticleID::Z0) swap(part2,part3);
  assert( part1->id() == ParticleID::h0 &&
 	  part2->id() == ParticleID::Z0 && part3->id() == ParticleID::gamma );
  int Qminloop = _minloop;
  int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) swap(Qmaxloop,Qminloop);
  double cw = sqrt(1.-sin2ThetaW()),sw=sqrt(sin2ThetaW()),tw = sw/cw;
  if(q2 != _q2last||_couplast==0.) {
    double g   = weakCoupling(q2);
    double e2 = sqr(electroMagneticCoupling(q2));
    _couplast = UnitRemoval::E * e2 * g / 16. / _mw/ sqr(Constants::pi);
    _q2last = q2;
  }
  norm(_couplast);
  Complex loop(0.);
  // quark loops
  for ( int i = Qminloop; i <= Qmaxloop; ++i ) {
    tcPDPtr qrk = getParticleData(i);
    Energy mass = (2 == _massopt) ? _theSM->mass(q2,qrk) : qrk->mass();
    double charge = i%2==0 ?
      generator()->standardModel()->eu() : generator()->standardModel()->ed();
    double gv = i%2==0 ?
      generator()->standardModel()->vu() : generator()->standardModel()->vd();
    double tau = 0.25*invariant(0,0)/sqr(mass), lambda(0.25*sqr(_mz/mass));
    loop += 3.*charge*gv *(I1(tau,lambda)-I2(tau,lambda))/(sw*cw);
  }
  // lepton loops
  int Lminloop = 3; // still fixed value
  int Lmaxloop = 3; // still fixed value
  for (int i = Lminloop; i <= Lmaxloop; ++i) {
    tcPDPtr lpt = getParticleData(9 + 2*i);
    Energy mass = (2 == _massopt) ? _theSM->mass(q2,lpt) : lpt->mass();
    double charge = generator()->standardModel()->ee();
    double gv     = generator()->standardModel()->ve();
    double tau = 0.25*invariant(0,0)/sqr(mass), lambda(0.25*sqr(_mz/mass));
    loop += charge*gv*(I1(tau,lambda)-I2(tau,lambda))/(sw*cw);
  }
  // W loop
  double tau = 0.25*invariant(0,0)/sqr(_mw), lambda(0.25*sqr(_mz/_mw));
  loop += ( 4.*(3.-sqr(tw))*I2(tau,lambda) +
	    ((1.+2.*tau)*sqr(tw)-(5.+2.*tau))*I1(tau,lambda))/tw;
  a00(loop);
  a11(0.0);
  a12(0.0);
  a21(-loop);
  a22(0.0);
  aEp(0.0);
  // test of the width calculation
  // Energy mh = getParticleData(25)->mass();
  // Energy pre = sqr(weakCoupling(q2))*pow(electroMagneticCoupling(q2),4)*mh*sqr(mh/_mw)
  //   /128./16./pow(Constants::pi,5)*pow(double(1.-sqr(_mz/mh)),3)*
  //   std::real(std::norm(loop));
}

Complex SMHZPVertex::I1(double tau,double lambda) const {
  return (-0.5+0.5/(tau-lambda)*(f(tau)-f(lambda))+
	  lambda/(tau-lambda)*(g(tau)-g(lambda)))/(tau-lambda);
}

Complex SMHZPVertex::I2(double tau,double lambda) const {
  return 0.5/(tau-lambda)*(f(tau)-f(lambda));
}

Complex SMHZPVertex::f(double tau) const {
  if(tau>0 && tau<= 1.) {
    return sqr(asin(sqrt(tau)));
  }
  else if(tau>1.) {
    double lx = log(sqrt(tau)+sqrt(tau-1));
    return -sqr(lx)+0.25*sqr(Constants::pi)+Complex(0.,1.)*Constants::pi*lx;
  }
  else {
    assert(false);
    return 0.;
  }
}

Complex SMHZPVertex::g(double tau) const {
  if(tau>0 && tau<= 1.) {
    return sqrt((1.-tau)/tau)*asin(sqrt(tau));
  }
  else if(tau>1.) {
    double lx = log(sqrt(tau)+sqrt(tau-1));
    double root = sqrt((tau-1.)/tau);
    return root*(lx-0.5*Complex(0,1)*Constants::pi);
  }
  else {
    assert(false);
    return 0.;
  }
}
#line 1 "./SMHPPVertex.cc"
// -*- C++ -*-
//
// SMHPPVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHPPVertex class.
//

#include "SMHPPVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "Herwig/Looptools/clooptools.h"

using namespace Herwig;
using namespace ThePEG;

void SMHPPVertex::persistentOutput(PersistentOStream & os) const {
  os << _theSM << ounit(_mw,GeV) << _massopt << _minloop << _maxloop 
     << _CoefRepresentation;
}

void SMHPPVertex::persistentInput(PersistentIStream & is, int) {
  is >> _theSM >> iunit(_mw,GeV) >> _massopt >> _minloop >> _maxloop 
     >> _CoefRepresentation;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHPPVertex,VVSLoopVertex>
describeHerwigSMHPPVertex("Herwig::SMHPPVertex", "Herwig.so");

void SMHPPVertex::Init() {
  
  static ClassDocumentation<SMHPPVertex> documentation
    ("This class implements the h0->gamma,gamma vertex.");

  static Parameter<SMHPPVertex,int> interfaceMinQuarkInLoop
    ("MinQuarkInLoop",
     "The minimum flavour of the quarks to include in the loops",
     &SMHPPVertex::_minloop, 6, 1, 6,
     false, false, Interface::limited);

  static Parameter<SMHPPVertex,int> interfaceMaxQuarkInLoop
    ("MaxQuarkInLoop",
     "The maximum flavour of the quarks to include in the loops",
     &SMHPPVertex::_maxloop, 6, 1, 6,
     false, false, Interface::limited);

  static Switch<SMHPPVertex,unsigned int> interfaceMassOption
    ("LoopMassScheme",
     "Switch for the treatment of the masses in the loops ",
     &SMHPPVertex::_massopt, 2, false, false);
  static SwitchOption interfaceHeavyMass
    (interfaceMassOption,
     "PoleMasses",
     "The loop is calculcated with the pole quark masses",
     1);
  static SwitchOption interfaceNormalMass
    (interfaceMassOption,
     "RunningMasses",
     "running quark masses are taken in the loop",
     2);

  static Switch<SMHPPVertex,unsigned int> interfaceScheme
    ("CoefficientScheme",
     "Which scheme for the tensor coefficients is applied",
     &SMHPPVertex::_CoefRepresentation, 1, false, false);
  static SwitchOption interfaceSchemeSimplified
    (interfaceScheme,
     "Simplified",
     "Represection suitable for the simplified the H-g-g and H-gamma-gamma vertices",
     1);
  static SwitchOption interfaceSchemeGeneral
    (interfaceScheme,
     "General",
     "Represection suitable for the Passarino-Veltman tensor reduction scheme",
     2);
}


void SMHPPVertex::setCoupling(Energy2 q2, tcPDPtr part2,
                              tcPDPtr part3, tcPDPtr part1) {
  assert( part1->id() == ParticleID::h0 &&
	  part2->id() == ParticleID::gamma && part3->id() == ParticleID::gamma );
  int Qminloop = _minloop;
  int Qmaxloop = _maxloop;
  if (_maxloop < _minloop) {
    Qmaxloop=_minloop;
    Qminloop=_maxloop;
  }
  switch (_CoefRepresentation) {
  case 1: {
    if(q2 != _q2last||_couplast==0.) {
      double g   = weakCoupling(q2);
      double e2 = sqr(electroMagneticCoupling(q2));
      _couplast = UnitRemoval::E * e2 * g / 8. / _mw/ sqr(Constants::pi);
      _q2last = q2;
    }
    norm(_couplast);
    Complex loop(0.);
    // quark loops
    for ( int i = Qminloop; i <= Qmaxloop; ++i ) {
      tcPDPtr qrk = getParticleData(i);
      Energy mass = (2 == _massopt) ? _theSM->mass(q2,qrk) : qrk->mass();
      Charge charge = qrk->charge();
      loop += Complex(3.*sqr(charge/ThePEG::Units::eplus) * Af(sqr(mass)/invariant(0,0)));
    }
    // lepton loops
    int Lminloop = 3; // still fixed value
    int Lmaxloop = 3; // still fixed value
    for (int i = Lminloop; i <= Lmaxloop; ++i) {
      tcPDPtr lpt = getParticleData(9 + 2*i);
      Energy mass = (2 == _massopt) ? _theSM->mass(q2,lpt) : lpt->mass();
      Charge charge = lpt->charge();
      loop += Complex(sqr(charge/ThePEG::Units::eplus) * Af(sqr(mass)/invariant(0,0)));
    }
    // W loop
    loop += Aw(sqr(_mw)/invariant(0,0));
    a00(loop);
    a11(0.0);
    a12(0.0);
    a21(-loop);
    a22(0.0);
    aEp(0.0);
    break;
  }
  case 2: {
    if(q2 != _q2last||_couplast==0.) {
      Looptools::clearcache();
      double e = electroMagneticCoupling(q2);
      _couplast = pow(e,3)/sqrt(sin2ThetaW());
      _q2last = q2;
    }
    norm(_couplast);
    // quarks
    int delta = Qmaxloop - Qminloop + 1;
    type.resize(delta,PDT::SpinUnknown);
    masses.resize(delta,ZERO);
    for (int i = 0; i < delta; ++i) {
      tcPDPtr q = getParticleData(_minloop+i);
      type[i] = PDT::Spin1Half;
      masses[i] = (2 == _massopt) ? _theSM->mass(q2,q) : q->mass();
      double copl = -masses[i]*3.*sqr(q->iCharge()/3.)/_mw/2.;
      couplings.push_back(make_pair(copl, copl));
    }
    // tau
    type.push_back(PDT::Spin1Half);
    tcPDPtr tau = getParticleData(ParticleID::tauminus);
    masses.push_back(_theSM->mass(q2,tau));
    double copl = -masses.back()*sqr(tau->iCharge()/3.)/_mw/2.;
    couplings.push_back(make_pair(copl, copl));
    // W
    type.push_back(PDT::Spin1);
    masses.push_back(_mw);
    const double mw = UnitRemoval::InvE*_mw;
    couplings.push_back(make_pair(mw,mw));
    setNParticles(delta+2);
    VVSLoopVertex::setCoupling(q2, part1, part2, part3);
    break;
  }
  }
}

Complex SMHPPVertex::Af(const double tau) const {
  return tau*(4. - W2(tau)*(1. - 4.*tau));
}

Complex SMHPPVertex::Aw(const double tau) const {
  return 0.5*(-3.*W2(tau)*tau*(4.*tau - 2.) - 12.*tau - 2.);
}

Complex SMHPPVertex::W2(double lambda) const {
  double pi = Constants::pi;

  if (0.0 == lambda) 
    return 0.0;

  if (lambda < 0.0) 
    return 4.*sqr(asinh(0.5*sqrt(-1./lambda)));

  double root(0.5*sqrt(1./lambda));
  Complex ac(0.);
  // formulae from NPB297,221
  if(root < 1.) {
    ac = -sqr(asin(root));
  } 
  else {
    double ex = acosh(root);
    ac = sqr(ex) - 0.25*sqr(pi) - pi*ex*Complex(0.,1.);
  }
  return 4.*ac;
}

SMHPPVertex::SMHPPVertex() 
  :_couplast(0.),_q2last(),_mw(),_massopt(1),
   _minloop(6),_maxloop(6),_CoefRepresentation(1) {
  orderInGs(0);
  orderInGem(3);
  colourStructure(ColourStructure::SINGLET);
}


// functions for loops for testing
// namespace {

// Complex F0(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return tau*(1.-tau*ft);
// }

// Complex FHalf(double tau,double eta) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return -2.*tau*(eta+(1.-tau*eta)*ft);
// }

// Complex F1(double tau) {
//   Complex ft;
//   if(tau>=1.)
//     ft = sqr(asin(1./sqrt(tau)));
//   else {
//     double etap = 1.+sqrt(1.-tau);
//     double etam = 1.-sqrt(1.-tau);
//     ft = -0.25*sqr(log(etap/etam)-Constants::pi*Complex(0.,1.));
//   }
//   return 2.+3.*tau+3.*tau*(2.-tau)*ft;
// }
// }


void SMHPPVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(22,22,25);
  _theSM = dynamic_ptr_cast<tcHwSMPtr>(generator()->standardModel());
  if( !_theSM ) 
    throw InitException() 
      << "SMHGGVertex::doinit() - The pointer to the SM object is null."
      << Exception::abortnow;
  _mw = getParticleData(ThePEG::ParticleID::Wplus)->mass();
  VVSLoopVertex::doinit();
//   // code to test the partial width
//   Energy mh = getParticleData(25)->mass();
//   Complex I(0.);
//   for(long ix=int(_minloop);ix<=int(_maxloop);++ix) {
//     tcPDPtr qrk = getParticleData(ix);
//     Energy mt = (2 == _massopt) ? _theSM->mass(sqr(mh),qrk) : qrk->mass();
//     double tau = sqr(2.*mt/mh);
//     I += 3.*sqr(double(qrk->iCharge())/3.)*FHalf(tau,1.);
//     cerr << "testing half " << FHalf(tau,1) << " " << Af(0.25*tau) << "\n";
//   }
//   for(long ix=15;ix<=15;++ix) {
//     tcPDPtr qrk = getParticleData(ix);
//     Energy mt = (2 == _massopt) ? _theSM->mass(sqr(mh),qrk) : qrk->mass();
//     double tau = sqr(2.*mt/mh);
//     I += sqr(double(qrk->iCharge())/3.)*FHalf(tau,1.);
//   }
//   I += F1(sqr(2.*_mw/mh));
//   Energy width = sqr(weakCoupling(sqr(mh))*sqr(electroMagneticCoupling(sqr(mh))))
//     /1024./pow(Constants::pi,5)/16.*sqr(mh/_mw)*mh*std::norm(I);
//   cerr << "testing anal " << width/GeV << "\n";
  if(loopToolsInitialized()) Looptools::ltexi();
}

#line 1 "./SMHHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMHHHVertex class.
//

#include "SMHHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMHHHVertex::SMHHHVertex() : ratio_(ZERO), couplast_(0.), q2last_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SMHHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMHHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SMHHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(ratio_,GeV);
}

void SMHHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(ratio_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMHHHVertex,SSSVertex>
describeHerwigSMHHHVertex("Herwig::SMHHHVertex", "Herwig.so");

void SMHHHVertex::Init() {

  static ClassDocumentation<SMHHHVertex> documentation
    ("The SMHHHVertex class implements the triple Higgs"
     " coupling in the Standard Model.");

}

void SMHHHVertex::doinit() {
  addToList(25,25,25);
  SSSVertex::doinit();
  ratio_ = -1.5*sqr(getParticleData(ParticleID::h0)->mass())/
    getParticleData(ParticleID::Wplus)->mass();
}

#ifndef NDEBUG
void SMHHHVertex::setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3) {
#else
void SMHHHVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,tcPDPtr) {
#endif
  assert(part1->id()==ParticleID::h0 &&
	 part2->id()==ParticleID::h0 &&
	 part3->id()==ParticleID::h0 );
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = weakCoupling(q2)*ratio_*UnitRemoval::InvE;
    q2last_=q2;
  }
  norm(couplast_);
}
#line 1 "./SMWWHHVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SMWWHHVertex class.
//

#include "SMWWHHVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SMWWHHVertex::SMWWHHVertex() : ratio_(0.), couplast_(0.), q2last_(ZERO) {
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

IBPtr SMWWHHVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SMWWHHVertex::fullclone() const {
  return new_ptr(*this);
}

void SMWWHHVertex::persistentOutput(PersistentOStream & os) const {
  os << ratio_;
}

void SMWWHHVertex::persistentInput(PersistentIStream & is, int) {
  is >> ratio_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SMWWHHVertex,Helicity::VVSSVertex>
describeHerwigSMWWHHVertex("Herwig::SMWWHHVertex", "Herwig.so");

void SMWWHHVertex::Init() {

  static ClassDocumentation<SMWWHHVertex> documentation
    ("The SMWWHHVertex class implements the coupling of two electroweeak"
     " gauge bosons and the Higgs boson in the Standard Model.");

}

void SMWWHHVertex::doinit() {
  addToList( 23, 23, 25, 25);
  addToList( 24,-24, 25, 25);
  VVSSVertex::doinit();
  ratio_ = 1./(1.-sin2ThetaW());
}

void SMWWHHVertex::setCoupling(Energy2 q2,
			       tcPDPtr part1,tcPDPtr,
#ifndef NDEBUG
			       tcPDPtr part3,tcPDPtr part4) {
#else
			       tcPDPtr,tcPDPtr) {
#endif
  assert(part3->id()==ParticleID::h0 && part4->id()==ParticleID::h0 );
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = sqr(weakCoupling(q2));
    q2last_=q2;
  }
  if(part1->id()==ParticleID::Z0) {
    norm(0.5*couplast_*ratio_);
  }
  else {
    norm(0.5*couplast_);
  }
}
#line 1 "./GenericSVVVertex.cc"
// -*- C++ -*-
//
// GenericSVVVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericSVVVertex class.
//

#include "GenericSVVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;  

GenericSVVVertex::GenericSVVVertex()
  :pids(ZERO),oas(0),oaew(0)
{}

void GenericSVVVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(pids[0],pids[1],pids[2]);
  orderInGs(oas);
  orderInGem(oaew);
  GeneralVVSVertex::doinit();
}



string GenericSVVVertex::dopids(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() != 3 )
    throw InitException() << "accepts only three particles.";

  for ( vector<string>::iterator p = process.begin();
	p != process.end(); ++p ) {
       int tmp;
       istringstream(*p) >> tmp;
       pids.push_back(tmp);
  }
  return "";
}





void GenericSVVVertex::persistentOutput(PersistentOStream & os) const {
  os << pids<<oas<<oaew;
}

void GenericSVVVertex::persistentInput(PersistentIStream & is, int) {
  is >> pids>>oas>>oaew;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GenericSVVVertex,Helicity::GeneralVVSVertex>
describeHerwigGenericSVVVertex("Herwig::GenericSVVVertex", "Herwig.so");

void GenericSVVVertex::Init() {
  
  static ClassDocumentation<GenericSVVVertex> documentation
    ("This class implements the s->v,v vertex");

  static Command<GenericSVVVertex> interfacepids
    ("pids",
     "Set the pids.",
     &GenericSVVVertex::dopids, false);

  static Parameter<GenericSVVVertex, int> interfaceOrderoas
    ("OrderInAlphaS",
     "The order in alpha_S",
     &GenericSVVVertex::oas, 2, 0, 0,
     false, false, Interface::lowerlim);
            
  static Parameter<GenericSVVVertex, int> interfaceOrderoaew
    ("OrderInAlphaEW",
     "The order in alpha_EW",
     &GenericSVVVertex::oaew, 2, 0, 0,
     false, false, Interface::lowerlim);
}

void GenericSVVVertex::setCoupling(Energy2,
#ifndef NDEBUG
				   tcPDPtr part2,
#else
				   tcPDPtr,
#endif
#ifndef NDEBUG
				   tcPDPtr part3,
#else
				   tcPDPtr,
#endif
#ifndef NDEBUG
				   tcPDPtr part1) {
#else
				   tcPDPtr) {
#endif
  assert(part1 && part2 && part3);
  assert(part1->id() == pids[0] &&
	 part2->id() == pids[1]  && part3->id() == pids[2] );
}

#line 1 "./GenericVVVVertex.cc"
// -*- C++ -*-
//
// GenericVVVVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GenericVVVVertex class.
//

#include "GenericVVVVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;  

GenericVVVVertex::GenericVVVVertex()
  :pids(ZERO),oas(0),oaew(0)
{}

void GenericVVVVertex::doinit() {
  //PDG codes for particles at vertices
  addToList(pids[0],pids[1],pids[2]);
  orderInGs(oas);
  orderInGem(oaew);
  VVVVertex::doinit();
}



string GenericVVVVertex::dopids(string in) {
  vector<string> process = StringUtils::split(in);
  if ( process.size() != 3 )
    throw InitException() << "accepts only three particles.";

  for ( vector<string>::iterator p = process.begin();
	p != process.end(); ++p ) {
       int tmp;
       istringstream(*p) >> tmp;
       pids.push_back(tmp);
  }
  return "";
}





void GenericVVVVertex::persistentOutput(PersistentOStream & os) const {
  os << pids<<oas<<oaew;
}

void GenericVVVVertex::persistentInput(PersistentIStream & is, int) {
  is >> pids>>oas>>oaew;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GenericVVVVertex,Helicity::VVVVertex>
describeHerwigGenericVVVVertex("Herwig::GenericVVVVertex", "Herwig.so");

void GenericVVVVertex::Init() {
  
  static ClassDocumentation<GenericVVVVertex> documentation
    ("This class implements the v->v,v vertex");

  static Command<GenericVVVVertex> interfacepids
    ("pids",
     "Set the pids.",
     &GenericVVVVertex::dopids, false);

  static Parameter<GenericVVVVertex, int> interfaceOrderoas
    ("OrderInAlphaS",
     "The order in alpha_S",
     &GenericVVVVertex::oas, 2, 0, 0,
     false, false, Interface::lowerlim);
            
  static Parameter<GenericVVVVertex, int> interfaceOrderoaew
    ("OrderInAlphaEW",
     "The order in alpha_EW",
     &GenericVVVVertex::oaew, 2, 0, 0,
     false, false, Interface::lowerlim);
}

void GenericVVVVertex::setCoupling(Energy2,
#ifndef NDEBUG
				   tcPDPtr part2,
#else
				   tcPDPtr,
#endif
#ifndef NDEBUG
				   tcPDPtr part3,
#else
				   tcPDPtr,
#endif
#ifndef NDEBUG
				   tcPDPtr part1) {
#else
				   tcPDPtr) {
#endif
  assert(part1 && part2 && part3);
  assert(part1->id() == pids[0] &&
	 part2->id() == pids[1]  && part3->id() == pids[2] );
}

