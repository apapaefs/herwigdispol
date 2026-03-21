#line 1 "./SextetModel.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetModel class.
//

#include "SextetModel.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Throw.h"

using namespace Herwig;

IBPtr SextetModel::clone() const {
  return new_ptr(*this);
}

IBPtr SextetModel::fullclone() const {
  return new_ptr(*this);
}

void SextetModel::persistentOutput(PersistentOStream & os) const {
  os << GVVVertex_ << PVVVertex_ << VVVVVertex_
     << GSSVertex_ << PSSVertex_ << VVSSVertex_
     << FFVVertex_ << FFSVertex_
     << g1L_ << g1R_ << g1pR_ << g1ppR_ << g2_ << g2p_ << g3L_
     << enableScalarSingletY43_ << enableScalarSingletY13_ 
     << enableScalarSingletY23_ << enableScalarTripletY13_ 
     << enableVectorDoubletY16_ << enableVectorDoubletY56_;
}

void SextetModel::persistentInput(PersistentIStream & is, int) {
  is >> GVVVertex_ >> PVVVertex_ >> VVVVVertex_
     >> GSSVertex_ >> PSSVertex_ >> VVSSVertex_
     >> FFVVertex_ >> FFSVertex_
     >> g1L_ >> g1R_ >> g1pR_ >> g1ppR_ >> g2_ >> g2p_ >> g3L_
     >> enableScalarSingletY43_ >> enableScalarSingletY13_ 
     >> enableScalarSingletY23_ >> enableScalarTripletY13_ 
     >> enableVectorDoubletY16_ >> enableVectorDoubletY56_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetModel,StandardModel>
  describeSextetModel("Herwig::SextetModel", "HwSextetModel.so");

void SextetModel::Init() {

  static ClassDocumentation<SextetModel> documentation
    ("The SextetModel class provides the Model class for models with new scalars"
     " or vectors in the sextet representation of SU(3)");

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVVVertex>
    interfaceVertexVDQVDQG
    ("Vertex/VDQVDQG",
     "The coupling of the gluon to two vector diquarks",
     &SextetModel::GVVVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVVVertex>
    interfaceVertexVDQVDQP
    ("Vertex/VDQVDQP",
     "The coupling of the photon to two vector diquarks",
     &SextetModel::PVVVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVVVVertex>
    interfaceVertexVDQVDQGG
    ("Vertex/VDQVDQGG",
     "The coupling of two gluons to two vector diquarks",
     &SextetModel::VVVVVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVSSVertex>
    interfaceVertexSDQSDQG
    ("Vertex/SDQSDQG",
     "The coupling of the gluon to two scalar diquarks",
     &SextetModel::GSSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVSSVertex>
    interfaceVertexSDQSDQP
    ("Vertex/SDQSDQP",
     "The coupling of the photon to two scalar diquarks",
     &SextetModel::PSSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractVVSSVertex>
    interfaceVertexSDQSDQGG
    ("Vertex/SDQSDQGG",
     "The coupling of two gluons to two scalar diquarks",
     &SextetModel::VVSSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractFFSVertex> 
    interfaceVertexFFSDQ
    ("Vertex/FFSDQ",
     "The coupling of two quarks to the scalar diquark",
     &SextetModel::FFSVertex_, false, false, true, false, false);

  static Reference<SextetModel,ThePEG::Helicity::AbstractFFVVertex> 
    interfaceVertexFFVDQ
    ("Vertex/FFVDQ",
     "The coupling of two quarks to the vector diquark",
     &SextetModel::FFVVertex_, false, false, true, false, false);

  static ParVector<SextetModel,double> interfaceg1L
    ("g1L",
     "The \\f$SU(2)\\f$ quark-doublet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1L_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg1R
    ("g1R",
     "The \\f$SU(2)\\f$ singlet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1R_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg1RPrime
    ("g1RPrime",
     "The \\f$SU(2)\\f$ singlet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1pR_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg1RDoublePrime
    ("g1RDoublePrime",
     "The \\f$SU(2)\\f$ singlet coupling to \\f$\\Phi_{6,1,1/3}\\f$.",
     &SextetModel::g1ppR_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg2
    ("g2",
     "The coupling to \\f$V^\\mu_{6,2,-1/6}\\f$.",
     &SextetModel::g2_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg2Prime
    ("g2Prime",
     "The coupling to \\f$V^\\mu_{6,2,5/6}\\f$.",
     &SextetModel::g2p_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static ParVector<SextetModel,double> interfaceg3L
    ("g3L",
     "Coupling to \\f$\\Phi_{6,3,1/3}\\f$.",
     &SextetModel::g3L_, 3, 0.0, 0, 0,
     false, false, Interface::nolimits);

  static Command<SextetModel> interfaceEnableParticles
    ("EnableParticles",
     "Enable specfic diquarks",
     &SextetModel::doEnable, false);

}

void SextetModel::doinit() {
  StandardModel::doinit();
  if ( !(enableScalarSingletY43_ || enableScalarSingletY13_ 
	 || enableScalarSingletY23_ || enableScalarTripletY13_ 
	 || enableVectorDoubletY16_ || enableVectorDoubletY56_ )) {
    Throw<Exception>() << "You have not enabled any Sextet diquarks. Use e.g.\n"
		       << "  do Model:EnableParticles Scalar Triplet Y=1/3\n"
		       << "to specify the spin, weak isospin and weak hypercharge." 
		       << Exception::runerror;
  }
  addVertex(GVVVertex_);
  addVertex(PVVVertex_);
  addVertex(VVVVVertex_);
  addVertex(GSSVertex_);
  addVertex(PSSVertex_);
  addVertex(VVSSVertex_);
  addVertex(FFVVertex_);
  addVertex(FFSVertex_);
}

string SextetModel::doEnable(string args) {
  int spin=-1;
  int weak=-1;
  int Y[2]={-1000000,-1000000};
  string orig=args;
  while ( !args.empty() ) {
    string arg = StringUtils::car(args);
    args = StringUtils::cdr(args);
    if      ( arg == "Scalar" ) spin=1;
    else if ( arg == "Vector" ) spin=3;
    else if ( arg == "Singlet" ) weak=1;
    else if ( arg == "Doublet" ) weak=2;
    else if ( arg == "Triplet" ) weak=3;
    else {
      if(arg.find("Y=")==string::npos) continue;
      arg = StringUtils::cdr(arg,"=");
      vector<string> split = StringUtils::split(arg,"/");
      if(split.size()!=2) continue;
      istringstream is1(split[0]);
      is1 >> Y[0];
      istringstream is2(split[1]);
      is2 >> Y[1];
    }
  }
  // check we read a value for all three quantum numbers
  if ( spin <0 || weak<0 || 0 || Y[0]== -1000000) {
    return string("SextetModel:EnableParticles couldn't termine spin, weak") + 
      string(" isospin or hypercharge for ") + orig + ".";
  }
  // check the values of Y
  if(!(Y[1]==3||Y[1]==6)) {
    return string("SextetModel:EnableParticles invalid weak") + 
      string(" hypercharge for ") + orig + ".";
  }
  // the various allowed combinations
  bool found = false;
  if(spin == 1 ) {
    found = true;
    if     ( weak == 1 && Y[0] ==  4 && Y[1] == 3) {
      enableScalarSingletY43_ = true;
    }
    else if( weak == 1 && Y[0] ==  1 && Y[1] == 3) {
      enableScalarSingletY13_ = true;
    }
    else if( weak == 1 && Y[0] == -2 && Y[1] == 3) {
      enableScalarSingletY23_ = true;
    }
    else if( weak == 3 && Y[0] ==  1 && Y[1] == 3) {
      enableScalarTripletY13_ = true;
    }
    else
      found = false;
  }
  else if(spin == 3 && weak == 2) {
    found = true;
    if     ( Y[0] == -1 && Y[1] == 6) {
      enableVectorDoubletY16_ = true;
    }
    else if( Y[0] ==  5 && Y[1] == 6) {
      enableVectorDoubletY56_ = true;
    }
    else
      found = false;
  }
  if(!found)
    return string("SextetModel:EnableParticles invalid combination") + 
      string(" of spin, weak isospin or hypercharge for ") + orig + ".";
  else
    return "";
}
#line 1 "./SextetGSSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGSSVertex class.
//

#include "SextetGSSVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr SextetGSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGSSVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetGSSVertex,Helicity::VSSVertex,false,true>
describeSextetGSSVertex("Herwig::SextetGSSVertex", "HwSextetModel.so");

void SextetGSSVertex::Init() {

  static ClassDocumentation<SextetGSSVertex> documentation
    ("The SextetGSSVertex class implements the coupling of the gluon"
     " to scalar diquarks.");

}

void SextetGSSVertex::doinit() {
  orderInGs (1);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // add the enabled particles
  if(model->ScalarSingletY43Enabled())
    addToList(21,ParticleID::ScalarDQSingletY43,
   	         ParticleID::ScalarDQSingletY43bar);
  if(model->ScalarSingletY13Enabled())
    addToList(21,ParticleID::ScalarDQSingletY13,
	         ParticleID::ScalarDQSingletY13bar);
  if(model->ScalarSingletY23Enabled())
    addToList(21,ParticleID::ScalarDQSingletY23,
	         ParticleID::ScalarDQSingletY23bar);
  if(model->ScalarTripletY13Enabled()) {
    addToList(21,ParticleID::ScalarDQTripletP,
	         ParticleID::ScalarDQTripletPbar);
    addToList(21,ParticleID::ScalarDQTriplet0,
	         ParticleID::ScalarDQTriplet0bar);
    addToList(21,ParticleID::ScalarDQTripletM,
	         ParticleID::ScalarDQTripletMbar);
  }
  Helicity::VSSVertex::doinit();
}

void SextetGSSVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
				  tcPDPtr part1,
#else
				  tcPDPtr ,
#endif
				  tcPDPtr part2, tcPDPtr ) {
  assert(part1->id()==ParticleID::g);
#ifndef NDEBUG
  long idq = abs(part2->id());
#endif
  assert(idq == ParticleID::ScalarDQSingletY43 ||
	 idq == ParticleID::ScalarDQSingletY13 ||
	 idq == ParticleID::ScalarDQSingletY23 ||
	 idq == ParticleID::ScalarDQTripletP   ||
	 idq == ParticleID::ScalarDQTriplet0   ||
	 idq == ParticleID::ScalarDQTripletM);
  if(q2 != q2Last_ || coupLast_ == 0.) {
    coupLast_ = strongCoupling(q2);
    q2Last_   = q2;
  }
  if(part2->id()>0) 
    norm(-coupLast_);
  else
    norm( coupLast_);
}
#line 1 "./SextetPSSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetPSSVertex class.
//

#include "SextetPSSVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr SextetPSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetPSSVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetPSSVertex,Helicity::VSSVertex,false,true>
describeSextetPSSVertex("Herwig::SextetPSSVertex", "HwSextetModel.so");

void SextetPSSVertex::Init() {

  static ClassDocumentation<SextetPSSVertex> documentation
    ("The SextetPSSVertex class implements the coupling of the gluon"
     " to scalar diquarks.");

}

void SextetPSSVertex::doinit() {
  orderInGs (0);
  orderInGem(1);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetPSSVertex::doinit()"
			       << Exception::runerror;
  // add the enabled particles
  if(model->ScalarSingletY43Enabled())
    addToList(22,ParticleID::ScalarDQSingletY43,
   	         ParticleID::ScalarDQSingletY43bar);
  if(model->ScalarSingletY13Enabled())
    addToList(22,ParticleID::ScalarDQSingletY13,
	         ParticleID::ScalarDQSingletY13bar);
  if(model->ScalarSingletY23Enabled())
    addToList(22,ParticleID::ScalarDQSingletY23,
	         ParticleID::ScalarDQSingletY23bar);
  if(model->ScalarTripletY13Enabled()) {
    addToList(22,ParticleID::ScalarDQTripletP,
	         ParticleID::ScalarDQTripletPbar);
    addToList(22,ParticleID::ScalarDQTriplet0,
	         ParticleID::ScalarDQTriplet0bar);
    addToList(22,ParticleID::ScalarDQTripletM,
	         ParticleID::ScalarDQTripletMbar);
  }
  Helicity::VSSVertex::doinit();
}

void SextetPSSVertex::setCoupling(Energy2 q2,
#ifndef NDEBUG
				  tcPDPtr part1,
#else
				  tcPDPtr ,
#endif
				  tcPDPtr part2, tcPDPtr ) {
  assert(part1->id()==ParticleID::gamma);
  tcPDPtr sca = part2->id()>0 ? part2 : tcPDPtr(part2->CC());
  assert(sca->id() == ParticleID::ScalarDQSingletY43 ||
	 sca->id() == ParticleID::ScalarDQSingletY13 ||
	 sca->id() == ParticleID::ScalarDQSingletY23 ||
	 sca->id() == ParticleID::ScalarDQTripletP   ||
	 sca->id() == ParticleID::ScalarDQTriplet0   ||
	 sca->id() == ParticleID::ScalarDQTripletM);
  if(q2 != q2Last_ || coupLast_ == 0.) {
    coupLast_ = electroMagneticCoupling(q2);
    q2Last_   = q2;
  }
  if(part2->id()>0) 
    norm(-sca->iCharge()/3.*coupLast_);
  else
    norm( sca->iCharge()/3.*coupLast_);
}
#line 1 "./SextetGGSSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGGSSVertex class.
//

#include "SextetGGSSVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "SextetModel.h"
#include "SextetParticles.h"

using namespace Herwig;

SextetGGSSVertex::SextetGGSSVertex() : q2last_(), couplast_() {
  colourStructure(ColourStructure::SU3TT6);
}

IBPtr SextetGGSSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGGSSVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetGGSSVertex,VVSSVertex>
describeSextetGGSSVertex("Herwig::SextetGGSSVertex", "HwSextetModel.so");

void SextetGGSSVertex::Init() {

  static ClassDocumentation<SextetGGSSVertex> documentation
    ("The SextetGGSSVertex class implements the coupling of two gluons to two"
     " scalar sextets");

}

void SextetGGSSVertex::setCoupling(Energy2 q2, tcPDPtr, tcPDPtr, tcPDPtr,
				   tcPDPtr) { 
  if(q2 != q2last_ || couplast_ == 0.) {
    couplast_ = sqr(strongCoupling(q2));
    q2last_ = q2;
  }
  norm(couplast_);
}

void SextetGGSSVertex::doinit() {
  orderInGs(2);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // add the enabled particles
  if(model->ScalarSingletY43Enabled())
    addToList(21,21,ParticleID::ScalarDQSingletY43,
	            ParticleID::ScalarDQSingletY43bar);
  if(model->ScalarSingletY13Enabled())
    addToList(21,21,ParticleID::ScalarDQSingletY13,
	            ParticleID::ScalarDQSingletY13bar);
  if(model->ScalarSingletY23Enabled())
    addToList(21,21,ParticleID::ScalarDQSingletY23,
	            ParticleID::ScalarDQSingletY23bar);
  if(model->ScalarTripletY13Enabled()) {
    addToList(21,21,ParticleID::ScalarDQTripletP,
	            ParticleID::ScalarDQTripletPbar);
    addToList(21,21,ParticleID::ScalarDQTriplet0,
	            ParticleID::ScalarDQTriplet0bar);
    addToList(21,21,ParticleID::ScalarDQTripletM,
	            ParticleID::ScalarDQTripletMbar);
  }
  VVSSVertex::doinit();
}
#line 1 "./SextetGVVVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGVVVertex class.
//

#include "SextetModel.h"
#include "SextetGVVVertex.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

IBPtr SextetGVVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGVVVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetGVVVertex,Helicity::VVVVertex>
describeSextetGVVVertex("Herwig::SextetGVVVertex", "HwSextetModel.so");

void SextetGVVVertex::Init() {

  static ClassDocumentation<SextetGVVVertex> documentation
    ("The SextetGVVVertex class implements the coupling of the gluon to two"
     " vector sextet particles");

}

void SextetGVVVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGVVVertex::doinit()"
			       << Exception::runerror;
  if(model->VectorDoubletY16Enabled()) {
    addToList(21,ParticleID::VectorDQY16P,
	         ParticleID::VectorDQY16Pbar);
    addToList(21,ParticleID::VectorDQY16M,
	         ParticleID::VectorDQY16Mbar);

  }
  if(model->VectorDoubletY56Enabled()) {
    addToList(21,ParticleID::VectorDQY56P,
	         ParticleID::VectorDQY56Pbar);
    addToList(21,ParticleID::VectorDQY56M,
	         ParticleID::VectorDQY56Mbar);
  }
  VVVVertex::doinit();
}

void SextetGVVVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2, 
				  tcPDPtr p3) {
  if(q2 != q2Last_ || coupLast_ == 0.) {
    q2Last_ = q2;
    coupLast_ = strongCoupling(q2);
  }
  if((p1->id()==ParticleID::g&&p2->id()>0&&p3->id()<0)||
     (p2->id()==ParticleID::g&&p3->id()>0&&p1->id()<0)||
     (p3->id()==ParticleID::g&&p1->id()>0&&p2->id()<0))
    norm(-coupLast_);
  else
    norm( coupLast_);
}
#line 1 "./SextetPVVVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetPVVVertex class.
//

#include "SextetModel.h"
#include "SextetPVVVertex.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

IBPtr SextetPVVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetPVVVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetPVVVertex,Helicity::VVVVertex>
describeSextetPVVVertex("Herwig::SextetPVVVertex", "HwSextetModel.so");

void SextetPVVVertex::Init() {

  static ClassDocumentation<SextetPVVVertex> documentation
    ("The SextetPVVVertex class implements the coupling of the gluon to two"
     " vector sextet particles");

}

void SextetPVVVertex::doinit() {
  orderInGs(1);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetPVVVertex::doinit()"
			       << Exception::runerror;
  if(model->VectorDoubletY16Enabled()) {
    addToList(22,ParticleID::VectorDQY16P,
	         ParticleID::VectorDQY16Pbar);
    addToList(22,ParticleID::VectorDQY16M,
	         ParticleID::VectorDQY16Mbar);

  }
  if(model->VectorDoubletY56Enabled()) {
    addToList(22,ParticleID::VectorDQY56P,
	         ParticleID::VectorDQY56Pbar);
    addToList(22,ParticleID::VectorDQY56M,
	         ParticleID::VectorDQY56Mbar);
  }
  VVVVertex::doinit();
}

void SextetPVVVertex::setCoupling(Energy2 q2, tcPDPtr p1, tcPDPtr p2, 
				  tcPDPtr p3) {
  if(q2 != q2Last_ || coupLast_ == 0.) {
    q2Last_ = q2;
    coupLast_ = electroMagneticCoupling(q2);
  }
  if(p1->id()==ParticleID::gamma) {
    norm(p3->iCharge()/3.*coupLast_);
  }
  else if(p2->id()==ParticleID::gamma) {
    norm(p1->iCharge()/3.*coupLast_);
  }
  else if(p3->id()==ParticleID::gamma) {
    norm(p2->iCharge()/3.*coupLast_);
  }
  else
    assert(false);
}
#line 1 "./SextetGGVVVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetGGVVVertex class.
//

#include "SextetGGVVVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Herwig;

IBPtr SextetGGVVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetGGVVVertex::fullclone() const {
  return new_ptr(*this);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeNoPIOClass<SextetGGVVVertex,Helicity::VVVVVertex>
describeSextetGGVVVertex("Herwig::SextetGGVVVertex", "HwSextetModel.so");

void SextetGGVVVertex::Init() {

  static ClassDocumentation<SextetGGVVVertex> documentation
    ("The SextetGGVVVertex class implements the coupling of two gluons to two vector"
     " sextet particles.");

}

void SextetGGVVVertex::doinit() {
  orderInGs(2);
  orderInGem(0);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGGVVVertex::doinit()"
			       << Exception::runerror;
  if(model->VectorDoubletY16Enabled()) {
    addToList(21,21,ParticleID::VectorDQY16P,
	            ParticleID::VectorDQY16Pbar);
    addToList(21,21,ParticleID::VectorDQY16M,
	            ParticleID::VectorDQY16Mbar);

  }
  if(model->VectorDoubletY56Enabled()) {
    addToList(21,21,ParticleID::VectorDQY56P,
	            ParticleID::VectorDQY56Pbar);
    addToList(21,21,ParticleID::VectorDQY56M,
	            ParticleID::VectorDQY56Mbar);
  }
  VVVVVertex::doinit();
}

void SextetGGVVVertex::setCoupling(Energy2 q2, tcPDPtr , tcPDPtr , 
				   tcPDPtr , tcPDPtr ) {
  if(q2 != q2Last_ || coupLast_ == 0. ) {
    q2Last_ = q2;
    coupLast_ = sqr(strongCoupling(q2));
  }
  norm(coupLast_);
  setType(1);
  setOrder(0,1,2,3);
}
#line 1 "./SextetFFSVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetFFSVertex class.
//

#include "SextetFFSVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr SextetFFSVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetFFSVertex::fullclone() const {
  return new_ptr(*this);
}

void SextetFFSVertex::persistentOutput(PersistentOStream & os) const {
  os << g1L_ << g1R_ << g1pR_ << g1ppR_ << g3L_;
}

void SextetFFSVertex::persistentInput(PersistentIStream & is, int) {
  is >> g1L_ >> g1R_ >> g1pR_ >> g1ppR_ >> g3L_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetFFSVertex,Helicity::FFSVertex>
describeSextetFFSVertex("Herwig::SextetFFSVertex", "HwSextetModel.so");

void SextetFFSVertex::Init() {

  static ClassDocumentation<SextetFFSVertex> documentation
    ("The SextetFFSVertex class implements the coupling of two "
     "fermions to a scalar sextet particle.");

}

void SextetFFSVertex::doinit() {
  orderInGs (0);
  orderInGem(1);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // extract the couplings
  g1L_     = model->g1L();
  g1R_     = model->g1R();
  g1pR_   = model->g1pR();
  g1ppR_ = model->g1ppR();
  g3L_     = model->g3L();
  // add the enabled particles
  if(model->ScalarSingletY43Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      if(g1ppR_[ix]!=0.) {
        addToList(  iu,  iu, ParticleID::ScalarDQSingletY43bar);
        addToList( -iu, -iu, ParticleID::ScalarDQSingletY43);
      }
    }
  }
  if(model->ScalarSingletY13Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g1L_[ix]!=0. || g1R_[ix]!=0.) {
        addToList(  id,  iu, ParticleID::ScalarDQSingletY13bar);
        addToList( -id, -iu, ParticleID::ScalarDQSingletY13);
      }
    }
  }
  if(model->ScalarSingletY23Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long id = 2*ix + 1;
      if(g1pR_[ix]!=0. ) {
        addToList( id, id,ParticleID::ScalarDQSingletY23bar);
        addToList(-id,-id,ParticleID::ScalarDQSingletY23);
      }
    }
  }
  if(model->ScalarTripletY13Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g3L_[ix]!=0. ) {
        addToList(  iu,  iu, ParticleID::ScalarDQTripletPbar);
        addToList( -iu, -iu, ParticleID::ScalarDQTripletP);
        addToList(  iu,  id, ParticleID::ScalarDQTriplet0bar);
        addToList( -iu, -id, ParticleID::ScalarDQTriplet0);
        addToList(  id,  id, ParticleID::ScalarDQTripletMbar);
        addToList( -id, -id, ParticleID::ScalarDQTripletM);
      }
    }
  }
  Helicity::FFSVertex::doinit();
}


void SextetFFSVertex::setCoupling(Energy2,tcPDPtr part1,
				  tcPDPtr part2,tcPDPtr part3) {
  long q1ID=(abs(part1->id())), q2ID=(abs(part2->id())),
    sDQID=(abs(part3->id()));
  //check scalar diquark
  assert( sDQID == ParticleID::ScalarDQSingletY43 || 
          sDQID == ParticleID::ScalarDQSingletY13 ||
          sDQID == ParticleID::ScalarDQSingletY23 ||
          sDQID == ParticleID::ScalarDQTripletP ||
          sDQID == ParticleID::ScalarDQTriplet0 ||
          sDQID == ParticleID::ScalarDQTripletM); 
  //check quarks
  assert(!(q1ID>6) && !(q2ID>6));
  bool part1Up = (q1ID==2 || q1ID==4 || q1ID==6);
#ifndef NDEBUG
  bool part2Up = (q2ID==2 || q2ID==4 || q2ID==6);
#endif
  Complex cRight, cLeft, prefactor(1.);
  if(sDQID==ParticleID::ScalarDQSingletY43){
    //should both be up type
    assert(part1Up && part2Up);
    if(q1ID==2)
      cRight=Complex(g1ppR_[0]);
    else if(q1ID==4)
      cRight=Complex(g1ppR_[1]);
    else
      cRight=Complex(g1ppR_[2]);
    cLeft=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQSingletY13){
    //should be one up one down type
    assert((part1Up && !part2Up) || (!part1Up && part2Up));
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;   
    if(upType==2){
      cRight=Complex(g1R_[0]);
      cLeft=Complex(2.*g1L_[0]);
    }
    else if(upType==4){
      cRight=Complex(g1R_[1]);
      cLeft=Complex(2.*g1L_[1]);
    }
    else {
      cRight=Complex(g1R_[2]);
      cLeft=Complex(2.*g1L_[2]);
    }
  }
  if(sDQID==ParticleID::ScalarDQSingletY23){
    //should both be down type
    assert(!part1Up && !part2Up);
    if(q1ID==1)
      cRight=Complex(g1pR_[0]);
    else if(q1ID==3)
      cRight=Complex(g1pR_[1]);
    else
      cRight=Complex(g1pR_[2]);
    cLeft=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQTripletP){
    //should both be up type
    assert(part1Up && part2Up);
    if(q1ID==2)
      cLeft=Complex(g3L_[0]);
    else if(q1ID==4)
      cLeft=Complex(g3L_[1]);
    else
      cLeft=Complex(g3L_[2]);  
    cRight=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQTriplet0){
    //should both one up and down type
    assert((part1Up && !part2Up) || (!part1Up && part2Up));
    //possibly doesn't couple
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;   
    if(upType==2)
      cLeft=Complex(g3L_[0]);
    else if(upType==4)
      cLeft=Complex(g3L_[1]);
    else
      cLeft=Complex(g3L_[2]);  
    cRight=Complex(0.);
  }
  if(sDQID==ParticleID::ScalarDQTripletM){
    //should one both be down type
    assert(!part1Up && !part2Up);
    if(q1ID==1)
      cLeft=Complex(g3L_[0]);
    else if(q1ID==3)
      cLeft=Complex(g3L_[1]);
    else
      cLeft=Complex(g3L_[2]);
    cRight=Complex(0.);
    prefactor=Complex(-1.); 
  }
  left(cLeft);
  right(cRight);
  norm(prefactor);
}
#line 1 "./SextetFFVVertex.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SextetFFVVertex class.
//

#include "SextetFFVVertex.h"
#include "SextetModel.h"
#include "SextetParticles.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

IBPtr SextetFFVVertex::clone() const {
  return new_ptr(*this);
}

IBPtr SextetFFVVertex::fullclone() const {
  return new_ptr(*this);
}

void SextetFFVVertex::persistentOutput(PersistentOStream & os) const {
  os << g2_ << g2p_ ;
}

void SextetFFVVertex::persistentInput(PersistentIStream & is, int) {
  is >> g2_ >> g2p_ ;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SextetFFVVertex,Helicity::FFVVertex>
describeSextetFFVVertex("Herwig::SextetFFVVertex", "HwSextetModel.so");

void SextetFFVVertex::Init() {

  static ClassDocumentation<SextetFFVVertex> documentation
    ("The SextetFFVVertex class implements the coupling of two "
     "fermions to a scalar sextet particle.");

}

void SextetFFVVertex::doinit() {
  orderInGs (0);
  orderInGem(1);
  SextetModelPtr model = 
    dynamic_ptr_cast<SextetModelPtr>(generator()->standardModel());
  if(!model) throw Exception() << "Must be using the SextetModel"
			       << " in SextetGSSVertex::doinit()"
			       << Exception::runerror;
  // extract the couplings
  g2_  = model->g2();
  g2p_= model->g2p();
  // add the enabled particles
  if(model->VectorDoubletY16Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g2_[ix]!=0.) {
        addToList(  -id, -iu, ParticleID::VectorDQY16P);
        addToList(  id,  iu, ParticleID::VectorDQY16Pbar);
        addToList( -id, -id, ParticleID::VectorDQY16M);
        addToList(  id,  id, ParticleID::VectorDQY16Mbar);
      }
    }
  }
  if(model->VectorDoubletY56Enabled()) {
    for(long ix=0;ix<3;++ix) {
      long iu = 2*ix + 2;
      long id = 2*ix + 1;
      if(g2p_[ix]!=0.) {
        addToList( -iu, -iu, ParticleID::VectorDQY56P);
        addToList(  iu,  iu, ParticleID::VectorDQY56Pbar);
        addToList( -id, -iu, ParticleID::VectorDQY56M);
        addToList(  id,  iu, ParticleID::VectorDQY56Mbar);
      }
    }
  }
  Helicity::FFVVertex::doinit();
}


void SextetFFVVertex::setCoupling(Energy2, tcPDPtr part1,
				  tcPDPtr part2,tcPDPtr part3) {
  long q1ID=(abs(part1->id())), q2ID=(abs(part2->id())),
    vDQID=(abs(part3->id()));
  //check scalar diquark
  assert( vDQID == ParticleID::VectorDQY16P || 
          vDQID == ParticleID::VectorDQY16M ||
          vDQID == ParticleID::VectorDQY56P ||
          vDQID == ParticleID::VectorDQY56M);
  //check quarks
  assert(!(q1ID>6) && !(q2ID>6));
  bool part1Up = (q1ID==2 || q1ID==4 || q1ID==6) ? true : false;
#ifndef NDEBUG
  bool part2Up = (q2ID==2 || q2ID==4 || q2ID==6) ? true : false;
#endif
  Complex cRight(1.,1.), cLeft(1.,1.), prefactor(1.,0.);

  if(vDQID==ParticleID::VectorDQY16P){
    //should be one up and down type
    assert((!part1Up && part2Up) || (part1Up && !part2Up));
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;
    if(upType==2)
      cRight=Complex(g2_[0]);
    else if(upType==4)
      cRight=Complex(g2_[1]);
    else
      cRight=Complex(g2_[2]);
    cLeft=Complex(0.);
  }
  if(vDQID==ParticleID::VectorDQY16M) {
    //should be both be down type
    assert(!part1Up && !part2Up);
    if(q1ID==1)
      cRight=Complex(g2_[0]);
    else if(q1ID==2)
      cRight=Complex(g2_[1]);
    else
      cRight=Complex(g2_[2]);
    cLeft=Complex(0.);
  }
  if(vDQID==ParticleID::VectorDQY56P) {
    //should both be up type
    assert(part1Up && part2Up);
    if(q1ID==2)
      cRight=Complex(g2p_[0]);
    else if(q1ID==4)
      cRight=Complex(g2p_[1]);
    else
      cRight=Complex(g2p_[2]);
    cLeft=Complex(0.);
  }
  if(vDQID==ParticleID::VectorDQY56M){
    //should be one up and down type
    assert((!part1Up && part2Up) || (part1Up && !part2Up));
    long upType;
    if(part1Up)
      upType=q1ID;
    else
      upType=q2ID;
    if(upType==2)
      cRight=Complex(g2p_[0]);
    else if(upType==4)
      cRight=Complex(g2p_[1]);
    else
      cRight=Complex(g2p_[2]);
    cLeft=Complex(0.);
  }
  left(cLeft);
  right(cRight);
  norm(prefactor);
}



  


