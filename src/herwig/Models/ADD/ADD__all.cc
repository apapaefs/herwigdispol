#line 1 "./ADDModel.cc"
// -*- C++ -*-
//
// ADDModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModel class.
//

#include "ADDModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void ADDModel::doinit() {
  addVertex(FFGRVertex_);
  addVertex(VVGRVertex_);
  addVertex(SSGRVertex_);
  addVertex(FFGGRVertex_);
  addVertex(FFWGRVertex_);
  addVertex(GGGGRVertex_);
  addVertex(WWWGRVertex_);
  BSMModel::doinit();
}

void ADDModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(mPlanckBar_,GeV) << ounit(md_,GeV) << delta_
     << ounit(lambdaT_,GeV)
     << FFGRVertex_ << VVGRVertex_ << SSGRVertex_ 
     << FFGGRVertex_ << FFWGRVertex_ 
     << GGGGRVertex_ << WWWGRVertex_;
}

void ADDModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mPlanckBar_,GeV) >> iunit(md_,GeV) >> delta_
     >> iunit(lambdaT_,GeV)
     >> FFGRVertex_ >> VVGRVertex_ >> SSGRVertex_
     >> FFGGRVertex_ >> FFWGRVertex_ 
     >> GGGGRVertex_ >> WWWGRVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModel,BSMModel>
describeHerwigADDModel("Herwig::ADDModel", "HwADDModel.so");

void ADDModel::Init() {
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractFFTVertex> interfaceVertexFFGR
    ("Vertex/FFGR",
     "Reference to the fermion-fermion-graviton vertex",
     &ADDModel::FFGRVertex_, false, false, true, false, false);

  static Reference<ADDModel,ThePEG::Helicity::AbstractVVTVertex> interfaceVertexVVGR
    ("Vertex/VVGR",
     "Reference to the vector-vector-graviton vertex",
     &ADDModel::VVGRVertex_, false, false, true, false, false);

  static Reference<ADDModel,ThePEG::Helicity::AbstractSSTVertex> interfaceVertexSSGR
    ("Vertex/SSGR",
     "Reference to the scalar-scalar-graviton vertex",
     &ADDModel::SSGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFGGR
    ("Vertex/FFGGR",
     "Reference to the fermion-antifermion-gluon graviton vertex",
     &ADDModel::FFGGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFWGR
    ("Vertex/FFWGR",
     "Reference to the fermion-antifermion-weak vector boson graviton vertex",
     &ADDModel::FFWGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexGGGGR
    ("Vertex/GGGGR",
     "Reference to the three gluon graviton vertex",
     &ADDModel::GGGGRVertex_, false, false, true, false, false);
  
  static Reference<ADDModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexWWWGR
    ("Vertex/WWWGR",
     "Reference to the three weak vector boson graviton vertex",
     &ADDModel::WWWGRVertex_, false, false, true, false, false);
  
  static ClassDocumentation<ADDModel> documentation
    ("The ADDModel class replaces the Standard Model class for the"
     " ADD model");
  
  static Parameter<ADDModel,unsigned int> interfaceDelta
    ("Delta",
     "Number of extra dimensions",
     &ADDModel::delta_, 2, 2, 1000,
     false, false, Interface::limited);

  static Parameter<ADDModel,Energy> interfaceReducedPlanckMass
    ("Reduced4dPlanckMass",
     "The reduced planck mass in 4 dimensions",
     &ADDModel::mPlanckBar_, GeV, 2.4e18*GeV, 1e17*GeV, 1e20*GeV,
     false, false, Interface::limited);

  static Parameter<ADDModel,Energy> interfaceDdPlanckMass
    ("DdPlanckMass",
     "The d dimension planck mass",
     &ADDModel::md_, GeV, 1000.*GeV, 100.0*GeV, 1e6*GeV,
     false, false, Interface::limited);

  static Parameter<ADDModel,Energy> interfaceLambdaT
    ("LambdaT",
     "The cut-off for virtual graviton processes",
     &ADDModel::lambdaT_, GeV, 1000.*GeV, 100.*GeV, 100000.0*GeV,
     false, false, Interface::limited);

}
#line 1 "./ADDModelFFGRVertex.cc"
// -*- C++ -*-
//
// ADDModelFFGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelFFGRVertex class.
//

#include "ADDModelFFGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void ADDModelFFGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelFFGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelFFGRVertex,FFTVertex>
describeHerwigADDModelFFGRVertex("Herwig::ADDModelFFGRVertex", "HwADDModel.so");

void ADDModelFFGRVertex::Init() {
  static ClassDocumentation<ADDModelFFGRVertex> documentation
    ("The ADDModelFFGRVertex class is the ADDModel calculation"
     " of the fermion-antifermion-graviton vertex");
  
}
  
void ADDModelFFGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(kappa_ * UnitRemoval::E));
}

ADDModelFFGRVertex::ADDModelFFGRVertex() : kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void ADDModelFFGRVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for (int ix=1;ix<7;++ix) addToList(-ix,ix,39);
  // the leptons
  for (int ix=11;ix<17;++ix) addToList(-ix,ix,39);
  FFTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD)
    throw Exception() << "Must have ADDModel in ADDModelFFGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

Complex ADDModelFFGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
				       Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
#line 1 "./ADDModelFFGGRVertex.cc"
// -*- C++ -*-
//
// ADDModelFFGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelFFGGRVertex class.
//

#include "ADDModelFFGGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"


using namespace Herwig;
using namespace ThePEG;

ADDModelFFGGRVertex::ADDModelFFGGRVertex() 
  : couplast_(0.), q2last_(ZERO), kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (1);
  colourStructure(ColourStructure::SU3TFUND);
}

void ADDModelFFGGRVertex::doinit() {
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21,39);
  }
  FFVTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) throw Exception() 
	      << "Must have ADDModel in ADDModelFFGGRVertex::doinit()"
	      << Exception::runerror;
  kappa_ = 2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelFFGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelFFGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelFFGGRVertex,FFVTVertex>
describeHerwigADDModelFFGGRVertex("Herwig::ADDModelFFGGRVertex", "HwADDModel.so");

void ADDModelFFGGRVertex::Init() {
  static ClassDocumentation<ADDModelFFGGRVertex> documentation
    ("The ADDModelFFGGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the ADD model.");
  
}

#ifndef NDEBUG
void ADDModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,
				      tcPDPtr cc, tcPDPtr) {
#else
void ADDModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				      tcPDPtr, tcPDPtr) {
#endif
  // work out the particles
  assert(cc->id()==ParticleID::g && abs(aa->id()) <= 6);
  // overall factor
  if(q2last_!=q2||couplast_==0.) {
    couplast_ = strongCoupling(q2);
    q2last_=q2;
  }
  left (1.);
  right(1.);
  // set the coupling
  norm(UnitRemoval::E * kappa_ * couplast_);
}

Complex ADDModelFFGGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}

#line 1 "./ADDModelFFWGRVertex.cc"
// -*- C++ -*-
//
// ADDModelFFWGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelFFWGRVertex class.
//

#include "ADDModelFFWGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelFFWGRVertex::ADDModelFFWGRVertex() 
  : charge_(17,0.), gl_(17,0.), gr_(17,0.),
    ckm_(3,vector<Complex>(3,0.0)), couplast_(0.),
    q2last_(ZERO), kappa_(ZERO), r_(ZERO) {
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void ADDModelFFWGRVertex::doinit() {
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,22,39);
    addToList(-ix,ix,23,39);
  }
  for(int ix=11;ix<17;++ix) {
    addToList(-ix,ix,22,39);
    addToList(-ix,ix,23,39);
  }
  // particles for outgoing W-
  // quarks
  for(int ix=1;ix<6;ix+=2) {
    for(int iy=2;iy<7;iy+=2) {
      addToList(-ix, iy, -24,39);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix, ix+1, -24,39);
  }
  // particles for outgoing W+
  // quarks
  for(int ix=2;ix<7;ix+=2) {
    for(int iy=1;iy<6;iy+=2) {
      addToList(-ix, iy, 24,39);
    }
  }
  // leptons
  for(int ix=11;ix<17;ix+=2) {
    addToList(-ix-1, ix, 24,39);
  }
  FFVTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) throw Exception() 
	      << "Must have ADDModel in ADDModelFFWGRVertex::doinit()"
	      << Exception::runerror;
  double sw2 = sin2ThetaW();
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    charge_[2*ix-1]  = hwADD->ed();
    charge_[2*ix ]   = hwADD->eu();
    charge_[2*ix+9 ] = hwADD->ee();
    charge_[2*ix+10] = hwADD->enu();
    gl_[2*ix-1]  = fact*(hwADD->vd()  + hwADD->ad() );
    gl_[2*ix ]   = fact*(hwADD->vu()  + hwADD->au() );
    gl_[2*ix+9 ] = fact*(hwADD->ve()  + hwADD->ae() );
    gl_[2*ix+10] = fact*(hwADD->vnu() + hwADD->anu());
    gr_[2*ix-1]  = fact*(hwADD->vd()  - hwADD->ad() );
    gr_[2*ix ]   = fact*(hwADD->vu()  - hwADD->au() );
    gr_[2*ix+9 ] = fact*(hwADD->ve()  - hwADD->ae() );
    gr_[2*ix+10] = fact*(hwADD->vnu() - hwADD->anu());
  }
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
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
    throw Exception() << "Must have access to the Herwig::StandardCKM object"
		      << "for the CKM matrix in SMFFWVertex::doinit()"
		      << Exception::runerror;
  }
}

void ADDModelFFWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << charge_ << gl_ << gr_ << ounit(kappa_,InvGeV) << ckm_ << ounit(r_,GeV);
}

void ADDModelFFWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> charge_ >> gl_ >> gr_ >> iunit(kappa_,InvGeV) >> ckm_ >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelFFWGRVertex,FFVTVertex>
describeHerwigADDModelFFWGRVertex("Herwig::ADDModelFFWGRVertex", "HwADDModel.so");

void ADDModelFFWGRVertex::Init() {
  static ClassDocumentation<ADDModelFFWGRVertex> documentation
    ("The ADDModelFFWGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the ADD model.");
  
}

void ADDModelFFWGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr bb,
				      tcPDPtr cc, tcPDPtr) {
  // work out the particles
  int iferm= abs(aa->id());
  int ibos = abs(cc->id());
  Complex coup;
  // overall factor
  assert( ibos >= 22 && ibos <= 24 );
  if( q2last_ != q2 || couplast_ == 0. ) {
    couplast_ = electroMagneticCoupling(q2);
    q2last_ = q2;
  }
  // photon
  if(ibos==22) {
    // alpha
    coup = Complex(UnitRemoval::E * kappa_ * couplast_);
    // _charge of particle
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    coup *= charge_[iferm];
    left (1.);
    right(1.);
  }
  // Z boson
  else if(ibos==23) {
    coup = Complex(UnitRemoval::E * kappa_ * couplast_);
    // _charge of particle
    assert((iferm>=1 && iferm<=6)||(iferm>=11 &&iferm<=16));
    left (gl_[iferm]);
    right(gr_[iferm]);
  }
  else if(ibos==24) {
    coup = Complex(UnitRemoval::E * kappa_ * couplast_) * 
      sqrt(0.5) / sqrt(sin2ThetaW());
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
      left(ckm_[iu-1][id-1]);
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
  // set the coupling
  norm(coup);
}

Complex ADDModelFFWGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
#line 1 "./ADDModelSSGRVertex.cc"
// -*- C++ -*-
//
// ADDModelSSGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelSSGRVertex class.
//

#include "ADDModelSSGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelSSGRVertex::ADDModelSSGRVertex() : kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

void ADDModelSSGRVertex::doinit() {
  addToList(25,25,39);
  SSTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() << "Must have ADDModel in ADDModelSSGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelSSGRVertex,SSTVertex>
describeHerwigADDModelSSGRVertex("Herwig::ADDModelSSGRVertex", "HwADDModel.so");

void ADDModelSSGRVertex::Init() {
  static ClassDocumentation<ADDModelSSGRVertex> documentation
    ("The ADDModelSSGRVertex class is the implementation of"
     " the ADDModel scalar-scalar-graviton vertex");
  
}

void ADDModelSSGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(kappa_ * UnitRemoval::E));
}

Complex ADDModelSSGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
				       Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
#line 1 "./ADDModelVVGRVertex.cc"
// -*- C++ -*-
//
// ADDModelVVGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelVVGRVertex class.
//

#include "ADDModelVVGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void ADDModelVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelVVGRVertex,VVTVertex>
describeHerwigADDModelVVGRVertex("Herwig::ADDModelVVGRVertex", "HwADDModel.so");

void ADDModelVVGRVertex::Init() {
 static ClassDocumentation<ADDModelVVGRVertex> documentation
    ("The ADDModelVVGRVertex class is the implementation"
     " of the ADDModel vector-vector-graviton vertex");
  
}
  
void ADDModelVVGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(UnitRemoval::E * kappa_));
}

ADDModelVVGRVertex::ADDModelVVGRVertex() : kappa_(ZERO), r_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void ADDModelVVGRVertex::doinit() {
  addToList(23,23,39);
  addToList(22,22,39);
  addToList(24,-24,39);
  addToList(21,21,39);
  VVTVertex::doinit();
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD)
    throw Exception() << "Must be ADDModel in ADDModelVVGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

Complex ADDModelVVGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
				       Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
#line 1 "./ADDModelGGGGRVertex.cc"
// -*- C++ -*-
//
// ADDModelGGGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelGGGGRVertex class.
//

#include "ADDModelGGGGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelGGGGRVertex::ADDModelGGGGRVertex() 
  : kappa_(ZERO), r_(ZERO), couplast_(0.), q2last_(ZERO) {
  orderInGem(1);
  orderInGs (1);
  colourStructure(ColourStructure::SU3F);
}

void ADDModelGGGGRVertex::doinit() {
  addToList(21, 21, 21, 39);
  VVVTVertex::doinit();
  // set the graviton coupling 
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() 
      << "Must have ADDModel in ADDModelGGGGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelGGGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << ounit(r_,GeV);
}

void ADDModelGGGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelGGGGRVertex,VVVTVertex>
describeHerwigADDModelGGGGRVertex("Herwig::ADDModelGGGGRVertex", "HwADDModel.so");

void ADDModelGGGGRVertex::Init() {
 static ClassDocumentation<ADDModelGGGGRVertex> documentation
    ("The ADDModelGGGGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}

#ifndef NDEBUG
void ADDModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
#else
void ADDModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				     tcPDPtr, tcPDPtr) {
#endif
  assert(a->id() == ParticleID::g && b->id() ==  ParticleID::g &&
	 c->id() == ParticleID::g);
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = strongCoupling(q2);
    q2last_=q2;
  }
  norm(Complex(couplast_*kappa_*UnitRemoval::E));
}

Complex ADDModelGGGGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
#line 1 "./ADDModelWWWGRVertex.cc"
// -*- C++ -*-
//
// ADDModelWWWGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ADDModelWWWGRVertex class.
//

#include "ADDModelWWWGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

ADDModelWWWGRVertex::ADDModelWWWGRVertex() 
  : kappa_(ZERO), r_(ZERO), couplast_(0.),
    q2last_(ZERO), zfact_(0.) {
  // order in the couplings
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

void ADDModelWWWGRVertex::doinit() {
  addToList(24,-24, 22, 39);
  addToList(24,-24, 23, 39);
  VVVTVertex::doinit();
  zfact_ = sqrt((1.-sin2ThetaW())/sin2ThetaW());
  // set the graviton coupling 
  tcHwADDPtr hwADD=dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() 
      << "Must have ADDModel in ADDModelWWWGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwADD->MPlanckBar();
  r_ = sqr(hwADD->LambdaT())/hwADD->MPlanckBar();
}

void ADDModelWWWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << zfact_ << ounit(r_,GeV);
}
void ADDModelWWWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> zfact_ >> iunit(r_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ADDModelWWWGRVertex,VVVTVertex>
describeHerwigADDModelWWWGRVertex("Herwig::ADDModelWWWGRVertex", "HwADDModel.so");

void ADDModelWWWGRVertex::Init() {
 static ClassDocumentation<ADDModelWWWGRVertex> documentation
    ("The ADDModelWWWGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the WWWGR vertex
void ADDModelWWWGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
  int ida=a->id();
  int idb=b->id();
  int idc=c->id();
  // first the overall normalisation
  if(q2!=q2last_||couplast_==0.) {
    couplast_ = electroMagneticCoupling(q2);
    q2last_=q2;
  }
  // W- W+ photon and cylic perms
  if((ida==-24 && idb== 24 && idc== 22) ||
     (ida== 22 && idb==-24 && idc== 24) || 
     (ida== 24 && idb== 22 && idc==-24) )
    norm(Complex(couplast_*kappa_*UnitRemoval::E));
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) ||
	  (ida== 22 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 22 && idc== 24) )
    norm(-Complex(couplast_*kappa_*UnitRemoval::E));
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) ||
	  (ida== 23 && idb==-24 && idc== 24) || 
	  (ida== 24 && idb== 23 && idc==-24) )
    norm(Complex(couplast_*zfact_*kappa_*UnitRemoval::E));
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) ||
	  (ida== 23 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 23 && idc== 24) )
    norm(-Complex(couplast_*zfact_*kappa_*UnitRemoval::E));
  else assert(false);
}

Complex ADDModelWWWGRVertex::propagator(int iopt, Energy2 q2,tcPDPtr part,
					Energy mass, Energy width) {
  if(part->id()!=ParticleID::Graviton)
    return VertexBase::propagator(iopt,q2,part,mass,width);
  else
    return Complex(4.*Constants::pi*UnitRemoval::E2/sqr(r_));
}
#line 1 "./GravitonMassGenerator.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GravitonMassGenerator class.
//

#include "GravitonMassGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/GenericWidthGenerator.h"
#include "ADDModel.h"

using namespace Herwig;

GravitonMassGenerator::GravitonMassGenerator()
  : prefactor_(0.), delta_(2), md_(1000.*GeV), mMin_(MeV) 
{}

IBPtr GravitonMassGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GravitonMassGenerator::fullclone() const {
  return new_ptr(*this);
}

void GravitonMassGenerator::persistentOutput(PersistentOStream & os) const {
  os << prefactor_ << delta_ << ounit(md_,GeV) << ounit(mMin_,GeV);
}

void GravitonMassGenerator::persistentInput(PersistentIStream & is, int) {
  is >> prefactor_ >> delta_ >> iunit(md_,GeV) >> iunit(mMin_,GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<GravitonMassGenerator,GenericMassGenerator>
describeHerwigGravitonMassGenerator("Herwig::GravitonMassGenerator", "HwADDModel.so");

void GravitonMassGenerator::Init() {

  static ClassDocumentation<GravitonMassGenerator> documentation
    ("The GravitonMassGenerator class generates the mass for external gravitions "
     "in the ADD model.");

  static Parameter<GravitonMassGenerator,Energy> interfaceMinimumMass
    ("MinimumMass",
     "Minimum gravition mass to avoid numerical problems",
     &GravitonMassGenerator::mMin_, GeV, MeV, MeV, GeV,
     false, false, Interface::limited);

}

void GravitonMassGenerator::doinit() {
  GenericMassGenerator::doinit();
  tcHwADDPtr hwADD = dynamic_ptr_cast<tcHwADDPtr>(generator()->standardModel());
  if(!hwADD) 
    throw Exception() << "Must have ADDModel in GravitonMassGenerator::doinit()"
		      << Exception::runerror;
  delta_ = hwADD->delta();
  md_ =  hwADD->MD();
  // calculate the prefactor
  prefactor_ = sqr(hwADD->MPlanckBar()/md_);
  // even no of dimensions
  if(delta_%2==0) {
    unsigned int n = delta_/2;
    prefactor_ *= 2.*pow(Constants::pi,int(n));
    for(unsigned int ix=1;ix<n;++ix) {
      prefactor_ /= double(ix);
    }
  }
  // odd number of dimensions
  else {
    unsigned int n = (delta_-1)/2;
    prefactor_ *= 2.*pow(Constants::pi,int(n));
    for(unsigned int ix=0;ix<n;++ix) {
      prefactor_ /= double(ix)+0.5;
    }
  }
}

Energy GravitonMassGenerator::mass(double & wgt, const ParticleData & ,
				   const Energy low,const Energy upp, int,
				   double r) const {
  Energy low2 = max(mMin_,low);
  double rlow = pow(double(low2/md_),int(delta_))/double(delta_);
  double rupp = pow(double(upp /md_),int(delta_))/double(delta_);
  double rho = rlow + (rupp-rlow)*r;
  wgt = (rupp-rlow)*prefactor_;
  return pow(double(delta_)*rho,1./double(delta_))*md_;
}
