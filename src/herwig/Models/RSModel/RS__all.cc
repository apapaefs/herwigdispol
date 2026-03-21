#line 1 "./RSModel.cc"
// -*- C++ -*-
//
// RSModel.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModel class.
//

#include "RSModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

void RSModel::doinit() {
  addVertex(FFGRVertex_);
  addVertex(VVGRVertex_);
  addVertex(SSGRVertex_);
  addVertex(FFGGRVertex_);
  addVertex(FFWGRVertex_);
  addVertex(GGGGRVertex_);
  addVertex(WWWGRVertex_);
  BSMModel::doinit();
}

void RSModel::persistentOutput(PersistentOStream & os) const {
  os << ounit(Lambda_pi_,GeV) 
     << FFGRVertex_ << VVGRVertex_ << SSGRVertex_ 
     << FFGGRVertex_ << FFWGRVertex_ 
     << GGGGRVertex_ << WWWGRVertex_;
}

void RSModel::persistentInput(PersistentIStream & is, int) {
  is >> iunit(Lambda_pi_,GeV) 
     >> FFGRVertex_ >> VVGRVertex_ >> SSGRVertex_
     >> FFGGRVertex_ >> FFWGRVertex_ 
     >> GGGGRVertex_ >> WWWGRVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModel,BSMModel>
describeHerwigRSModel("Herwig::RSModel", "HwRSModel.so");

void RSModel::Init() {
  

  static Reference<RSModel,ThePEG::Helicity::AbstractFFTVertex> interfaceVertexFFGR
    ("Vertex/FFGR",
     "Reference to the fermion-fermion-graviton vertex",
     &RSModel::FFGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractVVTVertex> interfaceVertexVVGR
    ("Vertex/VVGR",
     "Reference to the vector-vector-graviton vertex",
     &RSModel::VVGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractSSTVertex> interfaceVertexSSGR
    ("Vertex/SSGR",
     "Reference to the scalar-scalar-graviton vertex",
     &RSModel::SSGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFGGR
    ("Vertex/FFGGR",
     "Reference to the fermion-antifermion-gluon graviton vertex",
     &RSModel::FFGGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractFFVTVertex> interfaceVertexFFWGR
    ("Vertex/FFWGR",
     "Reference to the fermion-antifermion-weak vector boson graviton vertex",
     &RSModel::FFWGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexGGGGR
    ("Vertex/GGGGR",
     "Reference to the three gluon graviton vertex",
     &RSModel::GGGGRVertex_, false, false, true, false, false);
  
  static Reference<RSModel,ThePEG::Helicity::AbstractVVVTVertex> interfaceVertexWWWGR
    ("Vertex/WWWGR",
     "Reference to the three weak vector boson graviton vertex",
     &RSModel::WWWGRVertex_, false, false, true, false, false);
  
  static Parameter<RSModel,Energy> interfaceLambda_pi
    ("Lambda_pi",
     "The coupling of the graviton to matter",
     &RSModel::Lambda_pi_, GeV, 10000*GeV, ZERO, 1.0e12*GeV,
     false, false, false);
  
  static ClassDocumentation<RSModel> documentation
    ("The RSModel class replaces the Standard Model class for the"
     " RS model",
     "The Randall-Sundrum model was constructed from \\cite{Randall:1999ee}.",
     "%\\cite{Randall:1999ee}\n"
     "\\bibitem{Randall:1999ee}\n"
     "  L.~Randall and R.~Sundrum,\n"
     "  ``A large mass hierarchy from a small extra dimension,''\n"
     "  Phys.\\ Rev.\\ Lett.\\  {\\bf 83}, 3370 (1999)\n"
     "  [arXiv:hep-ph/9905221].\n"
     "  %%CITATION = PRLTA,83,3370;%%\n"
     );
  
}
#line 1 "./RSModelFFGRVertex.cc"
// -*- C++ -*-
//
// RSModelFFGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFGRVertex class.
//

#include "RSModelFFGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void RSModelFFGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelFFGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelFFGRVertex,FFTVertex>
describeHerwigRSModelFFGRVertex("Herwig::RSModelFFGRVertex", "HwRSModel.so");

void RSModelFFGRVertex::Init() {
  static ClassDocumentation<RSModelFFGRVertex> documentation
    ("The RSModelFFGRVertex class is the RSModel calculation"
     " of the fermion-antifermion-graviton vertex");
  
}
  
void RSModelFFGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(kappa_ * UnitRemoval::E));
}

RSModelFFGRVertex::RSModelFFGRVertex() : kappa_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void RSModelFFGRVertex::doinit() {
  // PDG codes for the particles
  // the quarks
  for (int ix=1;ix<7;++ix) addToList(-ix,ix,39);
  // the leptons
  for (int ix=11;ix<17;++ix) addToList(-ix,ix,39);
  FFTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS)
    throw Exception() << "Must have RSModel in RSModelFFGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}
#line 1 "./RSModelFFGGRVertex.cc"
// -*- C++ -*-
//
// RSModelFFGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFGGRVertex class.
//

#include "RSModelFFGGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelFFGGRVertex::RSModelFFGGRVertex() 
  : couplast_(0.), q2last_(ZERO), kappa_(ZERO) {
  orderInGem(1);
  orderInGs (1);
  colourStructure(ColourStructure::SU3TFUND);
}

void RSModelFFGGRVertex::doinit() {
  for(int ix=1;ix<7;++ix) {
    addToList(-ix,ix,21,39);
  }
  FFVTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) throw Exception() 
	      << "Must have RSModel in RSModelFFGGRVertex::doinit()"
	      << Exception::runerror;
  kappa_ = 2./hwRS->lambda_pi();
}

void RSModelFFGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelFFGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelFFGGRVertex,FFVTVertex>
describeHerwigRSModelFFGGRVertex("Herwig::RSModelFFGGRVertex", "HwRSModel.so");

void RSModelFFGGRVertex::Init() {
  static ClassDocumentation<RSModelFFGGRVertex> documentation
    ("The RSModelFFGGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the RS model.");
  
}

// FFGGR coupling
#ifndef NDEBUG
void RSModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr,
				     tcPDPtr cc, tcPDPtr) {
#else
void RSModelFFGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				      tcPDPtr, tcPDPtr) {
#endif
  // work out the particles
  assert(cc->id()==ParticleID::g && abs(aa->id()) <=6 );
  // overall factor
  if(q2last_ != q2 || couplast_ == 0. ) {
    couplast_ = strongCoupling(q2);
    q2last_ = q2;
  }
  left (1.);
  right(1.);
  // set the coupling
  norm( UnitRemoval::E * kappa_ * couplast_);
}
#line 1 "./RSModelFFWGRVertex.cc"
// -*- C++ -*-
//
// RSModelFFWGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelFFWGRVertex class.
//

#include "RSModelFFWGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Models/StandardModel/StandardCKM.h"

using namespace Herwig;
using namespace ThePEG;

RSModelFFWGRVertex::RSModelFFWGRVertex() 
  : charge_(17,0.), gl_(17,0.), gr_(17,0.),
    ckm_(3,vector<Complex>(3,0.0)),
    couplast_(0.), q2last_(ZERO), kappa_(ZERO) {
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void RSModelFFWGRVertex::doinit() {
  for(int ix=11;ix<17;++ix) {
    addToList(-ix,ix,22,39);
    addToList(-ix,ix,23,39);
  }
  for(int ix=1;ix<7;++ix) {
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
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) throw Exception() 
	      << "Must have RSModel in RSModelFFWGRVertex::doinit()"
	      << Exception::runerror;
  double sw2 = sin2ThetaW();
  double fact = 0.25/sqrt(sw2*(1.-sw2));
  for(int ix=1;ix<4;++ix) {
    charge_[2*ix-1]  = hwRS->ed();
    charge_[2*ix ]   = hwRS->eu();
    charge_[2*ix+9 ] = hwRS->ee();
    charge_[2*ix+10] = hwRS->enu();
    gl_[2*ix-1]  = fact*(hwRS->vd()  + hwRS->ad() );
    gl_[2*ix ]   = fact*(hwRS->vu()  + hwRS->au() );
    gl_[2*ix+9 ] = fact*(hwRS->ve()  + hwRS->ae() );
    gl_[2*ix+10] = fact*(hwRS->vnu() + hwRS->anu());
    gr_[2*ix-1]  = fact*(hwRS->vd()  - hwRS->ad() );
    gr_[2*ix ]   = fact*(hwRS->vu()  - hwRS->au() );
    gr_[2*ix+9 ] = fact*(hwRS->ve()  - hwRS->ae() );
    gr_[2*ix+10] = fact*(hwRS->vnu() - hwRS->anu());
  }
  kappa_ = 2./hwRS->lambda_pi();
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

void RSModelFFWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << charge_ << gl_ << gr_ << ounit(kappa_,InvGeV) << ckm_;
}

void RSModelFFWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> charge_ >> gl_ >> gr_ >> iunit(kappa_,InvGeV) >> ckm_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelFFWGRVertex,FFVTVertex>
describeHerwigRSModelFFWGRVertex("Herwig::RSModelFFWGRVertex", "HwRSModel.so");

void RSModelFFWGRVertex::Init() {
  static ClassDocumentation<RSModelFFWGRVertex> documentation
    ("The RSModelFFWGRVertexxs class is the implementation"
     " of the two fermion vector coupling for the RS model.");
  
}

// FFWGR coupling
void RSModelFFWGRVertex::setCoupling(Energy2 q2,tcPDPtr aa,tcPDPtr bb,
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
#line 1 "./RSModelSSGRVertex.cc"
// -*- C++ -*-
//
// RSModelSSGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelSSGRVertex class.
//

#include "RSModelSSGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelSSGRVertex::RSModelSSGRVertex() : kappa_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

void RSModelSSGRVertex::doinit() {
  addToList(25,25,39);
  SSTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) 
    throw Exception() << "Must have RSModel in RSModelSSGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}

void RSModelSSGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelSSGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelSSGRVertex,SSTVertex>
describeHerwigRSModelSSGRVertex("Herwig::RSModelSSGRVertex", "HwRSModel.so");

void RSModelSSGRVertex::Init() {
  static ClassDocumentation<RSModelSSGRVertex> documentation
    ("The RSModelSSGRVertex class is the implementation of"
     " the RSModel scalar-scalar-graviton vertex");
  
}

void RSModelSSGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
    norm(Complex(kappa_ * UnitRemoval::E));
}
#line 1 "./RSModelVVGRVertex.cc"
// -*- C++ -*-
//
// RSModelVVGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelVVGRVertex class.
//

#include "RSModelVVGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

void RSModelVVGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}

void RSModelVVGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelVVGRVertex,VVTVertex>
describeHerwigRSModelVVGRVertex("Herwig::RSModelVVGRVertex", "HwRSModel.so");

void RSModelVVGRVertex::Init() {
 static ClassDocumentation<RSModelVVGRVertex> documentation
    ("The RSModelVVGRVertex class is the implementation"
     " of the RSModel vector-vector-graviton vertex");
  
}
  
void RSModelVVGRVertex::setCoupling(Energy2,tcPDPtr,tcPDPtr, tcPDPtr) {
  norm(Complex(UnitRemoval::E * kappa_));
}

RSModelVVGRVertex::RSModelVVGRVertex() : kappa_(ZERO) {
  orderInGem(1);
  orderInGs (0);
  colourStructure(ColourStructure::DELTA);
}

void RSModelVVGRVertex::doinit() {
  addToList(23,23,39);
  addToList(22,22,39);
  addToList(24,-24,39);
  addToList(21,21,39);
  VVTVertex::doinit();
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS)
    throw Exception() << "Must be RSModel in RSModelVVGRVertex::doinit()"
		      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}
#line 1 "./RSModelWWWGRVertex.cc"
// -*- C++ -*-
//
// RSModelVVVGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelWWWGRVertex class.
//

#include "RSModelWWWGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelWWWGRVertex::RSModelWWWGRVertex() 
  : kappa_(ZERO), _couplast(0.), 
    _q2last(ZERO), _zfact(0.) {
  // order in the couplings
  orderInGem(2);
  orderInGs (0);
  colourStructure(ColourStructure::SINGLET);
}

void RSModelWWWGRVertex::doinit() {
  addToList(24,-24, 22, 39);
  addToList(24,-24, 23, 39);
  VVVTVertex::doinit();
  _zfact = sqrt((1.-sin2ThetaW())/sin2ThetaW());
  // set the graviton coupling 
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) 
    throw Exception() 
      << "Must have RSModel in RSModelWWWGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}

void RSModelWWWGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV) << _zfact;
}
void RSModelWWWGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV) >> _zfact;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelWWWGRVertex,VVVTVertex>
describeHerwigRSModelWWWGRVertex("Herwig::RSModelWWWGRVertex", "HwRSModel.so");

void RSModelWWWGRVertex::Init() {
 static ClassDocumentation<RSModelWWWGRVertex> documentation
    ("The RSModelWWWGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the WWWGR vertex
void RSModelWWWGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
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
     (ida== 24 && idb== 22 && idc==-24) )
    norm(Complex(_couplast*kappa_*UnitRemoval::E));
  // W+ W- photon (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 22) ||
	  (ida== 22 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 22 && idc== 24) )
    norm(-Complex(_couplast*kappa_*UnitRemoval::E));
  // W- W+ Z and cylic perms
  else if((ida==-24 && idb== 24 && idc== 23) ||
	  (ida== 23 && idb==-24 && idc== 24) || 
	  (ida== 24 && idb== 23 && idc==-24) )
    norm(Complex(_couplast*_zfact*kappa_*UnitRemoval::E));
  // W+ W- Z (anticylic perms of above)
  else if((ida== 24 && idb==-24 && idc== 23) ||
	  (ida== 23 && idb== 24 && idc==-24) || 
	  (ida==-24 && idb== 23 && idc== 24) )
    norm(-Complex(_couplast*_zfact*kappa_*UnitRemoval::E));
  else assert(false);
}
#line 1 "./RSModelGGGGRVertex.cc"
// -*- C++ -*-
//
// RSModelGGGGRVertex.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSModelGGGGRVertex class.
//

#include "RSModelGGGGRVertex.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

RSModelGGGGRVertex::RSModelGGGGRVertex() 
  : kappa_(ZERO), _couplast(0.), _q2last(ZERO) {
  orderInGem(1);
  orderInGs (1);
  colourStructure(ColourStructure::SU3FF);
}

void RSModelGGGGRVertex::doinit() {
  addToList(21, 21, 21, 39);
  VVVTVertex::doinit();
  // set the graviton coupling 
  tcHwRSPtr hwRS=dynamic_ptr_cast<tcHwRSPtr>(generator()->standardModel());
  if(!hwRS) 
    throw Exception() 
      << "Must have RSModel in RSModelGGGGRVertex::doinit()"
      << Exception::runerror;
  kappa_=2./hwRS->lambda_pi();
}

void RSModelGGGGRVertex::persistentOutput(PersistentOStream & os) const {
  os << ounit(kappa_,InvGeV);
}
void RSModelGGGGRVertex::persistentInput(PersistentIStream & is, int) {
  is >> iunit(kappa_,InvGeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RSModelGGGGRVertex,VVVTVertex>
describeHerwigRSModelGGGGRVertex("Herwig::RSModelGGGGRVertex", "HwRSModel.so");

void RSModelGGGGRVertex::Init() {
 static ClassDocumentation<RSModelGGGGRVertex> documentation
    ("The RSModelGGGGRVertex class is the four point coupling"
     " of three vector bosons and a graviton in the Randell-Sundrum model.");
  
}


// couplings for the GGGGR vertex
#ifndef NDEBUG
void RSModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr a,tcPDPtr b,
				     tcPDPtr c, tcPDPtr) {
#else
void RSModelGGGGRVertex::setCoupling(Energy2 q2,tcPDPtr,tcPDPtr,
				     tcPDPtr, tcPDPtr) {
#endif
  assert(a->id() == ParticleID::g && b->id() ==  ParticleID::g &&
	 c->id() == ParticleID::g);
  if(q2!=_q2last||_couplast==0.) {
    _couplast = strongCoupling(q2);
    _q2last=q2;
  }
  norm(Complex(_couplast*kappa_*UnitRemoval::E));
}
