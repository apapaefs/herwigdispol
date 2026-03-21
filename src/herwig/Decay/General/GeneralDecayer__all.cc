#line 1 "./FFSDecayer.cc"
// -*- C++ -*-
//
// FFSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFSDecayer class.
//

#include "FFSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FFSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractFFSVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<FFSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(inV.at(inter));
    outgoingVertexF_[inter] = AbstractFFVVertexPtr();
    outgoingVertexS_[inter] = AbstractVSSVertexPtr();
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexS_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
      else
	outgoingVertexS_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
    }
  }
}

void FFSDecayer::persistentOutput(PersistentOStream & os) const {
  os << perturbativeVertex_       << vertex_
     << incomingVertex_   << outgoingVertexF_
     << outgoingVertexS_;
}

void FFSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> perturbativeVertex_       >> vertex_
     >> incomingVertex_   >> outgoingVertexF_
     >> outgoingVertexS_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFSDecayer,GeneralTwoBodyDecayer>
describeHerwigFFSDecayer("Herwig::FFSDecayer", "Herwig.so");

void FFSDecayer::Init() {

  static ClassDocumentation<FFSDecayer> documentation
    ("The FFSDecayer class implements the decay of a fermion to "
     "a fermion and a scalar.");

}

void FFSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  int itype[2];
  if(part.dataPtr()->CC())        itype[0] = part.id()    > 0? 0:1;
  else                              itype[0] = 2;
  if(decay[0]->dataPtr()->CC())     itype[1] = decay[0]->id() > 0? 0:1;
  else                              itype[1] = 2;
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  // for the decaying particle
  if(ferm) {
    SpinorWaveFunction::
      constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorWaveFunction::constructSpinInfo(wave_,decay[0],outgoing,true);
  }
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

double FFSDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  int itype[2];
  if(part.dataPtr()->CC())  itype[0] = part.id()    > 0? 0:1;
  else                      itype[0] = 2;
  if(outgoing[0]->CC())     itype[1] = outgoing[0]->id() > 0? 0:1;
  else                              itype[1] = 2;
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));

  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,momenta[0],outgoing[0],Helicity::outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(wave_   ,momenta[0],outgoing[0],Helicity::outgoing);
  ScalarWaveFunction scal(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      if(ferm) (*ME())(if1, if2, 0) = 0.;
      else     (*ME())(if2, if1, 0) = 0.;
      for(auto vert : vertex_) {
	if(ferm) (*ME())(if1, if2, 0) += 
		   vert->evaluate(scale,wave_[if1],wavebar_[if2],scal);
	else     (*ME())(if2, if1, 0) += 
		   vert->evaluate(scale,wave_[if1],wavebar_[if2],scal);
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy FFSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    double mu1(0.),mu2(0.);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin() == PDT::Spin1Half) {
      mu1 = outa.second/inpart.second;
      mu2 = outb.second/inpart.second;
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, outa.first, outb.first);
    }
    else {
      mu1 = outb.second/inpart.second;
      mu2 = outa.second/inpart.second;
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, outb.first, outa.first);
      
    }
    double c2 = norm(perturbativeVertex_[0]->norm());
    Complex cl = perturbativeVertex_[0]->left();
    Complex cr = perturbativeVertex_[0]->right();
    double me2 = c2*( (norm(cl) + norm(cr))*(1. + sqr(mu1) - sqr(mu2))
		      + 2.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/16./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double FFSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  int iscal (0), iferm (1), iglu (2);
  // get location of outgoing fermion/scalar
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,iferm);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(false);
  if(itype[0] == itype[1] ) {
    ferm = itype[0]==0 || (itype[0]==2 && decay[iscal]->id() < 0);
  }
  else if(itype[0] == 2) {
    ferm = itype[1]==0;
  }
  else if(itype[1] == 2) {
    ferm = itype[0]==0;
  }
  else if((itype[0] == 1 && itype[1] == 0) ||
	  (itype[0] == 0 && itype[1] == 1)) {
    if(abs(inpart.id())<=16) {
      ferm = itype[0]==0;
    }
    else if(abs(decay[iferm]->id())<=16) {
      ferm = itype[1]==0;
    }
    else {
      ferm = true;
    }
  }
  else
    assert(false);
  if(meopt==Initialize) {
    // create spinor (bar) for decaying particle
    if(ferm) {
      SpinorWaveFunction::calculateWaveFunctions(wave3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
						 incoming);
      if(wave3_[0].wave().Type() != SpinorType::u)
   	for(unsigned int ix = 0; ix < 2; ++ix) wave3_[ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_,rho3_, const_ptr_cast<tPPtr>(&inpart), 
						    incoming);
      if(wavebar3_[0].wave().Type() != SpinorType::v)
   	for(unsigned int ix = 0; ix < 2; ++ix) wavebar3_[ix].conjugate();
    }
  }
  // setup spin information when needed 
  if(meopt==Terminate) {
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(wave3_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar3_,decay[iferm],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar3_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(wave3_,decay[iferm],outgoing,true);
    }
    ScalarWaveFunction::constructSpinInfo(        decay[iscal],outgoing,true);
    VectorWaveFunction::constructSpinInfo(gluon_, decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calulate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin0,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(wave3_   , decay[iferm],outgoing);
  
  ScalarWaveFunction swave3_(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_,   decay[iglu ],outgoing,true);

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),decay[iglu ]->dataPtr(),10,
					  outgoing));
    }
  }
#endif
  
  if (! ((incomingVertex_[inter]  && (outgoingVertexF_[inter] || outgoingVertexS_[inter])) ||
	 (outgoingVertexF_[inter] &&  outgoingVertexS_[inter])))
    throw Exception()
      << "Invalid vertices for radiation in FFS decay in FFSDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int F(1), S(2);
  if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 && 
      decay[iferm]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,S);
  else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar && 
	   decay[iscal]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,S);


  Complex diag;
  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
      for(unsigned int ig = 0; ig < 2; ++ig) {
   	// radiation from the incoming fermion
   	if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
   	  assert(incomingVertex_[inter]);
	  if (ferm) {
	    SpinorWaveFunction spinorInter =
	      incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wave3_[ifi],
					       gluon_[2*ig],inpart.mass());

	    assert(wave3_[ifi].particle()->id()==spinorInter.particle()->id());
	    diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,spinorInter,wavebar3_[ifo],swave3_);
	  }
	  else {
	    SpinorBarWaveFunction spinorBarInter = 
	      incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wavebar3_[ifi],
					       gluon_[2*ig],inpart.mass());

	    assert(wavebar3_[ifi].particle()->id()==spinorBarInter.particle()->id());
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag+= vertex->evaluate(scale,wave3_[ifo], spinorBarInter,swave3_);
	  }
	  if(!couplingSet) {
	    gs = abs(incomingVertex_[inter]->norm());
	    couplingSet = true;
	  }
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(ifi, 0, ifo, ig) += 
	       colourFlow[0][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
	}
	  
  	// radiation from outgoing fermion
   	if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
  	  assert(outgoingVertexF_[inter]);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[iferm]->dataPtr();
	  if(off->CC()) off = off->CC();  	  
	  if (ferm) {	    
	    SpinorBarWaveFunction spinorBarInter = 
	      outgoingVertexF_[inter]->evaluate(scale,3,off,wavebar3_[ifo],
						gluon_[2*ig],decay[iferm]->mass());
	    
	    assert(wavebar3_[ifo].particle()->id()==spinorBarInter.particle()->id());
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag+= vertex->evaluate(scale,wave3_[ifi],spinorBarInter,swave3_);
	  }
	  else {
	    SpinorWaveFunction spinorInter = 
	      outgoingVertexF_[inter]->evaluate(scale,3,off,wave3_[ifo],
						gluon_[2*ig],decay[iferm]->mass());
	      
	    assert(wave3_[ifo].particle()->id()==spinorInter.particle()->id());
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag+= vertex->evaluate(scale,spinorInter,wavebar3_[ifi],swave3_);
	  }
	  if(!couplingSet) {
	    gs = abs(outgoingVertexF_[inter]->norm());
	    couplingSet = true;
	  }
	  for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	    (*ME[colourFlow[F][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[F][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
  	}

  	// radiation from outgoing scalar
   	if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
  	  assert(outgoingVertexS_[inter]);
	  // ensure you get correct ougoing particle from first vertex
	  tcPDPtr off = decay[iscal]->dataPtr();
	  if(off->CC()) off = off->CC();
	  ScalarWaveFunction  scalarInter = 
	    outgoingVertexS_[inter]->evaluate(scale,3,off,gluon_[2*ig],
					      swave3_,decay[iscal]->mass());
	    
	  assert(swave3_.particle()->id()==scalarInter.particle()->id());
	  if (ferm){
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag += vertex->evaluate(scale,wave3_[ifi],wavebar3_[ifo],scalarInter);
	  }
	  else {
	    diag = 0.;
	    for(auto vertex :vertex_)
	      diag += vertex->evaluate(scale,wave3_[ifo],wavebar3_[ifi],scalarInter);
	  }
	  if(!couplingSet) {
	    gs = abs(outgoingVertexS_[inter]->norm());
	    couplingSet = true;
	  }
	  for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
  	    (*ME[colourFlow[S][ix].first])(ifi, 0, ifo, ig) += 
	      colourFlow[S][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
  	}
      }
    }
  }

  // contract matrices
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif

  // return the answer
  return output;
}
#line 1 "./FFVDecayer.cc"
// -*- C++ -*-
//
// FFVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVDecayer class.
//

#include "FFVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractFFVVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<FFVVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(inV.at(inter));
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexV_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::FFV)
	outgoingVertexF_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
      else
	outgoingVertexV_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    }
  }
}

void FFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_
     << incomingVertex_   << outgoingVertexF_
     << outgoingVertexV_;
}

void FFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertexF_
     >> outgoingVertexV_;
}

void FFVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // type of process
  int itype[2];
  if(part.dataPtr()->CC())      itype[0] = part.id() > 0 ? 0 : 1;
  else                          itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                          itype[1] = 2;  
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  // for the decaying particle
  if(ferm) {
    SpinorWaveFunction::
      constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorWaveFunction::constructSpinInfo(wave_,decay[0],outgoing,true);
  }
  VectorWaveFunction::
    constructSpinInfo(vector_,decay[1],outgoing,true,false);
}

double FFVDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // type of process
  int itype[2];
  if(part.dataPtr()->CC()) itype[0] = part.id() > 0 ? 0 : 1;
  else                     itype[0] = 2;
  if(outgoing[0]->CC())    itype[1] = outgoing[0]->id() > 0 ? 0 : 1;
  else                     itype[1] = 2;  
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  Energy2 scale(sqr(part.mass()));
  if(ferm)
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,momenta[0],outgoing[0],Helicity::outgoing);
  else
    SpinorWaveFunction::
      calculateWaveFunctions(wave_   ,momenta[0],outgoing[0],Helicity::outgoing);
  bool massless = outgoing[1]->mass()==ZERO;
  VectorWaveFunction::
    calculateWaveFunctions(vector_,momenta[1],outgoing[1],Helicity::outgoing,massless);
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 2; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(massless && vhel == 1) ++vhel;
	if(ferm)
	  (*ME())(if1, if2,vhel) = 0.;
	else
	  (*ME())(if2, if1, vhel) = 0.;
	for(auto vertex : vertex_) {
	  if(ferm)
	    (*ME())(if1, if2,vhel) += 
	      vertex->evaluate(scale,wave_[if1],wavebar_[if2],vector_[vhel]);
	  else
	    (*ME())(if2, if1, vhel) += 
	      vertex->evaluate(scale,wave_[if1],wavebar_[if2],vector_[vhel]);
	}
      }
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy FFVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outa.first->iSpin() == PDT::Spin1Half)
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in,
				       outa.first, outb.first);
    else {
      swap(mu1,mu2);
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second),in,
				       outb.first,outa.first);
    }
    Complex cl(perturbativeVertex_[0]->left()),cr(perturbativeVertex_[0]->right());
    double me2(0.);
    if( mu2 > 0. ) {
      me2 = (norm(cl) + norm(cr))*(1. + sqr(mu1*mu2) + sqr(mu2) 
				   - 2.*sqr(mu1) - 2.*sqr(mu2*mu2) 
				   +  sqr(mu1*mu1))
	- 6.*mu1*sqr(mu2)*(conj(cl)*cr + conj(cr)*cl).real();
      me2 /= sqr(mu2);
    }
    else
      me2 = 2.*( (norm(cl) + norm(cr))*(sqr(mu1) + 1.) 
		 - 4.*mu1*(conj(cl)*cr + conj(cr)*cl).real() );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm/16./Constants::pi; 
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer 
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFVDecayer,GeneralTwoBodyDecayer>
describeHerwigFFVDecayer("Herwig::FFVDecayer", "Herwig.so");

void FFVDecayer::Init() {

  static ClassDocumentation<FFVDecayer> documentation
    ("The FFVDecayer class implements the decay of a fermion to a fermion and a vector boson");

}

double  FFVDecayer::threeBodyME(const int , const Particle & inpart,
				const ParticleVector & decay,
				ShowerInteraction inter, MEOption meopt) {

  int iferm (0), ivect (1), iglu (2);
  // get location of outgoing lepton/vector
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin1Half) swap(iferm,ivect);
  // work out whether inpart is a fermion or antifermion
  int itype[2];
  if(inpart.dataPtr()->CC())        itype[0] = inpart.id() > 0 ? 0 : 1;
  else                              itype[0] = 2;
  if(decay[iferm]->dataPtr()->CC()) itype[1] = decay[iferm]->id() > 0 ? 0 : 1;
  else                              itype[1] = 2;

  bool ferm(itype[0] == 0 || itype[1] == 0 || 
	   (itype[0] == 2 && itype[1] == 2 && decay[ivect]->id() < 0));  

  // no emissions from massive vectors
  bool massless = decay[ivect]->dataPtr()->mass()==ZERO;

  if(meopt==Initialize) {
    // create spinor (bar) for decaying particle
    if(ferm) {
      SpinorWaveFunction::calculateWaveFunctions(wave3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
						 incoming);
      if(wave3_[0].wave().Type() != SpinorType::u)
   	for(unsigned int ix = 0; ix < 2; ++ix) wave3_[ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_,rho3_, const_ptr_cast<tPPtr>(&inpart), 
						    incoming);
      if(wavebar3_[0].wave().Type() != SpinorType::v)
   	for(unsigned int ix = 0; ix < 2; ++ix) wavebar3_[ix].conjugate();
    }
  }
  // setup spin information when needed 
  if(meopt==Terminate) {
    if(ferm) {
      SpinorWaveFunction::
	constructSpinInfo(wave3_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorBarWaveFunction::constructSpinInfo(wavebar3_,decay[iferm],outgoing,true);
    }
    else {
      SpinorBarWaveFunction::
	constructSpinInfo(wavebar3_,const_ptr_cast<tPPtr>(&inpart),incoming,true);
      SpinorWaveFunction::constructSpinInfo(wave3_,decay[iferm],outgoing,true);
    }
    VectorWaveFunction::constructSpinInfo(vector3_, decay[ivect],outgoing,true,massless);
    VectorWaveFunction::constructSpinInfo(gluon_,   decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calulate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, PDT::Spin1Half,
								       PDT::Spin1,     PDT::Spin1)));

  // create wavefunctions
  if (ferm)  SpinorBarWaveFunction::calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  else       SpinorWaveFunction::   calculateWaveFunctions(wave3_   , decay[iferm],outgoing);
  
  VectorWaveFunction::calculateWaveFunctions(vector3_, decay[ivect],outgoing,massless);
  VectorWaveFunction::calculateWaveFunctions(gluon_,   decay[iglu ],outgoing,true );

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  					  decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  if (! ((incomingVertex_[inter]  && (outgoingVertexF_[inter]  || outgoingVertexV_[inter])) ||
	 (outgoingVertexF_[inter] &&  outgoingVertexV_[inter])))
    throw Exception()
      << "Invalid vertices for QCD radiation in FFV decay in FFVDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int F(1), V(2);
  if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar && 
      decay[ivect]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,V);
  else if (decay[ivect]->dataPtr()->iColour()==PDT::Colour3 && 
	   decay[iferm]->dataPtr()->iColour()==PDT::Colour8)
    swap(F,V);

  Complex diag;
  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ifi = 0; ifi < 2; ++ifi) {
    for(unsigned int ifo = 0; ifo < 2; ++ifo) {
      for(unsigned int iv = 0; iv < 3; ++iv) {
	for(unsigned int ig = 0; ig < 2; ++ig) {
	  // radiation from the incoming fermion
	  if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(incomingVertex_[inter]);
	    if (ferm){
	      SpinorWaveFunction spinorInter =
		incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wave3_[ifi],
						  gluon_[2*ig],inpart.mass());
	      
	      assert(wave3_[ifi].particle()->id()==spinorInter.particle()->id());
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,spinorInter,wavebar3_[ifo],vector3_[iv]);
	    }
	    else {
	      SpinorBarWaveFunction spinorBarInter = 
		incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),wavebar3_[ifi],
						  gluon_[2*ig],inpart.mass());
	      
	      assert(wavebar3_[ifi].particle()->id()==spinorBarInter.particle()->id());
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifo], spinorBarInter,vector3_[iv]);
	    }
	    if(!couplingSet) {
	      gs = abs(incomingVertex_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	      (*ME[colourFlow[0][ix].first])(ifi, ifo, iv, ig) += 
		colourFlow[0][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	  // radiation from outgoing fermion
	  if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexF_[inter]);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[iferm]->dataPtr();
	    if(off->CC()) off = off->CC(); 
	    if (ferm) {	    
	      SpinorBarWaveFunction spinorBarInter = 
		outgoingVertexF_[inter]->evaluate(scale,3,off,wavebar3_[ifo],
						  gluon_[2*ig],decay[iferm]->mass());
	      
	      assert(wavebar3_[ifo].particle()->id()==spinorBarInter.particle()->id());
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifi],spinorBarInter,vector3_[iv]);
	    }
	    else {
	      SpinorWaveFunction spinorInter = 
		outgoingVertexF_[inter]->evaluate(scale,3,off,wave3_[ifo],
						  gluon_[2*ig],decay[iferm]->mass());
		
	      assert(wave3_[ifo].particle()->id()==spinorInter.particle()->id());
	      
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,spinorInter,wavebar3_[ifi],vector3_[iv]);
	    }
	    if(!couplingSet) {
	      gs = abs(outgoingVertexF_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[F].size();++ix) {
	      (*ME[colourFlow[F][ix].first])(ifi, ifo, iv, ig) += 
		 colourFlow[F][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	  
	  // radiation from outgoing vector
	  if((decay[ivect]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[ivect]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexV_[inter]);
	    // ensure you get correct ougoing particle from first vertex
	    tcPDPtr off = decay[ivect]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction  vectorInter = 
	      outgoingVertexV_[inter]->evaluate(scale,3,off,gluon_[2*ig],
						vector3_[iv],decay[ivect]->mass());
	    
	    assert(vector3_[iv].particle()->id()==vectorInter.particle()->id());
	    if (ferm) {
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifi],wavebar3_[ifo],vectorInter);
	    }
	    else {
	      diag = 0.;
	      for(auto vertex : vertex_)
		diag += vertex->evaluate(scale,wave3_[ifo],wavebar3_[ifi],vectorInter);
	    }
	    if(!couplingSet) {
	      gs = abs(outgoingVertexV_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[V].size();++ix) {
	      (*ME[colourFlow[V][ix].first])(ifi, ifo, iv, ig) += 
		colourFlow[V][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	}
	if(massless) ++iv;
      }
    }
  }

  // contract matrices
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha(S,eM)
  output *= (4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}
#line 1 "./SFFDecayer.cc"
// -*- C++ -*-
//
// SFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SFFDecayer class.
//

#include "SFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SFFDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractFFSVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<FFSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
  }
}

void SFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_ 
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void SFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_ 
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SFFDecayer,GeneralTwoBodyDecayer>
describeHerwigSFFDecayer("Herwig::SFFDecayer", "Herwig.so");

void SFFDecayer::Init() {

  static ClassDocumentation<SFFDecayer> documentation
    ("This class implements to decay of a scalar to 2 fermions");

}

void SFFDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  int iferm(1),ianti(0);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) swap(iferm,ianti);
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),incoming,true);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double SFFDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin1Half)));
  // work out which is the fermion and antifermion
  int iferm(1),ianti(0);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(outgoing[ix]->CC()) itype[ix] = outgoing[ix]->id()>0 ? 0:1;
    else                   itype[ix] = 2;
  }
  if(itype[0]==0||itype[1]==1||(itype[0]==2&&itype[1]==2)) swap(iferm,ianti);

  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar_,momenta[iferm],outgoing[iferm],Helicity::outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave_   ,momenta[ianti],outgoing[ianti],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int ifm = 0; ifm < 2; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(iferm > ianti) (*ME())(0, ia, ifm) = 0.;
      else              (*ME())(0, ifm, ia) = 0.;
      for(auto vert : vertex_) {
	if(iferm > ianti){
	  (*ME())(0, ia, ifm) += vert->evaluate(scale,wave_[ia],
						wavebar_[ifm],swave_);
	}
	else {
	  (*ME())(0, ifm, ia) += vert->evaluate(scale,wave_[ia],
						wavebar_[ifm],swave_);
	}
      }
    }
  }

  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy SFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), outb.first, outa.first,
				     in);
    double mu1(outa.second/inpart.second),mu2(outb.second/inpart.second);
    double c2 = norm(perturbativeVertex_[0]->norm());
    Complex al(perturbativeVertex_[0]->left()), ar(perturbativeVertex_[0]->right());
    double me2 = -c2*( (norm(al) + norm(ar))*( sqr(mu1) + sqr(mu2) - 1.)
		       + 2.*(ar*conj(al) + al*conj(ar)).real()*mu1*mu2 );
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					outb.second);
    Energy output = me2*pcm/(8*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SFFDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the fermion and antifermion
  int ianti(0), iferm(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(iferm, ianti);
  if(itype[0]==2 && itype[1]==1) swap(iferm, ianti);
  if(itype[0]==0 && itype[1]==0 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  if(itype[0]==1 && itype[1]==1 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  
  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::
      calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar3_ ,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave3_    ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_    ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,     PDT::Spin1Half,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave3_   , decay[ianti],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(gluon_   , decay[iglu ],outgoing,true);

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  				          decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr outgoingVertexF;
  AbstractFFVVertexPtr outgoingVertexA;
  identifyVertices(iferm, ianti, inpart, decay, outgoingVertexF, outgoingVertexA,
		   inter);

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);

  Energy2 scale(sqr(inpart.mass()));
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ifm = 0; ifm < 2; ++ifm) {
    for(unsigned int ia = 0; ia < 2; ++ia) {
      for(unsigned int ig = 0; ig < 2; ++ig) {
	// radiation from the incoming scalar
	if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	  assert(incomingVertex_[inter]);

	  ScalarWaveFunction scalarInter = 
	    incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),
					      gluon_[2*ig],swave3_,inpart.mass());

	  assert(swave3_.particle()->id()==scalarInter.particle()->id());

	  if(!couplingSet) {
	    gs = abs(incomingVertex_[inter]->norm());
	    couplingSet = true;
	  }
	  Complex diag = 0.;
	  for(auto vertex : vertex_)
	    diag += vertex->evaluate(scale,wave3_[ia],
				     wavebar3_[ifm],scalarInter);
	  for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	    (*ME[colourFlow[0][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[0][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
	}
	// radiation from outgoing fermion
	if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED)) {
	  assert(outgoingVertexF);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[iferm]->dataPtr();
	  if(off->CC()) off = off->CC();
	  SpinorBarWaveFunction interS = 
	    outgoingVertexF->evaluate(scale,3,off,wavebar3_[ifm],
				      gluon_[2*ig],decay[iferm]->mass());
	  
	  assert(wavebar3_[ifm].particle()->id()==interS.particle()->id());
	  
	  if(!couplingSet) {
	    gs = abs(outgoingVertexF->norm());
	    couplingSet = true;
	  }
	  Complex diag = 0.;
	  for(auto vertex : vertex_)
	    diag += vertex->evaluate(scale,wave3_[ia], interS,swave3_);
	  for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	    (*ME[colourFlow[1][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[1][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
	}
	
	// radiation from outgoing antifermion
	if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	   (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED)) {
	  assert(outgoingVertexA);
	  // ensure you get correct outgoing particle from first vertex
	  tcPDPtr off = decay[ianti]->dataPtr();
	  if(off->CC()) off = off->CC();
	  SpinorWaveFunction  interS = 
	    outgoingVertexA->evaluate(scale,3,off,wave3_[ia],
				      gluon_[2*ig],decay[ianti]->mass());
	  
	  assert(wave3_[ia].particle()->id()==interS.particle()->id());

	  if(!couplingSet) {
	    gs = abs(outgoingVertexA->norm());
	    couplingSet = true;
	  }
	  Complex diag = 0.;
	  for(auto vertex : vertex_)
	    diag += vertex->evaluate(scale,interS,wavebar3_[ifm],swave3_);
	  for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	    (*ME[colourFlow[2][ix].first])(0, ia, ifm, ig) += 
	      colourFlow[2][ix].second*diag;
	  }
#ifdef GAUGE_CHECK
	  total+=norm(diag);
#endif
	}
      }
    }
  }
  
  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha(S,EM)
  output *= (4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
    // return the answer
  return output;
}

void SFFDecayer::identifyVertices(const int iferm, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractFFVVertexPtr & outgoingVertexF, 
				  AbstractFFVVertexPtr & outgoingVertexA,
				  ShowerInteraction inter) {
  // QCD
  if(inter==ShowerInteraction::QCD) {
    // work out which fermion each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[iferm]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[iferm]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour8))) {
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) {
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
      if     (outgoingVertex1_[inter]) outgoingVertexF = outgoingVertex1_[inter];
      else if(outgoingVertex2_[inter]) outgoingVertexF = outgoingVertex2_[inter];
      }
      else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour8) {
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
      }
      else if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3bar &&
	      decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) {
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[iferm]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexA = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
      else if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	      decay[ianti]->dataPtr()->iColour()==PDT::Colour3) {
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour6 ||
	    inpart.dataPtr()->iColour()==PDT::Colour6bar) {
      if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
  
    if (! ((incomingVertex_[inter]  && (outgoingVertexF  || outgoingVertexA)) ||
	   ( outgoingVertexF &&  outgoingVertexA))) {
      throw Exception()
	<< "Invalid vertices for QCD radiation in SFF decay in SFFDecayer::identifyVertices"
	<< Exception::runerror;
    }
  }
  // QED
  else {
    if(decay[iferm]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr())))
	outgoingVertexF = outgoingVertex1_[inter];
      else
	outgoingVertexF = outgoingVertex2_[inter];
    }
    if(decay[ianti]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr())))
	outgoingVertexA = outgoingVertex1_[inter];
      else
	outgoingVertexA = outgoingVertex2_[inter];
    }
  }
}
#line 1 "./SSSDecayer.cc"
// -*- C++ -*-
//
// SSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSSDecayer class.
//

#include "SSSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SSSDecayer::fullclone() const {
  return new_ptr(*this);
}

void SSSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractSSSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<SSSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
  }
}

void SSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void SSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSSDecayer,GeneralTwoBodyDecayer>
describeHerwigSSSDecayer("Herwig::SSSDecayer", "Herwig.so");

void SSSDecayer::Init() {

  static ClassDocumentation<SSSDecayer> documentation
    ("This class implements the decay of a scalar to 2 scalars.");

}

void SSSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::
      constructSpinInfo(decay[ix],outgoing,true);
}

double SSSDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  ScalarWaveFunction s1(momenta[0],outgoing[0],Helicity::outgoing);
  ScalarWaveFunction s2(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  (*ME())(0,0,0) = 0.;
  for(auto vert : vertex_) {
    (*ME())(0,0,0) += vert->evaluate(scale,s1,s2,swave_);
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy SSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0] && !perturbativeVertex_[0]->kinematics()) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, in, outa.first, outb.first);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second, outa.second,
					       outb.second);
    double c2 = norm(perturbativeVertex_[0]->norm());
    Energy pWidth = c2*pcm/8./Constants::pi/scale*UnitRemoval::E2;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SSSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the scalar and anti scalar
  int ianti(0), iscal(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(ianti, iscal);
  if(itype[0]==2 && itype[1]==1) swap(ianti, iscal);
  if(itype[0]==0 && itype[1]==0 && abs(decay[0]->dataPtr()->id())>abs(decay[1]->dataPtr()->id())) 
    swap(iscal, ianti);
  if(itype[0]==1 && itype[1]==1 && abs(decay[0]->dataPtr()->id())<abs(decay[1]->dataPtr()->id())) 
    swap(iscal, ianti);

  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    ScalarWaveFunction::
      constructSpinInfo(decay[iscal],outgoing,true);
    ScalarWaveFunction::
      constructSpinInfo(decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_,decay[iglu ],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin0,
								       PDT::Spin0, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  ScalarWaveFunction anti(decay[ianti]->momentum(), decay[ianti]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_,decay[iglu ],outgoing,true);

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
					  decay[iglu ]->dataPtr(),10,
					  outgoing));
    }
  }
#endif

  AbstractVSSVertexPtr outgoingVertexS;
  AbstractVSSVertexPtr outgoingVertexA;
  identifyVertices(iscal, ianti, inpart, decay, outgoingVertexS, outgoingVertexA,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int ig = 0; ig < 2; ++ig) {
    // radiation from the incoming scalar
    if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
       (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
      assert(incomingVertex_[inter]);
      ScalarWaveFunction scalarInter = 
	incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),
					  gluon_[2*ig],swave3_,inpart.mass());

      assert(swave3_.particle()->id()==scalarInter.particle()->id());

      Complex diag = 0.;
      for(auto vertex : vertex_)
	diag += vertex->evaluate(scale,scal,anti,scalarInter);
      if(!couplingSet) {
	gs = abs(incomingVertex_[inter]->norm());
	couplingSet = true;
      }
      for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	(*ME[colourFlow[0][ix].first])(0, 0, 0, ig) += 
	   colourFlow[0][ix].second*diag; 
      }
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
    }
    // radiation from the outgoing scalar
    if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
       (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
      assert(outgoingVertexS);
      // ensure you get correct outgoing particle from first vertex
      tcPDPtr off = decay[iscal]->dataPtr();
      if(off->CC()) off = off->CC();
      ScalarWaveFunction scalarInter = 
	outgoingVertexS->evaluate(scale,3,off,gluon_[2*ig],scal,decay[iscal]->mass());

      assert(scal.particle()->id()==scalarInter.particle()->id());

      Complex diag = 0.;
      for(auto vertex : vertex_)
	diag += vertex->evaluate(scale,swave3_,anti,scalarInter);
      if(!couplingSet) {
	gs = abs(outgoingVertexS->norm());
	couplingSet = true;
      }
      for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	(*ME[colourFlow[1][ix].first])(0, 0, 0, ig) += 
	   colourFlow[1][ix].second*diag;
      }
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
    }
    // radiation from the outgoing anti scalar
    if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
       (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
      assert(outgoingVertexA);
      // ensure you get correct outgoing particle from first vertex
      tcPDPtr off = decay[ianti]->dataPtr();
      if(off->CC()) off = off->CC();
      ScalarWaveFunction scalarInter = 
	outgoingVertexA->evaluate(scale,3,off, gluon_[2*ig],anti,decay[ianti]->mass());

      assert(anti.particle()->id()==scalarInter.particle()->id());

      Complex diag = 0.;
      for(auto vertex : vertex_)
	diag += vertex->evaluate(scale,swave3_,scal,scalarInter);
      if(!couplingSet) {
	gs = abs(outgoingVertexA->norm());
	couplingSet = true;
      }
      for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	(*ME[colourFlow[2][ix].first])(0, 0, 0, ig) += 
	  colourFlow[2][ix].second*diag;
      }
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
    }
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}


void SSSDecayer::identifyVertices(const int iscal, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractVSSVertexPtr & outgoingVertexS, 
				  AbstractVSSVertexPtr & outgoingVertexA,
				  ShowerInteraction inter){
  // QCD
  if(inter==ShowerInteraction::QCD) {
    // work out which scalar each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[iscal]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[iscal]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    // one outgoing vertex
    else if(inpart.dataPtr()    ->iColour()==PDT::Colour3){ 
      if(decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexS = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexS = outgoingVertex2_[inter]; 
      }
      else if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&  
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[ianti]->dataPtr()->id()))){
	  outgoingVertexS = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertexS = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
      }
    }
    else if(inpart.dataPtr()    ->iColour()==PDT::Colour3bar){
      if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[iscal]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexA = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	       decay[iscal]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->dataPtr()->id()))){
	  outgoingVertexS = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexS = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertexS || outgoingVertexA)) ||
	   ( outgoingVertexS &&  outgoingVertexA)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in SSS decay in SSSDecayer::identifyVertices"
	<< Exception::runerror;
  }
  // QED
  else {
    if(decay[iscal]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iscal]->dataPtr())))
	outgoingVertexS = outgoingVertex1_[inter];
      else
	outgoingVertexS = outgoingVertex2_[inter];
    }
    if(decay[ianti]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr())))
	outgoingVertexA = outgoingVertex1_[inter];
      else
	outgoingVertexA = outgoingVertex2_[inter];
    }
  }
}
#line 1 "./SSVDecayer.cc"
// -*- C++ -*-
//
// SSVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SSVDecayer class.
//

#include "SSVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SSVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SSVDecayer::fullclone() const {
  return new_ptr(*this);
}

void SSVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex)
    vertex_            .push_back(dynamic_ptr_cast<AbstractVSSVertexPtr>(vert));
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter]  = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractVVSSVertexPtr>(fourV.at(inter));
    outgoingVertexS_[inter] = AbstractVSSVertexPtr();
    outgoingVertexV_[inter] = AbstractVVVVertexPtr();
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
      else
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    }
  }
}

void SSVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_
     << incomingVertex_   << outgoingVertexS_
     << outgoingVertexV_  << fourPointVertex_;
}

void SSVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_
     >> incomingVertex_   >> outgoingVertexS_
     >> outgoingVertexV_  >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SSVDecayer,GeneralTwoBodyDecayer>
describeHerwigSSVDecayer("Herwig::SSVDecayer", "Herwig.so");

void SSVDecayer::Init() {

  static ClassDocumentation<SSVDecayer> documentation
    ("This implements the decay of a scalar to a vector and a scalar");

}

void SSVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int isc(0),ivec(1);
  if(decay[0]->dataPtr()->iSpin() != PDT::Spin0) swap(isc,ivec);
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),incoming,true);
  ScalarWaveFunction::
    constructSpinInfo(decay[isc],outgoing,true);
  VectorWaveFunction::
    constructSpinInfo(vector_,decay[ivec],outgoing,true,false);
}

double SSVDecayer::me2(const int,const Particle & part,
		       const tPDVector & decay,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  unsigned int isc(0),ivec(1);
  if(decay[0]->iSpin() != PDT::Spin0) swap(isc,ivec);
  if(!ME()) {
    if(ivec==1)
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin1)));
    else
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin0)));
  }
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  bool massless = decay[1]->id()==ParticleID::gamma || decay[1]->id()==ParticleID::g;
  VectorWaveFunction::
    calculateWaveFunctions(vector_,momenta[ivec],decay[ivec],Helicity::outgoing,massless);
  ScalarWaveFunction sca(momenta[isc],decay[isc],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  //make sure decay matrix element is in the correct order
  double output(0.);
  if(ivec == 0) {
    for(unsigned int ix = 0; ix < 3; ++ix) {
      (*ME())(0, ix, 0) = 0.;
      for(auto vert : vertex_)
	(*ME())(0, ix, 0) += vert->evaluate(scale,vector_[ix],sca, swave_);
    }
  }
  else {
    for(unsigned int ix = 0; ix < 3; ++ix) {
      (*ME())(0, 0, ix) = 0.;
      for(auto vert : vertex_)
	(*ME())(0, 0, ix) += vert->evaluate(scale,vector_[ix],sca,swave_);
    }
  }
  output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),decay[0],decay[1]);
  // return the answer
  return output;
}

Energy SSVDecayer:: partialWidth(PMPair inpart, PMPair outa,
				 PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
}


double  SSVDecayer::threeBodyME(const int , const Particle & inpart,
				const ParticleVector & decay,
				ShowerInteraction inter, MEOption meopt) {
  int iscal (0), ivect (1), iglu (2);
  // get location of outgoing scalar/vector
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) swap(iscal,ivect);

  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    ScalarWaveFunction::
      constructSpinInfo(decay[iscal],outgoing,true);
     VectorWaveFunction::
      constructSpinInfo(vector3_,decay[ivect],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(gluon_,  decay[iglu ],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin0,
								       PDT::Spin1, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal(decay[iscal]->momentum(),  decay[iscal]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(vector3_,decay[ivect],outgoing,false);
  VectorWaveFunction::calculateWaveFunctions(gluon_,  decay[iglu ],outgoing,true );

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  					  decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  if (! ((incomingVertex_[inter]  && (outgoingVertexS_[inter] || outgoingVertexV_[inter])) ||
	 (outgoingVertexS_[inter] &&  outgoingVertexV_[inter])))
    throw Exception()
      << "Invalid vertices for radiation in SSV decay in SSVDecayer::threeBodyME"
      << Exception::runerror;


  // sort out colour flows
  int S(1), V(2);
  if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3bar &&
      decay[ivect]->dataPtr()->iColour()==PDT::Colour8)
    swap(S,V);
  else if (decay[ivect]->dataPtr()->iColour()==PDT::Colour3 &&
	   decay[iscal]->dataPtr()->iColour()==PDT::Colour8)
    swap(S,V);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ig = 0; ig < 2; ++ig) {
      // radiation from the incoming scalar
      if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(incomingVertex_[inter]);
	ScalarWaveFunction scalarInter =
	  incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),
					   gluon_[2*ig],swave3_,inpart.mass());

	assert(swave3_.particle()->id()==scalarInter.particle()->id());
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vector3_[iv],scal,scalarInter);
	if(!couplingSet) {
	  gs = abs(incomingVertex_[inter]->norm());
	  couplingSet = true;
	}
	for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	  (*ME[colourFlow[0][ix].first])(0, 0, iv, ig) +=
	    colourFlow[0][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }
      // radiation from the outgoing scalar
      if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(outgoingVertexS_[inter]);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[iscal]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter =
	  outgoingVertexS_[inter]->evaluate(scale,3,off,gluon_[2*ig],scal,decay[iscal]->mass());

	assert(scal.particle()->id()==scalarInter.particle()->id());

	if(!couplingSet) {
	  gs = abs(outgoingVertexS_[inter]->norm());
	  couplingSet = true;
	}
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vector3_[iv],scalarInter,swave3_);
	for(unsigned int ix=0;ix<colourFlow[S].size();++ix) {
	  (*ME[colourFlow[S][ix].first])(0, 0, iv, ig) +=
	    colourFlow[S][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }

      // radiation from outgoing vector
      if((decay[ivect]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (decay[ivect]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(outgoingVertexV_[inter]);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[ivect]->dataPtr();
	if(off->CC()) off = off->CC();
	VectorWaveFunction  vectorInter =
	  outgoingVertexV_[inter]->evaluate(scale,3,off,gluon_[2*ig],
					    vector3_[iv],decay[ivect]->mass());

	assert(vector3_[iv].particle()->id()==vectorInter.particle()->id());

	if(!couplingSet) {
	  gs = abs(outgoingVertexV_[inter]->norm());
	  couplingSet = true;
	}
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vectorInter,scal,swave3_);
	for(unsigned int ix=0;ix<colourFlow[V].size();++ix) {
	  (*ME[colourFlow[V][ix].first])(0, 0, iv, ig) +=
	    colourFlow[V][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }
      // radiation from 4 point vertex
      if (fourPointVertex_[inter]) {
	Complex diag =  fourPointVertex_[inter]->evaluate(scale, gluon_[2*ig], vector3_[iv],
							  scal, swave3_);
	for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	  (*ME[colourFlow[3][ix].first])(0, 0, iv, ig) +=
	     colourFlow[3][ix].second*diag;
	}
#ifdef GAUGE_CHECK
	total+=norm(diag);
#endif
      }
    }
  }

  // contract matrices
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}
#line 1 "./SVVDecayer.cc"
// -*- C++ -*-
//
// SVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SVVDecayer class.
//

#include "SVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/Vertex/Scalar/VVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void SVVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVVSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VVSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
  }
}

void SVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void SVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SVVDecayer,GeneralTwoBodyDecayer>
describeHerwigSVVDecayer("Herwig::SVVDecayer", "Herwig.so");

void SVVDecayer::Init() {

  static ClassDocumentation<SVVDecayer> documentation
    ("This implements the decay of a scalar to 2 vector bosons.");

}


void SVVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),incoming,true);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      constructSpinInfo(vectors_[ix],decay[ix],outgoing,true,photon[ix]);
}

double SVVDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = outgoing[ix]->mass()==ZERO;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix],momenta[ix],outgoing[ix],Helicity::outgoing,photon[ix]);
  
  
  Energy2 scale(sqr(part.mass()));
  unsigned int iv1,iv2;
  for(iv2 = 0; iv2 < 3; ++iv2) {
    if( photon[1] && iv2 == 1 ) ++iv2;
    for(iv1=0;iv1<3;++iv1) {
      if( photon[0] && iv1 == 1) ++iv1;
      (*ME())(0, iv1, iv2) = 0.;
      for(auto vert : vertex_)
	(*ME())(0, iv1, iv2) += vert->evaluate(scale,vectors_[0][iv1],
					       vectors_[1][iv2],swave_);
    }
  }
  double output = ME()->contract(rho_).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy SVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, outa.first , 
				    outb.first, in);
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    double m1pm2 = mu1sq + mu2sq;
    double me2(0.); 
    if( mu1sq > 0. && mu2sq > 0.)
      me2 = ( m1pm2*(m1pm2 - 2.) + 8.*mu1sq*mu2sq + 1.)/4./mu1sq/mu2sq;
    else if( mu1sq == 0. || mu2sq == 0. )
      me2 = 3.;
    else 
      me2 = 4.;
    
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*
      me2*pcm/(8*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double SVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  if(meopt==Initialize) {
    // create scalar wavefunction for decaying particle
    ScalarWaveFunction::
      calculateWaveFunctions(rho3_,const_ptr_cast<tPPtr>(&inpart),incoming);
    swave3_ = ScalarWaveFunction(inpart.momentum(),inpart.dataPtr(),incoming);
  }
  if(meopt==Terminate) {
    ScalarWaveFunction::
      constructSpinInfo(const_ptr_cast<tPPtr>(&inpart),incoming,true);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[0],decay[0],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[1],decay[1],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(gluon_      ,decay[2],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin0, PDT::Spin1,
								       PDT::Spin1, PDT::Spin1)));
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix)
    massless[ix] = decay[ix]->mass()!=ZERO;
  // create wavefunctions
  VectorWaveFunction::calculateWaveFunctions(vectors3_[0],decay[0],outgoing,massless[0]);
  VectorWaveFunction::calculateWaveFunctions(vectors3_[1],decay[1],outgoing,massless[1]);
  VectorWaveFunction::calculateWaveFunctions(gluon_      ,decay[2],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[2]->momentum(),
  					  decay[2]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  // get the outgoing vertices
  AbstractVVVVertexPtr outgoingVertex1;
  AbstractVVVVertexPtr outgoingVertex2;
  identifyVertices(inpart,decay, outgoingVertex1, outgoingVertex2,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif

   for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
     if(massless[0] && iv1==1) continue;
     for(unsigned int iv2 = 0; iv2 < 3; ++iv2) {
       if(massless[1] && iv2==1) continue;
       for(unsigned int ig = 0; ig < 2; ++ig) {
	 // radiation from the incoming vector
	 if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(incomingVertex_[inter]);	
	   ScalarWaveFunction scalarInter = 
	     incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),gluon_[2*ig],
					      swave3_,inpart.mass());
	
	   assert(swave3_.particle()->id()==scalarInter.particle()->id());
	
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv1],
				      vectors3_[1][iv2],scalarInter);
	   if(!couplingSet) {
	     gs = abs(incomingVertex_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	     (*ME[colourFlow[0][ix].first])(0, iv1, iv2, ig) += 
	       colourFlow[0][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 // radiation from the 1st outgoing vector
	 if((decay[0]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[0]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertex1);
	// ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[0]->dataPtr();
	   if(off->CC()) off = off->CC();
	   VectorWaveFunction vectorInter = 
	     outgoingVertex1->evaluate(scale,3,off,gluon_[2*ig],vectors3_[0][iv1],decay[0]->mass());
	   
	   assert(vectors3_[0][iv1].particle()->id()==vectorInter.particle()->id());
	 
	   Complex diag =0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectorInter,vectors3_[1][iv2],swave3_);
	   if(!couplingSet) {
	     gs = abs(outgoingVertex1->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	     (*ME[colourFlow[1][ix].first])(0, iv1, iv2, ig) += 
	       colourFlow[1][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 
	 if((decay[1]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[1]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertex2);
	   // ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[1]->dataPtr();
	   if(off->CC()) off = off->CC();
	   VectorWaveFunction vectorInter = 
	     outgoingVertex2->evaluate(scale,3,off, gluon_[2*ig],vectors3_[1][iv2],decay[1]->mass());
	
	   assert(vectors3_[1][iv2].particle()->id()==vectorInter.particle()->id());
	
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv1],vectorInter,swave3_);
	   if(!couplingSet) {
	     gs = abs(outgoingVertex2->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	     (*ME[colourFlow[2][ix].first])(0, iv1, iv2, ig) += 
	       colourFlow[2][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
       }
     }
   }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}

void SVVDecayer::identifyVertices(const Particle & inpart, const ParticleVector & decay, 
				  AbstractVVVVertexPtr & outgoingVertex1, 
				  AbstractVVVVertexPtr & outgoingVertex2,
				  ShowerInteraction inter) {
  if(inter==ShowerInteraction::QCD) {
    // work out which scalar each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[0]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[1]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[0]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[1]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex2_[inter];
	outgoingVertex2 = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[0]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[1]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex2_[inter];
	outgoingVertex2 = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[0]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[1]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertex1 = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertex1 = outgoingVertex2_[inter];
      }
      else if (decay[0]->dataPtr()->iColour()==PDT::Colour3 &&
	       decay[1]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[1]->dataPtr()->id()))){
	  outgoingVertex1 = outgoingVertex2_[inter];
	  outgoingVertex2 = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertex1 = outgoingVertex1_[inter];
	  outgoingVertex2 = outgoingVertex2_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[1]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[0]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertex2 = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (decay[0]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[1]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->dataPtr()->id()))){
	  outgoingVertex1 = outgoingVertex1_[inter];
	  outgoingVertex2 = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertex1 = outgoingVertex2_[inter];
	  outgoingVertex2 = outgoingVertex1_[inter];
	}
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertex1  || outgoingVertex2)) ||
	   ( outgoingVertex1 &&  outgoingVertex2)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in SVV decay in SVVDecayer::identifyVertices"
	<< Exception::runerror;
  }
  else {
    if(decay[0]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[0]->dataPtr())))
	outgoingVertex1 = outgoingVertex1_[inter];
      else
	outgoingVertex1 = outgoingVertex2_[inter];
    }
    if(decay[1]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[1]->dataPtr())))
	outgoingVertex2 = outgoingVertex1_[inter];
      else
	outgoingVertex2 = outgoingVertex2_[inter];
    }
  }
}
#line 1 "./TFFDecayer.cc"
// -*- C++ -*-
//
// TFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TFFDecayer class.
//

#include "TFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
IBPtr TFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void TFFDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractFFTVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<FFTVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractFFVTVertexPtr>(fourV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr> (outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr> (outV[1].at(inter));
  }
}

void TFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_          << perturbativeVertex_
     << outgoingVertex1_ << outgoingVertex2_
     << fourPointVertex_;
}

void TFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_          >> perturbativeVertex_
     >> outgoingVertex1_ >> outgoingVertex2_
     >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TFFDecayer,GeneralTwoBodyDecayer>
describeHerwigTFFDecayer("Herwig::TFFDecayer", "Herwig.so");

void TFFDecayer::Init() {

  static ClassDocumentation<TFFDecayer> documentation
    ("The TFFDecayer class implements the decay of a tensor particle "
     "to 2 fermions ");
  
}

void TFFDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int iferm(0),ianti(1);
  if(decay[0]->id()>=0) swap(iferm,ianti);
  TensorWaveFunction::
    constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
		      incoming,true,false);
  SpinorBarWaveFunction::
    constructSpinInfo(wavebar_,decay[iferm],outgoing,true);
  SpinorWaveFunction::
    constructSpinInfo(wave_   ,decay[ianti],outgoing,true);
}

double TFFDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  unsigned int iferm(0),ianti(1);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin1Half,PDT::Spin1Half)));
  if(outgoing[0]->id()>=0) swap(iferm,ianti);
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar_,momenta[iferm],outgoing[iferm],Helicity::outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave_   ,momenta[ianti],outgoing[ianti],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  unsigned int thel,fhel,ahel;
  for(thel=0;thel<5;++thel) {
    for(fhel=0;fhel<2;++fhel) {
      for(ahel=0;ahel<2;++ahel) {
	if(iferm > ianti) {
	  (*ME())(thel,fhel,ahel) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(thel,fhel,ahel) += 
	      vert->evaluate(scale,wave_[ahel],
			     wavebar_[fhel],tensors_[thel]);
	}
	else {
	  (*ME())(thel,ahel,fhel) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(thel,ahel,fhel) += 
	      vert->evaluate(scale,wave_[ahel],
			     wavebar_[fhel],tensors_[thel]);
	}
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy TFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale = sqr(inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, in, outa.first, outb.first);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1- 4.*musq);
    double me2 = b*b*(5-2*b*b)*scale/120.*UnitRemoval::InvE2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm/(8.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double TFFDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the fermion and antifermion
  int ianti(0), iferm(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(iferm, ianti);
  if(itype[0]==2 && itype[1]==1) swap(iferm, ianti);
  if(itype[0]==0 && itype[1]==0 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  if(itype[0]==1 && itype[1]==1 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);

  if(meopt==Initialize) {
    // create tensor wavefunction for decaying particle
    TensorWaveFunction::
      calculateWaveFunctions(tensors3_, rho3_, const_ptr_cast<tPPtr>(&inpart), incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(tensors3_, const_ptr_cast<tPPtr>(&inpart),incoming,true, false);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar3_ ,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave3_    ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_    ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin2,     PDT::Spin1Half,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave3_   , decay[ianti],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(gluon_   , decay[iglu ],outgoing,true);

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  				          decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif
  
  if (! (outgoingVertex1_[inter] && outgoingVertex2_[inter]))
    throw Exception()
      << "Invalid vertices for QCD radiation in TFF decay in TFFDecayer::threeBodyME"
      << Exception::runerror;

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr outgoingVertexF = outgoingVertex1_[inter];
  AbstractFFVVertexPtr outgoingVertexA = outgoingVertex2_[inter];

  if(outgoingVertex1_[inter]!=outgoingVertex2_[inter] &&
     outgoingVertex1_[inter]->isIncoming(getParticleData(decay[ianti]->id())))
    swap (outgoingVertexF, outgoingVertexA);  
  
  if(! (inpart.dataPtr()->iColour()==PDT::Colour0)){
    throw Exception()
      << "Invalid vertices for QCD radiation in TFF decay in TFFDecayer::threeBodyME"
      << Exception::runerror;
  }

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int it = 0; it < 5; ++it) {  
    for(unsigned int ifm = 0; ifm < 2; ++ifm) {
      for(unsigned int ia = 0; ia < 2; ++ia) {
	for(unsigned int ig = 0; ig < 2; ++ig) {

	  // radiation from outgoing fermion
	  if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexF);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[iferm]->dataPtr();
	    if(off->CC()) off = off->CC();
	    SpinorBarWaveFunction interS = 
	      outgoingVertexF->evaluate(scale,3,off,wavebar3_[ifm],
					gluon_[2*ig],decay[iferm]->mass());
	  
	    assert(wavebar3_[ifm].particle()->id()==interS.particle()->id());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,wave3_[ia], interS,tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertexF->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(it, ifm, ia, ig) += 
		colourFlow[1][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from outgoing antifermion
	  if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexA);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[ianti]->dataPtr();
	    if(off->CC()) off = off->CC();
	    SpinorWaveFunction  interS = 
	      outgoingVertexA->evaluate(scale,3,off,wave3_[ia],
					gluon_[2*ig],decay[ianti]->mass());
	    
	    assert(wave3_[ia].particle()->id()==interS.particle()->id());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,interS,wavebar3_[ifm],tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertexA->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(it, ifm, ia, ig) += 
		colourFlow[2][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from 4 point vertex
	  if (fourPointVertex_[inter]) {
	    Complex diag = fourPointVertex_[inter]->evaluate(scale, wave3_[ia], wavebar3_[ifm],
							     gluon_[2*ig], tensors3_[it]);
	    for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	      (*ME[colourFlow[3][ix].first])(it, ifm, ia, ig) += 
		colourFlow[3][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	}
      }
    }
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(s,em)
  output *= (4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}

#line 1 "./TSSDecayer.cc"
// -*- C++ -*-
//
// TSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TSSDecayer class.
//

#include "TSSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TSSDecayer::fullclone() const {
  return new_ptr(*this);
}


void TSSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & ,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & ,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractSSTVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<SSTVertexPtr>        (vert));
  }
}


void TSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_;
}

void TSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TSSDecayer,GeneralTwoBodyDecayer>
describeHerwigTSSDecayer("Herwig::TSSDecayer", "Herwig.so");

void TSSDecayer::Init() {

  static ClassDocumentation<TSSDecayer> documentation
    ("This class implements the decay of a tensor particle into "
     "2 scalars.");

}

void TSSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  TensorWaveFunction::
    constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
		      incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::
      constructSpinInfo(decay[ix],outgoing,true);
}

double TSSDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
			     incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  ScalarWaveFunction sca1(momenta[0],outgoing[0],Helicity::outgoing);
  ScalarWaveFunction sca2(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int thel=0;thel<5;++thel) {
    (*ME())(thel,0,0) =0.;
    for(auto vert : vertex_)
      (*ME())(thel,0,0) += vert->evaluate(scale,sca1,sca2,tensors_[thel]); 
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}


Energy TSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, outa.first, outb.first, in);
    double musq = sqr(outa.second/inpart.second);
    double b = sqrt(1. - 4.*musq);
    double me2 = scale*pow(b,4)/120*UnitRemoval::InvE2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm/(8.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

#line 1 "./TVVDecayer.cc"
// -*- C++ -*-
//
// TVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TVVDecayer class.
//

#include "TVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/LorentzTensor.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr TVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr TVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void TVVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVVTVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VVTVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractVVVTVertexPtr>(fourV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr> (outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr> (outV[1].at(inter));
  }
}

void TVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_          << perturbativeVertex_
     << outgoingVertex1_ << outgoingVertex2_
     << fourPointVertex_;
}

void TVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_          >> perturbativeVertex_
     >> outgoingVertex1_ >> outgoingVertex2_
     >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<TVVDecayer,GeneralTwoBodyDecayer>
describeHerwigTVVDecayer("Herwig::TVVDecayer", "Herwig.so");

void TVVDecayer::Init() {

  static ClassDocumentation<TVVDecayer> documentation
    ("This class implements the decay of a tensor to 2 vector bosons");

}

void TVVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = decay[ix]->mass()==ZERO;
  TensorWaveFunction::
    constructSpinInfo(tensors_,const_ptr_cast<tPPtr>(&part),
		      incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      constructSpinInfo(vectors_[ix],decay[ix],outgoing,true,photon[ix]);
}

double TVVDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin2,PDT::Spin1,PDT::Spin1)));
  bool photon[2];
  for(unsigned int ix=0;ix<2;++ix)
    photon[ix] = outgoing[ix]->mass()==ZERO;
  if(meopt==Initialize) {
    TensorWaveFunction::
      calculateWaveFunctions(tensors_,rho_,const_ptr_cast<tPPtr>(&part),
  			     incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix],momenta[ix],outgoing[ix],Helicity::outgoing,photon[ix]);
  Energy2 scale(sqr(part.mass()));
  unsigned int thel,v1hel,v2hel;
  for(thel=0;thel<5;++thel) {
    for(v1hel=0;v1hel<3;++v1hel) {
      for(v2hel=0;v2hel<3;++v2hel) {
  	(*ME())(thel,v1hel,v2hel) = 0.;
  	  for(auto vert : vertex_)
  	    (*ME())(thel,v1hel,v2hel) += vert->evaluate(scale,
  							vectors_[0][v1hel],
  							vectors_[1][v2hel],
  							tensors_[thel]);
  	if(photon[1]) ++v2hel;
      }
      if(photon[0]) ++v1hel;
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}
  
Energy TVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(scale, outa.first, outb.first, in);
    double mu2 = sqr(outa.second/inpart.second);
    double b = sqrt(1 - 4.*mu2);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy2 me2;
    if(outa.second > ZERO && outb.second > ZERO)
      me2 = scale*(30 - 20.*b*b + 3.*pow(b,4))/120.; 
    else 
      me2 = scale/10.;
    
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm
      /(8.*Constants::pi)*UnitRemoval::InvE2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double TVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix)
    massless[ix] = decay[ix]->mass()==ZERO;
  int iglu(2);  
  if(meopt==Initialize) {
    // create tensor wavefunction for decaying particle
    TensorWaveFunction::
      calculateWaveFunctions(tensors3_, rho3_, const_ptr_cast<tPPtr>(&inpart), incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    TensorWaveFunction::
      constructSpinInfo(tensors3_, const_ptr_cast<tPPtr>(&inpart),incoming,true, false);
    for(unsigned int ix=0;ix<2;++ix)
      VectorWaveFunction::
	constructSpinInfo(vectors3_[ix],decay[ix   ],outgoing,true, massless[ix]);
    VectorWaveFunction::
        constructSpinInfo(gluon_       ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin2, PDT::Spin1,
								       PDT::Spin1, PDT::Spin1)));
  // create wavefunctions
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors3_[ix],decay[ix   ],outgoing,massless[ix]);
  VectorWaveFunction::
      calculateWaveFunctions(gluon_       ,decay[iglu ],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  				          decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif
  
  // work out which vector each outgoing vertex corresponds to 
  if(outgoingVertex1_[inter]!=outgoingVertex2_[inter] &&
     outgoingVertex1_[inter]->isIncoming(getParticleData(decay[1]->id())))
    swap(outgoingVertex1_[inter], outgoingVertex2_[inter]);
  
  if (! (outgoingVertex1_[inter] && outgoingVertex2_[inter]))
    throw Exception()
      << "Invalid vertices for radiation in TVV decay in TVVDecayer::threeBodyME"
      << Exception::runerror;

  if( !(!inpart.dataPtr()->coloured() && inter ==ShowerInteraction::QCD) &&
      !(!inpart.dataPtr()->charged()  && inter ==ShowerInteraction::QED))
    throw Exception()
      << "Invalid vertices for radiation in TVV decay in TVVDecayer::threeBodyME"
      << Exception::runerror;


  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int it = 0; it < 5; ++it) {  
    for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
      for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
	for(unsigned int ig = 0; ig < 2; ++ig) {

	  // radiation from first outgoing vector
	  if((decay[0]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[0]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertex1_[inter]);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[0]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction vectInter = 
	      outgoingVertex1_[inter]->evaluate(scale,3,off,gluon_[2*ig],
						 vectors3_[0][iv0],decay[0]->mass());
	  
	    assert(vectors3_[0][iv0].particle()->PDGName()==vectInter.particle()->PDGName());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,vectors3_[1][iv1], 
				       vectInter,tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertex1_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[1][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from second outgoing vector
	  if((decay[1]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[1]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertex2_[inter]);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[1]->dataPtr();
	    if(off->CC()) off = off->CC();
	    VectorWaveFunction  vectInter = 
	      outgoingVertex2_[inter]->evaluate(scale,3,off,vectors3_[1][iv1],
						gluon_[2*ig],decay[1]->mass());
	    
	    assert(vectors3_[1][iv1].particle()->PDGName()==vectInter.particle()->PDGName());
	    
	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,vectInter,vectors3_[0][iv0],
					tensors3_[it]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertex2_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[2][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from 4 point vertex
	  if (fourPointVertex_[inter]) {
	    Complex diag = fourPointVertex_[inter]->evaluate(scale, vectors3_[0][iv0],
							     vectors3_[1][iv1],gluon_[2*ig], 
							     tensors3_[it]);
	    for(unsigned int ix=0;ix<colourFlow[3].size();++ix) {
	      (*ME[colourFlow[3][ix].first])(it, iv0, iv1, ig) += 
		colourFlow[3][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	}
	if(massless[1]) ++iv1;
      }
      if(massless[0]) ++iv0;
    }
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(s,em)
  output *= (4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}
#line 1 "./VFFDecayer.cc"
// -*- C++ -*-
//
// VFFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VFFDecayer class.
//

#include "VFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/Vertex/Vector/FFVVertex.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
using namespace Herwig;
using namespace ThePEG::Helicity;
IBPtr VFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void VFFDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractFFVVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<FFVVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractFFVVertexPtr>(outV[1].at(inter));
  }
}

void VFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void VFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VFFDecayer,GeneralTwoBodyDecayer>
describeHerwigVFFDecayer("Herwig::VFFDecayer", "Herwig.so");

void VFFDecayer::Init() {

  static ClassDocumentation<VFFDecayer> documentation
    ("The VFFDecayer implements the matrix element for the"
     " decay of a vector to fermion-antifermion pair");

}

void VFFDecayer::
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

double VFFDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  int iferm(1),ianti(0);
  if(outgoing[0]->id()>0) swap(iferm,ianti);
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1Half,PDT::Spin1Half)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
  					       const_ptr_cast<tPPtr>(&part),
  					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar_,momenta[iferm],outgoing[iferm],Helicity::outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave_   ,momenta[ianti],outgoing[ianti],Helicity::outgoing);
  // compute the matrix element
  Energy2 scale(part.mass()*part.mass());
  for(unsigned int ifm = 0; ifm < 2; ++ifm) { //loop over fermion helicities
    for(unsigned int ia = 0; ia < 2; ++ia) {// loop over antifermion helicities
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {//loop over vector helicities
  	if(iferm > ianti) {
  	  (*ME())(vhel, ia, ifm) = 0.;
  	  for(auto vert : vertex_)
  	    (*ME())(vhel, ia, ifm) += 
  	      vert->evaluate(scale,wave_[ia],
  			     wavebar_[ifm],vectors_[vhel]);
  	}
  	else {
  	  (*ME())(vhel,ifm,ia)= 0.;
  	  for(auto vert : vertex_)
  	    (*ME())(vhel,ifm,ia) +=
  	      vert->evaluate(scale,wave_[ia],
  			     wavebar_[ifm],vectors_[vhel]);
  	}
      }
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy VFFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    double mu1(outa.second/inpart.second), mu2(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), outa.first, outb.first,in);
    Complex cl(perturbativeVertex_[0]->left()), cr(perturbativeVertex_[0]->right());
    double me2 = (norm(cl) + norm(cr))*( sqr(sqr(mu1) - sqr(mu2)) 
					 + sqr(mu1) + sqr(mu2) - 2.)
      - 6.*(cl*conj(cr) + cr*conj(cl)).real()*mu1*mu2;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = -norm(perturbativeVertex_[0]->norm())*me2*pcm / 
      (24.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}


double VFFDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the fermion and antifermion
  int ianti(0), iferm(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(iferm, ianti);
  if(itype[0]==2 && itype[1]==1) swap(iferm, ianti);
  if(itype[0]==0 && itype[1]==0 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);
  if(itype[0]==1 && itype[1]==1 && decay[0]->dataPtr()->id()<decay[1]->dataPtr()->id()) 
    swap(iferm, ianti);

  if(meopt==Initialize) {
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(vector3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
					       incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(vector3_ ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,false);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar3_,decay[iferm],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave3_   ,decay[ianti],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_   ,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);

  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1,     PDT::Spin1Half,
								       PDT::Spin1Half, PDT::Spin1)));
  // create wavefunctions
  SpinorBarWaveFunction::
    calculateWaveFunctions(wavebar3_, decay[iferm],outgoing);
  SpinorWaveFunction::
    calculateWaveFunctions(wave3_   , decay[ianti],outgoing);
  VectorWaveFunction::
    calculateWaveFunctions(gluon_   , decay[iglu ],outgoing,true);

  // gauge invariance test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  				          decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  // identify fermion and/or anti-fermion vertex
  AbstractFFVVertexPtr outgoingVertexF;
  AbstractFFVVertexPtr outgoingVertexA;
  identifyVertices(iferm, ianti, inpart, decay, outgoingVertexF, outgoingVertexA,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ifm = 0; ifm < 2; ++ifm) {
      for(unsigned int ia = 0; ia < 2; ++ia) {
	for(unsigned int ig = 0; ig < 2; ++ig) {
	  // radiation from the incoming vector
	  if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(incomingVertex_[inter]);
	    
	    VectorWaveFunction vectorInter = 
	      incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),vector3_[iv],
						gluon_[2*ig],inpart.mass());
	    
	    assert(vector3_[iv].particle()->id()==vectorInter.particle()->id());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,wave3_[ia],wavebar3_[ifm],vectorInter);
	    if(!couplingSet) {
              gs = abs(incomingVertex_[inter]->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	      (*ME[colourFlow[0][ix].first])(iv, ia, ifm, ig) += 
		 colourFlow[0][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
  
	  // radiation from outgoing fermion
	  if((decay[iferm]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[iferm]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexF);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[iferm]->dataPtr();
	    if(off->CC()) off = off->CC();
	    SpinorBarWaveFunction interS = 
	      outgoingVertexF->evaluate(scale,3,off,wavebar3_[ifm],
						gluon_[2*ig],decay[iferm]->mass());
	    
	    assert(wavebar3_[ifm].particle()->id()==interS.particle()->id());
	    
	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,wave3_[ia], interS,vector3_[iv]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertexF->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	      (*ME[colourFlow[1][ix].first])(iv, ia, ifm, ig) += 
		colourFlow[1][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }

	  // radiation from outgoing antifermion
	  if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	     (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	    assert(outgoingVertexA);
	    // ensure you get correct outgoing particle from first vertex
	    tcPDPtr off = decay[ianti]->dataPtr();
	    if(off->CC()) off = off->CC();
	    SpinorWaveFunction  interS = 
	      outgoingVertexA->evaluate(scale,3,off,wave3_[ia],
						gluon_[2*ig],decay[ianti]->mass());

	    assert(wave3_[ia].particle()->id()==interS.particle()->id());

	    Complex diag = 0.;
	    for(auto vertex : vertex_)
	      diag += vertex->evaluate(scale,interS,wavebar3_[ifm],vector3_[iv]);
	    if(!couplingSet) {
	      gs = abs(outgoingVertexA->norm());
	      couplingSet = true;
	    }
	    for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	      (*ME[colourFlow[2][ix].first])(iv, ia, ifm, ig) += 
		colourFlow[2][ix].second*diag;
	    }
#ifdef GAUGE_CHECK
	    total+=norm(diag);
#endif
	  }
	}
      }
    }
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  //return output
  return output;
}


void VFFDecayer::identifyVertices(const int iferm, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractFFVVertexPtr & outgoingVertexF, 
				  AbstractFFVVertexPtr & outgoingVertexA,
				  ShowerInteraction inter){
  // QCD vertices
  if(inter==ShowerInteraction::QCD) {
    // work out which fermion each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[iferm]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[iferm]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))){
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))){
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))){
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))){
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexF = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexF = outgoingVertex2_[inter];
      }
      else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour3 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr()))){
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[iferm]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexA = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (decay[iferm]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))){
	  outgoingVertexF = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexF = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour6 ||
	    inpart.dataPtr()->iColour()==PDT::Colour6bar) {
      if (outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr()))) {
	outgoingVertexF = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else {
	outgoingVertexF = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertexF  || outgoingVertexA)) ||
	   ( outgoingVertexF &&  outgoingVertexA)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in VFF decay in VFFDecayer::identifyVertices"
	<< Exception::runerror;
  }
  // QED
  else {
    if(decay[iferm]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iferm]->dataPtr())))
	outgoingVertexF = outgoingVertex1_[inter];
      else
	outgoingVertexF = outgoingVertex2_[inter];
    }
    if(decay[ianti]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr())))
	outgoingVertexA = outgoingVertex1_[inter];
      else
	outgoingVertexA = outgoingVertex2_[inter];
    }
  }
}

#line 1 "./VSSDecayer.cc"
// -*- C++ -*-
//
// VSSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VSSDecayer class.
//

#include "VSSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VSSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VSSDecayer::fullclone() const {
  return new_ptr(*this);
}

void VSSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> ) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVSSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VSSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
  }
}

void VSSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_           << perturbativeVertex_
     << incomingVertex_   << outgoingVertex1_
     << outgoingVertex2_;
}

void VSSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_           >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VSSDecayer,GeneralTwoBodyDecayer>
describeHerwigVSSDecayer("Herwig::VSSDecayer", "Herwig.so");

void VSSDecayer::Init() {

  static ClassDocumentation<VSSDecayer> documentation
    ("This implements the decay of a vector to 2 scalars");

}

void VSSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    ScalarWaveFunction::
      constructSpinInfo(decay[ix],outgoing,true);
}

double VSSDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  ScalarWaveFunction sca1(momenta[0],outgoing[0],Helicity::outgoing);
  ScalarWaveFunction sca2(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(ix,0,0) = 0.;
    for(auto vert : vertex_)
      (*ME())(ix,0,0) += vert->evaluate(scale,vectors_[ix],sca1,sca2);
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy VSSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, outa.first,
				     outb.first);
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    double me2 = 4.*sqr(pcm/inpart.second);
    Energy output = norm(perturbativeVertex_[0]->norm())*me2*pcm /
      (24.*Constants::pi);
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double VSSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  // work out which is the scalar and anti-scalar
  int ianti(0), iscal(1), iglu(2);
  int itype[2];
  for(unsigned int ix=0;ix<2;++ix) {
    if(decay[ix]->dataPtr()->CC()) itype[ix] = decay[ix]->id()>0 ? 0:1;
    else                           itype[ix] = 2;
  }
  if(itype[0]==0 && itype[1]!=0) swap(ianti, iscal);
  if(itype[0]==2 && itype[1]==1) swap(ianti, iscal);
  if(itype[0]==0 && itype[1]==0 && abs(decay[0]->dataPtr()->id())>abs(decay[1]->dataPtr()->id())) 
    swap(iscal, ianti);
  if(itype[0]==1 && itype[1]==1 && abs(decay[0]->dataPtr()->id())<abs(decay[1]->dataPtr()->id())) 
    swap(iscal, ianti);

 if(meopt==Initialize) {
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(vector3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
					       incoming, false);
  }
  // setup spin information when needed
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(vector3_ ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,false);
    ScalarWaveFunction::constructSpinInfo(       decay[iscal],outgoing,true);
    ScalarWaveFunction::constructSpinInfo(       decay[ianti],outgoing,true);
    VectorWaveFunction::constructSpinInfo(gluon_,decay[iglu ],outgoing,true,false);
    return 0.;
  }

  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1, PDT::Spin0,
								       PDT::Spin0, PDT::Spin1)));

  // create wavefunctions
  ScalarWaveFunction scal(decay[iscal]->momentum(), decay[iscal]->dataPtr(),outgoing);
  ScalarWaveFunction anti(decay[ianti]->momentum(), decay[ianti]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_,decay[iglu ],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[iglu ]->momentum(),
  					  decay[iglu ]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  // identify scalar and/or anti-scalar vertex
  AbstractVSSVertexPtr outgoingVertexS;
  AbstractVSSVertexPtr outgoingVertexA;
  identifyVertices(iscal, ianti, inpart, decay, outgoingVertexS, outgoingVertexA,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif
  for(unsigned int iv = 0; iv < 3; ++iv) {
    for(unsigned int ig = 0; ig < 2; ++ig) {
      // radiation from the incoming vector
      if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(incomingVertex_[inter]);	
	VectorWaveFunction vectorInter = 
	  incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),vector3_[iv],
					    gluon_[2*ig],inpart.mass());
	
	assert(vector3_[iv].particle()->id()==vectorInter.particle()->id());
	
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vectorInter,scal,anti);
	if(!couplingSet) {
	  gs = abs(incomingVertex_[inter]->norm());
	  couplingSet = true;
	}
	for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	  (*ME[colourFlow[0][ix].first])(iv, 0, 0, ig) += 
	    colourFlow[0][ix].second*diag;
	}
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
      }
      // radiation from the outgoing scalar
      if((decay[iscal]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (decay[iscal]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(outgoingVertexS);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[iscal]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter = 
	  outgoingVertexS->evaluate(scale,3,off,gluon_[2*ig],scal,decay[iscal]->mass());
	
	assert(scal.particle()->id()==scalarInter.particle()->id());
	
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vector3_[iv],anti,scalarInter);
	if(!couplingSet) {
	  gs = abs(outgoingVertexS->norm());
	  couplingSet = true;
	}
	for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	  (*ME[colourFlow[1][ix].first])(iv, 0, 0, ig) += 
	    colourFlow[1][ix].second*diag;
	}
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
      }
      
      if((decay[ianti]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	 (decay[ianti]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	assert(outgoingVertexA);
	// ensure you get correct outgoing particle from first vertex
	tcPDPtr off = decay[ianti]->dataPtr();
	if(off->CC()) off = off->CC();
	ScalarWaveFunction scalarInter = 
	  outgoingVertexA->evaluate(scale,3,off, gluon_[2*ig],anti,decay[ianti]->mass());
	
	assert(anti.particle()->id()==scalarInter.particle()->id());
	
	Complex diag = 0.;
	for(auto vertex : vertex_)
	  diag += vertex->evaluate(scale,vector3_[iv],scal,scalarInter);
	if(!couplingSet) {
	  gs = abs(outgoingVertexA->norm());
	  couplingSet = true;
	}
	for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	  (*ME[colourFlow[2][ix].first])(iv, 0, 0, ig) += 
	    colourFlow[2][ix].second*diag;
	}
#ifdef GAUGE_CHECK
      total+=norm(diag);
#endif
      }
    }
  }

  // contract matrices 
  double output=0.;
  for(unsigned int ix=0; ix<nflow; ++ix){
    for(unsigned int iy=0; iy<nflow; ++iy){
      output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
    }
  }
  // divide by alpha_(S,EM)
  output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
  double ratio = output/total;
  if(abs(ratio)>1e-20) {
    generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
    for(unsigned int ix=0;ix<decay.size();++ix)
      generator()->log() << *decay[ix] << "\n";
    generator()->log() << "Test of gauge invariance " << ratio << "\n";
  }
#endif
  // return the answer
  return output;
}


void VSSDecayer::identifyVertices(const int iscal, const int ianti,
				  const Particle & inpart, const ParticleVector & decay, 
				  AbstractVSSVertexPtr & outgoingVertexS, 
				  AbstractVSSVertexPtr & outgoingVertexA,
				  ShowerInteraction inter){
  if(inter==ShowerInteraction::QCD) {
    // work out which scalar each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[iscal]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[iscal]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[ianti]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[iscal]->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[ianti]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexS = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexS = outgoingVertex2_[inter];
      }
    else if (decay[iscal]->dataPtr()->iColour()==PDT::Colour3 &&
	     decay[ianti]->dataPtr()->iColour()==PDT::Colour8){
      if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[ianti]->dataPtr()->id()))){
	outgoingVertexS = outgoingVertex2_[inter];
	outgoingVertexA = outgoingVertex1_[inter];
      }
      else {
	outgoingVertexS = outgoingVertex1_[inter];
	outgoingVertexA = outgoingVertex2_[inter];
      }
    }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[iscal]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertexA = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertexA = outgoingVertex2_[inter];
      }
      else if (decay[iscal]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[ianti]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[iscal]->dataPtr()->id()))){
	  outgoingVertexS = outgoingVertex1_[inter];
	  outgoingVertexA = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertexS = outgoingVertex2_[inter];
	  outgoingVertexA = outgoingVertex1_[inter];
	}
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertexS  || outgoingVertexA)) ||
	   ( outgoingVertexS &&  outgoingVertexA)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in VSS decay in VSSDecayer::identifyVertices"
	<< Exception::runerror;
  }
  else {
    if(decay[iscal]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[iscal]->dataPtr())))
	outgoingVertexS = outgoingVertex1_[inter];
      else
	outgoingVertexS = outgoingVertex2_[inter];
    }
    if(decay[ianti]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[ianti]->dataPtr())))
	outgoingVertexA = outgoingVertex1_[inter];
      else
	outgoingVertexA = outgoingVertex2_[inter];
    }
  }
}
#line 1 "./VVSDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVSDecayer class.
//

#include "VVSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;
IBPtr VVSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVSDecayer::fullclone() const {
  return new_ptr(*this);
}
void VVSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractVVSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<VVSVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter]  = dynamic_ptr_cast<AbstractVVVVertexPtr>(inV.at(inter));
    outgoingVertexS_[inter] = AbstractVSSVertexPtr();
    outgoingVertexV_[inter] = AbstractVVVVertexPtr();  
    if(outV[0].at(inter)) {
      if (outV[0].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[0].at(inter));
      else
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    }
    if(outV[1].at(inter)) {
      if (outV[1].at(inter)->getName()==VertexType::VSS)
	outgoingVertexS_[inter]   = dynamic_ptr_cast<AbstractVSSVertexPtr>(outV[1].at(inter));
      else 
	outgoingVertexV_[inter]   = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    }
  }
}

void VVSDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_
     << incomingVertex_ << outgoingVertexS_ << outgoingVertexV_;
}

void VVSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_
     >> incomingVertex_ >> outgoingVertexS_ >> outgoingVertexV_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVSDecayer,GeneralTwoBodyDecayer>
describeHerwigVVSDecayer("Herwig::VVSDecayer", "Herwig.so");

void VVSDecayer::Init() {

  static ClassDocumentation<VVSDecayer> documentation
    ("The VVSDecayer class implements the decay of a vector"
     " to a vector and a scalar");

}

void VVSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool massless = ( decay[0]->id()==ParticleID::gamma || 
		    decay[0]->id()==ParticleID::g );
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  VectorWaveFunction::
    constructSpinInfo(vectors_[1],decay[0],outgoing,true,massless);
  ScalarWaveFunction::
    constructSpinInfo(decay[1],outgoing,true);
}
  
double VVSDecayer::me2(const int,const Particle & part,
		      const tPDVector & outgoing,
		      const vector<Lorentz5Momentum> & momenta,
		      MEOption meopt) const {
  bool massless = ( outgoing[0]->id()==ParticleID::gamma || 
		    outgoing[0]->id()==ParticleID::g );
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin0)));
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  VectorWaveFunction::
    calculateWaveFunctions(vectors_[1],momenta[0],outgoing[0],
			   Helicity::outgoing,massless);
  ScalarWaveFunction sca(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int in=0;in<3;++in) {
    for(unsigned int out=0;out<3;++out) {
      if(massless&&out==1) ++out;
      (*ME())(in,out,0) = 0.;
      for(auto vert : vertex_)
	(*ME())(in,out,0) += 
	  vert->evaluate(scale,vectors_[0][in],vectors_[1][out],sca);
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy VVSDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy2 scale(sqr(inpart.second));
    double mu1sq = sqr(outa.second/inpart.second);
    double mu2sq = sqr(outb.second/inpart.second);
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if( outb.first->iSpin() == PDT::Spin0 )
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, 
				       outa.first, outb.first);
    else {
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in, 
				       outb.first, outa.first);
      swap(mu1sq, mu2sq);
    }
    double vn = norm(perturbativeVertex_[0]->norm());
    if(vn == ZERO || mu1sq == ZERO) return ZERO;
    double me2 = 2. + 0.25*sqr(1. + mu1sq - mu2sq)/mu1sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy output = vn*me2*pcm/(24.*Constants::pi)/scale*UnitRemoval::E2;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double VVSDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  unsigned int ivec(0),isca(1);
  if(decay[ivec]->dataPtr()->iSpin()!=PDT::Spin1) swap(ivec,isca);
  if(meopt==Initialize) {
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(vectors3_[0], rho3_,
					       const_ptr_cast<tPPtr>(&inpart), 
					       incoming, false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(vectors3_[0] ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[1],decay[ivec],outgoing,true,false); 
    ScalarWaveFunction::
      constructSpinInfo(decay[isca],outgoing,true);
    VectorWaveFunction::
      constructSpinInfo(gluon_      ,decay[2],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1, PDT::Spin1,
								       PDT::Spin0, PDT::Spin1)));
  bool massless= decay[ivec]->mass()!=ZERO;
  // create wavefunctions
  VectorWaveFunction::calculateWaveFunctions(vectors3_[1],decay[0],outgoing,massless);
  ScalarWaveFunction scal(decay[isca]->momentum(),  decay[isca]->dataPtr(),outgoing);
  VectorWaveFunction::calculateWaveFunctions(gluon_      ,decay[2],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[2]->momentum(),
  					  decay[2]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  Energy2 scale(sqr(inpart.mass()));
  
  const GeneralTwoBodyDecayer::CFlow & colourFlow
    = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif

   for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
     for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
       if(massless && iv1==1) continue;
       for(unsigned int ig = 0; ig < 2; ++ig) {
	 // radiation from the incoming vector
	 if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(incomingVertex_[inter]);
	   VectorWaveFunction vectorInter = 
	     incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),vectors3_[0][iv0],
					      gluon_[2*ig],inpart.mass());
	   
	   assert(vectors3_[0][iv0].particle()->id()==vectorInter.particle()->id());
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectorInter,vectors3_[1][iv1],scal);
	   if(!couplingSet) {
	     gs = abs(incomingVertex_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	     (*ME[colourFlow[0][ix].first])(iv0, iv1, 0, ig) += 
	       colourFlow[0][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 // radiation from the outgoing vector
	 if((decay[ivec]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[ivec]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertexV_[inter]);
	   // ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[ivec]->dataPtr();
	   if(off->CC()) off = off->CC();
	   VectorWaveFunction vectorInter = 
	     outgoingVertexV_[inter]->evaluate(scale,3,off,gluon_[2*ig],vectors3_[1][iv1],decay[ivec]->mass());
	   
	   assert(vectors3_[1][iv1].particle()->id()==vectorInter.particle()->id());
	   
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv0],vectorInter,scal);
	   if(!couplingSet) {
	     gs = abs(outgoingVertexV_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	     (*ME[colourFlow[1][ix].first])(iv0, iv1, 0, ig) += 
	       colourFlow[1][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
	 // radiation from the outgoing scalar
	 if((decay[isca]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	    (decay[isca]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	   assert(outgoingVertexS_[inter]);
	   // ensure you get correct outgoing particle from first vertex
	   tcPDPtr off = decay[isca]->dataPtr();
	   if(off->CC()) off = off->CC();
	   ScalarWaveFunction scalarInter = 
	     outgoingVertexS_[inter]->evaluate(scale,3,off,gluon_[2*ig],scal,decay[isca]->mass());
	
	   assert(scal.particle()->id()==scalarInter.particle()->id());
	 
	   Complex diag = 0.;
	   for(auto vertex : vertex_)
	     diag += vertex->evaluate(scale,vectors3_[0][iv0],vectors3_[1][iv1],scalarInter);
	   if(!couplingSet) {
	     gs = abs(outgoingVertexS_[inter]->norm());
	     couplingSet = true;
	   }
	   for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	     (*ME[colourFlow[2][ix].first])(iv0, iv1, 0, ig) += 
	       colourFlow[2][ix].second*diag;
	   }
#ifdef GAUGE_CHECK
	   total+=norm(diag);
#endif
	 }
       }
     }
   }
   
   // contract matrices 
   double output=0.;
   for(unsigned int ix=0; ix<nflow; ++ix){
     for(unsigned int iy=0; iy<nflow; ++iy){
       output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
     }
   }
   // divide by alpha_(S,EM)
   output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
   double ratio = output/total;
   if(abs(ratio)>1e-20) {
     generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
     for(unsigned int ix=0;ix<decay.size();++ix)
       generator()->log() << *decay[ix] << "\n";
     generator()->log() << "Test of gauge invariance " << ratio << "\n";
   }
#endif
   // return the answer
   return output;
}
#line 1 "./VVVDecayer.cc"
// -*- C++ -*-
//
// VVVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVDecayer class.
//

#include "VVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr VVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void VVVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> & inV,
			      const vector<map<ShowerInteraction,VertexBasePtr> > & outV,
			      map<ShowerInteraction,VertexBasePtr> fourV) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractVVVVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<VVVVertexPtr>        (vert));
  }
  vector<ShowerInteraction> itemp={ShowerInteraction::QCD,ShowerInteraction::QED};
  for(auto & inter : itemp) {
    incomingVertex_[inter]  = dynamic_ptr_cast<AbstractVVVVertexPtr>(inV.at(inter));
    outgoingVertex1_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[0].at(inter));
    outgoingVertex2_[inter] = dynamic_ptr_cast<AbstractVVVVertexPtr>(outV[1].at(inter));
    fourPointVertex_[inter] = dynamic_ptr_cast<AbstractVVVVVertexPtr>(fourV.at(inter));
  }
}

void VVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_
     << incomingVertex_  << outgoingVertex1_
     << outgoingVertex2_ << fourPointVertex_;
}

void VVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_
     >> incomingVertex_   >> outgoingVertex1_
     >> outgoingVertex2_ >> fourPointVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VVVDecayer,GeneralTwoBodyDecayer>
describeHerwigVVVDecayer("Herwig::VVVDecayer", "Herwig.so");

void VVVDecayer::Init() {

  static ClassDocumentation<VVVDecayer> documentation
    ("The VVVDecayer class implements the decay of a vector boson "
     "into 2 vector bosons");

}

void VVVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix) 
    massless[ix] = (decay[ix]->id()==ParticleID::gamma ||
		    decay[ix]->id()==ParticleID::g);
  VectorWaveFunction::constructSpinInfo(vectors_[0],const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      constructSpinInfo(vectors_[ix+1],decay[ix],outgoing,true,massless[ix]);
}

double VVVDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin1,PDT::Spin1)));
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix) 
    massless[ix] = (outgoing[ix]->id()==ParticleID::gamma ||
		    outgoing[ix]->id()==ParticleID::g);
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_[0],rho_,
					       const_ptr_cast<tPPtr>(&part),
					       incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  for(unsigned int ix=0;ix<2;++ix)
    VectorWaveFunction::
      calculateWaveFunctions(vectors_[ix+1],momenta[ix],outgoing[ix],Helicity::outgoing,massless[ix]);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int iv3=0;iv3<3;++iv3) {
    for(unsigned int iv2=0;iv2<3;++iv2) {
      for(unsigned int iv1=0;iv1<3;++iv1) {
	(*ME())(iv1,iv2,iv3) = 0.;
	for(auto vert : vertex_) {
	  (*ME())(iv1,iv2,iv3) += vert->
	    evaluate(scale,vectors_[1][iv2],vectors_[2][iv3],vectors_[0][iv1]);
	}
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy VVVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), in,
					outa.first, outb.first);
    double mu1(outa.second/inpart.second), mu1sq(sqr(mu1)),
      mu2(outb.second/inpart.second), mu2sq(sqr(mu2));
    double vn = norm(perturbativeVertex_[0]->norm());
    if(vn == ZERO || mu1sq == ZERO || mu2sq == ZERO) return ZERO;
    double me2 = 
      (mu1 - mu2 - 1.)*(mu1 - mu2 + 1.)*(mu1 + mu2 - 1.)*(mu1 + mu2 + 1.)
      * (sqr(mu1sq) + sqr(mu2sq) + 10.*(mu1sq*mu2sq + mu1sq + mu2sq) + 1.)
      /4./mu1sq/mu2sq;
    Energy pcm = Kinematics::pstarTwoBodyDecay(inpart.second,outa.second,
					outb.second);
    Energy pWidth = vn*me2*pcm/24./Constants::pi;
    // colour factor
    pWidth *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return pWidth;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

double VVVDecayer::threeBodyME(const int , const Particle & inpart,
			       const ParticleVector & decay,
			       ShowerInteraction inter, MEOption meopt) {
  if(meopt==Initialize) {
    // create vector wavefunction for decaying particle
    VectorWaveFunction::calculateWaveFunctions(vector3_, rho3_, const_ptr_cast<tPPtr>(&inpart), 
					       incoming, false);
  }
  if(meopt==Terminate) {
    VectorWaveFunction::
      constructSpinInfo(vector3_ ,const_ptr_cast<tPPtr>(&inpart),outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[0],decay[0],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(vectors3_[1],decay[1],outgoing,true,false);
    VectorWaveFunction::
      constructSpinInfo(gluon_      ,decay[2],outgoing,true,false);
    return 0.;
  }
  // calculate colour factors and number of colour flows
  unsigned int nflow;
  vector<DVector> cfactors = getColourFactors(inpart, decay, nflow);
  vector<GeneralDecayMEPtr> ME(nflow,new_ptr(GeneralDecayMatrixElement(PDT::Spin1, PDT::Spin1,
								       PDT::Spin1, PDT::Spin1)));
  bool massless[2];
  for(unsigned int ix=0;ix<2;++ix)
    massless[ix] = decay[ix]->mass()!=ZERO;
  // create wavefunctions
  VectorWaveFunction::calculateWaveFunctions(vectors3_[0],decay[0],outgoing,massless[0]);
  VectorWaveFunction::calculateWaveFunctions(vectors3_[1],decay[1],outgoing,massless[1]);
  VectorWaveFunction::calculateWaveFunctions(gluon_      ,decay[2],outgoing,true);

  // gauge test
#ifdef GAUGE_CHECK
  gluon_.clear();
  for(unsigned int ix=0;ix<3;++ix) {
    if(ix==1) gluon_.push_back(VectorWaveFunction());
    else {
      gluon_.push_back(VectorWaveFunction(decay[2]->momentum(),
  					  decay[2]->dataPtr(),10,
  					  outgoing));
    }
  }
#endif

  // get the outgoing vertices
  AbstractVVVVertexPtr outgoingVertex1;
  AbstractVVVVertexPtr outgoingVertex2;
  identifyVertices(inpart,decay, outgoingVertex1, outgoingVertex2,inter);

  Energy2 scale(sqr(inpart.mass()));

  const GeneralTwoBodyDecayer::CFlow & colourFlow
        = colourFlows(inpart, decay);
  double gs(0.);
  bool couplingSet(false);
#ifdef GAUGE_CHECK
  double total=0.;
#endif

   for(unsigned int iv0 = 0; iv0 < 3; ++iv0) {
     for(unsigned int iv1 = 0; iv1 < 3; ++iv1) {
       if(massless[0] && iv1==1) continue;
       for(unsigned int iv2 = 0; iv2 < 3; ++iv2) {
	 if(massless[1] && iv2==1) continue;
	 for(unsigned int ig = 0; ig < 2; ++ig) {
	   // radiation from the incoming vector
	   if((inpart.dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	      (inpart.dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	     assert(incomingVertex_[inter]);
	     VectorWaveFunction vectorInter = 
	       incomingVertex_[inter]->evaluate(scale,3,inpart.dataPtr(),vector3_[iv0],
						gluon_[2*ig],inpart.mass());
	     
	     assert(vector3_[iv0].particle()->id()==vectorInter.particle()->id());
	     
	     Complex diag = 0.;
	     for(auto vertex : vertex_)
	       diag += vertex->evaluate(scale,vectorInter,vectors3_[0][iv1],
					vectors3_[1][iv2]);
	     if(!couplingSet) {
	       gs = abs(incomingVertex_[inter]->norm());
	       couplingSet = true;
	     }
	     for(unsigned int ix=0;ix<colourFlow[0].size();++ix) {
	       (*ME[colourFlow[0][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[0][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
	   }
	   // radiation from the 1st outgoing vector
	   if((decay[0]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	      (decay[0]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	     assert(outgoingVertex1);
	     // ensure you get correct outgoing particle from first vertex
	     tcPDPtr off = decay[0]->dataPtr();
	     if(off->CC()) off = off->CC();
	     VectorWaveFunction vectorInter = 
	       outgoingVertex1->evaluate(scale,3,off,gluon_[2*ig],vectors3_[0][iv1],decay[0]->mass());
	     
	     assert(vectors3_[0][iv1].particle()->id()==vectorInter.particle()->id());
	     
	     Complex diag = 0.;
	     for(auto vertex : vertex_)
	       diag += vertex->evaluate(scale,vector3_[iv0],vectorInter,vectors3_[1][iv2]);
	     if(!couplingSet) {
	       gs = abs(outgoingVertex1->norm());
	       couplingSet = true;
	     }
	     for(unsigned int ix=0;ix<colourFlow[1].size();++ix) {
	       (*ME[colourFlow[1][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[1][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
	   }
	   // radiation from the second outgoing vector
	   if((decay[1]->dataPtr()->coloured() && inter==ShowerInteraction::QCD) ||
	      (decay[1]->dataPtr()->charged()  && inter==ShowerInteraction::QED) ) {
	     assert(outgoingVertex2);
	     // ensure you get correct outgoing particle from first vertex
	     tcPDPtr off = decay[1]->dataPtr();
	     if(off->CC()) off = off->CC();
	     VectorWaveFunction vectorInter = 
	       outgoingVertex2->evaluate(scale,3,off, gluon_[2*ig],vectors3_[1][iv2],decay[1]->mass());
	     
	     assert(vectors3_[1][iv2].particle()->id()==vectorInter.particle()->id());
	     
	     Complex diag = 0.;
	     for(auto vertex : vertex_)
	       diag += vertex->evaluate(scale,vector3_[iv0],vectors3_[0][iv1],vectorInter);
	     if(!couplingSet) {
	       gs = abs(outgoingVertex2->norm());
	       couplingSet = true;
	     }
	     for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	       (*ME[colourFlow[2][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[2][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
	   }
	   // 4-point vertex
	   if (fourPointVertex_[inter]) {
	     Complex diag = fourPointVertex_[inter]->evaluate(scale,0, vector3_[iv0],
							      vectors3_[0][iv1],
							      vectors3_[1][iv2],gluon_[2*ig]);
	     for(unsigned int ix=0;ix<colourFlow[2].size();++ix) {
	       (*ME[colourFlow[3][ix].first])(iv0, iv1, iv2, ig) += 
		 colourFlow[3][ix].second*diag;
	     }
#ifdef GAUGE_CHECK
	     total+=norm(diag);
#endif
	   }
	 }
       }
     }
   }
   
   // contract matrices 
   double output=0.;
   for(unsigned int ix=0; ix<nflow; ++ix){
     for(unsigned int iy=0; iy<nflow; ++iy){
       output+=cfactors[ix][iy]*(ME[ix]->contract(*ME[iy],rho3_)).real();
     }
   }
   // divide by alpha_(S,EM)
   output*=(4.*Constants::pi)/sqr(gs);
#ifdef GAUGE_CHECK
   double ratio = output/total;
   if(abs(ratio)>1e-20) {
     generator()->log() << "Test of gauge invariance in decay\n" << inpart << "\n";
     for(unsigned int ix=0;ix<decay.size();++ix)
       generator()->log() << *decay[ix] << "\n";
     generator()->log() << "Test of gauge invariance " << ratio << "\n";
   }
#endif
   // return the answer
   return output;
}

void VVVDecayer::identifyVertices(const Particle & inpart, const ParticleVector & decay, 
				  AbstractVVVVertexPtr & outgoingVertex1, 
				  AbstractVVVVertexPtr & outgoingVertex2,
				  ShowerInteraction inter) {
  if(inter==ShowerInteraction::QCD) {
    // work out which scalar each outgoing vertex corresponds to 
    // two outgoing vertices
    if( inpart.dataPtr()       ->iColour()==PDT::Colour0     &&
	((decay[0]->dataPtr()->iColour()==PDT::Colour3     &&
	  decay[1]->dataPtr()->iColour()==PDT::Colour3bar) ||
	 (decay[0]->dataPtr()->iColour()==PDT::Colour8     &&
	  decay[1]->dataPtr()->iColour()==PDT::Colour8))){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex2_[inter];
	outgoingVertex2 = outgoingVertex1_[inter];
      }
    }
    else if(inpart.dataPtr()       ->iColour()==PDT::Colour8 &&
	    decay[0]->dataPtr()->iColour()==PDT::Colour3 &&
	    decay[1]->dataPtr()->iColour()==PDT::Colour3bar){
      if(outgoingVertex1_[inter]==outgoingVertex2_[inter]){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex1_[inter];
	outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (outgoingVertex2_[inter]->isIncoming(getParticleData(decay[0]->id()))){
	outgoingVertex1 = outgoingVertex2_[inter];
	outgoingVertex2 = outgoingVertex1_[inter];
      }
    }
    
    // one outgoing vertex
    else if(inpart.dataPtr()->iColour()==PDT::Colour3){
      if(decay[0]->dataPtr()->iColour()==PDT::Colour3 &&  
	 decay[1]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertex1 = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertex1 = outgoingVertex2_[inter];
      }
      else if (decay[0]->dataPtr()->iColour()==PDT::Colour3 &&
	       decay[1]->dataPtr()->iColour()==PDT::Colour8){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[1]->dataPtr()->id()))){
	  outgoingVertex1 = outgoingVertex2_[inter];
	  outgoingVertex2 = outgoingVertex1_[inter];
	}
	else {
	  outgoingVertex1 = outgoingVertex1_[inter];
	  outgoingVertex2 = outgoingVertex2_[inter];
	}
      }
    }
    else if(inpart.dataPtr()->iColour()==PDT::Colour3bar){
      if(decay[1]->dataPtr()->iColour()==PDT::Colour3bar &&  
	 decay[0]->dataPtr()->iColour()==PDT::Colour0){
	if     (outgoingVertex1_[inter]) outgoingVertex2 = outgoingVertex1_[inter];
	else if(outgoingVertex2_[inter]) outgoingVertex2 = outgoingVertex2_[inter];
      }
      else if (decay[0]->dataPtr()->iColour()==PDT::Colour8 &&
	       decay[1]->dataPtr()->iColour()==PDT::Colour3bar){
	if (outgoingVertex1_[inter]->isIncoming(getParticleData(decay[0]->dataPtr()->id()))){
	  outgoingVertex1 = outgoingVertex1_[inter];
	  outgoingVertex2 = outgoingVertex2_[inter];
	}
	else {
	  outgoingVertex1 = outgoingVertex2_[inter];
	  outgoingVertex2 = outgoingVertex1_[inter];
	}
      }
    }
    
    if (! ((incomingVertex_[inter]  && (outgoingVertex1  || outgoingVertex2)) ||
	   ( outgoingVertex1 &&  outgoingVertex2)))
      throw Exception()
	<< "Invalid vertices for QCD radiation in VVV decay in VVVDecayer::identifyVertices"
	<< Exception::runerror;
  }
  else {
    if(decay[0]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[0]->dataPtr())))
	outgoingVertex1 = outgoingVertex1_[inter];
      else
	outgoingVertex1 = outgoingVertex2_[inter];
    }
    if(decay[1]->dataPtr()->charged()) {
      if (outgoingVertex1_[inter] &&
	  outgoingVertex1_[inter]->isIncoming(const_ptr_cast<tPDPtr>(decay[1]->dataPtr())))
	outgoingVertex2 = outgoingVertex1_[inter];
      else
	outgoingVertex2 = outgoingVertex2_[inter];
    }
  }
}
#line 1 "./SRFDecayer.cc"
// -*- C++ -*-
//
// SRFDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SRFDecayer class.
//

#include "SRFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr SRFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SRFDecayer::fullclone() const {
  return new_ptr(*this);
}

void SRFDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > &,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_             .push_back(dynamic_ptr_cast<AbstractRFSVertexPtr>(vert));
    perturbativeVertex_ .push_back(dynamic_ptr_cast<RFSVertexPtr>        (vert));
  }
}

void SRFDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_;
}

void SRFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SRFDecayer,GeneralTwoBodyDecayer>
describeHerwigSRFDecayer("Herwig::SRFDecayer", "Herwig.so");

void SRFDecayer::Init() {

  static ClassDocumentation<SRFDecayer> documentation
    ("This class implements to decay of a scalar to a spin-3/2 and"
     " spin-1/2 fermion");

}

void SRFDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  unsigned int irs=0,ifm=1;
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) swap(irs,ifm);
  bool ferm = decay[ifm]->id()<0;
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),incoming,true);
  if(ferm) {
    RSSpinorBarWaveFunction::
      constructSpinInfo(RSwavebar_,decay[irs],outgoing,true);
    SpinorWaveFunction::
      constructSpinInfo(wave_     ,decay[ifm],outgoing,true);
  }
  else {
    RSSpinorWaveFunction::
      constructSpinInfo(RSwave_ ,decay[irs],outgoing,true);
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,decay[ifm],outgoing,true);
  }
}

double SRFDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  unsigned int irs=0,ifm=1;
  if(outgoing[0]->iSpin()==PDT::Spin1Half) swap(irs,ifm);
  if(!ME()) {
    if(irs==0)
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin3Half,PDT::Spin1Half)));
    else
      ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin1Half,PDT::Spin3Half)));
  }
  bool ferm = outgoing[ifm]->id()<0;
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(ferm) {
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(RSwavebar_,momenta[irs],outgoing[irs],Helicity::outgoing);
    SpinorWaveFunction::
      calculateWaveFunctions(wave_     ,momenta[ifm],outgoing[ifm],Helicity::outgoing);
  }
  else {
    RSSpinorWaveFunction::
      calculateWaveFunctions(RSwave_ ,momenta[irs],outgoing[irs],Helicity::outgoing);
    SpinorBarWaveFunction::
      calculateWaveFunctions(wavebar_,momenta[ifm],outgoing[ifm],Helicity::outgoing);
  }
  Energy2 scale(sqr(part.mass()));
  for(unsigned int ifm = 0; ifm < 4; ++ifm){
    for(unsigned int ia = 0; ia < 2; ++ia) {
      if(irs==0) {
	if(ferm) {
	  (*ME())(0, ifm, ia) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(0, ifm, ia) += vert->evaluate(scale,wave_[ia],
						  RSwavebar_[ifm],swave_);
	}
	else {
	  (*ME())(0, ifm, ia) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(0, ifm, ia) += vert->evaluate(scale,RSwave_[ifm],
						  wavebar_[ia],swave_);
	}
      }
      else {
	if(ferm) {
	  (*ME())(0, ia, ifm) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(0, ia, ifm) += vert->evaluate(scale,wave_[ia],
						  RSwavebar_[ifm],swave_);
	}
	else {
	  (*ME())(0, ia, ifm) = 0.;
	  for(auto vert : vertex_)
	    (*ME())(0, ia, ifm) += vert->evaluate(scale,RSwave_[ifm],
						  wavebar_[ia],swave_);
	}
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[irs],outgoing[ifm]);
  // return the answer
  return output;
}

Energy SRFDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy q = inpart.second;
    Energy m1 = outa.second, m2 = outb.second;
    // couplings
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    if(outa.first->iSpin()==PDT::Spin1Half) {
      swap(m1,m2);
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second),outb.first,
				       outa.first, in);
    }
    else {
      perturbativeVertex_[0]->setCoupling(sqr(inpart.second),outa.first,
				       outb.first, in);
    }
    Complex left  = perturbativeVertex_[0]-> left()*perturbativeVertex_[0]-> norm();
    Complex right = perturbativeVertex_[0]->right()*perturbativeVertex_[0]-> norm();
    complex<InvEnergy> A1 = 0.5*(left+right)*UnitRemoval::InvE;
    complex<InvEnergy> B1 = 0.5*(right-left)*UnitRemoval::InvE;
    Energy2 q2(q*q),m12(m1*m1),m22(m2*m2);
    Energy2 pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
    Energy pcm(sqrt(pcm2));
    Energy Qp(sqrt(-sqr(m2+m1)+q2)),Qm(sqrt(-sqr(m2-m1)+q2));
    double r23(sqrt(2./3.));
    complex<Energy> h1(-2.*r23*pcm*q/m1*Qm*B1);
    complex<Energy> h2( 2.*r23*pcm*q/m1*Qp*A1);
    double me2 = real(h1*conj(h1)+h2*conj(h2))/2./sqr(inpart.second);
    Energy output = me2*pcm/8./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

#line 1 "./FRSDecayer.cc"
// -*- C++ -*-
//
// FRSDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FRSDecayer class.
//

#include "FRSDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FRSDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FRSDecayer::fullclone() const {
  return new_ptr(*this);
}

void FRSDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > &,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractRFSVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<RFSVertexPtr>        (vert));
  }
}

void FRSDecayer::persistentOutput(PersistentOStream & os) const {
  os << perturbativeVertex_ << vertex_;
}

void FRSDecayer::persistentInput(PersistentIStream & is, int) {
  is >> perturbativeVertex_ >> vertex_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FRSDecayer,GeneralTwoBodyDecayer>
describeHerwigFRSDecayer("Herwig::FRSDecayer", "Herwig.so");

void FRSDecayer::Init() {

  static ClassDocumentation<FRSDecayer> documentation
    ("The FRSDecayer class implements the decay of a fermion to "
     "a spin-3/2 fermion and a scalar.");

}

void FRSDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool ferm = part.id() > 0;
  // for the decaying particle
  if(ferm) {
    SpinorWaveFunction::
      constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&part),incoming,true);
    RSSpinorBarWaveFunction::constructSpinInfo(RSwavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&part),incoming,true);
    RSSpinorWaveFunction::constructSpinInfo(RSwave_,decay[0],outgoing,true);
  }
  ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
}

double FRSDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  bool ferm = part.id() > 0;
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin0)));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  if(ferm)
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(RSwavebar_,momenta[0],outgoing[1],Helicity::outgoing);
  else
    RSSpinorWaveFunction::
      calculateWaveFunctions(RSwave_   ,momenta[0],outgoing[1],Helicity::outgoing);
  ScalarWaveFunction scal(momenta[1],outgoing[1],Helicity::outgoing);
  Energy2 scale(sqr(part.mass()));
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 4; ++if2) {
      (*ME())(if1, if2, 0) = 0.;
      for(auto vert : vertex_) {
	if(ferm) (*ME())(if1, if2, 0) +=
		   vert->evaluate(scale,wave_[if1],RSwavebar_[if2],scal);
	else     (*ME())(if1, if2, 0) += 
		   vert->evaluate(scale,RSwave_[if2],wavebar_[if1],scal);
      }
    }
  }
  double output = (ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // test code
//   Energy q = part.mass();
//   Energy m1 = decay[0]->mass();
//   Energy m2 = decay[1]->mass();
//   Energy2 q2(q*q),m12(m1*m1),m22(m2*m2);
//   Energy2 pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
//   Energy pcm(sqrt(pcm2));
//   Energy Qp(sqrt((q+m1)*(q+m1)-m22)),Qm(sqrt((q-m1)*(q-m1)-m22));
//   double r23(sqrt(2./3.));
//   // couplings
//   Complex left  = perturbativeVertex_-> left()*perturbativeVertex_-> norm();
//   Complex right = perturbativeVertex_->right()*perturbativeVertex_-> norm();
//   complex<InvEnergy> A1 = 0.5*(left+right)*UnitRemoval::InvE;
//   complex<InvEnergy> B1 = 0.5*(right-left)*UnitRemoval::InvE;
//   complex<Energy> h1(-2.*r23*pcm*q/m1*Qm*B1);
//   complex<Energy> h2( 2.*r23*pcm*q/m1*Qp*A1);
//   cout << "testing 1/2->3/2 0 "
//        << output*scale/GeV2 << "   " 
//        << real(h1*conj(h1)+h2*conj(h2))/4./GeV2     << "   " 
//        << real(h1*conj(h1)+h2*conj(h2))/4./(output*scale) << endl;
  // return the answer
  return output;
}

Energy FRSDecayer::partialWidth(PMPair inpart, PMPair outa,
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy q = inpart.second;
    Energy m1 = outa.second;
    Energy m2 = outb.second;
    Energy2 q2(q*q),m12(m1*m1),m22(m2*m2);
    Energy2 pcm2(0.25*(q2*(q2-2.*m12-2.*m22)+(m12-m22)*(m12-m22))/q2);
    Energy pcm(sqrt(pcm2));
    Energy Qp(sqrt((q+m1)*(q+m1)-m22)),Qm(sqrt((q-m1)*(q-m1)-m22));
    double r23(sqrt(2./3.));
    // couplings
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), outa.first, 
				     in,  outb.first);
    Complex left  = perturbativeVertex_[0]-> left()*perturbativeVertex_[0]-> norm();
    Complex right = perturbativeVertex_[0]->right()*perturbativeVertex_[0]-> norm();
    complex<InvEnergy> A1 = 0.5*(left+right)*UnitRemoval::InvE;
    complex<InvEnergy> B1 = 0.5*(right-left)*UnitRemoval::InvE;
    complex<Energy> h1(-2.*r23*pcm*q/m1*Qm*B1);
    complex<Energy> h2( 2.*r23*pcm*q/m1*Qp*A1);
    double me2 = real(h1*conj(h1)+h2*conj(h2))/4./sqr(inpart.second);
    Energy output = me2*pcm/8./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    // return the answer
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

#line 1 "./FRVDecayer.cc"
// -*- C++ -*-
//
// FRVDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FRVDecayer class.
//

#include "FRVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "Herwig/Utilities/Kinematics.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr FRVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FRVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FRVDecayer::setDecayInfo(PDPtr incoming, PDPair outgoing,
			      vector<VertexBasePtr> vertex,
			      map<ShowerInteraction,VertexBasePtr> &,
			      const vector<map<ShowerInteraction,VertexBasePtr> > &,
			      map<ShowerInteraction,VertexBasePtr>) {
  decayInfo(incoming,outgoing);
  for(auto vert : vertex) {
    vertex_            .push_back(dynamic_ptr_cast<AbstractRFVVertexPtr>(vert));
    perturbativeVertex_.push_back(dynamic_ptr_cast<RFVVertexPtr>        (vert));
  }
}

void FRVDecayer::persistentOutput(PersistentOStream & os) const {
  os << vertex_ << perturbativeVertex_;
}

void FRVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> vertex_ >> perturbativeVertex_;
}

void FRVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool ferm = part.id() > 0;
  // for the decaying particle
  if(ferm) {
    SpinorWaveFunction::
      constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&part),incoming,true);
    RSSpinorBarWaveFunction::constructSpinInfo(RSwavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&part),incoming,true);
    RSSpinorWaveFunction::constructSpinInfo(RSwave_,decay[0],outgoing,true);
  }
  VectorWaveFunction::
    constructSpinInfo(vector_,decay[1],outgoing,true,false);
}

double FRVDecayer::me2(const int,const Particle & part,
		       const tPDVector & outgoing,
		       const vector<Lorentz5Momentum> & momenta,
		       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin1)));
  // decaying fermion or antifermion
  bool ferm = part.id() > 0;
  // initialize
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  Energy2 scale(sqr(part.mass()));
  if(ferm)
    RSSpinorBarWaveFunction::
      calculateWaveFunctions(RSwavebar_,momenta[0],outgoing[1],Helicity::outgoing);
  else
    RSSpinorWaveFunction::
      calculateWaveFunctions(RSwave_   ,momenta[0],outgoing[1],Helicity::outgoing);
  bool massless = outgoing[1]->mass()==ZERO;
  VectorWaveFunction::
    calculateWaveFunctions(vector_,momenta[1],outgoing[1],Helicity::outgoing,massless);
  // loop over helicities
  for(unsigned int if1 = 0; if1 < 2; ++if1) {
    for(unsigned int if2 = 0; if2 < 4; ++if2) {
      for(unsigned int vhel = 0; vhel < 3; ++vhel) {
	if(massless && vhel == 1) ++vhel;
	(*ME())(if1, if2,vhel) = 0.;
	for(auto vert : vertex_) {
	  if(ferm)
	    (*ME())(if1, if2,vhel) += 
	      vert->evaluate(scale,wave_[if1],
			     RSwavebar_[if2],vector_[vhel]);
	  else
	    (*ME())(if1, if2, vhel) += 
	      vert->evaluate(scale,RSwave_[if2],
			     wavebar_[if1],vector_[vhel]);
	}
      }
    }
  }
  double output=(ME()->contract(rho_)).real()/scale*UnitRemoval::E2;
  // test
//   Energy m1(part.mass()),m2(decay[0]->mass()),m3(decay[1]->mass());
//   Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
//   Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
//   double r2(sqrt(2.)),r3(sqrt(3.));
//   Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
//   vector<Complex> left  = perturbativeVertex_-> left();
//   vector<Complex> right = perturbativeVertex_->right();
//   Complex A1 = 0.5*(left [0]+right[0])*perturbativeVertex_-> norm();
//   Complex B1 = 0.5*(right[0]- left[0])*perturbativeVertex_-> norm();
//   complex<InvEnergy> A2 = 0.5*(left [1]+right[1])*perturbativeVertex_-> norm()*UnitRemoval::InvE;
//   complex<InvEnergy> B2 = 0.5*(right[1]- left[1])*perturbativeVertex_-> norm()*UnitRemoval::InvE;
//   complex<InvEnergy2> A3 = 0.5*(left [2]+right[2])*perturbativeVertex_-> norm()*UnitRemoval::InvE2;
//   complex<InvEnergy2> B3 = 0.5*(right[2]- left[2])*perturbativeVertex_-> norm()*UnitRemoval::InvE2;
//   complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
//   complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2));
//   complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2));
//   complex<Energy> h5(ZERO),h6(ZERO);
//   if(decay[1]->mass()>ZERO) {
//     h5 = -2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2
//     			   +m12*pcm*pcm*A3);
//     h6 =  2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2
// 					   +m12*pcm*pcm*B3);
//   }
//   cout << "testing 1/2->3/2 1 " << part.id() << " "
//        << output << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass()) << "   " 
//        << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
// 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass())/output << endl;
  // colour and identical particle factors
  output *= colourFactor(part.dataPtr(),outgoing[0],outgoing[1]);
  // return the answer
  return output;
}

Energy FRVDecayer::partialWidth(PMPair inpart, PMPair outa, 
				PMPair outb) const {
  if( inpart.second < outa.second + outb.second  ) return ZERO;
  if(perturbativeVertex_.size()==1 &&
     perturbativeVertex_[0]) {
    Energy m1(inpart.second),m2(outa.second),m3(outb.second);
    Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
    Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
    double r2(sqrt(2.)),r3(sqrt(3.));
    Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
    // couplings
    tcPDPtr in = inpart.first->CC() ? tcPDPtr(inpart.first->CC()) : inpart.first;
    perturbativeVertex_[0]->setCoupling(sqr(inpart.second), outa.first, 
				     in,  outb.first);
    vector<Complex> left  = perturbativeVertex_[0]-> left();
    vector<Complex> right = perturbativeVertex_[0]->right();
    Complex A1 = 0.5*(left [0]+right[0])*perturbativeVertex_[0]-> norm();
    Complex B1 = 0.5*(right[0]- left[0])*perturbativeVertex_[0]-> norm();
    complex<InvEnergy> A2 = 0.5*(left [1]+right[1])*perturbativeVertex_[0]-> norm()*UnitRemoval::InvE;
    complex<InvEnergy> B2 = 0.5*(right[1]- left[1])*perturbativeVertex_[0]-> norm()*UnitRemoval::InvE;
    complex<InvEnergy2> A3 = 0.5*(left [2]+right[2])*perturbativeVertex_[0]-> norm()*UnitRemoval::InvE2;
    complex<InvEnergy2> B3 = 0.5*(right[2]- left[2])*perturbativeVertex_[0]-> norm()*UnitRemoval::InvE2;
    complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
    complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2));
    complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2));
    complex<Energy> h5(ZERO),h6(ZERO);
    if(outb.second>ZERO) {
      h5 = -2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2
			       +m12*pcm*pcm*A3);
      h6 =  2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2
			       +m12*pcm*pcm*B3);
    }
    double me2 = 0.25*real(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
			   h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(inpart.second);
    Energy output = me2*pcm/8./Constants::pi;
    // colour factor
    output *= colourFactor(inpart.first,outa.first,outb.first);
    return output;
  }
  else {
    return GeneralTwoBodyDecayer::partialWidth(inpart,outa,outb);
  }
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FRVDecayer,GeneralTwoBodyDecayer>
describeHerwigFRVDecayer("Herwig::FRVDecayer", "Herwig.so");

void FRVDecayer::Init() {

  static ClassDocumentation<FRVDecayer> documentation
    ("The FRVDecayer class handles the decay of a fermion to "
     "a spin-3/2 particle and a vector boson.");

}
#line 1 "./FtoFFFDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFFFDecayer class.
//

#include "FtoFFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;

IBPtr FtoFFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FtoFFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void FtoFFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << vec_ << ten_;
}

void FtoFFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> vec_ >> ten_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FtoFFFDecayer,GeneralThreeBodyDecayer>
describeHerwigFtoFFFDecayer("Herwig::FtoFFFDecayer", "Herwig.so");

void FtoFFFDecayer::Init() {

  static ClassDocumentation<FtoFFFDecayer> documentation
    ("The FtoFFFDecayer class implements the general decay of a fermion to "
     "three fermions.");

}

void FtoFFFDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in FtoFFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in FtoFFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
  }
}

void FtoFFFDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()<0) 
    SpinorWaveFunction::constructSpinInfo(inwave_.first,
					  const_ptr_cast<tPPtr>(&part),
					  Helicity::incoming,true);
  else
    SpinorBarWaveFunction::constructSpinInfo(inwave_.second,
					     const_ptr_cast<tPPtr>(&part),
					     Helicity::incoming,true);
  // outgoing particles
  for(unsigned int ix = 0; ix < 3; ++ix) {
    SpinorWaveFunction::
      constructSpinInfo(outwave_[ix].first,decay[ix],Helicity::outgoing,true);
  }
}

double  FtoFFFDecayer::me2(const int ichan,const Particle & part,
			   const tPDVector & outgoing,
			   const vector<Lorentz5Momentum> & momenta,
			   MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != part.id();
  // special handling or first/last call
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  Energy2 scale(sqr(part.mass()));
  if(meopt==Initialize) {
    SpinorWaveFunction::
      calculateWaveFunctions(inwave_.first,rho_,const_ptr_cast<tPPtr>(&part),
			    Helicity::incoming);
    inwave_.second.resize(2);
    if(inwave_.first[0].wave().Type() == SpinorType::u) {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	inwave_.second[ix] = inwave_.first[ix].bar();
	inwave_.second[ix].conjugate();
      }
    }
    else {
      for(unsigned int ix = 0; ix < 2; ++ix) {
	inwave_.second[ix] = inwave_.first[ix].bar();
	inwave_.first[ix].conjugate();
      }
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  // outgoing particles
  for(unsigned int ix = 0; ix < 3; ++ix) {
    SpinorWaveFunction::
      calculateWaveFunctions(outwave_[ix].first,momenta[ix],outgoing[ix],Helicity::outgoing);
    outwave_[ix].second.resize(2);
    if(outwave_[ix].first[0].wave().Type() == SpinorType::u) {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave_[ix].second[iy] = outwave_[ix].first[iy].bar();
	outwave_[ix].first[iy].conjugate();
      }
    }
    else {
      for(unsigned int iy = 0; iy < 2; ++iy) {
	outwave_[ix].second[iy] = outwave_[ix].first[iy].bar();
	outwave_[ix].second[iy].conjugate();
      }
    }
  }
  bool ferm = part.id()>0; 
  vector<Complex> flows(ncf, Complex(0.)),largeflows(ncf, Complex(0.)); 
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  vector<GeneralDecayMEPtr> mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,
								      PDT::Spin1Half,PDT::Spin1Half)));
  unsigned int ihel[4];
  for(ihel[0] = 0; ihel[0] < 2; ++ihel[0]) {
    for(ihel[1] = 0; ihel[1] < 2; ++ihel[1]) {
      for(ihel[2] = 0; ihel[2] < 2; ++ihel[2]) {
	for(ihel[3] = 0; ihel[3] < 2; ++ihel[3]) {
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag=0;
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    if(ichan>=0&&diagramMap()[ichan]!=idiag) {
	      ++idiag;
	      continue;
	    }
	    // the sign from normal ordering
	    double sign = ferm ? 1. : -1;
	    // outgoing wavefunction and NO sign
	    if     (dit->channelType==TBDiagram::channel23) sign *= -1.;
	    else if(dit->channelType==TBDiagram::channel13) sign *=  1.;
	    else if(dit->channelType==TBDiagram::channel12) sign *= -1.;
	    else throw Exception()
	      << "Unknown diagram type in FtoFFFDecayer::me2()" << Exception::runerror;
	    // wavefunctions
	    SpinorWaveFunction    w0,w3;
	    SpinorBarWaveFunction w1,w2;
	    // incoming wavefunction
	    if(ferm) {
	      w0 = inwave_.first [ihel[0]];
	      w1 = outwave_[dit->channelType].second[ihel[dit->channelType+1]];
	    }
	    else {
	      w0 = outwave_[dit->channelType].first [ihel[dit->channelType+1]];
	      w1 = inwave_.second[ihel[0]];
	    }
	    if(outgoing[out2[dit->channelType]]->id()<0&&
	       outgoing[out3[dit->channelType]]->id()>0) {
	      w2 = outwave_[out3[dit->channelType]].second[ihel[out3[dit->channelType]+1]];
	      w3 = outwave_[out2[dit->channelType]].first [ihel[out2[dit->channelType]+1]];
	      sign *= -1.;
	    }
	    else {
	      w2 = outwave_[out2[dit->channelType]].second[ihel[out2[dit->channelType]+1]];
	      w3 = outwave_[out3[dit->channelType]].first [ihel[out3[dit->channelType]+1]];
	    }
	    tcPDPtr offshell = dit->intermediate;
	    if(cc&&offshell->CC()) offshell=offshell->CC();
	    Complex diag(0.);
	    // intermediate scalar
	    if     (offshell->iSpin() == PDT::Spin0) { 
	      ScalarWaveFunction inters = sca_[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = sca_[idiag].second->evaluate(scale,w3,w2,inters);
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = vec_[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = vec_[idiag].second->evaluate(scale,w3,w2,interv);
	    }
	    // intermediate tensor
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction intert = ten_[idiag].first->
		evaluate(scale, widthOption(), offshell, w0, w1);
	      diag = ten_[idiag].second->evaluate(scale,w3,w2,intert);
	    }
	    // unknown
	    else throw Exception()
	      << "Unknown intermediate in FtoFFFDecayer::me2()" 
	      << Exception::runerror;
	    // apply NO sign
	    diag *= sign;
	    // matrix element for the different colour flows
	    if(ichan<0) {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }
	    else {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }
	    ++idiag;
	  }
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    (*mes[ix])(ihel[0],ihel[1],ihel[2],ihel[3]) =      flows[ix];
	    (*mel[ix])(ihel[0],ihel[1],ihel[2],ihel[3]) = largeflows[ix];
	  }
	}
      }
    }
  }
  double me2(0.);
  if(ichan<0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *=UseRandom::rnd();
    for(unsigned int ix=0;ix<pflows.size();++ix) {
      if(ptotal<=pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal-=pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}

WidthCalculatorBasePtr FtoFFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<FtoFFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass(),
		 relativeError()));
}
#line 1 "./StoSFFDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoSFFDecayer class.
//

#include "StoSFFDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

IBPtr StoSFFDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoSFFDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoSFFDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << fer_ << vec_ << ten_ << RSfer_ << four_;
}

void StoSFFDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> fer_ >> vec_ >> ten_ >> RSfer_ >> four_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StoSFFDecayer,GeneralThreeBodyDecayer>
describeHerwigStoSFFDecayer("Herwig::StoSFFDecayer", "Herwig.so");

void StoSFFDecayer::Init() {

  static ClassDocumentation<StoSFFDecayer> documentation
    ("The StoSFFDecayer class implements the general decay of a scalar to "
     "a scalar and two fermions.");

}

WidthCalculatorBasePtr StoSFFDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<StoSFFDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),outgoing()[2]->mass(),
		  relativeError()));
}

void StoSFFDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  fer_.resize(ndiags);
  RSfer_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  four_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    // four point vertex
    if(!offshell) {
      four_[ix] = dynamic_ptr_cast<AbstractFFSSVertexPtr>(current.vertices.first);
      continue;
    }
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractSSSVertexPtr vert1 = dynamic_ptr_cast<AbstractSSSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      fer_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVSSVertexPtr vert1 = dynamic_ptr_cast<AbstractVSSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractSSTVertexPtr vert1 = dynamic_ptr_cast<AbstractSSTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin3Half) {
      AbstractRFSVertexPtr vert1 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	(current.vertices.first);
      AbstractRFSVertexPtr vert2 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a RS fermion diagram in StoSFFDecayer::setupDiagrams()"
	<< Exception::runerror;
      RSfer_[ix] = make_pair(vert1, vert2);
    }
  }
}

void StoSFFDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),
		      Helicity::incoming,true);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin()==PDT::Spin0) {
      ScalarWaveFunction::constructSpinInfo(decay[ix],Helicity::outgoing,true);
    }
    else {
      SpinorWaveFunction::
	constructSpinInfo(outspin_[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
}

double StoSFFDecayer::me2(const int ichan, const Particle & part,
			  const tPDVector & outgoing,
			  const vector<Lorentz5Momentum> & momenta,
			  MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != part.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),
			     Helicity::incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),
				Helicity::incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  // get the wavefunctions for all the particles
  ScalarWaveFunction outScalar;
  unsigned int isca(0);
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    if(outgoing[ix]->iSpin()==PDT::Spin0) {
      isca = ix;
      outScalar = ScalarWaveFunction(momenta[ix],outgoing[ix],Helicity::outgoing);
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(outspin_[ix].first,momenta[ix],outgoing[ix],Helicity::outgoing);
      outspin_[ix].second.resize(2);
      if(outspin_[ix].first[0].wave().Type() == SpinorType::u) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].second[iy].conjugate();
	}
      }
    }
  }
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  Energy2 scale(sqr(part.mass()));  
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.)); 
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      isca==0 ? PDT::Spin0 : PDT::Spin1Half,
					      isca==1 ? PDT::Spin0 : PDT::Spin1Half,
					      isca==2 ? PDT::Spin0 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      isca == 0 ? PDT::Spin0 : PDT::Spin1Half,
					      isca == 1 ? PDT::Spin0 : PDT::Spin1Half,
					      isca == 2 ? PDT::Spin0 : PDT::Spin1Half)));
  static const unsigned int out2[3]={1,0,0},out3[3]={2,2,1};
  for(unsigned int s1 = 0;s1 < 2; ++s1) {
    for(unsigned int s2 = 0;s2 < 2; ++s2) {
      flows = vector<Complex>(ncf, Complex(0.));
      largeflows = vector<Complex>(ncf, Complex(0.));
      unsigned int idiag(0);
      for(vector<TBDiagram>::const_iterator dit = getProcessInfo().begin();
	  dit != getProcessInfo().end(); ++dit) {
	// channels if selecting
	if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	  ++idiag;
	  continue;
	}
	tcPDPtr offshell = dit->intermediate;
	Complex diag;
	if(offshell) {
	  if(cc&&offshell->CC()) offshell=offshell->CC();
	  double sign = out3[dit->channelType] < out2[dit->channelType] ? 1. : -1.;
	  // intermediate scalar
	  if     (offshell->iSpin() == PDT::Spin0) {
	    ScalarWaveFunction inters = sca_[idiag].first->
	      evaluate(scale, widthOption(), offshell, swave_, outScalar);
	    unsigned int h1(s1),h2(s2);
	    if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	    if(outgoing[out2[dit->channelType]]->id()<0&&
	       outgoing[out3[dit->channelType]]->id()>0) {
	      diag =-sign*sca_[idiag].second->
		evaluate(scale,
			 outspin_[out2[dit->channelType]].first [h1],
			 outspin_[out3[dit->channelType]].second[h2],inters);
	    }
	    else {
	      diag = sign*sca_[idiag].second->
		evaluate(scale,
			 outspin_[out3[dit->channelType]].first [h2],
			 outspin_[out2[dit->channelType]].second[h1],inters);
	    }
	  }
	  // intermediate fermion
	  else if(offshell->iSpin() == PDT::Spin1Half) {
	    int iferm = 
	      outgoing[out2[dit->channelType]]->iSpin()==PDT::Spin1Half
	      ? out2[dit->channelType] : out3[dit->channelType];
	    unsigned int h1(s1),h2(s2);
	    if(dit->channelType>iferm) swap(h1,h2);
	    sign = iferm<dit->channelType ? 1. : -1.;
	    
	    
	    if((outgoing[dit->channelType]->id() < 0 &&outgoing[iferm]->id() > 0 ) ||
	       (outgoing[dit->channelType]->id()*offshell->id()>0)) {
	      SpinorWaveFunction    inters = fer_[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 outspin_[dit->channelType].first [h1],swave_);
	      diag = -sign*fer_[idiag].second->
		evaluate(scale,inters,outspin_[iferm].second[h2],outScalar);
	    }
	    else {
	    SpinorBarWaveFunction inters = fer_[idiag].first->
	      evaluate(scale,widthOption(),offshell,
		       outspin_[dit->channelType].second[h1],swave_);
	    diag =  sign*fer_[idiag].second->
	      evaluate(scale,outspin_[iferm].first [h2],inters,outScalar);
	    }
	  }
	  // intermediate vector
	  else if(offshell->iSpin() == PDT::Spin1) {
	    VectorWaveFunction interv = vec_[idiag].first->
	      evaluate(scale, widthOption(), offshell, swave_, outScalar);
	    unsigned int h1(s1),h2(s2);
	    if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	  if(outgoing[out2[dit->channelType]]->id()<0&&
	     outgoing[out3[dit->channelType]]->id()>0) {
	    diag =-sign*vec_[idiag].second->
	      evaluate(scale,
		       outspin_[out2[dit->channelType]].first [h1],
		       outspin_[out3[dit->channelType]].second[h2],interv);
	  }
	  else {
	    diag = sign*vec_[idiag].second->
	      evaluate(scale,
		       outspin_[out3[dit->channelType]].first [h2],
		       outspin_[out2[dit->channelType]].second[h1],interv);
	  }
	  }
	  // intermediate tensor
	  else if(offshell->iSpin() == PDT::Spin2) {
	    TensorWaveFunction intert = ten_[idiag].first->
	      evaluate(scale, widthOption(), offshell, swave_, outScalar);
	    unsigned int h1(s1),h2(s2);
	    if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	    if(outgoing[out2[dit->channelType]]->id()<0&&
	       outgoing[out3[dit->channelType]]->id()>0) {
	      diag =-sign*ten_[idiag].second->
		evaluate(scale,
			 outspin_[out2[dit->channelType]].first [h1],
			 outspin_[out3[dit->channelType]].second[h2],intert);
	    }
	    else {
	      diag = sign*ten_[idiag].second->
		evaluate(scale,
			 outspin_[out3[dit->channelType]].first [h2],
			 outspin_[out2[dit->channelType]].second[h1],intert);
	    }
	  }
	  // intermediate RS fermion
	  else if(offshell->iSpin() == PDT::Spin3Half) {
	    int iferm = 
	      outgoing[out2[dit->channelType]]->iSpin()==PDT::Spin1Half
	      ? out2[dit->channelType] : out3[dit->channelType];
	    unsigned int h1(s1),h2(s2);
	    if(dit->channelType>iferm) swap(h1,h2);
	    sign = iferm<dit->channelType ? 1. : -1.;
	    if((outgoing[dit->channelType]->id() < 0 &&outgoing[iferm]->id() > 0 ) ||
	       (outgoing[dit->channelType]->id()*offshell->id()>0)) {
	      RSSpinorWaveFunction    inters = RSfer_[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 outspin_[dit->channelType].first [h1],swave_);
	      diag = -sign*RSfer_[idiag].second->
		evaluate(scale,inters,outspin_[iferm].second[h2],outScalar);
	    }
	    else {
	      RSSpinorBarWaveFunction inters = RSfer_[idiag].first->
		evaluate(scale,widthOption(),offshell,
			 outspin_[dit->channelType].second[h1],swave_);
	      diag =  sign*RSfer_[idiag].second->
		evaluate(scale,outspin_[iferm].first [h2],inters,outScalar);
	    }
	  }
	  // unknown
	  else throw Exception()
		 << "Unknown intermediate in StoSFFDecayer::me2()" 
		 << Exception::runerror;
	}
	// four point diagram
	else {
	    unsigned int o2 = isca > 0 ? 0 : 1;
	    unsigned int o3 = isca < 2 ? 2 : 1;
	    if(outgoing[o2]->id() < 0 && outgoing[o3]->id() > 0) {
	      diag =-four_[idiag]->
	    	evaluate(scale, outspin_[o2].first[s1],
	    		 outspin_[o3].second[s2], outScalar, swave_);
	    }
	    else {
	      diag = four_[idiag]->
	    	evaluate(scale, outspin_[o3].first[s2],
	    		 outspin_[o2].second[s1], outScalar, swave_);
	    }
	}
	// matrix element for the different colour flows
	if(ichan < 0) {
	  for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	    flows[dit->colourFlow[iy].first - 1] += 
	      dit->colourFlow[iy].second * diag;
	  }
	  for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	    largeflows[dit->largeNcColourFlow[iy].first - 1] += 
	      dit->largeNcColourFlow[iy].second * diag;
	  }
	}
	else {
	  for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	    if(dit->colourFlow[iy].first - 1 != colourFlow()) continue;
	    flows[dit->colourFlow[iy].first - 1] += 
	      dit->colourFlow[iy].second * diag;
	  }
	  for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	    if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
	    largeflows[dit->largeNcColourFlow[iy].first - 1] += 
	      dit->largeNcColourFlow[iy].second * diag;
	  }
	}
	++idiag;
      }
      for(unsigned int ix = 0; ix < ncf; ++ix) {
	if(isca == 0) {
	  (*mes[ix])(0, 0, s1, s2) = flows[ix];
	  (*mel[ix])(0, 0, s1, s2) = largeflows[ix];
	}
	else if(isca == 1 ) { 
	  (*mes[ix])(0, s1, 0, s2) = flows[ix];
	  (*mel[ix])(0, s1, 0, s2) = largeflows[ix];
	}
	else if(isca == 2) { 
	  (*mes[ix])(0, s1,s2, 0) = flows[ix];
	  (*mel[ix])(0, s1,s2, 0) = largeflows[ix] ;
	}
      }
    }
  }
  double me2(0.);
  if(ichan < 0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix == iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix = 0;ix < pflows.size(); ++ix) {
      if(ptotal <= pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal -= pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}
#line 1 "./StoFFVDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StoFFVDecayer class.
//

#include "StoFFVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

IBPtr StoFFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr StoFFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void StoFFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << fer_ << vec_ << RSfer_ << four_;
}

void StoFFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> fer_ >> vec_ >> RSfer_ >> four_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StoFFVDecayer,GeneralThreeBodyDecayer>
describeHerwigStoFFVDecayer("Herwig::StoFFVDecayer", "Herwig.so");

void StoFFVDecayer::Init() {

  static ClassDocumentation<StoFFVDecayer> documentation
    ("The StoFFVDecayer class implements the general decay of a scalar to "
     "a two fermions and a vector.");

}

WidthCalculatorBasePtr StoFFVDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<StoFFVDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),
		  outgoing()[2]->mass(),relativeError()));
}

void StoFFVDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  fer_.resize(ndiags);
  RSfer_.resize(ndiags);
  vec_.resize(ndiags);
  four_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    // four point vertex
    if(!offshell) {
      four_[ix] = dynamic_ptr_cast<AbstractFFVSVertexPtr>(current.vertices.first);
      continue;
    }
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractVSSVertexPtr vert1 = dynamic_ptr_cast<AbstractVSSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in StoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in StoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      fer_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVSVertexPtr vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in StoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin3Half) {
      AbstractRFSVertexPtr vert1 = dynamic_ptr_cast<AbstractRFSVertexPtr>
	(current.vertices.first);
      AbstractRFVVertexPtr vert2 = dynamic_ptr_cast<AbstractRFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a RS fermion diagram in StoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      RSfer_[ix] = make_pair(vert1, vert2);
    }
  }
}

void StoFFVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  ScalarWaveFunction::
    constructSpinInfo(const_ptr_cast<tPPtr>(&part),
		      Helicity::incoming,true);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin()==PDT::Spin1) {
      VectorWaveFunction::constructSpinInfo(outVector_,decay[ix],
					    Helicity::outgoing,true,false);
    }
    else {
      SpinorWaveFunction::
	constructSpinInfo(outspin_[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
}

double StoFFVDecayer::me2(const int ichan, const Particle & part,
			  const tPDVector & outgoing,
			  const vector<Lorentz5Momentum> & momenta,
			  MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != part.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),
			     Helicity::incoming);
    swave_ = ScalarWaveFunction(part.momentum(),part.dataPtr(),
				Helicity::incoming);
    // fix rho if no correlations
    fixRho(rho_);
  }
  unsigned int ivec(0);
  bool massless(false);
  for(unsigned int ix = 0; ix < outgoing.size();++ix) {
    if(outgoing[ix]->iSpin() == PDT::Spin1) {
      ivec = ix;
      massless = outgoing[ivec]->mass()==ZERO;
      VectorWaveFunction::
	calculateWaveFunctions(outVector_, momenta[ix],outgoing[ix], Helicity::outgoing,massless);
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(outspin_[ix].first,momenta[ix],outgoing[ix],Helicity::outgoing);
      outspin_[ix].second.resize(2);
      // Need a ubar and a v spinor
      if(outspin_[ix].first[0].wave().Type() == SpinorType::u) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].second[iy].conjugate();
	}
      }
    }
  }
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  Energy2 scale(sqr(part.mass()));
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.));
  // setup the DecayMatrixElement
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin0,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  //the channel possiblities
  static const unsigned int out2[3] = {1,0,0}, out3[3] = {2,2,1};
  for(unsigned int s1 = 0; s1 < 2; ++s1) {
    for(unsigned int s2 = 0; s2 < 2; ++s2) {
      for(unsigned int v1 = 0; v1 < 3; ++v1) {
	if(massless&&v1==1) ++v1;
	flows = vector<Complex>(ncf, Complex(0.));
	largeflows = vector<Complex>(ncf, Complex(0.));
	unsigned int idiag(0);
	Complex diag;
	for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	    dit!=getProcessInfo().end();++dit) {
	  // channels if selecting
	  if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	    ++idiag;
	    continue;
	  }
	  tcPDPtr offshell = dit->intermediate;
	  if(offshell) {
	    if(cc&&offshell->CC()) offshell=offshell->CC();
	    unsigned int o2(out2[dit->channelType]), o3(out3[dit->channelType]);
	    double sign = (o3 < o2) ? 1. : -1.;
	    // intermediate scalar
	    if(offshell->iSpin() == PDT::Spin0) {
	      ScalarWaveFunction inters = sca_[idiag].first->
		evaluate(scale, widthOption(), offshell, outVector_[v1], swave_);
	      unsigned int h1(s1),h2(s2);
	      if(o2 > o3) swap(h1, h2);
	      if(outgoing[o2]->id() < 0 &&  outgoing[o3]->id() > 0) {
		diag = -sign*sca_[idiag].second->
		  evaluate(scale,outspin_[o2].first[h1],
			   outspin_[o3].second[h2],inters);
	      }
	      else {
		diag = sign*sca_[idiag].second->
		  evaluate(scale, outspin_[o3].first [h2],
			   outspin_[o2].second[h1],inters);
	      }
	    }
	    // intermediate fermion
	    else if(offshell->iSpin() == PDT::Spin1Half) {
	      int iferm = (outgoing[o2]->iSpin() == PDT::Spin1Half) 
		? o2 : o3;
	      unsigned int h1(s1),h2(s2);
	      if(dit->channelType > iferm) swap(h1, h2);
	      sign = iferm < dit->channelType ? 1. : -1.;
	      if((outgoing[dit->channelType]->id() < 0 && outgoing[iferm]->id() > 0 ) ||
		 (outgoing[dit->channelType]->id()*offshell->id()>0)) {
		SpinorWaveFunction inters = fer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].first[h1], swave_);
		diag = -sign*fer_[idiag].second->
		  evaluate(scale,inters,outspin_[iferm].second[h2], outVector_[v1]);
	      }
	      else {
		SpinorBarWaveFunction inters = fer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].second[h1],swave_);
		diag =  sign*fer_[idiag].second->
		  evaluate(scale,outspin_[iferm].first [h2],inters, outVector_[v1]);
	      }
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = vec_[idiag].first->
		evaluate(scale, widthOption(), offshell, outVector_[v1], swave_);
	      unsigned int h1(s1),h2(s2);
	      if(o2 > o3) swap(h1,h2);
	      if(outgoing[o2]->id() < 0 && outgoing[o3]->id() > 0) {
		diag =-sign*vec_[idiag].second->
		  evaluate(scale, outspin_[o2].first[h1],
			   outspin_[o3].second[h2], interv);
	      }
	      else {
		diag = sign*vec_[idiag].second->
		  evaluate(scale, outspin_[o3].first[h2],
			   outspin_[o2].second[h1], interv);
	      }
	    }
	    // intermediate RS fermion
	    else if(offshell->iSpin() == PDT::Spin3Half) {
	      int iferm = (outgoing[o2]->iSpin() == PDT::Spin1Half) 
		? o2 : o3;
	      unsigned int h1(s1),h2(s2);
	      if(dit->channelType > iferm) swap(h1, h2);
	      sign = iferm < dit->channelType ? 1. : -1.;
	      if((outgoing[dit->channelType]->id() < 0 && outgoing[iferm]->id() > 0 ) ||
		 (outgoing[dit->channelType]->id()*offshell->id()>0)) {
		RSSpinorWaveFunction inters = RSfer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].first[h1], swave_);
		diag = -sign*RSfer_[idiag].second->
		  evaluate(scale,inters,outspin_[iferm].second[h2], outVector_[v1]);
	      }
	      else {
		RSSpinorBarWaveFunction inters = RSfer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].second[h1],swave_);
		diag =  sign*RSfer_[idiag].second->
		  evaluate(scale,outspin_[iferm].first [h2],inters, outVector_[v1]);
	      }
	    }
	    // unknown
	    else throw Exception()
		   << "Unknown intermediate in StoFFVDecayer::me2()" 
		   << Exception::runerror;
	  }
	  else {
	    unsigned int o2 = ivec > 0 ? 0 : 1;
	    unsigned int o3 = ivec < 2 ? 2 : 1;
	    if(outgoing[o2]->id() < 0 && outgoing[o3]->id() > 0) {
	      diag =-four_[idiag]->
	    	evaluate(scale, outspin_[o2].first[s1],
	    		 outspin_[o3].second[s2], outVector_[v1], swave_);
	    }
	    else {
	      diag = four_[idiag]->
	    	evaluate(scale, outspin_[o3].first[s2],
	    		 outspin_[o2].second[s1], outVector_[v1], swave_);
	    }
	  }
	  // matrix element for the different colour flows
	  if(ichan < 0) {
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	      largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		dit->largeNcColourFlow[iy].second * diag;
	    }
	  }
	  else {
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      if(dit->colourFlow[iy].first - 1 != colourFlow()) continue;
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	      if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
	      largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		dit->largeNcColourFlow[iy].second * diag;
	    }
	  }
	  ++idiag;
	} //end of diagrams
	// now add the flows to the me2 with appropriate colour factors
	for(unsigned int ix = 0; ix < ncf; ++ix) {
	  if ( ivec == 0 ) { 
	    (*mes[ix])(0, v1, s1, s2) = flows[ix];
	    (*mel[ix])(0, v1, s1, s2) = largeflows[ix];
	  }
	  else if( ivec == 1 ) {
	    (*mes[ix])(0, s1, v1, s2) = flows[ix];
	    (*mel[ix])(0, s1, v1, s2) = largeflows[ix];
	  }
	  else if( ivec == 2 ) {
	    (*mes[ix])(0, s1, s2, v1) = flows[ix];
	    (*mel[ix])(0, s1, s2, v1) = largeflows[ix];
	  }
	}
      }
    }
  }
  double me2(0.);
  if(ichan < 0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix == iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix = 0;ix < pflows.size(); ++ix) {
      if(ptotal <= pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal -= pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}
#line 1 "./VtoFFVDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VtoFFVDecayer class.
//

#include "VtoFFVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;

IBPtr VtoFFVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr VtoFFVDecayer::fullclone() const {
  return new_ptr(*this);
}

void VtoFFVDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << fer_ << vec_ << ten_;
}

void VtoFFVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> fer_ >> vec_ >> ten_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VtoFFVDecayer,GeneralThreeBodyDecayer>
describeHerwigVtoFFVDecayer("Herwig::VtoFFVDecayer", "Herwig.so");

void VtoFFVDecayer::Init() {

  static ClassDocumentation<VtoFFVDecayer> documentation
    ("The VtoFFVDecayer class implements the general three-body "
     "decay of a vector to a two fermions and a vector.");

}

WidthCalculatorBasePtr VtoFFVDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<VtoFFVDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),
		  outgoing()[2]->mass(),relativeError()));
}

void VtoFFVDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  fer_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if( offshell->CC() ) offshell = offshell->CC();
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractVVSVertexPtr vert1 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(current.vertices.first);
      AbstractFFSVertexPtr vert2 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a fermion diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      fer_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractVVVVertexPtr vert1 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractVVTVertexPtr vert1 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(current.vertices.first);
      AbstractFFTVertexPtr vert2 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in VtoFFVDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
  }
}
void VtoFFVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::
    constructSpinInfo(inVector_,const_ptr_cast<tPPtr>(&part),
		      Helicity::incoming,true,false);
  for(unsigned int ix=0;ix<decay.size();++ix) {
    if(decay[ix]->dataPtr()->iSpin()==PDT::Spin1) {
      VectorWaveFunction::constructSpinInfo(outVector_,decay[ix],
					    Helicity::outgoing,true,false);
    }
    else {
      SpinorWaveFunction::
	constructSpinInfo(outspin_[ix].first,decay[ix],Helicity::outgoing,true);
    }
  }
}

double VtoFFVDecayer::me2(const int ichan,const Particle & part,
			  const tPDVector & outgoing,
			  const vector<Lorentz5Momentum> & momenta,
			  MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != part.id();
  // special handling or first/last call
  if(meopt==Initialize) {
    VectorWaveFunction::
      calculateWaveFunctions(inVector_,rho_,const_ptr_cast<tPPtr>(&part),
			     Helicity::incoming,false);
    // fix rho if no correlations
    fixRho(rho_);
  }
  unsigned int ivec(0);
  bool massless = false;
  for(unsigned int ix = 0; ix < outgoing.size();++ix) {
    if(outgoing[ix]->iSpin() == PDT::Spin1) {
      massless = outgoing[ix]->id()==ParticleID::g || outgoing[ix]->id()==ParticleID::gamma;
      ivec = ix;
      VectorWaveFunction::
	calculateWaveFunctions(outVector_, momenta[ix], outgoing[ix], Helicity::outgoing, massless );
    }
    else {
      SpinorWaveFunction::
	calculateWaveFunctions(outspin_[ix].first,momenta[ix], outgoing[ix],Helicity::outgoing);
      outspin_[ix].second.resize(2);
      // Need a ubar and a v spinor
      if(outspin_[ix].first[0].wave().Type() == SpinorType::u) {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].first[iy].conjugate();
	}
      }
      else {
	for(unsigned int iy = 0; iy < 2; ++iy) {
	  outspin_[ix].second[iy] = outspin_[ix].first[iy].bar();
	  outspin_[ix].second[iy].conjugate();
	}
      }
    }
  }
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  Energy2 scale(sqr(part.mass()));
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.));
  // setup the DecayMatrixElement
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1,
					      ivec == 0 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 1 ? PDT::Spin1 : PDT::Spin1Half,
					      ivec == 2 ? PDT::Spin1 : PDT::Spin1Half)));
  
  //the channel possiblities
  static const unsigned int out2[3] = {1,0,0}, out3[3] = {2,2,1};
  for(unsigned int vi = 0; vi < 3; ++vi) {
    for(unsigned int s1 = 0; s1 < 2; ++s1) {
      for(unsigned int s2 = 0; s2 < 2; ++s2) {
	for(unsigned int v1 = 0; v1 < 3; ++v1) {
	  if ( massless && v1 == 1 ) continue; 
	  flows = vector<Complex>(ncf, Complex(0.));
	  largeflows = vector<Complex>(ncf, Complex(0.));
	  unsigned int idiag(0);
	  for(vector<TBDiagram>::const_iterator dit=getProcessInfo().begin();
	      dit!=getProcessInfo().end();++dit) {
	    // channels if selecting
	    if( ichan >=0 && diagramMap()[ichan] != idiag ) {
	      ++idiag;
	      continue;
	    }
	    tcPDPtr offshell = dit->intermediate;
	    if(cc&&offshell->CC()) offshell=offshell->CC();
	    Complex diag;
	    unsigned int o2(out2[dit->channelType]), o3(out3[dit->channelType]);
	    double sign = (o3 < o2) ? 1. : -1.;
	    // intermediate scalar
	    if(offshell->iSpin() == PDT::Spin0) {
	      ScalarWaveFunction inters = sca_[idiag].first->
		evaluate(scale, widthOption(), offshell, outVector_[v1], 
			 inVector_[vi]);
	      unsigned int h1(s1),h2(s2);
	      if(o2 > o3) swap(h1, h2);
	      if(outgoing[o2]->id() < 0 &&  outgoing[o3]->id() > 0) {
		diag = -sign*sca_[idiag].second->
		  evaluate(scale,outspin_[o2].first[h1],
			   outspin_[o3].second[h2],inters);
	      }
	      else {
		diag = sign*sca_[idiag].second->
		  evaluate(scale, outspin_[o3].first [h2],
			   outspin_[o2].second[h1],inters);
	      }
	    }
	    // intermediate fermion
	    else if(offshell->iSpin() == PDT::Spin1Half) {
	      int iferm = (outgoing[o2]->iSpin() == PDT::Spin1Half) 
		? o2 : o3;
	      unsigned int h1(s1),h2(s2);
	      if(dit->channelType > iferm) swap(h1, h2);
	      sign = iferm < dit->channelType ? 1. : -1.;
	      if(outgoing[dit->channelType]->id() < 0 && outgoing[iferm]->id() > 0 ) {
		SpinorWaveFunction inters = fer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].first[h1], inVector_[vi]);
		diag = -sign*fer_[idiag].second->
		  evaluate(scale,inters,outspin_[iferm].second[h2], outVector_[v1]);
	      }
	      else {
		SpinorBarWaveFunction inters = fer_[idiag].first->
		  evaluate(scale,widthOption(),offshell,
			   outspin_[dit->channelType].second[h1],inVector_[vi]);
		diag =  sign*fer_[idiag].second->
		  evaluate(scale,outspin_[iferm].first [h2],inters, outVector_[v1]);
	      }
	    }
	    // intermediate vector
	    else if(offshell->iSpin() == PDT::Spin1) {
	      VectorWaveFunction interv = vec_[idiag].first->
		evaluate(scale, widthOption(), offshell, outVector_[v1], 
			 inVector_[vi]);
	      unsigned int h1(s1),h2(s2);
	      if(o2 > o3) swap(h1,h2);
	      if(outgoing[o2]->id() < 0 && outgoing[o3]->id() > 0) {
		diag =-sign*vec_[idiag].second->
		  evaluate(scale, outspin_[o2].first[h1],
			   outspin_[o3].second[h2], interv);
	      }
	      else {
		diag = sign*vec_[idiag].second->
		  evaluate(scale, outspin_[o3].first[h2],
			   outspin_[o2].second[h1], interv);
	      }
	    }
	    else if(offshell->iSpin() == PDT::Spin2) {
	      TensorWaveFunction intert = ten_[idiag].first->
		evaluate(scale, widthOption(), offshell, inVector_[vi], 
			 outVector_[v1]);
	      unsigned int h1(s1),h2(s2);
	      if(out2[dit->channelType]>out3[dit->channelType]) swap(h1,h2);
	      if(outgoing[out2[dit->channelType]]->id()<0&&
		 outgoing[out3[dit->channelType]]->id()>0) {
		diag =-sign*ten_[idiag].second->
		  evaluate(scale,
			   outspin_[out2[dit->channelType]].first [h1],
			   outspin_[out3[dit->channelType]].second[h2],intert);
	      }
	      else {
		diag = sign*ten_[idiag].second->
		  evaluate(scale,
			   outspin_[out3[dit->channelType]].first [h2],
			   outspin_[out2[dit->channelType]].second[h1],intert);
	      }
	    }
	    // unknown
	    else throw Exception()
	      << "Unknown intermediate in VtoFFVDecayer::me2()" 
	      << Exception::runerror;
	    
	    // matrix element for the different colour flows
	    if(ichan < 0) {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }
	    else {
	      for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1 != colourFlow()) continue;
		flows[dit->colourFlow[iy].first - 1] += 
		  dit->colourFlow[iy].second * diag;
	      }
	      for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
		if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
		largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		  dit->largeNcColourFlow[iy].second * diag;
	      }
	    }	    
	    ++idiag;
	  } //end of diagrams
	  // now add the flows to the me2 with appropriate colour factors
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    if (ivec == 0) {
	      (*mes[ix])(vi, v1, s1, s2) = flows[ix];
	      (*mel[ix])(vi, v1, s1, s2) = largeflows[ix];
	    }
	    else if(ivec == 1) {
	      (*mes[ix])(vi, s1, v1, s2) = flows[ix];
	      (*mel[ix])(vi, s1, v1, s2) = largeflows[ix];
	      
	    }
	    else if(ivec == 2) { 
	      (*mes[ix])(vi, s1, s2, v1) = flows[ix];
	      (*mel[ix])(vi, s1, s2, v1) = largeflows[ix];
	    }
	  }

	}
      }
    }
  }
  double me2(0.);
  if(ichan < 0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix == iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix = 0;ix < pflows.size(); ++ix) {
      if(ptotal <= pflows[ix]) {
	colourFlow(ix);
	ME(mes[ix]);
	break;
      }
      ptotal -= pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}
#line 1 "./FtoFVVDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FtoFVVDecayer class.
//

#include "FtoFVVDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include <numeric>

using namespace Herwig;

IBPtr FtoFVVDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FtoFVVDecayer::fullclone() const {
  return new_ptr(*this);
}

void FtoFVVDecayer::persistentOutput(PersistentOStream & os) const {
  os << sca_ << fer_ << vec_ << ten_;
}

void FtoFVVDecayer::persistentInput(PersistentIStream & is, int) {
  is >> sca_ >> fer_ >> vec_ >> ten_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FtoFVVDecayer,GeneralThreeBodyDecayer>
describeHerwigFtoFVVDecayer("Herwig::FtoFVVDecayer", "Herwig.so");

void FtoFVVDecayer::Init() {

  static ClassDocumentation<FtoFVVDecayer> documentation
    ("The FtoFVVDecayer class implements the general decay of a fermion to "
     "a fermion and a pair of vectors.");

}

void FtoFVVDecayer::setupDiagrams(bool kinCheck) {
  GeneralThreeBodyDecayer::setupDiagrams(kinCheck);
  if(outgoing().empty()) return;
  unsigned int ndiags = getProcessInfo().size();
  sca_.resize(ndiags);
  fer_.resize(ndiags);
  vec_.resize(ndiags);
  ten_.resize(ndiags);
  for(unsigned int ix = 0;ix < ndiags; ++ix) {
    TBDiagram current = getProcessInfo()[ix];
    tcPDPtr offshell = current.intermediate;
    if(offshell->iSpin() == PDT::Spin0) {
      AbstractFFSVertexPtr vert1 = dynamic_ptr_cast<AbstractFFSVertexPtr>
	(current.vertices.first);
      AbstractVVSVertexPtr vert2 = dynamic_ptr_cast<AbstractVVSVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFVVDecayer::setupDiagrams()"
	<< Exception::runerror;
      sca_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1Half) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractFFVVertexPtr vert2 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a scalar diagram in FtoFVVDecayer::setupDiagrams()"
	<< Exception::runerror;
      fer_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin1) {
      AbstractFFVVertexPtr vert1 = dynamic_ptr_cast<AbstractFFVVertexPtr>
	(current.vertices.first);
      AbstractVVVVertexPtr vert2 = dynamic_ptr_cast<AbstractVVVVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a vector diagram in FtoFVVDecayer::setupDiagrams()"
	<< Exception::runerror;
      vec_[ix] = make_pair(vert1, vert2);
    }
    else if(offshell->iSpin() == PDT::Spin2) {
      AbstractFFTVertexPtr vert1 = dynamic_ptr_cast<AbstractFFTVertexPtr>
	(current.vertices.first);
      AbstractVVTVertexPtr vert2 = dynamic_ptr_cast<AbstractVVTVertexPtr>
	(current.vertices.second);
      if(!vert1||!vert2) throw Exception() 
	<< "Invalid vertices for a tensor diagram in FtoFVVDecayer::setupDiagrams()"
	<< Exception::runerror;
      ten_[ix] = make_pair(vert1, vert2);
    }
  }
}
void FtoFVVDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  bool ferm( part.id() > 0 );
  // for the decaying particle
  if(ferm)
    SpinorWaveFunction::constructSpinInfo(fwave_,
					  const_ptr_cast<tPPtr>(&part),
					  Helicity::incoming,true);
  else
    SpinorBarWaveFunction::constructSpinInfo(fbwave_,
					     const_ptr_cast<tPPtr>(&part),
					     Helicity::incoming,true);
  int ivec(-1);
  // outgoing particles
  for(int ix = 0; ix < 3; ++ix) {
    tPPtr p = decay[ix];
    if( p->dataPtr()->iSpin() == PDT::Spin1Half ) {
      if( ferm ) {
	SpinorBarWaveFunction::
	  constructSpinInfo(fbwave_, p, Helicity::outgoing,true);
      }
      else {
	SpinorWaveFunction::
	  constructSpinInfo(fwave_ , p, Helicity::outgoing,true);
      }
    }
    else if( p->dataPtr()->iSpin() == PDT::Spin1 ) {
      if( ivec < 0 ) {
	ivec = ix;
	VectorWaveFunction::
	  constructSpinInfo(vwave_.first , p, Helicity::outgoing, true, false);
      }
      else {
	VectorWaveFunction::
	  constructSpinInfo(vwave_.second, p, Helicity::outgoing, true, false);
      }
    }
  }
}

double  FtoFVVDecayer::me2(const int ichan,const Particle & part,
			   const tPDVector & outgoing,
			   const vector<Lorentz5Momentum> & momenta,
			   MEOption meopt) const {
  // particle or CC of particle
  bool cc = (*getProcessInfo().begin()).incoming != part.id();
  // special handling or first/last call
  //Set up wave-functions
  bool ferm( part.id() > 0 );
  if(meopt==Initialize) {
    if( ferm ) {
      SpinorWaveFunction::
	calculateWaveFunctions(fwave_,rho_,const_ptr_cast<tPPtr>(&part),
			       Helicity::incoming);
      if( fwave_[0].wave().Type() != SpinorType::u )
	fwave_[0].conjugate();
      if( fwave_[1].wave().Type() != SpinorType::u )
	fwave_[1].conjugate();
    }
    else {
      SpinorBarWaveFunction::
	calculateWaveFunctions(fbwave_, rho_, const_ptr_cast<tPPtr>(&part),
			       Helicity::incoming);
      if( fbwave_[0].wave().Type() != SpinorType::v )
	fbwave_[0].conjugate();
      if( fbwave_[1].wave().Type() != SpinorType::v )
	fbwave_[1].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  // outgoing, keep track of fermion and first occurrence of vector positions
  int isp(-1), ivec(-1);
  // outgoing particles
  pair<bool,bool> mass = make_pair(false,false);
  for(int ix = 0; ix < 3; ++ix) {
    if( outgoing[ix]->iSpin() == PDT::Spin1Half ) {
      isp = ix;
      if( ferm ) {
	SpinorBarWaveFunction::
	  calculateWaveFunctions(fbwave_, momenta[ix],outgoing[ix], Helicity::outgoing);
	if( fbwave_[0].wave().Type() != SpinorType::u )
	  fbwave_[0].conjugate();
	if( fbwave_[1].wave().Type() != SpinorType::u )
	  fbwave_[1].conjugate();
      }
      else {
	SpinorWaveFunction::
	  calculateWaveFunctions(fwave_, momenta[ix],outgoing[ix], Helicity::outgoing);
	if( fwave_[0].wave().Type() != SpinorType::v )
	  fwave_[0].conjugate();
	if( fwave_[1].wave().Type() != SpinorType::v )
	  fwave_[1].conjugate();
      }
    }
    else if( outgoing[ix]->iSpin() == PDT::Spin1 ) {
      bool massless = outgoing[ix]->id() == ParticleID::gamma || outgoing[ix]->id() == ParticleID::g;
      if( ivec < 0 ) {
	ivec = ix;
	VectorWaveFunction::
	  calculateWaveFunctions(vwave_.first , momenta[ix],outgoing[ix], Helicity::outgoing, massless);
	mass.first = massless;
      }
      else {
	VectorWaveFunction::
	  calculateWaveFunctions(vwave_.second, momenta[ix],outgoing[ix], Helicity::outgoing, massless);
	mass.second = massless;
      }
    }
  }
  assert(isp >= 0 && ivec >= 0);
  Energy2 scale(sqr(part.mass()));
  const vector<vector<double> > cfactors(getColourFactors());
  const vector<vector<double> > nfactors(getLargeNcColourFactors());
  const size_t ncf(numberOfFlows());
  vector<Complex> flows(ncf, Complex(0.)), largeflows(ncf, Complex(0.)); 
  vector<GeneralDecayMEPtr> 
    mes(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, 
					      (isp == 0) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 1) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 2) ? PDT::Spin1Half : PDT::Spin1)));
  vector<GeneralDecayMEPtr> 
    mel(ncf,new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half, 
					      (isp == 0) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 1) ? PDT::Spin1Half : PDT::Spin1,
					      (isp == 2) ? PDT::Spin1Half : PDT::Spin1)));
  //Helicity calculation
  for( unsigned int if1 = 0; if1 < 2; ++if1 ) {
    for( unsigned int if2 = 0; if2 < 2; ++if2 ) {
      for( unsigned int iv1 = 0; iv1 < 3; ++iv1 ) {
	if ( mass.first && iv1 == 1 ) continue;
	for( unsigned int iv2 = 0; iv2 < 3; ++iv2 ) {
	  if ( mass.second && iv2 == 1 ) continue;
	  flows = vector<Complex>(ncf, Complex(0.));
	largeflows = vector<Complex>(ncf, Complex(0.));
	unsigned int idiag(0);
	for(vector<TBDiagram>::const_iterator dit = getProcessInfo().begin();
	    dit != getProcessInfo().end(); ++dit) {
	  // If we are selecting a particular channel
	  if( ichan >= 0 && diagramMap()[ichan] != idiag ) {
	    ++idiag;
	    continue;
	  }
	  tcPDPtr offshell = (*dit).intermediate;
	  if(cc&&offshell->CC()) offshell=offshell->CC();
	  Complex diag;
	  if( offshell->iSpin() == PDT::Spin1Half ) {
	    // Make sure we connect the correct particles 
	    VectorWaveFunction vw1, vw2;
	    if( (*dit).channelType == TBDiagram::channel23 ) {
	      vw1 = vwave_.first[iv1];
	      vw2 = vwave_.second[iv2];
	    }
	    else if( (*dit).channelType == TBDiagram::channel12 ) {
	      vw1 = vwave_.second[iv2];
	      vw2 = vwave_.first[iv1];
	    }
	    else {
	      if( ivec < isp ) {
		vw1 = vwave_.second[iv2];
		vw2 = vwave_.first[iv1];
	      }
	      else {
		vw1 = vwave_.first[iv1];
		vw2 = vwave_.second[iv2];
	      }
	    }
	    if( ferm ) {
	      SpinorWaveFunction inters = 
		fer_[idiag].first->evaluate(scale, widthOption(), offshell,
					    fwave_[if1], vw1);
	      diag = fer_[idiag].second->evaluate(scale, inters, fbwave_[if2],
						  vw2);
	    }
	    else {
	      SpinorBarWaveFunction inters = 
		fer_[idiag].first->evaluate(scale, widthOption(), offshell,
					    fbwave_[if2], vw1);
	      diag = fer_[idiag].second->evaluate(scale, fwave_[if1], inters, 
						  vw2);
	    }
	  }
	  else if( offshell->iSpin() == PDT::Spin0 ) {
	    ScalarWaveFunction inters = 
	      sca_[idiag].first->evaluate(scale, widthOption(), offshell, 
					  fwave_[if1], fbwave_[if2]);
	    diag = sca_[idiag].second->evaluate(scale, vwave_.first[iv1],
						vwave_.second[iv2], inters);
	  }
	  else if( offshell->iSpin() == PDT::Spin1 ) {
	    VectorWaveFunction interv = 
	      vec_[idiag].first->evaluate(scale, widthOption(), offshell, 
					  fwave_[if1], fbwave_[if2]);
	    diag = vec_[idiag].second->evaluate(scale, vwave_.first[iv1],
						vwave_.second[iv2], interv);
	  } 
	  else if( offshell->iSpin() == PDT::Spin2 ) {
	    TensorWaveFunction intert = 
	      ten_[idiag].first->evaluate(scale, widthOption(), offshell, 
					  fwave_[if1], fbwave_[if2]);
	    diag = ten_[idiag].second->evaluate(scale, vwave_.first[iv1],
						vwave_.second[iv2], intert);
	  }
	  else 
	    throw Exception()
	      << "Unknown intermediate in FtoFVVDecayer::me2()" 
	      << Exception::runerror;
	  //NO sign
	  if( !ferm ) diag *= -1;

	  if(ichan<0) {
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	      largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		dit->largeNcColourFlow[iy].second * diag;
	    }
	  }
	  else {
	    for(unsigned iy = 0; iy < dit->colourFlow.size(); ++iy) {
	      if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
	      flows[dit->colourFlow[iy].first - 1] += 
		dit->colourFlow[iy].second * diag;
	    }
	    for(unsigned iy = 0; iy < dit->largeNcColourFlow.size(); ++iy) {
	      if(dit->colourFlow[iy].first - 1!=colourFlow()) continue;
	      largeflows[dit->largeNcColourFlow[iy].first - 1] += 
		dit->largeNcColourFlow[iy].second * diag;
	    }
	  }
	  ++idiag;
	}// end diagram loop

	// now add the flows to the me2 
	unsigned int h1(if1), h2(if2);
	if( !ferm ) swap(h1,h2);
	  for(unsigned int ix = 0; ix < ncf; ++ix) {
	    if(isp == 0) {
	      (*mes[ix])(h1, h2, iv1, iv2) = flows[ix];
	      (*mel[ix])(h1, h2, iv1, iv2) = largeflows[ix];
	    }
	    else if(isp == 1) { 
	      (*mes[ix])(h1, iv1, h2, iv2) = flows[ix];
	      (*mel[ix])(h1, iv1, h2, iv2) = largeflows[ix];
	    }
	    else if(isp == 2) { 
	      (*mes[ix])(h1, iv1, iv2, h2) = flows[ix];
	      (*mel[ix])(h1, iv1, h2, iv2) = largeflows[ix];
	    }
	  }

	}//v2
      }//v1
    }//f2
  }//f1
  
  //Finally, work out me2. This depends on whether we are selecting channels
  //or not
  double me2(0.);
  if(ichan<0) {
    vector<double> pflows(ncf,0.);
    for(unsigned int ix = 0; ix < ncf; ++ix) {
      for(unsigned int iy = 0; iy < ncf; ++ iy) {
	double con = cfactors[ix][iy]*(mes[ix]->contract(*mes[iy],rho_)).real();
	me2 += con;
	if(ix==iy) {
	  con = nfactors[ix][iy]*(mel[ix]->contract(*mel[iy],rho_)).real();
	  pflows[ix] += con;
	}
      }
    }
    double ptotal(std::accumulate(pflows.begin(),pflows.end(),0.));
    ptotal *= UseRandom::rnd();
    for(unsigned int ix=0;ix<pflows.size();++ix) {
      if(ptotal<=pflows[ix]) {
        colourFlow(ix);
        ME(mes[ix]);
        break;
      }
      ptotal-=pflows[ix];
    }
  }
  else {
    unsigned int iflow = colourFlow();
    me2 = nfactors[iflow][iflow]*(mel[iflow]->contract(*mel[iflow],rho_)).real();
  }
  // return the matrix element squared
  return me2;
}

WidthCalculatorBasePtr FtoFVVDecayer::
threeBodyMEIntegrator(const DecayMode & ) const {
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  vector<double> inpow,inweights;
  constructIntegratorChannels(intype,inmass,inwidth,inpow,inweights);
  return new_ptr(ThreeBodyAllOnCalculator<FtoFVVDecayer>
		 (inweights,intype,inmass,inwidth,inpow,*this,0,
		  outgoing()[0]->mass(),outgoing()[1]->mass(),
		  outgoing()[2]->mass(),relativeError()));
}
#line 1 "./FFVCurrentDecayer.cc"
// -*- C++ -*-
//
// FFVCurrentDecayer.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2019 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FFVCurrentDecayer class.
//

#include "FFVCurrentDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

using ThePEG::Helicity::VectorWaveFunction;
using ThePEG::Helicity::SpinorWaveFunction;
using ThePEG::Helicity::SpinorBarWaveFunction;
using ThePEG::Helicity::Direction;
using ThePEG::Helicity::incoming;
using ThePEG::Helicity::outgoing;

IBPtr FFVCurrentDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr FFVCurrentDecayer::fullclone() const {
  return new_ptr(*this);
}

void FFVCurrentDecayer::doinit() {
  FFVPtr_ = dynamic_ptr_cast<FFVVertexPtr>(vertex());
  GeneralCurrentDecayer::doinit();
}


void FFVCurrentDecayer::rebind(const TranslationMap & trans)
  {
  FFVPtr_ = trans.translate(FFVPtr_);
  GeneralCurrentDecayer::rebind(trans);
}

IVector FFVCurrentDecayer::getReferences() {
  IVector ret = GeneralCurrentDecayer::getReferences();
  ret.push_back(FFVPtr_);
  return ret;
}

void FFVCurrentDecayer::persistentOutput(PersistentOStream & os) const {
  os << FFVPtr_;
}

void FFVCurrentDecayer::persistentInput(PersistentIStream & is, int) {
  is >> FFVPtr_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<FFVCurrentDecayer,GeneralCurrentDecayer>
describeHerwigFFVCurrentDecayer("Herwig::FFVCurrentDecayer", "Herwig.so");

void FFVCurrentDecayer::Init() {

  static ClassDocumentation<FFVCurrentDecayer> documentation
    ("There is no documentation for the FFVCurrentDecayer class");

}
void FFVCurrentDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // setup spin info when needed
  // fermion types
  int itype[2];
  if(part.dataPtr()->CC())    itype[0] = part.id() > 0 ? 0 : 1;
  else                          itype[0] = 2;
  if(decay[0]->dataPtr()->CC()) itype[1] = decay[0]->id() > 0 ? 0 : 1;
  else                          itype[1] = 2;
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  // for the decaying particle
  if(ferm) {
    SpinorWaveFunction::
      constructSpinInfo(wave_,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorBarWaveFunction::constructSpinInfo(wavebar_,decay[0],outgoing,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(wavebar_,const_ptr_cast<tPPtr>(&part),incoming,true);
    SpinorWaveFunction::constructSpinInfo(wave_,decay[0],outgoing,true);
  }
  weakCurrent()->constructSpinInfo(ParticleVector(decay.begin()+1,decay.end()));
}

double FFVCurrentDecayer::me2(const int ichan, const Particle & part,
			      const tPDVector & outgoing,
			      const vector<Lorentz5Momentum> & momenta,
			      MEOption meopt) const {
  // fermion types
  int itype[2];
  if(part.dataPtr()->CC()) itype[0] = part.id() > 0 ? 0 : 1;
  else                     itype[0] = 2;
  if(outgoing[0]->CC())    itype[1] = outgoing[0]->id() > 0 ? 0 : 1;
  else                     itype[1] = 2;
  //Need to use different barred or unbarred spinors depending on 
  //whether particle is cc or not.
  bool ferm(itype[0] == 0 || itype[1] == 0 || (itype[0] == 2 && itype[1] == 2));
  if(meopt==Initialize) {
    // spinors and rho
    if(ferm) {
      SpinorWaveFunction   ::calculateWaveFunctions(wave_,rho_,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
      if(wave_[0].wave().Type() != SpinorType::u)
  	for(unsigned int ix = 0; ix < 2; ++ix) wave_   [ix].conjugate();
    }
    else {
      SpinorBarWaveFunction::calculateWaveFunctions(wavebar_,rho_,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
      if(wavebar_[0].wave().Type() != SpinorType::v)
  	for(unsigned int ix = 0; ix < 2; ++ix) wavebar_[ix].conjugate();
    }
    // fix rho if no correlations
    fixRho(rho_);
  }
  Energy2 scale(sqr(part.mass()));
  // spinors for the outgoing ferimon
  if(ferm) {
    wavebar_.resize(2);
    SpinorBarWaveFunction wbar(momenta[0],outgoing[0],Helicity::outgoing);
    for(unsigned int ihel=0;ihel<2;++ihel) {
      wbar.reset(ihel);
      wavebar_[ihel] = wbar;
    }
  }
  else {
    wave_   .resize(2);
    SpinorWaveFunction    w   (momenta[0],outgoing[0],Helicity::outgoing);
    for(unsigned int ihel=0;ihel<2;++ihel) {
      w.reset(ihel);
      wave_   [ihel] = w;
    }
  }
  // calculate the hadron current
  Energy q;
  vector<LorentzPolarizationVectorE> 
    hadron(weakCurrent()->current(tcPDPtr(),FlavourInfo(),mode(),ichan,q,tPDVector(outgoing.begin()+1,outgoing.end()),
				  vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt));
  // prefactor
  double pre = sqr(pow(part.mass()/q,int(outgoing.size()-3)));
  // work out the mapping for the hadron vector
  vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
  vector<PDT::Spin> ispin(outgoing.size());
  int itemp(1);
  unsigned int hhel,ix(outgoing.size());
  do {
    --ix;
    ispin[ix]=outgoing[ix]->iSpin();
    itemp*=ispin[ix];
    constants[ix]=itemp;
  }
  while(ix>0);
  constants[outgoing.size()]=1;
  constants[0]=constants[1];
  // compute the matrix element
  GeneralDecayMEPtr newME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,ispin)));
  VectorWaveFunction vWave;
  tcPDPtr vec= part.dataPtr()->iCharge()-outgoing[0]->iCharge() > 0
    ? getParticleData(ParticleID::Wplus) : getParticleData(ParticleID::Wminus);
  Lorentz5Momentum vmom=part.momentum()-momenta[0];
  vmom.rescaleMass();
  for(hhel=0;hhel<hadron.size();++hhel) {
    // map the index for the hadrons to a helicity state
    for(ix=outgoing.size();ix>1;--ix) ihel[ix]=(hhel%constants[ix-1])/constants[ix];
    vWave=VectorWaveFunction(vmom,vec,hadron[hhel]*UnitRemoval::InvE,Helicity::outgoing);
    for(unsigned int if1 = 0; if1 < 2; ++if1) {
      for(unsigned int if2 = 0; if2 < 2; ++if2) {
  	ihel[0]=if1;
  	ihel[1]=if2;
  	if(!ferm) swap(ihel[0],ihel[1]);
  	(*newME)(ihel) = FFVPtr_->evaluate(scale,wave_[if1],wavebar_[if2],vWave);
      }
    }
  }
  // store the matrix element
  ME(newME);
  // multiply by the CKM element
  int iq,ia;
  weakCurrent()->decayModeInfo(mode(),iq,ia);
  double ckm(1.);
  if(iq<=6) {
    if(iq%2==0) ckm = SM().CKM(iq/2-1,(abs(ia)-1)/2);
    else        ckm = SM().CKM(abs(ia)/2-1,(iq-1)/2);
  }
  pre /= 0.125*sqr(FFVPtr_->weakCoupling(scale));
  double output(0.5*pre*ckm*(ME()->contract(rho_)).real()*
  		sqr(SM().fermiConstant()*UnitRemoval::E2));
  return output;
}
 
Energy FFVCurrentDecayer::partialWidth(tPDPtr part, tPDPtr outa,
				       vector<tPDPtr> currout) {
  vector<long> id;
  id.push_back(part->id());
  id.push_back(outa->id());
  for(unsigned int ix=0;ix<currout.size();++ix) id.push_back(currout[ix]->id());
  bool cc;
  int mode=modeNumber(cc,id);
  imode(mode);
  // return initializePhaseSpaceMode(mode,true,true);  
  assert(false);
}
