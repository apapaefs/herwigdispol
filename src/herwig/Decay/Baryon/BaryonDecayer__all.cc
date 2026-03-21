#line 1 "./Baryon1MesonDecayerBase.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Baryon1MesonDecayerBase class.
//

#include "Baryon1MesonDecayerBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/DecayMode.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/Decay/TwoBodyDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractNoPIOClass<Baryon1MesonDecayerBase,DecayIntegrator>
describeHerwigBaryon1MesonDecayerBase("Herwig::Baryon1MesonDecayerBase", "HwBaryonDecay.so");

void Baryon1MesonDecayerBase::Init() {

  static ClassDocumentation<Baryon1MesonDecayerBase> documentation
    ("The Baryon1MesonDecayerBase class is the base class for"
     " the decays of the baryons to a baryon and a pseudoscalar or vector meson.");

}

void Baryon1MesonDecayerBase::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()>0) {
    // incoming particle
    if(part.dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorWaveFunction::
	constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else if(part.dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorWaveFunction::
	constructSpinInfo(_inThreeHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else
      assert(false);
    // outgoing fermion
    if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
						 decay[0],outgoing,true);
    }
    else
      assert(false);
  }
  else {
    // incoming particle
    if(part.dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorBarWaveFunction::
	constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else if(part.dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorBarWaveFunction::
	constructSpinInfo(_inThreeHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
    }
    else
      assert(false);
    // outgoing fermion
    if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
    else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
					      decay[0],outgoing,true);
    }
    else
      assert(false);
  }
  // outgoing meson
  if(decay[1]->dataPtr()->iSpin()==PDT::Spin0) {
    ScalarWaveFunction::constructSpinInfo(decay[1],outgoing,true);
  }
  else if(decay[1]->dataPtr()->iSpin()==PDT::Spin1) {
    VectorWaveFunction::constructSpinInfo(_inVec,decay[1],outgoing,true,
					  decay[1]->id()==ParticleID::gamma);
  }
  else
    assert(false);
}

// return the matrix element squared for a given mode and phase-space channel
// (inherited from DecayIntegrator and implemented here)
double Baryon1MesonDecayerBase::me2(const int ichan,const Particle & part,
				    const tPDVector & outgoing,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  double me(0.);
  // decide which matrix element we are doing
  // incoming spin-1/2 particle
  if(part.dataPtr()->iSpin()==2) {
    // decay to spin-1/2 particle
    if(outgoing[0]->iSpin()==2) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=halfHalfScalar(ichan,part,outgoing,momenta,meopt);
      // vector meson
      else if(outgoing[1]->iSpin()==3)
   	me=halfHalfVector(ichan,part,outgoing,momenta,meopt);
      else
   	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				      << "Baryon1MesonDecayerBase::me2()" 
				      << Exception::abortnow;
    }
    // decay to spin-3/2 particle
    else if(outgoing[0]->iSpin()==4) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=halfThreeHalfScalar(ichan,part,outgoing,momenta,meopt);
      // vector meson
      else if(outgoing[1]->iSpin()==3)
   	me=halfThreeHalfVector(ichan,part,outgoing,momenta,meopt);
      else
   	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				      << "Baryon1MesonDecayerBase::me2()" 
				      << Exception::abortnow;
    }
    // unknown
    else
      throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
				    << "Baryon1MesonDecayerBase::me2()" 
				    << Exception::abortnow;
  }
  // incoming spin-3/2 particle
  else if(part.dataPtr()->iSpin()==4) {
    // decay to spin-1/2 particle
    if(outgoing[0]->iSpin()==2) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=threeHalfHalfScalar(ichan,part,outgoing,momenta,meopt);
      // vector meson
      else if(outgoing[1]->iSpin()==3)
   	me=threeHalfHalfVector(ichan,part,outgoing,momenta,meopt);
      else
   	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				      << "Baryon1MesonDecayerBase::me2()" 
				      << Exception::abortnow;
    }
    // decay to spin-3/2 particle
    else if(outgoing[0]->iSpin()==4) {
      // scalar meson
      if(outgoing[1]->iSpin()==1)
   	me=threeHalfThreeHalfScalar(ichan,part,outgoing,momenta,meopt);
      else
   	throw DecayIntegratorError() << "Unknown outgoing meson spin in "
				      << "Baryon1MesonDecayerBase::me2()" 
				      << Exception::abortnow;
    }
    // unknown
    else
      throw DecayIntegratorError() << "Unknown outgoing baryon spin in "
				    << "Baryon1MesonDecayerBase::me2()" 
				    << Exception::abortnow;
  }
  // unknown
  else
    throw DecayIntegratorError() << "Unknown incoming spin in "
				  << "Baryon1MesonDecayerBase::me2()" 
				  << Exception::abortnow;
  return me;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a pseudoscalar meson
double Baryon1MesonDecayerBase::
halfHalfScalar(const int,const Particle & part, const tPDVector & outgoing,
	       const vector<Lorentz5Momentum> & momenta,
	       MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A,B;
  halfHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),A,B);
  Complex left,right;
  // coupling for an incoming particle
  if(part.id()>0) {
    left  = (A-B);
    right = (A+B);
  }
  // coupling for an incoming antiparticle
  else {
    left  = conj(A+B);
    right = conj(A-B);
  }
  // calculate the matrix element
  vector<unsigned int> ispin(3,0);
  unsigned int ix,iy;
  // Complex output(0.);
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      if(outgoing[0]->id()>0){ispin[0]=iy;ispin[1]=ix;}
      else{ispin[0]=ix;ispin[1]=iy;}
      (*ME())(ispin)=Complex(_inHalf[iy].generalScalar(_inHalfBar[ix],left,right)/part.mass());
    }
  }
  // test of the matrix element
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // Complex h1(2.*Qp*A/part.mass()),h2(-2.*Qm*B/part.mass());
  // cout << "testing 1/2->1/2 0 " 
  //      << 0.5*output << "   " 
  //      << 0.25*(h1*conj(h1)+h2*conj(h2)) << "   " 
  //      << 0.5*(h1*conj(h1)+h2*conj(h2))/output << endl;
  // cout << "testing alpha " << 
  //   (norm(0.5*(h1+h2))-norm(0.5*(h1-h2)))/
  //   (norm(0.5*(h1+h2))+norm(0.5*(h1-h2))) << "\n";
  // store the matrix element
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-1/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfHalfVector(const int,const Particle & part, const tPDVector & outgoing,
	       const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin1Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=outgoing[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  _inVec.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    _inVec[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A1,A2,B1,B2;
  halfHalfVectorCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			 A1,A2,B1,B2);
  Complex lS,rS,lV,rV;
  complex<Energy> scalar;
  // couplings for an incoming particle
  if(part.id()>0) {
    lS = (A2-B2);
    rS = (A2+B2);
    lV = (A1-B1);
    rV = (A1+B1);
  }
  else {
    lS = -conj(A2+B2);
    rS = -conj(A2-B2);
    lV =  conj(A1-B1);
    rV =  conj(A1+B1);
  }
  // calculate the matrix element
  // decide which type of mode to do
  Energy msum(part.mass()+momenta[0].mass());
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  Complex prod;
  // Complex output(0.);
  unsigned int ix,iy;
  for(ix=0;ix<2;++ix) {
    for(iy=0;iy<2;++iy) {
      // scalar like piece
      scalar = _inHalf[iy].generalScalar(_inHalfBar[ix],lS,rS);
      // vector like piece
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      if(outgoing[0]->id()>0) {
	ispin[0] = iy;
	ispin[1] = ix;
      }
      else {
	ispin[0] = ix;
	ispin[1] = iy;
      }
      for(ispin[2]=0;ispin[2]<3;++ispin[2]) {
	ispin[2]=ispin[2];
	prod=_inVec[ispin[2]].dot(part.momentum())/msum;
	(*ME())(ispin)=(svec.dot(_inVec[ispin[2]])+prod*scalar)/part.mass();
 	// output += norm((*ME())(ispin));
      }
    }
  }
  // test of the matrix element
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // double r2(sqrt(2.));
  // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
  // if(m3!=ZERO) {
  //   Complex h1(2.*r2*Qp*B1/part.mass()),h2(-2.*r2*Qm*A1/part.mass()),
  //     h3(2./m3*(Qp*(m1-m2)*B1-Qm*m1*B2*pcm/(m1+m2))/part.mass()),
  //     h4(2./m3*(Qm*(m1+m2)*A1+Qp*m1*A2*pcm/(m1+m2))/part.mass());
  //   generator()->log() << "testing 1/2->1/2 1 " 
  // 		       << 0.5*output << "   " 
  // 		       << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4)) << "   " 
  // 		       << 0.50*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+h4*conj(h4))/output 
  // 		       << "\n";
  //   generator()->log() << "alpha = " << 2.*(norm(h3)+norm(h4))/(norm(h1)+norm(h2))-1.
  // 		       << "\n";
  //   generator()->log() << "masses " << m1/GeV << " " << m2/GeV << " " << m3/GeV << "\n";
  // }
  // return the answer
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
halfThreeHalfScalar(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  // get the couplings
  Complex A,B,left,right;
  Energy msum(part.mass()+momenta[0].mass());
  halfThreeHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A,B);
  // incoming particle
  if(part.id()>0) {
    left=(A-B);
    right=(A+B);
  }
  // incoming anti-particle
  else {
    left=conj(A+B);
    right=conj(A-B);
  }
  vector<unsigned int> ispin(3,0);
  for(unsigned ixa=0;ixa<2;++ixa) {
    for(unsigned int iya=0;iya<4;++iya) {
      unsigned int ix(iya),iy(ixa);
      if(outgoing[0]->id()<0) swap(ix,iy);
      ispin[0]=ixa;
      ispin[1]=iya;
      complex<double> value = _inHalf[iy].generalScalar(_inHalfBar[ix],left,right)
	*UnitRemoval::E/part.mass()/msum;
      (*ME())(ispin) = value;
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
   // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
   // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
   // double r23(sqrt(2./3.));
   // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
   // complex<Energy> h1(-2.*r23*pcm*m1/m2*Qm*B/(m1+m2)),h2( 2.*r23*pcm*m1/m2*Qp*A/(m1+m2));
   // cout << "testing 1/2->3/2 0 " << part.id() << " "
   //      << output << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass()) << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass())/output << endl;
  // return the answer
  return output;
}

// matrix element for the decay of a spin-1/2 fermion to a spin-3/2 fermion and
// a vector meson
double Baryon1MesonDecayerBase::
halfThreeHalfVector(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin1Half,PDT::Spin3Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=outgoing[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    // matrix element
  }
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  // wavefunctions for outgoing fermion
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  ME()->zero();
  _inVec.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    _inVec[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3;
  halfThreeHalfVectorCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(part.mass()+momenta[0].mass());
  Complex lS,rS,lV,rV,left,right;
  // incoming particle
  if(part.id()>0) {
    lS=(A3-B3);rS=(A3+B3);
    lV=(A2-B2);rV=(A2+B2);
    left=(A1-B1);right=(A1+B1);
  }
  // incoming anti-particle
  else {
    lS=conj(A3+B3);rS=conj(A3-B3);
    lV=-conj(A2-B2);rV=-conj(A2+B2);
    left=conj(A1+B1);right=conj(A1-B1);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  Complex prod;
  complex<Energy> scalar;
  LorentzSpinor<SqrtEnergy> stemp;
  LorentzSpinorBar<SqrtEnergy> sbtemp;
  for(unsigned iya=0;iya<4;++iya) {
    ispin[1]=iya;
    // piece where the vector-spinor is dotted with the momentum of the
    // incoming fermion
    for(unsigned ixa=0;ixa<2;++ixa) {
      unsigned int ix(iya),iy(ixa);
      if(outgoing[0]->id()<0) swap(ix,iy);
      scalar = _inHalf[iy].generalScalar (_inHalfBar[ix],lS,rS);
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      ispin[0]=ixa;
      for(unsigned int iz=0;iz<3;++iz) {
	ispin[2]=iz;
	prod=_inVec[iz].dot(part.momentum())/msum;
	(*ME())(ispin) += (svec.dot(_inVec[iz])+prod*scalar)*
	  UnitRemoval::E/msum/part.mass();
      }
    }
    // the piece where the vector spinor is dotted with the polarization vector
    for(unsigned int iz=0;iz<3;++iz) {
      ispin[2]=iz;
      if(outgoing[0]->id()>0) sbtemp = _inThreeHalfBar[iya].dot(_inVec[iz]);
      else                 stemp  = _inThreeHalf[iya].dot(_inVec[iz]);
      for(unsigned int ixa=0;ixa<2;++ixa) {
	ispin[0]=ixa;
	if(outgoing[0]->id()>0) stemp  = _inHalf[ixa];
	else                    sbtemp = _inHalfBar[ixa];
	(*ME())(ispin) += Complex(stemp.generalScalar(sbtemp,left,right)/part.mass());
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  // test of the matrix element
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // double r2(sqrt(2.)),r3(sqrt(3.));
  // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
  // complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
  // complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m2*A2/msum));
  // complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m2*B2/msum));
  // complex<Energy> h5(-2.*r2/r3/m2/m3*Qp*(0.5*(m12-m22-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2/msum
  // 					 +m12*pcm*pcm*A3/msum/msum));
  // complex<Energy> h6( 2.*r2/r3/m2/m3*Qm*(0.5*(m12-m22-m32)*B1-0.5*Qp*Qp*(m1-m2)*B2/msum
  // 					 +m12*pcm*pcm*B3/msum/msum));
  // cout << "testing 1/2->3/2 1 " << part.id() << " "
  //      << output << "   " 
  //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
  // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass()) << "   " 
  //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
  // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass())/output << endl;
  // return the answer
  return output;
}


// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::
threeHalfHalfScalar(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
  }
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  LorentzPolarizationVector out=UnitRemoval::InvE*momenta[0];
  if(part.id()>0) {
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(out);
  }
  else {
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(out);
  }
  // get the couplings
  Complex A,B;
  Energy msum=part.mass()+momenta[0].mass();
  threeHalfHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A,B);
  Complex left,right;
  // incoming particle
  if(part.id()>0) {
    left=(A-B);
    right=(A+B);
  }
  // incoming anti-particle
  else {
    left=conj(A+B);
    right=conj(A-B);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  for(unsigned ixa=0;ixa<2;++ixa) {
    for(unsigned iya=0;iya<4;++iya) {
      unsigned int iy=iya,ix=ixa;
      if(outgoing[0]->id()<0) swap(ix,iy);
      ispin[0]=iya;
      ispin[1]=ixa;
      (*ME())(ispin) = Complex(_inHalf[iy].generalScalar(_inHalfBar[ix],left,right)*
			       UnitRemoval::E/msum/part.mass());
    }
  }
  // test of the matrix element
  // double test = (ME()->contract(RhoDMatrix(PDT::Spin3Half))).real();
  // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
  // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
  // double r23(sqrt(2./3.));
  // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
  // complex<Energy> h1(-2.*r23*pcm*Qm*B/(m1+m2)),h2( 2.*r23*pcm*Qp*A/(m1+m2));
  // generator()->log() << "testing 3/2->1/2 0 " << part.id() << " "
  // 		     << test << "   " 
  // 		     << 0.125*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass()) << "   " 
  // 		     << 0.125*(h1*conj(h1)+h2*conj(h2))/sqr(part.mass())/test << endl;
  // return the answer
  return (ME()->contract(_rho)).real();;
}

// matrix element for the decay of a spin-3/2 fermion to a spin-3/2 fermion and
// a scalar meson
double Baryon1MesonDecayerBase::threeHalfThreeHalfScalar(const int, const Particle & part,
							 const tPDVector & outgoing,
							 const vector<Lorentz5Momentum> & momenta,
							 MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin3Half,PDT::Spin0)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    // matrix element
  }
  // spinors for the decay product
  // wavefunctions for outgoing fermion
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
    }
  }
  LorentzPolarizationVector in = UnitRemoval::InvE*part.momentum();
  _inHalf.resize(_inThreeHalf.size());
  _inHalfBar.resize(_inThreeHalfBar.size());
  for(unsigned int ix=0;ix<_inThreeHalf.size();++ix) {
    _inHalf[ix] = _inThreeHalf[ix].dot(in);
    _inHalfBar[ix] = _inThreeHalfBar[ix].dot(in);
  }
  // get the couplings
  Complex A1,B1,A2,B2;
  Energy msum(part.mass()+momenta[0].mass());
  threeHalfThreeHalfScalarCoupling(imode(),part.mass(),momenta[0].mass(),
				   momenta[1].mass(),A1,A2,B1,B2);
  Complex left1,right1,left2,right2;
  // incoming particle
  if(part.id()>0) {
    left1=(A1-B1); right1=(A1+B1);
    left2=(A2-B2); right2=(A2+B2);
  }
  // incoming anti-particle
  else {
    left1=(A1+B1); right1=(A1-B1);
    left2=(A2+B2); right2=(A2-B2);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3,0);
  for(unsigned ixa=0;ixa<4;++ixa) {
    for(unsigned iya=0;iya<4;++iya) {
      unsigned int iy=iya,ix=ixa;
      if(outgoing[0]->id()<0) swap(ix,iy);
      ispin[0]=iya;
      ispin[1]=ixa;
      (*ME())(ispin)=Complex((_inThreeHalf[iy].generalScalar(_inThreeHalfBar[ix],left1,right1)
			      +_inHalf[iy].generalScalar( _inHalfBar[ix],left2,right2)
			      *UnitRemoval::E2/sqr(msum))/part.mass());
    }
  }
  // return the answer
  return (ME()->contract(_rho)).real();
}

// matrix element for the decay of a spin-3/2 fermion to a spin-1/2 fermion and
// a vector meson

double Baryon1MesonDecayerBase::
threeHalfHalfVector(const int,const Particle & part, const tPDVector & outgoing,
		    const vector<Lorentz5Momentum> & momenta, MEOption meopt) const {
  if(!ME())
    ME(new_ptr(TwoBodyDecayMatrixElement(PDT::Spin3Half,PDT::Spin1Half,PDT::Spin1)));
  // check if the outgoing meson is really a photon
  bool photon=outgoing[1]->id()==ParticleID::gamma;
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0) {
      RSSpinorWaveFunction   ::calculateWaveFunctions(_inThreeHalf,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
    else {
      RSSpinorBarWaveFunction::calculateWaveFunctions(_inThreeHalfBar,_rho,
						      const_ptr_cast<tPPtr>(&part),
						      incoming);
    }
  }
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalfBar[ix] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ix,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ix=0;ix<2;++ix)
      _inHalf[ix] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ix,Helicity::outgoing);
  }
  LorentzPolarizationVector out=UnitRemoval::InvE*momenta[0];
  if(part.id()>0) {
    _inHalf.resize(_inThreeHalf.size());
    for(unsigned int ix=0;ix<_inThreeHalf.size();++ix)
      _inHalf[ix] = _inThreeHalf[ix].dot(out);
  }
  else {
    _inHalfBar.resize(_inThreeHalfBar.size());
    for(unsigned int ix=0;ix<_inThreeHalfBar.size();++ix)
      _inHalfBar[ix] = _inThreeHalfBar[ix].dot(out);
  }
  // wavefunctions for outgoing fermion
  ME()->zero();
  _inVec.resize(3);
  for(unsigned int ix=0;ix<3;++ix) {
    if(photon && ix==1) continue;
    _inVec[ix] = HelicityFunctions::polarizationVector(-momenta[1],ix,Helicity::outgoing);
  }
  // get the couplings
  Complex A1,A2,A3,B1,B2,B3,prod;
  threeHalfHalfVectorCoupling(imode(),part.mass(),momenta[0].mass(),momenta[1].mass(),
			      A1,A2,A3,B1,B2,B3);
  Energy msum(part.mass()+momenta[0].mass());
  Complex lS,rS,lV,rV,left,right;
  // incoming particle
  if(part.id()>0) {
    lS    = (A3-B3);
    rS    = (A3+B3);
    lV    = (A2-B2);
    rV    = (A2+B2);
    left  = (A1-B1);
    right = (A1+B1);
  }
  // incoming anti-particle
  else {
    lS    = conj(A3+B3);
    rS    = conj(A3-B3);
    lV    =-conj(A2-B2);
    rV    =-conj(A2+B2);
    left  = conj(A1+B1);
    right = conj(A1-B1);
  }
  // compute the matrix element
  vector<unsigned int> ispin(3);
  LorentzVector<complex<Energy> > svec;
  LorentzSpinor<SqrtEnergy> stemp;
  LorentzSpinorBar<SqrtEnergy> sbtemp;
  complex<Energy> scalar;
  for(unsigned iya=0;iya<4;++iya) {
    ispin[0]=iya;
    for(unsigned ixa=0;ixa<2;++ixa) {
      unsigned int iy=iya,ix=ixa;
      if(outgoing[0]->id()<0) swap(ix,iy);
      scalar = _inHalf[iy].generalScalar( _inHalfBar[ix],lS,rS);
      svec   = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV,rV);
      ispin[1]=ixa;
      for(unsigned int iz=0;iz<3;++iz) {
	ispin[2]=iz;
	prod=_inVec[iz].dot(momenta[0])/msum;
	(*ME())(ispin) += (svec.dot(_inVec[iz])+prod*scalar)*
	  UnitRemoval::E/msum/part.mass();
      }
    }
    // the piece where the vector spinor is dotted with the polarization vector
    for(unsigned iz=0;iz<3;++iz) {
      ispin[2]=iz;
      if(outgoing[0]->id()>0) stemp  = _inThreeHalf[iya].dot(_inVec[iz]);
      else                    sbtemp = _inThreeHalfBar[iya].dot(_inVec[iz]);
      for(unsigned int ixa=0;ixa<2;++ixa) {
	ispin[1]=ixa;
	if(outgoing[0]->id()>0) sbtemp = _inHalfBar[ixa];
	else                 stemp  = _inHalf[ixa];
	(*ME())(ispin) += Complex(stemp.generalScalar(sbtemp,left,right)/part.mass());
      }
    }
  }
  double output = (ME()->contract(_rho)).real();
  // testing code
   // Energy m1(part.mass()),m2(momenta[0].mass()),m3(momenta[1].mass());
   // Energy2 m12(m1*m1),m22(m2*m2),m32(m3*m3);
   // Energy Qp(sqrt(sqr(m1+m2)-sqr(m3))),Qm(sqrt(sqr(m1-m2)-sqr(m3)));
   // double r2(sqrt(2.)),r3(sqrt(3.));
   // Energy pcm(Kinematics::pstarTwoBodyDecay(m1,m2,m3));
   // complex<Energy> h1(-2.*Qp*A1),h2(2.*Qm*B1);
   // complex<Energy> h3(-2./r3*Qp*(A1-Qm*Qm/m1*A2/msum));
   // complex<Energy> h4( 2./r3*Qm*(B1-Qp*Qp/m1*B2/msum));
   // complex<Energy> h5(ZERO),h6(ZERO);
   // if(m3!=ZERO) {
   //   h5 = (-2.*r2/r3/m1/m3*Qp*(0.5*(m22-m12-m32)*A1+0.5*Qm*Qm*(m1+m2)*A2/msum
   // 					 +m12*pcm*pcm*A3/msum/msum));
   //   h6 = ( 2.*r2/r3/m1/m3*Qm*(0.5*(m22-m12-m32)*B1-0.5*Qp*Qp*(m2-m1)*B2/msum
   // 			       +m22*pcm*pcm*B3/msum/msum));
   // }
   // cout << "testing 3/2->1/2 1 " << part.id() << " "
   //      << output << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
   // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass()) << "   " 
   //      << 0.25*(h1*conj(h1)+h2*conj(h2)+h3*conj(h3)+
   // 		h4*conj(h4)+h5*conj(h5)+h6*conj(h6))/sqr(part.mass())/output << endl;
  // return the answer
  return output;
}


void Baryon1MesonDecayerBase::halfHalfScalarCoupling(int,Energy,Energy,Energy,
						     Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfScalarCoupling()"
			      << " called from base class this must be implemented"
			      << " in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::halfHalfVectorCoupling(int,Energy,Energy,Energy,
						     Complex&,Complex&,
						     Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfHalfVectorCoupling()" 
			       << " called from base class this must be implemented " 
			       << "in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfScalarCoupling"
			       << "() called from base class this must be implemented"
			       << " in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::halfThreeHalfVectorCoupling"
			       << "() called from base class this must be implemented "
			       << "in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::threeHalfHalfScalarCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfHalfScalarCoupling"
			       << "() called from base class this must be implemented"
			       << " in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::threeHalfHalfVectorCoupling(int,Energy,Energy,Energy,
							  Complex&,Complex&,
							  Complex&,Complex&,
							  Complex&,Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfHalfVectorCoupling"
			       << "() called from base class this must be implemented "
			       << "in the inheriting class" << Exception::abortnow;
}

void Baryon1MesonDecayerBase::threeHalfThreeHalfScalarCoupling(int,Energy,Energy,
							       Energy,Complex&,
							       Complex&,Complex&,
							       Complex&) const {
  throw DecayIntegratorError() << "Baryon1MesonDecayerBase::threeHalfThreeHalfScalar"
			       << "Coupling() called from base class this must be "
			       << "implemented in the inheriting class" 
			       << Exception::abortnow;
}

bool Baryon1MesonDecayerBase::twoBodyMEcode(const DecayMode & dm,int & mecode,
					    double & coupling) const {
  coupling=1.;
  unsigned int inspin(dm.parent()->iSpin()),outspin,outmes;
  ParticleMSet::const_iterator pit(dm.products().begin());
  bool order; 
  if((**pit).iSpin()%2==0) {
    order=true;
    outspin=(**pit).iSpin();
    ++pit;outmes=(**pit).iSpin();
  }
  else {
    order=false;
    outmes=(**pit).iSpin();++pit;
    outspin=(**pit).iSpin();
  }
  mecode=-1;
  if(inspin==2) {
    if(outspin==2){if(outmes==1){mecode=101;}else{mecode=102;}}
    else if(outspin==4){if(outmes==1){mecode=103;}else{mecode=104;}}
  }
  else if(inspin==4) {
    if(outspin==2){if(outmes==1){mecode=105;}else{mecode=106;}}
    else if(outspin==4){if(outmes==1){mecode=107;}else{mecode=108;}}
  }
  return order;
}

void Baryon1MesonDecayerBase::dataBaseOutput(ofstream & os,bool header) const {
  DecayIntegrator::dataBaseOutput(os,header);
}
#line 1 "./BaryonFactorizedDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BaryonFactorizedDecayer class.
//

#include "BaryonFactorizedDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "Herwig/Decay/DecayVertex.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

BaryonFactorizedDecayer::BaryonFactorizedDecayer() {
  // default values taken from PRD56, 2799
  _a1c= 1.1;
  _a2c=-0.5;
  _a1b= 1.0;
  _a2b= 0.28;
  // intermediates
  generateIntermediates(true);
}

void BaryonFactorizedDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator::doinitrun();
  _weights.clear();_wgtloc.clear();_wgtmax.clear();
  for(unsigned int ix=0;ix<numberModes();++ix) {
    _wgtmax.push_back(mode(ix)->maxWeight());
    _wgtloc.push_back(_weights.size());
    for(unsigned int iy=0;iy<mode(ix)->channels().size();++iy)
      _weights.push_back(mode(ix)->channels()[iy].weight());
  }
}

void BaryonFactorizedDecayer::doinit() {
  DecayIntegrator::doinit();
  // get the CKM matrix (unsquared for interference)
  Complex ckmmat[3][3];
  vector< vector<Complex > > CKM(_theCKM->getUnsquaredMatrix(SM().families()));
  for(unsigned int ix=0;ix<3;++ix) {
    for(unsigned int iy=0;iy<3;++iy) {
      ckmmat[ix][iy]=CKM[ix][iy];
    }
  }
  // make sure the current and form factor got initialised
  _current->init();
  _form->init();
  // find all the possible modes
  vector<unsigned int> tformmap,tcurrmap;
  vector<int> inquark,outquark,currq,curra;
  tPDVector incoming;
  vector<tPDVector> outgoing;
  for(unsigned int iform=0;iform<_form->numberOfFactors();++iform) {
    // particles from the form factor
    int id0,id1;
    _form->particleID (iform,id0,id1);
    int spect1,spect2,inq,outq,ispin,ospin;
    _form->formFactorInfo(iform,ispin,ospin,spect1,spect2,inq,outq);
    // particles from the form factor
    tPDPtr in  = getParticleData(id0);
    tPDPtr out = getParticleData(id1);
    // the charge of the decay products
    int Wcharge = in->iCharge()-out->iCharge();
    // max mass for the particles in the current
    Energy min = in->massMax()-out->massMin();
    for(unsigned int icurr=0;icurr<_current->numberOfModes();++icurr) {
      // get the particles from the current
      int iq,ia;
      _current->decayModeInfo(icurr,iq,ia);
      tPDVector ptemp=_current->particles(Wcharge,icurr,iq,ia);
      tPDVector outV = {out};
      outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
      Energy minb=ZERO;
      for(unsigned int iz=0;iz<ptemp.size();++iz) minb+=ptemp[iz]->massMin();
      // valid mode
      if(outV.size()>1&&minb<min&&
	 (Wcharge!=0||(Wcharge==0&&((inq>0&&inq%2!=iq%2)||
				    (inq<0&&abs(inq)%2!=abs(ia)%2))))) {
	tformmap.push_back(iform);tcurrmap.push_back(icurr);
	incoming.push_back(in);
	outgoing.push_back(outV);
	inquark.push_back(inq);outquark.push_back(outq);
	currq.push_back(iq);curra.push_back(ia);
      }
      // if the meson is neutral try the CC mode
      if(Wcharge==0&&iq!=-ia&&((inq>0&&inq%2!=iq%2)||
			       (inq<0&&abs(inq)%2!=abs(ia)%2))) {
	ptemp=_current->particles(Wcharge,icurr,-ia,-iq);
	tPDVector outV = {out};
	outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
	Energy minb=ZERO;
	for(unsigned int iz=0;iz<ptemp.size();++iz) minb+=ptemp[iz]->massMin();
	if(outV.size()>1&&minb<min) {
  	  tformmap.push_back(iform);tcurrmap.push_back(icurr);
	  incoming.push_back(in);
  	  outgoing.push_back(outV);
  	  inquark.push_back(inq);outquark.push_back(outq);
  	  currq.push_back(-ia);curra.push_back(-iq);
  	}
      }
    }
  }
  _formmap.clear();
  _currentmap.clear();
  // loop over the modes and find the dupliciates
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    while(true) {
      if ( outgoing[ix].empty() ) break;
      vector<unsigned int> modeloc;
      vector<bool> modecc;
      findModes(ix,incoming,outgoing,modeloc,modecc);
      // if more than two outgoing particles only allow one diagram
      if ( outgoing[ix].size() > 2 && !modeloc.empty() ) {break;}
      // create the mode and set the particles as for the first instance
      PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],1.));
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      Energy min = incoming[ix]->massMax()-outgoing[ix][0]->massMin();
      int Wcharge = incoming[ix]->iCharge()-outgoing[ix][0]->iCharge();
      bool done = _current->createMode(Wcharge,tcPDPtr(),FlavourInfo(),
				       tcurrmap[ix],mode,1,0,channel,min);
      if(!done){throw InitException() << "Failed to construct mode in "
				      << "BaryonFactorizedDecayer::doinit()." 
				      << Exception::abortnow;}
      // set the parameters for the additional modes
      vector<unsigned int>ttform,ttcurr;
      ttform.push_back(tformmap[ix]);ttcurr.push_back(tcurrmap[ix]);
      for(unsigned int iy=0;iy<modeloc.size();++iy) {
	ttform.push_back(tformmap[modeloc[iy]]);
	ttcurr.push_back(tcurrmap[modeloc[iy]]);
      }
      vector<Complex> tCKM;
      for(unsigned int iy=0;iy<ttcurr.size();++iy) {
	// get the quarks involved in the process
	int iq,ia,inq,outq;
	if(iy==0) {
	  iq=currq[ix];
	  ia=curra[ix];
	  inq=inquark[ix];
	  outq=outquark[ix];
	}
	else {
	  if(!modecc[iy-1]) {
	    iq=currq[modeloc[iy-1]];
	    ia=curra[modeloc[iy-1]];
	    inq=inquark[modeloc[iy-1]];
	    outq=outquark[modeloc[iy-1]];
	  }
	  else {
	    ia=-currq[modeloc[iy-1]];
	    iq=-curra[modeloc[iy-1]];
	    inq=-inquark[modeloc[iy-1]];
	    outq=-outquark[modeloc[iy-1]];
	  }
 	}
	int id0,id1;
	_form->particleID(ttform[iy],id0,id1);
	int Wcharge = getParticleData(id0)->iCharge()-getParticleData(id1)->iCharge();
	Complex ckm=1.;
	if(Wcharge!=0) {
	  if(abs(iq)%2==0){ckm *= conj(ckmmat[abs(iq)/2-1][(abs(ia)-1)/2]);}
	  else{ckm *= conj(ckmmat[abs(ia)/2-1][(abs(iq)-1)/2]);}
	  if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(outq)-1)/2];}
	  else{ckm *= ckmmat[abs(outq)/2-1][(abs(inq)-1)/2];}
	  if(abs(inq)==5){ckm*=_a1b;}
	  else{ckm*=_a1c;}
	}
	else {
	  if(inq>0) {
	    if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(iq)-1)/2];}
	    else{ckm *= ckmmat[abs(iq)/2-1][(abs(inq)-1)/2];}
	    if(abs(outq)%2==0)
	      {ckm *= conj(ckmmat[abs(outq)/2-1][(abs(ia)-1)/2]);}
	    else{ckm *= conj(ckmmat[abs(ia)/2-1][(abs(outq)-1)/2]);}
	  }
	  else {
	    if(abs(inq)%2==0){ckm *= ckmmat[abs(inq)/2-1][(abs(ia)-1)/2];}
	    else{ckm *= ckmmat[abs(ia)/2-1][(abs(inq)-1)/2];}
	    if(abs(outq)%2==0)
	      {ckm *= conj(ckmmat[abs(outq)/2-1][(abs(iq)-1)/2]);}
	    else{ckm *= conj(ckmmat[abs(iq)/2-1][(abs(outq)-1)/2]);}
	  }
	  if(abs(inq)==5){ckm*=_a2b;}
	  else{ckm*=_a2c;}
	}
	if((abs(inq)%2==0&&inq<0)||(abs(inq)%2!=0&&inq>0)){ckm=conj(ckm);}
	tCKM.push_back(ckm);
      }
      // add the parameters for the mode to the list
      _currentmap.push_back(ttcurr);
      _formmap.push_back(ttform);
      _factCKM.push_back(tCKM);
      double maxweight(0.);
      // add the mode to the list
      if(_wgtmax.size()>numberModes()) maxweight=_wgtmax[numberModes()];
      // the weights for the channel
      vector<double> channelwgts;
      if(_wgtloc.size()>numberModes()&&
	 _wgtloc[numberModes()]+mode->channels().size()<=_weights.size()) {
	vector<double>::iterator start=_weights.begin()+_wgtloc[numberModes()];
	vector<double>::iterator end  = start+mode->channels().size();
	channelwgts=vector<double>(start,end);
      }
      else {
	channelwgts.resize(mode->channels().size(),1./(mode->channels().size()));
      }
      // don't need channels for two body decays
      if(outgoing[ix].size()==2) {
	channelwgts.clear();
	mode = new_ptr(PhaseSpaceMode(incoming[ix],outgoing[ix],maxweight));
      }
      else {
       mode->maxWeight(maxweight);
       mode->setWeights(channelwgts);
      }
      addMode(mode);
      // resize the duplicate modes to remove them
      for(unsigned int iy=0;iy<modeloc.size();++iy) outgoing[modeloc[iy]]=tPDVector();
      break;
    }
  }
}

bool BaryonFactorizedDecayer::accept(tcPDPtr parent, const tPDVector & children) const {
  bool allowed=false;
  unsigned int iform(0),ix;
  int idin(parent->id()),ibaryon,foundb,id0,id1;
  vector<int> idall,idother;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit){idall.push_back((**pit).id());}
  // loop over the particles in the form factor
  do {
    _form->particleID(iform,id0,id1);
    ibaryon=0;
    if(id0==idin){ibaryon=id1;}
    else if(id0==-idin){ibaryon=-id1;}
    if(ibaryon!=0) {
      foundb=false;
      idother.clear();
      for(ix=0;ix<idall.size();++ix) {
	if(idall[ix]==ibaryon){foundb=true;}
	else{idother.push_back(idall[ix]);}
      }
      if(foundb){allowed=_current->accept(idother);}
    }
    ++iform;
  }
  while(!allowed&&iform<_form->numberOfFactors());
  return allowed;
}

int BaryonFactorizedDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  unsigned int ix,iy;
  int idin(parent->id()),ibaryon,foundb,id0,id1,icurr(-1),iform(0);
  vector<int> idall,idother;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit){idall.push_back((**pit).id());}
  // loop over the particles in the form factor
  do
    {
      _form->particleID(iform,id0,id1);
      ibaryon=0;
      if(id0==idin){ibaryon=id1;}
      else if(id0==-idin){ibaryon=-id1;}
      ++iform;
      foundb=false;
      idother.clear();
      for(ix=0;ix<idall.size();++ix)
	{
	  if(idall[ix]==ibaryon){foundb=true;}
	  else{idother.push_back(idall[ix]);}
	}
      if(foundb){icurr=_current->decayMode(idother);}
    }
  while(icurr<0&&iform<int(_form->numberOfFactors()));
  // now find the mode
  int imode=-1;
  ix=0;
  --iform;
  do
    {
      for(iy=0;iy<_currentmap[ix].size();++iy)
	{if(int(_currentmap[ix][iy])==icurr&&int(_formmap[ix][iy])==iform){imode=ix;}}
      ++ix;
    }
  while(imode<0&&ix<numberModes());
  if(imode<0){throw DecayIntegratorError() << "Unable to find the mode in " 
					   << "BaryonFactorizedDecayer::decay()" 
					   << Exception::abortnow;}
  // generate the mode
  cc=id0!=idin;
  return imode;
}


void BaryonFactorizedDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _a1b << _a2b <<_a1c <<_a2c 
     << _currentmap << _formmap << _factCKM << _wgtloc << _wgtmax << _weights 
     << _theCKM;
}

void BaryonFactorizedDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _a1b >> _a2b >>_a1c >>_a2c 
     >> _currentmap >> _formmap >> _factCKM >> _wgtloc >> _wgtmax >> _weights 
     >> _theCKM;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<BaryonFactorizedDecayer,DecayIntegrator>
describeHerwigBaryonFactorizedDecayer("Herwig::BaryonFactorizedDecayer", "HwBaryonDecay.so");

void BaryonFactorizedDecayer::Init() {

  static ClassDocumentation<BaryonFactorizedDecayer> documentation
    ("The BaryonFactorizedDecayer class combines the baryon form factor and a"
     " weak current to perform a decay in the naive factorization approximation.");

  static Reference<BaryonFactorizedDecayer,WeakCurrent> interfaceWeakCurrent
    ("Current",
     "The reference for the decay current to be used.",
     &BaryonFactorizedDecayer::_current, false, false, true, false, false);

  static ParVector<BaryonFactorizedDecayer,int> interfaceWeightLocation
    ("WeightLocation",
     "The locations of the weights for a given channel in the vector",
     &BaryonFactorizedDecayer::_wgtloc,
     0, 0, 0, 0, 10000, false, false, true);

  static ParVector<BaryonFactorizedDecayer,double> interfaceWeightMax
    ("MaximumWeight",
     "The maximum weight for a given channel.",
     &BaryonFactorizedDecayer::_wgtmax,
     0, 0, 0, 0., 100., false, false, true);

  static ParVector<BaryonFactorizedDecayer,double> interfaceWeights
    ("Weights",
     "The weights for the integration.",
     &BaryonFactorizedDecayer::_weights,
     0, 0, 0, 0., 1., false, false, true);

  static Reference<BaryonFactorizedDecayer,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form-factor",
     &BaryonFactorizedDecayer::_form, true, true, true, false, false);

  static Parameter<BaryonFactorizedDecayer,double> interfacea1Bottom
    ("a1Bottom",
     "The factorization paramter a_1 for decays of bottom baryons",
     &BaryonFactorizedDecayer::_a1b, 1., -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea2Bottom
    ("a2Bottom",
     "The factorization paramter a_2 for decays of bottom baryons",
     &BaryonFactorizedDecayer::_a2b, 0.28, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea1Charm
    ("a1Charm",
     "The factorization paramter a_1 for decays of charm baryons",
     &BaryonFactorizedDecayer::_a1c, 1.1, -10.0, 10.0,
     false, false, true);

  static Parameter<BaryonFactorizedDecayer,double> interfacea2Charm
    ("a2Charm",
     "The factorization paramter a_2 for decays of charm baryons",
     &BaryonFactorizedDecayer::_a2c, -0.5, -10.0, 10.0,
     false, false, true);

  static Reference<BaryonFactorizedDecayer,StandardCKM> interfaceCKM
    ("CKM",
     "Reference to the Standard Model object",
     &BaryonFactorizedDecayer::_theCKM, false, false, true, false);
}

void BaryonFactorizedDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()>0) {
    SpinorWaveFunction::
      constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
  }
  // decay product
  // spin 1/2
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
    if(part.id()>0) {
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else {
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
  }
  // spin 3/2
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
    if(part.id()>0) {
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
  						 decay[0],outgoing,true);
      
    }
    else {
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
  					      decay[0],outgoing,true);
    }
  }
  else
    assert(false);
  // and the stuff from the current
  _current->constructSpinInfo(ParticleVector(decay.begin()+1,decay.end()));
}

double BaryonFactorizedDecayer::me2(const int ichan, const Particle & part,
				    const tPDVector & outgoing,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  double me(0.);
  assert(part.dataPtr()->iSpin()==PDT::Spin1Half);
  if(outgoing[0]->iSpin()==PDT::Spin1Half)
    me=halfHalf(ichan,part,outgoing,momenta,meopt);
  else if(outgoing[0]->iSpin()==PDT::Spin3Half)
    me=halfThreeHalf(ichan,part,outgoing,momenta,meopt);
  else
    assert(false);
  return me;
}

// matrix element for a 1/2 -> 1/2 decay
double BaryonFactorizedDecayer::halfHalf(const int ichan, const Particle & part,
					 const tPDVector & outgoing,
					 const vector<Lorentz5Momentum> & momenta,
					 MEOption meopt) const {
  Energy scale;
  // extract the spins of the particles
  vector<PDT::Spin> spin;
  for(unsigned ix=0;ix<outgoing.size();++ix) 
    spin.push_back(outgoing[ix]->iSpin());
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,spin)));
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
						    const_ptr_cast<tPPtr>(&part),
						    incoming);
  }
  ME()->zero();
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel)
      _inHalfBar[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ihel,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel)
      _inHalf[ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ihel,Helicity::outgoing);
  }
  // get the information on the form-factor
  int id0(part.id()),id1(outgoing[0]->id());
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the baryonic part of the current for the decay
  double pre(0.);
  for(unsigned int mode=0;mode<_formmap[imode()].size();++mode) {
    Complex f1v,f2v,f3v,f1a,f2a,f3a;
    // calculate the form factor piece
    _form->SpinHalfSpinHalfFormFactor(q2,_formmap[imode()][mode],id0,id1,m0,m1,
				      f1v,f2v,f3v,f1a,f2a,f3a,FlavourInfo());
    Complex  left = f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
    Complex right = f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
    vector<LorentzPolarizationVectorE> baryon(4);
    for(unsigned int ix=0;ix<2;++ix) {
      for(unsigned int iy=0;iy<2;++iy) {
	LorentzPolarizationVectorE 
	  vtemp = _inHalf[ix].generalCurrent(_inHalfBar[iy],left,right);
	complex<Energy> vspin=_inHalf[ix].scalar(_inHalfBar[iy]);      
	complex<Energy> aspin=_inHalf[ix].pseudoScalar(_inHalfBar[iy]);
	// the momentum like pieces
 	if(part.id()>0) {
 	  vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
 	  vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
 	}
 	else {
 	  vtemp+= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
 	  vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
 	}
 	if(part.id()>0) baryon[2*ix+iy]=vtemp;
 	else            baryon[2*iy+ix]=vtemp;
      }
    }
    // construct the weak current
    vector<LorentzPolarizationVectorE> hadron =
      _current->current(tcPDPtr(),FlavourInfo(),
			_currentmap[imode()][mode],ichan,scale,
			tPDVector(outgoing.begin()+1,outgoing.end()),
			vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt);
    pre=pow(part.mass()/scale,int(outgoing.size()-3));pre*=pre;
    vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
    int itemp=1;
    unsigned int ibar=0;
    for(int iz=int(outgoing.size()-1);iz>=0;--iz) {
      if(abs(outgoing[iz]->id())!=id1) {
 	itemp *= outgoing[iz]->iSpin();
 	constants[iz]=itemp;
      }
      else ibar=iz;
      constants[outgoing.size()]=1;
      constants[ibar]=constants[ibar+1];
    }
    for(unsigned int mhel=0;mhel<baryon.size();++mhel) {
      ihel[0     ]=mhel/2;
      ihel[ibar+1]=mhel%2;
      for(unsigned int lhel=0;lhel<hadron.size();++lhel) {
 	// map the index for the hadrons to a helicity state
 	for(unsigned int ix=outgoing.size();ix>0;--ix) {
 	  if(ix-1!=ibar){ihel[ix]=(lhel%constants[ix-1])/constants[ix];}}
 	(*ME())(ihel) += Complex(hadron[lhel].dot(baryon[mhel])*
				 _factCKM[imode()][mode]*SM().fermiConstant());
      }
    }
  }
  // return the answer
  return 0.5*pre*(ME()->contract(_rho)).real();
}

// matrix element for a 1/2 -> 3/2 decay
double BaryonFactorizedDecayer::halfThreeHalf(const int ichan, const Particle & part,
					      const tPDVector & outgoing,
					      const vector<Lorentz5Momentum> & momenta,
					      MEOption meopt) const {
  // spins
  Energy scale;
  vector<PDT::Spin> spin(outgoing.size());
  for(unsigned int ix=0;ix<outgoing.size();++ix)
    spin[ix]=outgoing[ix]->iSpin();
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,spin)));
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
  }
  ME()->zero();
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  // get the information on the form-factor
  int id0(part.id()),id1(outgoing[0]->id());
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  InvEnergy ms(1./(m0+m1));
  InvEnergy2 ms2(ms*ms);
  double pre(0.);
  for(unsigned int mode=0;mode<_formmap[imode()].size();++mode) {
    // calculate the form factors
    Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
    _form->SpinHalfSpinThreeHalfFormFactor(q2,_formmap[imode()][mode],id0,id1,m0,m1,
  					   f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a,FlavourInfo());
    complex<InvEnergy2> lS1,lS2,rS1,rS2;
    Complex left,right;
    complex<InvEnergy> lV,rV;
    if(part.id()>0) {
      left  = f1a-f1v;
      right = f1a+f1v; 
      lS1   = ms2*(f3a-f4a-f3v+f4v);
      rS1   = ms2*(f3a-f4a+f3v-f4v);
      lS2   = ms2*(f4a-f4v);
      rS2   = ms2*(f4a+f4v);
      lV    = ms*(f2a-f2v);
      rV    = ms*(f2a+f2v);
    }
    else {
      left  = conj(f1a+f1v);
      right = conj(f1a-f1v); 
      lS1   = ms2*conj(f3a-f4a+f3v-f4v);
      rS1   = ms2*conj(f3a-f4a-f3v+f4v);
      lS2   = ms2*conj(f4a-f4v);
      rS2   = ms2*conj(f4a+f4v);
      lV    = ms *conj(f2a-f2v);
      rV    = ms *conj(f2a+f2v);
    }
    // construct the vectors for the decay
    LorentzPolarizationVectorE baryon[4][2];
    for(unsigned int iya=0;iya<4;++iya) {
      for(unsigned int ixa=0;ixa<2;++ixa) {
	unsigned int ix,iy;
  	if(outgoing[0]->id()>0) {
	  ix=iya;
	  iy=ixa;
	}
	else {
	  ix=ixa;
	  iy=iya;
	}
  	// scalar like terms
	complex<Energy> lfact = _inHalf[iy].leftScalar( _inHalfBar[ix]);
	complex<Energy> rfact = _inHalf[iy].rightScalar(_inHalfBar[ix]);
	Complex       scalar1 = (lS1*lfact+rS1*rfact)*UnitRemoval::E;
	Complex       scalar2 = (lS2*lfact+rS2*rfact)*UnitRemoval::E;
	LorentzPolarizationVector  svec = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV/ms,rV/ms)*ms;
	LorentzPolarizationVectorE tvec;
   	if(part.id()>0) {
   	  tvec=_inThreeHalfBar[ix].generalCurrent(_inHalf[iy],left,right);
   	}
   	else {
   	  tvec=_inThreeHalf[iy].generalCurrent(_inHalfBar[ix],left,right);
   	}
  	baryon[iya][ixa] = tvec+svec*UnitRemoval::E
  	  +scalar1*momenta[0]+scalar2*part.momentum();
      }
    }
    vector<LorentzPolarizationVectorE> hadron =
      _current->current(tcPDPtr(),FlavourInfo(),_currentmap[imode()][mode],ichan,scale,
			tPDVector(outgoing.begin()+1,outgoing.end()),
			vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt);
    // prefactor
    pre  = pow(part.mass()/scale,int(outgoing.size()-3));
    pre *= pre;
    // work out the mapping for the hadron vector
    vector<unsigned int> constants(outgoing.size()+1),ihel(outgoing.size()+1);
    int itemp = 1;
    int ibar  = 0;
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      if(abs(outgoing[ix]->id())!=id1) {
   	itemp*=outgoing[ix]->iSpin();
   	constants[ix]=itemp;
      }
      else{ibar=ix;}
    }
    constants[outgoing.size()]=1;
    constants[ibar]=constants[ibar+1];
    for(unsigned int iya=0;iya<4;++iya) {
      ihel[1]=iya;
      for(unsigned int ixa=0;ixa<2;++ixa) {
   	ihel[0]=ixa;
	for(unsigned int lhel=0;lhel<hadron.size();++lhel) {
	  // map the index for the hadrons to a helicity state
	  for(int ix=int(outgoing.size());ix>0;--ix)
	    {if(ix-1!=ibar){ihel[ix]=(lhel%constants[ix-1])/constants[ix];}}
	  (*ME())(ihel) += Complex(hadron[lhel].dot(baryon[iya][ixa])*
				   _factCKM[imode()][mode]*SM().fermiConstant());
	}
      }
    }
  }
  // return the answer
  return 0.5*pre*(ME()->contract(_rho)).real();
}

void BaryonFactorizedDecayer::findModes(unsigned int imode,
					tPDVector & incoming,
					vector<tPDVector> & outgoing,
					vector<unsigned int> & loc,
					vector<bool> & cc) {
  // get the id's for the mode
  // incoming
  int id_in    = incoming[imode]->id();
  int idbar_in = incoming[imode]->CC() ?
    incoming[imode]->CC()->id() : incoming[imode]->id();
  // outgoing
  vector<int> id_out,idbar_out;
  for(unsigned int ix=0;ix<outgoing[imode].size();++ix) {
    id_out.push_back(outgoing[imode][ix]->id());
    if(outgoing[imode][ix]->CC())
      idbar_out.push_back(outgoing[imode][ix]->CC()->id());
    else
      idbar_out.push_back(id_out[ix]);
  }
  // loop over the modes
  for(unsigned int ix=0;ix<outgoing.size();++ix) {
    if(ix==imode||outgoing[ix].empty()) continue;
    assert(!outgoing[ix].empty());
    assert(incoming[ix]);
    // the particle mode
    if(incoming[ix]->id()==id_in&&outgoing[ix].size()==id_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound = 0;
      for(unsigned int iy=0;iy<id_out.size();++iy) {
	int idtemp = outgoing[ix][iy]->id();
	unsigned int iz = 0;
	bool found = false;
	do {
	  if(idtemp==id_out[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<id_out.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==id_out.size()) {
	cc.push_back(false);
	loc.push_back(ix);
      }
    }
    // the charge conjugate mode
    if(incoming[ix]->id()==idbar_in&&outgoing[ix].size()==idbar_out.size()) {
      vector<bool> done(id_out.size(),false);
      unsigned int nfound=0;
      for(unsigned int iy=0;iy<idbar_out.size();++iy) {
	int idtemp=outgoing[ix][iy]->id();
	unsigned int iz=0;
	bool found = false;
	do {
	  if(idtemp==idbar_out[iz]&&!done[iz]) {
	    done[iz]=true;
	    found=true;
	  }
	  ++iz;
	}
	while(iz<idbar_out.size()&&!found);
	if(found) ++nfound;
      }
      if(nfound==idbar_out.size()) {
	cc.push_back(false);
	loc.push_back(ix);
      }
    }
  }
}

// output the setup information for the particle database
void BaryonFactorizedDecayer::dataBaseOutput(ofstream & output, bool header) const {
  unsigned int ix;
  if(header){output << "update decayers set parameters=\"";}
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":a1Bottom "  << _a1b << "\n";
  output << "newdef " << name() << ":a2Bottom "  << _a2b << "\n";
  output << "newdef " << name() << ":a1Charm "   << _a1c << "\n";
  output << "newdef " << name() << ":a2Charm "   << _a2c << "\n";
  output << "newdef " << name() << ":CKM "       << _theCKM->name() << " \n";
  for(ix=0;ix<_wgtloc.size();++ix)
    {output << "insert " << name() << ":WeightLocation " << ix << " " 
	    << _wgtloc[ix] << "\n";}
  for(ix=0;ix<_wgtmax.size();++ix)
    {output << "insert " << name() << ":MaximumWeight "  << ix << " " 
	    << _wgtmax[ix] << "\n";}
  for(ix=0;ix<_weights.size();++ix)
    {output << "insert " << name() << ":Weights "        << ix << " " 
	    << _weights[ix] << "\n";}
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":Current " << _current->name() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":FormFactor " << _form->name() << " \n";
  if(header){output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;}
}
#line 1 "./KornerKramerCharmDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the KornerKramerCharmDecayer class.
//

#include "KornerKramerCharmDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

KornerKramerCharmDecayer::KornerKramerCharmDecayer() {
  // one over the number of colours 
  oneNC_=0.;
  // pseudoscalar meson decay constants
  fpi_ = 131.7*MeV;
  fk_  = 160.6*MeV;
  // vector decay constants
  frho_   = 0.272;
  fKstar_ = 0.238;
  // masses for the form-factors for the factorizing diagrams
  mdcplus_  = 2.42*GeV;
  mdcminus_ = 2.01*GeV;
  mscplus_  = 2.54*GeV;
  mscminus_ = 2.11*GeV;
  // perturbative factors
  cplus_  = 0.73;
  cminus_ = 1.90;
  // factors for the non-factorizing diagrams
  H2_ = 0.119*GeV;
  H3_ =-0.011*GeV;
  // lambda_c to lambda pi+
  incoming_.push_back(4122);outgoingB_.push_back(3122);outgoingM_.push_back(211);
  maxweight_.push_back(0.0153611);
  I1_.push_back(-1./3.);I2_.push_back(-1./3.);
  I3_.push_back(-1./3.);I4_.push_back( 2./3.);I5_.push_back( 1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // lambda_c to sigma_0 pi+
  incoming_.push_back(4122);outgoingB_.push_back(3212);outgoingM_.push_back(211);
  maxweight_.push_back(0.00684877);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-1./sqrt(3.));I4_.push_back(0.);I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(1./sqrt(3.));Ihat4_.push_back(-2./sqrt(3.));
  // lambda_c to  sigma+ pi0
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(111);
  maxweight_.push_back(0.00686766);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(1./sqrt(3.));I4_.push_back(0.);I5_.push_back(1./sqrt(12.));
  Ihat3_.push_back(-1./sqrt(3.));Ihat4_.push_back(2./sqrt(3.));
  // lambda_c+ to sigma+ eta
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(221);
  maxweight_.push_back(0.00311807);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-0.169101981);I4_.push_back(0.577350262);I5_.push_back(0.204124141);
  Ihat3_.push_back(0.408248281);Ihat4_.push_back(-0.816496562);
  // lambda_c+ to sigma+ eta'
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(331);
  maxweight_.push_back(0.0267591);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.985598555);I4_.push_back(-0.577350261);I5_.push_back(0.204124147);
  Ihat3_.push_back( 0.408248294);Ihat4_.push_back(-0.816496587);
  // lambda_c to p k0
  incoming_.push_back(4122);outgoingB_.push_back(2212);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0444906);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c+ to xi K+
  incoming_.push_back(4122);outgoingB_.push_back(3322);outgoingM_.push_back(321);
  maxweight_.push_back(0.00533608);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c to sigma+ Kbar0
  incoming_.push_back(4232);outgoingB_.push_back(3222);outgoingM_.push_back(-311);
  maxweight_.push_back(0.135814);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(6.));Ihat4_.push_back(4./sqrt(6.));
  // xi_c to xi_0 pi+
  incoming_.push_back(4232);outgoingB_.push_back(3322);outgoingM_.push_back(211);
  maxweight_.push_back(0.0745085);
  I1_.push_back(1./sqrt(6.));I2_.push_back(1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2./sqrt(6.));Ihat4_.push_back(-4./sqrt(6.));
  // xi_c to lambda kbar0
  incoming_.push_back(4132);outgoingB_.push_back(3122);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0025821);
  I1_.push_back(1./6.);I2_.push_back(1./6.);
  I3_.push_back(2./3.);I4_.push_back(-1./3.);I5_.push_back(1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // xi_c to sigma k0
  incoming_.push_back(4132);outgoingB_.push_back(3212);outgoingM_.push_back(-311);
  maxweight_.push_back(0.024961);
  I1_.push_back(1./sqrt(12.));I2_.push_back(1./sqrt(12.));
  I3_.push_back(0.);I4_.push_back(-2./sqrt(12.));I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(2./sqrt(12.));Ihat4_.push_back(-4./sqrt(12.));
  // xi_c to sigma k-
  incoming_.push_back(4132);outgoingB_.push_back(3222);outgoingM_.push_back(-321);
  maxweight_.push_back(0.00250939);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi0 pi0
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(111);
  maxweight_.push_back(0.000739771);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(2./sqrt(12.));I4_.push_back(-2./sqrt(12.));I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(12.));Ihat4_.push_back(4./sqrt(12.));
  // xi_c0 to xi0 eta
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(221);
  maxweight_.push_back(0.00490303);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-0.169101981);I4_.push_back(-0.408248281);I5_.push_back(-0.288675131);
  Ihat3_.push_back(0.408248281);Ihat4_.push_back(-0.816496562);
  // xi_c0 to xi0 eta'
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(331);
  maxweight_.push_back(0.0178693);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.985598555);I4_.push_back(-0.408248294);I5_.push_back(0.288675131);
  Ihat3_.push_back(0.408248294);Ihat4_.push_back(-0.816496587);
  // xi_c0 to xi- pi+
  incoming_.push_back(4132);outgoingB_.push_back(3312);outgoingM_.push_back(211);
  maxweight_.push_back(0.0220038);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c to xi0 K0
  incoming_.push_back(4332);outgoingB_.push_back(3322);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0250268);
  I1_.push_back(1.);I2_.push_back(-1.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2.);Ihat4_.push_back(0.);
  // lambda_c to lambda rho+
  incoming_.push_back(4122);outgoingB_.push_back(3122);outgoingM_.push_back(213);
  maxweight_.push_back(0.467008);
  I1_.push_back(-1./3.);I2_.push_back(-1./3.);
  I3_.push_back(-1./3.);I4_.push_back( 2./3.);I5_.push_back( 1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // lambda_c to sigma_0 rho+
  incoming_.push_back(4122);outgoingB_.push_back(3212);outgoingM_.push_back(213);
  maxweight_.push_back(0.0702498);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(-1./sqrt(3.));I4_.push_back(0.);I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(1./sqrt(3.));Ihat4_.push_back(-2./sqrt(3.));
  // lambda_c to  sigma+ rho
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(113);
  maxweight_.push_back(0.0697802);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(1./sqrt(3.));I4_.push_back(0.);
  I5_.push_back(1./sqrt(12.));
  Ihat3_.push_back(-1./sqrt(3.));Ihat4_.push_back(2./sqrt(3.));
  // lambda_c to sigma+ omega
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(223);
  maxweight_.push_back(0.252355);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.57735033);I4_.push_back(0.);I5_.push_back( 0.288675151);
  Ihat3_.push_back(0.577350303);Ihat4_.push_back( -1.15470061);
  // lambda_c to simga+ phi
  incoming_.push_back(4122);outgoingB_.push_back(3222);outgoingM_.push_back(333);
  maxweight_.push_back(0.0063064);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.816496614);I4_.push_back(-0.816496624);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to p k*0bar
  incoming_.push_back(4122);outgoingB_.push_back(2212);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0996461);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c+ to xi K*+
  incoming_.push_back(4122);outgoingB_.push_back(3322);outgoingM_.push_back(323);
  maxweight_.push_back(0.00413946);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c to sigma+ K*bar0
  incoming_.push_back(4232);outgoingB_.push_back(3222);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0583248);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(6.));Ihat4_.push_back(4./sqrt(6.));
  // xi_c to xi_0 rho+
  incoming_.push_back(4232);outgoingB_.push_back(3322);outgoingM_.push_back(213);
  maxweight_.push_back(2.18754);
  I1_.push_back(1./sqrt(6.));I2_.push_back(1./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2./sqrt(6.));Ihat4_.push_back(-4./sqrt(6.));
  // xi_c to lambda k*bar0
  incoming_.push_back(4132);outgoingB_.push_back(3122);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0373243);
  I1_.push_back(1./6.);I2_.push_back(1./6.);
  I3_.push_back(2./3.);I4_.push_back(-1./3.);I5_.push_back(1./6.);
  Ihat3_.push_back(-1./3.);Ihat4_.push_back(2./3.);
  // xi_c to sigma k*0
  incoming_.push_back(4132);outgoingB_.push_back(3212);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0135937);
  I1_.push_back(1./sqrt(12.));I2_.push_back(1./sqrt(12.));
  I3_.push_back(0.);I4_.push_back(-2./sqrt(12.));I5_.push_back(-1./sqrt(12.));
  Ihat3_.push_back(2./sqrt(12.));Ihat4_.push_back(-4./sqrt(12.));
  // xi_c to sigma k*-
  incoming_.push_back(4132);outgoingB_.push_back(3222);outgoingM_.push_back(-323);
  maxweight_.push_back(0.00877507);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(2./sqrt(6.));I5_.push_back(1./sqrt(6.));
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi0 rho0
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(113);
  maxweight_.push_back(0.0364247);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(2./sqrt(12.));I4_.push_back(-2./sqrt(12.));I5_.push_back(0.);
  Ihat3_.push_back(-2./sqrt(12.));Ihat4_.push_back(4./sqrt(12.));
  // xi_c0 to xi0 omega
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(223);
  maxweight_.push_back(0.134395);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back( 0.57735033);I4_.push_back(-0.577350303);I5_.push_back(0.);
  Ihat3_.push_back(0.577350303);Ihat4_.push_back(-1.15470061);
  // xi_c0 to xi0 phi
  incoming_.push_back(4132);outgoingB_.push_back(3322);outgoingM_.push_back(333);
  maxweight_.push_back(0.00676745);
  I1_.push_back(0.);I2_.push_back(0.);
  I3_.push_back(0.816496614);I4_.push_back(0.);I5_.push_back(0.408248312);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi- rho+
  incoming_.push_back(4132);outgoingB_.push_back(3312);outgoingM_.push_back(213);
  maxweight_.push_back(0.309733);
  I1_.push_back(-1./sqrt(6.));I2_.push_back(-1./sqrt(6.));
  I3_.push_back(-2./sqrt(6.));I4_.push_back(2./sqrt(6.));I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c to xi0 K*0
  incoming_.push_back(4332);outgoingB_.push_back(3322);outgoingM_.push_back(-313);
  maxweight_.push_back(0.019967);
  I1_.push_back(1.);I2_.push_back(-1.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(2.);Ihat4_.push_back(0.);
  // lambda_c to sigma*0 pi+
  incoming_.push_back(4122);outgoingB_.push_back(3214);outgoingM_.push_back(211);
  maxweight_.push_back(0.0254397);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ pi0
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(111);
  maxweight_.push_back(0.0256531);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ eta
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(221);
  maxweight_.push_back(0.0299659);
  I1_.push_back(0.);I2_.push_back( 0.569036);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ eta'
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(331);
  maxweight_.push_back(4.4186e-07);
  I1_.push_back(0.);I2_.push_back(-0.097631);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta+ kbar0
  incoming_.push_back(4122);outgoingB_.push_back(2214);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0215801);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta++ k-
  incoming_.push_back(4122);outgoingB_.push_back(2224);outgoingM_.push_back(-321);
  maxweight_.push_back(0.0649794);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to xi*0 k+
  incoming_.push_back(4122);outgoingB_.push_back(3324);outgoingM_.push_back(321);
  maxweight_.push_back(0.0198121);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*0 Kbar0
  incoming_.push_back(4132);outgoingB_.push_back(3214);outgoingM_.push_back(-311);
  maxweight_.push_back(0.0118739);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*+ Kbar-
  incoming_.push_back(4132);outgoingB_.push_back(3224);outgoingM_.push_back(-321);
  maxweight_.push_back(0.0240231);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 pi0
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(111);
  maxweight_.push_back(0.0173433);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 eta
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(221);
  maxweight_.push_back(0.000939099);
  I1_.push_back(0.);I2_.push_back(0.097631);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 eta'
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(331);
  maxweight_.push_back(1.31513e-05);
  I1_.push_back(0.);I2_.push_back(-0.569036);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*- pi+
  incoming_.push_back(4132);outgoingB_.push_back(3314);outgoingM_.push_back(211);
  maxweight_.push_back(0.0338039);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to omega- K+
  incoming_.push_back(4132);outgoingB_.push_back(3334);outgoingM_.push_back(321);
  maxweight_.push_back(0.00715116);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to xi*0 kbar0
  incoming_.push_back(4332);outgoingB_.push_back(3324);outgoingM_.push_back(-311);
  maxweight_.push_back(0.00850994);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to omega pi+
  incoming_.push_back(4332);outgoingB_.push_back(3334);outgoingM_.push_back(211);
  maxweight_.push_back(0.012698);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*0 rho+
  incoming_.push_back(4122);outgoingB_.push_back(3214);outgoingM_.push_back(213);
  maxweight_.push_back(0.0478329);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ rho0
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(113);
  maxweight_.push_back(0.0476973);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ omega
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(223);
  maxweight_.push_back(0.0260017);
  I1_.push_back(0.);I2_.push_back(1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to sigma*+ phi
  incoming_.push_back(4122);outgoingB_.push_back(3224);outgoingM_.push_back(333);
  maxweight_.push_back(3.55986e-07);
  I1_.push_back(0.);I2_.push_back(-sqrt(2.)/3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta+ kbar0
  incoming_.push_back(4122);outgoingB_.push_back(2214);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0988993);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to delta++ k-
  incoming_.push_back(4122);outgoingB_.push_back(2224);outgoingM_.push_back(-323);
  maxweight_.push_back(0.303435);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(6.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // lambda_c to xi*0 k+
  incoming_.push_back(4122);outgoingB_.push_back(3324);outgoingM_.push_back(323);
  maxweight_.push_back(9.70866e-05);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*0 Kbar0
  incoming_.push_back(4132);outgoingB_.push_back(3214);outgoingM_.push_back(-313);
  maxweight_.push_back(0.0281627);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to sigma*- K-
  incoming_.push_back(4132);outgoingB_.push_back(3224);outgoingM_.push_back(-323);
  maxweight_.push_back(0.0575937);
  I1_.push_back(0.);I2_.push_back(-2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 pi0
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(113);
  maxweight_.push_back(0.0480608);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 omega
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(223);
  maxweight_.push_back(0.0223044);
  I1_.push_back(0.);I2_.push_back(-1./3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*0 phi
  incoming_.push_back(4132);outgoingB_.push_back(3324);outgoingM_.push_back(333);
  maxweight_.push_back(3.187e-09);
  I1_.push_back(0.);I2_.push_back(-sqrt(2.)/3.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to xi*- rho+
  incoming_.push_back(4132);outgoingB_.push_back(3314);outgoingM_.push_back(213);
  maxweight_.push_back(0.0958052);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // xi_c0 to omega- K*+
  incoming_.push_back(4132);outgoingB_.push_back(3334);outgoingM_.push_back(323);
  maxweight_.push_back(0.000205851);
  I1_.push_back(0.);I2_.push_back(2./sqrt(18.));
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to xi*0 k*bar0
  incoming_.push_back(4332);outgoingB_.push_back(3324);outgoingM_.push_back(-313);
  maxweight_.push_back(0.050904);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // omega_c0 to omega rho+
  incoming_.push_back(4332);outgoingB_.push_back(3334);outgoingM_.push_back(213);
  maxweight_.push_back(0.0814021);
  I1_.push_back(-1./sqrt(3.));I2_.push_back(0.);
  I3_.push_back(0.);I4_.push_back(0.);I5_.push_back(0.);
  Ihat3_.push_back(0.);Ihat4_.push_back(0.);
  // initial size of the vectors
  initsize_=incoming_.size();
  // intermediates
  generateIntermediates(false);
}

void KornerKramerCharmDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

void KornerKramerCharmDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the vectors have the same size
  unsigned int isize=incoming_.size();
  if(isize!=I1_.size()||isize!=I2_.size()||isize!=I3_.size()||isize!=I4_.size()||
     isize!=I5_.size()||isize!=Ihat3_.size()||isize!=Ihat4_.size()||
     isize!=incoming_.size()||isize!=outgoingB_.size()||isize!=outgoingM_.size()||
     isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in"
			  << " KornerKramerCharmDecayer::doinit()" 
			  << Exception::abortnow;
  // compute the various coefficients
  Energy m1,m2,m3,fmes(ZERO); 
  Energy2 P1P2,Qplus,gmes(ZERO);
  double Fnonfact,A3,B3,Ffact[2];
  Energy H2,H3,A2,B2;
  Energy2 A,B;
  int mspin,bspin;
  double chi,
    chiplus( 0.5*(cplus_*(1.+oneNC_)+cminus_*(1.-oneNC_))),
    chiminus(0.5*(cplus_*(1.+oneNC_)-cminus_*(1.-oneNC_))); 
  InvEnergy2 pre(SM().fermiConstant()*sqrt(SM().CKM(1,1)*SM().CKM(0,0)/2.)),mform2[2];
  // testing only
//   pre = SM().fermiConstant()*0.974/sqrt(2.);
  for(unsigned int ix=0;ix<isize;++ix) {
    // get the mass of the particles
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoingB_[ix]),
    		     getParticleData(outgoingM_[ix])};
    PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(in,out,maxweight_[ix]));
    addMode(mode);
    m1=in->mass();
    m2=out[0]->mass();
    m3=out[1]->mass();
    // dot product P1P2
    P1P2 = 0.5*(m1*m1+m2*m2-m3*m3);
    // formfactor for the non-factorizating diagrams and coefficients
    Fnonfact = 1./(1.-m3*m3/mscminus_/mscminus_)*
      (1.-(m1-m2)*(m1-m2)/mscminus_/mscminus_);
    Fnonfact *= Fnonfact*Fnonfact;
    H2 = H2_*Fnonfact;
    H3 = H3_*Fnonfact;
    if(abs(outgoingM_[ix])==211) {
      fmes=fpi_;
      chi=chiplus;
      mform2[0]=1./(mscminus_*mscminus_);
      mform2[1]=1./(mscplus_ *mscplus_);
    }
    else if(abs(outgoingM_[ix])==311) {
      fmes=fk_;
      chi=chiminus;
      mform2[0]=1./(mdcminus_*mdcminus_);
      mform2[1]=1./(mdcplus_ *mdcplus_);
    }
    else if(abs(outgoingM_[ix])==213) {
      gmes = frho_*m3*m3;
      chi=chiplus;
      mform2[0]=1./(mscminus_*mscminus_);
      mform2[1]=1./(mscplus_ *mscplus_);
    }
    else if(abs(outgoingM_[ix])==313) {
      gmes = fKstar_*m3*m3;
      chi=chiminus;
      mform2[0]=1./(mdcminus_*mdcminus_);
      mform2[1]=1./(mdcplus_ *mdcplus_);
    }
    else {
      fmes=ZERO;
      chi=0.;
      mform2[0]=1./(mscminus_*mscminus_);
      mform2[1]=1./(mscplus_ *mscplus_);
      gmes=ZERO;
    }
    // form factor for the factorising diagrams
    for(unsigned int iy=0;iy<2;++iy) {
      Ffact[iy]  = 1./(1.-m3*m3*mform2[iy])*(1.-(m1-m2)*(m1-m2)*mform2[iy]);
      Ffact[iy] *= Ffact[iy]*Ffact[iy];
    }
    // invariants
    Qplus  = (m1+m2)*(m1+m2)-m3*m3;
    // decide which type of decay
    mspin=getParticleData(outgoingM_[ix])->iSpin();
    bspin=getParticleData(outgoingB_[ix])->iSpin();
    if(bspin==2) {
      if(mspin==1) {
	// non-factorizing piece of the diagrams
	A = -0.75*H2/m1/m2*cminus_*( (m1*P1P2-m1*m2*(m2+m3))*I3_[ix]
				     -(m2*P1P2-m1*m2*(m1+m3))*Ihat3_[ix]);
	B = 0.25/m1/m2*cminus_*( 0.5*H2*Qplus*( m1*(   I3_[ix]+2.*    I4_[ix])
						+m2*(Ihat3_[ix]+2.*Ihat4_[ix]))
				 +H3*m1*m2*(m1+m2+m3)*12.*I5_[ix]);
	// the factorizing piece
	A +=0.25/m1/m2*chi*fmes*Qplus*(m1-m2)*(2.*I1_[ix]+I2_[ix])*Ffact[0];
	B +=0.25/m1/m2*chi*fmes*Qplus*(m1+m2)/3.*(4.*I1_[ix]+5.*I2_[ix])*Ffact[1];
	// add to vectors
	A1_.push_back(A*pre);B1_.push_back(B*pre);
	A2_.push_back(ZERO);B2_.push_back(ZERO);
	A3_.push_back(ZERO);B3_.push_back(ZERO);
      }
      else if(mspin==3) {
	// non-factorizing terms
	A  = -0.25*H2/m1/m2*cminus_*
	  ( (m1*P1P2-m1*m2*(m2+m3))*(   I3_[ix]+2.*   I4_[ix])
	    +(m2*P1P2-m1*m2*(m1+m3))*(Ihat3_[ix]+2.*Ihat4_[ix]));
	A2 = -0.25*H2/m1/m2*cminus_*(m1+m2+m3)*( m1*(   I3_[ix]+2.*   I4_[ix])
						 -m2*(Ihat3_[ix]+2.*Ihat4_[ix]));
	B  = +0.25/m1/m2*cminus_*(0.5*H2*Qplus*( m1*(   I3_[ix]+2.*   I4_[ix])
						 +m2*(Ihat3_[ix]+2.*Ihat4_[ix]))
				  +H3*m1*m2*(m1+m2+m3)*12.*I5_[ix]);
	B2 = -0.25/m1/m2*cminus_*
	  ( H2*( m1*(m1+m2)*(   I3_[ix]+2.*   I4_[ix])-3.*m1*m3*   I3_[ix]
		 +m2*(m1+m2)*(Ihat3_[ix]+2.*Ihat4_[ix])-3.*m2*m3*Ihat3_[ix])
	    +H3*m1*m2*24.*I5_[ix]);
	// the factorizing piece
	A  += -0.25/m1/m2*chi*gmes*Qplus/3.*(4.*I1_[ix]+5.*I2_[ix])*Ffact[1];
	B  +=  0.25/m1/m2*chi*gmes*Qplus/3.*(4.*I1_[ix]+5.*I2_[ix])*Ffact[0];
	B2 +=  0.25/m1/m2*chi*gmes*(m1+m2)*4./3.*(I1_[ix]-I2_[ix])*Ffact[0];
	// add to vectors
	A1_.push_back(A*pre);B1_.push_back(B*pre);
	A2_.push_back(A2*pre);B2_.push_back(B2*pre);
	A3_.push_back(ZERO);B3_.push_back(ZERO);
      }
      else
	throw InitException() << "Invalid outgoing meson spin in"
			      << " KornerKramerCharmDecayer::doinit()" 
			      << Exception::abortnow;
    }
    else if(bspin==4) {
      if(mspin==1) {
	// first the non-factorizing piece
	B2  = -1.5*cminus_*H2*I2_[ix];
	// then the factorizing piece
	B2  += chi*fmes*(1.+m2/m1)*I1_[ix]*Ffact[1];
	// add to vectors
	// make the coupling dimensionless
	A1_.push_back(0.);B1_.push_back(0.);
	A2_.push_back(ZERO);B2_.push_back(B2*pre);
	A3_.push_back(ZERO);B3_.push_back(ZERO);
      }
      else if(mspin==3) {
	InvEnergy norm(0.75/m1/m2*cminus_*H2*I2_[ix]);
	// first the non-factorizing piece
	A  = -norm*m1*(P1P2-m2*(m2+m3))*2.;
	A2 = ZERO;
	A3 =  norm*m1*2.;
	B  = -norm*m1*Qplus;
	B2 = -norm*m1*m2*2.;
	B3 =  norm*m1*2.;
	// then the factorizing piece
	double norm2 = 0.5/m1/m2*chi*gmes*I1_[ix];
	A  += norm2*Qplus*Ffact[1];
      //A2 += ZERO;
	A3 +=-norm2*2.*Ffact[1];
	B  += norm2*Qplus*Ffact[0];
	B2 += norm2*m2*2.*Ffact[0];
	B3 +=-norm2*2.*Ffact[0];
	// add to vectors
	A1_.push_back( A*pre);B1_.push_back( B*pre);
	A2_.push_back(2.*A2*pre);B2_.push_back(2.*B2*pre);
	A3_.push_back(A3*pre);B3_.push_back(B3*pre);
      }
      else
	throw InitException() << "Invalid outgoing meson spin in"
			      << " KornerKramerCharmDecayer::doinit()" 
			      << Exception::abortnow;
    }
    else
      throw InitException() << "Invalid outgoing baryon spin in"
			    << " KornerKramerCharmDecayer::doinit()" 
			    << Exception::abortnow;
  }
}

int KornerKramerCharmDecayer::modeNumber(bool & cc,tcPDPtr parent,
					 const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do {
    if(id0==incoming_[ix]) {
      if((id1==outgoingB_[ix]&&id2==outgoingM_[ix])||
	 (id2==outgoingB_[ix]&&id1==outgoingM_[ix])) imode=ix;
    }
    else if(id0==-incoming_[ix]) {
      if((id1==-outgoingB_[ix]&&id2==-outgoingM_[ix])||
	 (id2==-outgoingB_[ix]&&id1==-outgoingM_[ix])) imode=ix;
      if(((id1==-outgoingB_[ix]&&id2==outgoingM_[ix])||
	  (id2==-outgoingB_[ix]&&id1==outgoingM_[ix]))&&
	 (outgoingM_[ix]==111||outgoingM_[ix]==221||outgoingM_[ix]==331||
	  outgoingM_[ix]==113||outgoingM_[ix]==223||outgoingM_[ix]==333))
	imode=ix;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  // charge conjugation if needed
  cc = id0<0;
  // return mode number
  return imode;
}

void KornerKramerCharmDecayer::persistentOutput(PersistentOStream & os) const {
  os << oneNC_ << ounit(fpi_,GeV) << ounit(fk_,GeV) << frho_ 
     << fKstar_ << ounit(mdcplus_,GeV) << ounit(mdcminus_,GeV) << ounit(mscplus_,GeV) 
     << ounit(mscminus_,GeV) << cplus_ << cminus_ << ounit(H2_,GeV) << ounit(H3_,GeV) 
     << I1_ << I2_ << I3_ << I4_ << I5_ << Ihat3_ << Ihat4_ << incoming_ << outgoingB_ 
     << outgoingM_ << maxweight_ << A1_ << ounit(A2_,1./GeV) << ounit(A3_,1./GeV2) 
     << B1_ << ounit(B2_,1./GeV) << ounit(B3_,1./GeV2) << initsize_;
}

void KornerKramerCharmDecayer::persistentInput(PersistentIStream & is, int) {
  is >> oneNC_ >> iunit(fpi_,GeV) >> iunit(fk_,GeV) >> frho_ 
     >> fKstar_ >> iunit(mdcplus_,GeV) >> iunit(mdcminus_,GeV) >> iunit(mscplus_,GeV) 
     >> iunit(mscminus_,GeV) >> cplus_ >> cminus_ >> iunit(H2_,GeV) >> iunit(H3_,GeV) 
     >> I1_ >> I2_ >> I3_ >> I4_ >> I5_ >> Ihat3_ >> Ihat4_ >> incoming_ >> outgoingB_ 
     >> outgoingM_ >> maxweight_ >> A1_ >> iunit(A2_,1./GeV) >> iunit(A3_,1./GeV2) 
     >> B1_ >> iunit(B2_,1./GeV) >> iunit(B3_,1./GeV2) >> initsize_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<KornerKramerCharmDecayer,Baryon1MesonDecayerBase>
describeHerwigKornerKramerCharmDecayer("Herwig::KornerKramerCharmDecayer", "HwBaryonDecay.so");

void KornerKramerCharmDecayer::Init() {

  static ClassDocumentation<KornerKramerCharmDecayer> documentation
    ("The KornerKramerCharmDecayer class implements the"
     " non-leptonic weak decay of charm baryons using the results of Z.Phys.C55,659.",
     "The non-leptonic charm decays were simulated using the KornerKramerCharmDecayer"
     "class which implements the model of \\cite{Korner:1992wi}.",
     "\\bibitem{Korner:1992wi}\n"
     "J.~G.~Korner and M.~Kramer,\n"
     "Z.\\ Phys.\\  C {\\bf 55} (1992) 659.\n"
     "%%CITATION = ZEPYA,C55,659;%%\n");

  static Parameter<KornerKramerCharmDecayer,double> interfaceOneOverNc
    ("OneOverNc",
     "One divided by the number of colours, the default is to take N_c to infinity.",
     &KornerKramerCharmDecayer::oneNC_, 0.0, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceFpi
    ("Fpi",
     "The decay constant for the pi meson.",
     &KornerKramerCharmDecayer::fpi_, MeV, 131.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceFK
    ("FK",
     "The decay constant for the K meson",
     &KornerKramerCharmDecayer::fk_, MeV, 160.6*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceFrho
    ("Frho",
     "The decay constant for the rho meson",
     &KornerKramerCharmDecayer::frho_, 0.272, 0.0, 1.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfacefKstar
    ("fKstar",
     "The decay constant for the K star meson",
     &KornerKramerCharmDecayer::fKstar_, 0.238, 0.0, 1.0,
     false, false, false);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMdcplus
    ("Mdcplus",
     "The mass of the 1+ dc meson for the form-factors",
     &KornerKramerCharmDecayer::mdcplus_, GeV, 2.42*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMscplus
    ("Mscplus",
     "The mass of the 1+ sc meson for the form-factors",
     &KornerKramerCharmDecayer::mscplus_, GeV, 2.54*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMdcminus
    ("Mdcminus",
     "The mass of the 1- dc meson for the form-factors",
     &KornerKramerCharmDecayer::mdcminus_, GeV, 2.01*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceMscminus
    ("Mscminus",
     "The mass of the 1- sc meson for the form-factors",
     &KornerKramerCharmDecayer::mscminus_, GeV, 2.11*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceCplus
    ("Cplus",
     "The c+ perturvative coefficient",
     &KornerKramerCharmDecayer::cplus_, 0.73, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,double> interfaceCminus
    ("Cminus",
     "The c+ perturvative coefficient",
     &KornerKramerCharmDecayer::cminus_, 1.90, 0.0, 10.0,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceH2
    ("H2",
     "The H2 parameter",
     &KornerKramerCharmDecayer::H2_, GeV, 0.119*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<KornerKramerCharmDecayer,Energy> interfaceH3
    ("H3",
     "The H3 parameter",
     &KornerKramerCharmDecayer::H3_, GeV,-0.011*GeV,-10.0*GeV, 10.0*GeV,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI1
    ("I1",
     "The I_1 invariant for the decay modes",
     &KornerKramerCharmDecayer::I1_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI2
    ("I2",
     "The I_2 invariant for the decay modes",
     &KornerKramerCharmDecayer::I2_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI3
    ("I3",
     "The I_3 invariant for the decay modes",
     &KornerKramerCharmDecayer::I3_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI4
    ("I4",
     "The I_4 invariant for the decay modes",
     &KornerKramerCharmDecayer::I4_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceI5
    ("I5",
     "The I_5 invariant for the decay modes",
     &KornerKramerCharmDecayer::I5_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceIhat3
    ("Ihat3",
     "The Ihat_3 invariant for the decay modes",
     &KornerKramerCharmDecayer::Ihat3_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceIhat4
    ("Ihat4",
     "The Ihat_4 invariant for the decay modes",
     &KornerKramerCharmDecayer::Ihat4_, -1, 0., -10.0, 10.0,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDG code of the incoming baryon",
     &KornerKramerCharmDecayer::incoming_, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceOutgoingB
    ("OutgoingB",
     "The PDG code of the outgoing baryon",
     &KornerKramerCharmDecayer::outgoingB_, -1, 0, 0, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,int> interfaceOutgoingM
    ("OutgoingM",
     "The PDG code of the outgoing meson",
     &KornerKramerCharmDecayer::outgoingM_, -1, 0, -1000000, 1000000,
     false, false, true);

  static ParVector<KornerKramerCharmDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &KornerKramerCharmDecayer::maxweight_,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void KornerKramerCharmDecayer::
 halfHalfScalarCoupling(int imode,Energy,Energy,Energy,
			Complex& A,Complex&B) const {
  useMe();
  A = A1_[imode];
  B = B1_[imode];
}

// couplings for spin-1/2 to spin-1/2 spin-1
void KornerKramerCharmDecayer::halfHalfVectorCoupling(int imode,
						      Energy m0,Energy m1,Energy,
						      Complex& A1,Complex& A2,
						      Complex& B1,Complex& B2) const {
  useMe();
  // conventions changed so that A is the coefficient of the 
  // non-gamma_5 piece in the base class
  A1 = B1_[imode];
  B1 = A1_[imode];
  A2 = B2_[imode]*(m0+m1);
  B2 = A2_[imode]*(m0+m1);
}

// couplings for spin-1/2 to spin-3/2 spin-0
void KornerKramerCharmDecayer::halfThreeHalfScalarCoupling(int imode, Energy m0,
							   Energy m1,Energy,Complex& A,
							   Complex& B) const {
  useMe();
  // conventions changed so that A is the coefficient of the 
  // non-gamma_5 piece in the base class
  A = B2_[imode]*(m0+m1);
  B = A2_[imode]*(m0+m1);
}

// couplings for spin-1/2 to spin-3/2 spin-1
void KornerKramerCharmDecayer::
halfThreeHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex& A1, Complex& A2, Complex& A3,
			    Complex& B1, Complex& B2, Complex& B3) const {
  useMe();
  A1 = A1_[imode];
  A2 = A2_[imode]*(m0+m1);
  A3 = A3_[imode]*(m0+m1)*(m0+m1);
  B1 = B1_[imode];
  B2 = B2_[imode]*(m0+m1);
  B3 = B3_[imode]*(m0+m1)*(m0+m1);
}
 
void KornerKramerCharmDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":OneOverNc " <<  oneNC_ << "\n";
  output << "newdef " << name() << ":Fpi " << fpi_/MeV << "\n";
  output << "newdef " << name() << ":FK " << fk_/MeV  << "\n";
  output << "newdef " << name() << ":Frho " << frho_ << "\n";
  output << "newdef " << name() << ":fKstar " << fKstar_ << "\n";
  output << "newdef " << name() << ":Mdcplus " << mdcplus_/GeV << "\n";
  output << "newdef " << name() << ":Mscplus " << mscplus_/GeV << "\n";
  output << "newdef " << name() << ":Mdcminus " << mdcminus_/GeV << "\n";
  output << "newdef " << name() << ":Mscminus " << mscminus_/GeV << "\n";
  output << "newdef " << name() << ":Cplus " << cplus_ << "\n";
  output << "newdef " << name() << ":Cminus " << cminus_ << "\n";
  output << "newdef " << name() << ":H2 " << H2_/GeV << "\n";
  output << "newdef " << name() << ":H3 " << H3_/GeV << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    if(ix<initsize_) {
      output << "newdef " << name() << ":I1 " 
	     << ix << " " << I1_[ix] << "\n";
      output << "newdef " << name() << ":I2 " 
	     << ix << " " << I2_[ix] << "\n";
      output << "newdef " << name() << ":I3 " 
	     << ix << " " << I3_[ix] << "\n";
      output << "newdef " << name() << ":I4 " 
	     << ix << " " << I4_[ix] << "\n";
      output << "newdef " << name() << ":I5 " 
	     << ix << " " << I5_[ix] << "\n";
      output << "newdef " << name() << ":Ihat3 " 
	     << ix << " " << Ihat3_[ix] << "\n";
      output << "newdef " << name() << ":Ihat4 " 
	     << ix << " " << Ihat4_[ix] << "\n";
      output << "newdef " << name() << ":Incoming " 
	     << ix << " " << incoming_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingB " 
	     << ix << " " << outgoingB_[ix] << "\n";
      output << "newdef " << name() << ":OutgoingM " 
	     << ix << " " << outgoingM_[ix] << "\n";
      output << "newdef " << name() << ":MaxWeight " 
	     << ix << " " << maxweight_[ix] << "\n";
    }
    else {
      output << "insert " << name() << ":I1 " 
	     << ix << " " << I1_[ix] << "\n";
      output << "insert " << name() << ":I2 " 
	     << ix << " " << I2_[ix] << "\n";
      output << "insert " << name() << ":I3 " 
	     << ix << " " << I3_[ix] << "\n";
      output << "insert " << name() << ":I4 " 
	     << ix << " " << I4_[ix] << "\n";
      output << "insert " << name() << ":I5 " 
	     << ix << " " << I5_[ix] << "\n";
      output << "insert " << name() << ":Ihat3 " 
	     << ix << " " << Ihat3_[ix] << "\n";
      output << "insert " << name() << ":Ihat4 " 
	     << ix << " " << Ihat4_[ix] << "\n";
      output << "insert " << name() << ":Incoming " 
	     << ix << " " << incoming_[ix] << "\n";
      output << "insert " << name() << ":OutgoingB " 
	     << ix << " " << outgoingB_[ix] << "\n";
      output << "insert " << name() << ":OutgoingM " 
		 << ix << " " << outgoingM_[ix] << "\n";
      output << "insert " << name() << ":MaxWeight " 
		 << ix << " " << maxweight_[ix] << "\n";
    }
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./OmegaXiStarPionDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OmegaXiStarPionDecayer class.
//

#include "OmegaXiStarPionDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

OmegaXiStarPionDecayer::OmegaXiStarPionDecayer()  {
  // the ids of the particles
  idin_  = 3334;
  idout_ = 3324;
  // the couplings from the paper
  Acomm_ =  20.91e-8;
  AP_    =  -9.20e-8;
  AS_    =  -6.32e-8;
  BP_    =  230.1e-8;
  BS_    = -100.8e-8;
  // maximum weight for the decay
  wgtmax_=0.0032;
  // intermediates
  generateIntermediates(false);
}

void OmegaXiStarPionDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the phase space
  tPDPtr    in  =  getParticleData(idin_);
  tPDVector out = {getParticleData(idout_),
  		   getParticleData(-211)};
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,wgtmax_));
  addMode(mode);
}

void OmegaXiStarPionDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) wgtmax_=mode(0)->maxWeight();
}

int OmegaXiStarPionDecayer::modeNumber(bool & cc,tcPDPtr parent,
				       const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  if(id0==idin_) {
    if((id1==idout_&&id2==-211)||
       (id2==idout_&&id1==-211)) imode=0;
  }
  else if(id0==-idin_) {
    if((id1==-idout_&&id2==211)||
       (id2==-idout_&&id1==211)) imode=0;
  }
  // charge conjugation
  cc=id0<0;
  // return the answer
  return imode;
}

void OmegaXiStarPionDecayer::persistentOutput(PersistentOStream & os) const {
  os << Acomm_ << AP_ << AS_ << BP_ << BS_ << idin_ << idout_ <<  wgtmax_;
}

void OmegaXiStarPionDecayer::persistentInput(PersistentIStream & is, int) {
  is >> Acomm_ >> AP_ >> AS_ >> BP_ >> BS_ >> idin_ >> idout_ >>  wgtmax_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<OmegaXiStarPionDecayer,Baryon1MesonDecayerBase>
describeHerwigOmegaXiStarPionDecayer("Herwig::OmegaXiStarPionDecayer", "HwBaryonDecay.so");

void OmegaXiStarPionDecayer::Init() {
  
  static ClassDocumentation<OmegaXiStarPionDecayer> documentation
    ("The OmegaXiStarPionDecayer class performs the weak decay"
     " of the Omega to Xi*0 and pi-",
     "The decay of the $\\Omega^-$ to $\\Xi^{*0}\\pi^-$ was simulated"
     " using the model of \\cite{Duplancic:2004dy}.",
     "\\bibitem{Duplancic:2004dy}\n"
     "G.~Duplancic, H.~Pasagic and J.~Trampetic,\n"
     "Phys.\\ Rev.\\  D {\\bf 70} (2004) 077506 [arXiv:hep-ph/0405162].\n"
     "%%CITATION = PHRVA,D70,077506;%%\n");

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAcomm
    ("Acomm",
     "The Acomm coupling for the decay",
     &OmegaXiStarPionDecayer::Acomm_, 20.91e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAP
    ("AP",
     "The A_P coupling for the decay",
     &OmegaXiStarPionDecayer::AP_, -9.20e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceAS
    ("AS",
     "The A_S coupling for the decay",
     &OmegaXiStarPionDecayer::AS_, -6.32e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceBP
    ("BP",
     "The B_P coupling for the decay",
     &OmegaXiStarPionDecayer::BP_, 230.1e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceBS
    ("BS",
     "The B_S coupling for the decay",
     &OmegaXiStarPionDecayer::BS_, -100.8e-8, -1.e-5, 1.e-5,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the decay",
     &OmegaXiStarPionDecayer::wgtmax_, 0.0032, 0., 100.,
     false, false, false);

  static Parameter<OmegaXiStarPionDecayer,int> interfaceIncoming
    ("Incoming",
     "The PDF code for the incoming baryon",
     &OmegaXiStarPionDecayer::idin_, 3334, 0, 1000000,
     false, false, true);

  static Parameter<OmegaXiStarPionDecayer,int> interfaceOutgoing
    ("Outgoing",
     "The PDF code for the outgoing baryon",
     &OmegaXiStarPionDecayer::idout_, 3324, 0, 1000000,
     false, false, true);
}

// couplings for spin-3/2 to spin-3/2 spin-0
void OmegaXiStarPionDecayer::
threeHalfThreeHalfScalarCoupling(int,Energy,Energy,Energy,
				 Complex&A1,Complex&A2,Complex&B1,Complex&B2) const {
  useMe();
  A2=0.;
  B2=0.;
  A1=Acomm_+AP_+AS_;
  B1=BP_+BS_;
}

void OmegaXiStarPionDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Acomm " << Acomm_ << "\n";
  output << "newdef " << name() << ":AP " << AP_ << "\n";
  output << "newdef " << name() << ":AS " << AS_ << "\n";
  output << "newdef " << name() << ":BP " << BP_ << "\n";
  output << "newdef " << name() << ":BS " << BS_ << "\n";
  output << "newdef " << name() << ":MaximumWeight " << wgtmax_ << "\n";
  output << "newdef " << name() << ":Incoming " << idin_ << "\n";
  output << "newdef " << name() << ":Outgoing " << idout_ << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SemiLeptonicBaryonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiLeptonicBaryonDecayer class.
//

#include "SemiLeptonicBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/LorentzPolarizationVector.h"
#include "ThePEG/Helicity/FermionSpinInfo.h"
#include "ThePEG/Helicity/RSFermionSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

SemiLeptonicBaryonDecayer::SemiLeptonicBaryonDecayer() {
  // intermediates
  generateIntermediates(true);
}

void SemiLeptonicBaryonDecayer::doinitrun() {
  _current->initrun();
  _form->initrun();
  DecayIntegrator::doinitrun();
  if(initialize()) {
    _maxwgt.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxwgt.push_back(mode(ix)->maxWeight());
  }
}

void SemiLeptonicBaryonDecayer::doinit() {
  DecayIntegrator::doinit();
  // make sure the current got initialised
  _current->init();
  // and the form factors
  _form->init();
  // the channels
  _modemap.clear();
  for(unsigned int ix=0;ix<_form->numberOfFactors();++ix) {
    int id0(0),id1(0);
    // get the external particles for this mode
    _form->particleID(ix,id0,id1);
    int inspin,spect1,spect2,inquark,outquark,outspin;
    _form->formFactorInfo(ix,inspin,outspin,spect1,spect2,inquark,outquark);
    // incoming  and outgoing particles
    tPDPtr in  = getParticleData(id0);
    tPDPtr out = getParticleData(id1);
    // charge of the W
    int Wcharge =(in->iCharge()-out->iCharge());
    Energy min = in->mass()+in->widthUpCut()
      -out->mass()+out->widthLoCut();
    _modemap.push_back(numberModes());
    for(unsigned int iy=0;iy<_current->numberOfModes();++iy) {
      int iq(0),ia(0);
      tPDVector outV = {out};
      tPDVector ptemp=_current->particles(Wcharge,iy,iq,ia);
      outV.insert(std::end(outV), std::begin(ptemp), std::end(ptemp));
      // create the mode
      PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,outV,1.));
      // create the first piece of the channel
      PhaseSpaceChannel channel((PhaseSpaceChannel(mode),0,1));
      // and the rest
      bool done = _current->createMode(Wcharge,tcPDPtr(),FlavourInfo(),
				       iy,mode,1,0,channel,min);
      // check the result
      if(done&&abs(Wcharge)==3&&inspin==2&&(outspin==2||outspin==4)) {
	// the maximum weight
	double maxweight = _maxwgt.size()>numberModes() ? _maxwgt[numberModes()] : 2.;
	mode->maxWeight(maxweight);
	addMode(mode);
      }
    }
  }
}

bool SemiLeptonicBaryonDecayer::accept(tcPDPtr parent,
				       const tPDVector & children) const {
  // find the non-lepton
  int ibar(0),idtemp,idin(parent->id());
  vector<int> idother; bool dummy;
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) ibar=idtemp;
    else idother.push_back(idtemp);
  }
  // check that the form factor exists
  if(_form->formFactorNumber(idin,ibar,dummy)<0) return false;
  // and the current
  return _current->accept(idother);
}

int SemiLeptonicBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  // find the ids of the particles for the decay current
  tPDVector::const_iterator pit  = children.begin();
  tPDVector::const_iterator pend = children.end();
  int idtemp,ibar(0),idin(parent->id());
  vector<int> idother;
  cc=false;
  for( ; pit!=pend;++pit) {
    idtemp=(**pit).id();
    if(abs(idtemp)>16) ibar=idtemp;
    else idother.push_back(idtemp);
  }
  return _modemap[_form->formFactorNumber(idin,ibar,cc)]
    +_current->decayMode(idother);
}


void SemiLeptonicBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << _current << _form << _maxwgt << _modemap;
}

void SemiLeptonicBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _current >> _form >> _maxwgt >> _modemap;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SemiLeptonicBaryonDecayer,DecayIntegrator>
describeHerwigSemiLeptonicBaryonDecayer("Herwig::SemiLeptonicBaryonDecayer", "HwBaryonDecay.so");

void SemiLeptonicBaryonDecayer::Init() {

  static ClassDocumentation<SemiLeptonicBaryonDecayer> documentation
    ("The SemiLeptonicBaryonDecayer class is designed for"
     " the semi-leptonic decay of the baryons.");

  static Reference<SemiLeptonicBaryonDecayer,LeptonNeutrinoCurrent> interfaceCurrent
    ("Current",
     "The current for the leptons produced in the decay.",
     &SemiLeptonicBaryonDecayer::_current, true, true, true, false, false);


  static Reference<SemiLeptonicBaryonDecayer,BaryonFormFactor> interfaceFormFactor
    ("FormFactor",
     "The form factor",
     &SemiLeptonicBaryonDecayer::_form, true, true, true, false, false);

  static ParVector<SemiLeptonicBaryonDecayer,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weights for the decays",
     &SemiLeptonicBaryonDecayer::_maxwgt,
     0, 0, 0, 0, 10000, false, false, true);

}

void SemiLeptonicBaryonDecayer::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // for the decaying particle
  if(part.id()>0) {
    SpinorWaveFunction::
      constructSpinInfo(_inHalf,const_ptr_cast<tPPtr>(&part),incoming,true);
  }
  else {
    SpinorBarWaveFunction::
      constructSpinInfo(_inHalfBar,const_ptr_cast<tPPtr>(&part),incoming,true);
  }
  // decay product
  // spin 1/2
  if(decay[0]->dataPtr()->iSpin()==PDT::Spin1Half) {
    if(part.id()>0) {
      SpinorBarWaveFunction::constructSpinInfo(_inHalfBar,decay[0],outgoing,true);
    }
    else {
      SpinorWaveFunction::constructSpinInfo(_inHalf,decay[0],outgoing,true);
    }
  }
  // spin 3/2
  else if(decay[0]->dataPtr()->iSpin()==PDT::Spin3Half) {
    if(part.id()>0) {
      RSSpinorBarWaveFunction::constructSpinInfo(_inThreeHalfBar,
  						 decay[0],outgoing,true);
      
    }
    else {
      RSSpinorWaveFunction::constructSpinInfo(_inThreeHalf,
  					      decay[0],outgoing,true);
    }
  }
  else
    assert(false);
  // and the stuff from the current
  _current->constructSpinInfo(ParticleVector(decay.begin()+1,decay.end()));
}

// combine the currents and form-factors to give the matrix element
double SemiLeptonicBaryonDecayer::me2(const int , const Particle & part,
				      const tPDVector & outgoing,
				      const vector<Lorentz5Momentum> & momenta,
				      MEOption meopt) const {
  assert(part.dataPtr()->iSpin()==PDT::Spin1Half);
  double me(0.);
  if(outgoing[0]->iSpin()==PDT::Spin1Half)
    me =      halfHalf(part,outgoing,momenta,meopt);
  else if(outgoing[0]->iSpin()==PDT::Spin3Half)
    me = halfThreeHalf(part,outgoing,momenta,meopt);
  else
    assert(false);
  return me;
}

// matrix element for a 1/2 -> 1/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfHalf(const Particle & part,
					   const tPDVector & outgoing,
					   const vector<Lorentz5Momentum> & momenta,
					   MEOption meopt) const {
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    // work out the mapping for the lepton vector
    _constants.resize(outgoing.size()+1);
    _ispin.resize(outgoing.size());
    int itemp(1);
    _ibar=0;
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      _ispin[ix]=outgoing[ix]->iSpin();
      if(abs(outgoing[ix]->id())<=16) {
  	itemp*=_ispin[ix];
  	_constants[ix]=itemp;
      }
      else _ibar=ix;
    }
    _constants[outgoing.size()]=1;
    _constants[_ibar]=_constants[_ibar+1];
  }
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,_ispin)));
  // spinors for the decay product
  if(part.id()>0) {
    _inHalfBar.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel)
      _inHalfBar[ihel] = HelicityFunctions::dimensionedSpinorBar(-momenta[0],ihel,Helicity::outgoing);
  }
  else {
    _inHalf.resize(2);
    for(unsigned int ihel=0;ihel<2;++ihel)
      _inHalf[ihel] = HelicityFunctions::dimensionedSpinor   (-momenta[0],ihel,Helicity::outgoing);
  }
  // get the information on the form-factor
  int spinin(0),spinout(0),spect1,spect2,inquark,outquark;
  int id0(part.id()),id1(outgoing[0]->id());
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the form factors
  Complex f1v,f2v,f3v,f1a,f2a,f3a;
  _form->SpinHalfSpinHalfFormFactor(q2,iloc,id0,id1,m0,m1,
  				    f1v,f2v,f3v,f1a,f2a,f3a,FlavourInfo());
  // calculate the hadronic current for the decay
  vector<LorentzPolarizationVectorE> hadron(4);
  Complex left  =f1v-f1a-f2v-double((m0-m1)/(m0+m1))*f2a;
  Complex right =f1v+f1a-f2v+double((m0-m1)/(m0+m1))*f2a;
  LorentzPolarizationVectorE vtemp;
  for(unsigned int ix=0;ix<2;++ix) {
    for(unsigned int iy=0;iy<2;++iy) {
      vtemp = _inHalf[ix].generalCurrent(_inHalfBar[iy],left,right);
      complex<Energy> vspin = _inHalf[ix].scalar(_inHalfBar[iy]);
      complex<Energy> aspin = _inHalf[ix].pseudoScalar(_inHalfBar[iy]);
      // the momentum like pieces
      if(part.id()>0) {
  	vtemp+= (f2v*vspin+f2a*aspin)/(m0+m1)*sum;
  	vtemp+= (f3v*vspin+f3a*aspin)/(m0+m1)*q;
      }
      else {
  	vtemp-= (f2v*vspin-f2a*aspin)/(m0+m1)*sum;
  	vtemp+= (f3v*vspin-f3a*aspin)/(m0+m1)*q;
      }
      if(part.id()>0) hadron[2*ix+iy]=vtemp;
      else            hadron[2*iy+ix]=vtemp;
    }
  }
  // construct the lepton current
  Energy scale;
  int mode((abs(outgoing[1]->id())-11)/2);
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(tcPDPtr(),FlavourInfo(),
			     mode,-1,scale,tPDVector(outgoing.begin()+1,outgoing.end()),
			     vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),
			     meopt));
  // matrix element
  vector<unsigned int> ihel(outgoing.size()+1);
  unsigned int mhel,ix,lhel;
  for(mhel=0;mhel<hadron.size();++mhel) {
    ihel[0      ]=mhel/spinout;
    ihel[_ibar+1]=mhel%spinout;
    for(lhel=0;lhel<lepton.size();++lhel) {
      // map the index for the leptons to a helicity state
      for(ix=outgoing.size();ix>0;--ix) {
  	if(ix-1!=_ibar) ihel[ix]=(lhel%_constants[ix-1])/_constants[ix];
      }
      (*ME())(ihel) = Complex(lepton[lhel].dot(hadron[mhel])*SM().fermiConstant());
    }
  }
  // ckm factor
  double ckm(1.);
  if(inquark<=6) {
    if(inquark%2==0) ckm = SM().CKM(inquark/2-1,(abs(outquark)-1)/2);
    else             ckm = SM().CKM(abs(outquark)/2-1,(inquark-1)/2);
  }
  // return the answer
  return 0.5*(ME()->contract(_rho)).real()*ckm; 
}



// matrix element for a 1/2 -> 3/2 semi-leptonic decay
double SemiLeptonicBaryonDecayer::halfThreeHalf(const Particle & part,
						const tPDVector & outgoing,
						const vector<Lorentz5Momentum> & momenta,
						MEOption meopt) const {
  // spinors etc for the decaying particle
  if(meopt==Initialize) {
    // spinors and rho
    if(part.id()>0)
      SpinorWaveFunction   ::calculateWaveFunctions(_inHalf,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    else
      SpinorBarWaveFunction::calculateWaveFunctions(_inHalfBar,_rho,
  						    const_ptr_cast<tPPtr>(&part),
  						    incoming);
    // work out the mapping for the lepton vector
    _constants.resize(outgoing.size()+1);
    _ispin.resize(outgoing.size());
    int itemp(1);
    _ibar=0;
    for(int ix=int(outgoing.size()-1);ix>=0;--ix) {
      _ispin[ix]=outgoing[ix]->iSpin();
      if(abs(outgoing[ix]->id())<=16) {
  	itemp*=_ispin[ix];
  	_constants[ix]=itemp;
      }
      else _ibar=ix;
    }
    _constants[outgoing.size()]=1;
    _constants[_ibar]=_constants[_ibar+1];
  }
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1Half,_ispin)));
  // spinors for the decay product
  LorentzPolarizationVector in=UnitRemoval::InvE*part.momentum();
  if(part.id()>0) {
    RSSpinorBarWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalfBar.resize(4);
    _inHalfBar.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalfBar[ihel]=swave.dimensionedWf();
      _inHalfBar[ihel] = _inThreeHalfBar[ihel].dot(in);
    }
  }
  else {
    RSSpinorWaveFunction swave(momenta[0],outgoing[0],Helicity::outgoing);
    _inThreeHalf.resize(4);
    _inHalf.resize(4);
    for(unsigned int ihel=0;ihel<4;++ihel) {
      swave.reset(ihel);
      _inThreeHalf[ihel]=swave.dimensionedWf();
      _inHalf[ihel] = _inThreeHalf[ihel].dot(in);
    }
  }
  // get the information on the form-factor
  int spinin(0),spinout(0),inquark,outquark,spect1,spect2;
  int id0(part.id()),id1(outgoing[0]->id());
  bool cc; int iloc(_form->formFactorNumber(id0,id1,cc));
  _form->formFactorInfo(iloc,spinin,spinout,spect1,spect2,inquark,outquark);
  // work out the value of q and calculate the form factors
  Lorentz5Momentum q(part.momentum()-momenta[0]);
  q.rescaleMass();
  Energy m0(part.mass()),m1(momenta[0].mass());
  Energy2 q2(q.mass2());
  Lorentz5Momentum sum(part.momentum()+momenta[0]);
  // calculate the form factors
  Complex f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a;
  _form->SpinHalfSpinThreeHalfFormFactor(q2,iloc,id0,id1,m0,m1,
  					 f1v,f2v,f3v,f4v,f1a,f2a,f3a,f4a,FlavourInfo());
  LorentzPolarizationVector vtemp;
  complex<InvEnergy2> lS1,lS2,rS1,rS2;
  complex<InvEnergy> lV,rV;
  Complex left,right;
  InvEnergy ms(1./(m0+m1));
  InvEnergy2 ms2(ms*ms);
  if(part.id()>0) {
    left  = f1a-f1v;
    right = f1a+f1v; 
    lS1   = ms2*(f3a-f4a-f3v+f4v);
    rS1   = ms2*(f3a-f4a+f3v-f4v);
    lS2   = ms2*(f4a-f4v);
    rS2   = ms2*(f4a+f4v);
    lV    = ms*(f2a-f2v);
    rV    = ms*(f2a+f2v);
  }
  else {
    left  = conj(f1a+f1v);
    right = conj(f1a-f1v); 
    lS1   = ms2*conj(f3a-f4a+f3v-f4v);
    rS1   = ms2*conj(f3a-f4a-f3v+f4v);
    lS2   = ms2*conj(f4a-f4v);
    rS2   = ms2*conj(f4a+f4v);
    lV    = ms *conj(f2a-f2v);
    rV    = ms *conj(f2a+f2v);
  }
  // calculate the hadronic current for the decay
  LorentzPolarizationVectorE hadron[4][2];
  // construct the vectors for the decay
  Complex scalar1,scalar2;
  complex<Energy> lfact,rfact;
  LorentzPolarizationVectorE tvec;
  LorentzPolarizationVector svec;
  for(unsigned int iya=0;iya<4;++iya) {
    for(unsigned int ixa=0;ixa<2;++ixa) {
      unsigned int ix=iya,iy=ixa;
      if(outgoing[0]->id()<0) swap(ix,iy);
      // scalar like terms
      lfact = _inHalf[iy].leftScalar(_inHalfBar[ix]);
      rfact = _inHalf[iy].rightScalar(_inHalfBar[ix]);
      scalar1 = Complex((lS1*lfact+rS1*rfact)*UnitRemoval::E);
      scalar2 = Complex((lS2*lfact+rS2*rfact)*UnitRemoval::E);
      svec = _inHalf[iy].generalCurrent(_inHalfBar[ix],lV/ms,rV/ms)*ms;
      if(part.id()>0) tvec=_inThreeHalfBar[ix].generalCurrent(_inHalf[iy],left,right);
      else              tvec=_inThreeHalf[iy].generalCurrent(_inHalfBar[ix],left,right);
      hadron[iya][ixa] = tvec+svec*UnitRemoval::E+scalar1*momenta[0]
  	+scalar2*part.momentum();
    }
  }
  // construct the lepton current
  Energy scale;
  int mode((abs(outgoing[1]->id())-11)/12);
  vector<LorentzPolarizationVectorE> 
    lepton(_current->current(tcPDPtr(),FlavourInfo(),mode,-1,scale,tPDVector(outgoing.begin()+1,outgoing.end()),
			     vector<Lorentz5Momentum>(momenta.begin()+1,momenta.end()),meopt));
  vector<unsigned int> ihel(outgoing.size()+1);
  for(unsigned int iya=0;iya<4;++iya) {
    ihel[1]=iya;
    for(unsigned int ixa=0;ixa<2;++ixa) {
      ihel[0]=ixa;
      for(unsigned int lhel=0;lhel<lepton.size();++lhel) {
  	ihel[2] = lhel/2;
  	ihel[3] = lhel%2;
  	(*ME())(ihel) = Complex(lepton[lhel].dot(hadron[iya][ixa])*SM().fermiConstant());
      }
    }  
  }
  // ckm factor
  double ckm(1.);
  if(inquark<=6) {
    if(inquark%2==0){ckm = SM().CKM(inquark/2-1,(abs(outquark)-1)/2);}
    else{ckm = SM().CKM(abs(outquark)/2-1,(inquark-1)/2);}
  }
  // return the answer
  return 0.5*(ME()->contract(_rho)).real()*ckm;
}

// output the setup information for the particle database
void SemiLeptonicBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  DecayIntegrator::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<_maxwgt.size();++ix) {
    output << "insert " << name() << ":MaximumWeight " << ix << " " 
	   << _maxwgt[ix] << " \n";
  }
  _current->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":Current " << _current->name() << " \n";
  _form->dataBaseOutput(output,false,true);
  output << "newdef " << name() << ":FormFactor " << _form->name() << " \n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./StrongHeavyBaryonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StrongHeavyBaryonDecayer class.
//

#include "StrongHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void StrongHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxWeight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) maxWeight_.push_back(mode(ix)->maxWeight());
      else         maxWeight_.push_back(1.);
    }
  }
}

StrongHeavyBaryonDecayer::StrongHeavyBaryonDecayer() {
  // strong couplings of the baryons to pions
  // coupling of the sigma_c to lambda_c pi
  _gsigma_clambda_cpi = 8.88/GeV;
  // coupling of xi*_c to xi_c pi
  _gxistar_cxi_cpi    = 8.34/GeV;
  // strong coupling for lambda_c1 to sigma_c pi
  _flambda_c1sigma_cpi=0.52;
  // strong coupling for xi_c1 to xi_c'
  _fxi_c1xi_cpi = 0.36;
  // strong coupling for lambda_c1star to sigma_c pi
  _flambda_c1starsigma_cpi=21.5/GeV2;
  // strong coupling for xi_ci* to xi_c'
  _fxi_c1starxi_cpi = 20./GeV2;
  // coupling of the sigma_b to lambda_b pi
  _gsigma_blambda_bpi = 8.88/GeV;
  // coupling of xi*_b to xi_b pi
  _gxistar_bxi_bpi    = 8.34/GeV;
  // strong coupling for lambda_b1 to sigma_b pi
  _flambda_b1sigma_bpi=0.52;
  // strong coupling for xi_b1 to xi_b'
  _fxi_b1xi_bpi = 0.36;
  // strong coupling for lambda_b1star to sigma_b pi
  _flambda_b1starsigma_bpi=21.5/GeV2;
  // strong coupling for xi_bi* to xi_b'
  _fxi_b1starxi_bpi = 20./GeV2;
  // intermediates
  generateIntermediates(false);
}

void StrongHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size()||isize!=maxWeight_.size())
    throw InitException() << "Inconsistent parameters in StrongHeavyBaryonDecayer"
			  << "::doinit()" << Exception::abortnow;
  // add the various decay modes
  double or2(1./sqrt(2.)),or3(1./sqrt(3.)),or6(1./sqrt(6.));
  // the decay modes
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first),
    		     getParticleData(outgoing_[ix].second)};
    PhaseSpaceModePtr mode;
    if(in&&out[0]&&out[1]) {
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else {
      mode = PhaseSpaceModePtr();
    }
    addMode(mode);
    if(outgoing_[ix].first==4122&&((incoming_[ix]==4222&&outgoing_[ix].second==211)||
			      (incoming_[ix]==4212&&outgoing_[ix].second==111)||
			      (incoming_[ix]==4112&&outgoing_[ix].second==-211)))
      _prefactor.push_back(-_gsigma_clambda_cpi*GeV*or3);
    else if((incoming_[ix]==4322&&((outgoing_[ix].first==4232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==4312&&((outgoing_[ix].first==4132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   0.5*_gxistar_cxi_cpi*GeV*or3 :
			   _gxistar_cxi_cpi*GeV*or6);
    else if(outgoing_[ix].first==4122&&((incoming_[ix]==4224&&outgoing_[ix].second== 211)||
				   (incoming_[ix]==4214&&outgoing_[ix].second== 111)||
				   (incoming_[ix]==4114&&outgoing_[ix].second==-211)))
      _prefactor.push_back(_gsigma_clambda_cpi*GeV);
    else if((incoming_[ix]==4324&&((outgoing_[ix].first==4232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==4314&&((outgoing_[ix].first==4132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==4232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 0.5*_gxistar_cxi_cpi*GeV :
			   _gxistar_cxi_cpi*or2*GeV);
    else if(incoming_[ix]==101242&&((outgoing_[ix].first==4222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==4212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_c1sigma_cpi);
    else if(incoming_[ix]==101244&&((outgoing_[ix].first==4222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==4212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_c1starsigma_cpi*or3*GeV2);
    else if((incoming_[ix]==102344&&((outgoing_[ix].first==4322&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101344&&((outgoing_[ix].first==4312&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_c1starxi_cpi*0.5*or6*GeV2 :
			   _fxi_c1starxi_cpi*0.5*or3*GeV2);
    else if((incoming_[ix]==102344&&((outgoing_[ix].first==4324&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101344&&((outgoing_[ix].first==4314&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==4324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_c1xi_cpi*0.5*or2 :
			   _fxi_c1xi_cpi*0.5);
    else if((incoming_[ix]==102342&&((outgoing_[ix].first==4322&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101342&&((outgoing_[ix].first==4312&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? _fxi_c1xi_cpi*0.5*or2 :
			   _fxi_c1xi_cpi*0.5);
    else if((incoming_[ix]==102342&&((outgoing_[ix].first==4324&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101342&&((outgoing_[ix].first==4314&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==4324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_c1starxi_cpi*0.5*or6*GeV2 :
			   _fxi_c1starxi_cpi*0.5*or3*GeV2);
    else if(outgoing_[ix].first==5122&&((incoming_[ix]==5222&&outgoing_[ix].second==211)||
				   (incoming_[ix]==5212&&outgoing_[ix].second==111)||
				   (incoming_[ix]==5112&&outgoing_[ix].second==-211)))
      _prefactor.push_back(_gsigma_blambda_bpi*GeV*or3);
    else if((incoming_[ix]==5322&&((outgoing_[ix].first==5232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==5312&&((outgoing_[ix].first==5132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   0.5*_gxistar_bxi_bpi*GeV*or3 : 
			   _gxistar_bxi_bpi*GeV*or6);
    else if(outgoing_[ix].first==5122&&((incoming_[ix]==5224&&outgoing_[ix].second== 211)||
				   (incoming_[ix]==5214&&outgoing_[ix].second== 111)||
				   (incoming_[ix]==5114&&outgoing_[ix].second==-211)))
      _prefactor.push_back(-_gsigma_blambda_bpi*GeV);
    else if((incoming_[ix]==5324&&((outgoing_[ix].first==5232&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5132&&outgoing_[ix].second==211)))||
	    (incoming_[ix]==5314&&((outgoing_[ix].first==5132&&outgoing_[ix].second==111)||
				   (outgoing_[ix].first==5232&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   -0.5*_gxistar_bxi_bpi*GeV : 
			   _gxistar_bxi_bpi*or2*GeV);
    else if(incoming_[ix]==101252&&((outgoing_[ix].first==5222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==5212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_b1sigma_bpi);
    else if(incoming_[ix]==101254&&((outgoing_[ix].first==5222&&outgoing_[ix].second==-211)||
				    (outgoing_[ix].first==5212&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5112&&outgoing_[ix].second== 211)))
      _prefactor.push_back(_flambda_b1starsigma_bpi*or3*GeV*GeV);
    else if((incoming_[ix]==102354&&((outgoing_[ix].first==5322&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101354&&((outgoing_[ix].first==5312&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_b1starxi_bpi*0.5*or6*GeV2 :
			   _fxi_b1starxi_bpi*0.5*or3*GeV2);
    else if((incoming_[ix]==102354&&((outgoing_[ix].first==5324&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101354&&((outgoing_[ix].first==5314&&outgoing_[ix].second== 111)||
				    (outgoing_[ix].first==5324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? _fxi_b1xi_bpi*0.5*or2 :
			   _fxi_b1xi_bpi*0.5);
    else if((incoming_[ix]==102352&&((outgoing_[ix].first==5322&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5312&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101352&&((outgoing_[ix].first==5312&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5322&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_b1xi_bpi*0.5*or2 : _fxi_b1xi_bpi*0.5);
    else if((incoming_[ix]==102352&&((outgoing_[ix].first==5324&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5314&&outgoing_[ix].second== 211)))||
	    (incoming_[ix]==101352&&((outgoing_[ix].first==5314&&outgoing_[ix].second== 111)||
				     (outgoing_[ix].first==5324&&outgoing_[ix].second==-211))))
      _prefactor.push_back(outgoing_[ix].second==111 ? 
			   _fxi_b1starxi_bpi*0.5*or6*GeV2 :
			   _fxi_b1starxi_bpi*0.5*or3*GeV2);
    else
      throw InitException() << "Unknown mode in StrongHeavyBaryonDecayer::doinit()"
			    << Exception::abortnow;
  }
}

void StrongHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_gsigma_clambda_cpi,1./GeV) << ounit(_gxistar_cxi_cpi,1./GeV) 
     << _flambda_c1sigma_cpi << _fxi_c1xi_cpi << ounit(_flambda_c1starsigma_cpi,1./GeV2) 
     << ounit(_fxi_c1starxi_cpi,1./GeV2) << ounit(_gsigma_blambda_bpi,1./GeV) 
     << ounit(_gxistar_bxi_bpi,1./GeV) << _flambda_b1sigma_bpi 
     << _fxi_b1xi_bpi << ounit(_flambda_b1starsigma_bpi,1./GeV2) 
     << ounit(_fxi_b1starxi_bpi,1./GeV2) 
     << incoming_ << outgoing_ << maxWeight_ << _prefactor << modeType_;
}

void StrongHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_gsigma_clambda_cpi,1./GeV) >> iunit(_gxistar_cxi_cpi,1./GeV) 
     >> _flambda_c1sigma_cpi >> _fxi_c1xi_cpi >> iunit(_flambda_c1starsigma_cpi,1./GeV2) 
     >> iunit(_fxi_c1starxi_cpi,1./GeV2) >> iunit(_gsigma_blambda_bpi,1./GeV) 
     >> iunit(_gxistar_bxi_bpi,1./GeV) >> _flambda_b1sigma_bpi 
     >> _fxi_b1xi_bpi >> iunit(_flambda_b1starsigma_bpi,1./GeV2) 
     >> iunit(_fxi_b1starxi_bpi,1./GeV2) 
     >> incoming_ >> outgoing_ >> maxWeight_ >> _prefactor >> modeType_;
}

int StrongHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==incoming_[ix]) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-incoming_[ix]) {
      if((id1==-outgoing_[ix].first&&id2==-outgoing_[ix].second)||
	 (id2==-outgoing_[ix].first&&id1==-outgoing_[ix].second)) {
	imode=ix;
	cc=true;
      }
      if(((id1==-outgoing_[ix].first&&id2==outgoing_[ix].second)||
	  (id2==-outgoing_[ix].first&&id1==outgoing_[ix].second))&&
	 (outgoing_[ix].second==111||outgoing_[ix].second==221||outgoing_[ix].second==331||
	  outgoing_[ix].second==223||outgoing_[ix].second==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<StrongHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigStrongHeavyBaryonDecayer("Herwig::StrongHeavyBaryonDecayer", "HwBaryonDecay.so");

void StrongHeavyBaryonDecayer::Init() {

  static ClassDocumentation<StrongHeavyBaryonDecayer> documentation
    ("The StrongHeavyBaryonDecayer class performs the strong decays of"
     " baryons containing a heavy quark using the results of hep-ph/9904421.",
     "The strong decays of the heavy baryons were simulated using the results of"
     "\\cite{Ivanov:1999bk}.",
     "\\bibitem{Ivanov:1999bk}\n"
     "M.~A.~Ivanov, J.~G.~Korner, V.~E.~Lyubovitskij and A.~G.~Rusetsky,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 094002\n"
     "[arXiv:hep-ph/9904421].\n"
     "%%CITATION = PHRVA,D60,094002;%%\n");

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegSigma_cLambda_cPi
    ("gSigma_cLambda_cPi",
     "The coupling of the Sigma_c to Lambda_c pi",
     &StrongHeavyBaryonDecayer::_gsigma_clambda_cpi, 1./GeV, 8.8/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegXiStar_cXi_cPi
    ("gXiStar_cXi_cPi",
     "The coupling of the Xi*_c to Xi_c pi",
     &StrongHeavyBaryonDecayer::_gxistar_cxi_cpi, 1./GeV, 8.34/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefLambda_c1Sigma_cPi
    ("fLambda_c1Sigma_cPi",
     "The coupling of the Lambda_c1 to Sigma_c pi",
     &StrongHeavyBaryonDecayer::_flambda_c1sigma_cpi, 0.52, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefXi_c1Xi_cPi
    ("fXi_c1Xi_cPi",
     "The coupling of the Xi_c1 to Xi_c pi",
     &StrongHeavyBaryonDecayer::_fxi_c1xi_cpi, 0.36, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefLambda_c1starSigma_cPi
    ("fLambda_c1*Sigma_cPi",
     "The coupling of Lambda_c1* to Sigma_c and pi",
     &StrongHeavyBaryonDecayer::_flambda_c1starsigma_cpi, 1./GeV2, 21.5/GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegSigma_bLambda_bPi
    ("gSigma_bLambda_bPi",
     "The coupling of the Sigma_b to Lambda_b pi",
     &StrongHeavyBaryonDecayer::_gsigma_blambda_bpi, 1./GeV, 8.8/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefXi_c1starXi_cPi
    ("fXi_c1*Xi_cPi",
     "The coupling of Xi_c1* to Xi_c and pi",
     &StrongHeavyBaryonDecayer::_fxi_c1starxi_cpi, 1./GeV2, 20./GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy> interfacegXiStar_bXi_bPi
    ("gXiStar_bXi_bPi",
     "The coupling of the Xi*_b to Xi_b pi",
     &StrongHeavyBaryonDecayer::_gxistar_bxi_bpi, 1./GeV, 8.34/GeV, ZERO, 20.0/GeV,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefLambda_b1Sigma_bPi
    ("fLambda_b1Sigma_bPi",
     "The coupling of the Lambda_b1 to Sigma_b pi",
     &StrongHeavyBaryonDecayer::_flambda_b1sigma_bpi, 0.52, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,double> interfacefXi_b1Xi_bPi
    ("fXi_b1Xi_bPi",
     "The coupling of the Xi_b1 to Xi_b pi",
     &StrongHeavyBaryonDecayer::_fxi_b1xi_bpi, 0.36, 0, 10,
     false, false, false);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefLambda_b1starSigma_bPi
    ("fLambda_b1*Sigma_bPi",
     "The coupling of Lambda_b1* to Sigma_b and pi",
     &StrongHeavyBaryonDecayer::_flambda_b1starsigma_bpi, 1./GeV2, 21.5/GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Parameter<StrongHeavyBaryonDecayer,InvEnergy2> interfacefXi_b1starXi_bPi
    ("fXi_b1*Xi_bPi",
     "The coupling of Xi_b1* to Xi_b and pi",
     &StrongHeavyBaryonDecayer::_fxi_b1starxi_bpi, 1./GeV2, 20./GeV2,
     ZERO, 100.0/GeV2,
     false, false, true);

  static Command<StrongHeavyBaryonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing scalars, coupling(MeV) and max weight for a decay",
     &StrongHeavyBaryonDecayer::setUpDecayMode, false);
  
  static Deleted<StrongHeavyBaryonDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceOutgoingB
    ("OutgoingB","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceOutgoingM
    ("OutgoingM","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceModeType
    ("ModeType","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<StrongHeavyBaryonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in StrongHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

}

// couplings for spin-1/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
		       Complex & A,Complex & B) const {
  useMe();
  if(modeType_[imode]==0) {
    A = _prefactor[imode];
    B = 0.;
  }
  else if(modeType_[imode]==1) {
    A = 0.;
    B = 0.5*_prefactor[imode]*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV;
  }
  else
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "halfHalfScalarCoupling() " 
				 << Exception::abortnow;
}

// couplings for spin-1/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
halfThreeHalfScalarCoupling(int imode, Energy m0, Energy m1, Energy m2,
			    Complex& A,Complex& B) const {
  useMe();
  if(modeType_[imode]==1) {
    A = _prefactor[imode]*(m0+m1)/GeV;
    B = 0.;
  }
  else if(modeType_[imode]==2) {
    A = 0.;
    B = 0.5*_prefactor[imode]*(m0+m1)*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV2;
  }
  else {
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "halfThreeHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

// couplings for spin-3/2 to spin-3/2 spin-0
void StrongHeavyBaryonDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Energy,Energy,Energy,
				 Complex& A1,Complex& A2,Complex& B1,Complex& B2) const {
  useMe();
  if(modeType_[imode]==0) {
    A1 = _prefactor[imode];
    B1 = 0.;
    A2=0.;
    B2=0.;
  }
  else {
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "threeHalfThreeHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

// couplings for spin-3/2 to spin-1/2 spin-0
void StrongHeavyBaryonDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy m2,
			    Complex& A,Complex& B) const {
  useMe();
  if(modeType_[imode]==1) {
    A = _prefactor[imode]*(m0+m1)/GeV;
    B = 0.;
  }
  else if(modeType_[imode]==2) {
    A = 0.;
    B = 0.5*_prefactor[imode]*(m0+m1)*((m0+m1)*(m0+m1)-m2*m2)/m0/GeV2;
  }
  else {
    throw Exception() << "Unknown mode in  StrongHeavyBaryonDecayer::"
				 << "threeHalfHalfScalarCoupling() " 
				 << Exception::abortnow;
  }
}

void StrongHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":gSigma_cLambda_cPi " 
	 << _gsigma_clambda_cpi*GeV << "\n";
  output << "newdef " << name() << ":gXiStar_cXi_cPi " 
	 << _gxistar_cxi_cpi*GeV << "\n";
  output << "newdef " << name() << ":fLambda_c1Sigma_cPi " 
	 << _flambda_c1sigma_cpi << "\n";
  output << "newdef " << name() << ":fXi_c1Xi_cPi " 
	 << _fxi_c1xi_cpi << "\n";
  output << "newdef " << name() << ":fLambda_c1*Sigma_cPi " 
	 << _flambda_c1starsigma_cpi*GeV2 << "\n";
  output << "newdef " << name() << ":fXi_c1*Xi_cPi " 
	 << _fxi_c1starxi_cpi*GeV2 << "\n";
  output << "newdef " << name() << ":gSigma_bLambda_bPi " 
	 << _gsigma_blambda_bpi*GeV << "\n";
  output << "newdef " << name() << ":gXiStar_bXi_bPi " 
	 << _gxistar_bxi_bpi*GeV << "\n";
  output << "newdef " << name() << ":fLambda_b1Sigma_bPi " 
	 << _flambda_b1sigma_bpi << "\n";
  output << "newdef " << name() << ":fXi_b1Xi_bPi " 
	 << _fxi_b1xi_bpi << "\n";
  output << "newdef " << name() << ":fLambda_b1*Sigma_bPi " 
	 << _flambda_b1starsigma_bpi*GeV2 << "\n";
  output << "newdef " << name() << ":fXi_b1*Xi_bPi " 
	 << _fxi_b1starxi_bpi*GeV2 << "\n";
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode " << incoming_[ix] << " "
	   << outgoing_[ix].first << " " << outgoing_[ix].second << " " << modeType_[ix]
	   << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string StrongHeavyBaryonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2 or 3/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2 or 3/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the type
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int itype = stoi(stype);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  modeType_.push_back(itype);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./NonLeptonicHyperonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonLeptonicHyperonDecayer class.
//

#include "NonLeptonicHyperonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void NonLeptonicHyperonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

void NonLeptonicHyperonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size() || isize!=a_.size()||
     isize!=b_.size()        || isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "NonLeptonicHyperonDecayer::doinit()" 
			  << Exception::runerror;
  // set up the decay modes
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix].first ),
    		     getParticleData(outgoing_[ix].second)};
    double wgtmax = maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    PhaseSpaceModePtr mode(new_ptr(PhaseSpaceMode(in,out,wgtmax)));
    addMode(mode);
    // test of the asummetries
    // Energy e1 = (sqr(in->mass())+sqr(out[0]->mass())-
    // 		 sqr(out[1]->mass()))/2./in->mass();
    // double btemp = b_[ix]*sqrt((e1-out[0]->mass())/(e1+out[0]->mass()));
    // double alpha = -2.*(a_[ix]*btemp)/(sqr(a_[ix])+sqr(btemp));
    // generator()->log() << "Asymmetry parameter for " << in->PDGName() << "->"
    // 		       << out[0]->PDGName() << "," << out[1]->PDGName()
    // 		       << " = " << alpha << "\n";
  }
}

NonLeptonicHyperonDecayer::NonLeptonicHyperonDecayer() {
  // intermediates
  generateIntermediates(false);
}

int NonLeptonicHyperonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing pa4rticles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0=parent->id();
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do {
    if(id0==incoming_[ix]) {
      if((id1==outgoing_[ix].first&&id2==outgoing_[ix].second)||
	 (id2==outgoing_[ix].first&&id1==outgoing_[ix].second)) imode=ix;
    }
    else if(id0==-incoming_[ix]) {
      if((id1==-outgoing_[ix].first&&id2==-outgoing_[ix].second)||
	 (id2==-outgoing_[ix].first&&id1==-outgoing_[ix].second)) imode=ix;
      if(((id1==-outgoing_[ix].first&&id2==outgoing_[ix].second)||
	  (id2==-outgoing_[ix].first&&id1==outgoing_[ix].second))&&
	 (outgoing_[ix].second==111||outgoing_[ix].second==221||outgoing_[ix].second==331||
	  outgoing_[ix].second==223||outgoing_[ix].second==333)) imode=ix;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  // charge conjugate
  cc=id0<0;
  // return the answer
  return imode;
}

void NonLeptonicHyperonDecayer::persistentOutput(PersistentOStream & os) const {
  os << incoming_ << outgoing_ << a_ << b_ << maxweight_;
}

void NonLeptonicHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incoming_ >> outgoing_ >> a_ >> b_ >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NonLeptonicHyperonDecayer,Baryon1MesonDecayerBase>
describeHerwigNonLeptonicHyperonDecayer("Herwig::NonLeptonicHyperonDecayer", "HwBaryonDecay.so");

void NonLeptonicHyperonDecayer::Init() {

  static ClassDocumentation<NonLeptonicHyperonDecayer> documentation
    ("The NonLeptonicHyperonDecayer class performs the non-leptonic"
     " weak decay of the hyperons.",
     "The non-leptonic hyperon decays were simulated using the "
     "NonLeptonicHyperonDecayer class which implements the model of"
     "\\cite{Borasoy:1999md}",
     "\\bibitem{Borasoy:1999md}\n"
     "B.~Borasoy and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 59} (1999) 094025 [arXiv:hep-ph/9902351].\n"
     "%%CITATION = PHRVA,D59,094025;%%\n");

  static Command<NonLeptonicHyperonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryon & meson, A and B couplings and max weight for a decay",
     &NonLeptonicHyperonDecayer::setUpDecayMode, false);

  static Deleted<NonLeptonicHyperonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceIncomingBaryon
    ("IncomingBaryon","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceOutgoingBaryon
    ("OutgoingBaryon","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceOutgoingMeson
    ("OutgoingMeson","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceA
    ("CouplingA","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<NonLeptonicHyperonDecayer> interfaceB
    ("CouplingB","The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");
}

// couplings for spin-1/2 to spin-1/2 spin-0
void NonLeptonicHyperonDecayer::halfHalfScalarCoupling(int imode,
						       Energy,Energy,Energy,
						       Complex& A,Complex& B) const {
  useMe();
  A=a_[imode];
  B=b_[imode];
}

void NonLeptonicHyperonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix].first << " " << outgoing_[ix].second
	   << " " << a_[ix] << " " << b_[ix] << " " << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

string NonLeptonicHyperonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  pair<long,long> out;
  out.first = stoi(stype);
  pData = getParticleData(out.first);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out.first) + "does not have spin 1/2";
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  out.second = stoi(stype);
  pData = getParticleData(out.second);
  if(!pData)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not exist";
  if(pData->iSpin()!=PDT::Spin0)
    return "Second outgoing particle with id " + std::to_string(out.second) + "does not have spin 0";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  a_.push_back(g);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  g = stof(stype);
  b_.push_back(g);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./NonLeptonicOmegaDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NonLeptonicOmegaDecayer class.
//

#include "NonLeptonicOmegaDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

NonLeptonicOmegaDecayer::NonLeptonicOmegaDecayer() {
  // couplings for the decays
  _fstar  =  0.5e-7;
  _dstar  = -3.98e-7;
  _omegad = -1.16e-8/1.8;
  _omegaf =  1.53e-8/1.8;
  _cbstar =  1.35;
  _sc     = -0.85;
  _c      =  1.50;
  _fpi    =  92.4*MeV;
  _hc     =  0.39e-7*GeV;
  _hpi    =  3.2e-7;
  _d      =  0.44e-7*GeV;
  _f      = -0.50e-7*GeV;
  // massses of the particles
  _mlambda = 1115.683*MeV;
  _mxi     = 1314.830*MeV;
  _momega  = 1672.450*MeV;
  _mxistar = 1531.800*MeV;
  _mpip    =  139.570*MeV;
  _mpi0    =  134.977*MeV;
  _mkp     =  493.667*MeV;
  _mk0     =  497.648*MeV;
  _mbstar  =  1620   *MeV;
  _mr      =  1500   *MeV;
  // use local values for the masses for the couplings
  _localmasses=true;
  // the PDG codes for the modes
  _incomingB = 3334;
  _outgoingB.resize(3);_outgoingM.resize(3);_maxweight.resize(3);
  _outgoingB[0] = 3122;_outgoingM[0] =-321;_maxweight[0]=1.5;
  _outgoingB[1] = 3322;_outgoingM[1] =-211;_maxweight[1]=0.4;
  _outgoingB[2] = 3312;_outgoingM[2] = 111;_maxweight[2]=0.2;
  // intermediates
  generateIntermediates(false);
}

void NonLeptonicOmegaDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

void NonLeptonicOmegaDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // reset the masses if needed
  if(!_localmasses) {
    _mlambda = getParticleData(3122)->mass();
    _mxi     = getParticleData(3322)->mass();
    _momega  = getParticleData(3334)->mass();
    _mxistar = getParticleData(3324)->mass();
    _mpip    = getParticleData(211)->mass();
    _mpi0    = getParticleData(111)->mass();
    _mkp     = getParticleData(321)->mass();
    _mk0     = getParticleData(311)->mass();
  }
  // calculate the couplings
  _a.resize(3);_b.resize(3);
  // couplings for lambda K (N.B. sign of B due to gamma_5 defn)
  _a[0] =  0.5*_c/sqrt(3.)/_fpi*((_d-3.*_f)/(_mlambda-_mxi)
			       +_hc/(_momega-_mxistar))
    +0.5*_cbstar/sqrt(3.)/_fpi*(_dstar-3.*_fstar)/(_mlambda/_mbstar-1.);
  _b[0] =  0.5*_sc/sqrt(3.)/_fpi*(_omegad-3.*_omegaf)/(_mlambda/_mr-1.);
  // couplings for xi0 pi-
  _a[1] = _c/   sqrt(2.)/_fpi*(_hc/3./(_momega-_mxistar)
			       +_hpi*_mpip*_mpip/2./(_mkp*_mkp-_mpip*_mpip));
  _b[1] = ZERO;
  // couplings for xi- pi0
  _a[2] = _c/2./         _fpi*(_hc/3./(_momega-_mxistar)
			       +_hpi*_mpi0*_mpi0/2./(_mk0*_mk0-_mpi0*_mpi0));
  _b[2] = ZERO;
  // set up the decay modes
  for(unsigned int ix=0;ix<_outgoingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_incomingB);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(_outgoingM[ix])};
    double wgtmax = _maxweight.size()>numberModes() ?
      _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

int NonLeptonicOmegaDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  do {
    if(id0==_incomingB) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) imode=ix;
    }
    else if(id0==-_incomingB) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) imode=ix;
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) imode=ix;
    }
    ++ix;
  }
  while(ix<_outgoingB.size()&&imode<0);
  // charge conjugation
  cc=id0<0;
  // return the answer
  return imode;
}

void NonLeptonicOmegaDecayer::persistentOutput(PersistentOStream & os) const {
  os << _dstar << _fstar << _omegad << _omegaf << _cbstar << _sc << _c 
     << ounit(_fpi,GeV) 
     << ounit(_hc,GeV) << _hpi << ounit(_d,GeV) << ounit(_f,GeV) << ounit(_mlambda,GeV) 
     << ounit(_mxi,GeV) << ounit(_momega,GeV) << ounit(_mxistar,GeV) << ounit(_mpip,GeV) 
     << ounit(_mpi0,GeV) << ounit(_mkp,GeV) << ounit(_mk0,GeV) << ounit(_mbstar,GeV) 
     << ounit(_mr,GeV) << _localmasses << _incomingB << _outgoingB << _outgoingM 
     << ounit(_a,1./GeV) << ounit(_b,1./GeV) << _maxweight;
}

void NonLeptonicOmegaDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _dstar >> _fstar >> _omegad >> _omegaf >> _cbstar >> _sc >> _c 
     >> iunit(_fpi,GeV) 
     >> iunit(_hc,GeV) >> _hpi >> iunit(_d,GeV) >> iunit(_f,GeV) >> iunit(_mlambda,GeV) 
     >> iunit(_mxi,GeV) >> iunit(_momega,GeV) >> iunit(_mxistar,GeV) >> iunit(_mpip,GeV) 
     >> iunit(_mpi0,GeV) >> iunit(_mkp,GeV) >> iunit(_mk0,GeV) >> iunit(_mbstar,GeV) 
     >> iunit(_mr,GeV) >> _localmasses >> _incomingB >> _outgoingB >> _outgoingM 
     >> iunit(_a,1./GeV) >> iunit(_b,1./GeV) >> _maxweight;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<NonLeptonicOmegaDecayer,Baryon1MesonDecayerBase>
describeHerwigNonLeptonicOmegaDecayer("Herwig::NonLeptonicOmegaDecayer", "HwBaryonDecay.so");

void NonLeptonicOmegaDecayer::Init() {

  static ClassDocumentation<NonLeptonicOmegaDecayer> documentation
    ("The NonLeptonicOmegaDecayer class performs the non-leptonic decays"
     " of the omega based on the results of hep-ph/9905398.",
     "The  non-leptonic decays of the Omega baryon were simulated "
     "using the NonLeptonicOmegaDecayer class based on the results of"
     "\\cite{Borasoy:1999ip}.",
     "\\bibitem{Borasoy:1999ip}\n"
     "B.~Borasoy and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 054021 [arXiv:hep-ph/9905398].\n"
     "%%CITATION = PHRVA,D60,054021;%%\n");

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceDStar
    ("DStar",
     "The d* coupling from hep-ph/9905398 multiplied by MB*.",
     &NonLeptonicOmegaDecayer::_dstar, -3.98e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceFStar
    ("FStar",
     "The f* coupling from hep-ph/9905398 multiplied by MB*.",
     &NonLeptonicOmegaDecayer::_fstar, 0.5e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceomegad
    ("omegad",
     "The omega_d coupling from hep-ph/9905398 multiplied by MR.",
     &NonLeptonicOmegaDecayer::_omegad, -1.16e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceomegaf
    ("omegaf",
     "The omega_f coupling from hep-ph/9905398 multiplied by MR.",
     &NonLeptonicOmegaDecayer::_omegaf, 1.53e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceCBstar
    ("CBstar",
     "The C_b* coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_cbstar, 1.35, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacesc
    ("sc",
     "The sc coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_sc,-0.85, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfaceC
    ("C",
     "The C coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_c, 1.5, -10.0, 10.0,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &NonLeptonicOmegaDecayer::_fpi, MeV, 92.4*MeV, ZERO, 200.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfacehc
    ("hc",
     "The h_c coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_hc, GeV, 0.39e-7*GeV, -10.0e-7*GeV, 10.0e-7*GeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,double> interfacehpi
    ("hpi",
     "The h_pi coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_hpi, 3.2e-7, -10.0e-7, 10.0e-7,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaced
    ("d",
     "The d coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_d, GeV, 0.44e-7*GeV, -1e-6*GeV, 1e-6*GeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfacef
    ("f",
     "The f coupling from hep-ph/9905398",
     &NonLeptonicOmegaDecayer::_f, GeV, -0.50e-7*GeV, -1e-6*GeV, 1e-6*GeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMLambda
    ("MLambda",
     "The mass of the Lambda baryon",
     &NonLeptonicOmegaDecayer::_mlambda, MeV, 1115.683*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMXi
    ("MXi",
     "The mass of the Xi baryon",
     &NonLeptonicOmegaDecayer::_mxi, MeV, 1314.830*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMOmega
    ("MOmega",
     "The mass of the Omega baryon",
     &NonLeptonicOmegaDecayer::_momega, MeV, 1672.450*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMXiStar
    ("MXiStar",
     "The mass of the XiStar baryon",
     &NonLeptonicOmegaDecayer::_mxistar, MeV, 1531.800*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMpiplus
    ("Mpiplus",
     "The mass of the charged pion",
     &NonLeptonicOmegaDecayer::_mpip, MeV, 139.57*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMKplus
    ("MKplus",
     "The mass of the charged kaon",
     &NonLeptonicOmegaDecayer::_mkp, MeV, 493.667*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMpi0
    ("Mpi0",
     "The mass of the neutral pion",
     &NonLeptonicOmegaDecayer::_mpi0, MeV, 134.977*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMK0
    ("MK0",
     "The mass of the neutral kaon",
     &NonLeptonicOmegaDecayer::_mk0, MeV, 497.648*MeV, ZERO, 1000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMBstar
    ("MBstar",
     "The mass of the excited B* resonnaces",
     &NonLeptonicOmegaDecayer::_mbstar, MeV, 1620.0*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Parameter<NonLeptonicOmegaDecayer,Energy> interfaceMR
    ("MR",
     "The mass of the excited R resonnaces",
     &NonLeptonicOmegaDecayer::_mr, MeV, 1620.0*MeV, ZERO, 10000.0*MeV,
     false, false, true);

  static Switch<NonLeptonicOmegaDecayer,bool> interfaceLocalMasses
    ("LocalMasses",
     "Use local values of all the masses for the couplings.",
     &NonLeptonicOmegaDecayer::_localmasses, true, false, false);
  static SwitchOption interfaceLocalMassesLocal
    (interfaceLocalMasses,
     "Local",
     "Use local values",
     true);
  static SwitchOption interfaceLocalMassesNonLocal
    (interfaceLocalMasses,
     "ParticleData",
     "Use values from the particle data objects",
     false);

  static ParVector<NonLeptonicOmegaDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &NonLeptonicOmegaDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void NonLeptonicOmegaDecayer::threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,
							  Energy,Complex& A,
							  Complex& B) const {
  useMe();
  A = _a[imode]*(m0+m1);
  B = _b[imode]*(m0+m1);
}

void NonLeptonicOmegaDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  output << "newdef " << name() << ":DStar " << _dstar<< "\n";
  output << "newdef " << name() << ":FStar " << _fstar << "\n";
  output << "newdef " << name() << ":omegad " << _omegad<< "\n";
  output << "newdef " << name() << ":omegaf " << _omegaf<< "\n";
  output << "newdef " << name() << ":CBstar " << _cbstar<< "\n";
  output << "newdef " << name() << ":sc " << _sc << "\n";
  output << "newdef " << name() << ":C " << _c << "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":hc " << _hc/GeV << "\n";
  output << "newdef " << name() << ":hpi " << _hpi<< "\n";
  output << "newdef " << name() << ":d " << _d/GeV << "\n";
  output << "newdef " << name() << ":f " << _f/GeV << "\n";
  output << "newdef " << name() << ":MLambda " << _mlambda/MeV << "\n";
  output << "newdef " << name() << ":MXi " << _mxi/MeV << "\n";
  output << "newdef " << name() << ":MOmega " << _momega/MeV << "\n";
  output << "newdef " << name() << ":MXiStar " << _mxistar/MeV << "\n";
  output << "newdef " << name() << ":Mpiplus " << _mpip/MeV << "\n";
  output << "newdef " << name() << ":MKplus " << _mkp/MeV << "\n";
  output << "newdef " << name() << ":Mpi0 " << _mpi0/MeV << "\n";
  output << "newdef " << name() << ":MK0 " << _mk0/MeV << "\n";
  output << "newdef " << name() << ":MBstar " << _mbstar/MeV << "\n";
  output << "newdef " << name() << ":MR " << _mr/MeV << "\n";
  output << "newdef " << name() << ":LocalMasses " << _localmasses << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\""
		    << fullName() << "\";" << endl;
}
#line 1 "./RadiativeHyperonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeHyperonDecayer class.
//

#include "RadiativeHyperonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeHyperonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

void RadiativeHyperonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  unsigned int isize(incomingB_.size());
  if(isize!=outgoingB_.size()||isize!=A_.size()||
     isize!=B_.size()        ||isize!=maxweight_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeHyperonDecayer::doinit()" 
			  << Exception::runerror;
  // set up the decay modes
  tPDPtr photon = getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    tPDPtr in = getParticleData(incomingB_[ix]);
    tPDVector out={getParticleData(outgoingB_[ix]),photon};
    double wgtmax = maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

int RadiativeHyperonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing pa4rticles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0=parent->id();
  int id1(children[0]->id());
  int id2(children[1]->id());
  if(id1==ParticleID::gamma) swap(id1,id2);
  if(id2!=ParticleID::gamma) return imode;
  unsigned int ix(0);
  do {
    if(id0==incomingB_[ix] && id1==outgoingB_[ix]) imode=ix;
    else if(id0==-incomingB_[ix] && id1==-outgoingB_[ix]) imode=ix;
    ++ix;
  }
  while(ix<incomingB_.size()&&imode<0);
  // charge conjugate
  cc=id0<0;
  // return the answer
  return imode;
}

IBPtr RadiativeHyperonDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr RadiativeHyperonDecayer::fullclone() const {
  return new_ptr(*this);
}

void RadiativeHyperonDecayer::persistentOutput(PersistentOStream & os) const {
  os << incomingB_ << outgoingB_ << ounit(A_,1./GeV) << ounit(B_,1./GeV) 
     << maxweight_;
}

void RadiativeHyperonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> incomingB_ >> outgoingB_ >> iunit(A_,1./GeV) >> iunit(B_,1./GeV) 
     >> maxweight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeHyperonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeHyperonDecayer("Herwig::RadiativeHyperonDecayer", "HwBaryonDecay.so");

void RadiativeHyperonDecayer::Init() {

  static ClassDocumentation<RadiativeHyperonDecayer> documentation
    ("The RadiativeHyperonDecayer class performs the radiative decays of hyperons.",
     "The radiative hyperons decays were simulated using the RadiativeHyperonDecayer"
     " class which implements the results of \\cite{Borasoy:1999nt}.",
     "\\bibitem{Borasoy:1999nt}\n"
     "B.~Borasoy and B.~R.~Holstein,\n"
     "Phys.\\ Rev.\\  D {\\bf 59} (1999) 054019 [arXiv:hep-ph/9902431].\n"
     "%%CITATION = PHRVA,D59,054019;%%\n");

  static Command<RadiativeHyperonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryon , A and B couplings (1/GeV) and max weight for a decay",
     &RadiativeHyperonDecayer::setUpDecayMode, false);
  
  static Deleted<RadiativeHyperonDecayer> interfaceMaxWeight
    ("MaxWeight",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceIncomingBaryon
    ("IncomingBaryon",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceOutgoingBaryon
    ("OutgoingBaryon",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceCouplingA
    ("CouplingA",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHyperonDecayer> interfaceCouplingB
    ("CouplingB",
     "The old methods of setting up a decay in NonLeptonicHyperonDecayer have been deleted, please use SetUpDecayMode");

}

void RadiativeHyperonDecayer::halfHalfVectorCoupling(int imode, Energy m0, Energy m1,
						     Energy,Complex& A1,Complex& A2,
						     Complex& B1,Complex& B2) const {
  useMe();
  A1 = A_[imode]*(m0+m1);
  B1 = B_[imode]*(m0-m1);
  A2 = 2.*A_[imode]*(m0+m1);
  B2 = 2.*B_[imode]*(m0+m1);
}

void RadiativeHyperonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incomingB_[ix] << " " << outgoingB_[ix] << " "
	   << A_[ix]*GeV << " " << B_[ix]*GeV << " "
	   << maxweight_[ix] << "\n";
  }
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}


string RadiativeHyperonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out;
  out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "First outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "First outgoing particle with id " + std::to_string(out) + "does not have spin 1/2";
  // get the coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double g = stof(stype);
  A_.push_back(g/GeV);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  g = stof(stype);
  B_.push_back(g/GeV);
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incomingB_.push_back(in);
  outgoingB_.push_back(out);
  maxweight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./RadiativeHeavyBaryonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeHeavyBaryonDecayer class.
//

#include "RadiativeHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxWeight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) maxWeight_.push_back(mode(ix)->maxWeight());
      else         maxWeight_.push_back(1.);
    }
  }
}

void RadiativeHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the parameters are consistent
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size() ||isize!=maxWeight_.size()||isize!=M1Coupling_.size()||
     isize!=E1Coupling_.size()||isize!=modeType_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeHeavyBaryonDecayer::doinit()" 
			  << Exception::abortnow;
  // the decay modes
  tPDPtr photon = getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix]),photon};
    PhaseSpaceModePtr mode;
    if(in&&out[0]) {
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else {
      mode = PhaseSpaceModePtr();
    }
    addMode(mode);
  }
}

int RadiativeHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					    const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id()),ibaryon;
  if(id1==ParticleID::gamma){ibaryon=id2;}
  else if(id2==ParticleID::gamma){ibaryon=id1;}
  else {return imode;}
  unsigned int ix(0);
  do {
    if(     id0== incoming_[ix]&&ibaryon== outgoing_[ix]) {
      imode=ix;
      cc=false;
    }
    else if(id0==-incoming_[ix]&&ibaryon==-outgoing_[ix]) {
      imode=ix;
      cc=true;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void RadiativeHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(M1Coupling_,1./GeV) << ounit(E1Coupling_,1./GeV2) 
     << incoming_ << outgoing_ << modeType_ << maxWeight_;
}

void RadiativeHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(M1Coupling_,1./GeV) >> iunit(E1Coupling_,1./GeV2) 
     >> incoming_ >> outgoing_ >> modeType_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeHeavyBaryonDecayer("Herwig::RadiativeHeavyBaryonDecayer", "HwBaryonDecay.so");

void RadiativeHeavyBaryonDecayer::Init() {

  static ClassDocumentation<RadiativeHeavyBaryonDecayer> documentation
    ("The RadiativeHeavyBaryonDecayer class is designed for the radiative decays of"
     " heavy baryons.",
     "The radiative decays of the heavy baryons were simulated using the results of"
     "\\cite{Ivanov:1999bk,Ivanov:1998wj}.",
     "\\bibitem{Ivanov:1999bk}\n"
     "M.~A.~Ivanov, J.~G.~Korner, V.~E.~Lyubovitskij and A.~G.~Rusetsky,\n"
     "Phys.\\ Rev.\\  D {\\bf 60} (1999) 094002\n"
     "[arXiv:hep-ph/9904421].\n"
     "%%CITATION = PHRVA,D60,094002;%%\n"
     "\\bibitem{Ivanov:1998wj}\n"
     "M.~A.~Ivanov, J.~G.~Korner and V.~E.~Lyubovitskij,\n"
     "Phys.\\ Lett.\\  B {\\bf 448} (1999) 143 [arXiv:hep-ph/9811370].\n"
     "%%CITATION = PHLTA,B448,143;%%\n");

  static Command<RadiativeHeavyBaryonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryons, type of decay, M1(1/MeV) and E1(1/MeV2) couplings and max weight for a decay",
     &RadiativeHeavyBaryonDecayer::setUpDecayMode, false);

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceIncoming
    ("Incoming","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceOutgoingB
    ("OutgoingB","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceModeType
    ("ModeType","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");
  
  static Deleted<RadiativeHeavyBaryonDecayer> interfaceM1Coupling
    ("M1Coupling","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceE1Coupling
    ("E1Coupling","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");

  static Deleted<RadiativeHeavyBaryonDecayer> interfaceMaxWeight
    ("MaxWeight","The old methods of setting up a decay in RadiativeHeavyBaryonDecayer have been deleted, please use SetUpDecayMode");
}

void RadiativeHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix] << " " << " " << modeType_[ix]
	   << M1Coupling_[ix]*MeV << " " << E1Coupling_[ix]*MeV2 << " "
	   << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void RadiativeHeavyBaryonDecayer::halfHalfVectorCoupling(int imode,Energy m0,Energy m1,
							 Energy,Complex& A1,
							 Complex& A2,Complex& B1,
							 Complex& B2) const {
  useMe();
  if(modeType_[imode]==0) {
    B1=-0.5*(sqr(m0) - sqr(m1))*E1Coupling_[imode];
    B2=-    (m0+m1)*(m0+m1)*E1Coupling_[imode];
    A1=0.;
    A2=0.;
  }
  else if(modeType_[imode]==1) {
    A1=-2.*(m0+m1)*M1Coupling_[imode];
    A2= 4.*(m0+m1)*M1Coupling_[imode];
    B1=0.;
    B2=0.;
  }
  else {
    throw Exception() << "Unknown type of mode " << modeType_[imode] 
				 << " in RadiativeHeavyBaryonDecayer::"
				 << "halfHalfVectorCoupling()" << Exception::abortnow;
  }
}

void RadiativeHeavyBaryonDecayer::threeHalfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A1,Complex& A2,
							      Complex& A3,Complex& B1,
							      Complex& B2,
							      Complex& B3) const {
  useMe();
  if(modeType_[imode]==0) {
    A1=-0.5*(m0*m0 - m1*m1)*E1Coupling_[imode];
    A3=-    (m0+m1)*(m0+m1)*E1Coupling_[imode];    
    A2=0.;
    B1=0.;
    B2=0.;
    B3=0.;
  }
  else if(modeType_[imode]==1) {
    B1=-(m0+m1)*M1Coupling_[imode];
    B2= (m0+m1)*M1Coupling_[imode];
    B3=0.;
    A1=0.;
    A2=0.;
    A3=0.;
  }
  else {
    throw Exception() << "Unknown type of mode " << modeType_[imode] 
				 << " in RadiativeHeavyBaryonDecayer::"
				 << "threeHalfHalfVectorCoupling()" 
				 << Exception::abortnow;
  }
}

string RadiativeHeavyBaryonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2 or 3/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "Outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Outgoing particle with id " + std::to_string(out) + "does not have spin 1/2";
  // type of mode
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int itype = stoi(stype);
  // M1 coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy m1 = stof(stype)/MeV;
  // E1 coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy2 e1 = stof(stype)/MeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  modeType_.push_back(itype);
  M1Coupling_.push_back(m1);
  E1Coupling_.push_back(e1);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./RadiativeDoublyHeavyBaryonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RadiativeDoublyHeavyBaryonDecayer class.
//

#include "RadiativeDoublyHeavyBaryonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void RadiativeDoublyHeavyBaryonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxWeight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) {
      if(mode(ix)) maxWeight_.push_back(mode(ix)->maxWeight());
      else         maxWeight_.push_back(1.);
    }
  }
}

void RadiativeDoublyHeavyBaryonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // check the parameters are consistent
  unsigned int isize(incoming_.size());
  if(isize!=outgoing_.size() ||isize!=maxWeight_.size()||isize!=coupling_.size())
    throw InitException() << "Inconsistent parameters in "
			  << "RadiativeDoublyHeavyBaryonDecayer::doinit()" 
			  << Exception::abortnow;
  // the decay modes
  tPDPtr photon = getParticleData(ParticleID::gamma);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    tPDPtr    in  =  getParticleData(incoming_[ix]);
    tPDVector out = {getParticleData(outgoing_[ix]),photon};
    PhaseSpaceModePtr mode;
    if(in&&out[0]) {
      mode = new_ptr(PhaseSpaceMode(in,out,maxWeight_[ix]));
    }
    else {
      mode = PhaseSpaceModePtr();
    }
    addMode(mode);
  }
}

int RadiativeDoublyHeavyBaryonDecayer::modeNumber(bool & cc,tcPDPtr parent,
						  const tPDVector & children) const {
  int imode(-1);
  // must be two outgoing particles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id()),ibaryon;
  if(id1==ParticleID::gamma){ibaryon=id2;}
  else if(id2==ParticleID::gamma){ibaryon=id1;}
  else {return imode;}
  unsigned int ix(0);
  do {
    if(     id0== incoming_[ix]&&ibaryon== outgoing_[ix]) {
      imode=ix;
      cc=false;
    }
    else if(id0==-incoming_[ix]&&ibaryon==-outgoing_[ix]) {
      imode=ix;
      cc=true;
    }
    ++ix;
  }
  while(ix<incoming_.size()&&imode<0);
  return imode;
}

void RadiativeDoublyHeavyBaryonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1./GeV2) << incoming_ << outgoing_ << maxWeight_;
}

void RadiativeDoublyHeavyBaryonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1./GeV2) >> incoming_ >> outgoing_ >> maxWeight_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<RadiativeDoublyHeavyBaryonDecayer,Baryon1MesonDecayerBase>
describeHerwigRadiativeDoublyHeavyBaryonDecayer("Herwig::RadiativeDoublyHeavyBaryonDecayer", "HwBaryonDecay.so");

void RadiativeDoublyHeavyBaryonDecayer::Init() {

  static ClassDocumentation<RadiativeDoublyHeavyBaryonDecayer> documentation
    ("The RadiativeDoublyHeavyBaryonDecayer class is designed for the radiative decays of"
     " doubly heavy baryons.");

  static Command<RadiativeDoublyHeavyBaryonDecayer> interfaceSetUpDecayMode
    ("SetUpDecayMode",
     "Set up the particles (incoming, outgoing baryons, M1(1/GeV2) coupling and max weight for a decay",
     &RadiativeDoublyHeavyBaryonDecayer::setUpDecayMode, false);

}

void RadiativeDoublyHeavyBaryonDecayer::dataBaseOutput(ofstream & output,bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  for(unsigned int ix=0;ix<incoming_.size();++ix) {
    output << "do " << name() << ":SetUpDecayMode "
	   << incoming_[ix] << " " << outgoing_[ix] << " " << " " 
	   << coupling_[ix]*GeV2 << " " << maxWeight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}

void RadiativeDoublyHeavyBaryonDecayer::halfHalfVectorCoupling(int imode,Energy m0,Energy m1,
							       Energy,Complex& A1,
							       Complex& A2,Complex& B1,
							       Complex& B2) const {
  useMe();
  A1 = -0.5*sqr(m0+m1)*coupling_[imode];
  A2 =      sqr(m0+m1)*coupling_[imode];
  B1 = 0.;
  B2 = 0.;
}

void RadiativeDoublyHeavyBaryonDecayer::threeHalfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A1,Complex& A2,
							      Complex& A3,Complex& B1,
							      Complex& B2,
							      Complex& B3) const {
  useMe();
  A1 = 0.;
  A2 = 0.;
  A3 = 0.;
  B1 = 0.5*sqr(m0+m1)*coupling_[imode];
  B2 =- m0*   (m0+m1)*coupling_[imode];
  B3 =-    sqr(m0+m1)*coupling_[imode];
}

string RadiativeDoublyHeavyBaryonDecayer::setUpDecayMode(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long in = stoi(stype);
  tcPDPtr pData = getParticleData(in);
  if(!pData)
    return "Incoming particle with id " + std::to_string(in) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half && pData->iSpin()!=PDT::Spin3Half)
    return "Incoming particle with id " + std::to_string(in) + "does not have spin 1/2 or 3/2";
  // and outgoing particles
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  long out = stoi(stype);
  pData = getParticleData(out);
  if(!pData)
    return "Outgoing particle with id " + std::to_string(out) + "does not exist";
  if(pData->iSpin()!=PDT::Spin1Half)
    return "Outgoing particle with id " + std::to_string(out) + "does not have spin 1/2";
  // M1 coupling
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy2 g = stof(stype)/GeV2;
  // and the maximum weight
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double wgt = stof(stype);
  // store the information
  incoming_.push_back(in);
  outgoing_.push_back(out);
  coupling_.push_back(g);
  maxWeight_.push_back(wgt);
  // success
  return "";
}
#line 1 "./SU3BaryonDecupletOctetPhotonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonDecupletOctetPhotonDecayer class.
// 

#include "SU3BaryonDecupletOctetPhotonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void SU3BaryonDecupletOctetPhotonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<incomingB_.size();++ix) {
    tPDPtr    in  =  getParticleData(incomingB_[ix]);
    tPDVector out = {getParticleData(outgoingB_[ix]),
    		     getParticleData(ParticleID::gamma)};
    double wgtmax = maxweight_.size()>numberModes() ? 
      maxweight_[numberModes()] : 1.;
    PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

void SU3BaryonDecupletOctetPhotonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    maxweight_.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      maxweight_.push_back(mode(ix)->maxWeight());
  }
}

SU3BaryonDecupletOctetPhotonDecayer::SU3BaryonDecupletOctetPhotonDecayer() {
  // the coupling
  C_=1.0/GeV;
  // the relative parities of the two baryon multiplets
  parity_=true;
  // PDG codes for the various octet baryons
  proton_   = 2212;
  neutron_  = 2112;
  sigma0_   = 3212;
  sigmap_   = 3222;
  sigmam_   = 3112;
  lambda_   = 3122;
  xi0_      = 3322;
  xim_      = 3312;
  // PDG codes for the various decuplet baryons
  deltapp_  = 2224;
  deltap_   = 2214;
  delta0_   = 2114;
  deltam_   = 1114;
  sigmasp_  = 3224;
  sigmas0_  = 3214;
  sigmasm_  = 3114;
  omega_    = 3334;
  xism_     = 3314;
  xis0_     = 3324;
  // intermediates
  generateIntermediates(false);
}

int SU3BaryonDecupletOctetPhotonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  if(incomingB_.size()==0) setupModes(0); 
  // must be two outgoing particles
  if(children.size()!=2||(children[0]->id()!=ParticleID::gamma&&
			  children[1]->id()!=ParticleID::gamma)) return -1;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  int iout = id1==ParticleID::gamma ? id2 : id1;
  unsigned int ix(0);
  cc=false;
  int imode(-1);
  do {
    if(id0==incomingB_[ix]) {
      if(iout==outgoingB_[ix]) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-incomingB_[ix]) {
      if(iout==-outgoingB_[ix]) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<incomingB_.size()&&imode<0);
  return imode;
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(C_,1./GeV) << parity_ << proton_ << neutron_ << sigma0_ << sigmap_ 
     << sigmam_ << lambda_ << xi0_ << xim_ << deltapp_ << deltap_ << delta0_ << deltam_
     << sigmasp_ << sigmas0_ << sigmasm_ << omega_ << xism_ << xis0_ << incomingB_ 
     << outgoingB_ << maxweight_ << ounit(prefactor_,1./GeV);
}

void SU3BaryonDecupletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(C_,1./GeV) >> parity_ >> proton_ >> neutron_ >> sigma0_ >> sigmap_ 
     >> sigmam_ >> lambda_ >> xi0_ >> xim_ >> deltapp_ >> deltap_ >> delta0_ >> deltam_
     >> sigmasp_ >> sigmas0_ >> sigmasm_ >> omega_ >> xism_ >> xis0_ >> incomingB_ 
     >> outgoingB_ >> maxweight_ >> iunit(prefactor_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonDecupletOctetPhotonDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonDecupletOctetPhotonDecayer("Herwig::SU3BaryonDecupletOctetPhotonDecayer", "HwBaryonDecay.so");

void SU3BaryonDecupletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetPhotonDecayer> documentation
    ("The SU3BaryonDecupletOctetPhotonDecayer class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a photon.");

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,InvEnergy> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetPhotonDecayer::C_, 1.0/GeV, 1.0/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetPhotonDecayer::parity_, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "The multiplets have the same parity.",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "The multiplets have different parities.",
     false);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::proton_, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::neutron_, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmap_, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigma0_, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmam_, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::lambda_, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xi0_, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xim_, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::deltapp_, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::deltap_, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::delta0_, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::deltam_, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmasp_, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmas0_, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::sigmasm_, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::omega_, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xis0_, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetPhotonDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetPhotonDecayer::xism_, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetPhotonDecayer::maxweight_,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonDecupletOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const {
  A3=0.;
  B3=0.;
  if(parity_) {
    A1 = 0.;
    B1 = -prefactor_[imode]*(m0+m1);
    A2 = 0.;
    B2 = prefactor_[imode]*(m0+m1);
  }
  else {
    A1= prefactor_[imode]*(m0-m1);
    B1=0.;
    A2= prefactor_[imode]*(m0+m1);
    B2=0.;
  }
}

// set up the decay modes
void SU3BaryonDecupletOctetPhotonDecayer::setupModes(unsigned int iopt) const {
  if(incomingB_.size()!=0&&iopt==0) return;
  if(iopt==1) {
    outgoingB_.clear();
    incomingB_.clear();
  }
  vector<InvEnergy> factor;
  vector<int> intemp,outtemp;
  double ortw(1./sqrt(12.)),orr(1./sqrt(3.));
  // decays of the delta+
  intemp.push_back(deltap_);outtemp.push_back(proton_);
  factor.push_back(C_*orr);
  // decays of the delta0
  intemp.push_back(delta0_);outtemp.push_back(neutron_);
  factor.push_back(C_*orr);
  // sigma*+
  intemp.push_back(sigmasp_);outtemp.push_back(sigmap_);
  factor.push_back(-C_*orr);
  // sigma*0
  intemp.push_back(sigmas0_);outtemp.push_back(lambda_);
  factor.push_back(-C_*.5);
  intemp.push_back(sigmas0_);outtemp.push_back(sigma0_);
  factor.push_back(C_*ortw);
  // xi*0
  intemp.push_back(xis0_);outtemp.push_back(xi0_);
  factor.push_back(-C_*orr);
  // set up the modes
  tPDVector extpart(2);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()) {
	incomingB_.push_back(intemp[ix]);
	outgoingB_.push_back(outtemp[ix]);
	if(iopt==1) prefactor_.push_back(factor[ix]);
      }
    }
  }
}
void SU3BaryonDecupletOctetPhotonDecayer::dataBaseOutput(ofstream & output,
							 bool header) const {
  if(header) output << "update decayers set parameters=\""; 
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Ccoupling " << C_*GeV<< "\n";
  output << "newdef " << name() << ":Parity " << parity_<< "\n";
  output << "newdef " << name() << ":Proton " << proton_ << "\n";
  output << "newdef " << name() << ":Neutron " << neutron_ << "\n";
  output << "newdef " << name() << ":Sigma+ " << sigmap_ << "\n";
  output << "newdef " << name() << ":Sigma0 " << sigma0_ << "\n";
  output << "newdef " << name() << ":Sigma- " << sigmam_ << "\n";
  output << "newdef " << name() << ":Lambda " << lambda_ << "\n";
  output << "newdef " << name() << ":Xi0 " << xi0_ << "\n";
  output << "newdef " << name() << ":Xi- " << xim_ << "\n";
  output << "newdef " << name() << ":Delta++ " << deltapp_ << "\n";
  output << "newdef " << name() << ":Delta+ " << deltap_ << "\n";
  output << "newdef " << name() << ":Delta0 " << delta0_ << "\n";
  output << "newdef " << name() << ":Delta- " << deltam_ << "\n";
  output << "newdef " << name() << ":Sigma*+ " << sigmasp_ << "\n";
  output << "newdef " << name() << ":Sigma*0 " << sigmas0_ << "\n";
  output << "newdef " << name() << ":Sigma*- " << sigmasm_ << "\n";
  output << "newdef " << name() << ":Omega " << omega_ << "\n";
  output << "newdef " << name() << ":Xi*0 " << xis0_ << "\n";
  output << "newdef " << name() << ":Xi*- " << xism_ << "\n";
  for(unsigned int ix=0;ix<maxweight_.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << maxweight_[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
}
#line 1 "./SU3BaryonDecupletOctetScalarDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonDecupletOctetScalarDecayer class.
//

#include "SU3BaryonDecupletOctetScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonDecupletOctetScalarDecayer::SU3BaryonDecupletOctetScalarDecayer() {
  // couplings and off-shell parameter
  _c=1.5;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // the pion decay constant
  _fpi=92.4*MeV;
  // PDG codes for the various octet baryons
  _proton   = 2212;
  _neutron  = 2112;
  _sigma0   = 3212;
  _sigmap   = 3222;
  _sigmam   = 3112;
  _lambda   = 3122;
  _xi0      = 3322;
  _xim      = 3312;
  // PDG codes for the various decuplet baryons
  _deltapp  = 2224;
  _deltap   = 2214;
  _delta0   = 2114;
  _deltam   = 1114;
  _sigmasp  = 3224;
  _sigmas0  = 3214;
  _sigmasm  = 3114;
  _omega    = 3334;
  _xism     = 3314;
  _xis0     = 3324;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonDecupletOctetScalarDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_incomingB[ix]);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(_outgoingM[ix])};
    double wgtmax = _maxweight.size()>numberModes()? 
      _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode =
      new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

void SU3BaryonDecupletOctetScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonDecupletOctetScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc=false;
  do {
    if(id0==_incomingB[ix]) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incomingB[ix]) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) {
	imode=ix;
	cc=true;
      }
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonDecupletOctetScalarDecayer::
persistentOutput(PersistentOStream & os) const {
  os << _c << _parity << ounit(_fpi,GeV) << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonDecupletOctetScalarDecayer::
persistentInput(PersistentIStream & is, int) {
  is >> _c >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonDecupletOctetScalarDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonDecupletOctetScalarDecayer("Herwig::SU3BaryonDecupletOctetScalarDecayer", "HwBaryonDecay.so");

void SU3BaryonDecupletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonDecupletOctetScalarDecayer> documentation
    ("The SU3BaryonDecupletOctetScalarDecayer class is designed for the"
     " decay of an SU(3) decuplet baryon to an SU(3) octet baryon and a pseudoscalar"
     " meson from the lightest multiplet.");

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,double> interfaceCcoupling
    ("Ccoupling",
     "The C coupling for the decuplet decays.",
     &SU3BaryonDecupletOctetScalarDecayer::_c, 1.5, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonDecupletOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonDecupletOctetScalarDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonDecupletOctetScalarDecayer::_fpi, MeV, 92.4*MeV, 50.0*MeV, 150.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltapp, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltap, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_delta0, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_deltam, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasp, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmas0, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_sigmasm, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_omega, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xis0, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonDecupletOctetScalarDecayer,long> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonDecupletOctetScalarDecayer::_xism, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonDecupletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonDecupletOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonDecupletOctetScalarDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex& A, Complex& B) const {
  if(_parity) {
    A = _prefactor[imode]*(m0+m1);
    B = 0.;
  }
  else {
    A = 0.;
    B = _prefactor[imode]*(m0+m1);
  }
}

// set up the decay modes
void SU3BaryonDecupletOctetScalarDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0) return;
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
    _outgoingM.clear();
  }
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.)),orr(1./sqrt(3.));
  // decays of the delta++
  intemp.push_back(_deltapp);outtemp.push_back(_proton);mestemp.push_back(211);
  factor.push_back(_c);
  intemp.push_back(_deltapp);outtemp.push_back(_sigmap);mestemp.push_back(321);
  factor.push_back(-_c);
  // decays of the delta+
  intemp.push_back(_deltap);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back(_c*orr);
  intemp.push_back(_deltap);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_deltap);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_c*orr);
  // decays of the delta0
  intemp.push_back(_delta0);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back(_c*orr);
  intemp.push_back(_delta0);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(_c*rt*orr);
  intemp.push_back(_delta0);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_c*orr);
  // decays of the delta-
  intemp.push_back(_deltam);outtemp.push_back(_neutron);mestemp.push_back(-211);
  factor.push_back(-_c);
  intemp.push_back(_deltam);outtemp.push_back(_sigmam);mestemp.push_back(311);
  factor.push_back(_c);
  // sigma*+
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasp);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasp);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(_c*ort);
  intemp.push_back(_sigmasp);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(_c*ort);
  // sigma*0
  intemp.push_back(_sigmas0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmas0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(_c*ort);
  intemp.push_back(_sigmas0);outtemp.push_back(_sigma0);mestemp.push_back(221);
  factor.push_back(_c*ort);
  // sigma*-
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_c*ors);
  intemp.push_back(_sigmasm);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_c*orr);
  intemp.push_back(_sigmasm);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(_c*ort);
  intemp.push_back(_sigmasm);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(_c*ort);
  // xi*0
  intemp.push_back(_xis0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_c*orr);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_xis0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_c*orr);
  intemp.push_back(_xis0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(_c*ors);
  intemp.push_back(_xis0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(_c*ort);
  intemp.push_back(_xis0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(_c*ort);
  // xi*-
  intemp.push_back(_xism);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_c*orr);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(_c*ors);
  intemp.push_back(_xism);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_c*orr);
  intemp.push_back(_xism);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(_c*ors);
  intemp.push_back(_xism);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(_c*ort);
  intemp.push_back(_xism);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(_c*ort);
  // set up the modes
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0&&mestemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      extpart[2]=getParticleData(mestemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin()) {
	_incomingB.push_back(intemp[ix]);
	_outgoingB.push_back(outtemp[ix]);
	_outgoingM.push_back(mestemp[ix]);
	if(iopt==1) _prefactor.push_back(factor[ix]/_fpi);
      }
    }
  }
}

void SU3BaryonDecupletOctetScalarDecayer::dataBaseOutput(ofstream & output, 
							 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Ccoupling " << _c << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":Proton " << _proton << "\n";
  output << "newdef " << name() << ":Neutron " << _neutron << "\n";
  output << "newdef " << name() << ":Sigma+ " << _sigmap << "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Sigma- " << _sigmam << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":Xi0 " << _xi0 << "\n";
  output << "newdef " << name() << ":Xi- " << _xim << "\n";
  output << "newdef " << name() << ":Delta++ " << _deltapp << "\n";
  output << "newdef " << name() << ":Delta+ " << _deltap << "\n";
  output << "newdef " << name() << ":Delta0 " << _delta0 << "\n";
  output << "newdef " << name() << ":Delta- " << _deltam << "\n";
  output << "newdef " << name() << ":Sigma*+ " << _sigmasp << "\n";
  output << "newdef " << name() << ":Sigma*0 " << _sigmas0 << "\n";
  output << "newdef " << name() << ":Sigma*- " << _sigmasm << "\n";
  output << "newdef " << name() << ":Omega " << _omega << "\n";
  output << "newdef " << name() << ":Xi*0 " << _xis0 << "\n";
  output << "newdef " << name() << ":Xi*- " << _xism << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix)
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SU3BaryonOctetDecupletScalarDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonOctetDecupletScalarDecayer class.
//

#include "SU3BaryonOctetDecupletScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonOctetDecupletScalarDecayer::SU3BaryonOctetDecupletScalarDecayer() {
  // couplings and off-shell parameter
  _c=1.35;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // the pion decay constant
  _fpi=130.7*MeV;
  // PDG codes for the various octet baryons
  _proton   =  12212;
  _neutron  =  12112;
  _sigma0   =  13212;
  _sigmap   =  13222;
  _sigmam   =  13112;
  _lambda   =  23122;
  _xi0      =  13322;
  _xim      =  13312;
  // PDG codes for the various decuplet baryons
  _deltapp  = 2224;
  _deltap   = 2214;
  _delta0   = 2114;
  _deltam   = 1114;
  _sigmasp  = 3224;
  _sigmas0  = 3214;
  _sigmasm  = 3114;
  _omega    = 3334;
  _xism     = 3314;
  _xis0     = 3324;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonOctetDecupletScalarDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_incomingB[ix]);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(_outgoingM[ix])};
    double wgtmax = _maxweight.size()>numberModes() ? 
      _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode =
      new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

void SU3BaryonOctetDecupletScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonOctetDecupletScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_incomingB[ix]) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incomingB[ix]) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) {
	imode=ix;
	cc=true;
      }
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonOctetDecupletScalarDecayer::
persistentOutput(PersistentOStream & os) const {
  os << _c << _parity << ounit(_fpi,GeV) << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _deltapp << _deltap << _delta0 << _deltam
     << _sigmasp << _sigmas0 << _sigmasm << _omega << _xism << _xis0 << _incomingB 
     << _outgoingB << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonOctetDecupletScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _c >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _deltapp >> _deltap >> _delta0 >> _deltam
     >> _sigmasp >> _sigmas0 >> _sigmasm >> _omega >> _xism >> _xis0 >> _incomingB 
     >> _outgoingB >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonOctetDecupletScalarDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonOctetDecupletScalarDecayer("Herwig::SU3BaryonOctetDecupletScalarDecayer", "HwBaryonDecay.so");

void SU3BaryonOctetDecupletScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetDecupletScalarDecayer> documentation
    ("The SU3BaryonOctetDecupletScalarDecayer class is designed for the"
     " decay of an SU(3) octet baryon to an SU(3) decuplet baryon and a pseudoscalar"
     " meson from the lightest multiplet.");

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,double> interfaceCoupling
    ("Coupling",
     "The coupling for the decuplet decays.",
     &SU3BaryonOctetDecupletScalarDecayer::_c, 1.35, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonOctetDecupletScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetDecupletScalarDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonOctetDecupletScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the proton-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_proton, 12212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the neutron-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_neutron, 12112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the Sigma+-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmap, 13222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the Sigma0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigma0, 13212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the Sigma--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmam, 13112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the Lambda-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_lambda, 23122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the Xi0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xi0, 13322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the Xi--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xim, 13312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltapp
    ("Delta++",
     "The PDG code for the Delta++ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltapp, 2224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltap
    ("Delta+",
     "The PDG code for the Delta+ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltap, 2214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDelta0
    ("Delta0",
     "The PDG code for the Delta0 like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_delta0, 2114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceDeltam
    ("Delta-",
     "The PDG code for the Delta- like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_deltam, 1114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmasp
    ("Sigma*+",
     "The PDG code for the Sigma*+ like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmasp, 3224, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmas0
    ("Sigma*0",
     "The PDG code for the Sigma*0 like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmas0, 3214, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceSigmasm
    ("Sigma*-",
     "The PDG code for the Sigma*- like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_sigmasm, 3114, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceOmega
    ("Omega",
     "The PDG code for the Omega like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_omega, 3334, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXis0
    ("Xi*0",
     "The PDG code for the Xi*0-like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xis0, 3324, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetDecupletScalarDecayer,int> interfaceXism
    ("Xi*-",
     "The PDG code for the Xi*--like baryon.",
     &SU3BaryonOctetDecupletScalarDecayer::_xism, 3314, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonOctetDecupletScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetDecupletScalarDecayer::_maxweight,
     0, 0, 0, 0., 2000., false, false, true);

}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonOctetDecupletScalarDecayer
::halfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
			      Complex& A, Complex& B) const {
  if(_parity) {
    A=_prefactor[imode]*(m0+m1);
    B=0.;
  }
  else {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
}

// couplings for spin-3/2 to spin-3/2 spin-0
void SU3BaryonOctetDecupletScalarDecayer::
threeHalfThreeHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
				 Complex& A1,Complex& A2,
				 Complex& B1,Complex& B2) const {
  A2=0.;B2=0.;
  if(_parity) {
    A1=0.;
    B1=_prefactor[imode]*(m0+m1);
  }
  else {
    A1=_prefactor[imode]*(m0-m1);
    B1=0.;
  }
}

// set up the decay modes
void SU3BaryonOctetDecupletScalarDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
    _outgoingM.clear();
  }
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.)),orr(1./sqrt(3.));
  // decays to the delta++
  outtemp.push_back(_deltapp);intemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back(_c);
  outtemp.push_back(_deltapp);intemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(-_c);
  // decays to the delta+
  outtemp.push_back(_deltap);intemp.push_back(_neutron);mestemp.push_back(-211);
  factor.push_back(_c*orr);
  outtemp.push_back(_deltap);intemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(_c*rt*orr);
  outtemp.push_back(_deltap);intemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(_c*rt*orr);
  outtemp.push_back(_deltap);intemp.push_back(_sigmap);mestemp.push_back(-311);
  factor.push_back(_c*orr);
  // decays to the delta0
  outtemp.push_back(_delta0);intemp.push_back(_proton);mestemp.push_back(211);
  factor.push_back(_c*orr);
  outtemp.push_back(_delta0);intemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(_c*rt*orr);
  outtemp.push_back(_delta0);intemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(_c*rt*orr);
  outtemp.push_back(_delta0);intemp.push_back(_sigmam);mestemp.push_back(-321);
  factor.push_back(_c*orr);
  // decays to the delta-
  outtemp.push_back(_deltam);intemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back(-_c);
  outtemp.push_back(_deltam);intemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_c);
  // decays to the sigma*+
  outtemp.push_back(_sigmasp);intemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmasp);intemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmasp);intemp.push_back(_proton);mestemp.push_back(311);
  factor.push_back(_c*orr);
  outtemp.push_back(_sigmasp);intemp.push_back(_xi0);mestemp.push_back(-321);
  factor.push_back(_c*orr);
  outtemp.push_back(_sigmasp);intemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(_c*ort);
  outtemp.push_back(_sigmasp);intemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(_c*ort);
  // decays to the sigma*0
  outtemp.push_back(_sigmas0);intemp.push_back(_neutron);mestemp.push_back(311);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_proton);mestemp.push_back(321);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_xim);mestemp.push_back(-321);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_xi0);mestemp.push_back(-311);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_sigmam);mestemp.push_back(-211);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_sigmap);mestemp.push_back(211);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmas0);intemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(_c*ort);
  outtemp.push_back(_sigmas0);intemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(_c*ort);
  // decays to the sigma*-
  outtemp.push_back(_sigmasm);intemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmasm);intemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(_c*ors);
  outtemp.push_back(_sigmasm);intemp.push_back(_neutron);mestemp.push_back(321);
  factor.push_back(_c*orr);
  outtemp.push_back(_sigmasm);intemp.push_back(_xim);mestemp.push_back(-311);
  factor.push_back(_c*orr);
  outtemp.push_back(_sigmasm);intemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(_c*ort);
  outtemp.push_back(_sigmasm);intemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(_c*ort);
  // decays to the xi*0
  outtemp.push_back(_xis0);intemp.push_back(_xim);mestemp.push_back(-211);
  factor.push_back(_c*orr);
  outtemp.push_back(_xis0);intemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(_c*ors);
  outtemp.push_back(_xis0);intemp.push_back(_sigmap);mestemp.push_back(321);
  factor.push_back(_c*orr);
  outtemp.push_back(_xis0);intemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(_c*ors);
  outtemp.push_back(_xis0);intemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(_c*ort);
  outtemp.push_back(_xis0);intemp.push_back(_lambda);mestemp.push_back(311);
  factor.push_back(_c*ort);
  // decays to the xi*-
  outtemp.push_back(_xism);intemp.push_back(_xi0);mestemp.push_back(211);
  factor.push_back(_c*orr);
  outtemp.push_back(_xism);intemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(_c*ors);
  outtemp.push_back(_xism);intemp.push_back(_sigmam);mestemp.push_back(311);
  factor.push_back(_c*orr);
  outtemp.push_back(_xism);intemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(_c*ors);
  outtemp.push_back(_xism);intemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(_c*ort);
  outtemp.push_back(_xism);intemp.push_back(_lambda);mestemp.push_back(321);
  factor.push_back(_c*ort);
  // set up the modes
  int inspin,outspin;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<outtemp.size();++ix) {
    if(outtemp[ix]!=0&&intemp[ix]!=0&&mestemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      extpart[2]=getParticleData(mestemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin()) {
	_incomingB.push_back(intemp[ix]);
	_outgoingB.push_back(outtemp[ix]);
	_outgoingM.push_back(mestemp[ix]);
	if(iopt==1) {
	  inspin  = extpart[0]->iSpin();
	  outspin = extpart[1]->iSpin();
	  if(inspin==2&&outspin==4)      _prefactor.push_back(factor[ix]/_fpi);
	  else if(inspin==4&&outspin==4) _prefactor.push_back(factor[ix]/_fpi);
	}
      }
    }
  }
}

void SU3BaryonOctetDecupletScalarDecayer::dataBaseOutput(ofstream & output,
							 bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Coupling " << _c<< "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":Proton " << _proton << "\n";
  output << "newdef " << name() << ":Neutron " << _neutron << "\n";
  output << "newdef " << name() << ":Sigma+ " << _sigmap << "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Sigma- " << _sigmam << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":Xi0 " << _xi0 << "\n";
  output << "newdef " << name() << ":Xi- " << _xim << "\n";
  output << "newdef " << name() << ":Delta++ " << _deltapp << "\n";
  output << "newdef " << name() << ":Delta+ " << _deltap << "\n";
  output << "newdef " << name() << ":Delta0 " << _delta0 << "\n";
  output << "newdef " << name() << ":Delta- " << _deltam << "\n";
  output << "newdef " << name() << ":Sigma*+ " << _sigmasp << "\n";
  output << "newdef " << name() << ":Sigma*0 " << _sigmas0 << "\n";
  output << "newdef " << name() << ":Sigma*- " << _sigmasm << "\n";
  output << "newdef " << name() << ":Omega " << _omega << "\n";
  output << "newdef " << name() << ":Xi*0 " << _xis0 << "\n";
  output << "newdef " << name() << ":Xi*- " << _xism << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SU3BaryonOctetOctetPhotonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3OctetOctetPhotonDecayer class.
//

#include "SU3BaryonOctetOctetPhotonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonOctetOctetPhotonDecayer::SU3BaryonOctetOctetPhotonDecayer() {
  // default values of the parameters
  // these are the values for first excited multiplet
  // the couplings of the anticommutator and communtator terms
  _lf=-0.009/GeV;
  _ld=-0.024/GeV;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // PDG codes for the various ground state baryons
  _proton   = 2212;
  _neutron  = 2112;
  _sigma0   = 3212;
  _sigmap   = 3222;
  _sigmam   = 3112;
  _lambda   = 3122;
  _xi0      = 3322;
  _xim      = 3312;
  // PDG codes for the various excited baryons
  _eproton  = 12212;
  _eneutron = 12112;
  _esigma0  = 13212;
  _esigmap  = 13222;
  _esigmam  = 13112;
  _elambda  = 23122;
  _exi0     = 13322;
  _exim     = 13312;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonOctetOctetPhotonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_incomingB[ix]);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(ParticleID::gamma)};
    double wgtmax = _maxweight.size()>numberModes() ? 
      _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode =
      new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
    // testing code
//     Energy MR = extpart[0]->mass();
//     Energy MB = extpart[1]->mass();
//     Energy kg = 0.5*(sqr(MR)-sqr(MB))/MR;
//     Energy2 Tji = 128.*sqr((sqr(MR)-sqr(MB))/2.*_prefactor[ix]/4.);
//     Energy width = 1./8./Constants::pi/sqr(MR)*kg*Tji;
//     generator()->log() << "testing " << extpart[0]->PDGName() << "->"
// 		       << extpart[1]->PDGName() << " " 
// 		       << extpart[2]->PDGName() << " " << width/MeV
// 		       << " MeV\n";
  }
}

void SU3BaryonOctetOctetPhotonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonOctetOctetPhotonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id()),id1(children[0]->id()),id2(children[1]->id()),iout;
  if(id1==ParticleID::gamma){iout=id2;}
  else if(id2==ParticleID::gamma){iout=id1;}
  else{return imode;}
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_incomingB[ix]) {
      if(iout==_outgoingB[ix]) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incomingB[ix]) {
      if(iout==-_outgoingB[ix]) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonOctetOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_lf,1./GeV) << ounit(_ld,1./GeV) <<  _parity << _proton << _neutron 
     << _sigma0 << _sigmap << _sigmam << _lambda << _xi0 << _xim << _eproton 
     << _eneutron << _esigma0 << _esigmap << _esigmam << _elambda << _exi0 << _exim 
     << _incomingB << _outgoingB << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonOctetOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_lf,1./GeV) >> iunit(_ld,1./GeV) >>  _parity >> _proton >> _neutron 
     >> _sigma0 >> _sigmap >> _sigmam >> _lambda >> _xi0 >> _xim >> _eproton 
     >> _eneutron >> _esigma0 >> _esigmap >> _esigmam >> _elambda >> _exi0 >> _exim 
     >> _incomingB >> _outgoingB >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonOctetOctetPhotonDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonOctetOctetPhotonDecayer("Herwig::SU3BaryonOctetOctetPhotonDecayer", "HwBaryonDecay.so");

void SU3BaryonOctetOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetOctetPhotonDecayer> documentation
    ("The SU3BaryonOctetOctetPhotonDecayer class is designed for the "
     "radiative decay of an octet baryon to another octet baryon.");

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,InvEnergy> interfaceFcoupling
    ("Fcoupling",
     "The F coupling of the baryon resonances",
     &SU3BaryonOctetOctetPhotonDecayer::_lf, 1./GeV, -0.009/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,InvEnergy> interfaceDcoupling
    ("Dcoupling",
     "The D coupling of the baryon resonances",
     &SU3BaryonOctetOctetPhotonDecayer::_ld, 1./GeV, -0.024/GeV, -10.0/GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonOctetOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetOctetPhotonDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedProton
    ("ExcitedProton",
     "The PDG code for the heavier proton-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_eproton, 12212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedNeutron
    ("ExcitedNeutron",
     "The PDG code for the heavier neutron-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_eneutron, 12112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigmap
    ("ExcitedSigma+",
     "The PDG code for the heavier Sigma+-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigmap, 13222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigma0
    ("ExcitedSigma0",
     "The PDG code for the heavier Sigma0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigma0, 13212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedSigmam
    ("ExcitedSigma-",
     "The PDG code for the heavier Sigma--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_esigmam, 13112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_elambda, 23122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedXi0
    ("ExcitedXi0",
     "The PDG code for the heavier Xi0-like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_exi0, 13322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetPhotonDecayer,int> interfaceExcitedXim
    ("ExcitedXi-",
     "The PDG code for the heavier Xi--like baryon.",
     &SU3BaryonOctetOctetPhotonDecayer::_exim, 13312, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonOctetOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetOctetPhotonDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-1
void SU3BaryonOctetOctetPhotonDecayer::halfHalfVectorCoupling(int imode,Energy m0,
							      Energy m1, Energy,
							      Complex&A1,
							      Complex&A2,Complex&B1,
							      Complex&B2) const {
  if(_parity) {
    A1=    _prefactor[imode]*(m0+m1);
    B1=0.;
    A2=-2.*_prefactor[imode]*(m0+m1);
    B2=0.;
  }
  else {
    A1=0.;
    B1=    _prefactor[imode]*(m1-m0);
    A2=0.;
    B2=-2.*_prefactor[imode]*(m0+m1);
  }
}


// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonOctetOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1, Energy,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const {
  if(_parity) {
    A1=0.;
    B1=-_prefactor[imode]*(m0+m1);
    A2=0.;
    B2= _prefactor[imode]*(m0+m1);
  }
  else {
    A1= _prefactor[imode]*(m0-m1);
    B1=0.;
    A2= _prefactor[imode]*(m0+m1);
    B2=0.;
  }
  A3=0.;
  B3=0.;
}

// set up the decay modes
void SU3BaryonOctetOctetPhotonDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0) return;
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
  }
  // set up for the various different decay modes
  vector<InvEnergy> factor;
  vector<int> intemp,outtemp;
  // decays of the excited proton
  intemp.push_back(_eproton);outtemp.push_back(_proton);
  factor.push_back(_lf+_ld/3.);
  // decays of the excited neutron
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);
  factor.push_back(-2.*_ld/3.);
  // decays of the excited lambda
  intemp.push_back(_elambda);outtemp.push_back(_sigma0);
  factor.push_back(2.*_ld/sqrt(3.));
  intemp.push_back(_elambda);outtemp.push_back(_lambda);
  factor.push_back(-_ld/3.);
  // decays of the excited sigma+
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);
  factor.push_back(_lf+_ld/3.);
  // decays of the excited sigma0
  intemp.push_back(_esigma0);outtemp.push_back(_sigma0);
  factor.push_back(_lf/3.);
  intemp.push_back(_esigma0);outtemp.push_back(_lambda);
  factor.push_back(2.*_ld/sqrt(3.));
  // decays of the excited simga-
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);
  factor.push_back(-_ld-_lf/3.);
  // decays of the excited xi-
  intemp.push_back(_exim);outtemp.push_back(_xim);
  factor.push_back(_ld/3.-_lf);
  // decays of the excited xi0
  intemp.push_back(_exi0);outtemp.push_back(_xi0);
  factor.push_back(-2.*_ld/3.);
  int inspin,outspin;
  tPDVector extpart(2);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()) {
	_incomingB.push_back(intemp[ix]);
	_outgoingB.push_back(outtemp[ix]);
	if(iopt==1) {
	  inspin  = extpart[0]->iSpin();
	  outspin = extpart[1]->iSpin();
	  factor[ix] *=2.;
	  if(inspin==2&&outspin==2)      _prefactor.push_back(2.*factor[ix]);
	  else if(inspin==4&&outspin==2) _prefactor.push_back(factor[ix]);
	  else throw Exception() 
	    << "Invalid combination of spins in "
	    << "SU3BaryonOctetOctetPhotonDecayer::" 
	    << "setupModes()" << Exception::abortnow;}
      }
    }
  }
}

void SU3BaryonOctetOctetPhotonDecayer::dataBaseOutput(ofstream & output,
						      bool header) const {
  if(header) output << "update decayers set parameters=\""; 
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Fcoupling " << _lf*GeV << "\n";
  output << "newdef " << name() << ":Dcoupling " << _ld*GeV << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Proton " << _proton << "\n";
  output << "newdef " << name() << ":Neutron " << _neutron << "\n";
  output << "newdef " << name() << ":Sigma+ " << _sigmap << "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Sigma- " << _sigmam << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":Xi0 " << _xi0 << "\n";
  output << "newdef " << name() << ":Xi- " << _xim << "\n"; 
  output << "newdef " << name() << ":ExcitedProton " << _eproton << "\n";
  output << "newdef " << name() << ":ExcitedNeutron " << _eneutron << "\n";
  output << "newdef " << name() << ":ExcitedSigma+ " << _esigmap << "\n";
  output << "newdef " << name() << ":ExcitedSigma0 " << _esigma0 << "\n";
  output << "newdef " << name() << ":ExcitedSigma- " << _esigmam << "\n";
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  output << "newdef " << name() << ":ExcitedXi0 " << _exi0 << "\n";
  output << "newdef " << name() << ":ExcitedXi- " << _exim << "\n"; 
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SU3BaryonOctetOctetScalarDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonOctetOctetScalarDecayer class.
//

#include "SU3BaryonOctetOctetScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;
using namespace ThePEG;

SU3BaryonOctetOctetScalarDecayer::SU3BaryonOctetOctetScalarDecayer() {
  // default values of the parameters
  // these are the values for first excited multiplet
  // the couplings of the anticommutator and communtator terms
  _sf= 0.11;
  _sd= 0.60;
  // the relative parities of the two baryon multiplets
  _parity=true;
  // the pion decay constant
  _fpi=130.7*MeV;
  // PDG codes for the various ground state baryons
  _proton   = 2212;
  _neutron  = 2112;
  _sigma0   = 3212;
  _sigmap   = 3222;
  _sigmam   = 3112;
  _lambda   = 3122;
  _xi0      = 3322;
  _xim      = 3312;
  // PDG codes for the various excited baryons
  _eproton  = 12212;
  _eneutron = 12112;
  _esigma0  = 13212;
  _esigmap  = 13222;
  _esigmam  = 13112;
  _elambda  = 23122;
  _exi0     = 13322;
  _exim     = 13312;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonOctetOctetScalarDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<_incomingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_incomingB[ix]);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(_outgoingM[ix])};
    double wgtmax = 
      _maxweight.size()>numberModes() ? _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode = 
      new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
    // testing code
//     Energy MR=extpart[0]->mass();
//     Energy MB=extpart[1]->mass();
//     Energy Mp=extpart[2]->mass();
//     Energy kp = 0.5/MR*sqrt((sqr(MR)-sqr(MB+Mp))*(sqr(MR)-sqr(MB-Mp)));
//     Energy width;
//     if(_parity) {
//       width = sqr((MR+MB)/MR)*kp/(8.*Constants::pi)*
// 	sqr(_prefactor[ix])*(sqr(MR-MB)-sqr(Mp));
//     }
//     else {
//     width  = sqr((MR-MB)/MR)*kp/(8.*Constants::pi)*
// 	sqr(_prefactor[ix])*(sqr(MR+MB)-sqr(Mp));
//     } 
//     generator()->log() << "Partial Width for "
// 		       << extpart[0]->PDGName() << "->"
// 		       << extpart[1]->PDGName() << " "
// 		       << extpart[2]->PDGName() << " "
// 		       << width/GeV << "\n";
  }
}
void SU3BaryonOctetOctetScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix) 
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonOctetOctetScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_incomingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_incomingB[ix]) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_incomingB[ix]) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) {
	imode=ix;
	cc=true;
      }
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_incomingB.size()&&imode<0);
  return imode;
}

void SU3BaryonOctetOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _sf << _sd << _parity << ounit(_fpi,GeV) << _proton << _neutron 
     << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _eproton << _eneutron << _esigma0 
     << _esigmap << _esigmam << _elambda << _exi0 << _exim << _incomingB << _outgoingB 
     << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonOctetOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _sf >> _sd >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron 
     >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _eproton >> _eneutron >> _esigma0 
     >> _esigmap >> _esigmam >> _elambda >> _exi0 >> _exim >> _incomingB >> _outgoingB 
     >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonOctetOctetScalarDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonOctetOctetScalarDecayer("Herwig::SU3BaryonOctetOctetScalarDecayer", "HwBaryonDecay.so");

void SU3BaryonOctetOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonOctetOctetScalarDecayer> documentation
    ("The SU3BaryonOctetOctetScalarDecayer class is designed for the"
     " decay of excited baryon resonances assuming SU(3) symmetry");

  static Parameter<SU3BaryonOctetOctetScalarDecayer,double> interfaceFcoupling
    ("Fcoupling",
     "The F coupling of the baryon resonances",
     &SU3BaryonOctetOctetScalarDecayer::_sf, 0.60, -20.0, 20.0,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,double> interfaceDcoupling
    ("Dcoupling",
     "The D coupling of the baryon resonances",
     &SU3BaryonOctetOctetScalarDecayer::_sd, 0.11, -20.0, 20.0,
     false, false, true);

  static Switch<SU3BaryonOctetOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonOctetOctetScalarDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonOctetOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedProton
    ("ExcitedProton",
     "The PDG code for the heavier proton-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_eproton, 12212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedNeutron
    ("ExcitedNeutron",
     "The PDG code for the heavier neutron-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_eneutron, 12112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigmap
    ("ExcitedSigma+",
     "The PDG code for the heavier Sigma+-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigmap, 13222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigma0
    ("ExcitedSigma0",
     "The PDG code for the heavier Sigma0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigma0, 13212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedSigmam
    ("ExcitedSigma-",
     "The PDG code for the heavier Sigma--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_esigmam, 13112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_elambda, 23112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedXi0
    ("ExcitedXi0",
     "The PDG code for the heavier Xi0-like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_exi0, 13322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonOctetOctetScalarDecayer,int> interfaceExcitedXim
    ("ExcitedXi-",
     "The PDG code for the heavier Xi--like baryon.",
     &SU3BaryonOctetOctetScalarDecayer::_exim, 13312, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonOctetOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonOctetOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 10000., false, false, true);
}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonOctetOctetScalarDecayer::halfHalfScalarCoupling(int imode,Energy m0,
							      Energy m1,Energy,
							      Complex& A,
							      Complex& B) const {
  if(_parity) {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
  else {
    A=_prefactor[imode]*(m0-m1);
    B=0.;
  }
}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonOctetOctetScalarDecayer::threeHalfHalfScalarCoupling(int imode,Energy m0,
								   Energy m1,Energy,
								   Complex& A,
								   Complex& B) const {
  if(_parity) {
    A=_prefactor[imode]*(m0+m1);
    B=0.;
  }
  else {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
}

// set up the decay modes
void SU3BaryonOctetOctetScalarDecayer::setupModes(unsigned int iopt) const {
  if(_incomingB.size()!=0&&iopt==0){return;}
  if(iopt==1) {
    _outgoingB.clear();
    _incomingB.clear();
    _outgoingM.clear();
  }
  // set up for the various different decay modes
  vector<double> factor;
  vector<int> intemp,outtemp,mestemp;
  double ort(1./sqrt(2.)),ors(1./sqrt(6.)),rt(sqrt(2.));
  // decays of the excited proton
  intemp.push_back(_eproton);outtemp.push_back(_neutron);mestemp.push_back(211);
  factor.push_back((_sd+_sf));
  intemp.push_back(_eproton);outtemp.push_back(_proton);mestemp.push_back(111);
  factor.push_back(-ort*(_sd+_sf));
  intemp.push_back(_eproton);outtemp.push_back(_sigmap);mestemp.push_back(311);
  factor.push_back(_sd-_sf);
  intemp.push_back(_eproton);outtemp.push_back(_sigma0);mestemp.push_back(321);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_eproton);outtemp.push_back(_proton);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_eproton);outtemp.push_back(_lambda);mestemp.push_back(321);
  factor.push_back(-ors*(_sd+3.*_sf));
  // decays of the excited neutron
  intemp.push_back(_eneutron);outtemp.push_back(_proton);mestemp.push_back(-211);
  factor.push_back((_sd+_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);mestemp.push_back(111);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_sigmam);mestemp.push_back(321);
  factor.push_back(_sd-_sf);
  intemp.push_back(_eneutron);outtemp.push_back(_sigma0);mestemp.push_back(311);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_eneutron);outtemp.push_back(_neutron);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_eneutron);outtemp.push_back(_lambda);mestemp.push_back(311);
  factor.push_back(-ors*(_sd+3.*_sf));
  // decays of the excited lambda
  intemp.push_back(_elambda);outtemp.push_back(_sigma0);mestemp.push_back(111);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_lambda);mestemp.push_back(221);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_elambda);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_elambda);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_elambda);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(-ors*(3.*_sf+_sd));
  intemp.push_back(_elambda);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(-ors*(3.*_sf+_sd));
  // decays of the excited sigma+
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);mestemp.push_back(111);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_sigma0);mestemp.push_back(211);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_xi0);mestemp.push_back(321);
  factor.push_back(_sd+_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_proton);mestemp.push_back(-311);
  factor.push_back(_sd-_sf);
  intemp.push_back(_esigmap);outtemp.push_back(_lambda);mestemp.push_back(211);
  factor.push_back(2.*ors*_sd);
  intemp.push_back(_esigmap);outtemp.push_back(_sigmap);mestemp.push_back(221);
  factor.push_back(2.*ors*_sd);
  // decays of the excited sigma0
  intemp.push_back(_esigma0);outtemp.push_back(_sigmam);mestemp.push_back(211);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigma0);outtemp.push_back(_sigmap);mestemp.push_back(-211);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigma0);outtemp.push_back(_xim);mestemp.push_back(321);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_xi0);mestemp.push_back(311);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_proton);mestemp.push_back(-321);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_neutron);mestemp.push_back(-311);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_esigma0);outtemp.push_back(_sigma0);mestemp.push_back(221);
  factor.push_back(2.*ors*_sd);
  intemp.push_back(_esigma0);outtemp.push_back(_lambda);mestemp.push_back(111);
  factor.push_back(2.*ors*_sd);
  // decays of the excited simga-
  intemp.push_back(_esigmam);outtemp.push_back(_sigma0);mestemp.push_back(-211);
  factor.push_back(rt*_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);mestemp.push_back(111);
  factor.push_back(-rt*_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_xim);mestemp.push_back(311);
  factor.push_back(_sd+_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_neutron);mestemp.push_back(-321);
  factor.push_back(_sd-_sf);
  intemp.push_back(_esigmam);outtemp.push_back(_lambda);mestemp.push_back(-211);
  factor.push_back(2.*_sd*ors);
  intemp.push_back(_esigmam);outtemp.push_back(_sigmam);mestemp.push_back(221);
  factor.push_back(2.*_sd*ors);
  // decays of the excited xi-
  intemp.push_back(_exim);outtemp.push_back(_sigmam);mestemp.push_back(-311);
  factor.push_back(_sd+_sf);
  intemp.push_back(_exim);outtemp.push_back(_sigma0);mestemp.push_back(-321);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_exim);outtemp.push_back(_xi0);mestemp.push_back(-211);
  factor.push_back(_sd-_sf);
  intemp.push_back(_exim);outtemp.push_back(_xim);mestemp.push_back(111);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_exim);outtemp.push_back(_lambda);mestemp.push_back(-321);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_exim);outtemp.push_back(_xim);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf+_sd));
  // decays of the excited xi0
  intemp.push_back(_exi0);outtemp.push_back(_sigmap);mestemp.push_back(-321);
  factor.push_back(_sd+_sf);
  intemp.push_back(_exi0);outtemp.push_back(_sigma0);mestemp.push_back(-311);
  factor.push_back(ort*(_sd+_sf));
  intemp.push_back(_exi0);outtemp.push_back(_xim);mestemp.push_back(211);
  factor.push_back(_sd-_sf);
  intemp.push_back(_exi0);outtemp.push_back(_xi0);mestemp.push_back(111);
  factor.push_back(ort*(_sd-_sf));
  intemp.push_back(_exi0);outtemp.push_back(_lambda);mestemp.push_back(-311);
  factor.push_back(ors*(3.*_sf-_sd));
  intemp.push_back(_exi0);outtemp.push_back(_xi0);mestemp.push_back(221);
  factor.push_back(ors*(3.*_sf+_sd));
  int inspin,outspin;
  tPDVector extpart(3);
  for(unsigned int ix=0;ix<intemp.size();++ix) {
    if(intemp[ix]!=0&&outtemp[ix]!=0&&mestemp[ix]!=0) {
      extpart[0]=getParticleData(intemp[ix]);
      extpart[1]=getParticleData(outtemp[ix]);
      extpart[2]=getParticleData(mestemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin()) {
	_incomingB.push_back(intemp[ix]);
	_outgoingB.push_back(outtemp[ix]);
	_outgoingM.push_back(mestemp[ix]);
	if(iopt==1) {
	  inspin  = extpart[0]->iSpin();
	  outspin = extpart[1]->iSpin();
	  if(inspin==2&&outspin==2)
	    _prefactor.push_back(ort*factor[ix]/_fpi);
	  else if(inspin==4&&outspin==2)
	    _prefactor.push_back(ort*factor[ix]/_fpi);
	  else
	    throw Exception()<< "Invalid combination of spins in "
					<< "SU3BaryonOctetOctetScalarDecayer::" 
					<< "setupModes()" 
					<< Exception::abortnow;
	}
      }
    }
  }
}
 
void SU3BaryonOctetOctetScalarDecayer::dataBaseOutput(ofstream & output,
						      bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Fcoupling " << _sf << "\n";
  output << "newdef " << name() << ":Dcoupling " << _sd << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":Proton " << _proton << "\n";
  output << "newdef " << name() << ":Neutron " << _neutron << "\n";
  output << "newdef " << name() << ":Sigma+ " << _sigmap << "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Sigma- " << _sigmam << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":Xi0 " << _xi0 << "\n";
  output << "newdef " << name() << ":Xi- " << _xim << "\n"; 
  output << "newdef " << name() << ":ExcitedProton " << _eproton << "\n";
  output << "newdef " << name() << ":ExcitedNeutron " << _eneutron << "\n";
  output << "newdef " << name() << ":ExcitedSigma+ " << _esigmap << "\n";
  output << "newdef " << name() << ":ExcitedSigma0 " << _esigma0 << "\n";
  output << "newdef " << name() << ":ExcitedSigma- " << _esigmam << "\n";
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  output << "newdef " << name() << ":ExcitedXi0 " << _exi0 << "\n";
  output << "newdef " << name() << ":ExcitedXi- " << _exim << "\n"; 
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SU3BaryonSingletOctetPhotonDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonSingletOctetPhotonDecayer class.
//

#include "SU3BaryonSingletOctetPhotonDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonSingletOctetPhotonDecayer::SU3BaryonSingletOctetPhotonDecayer() {
  // the coupling
  _c=0.252/GeV;
  // the relative parities of the two baryon multiplets
  _parity=false;
  // PDG codes for the various ground state baryons
  _sigma0   = 3212;
  _lambda   = 3122;
  // PDG codes for the excited baryon
  _elambda  = 3124;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonSingletOctetPhotonDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<_outgoingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_elambda);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(ParticleID::gamma)};
    double wgtmax = _maxweight.size()>numberModes() ?
      _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode=new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

void SU3BaryonSingletOctetPhotonDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonSingletOctetPhotonDecayer::modeNumber(bool & cc,tcPDPtr parent,
					  const tPDVector & children) const {
  int imode(-1);
  if(_outgoingB.size()==0){setupModes(0);}
  // must be two outgoing particles
  if(children.size()!=2){return imode;}
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id()),iout;
  if(id1==ParticleID::gamma){iout=id2;}
  else if(id2==ParticleID::gamma){iout=id1;}
  else{return imode;}
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_elambda){
      if(iout==_outgoingB[ix]) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_elambda) {
      if(iout==-_outgoingB[ix]) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_outgoingB.size()&&imode<0);
  return imode;
}

void SU3BaryonSingletOctetPhotonDecayer::persistentOutput(PersistentOStream & os) const {
  os << ounit(_c,1./GeV) << _parity << _sigma0 << _lambda << _elambda << _outgoingB 
     << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonSingletOctetPhotonDecayer::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_c,1./GeV) >> _parity >> _sigma0 >> _lambda >> _elambda >> _outgoingB 
     >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonSingletOctetPhotonDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonSingletOctetPhotonDecayer("Herwig::SU3BaryonSingletOctetPhotonDecayer", "HwBaryonDecay.so");

void SU3BaryonSingletOctetPhotonDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetPhotonDecayer> documentation
    ("The SU3BaryonSingletOctetPhotonDecayer class performs the decay"
     " of a singlet baryon to an octet baryon and a photon.");

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,InvEnergy> interfaceCoupling
    ("Coupling",
     "The C coupling of the baryon resonances.",
     &SU3BaryonSingletOctetPhotonDecayer::_c, 1./GeV, 0.252/GeV, -10./GeV, 10.0/GeV,
     false, false, true);

  static Switch<SU3BaryonSingletOctetPhotonDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetPhotonDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetPhotonDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetPhotonDecayer::_elambda, 3124, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetPhotonDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetPhotonDecayer::_maxweight,
     0, 0, 0, -10000, 10000, false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-1
void SU3BaryonSingletOctetPhotonDecayer::
halfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
		       Complex&A1,Complex&A2,Complex&B1,Complex&B2) const {
  if(_parity) {
    A1=    _prefactor[imode]*(m0+m1);
    B1=0.;
    A2=-2.*_prefactor[imode]*(m0+m1);
    B2=0.;
  }
  else {
    A1=0.;
    B1=    _prefactor[imode]*(m1-m0);
    A2=0.;
    B2=-2.*_prefactor[imode]*(m0+m1);
  }
}

// couplings for spin-1/2 to spin-3/2 spin-1
void SU3BaryonSingletOctetPhotonDecayer::
threeHalfHalfVectorCoupling(int imode,Energy m0,Energy m1,Energy,
			    Complex&A1,Complex&A2,Complex&A3,
			    Complex&B1,Complex&B2,Complex&B3) const {
  A3=0.;B3=0.;
  if(_parity) {
    A1=0.;
    B1=-_prefactor[imode]*(m0+m1);
    A2=0.;
    B2= _prefactor[imode]*(m0+m1);
  }
  else {
    A1=_prefactor[imode]*(m0-m1);B1=0.;
    A2=_prefactor[imode]*(m0+m1);B2=0.;
  }
}

// set up the decay modes
void SU3BaryonSingletOctetPhotonDecayer::setupModes(unsigned int iopt) const {
  if(_outgoingB.size()!=0&&iopt==0){return;}
  if(iopt==1){_outgoingB.clear();}
  // set up for the various different decay modes
  vector<int> outtemp;
  vector<InvEnergy> factor;
  if(_elambda==0)
    throw Exception() << "Invalid incoming particle in "
				 << "SU3BaryonSingletOctetScalarDecayer::" 
				 << "setupModes()" << Exception::abortnow;
  // decays of the excited lambda
  outtemp.push_back(_sigma0);factor.push_back(_c/sqrt(2.));
  outtemp.push_back(_lambda);factor.push_back(_c/sqrt(6.));
  tPDVector extpart(2);extpart[0]=getParticleData(_elambda);
  int inspin(extpart[0]->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix) {
    if(outtemp[ix]!=0) {
      extpart[1]=getParticleData(outtemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()) {
	_outgoingB.push_back(outtemp[ix]);
	if(iopt==1) {
	  outspin = extpart[1]->iSpin();
	  factor[ix] *=2.;
	  if(inspin==2&&outspin==2)      _prefactor.push_back(2.*factor[ix]);
	  else if(inspin==4&&outspin==2) _prefactor.push_back(   factor[ix]);
	  else
	    throw Exception() 
	      << "Invalid combination of spins in "
	      << "SU3BaryonSingletOctetScalarDecayer::" 
	      << "setupModes()" << Exception::abortnow;
	}
      }
    }
  }
}

void SU3BaryonSingletOctetPhotonDecayer::dataBaseOutput(ofstream & output,
							bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Coupling " << _c*GeV << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) {
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  }
  if(header) output << "\n\" where BINARY ThePEGName=\"" 
		    << fullName() << "\";" << endl;
}
#line 1 "./SU3BaryonSingletOctetScalarDecayer.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SU3BaryonSingletOctetScalarDecayer class.
//

#include "SU3BaryonSingletOctetScalarDecayer.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

SU3BaryonSingletOctetScalarDecayer::SU3BaryonSingletOctetScalarDecayer() {
  // the coupling
  _c=0.39;
  // the relative parities of the two baryon multiplets
  _parity=false;
  // the pion decay constant
  _fpi=130.7*MeV;
  // PDG codes for the various ground state baryons
  _proton   = 2212;
  _neutron  = 2112;
  _sigma0   = 3212;
  _sigmap   = 3222;
  _sigmam   = 3112;
  _lambda   = 3122;
  _xi0      = 3322;
  _xim      = 3312;
  // PDG codes for the excited baryon
  _elambda  = 13122;
  // intermediates
  generateIntermediates(false);
}

void SU3BaryonSingletOctetScalarDecayer::doinit() {
  Baryon1MesonDecayerBase::doinit();
  // set up the decay modes
  setupModes(1);
  // set up the phase space and the couplings
  for(unsigned int ix=0;ix<_outgoingB.size();++ix) {
    tPDPtr    in  =  getParticleData(_elambda);
    tPDVector out = {getParticleData(_outgoingB[ix]),
    		     getParticleData(_outgoingM[ix])};
    double wgtmax = _maxweight.size()>numberModes() 
      ? _maxweight[numberModes()] : 1.;
    PhaseSpaceModePtr mode = 
      new_ptr(PhaseSpaceMode(in,out,wgtmax));
    addMode(mode);
  }
}

void SU3BaryonSingletOctetScalarDecayer::doinitrun() {
  Baryon1MesonDecayerBase::doinitrun();
  if(initialize()) {
    _maxweight.clear();
    for(unsigned int ix=0;ix<numberModes();++ix)
      _maxweight.push_back(mode(ix)->maxWeight());
  }
}

int SU3BaryonSingletOctetScalarDecayer::modeNumber(bool & cc,tcPDPtr parent,
						   const tPDVector & children) const {
  int imode(-1);
  if(_outgoingB.size()==0) setupModes(0);
  // must be two outgoing particles
  if(children.size()!=2) return imode;
  // ids of the particles
  int id0(parent->id());
  int id1(children[0]->id());
  int id2(children[1]->id());
  unsigned int ix(0);
  cc =false;
  do {
    if(id0==_elambda) {
      if((id1==_outgoingB[ix]&&id2==_outgoingM[ix])||
	 (id2==_outgoingB[ix]&&id1==_outgoingM[ix])) {
	imode=ix;
	cc=false;
      }
    }
    else if(id0==-_elambda) {
      if((id1==-_outgoingB[ix]&&id2==-_outgoingM[ix])||
	 (id2==-_outgoingB[ix]&&id1==-_outgoingM[ix])) {
	imode=ix;
	cc=true;
      }
      if(((id1==-_outgoingB[ix]&&id2==_outgoingM[ix])||
	  (id2==-_outgoingB[ix]&&id1==_outgoingM[ix]))&&
	 (_outgoingM[ix]==111||_outgoingM[ix]==221||_outgoingM[ix]==331||
	  _outgoingM[ix]==223||_outgoingM[ix]==333)) {
	imode=ix;
	cc=true;
      }
    }
    ++ix;
  }
  while(ix<_outgoingB.size()&&imode<0);
  return imode;
}

void SU3BaryonSingletOctetScalarDecayer::persistentOutput(PersistentOStream & os) const {
  os << _c << _parity << ounit(_fpi,GeV) << _proton << _neutron << _sigma0 << _sigmap 
     << _sigmam << _lambda << _xi0 << _xim << _elambda << _outgoingB 
     << _outgoingM << _maxweight << ounit(_prefactor,1./GeV);
}

void SU3BaryonSingletOctetScalarDecayer::persistentInput(PersistentIStream & is, int) {
  is >> _c >> _parity >> iunit(_fpi,GeV) >> _proton >> _neutron >> _sigma0 >> _sigmap 
     >> _sigmam >> _lambda >> _xi0 >> _xim >> _elambda >> _outgoingB 
     >> _outgoingM >> _maxweight >> iunit(_prefactor,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<SU3BaryonSingletOctetScalarDecayer,Baryon1MesonDecayerBase>
describeHerwigSU3BaryonSingletOctetScalarDecayer("Herwig::SU3BaryonSingletOctetScalarDecayer", "HwBaryonDecay.so");

void SU3BaryonSingletOctetScalarDecayer::Init() {

  static ClassDocumentation<SU3BaryonSingletOctetScalarDecayer> documentation
    ("The SU3BaryonSingletOctetScalarDecayer class is designed for"
     "the decay of an excited SU(3) singlet baryon");

  static Parameter<SU3BaryonSingletOctetScalarDecayer,double> interfaceCcoupling
    ("Coupling",
     "The C coupling of the baryon resonances",
     &SU3BaryonSingletOctetScalarDecayer::_c, 0.39, -10.0, 10.0,
     false, false, true);

  static Switch<SU3BaryonSingletOctetScalarDecayer,bool> interfaceParity
    ("Parity",
     "The relative parities of the two multiplets.",
     &SU3BaryonSingletOctetScalarDecayer::_parity, true, false, false);
  static SwitchOption interfaceParitySame
    (interfaceParity,
     "Same",
     "Same parity",
     true);
  static SwitchOption interfaceParityDifferent
    (interfaceParity,
     "Different",
     "Opposite parity",
     false);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,Energy> interfaceFpi
    ("Fpi",
     "The pion decay constant.",
     &SU3BaryonSingletOctetScalarDecayer::_fpi, MeV, 130.7*MeV, 100.0*MeV, 200.0*MeV,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceProton
    ("Proton",
     "The PDG code for the lighter proton-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_proton, 2212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceNeutron
    ("Neutron",
     "The PDG code for the lighter neutron-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_neutron, 2112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmap
    ("Sigma+",
     "The PDG code for the lighter Sigma+-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmap, 3222, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigma0
    ("Sigma0",
     "The PDG code for the lighter Sigma0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigma0, 3212, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceSigmam
    ("Sigma-",
     "The PDG code for the lighter Sigma--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_sigmam, 3112, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceLambda
    ("Lambda",
     "The PDG code for the lighter Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_lambda, 3122, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXi0
    ("Xi0",
     "The PDG code for the lighter Xi0-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xi0, 3322, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceXim
    ("Xi-",
     "The PDG code for the lighter Xi--like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_xim, 3312, 0, 1000000,
     false, false, true);

  static Parameter<SU3BaryonSingletOctetScalarDecayer,int> interfaceExcitedLambda
    ("ExcitedLambda",
     "The PDG code for the heavier Lambda-like baryon.",
     &SU3BaryonSingletOctetScalarDecayer::_elambda, 13122, 0, 1000000,
     false, false, true);

  static ParVector<SU3BaryonSingletOctetScalarDecayer,double> interfaceMaxWeight
    ("MaxWeight",
     "The maximum weight for the decay mode",
     &SU3BaryonSingletOctetScalarDecayer::_maxweight,
     0, 0, 0, 0., 100., false, false, true);

}

// couplings for spin-1/2 to spin-1/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::
halfHalfScalarCoupling(int imode,Energy m0,Energy m1,Energy,
		       Complex& A, Complex& B) const {
  if(_parity) {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
  else {
    A=_prefactor[imode]*(m0-m1);
    B=0.;
  }
}

// couplings for spin-1/2 to spin-3/2 spin-0
void SU3BaryonSingletOctetScalarDecayer::
threeHalfHalfScalarCoupling(int imode,Energy m0, Energy m1,Energy,
			    Complex& A, Complex& B) const {
  if(_parity) {
    A=_prefactor[imode]*(m0+m1);
    B=0.;
  }
  else {
    A=0.;
    B=_prefactor[imode]*(m0+m1);
  }
}

// set up the decay modes
void SU3BaryonSingletOctetScalarDecayer::setupModes(unsigned int iopt) const {
  if(_outgoingB.size()!=0&&iopt==0) return;
  if(iopt==1) {
    _outgoingB.clear();
    _outgoingM.clear();
  }
  // set up for the various different decay modes
  vector<int> outtemp,mestemp;
  double rt(sqrt(2.));
  if(_elambda==0)
    throw Exception() << "Invalid incoming particle in "
				 << "SU3BaryonSingletOctetScalarDecayer::" 
				 << "setupModes()" << Exception::abortnow;
  // decays of the excited lambda
  outtemp.push_back(_sigma0);mestemp.push_back(111);
  outtemp.push_back(_sigmap);mestemp.push_back(-211);
  outtemp.push_back(_sigmam);mestemp.push_back(211);
  outtemp.push_back(_lambda);mestemp.push_back(221);
  outtemp.push_back(_xim);mestemp.push_back(321);
  outtemp.push_back(_xi0);mestemp.push_back(311);
  outtemp.push_back(_proton);mestemp.push_back(-321);
  outtemp.push_back(_neutron);mestemp.push_back(-311);
  tPDVector extpart(3);
  extpart[0]=getParticleData(_elambda);
  int inspin(extpart[0]->iSpin()),outspin;
  for(unsigned int ix=0;ix<outtemp.size();++ix) {
    if(outtemp[ix]!=0&&mestemp[ix]!=0) {
      extpart[1]=getParticleData(outtemp[ix]);
      extpart[2]=getParticleData(mestemp[ix]);
      if(extpart[0]->massMax()>extpart[1]->massMin()+extpart[2]->massMin()) {
	_outgoingB.push_back(outtemp[ix]);
	_outgoingM.push_back(mestemp[ix]);
	if(iopt==1) {
	  outspin = extpart[1]->iSpin();
	  if(inspin==2&&outspin==2)      _prefactor.push_back(_c*rt/_fpi);
	  else if(inspin==4&&outspin==2) _prefactor.push_back(_c*rt/_fpi);
	  else throw Exception() 
	    << "Invalid combination of spins in "
	    << "SU3BaryonSingletOctetScalarDecayer::" 
	    << "setupModes()" << Exception::abortnow;
	}
      }
    }
  }
}

void SU3BaryonSingletOctetScalarDecayer::dataBaseOutput(ofstream & output,
							bool header) const {
  if(header) output << "update decayers set parameters=\"";
  Baryon1MesonDecayerBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":Coupling " << _c << "\n";
  output << "newdef " << name() << ":Parity " << _parity<< "\n";
  output << "newdef " << name() << ":Fpi " << _fpi/MeV << "\n";
  output << "newdef " << name() << ":Proton " << _proton << "\n";
  output << "newdef " << name() << ":Neutron " << _neutron << "\n";
  output << "newdef " << name() << ":Sigma+ " << _sigmap << "\n";
  output << "newdef " << name() << ":Sigma0 " << _sigma0 << "\n";
  output << "newdef " << name() << ":Sigma- " << _sigmam << "\n";
  output << "newdef " << name() << ":Lambda " << _lambda << "\n";
  output << "newdef " << name() << ":Xi0 " << _xi0 << "\n";
  output << "newdef " << name() << ":Xi- " << _xim << "\n"; 
  output << "newdef " << name() << ":ExcitedLambda " << _elambda << "\n";
  for(unsigned int ix=0;ix<_maxweight.size();++ix) 
    output << "insert " << name() << ":MaxWeight " << ix << " " 
	   << _maxweight[ix] << "\n";
  if(header) 
    output << "\n\" where BINARY ThePEGName=\"" 
	   << fullName() << "\";" << endl;
}
