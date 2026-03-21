#line 1 "./DalitzResonance.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzResonance class.
//

#include "DalitzResonance.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "FlatteResonance.h"
#include "MIPWA.h"
#include "PiPiI2.h"
#include "DalitzKMatrix.h"
#include "DalitzLASS.h"
#include "DalitzGS.h"
#include "DalitzSigma.h"

using namespace Herwig;

void DalitzResonance::persistentOutput(PersistentOStream & os) const {
  os << id << oenum(type) << ounit(mass,GeV) << ounit(width,GeV)
     << daughter1 << daughter2 << spectator
     << amp << ounit(R,1./GeV);
}

void DalitzResonance::persistentInput(PersistentIStream & is, int) {
  is >> id >> ienum(type) >> iunit(mass,GeV) >> iunit(width,GeV)
     >> daughter1 >> daughter2 >> spectator
     >> amp >> iunit(R,1./GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzResonance,Base>
  describeHerwigDalitzResonance("Herwig::DalitzResonance", "HwDalitzDecay.so");

void DalitzResonance::Init() {

  static ClassDocumentation<DalitzResonance> documentation
    ("The DalitzResonance class provides a container class for"
     " information on resonances in multi-body dalitz decays.");

}

Complex DalitzResonance::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  static const Complex ii = Complex(0.,1.);
  // non-resonant pieces
  if(abs(type)/10==10) return 1.;
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(mAB) -sqr(mA)-sqr(mB)) - sqr(mA*mB))/mAB;
  if(type==ResonanceType::BABARf0) {
    double rho = 2.*pAB/mAB;
    return GeV2/(sqr(mass)-sqr(mAB)-ii*mass*width*rho);
  }
  else if (type==ResonanceType::Spin0Complex) {
    complex<Energy> sR(mass,width);
    return GeV2/(sqr(sR)-sqr(mAB));
  }
  else if (type==ResonanceType::Flattef0  ||
	   type==ResonanceType::Flattea0 ||
	   type==ResonanceType::FlatteKstar0) {
    assert(false);
  }
  //  on-shell
  Energy  pR=sqrt(0.25*sqr( mass*mass - sqr(mA) - sqr(mB)) - sqr(mA*mB))/mass;
  // Blatt-Weisskopf factors
  double fR=1;
  unsigned int power(1);
  if(type!=ResonanceType::Spin0 &&
     type!=ResonanceType::Spin0E691) {
    double r1A(R*pR),r1B(R*pAB);
    // Blatt-Weisskopf factors and spin piece
    switch (type) {
    case ResonanceType::Spin0Gauss:
      fR = exp(-(r1B-r1A)/12.);
      // cerr << "testing scalar B " <<   exp(+r1A/12.) << "\n";
      break;
    case ResonanceType::Spin1: case ResonanceType::Spin1E691 :
      fR=sqrt( (1. + sqr(r1A)) / (1. + sqr(r1B)) );
      // cerr << "testing vector B " << sqrt(1. + sqr(r1A))   << "\n";
      power=3;
      break;
    case ResonanceType::Spin2: case ResonanceType::Spin2E691:
      fR = sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))) / (9. + sqr(r1B)*(3.+sqr(r1B))));
      // cerr << "testing tensor B " <<  sqrt( (9. + sqr(r1A)*(3.+sqr(r1A))))  << "\n";
      power=5;
      break;
    default :
      assert(false);
    }
  }
  // multiply by Breit-Wigner piece and return
  if (type/10 == 1 ) {
    return fR*sqrt(0.5*width/GeV/Constants::pi)*GeV/(mAB-mass-complex<Energy>(ZERO,0.5*width));
  }
  else {
    Energy gam = width*pow(pAB/pR,power)*(mass/mAB)*fR*fR;
    return fR*GeV2/(sqr(mass)-sqr(mAB)-mass*gam*ii);
  }
}

void DalitzResonance::dataBaseOutput(ofstream & output) {
  output << id << " " << oenum(type) << " "
	 << mass/GeV << " " << width/GeV << " "
	 << daughter1 << " " << daughter2 << " "
	 << spectator << " " 
	 << abs(amp) << " " << arg(amp) << " "
	 << R*GeV; 
}

DalitzResonancePtr DalitzResonance::readResonance(string arg, string & error) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long id = stoi(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  ResonanceType::Type type = static_cast<ResonanceType::Type>(stoi(stype));
  // mass and width
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy mass = stof(stype)*GeV;
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  Energy width = stof(stype)*GeV;
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int d1 = stoi(stype);
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int d2 = stoi(stype);
  // children
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  int sp = stoi(stype);
  if (sp==d1 || sp ==d2 || d1 == d2) {
    error =  "Daughters and spectator must all be different not " + std::to_string(d1) + ", " + std::to_string(d2) + ", " + std::to_string(sp);
    return DalitzResonancePtr();
  }
  // magnitude and phase
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double mag = stof(stype);
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  double phi = stof(stype);
  // radius
  stype = StringUtils::car(arg);
  arg   = StringUtils::cdr(arg);
  InvEnergy r = stof(stype)/GeV;
  // special for flate
  if( type==ResonanceType::Flattef0 ||
      type==ResonanceType::Flattea0 ||
      type==ResonanceType::FlatteKstar0) {
    // Flatte parameters
    // magnitude and phase
    vector<Energy> f;
    while (arg!="") {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      f.push_back( stof(stype)*GeV );
    }
    // add to list
    return new_ptr(FlatteResonance(id,type,mass,width,d1,d2,sp,mag,phi,r,f));
  }
  // MIPWA
  else if(type==ResonanceType::Spin0MIPWA) {
    // no of entries in table
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    int nn = stoi(stype);
    vector<Energy> en; en.reserve(nn);
    vector<double> mag2,phase2; mag2.reserve(nn); phase2.reserve(nn);
    for(int ix=0;ix<nn;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      Energy ee = stof(stype)*GeV;
      en.push_back(ee);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double mm=stof(stype);
      mag2.push_back(mm);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double pp=stof(stype);
      phase2.push_back(pp);
    }
    return new_ptr(MIPWA(id,type,mass,width,d1,d2,sp,mag,phi,r,en,mag2,phase2));
  }
  // I=2 pipi
  else if(type==ResonanceType::PiPiI2) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy a = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy2 b = stof(stype)/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy4 c = stof(stype)/GeV2/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy6 d = stof(stype)/GeV2/GeV2/GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mmin = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mmax = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double dEta=stof(stype);
    return new_ptr(PiPiI2(id,type,mass,width,d1,d2,sp,mag,phi,r,
			  a,b,c,d,mmin,mmax,dEta));
  }
  // K-matrix
  else if(type==ResonanceType::KMatrix) {
    // no of poles
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int npole= stoi(stype);
    // no of channels
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int nchannels= stoi(stype);
    // matrix location
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int imat = stoi(stype);
    // this channel
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int chan = stoi(stype);
    // expansion point for the constants terms
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy2 sc = GeV2*stof(stype);
    // type of expansion
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int itype= stoi(stype);
    vector<pair<double,double> > beta;
    // first loop over the coefficients of the poles
    for(unsigned int ix=0;ix<npole;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      double b = stof(stype);
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      beta.push_back(make_pair(b,stof(stype)));
    }
    // now over the power series for the different channels
    vector<pair<double,vector<double > > > coeffs(nchannels);
    for(unsigned int ix=0;ix<nchannels;++ix) {
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      unsigned int nterms = stoi(stype);
      for(unsigned int iy=0;iy<nterms;++iy) {
    	stype = StringUtils::car(arg);
    	arg   = StringUtils::cdr(arg);
    	coeffs[ix].second.push_back(stof(stype));
      }
      stype = StringUtils::car(arg);
      arg   = StringUtils::cdr(arg);
      coeffs[ix].first = stof(stype);
    }
    // finally make the channel
    return new_ptr(DalitzKMatrix(id,type,mass,width,d1,d2,sp,mag,phi,r,imat,chan,sc,itype,beta,coeffs));
  }
  // LASS
  else if(type==ResonanceType::LASS) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    unsigned int iopt = stoi(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double FNR = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double phiNR = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double FRes = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    double phiRes = stof(stype);
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy ascat = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy reff = stof(stype)/GeV;
    // finally make the channel
    return new_ptr(DalitzLASS(id,type,mass,width,d1,d2,sp,mag,phi,r,iopt,
			      FNR,phiNR,FRes,phiRes,ascat,reff));
  }
  // Bugg sigma form
  else if(type==ResonanceType::Sigma) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy2 a = stof(stype)*GeV2;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy b1 = stof(stype)*GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    InvEnergy b2 = stof(stype)/GeV;
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy g4pi = stof(stype)*GeV;
    return new_ptr(DalitzSigma(id,type,mass,width,d1,d2,sp,mag,phi,r,a,b1,b2,g4pi));
  }
  // GS form
  else if(type==ResonanceType::Spin1GS) {
    stype = StringUtils::car(arg);
    arg   = StringUtils::cdr(arg);
    Energy mpi = stof(stype)*GeV;
    return new_ptr(DalitzGS(id,type,mass,width,d1,d2,sp,mag,phi,r,mpi));
  }
  // otherwise add to list
  else {
    return new_ptr(DalitzResonance(id,type,mass,width,d1,d2,sp,mag,phi,r));
  }
}
#line 1 "./FlatteResonance.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlatteResonance class.
//

#include "FlatteResonance.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"

using namespace Herwig;
void FlatteResonance::persistentOutput(PersistentOStream & os) const {
  os << ounit(g_,GeV);
}

void FlatteResonance::persistentInput(PersistentIStream & is, int) {
  is >> iunit(g_,GeV);
}

//The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<FlatteResonance,DalitzResonance>
describeHerwigFlatteResonance("Herwig::FlatteResonance", "HwDalitzDecay.so");

void FlatteResonance::Init() {

  static ClassDocumentation<FlatteResonance> documentation
    ("The FlatteResonance class implements the Flatte lineshape for Dalitz decays.");

}

void FlatteResonance::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  for(const Energy & g : g_)
    output << " " << g/GeV; 
}

namespace {

double rho2(const Energy2 & q2, const Energy & m1, const Energy & m2) {
  return (1.-sqr(m1+m2)/q2)*(1.-sqr(m1-m2)/q2);
}

}

Complex FlatteResonance::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static const Complex ii = Complex(0.,1.);
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  Energy mK  = CurrentGenerator::current().getParticleData(321)->mass();
  Energy2 q2=sqr(mAB);
  if(type==ResonanceType::Flattef0) {
    assert(g_.size()==2);
    complex<Energy2> MGamma = sqr(g_[0])*sqrt(max(0.,rho2(q2,mpi,mpi)));
    double arg = rho2(q2,mK,mK);
    MGamma += mAB>2.*mK ? sqr(g_[1])*sqrt(arg) : sqr(g_[1])*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else if(type==ResonanceType::Flattea0) {
    assert(g_.size()==2 || g_.size()==3);
    Energy meta = CurrentGenerator::current().getParticleData(221)->mass();
    complex<Energy2> MGamma = sqr(g_[0])*sqrt(max(0.,rho2(q2,meta,mpi)));
    double arg = rho2(q2,mK,mK);
    MGamma += mAB>2.*mK ? sqr(g_[1])*sqrt(arg) : sqr(g_[1])*ii*sqrt(abs(arg));
    if(g_.size()==3 ) {
      Energy metap = CurrentGenerator::current().getParticleData(331)->mass();
      arg = rho2(q2,metap,mpi);
      MGamma += mAB>mpi+metap ? sqr(g_[2])*sqrt(arg) : sqr(g_[2])*ii*sqrt(abs(arg));
    }
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else if(type==ResonanceType::FlatteKstar0) {
    assert(g_.size()==2);
    Energy metaP = CurrentGenerator::current().getParticleData(331)->mass();
    complex<Energy2> MGamma = sqr(g_[0])*sqrt(max(0.,rho2(q2,mK,mpi)));
    double arg = rho2(q2,mK,metaP);
    MGamma += mAB>mK+metaP ? sqr(g_[1])*sqrt(arg) : sqr(g_[1])*ii*sqrt(abs(arg));
    return GeV2/(sqr(mass)-sqr(mAB)-ii*MGamma);
  }
  else {
    assert(false);
  }
}
#line 1 "./MIPWA.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MIPWA class.
//

#include "MIPWA.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void MIPWA::persistentOutput(PersistentOStream & os) const {
  os << ounit(energy_,GeV) << mag_ << phase_;
}

void MIPWA::persistentInput(PersistentIStream & is, int) {
  is >> iunit(energy_,GeV) >> mag_ >> phase_;
}

// The following static variable is needed for the type
// description system in ThePEG. 
DescribeClass<MIPWA,DalitzResonance>
  describeHerwigMIPWA("Herwig::MIPWA", "HwDalitzDecay.so");

void MIPWA::Init() {

  static ClassDocumentation<MIPWA> documentation
    ("The MIPWA class allows the use on experimental extractions from"
     " Model Independent Partial Wave Analyses. ");

}

Complex MIPWA::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static Complex ii(0.,1.);
  if(!iMag_) {
    iMag_   = make_InterpolatorPtr(mag_  ,energy_,3);
    iPhase_ = make_InterpolatorPtr(phase_,energy_,3);
  }
  return (*iMag_)(mAB)*exp(ii*(*iPhase_)(mAB));
}

void MIPWA::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << energy_.size();
  for(unsigned int ix=0;ix<energy_.size();++ix)
    output << " " << energy_[ix]/GeV << " " << mag_[ix] << " " << phase_[ix];
}
#line 1 "./PiPiI2.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PiPiI2 class.
//

#include "PiPiI2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void PiPiI2::persistentOutput(PersistentOStream & os) const {
  os << ounit(a_,1./GeV) << ounit(b_,1./GeV2) << ounit(c_,1./GeV2/GeV2) << ounit(d_,1./GeV2/GeV2/GeV2)
     << ounit(mmin_,GeV) << ounit(mmax_,GeV) << deltaEta_;
}

void PiPiI2::persistentInput(PersistentIStream & is, int) {
  is >> iunit(a_,1./GeV) >> iunit(b_,1./GeV2) >> iunit(c_,1./GeV2/GeV2) >> iunit(d_,1./GeV2/GeV2/GeV2)
     >> iunit(mmin_,GeV) >> iunit(mmax_,GeV) >> deltaEta_;
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<PiPiI2,DalitzResonance>
describeHerwigPiPiI2("Herwig::PiPiI2", "HwDalitzDecay.so");

void PiPiI2::Init() {

  static ClassDocumentation<PiPiI2> documentation
    ("The PiPiI2 class provides an implementation of the I=2 s-eave for pipi.");

}

Complex PiPiI2::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static Complex ii(0.,1.);
  Energy2 m2 = sqr(mAB);
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  double delta = -a_*sqrt(0.25*m2-sqr(mpi))/(1.+m2*(b_+m2*(c_+d_*m2)));
  double eta = 1.;
  if(mAB>mmax_) {
    eta = 1. - deltaEta_;
  }
  else if(mAB>mmin_) {
    eta = 1. - 0.5*deltaEta_*(1.-cos(Constants::pi*(mAB-mmin_)/(mmax_-mmin_)));
  }
  return -0.5*ii*(eta*exp(2.*ii*delta)-1.);
}

void PiPiI2::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << a_*GeV << " " << b_*GeV2 << " " << c_*GeV2*GeV2 << " " << d_*GeV2*GeV2*GeV2
	 << " " << mmin_/GeV << " " << mmax_/GeV << " " << deltaEta_;
}
#line 1 "./DalitzKMatrix.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzKMatrix class.
//

#include "DalitzKMatrix.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DalitzKMatrix::persistentOutput(PersistentOStream & os) const {
  os << kMatrix_ << channel_ << imat_ << ounit(sc_,GeV2) << expType_ << beta_ << coeffs_;
}

void DalitzKMatrix::persistentInput(PersistentIStream & is, int) {
  is >> kMatrix_ >> channel_ >> imat_ >> iunit(sc_,GeV2) >> expType_ >> beta_ >> coeffs_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzKMatrix,DalitzResonance>
describeHerwigDalitzKMatrix("Herwig::DalitzKMatrix", "HwDalitzDecay.so");

void DalitzKMatrix::Init() {

  static ClassDocumentation<DalitzKMatrix> documentation
    ("The DalitzKMatrix class allows the use of \f$K\f$-matrices in Dalitz decays");

}

Complex DalitzKMatrix::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  Energy2 s = sqr(mAB);
  double sHat = (s-sc_)/GeV2;
  // construct the p-vector
  ublas::vector<Complex> pVector(coeffs_.size());
  // compute the terms
  for(unsigned int ix=0;ix<coeffs_.size();++ix) {
    Complex val(0.);
    // first the pole piece
    for(unsigned int iy=0;iy<kMatrix_->poles().size();++iy) {
      Complex piece = GeV*beta_[iy]*kMatrix_->poleCouplings()[iy][ix]/kMatrix_->poles()[iy];
      for(unsigned int iz=0;iz<kMatrix_->poles().size();++iz) {
	if(iz==iy) continue;
	piece *= 1. - s/kMatrix_->poles()[iz];
      }
      val +=piece;
    }
    Complex fact=exp(Complex(0.,coeffs_[ix].first));
    for(unsigned int iz=0;iz<kMatrix_->poles().size();++iz)
      fact *= 1. - s/kMatrix_->poles()[iz];
    // then the polynomial piece
    double poly=coeffs_[ix].second[0];
    if(expType_==0) {
      for(unsigned int iz=1;iz<coeffs_[ix].second.size();++iz)
	poly += coeffs_[ix].second[iz]*pow(sHat,iz);
    }
    else {
      poly *= (GeV2-sc_)/(s-sc_);
    }
    // store the answer
    pVector[ix]=val+fact*poly;
  }
  ublas::vector<Complex> amps = kMatrix_->amplitudes(s,pVector,true);
  return amps[channel_];
}

void DalitzKMatrix::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << kMatrix_->poles().size()
	 << " " << kMatrix_->numberOfChannels() << " "
	 << imat_ << " " << channel_ << " " << sc_/GeV2 << " " << expType_;
  for(unsigned int ix=0;ix<beta_.size();++ix)
    output << " " << abs(beta_[ix]) << " " << arg(beta_[ix]);
  for(unsigned int ix=0;ix<coeffs_.size();++ix) {
    output << " " << coeffs_[ix].second.size();
    for(unsigned int iy=0;iy<coeffs_[ix].second.size();++iy)
      output << " " << coeffs_[ix].second[iy];
    output << " " << coeffs_[ix].first;
  }
}
#line 1 "./DalitzLASS.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzLASS class.
//

#include "DalitzLASS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DalitzLASS::persistentOutput(PersistentOStream & os) const {
  os << opt_ << FNR_ << phiNR_ << FRes_ << phiRes_
     << ounit(aScat_,1./GeV) << ounit(rEff_,1./GeV);
}
void DalitzLASS::persistentInput(PersistentIStream & is, int) {
  is >> opt_ >> FNR_ >> phiNR_ >> FRes_ >> phiRes_
     >> iunit(aScat_,1./GeV) >> iunit(rEff_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzLASS,DalitzResonance>
describeHerwigDalitzLASS("Herwig::DalitzLASS", "HwDalitzDecay.so");

void DalitzLASS::Init() {

  static ClassDocumentation<DalitzLASS> documentation
    ("The DalitzLASS class implements the LASS parameterization of the Kpi s-wave.");

}

Complex DalitzLASS::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  static const Complex ii = Complex(0.,1.);
  // momenta for the resonance decay
  // off-shell
  Energy pAB=sqrt(0.25*sqr(sqr(mAB) -sqr(mA)-sqr(mB)) - sqr(mA*mB))/mAB;
  // on-shell
  Energy  pR=sqrt(0.25*sqr( mass*mass - sqr(mA) - sqr(mB)) - sqr(mA*mB))/mass;
  // non-resonant phase
  double NRphase = phiNR_+atan(1./(1./(aScat_*pAB)+0.5*rEff_*pAB));
  // resonant phase
  Energy Gamma  = width*(pAB/pR)*mass/mAB;
  double Rphase = atan(mass*Gamma/(sqr(mass)-sqr(mAB)));
  // return the result
  // BABar/BES and hopefully right form
  if (opt_==0) {
    return double(mAB/pAB)*(FNR_*sin(NRphase)*exp(ii*NRphase) +FRes_*sin(Rphase)*exp(ii*(Rphase+phiRes_+2.*NRphase)));
  }
  // BELLE form
  else if(opt_==1) {
    return (FNR_*sin(NRphase)*exp(ii*NRphase) +FRes_*sin(Rphase)*exp(ii*(Rphase+phiRes_+2.*NRphase)));
  }
  // original LASS form
  else if(opt_==2) {
    double delta = Rphase+NRphase;
    return double(mAB/pAB)*sin(delta)*exp(ii*delta);
  }
  else
    assert(false);
}

void DalitzLASS::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << " " << opt_ << " " << FNR_ << " " << phiNR_ << " " << FRes_ << " " << phiRes_
	 << " " << aScat_*GeV << " " << rEff_*GeV;
}
#line 1 "./DalitzGS.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzGS class.
//

#include "DalitzGS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Herwig/Decay/ResonanceHelpers.h"

using namespace Herwig;

DalitzGS::DalitzGS(long pid, ResonanceType::Type rtype, Energy m, Energy w,
		   unsigned int d1, unsigned int d2, unsigned int s,
		   double mag, double phi, InvEnergy rr, Energy mpi)
  : DalitzResonance(pid,rtype,m,w,d1,d2,s,mag,phi,rr), mpi_(mpi) {
  hres_ = Resonance::Hhat(sqr(mass),mass,width,mpi_,mpi_);
  dh_   = Resonance::dHhatds(mass,width,mpi_,mpi_);
  h0_   = Resonance::H(ZERO,mass,width,mpi_,mpi_,dh_,hres_);
}

void DalitzGS::persistentOutput(PersistentOStream & os) const {
  os << ounit(mpi_,GeV) << dh_ << ounit(hres_,GeV2) << ounit(h0_,GeV2);
}

void DalitzGS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mpi_,GeV) >> dh_ >> iunit(hres_,GeV2) >> iunit(h0_,GeV2);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzGS,DalitzResonance>
describeHerwigDalitzGS("Herwig::DalitzGS", "HwDalitzDalitz.so");

void DalitzGS::Init() {

  static ClassDocumentation<DalitzGS> documentation
    ("The DalitzGS class implements the Gounaris and Sakurai Phys. Rev. "
     "Lett. 21, 244 (1968) form for the propagator.");

}

void DalitzGS::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << mpi_/GeV;
}

Complex DalitzGS::BreitWigner(const Energy & mAB, const Energy & mA, const Energy & mB) const {
  Energy pAB=sqrt(0.25*sqr(sqr(mAB) -sqr(mA)-sqr(mB)) - sqr(mA*mB))/mAB;
  Energy  pR=sqrt(0.25*sqr( mass*mass - sqr(mA) - sqr(mB)) - sqr(mA*mB))/mass;
  double r1A(R*pR),r1B(R*pAB);
  double fR=sqrt( (1. + sqr(r1A)) / (1. + sqr(r1B)) );
  return fR*GeV2/sqr(mass)*Resonance::BreitWignerGS(sqr(mAB),mass,width,mpi_,mpi_,h0_,dh_,hres_);
}
#line 1 "./DalitzSigma.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzSigma class.
//

#include "DalitzSigma.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

void DalitzSigma::persistentOutput(PersistentOStream & os) const {
  os << ounit(a_,GeV2) << ounit(b1_,GeV) << ounit(b2_,1./GeV) << ounit(g4Pi_,GeV);
}

void DalitzSigma::persistentInput(PersistentIStream & is, int) {
  is >> iunit(a_,GeV2) >> iunit(b1_,GeV) >> iunit(b2_,1./GeV) >> iunit(g4Pi_,GeV);
}


// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<DalitzSigma,DalitzResonance>
  describeHerwigDalitzSigma("Herwig::DalitzSigma", "HwDalitzDecay.so");

void DalitzSigma::Init() {

  static ClassDocumentation<DalitzSigma> documentation
    ("The DalitzSigma class implements the model of Bou and Zou for the sigma propagator");

}

void DalitzSigma::dataBaseOutput(ofstream & output) {
  DalitzResonance::dataBaseOutput(output);
  output << " " << a_/GeV2 << " " << b1_/GeV << " " << b2_*GeV << " " << g4Pi_/GeV;
}

Complex DalitzSigma::BreitWigner(const Energy & mAB, const Energy & , const Energy & ) const {
  static const Complex II(0.,1.);
  Energy2 s(sqr(mAB));
  Energy mpi = CurrentGenerator::current().getParticleData(111)->mass();
  Energy2 sA=0.5*sqr(mpi);
  // two pion width
  Energy gamma = (b1_+b2_*s)*exp(-(s-sqr(mass))/a_)*(s-sA)/(sqr(mass)-sA)*
    sqrt(1.-4.*sqr(mpi)/s)/sqrt(1.-4.*sqr(mpi/mass));
  if(mAB>4.*mpi) {
    gamma+= g4Pi_*rho4pi(s,mpi)/rho4pi(sqr(mass),mpi);
  }
  // return the propagator
  return GeV2/(sqr(mass)-s-II*mass*gamma);
}
#line 1 "./DalitzBase.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DalitzBase class.
//

#include "DalitzBase.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "Herwig/Decay/PhaseSpaceMode.h"

using namespace Herwig;

void DalitzBase::persistentOutput(PersistentOStream & os) const {
  os << ounit(rParent_,1./GeV) << resonances_ << maxWgt_ << weights_ 
     << channel1_ << channel2_ << incoming_ << outgoing_ << useAllK0_
     << kMatrix_;
}

void DalitzBase::persistentInput(PersistentIStream & is, int) {
  is >> iunit(rParent_,1./GeV) >> resonances_ >> maxWgt_ >> weights_ 
     >> channel1_ >> channel2_ >> incoming_ >> outgoing_ >> useAllK0_
     >> kMatrix_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeAbstractClass<DalitzBase,DecayIntegrator>
describeHerwigDalitzBase("Herwig::DalitzBase", "HwDalitzDecay.so");

void DalitzBase::Init() {

  static ClassDocumentation<DalitzBase> documentation
    ("The DalitzBase class provides a base class for the implementation of three-body Dalitz decays.");

  static Command<DalitzBase> interfaceSetExternal
    ("SetExternal",
     "Set the external particles for the decay mode",
     &DalitzBase::setExternal, false);
  
  static Command<DalitzBase> interfaceAddChannel
    ("AddChannel",
     "Add a channel for the description of the matrix element",
     &DalitzBase::addChannel, false);

  static Parameter<DalitzBase,InvEnergy> interfaceParentRadius
    ("ParentRadius",
     "The radius parameter for the Blatt-Weisskopf form-factor for the D",
     &DalitzBase::rParent_, 1./GeV, 5./GeV, ZERO, 10./GeV,
     false, false, Interface::limited);

  static Parameter<DalitzBase,double> interfaceMaximumWeight
    ("MaximumWeight",
     "The maximum weight for the phase-space sampling",
     &DalitzBase::maxWgt_, 1.0, 0.0, 1e20,
     false, false, Interface::limited);

  static ParVector<DalitzBase,double> interfaceWeights
    ("Weights",
     "The weights for the different channels for the phase-space integration",
     &DalitzBase::weights_, -1, 1.0, 0.0, 1.0,
     false, false, Interface::limited);
  
  static Parameter<DalitzBase,int> interfaceChannel1
    ("Channel1",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &DalitzBase::channel1_, -1, -1, 100,
     false, false, Interface::limited);
  
  static Parameter<DalitzBase,int> interfaceChannel2
    ("Channel2",
     "The first allowed channel, for debugging/calculation of fit fractions only",
     &DalitzBase::channel2_, -1, -1, 100,
     false, false, Interface::limited);
  
  static Switch<DalitzBase,bool> interfaceUseAllK0
    ("UseAllK0",
     "Use all K0 mesons when matching the mode",
     &DalitzBase::useAllK0_, false, false, false);
  static SwitchOption interfaceUseAllK0No
    (interfaceUseAllK0,
     "No",
     "Just use the identified state",
     false);
  static SwitchOption interfaceUseAllK0Yes
    (interfaceUseAllK0,
     "Yes",
     "Use all the states",
     true);

  static RefVector<DalitzBase,KMatrix> interfaceKMatrices
    ("KMatrices",
     "Any K-matrices needed to simulate the decay",
     &DalitzBase::kMatrix_, -1, false, false, true, false, false);

}

void DalitzBase::doinit() {
  if(incoming_!=0) {
    tPDPtr in = getParticleData(incoming_);
    vector<tPDPtr> out = {getParticleData(outgoing_[0]),
			  getParticleData(outgoing_[1]),
			  getParticleData(outgoing_[2])};
    createMode(in,out);
  }
  DecayIntegrator::doinit();
}

void DalitzBase::doinitrun() {
  if(!kMatrix_.empty()) {
    for(unsigned int ix=0;ix<resonances().size();++ix) {
      Ptr<Herwig::DalitzKMatrix>::transient_pointer mat =
	dynamic_ptr_cast<Ptr<Herwig::DalitzKMatrix>::transient_pointer>(resonances()[ix]);
      if(mat) {
	mat->setKMatrix(kMatrix_[mat->imatrix()]);
	// Energy2 s=3.*GeV2;
	// Complex amp = resonances()[ix]->BreitWigner(sqrt(s),0.139*GeV,0.139*GeV);
	// Energy2 s=0.5*GeV2;
	// while(s<3.5*GeV2) {
	//   Complex amp = resonances()[ix]->BreitWigner(sqrt(s),0.139*GeV,0.139*GeV);
	//   cerr << s/GeV2 << " " << abs(amp) << "\n";
	//   s+=0.01*GeV2;
	// }
      }
    }
  }
  DecayIntegrator::doinitrun();
  weights_.resize(mode(0)->channels().size());
  maxWgt_ = mode(0)->maxWeight();
  for(unsigned int iz=0;iz<mode(0)->channels().size();++iz) {
    weights_[iz]=mode(0)->channels()[iz].weight();
  }
}

void DalitzBase::createMode(tPDPtr in, tPDVector out) {
  PhaseSpaceModePtr mode = new_ptr(PhaseSpaceMode(in,out,maxWgt_));
  if(weights_.size()!=resonances().size()) {
    weights_=vector<double>(resonances_.size(),1./double(resonances_.size()));
  }
  unsigned int ix=0;
  for(DalitzResonancePtr res : resonances_) {
    tPDPtr resonance = getParticleData(res->id);
    if(resonance) {
      mode->addChannel((PhaseSpaceChannel(mode),0,resonance,0,res->spectator+1,1,res->daughter1+1,1,res->daughter2+1));
      resetIntermediate(resonance,res->mass,abs(res->width));
      ++ix;
    }
  }
  addMode(mode);
}

void DalitzBase::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DalitzBase base class
  DecayIntegrator::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ParentRadius " << rParent_*GeV << "\n";
  output << "newdef " << name() << ":UseAllK0 " << useAllK0_ << "\n";
  output << "newdef " << name() << ":MaximumWeight " << maxWgt_ << "\n";
  for(unsigned int ix=0;ix<kMatrix_.size();++ix) {
    output << "insert " << name() << ":KMatrices " << ix << " " << kMatrix_[ix]->fullName() << "\n";
  }
  for(unsigned int ix=0;ix<weights_.size();++ix) {
    output << "insert " << name() << ":Weights "
	   << ix << " " << weights_[ix] << "\n";
  }
  for(unsigned int ix=0;ix<resonances_.size();++ix) {
    output << "do " << name() << ":AddChannel ";
    resonances_[ix]->dataBaseOutput(output);
    output << "\n";
  }
  output << "do " << name() << ":SetExternal " << incoming_;
  for(unsigned int ix=0;ix<3;++ix) output << " " << outgoing_[ix];
  output << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}

string DalitzBase::addChannel(string arg) {
  string error;
  DalitzResonancePtr res = DalitzResonance::readResonance(arg,error);
  if (res)
    resonances_.push_back(res);
  else
    return error;
  // success
  return "";
}

string DalitzBase::setExternal(string arg) {
  // parse first bit of the string
  string stype = StringUtils::car(arg);
  arg          = StringUtils::cdr(arg);
  // extract PDG code for the incoming particle
  long id = stoi(stype);
  tPDPtr in = getParticleData(id);
  if(!in)
    return "Incoming particle with id " + std::to_string(id) + "does not exist";
  tPDVector out;
  for(unsigned int ix=0;ix<3;++ix) {
    string stype = StringUtils::car(arg);
    arg          = StringUtils::cdr(arg);
    long in = stoi(stype);
    tPDPtr pData = getParticleData(in);
    if(!pData)
      return "Outgoing particle with id " + std::to_string(in) + "does not exist";
    out.push_back(pData);
  }
  incoming_ = in->id();
  outgoing_ = {out[0]->id(),out[1]->id(),out[2]->id()};
  // success
  return "";
}

int DalitzBase::modeNumber(bool & cc,tcPDPtr parent,
				const tPDVector & children) const {
  // must be three decay products
  if(children.size()!=3) return -1;
  // ids of the outgoing particles
  map<long,unsigned int> ids;
  for(unsigned int iy=0;iy<3;++iy) {
    long id = outgoing_[iy];
    if(useAllK0_ && (id==130 || id==310 || id==311 || id==-311)) id=130; 
    if(ids.find(id)!=ids.end())
      ids[id]+=1;
    else
      ids[id]+=1;
  }
  if(incoming_==parent->id()) {
    bool found=true;
    map<long,unsigned int> ids2;
    for(unsigned int iy=0;iy<children.size();++iy) {
    long id = children[iy]->id();
    if(useAllK0_ && (id==130 || id==310 || id==311 || id==-311)) id=130; 
      if(ids2.find(id)!=ids2.end())
	ids2[id]+=1;
      else
	ids2[id] =1;
    }
    for (const auto& kv : ids2) {
      if(ids.find(kv.first)==ids.end() ||
	 ids[kv.first]!=kv.second) {
	found=false;
	break;
      }
    }
    if (found) {
      cc = false;
      return 0;
    }
  }
  if( (!parent->CC() && incoming_==parent->id()) ||
      ( parent->CC() && incoming_==parent->CC()->id()) ) {
    map<long,unsigned int> ids2;
    for(unsigned int iy=0;iy<3;++iy) {
      tPDPtr part = children[iy]->CC() ? children[iy]->CC() : children[iy];
      long id = part->id();
      if(useAllK0_ && (id==130 || id==310 || id==311 || id==-311)) id=130; 
      if(ids2.find(id)!=ids2.end())
	ids2[id]+=1;
      else
	ids2[id]+=1;
    }
    bool found=true;
    for (const auto& kv : ids2) {
      if(ids.find(kv.first)==ids.end() ||
	 ids[kv.first]!=kv.second) {
	found=false;
	break;
      }
    }
    if (found) {
      cc = true;
      return 0;
    }
  }
  return -1;
}
#line 1 "./ScalarTo3ScalarDalitz.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ScalarTo3ScalarDalitz class.
//

#include "ScalarTo3ScalarDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;
using namespace ThePEG::Helicity;

IBPtr ScalarTo3ScalarDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr ScalarTo3ScalarDalitz::fullclone() const {
  return new_ptr(*this);
}

void ScalarTo3ScalarDalitz::persistentOutput(PersistentOStream & os) const {
  os << useResonanceMass_ ;
}

void ScalarTo3ScalarDalitz::persistentInput(PersistentIStream & is, int) {
  is >> useResonanceMass_;
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<ScalarTo3ScalarDalitz,DalitzBase>
describeHerwigScalarTo3ScalarDalitz("Herwig::ScalarTo3ScalarDalitz", "HwDalitzDecay.so");

void ScalarTo3ScalarDalitz::Init() {

  static ClassDocumentation<ScalarTo3ScalarDalitz> documentation
    ("The ScalarTo3ScalarDalitz class provides a base class for "
     "weak three-body decays of bottom and charm mesons");

  static Switch<ScalarTo3ScalarDalitz,bool> interfaceResonanceMass
    ("ResonanceMass",
     "Whether to use the kinematic mass or the resonance pole mass for the denominator in kinematic expressions",
     &ScalarTo3ScalarDalitz::useResonanceMass_, false, false, false);
  static SwitchOption interfaceResonanceMassYes
    (interfaceResonanceMass,
     "Yes",
     "Use the resonance mass, to be avoided only use if do in experimental fit",
     true);
  static SwitchOption interfaceResonanceMassNo
    (interfaceResonanceMass,
     "No",
     "Use the correct kinematic mass",
     false);

}

void ScalarTo3ScalarDalitz::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  // set up the spin information for the decay products
  ScalarWaveFunction::constructSpinInfo(const_ptr_cast<tPPtr>(&part),
					incoming,true);
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double ScalarTo3ScalarDalitz::me2(const int ichan, const Particle & part,
			    const tPDVector & ,
			    const vector<Lorentz5Momentum> & momenta,
			    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin0,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    ScalarWaveFunction::
      calculateWaveFunctions(rho_,const_ptr_cast<tPPtr>(&part),incoming);
  }
  // set the kinematics
  mD_ = part.mass();
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    mOut_[ix]=momenta[ix].mass();
    for(unsigned int iy=ix+1;iy<momenta.size();++iy) {
      m2_[ix][iy]=(momenta[ix]+momenta[iy]).m();
      m2_[iy][ix]=m2_[ix][iy];
    }
  }
  // now compute the matrix element
  Complex amp = amplitude(ichan);
  (*ME())(0,0,0,0) = amp;
  return norm(amp);
}

Complex ScalarTo3ScalarDalitz::resAmp(unsigned int i) const {
  Complex output = resonances()[i]->amp;
  if (resonances()[i]->type==ResonanceType::NonResonant) return output;
  // mass of the resonance
  const Energy & mR = resonances()[i]->mass ;
  // locations of the outgoing particles
  const unsigned int &d1 = resonances()[i]->daughter1;
  const unsigned int &d2 = resonances()[i]->daughter2;
  const unsigned int &sp = resonances()[i]->spectator;
  // compute the Breit-Wigner times resonance form factor piece
  output *= resonances()[i]->BreitWigner(m2_[d1][d2],mOut_[d1],mOut_[d2]);
  // angular piece
  // mass for the denominator
  Energy mDen = useResonanceMass_ ? resonances()[i]->mass : m2_[d1][d2];
  // denominator for the older form of the amplitude
  Energy2 denom = GeV2;
  if (resonances()[i]->type/10 == 1 ) { 
    Energy2 pa2 = 0.25*(sqr(m2_[d1][d2])-2.*(sqr(mOut_[d1])+sqr(mOut_[d2])) + sqr(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(m2_[d1][d2]));
    Energy2 pc2 = 0.25*(sqr(m2_[d1][d2])-2.*(sqr(mD_      )+sqr(mOut_[sp])) + sqr(sqr(mD_      )-sqr(mOut_[sp]))/sqr(m2_[d1][d2]));
    denom = 4.*sqrt(pa2*pc2);
  }
  // vectors
  if(abs(resonances()[i]->type)%10==3) {
    output *= (sqr(m2_[d2][sp])-sqr(m2_[d1][sp])
		  + (  sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/sqr(mDen) )/denom;
  }
  else if(abs(resonances()[i]->type)%10==5) {
    output *= 1./sqr(denom)*( sqr( sqr(m2_[d2][sp]) - sqr(m2_[d1][sp]) +(sqr(mD_)-sqr(mOut_[sp]))*(sqr(mOut_[d1])-sqr(mOut_[d2]))/(sqr(mDen))) -
			      (sqr(m2_[d1][d2])-2*      sqr(mD_)-2*sqr(mOut_[sp]) + sqr((sqr(      mD_) - sqr(mOut_[sp]))/mDen))*
			      (sqr(m2_[d1][d2])-2*sqr(mOut_[d1])-2*sqr(mOut_[d2]) + sqr((sqr(mOut_[d1]) - sqr(mOut_[d2]))/mDen))/3.);
  }
  // spin zero and non-resonant is done now
  if((abs(resonances()[i]->type)%10==1 && resonances()[i]->type != ResonanceType::Spin0Gauss ) ||
     abs(resonances()[i]->type/10)==10) {
    return output;
  }
  // spin piece x Blatt-Weisskopf for parent 
  else {
    double fD=1.;
    // for the D decay
    Energy pD  = sqrt(max(ZERO,(0.25*sqr(sqr(mD_)-sqr(mR)-sqr(mOut_[sp])) - sqr(mR*mOut_[sp]))/sqr(mD_)));
    Energy pDAB= sqrt( 0.25*sqr(sqr(mD_)-sqr(m2_[d1][d2])-sqr(mOut_[sp])) - sqr(m2_[d1][d2]*mOut_[sp]))/mD_;
    double r2A(parentRadius()   *pD),r2B(parentRadius()   *pDAB);
    // Blatt-Weisskopf factors and spin piece
    switch (resonances()[i]->type) {
    case ResonanceType::Spin0Gauss:
      fD = exp(-(r2B-r2A)/12.);
      break;
    case ResonanceType::Spin1: case ResonanceType::Spin1E691 : case ResonanceType::Spin1GS : 
      fD=sqrt( (1. + sqr(r2A)) / (1. + sqr(r2B)) );
      break;
    case ResonanceType::Spin2: case ResonanceType::Spin2E691:
      fD = sqrt( (9. + sqr(r2A)*(3.+sqr(r2A))) / (9. + sqr(r2B)*(3.+sqr(r2B))));
      break;
    default :
      assert(false);
    }
    output *= fD;
    return output;
  }
}

void ScalarTo3ScalarDalitz::dataBaseOutput(ofstream & output, bool header) const {
  if(header) output << "update decayers set parameters=\"";
  // parameters for the DalitzBase base class
  DalitzBase::dataBaseOutput(output,false);
  output << "newdef " << name() << ":ResonanceMass " << useResonanceMass_ << "\n";
  output << "\n";
  if(header) {
    output << "\n\" where BINARY ThePEGName=\"" << fullName() << "\";" << endl;
  }
}
#line 1 "./VectorTo3PseudoScalarDalitz.cc"
// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorTo3PseudoScalarDalitz class.
//

#include "VectorTo3PseudoScalarDalitz.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Helicity/epsilon.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "Herwig/PDT/ThreeBodyAllOnCalculator.h"
#include "Herwig/Decay/GeneralDecayMatrixElement.h"

using namespace Herwig;

IBPtr VectorTo3PseudoScalarDalitz::clone() const {
  return new_ptr(*this);
}

IBPtr VectorTo3PseudoScalarDalitz::fullclone() const {
  return new_ptr(*this);
}

void VectorTo3PseudoScalarDalitz::persistentOutput(PersistentOStream & os) const {
  os << ounit(coupling_,1./GeV);
}

void VectorTo3PseudoScalarDalitz::persistentInput(PersistentIStream & is, int) {
  is >> iunit(coupling_,1./GeV);
}

// The following static variable is needed for the type
// description system in ThePEG.
DescribeClass<VectorTo3PseudoScalarDalitz,DalitzBase>
describeHerwigVectorTo3PseudoScalarDalitz("Herwig::VectorTo3PseudoScalarDalitz", "HwDalitzDecay.so");

void VectorTo3PseudoScalarDalitz::Init() {

  static ClassDocumentation<VectorTo3PseudoScalarDalitz> documentation
    ("The VectorTo3PseudoScalarDalitz class provides a base class "
     "for the decay of vector mesons to 3 pseudoscalar mesons");

  static Parameter<VectorTo3PseudoScalarDalitz,InvEnergy> interfaceCouopling
    ("Coupling",
     "The coupling for the normalisation of the mode",
     &VectorTo3PseudoScalarDalitz::coupling_, 1./GeV, 1./GeV, 0./GeV, 1000./GeV,
     false, false, Interface::limited);

}

void VectorTo3PseudoScalarDalitz::
constructSpinInfo(const Particle & part, ParticleVector decay) const {
  VectorWaveFunction::constructSpinInfo(vectors_,const_ptr_cast<tPPtr>(&part),
					incoming,true,false);
  // set up the spin information for the decay products
  for(unsigned int ix=0;ix<3;++ix)
    ScalarWaveFunction::constructSpinInfo(decay[ix],outgoing,true);
}

double VectorTo3PseudoScalarDalitz::me2(const int ichan, const Particle & part,
				    const tPDVector & ,
				    const vector<Lorentz5Momentum> & momenta,
				    MEOption meopt) const {
  if(!ME())
    ME(new_ptr(GeneralDecayMatrixElement(PDT::Spin1,PDT::Spin0,PDT::Spin0,PDT::Spin0)));
  useMe();
  if(meopt==Initialize) {
    VectorWaveFunction::calculateWaveFunctions(vectors_,rho_,
  						const_ptr_cast<tPPtr>(&part),
  						incoming,false);
  }
  // set the kinematics
  mD_ = part.mass();
  for(unsigned int ix=0;ix<momenta.size();++ix) {
    mOut_[ix]=momenta[ix].mass();
    for(unsigned int iy=ix+1;iy<momenta.size();++iy) {
      m2_[ix][iy]=(momenta[ix]+momenta[iy]).m();
      m2_[iy][ix]=m2_[ix][iy];
    }
  }
  // now compute the matrix element
  complex<InvEnergy2> amp = amplitude(ichan);
  // polarization vector piece
  LorentzVector<complex<Energy3> > scalar = epsilon(momenta[0],momenta[1],momenta[2]);
  // compute the matrix element
  for(unsigned int ix=0;ix<3;++ix) {
    (*ME())(ix,0,0,0) = Complex(coupling_*amp*scalar.dot(vectors_[ix]));
  }
  // return the answer
  return (ME()->contract(rho_)).real();
}

double VectorTo3PseudoScalarDalitz::
threeBodyMatrixElement(const int , const Energy2 q2,
		       const  Energy2 s3, const Energy2 s2, const Energy2 s1, const 
		       Energy m1, const Energy m2, const Energy m3) const {
  mD_ = sqrt(q2);
  mOut_[0] = m1;
  mOut_[1] = m2;
  mOut_[2] = m3;
  m2_[0][1]=m2_[1][0]=sqrt(s3);
  m2_[0][2]=m2_[2][0]=sqrt(s2);
  m2_[1][2]=m2_[2][1]=sqrt(s1);
  // now compute the matrix element
  // amplitide
  complex<InvEnergy2> amp = amplitude(-1);
  // epsilon piece
  Energy6 kin = (pow<4,1>(m1)*(-2*(sqr(m2) + sqr(m3)) + s1) + pow<4,1>(m2)*(-2*sqr(m3) + s2) +
		 s3*(pow<4,1>(m3) + s1*s2 - sqr(m3)*(s1 + s2 + s3)) - 
     sqr(m1)*(2*pow<4,1>(m2) + 2*pow<4,1>(m3) + sqr(m2)*(4*sqr(m3) - 3*(s1 + s2) - s3) + s1*(s1 + s2 + s3) - 
        sqr(m3)*(3*s1 + s2 + 3*s3)) - sqr(m2)*(2*pow<4,1>(m3) + s2*(s1 + s2 + s3) - sqr(m3)*(s1 + 3*(s2 + s3))))/12.;
  return norm(amp*coupling_*GeV*GeV2)*kin/GeV2/GeV2/GeV2;
} 

WidthCalculatorBasePtr 
VectorTo3PseudoScalarDalitz::threeBodyMEIntegrator(const DecayMode & ) const {
  int imode=0;
  // construct the integrator
  vector<double> inweights;
  vector<int> intype;
  vector<Energy> inmass,inwidth;
  inweights.reserve(resonances().size());
  intype.reserve(resonances().size());
  inmass.reserve(resonances().size());
  inwidth.reserve(resonances().size());
  int iloc=-1;
  vector<double> inpow(2,0.0);
  for(unsigned int ix=0;ix<resonances().size();++ix) {
    tPDPtr resonance = getParticleData(resonances()[ix]->id);
    if(resonance) {
      ++iloc;
      inweights.push_back(weights()[iloc]);
      inmass.push_back(resonances()[ix]->mass);
      inwidth.push_back(abs(resonances()[ix]->width));
      intype.push_back(resonances()[ix]->spectator+1);
    }
  }
  tcDecayIntegratorPtr decayer(this);
  WidthCalculatorBasePtr output(
    new_ptr(ThreeBodyAllOnCalculator<VectorTo3PseudoScalarDalitz>
  	    (inweights,intype,inmass,inwidth,inpow,
  	     *this,imode,
	     mode(0)->outgoing()[0]->mass(),
	     mode(0)->outgoing()[1]->mass(),
	     mode(0)->outgoing()[2]->mass())));
  return output;
}

void VectorTo3PseudoScalarDalitz::dataBaseOutput(ofstream & output,
					     bool header) const {
  if(header){output << "update decayers set parameters=\"";}
  output << "newdef " << name() << ":Coupling " << coupling_*GeV << "\n";
  // parameters for the DalitzBase base class
  DalitzBase::dataBaseOutput(output,false);
  if(header){output << "\n\" where BINARY ThePEGName=\"" 
   		    << fullName() << "\";" << endl;}
}

complex<InvEnergy2> VectorTo3PseudoScalarDalitz::resAmp(unsigned int i) const {
  // can't have a scalar here on spin/parity grounds
  assert(resonances()[i]->type%10!=1);
  // shouldn't have E691 stuff either
  assert(resonances()[i]->type%10!=1);
  // amplitude
  Complex output = resonances()[i]->amp;
  if (resonances()[i]->type==ResonanceType::NonResonant) return output/GeV2;
  // mass of the resonance
  const Energy & mR = resonances()[i]->mass ;
  // locations of the outgoing particles
  const int &d1 = resonances()[i]->daughter1;
  const int &d2 = resonances()[i]->daughter2;
  const int &sp = resonances()[i]->spectator;
  // epsilon piece  =eps(d1,d2,sp)
  double sign = (sp-d1)*(d1-d2)*(d2-sp)/2.;
  // compute the Breit-Wigner times resonance form factor piece
  output *= sign*resonances()[i]->BreitWigner(m2_[d1][d2],mOut_[d1],mOut_[d2]);
  // Blatt-Weisskopf factors
  // for the D decay
  Energy pD  = sqrt(max(ZERO,(0.25*sqr(sqr(mD_)-sqr(mR)-sqr(mOut_[sp])) - sqr(mR*mOut_[sp]))/sqr(mD_)));
  Energy pDAB= sqrt( 0.25*sqr(sqr(mD_)-sqr(m2_[d1][d2])-sqr(mOut_[sp])) - sqr(m2_[d1][d2]*mOut_[sp]))/mD_;
  double r2A(parentRadius()   *pD),r2B(parentRadius()   *pDAB);
  // Blatt-Weisskopf factors and spin piece
  switch (resonances()[i]->type) {
  case ResonanceType::Spin1: case ResonanceType::Spin1GS : 
    output *= sqrt( (1. + sqr(r2A)) / (1. + sqr(r2B)) );
    break;
  case ResonanceType::Spin2:
    output *= sqrt( (9. + sqr(r2A)*(3.+sqr(r2A))) / (9. + sqr(r2B)*(3.+sqr(r2B))));
    // spin piece
    output *= (sqr(mD_) - sqr(mOut_[sp]) + sqr(m2_[d1][d2]))/(sqr(mD_)*sqr(m2_[d1][d2]))/GeV2*
	       ((-sqr(m2_[sp][d1]) + sqr(m2_[sp][d2]))*sqr(m2_[d1][d2]) +
		(mD_ - mOut_[sp])*(mD_ + mOut_[sp])*(mOut_[d1] - mOut_[d2])*(mOut_[d1] + mOut_[d2]));
    break;
  default :
    assert(false);
  }
  return output/GeV2;
}
