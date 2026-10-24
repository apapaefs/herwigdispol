// Loaded by Herwig read: exercise the actual HwShower methods, not copies.
#include "Herwig/Shower/QTilde/QTildeShowerHandler.h"
#include "Herwig/Shower/QTilde/Kinematics/KinematicsReconstructor.h"
#include "Herwig/Shower/QTilde/Base/PartnerFinder.h"
#include "Herwig/PDF/HwRemDecayer.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/HalfHalfOneSplitFn.h"
#include "Herwig/Shower/QTilde/SplittingFunctions/OneHalfHalfSplitFn.h"
#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/PDF/LHAPDF6.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ThePEG;
using namespace Herwig;

namespace {
int checks = 0, failures = 0;
void check(bool good, const std::string & label) {
  ++checks;
  if (!good) { ++failures; std::cout << "CHECK_FAIL " << label << '\n'; }
}
void near(double got, double expected, const std::string & label,
          double tolerance=2.e-12) {
  check(std::isfinite(got) && std::abs(got-expected) <=
        tolerance*std::max(1.,std::abs(expected)), label);
  if (!std::isfinite(got) || std::abs(got-expected) >
      tolerance*std::max(1.,std::abs(expected)))
    std::cout << "  got=" << got << " expected=" << expected << '\n';
}

struct Query { long id; double x, q2; };
class RecordingPDF : public PDFBase {
public:
  int kind = 0; // 0: xf, 1: longitudinal difference, 2: transverse difference
  bool constant = false;
  mutable std::vector<Query> queries;
  bool canHandleParticle(tcPDPtr) const override { return true; }
  cPDVector partons(tcPDPtr) const override { return {}; }
  IBPtr clone() const override { return new_ptr(*this); }
  IBPtr fullclone() const override { return new_ptr(*this); }
  static double shape(long id, double x) {
    return (id==21 ? 0.7 : 1.+0.02*std::abs(id))*std::pow(1.-x,3);
  }
  double ratio(long id, double x) const {
    if (!kind) return 1.;
    if (constant) return kind==1 ? 0.23 : 0.11;
    if (kind==2) return id==1 ? -0.08+0.10*x : 0.10+0.15*x;
    if (id==21) return 0.12+0.30*x;
    return id==1 ? -0.18+0.10*x : 0.24+0.20*x;
  }
  // Deliberately override only xfx. The real inherited xfl performs exp(-l).
  double xfx(tcPDPtr, tcPDPtr parton, Energy2 q2, double x,
             double=0., Energy2=ZERO) const override {
    queries.push_back({parton->id(),x,q2/GeV2});
    return shape(parton->id(),x)*ratio(parton->id(),x);
  }
};
class FixtureHandler : public QTildeShowerHandler {
public:
  using ShowerHandler::setCurrentHandler;
  using ShowerHandler::unSetCurrentHandler;
};
struct HandlerScope {
  Ptr<FixtureHandler>::pointer handler;
  HandlerScope() : handler(new_ptr(FixtureHandler())) {
    if (ShowerHandler::currentHandlerIsSet())
      throw std::runtime_error("Run the fixture in its own Herwig read process");
    near(handler->factorizationScaleFactor(),1.,"handler scale default");
    handler->setCurrentHandler();
  }
  ~HandlerScope() { handler->unSetCurrentHandler(); }
};
template<class Base> class ExposedSplit : public Base {
public:
  using SudakovFormFactor::setPDF;
  using SudakovFormFactor::PDFVetoRatio;
  using SudakovFormFactor::pdfMax;
  using SudakovFormFactor::pdfFactor;
};
class RealPDF : public LHAPDF {
public:
  void select(const std::string & name) { setPDFName(name); initPDFptr(); }
};

void expectQueries(const RecordingPDF & pdf, const std::vector<Query> & expected,
                   const std::string & name) {
  check(pdf.queries.size()==expected.size(),name+" query count");
  for (size_t i=0; i<std::min(pdf.queries.size(),expected.size()); ++i) {
    const auto & got=pdf.queries[i]; const auto & want=expected[i];
    check(got.id==want.id,name+" flavour "+std::to_string(i));
    near(got.x,want.x,name+" x "+std::to_string(i));
    near(got.q2,want.q2,name+" scale "+std::to_string(i));
  }
}

template<class Split>
void vetoCase(const std::string & label, tcBeamPtr beam, tcPDPtr parent,
              tcPDPtr child, bool flavourChange, double x, double z,
              const ThreeVector<double> & polarization, const RhoDMatrix & rho,
              bool constant=false, bool enabled=true, double t=100.,
              double scaleFactor=1.25) {
  auto u=new_ptr(RecordingPDF());
  auto l=new_ptr(RecordingPDF()); l->kind=1; l->constant=constant;
  auto tr=new_ptr(RecordingPDF()); tr->kind=2; tr->constant=constant;
  auto split=new_ptr(ExposedSplit<Split>());
  split->setPDF(u,2.*GeV,enabled,polarization,l,tr);
  near(split->pdfMax(),35.,label+" PDFmax");
  check(split->pdfFactor()==0,label+" PDFFactor");
  const double q2=std::max(4.,t*scaleFactor*scaleFactor);
  const double xp=x/z;
  const double base=RecordingPDF::shape(parent->id(),xp)/RecordingPDF::shape(child->id(),x);
  const double d=(rho(1,1)-rho(0,0)).real();
  const Complex pt(polarization.x(),polarization.y());
  const double coherence=2.*(pt*rho(0,1)).real();
  const double numerator=1.+polarization.z()*d*l->ratio(parent->id(),xp)*
    (flavourChange ? (2.*z-1.)/(z*z+(1.-z)*(1.-z)) : 1.)+
    (flavourChange ? 0. : 2.*z/(1.+z*z)*coherence*tr->ratio(parent->id(),xp));
  const double denominator=1.+polarization.z()*d*l->ratio(child->id(),x)+
    coherence*tr->ratio(child->id(),x);
  const double expected=base*(enabled ? numerator/denominator : 1.)/35.;
  check(expected>0. && expected<1.,label+" controlled acceptance");
  try {
    const double got=split->PDFVetoRatio(t*GeV2,x,z,parent,child,beam,rho,scaleFactor);
    near(got,expected,label+" veto value");
  } catch (const std::exception & e) {
    check(false,label+" unexpected exception: "+e.what());
  }
  expectQueries(*u,{{parent->id(),xp,q2},{child->id(),x,q2}},label+" sum");
  const std::vector<Query> pair={{child->id(),x,q2},{parent->id(),xp,q2}};
  expectQueries(*l,enabled && polarization.z()!=0. ? pair : std::vector<Query>(),label+" longitudinal");
  expectQueries(*tr,enabled && std::norm(pt)!=0. ? pair : std::vector<Query>(),label+" transverse");
}

void matrixCase(const std::string & label, tcBeamPtr beam, tcPDPtr parton,
                double x, const ThreeVector<double> & p, bool constant=false) {
  auto u=new_ptr(RecordingPDF());
  auto l=new_ptr(RecordingPDF()); l->kind=1; l->constant=constant;
  auto tr=new_ptr(RecordingPDF()); tr->kind=2; tr->constant=constant;
  auto split=new_ptr(ExposedSplit<HalfHalfOneSplitFn>());
  split->setPDF(u,2.*GeV,true,p,l,tr);
  const auto h=split->calculateHMatrix(beam,parton,x,169.*GeV2);
  const unsigned int hi=parton->iSpin()-1;
  near(h(0,0).real(),.5*(1.-p.z()*l->ratio(parton->id(),x)),label+" H00");
  near(h(hi,hi).real(),.5*(1.+p.z()*l->ratio(parton->id(),x)),label+" Hhh");
  const Complex off=.5*Complex(p.x(),p.y())*tr->ratio(parton->id(),x);
  near(h(0,hi).real(),off.real(),label+" offdiag real");
  near(h(0,hi).imag(),off.imag(),label+" offdiag imag");
  near(h(hi,0).imag(),-off.imag(),label+" conjugate");
  if (hi==2) near(h(1,1).real(),0.,label+" massless vector middle");
  const std::vector<Query> query={{parton->id(),x,169.}};
  expectQueries(*u,query,label+" sum");
  expectQueries(*l,p.z()!=0. ? query : std::vector<Query>(),label+" longitudinal");
  expectQueries(*tr,p.x()!=0. || p.y()!=0. ? query : std::vector<Query>(),label+" transverse");
}

void savedState(tcBeamPtr beam, tcPDPtr parton) {
  auto u=new_ptr(RealPDF()); u->select("NNPDF40_nlo_pch_as_01180");
  auto l=new_ptr(RealPDF()); l->select("NNPDFpol20_nlo_as_01180");
  auto split=new_ptr(ExposedSplit<HalfHalfOneSplitFn>());
  split->setPDF(u,1.*GeV,true,ThreeVector<double>(0.,0.,-1.),l,tcPDFPtr());
  const double x=.81157378587876683, z=.9956721837705923;
  const Energy2 scale=60153.630625064412*GeV2;
  const double old=u->xfx(beam,parton,scale,x);
  const double next=u->xfx(beam,parton,scale,x/z);
  const double rc=l->xfx(beam,parton,scale,x)/old;
  const double rp=l->xfx(beam,parton,scale,x/z)/next;
  near(rc,.0322373798549054,"saved child ratio",1.e-10);
  near(rp,.0363847351984158,"saved parent ratio",1.e-10);
  near(l->xfl(beam,parton,scale,-std::log(x))/old,rc,"API coordinate equivalence");
  const auto h=split->calculateHMatrix(beam,parton,x,scale);
  near(h(1,1).real(),.483881310072547,"saved old density",1.e-10);
  RhoDMatrix rho(PDT::Spin1Half); rho(0,0)=0.; rho(1,1)=1.;
  const double spin=(1.-rp)/(1.-rc);
  near(spin,.995714491077483,"saved analytic spin factor",1.e-10);
  try {
    const double got=split->PDFVetoRatio(scale,x,z,parton,parton,beam,rho,1.);
    near(got/(next/old/35.),spin,"saved real PDF veto",1.e-10);
  } catch (const std::exception & e) {
    check(false,std::string("saved state threw: ")+e.what());
  }
}

void run() {
  HandlerScope scope;
  auto beam=dynamic_ptr_cast<Ptr<BeamParticleData>::pointer>(BeamParticleData::Create(2212,"p+"));
  auto quark=ParticleData::Create(2,"u"); quark->iSpin(PDT::Spin1Half);
  auto down=ParticleData::Create(1,"d"); down->iSpin(PDT::Spin1Half);
  auto gluon=ParticleData::Create(21,"g"); gluon->iSpin(PDT::Spin1);
  auto anti=ParticleData::Create(-2,"ubar"); anti->iSpin(PDT::Spin1Half);
  RhoDMatrix rho(PDT::Spin1Half);
  rho(0,0)=.65; rho(1,1)=.35; rho(0,1)=Complex(.09,.04); rho(1,0)=std::conj(rho(0,1));
  const ThreeVector<double> both(.24,-.17,.6), longitudinal(0.,0.,.6), transverse(.24,-.17,0.);
  for (double x : {.08,.27,.51}) {
    vetoCase<HalfHalfOneSplitFn>("combined "+std::to_string(x),beam,quark,quark,false,x,.77,both,rho);
    matrixCase("terminal "+std::to_string(x),beam,quark,x,both);
    matrixCase("parent "+std::to_string(x),beam,quark,x/.77,both);
    matrixCase("constant "+std::to_string(x),beam,quark,x,both,true);
  }
  vetoCase<HalfHalfOneSplitFn>("constant longitudinal",beam,quark,quark,false,.21,.73,longitudinal,rho,true);
  vetoCase<HalfHalfOneSplitFn>("negative difference",beam,down,down,false,.21,.73,longitudinal,rho);
  vetoCase<HalfHalfOneSplitFn>("transverse only",beam,quark,quark,false,.21,.73,transverse,rho);
  vetoCase<HalfHalfOneSplitFn>("real transverse",beam,quark,quark,false,.21,.73,ThreeVector<double>(.24,0.,0.),rho);
  matrixCase("real transverse density",beam,quark,.21,ThreeVector<double>(.24,0.,0.));
  vetoCase<HalfHalfOneSplitFn>("frozen scale",beam,quark,quark,false,.21,.73,both,rho,false,true,.25);
  vetoCase<OneHalfHalfSplitFn>("flavour change",beam,gluon,down,true,.21,.73,longitudinal,rho);
  vetoCase<HalfHalfOneSplitFn>("unpolarized bypass",beam,quark,quark,false,.21,.73,both,rho,false,false);
  vetoCase<HalfHalfOneSplitFn>("zero polarization",beam,quark,quark,false,.21,.73,ThreeVector<double>(0.,0.,0.),rho);
  matrixCase("gluon density",beam,gluon,.19,longitudinal);
  if (std::getenv("HERWIG_COORDINATE_REAL_PDFS")) savedState(beam,anti);
}
}

// Loaded after the ordinary Herwig repository and its dependencies.
__attribute__((constructor)) static void coordinateRegression() {
  std::cout << std::setprecision(17);
  try { run(); }
  catch (const std::exception & e) { check(false,std::string("fixture exception: ")+e.what()); }
  catch (...) { check(false,"fixture raised a nonstandard exception"); }
  std::cout << "PDF_COORDINATE_REGRESSION checks=" << checks
            << " failures=" << failures << std::endl;
}
