// -*- C++ -*-
// Native gg -> gamma g quark-box process. GPL version 3; see COPYING.
#include "MEPP2GammaJetLoop.h"
#include "GammaJetLoopAmplitudes.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "Herwig/Models/StandardModel/StandardModel.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/MatrixElement/Tree2toNDiagram.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/SubProcess.h"

using namespace Herwig;

MEPP2GammaJetLoop::MEPP2GammaJetLoop(): amplitudeMode_(2),scalePreFactor_(1.) {
  massOption(vector<unsigned int>(2,0));
}

Energy2 MEPP2GammaJetLoop::scale() const {
  const Energy2 s=sHat(),t=tHat(),u=uHat();
  return scalePreFactor_*2.*s*t*u/(s*s+t*t+u*u);
}

void MEPP2GammaJetLoop::getDiagrams() const {
  tcPDPtr g=getParticleData(ParticleID::g),a=getParticleData(ParticleID::gamma);
  // Dummy s-channel topology, as in the native gg -> Higgs g loop ME.
  add(new_ptr((Tree2toNDiagram(2),g,g,1,g,3,a,3,g,-1)));
}

Selector<MEBase::DiagramIndex>
MEPP2GammaJetLoop::diagrams(const DiagramVector &diags) const {
  Selector<DiagramIndex> sel;
  for(DiagramIndex i=0;i<diags.size();++i) sel.insert(1.,i);
  return sel;
}

Selector<const ColourLines *>
MEPP2GammaJetLoop::colourGeometries(tcDiagPtr) const {
  static const ColourLines forward("1 3 5, -5 -3 -2, 2 -1");
  static const ColourLines backward("-1 -3 -5, 5 3 2, -2 1");
  Selector<const ColourLines *> sel;
  sel.insert(0.5,&forward);sel.insert(0.5,&backward);
  return sel;
}

pair<RhoDMatrix,RhoDMatrix> MEPP2GammaJetLoop::incomingDensities() const {
  auto rho=getRhoMatrices();
  for(RhoDMatrix *r:{&rho.first,&rho.second}) {
    RhoDMatrix transverse(PDT::Spin1,false);
    const double trace=real((*r)(0,0)+(*r)(2,2));
    if(trace>0.) for(unsigned i:{0u,2u}) for(unsigned j:{0u,2u})
      transverse(i,j)=(*r)(i,j)/trace;
    *r=transverse;
  }
  return rho;
}

ProductionMatrixElement MEPP2GammaJetLoop::amplitudes(
  const std::array<Lorentz5Momentum,4>&mom,
  const std::array<vector<VectorWaveFunction>,4>&waves) const {
  ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin1,PDT::Spin1);
  GammaJetLoop::Momenta p;
  for(unsigned i=0;i<4;++i) {
    const double sign=i<2?1.:-1.;
    p[i]={{sign*double(mom[i].e()/GeV),sign*double(mom[i].x()/GeV),
           sign*double(mom[i].y()/GeV),sign*double(mom[i].z()/GeV)}};
  }
  GammaJetLoop::HelicityFunctions h{};
  if(amplitudeMode_!=0) {
    const auto s=real(2.L*GammaJetLoop::dot(p[0],p[1]));
    const auto t=real(2.L*GammaJetLoop::dot(p[1],p[2]));
    const auto u=real(2.L*GammaJetLoop::dot(p[0],p[2]));
    if(amplitudeMode_==1) {
      h=GammaJetLoop::massless(s,t,u);
      for(auto &a:h) a/=3.L; // sum Q(u,d,s,c,b)=1/3; no top
    } else {
      // The exactly massless u,d,s contribution cancels coherently.
      for(long flavour=4;flavour<=6;++flavour) {
        tcPDPtr q=getParticleData(flavour);
        auto f=GammaJetLoop::massive(s,t,u,double(q->mass()/GeV));
        for(unsigned j=0;j<5;++j) h[j]+=double(q->iCharge())/3.L*f[j];
      }
    }
  }
  // M = 2 sqrt(alphaEM alphaS^3) d^{abc} sum Qq hq. Colour averaging
  // is applied once in me2(); the hard vertex needs only the spin tensor.
  const double normalization=2.*sqrt(SM().alphaEM(ZERO)*pow(SM().alphaS(scale()),3));
  for(unsigned a:{0u,2u}) for(unsigned b:{0u,2u})
    for(unsigned c:{0u,2u}) for(unsigned d:{0u,2u}) {
      GammaJetLoop::Complex value=1.L;
      if(amplitudeMode_!=0) {
        std::array<GammaJetLoop::Vector,4> eps;
        unsigned hel[4]={a,b,c,d};
        for(unsigned i=0;i<4;++i) {
          const auto &w=waves[i][hel[i]];
          eps[i]={{w.t(),w.x(),w.y(),w.z()}};
        }
        value=GammaJetLoop::contract(p,eps,h);
      }
      me(a,b,c,d)=normalization*Complex(double(value.real()),double(value.imag()));
    }
  return me;
}

double MEPP2GammaJetLoop::me2() const {
  std::array<Lorentz5Momentum,4> mom;
  std::array<vector<VectorWaveFunction>,4> waves;
  for(unsigned i=0;i<4;++i) {
    mom[i]=meMomenta()[i];
    VectorWaveFunction::calculateWaveFunctions(waves[i],mom[i],mePartonData()[i],
                                               i<2?incoming:outgoing,true);
  }
  auto me=amplitudes(mom,waves);
  const auto rho=incomingDensities();
  // sum_abc d^{abc}d^{abc}/(8*8) = (40/3)/64 = 5/24.
  return 5./24.*me.average(rho.first,rho.second);
}

void MEPP2GammaJetLoop::constructVertex(tSubProPtr sub) {
  ParticleVector hard{sub->incoming().first,sub->incoming().second,
                      sub->outgoing()[0],sub->outgoing()[1]};
  if(hard[2]->id()!=ParticleID::gamma) swap(hard[2],hard[3]);
  std::array<Lorentz5Momentum,4> mom;
  std::array<vector<VectorWaveFunction>,4> waves;
  for(unsigned i=0;i<4;++i) {
    mom[i]=hard[i]->momentum();
    VectorWaveFunction(waves[i],hard[i],i<2?incoming:outgoing,i>=2,true,true);
  }
  auto rho=incomingDensities();
  hard[0]->spinInfo()->rhoMatrix(rho.first);
  hard[1]->spinInfo()->rhoMatrix(rho.second);
  for(unsigned i=2;i<4;++i) {
    RhoDMatrix d(PDT::Spin1,false);d(0,0)=d(2,2)=0.5;
    hard[i]->spinInfo()->DMatrix(d);
  }
  HardVertexPtr vertex=new_ptr(HardVertex());
  vertex->ME(amplitudes(mom,waves));
  for(auto part:hard) part->spinInfo()->productionVertex(vertex);
}

void MEPP2GammaJetLoop::persistentOutput(PersistentOStream&os) const {
  os<<amplitudeMode_<<scalePreFactor_;
}
void MEPP2GammaJetLoop::persistentInput(PersistentIStream&is,int) {
  is>>amplitudeMode_>>scalePreFactor_;
}
DescribeClass<MEPP2GammaJetLoop,HwMEBase>
 describeHerwigMEPP2GammaJetLoop("Herwig::MEPP2GammaJetLoop","HwMEHadron.so");
void MEPP2GammaJetLoop::Init() {
  static ClassDocumentation<MEPP2GammaJetLoop> doc
    ("Native finite quark-loop gg -> gamma g with spin correlations. "
     "Massive includes massless u,d,s and massive c,b,t. Use a QCD-only shower.");
  static Switch<MEPP2GammaJetLoop,unsigned int> mode
    ("AmplitudeMode","Quark-box amplitudes; Flat is a diagnostic only.",
     &MEPP2GammaJetLoop::amplitudeMode_,2,false,false);
  static SwitchOption flat(mode,"Flat","Unit reduced transverse amplitudes (diagnostic).",0);
  static SwitchOption massless(mode,"Massless5","Massless u,d,s,c,b; top omitted.",1);
  static SwitchOption massive(mode,"Massive","Massless u,d,s; particle-data c,b,t masses.",2);
  static Parameter<MEPP2GammaJetLoop,double> scale
    ("ScalePreFactor","Prefactor of 2 s t u / (s^2+t^2+u^2).",
     &MEPP2GammaJetLoop::scalePreFactor_,1.,0.,10.,false,false,Interface::limited);
}
