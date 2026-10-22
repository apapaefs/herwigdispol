// -*- C++ -*-
// Diagnostic analysis only: read the native hard vertex without rebuilding it.
// GPL version 3; see ../../src/herwig/COPYING.
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "MEPP2GammaJetLoop.h"
#include <fstream>
#include <iomanip>
#include <set>
#include <stdexcept>

namespace Herwig {
using namespace ThePEG;
using namespace ThePEG::Helicity;

class NativeHelicityAudit: public AnalysisHandler {
public:
  static void Init();
  void analyze(tEventPtr event, long, int loop, int state) override {
    if(!event || loop>0 || state) return;
    if(!seen_.insert(event->number()).second) return;
    auto sub=event->primarySubProcess();
    if(!sub || sub->outgoing().size()!=2)
      throw std::runtime_error("NativeHelicityAudit requires Born 2->2 events");
    auto native=dynamic_ptr_cast<Ptr<MEPP2GammaJetLoop>::transient_const_pointer>(sub->handler());
    if(!native) throw std::runtime_error("NativeHelicityAudit requires MEPP2GammaJetLoop");
    ParticleVector hard{sub->incoming().first,sub->incoming().second,
                        sub->outgoing()[0],sub->outgoing()[1]};
    if(hard[2]->id()!=22) std::swap(hard[2],hard[3]);
    if(hard[0]->id()!=21 || hard[1]->id()!=21 || hard[2]->id()!=22 || hard[3]->id()!=21)
      throw std::runtime_error("Unexpected native external-particle ordering");
    for(auto p:hard) if(!p->spinInfo())
      throw std::runtime_error("Missing native SpinInfo");
    auto vertex=dynamic_ptr_cast<Ptr<HardVertex>::transient_const_pointer>(
      hard[2]->spinInfo()->productionVertex());
    if(!vertex) throw std::runtime_error("Missing native HardVertex");
    for(auto p:hard) if(p->spinInfo()->productionVertex()!=vertex)
      throw std::runtime_error("The four particles do not share a HardVertex");

    std::ofstream out(generator()->filename()+".helicities",std::ios::app);
    if(!out) throw std::runtime_error("Cannot append helicity audit output");
    out<<std::setprecision(17)<<"EVENT "<<event->number()<<'\n';
    auto complex=[&out](Complex z){out<<' '<<z.real()<<' '<<z.imag();};
    const Energy2 mu2=native->scale();
    const auto sm=generator()->standardModel();
    out<<"PARAM "<<double(mu2/GeV2)<<' '<<sm->alphaS(mu2)<<' '<<sm->alphaEM(ZERO);
    for(long q=4;q<=6;++q) out<<' '<<double(getParticleData(q)->mass()/GeV);
    out<<'\n';
    const auto &amp=vertex->ME();
    out<<"ME2 "<<native->me2()<<' '<<native->lastXComb().lastME2()<<' '
       <<5./24.*amp.average(hard[0]->spinInfo()->rhoMatrix(),
                            hard[1]->spinInfo()->rhoMatrix())<<'\n';
    for(unsigned leg=0;leg<4;++leg) {
      const auto &p=hard[leg]->momentum();
      out<<"MOM "<<leg<<' '<<double(p.e()/GeV)<<' '<<double(p.x()/GeV)
         <<' '<<double(p.y()/GeV)<<' '<<double(p.z()/GeV)<<'\n';
      auto spin=dynamic_ptr_cast<Ptr<VectorSpinInfo>::transient_const_pointer>(hard[leg]->spinInfo());
      if(!spin) throw std::runtime_error("Missing native vector spin basis");
      // These are the actual stored production states, transformed with the
      // event into its current frame. Do not regenerate wavefunctions here.
      for(unsigned hel:{0u,2u}) {
        const auto &eps=spin->getCurrentBasisState(hel);
        out<<"POL "<<leg<<' '<<hel;
        complex(eps.t());complex(eps.x());complex(eps.y());complex(eps.z());out<<'\n';
      }
      if(leg<2) for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) {
        out<<"RHO "<<leg<<' '<<i<<' '<<j;
        complex(spin->rhoMatrix()(i,j));out<<'\n';
      }
    }
    for(unsigned a=0;a<3;++a) for(unsigned b=0;b<3;++b)
      for(unsigned c=0;c<3;++c) for(unsigned d=0;d<3;++d) {
        out<<"AMP "<<a<<' '<<b<<' '<<c<<' '<<d;
        complex(amp(a,b,c,d));out<<'\n';
      }
    for(unsigned leg=0;leg<2;++leg) {
      const auto rho=vertex->getRhoMatrix(leg,true);
      for(unsigned i=0;i<3;++i) for(unsigned j=0;j<3;++j) {
        out<<"OUTRHO "<<leg<<' '<<i<<' '<<j;complex(rho(i,j));out<<'\n';
      }
    }
    out<<"END\n";
  }
protected:
  void doinitrun() override {
    AnalysisHandler::doinitrun();seen_.clear();
    std::ofstream out(generator()->filename()+".helicities");
    if(!out) throw std::runtime_error("Cannot create helicity audit output");
  }
  IBPtr clone() const override {return new_ptr(*this);}
  IBPtr fullclone() const override {return new_ptr(*this);}
private:
  std::set<long> seen_;
};

DescribeNoPIOClass<NativeHelicityAudit,AnalysisHandler>
 describeNativeHelicityAudit("Herwig::NativeHelicityAudit","NativeHelicityAudit.so");
void NativeHelicityAudit::Init() {
  static ClassDocumentation<NativeHelicityAudit> doc(
    "Diagnostic export of the actual MEPP2GammaJetLoop HardVertex amplitudes, "
    "stored polarization bases, couplings and density matrices. Born events only.");
}
}
