#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "Herwig/MatrixElement/HardVertex.h"
#include "ThePEG/Helicity/VectorSpinInfo.h"
#include <cassert>
#include <iostream>
using namespace Herwig;
using namespace ThePEG;
using namespace ThePEG::Helicity;
int main() {
  for(bool flat:{true,false}) {
    ProductionMatrixElement me(PDT::Spin1,PDT::Spin1,PDT::Spin1,PDT::Spin1);
    for(unsigned a:{0u,2u}) for(unsigned b:{0u,2u})
      for(unsigned c:{0u,2u}) for(unsigned d:{0u,2u})
        me(a,b,c,d)=flat?Complex(1.):Complex(1.+a+2.*b+3.*c+5.*d,.7+a*b-.2*c*d);
    RhoDMatrix unpol(PDT::Spin1,false);unpol(0,0)=unpol(2,2)=.5;
    double average=0.;
    for(unsigned a:{0u,2u}) for(unsigned b:{0u,2u}) {
      RhoDMatrix first(PDT::Spin1,false),second(PDT::Spin1,false);
      first(a,a)=second(b,b)=1.;average+=me.average(first,second)/4.;
    }
    assert(abs(average-me.average(unpol,unpol))<1.e-10);
    if(flat) assert(abs(average-4.)<1.e-12);
    for(bool coherent:{false,true}) {
      RhoDMatrix first=unpol,second=unpol;
      if(coherent){first(0,2)=Complex(.2,.1);first(2,0)=conj(first(0,2));
                   second(0,0)=.7;second(2,2)=.3;second(0,2)=.15;second(2,0)=.15;}
      HardVertexPtr v=new_ptr(HardVertex());v->ME(me);
      vector<Ptr<SpinInfo>::pointer> spin;
      for(unsigned i=0;i<4;++i) {
        spin.push_back(new_ptr(VectorSpinInfo(Lorentz5Momentum(),i>=2)));
        spin.back()->productionVertex(v);
        spin.back()->DMatrix(unpol);
      }
      spin[0]->rhoMatrix(first);spin[1]->rhoMatrix(second);
      auto rho=v->getRhoMatrix(0,false);
      Complex direct[3][3]={};
      for(unsigned c:{0u,2u}) for(unsigned cp:{0u,2u})
        for(unsigned a:{0u,2u}) for(unsigned ap:{0u,2u})
          for(unsigned b:{0u,2u}) for(unsigned bp:{0u,2u})
            for(unsigned d:{0u,2u}) direct[c][cp]+=
              me(a,b,c,d)*conj(me(ap,bp,cp,d))*first(a,ap)*second(b,bp)*.5;
      double tr=real(direct[0][0]+direct[2][2]);
      for(unsigned c=0;c<3;++c) for(unsigned cp=0;cp<3;++cp)
        assert(abs(rho(c,cp)-direct[c][cp]/tr)<1.e-11);
      assert(abs(rho(0,0)+rho(2,2)-1.)<1.e-12);
      assert(real(rho(0,0)*rho(2,2))-norm(rho(0,2))>-1.e-12);
      assert(abs(rho(0,2)-conj(rho(2,0)))<1.e-12);
    }
  }
  std::cout<<"PASS: flat normalization, four-helicity closure and HardVertex "
              "density matrices including off-diagonal entries\n";
}
