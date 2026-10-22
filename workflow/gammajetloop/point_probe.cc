#include "GammaJetLoopAmplitudes.h"
#include <iostream>
#include <iomanip>
using namespace Herwig::GammaJetLoop;
// Input: mc mb mt (-1 omits that flavour), four positive-energy momenta, then sixteen
// groups of four complex polarization vectors in (E,x,y,z) order.
int main() {
  Real masses[3];for(auto &m:masses) std::cin>>m;
  Momenta p;
  for(unsigned i=0;i<4;++i) for(auto &v:p[i]) {Real x;std::cin>>x;v=(i<2?x:-x);}
  Real s=real(2.L*dot(p[0],p[1])),t=real(2.L*dot(p[1],p[2])),u=real(2.L*dot(p[0],p[2]));
  HelicityFunctions h{};
  for(unsigned i=0;i<3;++i) {
    if(masses[i]<0.L) continue;
    auto f=massive(s,t,u,masses[i]);
    for(unsigned j=0;j<5;++j) h[j]+=(i==1?-1.L:2.L)/3.L*f[j];
  }
  std::cout<<std::setprecision(18);
  for(unsigned n=0;n<16;++n) {
    std::array<Vector,4> e;
    for(auto &v:e) for(auto &z:v) {Real a,b;std::cin>>a>>b;z=Complex(a,b);}
    Complex a=contract(p,e,h);
    std::cout<<a.real()<<" "<<a.imag()<<"\n";
  }
}
