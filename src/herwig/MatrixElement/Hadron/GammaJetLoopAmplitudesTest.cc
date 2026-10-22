#include "GammaJetLoopAmplitudes.h"
#include <cassert>
#include <iostream>
#include <cmath>
using namespace Herwig::GammaJetLoop;

void close(Complex a,Complex b,Real tolerance=2.e-9L) {
  if(std::abs(a-b)>tolerance*std::max(1.L,std::abs(b))) {
    std::cerr<<"Mismatch "<<a<<" versus "<<b<<"\n";std::abort();
  }
}
Vector transform(Vector a) {
  const Real c=std::cos(.73L),s=std::sin(.73L),beta=.37L,g=1.L/std::sqrt(1-beta*beta);
  Vector r=a;r[1]=c*a[1]-s*a[2];r[2]=s*a[1]+c*a[2];
  Vector b=r;b[0]=g*(r[0]-beta*r[1]);b[1]=g*(r[1]-beta*r[0]);return b;
}
int main() {
  // Independent epsilon^0 BCG ancillary values at x=-s13/s12=0.3,
  // divided by two. The ACS ordering uses t=s23, u=s13.
  auto ml=massless(1.L,-.7L,-.3L);
  close(ml[0],1.L);close(ml[1],1.L);
  close(ml[2],-3.731461094783885178124L);
  close(ml[3],Complex(-.3763056692891186955926L,2.579506312189497750307L));
  close(ml[4],Complex(-.03191616931501568166703L,.7485949803271134380895L));
  auto h=massive(10.L,-3.L,-7.L,1.L);
  const HelicityFunctions ref{{{.814357117838985L,-.676871940346723L},
    {.000445781982268L,-.0893950517801L},{-.956745775787L,2.70748776139L},
    {.1035196588L,.193361009454L},{.461575647268L,.99421021659L}}};
  for(unsigned i=0;i<5;++i) close(h[i],ref[i]);
  h=massive(.001L,-.0003L,-.0007L,1.L);
  close(h[0]/1.e-8L,-1.3167500330254L);
  close(h[1]/1.e-14L,-8.3333334330857L);
  close(h[2]/1.e-8L,3.055714338962L);
  for(Real s:{.499999L,.500001L,3.999999L,4.L,4.000001L,10.L,100.L}) {
    auto f=massive(s,-s*.3L,-s*.7L,1.L);
    for(auto z:f) assert(std::isfinite(z.real())&&std::isfinite(z.imag()));
  }
  for(auto mode:{0,1}) for(Real c:{-.9L,-.4L,0.L,.4L,.9L}) {
    Real st=std::sqrt(1.L-c*c),r=std::sqrt(.5L);
    Momenta p{{{{1.L,0.L,0.L,1.L}},{{1.L,0.L,0.L,-1.L}},
               {{-1.L,-st,0.L,-c}},{{-1.L,st,0.L,c}}}};
    auto f=mode?massive(4.L,-2.L*(1+c),-2.L*(1-c),1.2L):
                massless(4.L,-2.L*(1+c),-2.L*(1-c));
    Real sum=0.L;
    for(unsigned n=0;n<16;++n) {
      std::array<Vector,4> e;
      Real dirs[4][3]={{0,0,1},{0,0,-1},{st,0,c},{-st,0,-c}};
      for(unsigned i=0;i<4;++i) {
        Real sign=(n>>i&1)?1.L:-1.L;
        e[i]={{0.L,r*dirs[i][2],Complex(0.L,r*sign),-r*dirs[i][0]}};
      }
      Complex a=contract(p,e,f);sum+=std::norm(a);
      for(unsigned i=0;i<4;++i) {
        auto gauge=e;gauge[i]=p[i];close(contract(p,gauge,f),0.L);
        gauge=e;for(unsigned j=0;j<4;++j) gauge[i][j]+=Complex(.17L,.3L)*p[i][j];
        close(contract(p,gauge,f),a);
      }
      auto pp=p,ee=e;
      for(unsigned i=0;i<4;++i){pp[i]=transform(p[i]);ee[i]=transform(e[i]);}
      close(contract(pp,ee,f),a);
      for(unsigned i=0;i<4;++i) for(unsigned j=1;j<4;++j){pp[i][j]*=-1.L;ee[i][j]*=-1.L;}
      close(contract(pp,ee,f),a);
      pp=p;ee=e;std::swap(pp[0],pp[1]);std::swap(ee[0],ee[1]);
      auto swapped=f;std::swap(swapped[3],swapped[4]);
      close(contract(pp,ee,swapped),a);
      auto phased=e;Complex phase=std::exp(Complex(0.L,.31L));
      for(auto &z:phased[2]) z*=phase;
      close(contract(p,phased,f),a*phase);
    }
    close(sum,2.L*std::norm(f[0])+8.L*std::norm(f[1])+
                2.L*(std::norm(f[2])+std::norm(f[3])+std::norm(f[4])));
  }
  std::cout<<"PASS: massive reference, low-energy series, threshold finiteness, "
              "Ward identities, gauge shifts, rotations/boosts, parity, Bose "
              "exchange, rephasing and helicity sum\n";
}
