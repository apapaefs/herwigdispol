// -*- C++ -*-
// Analytic one-loop quark boxes for gg -> gamma g. GPL version 3.
// Sources: arXiv:2212.14069 (epsilon^0 ancillary coefficients),
// arXiv:2312.16966, sections 2, 4.1 and appendices A,C,F.
#include "GammaJetLoopAmplitudes.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Herwig { namespace GammaJetLoop {
namespace {
const Real pi = std::acos(-1.L);
const Complex I(0.L,1.L);

// Long-double version of the reflection/inversion/series strategy used by
// Herwig::Math::Li2. Local precision is needed for the massive-box cancellations.
Complex li2(Complex z) {
  if (std::abs(z)==0.L) return 0.L;
  if (std::abs(z-Complex(1.L))<1.e-19L) return pi*pi/6.L;
  if (std::abs(z)>1.L)
    return -li2(1.L/z)-pi*pi/6.L-std::log(-z)*std::log(-z)/2.L;
  if (z.real()>0.5L)
    return pi*pi/6.L-std::log(z)*std::log(1.L-z)-li2(1.L-z);
  if (z.real()<-0.3L)
    return -li2(z/(z-1.L))-std::log(1.L-z)*std::log(1.L-z)/2.L;
  Complex result=0.L, power=z;
  for (unsigned n=1;n<100000;++n) {
    Complex term=power/(Real(n)*n);
    result+=term;
    if (std::abs(term)<2.e-20L*std::max(1.L,std::abs(result))) return result;
    power*=z;
  }
  throw std::runtime_error("GammaJetLoop: dilogarithm did not converge");
}

Complex logarithm(Real x) {
  // x -> x+i0: log((beta-1)/(beta+1)), with stable spacelike ratio.
  Complex beta=std::sqrt(1.L-4.L/Complex(x,1.e-30L*std::max(1.L,std::abs(x))));
  return std::log(-4.L/(Complex(x,1.e-30L*std::max(1.L,std::abs(x)))*
                         (1.L+beta)*(1.L+beta)));
}

Complex box(Real x, Real y) {
  // f6^(1,2)/sqrt(x*y*(x*y-4*(x+y))) = m^4 D0/2.
  // Symmetry allows regions I/II only, avoiding a separate region III path.
  if(y>0.L) std::swap(x,y);
  const Real r=(2.L-y+std::sqrt(y*(y-4.L)))/2.L;
  const Real b=4.L*(r-1.L)*(r-1.L)-x*(1.L+r*r);
  const Real discr=b*b-4.L*x*x*r*r;
  Complex root=std::sqrt(Complex(discr,-1.e-30L*std::max(1.L,std::abs(discr))));
  Complex q;
  if(b>=0.L) q=(-b-root)/(2.L*x*r*r);
  else q=2.L*x/(-b+root);
  Complex w=std::sqrt(q), z=r*w;
  if(discr>=0.L) {w=Complex(w.real(),-1.e-30L);z=Complex(z.real(),-1.e-30L);}
  auto L=[](Complex a) {return std::log(1.L-a)-std::log(1.L+a);};
  auto D=[](Complex a) {return li2(a)-li2(-a);};
  Complex f=(std::log(w)-std::log(z))*(L(w)+L(z))+2.L*(D(w)-D(z));
  if(x>0.L || y>0.L) f=-f;
  return f/std::sqrt(Complex(x*y*(x*y-4.L*(x+y))));
}

// Convergent low-energy series of the same scalar integrals, with the
// analytically vanishing constant and linear amplitude terms removed.
// This also avoids losing heavy-quark decoupling through cancellation.
HelicityFunctions lowEnergy(const std::array<Real,3> &x) {
  constexpr unsigned N=18;
  const Real scale=std::max({std::abs(x[0]),std::abs(x[1]),std::abs(x[2])});
  Real c[3]={x[0]/scale,x[1]/scale,x[2]/scale};
  Real B[3][N+2]={},T[3][N+2]={},D[3][N+2]={};
  Real fac[2*N+6];fac[0]=1.L;
  for(unsigned n=1;n<2*N+6;++n) fac[n]=n*fac[n-1];
  for(unsigned i=0;i<3;++i) for(unsigned n=0;n<=N+1;++n) {
    B[i][n]=fac[n]*fac[n]/(2.L*fac[2*n+1])*std::pow(c[i],n);
    if(n) T[i][n]=-fac[n-1]*fac[n-1]/fac[2*n]*std::pow(c[i],n);
    for(unsigned a=0;a<=n;++a)
      D[i][n]+=(n+1.L)*fac[n]*fac[a]*fac[n-a]/(2.L*fac[2*n+3])*
        std::pow(c[i],a)*std::pow(c[(i+1)%3],n-a);
  }
  HelicityFunctions h{};
  for(unsigned n=2;n<=N;++n) {
    Real all=0.L,one=0.L,tSum=0.L,inv=0.L;
    for(unsigned i=0;i<3;++i) {
      all-=4.L*D[i][n];
      one-=4.L*D[i][n]+2.L*c[i]*c[(i+1)%3]/c[(i+2)%3]*D[i][n-1];
      tSum+=T[i][n+1];inv+=1.L/c[i];
    }
    one+=2.L*inv*tSum;
    const Real power=std::pow(scale,n);
    h[0]+=power*all;
    if(n>=3) h[1]+=power*one;
    for(unsigned i=0;i<3;++i) {
      const unsigned j=(i+1)%3,k=(i+2)%3;
      const Real a=c[i],b=c[j],d=c[k];
      Real two=(b-d)/a*(b*B[j][n-1]-4.L*B[j][n])+
        (d-b)/a*(d*B[k][n-1]-4.L*B[k][n])-
        (1.L-2.L*b*d/(a*a))*(T[j][n]+T[k][n])+
        4.L/a*(T[j][n+1]+T[k][n+1])+
        2.L*a*(D[i][n-1]+D[k][n-1])-4.L*(D[i][n]+D[k][n])-
        4.L*D[j][n]-(-2.L*a+8.L*b*d/a)*D[j][n-1]-
        (-b*d+2.L*b*b*d*d/(a*a))*D[j][n-2];
      h[i+2]+=power*two;
    }
  }
  return h;
}
}

HelicityFunctions massless(Real s, Real t, Real u) {
  // H/2 from the finite epsilon^0 ancillary coefficients; G(0,x)=log x,
  // G(1,x)=log(1-x). Analytic continuation implements the physical +i0.
  HelicityFunctions h{{1.L,1.L,0.L,0.L,0.L}};
  Real x[3]={s,t,u};
  for(unsigned i=0;i<3;++i) {
    Real a=x[i],b=x[(i+1)%3],c=x[(i+2)%3];
    Complex l=std::log(Complex(-b,-0.L))-std::log(Complex(-c,-0.L));
    h[i+2]=-1.L-(b-c)/a*l-(b*b+c*c)/(2.L*a*a)*(l*l+pi*pi);
  }
  return h;
}

HelicityFunctions massive(Real s, Real t, Real u, Real mass) {
  if(mass<0.L || !std::isfinite(mass))
    throw std::invalid_argument("GammaJetLoop: invalid quark mass");
  if(mass==0.L) return massless(s,t,u);
  const std::array<Real,3> x{{s/(mass*mass),t/(mass*mass),u/(mass*mass)}};
  if(std::min({std::abs(x[0]),std::abs(x[1]),std::abs(x[2])})>1.e12L)
    return massless(s,t,u);
  if(std::max({std::abs(x[0]),std::abs(x[1]),std::abs(x[2])})<0.5L)
    return lowEnergy(x);
  Complex B[3],T[3],D[3];
  for(unsigned i=0;i<3;++i) {
    const Complex l=logarithm(x[i]);
    B[i]=x[i]==4.L ? Complex(-0.5L) :
      l/(x[i]*std::sqrt(1.L-4.L/Complex(x[i],1.e-30L*std::max(1.L,std::abs(x[i])))));
    T[i]=l*l/2.L;
    D[i]=box(x[i],x[(i+1)%3]);
  }
  HelicityFunctions h{{1.L,1.L,0.L,0.L,0.L}};
  Complex tsum=T[0]+T[1]+T[2];
  h[1]+=2.L*(1.L/x[0]+1.L/x[1]+1.L/x[2])*tsum;
  for(unsigned i=0;i<3;++i) {
    const unsigned j=(i+1)%3,k=(i+2)%3;
    Real a=x[i],b=x[j],c=x[k];
    h[0]-=4.L*D[i];
    h[1]-=2.L*(a*b+2.L*c)/c*D[i];
    h[i+2]=-1.L+(b-4.L)*(b-c)/a*B[j]+(c-4.L)*(c-b)/a*B[k]-
      (1.L-4.L/a-2.L*b*c/(a*a))*(T[j]+T[k])+
      2.L*(a-2.L)*(D[i]+D[k])-
      (4.L-2.L*a-b*c+2.L*b*c*(b*c+4.L*a)/(a*a))*D[j];
  }
  for(auto a:h) if(!std::isfinite(a.real()) || !std::isfinite(a.imag()))
    throw std::runtime_error("GammaJetLoop: non-finite massive amplitude");
  return h;
}

Complex dot(const Vector&a,const Vector&b) {
  return a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
}

Complex contract(const Momenta&p,const std::array<Vector,4>&e,
                 const HelicityFunctions&h) {
  const Complex s=2.L*dot(p[0],p[1]), t=2.L*dot(p[1],p[2]), u=2.L*dot(p[0],p[2]);
  const Complex AS=(h[0]+h[2]+h[3]+h[4]+4.L*h[1])/2.L;
  const Complex C=4.L*h[1]/(s*u);
  // In four dimensions only A1+A2+A3 is physical: choose Ai=AS/3,
  // and determine the corresponding B_i from ACS eq. (2.27).
  Complex result=0.L;
  for(unsigned k=0;k<3;++k) {
    unsigned a=0,b=1,c=2;
    Complex S=s,T=t,U=u,H=h[2];
    if(k==1) {a=2;b=1;c=0;S=t;T=s;U=u;H=h[3];}
    if(k==2) {a=2;b=0;c=1;S=u;T=s;U=t;H=h[4];}
    auto ep=[&](unsigned i,unsigned j){return dot(e[i],p[j]);};
    const Complex g=dot(e[a],e[b])-2.L*ep(b,a)*ep(a,b)/S;
    const Complex A=g*(dot(e[c],e[3])+2.L*(ep(c,a)+ep(c,b))*ep(3,c)/S);
    const Complex B=-g*(ep(c,a)-U/T*ep(c,b))*(ep(3,b)-U/S*ep(3,c))+
      (ep(a,c)-U/S*ep(a,b))*(ep(b,a)-S/T*ep(b,c))*
      (dot(e[c],e[3])+2.L*ep(c,b)*ep(3,a)/T+2.L*ep(c,a)*ep(3,b)/U);
    result+=AS/3.L*A+(h[0]+H+2.L*h[1]-4.L*AS/3.L)/U*B;
  }
  result+=C/(s*t*t*u)*(s*dot(e[0],p[2])-u*dot(e[0],p[1]))*
    (s*dot(e[1],p[2])-t*dot(e[1],p[0]))*
    (u*dot(e[2],p[1])-t*dot(e[2],p[0]))*
    (s*dot(e[3],p[1])-u*dot(e[3],p[2]));
  return result;
}
}}
