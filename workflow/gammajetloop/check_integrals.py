"""Independent high-precision reference for 2312.16966, section 4.1."""
import mpmath as m
m.mp.dps=50

def bubble(x):
 b=m.sqrt(1-4/m.mpc(x,'1e-45'))
 return m.log((b-1)/(b+1))/(x*b)
def triangle(x):
 b=m.sqrt(1-4/m.mpc(x,'1e-45'))
 return m.log((b-1)/(b+1))**2/2

def box(x,y):
 if y>0: x,y=y,x
 r=(2-y+m.sqrt(y*(y-4)))/2
 b=4*(r-1)**2-x*(1+r*r)
 disc=b*b-4*x*x*r*r
 # Region II below threshold: choose positive imaginary w.
 root=m.sqrt(m.mpc(disc,'-1e-45'))
 q=2*x/(-b+root)
 w=m.sqrt(q); z=r*w
 if abs(w.imag)<m.mpf('1e-35'): w=m.mpc(w.real,'-1e-45');z=m.mpc(z.real,'-1e-45')
 L=lambda a:m.log(1-a)-m.log(1+a)
 D=lambda a:m.polylog(2,a)-m.polylog(2,-a)
 f=(m.log(w)-m.log(z))*(L(w)+L(z))+2*(D(w)-D(z))
 if x>0 or y>0:f=-f
 return f/m.sqrt(m.mpc(x*y*(x*y-4*(x+y)),0))

def amp(s,t,u):
 xs=[m.mpf(s),m.mpf(t),m.mpf(u)]
 B=list(map(bubble,xs));T=list(map(triangle,xs))
 D=[box(xs[i],xs[(i+1)%3]) for i in range(3)]
 allp=1-4*sum(D)
 one=1+2*sum(1/x for x in xs)*sum(T)-sum(2*(xs[i]*xs[(i+1)%3]+2*xs[(i+2)%3])/xs[(i+2)%3]*D[i] for i in range(3))
 two=[]
 for i in range(3):
  x,y,z=xs[i],xs[(i+1)%3],xs[(i+2)%3]
  b,c=B[(i+1)%3],B[(i+2)%3]
  f,g=T[(i+1)%3],T[(i+2)%3]
  two.append(-1+(y-4)*(y-z)/x*b+(z-4)*(z-y)/x*c-(1-4/x-2*y*z/x**2)*(f+g)+2*(x-2)*(D[i]+D[(i+2)%3])-(4-2*x-y*z+2*y*z*(y*z+4*x)/x**2)*D[(i+1)%3])
 return [allp,one]+two
if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--quadrature',action='store_true',help='Also integrate two subthreshold boxes (slow)')
 args=parser.parse_args()
 for x,y in [(-.2,-.3),(1.,-.3),(3.9,-1.3),(4.1,-1.3),(10.,-3.),(100.,-30.)]:
  print(x,y,'D/2',m.nstr(box(x,y),15),'amp',*[m.nstr(a,12) for a in amp(x,y,-x-y)] if x>0 else [])
  if args.quadrature and x<2 and y<0:
   ref=m.quad(lambda a:m.quad(lambda b:m.quad(lambda c:1/(1-x*a*c-y*b*(1-a-b-c))**2,[0,1-a-b]),[0,1-a]),[0,1])/2
   print('integral',m.nstr(ref,15))
