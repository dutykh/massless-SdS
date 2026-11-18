#!/usr/bin/env python3
"""Targeted high-precision eigenpair refinement of the selected imaginary mode.

Every collocation node, horizon, metric coefficient and polynomial derivative
is rebuilt in mpmath. No binary64 matrix is promoted to high precision. Initial
eigenpairs come solely from the preceding binary64 spectral continuation. The
independent Frobenius frequency and eigenfunction are used only after Newton
convergence for error assessment. The pure-imaginary QEP is real in alpha=-Im(w).
"""
from pathlib import Path
import argparse,csv,json,time
import mpmath as mp

HERE=Path(__file__).resolve().parent
LAM='0.06';YSTAR='-0.3'

def bg():
 lam=mp.mpf(LAM);theta=mp.acos(3*mp.sqrt(lam))/3
 a=mp.cos(theta)/(2*mp.sqrt(lam));d=mp.sqrt(3)*mp.sin(theta)/mp.sqrt(lam)
 return lam,a,d,a-d/2,a+d/2,-2*a

def coeffs(y):
 lam,a,d,h,c,x0=bg();x=a+d*y/2
 F=(x-x0)/(3*x);Fy=x0*d/(6*x*x)
 g=(1-y*y)*F;gy=-2*y*F+(1-y*y)*Fy
 C=3*a*a*d/((2*h+c)*(h+2*c))
 G=y+C*(1-y*y)/x;Gy=1+C*(-2*y/x-d*(1-y*y)/(2*x*x))
 R=(1-2*y*C/x-(C/x)**2*(1-y*y))/F
 V=(6/x**2+2*lam*d*gy/x)/(4*lam)
 return g,gy,G,Gy,R,V

def basis(y,N):
 T=[mp.mpf(1),y];D=[mp.mpf(0),mp.mpf(1)];DD=[mp.mpf(0),mp.mpf(0)]
 for k in range(1,N-1):
  T.append(2*y*T[-1]-T[-2])
  D.append(2*T[-2]+2*y*D[-1]-D[-2])
  DD.append(4*D[-2]+2*y*DD[-1]-DD[-2])
 return T[:N],D[:N],DD[:N]

def matrices(N):
 A0=mp.matrix(N);B1=mp.matrix(N);A2=mp.matrix(N)
 for j in range(N):
  y=mp.cos(mp.pi*(mp.mpf(j)+mp.mpf('.5'))/N)
  T,T1,T2=basis(y,N);g,gy,G,Gy,R,V=coeffs(y)
  r0=[g*T2[k]+gy*T1[k]-V*T[k] for k in range(N)]
  r1=[2*G*T1[k]+Gy*T[k] for k in range(N)]
  r2=[R*T[k] for k in range(N)]
  scale=max(abs(z) for row in [r0,r1,r2] for z in row)
  for k in range(N):A0[j,k]=r0[k]/scale;B1[j,k]=r1[k]/scale;A2[j,k]=r2[k]/scale
 return A0,B1,A2,mp.matrix(basis(mp.mpf(YSTAR),N)[0])

def chebval(y,c):
 b1=mp.mpf(0);b2=mp.mpf(0)
 for a in reversed(c[1:]):b0=2*y*b1-b2+a;b2=b1;b1=b0
 return c[0]+y*b1-b2

def chebder(c):
 n=len(c)-1
 if n<1:return [mp.mpf(0)]
 b=[mp.mpf(0)]*(n+2)
 for j in range(n-1,-1,-1):b[j]=b[j+2]+2*(j+1)*c[j+1]
 b[0]/=2
 return b[:n]

def diagnostics(alpha,c,A0,B1,A2,reference):
 lam,a,d,h,cx,x0=bg();sigma=lam*d
 p1=chebder(c);p2=chebder(p1)
 max_raw=mp.mpf(0);max_rawscale=mp.mpf(0);max_reduced=mp.mpf(0);max_redscale=mp.mpf(0)
 endpoint=[]
 for j in range(1025):
  y=mp.cos(mp.pi*j/1024);p=chebval(y,c);dp=chebval(y,p1);ddp=chebval(y,p2)
  g,gy,G,Gy,R,V=coeffs(y)
  terms=[g*ddp,(gy+2*alpha*G)*dp,(-V+alpha*Gy-alpha*alpha*R)*p]
  max_reduced=max(max_reduced,abs(sum(terms)));max_redscale=max(max_redscale,sum(abs(t) for t in terms))
  if j in [0,1024]:endpoint.append(abs(sum(terms)))
  else:
   x=a+d*y/2;f=lam*d*d*(1-y*y)*(x-x0)/(3*x);fx=1/x**2-8*lam*x/3
   gg=f/(lam*d*d);ggp=fx/(2*lam*d)
   ph=3*h/(2*(h-x0));pc=3*cx/(2*(cx-x0))
   H=-ph/(1+y)+pc/(1-y);Hp=ph/(1+y)**2+pc/(1-y)**2
   U=f*(6/x**2+fx/x)/(4*sigma*sigma)
   raw=[gg*gg*ddp,2*alpha*gg*gg*H*dp,alpha*gg*gg*Hp*p,
        alpha*alpha*gg*gg*H*H*p,gg*ggp*dp,alpha*gg*ggp*H*p,-alpha*alpha*p,-U*p]
   max_raw=max(max_raw,abs(sum(raw)));max_rawscale=max(max_rawscale,sum(abs(t) for t in raw))
 C=mp.matrix(c);Q=A0+alpha*B1-alpha*alpha*A2
 p_res=mp.norm(Q*C)/(mp.norm(C)*(mp.norm(A0)+abs(alpha)*mp.norm(B1)+alpha*alpha*mp.norm(A2)))
 phi=[];diff=[];refnorm=mp.mpf(0)
 for point in reference['points']:
  y=mp.mpf(point['y']);rr=mp.mpf(point['phi']);value=chebval(y,c)
  phi.append(value);diff.append(abs(value-rr));refnorm=max(refnorm,abs(rr))
 alpha_phys_ref=mp.mpf(reference['alpha']);wref=alpha_phys_ref/sigma
 return {
  'varpi_im':str(-alpha),'Omega_im':str(-alpha*sigma),
  'error_to_independent_varpi_reference':str(abs(alpha-wref)),
  'error_to_independent_Omega_reference':str(abs(alpha*sigma-alpha_phys_ref)),
  'relative_pencil_residual':str(p_res),
  'independent_original_ODE_residual':str(max_raw/max_rawscale),
  'relative_reduced_ODE_residual_including_endpoints':str(max_reduced/max_redscale),
  'relative_reduced_endpoint_residual':str(max(endpoint)/max_redscale),
  'last_eight_coefficient_norm_fraction':str(mp.sqrt(sum(t*t for t in c[-8:]))/mp.norm(C)),
  'eigenfunction_error_to_independent_Frobenius_Linf':str(max(diff)/refnorm),
  'coefficients':[str(t) for t in c],'phi':[str(t) for t in phi],
 }

def solve(N,dps,seed,reference):
 started=time.time();mp.mp.dps=dps
 A0,B1,A2,tstar=matrices(N)
 alpha=mp.mpf(str(-seed['varpi_im']))
 c=mp.matrix([mp.mpf(str(x)) for x in seed['coeff_re']])
 if max(abs(x) for x in seed['coeff_im'])>1e-9*max(abs(x) for x in seed['coeff_re']):
  raise RuntimeError('Seed is not a nearly real normalized imaginary-mode eigenvector')
 c/=mp.fdot(tstar,c)
 iterations=[];tol=mp.mpf(10)**(-min(50,dps-15))
 for iteration in range(8):
  Q=A0+alpha*B1-alpha*alpha*A2
  r=Q*c;normalization=mp.fdot(tstar,c)-1
  J=mp.matrix(N+1);rhs=mp.matrix(N+1,1)
  qa=(B1-2*alpha*A2)*c
  for j in range(N):
   for k in range(N):J[j,k]=Q[j,k]
   J[j,N]=qa[j];J[N,j]=tstar[j];rhs[j]=-r[j]
  rhs[N]=-normalization
  step=mp.lu_solve(J,rhs)
  for j in range(N):c[j]+=step[j]
  alpha+=step[N]
  relstep=mp.norm(step[:N])/mp.norm(c)
  record={'iteration':iteration,'frequency_step':str(step[N]),'relative_vector_step':str(relstep),'elapsed_seconds':time.time()-started}
  iterations.append(record)
  print('MP',N,'iteration',iteration,'dw',mp.nstr(abs(step[N]),5),'da',mp.nstr(relstep,5),'elapsed',round(time.time()-started,2),flush=True)
  if abs(step[N])<tol and relstep<tol:break
 else:raise RuntimeError('Newton failed to converge')
 if abs(alpha-mp.mpf(str(-seed['varpi_im'])))>mp.mpf('.05'):
  raise RuntimeError('Newton escaped original spectral candidate neighborhood')
 result={'N':N,'dps':dps,'iterations':iterations,'elapsed_seconds':time.time()-started,
         'seed':'binary64 spectral candidate, no independent reference in solve',
         **diagnostics(alpha,list(c),A0,B1,A2,reference)}
 print('MP COMPLETE',N,'reference error',mp.nstr(mp.mpf(result['error_to_independent_varpi_reference']),6),'field error',mp.nstr(mp.mpf(result['eigenfunction_error_to_independent_Frobenius_Linf']),6),flush=True)
 return result

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--orders',type=int,nargs='+',default=[24,48,96,128]);parser.add_argument('--dps',type=int,default=70);parser.add_argument('--diagnostics-only',action='store_true')
 args=parser.parse_args()
 source=json.loads((HERE/'broad_resolution_results.json').read_text())
 binary=next(c for c in source['cases'] if c['definition']['name']=='scalar_imaginary')['results']
 reference=json.loads((HERE/'ds_reference_eigenfunction.json').read_text())
 path=HERE/'imaginary_multiprecision_results.json'
 report={'scope':__doc__,'method':'Fully rebuilt real imaginary-frequency QEP; bordered Newton for eigenvalue and normalized eigenvector, with Phi(-.3)=1. No reference root in selection or Newton.','mpmath':mp.__version__,'results':[]}
 if args.diagnostics_only:
  report=json.loads(path.read_text())
  for result in report['results']:
   mp.mp.dps=result['dps'];N=result['N']
   A0,B1,A2,_=matrices(N)
   result.update(diagnostics(-mp.mpf(result['varpi_im']),[mp.mpf(t) for t in result['coefficients']],A0,B1,A2,reference))
   print('Refreshed diagnostics',N,flush=True)
  report['offgrid_test_points']=1025
  path.write_text(json.dumps(report,indent=2)+'\n')
 elif path.exists():
  old=json.loads(path.read_text());report['results']=[r for r in old['results'] if r['N'] not in args.orders]
 for N in ([] if args.diagnostics_only else args.orders):
  seed=next(r for r in binary if r['N']==N)
  report['results'].append(solve(N,args.dps,seed,reference))
  report['results'].sort(key=lambda r:r['N'])
  report['offgrid_test_points']=1025
  path.write_text(json.dumps(report,indent=2)+'\n')
 compact=[{k:v for k,v in r.items() if k not in ['coefficients','phi','iterations']} for r in report['results']]
 with (HERE/'imaginary_multiprecision_results.csv').open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(compact[0]));writer.writeheader();writer.writerows(compact)

if __name__=='__main__':main()
