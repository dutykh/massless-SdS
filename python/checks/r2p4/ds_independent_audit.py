#!/usr/bin/env python3
"""Independent fixed-root cancellation and perturbation audit of selected SdS roots.

No root search is performed. At220decimal digits, recompute local series term
sums independently from the polynomial coefficients, quantify cancellation,
and evaluate the connection determinant at the stored root and alpha+/-1e-9.
Run after ds_extend_continuation.py. This verifies sampled points only, not
interval-certified continuation across intermediate exceptional parameters.
"""
import sys,json
from pathlib import Path
import mpmath as mp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from ds_branch_check import coefficients,geometry,connection
mp.mp.dps=220
r=json.loads((HERE/'ds_extended_results.json').read_text())
print('samples',len(r['path']),'selected',len(r['selected']))
print('minimum sampled integer distance',mp.nstr(min(mp.mpf(p['minimum_integer_index_distance']) for p in r['path']),15))
for entry in r['selected']:
 lam=mp.mpf(entry['Lambda_tilde'])
 if lam not in (mp.mpf('.06'),mp.mpf('.111'),mp.mpf('.1111')):continue
 alpha=mp.mpf(entry['validation'][-1]['alpha']);N=entry['validation'][-1]['N'];rh,rc=geometry(lam);D=rc-rh
 print('LAM',mp.nstr(lam,8),'alpha',mp.nstr(alpha,24),'N',N)
 for t in [mp.mpf('.35') if lam==mp.mpf('.06') else mp.mpf('.45'),mp.mpf('.5') if rh/D>mp.mpf('.5') else mp.mpf('.45')*rh/D]:
  out=[]
  for R,sigma,q in [(rh,1,t),(rc,-1,1-t)]:
   A,B,C,rho=coefficients(lam,R,D,sigma,alpha);cc=[mp.mpf(1)]
   for n in range(1,N+1):
    tot=mp.mpf(0)
    for k in range(1,min(n,max(len(A)-3,len(B)-2,len(C)-1))+1):
     p=rho+n-k;v=(A[k+2]*p*(p-1) if k+2<len(A) else 0)+(B[k+1]*p if k+1<len(B) else 0)+(C[k] if k<len(C) else 0)
     tot+=v*cc[n-k]
    cc.append(-tot/(A[2]*n*(n+2*rho)))
   terms=[cc[k]*q**k for k in range(N+1)];dt=[(k+rho)*cc[k]*q**(k-1) for k in range(N+1)]
   F=mp.fsum(terms);der=mp.fsum(dt)
   condF=mp.fsum(abs(x) for x in terms)/abs(F);condD=mp.fsum(abs(x) for x in dt)/abs(der)
   out.append((mp.log10(condF),mp.log10(condD)))
  w=connection(lam,alpha,N,match_fraction=t)
  wm=connection(lam,alpha-mp.mpf('1e-9'),N,match_fraction=t);wp=connection(lam,alpha+mp.mpf('1e-9'),N,match_fraction=t)
  print('t',mp.nstr(t,9),'lost decimal digits(left F,D;right F,D)',[[mp.nstr(v,8) for v in pair] for pair in out])
  print('W(-delta),W(root),W(+delta)',mp.nstr(wm,9),mp.nstr(w,9),mp.nstr(wp,9))
