"""Independent binary64 diagnostic only; no supplied Maple or Julia execution.
Direct radiative-factor radial QEP collocated on Chebyshev roots.
Limited to the first three positive-real damped roots, N=30,40,50.
"""
import json
from pathlib import Path
import numpy as np
from numpy.polynomial.chebyshev import chebvander,chebder,chebval
from scipy.linalg import eig
OUT=Path(__file__).resolve().parent

def pencil(N,sector,ell,Lambda):
    roots=np.sort(np.roots([4*Lambda/3,0,-1,1]).real)
    x0,xh,xc=roots; zc=xc/xh; z0=x0/xh; d=zc-1
    ah=2*xh*xh/(3-2*xh); ac=2*xc*xc/(3-2*xc)
    y=np.cos(np.pi*(np.arange(N)+.5)/N); z=(zc+1)/2+d*y/2
    a=z-1; b=z-zc; c=z-z0
    f=-4*Lambda*xh*xh/3*a*b*c/z
    fp=f*(1/a+1/b+1/c-1/z)
    H=-ah/a+ac/b
    F2HprimeplusFFprimeH=f*f*((ac-ah)/(a*b)+H*(1/c-1/z))
    if sector=='em': V=f*ell*(ell+1)/z**2
    elif sector=='odd': V=f*(ell*(ell+1)/z**2-3/(xh*z**3))
    elif sector=='even':
        cc=ell*(ell+1)/2-1; R=2*xh*z
        V=8*xh*xh*f*(9+3*cc*cc*R*R+cc*cc*(1+cc)*R**3+3*(3*cc*R-Lambda*R**3))/(R**3*(3+cc*R)**2)
    else:raise ValueError(sector)
    T=chebvander(y,N-1); T1=np.zeros_like(T); T2=np.zeros_like(T)
    for j in range(1,N):
        coeff=np.zeros(N); coeff[j]=1
        T1[:,j]=chebval(y,chebder(coeff))
        T2[:,j]=chebval(y,chebder(coeff,2))
    A0=(4*f*f/d**2)[:,None]*T2+(2*f*fp/d)[:,None]*T1-V[:,None]*T
    A1=1j*((4*f*f*H/d)[:,None]*T1+F2HprimeplusFFprimeH[:,None]*T)
    A2=(4*xh*xh-(f*H)**2)[:,None]*T
    scale=np.max(np.abs(np.concatenate([A0,A1,A2],axis=1)),axis=1)
    A0/=scale[:,None];A1/=scale[:,None];A2/=scale[:,None]
    Z=np.zeros((N,N));I=np.eye(N)
    left=np.block([[Z,I],[-A0,-A1]]);right=np.block([[I,Z],[Z,A2]])
    vals,vecs=eig(left,right)
    ix=np.flatnonzero(np.isfinite(vals)&(vals.real>1e-6)&(vals.imag<0))
    ix=ix[np.argsort(-vals[ix].imag)][:3]
    ans=[]
    for k in ix:
        w=vals[k];v=vecs[:N,k]
        residual=np.linalg.norm((A0+w*A1+w*w*A2)@v)/(np.linalg.norm(v)*(np.linalg.norm(A0)+abs(w)*np.linalg.norm(A1)+abs(w*w)*np.linalg.norm(A2)))
        ans.append({'re':w.real,'im':w.imag,'relative_pencil_residual':float(residual)})
    return ans

cases=[('em',2,.11,[.04614-.00962j,.04612-.02886j,.04608-.04810j]),('odd',3,.11,[.06009-.00961j,.06007-.02886j,.06005-.04809j]),('odd',2,.09,[.1626-.0414j,.1608-.1242j,.1570-.2071j]),('even',2,.09,[.1626-.0414j,.1608-.1242j,.1570-.2071j])]
all_results=[]
for sector,ell,lam,reported in cases:
    item={'sector':sector,'ell':ell,'Lambda_tilde':lam,'orders':{},'manuscript_SM':[{'re':w.real,'im':w.imag} for w in reported]}
    for N in [30,40,50]:item['orders'][str(N)]=pencil(N,sector,ell,lam)
    for n in range(3):
        r=item['orders']['50'][n];w=complex(r['re'],r['im']);old=item['orders']['40'][n]
        r['difference_N40']=abs(w-complex(old['re'],old['im']))
        r['difference_printed_SM']=abs(w-reported[n])
    all_results.append(item)
    print(sector,ell,lam)
    for N,rows in item['orders'].items():print(N,[(round(x['re'],12),round(x['im'],12)) for x in rows])
(OUT/'low_mode_diagnostic.json').write_text(json.dumps({'precision':'IEEE binary64','description':__doc__,'results':all_results},indent=2))
