#!/usr/bin/env python3
"""New bounded, independent axial/polar SdS spectral and radiation check.

Uses the dimensional Regge--Wheeler and Zerilli potentials directly; no
Darboux transform, exported Maple matrices, or manuscript frequency values
are used to construct/initialize either parity. Each parity is discretized
separately. A separate original-equation Frobenius shooting solver checks
both parities, initialized by Nariai frequencies and its own continuation.

Scope: ell=2, two positive-real low modes, Lambda_tilde=.06,.09,.11 only.
No claim about purely imaginary roots, all overtones, or completeness.
Run: OPENBLAS_NUM_THREADS=1 python parity_check.py
"""
from pathlib import Path
from functools import lru_cache
import csv,itertools,json,platform,time
import numpy as np
import scipy
from numpy.polynomial.chebyshev import chebvander,chebval,chebder
from scipy.linalg import eig
from scipy.integrate import solve_ivp
from scipy.optimize import root
from threadpoolctl import threadpool_limits

OUT=Path(__file__).resolve().parent
LAMBDAS=(.06,.09,.11)
ORDERS=(40,60,80,100)
PARITIES=('odd','even')
ELL=2

def geometry(lam):
    theta=np.arccos(3*np.sqrt(lam))/3
    a=np.cos(theta)/(2*np.sqrt(lam))
    d=np.sqrt(3)*np.sin(theta)/np.sqrt(lam)
    xh,xc=a-d/2,a+d/2
    return dict(lam=lam,a=a,d=d,xh=xh,xc=xc,x0=-2*a,
                rh=2*xh,rc=2*xc,R=2*d,S=lam*d)

def potential_over_f(r,lam,parity,ell=ELL):
    """Physical V/f at M=1, separately specified for each parity."""
    if parity=='odd':return ell*(ell+1)/r**2-6/r**3
    if parity!='even':raise ValueError(parity)
    c=ell*(ell+1)/2-1
    return 2*(9+3*c*c*r*r+c*c*(1+c)*r**3+3*(3*c*r-lam*r**3))/(r**3*(3+c*r)**2)

def coefficients(y,lam,parity):
    """Exact regular y equation: g Phi''+(gy+2iwG)Phi'
       +(-Vbar+iwGy+w²R)Phi=0, w=Omega/[lam*(xc-xh)]."""
    b=geometry(lam);a,d,xh,xc,x0=(b[k] for k in ('a','d','xh','xc','x0'))
    x=a+d*y/2;F=(x-x0)/(3*x);Fy=x0*d/(6*x*x)
    g=(1-y*y)*F;gy=-2*y*F+(1-y*y)*Fy
    C=3*a*a*d/((2*xh+xc)*(xh+2*xc))
    G=y+C*(1-y*y)/x
    Gy=1+C*(-2*y/x-d*(1-y*y)/(2*x*x))
    R=(1-2*y*C/x-(C/x)**2*(1-y*y))/F
    # Since M=1 and r=2x, U=4*(V/f), hence Vbar=U/(4lam).
    Vbar=potential_over_f(2*x,lam,parity)/lam
    return g,gy,G,Gy,R,Vbar

@lru_cache(None)
def basis(N):
    y=np.cos(np.pi*(np.arange(N)+.5)/N)
    T=chebvander(y,N-1);T1=np.zeros_like(T);T2=np.zeros_like(T)
    for j in range(1,N):
        c=np.zeros(N);c[j]=1
        T1[:,j]=chebval(y,chebder(c));T2[:,j]=chebval(y,chebder(c,2))
    return y,T,T1,T2

def pencil(N,lam,parity):
    y,T,T1,T2=basis(N);g,gy,G,Gy,R,V=coefficients(y,lam,parity)
    A0=g[:,None]*T2+gy[:,None]*T1-V[:,None]*T
    A1=1j*(2*G[:,None]*T1+Gy[:,None]*T);A2=R[:,None]*T
    row=np.max(np.abs(np.hstack((A0,A1,A2))),axis=1)
    A0/=row[:,None];A1/=row[:,None];A2/=row[:,None]
    I,Z=np.eye(N),np.zeros((N,N))
    vals,vec=eig(np.block([[Z,I],[-A0,-A1]]),np.block([[I,Z],[Z,A2]]))
    # Same explicit selection rectangle at every parity/resolution/background.
    keep=np.isfinite(vals)&(vals.real>.5)&(vals.real<4)&(vals.imag> -3.2)&(vals.imag<-.05)
    ix=np.flatnonzero(keep);ix=ix[np.argsort(-vals[ix].imag)]
    if len(ix)<2:raise RuntimeError('Fewer than two roots in window')
    return vals[ix[:2]],vec[:N,ix[:2]],(A0,A1,A2),vals[ix]

def diagnostics(w,c,matrices,lam,parity):
    c=c/np.max(np.abs(c));A0,A1,A2=matrices
    pencil=np.linalg.norm((A0+w*A1+w*w*A2)@c)/(np.linalg.norm(c)*(np.linalg.norm(A0)+abs(w)*np.linalg.norm(A1)+abs(w*w)*np.linalg.norm(A2)))
    yy=np.cos(np.pi*np.arange(8*len(c)+1)/(8*len(c)))
    p,p1,p2=chebval(yy,c),chebval(yy,chebder(c)),chebval(yy,chebder(c,2))
    g,gy,G,Gy,R,V=coefficients(yy,lam,parity)
    terms=np.array([g*p2,(gy+2j*w*G)*p1,(-V+1j*w*Gy+w*w*R)*p])
    resid=np.abs(terms.sum(axis=0));scale=np.max(np.sum(np.abs(terms),axis=0))
    return dict(pencil_residual=float(pencil),offgrid_residual=float(np.max(resid)/scale),
                endpoint_residual=float(np.max(resid[[0,-1]])/scale),tail_fraction=float(np.linalg.norm(c[-8:])/np.linalg.norm(c)))

def series_data(bg,parity,side,order):
    """Taylor coefficients of g²/q², gg'/q, V in original q equation.

    g=f/[Lambda*R²/2], V=2*g*(Vphysical/f)/Lambda.
    The polar rational potential is expanded directly, independently of odd.
    """
    lam,rh,rc,R=(bg[k] for k in ('lam','rh','rc','R'))
    a,b=(rh,R) if side=='left' else (rc,-R)
    N=order+4;k=np.arange(N)
    inv1=(-b/a)**k/a;inv2=(k+1)*(-b/a)**k/a**2
    inv3=(k+1)*(k+2)/2*(-b/a)**k/a**3
    ratio=(rh+rc)*inv1;ratio[0]+=1
    gs=np.zeros(N);gs[1:]+=2/3*ratio[:-1];gs[2:]-=2/3*ratio[:-2]
    gp=np.arange(1,N)*gs[1:]
    if parity=='odd':B=ELL*(ELL+1)*inv2-6*inv3
    elif parity=='even':
        c=ELL*(ELL+1)/2-1;r=np.array([a,b])
        r2=np.convolve(r,r);r3=np.convolve(r2,r)
        numerator=np.zeros(4);numerator[0]=9
        numerator[:2]+=9*c*r;numerator[:3]+=3*c*c*r2
        numerator+=(c*c*(1+c)-3*lam)*r3
        invc2=(k+1)*(-c*b/(3+c*a))**k/(3+c*a)**2
        B=2*np.convolve(numerator,np.convolve(inv3,invc2))[:N]
    else:raise ValueError(parity)
    Vs=2/lam*np.convolve(gs,B)[:N]
    P=np.convolve(gs[1:],gs[1:])[:order+1]
    Q=np.convolve(gs[1:],gp)[:order+1]
    return gs[1],P,Q,Vs[:order+1]

class Matcher:
    def __init__(self,lam,parity,order=24,eps=.04,mid=.43,rtol=2e-12,atol=2e-13):
        self.bg=geometry(lam);self.parity=parity
        self.order,self.eps,self.mid,self.rtol,self.atol=order,eps,mid,rtol,atol
        self.series={s:series_data(self.bg,parity,s,order) for s in ('left','right')}
    def coefficients(self,q,side):
        bg=self.bg;rh,rc,R,lam=(bg[k] for k in ('rh','rc','R','lam'))
        a,b=(rh,R) if side=='left' else (rc,-R)
        r=a+b*q;v=1+(rh+rc)/r;vp=-(rh+rc)*b/r**2
        g=2/3*q*(1-q)*v;gp=2/3*((1-2*q)*v+q*(1-q)*vp)
        V=2*g/lam*potential_over_f(r,lam,self.parity)
        return g,gp,V
    def local_data(self,w,side):
        slope,P,Q,V=self.series[side];p=-1j*w/slope
        a=np.zeros(self.order+1,dtype=complex);a[0]=1
        for n in range(1,self.order+1):
            total=0j
            for k in range(1,n+1):
                h=p+n-k;total+=(P[k]*h*(h-1)+Q[k]*h-V[k])*a[n-k]
            a[n]=-total/(slope*slope*n*(n+2*p))
        q=self.eps;F=np.polynomial.polynomial.polyval(q,a)
        Fq=np.polynomial.polynomial.polyval(q,np.arange(1,len(a))*a[1:])
        return np.array([q**p*F,q**(p-1)*(p*F+q*Fq)]),a
    def propagate(self,w,side):
        initial,_=self.local_data(w,side);end=self.mid if side=='left' else 1-self.mid
        def rhs(q,Y):
            g,gp,V=self.coefficients(q,side)
            return np.array([Y[1],-gp/g*Y[1]-(w*w-V)/(g*g)*Y[0]])
        sol=solve_ivp(rhs,(self.eps,end),initial,method='DOP853',rtol=self.rtol,atol=self.atol)
        if not sol.success:raise RuntimeError(sol.message)
        return sol.y[:,-1]
    def determinant(self,w):
        L,R=self.propagate(w,'left'),self.propagate(w,'right')
        # R[1] is d/dq=-d/dt, q=1-t, t=(r-rh)/R.
        return -(L[0]*R[1]+L[1]*R[0])/(np.linalg.norm(L)*np.linalg.norm(R))
    def solve(self,seed):
        def fun(v):
            w=self.determinant(complex(*v));return [w.real,w.imag]
        sol=root(fun,[seed.real,seed.imag],tol=1e-10)
        w=complex(*sol.x);res=abs(self.determinant(w))
        if not sol.success or res>1e-8:raise RuntimeError((self.parity,self.bg['lam'],sol.message,res))
        return w,res

def endpoint_check(matcher,w,side):
    _,a=matcher.local_data(w,side);q=matcher.eps
    slope=matcher.series[side][0];p=-1j*w/slope
    F=np.polynomial.polynomial.polyval(q,a)
    F1=np.polynomial.polynomial.polyval(q,np.arange(1,len(a))*a[1:])
    F2=np.polynomial.polynomial.polyval(q,np.arange(2,len(a))*np.arange(1,len(a)-1)*a[2:])
    # Drop the nonzero common q^p from the original ODE residual.
    g,gp,V=matcher.coefficients(q,side)
    terms=np.array([g*g*(p*(p-1)*F/q**2+2*p*F1/q+F2),g*gp*(p*F/q+F1),(w*w-V)*F])
    return float(abs(terms.sum())/np.abs(terms).sum())

def cdict(z):return {'real':float(z.real),'imag':float(z.imag)}

def main():
    start=time.time();spectral=[];jost=[]
    with threadpool_limits(limits=1):
      for lam,parity,N in itertools.product(LAMBDAS,PARITIES,ORDERS):
        vals,vecs,mats,allvals=pencil(N,lam,parity)
        for n,w in enumerate(vals):
            other=np.delete(allvals,n)
            spectral.append(dict(Lambda_tilde=lam,parity=parity,N=N,n=n,
                varpi=cdict(w),Omega=cdict(geometry(lam)['S']*w),
                nearest_other_root=float(np.min(np.abs(w-other))),
                **diagnostics(w,vecs[:,n],mats,lam,parity)))
      for parity,n in itertools.product(PARITIES,(0,1)):
        seed=complex(np.sqrt(15)/2,-n-.5)
        # Each parity independently initialized at .11 and followed downward.
        for lam in reversed(LAMBDAS):
            coarse=Matcher(lam,parity);w,res=coarse.solve(seed)
            fine=Matcher(lam,parity,32,.025,.57,8e-14,8e-15);wf,resf=fine.solve(w)
            third=Matcher(lam,parity,36,.06,.37,8e-14,8e-15)
            rest=abs(third.determinant(wf));scale=geometry(lam)['S']
            record=dict(Lambda_tilde=lam,parity=parity,n=n,varpi=cdict(wf),Omega=cdict(scale*wf),
                coarse_to_refined_varpi_difference=float(abs(wf-w)),
                coarse_to_refined_Omega_difference=float(scale*abs(wf-w)),
                refined_matching_residual=float(resf),third_matching_residual=float(rest),
                endpoint_series_residual=max(endpoint_check(fine,wf,s) for s in ('left','right')))
            jost.append(record);seed=wf
            print('Jost',parity,lam,n,record['Omega'],rest,flush=True)
    rows=[]
    for lam,n in itertools.product(LAMBDAS,(0,1)):
        item={'Lambda_tilde':lam,'ell':ELL,'n':n}
        last={};jw={}
        for parity in PARITIES:
            ss=[r for r in spectral if r['Lambda_tilde']==lam and r['n']==n and r['parity']==parity]
            values=[complex(r['Omega']['real'],r['Omega']['imag']) for r in ss]
            last[parity]=values[-1]
            item[parity+'_Omega_real']=values[-1].real;item[parity+'_Omega_imag']=values[-1].imag
            item[parity+'_resolution_spread']=float(max(abs(a-b) for a,b in itertools.combinations(values,2)))
            jj=next(r for r in jost if r['Lambda_tilde']==lam and r['n']==n and r['parity']==parity)
            jw[parity]=complex(jj['Omega']['real'],jj['Omega']['imag'])
            item[parity+'_spectral_Jost_difference']=float(abs(values[-1]-jw[parity]))
            item[parity+'_max_offgrid_residual']=max(r['offgrid_residual'] for r in ss)
            item[parity+'_Jost_third_residual']=jj['third_matching_residual']
        item['absolute_parity_difference']=float(abs(last['odd']-last['even']))
        item['Jost_parity_difference']=float(abs(jw['odd']-jw['even']))
        item['passed']=bool(max(item[p+'_resolution_spread'] for p in PARITIES)<1e-8
           and max(item[p+'_spectral_Jost_difference'] for p in PARITIES)<1e-9
           and max(item[p+'_max_offgrid_residual'] for p in PARITIES)<1e-7
           and max(item[p+'_Jost_third_residual'] for p in PARITIES)<1e-8)
        rows.append(item)
    maxima={key:max(r[key] for r in rows) for key in rows[0] if key not in ('Lambda_tilde','ell','n','passed') and 'Omega_' not in key}
    report={'scope':__doc__,'runtime':dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
      'normalization':'M=1; Omega=M*omega; varpi=Omega/[Lambda_tilde*(xc-xh)], exp(-i*omega*t).',
      'spectral_method':'Interior Chebyshev roots, N coefficients (degree N-1); direct independent parity potentials; same common row scaling of pencil coefficients; generalized companion pencil.',
      'selection':'For each parity, N and Lambda, select two least-damped positive-real roots in .5<Re(varpi)<4, -3.2<Im(varpi)<-.05; compare four resolutions, do not use odd roots to select even roots.',
      'radiation_check':'Original master equation, independent Jost roots for each parity: Nariai initializes .11; each parity separately continues .11 -> .09 -> .06. Unit leading Frobenius coefficient and p=-i*varpi/g_q(0) at each horizon.',
      'settings':{'N':ORDERS,'Lambda_tilde':LAMBDAS,'ell':ELL,'n':[0,1],
       'Jost_coarse':dict(order=24,eps=.04,match_t=.43,rtol=2e-12,atol=2e-13),
       'Jost_refined':dict(order=32,eps=.025,match_t=.57,rtol=8e-14,atol=8e-15),
       'Jost_third_evaluation':dict(order=36,eps=.06,match_t=.37,rtol=8e-14,atol=8e-15)},
      'limits':'Binary64 numerical differences/residuals are diagnostics, not rigorous error bounds. This is not historical Maple/Julia reproduction. No imaginary-mode or all-sector/all-overtone assertion.',
      'rows':rows,'maxima':maxima,'spectral':spectral,'jost':jost,'elapsed_seconds':time.time()-start}
    (OUT/'parity_results.json').write_text(json.dumps(report,indent=2)+'\n')
    with (OUT/'parity_results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    print(json.dumps({'passed':sum(r['passed'] for r in rows),'total':len(rows),'maxima':maxima,'elapsed':report['elapsed_seconds']},indent=2),flush=True)
    if not all(r['passed'] for r in rows):raise RuntimeError('One or more checks failed')

if __name__=='__main__':main()
