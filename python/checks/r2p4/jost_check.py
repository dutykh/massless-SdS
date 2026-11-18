"""Independent two-horizon Frobenius/Jost shooting for low SdS modes.

No imported spectral matrices or eigenvalues are used. The radial master
functions themselves are integrated in t=(r-rh)/(rc-rh). Endpoint data select
the actual ingoing/outgoing Frobenius sectors; the regularity of an unspecified
factorized remainder is not used as a boundary criterion.

M=1, varpi=omega/(Lambda*(rc-rh)/2).
"""
from pathlib import Path
import json,math,time,csv,platform
import scipy
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

BASE=Path(__file__).resolve().parent

def geom(lam):
    a=np.arccos(3*np.sqrt(lam))
    xh=np.cos((np.pi+a)/3)/np.sqrt(lam)
    xc=np.cos((np.pi-a)/3)/np.sqrt(lam)
    rh,rc=2*xh,2*xc
    return dict(lam=float(lam),rh=rh,rc=rc,R=rc-rh,S=lam*(rc-rh)/2,xh=xh,xc=xc)

def series_data(G,spin,ell,side,order):
    lam,rh,rc,R,S=[G[k] for k in ['lam','rh','rc','R','S']]
    a,b,sgn=(rh,R,1) if side=='left' else (rc,-R,-1)
    N=order+4
    k=np.arange(N)
    inv1=(-b/a)**k/a
    inv2=(k+1)*(-b/a)**k/a**2
    inv3=(k+1)*(k+2)/2*(-b/a)**k/a**3
    ratio=(rh+rc)*inv1;ratio[0]+=1
    gs=np.zeros(N)
    gs[1:]+=2/3*ratio[:-1]
    gs[2:]-=2/3*ratio[:-2]
    gp=np.arange(1,N)*gs[1:]
    B=ell*(ell+1)*inv2
    if spin==0:
        prod=np.convolve(gp,inv1)[:N]
        B+=S*sgn*prod
    if spin==2:B-=6*inv3
    Vs=2/lam*np.convolve(gs,B)[:N]
    P=np.convolve(gs[1:],gs[1:])[:order+1]
    Q=np.convolve(gs[1:],gp)[:order+1]
    return gs[1],P,Q,Vs[:order+1]

class Matcher:
    def __init__(self,lam,spin,ell,order=24,eps=.04,mid=.43,rtol=2e-12,atol=2e-13):
        self.G=geom(lam);self.spin=spin;self.ell=ell
        self.order=order;self.eps=eps;self.mid=mid;self.rtol=rtol;self.atol=atol
        self.series={side:series_data(self.G,spin,ell,side,order) for side in ['left','right']}
        self.calls=0
    def coefficients(self,q,side):
        G=self.G;rh,rc,R,S,lam=[G[k] for k in ['rh','rc','R','S','lam']]
        a,b,sgn=(rh,R,1) if side=='left' else (rc,-R,-1)
        r=a+b*q
        v=1+(rh+rc)/r
        vp=-(rh+rc)*b/r**2
        g=2/3*q*(1-q)*v
        gp=2/3*((1-2*q)*v+q*(1-q)*vp)
        B=self.ell*(self.ell+1)/r**2
        if self.spin==0:B+=S*sgn*gp/r
        if self.spin==2:B-=6/r**3
        V=2*g/lam*B
        return g,gp,V
    def local_data(self,w,side):
        slope,P,Q,V=self.series[side]
        p=-1j*w/slope
        coeff=np.zeros(self.order+1,dtype=complex);coeff[0]=1
        for n in range(1,self.order+1):
            total=0j
            for k in range(1,n+1):
                h=p+n-k
                total+=(P[k]*h*(h-1)+Q[k]*h-V[k])*coeff[n-k]
            den=slope*slope*n*(n+2*p)
            coeff[n]=-total/den
        q=self.eps
        poly=np.polynomial.polynomial.polyval(q,coeff)
        deriv=np.polynomial.polynomial.polyval(q,np.arange(1,len(coeff))*coeff[1:])
        y=q**p*poly
        yp=q**(p-1)*(p*poly+q*deriv)
        return np.array([y,yp]),coeff
    def propagate(self,w,side):
        initial,coeff=self.local_data(w,side)
        end=self.mid if side=='left' else 1-self.mid
        def rhs(q,Y):
            g,gp,V=self.coefficients(q,side)
            return np.array([Y[1],-gp/g*Y[1]-(w*w-V)/(g*g)*Y[0]])
        sol=solve_ivp(rhs,(self.eps,end),initial,method='DOP853',rtol=self.rtol,atol=self.atol)
        if not sol.success:raise RuntimeError(sol.message)
        return sol.y[:,-1],sol.nfev
    def determinant(self,w):
        self.calls+=1
        L,_=self.propagate(w,'left');R,_=self.propagate(w,'right')
        W=-(L[0]*R[1]+L[1]*R[0])
        return W/(np.linalg.norm(L)*np.linalg.norm(R))
    def solve(self,seed):
        def realfun(z):
            v=self.determinant(z[0]+1j*z[1]);return [v.real,v.imag]
        sol=root(realfun,[seed.real,seed.imag],tol=1e-10)
        w=sol.x[0]+1j*sol.x[1];res=abs(self.determinant(w))
        if res>1e-8:raise RuntimeError((self.G['lam'],self.spin,w,res,sol.message))
        return w,res,sol.success,str(sol.message)

def run(lams=(.11,.111,.1111,.11111),spins=(0,1,2),ells=(2,),ns=(0,1)):
    rows=[]
    for ell in ells:
      for spin in spins:
        A=ell*(ell+1)-(2 if spin==2 else 0)
        for n in ns:
          exact=complex(np.sqrt(A-.25),-n-.5)
          seed=exact
          for lam in reversed(lams):
            start=time.time()
            matcher=Matcher(lam,spin,ell)
            w,res,ok,msg=matcher.solve(seed)
            fine=Matcher(lam,spin,ell,order=32,eps=.025,mid=.57,rtol=8e-14,atol=8e-15)
            wf,resf,okf,msgf=fine.solve(w)
            # Evaluate refined root independently with a third starting point,
            # series order, and match location; no additional root search.
            third=Matcher(lam,spin,ell,order=36,eps=.06,mid=.37,rtol=8e-14,atol=8e-15)
            rest=abs(third.determinant(wf))
            G=matcher.G
            row=dict(lam=lam,spin=spin,ell=ell,n=n,varpi_re=wf.real,varpi_im=wf.imag,Omega_re=G['S']*wf.real,Omega_im=G['S']*wf.imag,scale=G['S'],nariai_re=exact.real,nariai_im=exact.imag,abs_difference=abs(wf-exact),relative_difference=abs(wf-exact)/abs(exact),refinement_difference=abs(wf-w),refined_normalized_residual=resf,independent_match_residual=rest,coarse_root_success=ok,refined_root_success=okf,elapsed_seconds=time.time()-start)
            rows.append(row);seed=wf
            print(json.dumps(row),flush=True)
    result={'method':'Independent original-master-equation Frobenius shooting with DOP853 integration and two-horizon Wronskian matching. No spectral eigenvalues used as seeds; Nariai analytic frequencies initialize the nearest-extremal case, then parameter continuation supplies seeds.','geometry':'M=1, t=(r-r_h)/(r_c-r_h), g=f/[Lambda (r_c-r_h)^2/2], frequency varpi=omega/[Lambda (r_c-r_h)/2].','endpoint_choice':'At either endpoint q=0 the explicitly selected radiation sector is psi=q^p sum a_n q^n, p=-i varpi/gprime(0); a0=1.','coarse':{'series_order':24,'endpoint_offset':.04,'match_t':.43,'rtol':2e-12,'atol':2e-13},'refined':{'series_order':32,'endpoint_offset':.025,'match_t':.57,'rtol':8e-14,'atol':8e-15},'third_check':{'series_order':36,'endpoint_offset':.06,'match_t':.37,'rtol':8e-14,'atol':8e-15},'frequency_convention':'The code variable w is the manuscript varpi=Omega/[Lambda_tilde(xc-xh)] with Omega=M*omega, M=1, and harmonic dependence exp(-i*omega*t_physical). Negative imaginary w is damping.', 'matching_residual':'Absolute local Wronskian divided by the product of Euclidean norms of (psi,dpsi/dt) for the two solutions, where t=(r-rh)/(rc-rh): abs(psiH*psiC_t-psiH_t*psiC)/(norm([psiH,psiH_t])*norm([psiC,psiC_t])). The right integration uses q=1-t, so psiC_t=-psiC_q. This is a normalized matching determinant, not an absolute Green-function-residue diagnostic.', 'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__}, 'limitations':'Double-precision numerical matching, not interval certification. Only the stated low-lying positive-real photon-branch continuations are checked; no de Sitter-family or spectral-completeness claim.','rows':rows}
    (BASE/'jost_check_results.json').write_text(json.dumps(result,indent=2))
    with (BASE/'jost_check_results.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    return rows

if __name__=='__main__':run()
