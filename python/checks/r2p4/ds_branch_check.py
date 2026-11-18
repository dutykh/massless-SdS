"""Independent scalar ell=2 SdS imaginary-branch Frobenius connection check.

Uses the original radial equation, not a collocation pencil. Each horizon
solution is the frequency-continued radiation Frobenius branch. Matching is
performed inside both series' convergence disks. Integer exponent differences
are excluded explicitly; no claim is made at those exceptional frequencies.
"""
from pathlib import Path
import json, math, time
import mpmath as mp

HERE=Path(__file__).resolve().parent

def add(a,b):
    ans=[mp.mpf(0)]*max(len(a),len(b))
    for k,v in enumerate(a): ans[k]+=v
    for k,v in enumerate(b): ans[k]+=v
    return ans

def mul(a,b):
    ans=[mp.mpf(0)]*(len(a)+len(b)-1)
    for j,x in enumerate(a):
        for k,y in enumerate(b): ans[j+k]+=x*y
    return ans

def scale(a,c):return [c*x for x in a]

def power(a,n):
    ans=[mp.mpf(1)]
    for j in range(n):ans=mul(ans,a)
    return ans

def geometry(lam):
    alpha=mp.acos(3*mp.sqrt(lam))
    rh=2*mp.cos((mp.pi+alpha)/3)/mp.sqrt(lam)
    rc=2*mp.cos((mp.pi-alpha)/3)/mp.sqrt(lam)
    return rh,rc

def coefficients(lam,R,D,sigma,alpha,ell=2):
    r=[R,sigma*D]
    # Exact horizon identity imposed before arithmetic expansion.
    P=[mp.mpf(0),sigma*D*(1-lam*R**2),-lam*R*D**2,-lam*sigma*D**3/3]
    Pr=add([1],scale(power(r,2),-lam))
    A=mul(power(r,2),power(P,2))
    B=scale(mul(mul(r,P),add(mul(r,Pr),scale(P,-1))),sigma*D)
    U=add(add(scale(r,ell*(ell+1)),[2]),scale(power(r,3),-2*lam/3))
    C=scale(add(scale(power(r,4),-alpha**2),scale(mul(P,U),-1)),D**2)
    # rho is the chosen radiation exponent at either endpoint.
    rho=-alpha/abs((1-lam*R**2)/R)
    return A,B,C,rho

def series(lam,R,D,sigma,alpha,q,N,ell=2):
    A,B,C,rho=coefficients(lam,R,D,sigma,alpha,ell)
    coeff=[mp.mpf(1)]
    for n in range(1,N+1):
        total=mp.mpf(0)
        for k in range(1,min(n,max(len(A)-3,len(B)-2,len(C)-1))+1):
            p=rho+n-k
            value=(A[k+2]*p*(p-1) if k+2<len(A) else 0)
            value+=(B[k+1]*p if k+1<len(B) else 0)
            value+=(C[k] if k<len(C) else 0)
            total+=value*coeff[n-k]
        den=A[2]*n*(n+2*rho)
        coeff.append(-total/den)
    F=mp.polyval(list(reversed(coeff)),q)
    Fq=mp.polyval([n*coeff[n] for n in range(N,0,-1)],q)
    # Remove the nonzero q^rho normalization factor before W calculation.
    derivative=sigma*(Fq+rho*F/q)/D
    tail=max(abs(coeff[-j]*q**(N+1-j)) for j in range(1,5))
    return F,derivative,rho,tail/(abs(F)+mp.mpf('1e-100'))

def connection(lam,alpha,N=300,match_factor=mp.mpf('.45'),details=False,match_fraction=None):
    rh,rc=geometry(lam);D=rc-rh
    t=min(mp.mpf('.5'),match_factor*rh/D) if match_fraction is None else mp.mpf(match_fraction)
    assert 0<t<1 and t<min(1,rh/D) and 1-t<1
    left=series(lam,rh,D,1,alpha,t,N)
    right=series(lam,rc,D,-1,alpha,1-t,N)
    W=left[0]*right[1]-left[1]*right[0]
    # For real alpha all entries and W are real. Normalization is nonzero.
    norm=mp.sqrt(left[0]**2+left[1]**2)*mp.sqrt(right[0]**2+right[1]**2)
    val=W/norm
    if details:return val, {'rh':str(rh),'rc':str(rc),'match_t':str(t),
        'left_series_radius_t':str(min(mp.mpf(1),rh/D)),
        'right_series_radius_1_minus_t':'1',
        'left_value_without_nonzero_power':str(left[0]),
        'left_r_derivative_without_nonzero_power':str(left[1]),
        'right_value_without_nonzero_power':str(right[0]),
        'right_r_derivative_without_nonzero_power':str(right[1]),
        'nu_h':str(-2*left[2]),'nu_c':str(-2*right[2]),
        'tail_left':str(left[3]),'tail_right':str(right[3])}
    return val

def locate(lam,guess,N=300):
    guess=mp.mpf(guess)
    return mp.findroot(lambda a:connection(lam,a,N), (guess*(1-mp.mpf('.001')),guess*(1+mp.mpf('.001'))),tol=mp.mpf('1e-35'),maxsteps=40,verify=True)

def locate_local(lam,guess,N=320,delta='1e-6',match_fraction=None):
    guess=mp.mpf(guess);delta=mp.mpf(delta)
    fun=lambda a:connection(lam,a,N,match_fraction=match_fraction)
    root=mp.findroot(fun,(guess-delta,guess+delta),tol=mp.mpf('1e-35'),maxsteps=30,verify=True)
    if abs(root-guess)>mp.mpf('.008'):raise RuntimeError('Local root search escaped predictor neighborhood')
    return root

if __name__=='__main__':
    # A bounded smoke test. The reproducible campaigns are the separate
    # ds_continuation.py and ds_extend_continuation.py entry points.
    mp.mp.dps=80
    lam=mp.mpf('.06');start=time.time()
    root=locate_local(lam,2*mp.sqrt(lam/3),N=320)
    value,details=connection(lam,root,400,details=True)
    print('Lambda=.06',mp.nstr(root,25),mp.nstr(value,5),details,
          'seconds',time.time()-start,flush=True)
