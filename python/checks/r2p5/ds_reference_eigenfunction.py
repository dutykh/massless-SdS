"""Independent Frobenius reference for scalar ell=2, Lambda_tilde=.06.

Phi is defined by psi=(1+y)^rho_h*(1-y)^rho_c*Phi, with
rho_h,c=-alpha/(2*kappa_h,c), Omega=-i*alpha, M=1, and Phi(-.3)=1.
This uses radiation-series solutions of the original radial equation, not
spectral eigenvectors. See ../r2p4/ds_method_summary.md for the recurrence.
"""
from pathlib import Path
import sys,json,time
import mpmath as mp
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'r2p4'))
from ds_branch_check import coefficients,geometry

def radiation_coefficients(lam,R,D,sigma,alpha,N):
    A,B,C,rho=coefficients(lam,R,D,sigma,alpha,ell=2)
    ans=[mp.mpf(1)]
    for n in range(1,N+1):
        total=mp.mpf(0)
        for k in range(1,min(n,max(len(A)-3,len(B)-2,len(C)-1))+1):
            p=rho+n-k
            value=(A[k+2]*p*(p-1) if k+2<len(A) else 0)
            value+=(B[k+1]*p if k+1<len(B) else 0)
            value+=(C[k] if k<len(C) else 0)
            total+=value*ans[n-k]
        ans.append(-total/(A[2]*n*(n+2*rho)))
    return ans,rho

class ScalarImaginaryReference:
    def __init__(self,N=500):
        data=json.loads((HERE.parent/'r2p4'/'ds_extended_results.json').read_text())
        row=min(data['selected'],key=lambda r:abs(mp.mpf(r['Lambda_tilde'])-mp.mpf('.06')))
        self.lam=mp.mpf('.06')
        self.alpha=mp.mpf(row['validation'][-1]['alpha'])
        self.rh,self.rc=geometry(self.lam);self.D=self.rc-self.rh
        self.ah,self.rho_h=radiation_coefficients(self.lam,self.rh,self.D,1,self.alpha,N)
        self.ac,self.rho_c=radiation_coefficients(self.lam,self.rc,self.D,-1,self.alpha,N)
        self.ah=list(reversed(self.ah));self.ac=list(reversed(self.ac))
        self.match=mp.mpf('.35')
        Fh=mp.polyval(self.ah,self.match)
        Fc=mp.polyval(self.ac,1-self.match)
        # Physical psi_h = C psi_c at the interior match point.
        self.C=(self.match**self.rho_h*Fh)/((1-self.match)**self.rho_c*Fc)
        # Common factor 2^(-rho_h-rho_c) cancels in normalization.
        self.normalization=(1-self.match)**(-self.rho_c)*Fh

    def phi(self,y):
        y=mp.mpf(y);t=(1+y)/2
        if not 0<=t<=1:raise ValueError('Require -1<=y<=1')
        if t<=self.match:
            val=(1-t)**(-self.rho_c)*mp.polyval(self.ah,t)
        else:
            val=self.C*t**(-self.rho_h)*mp.polyval(self.ac,1-t)
        return val/self.normalization

if __name__=='__main__':
    start=time.time()
    with mp.workdps(100):
        first=ScalarImaginaryReference(500)
        grid=[mp.cos(mp.pi*j/1024) for j in range(1025)]
        first_values=[first.phi(y) for y in grid]
        report={'method':__doc__,'normalization_y':'-.3','normalization_value':'1',
                'Lambda_tilde':str(first.lam),'alpha':str(first.alpha),
                'rho_h':str(first.rho_h),'rho_c':str(first.rho_c),
                'grid':'y_j=cos(pi*j/1024), j=0,...,1024',
                'matching_t':str(first.match),
                'settings':[{'dps':100,'N':500},{'dps':130,'N':650}]}
    with mp.workdps(130):
        second=ScalarImaginaryReference(650)
        grid=[mp.cos(mp.pi*j/1024) for j in range(1025)]
        values=[second.phi(y) for y in grid]
        report['max_absolute_difference']=str(max(abs(a-b) for a,b in zip(first_values,values)))
        report['max_relative_sup_difference']=str(max(abs(a-b) for a,b in zip(first_values,values))/max(abs(v) for v in values))
        report['phi_normalization_check']=str(second.phi(mp.mpf('-.3')))
        report['points']=[{'j':j,'y':str(y),'phi':str(v)} for j,(y,v) in enumerate(zip(grid,values))]
    report['elapsed_seconds']=time.time()-start
    out=HERE/'ds_reference_eigenfunction.json'
    out.write_text(json.dumps(report,indent=2))
    print('Saved',out,'max absolute difference',report['max_absolute_difference'],flush=True)
