"""Independent endpoint-complete Phi references from physical Jost shooting.

Uses the already checked off-axis frequencies only; no spectral coefficients
or eigenvectors enter the reference. Two endpoint Frobenius solutions of the
original master equation are continued to t=1/2 and matched there. Their
known endpoint series, rather than singular psi evaluation, supply Phi at
the endpoints. Every returned Phi is normalized to Phi(y=-0.3)=1.
"""
from pathlib import Path
import sys,json
import numpy as np
from scipy.integrate import solve_ivp

BASE=Path(__file__).resolve().parent
sys.path.insert(0,str(BASE.parent/'r2p4'))
from jost_check import Matcher

def reference(lam,spin,ell,w,y,order=32,eps=.025,rtol=8e-14,atol=8e-15):
    m=Matcher(lam,spin,ell,order=order,eps=eps,mid=.5,rtol=rtol,atol=atol)
    data={}
    for side in ('left','right'):
        initial,coeff=m.local_data(w,side)
        p=-1j*w/m.series[side][0]
        def rhs(q,Y):
            g,gp,V=m.coefficients(q,side)
            return [Y[1],-gp/g*Y[1]-(w*w-V)/(g*g)*Y[0]]
        sol=solve_ivp(rhs,(eps,.5),initial,method='DOP853',dense_output=True,rtol=rtol,atol=atol)
        if not sol.success:raise RuntimeError(sol.message)
        data[side]=dict(p=p,coeff=coeff,sol=sol)
    L=data['left']['sol'].y[:,-1]
    R=data['right']['sol'].y[:,-1]*np.array([1,-1])
    factor=np.vdot(R,L)/np.vdot(R,R)
    mismatch=np.linalg.norm(L-factor*R)/np.linalg.norm(L)
    ph,pc=data['left']['p'],data['right']['p']
    def value(t):
        side='left' if t<=.5 else 'right'
        q=t if side=='left' else 1-t
        rec=data[side]
        if q<=eps:
            otherp=pc if side=='left' else ph
            phi=np.polynomial.polynomial.polyval(q,rec['coeff'])/(1-q)**otherp
        else:
            psi=rec['sol'].sol(q)[0]
            phi=psi/(t**ph*(1-t)**pc)
        return phi if side=='left' else factor*phi
    normalization=value(.35)
    phi=np.array([value((1+v)/2)/normalization for v in y])
    return phi,dict(series_order=order,endpoint_offset=eps,match_t=.5,rtol=rtol,atol=atol,relative_match_vector_error=float(mismatch),left_power_re=ph.real,left_power_im=ph.imag,right_power_re=pc.real,right_power_im=pc.imag,normalization_t=.35)

def run():
    source=json.loads((BASE/'broad_resolution_results.json').read_text())
    y=np.array(source['grid_y']);rows=[]
    for case in source['cases']:
        d=case['definition']
        if d['kind']!='oscillatory':continue
        ref=case['independent_frequency_reference']
        w=complex(ref['re'],ref['im'])
        spin=0 if d['sector']=='scalar' else 2
        a,sa=reference(d['Lambda_tilde'],spin,d['ell'],w,y)
        b,sb=reference(d['Lambda_tilde'],spin,d['ell'],w,y,order=40,eps=.04)
        reference_norm=float(np.max(np.abs(b)))
        configuration_absolute=float(np.max(np.abs(a-b)))
        errors=[]
        for r in case['results']:
            numeric=np.array(r['phi_re'])+1j*np.array(r['phi_im'])
            absolute=float(np.max(np.abs(numeric-b)))
            errors.append(dict(N=r['N'],eigenfunction_error_to_independent_Frobenius_Linf=absolute,eigenfunction_difference_absolute_Linf=absolute,eigenfunction_difference_relative_Linf=absolute/reference_norm))
        rows.append(dict(definition=d,frequency=ref,configuration_a=sa,configuration_b=sb,reference_function_Linf_norm=reference_norm,reference_configuration_difference_Linf=configuration_absolute,reference_configuration_difference_absolute_Linf=configuration_absolute,reference_configuration_difference_relative_Linf=configuration_absolute/reference_norm,phi_re=b.real.tolist(),phi_im=b.imag.tolist(),spectral_comparison=errors))
    result=dict(method='Original-master-equation two-horizon Jost Frobenius series plus dense DOP853 integration, using independently checked fixed frequencies and no spectral eigenvectors.',factoring='Phi=psi/[t^p_h(1-t)^p_c], p_h,p_c the physical radiation Frobenius powers; this differs from factoring (1+y),(1-y) only by a constant removed by normalization.',normalization='Phi(y=-.3)=1, performed by direct reference evaluation, not interpolation.',norm_definitions='Absolute Linf is max_j abs(Phi_a-Phi_b). Relative Linf divides this by max_j abs(Phi_b), configuration_b. The legacy fields reference_configuration_difference_Linf and nested eigenfunction_error_to_independent_Frobenius_Linf retain ABSOLUTE values; use the explicit absolute/relative fields.',grid_y=y.tolist(),scope='Binary64 numerical references, with two endpoint offsets and series orders. Configuration difference is a practical consistency diagnostic, not a rigorous error bound; fixed-frequency uncertainty remains.',cases=rows)
    (BASE/'offaxis_reference_eigenfunctions.json').write_text(json.dumps(result,indent=2))
    for r in rows:print(r['definition']['name'],'absolute',r['reference_configuration_difference_absolute_Linf'],'relative',r['reference_configuration_difference_relative_Linf'],r['configuration_b']['relative_match_vector_error'],r['spectral_comparison'])

if __name__=='__main__':run()
