#!/usr/bin/env python3
"""New binary64 SdS convergence experiment at eight substantially separated N.

Candidate selection is spectral and resolution-continuous; horizon reference
frequencies enter only error assessment after selection. No radiation claim is
deduced from convergence, coefficient decay, or a differential residual alone.
"""
from pathlib import Path
import csv
import json
import sys
import numpy as np
import scipy
from numpy.polynomial.chebyshev import chebval,chebder
from scipy.linalg import eig
from threadpoolctl import threadpool_limits

HERE=Path(__file__).resolve().parent
PRIOR=HERE.parent/'r2p4'
sys.path.insert(0,str(PRIOR))
from continue_near_nariai import basis,coefficients,geometry

ORDERS=[12,16,24,32,48,64,96,128]
Y=np.cos(np.pi*np.arange(1025)/1024)
Y_ANCHOR=-.3
CASES=[
 {'name':'scalar_oscillatory','sector':'scalar','Lambda_tilde':.09,'ell':2,'n':1,'kind':'oscillatory'},
 {'name':'axial_oscillatory','sector':'odd','Lambda_tilde':.11,'ell':2,'n':0,'kind':'oscillatory'},
 {'name':'scalar_imaginary','sector':'scalar','Lambda_tilde':.06,'ell':2,'n':0,'kind':'imaginary'},
]

def references():
 extra=json.loads((PRIOR/'jost_extra_cell_check.json').read_text())['target']
 anchor=json.loads((PRIOR/'jost_check_results.json').read_text())['rows']
 axial=next(x for x in anchor if x['spin']==2 and x['n']==0 and x['lam']==.11)
 ds=json.loads((PRIOR/'ds_extended_results.json').read_text())['selected']
 pure=next(x for x in ds if float(x['Lambda_tilde'])==.06)
 return {
  'scalar_oscillatory':{'re':extra['varpi_re'],'im':extra['varpi_im'],'source':'r2p4/jost_extra_cell_check.json, target','reference_refinement_change':extra['coarse_refined_difference']},
  'axial_oscillatory':{'re':axial['varpi_re'],'im':axial['varpi_im'],'source':'r2p4/jost_check_results.json, spin2 ell2 n0 lam.11'},
  'scalar_imaginary':{'re':0.,'im':float(pure['varpi_im']),'source':'r2p4/ds_extended_results.json, selected lam.06','Omega_im_exact_record':str(-float(pure['alpha'])),'alpha_decimal':pure['alpha'],'varpi_im_decimal':pure['varpi_im']},
 }

def all_eigenpairs(N,case):
 y,T,T1,T2=basis(N)
 g,gy,G,Gy,R,V=coefficients(y,case['Lambda_tilde'],case['sector'])
 A0=g[:,None]*T2+gy[:,None]*T1-V[:,None]*T
 A1=1j*(2*G[:,None]*T1+Gy[:,None]*T)
 A2=R[:,None]*T
 row=np.max(abs(np.hstack([A0,A1,A2])),axis=1)
 A0/=row[:,None];A1/=row[:,None];A2/=row[:,None]
 I=np.eye(N);Z=np.zeros((N,N))
 val,vec=eig(np.block([[Z,I],[-A0,-A1]]),np.block([[I,Z],[Z,A2]]))
 return val,vec[:N,:],(A0,A1,A2)

def candidate_indices(vals,case):
 if case['kind']=='oscillatory':
  window=(vals.real>.5)&(vals.real<4)&(vals.imag> -3.2)&(vals.imag<-.05)
 else:
  window=(abs(vals.real)<1e-6)&(vals.imag> -8)&(vals.imag<-.05)
 return np.flatnonzero(np.isfinite(vals)&window)

def independent_master_residual(w,c,case):
 # Compute the original unfactored SdS equation divided by the nonzero
 # radiation prefactor at interior test points. The assembly's cancellation-
 # reduced R,G coefficient formulas are not used in this residual.
 yy=Y[1:-1];lam=case['Lambda_tilde'];bg=geometry(lam)
 h,cx,x0,d=(bg[k] for k in ['xh','xc','x0','d'])
 x=(h+cx)/2+d*yy/2
 f=lam*d*d*(1-yy*yy)*(x-x0)/(3*x)
 fx=1/(x*x)-8*lam*x/3
 g=f/(lam*d*d);gy=fx/(2*lam*d)
 ph=3*h/(2*(h-x0));pc=3*cx/(2*(cx-x0))
 H=-ph/(1+yy)+pc/(1-yy)
 Hp=ph/(1+yy)**2+pc/(1-yy)**2
 Vx=f*(6/x**2+(fx/x if case['sector']=='scalar' else 0)
       -(3/x**3 if case['sector']=='odd' else 0))
 U=Vx/(4*(lam*d)**2)
 p=chebval(yy,c);p1=chebval(yy,chebder(c));p2=chebval(yy,chebder(c,2))
 terms=np.array([g*g*p2,2j*w*g*g*H*p1,1j*w*g*g*Hp*p,
                 -w*w*g*g*H*H*p,g*gy*p1,1j*w*g*gy*H*p,w*w*p,-U*p])
 scale=np.max(np.sum(abs(terms),axis=0))
 return float(np.max(abs(np.sum(terms,axis=0)))/scale)

def eigenpair_diagnostics(w,coeff,mats,case):
 norm=chebval(Y_ANCHOR,coeff)
 if abs(norm)<1e-8*np.max(abs(coeff)):
  raise RuntimeError('Chosen normalization point is too close to a node')
 c=coeff/norm
 A0,A1,A2=mats
 pr=np.linalg.norm((A0+w*A1+w*w*A2)@c)/(np.linalg.norm(c)*(
     np.linalg.norm(A0)+abs(w)*np.linalg.norm(A1)+abs(w*w)*np.linalg.norm(A2)))
 phi=chebval(Y,c)
 p1=chebval(Y,chebder(c));p2=chebval(Y,chebder(c,2))
 g,gy,G,Gy,R,V=coefficients(Y,case['Lambda_tilde'],case['sector'])
 terms=np.array([g*p2,(gy+2j*w*G)*p1,(-V+1j*w*Gy+w*w*R)*phi])
 scale=np.max(np.sum(abs(terms),axis=0))
 residual=np.sum(terms,axis=0)
 metrics={
  'relative_pencil_residual':float(pr),
  'independent_original_ODE_residual':independent_master_residual(w,c,case),
  'relative_reduced_ODE_residual_including_endpoints':float(np.max(abs(residual))/scale),
  'relative_reduced_endpoint_residual':float(np.max(abs(residual[[0,-1]]))/scale),
  'last_eight_coefficient_norm_fraction':float(np.linalg.norm(c[-8:])/np.linalg.norm(c)),
 }
 return c,phi,metrics

def main():
 refs=references();allrows=[];cases=[]
 field_file=HERE/'ds_reference_eigenfunction.json'
 field_reference=json.loads(field_file.read_text()) if field_file.exists() else None
 offaxis_file=HERE/'offaxis_reference_eigenfunctions.json'
 offaxis_reference=json.loads(offaxis_file.read_text()) if offaxis_file.exists() else None
 with threadpool_limits(limits=1):
  for case in CASES:
   rows=[];previous=None;prev_phi=None
   for N in ORDERS:
    vals,vec,mats=all_eigenpairs(N,case);idx=candidate_indices(vals,case)
    if not len(idx):raise RuntimeError(f'No spectral candidate at {case["name"]},N{N}')
    if previous is None:
     ordered=idx[np.argsort(-vals[idx].imag)]
     chosen=ordered[case['n']]
    else:
     chosen=idx[np.argmin(abs(vals[idx]-previous))]
    w=vals[chosen];alternative=np.delete(vals[idx],np.flatnonzero(idx==chosen))
    gap=float(np.min(abs(w-alternative))) if len(alternative) else None
    step=None if previous is None else float(abs(w-previous))
    if step is not None and (step>.1 or (gap is not None and step>.25*gap)):
     raise RuntimeError('Ambiguous resolution association')
    coeff,phi,metrics=eigenpair_diagnostics(w,vec[:,chosen],mats,case)
    ref=refs[case['name']];rw=complex(ref['re'],ref['im']);sf=geometry(case['Lambda_tilde'])['scale']
    field_error=None;field_ref_variation=None;field_ref_variation_abs=None
    if case['kind']=='imaginary' and field_reference:
     independent_phi=np.array([float(p['phi']) for p in field_reference['points']])
     field_error=float(np.max(abs(phi-independent_phi))/np.max(abs(independent_phi)))
     field_ref_variation=float(field_reference['max_relative_sup_difference'])
     field_ref_variation_abs=float(field_reference['max_absolute_difference'])
    elif case['kind']=='oscillatory' and offaxis_reference:
     refcase=next(x for x in offaxis_reference['cases'] if x['definition']['name']==case['name'])
     independent_phi=np.array(refcase['phi_re'])+1j*np.array(refcase['phi_im'])
     field_error=float(np.max(abs(phi-independent_phi))/np.max(abs(independent_phi)))
     field_ref_variation=refcase['reference_configuration_difference_relative_Linf']
     field_ref_variation_abs=refcase['reference_configuration_difference_absolute_Linf']
    row={'case':case['name'],'N':N,'Lambda_tilde':case['Lambda_tilde'],
         'varpi_re':float(w.real),'varpi_im':float(w.imag),
         'Omega_re':float((w*sf).real),'Omega_im':float((w*sf).imag),
         'error_to_independent_varpi_reference':float(abs(w-rw)),
         'error_to_independent_Omega_reference':float(abs(w-rw)*sf),
         'change_from_previous_N':step,'nearest_candidate_separation_in_window':gap,
         'eigenfunction_change_from_previous_N_Linf':None if prev_phi is None else float(np.max(abs(phi-prev_phi))/np.max(abs(prev_phi))),
         'eigenfunction_error_to_independent_Frobenius_Linf':field_error,
         'independent_eigenfunction_configuration_difference_Linf':field_ref_variation,
         'independent_eigenfunction_configuration_difference_absolute_Linf':field_ref_variation_abs,
         'independent_eigenfunction_configuration_difference_relative_Linf':field_ref_variation,
         **metrics,
         'coeff_re':coeff.real.tolist(),'coeff_im':coeff.imag.tolist(),
         'coefficient_absolute_values':abs(coeff).tolist(),
         'phi_re':phi.real.tolist(),'phi_im':phi.imag.tolist()}
    rows.append(row);previous=w;prev_phi=phi
    print(case['name'],N,'w=',w,'reference_error=',row['error_to_independent_varpi_reference'],'rawres=',row['independent_original_ODE_residual'],flush=True)
   highest=complex(1)*np.array(rows[-1]['phi_re'])+1j*np.array(rows[-1]['phi_im'])
   for row in rows:
    pp=np.array(row['phi_re'])+1j*np.array(row['phi_im'])
    row['eigenfunction_difference_from_N128_Linf']=float(np.max(abs(pp-highest))/np.max(abs(highest)))
   cases.append({'definition':case,'independent_frequency_reference':refs[case['name']],
                 'results':rows})
   allrows.extend(rows)
 report={'scope':__doc__,'precision':'IEEE binary64',
         'versions':{'numpy':np.__version__,'scipy':scipy.__version__},
         'N':ORDERS,'frequency':'varpi=Omega/[Lambda_tilde*(xc-xh)]',
         'selection':'At N12 choose least damped imaginary candidate, or positive-real oscillatory overtone n by damping order, in fixed windows. At every subsequent N choose nearest previous spectral root; require displacement<.1 and <one quarter of nearest competing-candidate separation when a competitor exists. Independent horizon references never enter selection.',
         'windows':{'oscillatory':'0.5<Re(varpi)<4,-3.2<Im(varpi)<-.05',
                    'imaginary':'abs(Re(varpi))<1e-6,-8<Im(varpi)<-.05'},
         'normalization':'Phi(y=-.3)=1; complex division fixes amplitude and phase. Differences are relative common-grid maximum norms.',
         'reference_variation_norms':'Explicit absolute and relative reference-configuration differences are recorded. Relative differences divide by the reference function maximum norm. The legacy independent_eigenfunction_configuration_difference_Linf field is an alias of the RELATIVE value.',
         'residual':'Original unfactored master equation evaluated after dividing out the nonzero radiation factor at 1023 interior test points, using metric and derivative; denominator is maximum sum of absolute values of the eight unfactored terms. Separate reduced-ODE endpoint check is also recorded.',
         'grid_y':Y.tolist(),'cases':cases,
         'independent_imaginary_eigenfunction_reference':'ds_reference_eigenfunction.json; sampled physical horizon Frobenius solution, same normalization; high-precision grid rounded to binary64 for comparison' if field_reference else None,
         'independent_offaxis_eigenfunction_reference':'offaxis_reference_eigenfunctions.json; original-equation Frobenius and DOP853 solutions; reference variation is retained and is not an exact error bound' if offaxis_reference else None,
         'limitations':'References are numerical, not rigorous error enclosures. N128 eigenfunction differences are interresolution diagnostics, not error relative to an exact function. The scalar imaginary eigenvalue deteriorates at high N in binary64; no monotonic convergence or boundary selection is inferred from small residuals or coefficient tails.'}
 (HERE/'broad_resolution_results.json').write_text(json.dumps(report,indent=2)+'\n')
 compact=[{k:v for k,v in row.items() if not isinstance(v,list)} for row in allrows]
 with (HERE/'broad_resolution_results.csv').open('w',newline='') as stream:
  writer=csv.DictWriter(stream,fieldnames=list(compact[0]));writer.writeheader();writer.writerows(compact)

if __name__=='__main__':main()
