"""Bounded dense scalar ell=2 purely imaginary SdS QNM continuation.

No spectral/collocation eigenvalues are used. Radiation Frobenius series from
the original radial equation are matched. The leading dS prediction supplies
only the initial guess. Numerical, not interval-certified, accuracy evidence.
"""
from pathlib import Path
import json,time,sys,platform
import mpmath as mp
from ds_branch_check import connection,locate_local,geometry

HERE=Path(__file__).resolve().parent
mp.mp.dps=80
selected={mp.mpf(x) for x in ['.06','.08','.09','.10','.11']}
path=[mp.mpf('.06')+j*mp.mpf('.0005') for j in range(101)]
report={'method':__doc__,'ell':2,'sector':'massless scalar','M':'1',
        'runtime':{'python':sys.version,'mpmath':mp.__version__,'platform':platform.platform()},
        'path_step_Lambda_tilde':'0.0005','path_dps':80,'path_series_order':320,
        'search':'Secant solve of physical two-horizon Frobenius connection determinant in alpha=-Im(Omega); predictor +/-1e-6, escape limit0.008. No frequency-domain completeness claim.',
        'acceptance':'Normalized W<1e-30 at order360; Frobenius index distances>1e-7; comparison of frequency displacement with previous step; selected checks vary order, decimal precision and matching point.',
        'path':[],'selected':[]}
previous=[]
start=time.time()
for j,lam in enumerate(path):
    if not previous:guess=2*mp.sqrt(lam/3)
    elif len(previous)==1:guess=previous[-1][1]+(lam-previous[-1][0])*previous[-1][1]/(2*previous[-1][0])
    else:guess=previous[-1][1]+(lam-previous[-1][0])*(previous[-1][1]-previous[-2][1])/(previous[-1][0]-previous[-2][0])
    alpha=locate_local(lam,guess,320)
    val,detail=connection(lam,alpha,360,details=True)
    dist=min(abs(mp.mpf(detail[k])-mp.nint(mp.mpf(detail[k]))) for k in ['nu_h','nu_c'])
    if abs(val)>mp.mpf('1e-30') or dist<mp.mpf('1e-7'):raise RuntimeError(('Acceptance failed',lam,alpha,val,dist))
    rh,rc=geometry(lam);sf=lam*(rc-rh)/2
    item={'Lambda_tilde':str(lam),'alpha':str(alpha),'varpi_im':str(-alpha/sf),
          'rh_omega_im':str(-rh*alpha),'normalized_W_order360':str(val),
          'minimum_integer_index_distance':str(dist),'prediction_error':str(abs(alpha-guess)),
          'step_frequency_change':str(abs(alpha-previous[-1][1])) if previous else None,
          'geometry':detail}
    report['path'].append(item);previous.append((lam,alpha))
    if lam in selected:
        checks=[]
        for dps,N,where in [(90,400,None),(100,480,None),(100,480,'0.35' if lam<=mp.mpf('.06') else '0.45')]:
            with mp.workdps(dps):
                rr=locate_local(lam,alpha,N,delta='1e-12',match_fraction=where)
                ww,dd=connection(lam,rr,N,details=True,match_fraction=where)
                slope=mp.diff(lambda a:connection(lam,a,N,match_fraction=where),rr)
                checks.append({'dps':dps,'N':N,'matching_t':dd['match_t'],
                               'alpha':str(rr),'difference_path':str(abs(rr-alpha)),
                               'normalized_W':str(ww),'d_normalized_W_dalpha':str(slope),
                               'geometry_and_series':dd})
        selected_item=dict(item);selected_item['validation']=checks
        with mp.workdps(120):
            checked_roots=[mp.mpf(x['alpha']) for x in checks]
            selected_item['maximum_alpha_difference']=str(max(checked_roots)-min(checked_roots))
        report['selected'].append(selected_item)
        print('selected',lam,mp.nstr(alpha,23),'scaled',mp.nstr(alpha/sf,15),'elapsed',round(time.time()-start,1),flush=True)
    elif j%10==0:print('path',j,lam,mp.nstr(alpha,15),'elapsed',round(time.time()-start,1),flush=True)
    (HERE/'ds_continuation_results.json').write_text(json.dumps(report,indent=2))

# Bounded further approach: independent root at .111, seeded by smooth alpha.
# Included only if the same checks succeed and no resonance is encountered.
report['complete']=True
report['elapsed_seconds']=time.time()-start
(HERE/'ds_continuation_results.json').write_text(json.dumps(report,indent=2))
print('COMPLETE',len(report['path']),'path points',len(report['selected']),'validation points',flush=True)
