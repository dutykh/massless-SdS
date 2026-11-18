"""Backward continuation to check the de Sitter identification of the branch.

This auxiliary run uses exactly the radiation connection determinant and
local root search of the main campaign, with longer series when the disks'
overlap shrinks. It makes no spectral-completeness claim.
"""
from pathlib import Path
import json,time
import mpmath as mp
from ds_branch_check import connection,locate_local,geometry

HERE=Path(__file__).resolve().parent
mp.mp.dps=100
main=json.loads((HERE/'ds_extended_results.json').read_text())
base=main['path'][:2]
previous=[(mp.mpf(x['Lambda_tilde']),mp.mpf(x['alpha'])) for x in reversed(base)]
path=[mp.mpf('.06')-j*mp.mpf('.002') for j in range(1,21)]
report={'method':__doc__,'path':[],'start_Lambda':'.06','step':'-.002'}
start=time.time()
def bracket_root(lam,guess,N):
    """Safeguarded local real-axis search; never bracket across an index pole."""
    rh,rc=geometry(lam)
    kappas=[abs((1-lam*R**2)/R)/2 for R in [rh,rc]]
    width=mp.mpf('.001')
    lo,hi=guess-width,guess+width
    boundaries=[lo,hi]
    for kappa in kappas:
        for n in range(max(1,int(lo/kappa)),int(hi/kappa)+2):
            pole=n*kappa
            if lo<pole<hi:boundaries.append(pole)
    boundaries=sorted(boundaries)
    candidates=[]
    for a,b in zip(boundaries[:-1],boundaries[1:]):
        inset=(b-a)*mp.mpf('1e-8');a+=inset;b-=inset
        def fun(x):
            _,g=connection(lam,x,N,details=True)
            return (mp.mpf(g['left_value_without_nonzero_power'])*mp.mpf(g['right_r_derivative_without_nonzero_power'])
                    -mp.mpf(g['left_r_derivative_without_nonzero_power'])*mp.mpf(g['right_value_without_nonzero_power']))
        if fun(a)*fun(b)<0:
            rr=mp.findroot(fun,(a,b),solver='ridder',tol=mp.mpf('1e-50'),maxsteps=100,verify=True)
            if a<rr<b and abs(connection(lam,rr,N))<mp.mpf('1e-30'):candidates.append(rr)
    if len(candidates)!=1:raise RuntimeError(('Local bracket not unique',lam,candidates))
    return candidates[0]

for lam in path:
    rh,rc=geometry(lam);D=rc-rh
    t=min(mp.mpf('.5'),mp.mpf('.45')*rh/D)
    N=max(400,int(mp.ceil(150/(-mp.log(1-t)))))
    guess=previous[-1][1]+(lam-previous[-1][0])*(previous[-1][1]-previous[-2][1])/(previous[-1][0]-previous[-2][0])
    alpha=bracket_root(lam,guess,N)
    W,detail=connection(lam,alpha,N+200,details=True)
    distance=min(abs(mp.mpf(detail[k])-mp.nint(mp.mpf(detail[k]))) for k in ['nu_h','nu_c'])
    if abs(W)>mp.mpf('1e-30') or distance<mp.mpf('1e-7'):
        raise RuntimeError((lam,alpha,W,distance))
    item={'Lambda_tilde':str(lam),'alpha':str(alpha),'ratio_to_2sqrtLambdaOver3':str(alpha/(2*mp.sqrt(lam/3))),
          'N':N,'check_N':N+200,'normalized_W_check':str(W),
          'index_integer_distance':str(distance),'predictor_correction':str(abs(alpha-guess)),
          'geometry':detail}
    report['path'].append(item);previous.append((lam,alpha))
    print(mp.nstr(lam,6),mp.nstr(alpha,24),'ratio',mp.nstr(alpha/(2*mp.sqrt(lam/3)),15),'N',N,'elapsed',round(time.time()-start,1),flush=True)
    (HERE/'ds_low_lambda_results.json').write_text(json.dumps(report,indent=2))
report['complete']=True
report['elapsed_seconds']=time.time()-start
(HERE/'ds_low_lambda_results.json').write_text(json.dumps(report,indent=2))
