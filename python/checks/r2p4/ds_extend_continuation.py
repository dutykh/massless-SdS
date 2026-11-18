"""Continue the directly matched scalar de Sitter branch closer to Nariai.

Run ds_continuation.py first. This extends its physical-connection path and
does not use the collocation matrix or its eigenvalue output.
"""
from pathlib import Path
import json,time,csv
import mpmath as mp
from ds_branch_check import connection,locate_local,geometry

HERE=Path(__file__).resolve().parent
report=json.loads((HERE/'ds_continuation_results.json').read_text())
assert report['complete'] and len(report['path'])==101
mp.mp.dps=180
steps=[(mp.mpf('.11')+j*mp.mpf('.00005'),120,600,720) for j in range(1,21)]
steps += [(mp.mpf('.111')+j*mp.mpf('.000005'),180,900,1050) for j in range(1,21)]
previous=[(mp.mpf(q['Lambda_tilde']),mp.mpf(q['alpha'])) for q in report['path'][-2:]]
report['extension_method']=__doc__
report['extension_steps']=['.11 to .111: step .00005, 120digits, order600/720',
                           '.111 to .1111: step .000005, 180digits, order900/1050']
start=time.time()
for j,(lam,dps,N,Ncheck) in enumerate(steps):
    with mp.workdps(dps):
        guess=previous[-1][1]+(lam-previous[-1][0])*(previous[-1][1]-previous[-2][1])/(previous[-1][0]-previous[-2][0])
        alpha=locate_local(lam,guess,N,delta='1e-8')
        val,detail=connection(lam,alpha,Ncheck,details=True)
        dist=min(abs(mp.mpf(detail[k])-mp.nint(mp.mpf(detail[k]))) for k in ['nu_h','nu_c'])
        if abs(val)>mp.mpf('1e-30') or dist<mp.mpf('1e-7'):raise RuntimeError(('Acceptance failed',lam,alpha,val,dist))
        if abs(alpha-previous[-1][1])>mp.mpf('.0002'):raise RuntimeError('Continuation jump')
        rh,rc=geometry(lam);sf=lam*(rc-rh)/2
        item={'Lambda_tilde':str(lam),'alpha':str(alpha),'varpi_im':str(-alpha/sf),
              'rh_omega_im':str(-rh*alpha),'normalized_W_higher_order':str(val),
              'minimum_integer_index_distance':str(dist),'prediction_error':str(abs(alpha-guess)),
              'step_frequency_change':str(abs(alpha-previous[-1][1])),
              'dps':dps,'series_order':N,'check_series_order':Ncheck,'geometry':detail}
        previous.append((lam,alpha));report['path'].append(item)
        if j in [19,39]:
            checks=[]
            for dd,nn,where in [(dps,N,None),(dps+20,Ncheck,None),(dps+20,Ncheck,'0.45')]:
                with mp.workdps(dd):
                    rr=locate_local(lam,alpha,nn,delta='1e-12',match_fraction=where)
                    ww,details=connection(lam,rr,nn,details=True,match_fraction=where)
                    slope=mp.diff(lambda a:connection(lam,a,nn,match_fraction=where),rr)
                    checks.append({'dps':dd,'N':nn,'matching_t':details['match_t'],
                                  'alpha':str(rr),'difference_path':str(abs(rr-alpha)),
                                  'normalized_W':str(ww),'d_normalized_W_dalpha':str(slope),
                                  'geometry_and_series':details})
            selected=dict(item);selected['validation']=checks
            with mp.workdps(dps+40):
                checked_roots=[mp.mpf(q['alpha']) for q in checks]
                selected['maximum_alpha_difference']=str(max(checked_roots)-min(checked_roots))
            report['selected'].append(selected)
            print('selected',mp.nstr(lam,8),mp.nstr(alpha,25),'scaled',mp.nstr(alpha/sf,18),'elapsed',round(time.time()-start,1),flush=True)
        elif j%5==0:print('extended',j,mp.nstr(lam,8),mp.nstr(alpha,15),'elapsed',round(time.time()-start,1),flush=True)
        (HERE/'ds_extended_results.json').write_text(json.dumps(report,indent=2))
report['extended_complete']=True
report['extension_elapsed_seconds']=time.time()-start
(HERE/'ds_extended_results.json').write_text(json.dumps(report,indent=2))
with (HERE/'ds_selected.csv').open('w') as f:
    out=csv.writer(f)
    out.writerow(['Lambda_tilde','minus_Im_Omega','minus_Im_varpi','minus_Im_rh_omega','max_alpha_variation','minimum_integer_index_distance'])
    for p in report['selected']:
        out.writerow([mp.nstr(mp.mpf(p['Lambda_tilde']),10),mp.nstr(mp.mpf(p['alpha']),24),
                      mp.nstr(-mp.mpf(p['varpi_im']),20),mp.nstr(-mp.mpf(p['rh_omega_im']),20),
                      p['maximum_alpha_difference'],p['minimum_integer_index_distance']])
print('COMPLETE',len(report['path']),'physical continuation points',len(report['selected']),'validation points',flush=True)
