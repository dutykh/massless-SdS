#!/usr/bin/env python3
"""Clarify stored absolute/relative norms without rerunning any eigensolve.

The original reference-configuration differences are absolute. Explicit
relative fields are calculated by division by the stored reference sup norm.
Existing broad-study frequencies and relative eigenfunction errors are checked
for exact preservation. This script is idempotent.
"""
from pathlib import Path
import csv,json
import numpy as np

HERE=Path(__file__).resolve().parent
offpath=HERE/'offaxis_reference_eigenfunctions.json'
broadpath=HERE/'broad_resolution_results.json'
off=json.loads(offpath.read_text());broad=json.loads(broadpath.read_text())
ds=json.loads((HERE/'ds_reference_eigenfunction.json').read_text())
preserved={(c['definition']['name'],r['N']):tuple(r[k] for k in
 ('varpi_re','varpi_im','Omega_re','Omega_im','eigenfunction_error_to_independent_Frobenius_Linf'))
 for c in broad['cases'] for r in c['results']}

for c in off['cases']:
 reference=np.array(c['phi_re'])+1j*np.array(c['phi_im'])
 norm=float(np.max(np.abs(reference)))
 absolute=c.get('reference_configuration_difference_absolute_Linf',c['reference_configuration_difference_Linf'])
 c['reference_function_Linf_norm']=norm
 c['reference_configuration_difference_absolute_Linf']=absolute
 c['reference_configuration_difference_relative_Linf']=absolute/norm
 original=next(b for b in broad['cases'] if b['definition']['name']==c['definition']['name'])
 for row in c['spectral_comparison']:
  numeric=next(r for r in original['results'] if r['N']==row['N'])
  phi=np.array(numeric['phi_re'])+1j*np.array(numeric['phi_im'])
  delta=float(np.max(np.abs(phi-reference)))
  row['eigenfunction_difference_absolute_Linf']=delta
  row['eigenfunction_difference_relative_Linf']=delta/norm
  assert row['eigenfunction_difference_relative_Linf']==numeric['eigenfunction_error_to_independent_Frobenius_Linf']
off['norm_definitions']='Absolute Linf is max_j abs(Phi_a-Phi_b). Relative Linf divides this by max_j abs(Phi_b), configuration_b. The legacy fields reference_configuration_difference_Linf and nested eigenfunction_error_to_independent_Frobenius_Linf retain ABSOLUTE values; use the explicit absolute/relative fields.'

for c in broad['cases']:
 if c['definition']['kind']=='oscillatory':
  ref=next(r for r in off['cases'] if r['definition']['name']==c['definition']['name'])
  absolute=ref['reference_configuration_difference_absolute_Linf']
  relative=ref['reference_configuration_difference_relative_Linf']
 else:
  absolute=float(ds['max_absolute_difference']);relative=float(ds['max_relative_sup_difference'])
 for row in c['results']:
  row['independent_eigenfunction_configuration_difference_Linf']=relative
  row['independent_eigenfunction_configuration_difference_absolute_Linf']=absolute
  row['independent_eigenfunction_configuration_difference_relative_Linf']=relative
broad['reference_variation_norms']='Explicit absolute and relative reference-configuration differences are recorded. Relative differences divide by the reference function maximum norm. The legacy independent_eigenfunction_configuration_difference_Linf field is an alias of the RELATIVE value.'
for c in broad['cases']:
 for row in c['results']:
  assert tuple(row[k] for k in ('varpi_re','varpi_im','Omega_re','Omega_im','eigenfunction_error_to_independent_Frobenius_Linf'))==preserved[c['definition']['name'],row['N']]
offpath.write_text(json.dumps(off,indent=2)+'\n')
broadpath.write_text(json.dumps(broad,indent=2)+'\n')
compact=[{k:v for k,v in r.items() if not isinstance(v,list)} for c in broad['cases'] for r in c['results']]
with (HERE/'broad_resolution_results.csv').open('w',newline='') as stream:
 writer=csv.DictWriter(stream,fieldnames=list(compact[0]));writer.writeheader();writer.writerows(compact)
for c in off['cases']:
 print(c['definition']['name'],'absolute',c['reference_configuration_difference_absolute_Linf'],'relative',c['reference_configuration_difference_relative_Linf'])
print('Preserved all 24 broad-study frequencies and relative eigenfunction errors exactly; no eigensolve rerun.')
