"""Summarize and structurally verify saved N240/N300 runs; no new eigensolve."""
from pathlib import Path
import csv,hashlib,json
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
KEYS=['N','dps','varpi_im','Omega_im','error_to_independent_varpi_reference',
 'error_to_independent_Omega_reference','relative_pencil_residual',
 'independent_original_ODE_residual','relative_reduced_ODE_residual_including_endpoints',
 'relative_reduced_endpoint_residual','last_eight_coefficient_norm_fraction',
 'eigenfunction_error_to_independent_Frobenius_Linf']
old=json.loads((ROOT/'checks/r2p5/imaginary_multiprecision_results.json').read_text())
seed=next(r for r in old['results'] if r['N']==128 and r['dps']==70)
rows=[{**{k:seed[k] for k in KEYS},'seed_type':'preceding binary64 spectral seed',
       'Newton_corrections':len(seed['iterations'])}]
checks={}
for N in [240,300]:
 report=json.loads((OUT/f'high_order_multiprecision_N{N}.json').read_text())
 r=report['result']
 assertions={
  'order_precision':r['N']==N and r['dps']==70,
  'complete_coefficient_vector':len(r['coefficients'])==N,
  'complete_1025_point_field':len(r['phi'])==1025,
  'actual_Newton_correction':report['actual_Newton_corrections']==len(r['iterations'])>=1,
  'reference_loaded_after_Newton':report['reference']['loaded_at_elapsed_seconds']>r['iterations'][-1]['elapsed_seconds'],
  'source_seed_frequency':report['seed']['source_varpi_im']==seed['varpi_im'],
  'exact_padding_length':report['seed']['zero_coefficients_added']==N-128,
  'original_solver_unchanged':report['unchanged_solver']['sha256']==hashlib.sha256((ROOT/report['unchanged_solver']['path']).read_bytes()).hexdigest(),
  'original_seed_file_unchanged':report['seed']['source_sha256']==hashlib.sha256((ROOT/report['seed']['source']).read_bytes()).hexdigest(),
  'independent_reference_unchanged':report['reference']['sha256']==hashlib.sha256((ROOT/report['reference']['path']).read_bytes()).hexdigest(),
 }
 if not all(assertions.values()):raise AssertionError(assertions)
 checks[str(N)]=assertions
 rows.append({**{k:r[k] for k in KEYS},'seed_type':'zero-padded N128 70-digit spectral eigenpair',
              'Newton_corrections':len(r['iterations'])})
 print(f'N={N}: {len(assertions)} checks pass; varpi difference={float(r["error_to_independent_varpi_reference"]):.8e}, field difference={float(r["eigenfunction_error_to_independent_Frobenius_Linf"]):.8e}')
summary={'scope':'New fully rebuilt 70-digit N240/N300 targeted imaginary spectral QEP checks versus an existing independent reference; N128 is the preceding spectral result for comparison. No old result modified.',
 'rows':rows,'checks':checks,'all_checks_passed':True,
 'limits':'Numerical diagnostics only. Fixed70-digit higher-order differences increase from the N128 minimum; this is not monotone truncation convergence or a rigorous accuracy bound.'}
(OUT/'high_order_multiprecision_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
with (OUT/'high_order_multiprecision_summary.csv').open('w',newline='') as stream:
 writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
