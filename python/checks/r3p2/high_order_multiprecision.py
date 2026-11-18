#!/usr/bin/env python3
"""Targeted N240/N300,70-digit imaginary spectral QEP extension.

Imports the unchanged Referee2 Newton implementation. The sole eigenpair
seed is the stored converged N128,70-digit spectral eigenpair, padded with
zero Chebyshev coefficients. No independent Jost/Frobenius frequency or
field enters the matrix construction, initialization or Newton iteration.
The reference is loaded lazily when the unchanged post-Newton diagnostics
first access it. The solver always computes a bordered Newton correction
before checking its stopping criterion. This is a new bounded calculation,
not a rerun of an original historical Maple/Julia campaign.
"""
from pathlib import Path
import argparse,hashlib,json,platform,sys,time
import mpmath as mp

ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
PREVIOUS=ROOT/'checks/r2p5'
sys.path.insert(0,str(PREVIOUS))
import imaginary_multiprecision_newton as implementation

class PostNewtonReference:
    """Loaded only at the existing diagnostics() function's first access."""
    def __init__(self,path):
        self.path=path
        self.data=None
        self.loaded_at_elapsed_seconds=None
        self.started=time.monotonic()
    def __getitem__(self,key):
        if self.data is None:
            self.data=json.loads(self.path.read_text())
            self.loaded_at_elapsed_seconds=time.monotonic()-self.started
            print('Independent reference first loaded in post-Newton diagnostics',
                  f'elapsed={self.loaded_at_elapsed_seconds:.3f}s',flush=True)
        return self.data[key]

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--order',type=int,choices=[240,300],required=True)
    parser.add_argument('--dps',type=int,default=70)
    args=parser.parse_args()
    if args.dps!=70:raise RuntimeError('This targeted reproducibility run is specified at70decimaldigits')
    started=time.monotonic();mp.mp.dps=args.dps
    seed_path=PREVIOUS/'imaginary_multiprecision_results.json'
    stored=json.loads(seed_path.read_text())
    source=next(r for r in stored['results'] if r['N']==128 and r['dps']==70)
    if len(source['coefficients'])!=128:raise RuntimeError('Unexpected seed coefficient count')
    # All entries are constructed directly as70-digit mpf values. No binary64
    # promotion or independent physical reference is involved in this seed.
    coefficients=[mp.mpf(c) for c in source['coefficients']]
    seed={'varpi_im':mp.mpf(source['varpi_im']),
          'coeff_re':coefficients+[mp.mpf(0)]*(args.order-128),
          'coeff_im':[mp.mpf(0)]*args.order}
    reference_path=PREVIOUS/'ds_reference_eigenfunction.json'
    reference=PostNewtonReference(reference_path)
    method=('Unchanged imaginary_multiprecision_newton.solve(); fully rebuilt '
            '70-digit real QEP; bordered Newton for alpha=-Im(varpi) and all '
            'N Chebyshev coefficients, normalized by Phi(-0.3)=1.')
    seed_description=('Converged N128,70-digit SPECTRAL eigenpair from the '
                      'preceding precision study, converted directly to70-digit '
                      'mpf arrays and zero-padded to the new N. No Jost or '
                      'independent reference value in the solve.')
    print('START',json.dumps({'N':args.order,'dps':args.dps,'Lambda_tilde':'.06',
       'spin':0,'ell':2,'seed_N':128,'seed_dps':70,'zero_padding':args.order-128,
       'seed_varpi_im':str(seed['varpi_im']),'reference_loaded_before_solve':reference.data is not None,
       'method':method}),flush=True)
    solver_hash_before=digest(Path(implementation.__file__))
    result=implementation.solve(args.order,args.dps,seed,reference)
    if not result['iterations']:raise RuntimeError('No actual Newton correction recorded')
    if digest(Path(implementation.__file__))!=solver_hash_before:raise RuntimeError('Original solver changed during the run')
    inherited_label=result['seed']
    result['seed']=seed_description
    result['frequency_change_from_N128_spectral_seed']=str(abs(mp.mpf(result['varpi_im'])-mp.mpf(source['varpi_im'])))
    report={'scope':__doc__,'method':method,'runtime':{'python':platform.python_version(),
      'mpmath':mp.__version__,'decimal_digits':args.dps},'sector':{'spin':0,'ell':2,'Lambda_tilde':'.06'},
      'seed':{'N':128,'dps':70,'source':str(seed_path.relative_to(ROOT)),
       'source_sha256':digest(seed_path),'source_varpi_im':source['varpi_im'],
       'zero_coefficients_added':args.order-128,'coeff_im':'all exactly zero',
       'description':seed_description,'inherited_solver_seed_label_overridden':inherited_label},
      'unchanged_solver':{'path':str(Path(implementation.__file__).relative_to(ROOT)),
                         'sha256':solver_hash_before},
      'reference':{'path':str(reference_path.relative_to(ROOT)),
        'sha256':digest(reference_path),'role':'frequency and1025-point field comparison only after Newton convergence',
        'loaded_at_elapsed_seconds':reference.loaded_at_elapsed_seconds},
      'offgrid_test_points':1025,'actual_Newton_corrections':len(result['iterations']),
      'post_Newton_total_elapsed_seconds':time.monotonic()-started,
      'result':result,'limits':'Numerical residuals and cross-method differences are diagnostics, not rigorous error bounds or a completeness proof.'}
    path=OUT/f'high_order_multiprecision_N{args.order}.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    summary={k:result[k] for k in ['N','dps','error_to_independent_varpi_reference',
      'error_to_independent_Omega_reference','relative_pencil_residual',
      'independent_original_ODE_residual','relative_reduced_ODE_residual_including_endpoints',
      'relative_reduced_endpoint_residual','last_eight_coefficient_norm_fraction',
      'eigenfunction_error_to_independent_Frobenius_Linf',
      'frequency_change_from_N128_spectral_seed']}
    summary.update(actual_Newton_corrections=len(result['iterations']),
                   total_elapsed_seconds=report['post_Newton_total_elapsed_seconds'])
    print('FINAL',json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
