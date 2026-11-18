#!/usr/bin/env python3
"""Bounded extension of the three R2p5 spectral-resolution tracks to N=300.

Imports the unchanged R2p5 assembly, window, and diagnostic functions. Each
track starts at its archived N128 spectral candidate and uses nearest-previous
association in the original fixed window. Independent references enter only
after candidate selection. Failed associations are recorded and do not advance
the track. This is a new binary64 experiment, not a historical Maple/Julia rerun.
"""
from pathlib import Path
import csv
import hashlib
import json
import platform
import sys
import time

import numpy as np
import scipy
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
PRIOR = HERE.parent / 'r2p5'
sys.path.insert(0, str(PRIOR))
import broad_resolution_study as base

ORDERS = [160, 190, 195, 200, 240, 300]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def raw_frequency(z):
    return {'re': float(z.real) if np.isfinite(z.real) else str(z.real),
            'im': float(z.imag) if np.isfinite(z.imag) else str(z.imag),
            'finite': bool(np.isfinite(z))}


def spread(rows, orders, scale):
    selected = [r for r in rows if r['N'] in orders and r['status'] == 'associated']
    result = {'requested_N': orders, 'available_associated_N': [r['N'] for r in selected],
              'complete': len(selected) == len(orders)}
    if len(selected) < 2:
        return {**result, 'varpi': None, 'Omega': None, 'maximizing_pair_N': None}
    best = max((abs(complex(a['varpi_re'], a['varpi_im']) -
                        complex(b['varpi_re'], b['varpi_im'])), a['N'], b['N'])
               for j, a in enumerate(selected) for b in selected[j + 1:])
    return {**result, 'varpi': best[0], 'Omega': best[0] * scale,
            'maximizing_pair_N': list(best[1:])}


def main():
    start = time.monotonic()
    original = json.loads((PRIOR / 'broad_resolution_results.json').read_text())
    imaginary_field = json.loads((PRIOR / 'ds_reference_eigenfunction.json').read_text())
    oscillatory_fields = json.loads((PRIOR / 'offaxis_reference_eigenfunctions.json').read_text())
    frequency_references = base.references()
    dependencies = [PRIOR / 'broad_resolution_study.py',
                    PRIOR / 'broad_resolution_results.json',
                    PRIOR / 'ds_reference_eigenfunction.json',
                    PRIOR / 'offaxis_reference_eigenfunctions.json',
                    HERE.parent / 'r2p4' / 'continue_near_nariai.py',
                    HERE.parent / 'r2p4' / 'jost_extra_cell_check.json',
                    HERE.parent / 'r2p4' / 'jost_check_results.json',
                    HERE.parent / 'r2p4' / 'ds_extended_results.json']
    report = {
        'scope': __doc__, 'precision': 'IEEE binary64', 'orders': ORDERS,
        'versions': {'python': platform.python_version(), 'numpy': np.__version__,
                     'scipy': scipy.__version__},
        'dependency_sha256': {str(p.relative_to(HERE.parent.parent)): digest(p) for p in dependencies},
        'selection': 'Start from each archived N128 spectral candidate, then nearest previously accepted spectral candidate in unchanged fixed windows. Require displacement <0.1 and, if a competitor exists, displacement <0.25 of its separation. No independent reference enters selection. Failed association is recorded and does not advance the previous accepted root.',
        'windows': original['windows'], 'grid_y': base.Y.tolist(),
        'normalization': original['normalization'], 'residual': original['residual'],
        'limitations': 'Only three tracked representatives. No claim about the full spectrum, historical matrix provenance, completeness, monotone convergence, or radiation admissibility from residuals/association alone. Reference frequency and field variations retain the limitations of the R2p5 study.',
        'cases': []}
    compact = []
    with threadpool_limits(limits=1):
        for case in base.CASES:
            oldcase = next(c for c in original['cases'] if c['definition']['name'] == case['name'])
            seed = next(r for r in oldcase['results'] if r['N'] == 128)
            previous = complex(seed['varpi_re'], seed['varpi_im'])
            previous_N = 128
            previous_phi = np.asarray(seed['phi_re']) + 1j * np.asarray(seed['phi_im'])
            scale = base.geometry(case['Lambda_tilde'])['scale']
            rows = []
            reference = frequency_references[case['name']]
            reference_w = complex(reference['re'], reference['im'])
            if case['kind'] == 'imaginary':
                reference_phi = np.asarray([float(p['phi']) for p in imaginary_field['points']])
                field_variation = float(imaginary_field['max_relative_sup_difference'])
            else:
                field_case = next(c for c in oscillatory_fields['cases'] if c['definition']['name'] == case['name'])
                reference_phi = np.asarray(field_case['phi_re']) + 1j * np.asarray(field_case['phi_im'])
                field_variation = field_case['reference_configuration_difference_relative_Linf']
            for N in ORDERS:
                t0 = time.monotonic()
                vals, vectors, matrices = base.all_eigenpairs(N, case)
                indices = base.candidate_indices(vals, case)
                row = {'case': case['name'], 'N': N, 'Lambda_tilde': case['Lambda_tilde'],
                       'previous_accepted_N': previous_N, 'raw_count': len(vals),
                       'finite_raw_count': int(np.count_nonzero(np.isfinite(vals))),
                       'window_candidate_count': len(indices),
                       'raw_roots': [raw_frequency(v) for v in vals],
                       'window_indices': indices.tolist()}
                if not len(indices):
                    row.update(status='no_candidate_in_fixed_window', elapsed_seconds=time.monotonic() - t0)
                else:
                    chosen = int(indices[np.argmin(abs(vals[indices] - previous))])
                    w = vals[chosen]
                    competitors = vals[indices[indices != chosen]]
                    gap = float(np.min(abs(w - competitors))) if len(competitors) else None
                    displacement = float(abs(w - previous))
                    accepted = displacement < .1 and (gap is None or displacement < .25 * gap)
                    coeff, phi, metrics = base.eigenpair_diagnostics(w, vectors[:, chosen], matrices, case)
                    row.update(
                        status='associated' if accepted else 'association_failed', selected_raw_index=chosen,
                        varpi_re=float(w.real), varpi_im=float(w.imag),
                        Omega_re=float((w * scale).real), Omega_im=float((w * scale).imag),
                        change_from_previous_accepted_N=displacement,
                        nearest_candidate_separation_in_window=gap,
                        association_step_pass=bool(displacement < .1),
                        association_gap_pass=bool(gap is None or displacement < .25 * gap),
                        error_to_independent_varpi_reference=float(abs(w - reference_w)),
                        error_to_independent_Omega_reference=float(abs(w - reference_w) * scale),
                        eigenfunction_change_from_previous_accepted_N_Linf=float(np.max(abs(phi - previous_phi)) / np.max(abs(previous_phi))),
                        eigenfunction_error_to_independent_Frobenius_Linf=float(np.max(abs(phi - reference_phi)) / np.max(abs(reference_phi))),
                        independent_eigenfunction_configuration_difference_relative_Linf=field_variation,
                        coeff_re=coeff.real.tolist(), coeff_im=coeff.imag.tolist(),
                        phi_re=phi.real.tolist(), phi_im=phi.imag.tolist(), **metrics,
                        elapsed_seconds=time.monotonic() - t0)
                    if accepted:
                        previous, previous_phi, previous_N = w, phi, N
                rows.append(row)
                compact.append({k: v for k, v in row.items() if not isinstance(v, (list, dict))})
                print(json.dumps(compact[-1], allow_nan=False), flush=True)
            spreads = {
                'N190_195_200': spread(rows, [190, 195, 200], scale),
                'N190_to300': spread(rows, [190, 195, 200, 240, 300], scale),
                'N240_300': spread(rows, [240, 300], scale),
                'N128_to300': spread([seed | {'status': 'associated'}] + rows, [128] + ORDERS, scale)}
            report['cases'].append({'definition': case, 'N128_spectral_seed': {
                'varpi_re': seed['varpi_re'], 'varpi_im': seed['varpi_im'],
                'source': 'checks/r2p5/broad_resolution_results.json'},
                'independent_frequency_reference': reference, 'results': rows, 'spreads': spreads})
            print('SPREADS', case['name'], json.dumps(spreads), flush=True)
            report['elapsed_seconds'] = time.monotonic() - start
            (HERE / 'high_order_binary64.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    fields = list(dict.fromkeys(k for row in compact for k in row))
    with (HERE / 'high_order_binary64.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(compact)
    print('COMPLETE', report['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
