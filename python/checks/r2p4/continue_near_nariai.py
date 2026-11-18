#!/usr/bin/env python3
"""Independent binary64 continuation of six low SdS photon-sphere branches.

This does not run the supplied Maple/Julia calculation.  It collocates a radial
equation derived directly from the SdS metric, using varpi=Omega/[lambda*(xc-xh)].
All six branches are followed in small parameter steps, independently at each
resolution.  No manuscript frequency table is read, fitted, or used as a seed.
Run: OPENBLAS_NUM_THREADS=1 python checks/r2p4/continue_near_nariai.py
"""
from __future__ import annotations

import csv
import itertools
import json
import platform
from functools import lru_cache
from pathlib import Path

import numpy as np
import scipy
from numpy.polynomial.chebyshev import chebder, chebval, chebvander
from scipy.linalg import eig
from scipy.optimize import linear_sum_assignment
from threadpoolctl import threadpool_limits

OUT = Path(__file__).resolve().parent
ANCHORS = np.array([.09, .10, .105, .11, .111, .1111, .11111])
STEPS_PER_INTERVAL = 12
GRID = np.concatenate([
    np.linspace(a, b, STEPS_PER_INTERVAL, endpoint=False)
    for a, b in zip(ANCHORS[:-1], ANCHORS[1:])
] + [ANCHORS[-1:]])
ORDERS = (40, 60, 80, 100)
SECTORS = ('scalar', 'em', 'odd')
ELL = 2
MAX_STEP = .08
CONVERGENCE_TOLERANCE = 1e-7
OFFGRID_TOLERANCE = 1e-7
PENCIL_TOLERANCE = 1e-12
TAIL_TOLERANCE = 1e-7


def geometry(lam):
    # Symmetric trigonometric roots; d is evaluated without subtracting roots.
    theta = np.arccos(3*np.sqrt(lam))/3
    a = np.cos(theta)/(2*np.sqrt(lam))
    d = np.sqrt(3)*np.sin(theta)/np.sqrt(lam)
    return {'a': a, 'd': d, 'xh': a-d/2, 'xc': a+d/2,
            'x0': -2*a, 'scale': lam*d}


def coefficients(y, lam, sector):
    """Returns g,g',G,G',R,Vbar in the regular residual equation.

    f=lambda*d^2*g, g=(1-y^2)F, F=(x-x0)/(3x), x=a+d*y/2.
    Psi=(1+y)^(-i*w*ph)*(1-y)^(-i*w*pc)*Phi,
    ph=3*xh/[2*(xh-x0)], pc=3*xc/[2*(xc-x0)].
    Let G=g*d(log(Psi/Phi))/dy/(i*w). Exact cancellation gives
        G=y+C*(1-y^2)/x, C=3*a^2*d/[(2*xh+xc)*(xh+2*xc)],
        R=(1-G^2)/g=[1-2*y*C/x-(C/x)^2*(1-y^2)]/F.
    The nonsingular residual equation is
      g*Phi''+(g'+2*i*w*G)*Phi'+(-Vbar+i*w*G'+w^2*R)*Phi=0.
    """
    bg = geometry(lam)
    a, d, xh, xc, x0 = (bg[k] for k in ('a','d','xh','xc','x0'))
    x = a+d*y/2
    F = (x-x0)/(3*x)
    Fy = x0*d/(6*x*x)
    g = (1-y*y)*F
    gy = -2*y*F+(1-y*y)*Fy
    C = 3*a*a*d/((2*xh+xc)*(xh+2*xc))
    G = y+C*(1-y*y)/x
    Gy = 1+C*(-2*y/x-d*(1-y*y)/(2*x*x))
    R = (1-2*y*C/x-(C/x)**2*(1-y*y))/F
    Vbar = ELL*(ELL+1)/x**2
    if sector == 'scalar':
        Vbar = Vbar+2*lam*d*gy/x
    elif sector == 'odd':
        Vbar = Vbar-3/x**3
    elif sector != 'em':
        raise ValueError(sector)
    Vbar = Vbar/(4*lam)
    return g, gy, G, Gy, R, Vbar


@lru_cache(None)
def basis(order):
    y = np.cos(np.pi*(np.arange(order)+.5)/order)
    T = chebvander(y, order-1)
    T1, T2 = np.zeros_like(T), np.zeros_like(T)
    for j in range(1, order):
        c = np.zeros(order); c[j] = 1
        T1[:,j] = chebval(y, chebder(c))
        T2[:,j] = chebval(y, chebder(c, 2))
    return y, T, T1, T2


def pencil(order, lam, sector):
    y, T, T1, T2 = basis(order)
    g, gy, G, Gy, R, Vbar = coefficients(y, lam, sector)
    A0 = g[:,None]*T2+gy[:,None]*T1-Vbar[:,None]*T
    A1 = 1j*(2*G[:,None]*T1+Gy[:,None]*T)
    A2 = R[:,None]*T
    # Identical nonzero scaling of every polynomial coefficient matrix.
    row = np.max(np.abs(np.hstack([A0,A1,A2])), axis=1)
    A0 /= row[:,None]; A1 /= row[:,None]; A2 /= row[:,None]
    I, Z = np.eye(order), np.zeros((order,order))
    vals, vecs = eig(np.block([[Z,I],[-A0,-A1]]),
                     np.block([[I,Z],[Z,A2]]))
    window = (np.isfinite(vals) & (vals.real > .5) & (vals.real < 4)
              & (vals.imag < -.05) & (vals.imag > -3.2))
    idx = np.flatnonzero(window)
    return vals[idx], vecs[:order,idx], (A0,A1,A2)


def diagnostics(w, c, matrices, lam, sector):
    c = c/np.max(np.abs(c))
    A0,A1,A2 = matrices
    pencil_residual = np.linalg.norm((A0+w*A1+w*w*A2)@c)/(
        np.linalg.norm(c)*(np.linalg.norm(A0)+abs(w)*np.linalg.norm(A1)
                           +abs(w*w)*np.linalg.norm(A2)))
    # A separate, overresolved grid is never used in the eigenvalue solve.
    # It includes both endpoints and tests the limiting first-order equations.
    yy = np.cos(np.pi*np.arange(8*len(c)+1)/(8*len(c)))
    p = chebval(yy,c); p1 = chebval(yy,chebder(c)); p2 = chebval(yy,chebder(c,2))
    g,gy,G,Gy,R,Vbar = coefficients(yy,lam,sector)
    terms = np.array([g*p2, (gy+2j*w*G)*p1, (-Vbar+1j*w*Gy+w*w*R)*p])
    balance = np.max(np.sum(np.abs(terms), axis=0))
    residual = np.sum(terms,axis=0)
    offgrid = np.max(np.abs(residual))/balance
    endpoint = np.max(np.abs(residual[[0,-1]]))/balance
    # Highest eight coefficients, normalized to the complete coefficient norm.
    tail = np.linalg.norm(c[-8:])/np.linalg.norm(c)
    return {'relative_pencil_residual':float(pencil_residual),
            'relative_offgrid_residual':float(offgrid),
            'relative_endpoint_residual':float(endpoint),
            'coefficient_tail_fraction':float(tail)}


def track(order, sector):
    rows = []
    prev = None
    for j, lam in enumerate(GRID):
        vals, vectors, matrices = pencil(order,float(lam),sector)
        if len(vals) < 2:
            raise RuntimeError(f'Too few candidates: {order},{sector},{lam}')
        if j == 0:
            # Independently initialize from the two least damped roots in the
            # stated window; do not seed with the analytical Nariai spectrum.
            selected = np.argsort(-vals.imag)[:2]
            predictor = vals[selected]
        else:
            predictor = prev.copy()
            if j > 1:
                predictor += (prev-prevprev)*(lam-GRID[j-1])/(GRID[j-1]-GRID[j-2])
            cost = np.abs(predictor[:,None]-vals[None,:])
            rr, selected = linear_sum_assignment(cost)
            assert np.array_equal(rr,np.arange(2))
        matched = vals[selected]
        for n,k in enumerate(selected):
            w = vals[k]
            step = 0. if j == 0 else float(abs(w-prev[n]))
            prediction_error = float(abs(w-predictor[n]))
            other = np.delete(vals,k)
            candidate_separation = float(np.min(np.abs(w-other)))
            alternative_distance = float(np.min(np.abs(predictor[n]-other)))
            unambiguous = prediction_error < .25*alternative_distance
            if step >= MAX_STEP or not unambiguous:
                raise RuntimeError(f'Ambiguous continuation: {order},{sector},{lam},n{n}')
            row = {'sector':sector,'ell':ELL,'n':n,'N':order,'step_index':j,
                   'Lambda_tilde':float(lam),'varpi_re':float(w.real),
                   'varpi_im':float(w.imag),'step_distance':step,
                   'predictor_error':prediction_error,
                   'nearest_candidate_separation':candidate_separation,
                   'alternative_predictor_distance':alternative_distance,
                   'candidate_count_in_window':len(vals),
                   **diagnostics(w,vectors[:,k],matrices,float(lam),sector)}
            rows.append(row)
        if j:
            prevprev = prev.copy()
        prev = matched.copy()
    print(f'{sector}, N={order}: {len(GRID)} continuation points; '
          f'largest step {max(r["step_distance"] for r in rows):.3g}',flush=True)
    return rows


def write_csv(path,rows):
    with path.open('w',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def main():
    all_rows = []
    with threadpool_limits(limits=1):
        for sector,order in itertools.product(SECTORS,ORDERS):
            all_rows.extend(track(order,sector))
    lookup = {(r['sector'],r['N'],r['step_index'],r['n']):r for r in all_rows}
    summaries = []
    for sector,j,n in itertools.product(SECTORS,range(len(GRID)),range(2)):
        rows = [lookup[sector,N,j,n] for N in ORDERS]
        ws = [complex(r['varpi_re'],r['varpi_im']) for r in rows]
        spread = max(abs(a-b) for a,b in itertools.combinations(ws,2))
        worst_offgrid = max(r['relative_offgrid_residual'] for r in rows)
        worst_pencil = max(r['relative_pencil_residual'] for r in rows)
        worst_tail = max(r['coefficient_tail_fraction'] for r in rows)
        passed = (spread < CONVERGENCE_TOLERANCE and
                  worst_offgrid < OFFGRID_TOLERANCE and
                  worst_pencil < PENCIL_TOLERANCE and worst_tail < TAIL_TOLERANCE)
        w = ws[-1]
        A = ELL*(ELL+1)-(2 if sector == 'odd' else 0)
        exact = complex(np.sqrt(A-.25),-(n+.5))
        bg = geometry(float(GRID[j]))
        summary = {'sector':sector,'ell':ELL,'n':n,'step_index':j,
                   'Lambda_tilde':float(GRID[j]),'xh':float(bg['xh']),
                   'xc':float(bg['xc']),'frequency_scale':float(bg['scale']),
                   'varpi_re':float(w.real),'varpi_im':float(w.imag),
                   'Omega_re':float((bg['scale']*w).real),
                   'Omega_im':float((bg['scale']*w).imag),
                   'Nariai_re':float(exact.real),'Nariai_im':float(exact.imag),
                   'distance_to_Nariai':float(abs(w-exact)),
                   'max_resolution_spread':float(spread),
                   'max_offgrid_residual':worst_offgrid,
                   'max_endpoint_residual':max(r['relative_endpoint_residual'] for r in rows),
                   'max_pencil_residual':worst_pencil,'max_tail_fraction':worst_tail,
                   'passed_acceptance':passed}
        for N,wn in zip(ORDERS,ws):
            summary[f'N{N}_re'] = float(wn.real); summary[f'N{N}_im'] = float(wn.imag)
        summaries.append(summary)
    anchor_indices = set(STEPS_PER_INTERVAL*k for k in range(len(ANCHORS)))
    anchors = [r for r in summaries if r['step_index'] in anchor_indices]
    failed = [r for r in summaries if not r['passed_acceptance']]
    report = {
        'description':__doc__,'precision':'IEEE 754 binary64; independent SciPy diagnostic',
        'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__,
        'frequency_definition':'varpi=Omega/[Lambda_tilde*(xc-xh)], Omega=M*omega, x=r/(2M)',
        'method':'Chebyshev polynomial coefficient collocation on interior roots; regular equation derived directly from metric; common row normalization; generalized companion pencil',
        'continuation':'Each N tracks separately. Start at .09 with two least damped positive-real roots in stated window. Subsequent roots assigned jointly to linear predictors by minimum total distance. Require step<.08 and predictor error<one quarter of distance to any alternative candidate.',
        'window':{'real':[.5,4.],'imag':[-3.2,-.05]},
        'limits':'Only six low positive-real branches are tested. This is not a reproduction of the supplied Maple/Julia run, a complete spectrum, or a test of purely imaginary mode counts. Off-grid and pencil residuals are numerical diagnostics; independent two-horizon Jost checks are separate.',
        'anchors':ANCHORS.tolist(),'steps_per_interval':STEPS_PER_INTERVAL,
        'grid_points':len(GRID),'resolutions':list(ORDERS),
        'acceptance':{'max_resolution_spread':CONVERGENCE_TOLERANCE,
                      'max_relative_offgrid_residual':OFFGRID_TOLERANCE,
                      'max_relative_pencil_residual':PENCIL_TOLERANCE,
                      'max_tail_fraction':TAIL_TOLERANCE},
        'passed':len(summaries)-len(failed),'total':len(summaries),
        'maxima':{key:max(r[key] for r in summaries) for key in
                  ['max_resolution_spread','max_offgrid_residual',
                   'max_endpoint_residual','max_pencil_residual','max_tail_fraction']},
        'max_continuation_step':max(r['step_distance'] for r in all_rows),
        'minimum_candidate_separation':min(r['nearest_candidate_separation'] for r in all_rows),
        'anchor_results':anchors,'all_summary_results':summaries,'failed_results':failed}
    write_csv(OUT/'near_nariai_continuation_all_resolutions.csv',all_rows)
    write_csv(OUT/'near_nariai_continuation_summary.csv',summaries)
    write_csv(OUT/'near_nariai_anchor_results.csv',anchors)
    (OUT/'near_nariai_continuation_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:report[k] for k in ['passed','total','maxima','max_continuation_step','minimum_candidate_separation']},indent=2))
    if failed:
        raise SystemExit('One or more diagnostic acceptance checks failed; inspect JSON.')


if __name__ == '__main__':
    main()
