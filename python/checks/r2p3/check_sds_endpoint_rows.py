#!/usr/bin/env python3
"""Exact source-level checks for the SdS collocation endpoint rows.

Reads decoded active Maple assembler inputs without executing Maple or Julia.
Uses exact horizon relations, parameterized by zc>1, to verify that the
unsimplified endpoint rows vanish and that division by 1-y^2 gives finite,
nontrivial first-order limiting rows. This is not a spectral recomputation.
"""
from pathlib import Path
import hashlib
import json
import re
import sympy as sp

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent
names = 'y zc z0 xh Lambda alphac alphah epsilon L c'
y, zc, z0, xh, lam, ac, ah, epsilon, ell, c = sp.symbols(names)
namespace = dict(zip(names.split(), (y, zc, z0, xh, lam, ac, ah, epsilon, ell, c)))
h = (zc**2 + zc + 1)/(zc*(zc + 1))
cosmological = 3/(4*h**2*(zc**2 + zc + 1))
substitutions = {
    xh: h, lam: cosmological, z0: -(1 + zc),
    ah: 2*h**2/(3 - 2*h),
    ac: 2*(zc*h)**2/(3 - 2*zc*h),
}
q = 1 - y**2
report = {
    'scope': 'Exact source-level endpoint identities only; no Maple/Julia execution or spectral recomputation.',
    'assumption': 'zc>1, with physical horizon relations imposed exactly before limits.',
    'sympy_version': sp.__version__,
    'sectors': [],
}
print(report['scope'])
print(report['assumption'])
for sector in ('s0_s1', 's2odd', 's2even'):
    path = ROOT/'code_review/extracted'/f'Matrix_Assembler_{sector}-Final Version.decoded.txt'
    text = path.read_text()
    active = text.split('===== Group 0 Input Equation 0 =====', 1)[1].split('===== Group 0 Output', 1)[0]
    assert 'xi := evalf(cos((2*i - 1)*Pi/(2*n))):' in active
    assert 'Lobatto' not in active
    definitions = {}
    for line in active.splitlines():
        match = re.fullmatch(r'\s*(L\d\d)\s*:=\s*y\s*->\s*(.*):', line)
        if match:
            name, expr = match.groups()
            definitions[name] = sp.sympify(expr.replace('^','**'), locals=namespace)
    assert len(definitions)==9
    coeffs = {name:sp.cancel(expr.subs(substitutions, simultaneous=True)) for name,expr in definitions.items()}
    reduced = {name:sp.cancel(expr/q) for name,expr in coeffs.items()}
    endpoint_results = {}
    for endpoint in (-1,1):
        for name,expr in coeffs.items():
            assert sp.cancel(expr.subs(y,endpoint)) == 0, (sector,name,endpoint,'raw')
        limits = {name:sp.factor(expr.subs(y,endpoint)) for name,expr in reduced.items()}
        for name,value in limits.items():
            assert not value.has(sp.nan,sp.zoo,sp.oo,-sp.oo),(sector,name,endpoint,'finite')
            derivative_limit = sp.cancel(-sp.diff(coeffs[name],y).subs(y,endpoint)/(2*endpoint))
            assert sp.cancel(value-derivative_limit)==0,(sector,name,endpoint,'derivative identity')
        assert limits['L02']==0
        assert limits['L01']!=0
        endpoint_results[str(endpoint)] = {name:str(value) for name,value in limits.items()}
    report['sectors'].append({
        'sector':sector, 'source':str(path),
        'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
        'active_grid':'Chebyshev-Gauss roots only; no optional Lobatto branch',
        'raw_endpoint_coefficient_zeros':18,
        'finite_reduced_limits_and_derivative_identities':18,
        'reduced_endpoint_coefficients':endpoint_results,
    })
    print(f'{sector}: PASS roots grid, 18 raw endpoint zeros, 18 finite quotient limits, 18 derivative-limit identities.')
    print(f'{sector}: PASS reduced L02 vanishes, reduced L01 is nonzero at both endpoints (generic first-order row).')

report['raw_endpoint_coefficient_zeros']=54
report['finite_reduced_limits_and_derivative_identities']=54
report['all_passed']=True
(OUT/'sds_endpoint_rows_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print('PASS: all exact endpoint tests. No code or manuscript was modified.')
print('Limit-row rule: Lbar_ij(sigma) = -dL_ij/dy(sigma)/(2 sigma), sigma = +/-1, after exact horizon identities.')
print('QNM branch admissibility is a separate question and is not established by these endpoint equations.')
