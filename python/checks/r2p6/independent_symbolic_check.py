"""Independent algebra checks for the scalar SdS stability addition.

This checks exact identities, not a numerical QNM computation or a PDE theorem.
The domain and boundary conclusions are assessed separately in the report.
"""
from pathlib import Path
import hashlib
import json
import sympy as s

out = Path(__file__).resolve().parent
root = out.parents[1]
r, M, Lam, ell = s.symbols('r M Lambda ell', positive=True)
f = 1 - 2*M/r - Lam*r**2/3
psi = s.Function('psi')(r)
dx = lambda u: f*s.diff(u, r)
D = lambda u: dx(u) - f*u/r
Dadj = lambda u: -dx(u) - f*u/r
V = f*(s.diff(f, r)/r + ell*(ell+1)/r**2)
H = lambda u: -dx(dx(u)) + V*u
H0 = lambda u: -dx(dx(u)) + f*s.diff(f, r)*u/r
zero = {
    'factorization_coefficient': s.simplify(dx(f/r)+(f/r)**2-f*s.diff(f,r)/r),
    'full_operator_identity': s.simplify(H(psi)-Dadj(D(psi))-f*ell*(ell+1)*psi/r**2),
    'static_first_order_residual': s.simplify(D(r)),
    'static_master_residual': s.simplify(H0(r)),
    'second_static_solution_residual': s.simplify(H0(r*s.Integral(1/(r**2*f),r))),
}
x, w, a, b, c, d = s.symbols('x omega a b c d', real=True)
R = a+s.I*b
T = c+s.I*d
inc = s.exp(-s.I*w*x) + R*s.exp(s.I*w*x)
trans = T*s.exp(-s.I*w*x)
current = lambda u: s.simplify((s.conjugate(u)*s.diff(u,x)-u*s.conjugate(s.diff(u,x)))/(2*s.I))
zero['incident_current_identity'] = s.simplify(current(inc)-w*(a*a+b*b-1))
zero['transmitted_current_identity'] = s.simplify(current(trans)+w*(c*c+d*d))

assert all(v == 0 for v in zero.values()), zero
manuscript = root/'revision/main_referee2_point6.tex'
response = root/'revision/response_referee2_point6.tex'
report = {
    'scope': __doc__,
    'sympy_version': s.__version__,
    'identities_passed': len(zero),
    'exact_residuals': {k: str(v) for k,v in zero.items()},
    'domain_assessment': {
        'conditions': 'M>0 and 0<M^2 Lambda<1/9, minimally coupled massless scalar',
        'form_domain': 'H^1(R), with closure taken in the form norm (q+L2 norm); bounded f/r makes this norm equivalent to H1',
        'operator_domain': 'H^2(R), since the real potential is smooth and bounded',
        'upper_half_plane': 'The prescribed branches and their derivatives decay exponentially; integration by parts has zero endpoint terms',
        'bounded_static_monopole': 'Unique up to scale: psi=C r; the other independent solution has logarithmic singularities at both simple horizons',
        'zero_not_L2': 'For C!=0, psi tends to the nonzero constants C r_h and C r_c at the two infinite tortoise-coordinate ends',
        'limitation': 'Mode stability only; no claim of spectral completeness, nonlinear stability, quantitative decay, or a damped-QNM L2 identity',
    },
    'source_hashes': {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in [manuscript,response] if p.exists()},
    'manuscript_modified': False,
    'numerical_QNM_run': False,
}
(out/'independent_symbolic_check.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
