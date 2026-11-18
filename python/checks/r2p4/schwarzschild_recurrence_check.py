#!/usr/bin/env python3
"""Exact independent derivation of the CF recurrence from the RW master ODE.
No spectral solve and no use of the uploaded Maple/Julia coefficients.
"""
from pathlib import Path
import json
import sympy as sp
z,rho,ell,spin,n=sp.symbols('z rho ell spin n')
x=1/(1-z)
angular=ell*(ell+1)
# z=1-1/x; logarithmic derivative of Leaver's prefactor with respect to z.
p = rho/z-rho/(1-z)-rho/(1-z)**2
# f d_x(f d_x psi) - [rho^2 + V]psi = 0, because w=i*rho.
a2=z**2*(1-z)**4
a1=z*(1-z)**3*(1-3*z)
v=z*(angular*(1-z)**2+(1-spin**2)*(1-z)**3)
scale=z*(1-z)**2
q2=sp.cancel(a2/scale)
q1=sp.cancel((a1+2*a2*p)/scale)
q0=sp.cancel((a2*(sp.diff(p,z)+p*p)+a1*p-rho**2-v)/scale)
q2,q1,q0=(sp.Poly(q,z) for q in (q2,q1,q0))
alpha=q2.coeff_monomial(z)*(n+1)*n+q1.coeff_monomial(1)*(n+1)
beta=q2.coeff_monomial(z**2)*n*(n-1)+q1.coeff_monomial(z)*n+q0.coeff_monomial(1)
gamma=q2.coeff_monomial(z**3)*(n-1)*(n-2)+q1.coeff_monomial(z**2)*(n-1)+q0.coeff_monomial(z)
expected=[n*n+(2*rho+2)*n+2*rho+1,-(2*n*n+(8*rho+2)*n+8*rho*rho+4*rho+angular-(spin*spin-1)),n*n+4*rho*n+4*rho*rho-spin*spin]
residuals=[sp.factor(a-b) for a,b in zip((alpha,beta,gamma),expected)]
assert residuals==[0,0,0]
print('PASS: exact alpha_n, beta_n, gamma_n recurrence identities from the original Schwarzschild master equation.')
print('Leaver prefactor: (x-1)^rho*x^(-2rho)*exp[-rho*(x-1)], rho=-2i*Omega.')
print('Horizon: prefactor/(x-1)^rho tends to1. Infinity: prefactor/[exp(-rho*x)*x^(-rho)] tends toexp(rho).')
print('The infinity factor is outgoing for exp(-i omega t); the CF minimal-sequence condition supplies the global selection.')
print('Scope: low oscillatory frequencies only; no statement about the negative imaginary axis.')
report={'scope':__doc__,'coefficient_residuals':[str(a) for a in residuals],'residual_ode_polynomials':[str(q.as_expr()) for q in(q2,q1,q0)],'sympy_version':sp.__version__}
Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
