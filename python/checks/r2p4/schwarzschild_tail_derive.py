"""Derive and verify the minimal asymptotic tail used by schwarzschild_cf.py."""
from pathlib import Path
import json
import importlib.util
import sympy as s
u,rho,L,eps=s.symbols('u rho L eps')
a=s.symbols('c1:6')
r=1+sum(a[k-1]*u**k for k in range(1,6))
rprev=1+sum(a[k-1]*u**k*(1-u*u)**(-s.Rational(k,2)) for k in range(1,6))
alpha=1+(2*rho+2)*u*u+(2*rho+1)*u**4
beta=-2-(8*rho+2)*u*u-(8*rho*rho+4*rho+L-eps)*u**4
gamma=1+4*rho*u*u+(4*rho*rho-eps-1)*u**4
expr=s.series(alpha*r+beta+gamma/rprev,u,0,8).removeO().expand()
sol={a[0]:-s.sqrt(2*rho)}
for power,var in zip(range(3,7),a[1:]):
 c=s.simplify(expr.coeff(u,power).subs(sol))
 v=s.simplify(s.solve(c,var)[0]); sol[var]=v
 print(var,'=',v)
residuals=[s.simplify(expr.coeff(u,k).subs(sol)) for k in range(2,7)]
assert residuals == [0]*5
# Check the actual numerical implementation's tail formula symbolically.
path=Path(__file__).with_name('schwarzschild_cf.py')
spec=importlib.util.spec_from_file_location('schwarzschild_cf',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
class SymbolicMP:
    sqrt=staticmethod(s.sqrt)
    mpf=staticmethod(s.Rational)
module.mp=SymbolicMP
N=s.symbols('N',positive=True)
coded=module.tail_ratio(rho,L,eps,N,5)
expected=1+sum(sol[a[k-1]]*N**(-s.Rational(k,2)) for k in range(1,6))
assert s.simplify(coded-expected)==0
print('PASS: recurrence coefficients of u^2 through u^6 vanish.')
print('PASS: actual solver tail matches all five derived coefficients.')
Path(__file__).with_suffix('.json').write_text(json.dumps({
    'scope':'Exact asymptotic-tail derivation and comparison with actual numerical source; no new spectral solve.',
    'coefficients':{str(k):str(v) for k,v in sol.items()},
    'recurrence_residuals':[str(v) for v in residuals],
    'coded_tail_residual':'0',
},indent=2)+'\n')
