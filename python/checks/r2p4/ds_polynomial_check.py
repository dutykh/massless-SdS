"""Exact check of original scalar ODE used in ds_branch_check.py."""
import sympy as s
r,Lam,w,ell=s.symbols('r Lambda omega ell')
P=r-2-Lam*r**3/3
f=P/r
V=f*(ell*(ell+1)/r**2+s.diff(f,r)/r)
A=r**2*P**2
B=r*P*(r*s.diff(P,r)-P)
C=w**2*r**4-P*(ell*(ell+1)*r+2-2*Lam*r**3/3)
assert s.factor(A-r**4*f**2)==0
assert s.factor(B-r**4*f*s.diff(f,r))==0
assert s.factor(C-r**4*(w**2-V))==0
rho,n,a2=s.symbols('rho n a2')
assert s.expand(a2*(rho+n)*(rho+n-1)+a2*(rho+n)-a2*rho**2-a2*n*(n+2*rho))==0
print('PASS: polynomial coefficients are exactly r^4 times original scalar master equation.')
print('PASS: Frobenius recurrence denominator a2*n*(n+2*rho).')
print('r^4 and the coordinate-derivative factors are nonzero in the positive interhorizon interval; they introduce no spectral zeros.')
