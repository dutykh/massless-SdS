"""Independent symbolic check of the SdS endpoint cancellations and Lobatto rows.

Temporary verification only; does not edit the manuscript or reproduce numerics.
"""
import sympy as s

y, c, Om, eps, L = s.symbols('y c Omega epsilon L', real=True)
d = c - 1
D = c*c+c+1
xh = D/(c*(c+1))
La = 3*c*c*(c+1)**2/(4*D**3)
xc, z0 = xh*c, -(1+c)
ah, ac = 2*xh*xh/(3-2*xh), 2*xc*xc/(3-2*xc)
zt, ze = (c+1)/d, (c+1-2*z0)/d
p = 1-y*y
g = s.factor(La*xh*xh*d*d/3*(y+ze)/(y+zt))
f = p*g
B = (ac-ah)*y+ac+ah
a = ac-ah
U = 4/(d*d*(y+zt))*(L/(y+zt)+eps*s.diff(f,y))
V = f*U

# Compact coefficient expressions, independently obtained from the factored ODE.
M2 = d*d*f*f/4
N0 = d*d*f*s.diff(f,y)/4
N1 = -d*d*f*f*B/(2*p)
C2 = d*d/4*(xh*xh*d*d-g*g*B*B)
C1 = -d*d*p/4*(a*g*g+g*s.diff(g,y)*B)
C0 = -d**4*V/16

# Verify the corrected expanded expressions occurring in the manuscript.
M2expanded = La**2*xh**4*d**6/36*p*p*(y+ze)**2/(y+zt)**2
N0expanded = La**2*xh**4*d**6/36*p*(y+ze)/(y+zt)**2*(2*z0*p/(d*(y+zt))-2*y*(y+ze))
C1expanded = -La**2*xh**4*d**6/36*(y+ze)/(y+zt)**2*(B*(2*z0*p/(d*(y+zt))-2*y*(y+ze))+(y+ze)*((ac-ah)*y*y+2*(ac+ah)*y+ac-ah))

def check(name, expr):
    result = s.factor(expr)
    assert result == 0, (name, result)
    print(name + ': 0')

check('M2 corrected expanded identity', M2-M2expanded)
check('N0 corrected expanded identity', N0-N0expanded)
check('C1 compact identity', C1-C1expanded)
A2 = xh*xh*d**4/(36*c*c)*(9*c*c-4*La*La*xh*xh*ac*ac*d*d*(c-z0)**2)
B2 = xh*xh*d**4/36*(9-4*La*La*xh*xh*ah*ah*d*d*(z0-1)**2)
check('A2 physical identity', A2)
check('B2 physical identity', B2)

for sig in (-1,1):
    gv = s.factor(g.subs(y,sig))
    Bv = s.factor(B.subs(y,sig))
    Kv = s.factor((a+B*s.diff(g,y)/g).subs(y,sig))
    Uv = s.factor(U.subs(y,sig))
    check(f'gB relation {sig}', gv*Bv+sig*xh*d)
    for name, expr in [('M2',M2),('N0',N0),('N1',N1),('C2',C2),('C1',C1),('C0',C0)]:
        check(f'endpoint zero {name} {sig}', s.cancel(expr).subs(y,sig))
    reduced = {name:s.factor(s.cancel(expr/p).subs(y,sig)) for name,expr in [('M2',M2),('N0',N0),('N1',N1),('C2',C2),('C1',C1),('C0',C0)]}
    factor = d*d*gv*gv/4
    expected = dict(M2=0,N0=-2*sig*factor,N1=-2*Bv*factor,C2=sig*Bv*Kv*factor,C1=-Kv*factor,C0=-d*d*Uv/(4*gv)*factor)
    for name in reduced:
        check(f'reduced coefficient {name} {sig}', reduced[name]-expected[name])
    # Full limiting equation divided by nonzero factor d^2*g_sig^2/4:
    # -2(sig+i Om Bsig) Phi' + [(sig Om^2 Bsig-i Om)Ksig
    #                                  -d^2 U_sig/(4g_sig)] Phi = 0.
    print(f'g({sig}) = {gv}')
    print(f'B({sig}) = {Bv}')
    print(f'K({sig}) = {Kv}')

print('All symbolic identities verified for c=zc>1 (nonextremal SdS).')
