"""Independent exact algebra checks; no numerical QNM computation."""
import sympy as s

r, x, M, c, L, a, b, q = s.symbols('r x M c L a b q', nonzero=True)
Hdim = 2*(9*M**3 + 3*c**2*M*r**2 + c**2*(1+c)*r**3
          + 3*M**2*(3*c*r-L*r**3))/(r**3*(3*M+c*r)**2)
U = (8*x**3*(c**3+c**2-3*L)+12*c**2*x**2+18*c*x+9)/(x**3*(2*c*x+3)**2)
conversion = s.factor(Hdim.subs({M:1, r:2*x}) - U/4)
assert conversion == 0

D = c**2*(1+c)-3*L
Nq = s.Poly(s.expand(18+18*c*(a+b*q)+6*c**2*(a+b*q)**2
                    +2*D*(a+b*q)**3),q)
Qq = s.Poly(s.expand(L*(a+b*q)**3*(3+c*(a+b*q))**2),q)
assert s.factor(Hdim.subs({M:1,r:a+b*q})/L-Nq.as_expr()/Qq.as_expr()) == 0

# Formal division in q, avoiding differentiation of the rational expression.
coeff=[]
for j in range(5):
    coeff.append(s.factor((Nq.nth(j)-sum(Qq.nth(k)*coeff[j-k]
                             for k in range(1,min(j,Qq.degree())+1)))/Qq.nth(0)))
residual = s.Poly(s.expand(Qq.as_expr()*sum(coeff[j]*q**j for j in range(5))
                            -Nq.as_expr()),q)
for j in range(5):
    assert s.factor(residual.nth(j)) == 0

print('Dimensional H(r), r=2x: H=U/4 verified exactly.')
print('Vbar=H/Lambda verified exactly.')
print('Rational q-series recurrence verified through q^4 exactly.')
print('q0 numerator:', s.factor(Nq.nth(0)))
print('q1 numerator:', s.factor(Nq.nth(1)))
print('q2 numerator:', s.factor(Nq.nth(2)))
print('q3 numerator:', s.factor(Nq.nth(3)))
print('Vbar coefficient0:', coeff[0])
print('Vbar coefficient1:', coeff[1])

# Both endpoints can use the same formulas.  At the right endpoint b must be -d.
for cv in (2,5):
    for lv in (.06,.08,.09):
        roots = sorted(float(s.re(z)) for z in s.nroots(r-2-lv*r**3/3)
                       if abs(float(s.im(z)))<1e-12 and s.re(z)>0)
        rh,rc=roots
        vals=[]
        for av,bv in ((rh,rc-rh),(rc,rh-rc)):
            vals.append([float(v.subs({a:av,b:bv,c:cv,L:lv})) for v in coeff[:3]])
        print({'c':cv,'Lambda':lv,'rh':rh,'rc':rc,
               'left_q_coefficients':vals[0],'right_q_coefficients':vals[1]})
