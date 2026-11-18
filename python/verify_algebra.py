import sympy as s

z, y, c, M, Om, eps, ell = s.symbols('z y c M Om eps ell', real=True)
H = (c*c+c+1)/(c*(c+1))
L = 3/(4*H**2*(c*c+c+1))
C = H*c
Z0 = -(1+c)
checks = []

def check(name, expr):
    residual = s.cancel(expr)
    assert residual == 0, (name, residual)
    checks.append(name)
    print('PASS:',name,flush=True)

f = 1-1/(H*z)-s.Rational(4,3)*L*H**2*z*z
fac = -4*L*H**2/(3*z)*(z-1)*(z-c)*(z-Z0)
check('Exact factorization in z', f-fac)
check('Event horizon root', f.subs(z,1))
check('Cosmological horizon root', f.subs(z,c))
check('Event horizon derivative', s.diff(f,z).subs(z,1)-(3-2*H)/H)
check('Cosmological horizon derivative', s.diff(f,z).subs(z,c)-H*(3-2*C)/C**2)

r = s.symbols('r', positive=True)
fr = 1-2*M/r-(L/M**2)*r*r/3
Ur = fr*(eps*s.diff(fr,r)/r+ell*(ell+1)/r**2)
Vz = f*(eps*s.diff(f,z)/z+ell*(ell+1)/z**2)
check('Potential scaling in Eq13', (2*M*H)**2*Ur.subs(r,2*M*H*z)-Vz)
check('Frequency scaling in Eq13', (2*M*H)**2*(Om/M)**2-4*H**2*Om**2)

ah = 2*H**2/(3-2*H)
ac = 2*C**2/(3-2*C)
check('Event indicial exponent', (-s.I*ah*Om)**2+(2*H*Om/s.diff(f,z).subs(z,1))**2)
check('Cosmological indicial exponent', (s.I*ac*Om)**2+(2*H*Om/s.diff(f,z).subs(z,c))**2)

zy = (c+1)/2+(c-1)*y/2
fy = s.factor(f.subs(z,zy))
ze = (c+1-2*Z0)/(c-1)
zt = (c+1)/(c-1)
fyp = L*H**2*(c-1)**2/3*(1-y*y)*(y+ze)/(y+zt)
check('Exact factorization in y', fy-fyp)
m2 = L**2*H**4/36*(c-1)**6*(1-y*y)**2*(y+ze)**2/(y+zt)**2
n0 = L**2*H**4/36*(c-1)**6*(1-y*y)*(y+ze)/(y+zt)**2*(2*Z0*(1-y*y)/((c-1)*(y+zt))-2*y*(y+ze))
check('Corrected M2', m2-(c-1)**2*fy**2/4)
check('Corrected N0', n0-(c-1)**2*fy*s.diff(fy,y)/4)
A2 = H**2*(c-1)**4/(36*c**2)*(9*c**2-4*L**2*H**2*ac**2*(c-1)**2*(c-Z0)**2)
B2 = H**2*(c-1)**4/36*(9-4*L**2*H**2*ah**2*(c-1)**2*(Z0-1)**2)
check('A2 vanishes identically', A2)
check('B2 vanishes identically', B2)

# Check the P_j transformation for an arbitrary lapse through the logarithmic
# derivative, without expanding the unrelated rational horizon relations.
a,b,F,Fp,V,Hsym=s.symbols('a b F Fp V Hsym')
D=(z-1)*(z-c)
B=(b-a)*z+c*a-b
E=(b-a)*z*z+2*(c*a-b)*z-c*c*a+b
logp=-s.I*a*Om/(z-1)+s.I*b*Om/(z-c)
p1=D*F*(D*Fp+2*s.I*Om*F*B)
p0=-Om**2*(F*F*B*B-4*Hsym**2*D*D)+s.I*Om*(D*F*Fp*B-F*F*E)-D*D*V
check('P1 from logarithmic derivative',p1-D*D*(2*F*F*logp+F*Fp))
check('P0 from logarithmic derivative',p0-D*D*(F*F*(s.diff(logp,z)+logp**2)+F*Fp*logp+4*Hsym**2*Om**2-V))

# The quadratic endpoint coefficient is the indicial polynomial times the
# nonzero squared separation of the horizons. These identities verify its
# cancellation exactly and leave a cubic or higher zero.
kh=(3-2*H)/H
kc=H*(3-2*C)/C**2
check('S0 quadratic cancellation at left endpoint',(c-1)**2*(4*H**2*Om**2-kh**2*ah**2*Om**2))
check('S0 quadratic cancellation at right endpoint',(c-1)**2*(4*H**2*Om**2-kc**2*ac**2*Om**2))

print('All',len(checks),'symbolic checks passed for arbitrary c=z_c>1.')
