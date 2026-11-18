#!/usr/bin/env python3
"""Exact checks of the independently derived rescaled SdS collocation equation.

Uses q=xc/xh>1 to impose all physical horizon identities exactly.  No numerical
spectral data or manuscript coefficient formulas are used.
"""
from pathlib import Path
import json
import sympy as s

q,y,w,L=s.symbols('q y w L', positive=True)
D=q*q+q+1
h=D/(q*(q+1)); c=q*h; x0=-h-c; d=c-h; a=(h+c)/2
lam=3/(4*h*h*D); sigma=lam*d; x=a+d*y/2
F=(x-x0)/(3*x); g=(1-y*y)*F; f=lam*d*d*g
ph=3*h/(2*(h-x0)); pc=3*c/(2*(c-x0))
H=-ph/(1+y)+pc/(1-y)
C=3*a*a*d/((2*h+c)*(h+2*c))
G=y+C*(1-y*y)/x
R=(1-2*y*C/x-(C/x)**2*(1-y*y))/F
checks=[]

def check(name, actual, expected):
    residual=s.factor(actual-expected)
    checks.append({'name':name,'residual':str(residual),'passed':residual==0})
    print(name, 'PASS' if residual==0 else f'FAIL: {residual}', flush=True)

check('scaled_metric',f,1-1/x-s.Rational(4,3)*lam*x*x)
check('black_hole_phase',ph,sigma*2*h*h/(3-2*h))
check('cosmological_phase',pc,-sigma*2*c*c/(3-2*c))
check('G_is_gH',G,g*H)
check('cancellation_free_R',R,(1-G*G)/g)
check('rescaled_second_derivative',4*f*f/d**2/(4*sigma*sigma),g*g)
check('rescaled_first_derivative',4*f*s.diff(f,y)/d**2/(4*sigma*sigma),g*s.diff(g,y))
check('phase_linear_term',g*s.diff(H,y)+s.diff(g,y)*H,s.diff(G,y))
check('phase_quadratic_term',(1-g*g*H*H)/g,R)
check('event_horizon_G',G.subs(y,-1),-1)
check('cosmological_horizon_G',G.subs(y,1),1)

for sector in ['scalar','em','odd']:
    Vraw=f*(L/x**2+(2*s.diff(f,y)/(d*x) if sector=='scalar' else 0)
            -(3/x**3 if sector=='odd' else 0))
    Vbar=(L/x**2+(2*lam*d*s.diff(g,y)/x if sector=='scalar' else 0)
          -(3/x**3 if sector=='odd' else 0))/(4*lam)
    check(f'{sector}_rescaled_potential',Vraw/(4*sigma*sigma),g*Vbar)
    check(f'{sector}_Nariai_potential',s.limit(s.cancel(Vbar),q,1),
          L-(2 if sector=='odd' else 0))

for name,expr,limit in [('g',g,1-y*y),('gy',s.diff(g,y),-2*y),
                        ('G',G,y),('Gy',s.diff(G,y),1),('R',R,1),
                        ('ph',ph,s.Rational(1,2)),('pc',pc,s.Rational(1,2))]:
    check(f'Nariai_limit_{name}',s.limit(s.cancel(expr),q,1),limit)

report={'scope':__doc__,'sympy':s.__version__,'passed':sum(c['passed'] for c in checks),
        'total':len(checks),'checks':checks}
Path(__file__).with_name('rescaled_equation_checks.json').write_text(json.dumps(report,indent=2)+'\n')
print(f"{report['passed']}/{report['total']} exact checks passed.")
raise SystemExit(0 if report['passed']==report['total'] else 1)
