#!/usr/bin/env python3
"""Symbolic checks for physical scalar/Maxwell master definitions.
No manuscript or numerical spectrum is changed.
"""
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent
r,w,L,Lam,delta,rm,chi=s.symbols('r omega L Lambda delta rm chi', nonzero=True)
f=s.Function('f')(r); p=s.Function('psi')(r)
checks=[]
def check(name,expr):
    residual=s.simplify(s.expand(expr))
    checks.append({'name':name,'passed':residual==0,'residual':str(residual)})
# Scalar field phi=e^{-iwt}psi(r)Y/r, Delta_S2 Y=-L Y.
scalar=(s.diff(r*r*f*s.diff(p/r,r),r)/(r*r)+(w*w/f-L/(r*r))*p/r)*r*f
scalar_master=f*s.diff(f*s.diff(p,r),r)+(w*w-f*(s.diff(f,r)/r+L/r**2))*p
check('minimally_coupled_scalar_reduction',scalar-scalar_master)
# Maxwell even gauge-invariant p_a=a_a-d_a a_e.
# Psi_e=r^2 F_tr/L; the two Maxwell constraints give these components.
pt=-f*s.diff(p,r);pr=s.I*w*p/f
Ftr=-s.I*w*pr-s.diff(pt,r)
em_master=f*s.diff(f*s.diff(p,r),r)+(w*w-f*L/r**2)*p
check('even_Maxwell_Ftr_identity_gives_master',f*(Ftr-L*p/r**2)-em_master)
check('even_Maxwell_t_constraint',s.diff(L*p,r)+L*pt/f)
check('even_Maxwell_r_constraint',s.I*w*L*p-f*L*pr)
check('odd_Maxwell_orbit_wave_operator',f*(s.diff(f*s.diff(p,r),r)+w*w*p/f-L*p/r**2)-em_master)
# X_A=epsilon_A{}^B D_B Y, epsilon_theta_phi=sin(theta).
th,ph=s.symbols('theta phi',real=True);Y=s.Function('Y')(th,ph)
Xt=s.diff(Y,ph)/s.sin(th);Xp=-s.sin(th)*s.diff(Y,th)
DeltaY=s.diff(s.sin(th)*s.diff(Y,th),th)/s.sin(th)+s.diff(Y,ph,2)/s.sin(th)**2
check('odd_vector_curl_sign',s.diff(Xp,th)-s.diff(Xt,ph)+s.sin(th)*DeltaY)
# On Nariai, r0^2=1/Lambda and epsilon_tau_chi=1/Lambda.
# General invariant master -r0^2 epsilon^{ab} F_ab/(2L)=F_tau_chi/L.
Ftc=s.symbols('F_tau_chi')
check('Nariai_even_field_strength_normalization',-(1/Lam)/(2*L)*(-2*Lam*Ftc)-Ftc/L)
# t=tau/(Lambda delta), r=rm+delta chi, so pulled-back F_tau_chi=F_tr/Lambda.
amp=s.symbols('amplitude')
pullback=(1/(Lam*delta))*delta*(L*amp/r**2)
check('SdS_to_Nariai_Maxwell_master_pullback',pullback.subs(r**2,1/Lam)-L*amp)
# Intrinsic product-space scalar and Maxwell orbit equation.
F=1-chi**2;P=s.Function('Psi')(chi);varpi=s.symbols('varpi')
box2=Lam*(varpi**2*P/F+s.diff(F*s.diff(P,chi),chi))
product=F/Lam*(box2-L*Lam*P)
nariai=F*s.diff(F*s.diff(P,chi),chi)+(varpi**2-L*F)*P
check('product_space_scalar_Maxwell_radial_equation',product-nariai)
# RW potential check of the distinct gravitational -2 dynamical term.
r0,M0=s.symbols('r0 M0',positive=True)
check('Nariai_gravitational_curvature_term',(L/r0**2-6*M0/r0**3).subs(M0,r0/3)*r0**2-(L-2))
report={'scope':'Physical variable definitions and exact reduction identities only; no new mode classification or spectral calculation.',
        'sympy_version':s.__version__,'passed':sum(c['passed'] for c in checks),'total':len(checks),'checks':checks,
        'conventions':{'time':'exp(-i omega t)','sphere':'unit sphere q_AB, epsilon_theta_phi=sin(theta), Delta Y=-L Y; L=ell(ell+1)',
        'odd_vector':'X_A=epsilon_A{}^B D_B Y, so D_A X_B-D_B X_A=L epsilon_AB Y',
        'even_master':'psi_e=r^2 [F_tr]_lm/L in SdS; Psi_e=[F_tau_chi]_lm/L in dimensionless Nariai coordinates',
        'odd_master':'psi_o=[epsilon^{AB}F_AB]_lm/(2L), with unit-sphere epsilon',
        'scalar_master':'phi=sum e^{-i omega t}psi_0(r)Y/r; Nariai phi=sum e^{-i varpi tau}Psi_0(chi)Y/r0'}}
(HERE/'scalar_maxwell_definitions_check.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
assert report['passed']==report['total']
