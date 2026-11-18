#!/usr/bin/env python3
"""New, independent low-overtone Schwarzschild continued-fraction calculation.

Units: x=r/(2M), w=2M*omega=2*Omega, rho=-i*w. The Regge-Wheeler
potential is (1-1/x)[ell(ell+1)/x^2+(1-s^2)/x^3]. Leaver's series is
psi=(x-1)^rho*x^(-2rho)*exp[-rho*(x-1)]*sum a_k[(x-1)/x]^k.
Its prefactor is event-ingoing and infinity-outgoing; the minimal solution
of the recurrence selects the outgoing branch by analytic continuation.
Only positive-real, low-overtone roots are computed here, NOT the NIA.

Primary sources:
E.W. Leaver, Proc. R. Soc. Lond. A402,285-298(1985), Eqs.(3)-(14),
https://blackholes.tecnico.ulisboa.pt/gritting/pdf/black_holes/Leaver_An-Analytic-representation-for-the-quasi-normal-modes-of-Kerr-black-holes.pdf
E.W. Leaver, 'Remarks on the continued-fraction method ...' (1991),
Eqs.(1)-(9), https://www.edleaver.com/Misc/EdLeaver/Publications/RemarksOnContinuedFractionMethodForComputingBlackHoleQuasinormalFrequenciesAndModes.pdf
The asymptotic tail coefficients are derived from that recurrence by
schwarzschild_tail_derive.py. This is not the uploaded Maple/Julia pipeline
and does not establish historical provenance for the manuscript tables.
"""
from pathlib import Path
import argparse, json, time
import mpmath as mp

OUT = Path(__file__).resolve().parent
# Seeds use only rough locations of the existing oscillatory table entries.
SEEDS = [
(0,0,0,'0.1105','-0.1049'),(0,0,1,'0.0861','-0.3481'),(0,0,2,'0.0757','-0.6011'),
(0,1,0,'0.2929','-0.0977'),(0,1,1,'0.2644','-0.3063'),(0,1,2,'0.2295','-0.5401'),(0,1,3,'0.2033','-0.7883'),
(0,2,0,'0.4836','-0.0968'),(0,2,1,'0.4639','-0.2956'),(0,2,2,'0.4305','-0.5086'),(0,2,3,'0.3939','-0.7381'),
(1,1,0,'0.2483','-0.0925'),(1,1,1,'0.2145','-0.2937'),
(1,2,0,'0.4576','-0.0950'),(1,2,1,'0.4365','-0.2907'),
(2,2,0,'0.3736','-0.0889'),(2,2,1,'0.3463','-0.2735'),
(2,3,0,'0.5994','-0.0927'),(2,3,1,'0.5826','-0.2813'),
]

def tail_ratio(rho, angular, eps, depth, order=5):
    """Minimal a_(N+1)/a_N asymptotic tail, Re(sqrt(2rho))>0."""
    if order == 0:
        return mp.mpc(0)
    if order not in (2,3,4,5):
        raise ValueError('tail order must be 0,2,3,4,5')
    cs = [
        -mp.sqrt(2*rho),
        2*rho-mp.mpf(3)/4,
        mp.sqrt(2)*(-16*angular-64*rho**2+80*rho-3)/(64*mp.sqrt(rho)),
        (16*angular+16*rho*(4*angular-4*eps-16*rho+5)+3)/(128*rho),
        mp.sqrt(2)*(256*angular**2-288*angular+128*rho**2*(16*angular+32*eps+32*rho**2+48*rho-79)+224*rho*(16*angular+3)-63)/(8192*rho**mp.mpf('1.5')),
    ]
    u=1/mp.sqrt(depth)
    return 1+sum(cs[k-1]*u**k for k in range(1,order+1))

def cf_residual(omega, spin, ell, inversion, depth, tail_order=5):
    """Inverted Leaver fraction. omega is manuscript Omega=M*omega_phys."""
    rho=-2j*omega
    angular=ell*(ell+1)
    eps=spin*spin-1
    aa=2*rho+2; ab=2*rho+1
    ba=8*rho+2; bb=8*rho*rho+4*rho+angular-eps
    ga=4*rho; gb=4*rho*rho-eps-1
    def alpha(k):return k*k+aa*k+ab
    def beta(k):return -(2*k*k+ba*k+bb)
    def gamma(k):return k*k+ga*k+gb
    ratio=tail_ratio(rho,angular,eps,depth,tail_order)
    for k in range(depth,inversion,-1):
        ratio=-gamma(k)/(beta(k)+alpha(k)*ratio)
    left=beta(0)
    for k in range(1,inversion+1):
        left=beta(k)-alpha(k-1)*gamma(k)/left
    return left+alpha(inversion)*ratio

def solve(spin,ell,n,seed,depth,dps,tail_order=5,inversion=None):
    with mp.workdps(dps):
        inv=n if inversion is None else inversion
        z=mp.mpc(seed)
        dz=mp.mpc('0.000001','0.0000007')
        f=lambda w:cf_residual(w,spin,ell,inv,depth,tail_order)
        root=mp.findroot(f,(z,z+dz),tol=mp.mpf(10)**(-min(dps-10,45)),maxsteps=30)
        residual=abs(f(root))
        if not (mp.re(root)>0 and mp.im(root)<0):raise AssertionError('Wrong quadrant')
        return root,residual

def string(z,digits=45):return mp.nstr(z,digits)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode',type=int,help='Only the zero-based index in SEEDS, for bounded diagnostics')
    args=parser.parse_args()
    started=time.time()
    rows=[]
    for idx,(spin,ell,n,re,im) in enumerate(SEEDS):
        if args.mode is not None and args.mode!=idx:continue
        seed=mp.mpc(re,im)
        steps=[];previous=None
        for depth,dps in [(256,40),(512,40),(1024,40),(1024,70)]:
            root,residual=solve(spin,ell,n,seed,depth,dps)
            with mp.workdps(75):
                difference=None if previous is None else abs(root-previous)
                steps.append({'depth':depth,'decimal_digits':dps,'tail_order':5,
                    'real':string(mp.re(root)),'imag':string(mp.im(root)),
                    'absolute_change':None if difference is None else string(difference),
                    'fraction_residual':string(residual)})
            seed=root;previous=root
        # Different tail order at the final depth independently probes truncation.
        root_tail3,res_tail3=solve(spin,ell,n,root,1024,70,tail_order=3)
        # Different inversion checks equivalent defining fractions at the same depth.
        root_inv,res_inv=solve(spin,ell,n,root,1024,70,inversion=n+1)
        with mp.workdps(75):
            row={'s':spin,'ell':ell,'n':n,'steps':steps,
                 'tail_order_3_change':string(abs(root_tail3-root)),
                 'inversion_plus_1_change':string(abs(root_inv-root)),
                 'final_real':string(mp.re(root)), 'final_imag':string(mp.im(root)),
                 'rounded_4dp':f'{float(mp.re(root)):.4f}{float(mp.im(root)):+.4f}i'}
            rows.append(row)
            print(f's={spin} ell={ell} n={n}: {string(root,20)}; depth change={steps[-2]["absolute_change"]}; precision change={steps[-1]["absolute_change"]}; tail3 change={row["tail_order_3_change"]}',flush=True)
        report={'scope':'New independent low-overtone Schwarzschild CF calculation, not historical Maple/Julia reproduction; no NIA classification.',
                'mpmath_version':mp.__version__,'normalization':'Omega=M*omega_phys; rho=-2i*Omega',
                'rows':rows,'elapsed_seconds':time.time()-started}
        suffix='' if args.mode is None else f'_mode{args.mode}'
        (OUT/f'schwarzschild_cf_results{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'Completed {len(rows)} modes in {time.time()-started:.1f}s.',flush=True)

if __name__=='__main__':main()
