# Scalar ell=2 imaginary-branch check

These are independent numerical radiation-boundary calculations for one
representative massless scalar branch. They are not interval-certified results,
an exhaustive search, or a theorem about every SdS branch. They use no
collocation eigenvalues as seeds or inputs.

## Equation and conventions

Set M=1, lambda=Lambda_tilde, Omega=omega=-i alpha, alpha>0. The physical time
dependence is exp(-i omega t). Put P(r)=r-2-lambda*r^3/3 and f=P/r. For ell=2,

    f^2 psi'' + f f' psi' + [omega^2-f(ell(ell+1)/r^2+f'/r)] psi = 0.

Multiplying by r^4 gives the polynomial equation A psi''+B psi'+C psi=0 with

    A=r^2 P^2,
    B=r P(r P'-P),
    C=-alpha^2 r^4-P[ell(ell+1)r+2-2lambda*r^3/3].

This identity and the indicial denominator are checked exactly with SymPy by
`ds_polynomial_check.py`. Multiplication by r^4 introduces no zeros in the
physical interval. No frequency-dependent denominator clearing is used.

At R=r_h choose sigma=+1 and q=(r-r_h)/D; at R=r_c choose sigma=-1 and
q=(r_c-r)/D, where D=r_c-r_h. Both radiation solutions have

    psi=q^rho F(q), rho=-alpha/|f'(R)|=-alpha/(2 kappa_R), F(0)=1.

The left solution is ingoing at the event horizon and the right solution is
outgoing at the cosmological horizon. These are the Frobenius sectors selected
by analytic continuation of the physical scattering solutions. Regularity of
an arbitrary factored solution is not the admissibility test.

Expanding the polynomial equation in q gives a finite-band coefficient
recurrence whose n-th denominator is A_2*n*(n+2rho). At all retained points,
nu=-2rho=alpha/kappa_R stays away from integers. No result is claimed at an
exceptional integer-index point.

## Series domains and determinant

The finite singularities of the equation are 0, r_h, r_c, and r_0=-r_h-r_c.
The event-horizon Taylor radius in q is min(1,r_h/D); the cosmological radius
is 1. Thus the default matching location

    t=(r-r_h)/D=min(1/2,0.45*r_h/D)

lies strictly inside both convergence disks. For the representative .06 case,
the event radius is 0.637889718523 and default t=0.287050373335; other selected
cases have event radius 1. Independent matching points are t=.35 at .06,
and t=.45 at all other selected parameters.

At the match point the nonzero q^rho factors are removed, but their logarithmic
derivatives are retained: the vector is (F, sigma*(F_q+rho*F/q)/D). Let v_h and
v_c be the two vectors. The real-axis matching function is

    Dcal=det(v_h,v_c)/(||v_h||_2 ||v_c||_2).

Its zeros coincide with the physical Jost Wronskian zeros; all removed factors
and norms are nonzero. A nonzero derivative with respect to alpha at a zero
implies a simple Jost zero. Each horizon solution has unit leading Frobenius
amplitude and is nonzero, so their Green-function numerator is not identically
zero. At the displayed representatives their matching-point values are also
nonzero. This excludes the zero-potential cancellation relevant to the exact
Nariai monopole and the endpoint-regularity ambiguity at exceptional indices.

## Continuation and numerical acceptance

The near-extremal campaign comprises 141 points:

| Interval | Lambda step | Decimal precision | Series order | Check order |
|---|---:|---:|---:|---:|
| .06 to .11 | .0005 | 80 | 320 | 360 |
| .11 to .111 | .00005 | 120 | 600 | 720 |
| .111 to .1111 | .000005 | 180 | 900 | 1050 |

The first guess is alpha=2 sqrt(lambda/3), the leading scalar ell=2 de Sitter
frequency. The second uses its local square-root slope. Subsequent predictors
use the previous two accepted roots and parameter values. A local secant search
uses predictor +/-1e-6 on the first interval and +/-1e-8 on the extensions,
with root tolerance 1e-35 and at most 30 iterations. A root is rejected if it
moves more than .008 from the predictor; the observed maximum predictor
correction after the initial seed is 4.86e-6. Extension frequency changes must
also be below .0002. This is an explicitly local search, not a global
frequency-window census.

At every accepted point the normalized determinant at the higher check order
must be below 1e-30, and both Frobenius-index distances from integers must be
greater than 1e-7. Observed maxima/minima were respectively 3.66e-49 and
5.81e-4. No exceptional sampled point required displacement. The index
screen applies to sampled points, not all unsampled parameters between them.

Each of the seven displayed representatives is re-solved at higher order and
precision, and at a different matching point. At .06 through .11 the checks
use (90 digits,N=400), (100 digits,N=480), and (100 digits,N=480,new point).
At .111 they use (120,600), (140,720), and (140,720,new point); at .1111 they
use (180,900), (200,1050), and (200,1050,new point). The largest root spread
is 2.47e-64. This many-digit agreement is numerical evidence, not a rigorous
error enclosure; report only conservative displayed digits in the paper.

The derivative of Dcal at the default match point is nonzero at every selected
parameter. Its values, in increasing lambda order, are approximately
3161.5968, 107.29768, 50.326494, -50.312328, -0.16990104, -0.010439966,
and 0.000039957744. Moving the match point changes this normalization and
its derivative, but not the computed root.

An independent evaluation in `ds_independent_audit.py` uses 220 decimal
digits, fixed refined roots at .06, .111, and .1111, and two matching points
per case. It sums the Frobenius terms directly, estimates cancellation, and
evaluates normalized W at alpha +/-1e-9. Both matching points exhibit a sign
change in each case. The worst lost-digit estimates are approximately 1.56,
16.71, and 47.49 respectively, leaving ample precision margin relative to the
production precisions. This audit does not perform a root search.

## Results and scope

| Lambda_tilde | -Im Omega | -Im varpi | -Im(r_h omega) |
|---:|---:|---:|---:|
| .06 | .2826623545036564 | 2.709360907059439 | .6270373804387558 |
| .08 | .3263244233032902 | 3.604514610155184 | .7645742076062627 |
| .09 | .3460849763778670 | 4.446896681832636 | .8413317419603011 |
| .10 | .3647700102588595 | 6.221474915619097 | .9330087111602768 |
| .11 | .3825369880767245 | 19.94765844737825 | 1.086059318221649 |
| .111 | .3842680392346436 | 63.16391472558721 | 1.132255975038725 |
| .1111 | .3844405936490753 | 199.7682213556714 | 1.146713905466271 |

Here varpi=Omega/[lambda(r_c-r_h)/2], exactly the manuscript convention
Omega/[Lambda_tilde(x_c-x_h)] with x=r/(2M). The calculation demonstrates
physical purely imaginary poles at .11, .111, and .1111 despite their absence
from the original fixed-resolution retained-eigenvalue table. The imaginary
branch has growing |varpi| over these increasingly near-extremal samples and
does not follow the finite-varpi light-ring continuation. No theorem about
its exact limiting Omega or globally unique branch topology follows from
these finite samples. Small modulations near integer horizon-index crossings
are retained under precision/order/matching checks.

## Reproduction

Run, from this directory:

    python ds_polynomial_check.py
    python ds_continuation.py
    python ds_extend_continuation.py

The outputs are `ds_continuation_results.json`, `ds_extended_results.json`,
and `ds_selected.csv`, with full roots, matching values, derivatives, geometry,
index distances, and numerical settings. `ds_low_lambda_continuation.py`
provides an auxiliary backward check of the de Sitter identification; its
explicit local brackets are split at all Frobenius-index poles.

The auxiliary backward calculation completed 20 points from .058 to .02 in
steps of -.002, starting from the .06 branch. It uses 100 decimal digits,
series orders 468 to 1387 (increased to keep the shrinking overlap resolved),
and check orders 200 higher. The raw radiation determinant is solved by a
safeguarded Ridder search in predictor +/- .001, splitting at every
Frobenius-index pole and requiring exactly one accepted root in the local
bracket. At .02 it gives alpha=.16326361316453362635981838,
alpha/[2sqrt(lambda/3)]=.999781364540614, and higher-order normalized W
=-1.19e-67. This is additional evidence for the de Sitter association, not a
proof of a unique global family label or of the Lambda=0 limiting problem.

`ds_fixed_root_shooting_check.json` records a separate original-master ODE
shooting implementation at the fixed .06 root. Its refined double-precision
normalized determinants are 4.48e-11 and 4.90e-11 in distinct endpoint and
matching configurations. Perturbing alpha by +/-1e-6 gives opposite signs,
approximately -1.4870e-5 and +1.4871e-5. This independently corroborates the
root without adjusting it, but is weaker precision evidence than the
arbitrary-precision series checks.
