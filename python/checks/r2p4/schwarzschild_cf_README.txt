Independent Schwarzschild low-overtone continued-fraction computation

Scope
This is a new Python calculation prepared during the referee revision. It is
not a reconstruction of the historical Maple/Julia calculations or provenance
for their original Schwarzschild entries. It computes only the 19 specified
positive-real, damped low-overtone frequencies. It makes no classification of
the scalar negative-imaginary-axis candidates formerly in the manuscript.

Primary formulation
E. W. Leaver, "An analytic representation for the quasi-normal modes of Kerr
black holes", Proc. R. Soc. Lond. A 402 (1985), 285-298,
DOI: 10.1098/rspa.1985.0119, especially Eqs. (3)-(14).
Primary full-text PDF:
https://blackholes.tecnico.ulisboa.pt/gritting/pdf/black_holes/Leaver_An-Analytic-representation-for-the-quasi-normal-modes-of-Kerr-black-holes.pdf
The DOI metadata was cross-checked against bibliographic records; the publisher
and DOI landing pages returned HTTP403 during this review.

Leaver's own clarification, "Remarks on the continued-fraction method for
computing black-hole quasinormal frequencies and modes" (received 17 October 1991),
Eqs. (1)-(9), distinguishes convergence of the continued fraction from its root
condition and discusses the minimal solution and radiation conditions:
https://www.edleaver.com/Misc/EdLeaver/Publications/RemarksOnContinuedFractionMethodForComputingBlackHoleQuasinormalFrequenciesAndModes.pdf

Definitions
x=r/(2M), Omega=M*omega_phys, rho=-2i*Omega, z=1-1/x.
The potential in dimensionless tortoise coordinate x_*=x+ln(x-1) is
V=(1-1/x)[ell(ell+1)/x^2+(1-s^2)/x^3].
The solver uses
psi=(x-1)^rho*x^(-2rho)*exp[-rho*(x-1)] sum_k a_k z^k.
Equivalently, multiply by a nonzero frequency-dependent constant to write
psi=exp(2i*Omega*x)*(x-1)^(-2i*Omega)*x^(4i*Omega) sum_k a_k z^k.
The prefactor supplies the event-ingoing and infinity-outgoing asymptotics.
The global minimal-sequence condition, and not the prefactor alone, selects
the infinity-outgoing solution by analytic continuation for these off-axis roots.

Recurrence in Omega
alpha_k=(k+1)(k+1-4i*Omega)
beta_k=-2(k-4i*Omega)(k+1-4i*Omega)-ell(ell+1)+s^2-1
gamma_k=(k-4i*Omega)^2-s^2
alpha_0*a_1+beta_0*a_0=0;
alpha_k*a_(k+1)+beta_k*a_k+gamma_k*a_(k-1)=0 for k>=1.
The defining fraction is beta_0-alpha_0*gamma_1/(beta_1-alpha_1*gamma_2/(...))=0.
For root n the solver uses the equivalent n-times-inverted fraction, then
checks the (n+1)-times-inverted fraction. Overtone labels follow the requested
low-mode branches and rough starting locations; no exhaustive search is claimed.

Execution
Dependencies available in this run: Python 3.12, mpmath 1.3.0, SymPy 1.14.0.
No qnm package was installed or used.
Run, from any directory:
  python /path/to/schwarzschild_recurrence_check.py
  python /path/to/schwarzschild_tail_derive.py
  python /path/to/schwarzschild_cf.py
The solver writes results next to its source. An optional --mode INDEX selects
one zero-based member of the 19-seed list. Full run elapsed time was 45.0 seconds.

Verification
The recurrence was derived independently by exact symbolic substitution into
the original Schwarzschild master equation. Three coefficient residuals vanish.
The asymptotic minimal tail is expanded through N^(-5/2); its five coefficients
are derived directly from the recurrence and compared with the solver's code.
Depths 256, 512, 1024 were evaluated at 40 decimal digits. Depth 1024 was repeated
at 70 digits, with tail order 3 instead of 5, and with the adjacent inversion.
Largest differences over all 19 modes:
  depth 512 -> 1024: 2.3409793649490e-12
  precision 40 -> 70: 1.1249982375737e-47
  tail order 5 -> 3 at 1024: 1.2727018954196e-11
  inversion n -> n+1 at 1024: 9.6555299771687e-78
Largest final finite-fraction residual: 4.7834647915028e-70.
These are numerical convergence diagnostics, not rigorous bounds on the
infinite-fraction roots. Four-decimal rounding is stable in all comparisons.
The extensive stored digits should not be interpreted as an accuracy claim.

Outputs
schwarzschild_cf_results.json: complete per-mode convergence records.
schwarzschild_cf_table_cells.json: table-cell mapping and maxima.
schwarzschild_cf.log: full numerical run transcript.
schwarzschild_recurrence_check.{py,json,log}: exact master-equation check.
schwarzschild_tail_derive.{py,json,log}: exact tail derivation/code check.

At four decimals 17 values agree with the original oscillatory table entries.
The two changed values are gravitational ell=2:
  n=0: 0.3737-0.0890i
  n=1: 0.3467-0.2739i
