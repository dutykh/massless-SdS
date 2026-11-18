Representative SdS resolution, eigenfunction and precision checks
================================================================

Scope
-----
These are new, bounded calculations for the fifth comment of Referee 2. They
do not rerun the historical Maple/Julia campaign or classify every spectral
root. Three previously radiation-checked representatives are examined:

  scalar ell=2, Lambda_tilde=.09, positive-real oscillatory n=1;
  axial  ell=2, Lambda_tilde=.11, positive-real oscillatory n=0;
  scalar ell=2, Lambda_tilde=.06, least-damped imaginary spectral candidate.

All three are computed at N=12,16,24,32,48,64,96,128 in binary64. Only the last
case is additionally solved at 70 decimal digits for N=24,48,96,128. This is
not a claim that 70 digits or any fixed precision cures boundary nonselection.

Equation and normalization
--------------------------
The direct metric-derived rescaled equation from ../r2p4 is used, with
  varpi=Omega/[Lambda_tilde*(xc-xh)], Omega=M*omega, x=r/(2M).
It factors the event-ingoing and cosmological-outgoing powers and analytically
cancels their removable coefficient factors. The residual equation is
  g Phi''+(g'+2i varpi G)Phi'+(-Vbar+i varpi G'+varpi^2 R)Phi=0.
Its derivation has 24 exact checks in ../r2p4/rescaled_equation_checks.json.

Every residual eigenfunction is normalized by its complex value at y=-.3,
so Phi(-.3)=1 fixes amplitude and phase. This point is safely away from a zero
for all three chosen representatives. Eigenfunctions are compared on the same
1025-point grid y_j=cos(pi*j/1024), including both endpoints. The relative
maximum-norm differences divide by the maximum reference-function magnitude.

Selection independent of reference frequencies
---------------------------------------------
At N=12 the oscillatory cases are selected by damping order in the fixed
window .5<Re(varpi)<4, -3.2<Im(varpi)<-.05. The imaginary case is the least
damped candidate with abs(Re(varpi))<1e-6, -8<Im(varpi)<-.05. At subsequent N,
the selected candidate is the nearest previous spectral root in the same
window. Its displacement must be below .1 and, when a competing candidate
exists, below one quarter of the separation from the nearest competitor.
The independently calculated horizon frequencies are never used for this
selection or as initial seeds. Their differences are calculated afterward.
This is an association of sampled numerical roots, not a completeness proof.

Recorded diagnostics
--------------------
The CSV and JSON record both Omega and varpi, the change from the previous N,
the difference from an independent radiation-selected frequency, the relative
norm of the highest eight Chebyshev coefficients, and normalized eigenfunction
differences from the preceding N, N=128, and independent physical references.
The complete coefficient arrays and normalized eigenfunctions are in JSON.
The N=128 binary64 eigenfunction is explicitly only an interresolution
comparator; it is not presumed to be the most accurate function.

The independent differential residual substitutes the reconstructed function
into the original unfactored SdS master equation after dividing out the
nonzero radiation prefactor. It uses the metric and its original derivative,
not the assembly's reduced G and R formulas. It is evaluated at the 1023
interior common-grid points. The maximum absolute residual is divided by the
maximum sum of the absolute values of its eight unfactored terms. Separate
reduced-equation and endpoint residuals, and a relative algebraic pencil
residual, are also supplied. Each residual tests the equation or discretized
eigenpair; none independently establishes a radiation condition.

Independent references and their limitations
--------------------------------------------
Frequency references are the original-equation Jost calculations saved in
../r2p4/jost_extra_cell_check.json, jost_check_results.json and
ds_extended_results.json. They are numerical references, not exact values.
The scalar oscillatory reference changed by 2.93e-11 in varpi when its shooting
settings were refined; very small spectral-reference differences therefore
cannot be interpreted as certified accuracy below that scale.

offaxis_reference_eigenfunctions.py evaluates the two oscillatory reference
functions using original-equation horizon Frobenius data and DOP853 integration,
with explicit endpoint-series values and direct complex normalization. The
two reference configurations differ in absolute function sup norm by 4.78e-11
(scalar) and 8.24e-14 (axial). Their reference-function sup norms are 18.5733
and 1.03790, respectively, giving RELATIVE configuration differences of
2.57611e-12 (scalar) and 7.93658e-14 (axial). Frequency-reference uncertainty
is additional. The observed RELATIVE scalar function-comparison plateau is
approximately 8e-12 at N24 through N64; it is not an exact spectral-error
measurement. Likewise, axial differences below its reference-variation scale
are not a claim of greater accuracy. All broad-study spectral function errors
are relative and have not changed in this metadata clarification.

Both reference JSON and broad-study JSON/CSV now expose explicitly named
absolute and relative reference-variation fields. In the off-axis reference
JSON, legacy unqualified configuration and nested spectral-comparison fields
retain their original ABSOLUTE values and are labeled accordingly. In the
broad-study JSON/CSV, the legacy unqualified reference-variation field is an
alias of the corrected RELATIVE value. Prefer the explicitly named fields.

ds_reference_eigenfunction.py constructs the imaginary-mode reference directly
from matched physical Frobenius series, without collocation. At 100 digits and
order500 versus 130 digits and order650, reference functions differ by at most
7.62e-96 absolutely, or 2.97e-97 in the normalized relative sup norm. These are
reference-stability checks, not rigorous bounds. Binary64 comparisons round
the common grid and reference values to the arithmetic being tested.

Main observed limitation and precision check
--------------------------------------------
For the imaginary representative, binary64 varpi-reference differences are
approximately 6.46e-6 at N12, 3.86e-8 at N16 and 2.34e-11 at N24. They then
deteriorate to 2.52e-8 at N64, 1.08e-6 at N96 and 1.11e-5 at N128, despite
small coefficient tails and residuals. Its independent relative eigenfunction
differences similarly deteriorate, reaching 3.70e-6 at N128. This is a direct
counterexample to interpreting a tiny residual, a tiny coefficient tail, or
increased N alone as sufficient evidence of frequency accuracy.

imaginary_multiprecision_newton.py rebuilds the entire pencil at 70 digits:
nodes, geometry, metric coefficients, Chebyshev values and both derivatives.
It does not promote binary64 matrices. A bordered Newton solve uses only the
selected binary64 eigenpair as its initial iterate and enforces Phi(-.3)=1.
For a pure-imaginary varpi=-i alpha, the real pencil is
  (A0+alpha B1-alpha^2 A2)a=0.
The Newton Jacobian borders this matrix by (B1-2alpha A2)a and the normalization
row. The independent root and field enter only after Newton has converged.
The Newton correction and vector update must both fall below 1e-50.

At N24,48,96,128 the observed varpi-reference differences are respectively
1.27e-12,1.18e-26,8.10e-56,4.25e-63; the independent relative eigenfunction
differences are 4.24e-13,3.96e-27,2.71e-56,1.42e-63. These results separate
spectral truncation from the binary64 deterioration for this particular
physical mode. They do not make high precision a substitute for branch
selection, an outgoing scattering test, or a rigorous error analysis.

Reproduction and files
----------------------
Dependencies: Python, numpy, scipy, mpmath, threadpoolctl; previous r2p4 scripts
and reference JSON files are intentionally reused, never modified. From the
workspace root:

  OPENBLAS_NUM_THREADS=1 python checks/r2p5/broad_resolution_study.py
  python checks/r2p5/ds_reference_eigenfunction.py
  python checks/r2p5/offaxis_reference_eigenfunctions.py
  OPENBLAS_NUM_THREADS=1 python checks/r2p5/broad_resolution_study.py
  python checks/r2p5/imaginary_multiprecision_newton.py

The first spectral run supplies the common grid for the independent references;
the second incorporates their eigenfunction comparisons. The Newton script
accepts --orders 24 48 96 128 and --dps 70. Existing computed eigenpairs can
have diagnostics recomputed with --diagnostics-only, without repeating Newton.
All final differential and function comparisons use the common 1025-point
grid. JSON retains the precision and full Newton iteration history.

Main outputs:
  broad_resolution_results.csv/.json
  imaginary_multiprecision_results.csv/.json
  ds_reference_eigenfunction.json
  offaxis_reference_eigenfunctions.json
  corresponding Python scripts and run logs.

No manuscript file, supplied Maple worksheet, or supplied Julia program was
modified as part of these checks.
