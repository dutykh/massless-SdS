Targeted higher-order multiprecision checks for Referee 3, point 2

This is a new bounded calculation for the scalar ell=2 purely imaginary
branch at Lambda_tilde=0.06. It uses the same fully rebuilt 70-decimal-digit
quadratic pencil and bordered Newton solver as the preceding precision study,
now at N=240 and 300. The imported original file
checks/r2p5/imaginary_multiprecision_newton.py is unchanged.

Run from the workspace/archive root in separate processes if desired:

  python -u checks/r3p2/high_order_multiprecision.py --order 240 --dps 70
  python -u checks/r3p2/high_order_multiprecision.py --order 300 --dps 70

Requires mpmath 1.3.0 (the recorded environment uses Python 3.12.14). The preceding
N128 spectral result and independent reference file must remain available at
checks/r2p5/imaginary_multiprecision_results.json and
checks/r2p5/ds_reference_eigenfunction.json. The wrapper writes only new results
under checks/r3p2; old results, the original solver, and manuscript are untouched.

Explicit change of initializer
The original smaller-N runs used binary64 spectral seeds. These two higher-N
runs use the previously converged N128, 70-digit SPECTRAL frequency and its 128
Chebyshev coefficients, converted directly to 70-digit mpmath mpf values.
The eigenvector is extended with 112 or 172 exactly zero coefficients, and all
imaginary coefficient entries in the seed are exactly zero. No independent
Jost/Frobenius frequency or field is used to initialize, select or solve the
new eigenpair. The inherited hard-coded 'binary64 spectral candidate' output
label is explicitly replaced in the new result and retained as an overridden
historical label in metadata. No original solver source is silently modified.

Full rebuild and actual correction
Every collocation node, horizon, geometry coefficient, Chebyshev polynomial
and derivative, and pencil entry is regenerated at 70 decimal digits for the new
N. The real imaginary-frequency pencil is Q(alpha)=A0+alpha*B1-alpha^2*A2,
where varpi=-i*alpha. The normalized eigenvector satisfies Phi(-0.3)=1.
The original bordered Newton system updates all N coefficients and alpha.
It always computes and applies an actual Newton correction before testing
whether the absolute frequency correction and relative vector correction are
both below 1e-50. The existing neighborhood guard of 0.05 in varpi is retained.
A tiny seed defect alone is not used as a stopping test.

Post-solve diagnostics
The independent Frobenius reference is lazily read only when the unchanged
diagnostics function first accesses it after Newton has converged. That read
time and each actual correction time are recorded. The original 1025-point
common Chebyshev-grid diagnostics, endpoint tests, independently reconstructed
unfactored-original-ODE residual, pencil residual, last-eight coefficient
fraction, and independent frequency/eigenfunction comparisons are retained.
Full coefficients, 1025 field values, runtime, source hashes and Newton history
are saved for each N. The independent reference is a numerical comparator,
not an exact solution or a rigorous error bound.

Interpretation
The N128 spectral solution is already much more accurate than the Newton
stopping threshold. The purpose of the higher-N rebuild is to test stability
at these larger orders and fixed 70-digit arithmetic; a small deterioration
from the best smaller-N result must be reported if observed. Small residuals
or coefficient tails alone do not certify frequency accuracy. These two
checks do not constitute a search for every mode or establish completeness.
They do not reproduce the historical Maple/Julia campaign.

Observed completed runs
Both larger orders performed one actual Newton correction before satisfying
the original stopping rule. N240 took 43.55 seconds and N300 took 77.05 seconds,
including all post-solve diagnostics. The independent varpi differences were
4.4638843e-61 and 2.5205223e-60, respectively; relative eigenfunction differences
were 1.4932515e-61 and 8.4316117e-61. The preceding N128 differences were
4.2451793e-63 and 1.4200928e-63. Thus increasing N at fixed 70-digit arithmetic
slightly degrades this already very accurate result; monotone improvement is
not claimed. The observed original-ODE residual is about 3.354e-67 at both
larger orders. All results remain numerical diagnostics, not certified bounds.

To regenerate only the summary and verify saved-run structure/source hashes:
  python checks/r3p2/high_order_multiprecision_summarize.py
