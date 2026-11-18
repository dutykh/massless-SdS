Independent axial/polar comparison for Referee 2, point 5

Scope
This is a NEW Python calculation made for the revision. It is not a rerun of
the historical Maple/Julia pipeline. It tests ell=2, n=0,1 at
Lambda_tilde=0.06,0.09,0.11. It makes no assertion about purely imaginary
roots, the algebraically special candidate, all overtones, or completeness.

Direct potentials and normalization
The physical master equation is
  f d_r(f d_r Psi)+(omega^2-V)Psi=0, exp(-i omega t).
M=1 and Omega=M omega. Each parity uses its own direct potential:
  Vodd/f = ell(ell+1)/r^2 - 6/r^3;
  Veven/f = 2[9+3c^2 r^2+c^2(1+c)r^3+3(3cr-Lambda r^3)]
            /[r^3(3+cr)^2], c=ell(ell+1)/2-1.
No Darboux transformation or axial eigenfunction is used to build the polar
potential, matrices, eigenfunctions or shooting data.

Set x=r/2, a=(xh+xc)/2, d=xc-xh, x=a+d y/2 and x0=-2a.
Then f=Lambda*d^2*g, g=(1-y^2)(x-x0)/(3x), and
  varpi=Omega/(Lambda*d).
After extracting the two radiation powers, the regular equation is
  g Phi''+(g'+2i varpi G)Phi'
    +[-Vbar+i varpi G'+varpi^2 R]Phi=0,
where Vbar=(V/f)/Lambda and the geometry-only G,R are implemented with
their removable endpoint cancellations performed analytically.
Equivalently Ueven=4(Veven/f) in the x equation, so
  Vbar_even=Ueven/(4Lambda).
parity_symbolic_check.py independently verifies this conversion exactly and
verifies a rational horizon-series division recurrence through fourth order.

Spectral calculation
Each parity is separately collocated at N=40,60,80,100 interior Chebyshev
roots with N coefficients (polynomial degree N-1). The same nonzero row
scaling is applied to all three coefficient matrices. No cross-parity
matrix construction, frequency seed or parity-enforcing matching is used.
For each parity, background and resolution, select the two least-damped
positive-real roots inside
  0.5<Re(varpi)<4, -3.2<Im(varpi)<-0.05.
Compare these roots across resolutions. The separate Jost check identifies
them as the desired off-axis physical roots in the stated finite set.

Radiation check
Both parities are separately solved in the original physical master
equation, using t=(r-rh)/(rc-rh). The endpoint coordinates are q=t at the
event horizon and q=1-t at the cosmological horizon. In the corresponding
g=f/[Lambda*(rc-rh)^2/2] equation, the radiation exponent is
  p=-i varpi/g_q(0)
at both endpoints. Each local series has leading coefficient one. The polar
rational potential is expanded directly in q. Its series uses the actual
right-endpoint sign r=rc-(rc-rh)q. The matching determinant includes
d/dt=-d/dq for the right solution.
Each parity is independently initialized by sqrt(15)/2-i(n+1/2) at .11,
then its own shooting roots initialize .09 and .06. No spectral root or
root from the other parity is used as a shooting seed.

Three settings probe shooting accuracy:
  coarse: order24, endpoint offset .04, matching t=.43,
          DOP853 rtol2e-12, atol2e-13;
  refined: order32, offset .025, matching t=.57,
           rtol8e-14, atol8e-15;
  third: order36, offset .06, matching t=.37, refined tolerances;
         evaluates the refined root without fitting another root.
All twelve coarse and twelve refined root searches returned success.
The horizon exponent differences are noninteger because these frequencies
have nonzero real parts. No exceptional Frobenius normalization is involved.

Results, in Omega units
All six axial/polar pairs pass the stored acceptance tests.
Maximum |Omega_odd-Omega_even| at N=100: 8.3352823e-11.
Maximum odd/even resolution spreads: 4.0836749e-11 / 4.5524561e-11.
Maximum spectral-to-Jost frequency difference: 4.3178244e-11.
Maximum independent Jost odd/even frequency difference: 1.1406550e-12.
Maximum coarse/refined Jost frequency difference: 2.0811299e-11.
The maximum parity discrepancy occurs for Lambda=.06, n=1; the separate
resolution and shooting differences show the same binary64 sensitivity.
All eight-decimal table frequencies agree between parities and methods.

Other diagnostics
Maximum normalized matrix-pencil residual: 1.8275727e-16.
Maximum globally term-normalized off-grid/endpoint ODE residual: 9.3666766e-11.
Maximum highest-eight-coefficient norm fraction: 1.8319054e-16.
Maximum refined normalized Jost determinant: 1.8960336e-13.
Maximum third-setting normalized determinant: 7.1279929e-12.
Maximum original-ODE residual for endpoint series data: 2.7477922e-16.
Normalized Wronskians divide by the product of the Euclidean norms of
(Psi,dPsi/dt) at the matching point. These diagnostics are not absolute
Green-function residues, rigorous frequency-error bounds or proof of
spectral completeness. Resolution agreement does not establish measured
exponential convergence because roundoff is already significant.

Reproduction
Python, NumPy, SciPy, threadpoolctl and SymPy are required.
  OPENBLAS_NUM_THREADS=1 python parity_check.py
  python parity_symbolic_check.py
The numerical program is self-contained; it imports no neighboring scripts
and reads no manuscript frequencies, historical matrices or stored results.
It writes parity_results.json and parity_results.csv next to itself.
The executed run took approximately 11 seconds in this environment.

Files
parity_check.py: complete spectral and original-equation shooting calculation.
parity_results.json: full roots, settings, diagnostics and six-row comparison.
parity_results.csv: comparison table in Omega units.
parity_run.log: actual numerical run transcript.
parity_symbolic_check.py/.log: independent exact normalization/series checks.
parity_table.tex: proposed small table fragment; manuscript was not edited.
