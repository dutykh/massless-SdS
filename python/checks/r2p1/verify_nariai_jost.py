#!/usr/bin/env python3
"""Exact and independent numerical checks for the Nariai boundary issue.

This script verifies the radiation-selected Jost determinant separately from
the polynomial regularity pencil. It does not modify manuscript or production
numerical files and does not rerun the supplied Maple/Julia calculation.
"""

from pathlib import Path
import json
import sympy as sp
import mpmath as mp


HERE = Path(__file__).resolve().parent
chi, z, u = sp.symbols("chi z u", real=True)
w, A, q = sp.symbols("w A q")
checks = []


def check(name, expression):
    residual = sp.simplify(expression)
    if residual != 0:
        residual = sp.factor(sp.cancel(residual))
    passed = residual == 0
    checks.append({"name": name, "passed": bool(passed), "residual": str(residual)})
    print(f"{name}: {'PASS' if passed else 'FAIL'}", flush=True)
    if not passed:
        print(residual, flush=True)


def regularized(phi, frequency=w, coefficient=A):
    return (
        (1 - chi**2) * sp.diff(phi, chi, 2)
        - 2 * (1 - sp.I * frequency) * chi * sp.diff(phi, chi)
        + (frequency**2 + sp.I * frequency - coefficient) * phi
    )


# The A=0 Jost problem is free on the complete tortoise line.
f_minus = sp.exp(-sp.I * w * u)
f_plus = sp.exp(sp.I * w * u)
W_free = f_minus * sp.diff(f_plus, u) - sp.diff(f_minus, u) * f_plus
check("free_Jost_W_equals_2iw", W_free - 2 * sp.I * w)
check("zero_frequency_constant_solution", regularized(sp.S.One, 0, 0))
check("zero_frequency_second_solution", regularized(sp.atanh(chi), 0, 0))
check("outgoing_reconstruction_plus_base", sp.cosh(u) + sp.sinh(u) - sp.exp(u))
check("outgoing_reconstruction_minus_base", sp.cosh(u) - sp.sinh(u) - sp.exp(-u))

for n in range(1, 7):
    for sign in (1, -1):
        phi = (1 + sign * chi)**n
        check(f"monopole_polynomial_n{n}_sign{sign}", regularized(phi, -sp.I * n, 0))
    check(f"free_W_at_minus_i_{n}", W_free.subs(w, -sp.I * n) - 2 * n)

# Hypergeometric parameters after z=(1+chi)/2. Here q^2=1/4-A.
a = sp.Rational(1, 2) - sp.I * w + q
b = sp.Rational(1, 2) - sp.I * w - q
c = 1 - sp.I * w
check("hypergeom_derivative_coefficient", a + b + 1 - 2 * (1 - sp.I * w))
check("hypergeom_potential_coefficient", a * b - (sp.Rational(1, 4) - q**2 - w**2 - sp.I * w))
check("hypergeom_connection_exponent", c - a - b - sp.I * w)
check("gamma_recurrence", sp.gammasimp(sp.gamma(c) + sp.I * w * sp.gamma(-sp.I * w)))

C_in = sp.gamma(c) * sp.gamma(-sp.I * w) / (sp.gamma(a) * sp.gamma(b))
W_gamma = -2 * sp.gamma(c)**2 / (sp.gamma(a) * sp.gamma(b))
check("Jost_W_equals_2iw_Cin", sp.gammasimp(W_gamma - 2 * sp.I * w * C_in))
check("A0_Cin_cancellation_equals_one", C_in.subs(q, sp.Rational(1, 2)) - 1)
check("A0_gamma_W_cancellation_equals_2iw", sp.gammasimp(W_gamma.subs(q, sp.Rational(1, 2))) - 2 * sp.I * w)

# Radiation-selected solutions for A>0. For physical sectors A>=2, q is
# imaginary. For 0<A<1/4, q lies in (0,1/2); at A=1/4 it vanishes.
# The parameters c below never become nonpositive integers in those domains.
for n in range(0, 6):
    frequency = -sp.I * (n + sp.Rational(1, 2) + q)
    c_n = sp.Rational(1, 2) - n - q
    check(f"positive_A_pole_condition_n{n}", a.subs(w, frequency) + n)
    polynomial_z = sum(
        sp.rf(-n, j) * sp.rf(-n - 2 * q, j)
        / (sp.rf(c_n, j) * sp.factorial(j)) * z**j
        for j in range(n + 1)
    )
    polynomial_chi = polynomial_z.subs(z, (1 + chi) / 2)
    check(
        f"positive_A_exact_outgoing_polynomial_n{n}",
        regularized(polynomial_chi, frequency, sp.Rational(1, 4) - q**2),
    )
    check(f"positive_A_polynomial_left_amplitude_n{n}", polynomial_z.subs(z, 0) - 1)
    endpoint = sp.rf(sp.Rational(1, 2) + q, n) / sp.rf(c_n, n)
    check(f"positive_A_polynomial_right_amplitude_n{n}", polynomial_z.subs(z, 1) - endpoint)

# The exact finite polynomial pencil explains why numerical convergence does
# not cure the scalar-monopole problem. Columns are operator images of chi^k.
finite_pencil_determinants = {}
for N in range(2, 7):
    matrix = sp.zeros(N)
    for k in range(N):
        image = sp.Poly(sp.expand(regularized(chi**k)), chi)
        for row in range(N):
            matrix[row, k] = image.coeff_monomial(chi**row)
    expected_general = sp.prod(w**2 + sp.I * (2 * k + 1) * w - k * (k + 1) - A for k in range(N))
    check(f"finite_pencil_general_determinant_N{N}", matrix.det() - expected_general)
    exact_free = sp.factor(matrix.subs(A, 0).det())
    expected_free = w * (w + sp.I * N) * sp.prod((w + sp.I * n)**2 for n in range(1, N))
    check(f"finite_pencil_spurious_monopole_ladder_N{N}", exact_free - expected_free)
    finite_pencil_determinants[str(N)] = str(exact_free)

# Independent evaluation of the Jost Wronskian from hypergeometric functions
# at a finite tortoise coordinate, compared with the gamma connection formula.
mp.mp.dps = 75
numerical = []
for coefficient, real_w, imaginary_w in [(2, "0.7", "-0.3"), (6, "2.0", "-2.7"), (0, "0.42", "-1.2")]:
    frequency = mp.mpc(real_w, imaginary_w)
    qq = mp.sqrt(mp.mpf("0.25") - coefficient)
    aa = mp.mpf("0.5") - 1j * frequency + qq
    bb = mp.mpf("0.5") - 1j * frequency - qq
    cc = 1 - 1j * frequency

    def jost(tortoise, left=True):
        zz = (1 + mp.tanh(tortoise)) / 2
        arg = zz if left else 1 - zz
        return (zz * (1 - zz))**(-1j * frequency / 2) * mp.hyp2f1(aa, bb, cc, arg)

    point = mp.mpf("0.173")
    direct = jost(point, True) * mp.diff(lambda t: jost(t, False), point)
    direct -= mp.diff(lambda t: jost(t, True), point) * jost(point, False)
    connection = -2 * mp.gamma(cc)**2 * mp.rgamma(aa) * mp.rgamma(bb)
    error = abs(direct - connection) / max(1, abs(connection))
    passed = error < mp.mpf("1e-65")
    numerical.append({
        "A": coefficient, "w": str(frequency),
        "relative_error": mp.nstr(error, 12), "passed": bool(passed),
    })
    print(f"numerical_Jost_A{coefficient}: {'PASS' if passed else 'FAIL'} ({mp.nstr(error, 5)})", flush=True)

report = {
    "sympy_version": sp.__version__,
    "mpmath_version": mp.__version__,
    "symbolic_passed": sum(row["passed"] for row in checks),
    "symbolic_total": len(checks),
    "symbolic_checks": checks,
    "finite_monopole_pencil_determinants": finite_pencil_determinants,
    "numerical_jost_checks": numerical,
    "interpretation": {
        "A_positive": "Gamma numerator remains finite at the denominator zeros; they are outgoing QNMs. A=1/4 has coalescent branches.",
        "A_zero": "Cancel coincident gamma factors first: C_in=1 and W=2iw. No nonzero damped outgoing resonance exists.",
        "zero_frequency": "Bounded constant static/threshold solution; do not classify by nonzero-frequency incoming/outgoing exponents.",
        "polynomial_ladder": "Exact finite-dimensional regularity eigenvalues persist under increasing N but are not free-wave outgoing resonances.",
    },
}
target = HERE / "nariai_jost_checks.json"
target.write_text(json.dumps(report, indent=2) + "\n")
print(f"{report['symbolic_passed']}/{report['symbolic_total']} symbolic checks passed.", flush=True)
print(target, flush=True)
raise SystemExit(0 if all(r["passed"] for r in checks + numerical) else 1)
