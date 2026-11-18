#!/usr/bin/env python3
"""Exact SdS Regge--Wheeler/Zerilli factorization and special-solution checks.

No numerical eigenvector is used or inferred. The script does not classify
frequency-dependent Jost limits at integer Frobenius differences; that
exception is explicitly retained in the output.
"""

from pathlib import Path
import json
import sympy as sp
import mpmath as mp

HERE = Path(__file__).resolve().parent
r, M, lam, c = sp.symbols("r M Lambda c", positive=True)
ell = sp.symbols("ell", integer=True, positive=True)
omega = sp.symbols("omega")
rh, rc = sp.symbols("rh rc", positive=True)
f = 1 - 2 * M / r - lam * r**2 / 3
k = c * (c + 1) / (3 * M)
q = c + 3 * M / r
W = k + 3 * M * f / (r * (c * r + 3 * M))
Vodd = f * (2 * (c + 1) / r**2 - 6 * M / r**3)
Veven = (
    2 * f
    * (9 * M**3 + 3 * c**2 * M * r**2 + c**2 * (1 + c) * r**3
       + 3 * M**2 * (3 * c * r - lam * r**3))
    / (r**3 * (3 * M + c * r)**2)
)


def D(expr):
    return f * sp.diff(expr, r)


checks = []


def check(name, expression):
    residual = sp.factor(sp.cancel(expression))
    passed = residual == 0
    checks.append({"name": name, "passed": passed, "residual": str(residual)})
    print(f"{name}: {'PASS' if passed else 'FAIL'}", flush=True)
    if not passed:
        print(residual, flush=True)


check("ell_form_of_special_frequency", k.subs(c, (ell - 1) * (ell + 2) / 2)
      - (ell - 1) * ell * (ell + 1) * (ell + 2) / (12 * M))
check("odd_superpotential_identity", Vodd - (W**2 - D(W) - k**2))
check("even_superpotential_identity", Veven - (W**2 + D(W) - k**2))
check("axial_kernel_D_plus_W", D(q) - k * q + W * q)
check("polar_kernel_minus_D_plus_W", -D(1 / q) - k / q + W / q)
check("axial_AS_reduced_wave_equation", D(D(q)) - 2 * k * D(q) - Vodd * q)
check("polar_AS_reduced_wave_equation", D(D(1 / q)) + 2 * k * D(1 / q) - Veven / q)

phi = sp.Function("phi")(r)
opA = lambda expression: D(expression) + W * expression
opB = lambda expression: -D(expression) + W * expression
Hodd = lambda expression: -D(D(expression)) + Vodd * expression
Heven = lambda expression: -D(D(expression)) + Veven * expression
check("factorization_BA", opB(opA(phi)) - Hodd(phi) - k**2 * phi)
check("factorization_AB", opA(opB(phi)) - Heven(phi) - k**2 * phi)
check("Darboux_intertwining", opA(Hodd(phi)) - Heven(opA(phi)))
check("inverse_amplitude_product", (k - sp.I * omega) * (k + sp.I * omega) - (omega**2 + k**2))

# Wronskian transformation on any two odd-sector solutions of frequency omega.
p, p1, s, s1, ww, kk = sp.symbols("p p1 s s1 W k")
Ap = p1 + ww * p
As = s1 + ww * s
DAp = (ww**2 - kk**2 - omega**2) * p + ww * p1
DAs = (ww**2 - kk**2 - omega**2) * s + ww * s1
check("Darboux_Wronskian_factor", Ap * DAs - DAp * As - (omega**2 + kk**2) * (p * s1 - p1 * s))

# Exact horizon derivatives and exponent differences.
lambda_h = 3 * (rh - 2 * M) / rh**3
lambda_c = 3 * (rc - 2 * M) / rc**3
kappa_h = (3 * M - rh) / rh**2
kappa_c = (rc - 3 * M) / rc**2
check("event_surface_gravity", sp.diff(f, r).subs({r: rh, lam: lambda_h}) / 2 - kappa_h)
check("cosmological_surface_gravity", -sp.diff(f, r).subs({r: rc, lam: lambda_c}) / 2 - kappa_c)
check("event_W_limit", W.subs({r: rh, lam: lambda_h}) - k)
check("cosmological_W_limit", W.subs({r: rc, lam: lambda_c}) - k)

mass_horizons = rh * rc * (rh + rc) / (2 * (rh**2 + rh * rc + rc**2))
polyell = 4 * c * (c + 1)
nu_h_radii = polyell * (rh**2 + rh * rc + rc**2)**2 / (3 * rc * (rc - rh) * (rc + rh) * (rc + 2 * rh))
nu_c_radii = polyell * (rh**2 + rh * rc + rc**2)**2 / (3 * rh * (rc - rh) * (rc + rh) * (2 * rc + rh))
check("event_integer_index_formula", (k / kappa_h).subs(M, mass_horizons) - nu_h_radii)
check("cosmological_integer_index_formula", (k / kappa_c).subs(M, mass_horizons) - nu_c_radii)

# A genuine exceptional example, not an assumption about the sampled cases.
h_example = 3 * M * (sp.sqrt(3) - 1)
c_example = 6 * M
lambda_example = 1 / (18 * M**2)
check("exceptional_event_horizon", f.subs({r: h_example, lam: lambda_example}))
check("exceptional_cosmological_horizon", f.subs({r: c_example, lam: lambda_example}))
check("exceptional_ell2_event_index_12", (k / kappa_h).subs({c: 2, rh: h_example}) - 12)
check("exceptional_ell2_cosmological_index_24", (k / kappa_c).subs({c: 2, rc: c_example}) - 24)

# Leading-amplitude Wronskian for nonresonant endpoints. This is explicitly
# conditional on the usual finite, analytically continued Jost branches.
u = sp.symbols("u", real=True)
qh = c + 3 * M / rh
qc = c + 3 * M / rc
left_at_right = qc / qh * sp.exp(-k * u)
right_out = sp.exp(k * u)
check("nonresonant_AS_Jost_Wronskian", left_at_right * sp.diff(right_out, u)
      - sp.diff(left_at_right, u) * right_out - 2 * k * qc / qh)

# Parameter checks do not reconstruct, inspect, or classify an eigenvector.
mp.mp.dps = 70
parameter_rows = []
for angular in (2, 3):
    cc = mp.mpf((angular - 1) * (angular + 2)) / 2
    special_k = cc * (cc + 1) / 3  # M=1
    for lambda_text in ("0.02", "0.04", "0.06", "0.08", "0.09", "0.10", "0.11"):
        value = mp.mpf(lambda_text)
        alpha = mp.acos(3 * mp.sqrt(value))
        horizon_h = 2 * mp.cos((mp.pi + alpha) / 3) / mp.sqrt(value)
        horizon_c = 2 * mp.cos((mp.pi - alpha) / 3) / mp.sqrt(value)
        gravity_h = (3 - horizon_h) / horizon_h**2
        gravity_c = (horizon_c - 3) / horizon_c**2
        nh, nc = special_k / gravity_h, special_k / gravity_c
        det = 2 * special_k * (cc + 3 / horizon_c) / (cc + 3 / horizon_h)
        parameter_rows.append({
            "ell": angular, "M": "1", "M2Lambda": lambda_text,
            "nu_h": mp.nstr(nh, 40), "nu_c": mp.nstr(nc, 40),
            "nu_h_distance_to_nearest_integer": mp.nstr(abs(nh - mp.nint(nh)), 25),
            "nu_c_distance_to_nearest_integer": mp.nstr(abs(nc - mp.nint(nc)), 25),
            "conditional_nonresonant_Jost_W_at_AS": mp.nstr(det, 40),
        })

report = {
    "sympy_version": sp.__version__, "mpmath_version": mp.__version__,
    "passed": sum(row["passed"] for row in checks), "total": len(checks),
    "checks": checks,
    "parameter_checks": parameter_rows,
    "scope": {
        "exact_solution": "Axial kernel (c+3M/r) exp(-k rstar), polar reciprocal kernel.",
        "general_second_solution": "psi2=psiAS*integral(psiAS^(-2) drstar).",
        "generic_classification": "For noninteger endpoint exponent differences, axial kernel is BH ingoing/cosmological incoming, and Wodd(-ik)=2k qc/qh is nonzero.",
        "integer_index_exception": "At k/kappa_h or k/kappa_c integer, analytic continuation of the frequency-dependent Jost solutions is required. Leading exponent or kernel existence alone does not establish its classification.",
        "numerical_eigenvectors": "Not supplied, constructed, or inferred in this verification.",
    },
}
target = HERE / "sds_algebraically_special_checks.json"
target.write_text(json.dumps(report, indent=2) + "\n")
print(f"{report['passed']}/{report['total']} exact checks passed.", flush=True)
print(target, flush=True)
raise SystemExit(0 if report["passed"] == report["total"] else 1)
