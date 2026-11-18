#!/usr/bin/env python3
"""Exact algebraic checks of the manuscript's Nariai potential amplitudes.

The common lapse is cancelled symbolically BEFORE imposing extremality.
These are checks of the finite coefficient r^2 V/f, not evaluations of 0/0.
No numerical frequencies, collocation approximations, or manuscript edits.
"""
from pathlib import Path
import hashlib
import json
import sympy as sp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
SOURCE = ROOT / "revision/main_referee3_point3.tex"
M, r, Lambda, L = sp.symbols("M r Lambda L", positive=True)
c = L / 2 - 1
f = 1 - 2 * M / r - Lambda * r**2 / 3

# Vodd and Veven, as printed in the inspected manuscript.
Vodd = f * (L / r**2 - 6 * M / r**3)
Veven = (
    2 * f * (
        9 * M**3 + 3 * c**2 * M * r**2 + c**2 * (1+c) * r**3
        + 3 * M**2 * (3*c*r - Lambda*r**3)
    ) / (r**3 * (3*M+c*r)**2)
)
Uscalar = f * (sp.diff(f, r)/r + L/r**2)
Umaxwell = f * L/r**2
extremal = {r: 3*M, Lambda: 1/(9*M**2)}

checks = []
for name, potential, expected in (
    ("gravitational_odd", Vodd, L-2),
    ("gravitational_even", Veven, L-2),
    ("massless_minimally_coupled_scalar", Uscalar, L),
    ("source_free_maxwell", Umaxwell, L),
):
    coefficient = sp.cancel(r**2 * potential/f)
    amplitude = sp.factor(coefficient.subs(extremal))
    residual = sp.factor(amplitude-expected)
    checks.append({
        "name": name,
        "coefficient_before_limit": str(coefficient),
        "r0_squared_times_V_over_f": str(amplitude),
        "expected": str(expected),
        "exact_residual": str(residual),
        "passed": residual == 0,
    })
    assert residual == 0, (name, residual)

# An independent expansion of the even numerator at r=3M also gives 2c.
even_at_r0 = sp.factor(sp.cancel(r**2*Veven/f).subs(extremal))
assert sp.factor(even_at_r0-2*c) == 0

output = {
    "method": "SymPy exact rational identities; no floating-point arithmetic",
    "source_inspected": str(SOURCE.relative_to(ROOT)),
    "source_sha256": hashlib.sha256(SOURCE.read_bytes()).hexdigest() if SOURCE.exists() else None,
    "definitions": {
        "L": "ell*(ell+1)", "c": "L/2-1", "r0": "3*M",
        "Lambda_at_extremality": "1/(9*M**2)",
        "physical_sphere_radius": "r0=Lambda**(-1/2)",
    },
    "important_order": "Cancel f before substituting the extremal parameters.",
    "domain": "M>0; radiative gravity ell>=2, hence c>0 and denominators nonzero.",
    "all_passed": all(item["passed"] for item in checks),
    "checks": checks,
}
(HERE / "nariai_potential_limits_results.json").write_text(
    json.dumps(output, indent=2) + "\n"
)
for item in checks:
    print(f'{item["name"]}: amplitude={item["r0_squared_times_V_over_f"]}; '
          f'residual={item["exact_residual"]}; passed={item["passed"]}')
