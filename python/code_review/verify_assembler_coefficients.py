#!/usr/bin/env python3
"""Exact source-level checks of the supplied Maple matrix assemblers.

Consumes the decoded active inputs produced by parser/decode_maple.py.
Does not execute Maple, run a spectral calculation, or modify user inputs.
All comparisons use symbolic rational expressions, with no floating-point
sampling. The dimensional gravitational potentials are represented by their
dimensionless y expressions in the manuscript.
"""

from pathlib import Path
import hashlib
import json
import re
import sympy as sp


HERE = Path(__file__).resolve().parent
EXTRACTED = HERE / "extracted"
names = "y zc z0 xh Lambda alphac alphah epsilon L c As"
symbols = sp.symbols(names)
namespace = dict(zip(names.split(), symbols))
y, zc, z0, xh, lam, ac, ah, epsilon, ell, c, As = symbols

delta = zc - 1
ze = (zc + 1 - 2 * z0) / delta
zt = (zc + 1) / delta
f_expected = lam * xh**2 * delta**2 / 3 * (1 - y**2) * (y + ze) / (y + zt)
fp_expected = sp.diff(f_expected, y)
bracket = (ac - ah) * y + ac + ah
quadratic = (ac - ah) * y**2 + 2 * (ac + ah) * y + ac - ah

v_scalar = (
    4 * f_expected / (delta**2 * (y + zt))
    * (ell * (ell + 1) / (y + zt) + epsilon * fp_expected)
)
v_odd = (
    4 * f_expected / (delta**2 * (y + zt)**2)
    * (ell * (ell + 1) - 6 / (xh * delta * (y + zt)))
)
v_even = (
    8 * f_expected
    * (
        xh**3 * delta**3 * (c**3 + c**2 - 3 * lam) * (y + zt)**3
        + 3 * c**2 * xh**2 * delta**2 * (y + zt)**2
        + 9 * c * xh * delta * (y + zt) + 9
    )
    / (xh * delta**3 * (y + zt)**3 * (c * xh * delta * (y + zt) + 3)**2)
)

common_expected = {
    "L01": delta**2 / 4 * f_expected * fp_expected,
    "L02": delta**2 / 4 * f_expected**2,
    "L10": -delta**2 / 4 * (
        f_expected * fp_expected / (1 - y**2) * bracket
        + f_expected**2 / (1 - y**2)**2 * quadratic
    ),
    "L11": -delta**2 / 2 * f_expected**2 / (1 - y**2) * bracket,
    "L12": sp.S.Zero,
    "L20": xh**2 * delta**4 / 4
    - delta**2 / 4 * f_expected**2 / (1 - y**2)**2 * bracket**2,
    "L21": sp.S.Zero,
    "L22": sp.S.Zero,
}
nariai_expected = {
    "L00": -As, "L01": -2 * y, "L02": 1 - y**2,
    "L10": sp.S.One, "L11": 2 * y, "L12": sp.S.Zero,
    "L20": sp.S.One, "L21": sp.S.Zero, "L22": sp.S.Zero,
}

checks = []
sources = []


def parse_definitions(path):
    raw = path.read_text()
    # Cached output is evidence of saved output only, never executable input.
    code = raw.split("===== Group 0 Input Equation 0 =====", 1)[1]
    code = code.split("===== Group 0 Output", 1)[0]
    values = {}
    pattern = re.compile(r"\s*(L\d\d|f|fp)\s*:=\s*y\s*->\s*(.*):$")
    for line in code.splitlines():
        match = pattern.fullmatch(line)
        if match:
            name, rhs = match.groups()
            if name in values:
                raise AssertionError(f"Duplicate active definition: {path.name}: {name}")
            # The decoder preserves the Maple display entity for Lambda.
            rhs = rhs.replace("&Lambda;", "Lambda").replace("^", "**")
            values[name] = sp.sympify(rhs, locals=namespace)
    sources.append({
        "file": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    })
    return values


def check(sector, name, actual, expected):
    residual = sp.factor(actual - expected)
    passed = residual == 0
    result = {
        "sector": sector, "check": name,
        "passed": passed, "exact_residual": str(residual),
    }
    checks.append(result)
    print(f"{sector}: {name}: {'PASS' if passed else 'FAIL'}", flush=True)
    if not passed:
        print(f"  residual = {residual}", flush=True)


for sector, potential in [("s0_s1", v_scalar), ("s2odd", v_odd), ("s2even", v_even)]:
    path = EXTRACTED / f"Matrix_Assembler_{sector}-Final Version.decoded.txt"
    actual = parse_definitions(path)
    assert len([name for name in actual if name.startswith("L")]) == 9
    expected = dict(common_expected)
    expected["L00"] = -delta**4 / 16 * potential
    for name in sorted(expected):
        check(sector, name, actual[name], expected[name])
    check(sector, "f_y_factorization", actual["f"], f_expected)
    check(sector, "fp_y_derivative", actual["fp"], fp_expected)

    # Impose exact horizon relations using zc>1 as the parameter. This
    # establishes that the coded f(y) is the original SdS metric function
    # evaluated at z(y), rather than merely an internally consistent factor.
    horizon_poly = zc**2 + zc + 1
    h = horizon_poly / (zc * (zc + 1))
    cosmological = 3 / (4 * h**2 * horizon_poly)
    z = (zc + 1 + delta * y) / 2
    substitutions = {xh: h, lam: cosmological, z0: -(1 + zc)}
    exact_f = actual["f"].subs(substitutions, simultaneous=True)
    original_f = 1 - 1 / (h * z) - sp.Rational(4, 3) * cosmological * h**2 * z**2
    check(sector, "f_y_equals_original_SdS", exact_f, original_f)

actual_nariai = parse_definitions(
    EXTRACTED / "Matrix_Assembler_all_s_Nariai-Final Version.decoded.txt"
)
assert len(actual_nariai) == 9
for name in sorted(nariai_expected):
    check("Nariai", name, actual_nariai[name], nariai_expected[name])

report = {
    "scope": "Exact source-level identities only; no Maple or spectral numerical execution.",
    "sympy_version": sp.__version__,
    "source_files": sources,
    "passed": sum(row["passed"] for row in checks),
    "total": len(checks),
    "coefficient_checks": 36,
    "metric_and_derivative_checks": 9,
    "checks": checks,
}
result_path = HERE / "assembler_coefficient_checks.json"
result_path.write_text(json.dumps(report, indent=2) + "\n")
print(f"{report['passed']}/{report['total']} exact checks passed.", flush=True)
print(result_path, flush=True)
raise SystemExit(0 if report["passed"] == report["total"] else 1)
