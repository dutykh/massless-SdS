"""Refine two existing scalar ell=2 imaginary SdS connection roots.

No branch extension or collocation computation is performed.  The existing
original-master Frobenius recurrence is re-evaluated with increased precision,
series order, and two new matching locations at each fixed background.
Run from any directory: python checks/r3p3/low_lambda_refinement.py
"""
from pathlib import Path
import hashlib
import json
import platform
import sys
import time

import mpmath as mp

HERE = Path(__file__).resolve().parent
DEPENDENCY = HERE.parent / "r2p4" / "ds_branch_check.py"
INPUT = HERE.parent / "r2p4" / "ds_low_lambda_results.json"
sys.path.insert(0, str(DEPENDENCY.parent))
from ds_branch_check import connection, geometry


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def vector_diagnostics(detail):
    a = mp.mpf(detail["left_value_without_nonzero_power"])
    b = mp.mpf(detail["left_r_derivative_without_nonzero_power"])
    c = mp.mpf(detail["right_value_without_nonzero_power"])
    d = mp.mpf(detail["right_r_derivative_without_nonzero_power"])
    norm = mp.sqrt(a*a+b*b)*mp.sqrt(c*c+d*d)
    return {"raw_determinant_without_nonzero_powers": str(a*d-b*c),
            "vector_norm_product": str(norm),
            "matching_value_product_without_nonzero_powers": str(a*c),
            "both_matching_values_nonzero": bool(a != 0 and c != 0),
            "both_solution_vectors_nonzero": bool(norm > 0)}


def main():
    old = json.loads(INPUT.read_text())
    anchors = {x["Lambda_tilde"]: x for x in old["path"]}
    report = {
        "method": __doc__, "sector": "minimally coupled massless scalar",
        "ell": 2, "M": "1", "time_convention": "exp(-i omega t)",
        "frequency_conventions": {
            "Omega": "M*omega=-i*alpha; M=1",
            "varpi": "Omega/[Lambda_tilde*(r_c-r_h)/2]",
            "scalar_fundamental_dS_ratio": "alpha/[2*sqrt(Lambda_tilde/3)]"},
        "search": "At each fixed Lambda use the archived root as seed; secant seeds alpha +/- 1e-14, tolerance 1e-90, at most 30 iterations; reject displacement >1e-8.",
        "checks": "Re-evaluate at N+300; require |normalized W|<1e-80 and integer-index distances>1e-7. Compare two root reconstructions, matching values and solution norms. Compare central determinant derivatives at h=1e-14 and h/2; require opposite perturbation signs and derivative relative disagreement <1e-14.",
        "diagnostic_step_refinement": "An initial h=1e-10 derivative comparison at Lambda=.02 passed the root/order/index checks but differed by 1.8336e-10 relatively because of finite-difference truncation. The derivative increments were reduced to 1e-14 and 5e-15 before the final complete run; initial failure log is retained separately.",
        "normalization": "Each physical radiation Frobenius solution has unit leading amplitude. The q^rho factors are nonzero at the matching point and are removed only after retaining their logarithmic derivatives. The determinant uses r derivatives; multiplication by f at the interior matching point converts to tortoise derivatives without altering its zero set. No frequency-dependent zero is multiplied into the determinant.",
        "provenance": {
            "source": str(Path(__file__).resolve()), "source_sha256": sha256(Path(__file__)),
            "recurrence_source": str(DEPENDENCY), "recurrence_sha256": sha256(DEPENDENCY),
            "input": str(INPUT), "input_sha256": sha256(INPUT),
            "python": platform.python_version(), "mpmath": mp.__version__},
        "limitations": [
            "Same original-ODE recurrence, independently varied numerical settings; not a second mathematical formulation.",
            "Observed numerical convergence is not an interval-certified error bound.",
            "Only Lambda_tilde=.02 and .04 are refined; no extension below .02.",
            "The de Sitter ratio supports a bounded numerical association, not a limiting theorem or globally unique family assignment.",
            "This revalidates newly constructed radiation solutions, not archived collocation eigenvectors."],
        "anchors": [], "complete": False}
    start = time.time()
    configurations = {
        "0.02": [(120, 3400, ".4"), (140, 3800, ".6")],
        "0.04": [(120, 2200, ".4"), (140, 2600, ".6")]}
    for lam_text, configs in configurations.items():
        original = anchors[lam_text]
        anchor = {"Lambda_tilde": lam_text, "original_auxiliary": original,
                  "refinements": []}
        report["anchors"].append(anchor)
        for dps, order, factor_text in configs:
            config_start = time.time()
            with mp.workdps(dps):
                lam = mp.mpf(lam_text)
                rh, rc = geometry(lam)
                t = min(mp.mpf(".45"), mp.mpf(factor_text)*rh/(rc-rh))
                seed = mp.mpf(original["alpha"])
                fn = lambda alpha: connection(lam, alpha, order, match_fraction=t)
                eps = mp.mpf("1e-14")
                root = mp.findroot(fn, (seed-eps, seed+eps),
                                   tol=mp.mpf("1e-90"), maxsteps=30, verify=True)
                if abs(root-seed) > mp.mpf("1e-8"):
                    raise RuntimeError("Refinement escaped original root neighborhood")
                value, detail = connection(lam, root, order, details=True, match_fraction=t)
                check, check_detail = connection(lam, root, order+300, details=True, match_fraction=t)
                index_distance = min(abs(mp.mpf(detail[key])-mp.nint(mp.mpf(detail[key])))
                                     for key in ("nu_h", "nu_c"))
                derivs, perturbations = [], []
                for step in (mp.mpf("1e-14"), mp.mpf("5e-15")):
                    minus = connection(lam, root-step, order+300, match_fraction=t)
                    plus = connection(lam, root+step, order+300, match_fraction=t)
                    derivs.append((plus-minus)/(2*step))
                    perturbations.append({"step": str(step), "D_minus": str(minus),
                                          "D_plus": str(plus),
                                          "opposite_signs": bool(minus*plus < 0)})
                derivative_relative_change = abs(derivs[0]-derivs[1])/abs(derivs[1])
                diagnostics = vector_diagnostics(check_detail)
                if not (abs(check) < mp.mpf("1e-80") and index_distance > mp.mpf("1e-7")
                        and all(x["opposite_signs"] for x in perturbations)
                        and derivative_relative_change < mp.mpf("1e-14")
                        and diagnostics["both_matching_values_nonzero"]):
                    raise RuntimeError(("Acceptance failed", lam_text, dps, str(check),
                                        str(index_distance), str(derivative_relative_change)))
                item = {"dps": dps, "N": order, "check_N": order+300,
                        "match_factor": factor_text, "matching_t": str(t),
                        "alpha": str(root), "difference_from_auxiliary": str(root-seed),
                        "varpi_im": str(-root/(lam*(rc-rh)/2)),
                        "ratio_to_2sqrtLambdaOver3": str(root/(2*mp.sqrt(lam/3))),
                        "normalized_W": str(value), "normalized_W_check": str(check),
                        "minimum_integer_index_distance": str(index_distance),
                        "derivative_estimates": list(map(str, derivs)),
                        "derivative_relative_change": str(derivative_relative_change),
                        "symmetric_perturbation_checks": perturbations,
                        "geometry_and_series": detail, "higher_order_series": check_detail,
                        "nonzero_solution_checks": diagnostics,
                        "elapsed_seconds": time.time()-config_start}
                anchor["refinements"].append(item)
                print("Lambda", lam_text, "dps", dps, "N", order,
                      "t", mp.nstr(t, 9), "alpha", mp.nstr(root, 28),
                      "delta_aux", mp.nstr(root-seed, 6), "Dcheck", mp.nstr(check, 6),
                      "dD", mp.nstr(derivs[-1], 10),
                      "seconds", round(time.time()-config_start, 1), flush=True)
            (HERE/"low_lambda_refinement_results.json").write_text(json.dumps(report, indent=2)+"\n")
        with mp.workdps(160):
            roots = [mp.mpf(x["alpha"]) for x in anchor["refinements"]]
            anchor["interconfiguration_alpha_difference"] = str(abs(roots[1]-roots[0]))
            anchor["maximum_absolute_auxiliary_difference"] = str(max(abs(x-mp.mpf(original["alpha"])) for x in roots))
            anchor["refined_ratio"] = anchor["refinements"][-1]["ratio_to_2sqrtLambdaOver3"]
            anchor["refined_alpha"] = anchor["refinements"][-1]["alpha"]
    report["complete"] = True
    report["elapsed_seconds"] = time.time()-start
    (HERE/"low_lambda_refinement_results.json").write_text(json.dumps(report, indent=2)+"\n")
    print("COMPLETE", round(time.time()-start, 1), "seconds", flush=True)


if __name__ == "__main__":
    main()
