"""Create a compact review summary from the executed two-anchor refinements."""
from pathlib import Path
import hashlib
import json
import mpmath as mp

HERE = Path(__file__).resolve().parent
mp.mp.dps = 180
report = json.loads((HERE / "low_lambda_refinement_results.json").read_text())
assert report["complete"]
rows = [r for a in report["anchors"] for r in a["refinements"]]
source = HERE / "low_lambda_refinement.py"
assert hashlib.sha256(source.read_bytes()).hexdigest() == report["provenance"]["source_sha256"]
summary = {
    "complete": True,
    "source_hash_matches_executed_run": True,
    "maximum_refinement_configuration_difference_in_alpha": str(max(
        mp.mpf(a["interconfiguration_alpha_difference"]) for a in report["anchors"])),
    "maximum_difference_from_previous_auxiliary_alpha": str(max(
        mp.mpf(a["maximum_absolute_auxiliary_difference"]) for a in report["anchors"])),
    "maximum_absolute_higher_order_normalized_determinant": str(max(
        abs(mp.mpf(r["normalized_W_check"])) for r in rows)),
    "minimum_integer_index_distance": str(min(
        mp.mpf(r["minimum_integer_index_distance"]) for r in rows)),
    "maximum_derivative_relative_change": str(max(
        mp.mpf(r["derivative_relative_change"]) for r in rows)),
    "all_symmetric_perturbations_change_sign": all(
        p["opposite_signs"] for r in rows for p in r["symmetric_perturbation_checks"]),
    "all_matching_values_and_solution_vectors_nonzero": all(
        r["nonzero_solution_checks"]["both_matching_values_nonzero"]
        and r["nonzero_solution_checks"]["both_solution_vectors_nonzero"] for r in rows),
    "normalization": report["frequency_conventions"],
    "anchors": [],
    "limitations": report["limitations"],
    "provenance": report["provenance"],
    "elapsed_seconds": report["elapsed_seconds"]}
for anchor in report["anchors"]:
    ratio = mp.mpf(anchor["refined_ratio"])
    last = anchor["refinements"][-1]
    rh = mp.mpf(last["geometry_and_series"]["rh"])
    rc = mp.mpf(last["geometry_and_series"]["rc"])
    summary["anchors"].append({
        "Lambda_tilde": anchor["Lambda_tilde"],
        "alpha": anchor["refined_alpha"],
        "varpi_im": last["varpi_im"],
        "ratio_to_2sqrtLambdaOver3": str(ratio),
        "ratio_minus_one": str(ratio-1),
        "ratio_change_from_auxiliary": str(ratio-mp.mpf(anchor["original_auxiliary"]["ratio_to_2sqrtLambdaOver3"])),
        "r_h_over_r_c": str(rh/rc),
        "refinement_alpha_difference": anchor["interconfiguration_alpha_difference"],
        "maximum_auxiliary_alpha_difference": anchor["maximum_absolute_auxiliary_difference"],
        "configurations": [{
            "dps": r["dps"], "N": r["N"], "check_N": r["check_N"],
            "matching_t": r["matching_t"],
            "d_normalized_D_dalpha": r["derivative_estimates"][-1],
            "derivative_relative_change": r["derivative_relative_change"],
            "normalized_D_check": r["normalized_W_check"],
            "integer_index_distance": r["minimum_integer_index_distance"],
            "matching_value_product_without_nonzero_powers": r["nonzero_solution_checks"]["matching_value_product_without_nonzero_powers"]
        } for r in anchor["refinements"]]})
(HERE / "low_lambda_refinement_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
for key, value in summary.items():
    if key.startswith("maximum") or key.startswith("minimum"):
        print(key, mp.nstr(mp.mpf(value), 16))
for anchor in summary["anchors"]:
    print("Lambda", anchor["Lambda_tilde"], "alpha", mp.nstr(mp.mpf(anchor["alpha"]), 32),
          "dSratio", mp.nstr(mp.mpf(anchor["ratio_to_2sqrtLambdaOver3"]), 25),
          "rh/rc", mp.nstr(mp.mpf(anchor["r_h_over_r_c"]), 15))
    for row in anchor["configurations"]:
        print("  settings", row["dps"], row["N"], row["check_N"],
              "t", mp.nstr(mp.mpf(row["matching_t"]), 18),
              "derivative", mp.nstr(mp.mpf(row["d_normalized_D_dalpha"]), 18),
              "nonzero_value_product", mp.nstr(mp.mpf(row["matching_value_product_without_nonzero_powers"]), 12))
