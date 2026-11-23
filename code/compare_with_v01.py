#!/usr/bin/env python3
"""Core-value compatibility comparison of V0.2 against the unchanged V0.1 baseline.

Usage: python3 compare_with_v01.py V01_ROOT V02_ROOT OUT_JSON
Compares: source identities, representative selection, per-curve final 45-digit G, G_red, regulator,
x, y, s (string equality) plus U, status, final precision; every population/block summary field.
Diagnostics whose definition changed in V0.2 (sign-invariant shape error) are reported side by side
under their own definitions, not as a compatibility target. Exit 1 on any core mismatch.
"""
import sys, json, hashlib, os

V1, V2, OUTJ = sys.argv[1:4]
res, bad = {}, []


def h(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


f1, f2 = json.load(open(os.path.join(V1, "results/frozen_populations.json"))), json.load(open(os.path.join(V2, "results/frozen_populations.json")))
res["sources_identical"] = [(a["block"], a["sha256"], a["sha256"] == b["sha256"]) for a, b in zip(f1["sources"], f2["sources"])]
res["population_counts"] = {"V0.1": [f1["population_A_records"], f1["population_B_classes"]],
                            "V0.2": [f2["population_A_records"], f2["population_B_classes"]]}
sel1, sel2 = h(os.path.join(V1, "results/population_B_selection.csv")), h(os.path.join(V2, "results/population_B_selection.csv"))
res["population_B_selection_sha256"] = {"V0.1": sel1, "V0.2": sel2, "identical": sel1 == sel2}
if not all(x[2] for x in res["sources_identical"]) or sel1 != sel2 or res["population_counts"]["V0.1"] != res["population_counts"]["V0.2"]:
    bad.append("source/selection identity")

R1 = {r["label"]: r for r in map(json.loads, open(os.path.join(V1, "results/per_curve_results.jsonl")))}
R2 = {r["label"]: r for r in map(json.loads, open(os.path.join(V2, "results/per_curve_results.jsonl")))}
CORE = ["G", "G_red", "reg", "x", "y", "s"]
OTHER = ["U", "det_U", "effective_bits", "margin_2b", "margin_h"]
mism = {k: [] for k in CORE + OTHER + ["status", "final_bits", "label_set"]}
if set(R1) != set(R2):
    mism["label_set"] = sorted(set(R1) ^ set(R2))
for lab in R1.keys() & R2.keys():
    a, b = R1[lab], R2[lab]
    for k in CORE + OTHER:
        if a["final"][k] != b["final"][k]:
            mism[k].append(lab)
    if a["status"] != b["status"]:
        mism["status"].append(lab)
    if a["final_bits"] != b["final_bits"]:
        mism["final_bits"].append(lab)
res["per_curve"] = {"records_V0.1": len(R1), "records_V0.2": len(R2),
                    "mismatch_counts": {k: len(v) for k, v in mism.items()},
                    "mismatch_examples": {k: v[:10] for k, v in mism.items() if v}}
if any(mism[k] for k in CORE + ["status", "final_bits", "label_set", "U"]):
    bad.append("per-curve core values")

S1, S2 = json.load(open(os.path.join(V1, "results/summary.json"))), json.load(open(os.path.join(V2, "results/summary.json")))
diffs = []
def walk(a, b, path):
    if isinstance(a, dict):
        for k in a:
            if k not in b:
                diffs.append((path + "/" + k, "missing in V0.2")); continue
            walk(a[k], b[k], path + "/" + k)
    elif a != b:
        diffs.append((path, "V0.1=%s V0.2=%s" % (a, b)))
for pop in ("A_all_rank2_records", "B_selected_class_representatives"):
    walk(S1[pop], S2[pop], pop)
walk(S1["reference"], S2["reference"], "reference")
res["summary_field_differences"] = diffs
res["summary_new_fields_in_V0.2"] = sorted({k for pop in ("A_all_rank2_records", "B_selected_class_representatives")
                                            for k in S2[pop]["overall"] if k not in S1[pop]["overall"]}
                                           | {"reference/" + k for k in S2["reference"] if k not in S1["reference"]}
                                           | {k for k in S2 if k not in S1})
if diffs:
    bad.append("summary values")

# diagnostics under their own definitions (not a compatibility target)
def mx(R, key):
    v = []
    for r in R.values():
        for p in r["pairs"]:
            x = p.get("diffs", {}).get(key)
            if x is not None:
                v.append(0.0 if "E" in x and x.startswith("0.") else float(x))
    return max(v) if v else None
def mx_shear(R, key):
    return max(float(c[key]) for r in R.values() for c in r["shear_controls"] if key in c)
res["diagnostics_by_definition"] = {
    "V0.1 inv_rel_max (componentwise relative, divides by |b|) max over all pairs": mx(R1, "inv_rel_max"),
    "V0.2 inv_scaled_max (diagonal-scale normalized) max over all pairs": mx(R2, "inv_scaled_max"),
    "V0.1 shear max_rel_diff_invariants max": mx_shear(R1, "max_rel_diff_invariants"),
    "V0.2 shear max_scaled_diff_invariants max": mx_shear(R2, "max_scaled_diff_invariants"),
    "unchanged definitions, max over all pairs V0.1 / V0.2": {k: [mx(R1, k), mx(R2, k)] for k in
        ("gram_rel_frobenius", "det_rel", "x_abs", "s_abs", "y_rel", "gred_rel_frobenius")},
}
res["status"] = "COMPATIBLE" if not bad else "MISMATCH: " + ", ".join(bad)
json.dump(res, open(OUTJ, "w"), indent=1)
print(json.dumps({k: res[k] for k in ("population_counts", "population_B_selection_sha256", "status")}, indent=1))
print("per-curve mismatches:", res["per_curve"]["mismatch_counts"])
print("summary field differences:", len(diffs), diffs[:5])
print(json.dumps(res["diagnostics_by_definition"], indent=1))
sys.exit(0 if not bad else 1)
