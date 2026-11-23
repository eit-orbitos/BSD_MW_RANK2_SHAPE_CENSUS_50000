#!/usr/bin/env python3
"""Analysis (V0.2) of the rank-two lattice-shape study: populations A and B, overall and per conductor block.

Usage: python3 analyze.py OUTDIR   (reads OUTDIR/per_curve_results.jsonl and population_B_selection.csv)

Numerical summaries are computed in PARI at 256-bit working precision from the 45-significant-digit
serialized PARI values (parsing at 256 bits; no float conversion). Floats are used only for plots.
Reference model: normalized hyperbolic (Haar) measure (3/pi) dx dy / y^2 on the SL2(Z) fundamental
domain; with x = |b|/h1 >= 0 the domain is folded to 0 <= x <= 1/2, density (6/pi) dx dy / y^2.
  E_mu[s]   = (3/pi) Cl_2(pi/3) = (6/pi) int_0^{1/2} arcsin(x)/x dx = 0.969197841658351...
  E_mu[s^2] = E_mu[s]/2 + 1/4 + 3 sqrt(3)/(8 pi)
  P_mu(y >= Y) = 3/(pi Y) for Y >= 1
  F_s(t) = (6/pi)[pi/6 - arcsin sqrt(1-t) - sqrt((1-t)/t) log(1/(2 sqrt(1-t)))], 3/4 <= t < 1
V0.2: the retry/failure ledger is written FIRST; unresolved records are retained in every table;
expected counts come from the frozen populations (frozen_populations.json, population_B_selection.csv
and the pinned source lines), not from the records that happened to be written; moments use resolved
records only; status INCOMPLETE and exit code 1 when anything is unresolved or missing.
x-bin selection uses PARI arithmetic: bin = min(9, floor(20*x)).
No convergence to this model is assumed. Distances below are descriptive, not test statistics:
isogenous curves (population A) and arithmetic populations are not IID samples.
"""
import sys, os, json, csv, math
import mwshape as m

p = m.pari
p.set_real_precision_bits(256)   # affects GP-string parsing and GP-level arithmetic defaults
BITS = 256
OUT = sys.argv[1]
BLOCKS = [(0, 9999), (10000, 19999), (20000, 29999), (30000, 39999), (40000, 49999)]
YS = [1, 2, 4, 8, 16]
WATKINS = {1: "93.5", 2: "35.6", 4: "9.7", 8: "1.97", 16: "0.35"}   # percent, Watkins 2008 §6.4 (different population)

PI = p.Pi(precision=BITS)
CL2 = p("imag(polylog(2, exp(I*Pi/3)))")            # Cl_2(pi/3) at 256 bits
REF_MEAN = 3 / PI * CL2
REF_MEAN_INT = 6 / PI * p("intnum(x=0,1/2,if(x,asin(x)/x,1))")
REF_M2 = REF_MEAN / 2 + p("1/4") + 3 * p.sqrt(p(3), precision=BITS) / (8 * PI)
REF_MED = p("solve(t=3/4+1e-40,1-1e-40,my(a=sqrt(1-t));(6/Pi)*(Pi/6-asin(a)-sqrt((1-t)/t)*log(1/(2*a)))-1/2)")
REF_F = p("(t)->if(t<=3/4,0,if(t>=1,1,my(a=sqrt(1-t));(6/Pi)*(Pi/6-asin(a)-sqrt((1-t)/t)*log(1/(2*a)))))")


def load():
    R = [json.loads(l) for l in open(os.path.join(OUT, "per_curve_results.jsonl"))]
    B = {r["label"] for r in csv.DictReader(open(os.path.join(OUT, "population_B_selection.csv")))}
    fz = json.load(open(os.path.join(OUT, "frozen_populations.json")))
    # frozen labels of population A, from the pinned source lines (independent of what was written)
    A = []
    for src in fz["sources"]:
        for line in open(src["path"]):
            t = line.split()
            if len(t) > 4 and t[4] == "2":
                A.append((t[0] + t[1] + t[2], int(t[0])))
    return R, B, fz, A


def is_resolved(r):
    return str(r.get("status", "")).startswith("RESOLVED")


def summarize(rows, expected):
    res = {"expected": expected, "completed": 0, "unresolved": 0}
    ok = [r for r in rows if is_resolved(r)]
    res["completed"] = len(ok)
    res["unresolved_in_records"] = len(rows) - len(ok)
    res["missing_records"] = expected - len(rows)
    res["unresolved"] = res["unresolved_in_records"] + res["missing_records"]
    res["boundary_flagged"] = sum(r["status"] == "RESOLVED_BOUNDARY_FLAGGED" for r in ok)
    if not ok:
        return res
    S = sorted(((p(r["final"]["s"]), r["label"]) for r in ok), key=lambda t: t[0])
    s = [t[0] for t in S]
    y = [p(r["final"]["y"]) for r in ok]
    n = len(s)
    mean = sum(s) / n
    m2 = sum(v * v for v in s) / n
    med = s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2
    # descriptive sup-distance between the empirical CDF and the reference CDF (evaluated at both sides of each jump)
    D = p(0)
    for i, v in enumerate(s):
        F = REF_F(v)
        D = max(D, abs(p(i + 1) / n - F), abs(p(i) / n - F))
    res.update({
        "n": n,
        "mean_s": m.sig(mean), "median_s": m.sig(med), "second_moment_s": m.sig(m2),
        "variance_s": m.sig(m2 - mean * mean),
        "min_s": m.sig(s[0]), "min_s_label": S[0][1], "max_s": m.sig(s[-1]), "max_s_label": S[-1][1],
        "mean_minus_reference": m.sig(mean - REF_MEAN, 12),
        "second_moment_minus_reference": m.sig(m2 - REF_M2, 12),
        "ecdf_sup_distance_to_reference_descriptive": m.sig(D, 12),
        "prop_y_ge": {str(Y): m.sig(p(sum(v >= Y for v in y)) / n, 12) for Y in YS},
        "count_y_ge": {str(Y): int(sum(v >= Y for v in y)) for Y in YS},
        "reference_prop_y_ge": {str(Y): m.sig(3 / (PI * Y), 12) for Y in YS},
        "max_y": m.sig(max(y), 12),
    })
    # x-distribution for y >= 1 (10 bins on [0, 1/2]; reference: uniform 10% each), cf. Watkins Table 1
    xs = [p(r["final"]["x"]) for r, yy in zip(ok, y) if yy >= 1]
    bins = [0] * 10
    for v in xs:
        bins[min(9, int(p.floor(v * 20)))] += 1          # PARI arithmetic (V0.1 used float(v*20))
    res["x_bins_given_y_ge_1_percent"] = [round(100.0 * c / len(xs), 3) for c in bins] if xs else []
    return res


def write_ledger(R, Aset):
    """Written first, so it exists even if a later step fails."""
    seen = {r.get("label") for r in R}
    rows = []
    for r in R:
        if r.get("retry_events") or not is_resolved(r) or r.get("precision_failures"):
            rows.append([r.get("label", "(unparsed)"), r.get("status"), r.get("final_bits"),
                         " | ".join(r.get("retry_events", [])),
                         json.dumps(r.get("precision_failures", {})), r.get("reason", r.get("flag_reason", ""))])
    for lab in sorted(Aset - seen):
        rows.append([lab, "MISSING", "", "", "", "record absent from per_curve_results.jsonl"])
    with open(os.path.join(OUT, "retry_failure_ledger.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "status", "final_bits", "retry_events", "precision_failures", "reason"])
        if rows:
            w.writerows(rows)
        else:
            w.writerow(["(none)", "all records RESOLVED at the first 192/256-bit pair; no retries; no failures", "", "", "", ""])
    return rows


def main():
    R, Bset, fz, Alist = load()
    Aset = {a for a, _ in Alist}
    ledger = write_ledger(R, Aset)
    pops = {"A_all_rank2_records": R, "B_selected_class_representatives": [r for r in R if r.get("label") in Bset]}
    expected = {"A_all_rank2_records": len(Aset), "B_selected_class_representatives": len(Bset)}
    frozenN = {"A_all_rank2_records": dict(Alist), "B_selected_class_representatives": {a: n for a, n in Alist if a in Bset}}
    out = {"reference": {
        "mean_s_closed_form_(3/pi)Cl2(pi/3)": m.sig(REF_MEAN, 40),
        "mean_s_integral_(6/pi)int_0^1/2_asin(x)/x": m.sig(REF_MEAN_INT, 40),
        "closed_form_minus_integral": m.sig(REF_MEAN - REF_MEAN_INT, 4),
        "second_moment_s": m.sig(REF_M2, 40),
        "median_s_solve_F(t)=1/2": m.sig(REF_MED, 40),
        "P(y>=Y)=3/(pi*Y)": {str(Y): m.sig(3 / (PI * Y), 20) for Y in YS},
        "Watkins_2008_section_6.4_percent_HISTORICAL_DIFFERENT_POPULATION": WATKINS,
        "note": "Watkins: 2,143,079 rank-2 curves of prime conductor < 1e10 (Stein-Watkins database); not an expected value for this sample."}}
    for name, rows in pops.items():
        o = {"overall": summarize(rows, expected[name]), "blocks": {}}
        for lo, hi in BLOCKS:
            sub = [r for r in rows if lo <= r.get("N", -1) <= hi]
            exp_blk = sum(lo <= n <= hi for n in frozenN[name].values())   # from the frozen population
            o["blocks"]["%05d-%05d" % (lo, hi)] = summarize(sub, exp_blk)
        out[name] = o
    json.dump(out, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    # flat summary table
    with open(os.path.join(OUT, "summary_table.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["population", "conductor_block", "expected", "completed", "unresolved", "mean_s", "median_s",
                    "second_moment_s", "min_s", "min_label", "max_s", "max_label", "mean_minus_ref",
                    "ecdf_sup_distance_descriptive"] + ["P(y>=%d)" % Y for Y in YS])
        for name in pops:
            for blk, sres in [("all", out[name]["overall"])] + list(out[name]["blocks"].items()):
                g = lambda k: sres.get(k, "")
                w.writerow([name, blk, sres["expected"], sres["completed"], sres["unresolved"], g("mean_s"),
                            g("median_s"), g("second_moment_s"), g("min_s"), g("min_s_label"),
                            g("max_s"), g("max_s_label"), g("mean_minus_reference"),
                            g("ecdf_sup_distance_to_reference_descriptive")] + [sres.get("prop_y_ge", {}).get(str(Y), "") for Y in YS])
        w.writerow(["reference_Haar", "-", "", "", "", out["reference"]["mean_s_closed_form_(3/pi)Cl2(pi/3)"], "",
                    out["reference"]["second_moment_s"], "0.75", "", "1", "", "0", "0"] +
                   [out["reference"]["P(y>=Y)=3/(pi*Y)"][str(Y)] for Y in YS])
    # compact per-curve table (final values)
    with open(os.path.join(OUT, "per_curve_shapes.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["label", "N", "class", "curve_number", "in_population_B", "status", "final_requested_bits",
                    "final_effective_bits", "h1", "b", "h2", "regulator", "x", "y", "s", "margin_2b", "margin_h"])
        for r in R:
            fr = r.get("final")
            base = [r.get("label"), r.get("N"), r.get("cls"), r.get("num"), r.get("label") in Bset, r.get("status"), r.get("final_bits")]
            if fr is None:   # unresolved: retained with empty numerical fields
                w.writerow(base + [""] * 10)
            else:
                w.writerow(base + [fr["effective_bits"], fr["G_red"][0], fr["G_red"][1], fr["G_red"][2], fr["reg"],
                                   fr["x"], fr["y"], fr["s"], fr["margin_2b"], fr["margin_h"]])
    # precision ledger: every curve, every pair (difference ledger)
    with open(os.path.join(OUT, "precision_difference_ledger.csv"), "w", newline="") as f:
        w = csv.writer(f)
        keys = ["gram_rel_frobenius", "det_rel", "x_abs", "s_abs", "y_rel", "inv_scaled_max", "gred_rel_frobenius"]
        w.writerow(["label", "low_bits", "high_bits", "low_effective_bits", "high_effective_bits", "usable", "same_U",
                    "stable", "boundary_sensitive"] + keys)
        for r in R:
            for pr in r.get("pairs", []):
                eb = lambda b: r["runs"].get(str(b), {}).get("effective_bits", "")
                w.writerow([r.get("label"), pr["low_bits"], pr["high_bits"], eb(pr["low_bits"]), eb(pr["high_bits"]),
                            pr.get("usable", True), pr.get("same_U", ""), pr.get("stable", ""), pr.get("boundary_sensitive", "")]
                           + [pr.get("diffs", {}).get(k, "") for k in keys])
    unres = sum(out[k]["overall"]["unresolved"] for k in pops)
    out["status"] = "COMPLETE" if unres == 0 else "INCOMPLETE"
    out["ledger_rows"] = len(ledger)
    json.dump(out, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(json.dumps({k: {kk: out[k]["overall"][kk] for kk in ("expected", "completed", "unresolved", "mean_s") if kk in out[k]["overall"]} for k in pops}, indent=1))
    print("STATUS:", out["status"], "| ledger rows:", len(ledger))
    sys.exit(0 if out["status"] == "COMPLETE" else 1)


if __name__ == "__main__":
    main()
