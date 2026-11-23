#!/usr/bin/env python3
"""Full run (V0.2): all rank-two records of the five pinned ecdata allgens blocks (N < 50,000).

Usage:
  python3 run_study.py ALLGENS_DIR SOURCE_COUNTS_JSON OUTDIR [--workers K] [--limit N]

Protocol (fixed before seeing results):
  * Sources: five allgens blocks, SHA-256 verified against source_counts.json before parsing.
  * Population A: every rank-2 record. Population B: per (conductor, class code), the record
    with the smallest Cremona curve number ("selected representative"; no optimality claim).
  * Rank and full-generator-basis provenance are ACCEPTED from the pinned source.
  * Exact: (X/Z, Y/Z) parsing, torsion parsed separately, both free generators (and torsion
    generators) checked on the curve with exact rationals.
  * Numerical: Gram matrix at requested 192 bits and again at 256 bits, curve object rebuilt per
    precision, explicit `precision=` on every precision-dependent PARI call; Gauss reduction with an
    exact integer U; checks det(U) = +-1, G_red = U~ G U, det(G_red) = det(G), positive definiteness,
    2|b| <= h1 <= h2, and the Gram matrix recomputed from the reduced points themselves.
  * Shear control (added because the supplied Cremona bases turned out to be already reduced, U = I):
    three deliberately non-reduced bases M*(P,Q), det M = 1, are rebuilt from points, reduced, and
    must give the same sign-invariant shape (h1, h2, |b|) to relative 1e-40.
  * Stability threshold 1e-25 (see compare() for the formulae). Boundary-sensitive: margin
    (h1-2|b|)/h1 or (h2-h1)/h2 below 1e-20, or different U between precisions.
    Unstable or boundary-sensitive -> retry 384 bits, then 512 bits. Never skip, clamp or drop.
  V0.2 control flow: a changed valid U between precisions is a retry trigger; a failed or invalid run at
  one precision is retained in the record and the ladder continues; acceptance needs a consecutive
  usable pair; at the top of the ladder a stable pair that is only boundary-sensitive or has different
  valid U is RESOLVED_BOUNDARY_FLAGGED. Exit status 1 if any record is UNRESOLVED (outputs still written).
"""
import sys, os, json, time, hashlib, argparse, traceback
from multiprocessing import Pool

STAB = "1e-25"
BOUNDARY = "1e-20"
INTERNAL = "1e-40"      # internal consistency at a single precision (all >= 256 effective bits)
LADDER = (192, 256, 384, 512)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def valid_single(o):
    """Internal-consistency validation at one precision. Returns list of failure strings."""
    import mwshape as m
    p = m.pari
    bad = []
    if not o.get("raw_positive_definite"):
        return ["raw Gram not positive definite"]
    if abs(o["det_U"]) != 1:
        bad.append("det(U) = %d" % o["det_U"])
    for k in ("rel_err_Gred_vs_UtGU", "rel_err_detGred_vs_detG", "rel_err_Gred_vs_heights_of_reduced_points"):
        if str(o[k]).startswith("not-applicable"):
            continue
        if p(o[k]) >= p(INTERNAL):
            bad.append("%s = %s" % (k, o[k]))
    if not o["red_positive_definite"]:
        bad.append("reduced Gram not positive definite")
    if o.get("effective_bits_below_requested"):
        bad.append("effective precision %d < requested %d" % (o["effective_bits"], o["requested_bits"]))
    if not o["red_inequalities"]:
        bad.append("reduction inequalities violated")
    return bad


def process(line):
    import mwshape as m
    p = m.pari
    t0 = time.time()
    rec = {"source_line": line.rstrip("\n")}
    try:
        row = m.parse_row(line)
    except Exception as e:
        rec.update(status="UNRESOLVED", reason="parse error: %s" % e)
        return rec
    rec.update(label=row["label"], N=row["N"], cls=row["cls"], num=row["num"], ainvs=row["ainvs"],
               torsion=row["torsion"],
               free_generators=[[m.frac_str(q["x"]), m.frac_str(q["y"])] for q in row["free"]],
               free_on_curve=[m.on_curve(row["ainvs"], q) for q in row["free"]],
               torsion_gens_on_curve=[m.on_curve(row["ainvs"], q) for q in row["torsion_gens"]])
    if not all(rec["free_on_curve"]) or not all(rec["torsion_gens_on_curve"]):
        rec.update(status="UNRESOLVED", reason="exact point-membership failure")
        return rec
    runs, nums, pairs, events, failures = {}, {}, [], [], {}

    def attempt(bits):
        """Compute one precision; numerical failures/exceptions are retained, never fatal."""
        try:
            o, n = m.compute_at(row, bits)
        except Exception as e:
            o, n = {"requested_bits": bits, "error": "exception: %s" % e,
                    "traceback": traceback.format_exc()}, None
        runs[bits] = o
        bad = [o["error"]] if n is None else valid_single(o)
        if bad:
            failures[bits] = bad
            n = None            # unusable for pairing
        nums[bits] = n

    try:
        attempt(LADDER[0]); attempt(LADDER[1])
        idx = 1                 # index of the highest precision computed so far
        final_pair = None       # (lo, hi) of the accepted consecutive usable pair
        top_flag = None
        while True:
            lo, hi = LADDER[idx - 1], LADDER[idx]
            reason = []
            if nums[lo] is None or nums[hi] is None:
                reason.append("unusable run at %s bits" % ",".join(str(b) for b in (lo, hi) if nums[b] is None))
                pair = {"low_bits": lo, "high_bits": hi, "usable": False,
                        "validation_failures": {str(b): failures[b] for b in (lo, hi) if b in failures}}
            else:
                same_U, d = m.compare(nums[lo], nums[hi])
                keys = ["gram_rel_frobenius", "det_rel", "x_abs", "s_abs", "y_rel", "inv_scaled_max"] + \
                       (["gred_rel_frobenius"] if same_U else [])
                stable = all(d[k] < p(STAB) for k in keys)
                boundary = (p(runs[hi]["margin_2b"]) < p(BOUNDARY)) or (p(runs[hi]["margin_h"]) < p(BOUNDARY))
                pair = {"low_bits": lo, "high_bits": hi, "usable": True, "same_U": same_U,
                        "diffs": {k: m.sig(v, 6) for k, v in d.items()},
                        "stable": stable, "boundary_sensitive": boundary}
                if not stable: reason.append("unstable")
                if boundary: reason.append("boundary-sensitive")
                if not same_U: reason.append("U differs (both valid)")
                if not reason:
                    final_pair = (lo, hi)
            pairs.append(pair)
            if final_pair:
                break
            if idx + 1 >= len(LADDER):
                # top of the ladder: a usable, stable pair that is only boundary-sensitive and/or has
                # different valid U is accepted with a flag (sign-invariant shape agrees); otherwise unresolved
                if pair.get("usable") and pair.get("stable"):
                    final_pair = (lo, hi)
                    top_flag = [r for r in reason]
                break
            idx += 1
            events.append("retry at %d bits: %s (pair %d/%d)" % (LADDER[idx], ", ".join(reason), lo, hi))
            attempt(LADDER[idx])
        if final_pair:
            hb = final_pair[1]
            rec["shear_controls"] = m.shear_control(row, hb, nums[hb])
            shear_ok = all(c["pass"] is True or str(c["pass"]).startswith("not-applicable") for c in rec["shear_controls"])
            rec["final_bits"] = hb
            rec["final"] = {k: runs[hb][k] for k in
                            ("effective_bits", "U", "det_U", "G", "G_red", "reg", "x", "y", "s", "margin_2b", "margin_h")}
            rec["boundary_sensitive_final"] = bool(pairs[-1].get("boundary_sensitive"))
            rec["U_changed_between_final_pair"] = not pairs[-1].get("same_U", True)
            rec["x_near_zero"] = bool(p(rec["final"]["x"]) < p("1e-30"))
            if not shear_ok:
                events.append("shear control failed")
                rec["status"] = "UNRESOLVED"; rec["reason"] = "shear control failed at %d bits" % hb
            elif top_flag:
                rec["status"] = "RESOLVED_BOUNDARY_FLAGGED"; rec["flag_reason"] = ", ".join(top_flag)
            else:
                rec["status"] = "RESOLVED"
        else:
            rec["status"] = "UNRESOLVED"
            rec["reason"] = "no consecutive usable stable pair up to %d bits" % LADDER[idx]
    except Exception as e:
        rec.update(status="UNRESOLVED", reason="exception: %s" % e, traceback=traceback.format_exc())
    rec["runs"] = {str(b): r for b, r in runs.items()}
    rec["pairs"] = pairs
    rec["retry_events"] = events
    rec["precision_failures"] = {str(b): f for b, f in failures.items()}
    rec["elapsed_s"] = round(time.time() - t0, 4)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("allgens_dir"); ap.add_argument("source_counts"); ap.add_argument("outdir")
    ap.add_argument("--workers", type=int, default=2); ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    t_start = time.time()
    sc = json.load(open(a.source_counts))
    lines, src = [], []
    for b in sc["blocks"]:
        path = os.path.join(a.allgens_dir, "allgens." + b["block"])
        h = sha256(path)
        src.append({"block": b["block"], "path": path, "sha256": h, "expected_sha256": b["sha256"], "match": h == b["sha256"]})
        if h != b["sha256"]:
            sys.exit("SHA-256 mismatch for %s" % path)
        cnt = 0
        for line in open(path):
            t = line.split()
            if len(t) > 4 and t[4] == "2":
                lines.append(line); cnt += 1
        src[-1]["rank2_records"] = cnt
        src[-1]["expected_rank2_records"] = b["rank2curves"]
    # Population B selection (frozen before any computation)
    best = {}
    for line in lines:
        t = line.split()
        key = (int(t[0]), t[1]); num = int(t[2])
        if key not in best or num < best[key][0]:
            best[key] = (num, t[0] + t[1] + t[2])
    popB = sorted((k[0], k[1], v[0], v[1]) for k, v in best.items())
    with open(os.path.join(a.outdir, "population_B_selection.csv"), "w") as f:
        f.write("conductor,class_code,selected_curve_number,label\n")
        for r in popB:
            f.write("%d,%s,%d,%s\n" % r)
    json.dump({"sources": src, "population_A_records": len(lines), "population_B_classes": len(popB),
               "commit": sc["commit"]}, open(os.path.join(a.outdir, "frozen_populations.json"), "w"), indent=1)
    print("Population A:", len(lines), " Population B:", len(popB), flush=True)
    if a.limit:
        step = max(1, len(lines) // a.limit)
        lines = lines[::step][:a.limit]
    out = open(os.path.join(a.outdir, "per_curve_results.jsonl"), "w")
    n = n_unres = 0
    with Pool(a.workers) as pool:
        for rec in pool.imap(process, lines, chunksize=8):
            n_unres += rec["status"] == "UNRESOLVED"
            out.write(json.dumps(rec, separators=(",", ":")) + "\n")
            n += 1
            if n % 1000 == 0:
                print("%d/%d  %.0fs" % (n, len(lines), time.time() - t_start), flush=True)
    out.close()
    json.dump({"records_written": n, "unresolved": n_unres, "elapsed_s": round(time.time() - t_start, 1), "workers": a.workers,
               "start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t_start)),
               "end_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())},
              open(os.path.join(a.outdir, "run_timing.json"), "w"), indent=1)
    print("done", n, "records in %.1fs" % (time.time() - t_start), "unresolved:", n_unres)
    sys.exit(1 if n_unres else 0)


if __name__ == "__main__":
    main()
