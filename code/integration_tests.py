#!/usr/bin/env python3
"""V0.2 integration tests of the repaired control paths, executed with PARI/cypari2.

All cases marked INJECTED are fault-injection controls on the real 389a1 row (or a small set of
real rows). They are NOT naturally occurring curve results.

Usage: python3 integration_tests.py ALLGENS_DIR WORKDIR OUT_JSON
"""
import sys, os, json, subprocess, shutil, hashlib
import mwshape as m
import run_study as rs

p = m.pari
p.set_real_precision_bits(256)
HERE = os.path.dirname(os.path.abspath(__file__))
ALLGENS, WORK, OUTJ = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(WORK, exist_ok=True)
L389 = "389 a 1 [0,1,1,-2,0] 2 [] [0:0:1] [1:0:1]\n"
T, fails = [], []


def record(name, ok, detail, injected):
    T.append({"name": name, "pass": bool(ok), "injected_control": injected, "detail": detail})
    if not ok:
        fails.append(name)
    print("%-4s %s%s" % ("OK" if ok else "FAIL", "[INJECTED] " if injected else "", name))


def old_inv_rel_max(a, b):
    """V0.1 definition, reproduced for comparison only."""
    return max(abs(a[k] - b[k]) / abs(b[k]) for k in range(3))


def run_case(spec_by_bits):
    m.INJECT.clear()
    for bits, spec in spec_by_bits.items():
        m.INJECT[("389a1", bits)] = spec
    rec = rs.process(L389)
    m.INJECT.clear()
    return rec


# --- 1. square / rectangular Gram forms (exact b = 0), arithmetic coverage
for name, G0 in (("square", "[1,0;0,1]"), ("rectangular", "[1,0;0,3]")):
    ga = p(G0) * p.sqrt(p(1), precision=192)
    gb = p(G0) * p.sqrt(p(1), precision=256)
    ra, Ua, _ = m.gauss_reduce(ga); rb, Ub, _ = m.gauss_reduce(gb)
    ia, ib = m.shape_from_reduced(ra, 192)[:3], m.shape_from_reduced(rb, 256)[:3]
    new = m.scaled_shape_error(ia, ib)
    try:
        old = old_inv_rel_max(ia, ib); old_s = str(old)
    except Exception as e:
        old_s = "V0.1 definition fails: %s" % type(e).__name__
    record("%s Gram (b = 0 exactly): scaled shape error defined and < 1e-25" % name,
           new < p("1e-25"), {"new_scaled_error": m.sig(new, 4), "old_V0.1_inv_rel_max": old_s}, False)

# exact b = 0 through the whole pipeline on a real row (INJECTED at every precision)
rec = run_case({192: "set_b:0", 256: "set_b:0"})
record("rectangular via pipeline: b set to exact 0 at 192 and 256 bits -> RESOLVED, x = 0, no division error",
       rec["status"] == "RESOLVED" and rec["final"]["x"] == "0",
       {"status": rec["status"], "x": rec["final"]["x"] if "final" in rec else None,
        "inv_scaled_max": rec["pairs"][0].get("diffs", {}).get("inv_scaled_max")}, True)

# --- 2. near-zero pairing under scaling: b = 1e-50*h1 at 192 bits, 2e-50*h1 at 256 bits (INJECTED)
rec = run_case({192: "b_over_h1:1e-50", 256: "b_over_h1:2e-50"})
ia = (p(rec["runs"]["192"]["G_red"][0]), p(rec["runs"]["192"]["G_red"][2]), abs(p(rec["runs"]["192"]["G_red"][1])))
ib = (p(rec["runs"]["256"]["G_red"][0]), p(rec["runs"]["256"]["G_red"][2]), abs(p(rec["runs"]["256"]["G_red"][1])))
record("near-zero pairing: V0.1 relative error would be 0.5 (spurious instability); V0.2 scaled error ~1e-50 -> RESOLVED",
       rec["status"] == "RESOLVED" and p(rec["pairs"][0]["diffs"]["inv_scaled_max"]) < p("1e-25"),
       {"status": rec["status"], "old_V0.1_inv_rel_max": m.sig(old_inv_rel_max(ia, ib), 4),
        "new_inv_scaled_max": rec["pairs"][0]["diffs"]["inv_scaled_max"]}, True)

# --- 3. valid boundary change of U: b straddles h1/2 alternately at every rung (INJECTED)
# 3a. guard control: injected constants parsed at only 256 bits must be caught by the effective-precision
#     guard at the 384/512-bit rungs (the run is then unusable and the record UNRESOLVED)
STRADDLE = {192: "b_over_h1:0.5000000000000000000000000000000000000000001",
            256: "b_over_h1:0.4999999999999999999999999999999999999999999",
            384: "b_over_h1:0.5000000000000000000000000000000000000000001",
            512: "b_over_h1:0.4999999999999999999999999999999999999999999"}
p.set_real_precision_bits(256)
rec = run_case(STRADDLE)
record("effective-precision guard: 256-bit injected constant at 384/512 bits -> runs rejected, record UNRESOLVED",
       rec["status"] == "UNRESOLVED" and all("effective precision" in f[0] for b, f in rec["precision_failures"].items())
       and set(rec["precision_failures"]) == {"384", "512"},
       {"status": rec["status"], "precision_failures": rec["precision_failures"]}, True)
# 3b. the boundary case itself, with injected constants parsed at 576 bits
p.set_real_precision_bits(576)
rec = run_case(STRADDLE)
p.set_real_precision_bits(256)
Us = {b: rec["runs"][str(b)]["U"] for b in (192, 256, 384, 512)}
record("valid boundary change of U: retried to 512 bits, U differs between rungs, accepted as RESOLVED_BOUNDARY_FLAGGED",
       rec["status"] == "RESOLVED_BOUNDARY_FLAGGED" and rec["final_bits"] == 512 and rec["U_changed_between_final_pair"]
       and all(rec["runs"][str(b)]["red_inequalities"] for b in Us),
       {"status": rec["status"], "final_bits": rec.get("final_bits"), "U_by_bits": {str(k): v for k, v in Us.items()},
        "retry_events": rec["retry_events"], "flag_reason": rec.get("flag_reason"),
        "final_pair_inv_scaled_max": rec["pairs"][-1]["diffs"]["inv_scaled_max"]}, True)

# --- 4. retry case: instability at 192 bits only (INJECTED perturbation 1e-20 * h2 on b)
rec = run_case({192: "perturb:1e-20"})
record("retry path: unstable 192/256 pair -> retry at 384 -> consecutive stable pair 256/384 -> RESOLVED at 384",
       rec["status"] == "RESOLVED" and rec["final_bits"] == 384 and len(rec["retry_events"]) == 1,
       {"status": rec["status"], "final_bits": rec.get("final_bits"), "retry_events": rec["retry_events"]}, True)

# --- 5a. exception at 192 only: retained, ladder continues (INJECTED)
rec = run_case({192: "exception"})
record("numerical failure retained: exception at 192 bits recorded, ladder continues, RESOLVED from pair 256/384",
       rec["status"] == "RESOLVED" and rec["final_bits"] == 384 and "192" in rec["precision_failures"],
       {"status": rec["status"], "final_bits": rec.get("final_bits"), "precision_failures": rec["precision_failures"],
        "retry_events": rec["retry_events"]}, True)
# --- 5b. exception at every precision -> UNRESOLVED with retained failures (INJECTED)
rec = run_case({b: "exception" for b in (192, 256, 384, 512)})
record("exception at all precisions -> UNRESOLVED, all four failures retained",
       rec["status"] == "UNRESOLVED" and len(rec["precision_failures"]) == 4,
       {"status": rec["status"], "reason": rec.get("reason"), "failed_bits": sorted(rec["precision_failures"])}, True)

# --- 6. unresolved export + process exit codes through the real scripts on a small real source block
mini = os.path.join(WORK, "mini_allgens"); out = os.path.join(WORK, "mini_out")
shutil.rmtree(mini, ignore_errors=True); shutil.rmtree(out, ignore_errors=True); os.makedirs(mini)
src_lines = []
for line in open(os.path.join(ALLGENS, "allgens.00000-09999")):
    t = line.split()
    if t[4] == "2" and len(src_lines) < 12:
        src_lines.append(line)
with open(os.path.join(mini, "allgens.mini"), "w") as f:
    f.writelines(src_lines)
h = hashlib.sha256(open(os.path.join(mini, "allgens.mini"), "rb").read()).hexdigest()
classes = {(l.split()[0], l.split()[1]) for l in src_lines}
json.dump({"commit": "mini-test (subset of pinned 00000-09999)", "blocks": [{"block": "mini", "sha256": h,
           "rank2curves": len(src_lines), "rank2classes": len(classes)}]},
          open(os.path.join(mini, "source_counts.json"), "w"))
lab0 = "".join(src_lines[0].split()[:3])
lab1 = "".join(src_lines[1].split()[:3])
env = dict(os.environ, MWSHAPE_INJECT_JSON=json.dumps({"%s@%d" % (lab0, b): "exception" for b in (192, 256, 384, 512)}))
r1 = subprocess.run([sys.executable, os.path.join(HERE, "run_study.py"), mini, os.path.join(mini, "source_counts.json"), out,
                     "--workers", "1"], capture_output=True, text=True, env=env)
# also delete one record from the export to test MISSING handling
recs = [json.loads(l) for l in open(os.path.join(out, "per_curve_results.jsonl"))]
with open(os.path.join(out, "per_curve_results.jsonl"), "w") as f:
    for rr in recs:
        if rr.get("label") != lab1:
            f.write(json.dumps(rr) + "\n")
r2 = subprocess.run([sys.executable, os.path.join(HERE, "analyze.py"), out], capture_output=True, text=True)
summ = json.load(open(os.path.join(out, "summary.json")))
ledger = open(os.path.join(out, "retry_failure_ledger.csv")).read()
A = summ["A_all_rank2_records"]["overall"]
ok = (r1.returncode == 1 and r2.returncode == 1 and summ["status"] == "INCOMPLETE" and lab0 in ledger and lab1 in ledger
      and "MISSING" in ledger and A["expected"] == len(src_lines) and A["completed"] == len(src_lines) - 2
      and A["unresolved"] == 2 and "mean_s" in A and os.path.exists(os.path.join(out, "per_curve_shapes.csv")))
record("unresolved export: run_study exits 1, analyze exits 1 with INCOMPLETE, ledger lists UNRESOLVED and MISSING, "
       "resolved-only summaries still written", ok,
       {"run_study_exit": r1.returncode, "analyze_exit": r2.returncode, "status": summ["status"],
        "A_expected_completed_unresolved": [A["expected"], A["completed"], A["unresolved"]],
        "unresolved_label_injected": lab0, "removed_record_label": lab1,
        "ledger": ledger.strip().splitlines()}, True)
# clean source block without injection -> both scripts exit 0
out2 = os.path.join(WORK, "mini_out_clean"); shutil.rmtree(out2, ignore_errors=True)
r3 = subprocess.run([sys.executable, os.path.join(HERE, "run_study.py"), mini, os.path.join(mini, "source_counts.json"), out2,
                     "--workers", "1"], capture_output=True, text=True)
r4 = subprocess.run([sys.executable, os.path.join(HERE, "analyze.py"), out2], capture_output=True, text=True)
s2 = json.load(open(os.path.join(out2, "summary.json")))
record("positive control: same 12-record block without injection -> run_study exit 0, analyze exit 0, COMPLETE",
       r3.returncode == 0 and r4.returncode == 0 and s2["status"] == "COMPLETE",
       {"run_study_exit": r3.returncode, "analyze_exit": r4.returncode, "status": s2["status"]}, False)

# --- 7. negative_controls.py and selftests.py exit codes are checked by run_all.sh (set -e)
json.dump({"tests": T, "failed": fails, "status": "PASS" if not fails else "FAIL"}, open(OUTJ, "w"), indent=1)
print("integration tests:", len(T), "run,", len(fails), "failed")
sys.exit(1 if fails else 0)
