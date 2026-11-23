#!/usr/bin/env python3
"""Negative controls for the per-curve pipeline (run_study.process). Each must NOT be RESOLVED.
Also one positive control (unmodified 389a1 row) that must be RESOLVED.
Usage: python3 negative_controls.py OUT_JSON"""
import sys, json
import run_study as rs

OK389 = "389 a 1 [0,1,1,-2,0] 2 [] [0:0:1] [1:0:1]"
cases = {
    "positive control 389a1 (unmodified)": (OK389, "RESOLVED"),
    "off-curve generator (y+1)": ("389 a 1 [0,1,1,-2,0] 2 [] [0:1:1] [1:0:1]", "UNRESOLVED"),
    "dependent generators (P, 2P)": ("389 a 1 [0,1,1,-2,0] 2 [] [0:0:1] [3:5:1]", "UNRESOLVED"),
    "repeated generator (P, P)": ("389 a 1 [0,1,1,-2,0] 2 [] [0:0:1] [0:0:1]", "UNRESOLVED"),
    "torsion point given as free generator (1088j1: (3,0) has order 2)": ("1088 j 1 [0,1,0,-25,39] 2 [] [1:4:1] [3:0:1]", "UNRESOLVED"),
    "torsion structure parsed as point (missing torsion generator)": ("1088 j 1 [0,1,0,-25,39] 2 [2] [1:4:1] [5:8:1]", "UNRESOLVED"),
    "wrong projective convention (817a1: intended point (-3/32,15/512) = ordinary triple [-48:15:512])": ("817 a 1 [0,1,1,1,6] 2 [] [4:9:1] [-48:15:512]", "UNRESOLVED"),
}
out, bad = [], 0
for name, (line, expect) in cases.items():
    rec = rs.process(line + "\n")
    got = rec.get("status")
    ok = (got == expect) if expect == "RESOLVED" else (got != "RESOLVED" and not str(got).startswith("RESOLVED"))
    bad += not ok
    out.append({"case": name, "line": line, "expected": expect, "got": got, "reason": rec.get("reason"), "pass": ok})
    print("%-75s expected %-10s got %-26s %s" % (name, expect, got, "OK" if ok else "*** FAIL ***"))
json.dump({"cases": out, "status": "PASS" if not bad else "FAIL"}, open(sys.argv[1], "w"), indent=1)
sys.exit(1 if bad else 0)
