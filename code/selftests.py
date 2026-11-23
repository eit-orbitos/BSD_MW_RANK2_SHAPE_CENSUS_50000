#!/usr/bin/env python3
"""Pre-run self-tests (must all pass before the full run).

1. ecdata format regressions: 817a1 (non-integral generator, Z=8) and 1088j1 (torsion
   structure [2] precedes points; torsion generator follows free generators).
   Negative controls: the wrong (X/Z^2, Y/Z^3) convention must fail; malformed rows must raise.
2. Reference lattices: square, hexagonal, and unimodularly sheared versions (exact rationals
   and 256-bit reals): reduction must recover the same sign-invariant shape (h1, h2, |b|).
3. Historical 19-row comparison against reference/reduced_basis_shape_extract.csv
   (rounded inputs, 15 significant digits): |s_new - s_hist| < 1e-12 and relative
   agreement of h1, h2, |b|, Reg < 1e-12. The historical CSV is read, never modified.

Usage: python3 selftests.py ALLGENS_DIR HISTORICAL_CSV OUT_JSON
"""
import sys, json, csv
from fractions import Fraction
import mwshape as m

pari = m.pari
results = {"tests": [], "status": None}
fail = []


def check(name, ok, detail=None):
    results["tests"].append({"name": name, "pass": bool(ok), "detail": detail})
    if not ok:
        fail.append(name)


# ---- 1. format regressions
R817 = "817 a 1 [0,1,1,1,6] 2 [] [4:9:1] [-6:15:8]"
R1088 = "1088 j 1 [0,1,0,-25,39] 2 [2] [1:4:1] [5:8:1] [3:0:1]"
r = m.parse_row(R817)
check("817a1 generator 2 parsed as (-3/4, 15/8)", (r["free"][1]["x"], r["free"][1]["y"]) == (Fraction(-3, 4), Fraction(15, 8)))
check("817a1 both free generators on curve (exact)", all(m.on_curve(r["ainvs"], p) for p in r["free"]))
wrong = {"x": Fraction(-6, 8 ** 2), "y": Fraction(15, 8 ** 3)}
check("817a1 NEGATIVE: (X/Z^2, Y/Z^3) convention gives off-curve point", not m.on_curve(r["ainvs"], wrong),
      "(-3/32, 15/512)")
r = m.parse_row(R1088)
check("1088j1 torsion structure parsed as [2], not as a point", r["torsion"] == [2])
check("1088j1 free generators are [1:4:1],[5:8:1]", [p["token"] for p in r["free"]] == ["[1:4:1]", "[5:8:1]"])
check("1088j1 torsion generator [3:0:1] separate and on curve", [p["token"] for p in r["torsion_gens"]] == ["[3:0:1]"]
      and m.on_curve(r["ainvs"], r["torsion_gens"][0]))
check("1088j1 free generators on curve (exact)", all(m.on_curve(r["ainvs"], p) for p in r["free"]))
for bad, why in [("1088 j 1 [0,1,0,-25,39] 2 [2] [1:4:1] [5:8:1]", "missing torsion generator"),
                 ("817 a 1 [0,1,1,1,6] 2 [] [4:9:1]", "missing free generator"),
                 ("817 a 1 [0,1,1,1] 2 [] [4:9:1] [-6:15:8]", "4 a-invariants"),
                 ("817 a 1 [0,1,1,1,6] 2 [] [4:9:0] [-6:15:8]", "Z = 0")]:
    try:
        m.parse_row(bad)
        check("NEGATIVE parse rejects: " + why, False)
    except m.ParseError:
        check("NEGATIVE parse rejects: " + why, True)

# ---- 2. reference lattices
def shape_inv(G):
    Gr, U, _ = m.gauss_reduce(G)
    h1, h2, ab, det, x, y, s = m.shape_from_reduced(Gr)
    return Gr, U, (h1, h2, ab), (x, y, s)

refs = {
    "square": (pari("[1,0;0,1]"), {"s": 1, "x": 0, "y": 1}),
    "hexagonal": (pari("[1,1/2;1/2,1]"), {"s": pari("3/4"), "x": pari("1/2"), "y": pari.sqrt(pari(3), precision=256) / 2}),
    "generic_reduced": (pari("[2,1/3;1/3,5]"), None),
}
shears = [pari("[1,7;0,1]"), pari("[2,1;1,1]"), pari("[5,2;2,1]"), pari("[0,1;1,0]"), pari("[-3,8;1,-3]")]
for name, (G0, expect) in refs.items():
    for exact in (True, False):
        Gbase = G0 if exact else G0 * pari.sqrt(pari(1), precision=256)   # 256-bit real copy
        Gr0, U0, inv0, shp0 = shape_inv(Gbase)
        if expect:
            ok = all(abs(shp0[i] - expect[k]) < pari("1e-60") for i, k in enumerate(["x", "y", "s"]))
            check("reference %s (%s): x, y, s as expected" % (name, "exact" if exact else "256-bit"), ok,
                  {k: str(v) for k, v in zip(["x", "y", "s"], shp0)})
        for M in shears:
            Gs = M.mattranspose() * Gbase * M
            Gr, U, inv, shp = shape_inv(Gs)
            okU = abs(pari.matdet(U)) == 1
            okT = (U.mattranspose() * Gs * U - Gr) == 0 if exact else m.frob(U.mattranspose() * Gs * U - Gr) < pari("1e-70")
            okI = all(abs(inv[k] - inv0[k]) <= pari("1e-70") * (1 + abs(inv0[k])) for k in range(3))
            okR = bool(2 * inv[2] <= inv[0] and inv[0] <= inv[1])
            check("reference %s (%s) sheared by %s: |detU|=1, Gred=U~GU, same (h1,h2,|b|), inequalities"
                  % (name, "exact" if exact else "256-bit", str(M)), okU and okT and okI and okR)
# boundary flags on the exact boundary cases
Gr, U, inv, shp = shape_inv(pari("[1,1/2;1/2,1]"))
check("hexagonal: both margins exactly 0 (boundary flagged)", (inv[0] - 2 * inv[2]) == 0 and (inv[1] - inv[0]) == 0)
Gr, U, inv, shp = shape_inv(pari("[1,0;0,1]"))
check("square: h1 == h2 margin exactly 0 and x exactly 0", (inv[1] - inv[0]) == 0 and shp[0] == 0)

# ---- 3. historical 19-row comparison (before full run)
allgens_dir, hist_csv, out_json = sys.argv[1], sys.argv[2], sys.argv[3]
hist = {r["label"]: r for r in csv.DictReader(open(hist_csv, newline=""))}
rows = {}
for blk in ("00000-09999",):
    for line in open("%s/allgens.%s" % (allgens_dir, blk)):
        t = line.split()
        lab = t[0] + t[1] + t[2]
        if lab in hist:
            rows[lab] = m.parse_row(line)
check("all 19 historical labels present in pinned allgens 00000-09999", set(rows) == set(hist), sorted(set(hist) - set(rows)))
hist_detail = {}
for lab, hr in sorted(hist.items()):
    o, n = m.compute_at(rows[lab], 192)
    h1, h2, ab = n["inv"]
    d = {
        "s_abs": abs(n["s"] - pari(hr["sin2_reduced"])),
        "h1_rel": abs(h1 / pari(hr["h1"]) - 1),
        "h2_rel": abs(h2 / pari(hr["h2"]) - 1),
        "absb_rel": abs(ab / abs(pari(hr["pairing"])) - 1),
        "reg_rel": abs(n["det"] / pari(hr["regulator"]) - 1),
    }
    ok = d["s_abs"] < pari("1e-12") and all(d[k] < pari("1e-12") for k in ("h1_rel", "h2_rel", "absb_rel", "reg_rel"))
    hist_detail[lab] = {k: m.sig(v, 4) for k, v in d.items()}
    hist_detail[lab]["U"] = o["U"]
    check("historical %s agrees with V3.15.4 CSV at rounded-input tolerance 1e-12" % lab, ok, hist_detail[lab])

results["historical_19"] = hist_detail
results["status"] = "PASS" if not fail else "FAIL"
results["failed"] = fail
json.dump(results, open(out_json, "w"), indent=1)
print("self-tests:", len(results["tests"]), "run,", len(fail), "failed ->", results["status"])
for f in fail:
    print("  FAILED:", f)
sys.exit(0 if not fail else 1)
