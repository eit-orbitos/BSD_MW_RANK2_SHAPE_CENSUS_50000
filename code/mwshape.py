#!/usr/bin/env python3
"""Core library: rank-two Mordell-Weil lattice-shape study (BSD L2, study candidate V0.2).

V0.2 repairs (Claude-authored from the written follow-up specification of 2026-11-23; the
supplied repair-candidate files were not received):
  * sign-invariant shape errors are normalized by the positive diagonal scale D = diag(sqrt h1, sqrt h2)
    of the higher-precision reduced Gram matrix (scaled_shape_error), also in the shear controls,
    so an exact or tiny pairing b no longer divides by zero or inflates a relative error;
  * test-only fault injection (INJECT) for labelled injected controls; empty in production.

Parsing of Cremona ecdata `allgens` rows, exact point membership, and Gauss
(Lagrange) reduction of the canonical-height Gram matrix with an exact integer
transformation U, at an explicitly requested PARI working precision.

IMPORTANT precision note (verified in this environment, cypari2 2.2.9 / PARI 2.19.1):
cypari2 method calls such as pari.ellheightmatrix(E, P) do NOT follow
pari.set_real_precision_bits(); they keep the precision of the first call.
Every precision-dependent call below therefore passes `precision=bits`
explicitly, and the effective precision of each result is read back with
pari.bitprecision() and recorded.
"""
import re
import json
from fractions import Fraction

import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10**9, silent=True)

_FMT = pari('(d,x)->strprintf(Str("%.",d,"e"),x)')  # PARI-side formatting: no float conversion
TRIPLE = re.compile(r"\[(-?\d+):(-?\d+):(-?\d+)\]")

# Height normalization: PARI >= 2.7 ellheight / ellheightmatrix, which agrees with the
# BSD regulator convention used by Cremona, Sage and LMFDB (e.g. Reg(389a1) = 0.152460177943...).
HEIGHT_NORMALIZATION = "PARI ellheightmatrix (PARI >= 2.7 normalization, = BSD/Cremona/LMFDB regulator convention)"


class ParseError(ValueError):
    pass


def parse_row(line):
    """Parse one allgens row. Returns a dict with exact data; raises ParseError."""
    parts = line.split()
    if len(parts) < 6:
        raise ParseError("too few fields")
    N, cls, num, rank = int(parts[0]), parts[1], int(parts[2]), int(parts[4])
    if not re.fullmatch(r"[a-z]+", cls):
        raise ParseError("bad class code")
    ainvs = json.loads(parts[3])
    if not (isinstance(ainvs, list) and len(ainvs) == 5 and all(type(a) is int for a in ainvs)):
        raise ParseError("bad a-invariants")
    torsion = json.loads(parts[5])  # torsion STRUCTURE, e.g. [] or [2] or [2,4]; never a point
    if not (isinstance(torsion, list) and all(type(t) is int and t >= 2 for t in torsion) and len(torsion) <= 2):
        raise ParseError("bad torsion structure")
    tokens = parts[6:]
    if len(tokens) != rank + len(torsion):
        raise ParseError("point count %d != rank %d + torsion gens %d" % (len(tokens), rank, len(torsion)))
    pts = []
    for tok in tokens:
        m = TRIPLE.fullmatch(tok)
        if not m:
            raise ParseError("bad projective triple %r" % tok)
        X, Y, Z = (int(v) for v in m.groups())
        if Z == 0:
            raise ParseError("point at infinity listed as generator")
        pts.append({"token": tok, "x": Fraction(X, Z), "y": Fraction(Y, Z)})  # ordinary projective: (X/Z, Y/Z)
    return {"N": N, "cls": cls, "num": num, "label": "%d%s%d" % (N, cls, num), "rank": rank,
            "ainvs": ainvs, "torsion": torsion, "free": pts[:rank], "torsion_gens": pts[rank:]}


def residual(ainvs, x, y):
    a1, a2, a3, a4, a6 = ainvs
    return y * y + a1 * x * y + a3 * y - (x ** 3 + a2 * x * x + a4 * x + a6)


def on_curve(ainvs, pt):
    return residual(ainvs, pt["x"], pt["y"]) == 0


def frac_str(q):
    return str(q.numerator) if q.denominator == 1 else "%d/%d" % (q.numerator, q.denominator)


def pari_point(pt):
    return pari("[%s,%s]" % (frac_str(pt["x"]), frac_str(pt["y"])))


def sig(x, digits=45):
    """Serialize a PARI real/exact number with `digits` significant digits (exact zero kept as '0')."""
    if x == 0:
        # exact zero -> "0"; a PARI real zero keeps its exponent (e.g. "0.E-58"), i.e. "below working precision"
        return "0" if str(pari.type(x)) != "t_REAL" else str(x).replace(" ", "")
    return str(_FMT(digits - 1, x)).replace(" ", "")


def gauss_reduce(G, max_steps=100000):
    """Lagrange/Gauss reduction of a positive-definite 2x2 Gram matrix G (PARI matrix).

    Returns (Gred, U, steps) with U an exact integer matrix, Gred = U~ * G * U computed
    incrementally, satisfying h1 <= h2 and 2|b| <= h1. Ties: no swap when h1 == h2;
    no translation when 2|b| <= h1 (so a pairing exactly at +-h1/2 is kept as is).
    """
    A = G
    U = pari.matid(2)
    SW = pari("[0,1;1,0]")
    for step in range(max_steps):
        if A[0, 0] > A[1, 1]:
            A = SW * A * SW
            U = U * SW
        if 2 * abs(A[0, 1]) <= A[0, 0]:
            return A, U, step
        m = (A[0, 1] / A[0, 0]).round()
        T = pari("[1,%s;0,1]" % (-m))
        A = T.mattranspose() * A * T
        U = U * T
    raise RuntimeError("gauss_reduce: no convergence")


def frob(M):
    return pari.sqrt(sum(M[i, j] ** 2 for i in range(2) for j in range(2)))


def shape_from_reduced(Gred, bits=256):
    h1, h2, b = Gred[0, 0], Gred[1, 1], Gred[0, 1]
    det = h1 * h2 - b * b
    x = abs(b) / h1
    # For t_REAL det, PARI uses det's own precision; `bits` only matters for exact rational input.
    y = pari.sqrt(det, precision=bits) / h1
    s = det / (h1 * h2)            # = y^2/(x^2+y^2)
    return h1, h2, abs(b), det, x, y, s


# Test-only fault injection: {(label, bits): spec}. Never populated by run_study.py in production.
#   "exception"        -> raise inside the computation
#   "perturb:<eps>"    -> add eps*G[1,1]-relative noise to G[0,1] (and its mirror) before reduction
#   "set_b:<value>"    -> replace the off-diagonal entry by an exact value (e.g. 0)
#   "b_over_h1:<r>"    -> replace the off-diagonal entry by r*G[0,0]
# Tests may also pass MWSHAPE_INJECT_JSON='{"<label>@<bits>": "<spec>", ...}' in the environment.
INJECT = {}
import os as _os
for _k, _v in json.loads(_os.environ.get("MWSHAPE_INJECT_JSON", "{}")).items():   # tests only: {"label@bits": spec}
    _lab, _b = _k.rsplit("@", 1)
    INJECT[(_lab, int(_b))] = _v


def _apply_injection(row, bits, G):
    spec = INJECT.get((row["label"], bits))
    if spec is None:
        return G, None
    if spec == "exception":
        raise RuntimeError("INJECTED exception at %d bits" % bits)
    G = pari.matrix(2, 2, [G[0, 0], G[0, 1], G[1, 0], G[1, 1]])
    if spec.startswith("perturb:"):
        e = pari(spec.split(":", 1)[1]) * G[1, 1]
        G[0, 1] = G[0, 1] + e; G[1, 0] = G[0, 1]
    elif spec.startswith("b_over_h1:"):
        v = pari(spec.split(":", 1)[1]) * G[0, 0]
        G[0, 1] = v; G[1, 0] = v
    elif spec.startswith("set_b:"):
        v = pari(spec.split(":", 1)[1])
        G[0, 1] = v; G[1, 0] = v
    return G, spec


def compute_at(row, bits, check_reduced_points=True):
    """Everything precision-dependent for one curve at requested `bits`. Rebuilds E each time."""
    E = pari.ellinit(row["ainvs"], precision=bits)   # fresh curve object: no reuse of cached lower-precision data
    P = [pari_point(p) for p in row["free"]]
    G = pari.ellheightmatrix(E, P, precision=bits)
    G, injected = _apply_injection(row, bits, G)
    out = {"requested_bits": bits}
    if injected:
        out["INJECTED_CONTROL"] = injected
    # positive definiteness with a relative tolerance: det(G) > 1e-30 * G11 * G22
    # (a dependent pair such as (P, 2P) has det(G) = 0 up to rounding and must be rejected)
    pd = bool(G[0, 0] > 0 and G[1, 1] > 0 and pari.matdet(G) > pari("1e-30") * G[0, 0] * G[1, 1])
    out["raw_positive_definite"] = pd
    if not pd:
        out["error"] = "raw Gram matrix not positive definite"
        return out, None
    Gred, U, steps = gauss_reduce(G)
    detU = pari.matdet(U)
    Gchk = U.mattranspose() * G * U
    detG = pari.matdet(G)
    h1, h2, ab, det, x, y, s = shape_from_reduced(Gred, bits)
    # actual effective precision: minimum over the raw/reduced Gram entries and derived reals
    eff = [int(pari.bitprecision(v)) for v in (G[0, 0], G[0, 1], G[1, 1], Gred[0, 0], Gred[0, 1], Gred[1, 1], det, y, s)
           if str(pari.type(v)) == "t_REAL"]
    out["effective_bits"] = min(eff)
    out["effective_bits_below_requested"] = bool(min(eff) < bits)
    out.update({
        "G": [sig(G[0, 0]), sig(G[0, 1]), sig(G[1, 1])],
        "U": [[int(U[0, 0]), int(U[0, 1])], [int(U[1, 0]), int(U[1, 1])]],
        "det_U": int(detU), "reduction_steps": int(steps),
        "G_red": [sig(Gred[0, 0]), sig(Gred[0, 1]), sig(Gred[1, 1])],
        "rel_err_Gred_vs_UtGU": sig(frob(Gred - Gchk) / frob(Gchk), 6),
        "rel_err_detGred_vs_detG": sig(abs(pari.matdet(Gred) - detG) / detG, 6),
        "reg": sig(detG), "x": sig(x), "y": sig(y), "s": sig(s),
        "margin_2b": sig((h1 - 2 * ab) / h1, 12),     # >= 0; 0 = on the |x| = 1/2 boundary
        "margin_h": sig((h2 - h1) / h2, 12),          # >= 0; 0 = on the unit-circle boundary
        "red_positive_definite": bool(h1 > 0 and det > pari("1e-30") * h1 * h2),
        "red_inequalities": bool(2 * ab <= h1 and h1 <= h2),
    })
    if check_reduced_points and not injected:
        # independent check of U: build reduced points exactly and recompute their Gram matrix
        Q = [pari.elladd(E, pari.ellmul(E, P[0], U[0, j]), pari.ellmul(E, P[1], U[1, j])) for j in range(2)]
        Gq = pari.ellheightmatrix(E, Q, precision=bits)
        out["rel_err_Gred_vs_heights_of_reduced_points"] = sig(frob(Gq - Gred) / frob(Gred), 6)
        out["reduced_points"] = [str(q) for q in Q]
    elif injected:
        out["rel_err_Gred_vs_heights_of_reduced_points"] = "not-applicable(injected)"
    num = {"G": G, "Gred": Gred, "U": U, "det": detG, "x": x, "y": y, "s": s,
           "inv": (h1, h2, ab), "injected": injected}
    return out, num


def scaled_shape_error(inv_a, inv_b):
    """Sign-invariant shape error between reduced (h1, h2, |b|) triples, normalized by the positive
    diagonal scale D = diag(sqrt h1, sqrt h2) of the reference triple inv_b:
        max(|dh1|/h1, |dh2|/h2, |d|b||/sqrt(h1*h2)),
    i.e. the max-entry difference of D^-1 G_red D^-1. Well defined when b = 0 exactly."""
    h1, h2, _ = inv_b
    return max(abs(inv_a[0] - h1) / h1, abs(inv_a[1] - h2) / h2,
               abs(inv_a[2] - inv_b[2]) / pari.sqrt(h1 * h2, precision=256))


def compare(a, b):
    """Stability comparison between two precisions (b = higher). Returns dict of PARI reals."""
    rel = lambda u, v: abs(u - v) / abs(v)
    same_U = str(a["U"]) == str(b["U"])
    d = {
        "gram_rel_frobenius": frob(a["G"] - b["G"]) / frob(b["G"]),
        "det_rel": rel(a["det"], b["det"]),
        "x_abs": abs(a["x"] - b["x"]),
        "s_abs": abs(a["s"] - b["s"]),
        "y_rel": rel(a["y"], b["y"]),
        # sign-invariant reduced shape (h1, h2, |b|), diagonal-scale normalized (V0.2 definition);
        # valid even when U differs near a boundary and when b = 0
        "inv_scaled_max": scaled_shape_error(a["inv"], b["inv"]),
    }
    if same_U:
        d["gred_rel_frobenius"] = frob(a["Gred"] - b["Gred"]) / frob(b["Gred"])
    return same_U, d


SHEAR_MATRICES = ("[2,1;1,1]", "[1,5;0,1]", "[3,-7;-2,5]")   # all det +1


def shear_control(row, bits, num_ref):
    """Control on real data: feed deliberately non-reduced bases M-combinations of the supplied
    generators, recompute their Gram matrix from the points, reduce, and require the same
    sign-invariant shape (h1, h2, |b|) as the reference reduction (relative difference < 1e-40)."""
    if num_ref.get("injected"):
        # the reference Gram matrix was artificially modified: comparing real sheared bases with it is meaningless
        return [{"M": Ms, "pass": "not-applicable(injected reference)"} for Ms in SHEAR_MATRICES]
    E = pari.ellinit(row["ainvs"], precision=bits)
    P = [pari_point(p) for p in row["free"]]
    res = []
    for Ms in SHEAR_MATRICES:
        M = pari(Ms)
        Q = [pari.elladd(E, pari.ellmul(E, P[0], M[0, j]), pari.ellmul(E, P[1], M[1, j])) for j in range(2)]
        Gs = pari.ellheightmatrix(E, Q, precision=bits)
        Gr, U, steps = gauss_reduce(Gs)
        h1, h2, ab, det, x, y, s = shape_from_reduced(Gr, bits)
        d = scaled_shape_error((h1, h2, ab), num_ref["inv"])
        res.append({"M": Ms, "det_U": int(pari.matdet(U)), "steps": int(steps),
                    "U_nontrivial": str(U) != str(pari.matid(2)),
                    "max_scaled_diff_invariants": sig(d, 6), "pass": bool(d < pari("1e-40") and abs(int(pari.matdet(U))) == 1)})
    return res
