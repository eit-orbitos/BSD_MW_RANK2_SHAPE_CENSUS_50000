#!/usr/bin/env python3
"""Plots for the rank-two lattice-shape study (floats used only here, for drawing).
Usage: python3 plots.py OUTDIR"""
import sys, os, json, csv, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import mpmath as mp

OUT = sys.argv[1]
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)
# reference palette (dataviz skill, slots 1-3 validated all-pairs; aqua < 3:1 -> direct labels + tables)
SURF, INK, INK2, GRID, REF = "#fcfcfb", "#0b0b0b", "#52514e", "#e4e3df", "#52514e"
C_A, C_B, C_W = "#2a78d6", "#eb6834", "#1baf7a"
plt.rcParams.update({"figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
                     "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2,
                     "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "lines.linewidth": 2})

rows = [r for r in csv.DictReader(open(os.path.join(OUT, "per_curve_shapes.csv"))) if r["status"].startswith("RESOLVED")]
A = {"s": np.array([float(r["s"]) for r in rows]), "x": np.array([float(r["x"]) for r in rows]),
     "y": np.array([float(r["y"]) for r in rows]), "N": np.array([int(r["N"]) for r in rows])}
mask_B = np.array([r["in_population_B"] == "True" for r in rows])
B = {k: v[mask_B] for k, v in A.items()}
summ = json.load(open(os.path.join(OUT, "summary.json")))
REF_MEAN = float(summ["reference"]["mean_s_closed_form_(3/pi)Cl2(pi/3)"])
WAT = {int(k): float(v) / 100 for k, v in summ["reference"]["Watkins_2008_section_6.4_percent_HISTORICAL_DIFFERENT_POPULATION"].items()}

mp.mp.dps = 30
def Fref(t):
    if t <= 0.75: return 0.0
    if t >= 1: return 1.0
    a = mp.sqrt(1 - t)
    return float((6 / mp.pi) * (mp.pi / 6 - mp.asin(a) - mp.sqrt((1 - t) / t) * mp.log(1 / (2 * a))))

def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIG, name + "." + ext), dpi=160, bbox_inches="tight")
    plt.close(fig)

def ecdf(v):
    v = np.sort(v); return v, np.arange(1, len(v) + 1) / len(v)

# 1. ECDF of s vs reference CDF
t = np.linspace(0.75, 1.0, 600); Ft = np.array([Fref(u) for u in t])
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.plot(t, Ft, color=REF, ls="--", lw=2, label="Haar reference $F_s$")
for d, c, lab in ((A, C_A, "A: all rank-2 records (n = %d)" % len(A["s"])), (B, C_B, "B: class representatives (n = %d)" % len(B["s"]))):
    xs, ys = ecdf(d["s"]); ax.step(xs, ys, where="post", color=c, lw=2, label=lab)
ax.set_xlim(0.75, 1.0); ax.set_ylim(0, 1)
ax.set_xlabel(r"$s=\sin^2\theta_{\mathrm{red}}$"); ax.set_ylabel("cumulative proportion")
ax.set_title("Empirical CDF of the reduced angle factor, N < 50,000", loc="left", fontsize=11)
ax.legend(frameon=False, loc="upper left")
ax.annotate("reference mean %.6f" % REF_MEAN, xy=(0.80, 0.55), color=INK2, fontsize=9)
ax.annotate("A mean %.6f" % A["s"].mean(), xy=(0.80, 0.49), color=INK2, fontsize=9)
save(fig, "fig1_ecdf_s_vs_reference")

# 2. ECDF minus reference
fig, ax = plt.subplots(figsize=(7.2, 3.6))
ax.axhline(0, color=REF, lw=1)
for d, c, lab in ((A, C_A, "A"), (B, C_B, "B")):
    xs, ys = ecdf(d["s"]); ax.plot(xs, ys - np.array([Fref(u) for u in xs]), color=c, lw=2, label=lab)
    i = np.argmax(np.abs(ys - np.array([Fref(u) for u in xs])))
ax.set_xlim(0.75, 1.0); ax.set_xlabel(r"$s$"); ax.set_ylabel(r"$F_{\mathrm{emp}}(s)-F_s(s)$")
ax.set_title("Empirical minus reference CDF (descriptive; not a test statistic)", loc="left", fontsize=11)
ax.legend(frameon=False)
save(fig, "fig2_ecdf_minus_reference")

# 3. folded fundamental domain, small multiples A | B
fig, axs = plt.subplots(1, 2, figsize=(9.6, 5.2), sharey=True)
xx = np.linspace(0, 0.5, 200)
for ax, d, c, lab in ((axs[0], A, C_A, "A: all records"), (axs[1], B, C_B, "B: class representatives")):
    ax.scatter(d["x"], d["y"], s=2, color=c, alpha=0.35, linewidths=0, rasterized=True)
    ax.plot(xx, np.sqrt(1 - xx ** 2), color=INK2, lw=1.5)
    ax.axhline(1, color=REF, ls="--", lw=1)
    ax.set_xlim(0, 0.5); ax.set_ylim(0.85, 9)
    ax.set_xlabel(r"$x=|b|/h_1$"); ax.set_title("%s (n = %d)" % (lab, len(d["x"])), loc="left", fontsize=10)
axs[0].set_ylabel(r"$y=\sqrt{\det G}/h_1$")
axs[0].annotate("y = 1", xy=(0.43, 1.05), color=INK2, fontsize=8)
axs[0].annotate(r"$x^2+y^2=1$", xy=(0.30, 0.89), color=INK2, fontsize=8)
fig.suptitle("Reduced rank-2 Mordell–Weil lattice shapes in the folded fundamental domain", x=0.06, ha="left", fontsize=11)
save(fig, "fig3_fundamental_domain")

# 4. tail of y: P(y >= Y)
Y = np.logspace(0, math.log10(16), 200)
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.loglog(Y, 3 / (math.pi * Y), color=REF, ls="--", lw=2, label=r"Haar reference $3/(\pi Y)$")
for d, c, lab in ((A, C_A, "A: all records"), (B, C_B, "B: class representatives")):
    pr = np.array([(d["y"] >= u).mean() for u in Y]); keep = pr > 0
    ax.loglog(Y[keep], pr[keep], color=c, lw=2, label=lab)
wy = sorted(WAT); ax.loglog(wy, [WAT[k] for k in wy], ls="none", marker="D", ms=8, color=C_W,
                            markeredgecolor=SURF, markeredgewidth=2, label="Watkins 2008 §6.4 (historical, different population)")
for k in wy:
    ax.annotate("%.3g%%" % (100 * WAT[k]), xy=(k, WAT[k]), xytext=(6, 4), textcoords="offset points", fontsize=8, color=INK2)
ax.set_xlabel("Y"); ax.set_ylabel(r"proportion with $y \geq Y$")
ax.set_title("Tail of the shape height y", loc="left", fontsize=11)
ax.legend(frameon=False, loc="lower left", fontsize=8.5)
save(fig, "fig4_y_tail")

# 5. per-block mean s and 6. per-block P(y>=1)  (two charts, no dual axis)
blocks = list(summ["A_all_rank2_records"]["blocks"])
mids = np.arange(len(blocks))
for key, fname, ylab, ref in (("mean_s", "fig5_block_mean_s", r"mean of $s$", REF_MEAN),
                              ("p1", "fig6_block_P_y_ge_1", r"proportion with $y\geq 1$", 3 / math.pi)):
    fig, ax = plt.subplots(figsize=(7.2, 3.8))
    ax.axhline(ref, color=REF, ls="--", lw=1.5)
    ax.annotate("Haar reference %.6f" % ref, xy=(len(blocks) - 1.3, ref), xytext=(0, 4), textcoords="offset points", fontsize=8, color=INK2)
    for pop, c, lab, off in (("A_all_rank2_records", C_A, "A", -0.08), ("B_selected_class_representatives", C_B, "B", 0.08)):
        v = [float(summ[pop]["blocks"][b]["mean_s"]) if key == "mean_s" else float(summ[pop]["blocks"][b]["prop_y_ge"]["1"]) for b in blocks]
        ax.plot(mids + off, v, ls="none", marker="o", ms=8, color=c, markeredgecolor=SURF, markeredgewidth=2, label=lab)
    if key == "p1":
        ax.axhline(WAT[1], color=C_W, ls=":", lw=1.5)
        ax.annotate("Watkins 93.5% (historical, different population)", xy=(0, WAT[1]), xytext=(0, -12), textcoords="offset points", fontsize=8, color=INK2)
    ax.set_xticks(mids); ax.set_xticklabels([b.replace("-", "–") for b in blocks], fontsize=8.5)
    ax.set_xlabel("conductor block"); ax.set_ylabel(ylab)
    ax.set_title(("Mean reduced angle factor by conductor block" if key == "mean_s" else "Proportion with y ≥ 1 by conductor block"), loc="left", fontsize=11)
    ax.legend(frameon=False, loc="center right")
    save(fig, fname)
print("figures written to", FIG)
