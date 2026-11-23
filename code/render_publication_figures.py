#!/usr/bin/env python3
"""Regenerate presentation figures from unchanged V0.2 CSV records.

Usage: python3 code/render_publication_figures.py STUDY_ROOT [OUTPUT_ROOT]
Floating-point conversion is used for plotting only, never to revise data.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import mpmath as mp
import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('study_root', type=Path)
parser.add_argument('output_root', type=Path, nargs='?', default=Path(__file__).resolve().parents[1])
args = parser.parse_args()
source = args.study_root / 'results/per_curve_shapes.csv'
before = hashlib.sha256(source.read_bytes()).hexdigest()
with source.open(newline='') as f:
    rows = list(csv.DictReader(f))
if len(rows) != 23612 or any(not r['status'].startswith('RESOLVED') for r in rows):
    raise ValueError('Expected complete, resolved V0.2 population')

def population(selected):
    return {key: np.array([float(r[key]) for r in selected]) for key in ['s', 'x', 'y']}

A = population(rows)
B = population([r for r in rows if r['in_population_B'] == 'True'])
if len(B['s']) != 18886:
    raise ValueError('Unexpected population B count')

figdir = args.output_root / 'figures'
figdir.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({
    'font.family': 'DejaVu Sans', 'font.size': 10,
    'axes.spines.top': False, 'axes.spines.right': False,
    'axes.grid': True, 'grid.alpha': 0.25,
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
    'pdf.fonttype': 42, 'ps.fonttype': 42,
})
blue, orange, green, gray = '#1763aa', '#c85811', '#18724f', '#454545'

def save(fig, name):
    for ext in ['pdf', 'png']:
        kw = {'metadata': {'CreationDate': None, 'ModDate': None}} if ext == 'pdf' else {}
        fig.savefig(figdir / f'{name}.{ext}', dpi=200, bbox_inches='tight', pad_inches=0.12, **kw)
    plt.close(fig)

mp.mp.dps = 40
def cdf(t):
    if t <= 0.75:
        return 0.
    if t >= 1:
        return 1.
    t = mp.mpf(float(t))
    a = mp.sqrt(1-t)
    return float(6/mp.pi * (mp.pi/6 - mp.asin(a) - mp.sqrt((1-t)/t)*mp.log(1/(2*a))))

fig, ax = plt.subplots(figsize=(7.2, 4.5), layout='constrained')
t = np.linspace(.75, 1, 900)
ax.plot(t, [cdf(v) for v in t], '--', color=gray, lw=1.8, label='Hyperbolic reference')
for p, color, ls, label in [(A, blue, '-', 'A: all curves (23,612)'), (B, orange, ':', 'B: class representatives (18,886)')]:
    xs = np.sort(p['s'])
    ax.step(xs, np.arange(1,len(xs)+1)/len(xs), where='post', lw=1.8, color=color, ls=ls, label=label)
ax.set(xlim=(.75,1), ylim=(0,1.02), xlabel=r'$s=\sin^2\theta_{\mathrm{red}}$', ylabel='Cumulative proportion')
ax.legend(loc='upper left', frameon=False, fontsize=9)
ax.text(.03,.56,'Mean s\nA: 0.953759\nB: 0.954750\nReference: 0.969198', transform=ax.transAxes,
        fontsize=9, va='top', bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
save(fig,'fig1_ecdf_s_vs_reference')

# The historical markers are pinned to the author-hosted version, not arXiv v2.
historical = {'1': '93.5', '2': '35.6', '4': '9.7', '8': '1.97', '16': '0.35'}
fig, ax = plt.subplots(figsize=(7.2,4.8), layout='constrained')
grid = np.geomspace(1,16,500)
ax.loglog(grid,3/(math.pi*grid),'--',color=gray,lw=1.8,label=r'Hyperbolic reference $3/(\pi Y)$')
for p,color,ls,label in [(A,blue,'-','A: all curves'),(B,orange,':','B: class representatives')]:
    # Add observed jump locations to the plotting grid.
    ys=np.sort(p['y'])
    thresholds=np.unique(np.r_[grid, ys[ys>=1]])
    probabilities=(len(ys)-np.searchsorted(ys,thresholds,side='left'))/len(ys)
    pos=probabilities>0
    ax.step(thresholds[pos],probabilities[pos],where='pre',color=color,lw=1.7,ls=ls,label=label)
wx=np.array([int(v) for v in historical])
wy=np.array([float(v)/100 for v in historical.values()])
ax.plot(wx,wy,'D',color=green,ms=6,label='Watkins, author-hosted version, §6.4')
for x,y in zip(wx,wy):
    ax.annotate(f'{100*y:g}%',(x,y),xytext=(5,4),textcoords='offset points',fontsize=8,color=green)
ax.set(xlim=(.95,19),ylim=(2e-5,1.6),xlabel=r'Threshold $Y$',ylabel=r'Proportion with $y\geq Y$')
ax.set_xticks([1,2,4,8,16],['1','2','4','8','16'])
ax.text(.98,.94,'A and B: zero at Y = 16\n(zero is omitted on the log scale)',transform=ax.transAxes,
        ha='right',va='top',fontsize=8.5,bbox={'facecolor':'white','edgecolor':'none','alpha':.9})
ax.legend(loc='lower left',bbox_to_anchor=(0,.16),frameon=False,fontsize=8)
save(fig,'fig2_y_tail')

fig,axs=plt.subplots(1,2,figsize=(8.6,4.5),sharey=True,layout='constrained')
xx=np.linspace(0,.5,300)
for ax,p,color,title in [(axs[0],A,blue,'A: all curves (23,612)'),(axs[1],B,orange,'B: class representatives (18,886)')]:
    ax.scatter(p['x'],p['y'],s=2.4,color=color,alpha=.3,linewidths=0,rasterized=True)
    ax.plot(xx,np.sqrt(1-xx**2),color=gray,lw=1.3)
    ax.axhline(1,ls='--',lw=.8,color=gray)
    ax.set(xlim=(0,.5),ylim=(.82,8.9),xlabel=r'$x=|b|/h_1$')
    ax.set_title(title,fontsize=9)
axs[0].set_ylabel(r'$y=\sqrt{\det G}/h_1$',labelpad=8)
save(fig,'fig3_fundamental_domain')

after=hashlib.sha256(source.read_bytes()).hexdigest()
if before != after:
    raise RuntimeError('Source changed during plotting')
evidence={
    'source_file':'results/per_curve_shapes.csv','source_sha256':before,
    'source_unchanged':True,'population_A':len(A['s']),'population_B':len(B['s']),
    'float_conversion':'plotting only','new_height_computation':False,
    'watkins_percentages':historical,
    'watkins_source_url':'https://magma.maths.usyd.edu.au/~watkins/papers/heur.pdf',
    'watkins_source_sha256':'0b7e6850026b51c160c38932c9731e4caaab795847429aff17c07efcde11709e',
    'figures':[str(p.relative_to(args.output_root)) for p in sorted(figdir.glob('*'))],
}
(args.output_root/'verification/FIGURE_PROVENANCE.json').write_text(json.dumps(evidence,indent=2)+'\n')
print(json.dumps({'status':'PASS','source_unchanged':True,'figures':len(evidence['figures'])}))
