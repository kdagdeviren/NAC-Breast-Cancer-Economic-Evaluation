# -*- coding: utf-8 -*-
"""Fig. 4 (calibration, eligible population) and Fig. 5 (decision curve, eligible population).

Requires oof_final.npy (main pipeline, from bootstrap_dca.py) and
oof_tam_icice.npy (fully nested, from tam_icice.py).
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

src = open('tam_icice.py', encoding='utf-8').read()
exec(src[:src.index('S=list(pd.read_csv')])
base, nested = np.load('oof_final.npy'), np.load('oof_tam_icice.npy')
ye, pe = y[elig], base[elig]

def bins(p, yy, k=8):
    g = pd.DataFrame({'p': p, 'y': yy, 'q': pd.qcut(p, k, duplicates='drop')}).groupby('q', observed=True)
    return g.p.mean().values, g.y.mean().values

fig, ax = plt.subplots(figsize=(4.4, 4.0))
ax.plot([0, 1], [0, 1], ls=':', color='#6c757d', lw=1, label='Perfect calibration')
for p, lab, col, mk in [(pe, 'Main pipeline', '#3d5a80', 'o'), (nested[elig], 'Fully nested', '#e76f51', 's')]:
    mp, mo = bins(p, ye); ax.plot(mp, mo, marker=mk, color=col, lw=1.6, ms=5, label=lab)
ax.set_xlabel('Predicted probability'); ax.set_ylabel('Observed proportion with RCB-II/III')
ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.legend(fontsize=7.5, frameon=False, loc='upper left')
plt.tight_layout(); plt.savefig('../figures/Fig4.png', dpi=600, bbox_inches='tight'); plt.close()

PT = np.linspace(0.30, 0.95, 66); n = len(ye); s1 = (d.Molekuler_i.values[elig] == 'Luminal A')
nb = lambda flag: np.array([((flag & (ye == 1)).sum()/n - (flag & (ye == 0)).sum()/n * (pt/(1-pt))) for pt in PT])
rule = np.array([(((pe >= pt) & (ye == 1)).sum()/n - ((pe >= pt) & (ye == 0)).sum()/n * (pt/(1-pt))) for pt in PT])
fig, ax = plt.subplots(figsize=(5.0, 3.4))
ax.plot(PT, rule, lw=2, color='#3d5a80', label='Prediction rule (S2)')
ax.plot(PT, nb(s1), lw=1.8, color='#2a9d8f', ls='--', label='Subtype rule (S1)')
ax.plot(PT, nb(np.ones(n, bool)), lw=1.4, color='#adb5bd', ls='-.', label='Flag all eligible')
ax.axhline(0, color='#22333b', lw=1, label='Flag none'); ax.axvline(0.80, ls=':', lw=1, color='#6c757d')
ax.set_xlabel('Threshold probability'); ax.set_ylabel('Net benefit')
ax.legend(fontsize=7.5, frameon=False, loc='lower left'); plt.tight_layout()
plt.savefig('../figures/Fig5.png', dpi=600, bbox_inches='tight')
