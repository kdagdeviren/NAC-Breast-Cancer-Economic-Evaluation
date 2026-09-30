# -*- coding: utf-8 -*-
"""Sensitivity of the economic result to the cost of a wrong flag (Section 3.18, Table 13, Fig. 9).

Net saving per patient vs S1 = base saving − (additional false flags × H) / n
where n = 320 (base-case cohort), and the two validation layers differ in both terms.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

N = 320                      # base-case cohort
LAM = 600_000                # illustrative willingness-to-pay, TRY per QALY
PPP = 15.351                 # TRY per international dollar (OECD 2025)
LAYERS = {'Base case':        dict(saving=3131, extra_false=4),
          'Fully nested':     dict(saving=3480, extra_false=10)}

def net_saving(H, saving, extra_false, n=N):
    return saving - extra_false * H / n

def break_even(saving, extra_false, n=N):
    return saving * n / extra_false

rows = []
for H in [0, 25_000, 50_000, 75_000, 100_000, 111_360, 150_000, 200_000, 250_480]:
    rows.append({'H_TRY': H,
                 **{k: round(net_saving(H, **v)) for k, v in LAYERS.items()}})
tab = pd.DataFrame(rows)
print(tab.to_string(index=False))
tab.to_csv('../outputs/false_flag_table.csv', index=False)

for k, v in LAYERS.items():
    be = break_even(**v)
    print(f'{k}: break-even H = {be:,.0f} TRY ({be/PPP:,.0f} int$; {be/LAM:.2f} QALY)')

H = np.linspace(0, 300_000, 601)
fig, ax = plt.subplots(figsize=(5.0, 3.4))
for (k, v), col, ls in zip(LAYERS.items(), ['#3d5a80', '#e76f51'], ['-', '--']):
    ax.plot(H/1000, net_saving(H, **v), lw=2, color=col, ls=ls,
            label=f"{k} ({v['extra_false']} additional false flags)")
    ax.plot([break_even(**v)/1000], [0], 'o', color=col, ms=6, zorder=5)
ax.axhline(0, color='#22333b', lw=1)
ax.set_xlabel('Harm per additional false flag, $H$ (thousand TRY)')
ax.set_ylabel('Net saving per patient vs S1 (TRY)')
sec = ax.secondary_xaxis('top', functions=(lambda x: x*1000/LAM, lambda q: q*LAM/1000))
sec.set_xlabel('Equivalent QALY loss per false flag')
ax.legend(fontsize=7.5, frameon=False); plt.tight_layout()
plt.savefig('../figures/Fig9.png', dpi=600, bbox_inches='tight')
