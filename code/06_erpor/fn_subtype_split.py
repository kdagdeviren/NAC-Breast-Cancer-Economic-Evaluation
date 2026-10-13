# -*- coding: utf-8 -*-
"""Fully nested re-run: the subtype breakdown and the outer-partition range.

Runs part B of ../03_validation/tam_icice.py three times and keeps every repeat
instead of overwriting them, which is what that script does. The reported fully
nested figures are those of the third partition (random_state 702):
  AUC 0.761 whole cohort, 0.694 eligible, scaled Brier 0.051,
  113 flagged, 97 concordant, 16 discordant, PPV 0.858, +42 over the subtype rule.
Produces, for the manuscript:
  Table 2, fully nested column
  Table 3, fully nested column, including the per-subtype split
  Section 3.2, the range across outer partitions

Run from this directory. Expects the cohort file at ../../data/Kohort_v17.xlsx,
as the other scripts in this repository do.
"""

import sys, numpy as np, pandas as pd, collections, warnings
warnings.filterwarnings("ignore")

SRC  = "../03_validation/tam_icice.py"
DATA = "../../data/Kohort_v17.xlsx"

src = open(SRC, encoding="utf-8").read()
src = src.replace("pd.read_excel('../../data/Kohort_v17.xlsx')", f"pd.read_excel(r'{DATA}')")
head = src[:src.index("S=list(pd.read_csv(")]

ns = {}
exec(compile(head, "tam_icice_head", "exec"), ns)

for k in ("X0", "y", "grp", "elig", "lumA", "VAR", "BLOK", "imp_modal",
          "fit_pred", "fit_pred_ic", "d"):
    globals()[k] = ns[k]

from sklearn.model_selection import StratifiedGroupKFold
from sklearn.metrics import roc_auc_score, brier_score_loss

KOMB = [('P',), ('P', 'O'), ('P', 'R'), ('P', 'O', 'R')]

# The author's loop overwrites oof on each repeat, so the reported figures come
# from the last repeat (r = 2). Keep every repeat so that can be verified.
oof_r = {}
for r in range(3):
    oof = np.full(len(y), np.nan)
    secilen = []
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=700 + r).split(X0, y, groups=grp):
        Xi = imp_modal(X0, tr)
        en, ea = None, -1
        for k in KOMB:
            vs = [v for v in sum([BLOK[b] for b in k], []) if v in VAR]
            ip = np.full(len(tr), np.nan)
            Xtr = X0.iloc[tr].reset_index(drop=True); ytr = y[tr]; gtr = grp[tr]
            for itr, ite in StratifiedGroupKFold(3, shuffle=True, random_state=11).split(Xtr, ytr, groups=gtr):
                Xii = imp_modal(Xtr, itr)
                ip[ite] = fit_pred_ic(Xii, vs, itr, ite, ytr, r)
            try:    a = roc_auc_score(y[tr], ip)
            except Exception: a = 0
            if a > ea: ea, en = a, k
        secilen.append(en)
        vs = [v for v in sum([BLOK[b] for b in en], []) if v in VAR]
        oof[te] = fit_pred(Xi, vs, tr, te, r, dummy=True)
    oof_r[r] = oof
    print(f"  repeat {r+1}/3 done | blocks chosen: {collections.Counter(secilen).most_common()}",
          flush=True)

np.save("../../outputs/oof_tam_icice_repeats.npy",
        np.column_stack([oof_r[r] for r in range(3)]))

lumB = (d.Molekuler_i.values == 'Luminal B (HER2 Negatif)')

def report(K, tag):
    m   = elig & (K >= 0.80)
    n   = int(m.sum()); tp = int(y[m].sum())
    print(f"\n--- {tag}")
    print(f"  AUC whole {roc_auc_score(y, K):.3f} | AUC eligible {roc_auc_score(y[elig], K[elig]):.3f}"
          f" | scaled Brier {1 - brier_score_loss(y[elig], K[elig]) / brier_score_loss(y[elig], np.full(elig.sum(), y[elig].mean())):.3f}")
    print(f"  flagged {n} | concordant {tp} | discordant {n - tp} | PPV {tp/n if n else float('nan'):.3f}")
    fa = int((m & lumA).sum()); fb = int((m & lumB).sum())
    print(f"  Luminal A flagged     {fa} of {int(lumA.sum())}  ({int(y[m & lumA].sum())} concordant)")
    print(f"  Luminal B HER2- flagged {fb} of {int(lumB.sum())}  ({int(y[m & lumB].sum())} concordant)")
    print(f"  Luminal A flagged by S1 but NOT S2: {int(lumA.sum()) - fa}"
          f"  (of which RCB-II/III: {int(y[lumA & ~m].sum())})")
    print(f"  additional concordant vs S1: +{tp - int(y[lumA].sum())}")

for r in range(3):
    report(oof_r[r], f"repeat {r} alone")
report(np.nanmean(np.column_stack([oof_r[r] for r in range(3)]), axis=1), "mean of 3 repeats")
