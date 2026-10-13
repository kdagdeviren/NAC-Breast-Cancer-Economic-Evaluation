# -*- coding: utf-8 -*-
"""Supplementary Figures S3 and S4.

S3, pathological complete response by molecular subtype, is redrawn in English
on the whole cohort (n = 328). The previously archived version was in Turkish,
had overlapping labels and carried obsolete rates (1.9% and 18.6% against the
1.6% and 20.3% the dataset gives). The three records with an equivocal HER2
result and no in situ hybridization confirmation are excluded, as in Table S4.4.

S4, the calibration curve of the base-case model, is redrawn with a decimal
point and the Brier score the recomputation gives (0.166), which is the value
the text reports; the previously archived version printed 0,168.

The base-case out-of-fold probabilities are recomputed with the loop of
../03_validation/bootstrap_dca.py and checked against the published
main-pipeline figures (108 flagged, 98 concordant, 10 discordant, PPV 0.907,
57 flagged in Luminal B with 52 concordant) before anything is drawn.

Run from this directory. Expects the cohort file at ../../data/Kohort_v17.xlsx,
as the other scripts in this repository do.
"""

import os, re, sys, numpy as np, pandas as pd, warnings
warnings.filterwarnings("ignore")
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC  = "../03_validation/bootstrap_dca.py"
DATA = "../../data/Kohort_v17.xlsx"
OUT  = "../../figures"
os.makedirs(OUT, exist_ok=True)

src = open(SRC, encoding="utf-8").read()
src = src.replace("pd.read_excel('../../data/Kohort_v17.xlsx')", f"pd.read_excel(r'{DATA}')")
head = src[:src.index("p = (pi.mean(1) + pm.mean(1)) / 2")] + "p = (pi.mean(1) + pm.mean(1)) / 2\n"
ns = {}
exec(compile(head, "bootstrap_dca_head", "exec"), ns)

d, y, elig, s1, lumB, p = (ns[k] for k in ("d", "y", "elig", "s1", "lumB", "p"))
from sklearn.metrics import roc_auc_score, brier_score_loss

m   = elig & (p >= 0.80)
n   = int(m.sum()); tp = int(y[m].sum())
mb  = lumB & (p >= 0.80)
print("=== base case (main pipeline) reproduction check")
print(f"  AUC whole {roc_auc_score(y, p):.3f}  AUC eligible {roc_auc_score(y[elig], p[elig]):.3f}")
print(f"  Brier (whole cohort) {brier_score_loss(y, p):.3f}")
print(f"  flagged {n} (published 108) | concordant {tp} (98) | discordant {n-tp} (10)")
print(f"  PPV {tp/n:.3f} (0.907)")
print(f"  Luminal B HER2- flagged {int(mb.sum())} (57), concordant {int(y[mb].sum())} (52)")
ok = (n == 108 and tp == 98 and int(mb.sum()) == 57 and int(y[mb].sum()) == 52)
print("  ->", "REPRODUCED" if ok else "DOES NOT MATCH - figures not regenerated")
if not ok:
    sys.exit(1)

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.edgecolor": "#22333b", "axes.labelcolor": "#22333b",
                     "text.color": "#22333b", "xtick.color": "#22333b",
                     "ytick.color": "#22333b"})

# ---------------------------------------------------------------- S3
d0 = pd.read_excel(DATA)          # whole cohort, as in Table 1
LABEL = {"Luminal A": "Luminal A",
         "Luminal B (HER2 Negatif)": "Luminal B\n(HER2−)",
         "Luminal B (HER2 Pozitif)": "Luminal B\n(HER2+)",
         "HER2-Zengin": "HER2-enriched",
         "Triple Negatif": "Triple negative"}
ELIG = {"Luminal A", "Luminal B (HER2 Negatif)"}

rows = []
for key, lab in LABEL.items():
    sub = d0[d0.Molekuler_i == key]
    if not len(sub):
        print(f"  !! subtype absent from the data: {key}"); continue
    rows.append((lab, 100 * (sub.RCB_Kategorize == 0).mean(), len(sub), key in ELIG))
rows.sort(key=lambda r: r[1])
print("\n=== Supplementary Figure S3, whole cohort (n = %d)" % len(d0))
for lab, pct, nn, el in rows:
    print(f"  {lab.replace(chr(10),' '):22s} {pct:5.1f}%  n = {nn:3d}{'   [eligible]' if el else ''}")

fig, ax = plt.subplots(figsize=(6.6, 3.4))
xs = np.arange(len(rows))
cols = ["#2a9d8f" if r[3] else "#b8c0c8" for r in rows]
ax.bar(xs, [r[1] for r in rows], color=cols, width=0.62)
for x, r in zip(xs, rows):
    ax.text(x, r[1] + 1.6, f"{r[1]:.1f}%", ha="center", va="bottom",
            fontsize=9, fontweight="bold")
ax.set_xticks(xs)
ax.set_xticklabels([f"{r[0]}\nn = {r[2]}" for r in rows], fontsize=8.5)
ax.set_ylabel("Pathological complete response rate (%)")
ax.set_ylim(0, max(r[1] for r in rows) * 1.28)
ax.tick_params(axis="x", length=0)
ymax = ax.get_ylim()[1]
el_x = [x for x, r in zip(xs, rows) if r[3]]
ax.annotate("", xy=(min(el_x) - 0.38, ymax * 0.93), xytext=(max(el_x) + 0.38, ymax * 0.93),
            arrowprops=dict(arrowstyle="-", color="#2a9d8f", lw=1.6))
ax.text((min(el_x) + max(el_x)) / 2, ymax * 0.95,
        "biology in which the rule operates", ha="center", va="bottom",
        fontsize=8.5, style="italic", color="#2a9d8f")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Supplementary_Figure_S3.png"), dpi=600, bbox_inches="tight")
plt.close()

# ---------------------------------------------------------------- S4
def bins(pp, yy, k=8):
    g = pd.DataFrame({"p": pp, "y": yy, "q": pd.qcut(pp, k, duplicates="drop")}) \
          .groupby("q", observed=True)
    return g.p.mean().values, g.y.mean().values

br = brier_score_loss(y, p)
mp, mo = bins(p, y)
print(f"\n=== Supplementary Figure S4: Brier = {br:.3f} (published 0.168), {len(mp)} bins")
fig, ax = plt.subplots(figsize=(4.2, 3.9))
ax.plot([0, 1], [0, 1], ls="--", color="#adb5bd", lw=1.2)
ax.plot(mp, mo, marker="o", color="#3d5a80", lw=1.8, ms=6)
ax.text(0.04, 0.95, f"Brier score = {br:.3f}", transform=ax.transAxes,
        ha="left", va="top", fontsize=9.5)
ax.set_xlabel("Predicted probability")
ax.set_ylabel("Observed proportion")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.set_xticks(np.arange(0, 1.01, 0.2)); ax.set_yticks(np.arange(0, 1.01, 0.2))
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Supplementary_Figure_S4.png"), dpi=600, bbox_inches="tight")
plt.close()
print("\nwritten to", OUT)
