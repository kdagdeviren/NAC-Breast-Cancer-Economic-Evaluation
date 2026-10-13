# -*- coding: utf-8 -*-
"""Supplementary Figures S1 and S2.

S1, the strategy diagram, is redrawn from ../05_figures/sekil_strateji.py with
the fully nested classification counts (113 / 97 / 16, PPV 0.858) in place of
the main-pipeline counts, and with the bootstrap intervals labeled as
main-pipeline intervals, which is what they are.

S2, the cost-effectiveness plane against the subtype rule, is drawn in English
at 600 dpi from the fully nested probabilistic draws. Its reference point
reproduces -3,480 TRY and +0.00140 QALY.

Run from this directory. Expects the cohort file at ../../data/Kohort_v17.xlsx,
as the other scripts in this repository do.
"""

import os, re, sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "../../figures"
os.makedirs(OUT, exist_ok=True)

# ------------------------------------------------------------------ S1
SRC = "../05_figures/sekil_strateji.py"
src = open(SRC, encoding="utf-8").read()

EDITS = [
 ("[('Flagged','108'),('Concordant / discordant','98 / 10'),('Predictive value','0.907')]",
  "[('Flagged','113'),('Concordant / discordant','97 / 16'),('Predictive value','0.858')]"),
 # the published intervals are the main-pipeline bootstrap intervals, so they
 # cannot be shown beside a fully nested point estimate; give both lines instead
 ("sat=[('Additional concordant','+43','+28 – +58',C['iyi']),\n"
  "     ('Additional discordant','+4','−1 – +9',C['kotu'])]",
  "sat=[('Additional concordant','+42','+43 (+28 – +58)',C['iyi']),\n"
  "     ('Additional discordant','+10','+4 (−1 – +9)',C['kotu'])]"),
 ("yaz(6.6,2.62,'Observed',8.4,'bold',c='#4a5568')",
  "yaz(6.6,2.62,'Fully nested',8.4,'bold',c='#4a5568')"),
 ("yaz(9.2,2.62,'95% confidence interval',8.4,'bold',c='#4a5568')",
  "yaz(9.2,2.62,'Main pipeline (95% CI)',8.4,'bold',c='#4a5568')"),
 ("yaz(6,0.85,'The gain excludes zero; the additional discordant classifications do not',\n"
  "    8.3,c='#4a5568')",
  "yaz(6,0.85,'Concordant gain is unchanged; discordant classifications more than double',\n"
  "    8.3,c='#4a5568')"),
]
for o, n in EDITS:
    assert src.count(o) == 1, f"S1 edit not found: {o[:50]}"
    src = src.replace(o, n)

# drop the script's own save chain and write one PNG where we want it
S1OUT = os.path.join(OUT, "Supplementary_Figure_S1.png")
src = re.sub(r"plt\.savefig\([^\)]*\)", "pass", src)
src += f"\nplt.savefig(r'{S1OUT}', dpi=600, bbox_inches='tight', facecolor='white')\n"
ns = {"__name__": "__main__"}
exec(compile(src, "sekil_strateji_patched", "exec"), ns)
plt.close("all")
assert os.path.exists(S1OUT), "S1 was not written"
print("  S1 regenerated with the fully nested counts")

# ------------------------------------------------------------------ S2
MB = "../04_economics/markov_base.py"
msrc = open(MB, encoding="utf-8").read()
mns = {}
exec(compile(msrc[:msrc.index("# ================================================================= 1. TABAN")],
             "markov_base_head", "exec"), mns)
BAZ, strateji = mns["BAZ"], mns["strateji"]
N = BAZ["kohort"]

gamma_par = lambda m, s: ((m / s) ** 2, s ** 2 / m)
def beta_par(m, s):
    v = s ** 2; a = m * (m * (1 - m) / v - 1)
    return max(a, 0.1), max(a * (1 - m) / m, 0.1)

rng = np.random.default_rng(2026)
d1 = np.empty((10000, 2)); d2 = np.empty((10000, 2))
for i in range(10000):
    P = dict(BAZ)
    for k, so in [("c_nak_sistemik", .15), ("c_adjuvan_sistemik", .15), ("c_cerrahi_nak", .20),
                  ("c_cerrahi_once", .20), ("c_rt", .15), ("c_lrr", .35), ("c_met_yillik", .40),
                  ("c_lenfodem_yillik", .40), ("c_izlem", .30), ("c_karar_destegi", .50)]:
        a, b = gamma_par(BAZ[k], BAZ[k] * so); P[k] = rng.gamma(a, b)
    for k, s in [("u_dfs", .03), ("u_kt", .04), ("u_ilk_yil", .04), ("u_lrr", .05), ("u_met", .05)]:
        a, b = beta_par(BAZ[k], s); P[k] = rng.beta(a, b)
    for k, s in [("du_lenfodem", .015), ("du_mastektomi", .015), ("du_fn", .03)]:
        P[k] = max(rng.normal(BAZ[k], s), 0)
    for k, s in [("adjuvan_kt_orani", .06), ("p_lenfodem_alnd", .05), ("p_lenfodem_slnb", .02),
                 ("mkc_nak", .05), ("mkc_once", .05), ("fn_insidans", .04), ("alnd_nak", .06),
                 ("alnd_once", .06), ("efs10_rcb01", .03), ("efs10_rcb23", .04)]:
        a, b = beta_par(BAZ[k], s); P[k] = rng.beta(a, b)
    P["lrr_ek_neoadjuvan"] = max(rng.normal(0.055, 0.016), 0)
    P["s2_n"] = int(np.clip(rng.normal(113, 8.4), 60, 188))
    d1[i] = strateji(P, P["s1_n"], 0)
    d2[i] = strateji(P, P["s2_n"], P["c_karar_destegi"])

dq = (d2[:, 1] - d1[:, 1]) / N
dc = (d2[:, 0] - d1[:, 0]) / N
b = {**BAZ, "s2_n": 113}
bq = (strateji(b, 113, b["c_karar_destegi"])[1] - strateji(b, b["s1_n"], 0)[1]) / N
bc = (strateji(b, 113, b["c_karar_destegi"])[0] - strateji(b, b["s1_n"], 0)[0]) / N
print(f"  S2 plane: base case {bc:+,.0f} TRY, {bq:+.5f} QALY | "
      f"SW quadrant (less costly, more QALYs) {100*((dc<0)&(dq>0)).mean():.1f}%")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(5.0, 4.0))
ax.scatter(dq, dc, s=3, alpha=0.18, color="#3d5a80", linewidths=0, rasterized=False)
lim = np.array([dq.min() * 1.05, dq.max() * 1.05])
ax.plot(lim, 600_000 * lim, ls="--", color="#e76f51", lw=1.5,
        label="600,000 TRY per QALY")
ax.axhline(0, color="#6c757d", lw=0.9); ax.axvline(0, color="#6c757d", lw=0.9)
ax.scatter([bq], [bc], s=70, color="#f4a261", edgecolor="#22333b", zorder=5,
           label="Reference analysis")
ax.set_xlabel("Incremental QALYs per patient, prediction rule versus subtype rule")
ax.set_ylabel("Incremental cost per patient (TRY)")
ax.set_xlim(*lim)
ax.legend(fontsize=8, frameon=False, loc="upper left")
ax.ticklabel_format(axis="y", style="plain")
plt.tight_layout()
plt.savefig(os.path.join(OUT, "Supplementary_Figure_S2.png"), dpi=600, bbox_inches="tight")
plt.close()
print("  S2 regenerated in English at 600 dpi")
