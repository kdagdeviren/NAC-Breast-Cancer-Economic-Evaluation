# -*- coding: utf-8 -*-
"""Economic model re-run on the fully nested flag count of 113.

Loads the parameters and model functions of ../04_economics/markov_base.py
without executing its reporting sections, then recomputes the base case, the
probabilistic analysis (10,000 draws, seed 2026), the three-strategy
acceptability curve, the expected value of perfect information and the
adherence analysis.

The adherence analysis has no script of its own in this repository. It is
derived here from the cost structure of markov_base.py:

    S1 - S2(a) per patient = (a * F - 61) * dC / 320 - c_karar_destegi

where F is the number flagged, dC the saving per patient switched to a
surgery-first pathway implied by the model, and the decision-support cost is
charged to all 320 patients, exactly as strateji() does. At 108 flagged this
reproduces the previously published break-even of 58.5% and the advantage over
S0 at 50% adherence of 3,620 TRY.

Produces, for the manuscript: Table 4; the CDK4/6 and adherence rows of
Table 5; the probabilistic figures in Section 3.5.

Run from this directory. Expects the cohort file at ../../data/Kohort_v17.xlsx,
as the other scripts in this repository do.
"""

import sys, numpy as np, pandas as pd

SRC = "../04_economics/markov_base.py"
src = open(SRC, encoding="utf-8").read()
head = src[:src.index("# ================================================================= 1. TABAN")]
ns = {}
exec(compile(head, "markov_base_head", "exec"), ns)

BAZ        = ns["BAZ"]
markov     = ns["markov"]
kisa_donem = ns["kisa_donem"]
strateji   = ns["strateji"]
gamma_par  = lambda ort, sd: ((ort/sd)**2, sd**2/ort)
def beta_par(ort, sd):
    v = sd**2; a = ort*(ort*(1-ort)/v - 1); b = a*(1-ort)/ort
    return max(a, 0.1), max(b, 0.1)

N_COH = BAZ["kohort"]

def run(P):
    return {"S0": strateji(P, 0, 0),
            "S1": strateji(P, P["s1_n"], 0),
            "S2": strateji(P, P["s2_n"], P["c_karar_destegi"])}

print("="*78)
print("1. BASE CASE, per patient")
print("="*78)
base = {}
for label, n2 in [("main pipeline (108)", 108), ("fully nested (113)", 113)]:
    P = dict(BAZ); P["s2_n"] = n2
    o = run(P)
    dm = (o["S2"][0]-o["S1"][0])/N_COH
    dq = (o["S2"][1]-o["S1"][1])/N_COH
    dm0 = (o["S2"][0]-o["S0"][0])/N_COH
    base[n2] = o
    print(f"{label}")
    for k in ["S0","S1","S2"]:
        print(f"   {k}: cost {o[k][0]/N_COH:12,.0f} TRY   QALY {o[k][1]/N_COH:.4f}")
    print(f"   S2-S1: {dm:+,.1f} TRY   {dq:+.5f} QALY")
    print(f"   S2-S0: {dm0:+,.1f} TRY")
    print(f"   share of S2-S0 obtainable from S1 alone: "
          f"{100*(o['S0'][0]-o['S1'][0])/(o['S0'][0]-o['S2'][0]):.1f}%\n")

# cost difference per switched patient, implied by the model
dC = (base[108]["S0"][0]-base[108]["S1"][0])/BAZ["s1_n"]
print(f"implied saving per patient switched to surgery-first: {dC:,.1f} TRY\n")

print("="*78)
print("2. PSA, 10,000 draws, seed 2026, fully nested (113)")
print("="*78)
def psa(n2_mean, n_iter=10000, seed=2026):
    rng = np.random.default_rng(seed)
    res = {"S0": [], "S1": [], "S2": []}
    for _ in range(n_iter):
        P = dict(BAZ)
        for k, sd_o in [("c_nak_sistemik",.15),("c_adjuvan_sistemik",.15),
                        ("c_cerrahi_nak",.20),("c_cerrahi_once",.20),("c_rt",.15),
                        ("c_lrr",.35),("c_met_yillik",.40),("c_lenfodem_yillik",.40),
                        ("c_izlem",.30),("c_karar_destegi",.50)]:
            a,b = gamma_par(BAZ[k], BAZ[k]*sd_o); P[k] = rng.gamma(a,b)
        for k, sd in [("u_dfs",.03),("u_kt",.04),("u_ilk_yil",.04),("u_lrr",.05),("u_met",.05)]:
            a,b = beta_par(BAZ[k], sd); P[k] = rng.beta(a,b)
        for k, sd in [("du_lenfodem",.015),("du_mastektomi",.015),("du_fn",.03)]:
            P[k] = max(rng.normal(BAZ[k], sd), 0)
        for k, sd in [("adjuvan_kt_orani",.06),("p_lenfodem_alnd",.05),("p_lenfodem_slnb",.02),
                      ("mkc_nak",.05),("mkc_once",.05),("fn_insidans",.04),
                      ("alnd_nak",.06),("alnd_once",.06),
                      ("efs10_rcb01",.03),("efs10_rcb23",.04)]:
            a,b = beta_par(BAZ[k], sd); P[k] = rng.beta(a,b)
        P["lrr_ek_neoadjuvan"] = max(rng.normal(0.055, 0.016), 0)
        P["s2_n"] = int(np.clip(rng.normal(n2_mean, 8.4), 60, 188))
        o = run(P)
        for k in res: res[k].append(o[k])
    return {k: np.array(v) for k, v in res.items()}

df = psa(113)
dmv = (df["S2"][:,0]-df["S1"][:,0])/N_COH
dqv = (df["S2"][:,1]-df["S1"][:,1])/N_COH
print(f"P(S2 less costly than S1)  = {100*(dmv<0).mean():.1f}%")
print(f"P(S2 more QALYs than S1)   = {100*(dqv>0).mean():.1f}%")
print(f"cost difference per patient: mean {dmv.mean():+,.0f}  "
      f"95% CI {np.percentile(dmv,2.5):+,.0f} to {np.percentile(dmv,97.5):+,.0f} TRY")
print(f"QALY difference per patient: mean {dqv.mean():+.5f}  "
      f"95% CI {np.percentile(dqv,2.5):+.5f} to {np.percentile(dqv,97.5):+.5f}")

print("\n" + "="*78)
print("3. CEAC over the THREE pre-specified strategies (S0, S1, S2)")
print("="*78)
lams = np.linspace(0, 3_000_000, 121)
ceac = {k: [] for k in ["S0","S1","S2"]}
for lam in lams:
    nmb = np.column_stack([lam*df[k][:,1]-df[k][:,0] for k in ["S0","S1","S2"]])
    win = nmb.argmax(axis=1)
    for i,k in enumerate(["S0","S1","S2"]): ceac[k].append((win==i).mean())
for lam in [0, 300_000, 600_000, 900_000, 1_800_000, 3_000_000]:
    i = int(np.argmin(abs(lams-lam)))
    print(f"  lambda={lam:>9,}  S2 {100*ceac['S2'][i]:5.1f}%   "
          f"S1 {100*ceac['S1'][i]:5.1f}%   S0 {100*ceac['S0'][i]:5.1f}%")
for k in ["S0","S1","S2"]:
    a = np.array(ceac[k]); print(f"  {k} range across lambda: {100*a.min():.1f}% to {100*a.max():.1f}%")

print("\n" + "="*78)
print("4. EVPI, three strategies")
print("="*78)
for lam in [0, 600_000, 1_800_000]:
    nmb = np.column_stack([lam*df[k][:,1]-df[k][:,0] for k in ["S0","S1","S2"]])
    evpi = (nmb.max(axis=1).mean() - nmb.mean(axis=0).max())/N_COH
    print(f"  lambda={lam:>9,}  EVPI = {evpi:,.0f} TRY per patient")

print("\n" + "="*78)
print("5. ADHERENCE, derived from the model's own cost structure")
print("="*78)
print("  S1 - S2(a) per patient = (a*F - 61) * dC / 320 - 155")
print(f"  with dC = {dC:,.1f} TRY and the decision-support cost charged to all 320\n")
def adh(a, F):
    return (a*F - BAZ["s1_n"]) * dC / N_COH - BAZ["c_karar_destegi"]
for F, lab in [(108,"main pipeline"), (113,"fully nested")]:
    be = (BAZ["c_karar_destegi"]*N_COH/dC + BAZ["s1_n"])/F
    print(f"  {lab} (F={F}): break-even adherence = {100*be:.1f}%")
    for a in [1.0, 0.8, 0.7, 0.6, 0.5]:
        print(f"      a={a:.0%}:  S1-S2 = {adh(a,F):+,.0f} TRY" , end="")
        o = run({**BAZ, "s2_n": F})
        s0 = ((BAZ['s1_n']*0 + 0) )  # placeholder
        # advantage over S0 at this adherence
        advS0 = a*F*dC/N_COH - BAZ["c_karar_destegi"]
        print(f"   |  S0-S2 = {advS0:+,.0f} TRY")
    print()
