# -*- coding: utf-8 -*-
"""Cost-effectiveness acceptability curve over the three pre-specified strategies.

The curve in markov_base.py includes the post-hoc combined strategy S3, which
is reported separately in Supplementary material 8. The main text reports the
curve over S0, S1 and S2 only, computed on the fully nested flag count of 113.

Produces, for the manuscript: Figure 4b.

Run from this directory. Expects the cohort file at ../../data/Kohort_v17.xlsx,
as the other scripts in this repository do.
"""
import numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','figure.dpi':600,
                     'axes.spines.top':False,'axes.spines.right':False})
SRC="../04_economics/markov_base.py"
src=open(SRC,encoding='utf-8').read()
ns={}; exec(compile(src[:src.index("# ================================================================= 1. TABAN")],"h","exec"),ns)
BAZ=ns['BAZ']; strateji=ns['strateji']
gamma_par=lambda o,s:((o/s)**2, s**2/o)
def beta_par(o,s):
    v=s**2; a=o*(o*(1-o)/v-1); b=a*(1-o)/o; return max(a,.1),max(b,.1)
rng=np.random.default_rng(2026)
res={'S0':[],'S1':[],'S2':[]}
for _ in range(10000):
    P=dict(BAZ)
    for k,so in [('c_nak_sistemik',.15),('c_adjuvan_sistemik',.15),('c_cerrahi_nak',.20),
                 ('c_cerrahi_once',.20),('c_rt',.15),('c_lrr',.35),('c_met_yillik',.40),
                 ('c_lenfodem_yillik',.40),('c_izlem',.30),('c_karar_destegi',.50)]:
        a,b=gamma_par(BAZ[k],BAZ[k]*so); P[k]=rng.gamma(a,b)
    for k,s in [('u_dfs',.03),('u_kt',.04),('u_ilk_yil',.04),('u_lrr',.05),('u_met',.05)]:
        a,b=beta_par(BAZ[k],s); P[k]=rng.beta(a,b)
    for k,s in [('du_lenfodem',.015),('du_mastektomi',.015),('du_fn',.03)]:
        P[k]=max(rng.normal(BAZ[k],s),0)
    for k,s in [('adjuvan_kt_orani',.06),('p_lenfodem_alnd',.05),('p_lenfodem_slnb',.02),
                ('mkc_nak',.05),('mkc_once',.05),('fn_insidans',.04),('alnd_nak',.06),
                ('alnd_once',.06),('efs10_rcb01',.03),('efs10_rcb23',.04)]:
        a,b=beta_par(BAZ[k],s); P[k]=rng.beta(a,b)
    P['lrr_ek_neoadjuvan']=max(rng.normal(0.055,0.016),0)
    P['s2_n']=int(np.clip(rng.normal(113,8.4),60,188))
    res['S0'].append(strateji(P,0,0)); res['S1'].append(strateji(P,P['s1_n'],0))
    res['S2'].append(strateji(P,P['s2_n'],P['c_karar_destegi']))
df={k:np.array(v) for k,v in res.items()}
lams=np.linspace(0,3_000_000,121); ceac={k:[] for k in ['S0','S1','S2']}
for lam in lams:
    nmb=np.column_stack([lam*df[k][:,1]-df[k][:,0] for k in ['S0','S1','S2']])
    w=nmb.argmax(axis=1)
    for i,k in enumerate(['S0','S1','S2']): ceac[k].append((w==i).mean())
C={'iyi':'#2a9d8f','notr':'#adb5bd','vurgu':'#f4a261'}
fig,ax=plt.subplots(figsize=(5.6,3.6))
ax.plot(lams/1000,ceac['S2'],lw=1.8,color=C['vurgu'],label='S2 prediction rule')
ax.plot(lams/1000,ceac['S1'],lw=1.8,color=C['iyi'],ls='--',label='S1 subtype rule')
ax.plot(lams/1000,ceac['S0'],lw=1.8,color=C['notr'],ls='-.',label='S0 cohort-observed pathway')
ax.axvline(600,ls=':',lw=1.0,color='#6c757d')
ax.set_xlabel('Threshold (thousand TRY per QALY)'); ax.set_ylabel('Probability of being optimal')
ax.set_ylim(0,1); ax.legend(fontsize=7.5,frameon=False)
plt.tight_layout()
out='../../figures/FN_CEAC3.png'
plt.savefig(out,bbox_inches='tight',dpi=600); print('saved',out)
for lam in [0,600_000,1_800_000,3_000_000]:
    i=int(np.argmin(abs(lams-lam)))
    print(f'  {lam:>9,}: S2 {100*ceac["S2"][i]:.1f}%  S1 {100*ceac["S1"][i]:.1f}%  S0 {100*ceac["S0"][i]:.1f}%')
