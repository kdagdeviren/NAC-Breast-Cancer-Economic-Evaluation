# -*- coding: utf-8 -*-
"""Harrell iyimserlik düzeltmesi — tam boru hattı her önyükleme örneğinde yeniden kurulur."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2

P=['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','TIL_i','HistolojikG_i']
O=['Metastaz_Yeri_i','NAC_Once_Evre_i']
R=['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Kitle_Sekli_i','Kitle_Konturu_i',
   'Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i','Kalsifikasyon_Dagilimi_i','Asimetri_i',
   'Multifokalite_Durumu_59','Cilt_Cekintisi_i','Meme_Basi_Retraksiyonu_i']
VARS=P+O+R
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values; grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
Xc=d[VARS].apply(lambda s:s.astype('category').cat.codes).values
du={v:pd.get_dummies(d[v].astype(str),drop_first=True).values.astype(float) for v in VARS}

def llr_p(X,yy):
    if X.shape[1]==0 or len(np.unique(yy))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(X,yy)
    q=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
    ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),X.shape[1])

def boru(idx_fit, idx_tahmin, seed=0):
    """Tam boru hattı: değişken seçimi + LR + HGB + ortalama"""
    sel=[v for v in VARS if llr_p(du[v][idx_fit], y[idx_fit])<0.10]
    if not sel: sel=VARS
    Xl=np.hstack([du[v] for v in sel])
    lr=LogisticRegression(C=0.05,max_iter=5000).fit(Xl[idx_fit],y[idx_fit])
    p1=lr.predict_proba(Xl[idx_tahmin])[:,1]
    hgb=HistGradientBoostingClassifier(max_iter=300,max_leaf_nodes=15,learning_rate=0.05,
        min_samples_leaf=15,l2_regularization=1.0,random_state=seed).fit(Xc[idx_fit],y[idx_fit])
    p2=hgb.predict_proba(Xc[idx_tahmin])[:,1]
    return (p1+p2)/2

tum=np.arange(len(y))
p_app=boru(tum,tum)
app_auc=roc_auc_score(y,p_app); app_auc_e=roc_auc_score(y[elig],p_app[elig])
app_br=brier_score_loss(y,p_app); app_br_e=brier_score_loss(y[elig],p_app[elig])
print(f'Görünürde (apparent): AUC tüm {app_auc:.3f} · uygun {app_auc_e:.3f} · Brier {app_br:.3f}/{app_br_e:.3f}')

B=200; rng=np.random.default_rng(2026)
hastalar=np.unique(grp); opt=[]
for b in range(B):
    hs=rng.choice(hastalar,len(hastalar),replace=True)
    idx=np.concatenate([np.where(grp==h)[0] for h in hs])
    if len(np.unique(y[idx]))<2: continue
    pb_b=boru(idx,idx,seed=b); pb_o=boru(idx,tum,seed=b)
    eb=elig[idx]
    try:
        o=dict(auc=roc_auc_score(y[idx],pb_b)-roc_auc_score(y,pb_o),
               auc_e=roc_auc_score(y[idx][eb],pb_b[eb])-roc_auc_score(y[elig],pb_o[elig]),
               br=brier_score_loss(y[idx],pb_b)-brier_score_loss(y,pb_o),
               br_e=brier_score_loss(y[idx][eb],pb_b[eb])-brier_score_loss(y[elig],pb_o[elig]))
        opt.append(o)
    except Exception: pass
    if (b+1)%50==0: print(f'  {b+1}/{B}')
O_=pd.DataFrame(opt).mean()
print(f'\nOrtalama iyimserlik ({len(opt)} örnek):')
print(f'  AUC tüm kohort   {O_.auc:+.4f}  → düzeltilmiş {app_auc-O_.auc:.3f}')
print(f'  AUC uygun pop.   {O_.auc_e:+.4f}  → düzeltilmiş {app_auc_e-O_.auc_e:.3f}')
print(f'  Brier tüm        {O_.br:+.4f}  → düzeltilmiş {app_br-O_.br:.3f}')
print(f'  Brier uygun      {O_.br_e:+.4f}  → düzeltilmiş {app_br_e-O_.br_e:.3f}')
pd.DataFrame([dict(gorunurde_auc=app_auc, iyimserlik_auc=O_.auc, duzeltilmis_auc=app_auc-O_.auc,
   gorunurde_auc_uygun=app_auc_e, iyimserlik_auc_uygun=O_.auc_e, duzeltilmis_auc_uygun=app_auc_e-O_.auc_e,
   B=len(opt))]).to_csv('../../outputs/optimizm_duzeltmesi.csv',index=False)
