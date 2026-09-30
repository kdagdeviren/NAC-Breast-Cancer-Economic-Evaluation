# -*- coding: utf-8 -*-
"""Hakem talebi 1: imputasyon/radyoloji bağımlılığı — alternatif model spesifikasyonları."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2

P=['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','TIL_i','HistolojikG_i']
O=['Metastaz_Yeri_i','NAC_Once_Evre_i']
R=['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Kitle_Sekli_i','Kitle_Konturu_i',
   'Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i','Kalsifikasyon_Dagilimi_i','Asimetri_i',
   'Multifokalite_Durumu_59','Cilt_Cekintisi_i','Meme_Basi_Retraksiyonu_i']
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values; grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
lumA=(d.Molekuler_i.values=='Luminal A')

# ham (imputasyonsuz) karsiliklarini bul → eksiklik oranlari
def ham_ad(v):
    for c in [v.replace('_i',''), v.replace('_i','_ham'), v]:
        if c in d0.columns: return c
    return None
eks={}
for v in P+O+R:
    h=ham_ad(v)
    eks[v]= d[h].isna().mean() if h and h in d.columns else np.nan
E=pd.Series(eks)
dusuk=[v for v in P+O+R if not np.isnan(E.get(v,np.nan)) and E[v]<=0.20]
print('Eksiklik oranları (ham kolonlardan):')
for v in P+O+R:
    print(f'  {v:32s} {E[v]*100 if not np.isnan(E[v]) else float("nan"):5.1f}%')

def calistir(VARS, ad):
    Xc=d[VARS].apply(lambda s:s.astype('category').cat.codes).values
    du={v:pd.get_dummies(d[v].astype(str),drop_first=True).values.astype(float) for v in VARS}
    def llr_p(X,yy):
        if X.shape[1]==0 or len(np.unique(yy))<2: return 1.0
        m=LogisticRegression(C=1e6,max_iter=2000).fit(X,yy)
        q=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
        ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
        return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),X.shape[1])
    oof=np.zeros(len(y)); k=0
    for r in range(5):
        pr=np.full(len(y),np.nan)
        for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=300+r).split(Xc,y,groups=grp):
            sel=[v for v in VARS if llr_p(du[v][tr],y[tr])<0.10] or VARS
            Xl=np.hstack([du[v] for v in sel])
            lr=LogisticRegression(C=0.05,max_iter=5000).fit(Xl[tr],y[tr])
            hgb=HistGradientBoostingClassifier(max_iter=300,max_leaf_nodes=15,learning_rate=0.05,
                min_samples_leaf=15,l2_regularization=1.0,random_state=r).fit(Xc[tr],y[tr])
            pr[te]=(lr.predict_proba(Xl[te])[:,1]+hgb.predict_proba(Xc[te])[:,1])/2
        oof+=pr; k+=1
    oof/=k
    m=elig&(oof>=0.80); tp=int(y[m].sum()); n=int(m.sum())
    tp1=int(y[lumA].sum())
    return dict(model=ad, k=len(VARS),
        AUC_tum=roc_auc_score(y,oof), AUC_uygun=roc_auc_score(y[elig],oof[elig]),
        Brier_uygun=brier_score_loss(y[elig],oof[elig]),
        isaretlenen=n, uyumlu=tp, uyumsuz=n-tp, PPV=tp/n if n else np.nan,
        ek_uyumlu=tp-tp1, ek_uyumsuz=(n-tp)-(int(lumA.sum())-tp1))

SPEC=[(P+O+R,'Tam model (23 değişken)'),
      (P+O,'Radyoloji bloğu çıkarıldı (10 değişken)'),
      (P,'Yalnızca patoloji (8 değişken)'),
      ([v for v in P+O+R if v in dusuk] or P+O, f'Yalnızca düşük eksiklikli değişkenler')]
sat=[calistir(v,a) for v,a in SPEC]
D=pd.DataFrame(sat)
print('\n═══ ALTERNATİF MODEL SPESİFİKASYONLARI\n')
print(D.round(3).to_string(index=False))
D.to_csv('../../outputs/alternatif_modeller.csv',index=False)
