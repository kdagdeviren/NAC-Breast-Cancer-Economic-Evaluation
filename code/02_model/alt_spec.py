# -*- coding: utf-8 -*-
"""Alternatif model spesifikasyonları — ANA BORU HATTIYLA BİREBİR aynı kod yolu."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
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
DUS=['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Cilt_Cekintisi_i',
     'Meme_Basi_Retraksiyonu_i','Metastaz_Yeri_i','NAC_Once_Evre_i']

def llr_p(X,yy):
    if X.shape[1]==0 or len(np.unique(yy))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(X,yy)
    q=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
    ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),X.shape[1])

def calistir(VARS, ad):
    Xc=d[VARS].apply(lambda s:s.astype('category').cat.codes).values
    du={v:pd.get_dummies(d[v].astype(str),drop_first=True).values.astype(float) for v in VARS}
    pi=np.zeros((len(y),5)); pm=np.zeros((len(y),5))
    for r in range(5):
        for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=300+r).split(Xc,y,groups=grp):
            sel=[v for v in VARS if llr_p(du[v][tr],y[tr])<0.10] or list(VARS)
            Xs=d[sel].astype(str)
            m1=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
                ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xs.iloc[tr],y[tr])
            pi[te,r]=m1.predict_proba(Xs.iloc[te])[:,1]
            m2=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
                min_samples_leaf=15,random_state=0).fit(Xc[tr],y[tr])
            pm[te,r]=m2.predict_proba(Xc[te])[:,1]
    K=(pi.mean(1)+pm.mean(1))/2
    m=elig&(K>=0.80); tp=int(y[m].sum()); n=int(m.sum()); tp1=int(y[lumA].sum())
    return dict(model=ad, k=len(VARS),
        AUC_tum=roc_auc_score(y,K), AUC_uygun=roc_auc_score(y[elig],K[elig]),
        Brier=brier_score_loss(y[elig],K[elig]),
        isaretlenen=n, uyumlu=tp, uyumsuz=n-tp, PPV=tp/n if n else np.nan,
        ek_uyumlu=tp-tp1, ek_uyumsuz=(n-tp)-(int(lumA.sum())-tp1))

SPEC=[(P+O+R,'Full model'),(P+O,'Radiology block removed'),(P,'Pathology only'),(DUS,'Low-missingness only')]
D=pd.DataFrame([calistir(v,a) for v,a in SPEC])
print(D.round(3).to_string(index=False))
D.to_csv('../../outputs/alternatif_modeller.csv',index=False)
print('\nAna metin referansı: AUC 0.820 / 0.772 · işaretlenen 108 · doğru 98 · PPV 0.907 · +43')
