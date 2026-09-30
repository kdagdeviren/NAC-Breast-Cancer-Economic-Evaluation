# -*- coding: utf-8 -*-
"""İç içe (nested) gruplandırılmış çapraz doğrulama:
omurga blok seçimi DIŞ katlamanın İÇİNDE yapılır. Hakem sorusu: omurga seçimi sızıntı yaratıyor mu?"""
import numpy as np, pandas as pd, warnings, itertools; warnings.filterwarnings('ignore')
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
BLOK={'P':P,'O':O,'R':R}
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values; grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
TUM=P+O+R
Xc_all=d[TUM].apply(lambda s:s.astype('category').cat.codes)
du={v:pd.get_dummies(d[v].astype(str),drop_first=True).values.astype(float) for v in TUM}

def llr_p(X,yy):
    if X.shape[1]==0 or len(np.unique(yy))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(X,yy)
    q=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
    ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),X.shape[1])

def fit_predict(vars_, tr, te, seed):
    sel=[v for v in vars_ if llr_p(du[v][tr], y[tr])<0.10] or list(vars_)
    Xl=np.hstack([du[v] for v in sel])
    lr=LogisticRegression(C=0.05,max_iter=5000).fit(Xl[tr],y[tr])
    p1=lr.predict_proba(Xl[te])[:,1]
    Xc=Xc_all[list(vars_)].values
    hgb=HistGradientBoostingClassifier(max_iter=300,max_leaf_nodes=15,learning_rate=0.05,
        min_samples_leaf=15,l2_regularization=1.0,random_state=seed).fit(Xc[tr],y[tr])
    return (p1+hgb.predict_proba(Xc[te])[:,1])/2

KOMB=[('P',),('P','O'),('P','R'),('P','O','R')]
oof=np.full(len(y),np.nan); secilen=[]
for r in range(3):
    for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=500+r).split(Xc_all,y,groups=grp):
        # IC dongu: blok kombinasyonunu YALNIZCA egitim verisiyle sec
        en_iyi=None; en_auc=-1
        for k in KOMB:
            vs=sum([BLOK[b] for b in k],[])
            ip=np.full(len(tr),np.nan)
            for itr,ite in StratifiedGroupKFold(3,shuffle=True,random_state=7).split(
                    Xc_all.iloc[tr],y[tr],groups=grp[tr]):
                ip[ite]=fit_predict(vs,tr[itr],tr[ite],seed=r)
            a=roc_auc_score(y[tr],ip)
            if a>en_auc: en_auc,en_iyi=a,k
        secilen.append(en_iyi)
        vs=sum([BLOK[b] for b in en_iyi],[])
        oof[te]=fit_predict(vs,tr,te,seed=r)
    print(f'  tekrar {r+1}/3 bitti')

import collections
print('\nDış katlamalarda seçilen blok kombinasyonu:', collections.Counter(secilen).most_common())
print(f'\n{"":26s} {"AUC":>7s} {"Brier":>7s}')
for ad,m in [('İÇ İÇE — tüm kohort',np.ones(len(y),bool)),('İÇ İÇE — uygun popülasyon',elig)]:
    print(f'{ad:26s} {roc_auc_score(y[m],oof[m]):7.3f} {brier_score_loss(y[m],oof[m]):7.3f}')
print()
p_bildirilen=np.load('../../outputs/oof_final.npy')
for ad,m in [('BİLDİRİLEN — tüm kohort',np.ones(len(y),bool)),('BİLDİRİLEN — uygun popülasyon',elig)]:
    print(f'{ad:26s} {roc_auc_score(y[m],p_bildirilen[m]):7.3f} {brier_score_loss(y[m],p_bildirilen[m]):7.3f}')
np.save('../../outputs/oof_nested.npy',oof)
pd.DataFrame({'oof_nested':oof}).to_csv('../../outputs/icice_cv_tahminleri.csv',index=False)
