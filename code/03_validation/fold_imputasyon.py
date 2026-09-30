# -*- coding: utf-8 -*-
"""Hakem talebi: katlama-içi imputasyon + kukla kodlama duyarlılığı + katmanlı AUC güven aralıkları."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2

HAM={'Histolojik_Yeni':'Histolojik_Yeni','ER_i':'ER','PR_i':'PR','HER2_i':'HER2','Molekuler_i':'Molekuler',
'Ki_67_i':'Ki_67','TIL_i':'TIL','HistolojikG_i':'HistolojikG','Metastaz_Yeri_i':'Metastaz_Yeri',
'NAC_Once_Evre_i':'NAC_Once_Evre','BI-RADS_i':'BI-RADS','Meme_Dansite_i':'Meme_Dansite',
'Lokalizasyon_i':'Lokalizasyon','Lezyon_Turu_i':'Lezyon_Turu','Kitle_Sekli_i':'Kitle_Sekli',
'Kitle_Konturu_i':'Kitle_Konturu','Kitle_Dansitesi_i':'Kitle_Dansitesi',
'Kalsifikasyon_Morfolojisi_i':'Kalsifikasyon_Morfolojisi','Kalsifikasyon_Dagilimi_i':'Kalsifikasyon_Dagilimi',
'Asimetri_i':'Asimetri','Multifokalite_Durumu_59':'Multifokalite_Durumu','Cilt_Cekintisi_i':'Cilt_Cekintisi',
'Meme_Basi_Retraksiyonu_i':'Meme_Basi_Retraksiyonu'}
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values; grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
lumA=(d.Molekuler_i.values=='Luminal A')
IMP=list(HAM.keys())
RAW=[HAM[v] if HAM[v] in d.columns else v for v in IMP]
EKSIK={'Bilinmiyor','nan','','Bilinmeyen','None'}

def llr_p(X,yy):
    if X.shape[1]==0 or len(np.unique(yy))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(X,yy)
    q=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
    ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),X.shape[1])

def kos(kolonlar, fold_imp=False, dummy_boost=False, etiket=''):
    Xdf=d[kolonlar].astype(str)
    pi=np.zeros((len(y),5)); pm=np.zeros((len(y),5))
    for r in range(5):
        for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=300+r).split(Xdf,y,groups=grp):
            X=Xdf.copy()
            if fold_imp:
                # EKSİK değerler YALNIZCA eğitim katlamasının modu ile doldurulur
                for c in kolonlar:
                    egit=X[c].iloc[tr]
                    gecerli=egit[~egit.isin(EKSIK)]
                    mod=gecerli.mode().iloc[0] if len(gecerli) else 'Bilinmiyor'
                    X[c]=X[c].where(~X[c].isin(EKSIK), mod)
            du={v:pd.get_dummies(X[v],drop_first=True).values.astype(float) for v in kolonlar}
            sel=[v for v in kolonlar if llr_p(du[v][tr],y[tr])<0.10] or list(kolonlar)
            m1=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
                ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(X[sel].iloc[tr],y[tr])
            pi[te,r]=m1.predict_proba(X[sel].iloc[te])[:,1]
            if dummy_boost:
                Xb=pd.get_dummies(X[kolonlar],drop_first=True).values.astype(float)
            elif fold_imp:
                Xb=X[kolonlar].apply(lambda s:s.astype('category').cat.codes).values
            else:
                Xb=d[kolonlar].apply(lambda s:s.astype('category').cat.codes).values
            m2=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
                min_samples_leaf=15,random_state=0).fit(Xb[tr],y[tr])
            pm[te,r]=m2.predict_proba(Xb[te])[:,1]
    K=(pi.mean(1)+pm.mean(1))/2
    m=elig&(K>=0.80); tp=int(y[m].sum()); n=int(m.sum()); tp1=int(y[lumA].sum())
    return dict(analiz=etiket, AUC_tum=roc_auc_score(y,K), AUC_uygun=roc_auc_score(y[elig],K[elig]),
        Brier=brier_score_loss(y[elig],K[elig]), isaretlenen=n, uyumlu=tp, uyumsuz=n-tp,
        PPV=tp/n if n else np.nan, ek_uyumlu=tp-tp1), K

S=[]
r1,K1=kos(IMP, False, False, 'Base case (pre-imputed data, integer-coded boosting)'); S.append(r1)
r2,_ =kos(RAW,  True,  False, 'Fold-internal imputation'); S.append(r2)
r3,_ =kos(IMP,  False, True,  'Dummy-coded boosting'); S.append(r3)
r4,_ =kos(RAW,  True,  True,  'Fold-internal imputation + dummy-coded boosting'); S.append(r4)
D=pd.DataFrame(S)
print(D.round(3).to_string(index=False))
D.to_csv('../../outputs/duyarlilik_imputasyon_kodlama.csv',index=False)
