# -*- coding: utf-8 -*-
"""İki talep:
(A) Model-tabanlı katlama-içi imputasyon — modal imputasyonun etkisini sızıntıdan ayırmak
(B) Tam iç içe boru hattı: imputasyon + kodlama + omurga seçimi + değişken seçimi hepsi dış katlamanın içinde
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2

P=['Histolojik_Yeni','ER','PR','HER2','Molekuler','Ki_67','TIL','HistolojikG']
O=['Metastaz_Yeri','NAC_Once_Evre']
R=['BI-RADS','Meme_Dansite','Lokalizasyon','Lezyon_Turu','Kitle_Sekli','Kitle_Konturu',
   'Kitle_Dansitesi','Kalsifikasyon_Morfolojisi','Kalsifikasyon_Dagilimi','Asimetri',
   'Multifokalite_Durumu','Cilt_Cekintisi','Meme_Basi_Retraksiyonu']
BLOK={'P':P,'O':O,'R':R}
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values; grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
lumA=(d.Molekuler_i.values=='Luminal A')
TUM=P+O+R
VAR=[v for v in TUM if v in d.columns]
EKSIK={'Bilinmiyor','nan','','Bilinmeyen','None','NaN'}
X0=d[VAR].astype(str)
M=X0.isin(EKSIK)   # eksiklik maskesi

def imp_modal(X, tr):
    Xi=X.copy()
    for c in X.columns:
        g=X[c].iloc[tr]; g=g[~g.isin(EKSIK)]
        mod=g.mode().iloc[0] if len(g) else 'Bilinmiyor'
        Xi[c]=Xi[c].where(~Xi[c].isin(EKSIK), mod)
    return Xi

def imp_model(X, tr, seed=0):
    """Model-tabanlı: her değişken için diğerlerinden sınıflandırıcı; YALNIZCA eğitim katlamasında kurulur."""
    Xi=imp_modal(X, tr)           # başlangıç doldurma
    kod={c:{v:i for i,v in enumerate(sorted(Xi[c].unique()))} for c in X.columns}
    for tur in range(2):          # iki geçiş
        for c in X.columns:
            eks=X[c].isin(EKSIK).values
            if eks.sum()==0: continue
            dig=[k for k in X.columns if k!=c]
            Z=np.column_stack([Xi[k].map(kod[k]).values for k in dig])
            eg=np.zeros(len(X),bool); eg[tr]=True; eg &= ~eks
            if eg.sum()<20 or Xi[c][eg].nunique()<2: continue
            yy=Xi[c][eg].map(kod[c]).values
            m=HistGradientBoostingClassifier(max_iter=60,max_leaf_nodes=7,learning_rate=0.15,
                min_samples_leaf=10,random_state=seed).fit(Z[eg],yy)
            ters={i:v for v,i in kod[c].items()}
            tah=m.predict(Z[eks])
            Xi.loc[eks,c]=[ters[int(t)] for t in tah]
    return Xi

def llr_p(Xd,yy):
    if Xd.shape[1]==0 or len(np.unique(yy))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(Xd,yy)
    q=np.clip(m.predict_proba(Xd)[:,1],1e-9,1-1e-9)
    ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),Xd.shape[1])

def fit_pred(Xi, vars_, tr, te, seed, dummy=True):
    du={v:pd.get_dummies(Xi[v],drop_first=True).values.astype(float) for v in vars_}
    sel=[v for v in vars_ if llr_p(du[v][tr],y[tr])<0.10] or list(vars_)
    m1=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
        ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xi[sel].iloc[tr],y[tr])
    p1=m1.predict_proba(Xi[sel].iloc[te])[:,1]
    Xb=(pd.get_dummies(Xi[vars_],drop_first=True).values.astype(float) if dummy
        else Xi[vars_].apply(lambda s:s.astype('category').cat.codes).values)
    m2=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
        min_samples_leaf=15,random_state=seed).fit(Xb[tr],y[tr])
    return (p1+m2.predict_proba(Xb[te])[:,1])/2

def fit_pred_ic(Xi, vars_, tr, te, yy_all, seed):
    du={v:pd.get_dummies(Xi[v],drop_first=True).values.astype(float) for v in vars_}
    sel=[v for v in vars_ if llr_p(du[v][tr],yy_all[tr])<0.10] or list(vars_)
    m1=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
        ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xi[sel].iloc[tr],yy_all[tr])
    p1=m1.predict_proba(Xi[sel].iloc[te])[:,1]
    Xb=pd.get_dummies(Xi[vars_],drop_first=True).values.astype(float)
    m2=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
        min_samples_leaf=15,random_state=seed).fit(Xb[tr],yy_all[tr])
    return (p1+m2.predict_proba(Xb[te])[:,1])/2

def ozet(K, ad):
    m=elig&(K>=0.80); tp=int(y[m].sum()); n=int(m.sum()); tp1=int(y[lumA].sum())
    return dict(analiz=ad, AUC_tum=roc_auc_score(y,K), AUC_uygun=roc_auc_score(y[elig],K[elig]),
        Brier=brier_score_loss(y[elig],K[elig]), isaretlenen=n, uyumlu=tp, uyumsuz=n-tp,
        PPV=tp/n if n else np.nan, ek_uyumlu=tp-tp1)

S=list(pd.read_csv('../../outputs/imputasyon_yontem.csv').to_dict('records'))
# ═══ A: zaten hesaplandı
for imp,ad in []:
    pi=np.zeros((len(y),5)); pm=np.zeros((len(y),5))
    for r in range(5):
        for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=300+r).split(X0,y,groups=grp):
            Xi=imp(X0,tr) if imp is imp_modal else imp(X0,tr,r)
            pr=fit_pred(Xi,VAR,tr,te,r,dummy=False)
            pi[te,r]=pr; pm[te,r]=pr
    S.append(ozet(pi.mean(1), ad)); print(' ',ad,'bitti')
D=pd.DataFrame(S); print(D.round(3).to_string(index=False))
D.to_csv('../../outputs/imputasyon_yontem.csv',index=False)

# ═══ B: TAM İÇ İÇE — imputasyon + kodlama + omurga seçimi + değişken seçimi hepsi dış katlamanın içinde
KOMB=[('P',),('P','O'),('P','R'),('P','O','R')]
oof=np.full(len(y),np.nan); secilen=[]
for r in range(3):
    for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=700+r).split(X0,y,groups=grp):
        Xi=imp_modal(X0,tr)                       # imputasyon YALNIZCA eğitim katlamasından
        en,ea=None,-1
        for k in KOMB:
            vs=[v for v in sum([BLOK[b] for b in k],[]) if v in VAR]
            ip=np.full(len(tr),np.nan)
            Xtr=X0.iloc[tr].reset_index(drop=True); ytr=y[tr]; gtr=grp[tr]
            for itr,ite in StratifiedGroupKFold(3,shuffle=True,random_state=11).split(Xtr,ytr,groups=gtr):
                Xii=imp_modal(Xtr, itr)
                ip[ite]=fit_pred_ic(Xii,vs,itr,ite,ytr,r)
            try: a=roc_auc_score(y[tr],ip)
            except Exception: a=0
            if a>ea: ea,en=a,k
        secilen.append(en)
        vs=[v for v in sum([BLOK[b] for b in en],[]) if v in VAR]
        oof[te]=fit_pred(Xi,vs,tr,te,r,dummy=True)
    print(f'  tam iç içe tekrar {r+1}/3')
import collections
print('\nSeçilen blok:',collections.Counter(secilen).most_common())
S.append(ozet(oof,'Fully nested (imputation + coding + backbone selection inside folds)'))
D=pd.DataFrame(S)
print(D.round(3).to_string(index=False))
D.to_csv('../../outputs/imputasyon_yontem.csv',index=False)
np.save('../../outputs/oof_tam_icice.npy',oof)
