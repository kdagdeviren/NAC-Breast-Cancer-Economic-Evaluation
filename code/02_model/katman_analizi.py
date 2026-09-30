"""LR ve GRADIENT BOOSTING KATMANLARININ AYRINTILI KARSILASTIRMASI — Kohort v17"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2
import collections

d0=pd.read_excel('../../data/Kohort_v17.xlsx')
P=['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','TIL_i','HistolojikG_i']
O=['Metastaz_Yeri_i','NAC_Once_Evre_i']
R=['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Kitle_Sekli_i','Kitle_Konturu_i',
   'Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i','Kalsifikasyon_Dagilimi_i','Asimetri_i',
   'Multifokalite_Durumu_59','Cilt_Cekintisi_i','Meme_Basi_Retraksiyonu_i']
VARS=P+O+R
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values; grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
s1=(d.Molekuler_i=='Luminal A').values
Xc=d[VARS].apply(lambda s:s.astype('category').cat.codes).values
du={v:pd.get_dummies(d[v].astype(str),drop_first=True).values.astype(float) for v in VARS}
def llr_p(X,yy):
    if X.shape[1]==0 or len(np.unique(yy))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(X,yy)
    q=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
    ll=(yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0=np.clip(yy.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0),X.shape[1])

pi=np.zeros((len(y),5)); pm=np.zeros((len(y),5)); pst=np.zeros((len(y),5))
nsel=[]; secim=collections.Counter()
for r in range(5):
    for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=300+r).split(Xc,y,groups=grp):
        sel=[v for v in VARS if llr_p(du[v][tr],y[tr])<0.10]; nsel.append(len(sel)); secim.update(sel)
        Xs=d[sel].astype(str)
        m1=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
            ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xs.iloc[tr],y[tr])
        pi[te,r]=m1.predict_proba(Xs.iloc[te])[:,1]
        m2=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
            min_samples_leaf=15,random_state=0).fit(Xc[tr],y[tr])
        pm[te,r]=m2.predict_proba(Xc[te])[:,1]
        # istifleme: ILR skoru ML ye ozellik olarak
        p1tr=np.zeros(len(tr))
        for itr,ite in StratifiedGroupKFold(5,shuffle=True,random_state=1).split(Xc[tr],y[tr],groups=grp[tr]):
            p1tr[ite]=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
                ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xs.iloc[tr[itr]],y[tr[itr]]
                ).predict_proba(Xs.iloc[tr[ite]])[:,1]
        m3=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
            min_samples_leaf=15,random_state=0).fit(np.column_stack([Xc[tr],p1tr]),y[tr])
        pst[te,r]=m3.predict_proba(np.column_stack([Xc[te],pi[te,r]]))[:,1]
A=pi.mean(1); B=pm.mean(1); ST=pst.mean(1); K=(A+B)/2
print('='*92); print('KATMAN KARŞILAŞTIRMASI — Kohort v17, hasta-gruplu 5×5 CV'); print('='*92)
print(f'{"Yapı":42s} {"AUC":>7s} {"Brier":>7s} {"eşik 0,80":>26s}')
for ad,pr in [('Katman 1 — lojistik regresyon (ILR)',A),('Katman 2 — gradyan artırma (ML)',B),
              ('İstifleme (ILR skoru ML’ye özellik)',ST),('KASKAD — olasılık ortalaması',K)]:
    m=elig&(pr>=0.80); n=int(m.sum()); tp=int(y[m].sum())
    print(f'{ad:42s} {roc_auc_score(y,pr):7.3f} {brier_score_loss(y,pr):7.3f}   n={n:3d} doğru={tp:3d} yanlış={n-tp:2d} PPV={tp/max(n,1):.3f}')
print(f'\nkatlama başına seçilen değişken: ort {np.mean(nsel):.1f} (min {min(nsel)}, maks {max(nsel)})')
print('\n25 katlamanın kaçında seçildi:')
for v,c in secim.most_common():
    print(f'   {v:32s} {c:2d}/25')
print('\nKatmanlar arası korelasyon: %.3f'%np.corrcoef(A,B)[0,1])
np.save('../../outputs/oof_final.npy',K)
