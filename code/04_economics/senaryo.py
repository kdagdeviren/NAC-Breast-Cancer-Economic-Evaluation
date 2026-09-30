import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2

d0 = pd.read_excel('../../data/Kohort_v17.xlsx')
P=['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','Tubul_i','Nukleer_i',
   'Mitotik_i','HistolojikG_i','TIL_i','Ecaderin_i']
O=['Metasdaz_Durumu_i','Metastaz_Yeri_i','NAC_Once_Evre_i']
R=['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Mimari_i','Kitle_Sekli_i',
   'Kitle_Konturu_i','Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i','Kalsifikasyon_Dagilimi_i',
   'Asimetri_i','Multifokalite_Durumu_59','2_Yildir_Stabil_i','Cilt_Cekintisi_i',
   'Meme_Basi_Retraksiyonu_i','Ameliyat_Oykusu_i','Kozmetik_Implant']
POR=P+O+R
ELIG=['Luminal A','Luminal B (HER2 Negatif)']

def llr_p(X,y):
    if X.shape[1]==0 or len(np.unique(y))<2: return 1.0
    m=LogisticRegression(C=1e6,max_iter=2000).fit(X,y)
    p=np.clip(m.predict_proba(X)[:,1],1e-9,1-1e-9)
    ll=(y*np.log(p)+(1-y)*np.log(1-p)).sum(); p0=np.clip(y.mean(),1e-9,1-1e-9)
    return chi2.sf(max(2*(ll-(y*np.log(p0)+(1-y)*np.log(1-p0)).sum()),0), X.shape[1])

d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values
grp=d.Hasta_ID.values
elig=np.isin(d.Molekuler_i.values,ELIG)
s1=(d.Molekuler_i=='Luminal A').values
print(f'cM0 kohort: {len(d)} kayit / {pd.Series(grp).nunique()} hasta | pozitif {y.sum()} ({100*y.mean():.1f}%)')
print(f'uygun (HR+/HER2-): {elig.sum()} | S1 (Luminal A): {s1.sum()}, dogru {int(y[s1].sum())}, '
      f'yanlis {int((~y[s1].astype(bool)).sum())}, PPV {y[s1].mean():.3f}')

# Madde 4: E-kaderin cikarildi. Madde 9: grade bilesenleri (tubul/nukleer/mitotik)
# grade ile esdogrusal oldugu icin cikarildi -> 28 degisken
# Madde 4  : E-kaderin (guvenilmez, %65,5 eksik)
# Madde 9  : grade bilesenleri (grade ile esdogrusal)
# Madde 15 : Metastaz Durumu (Metastaz Yeri'nin ikili turevi)
# Madde 7b : Mimari (%62,8 eksik) ve 2_Yildir_Stabil (%97 eksik) performansi dusuruyor;
#            Kozmetik_Implant ve Ameliyat_Oykusu sifir katki
DISLA=['Ecaderin_i','Tubul_i','Nukleer_i','Mitotik_i','Metasdaz_Durumu_i',
       'Mimari_i','2_Yildir_Stabil_i','Kozmetik_Implant','Ameliyat_Oykusu_i']
for ad,VARS in [('23 degisken (taban omurga)',[v for v in POR if v not in DISLA])]:
    Xc=d[VARS].apply(lambda s:s.astype('category').cat.codes).values
    dums={v:pd.get_dummies(d[v].astype(str),drop_first=True).values.astype(float) for v in VARS}
    pi=np.zeros((len(y),5)); pm=np.zeros((len(y),5))
    for r in range(5):
        for tr,te in StratifiedGroupKFold(5,shuffle=True,random_state=300+r).split(Xc,y,groups=grp):
            sel=[v for v in VARS if llr_p(dums[v][tr],y[tr])<0.10]; Xs=d[sel].astype(str)
            pi[te,r]=Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
                ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xs.iloc[tr],y[tr]).predict_proba(Xs.iloc[te])[:,1]
            pm[te,r]=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,
                min_samples_leaf=15,random_state=0).fit(Xc[tr],y[tr]).predict_proba(Xc[te])[:,1]
    p=(pi.mean(1)+pm.mean(1))/2
    print(f'\n=== {ad}  AUC={roc_auc_score(y,p):.3f}  Brier={brier_score_loss(y,p):.3f}')
    print(f'   {"esik":>5} {"n":>4} {"dogru":>6} {"yanlis":>7} {"PPV":>7}  S1 uzerine')
    for t in [0.70,0.75,0.80,0.85]:
        s=elig&(p>=t); n=int(s.sum()); tp=int(y[s].sum())
        print(f'   {t:5.2f} {n:4d} {tp:6d} {n-tp:7d} {tp/max(n,1):7.3f}  '
              f'+{tp-int(y[s1].sum())} dogru {(n-tp)-int((s1&(y==0)).sum()):+d} yanlis')
    lb=(d.Molekuler_i=='Luminal B (HER2 Negatif)').values
    s=lb&(p>=0.80); n=int(s.sum()); tp=int(y[s].sum())
    print(f'   Luminal B (HER2-) alt grubu: {int(lb.sum())} hastadan {n} isaretlenir, '
          f'{tp} dogru, {n-tp} yanlis, PPV={tp/max(n,1):.3f}')
