# -*- coding: utf-8 -*-
"""Hakem taleplerine yanıt: hedef popülasyonda performans, kalibrasyon,
iyimserlik düzeltmesi, hasta düzeyi çapraz tablo, prevalans duyarlılığı."""
import numpy as np, pandas as pd, warnings, collections
warnings.filterwarnings('ignore')
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
VARS=P+O+R
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values
grp=d.Hasta_ID.values
p=np.load('../../outputs/oof_final.npy')
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
lumA=(d.Molekuler_i.values=='Luminal A'); lumB=(d.Molekuler_i.values=='Luminal B (HER2 Negatif)')

def _irls(X, yy, offset=None):
    """Binom GLM (logit) — IRLS; katsayı ve standart hata döndürür."""
    X=np.asarray(X,float); n,k=X.shape
    off=np.zeros(n) if offset is None else np.asarray(offset,float)
    b=np.zeros(k)
    for _ in range(60):
        eta=X@b+off; mu=1/(1+np.exp(-eta)); w=np.clip(mu*(1-mu),1e-9,None)
        z=eta-off+(yy-mu)/w
        WX=X*w[:,None]
        b_new=np.linalg.solve(X.T@WX, WX.T@z)
        if np.max(np.abs(b_new-b))<1e-10: b=b_new; break
        b=b_new
    eta=X@b+off; mu=1/(1+np.exp(-eta)); w=np.clip(mu*(1-mu),1e-9,None)
    cov=np.linalg.inv(X.T@(X*w[:,None]))
    return b, np.sqrt(np.diag(cov))

def kalib(yy,pp):
    """kalibrasyon kesişimi, eğimi, E:O, Brier, ölçekli Brier"""
    yy=np.asarray(yy,float); q=np.clip(pp,1e-6,1-1e-6); lo=np.log(q/(1-q))
    X=np.column_stack([np.ones(len(yy)), lo])
    b,se=_irls(X,yy); egim=b[1]; egim_se=se[1]
    b0,se0=_irls(np.ones((len(yy),1)), yy, offset=lo)
    kesisim=b0[0]; kes_se=se0[0]
    br=brier_score_loss(yy,pp); pr=yy.mean(); br0=pr*(1-pr)
    return dict(n=len(yy), olay=int(yy.sum()), prevalans=pr,
                AUC=roc_auc_score(yy,pp), Brier=br, Brier_null=br0,
                Brier_olcekli=1-br/br0, E_O=pp.mean()/pr,
                kal_kesisim=kesisim, kk_alt=kesisim-1.96*kes_se, kk_ust=kesisim+1.96*kes_se,
                kal_egim=egim, ke_alt=egim-1.96*egim_se, ke_ust=egim+1.96*egim_se)

print('═══ 1. PERFORMANS — POPÜLASYON KATMANLARINA GÖRE (hakem talebi)\n')
sat=[]
for ad,m in [('Tüm cM0 kohortu',np.ones(len(y),bool)),
             ('Uygun popülasyon (HR+/HER2−)',elig),
             ('Luminal A',lumA),('Luminal B (HER2−)',lumB)]:
    r=kalib(y[m],p[m]); r['katman']=ad; sat.append(r)
K=pd.DataFrame(sat)[['katman','n','olay','prevalans','AUC','Brier','Brier_null','Brier_olcekli',
                     'E_O','kal_kesisim','kk_alt','kk_ust','kal_egim','ke_alt','ke_ust']]
print(K.round(3).to_string(index=False))
K.to_csv('../../outputs/katman_performansi.csv',index=False)

print('\n═══ 2. HASTA DÜZEYİ S1 × S2 ÇAPRAZ TABLO (hakem talebi)\n')
s1=lumA; s2=elig&(p>=0.80)
tab=pd.crosstab([s1[elig],s2[elig]], y[elig], rownames=['S1','S2'], colnames=['RCB-II/III'])
tab.columns=['RCB-0/I','RCB-II/III']; tab['toplam']=tab.sum(1)
print(tab.to_string())
print()
for (a,b),ad in [((True,True),'Her ikisi işaretledi'),((True,False),'Yalnız S1 işaretledi (model eledi)'),
                 ((False,True),'Yalnız S2 işaretledi (modelin net katkısı)'),((False,False),'Hiçbiri')]:
    m=elig&(s1==a)&(s2==b)
    print(f'  {ad:38s} n={m.sum():3d}  RCB-II/III={int(y[m].sum()):3d}  oran={y[m].mean() if m.sum() else 0:.3f}')
tab.to_csv('../../outputs/S1_S2_caprazTablo.csv')

print('\n═══ 3. PREVALANS DUYARLILIĞI — PPV (hakem talebi)\n')
# esikteki duyarlilik ve ozgulluk sabit tutulup prevalans degistirilir
m=elig&(p>=0.80)
TP=int(y[m].sum()); FP=int(m.sum()-TP)
FN=int(y[elig].sum()-TP); TN=int((~y[elig].astype(bool)&elig[elig]).sum()-FP) if False else int((elig.sum()-y[elig].sum())-FP)
duy=TP/(TP+FN); ozg=TN/(TN+FP)
print(f'  Gözlenen: duyarlılık {duy:.3f} · özgüllük {ozg:.3f} · prevalans {y[elig].mean():.3f} · PPV {TP/(TP+FP):.3f}\n')
print(f'  {"prevalans":>10s} {"PPV":>7s} {"NPV":>7s}')
pv=[]
for pre in [0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.766,0.80]:
    ppv=duy*pre/(duy*pre+(1-ozg)*(1-pre)); npv=ozg*(1-pre)/(ozg*(1-pre)+(1-duy)*pre)
    pv.append(dict(prevalans=pre,PPV=ppv,NPV=npv))
    print(f'  {pre:10.3f} {ppv:7.3f} {npv:7.3f}')
pd.DataFrame(pv).to_csv('../../outputs/prevalans_duyarliligi.csv',index=False)
