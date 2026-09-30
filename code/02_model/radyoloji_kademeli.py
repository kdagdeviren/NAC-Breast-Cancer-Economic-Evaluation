"""MADDE 7b — radyoloji blogunun kademeli eklenmesi.
17 radyoloji degiskeni ham veri eksikligine gore siralanir; en az eksik olandan baslayarak
teker teker eklenir. Amac: imputasyona en az bagimli, en iyi performansli alt kumeyi bulmak.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.metrics import roc_auc_score, brier_score_loss
from scipy.stats import chi2

d0 = pd.read_excel('../../outputs/Kohort_v16.xlsx')
ham = pd.read_excel('../../data/Kagan_TEZ.xlsx', sheet_name=0)

P = ['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','TIL_i','HistolojikG_i']
O = ['Metastaz_Yeri_i','NAC_Once_Evre_i']
R = ['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Mimari_i','Kitle_Sekli_i',
     'Kitle_Konturu_i','Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i',
     'Kalsifikasyon_Dagilimi_i','Asimetri_i','Multifokalite_Durumu_59','2_Yildir_Stabil_i',
     'Cilt_Cekintisi_i','Meme_Basi_Retraksiyonu_i','Ameliyat_Oykusu_i','Kozmetik_Implant']
HAM_AD = {'BI-RADS_i':'BI-RADS','Meme_Dansite_i':'Meme_Dansite','Lokalizasyon_i':'Lokalizasyon',
  'Lezyon_Turu_i':'Lezyon_Turu','Mimari_i':'Mimari','Kitle_Sekli_i':'Kitle_Sekli',
  'Kitle_Konturu_i':'Kitle_Konturu','Kitle_Dansitesi_i':'Kitle_Dansitesi',
  'Kalsifikasyon_Morfolojisi_i':'Kalsifikasyon_Morfolojisi',
  'Kalsifikasyon_Dagilimi_i':'Kalsifikasyon_Dagilimi','Asimetri_i':'Asimetri',
  'Multifokalite_Durumu_59':'Multifokalite_Durumu','2_Yildir_Stabil_i':'2_Yildir_Stabil',
  'Cilt_Cekintisi_i':'Cilt_Cekintisi','Meme_Basi_Retraksiyonu_i':'Meme_Basi_Retraksiyonu',
  'Ameliyat_Oykusu_i':'Ameliyat_Oykusu','Kozmetik_Implant':'Kozmetik_Implant'}

eksik = {}
for v in R:
    h = HAM_AD.get(v, v)
    if h in ham.columns:
        s = ham[h]
        eksik[v] = int((s.isna() | s.astype(str).str.contains('Bilinmiyor|Değerlendirilemedi',
                                                              na=False)).sum())
    else:
        eksik[v] = 0
SIRA = sorted(R, key=lambda v: (eksik[v], v))

print('RADYOLOJI DEGISKENLERI — ham veri eksikligine gore sirali')
for i, v in enumerate(SIRA, 1):
    print(f'  {i:2d}. {v:32s} eksik {eksik[v]:3d} (%{100*eksik[v]/len(ham):4.1f})')

d = d0[d0.Analiz_Kohortu == 'taban senaryo (cM0)'].reset_index(drop=True)
y = (d.RCB_Kategorize >= 2).astype(int).values
grp = d.Hasta_ID.values
elig = np.isin(d.Molekuler_i.values, ['Luminal A', 'Luminal B (HER2 Negatif)'])
s1 = (d.Molekuler_i == 'Luminal A').values
TUM = P + O + R
dums = {v: pd.get_dummies(d[v].astype(str), drop_first=True).values.astype(float) for v in TUM}

def llr_p(X, yy):
    if X.shape[1] == 0 or len(np.unique(yy)) < 2: return 1.0
    m = LogisticRegression(C=1e6, max_iter=2000).fit(X, yy)
    p = np.clip(m.predict_proba(X)[:, 1], 1e-9, 1 - 1e-9)
    ll = (yy*np.log(p) + (1-yy)*np.log(1-p)).sum()
    p0 = np.clip(yy.mean(), 1e-9, 1-1e-9)
    return chi2.sf(max(2*(ll - (yy*np.log(p0)+(1-yy)*np.log(1-p0)).sum()), 0), X.shape[1])

# katlamalar tum konfigurasyonlarda ayni -> tarama p-degerlerini bir kez hesapla
X0 = d[TUM].apply(lambda s: s.astype('category').cat.codes).values
katlar, pdeg = [], []
for r in range(5):
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=300+r).split(X0, y, groups=grp):
        katlar.append((tr, te))
        pdeg.append({v: llr_p(dums[v][tr], y[tr]) for v in TUM})
print(f'\n{len(katlar)} katlama icin tarama tamamlandi\n')

print(f'{"k":>3} {"degisken":>9} {"eklenen":32s} {"AUC":>6} {"Brier":>6} '
      f'{"n":>4} {"dogru":>6} {"yanlis":>7} {"PPV":>6} {"S1+":>5}')
sonuc = []
for k in range(0, 18):
    VARS = P + O + SIRA[:k]
    Xc = d[VARS].apply(lambda s: s.astype('category').cat.codes).values
    pi = np.zeros((len(y), 5)); pm = np.zeros((len(y), 5))
    for j, (tr, te) in enumerate(katlar):
        r = j // 5
        sel = [v for v in VARS if pdeg[j][v] < 0.10]
        Xs = d[sel].astype(str)
        pi[te, r] = Pipeline([('oh', OneHotEncoder(handle_unknown='ignore', drop='first',
                                                   sparse_output=False)),
                              ('m', LogisticRegression(penalty='l2', C=0.05, max_iter=5000))
                              ]).fit(Xs.iloc[tr], y[tr]).predict_proba(Xs.iloc[te])[:, 1]
        pm[te, r] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05,
                        max_leaf_nodes=15, min_samples_leaf=15, random_state=0
                        ).fit(Xc[tr], y[tr]).predict_proba(Xc[te])[:, 1]
    p = (pi.mean(1) + pm.mean(1)) / 2
    s = elig & (p >= 0.80); n = int(s.sum()); tp = int(y[s].sum())
    auc = roc_auc_score(y, p); br = brier_score_loss(y, p)
    ek = SIRA[k-1] if k > 0 else '— (radyoloji yok)'
    print(f'{k:3d} {len(VARS):9d} {ek:32s} {auc:6.3f} {br:6.3f} {n:4d} {tp:6d} {n-tp:7d} '
          f'{tp/max(n,1):6.3f} {tp-int(y[s1].sum()):+5d}')
    sonuc.append(dict(k=k, degisken=len(VARS), eklenen=ek, AUC=round(auc,3),
                      Brier=round(br,3), n=n, dogru=tp, yanlis=n-tp,
                      PPV=round(tp/max(n,1),3), S1_uzerine=tp-int(y[s1].sum()),
                      kumulatif_eksik=sum(eksik[v] for v in SIRA[:k])))

t = pd.DataFrame(sonuc)
t.to_csv('../../outputs/radyoloji_kademeli.csv', index=False)
print('\nEN IYI AUC:', t.loc[t.AUC.idxmax(), ['k','degisken','AUC','PPV','n']].to_dict())
print('EN IYI PPV:', t.loc[t.PPV.idxmax(), ['k','degisken','AUC','PPV','n']].to_dict())
