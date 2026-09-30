"""BOOTSTRAP GA + KARAR EGRISI ANALIZI — nihai omurga (23 degisken, Kohort v15).
Bootstrap hasta duzeyinde (kume) yapilir; bilateral olgular birlikte cekilir.
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

d0 = pd.read_excel('../../data/Kohort_v17.xlsx')
P = ['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','TIL_i','HistolojikG_i']
O = ['Metastaz_Yeri_i','NAC_Once_Evre_i']
R = ['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Kitle_Sekli_i',
     'Kitle_Konturu_i','Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i',
     'Kalsifikasyon_Dagilimi_i','Asimetri_i','Multifokalite_Durumu_59','Cilt_Cekintisi_i',
     'Meme_Basi_Retraksiyonu_i']
VARS = P + O + R
ESIK = 0.80

d = d0[d0.Analiz_Kohortu == 'taban senaryo (cM0)'].reset_index(drop=True)
y = (d.RCB_Kategorize >= 2).astype(int).values
grp = d.Hasta_ID.values
elig = np.isin(d.Molekuler_i.values, ['Luminal A', 'Luminal B (HER2 Negatif)'])
s1 = (d.Molekuler_i == 'Luminal A').values
lumB = (d.Molekuler_i == 'Luminal B (HER2 Negatif)').values
Xc = d[VARS].apply(lambda s: s.astype('category').cat.codes).values
dums = {v: pd.get_dummies(d[v].astype(str), drop_first=True).values.astype(float) for v in VARS}

def llr_p(X, yy):
    if X.shape[1] == 0 or len(np.unique(yy)) < 2: return 1.0
    m = LogisticRegression(C=1e6, max_iter=2000).fit(X, yy)
    p = np.clip(m.predict_proba(X)[:, 1], 1e-9, 1-1e-9)
    ll = (yy*np.log(p)+(1-yy)*np.log(1-p)).sum(); p0 = np.clip(yy.mean(), 1e-9, 1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(p0)+(1-yy)*np.log(1-p0)).sum()), 0), X.shape[1])

pi = np.zeros((len(y), 5)); pm = np.zeros((len(y), 5))
for r in range(5):
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=300+r).split(Xc, y, groups=grp):
        sel = [v for v in VARS if llr_p(dums[v][tr], y[tr]) < 0.10]
        Xs = d[sel].astype(str)
        pi[te, r] = Pipeline([('oh', OneHotEncoder(handle_unknown='ignore', drop='first',
                              sparse_output=False)),
                              ('m', LogisticRegression(penalty='l2', C=0.05, max_iter=5000))
                              ]).fit(Xs.iloc[tr], y[tr]).predict_proba(Xs.iloc[te])[:, 1]
        pm[te, r] = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05,
                        max_leaf_nodes=15, min_samples_leaf=15, random_state=0
                        ).fit(Xc[tr], y[tr]).predict_proba(Xc[te])[:, 1]
p = (pi.mean(1) + pm.mean(1)) / 2
np.save('../../outputs/oof_final.npy', p)

def olcumler(idx):
    yy, pp = y[idx], p[idx]
    el, ss1, lb = elig[idx], s1[idx], lumB[idx]
    m = el & (pp >= ESIK); n = int(m.sum()); tp = int(yy[m].sum())
    m1 = lb & (pp >= ESIK); n1 = int(m1.sum())
    return dict(
        AUC=roc_auc_score(yy, pp) if len(np.unique(yy)) > 1 else np.nan,
        Brier=brier_score_loss(yy, pp),
        uygun=int(el.sum()), n=n, dogru=tp, yanlis=n-tp, PPV=tp/max(n, 1),
        kohort_pay=100*n/len(yy),
        S1_n=int(ss1.sum()), S1_dogru=int(yy[ss1].sum()), S1_PPV=yy[ss1].mean() if ss1.sum() else np.nan,
        net_dogru=tp-int(yy[ss1].sum()), net_yanlis=(n-tp)-int((ss1 & (yy == 0)).sum()),
        lumB_n=n1, lumB_dogru=int(yy[m1].sum()), lumB_PPV=int(yy[m1].sum())/max(n1, 1))

goz = olcumler(np.arange(len(y)))
print('=== NIHAI TABAN SENARYO (23 degisken, esik %.2f)' % ESIK)
for k, v in goz.items(): print(f'  {k:12s} {v:.3f}' if isinstance(v, float) else f'  {k:12s} {v}')

# ---- bootstrap: hasta duzeyinde kume cekimi
B = 2000
hastalar = pd.unique(grp)
idx_by = {h: np.where(grp == h)[0] for h in hastalar}
rng = np.random.default_rng(42)
kayit = []
for b in range(B):
    sec = rng.choice(hastalar, size=len(hastalar), replace=True)
    idx = np.concatenate([idx_by[h] for h in sec])
    if len(np.unique(y[idx])) < 2: continue
    kayit.append(olcumler(idx))
bs = pd.DataFrame(kayit)
print(f'\n=== BOOTSTRAP %95 GA (B={len(bs)}, hasta duzeyinde kume cekimi)')
satir = [('Ayrım gücü (AUC)','AUC',3), ('Kalibrasyon (Brier)','Brier',3),
         ('Uygun popülasyon (n)','uygun',0), ('İşaretlenen olgu','n',0),
         ('Kohortun yüzdesi','kohort_pay',1), ('Doğru yönlendirme','dogru',0),
         ('Yanlış yönlendirme','yanlis',0), ('Kestirim değeri (PPV)','PPV',3),
         ('S1 kestirim değeri','S1_PPV',3), ('S1 üzerine ek doğru','net_dogru',0),
         ('S1 üzerine ek yanlış','net_yanlis',0),
         ('Luminal B (HER2−) işaretlenen','lumB_n',0),
         ('Luminal B (HER2−) PPV','lumB_PPV',3)]
print(f'{"Parametre":32s} {"Gözlenen":>10s} {"%95 GA":>20s}')
tab = []
for ad, k, nd in satir:
    lo, hi = np.nanpercentile(bs[k], [2.5, 97.5])
    f = f'%.{nd}f'
    print(f'{ad:32s} {f%goz[k]:>10s} {f%lo} – {f%hi:>10s}')
    tab.append(dict(parametre=ad, gozlenen=round(goz[k], nd), alt=round(lo, nd), ust=round(hi, nd)))
pd.DataFrame(tab).to_csv('../../outputs/bootstrap_GA.csv', index=False)

# ---- karar egrisi analizi
print('\n=== KARAR EGRISI ANALIZI (net fayda)')
print(f'{"pt":>5} {"karar kurali":>13} {"S1 alt tip":>11} {"hepsini isaretle":>17} {"hicbiri":>8}')
dca = []
N = len(y)
for pt in [0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90]:
    w = pt/(1-pt)
    def nf(mask):
        tp = int(y[mask].sum()); fp = int(mask.sum())-tp
        return tp/N - (fp/N)*w
    kk = nf(elig & (p >= pt)); ss = nf(s1); hh = nf(np.ones(N, bool))
    print(f'{pt:5.2f} {kk:13.3f} {ss:11.3f} {hh:17.3f} {0.0:8.3f}')
    dca.append(dict(pt=pt, karar_kurali=round(kk, 3), S1=round(ss, 3),
                    hepsi=round(hh, 3), hicbiri=0.0))
pd.DataFrame(dca).to_csv('../../outputs/karar_egrisi.csv', index=False)
