"""IKI PERDELI SUZGEC — ikinci katman, birinci katmanin gecirdigi alt grup icin
ozel egitilmis ayri bir model olarak sinanir.

Perde 1: p >= 0.70  -> geniş havuz
Perde 2: bu alt grup uzerinde egitilmis ikinci model, yanlis pozitifleri ayiklar

Sizinti onlemi: ikinci katman her katlamanin YALNIZCA egitim verisindeki alt grup
uzerinde egitilir, test verisindeki alt gruba uygulanir.
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import HistGradientBoostingClassifier
from scipy.stats import chi2

d0 = pd.read_excel('../../data/Kohort_v17.xlsx')
P = ['Histolojik_Yeni','ER_i','PR_i','HER2_i','Molekuler_i','Ki_67_i','TIL_i','HistolojikG_i']
O = ['Metastaz_Yeri_i','NAC_Once_Evre_i']
R = ['BI-RADS_i','Meme_Dansite_i','Lokalizasyon_i','Lezyon_Turu_i','Kitle_Sekli_i',
     'Kitle_Konturu_i','Kitle_Dansitesi_i','Kalsifikasyon_Morfolojisi_i',
     'Kalsifikasyon_Dagilimi_i','Asimetri_i','Multifokalite_Durumu_59','Cilt_Cekintisi_i',
     'Meme_Basi_Retraksiyonu_i']
VARS = P + O + R
d = d0[d0.Analiz_Kohortu == 'taban senaryo (cM0)'].reset_index(drop=True)
y = (d.RCB_Kategorize >= 2).astype(int).values
grp = d.Hasta_ID.values
elig = np.isin(d.Molekuler_i.values, ['Luminal A','Luminal B (HER2 Negatif)'])
s1 = (d.Molekuler_i == 'Luminal A').values
Xc = d[VARS].apply(lambda s: s.astype('category').cat.codes).values
dums = {v: pd.get_dummies(d[v].astype(str), drop_first=True).values.astype(float) for v in VARS}
p1 = np.load('../../outputs/oof_final.npy')     # birinci katman OOF olasiliklari

def llr_p(X, yy):
    if X.shape[1] == 0 or len(np.unique(yy)) < 2: return 1.0
    m = LogisticRegression(C=1e6, max_iter=2000).fit(X, yy)
    q = np.clip(m.predict_proba(X)[:,1], 1e-9, 1-1e-9)
    ll = (yy*np.log(q)+(1-yy)*np.log(1-q)).sum(); q0 = np.clip(yy.mean(), 1e-9, 1-1e-9)
    return chi2.sf(max(2*(ll-(yy*np.log(q0)+(1-yy)*np.log(1-q0)).sum()),0), X.shape[1])

def ozet(mask, ad):
    n = int(mask.sum()); tp = int(y[mask].sum()); fp = n-tp
    kacan = int((elig & (y==1)).sum()) - tp
    print(f'  {ad:44s} n={n:4d} doğru={tp:4d} yanlış={fp:3d} PPV={tp/max(n,1):.3f} '
          f'kaçırılan={kacan:3d} S1+{tp-int(y[s1].sum()):+3d}/{fp-int((s1&(y==0)).sum()):+3d}')
    return n, tp, fp, kacan

print('='*104); print('REFERANS — tek eşik'); print('='*104)
for t in [0.70, 0.75, 0.80, 0.85]:
    ozet(elig & (p1>=t), f'tek eşik {t:.2f}')

# ---------------------------------------------------------------- IKINCI PERDE
print('\n'+'='*104)
print('İKİ PERDELİ SÜZGEÇ — perde 1: p1 ≥ eşik1 · perde 2: alt grup modeli ≥ eşik2')
print('='*104)
for ESIK1 in [0.65, 0.70, 0.75]:
    alt = elig & (p1 >= ESIK1)
    p2 = np.full(len(y), np.nan)
    for r in range(5):
        for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=300+r).split(Xc, y, groups=grp):
            tr_alt = tr[alt[tr]]; te_alt = te[alt[te]]
            if len(tr_alt) < 30 or len(np.unique(y[tr_alt])) < 2 or len(te_alt) == 0: continue
            sel = [v for v in VARS if llr_p(dums[v][tr_alt], y[tr_alt]) < 0.25]
            if not sel: sel = VARS[:5]
            Xs = d[sel].astype(str)
            m1 = Pipeline([('oh',OneHotEncoder(handle_unknown='ignore',drop='first',sparse_output=False)),
                           ('m',LogisticRegression(penalty='l2',C=0.05,max_iter=5000))]).fit(Xs.iloc[tr_alt], y[tr_alt])
            m2 = HistGradientBoostingClassifier(max_iter=200, learning_rate=0.05, max_leaf_nodes=7,
                    min_samples_leaf=10, random_state=0).fit(Xc[tr_alt], y[tr_alt])
            pr = (m1.predict_proba(Xs.iloc[te_alt])[:,1] + m2.predict_proba(Xc[te_alt])[:,1]) / 2
            p2[te_alt] = pr if np.isnan(p2[te_alt]).all() else np.nanmean([p2[te_alt], pr], axis=0)
    print(f'\n--- Perde 1 eşiği {ESIK1:.2f} (havuz {int(alt.sum())} hasta, '
          f'{int(y[alt].sum())} doğru / {int((~y[alt].astype(bool)).sum())} yanlış)')
    for E2 in [0.60, 0.70, 0.80, 0.85, 0.90]:
        m = alt & (p2 >= E2) & ~np.isnan(p2)
        if m.sum() < 20: continue
        ozet(m, f'  + perde 2 eşiği {E2:.2f}')

# ---------------------------------------------------------------- KLINIK KURAL
print('\n'+'='*104); print('ALTERNATİF — perde 2 olarak klinik kural'); print('='*104)
alt = elig & (p1 >= 0.70)
print(f'  havuz: {int(alt.sum())} hasta ({int(y[alt].sum())} doğru / {int((~y[alt].astype(bool)).sum())} yanlış)')
KURAL = {
 'nodal tutulum var': (d.Metastaz_Yeri_i != 'Yok').values,
 'Ki-67 düşük değil': (d.Ki_67_i != 'Düşük').values,
 'grade 2-3': d.HistolojikG_i.astype(str).str.contains('2|3', na=False).values,
 'nodal + Ki-67': ((d.Metastaz_Yeri_i != 'Yok') & (d.Ki_67_i != 'Düşük')).values,
 'nodal + grade 2-3': ((d.Metastaz_Yeri_i != 'Yok') &
                       d.HistolojikG_i.astype(str).str.contains('2|3', na=False)).values,
}
for ad, k in KURAL.items():
    ozet(alt & k, f'  + {ad}')
