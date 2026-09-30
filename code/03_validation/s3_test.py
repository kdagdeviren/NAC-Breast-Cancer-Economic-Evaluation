"""S3 BIRLESIK STRATEJI — bootstrap, esik taramasi ve karar egrisi
S3 = S1 (Luminal A'nin tumu) ∪ S2 (model esigi gecenler)
"""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')

p = np.load('../../outputs/oof_final.npy')
d0 = pd.read_excel('../../data/Kohort_v17.xlsx')
d = d0[d0.Analiz_Kohortu == 'taban senaryo (cM0)'].reset_index(drop=True)
y = (d.RCB_Kategorize >= 2).astype(int).values
grp = d.Hasta_ID.values
elig = np.isin(d.Molekuler_i.values, ['Luminal A', 'Luminal B (HER2 Negatif)'])
la = (d.Molekuler_i == 'Luminal A').values
ESIK = 0.80

def olc(idx, esik=ESIK):
    yy, pp = y[idx], p[idx]
    el, l = elig[idx], la[idx]
    top = int(yy[el & (yy == 1)].sum())
    out = {}
    for ad, m in [('S1', l), ('S2', el & (pp >= esik)), ('S3', l | (el & (pp >= esik)))]:
        n = int(m.sum()); tp = int(yy[m].sum())
        out[ad] = dict(n=n, dogru=tp, yanlis=n-tp, PPV=tp/max(n,1), kacirilan=top-tp,
                       duyarlilik=tp/max(top,1))
    return out

g = olc(np.arange(len(y)))
print('=' * 92); print('TABAN — üç strateji'); print('=' * 92)
print(f'{"":6s} {"işaret":>7s} {"doğru":>7s} {"yanlış":>7s} {"kaçırılan":>10s} {"PPV":>7s} {"duyarlılık":>11s}')
for k in ['S1','S2','S3']:
    v = g[k]
    print(f'{k:6s} {v["n"]:7d} {v["dogru"]:7d} {v["yanlis"]:7d} {v["kacirilan"]:10d} '
          f'{v["PPV"]:7.3f} {v["duyarlilik"]:11.3f}')

# ---- bootstrap
hastalar = pd.unique(grp); idx_by = {h: np.where(grp == h)[0] for h in hastalar}
rng = np.random.default_rng(2026); B = 2000
kayit = []
for b in range(B):
    sec = rng.choice(hastalar, size=len(hastalar), replace=True)
    idx = np.concatenate([idx_by[h] for h in sec])
    if len(np.unique(y[idx])) < 2: continue
    o = olc(idx)
    kayit.append({f'{k}_{m}': o[k][m] for k in ['S1','S2','S3'] for m in
                  ['n','dogru','yanlis','PPV','kacirilan','duyarlilik']})
bs = pd.DataFrame(kayit)

print('\n' + '=' * 92); print(f'BOOTSTRAP %95 GA (B={len(bs)})'); print('=' * 92)
for k in ['S1','S2','S3']:
    print(f'\n  {k}')
    for m, nd in [('n',0),('dogru',0),('yanlis',0),('kacirilan',0),('PPV',3),('duyarlilik',3)]:
        v = bs[f'{k}_{m}']
        f = f'%.{nd}f'
        print(f'    {m:12s} {f % g[k][m]:>8s}  ({f % np.percentile(v,2.5)} – {f % np.percentile(v,97.5)})')

print('\n' + '=' * 92); print('S3 − S2 FARKI (asıl soru)'); print('=' * 92)
for m, ad in [('dogru','ek doğru'), ('yanlis','ek yanlış'), ('kacirilan','kaçırılan değişimi'),
              ('n','ek işaretlenen')]:
    f = bs[f'S3_{m}'] - bs[f'S2_{m}']
    goz = g['S3'][m] - g['S2'][m]
    print(f'  {ad:22s} {goz:+4d}  (%95 GA {np.percentile(f,2.5):+.0f} – {np.percentile(f,97.5):+.0f})'
          f'   sıfırı içeriyor mu: {"EVET" if np.percentile(f,2.5)<=0<=np.percentile(f,97.5) else "HAYIR"}')

print('\n' + '=' * 92); print('S3 − S1 FARKI'); print('=' * 92)
for m, ad in [('dogru','ek doğru'), ('yanlis','ek yanlış')]:
    f = bs[f'S3_{m}'] - bs[f'S1_{m}']
    goz = g['S3'][m] - g['S1'][m]
    print(f'  {ad:22s} {goz:+4d}  (%95 GA {np.percentile(f,2.5):+.0f} – {np.percentile(f,97.5):+.0f})')

# ---- esik taramasi
print('\n' + '=' * 92); print('S3 EŞİK TARAMASI'); print('=' * 92)
print(f'{"eşik":>6s} {"işaret":>7s} {"doğru":>7s} {"yanlış":>7s} {"kaçırılan":>10s} {"PPV":>7s} {"S1 üzerine":>22s}')
for t in [0.70, 0.75, 0.80, 0.85, 0.90, 1.01]:
    o = olc(np.arange(len(y)), t)['S3']
    ad = 'sadece S1' if t > 1 else f'{t:.2f}'
    print(f'{ad:>6s} {o["n"]:7d} {o["dogru"]:7d} {o["yanlis"]:7d} {o["kacirilan"]:10d} '
          f'{o["PPV"]:7.3f} {o["dogru"]-g["S1"]["dogru"]:+8d} doğru {o["yanlis"]-g["S1"]["yanlis"]:+5d} yanlış')

# ---- karar egrisi
print('\n' + '=' * 92); print('KARAR EĞRİSİ — net fayda'); print('=' * 92)
N = len(y)
print(f'{"pt":>6s} {"S3":>9s} {"S2":>9s} {"S1":>9s} {"tümü":>9s}')
for pt in [0.5,0.6,0.7,0.75,0.8,0.85,0.9]:
    w = pt/(1-pt)
    def nf(m):
        tp = int(y[m].sum()); fp = int(m.sum())-tp
        return tp/N - (fp/N)*w
    print(f'{pt:6.2f} {nf(la|(elig&(p>=pt))):9.3f} {nf(elig&(p>=pt)):9.3f} {nf(la):9.3f} {nf(elig):9.3f}')
