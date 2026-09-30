"""ESIK EGRISI — bootstrap guven araliklariyla (nihai omurga, Kohort v15)."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
p = np.load('../../outputs/oof_final.npy')
d0 = pd.read_excel('../../data/Kohort_v17.xlsx')
d = d0[d0.Analiz_Kohortu == 'taban senaryo (cM0)'].reset_index(drop=True)
y = (d.RCB_Kategorize >= 2).astype(int).values
grp = d.Hasta_ID.values
elig = np.isin(d.Molekuler_i.values, ['Luminal A', 'Luminal B (HER2 Negatif)'])
s1 = (d.Molekuler_i == 'Luminal A').values
ESIKLER = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90]

def olc(idx, t):
    yy, pp, el, ss1 = y[idx], p[idx], elig[idx], s1[idx]
    m = el & (pp >= t); n = int(m.sum()); tp = int(yy[m].sum()); fp = n - tp
    s1d = int(yy[ss1].sum()); s1y = int((ss1 & (yy == 0)).sum())
    w = t / (1 - t)
    return dict(n=n, dogru=tp, yanlis=fp, PPV=tp/max(n,1), pay=100*n/len(yy),
                duyarlilik=tp/max(yy.sum(),1), net_dogru=tp-s1d, net_yanlis=fp-s1y,
                net_fayda=tp/len(yy)-(fp/len(yy))*w)

hastalar = pd.unique(grp); idx_by = {h: np.where(grp == h)[0] for h in hastalar}
rng = np.random.default_rng(42); B = 2000
boot = {t: [] for t in ESIKLER}
for b in range(B):
    sec = rng.choice(hastalar, size=len(hastalar), replace=True)
    idx = np.concatenate([idx_by[h] for h in sec])
    if len(np.unique(y[idx])) < 2: continue
    for t in ESIKLER: boot[t].append(olc(idx, t))

print(f'{"Threshold":>5} {"Flagged":>18} {"Concordant":>16} {"Discordant":>14} {"PPV":>20} '
      f'{"vs S1 concordant":>14} {"vs S1 discordant":>14} {"Net benefit":>10}')
sat = []
for t in ESIKLER:
    g = olc(np.arange(len(y)), t); bs = pd.DataFrame(boot[t])
    def ga(k, nd=0):
        lo, hi = np.nanpercentile(bs[k], [2.5, 97.5])
        f = f'%.{nd}f'
        return f'{f%g[k]} ({f%lo}–{f%hi})'
    print(f'{t:5.2f} {ga("n"):>18} {ga("dogru"):>16} {ga("yanlis"):>14} {ga("PPV",3):>20} '
          f'{ga("net_dogru"):>14} {ga("net_yanlis"):>14} {g["net_fayda"]:10.3f}')
    r = dict(esik=t, net_fayda=round(g['net_fayda'], 3),
             kohort_yuzdesi=round(g['pay'], 1), duyarlilik=round(g['duyarlilik'], 3))
    for k, nd in [('n',0),('dogru',0),('yanlis',0),('PPV',3),('net_dogru',0),('net_yanlis',0)]:
        lo, hi = np.nanpercentile(bs[k], [2.5, 97.5])
        r[k] = round(g[k], nd); r[k+'_alt'] = round(lo, nd); r[k+'_ust'] = round(hi, nd)
    sat.append(r)
pd.DataFrame(sat).to_csv('../../outputs/esik_egrisi.csv', index=False)
print('\nS1 karşılaştırıcı: 59 işaretlenen, 53 doğru, 6 yanlış, PPV 0,898')
