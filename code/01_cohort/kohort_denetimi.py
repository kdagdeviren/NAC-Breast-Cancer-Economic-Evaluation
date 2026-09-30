"""KOHORT DENETIMI — yapisal, klinik ve capraz tutarlilik kontrolleri."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
pd.set_option('display.width', 200)

d = pd.read_excel('../../data/Kagan_TEZ.xlsx', sheet_name=0)
B = []
def bul(kod, baslik, detay=''):
    B.append((kod, baslik, detay)); print(f'[{kod}] {baslik}')
    if detay: print(detay)

print('=' * 90); print('1. YAPISAL'); print('=' * 90)
print(f'satir {len(d)} | kolon {d.shape[1]} | tekil Hasta_ID {d.Hasta_ID.nunique()}')
dup = d.Hasta_ID[d.Hasta_ID.duplicated()].tolist()
if dup: bul('Y1', f'tekrarlayan Hasta_ID: {dup}')
# ayni hastanin iki memesi?
kok = d.Hasta_ID.astype(str).str[:-1]
cift = kok[kok.duplicated()].tolist()
if cift: bul('Y2', f'ayni numara iki memede: {sorted(set(cift))}')
tam_bos = [c for c in d.columns if d[c].isna().all()]
if tam_bos: bul('Y3', f'tamamen bos kolon: {tam_bos}')
tekil = [c for c in d.columns if d[c].nunique(dropna=True) == 1]
if tekil: bul('Y4', f'tek degerli (bilgi tasimayan) kolon: {tekil}')
print('\nTani yili dagilimi:'); print(d['Tanı Yılı'].astype(str).value_counts().sort_index().to_string())

print('\n' + '=' * 90); print('2. RCB IC TUTARLILIGI'); print('=' * 90)
print(pd.crosstab(d.RCB_Sinif, d.RCB_Kategorize, dropna=False))
# Symmans esikleri
esik = {0: (0, 0), 1: (0, 1.36), 2: (1.36, 3.28), 3: (3.28, 99)}
kot = []
for k, (lo, hi) in esik.items():
    m = d.RCB_Kategorize == k
    if k == 0:
        bad = d[m & (d.RCB_Skor.fillna(0) > 0)]
    else:
        bad = d[m & ~d.RCB_Skor.between(lo, hi)]
    if len(bad): kot.append((k, bad.Hasta_ID.tolist(), bad.RCB_Skor.tolist()))
if kot:
    bul('R1', 'RCB skoru sinif esigiyle uyusmayan olgular',
        '\n'.join(f'  RCB-{k}: {h} skor={np.round(s,3).tolist()}' for k, h, s in kot))
else:
    print('RCB_Skor -> RCB_Kategorize esikleri tutarli')
b = d[(d.RCB_Kategorize == 0) & (d['%CA'].fillna(0) != 0)]
if len(b): bul('R2', f'RCB-0 oldugu halde %CA != 0: {b.Hasta_ID.tolist()} '
                     f'(%CA={b["%CA"].tolist()})')
b = d[(d.RCB_Kategorize == 0) & ((d.N.fillna(0) != 0) | (d.M.fillna(0) != 0))]
if len(b): bul('R3', f'RCB-0 oldugu halde N veya M != 0: {b.Hasta_ID.tolist()}')
b = d[(d.RCB_Kategorize > 0) & (d['%CA'].fillna(-1) == 0) & (d.N.fillna(0) == 0)]
if len(b): bul('R4', f'RCB>0 ama %CA=0 ve N=0: {b.Hasta_ID.tolist()}')

print('\n' + '=' * 90); print('3. MOLEKULER ALT TIP vs ER/PR/HER2'); print('=' * 90)
print(pd.crosstab([d.ER_i, d.PR_i, d.HER2_i], d.Molekuler_i, dropna=False).to_string())

print('\n' + '=' * 90); print('4. NOTTINGHAM GRADE = TUBUL + NUKLEER + MITOTIK'); print('=' * 90)
def skor(s):
    return pd.to_numeric(s.astype(str).str.extract(r'(\d)')[0], errors='coerce')
tp, nk, mt = skor(d.Tubul_i), skor(d.Nukleer_i), skor(d.Mitotik_i)
tot = tp + nk + mt
bek = pd.cut(tot, [2, 5, 7, 9], labels=['1', '2', '3'])
gz = d.HistolojikG_i.astype(str).str.extract(r'(\d)')[0]
uy = pd.crosstab(bek, gz, dropna=False)
print(uy)
ce = d[(bek.astype(str) != gz) & bek.notna() & gz.notna()]
if len(ce):
    bul('G1', f'Nottingham bilesen toplami ile bildirilen grade uyusmuyor: {len(ce)} olgu',
        '  ornek: ' + ', '.join(f'{r.Hasta_ID}(t{int(tp[i])}+n{int(nk[i])}+m{int(mt[i])}'
                                f'={int(tot[i])} -> beklenen G{bek[i]}, kayitli G{gz[i]})'
                                for i, r in list(ce.iterrows())[:8]))

print('\n' + '=' * 90); print('5. MENOPOZ vs YAS'); print('=' * 90)
print(pd.crosstab(d.Yas_Grubu, d.Menapoz_i, dropna=False))
g = d[(d.Menapoz_i.astype(str).str.contains('Post', na=False)) & (d.Tani_Yasi < 40)]
if len(g): bul('M1', f'40 yas alti postmenopozal: {g.Hasta_ID.tolist()} '
                     f'yas={g.Tani_Yasi.tolist()}')
g = d[(d.Menapoz_i.astype(str).str.contains('Pre', na=False)) & (d.Tani_Yasi > 58)]
if len(g): bul('M2', f'58 yas ustu premenopozal: {g.Hasta_ID.tolist()} '
                     f'yas={g.Tani_Yasi.tolist()}')

print('\n' + '=' * 90); print('6. VKI = KILO / BOY^2'); print('=' * 90)
boy = pd.to_numeric(d.Boy_imp, errors='coerce'); kilo = pd.to_numeric(d.Kilo_imp, errors='coerce')
boy_m = np.where(boy > 3, boy / 100, boy)
vki_h = kilo / boy_m ** 2
fark = (vki_h - pd.to_numeric(d.VKI, errors='coerce')).abs()
kt = d[fark > 1]
print(f'VKI yeniden hesabi ile |fark| > 1 olan: {len(kt)}')
if len(kt): bul('V1', f'VKI tutarsizligi: {kt.Hasta_ID.head(10).tolist()}',
                f'  ornek fark: {np.round(fark[fark>1].head(5).values,2).tolist()}')
for k, lo, hi in [('Boy(cm)', 130, 200), ('Kilo(kg)', 35, 180), ('Yas', 18, 90)]:
    s = {'Boy(cm)': boy, 'Kilo(kg)': kilo, 'Yas': d.Tani_Yasi}[k]
    ay = d[(s < lo) | (s > hi)]
    if len(ay): bul('V2', f'{k} sinir disi: {ay.Hasta_ID.tolist()} -> {s[ay.index].tolist()}')

print('\n' + '=' * 90); print('7. BIYOKIMYA — deger ile durum etiketi'); print('=' * 90)
cif = [('ALP', 'ALP_Değer', 'ALP_Durumu'), ('ALT', 'ALT_Değer', 'ALT_Durumu'),
       ('AST', 'AST_Değer', 'AST_Durumu'), ('BUN', 'BUN_Değer', 'BUN_Durumu'),
       ('CA15-3', 'CA15-3_Değer', 'CA15-3_Durumu'), ('CEA', 'CEA_Değer', 'CEA_Durumu'),
       ('CRP', 'CRP_Değer', 'CRP_Durumu'), ('GGT', 'GGT_Değer', 'GGT_Durumu'),
       ('Glukoz', 'Glukoz_Değer', 'Glukoz_Durumu'), ('HbA1c', 'HbA1c_Değer', 'HbA1c_Durumu'),
       ('Kreatinin', 'Kreatinin_Değer', 'Kreatinin_Durumu'), ('LDH', 'LDH_Değer', 'LDH_Durumu'),
       ('TSH', 'TSH_Değer', 'TSH_Durumu'), ('eGFR', 'e-GFR (CKD-EPI)_Değer', 'e-GFR_Durumu')]
for ad, dc, sc in cif:
    if dc not in d.columns: continue
    v = pd.to_numeric(d[dc], errors='coerce')
    norm = d[sc].astype(str).str.contains('Normal', na=False)
    if v.notna().sum() < 10: continue
    ov = v[norm].agg(['min', 'max']); yv = v[~norm & v.notna()].agg(['min', 'max'])
    print(f'  {ad:10s} normal aralik {ov.min():8.2f}-{ov.max():8.2f} | '
          f'anormal {yv.min():8.2f}-{yv.max():8.2f} | eksik deger {int(v.isna().sum())}')
dm = d.DM.astype(str).str.contains('Var|Evet|1', na=False)
h = pd.to_numeric(d.HbA1c_Değer, errors='coerce')
g = d[(~dm) & (h > 6.5)]
if len(g): bul('B1', f'DM kaydi yok ama HbA1c > 6.5: {len(g)} olgu '
                     f'({g.Hasta_ID.head(8).tolist()})')
g = d[dm & (h < 5.7) & h.notna()]
if len(g): bul('B2', f'DM var ama HbA1c < 5.7: {len(g)} olgu')

print('\n' + '=' * 90); print('8. RADYOLOJI IC TUTARLILIGI'); print('=' * 90)
print(pd.crosstab(d['BI-RADS_i'], d.Lezyon_Turu_i, dropna=False).to_string())
g = d[(d['BI-RADS_i'].astype(str).isin(['0', '1', '2'])) & (d.RCB_Kategorize.notna())]
if len(g): bul('D1', f'kanser tanisi kesin oldugu halde BI-RADS 0/1/2: {len(g)} olgu '
                     f'({g.Hasta_ID.head(8).tolist()})')
g = d[(d['2_Yildir_Stabil_i'].astype(str).str.contains('Var|Evet|Stabil', na=False))]
print(f'2 yildir stabil = var: {len(g)} olgu')
kal = d.Kalsifikasyon_Morfolojisi_i.astype(str); dag = d.Kalsifikasyon_Dagilimi_i.astype(str)
g = d[(kal.str.contains('Yok', na=False)) & (~dag.str.contains('Yok|Değerlend', na=False))]
if len(g): bul('D2', f'kalsifikasyon morfolojisi yok ama dagilimi tanimli: {len(g)} olgu')
sek = d.Kitle_Sekli_i.astype(str); kon = d.Kitle_Konturu_i.astype(str)
g = d[(sek.str.contains('Yok', na=False)) & (~kon.str.contains('Yok|Değerlend', na=False))]
if len(g): bul('D3', f'kitle sekli yok ama konturu tanimli: {len(g)} olgu')

print('\n' + '=' * 90); print('9. IMPUTASYON: ham deger varken degisen hucreler'); print('=' * 90)
ciftler = [(c.replace('_i', ''), c) for c in d.columns if c.endswith('_i')]
tab = []
for ham, imp in ciftler:
    if ham not in d.columns: continue
    m = d[ham].notna() & (~d[ham].astype(str).str.contains('Bilinmiyor|Değerlendirilemedi', na=False))
    if m.sum() < 20: continue
    fark = (d.loc[m, ham].astype(str) != d.loc[m, imp].astype(str)).sum()
    if fark: tab.append((ham, int(m.sum()), int(fark), round(100 * fark / m.sum(), 1)))
if tab:
    t = pd.DataFrame(tab, columns=['degisken', 'ham_dolu', 'degisen', '%'])
    bul('I1', 'imputasyon dolu hucreyi degistirmis gorunuyor',
        t.sort_values('%', ascending=False).to_string(index=False))
else:
    print('dolu hucrelerde ham ile imputasyonlu surum ayni')

print('\n' + '=' * 90); print('10. i48-i64 KONTAMINASYONU (yeniden dogrulama)'); print('=' * 90)
sat = np.zeros(len(d), dtype=bool); ozet = []
for n in range(48, 65):
    a, b_ = f'i{n}', f'i{n}_e'
    if a in d.columns and b_ in d.columns:
        fk = d[a].astype(str) != d[b_].astype(str)
        sat |= fk.values
        if fk.sum(): ozet.append((a, int(fk.sum())))
bul('K1', f'i48-i64 ile i*_e ayrisan satir sayisi: {int(sat.sum())}/{len(d)}',
    '  kolon basina: ' + ', '.join(f'{k}={v}' for k, v in ozet))
print('  ayrisan satirlarin RCB dagilimi:')
print(pd.crosstab(sat, d.RCB_Kategorize).to_string())

print('\n' + '=' * 90); print('11. EKSIK VERI YOGUNLUGU (ham degiskenler)'); print('=' * 90)
ham = [c for c in d.columns if not c.startswith('i') and not c.endswith('_i')
       and not c.endswith('_e') and d[c].dtype == object]
ek = []
for c in ham:
    s = d[c].astype(str)
    n = s.str.contains('Bilinmiyor|Değerlendirilemedi', na=False).sum() + d[c].isna().sum()
    if n: ek.append((c, int(n), round(100 * n / len(d), 1)))
t = pd.DataFrame(ek, columns=['degisken', 'eksik', '%']).sort_values('%', ascending=False)
print(t.head(25).to_string(index=False))

print('\n' + '=' * 90)
print(f'TOPLAM ISARETLENEN BULGU: {len(B)}')
for k, h, _ in B: print(f'  [{k}] {h}')
