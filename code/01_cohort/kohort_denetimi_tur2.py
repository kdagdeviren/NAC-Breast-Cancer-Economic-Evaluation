"""KOHORT DENETIMI — 2. TUR: ilk turda kapsanmayan tum kolonlar."""
import numpy as np, pandas as pd, warnings
warnings.filterwarnings('ignore')
pd.set_option('display.width', 220)
d = pd.read_excel('../../data/Kagan_TEZ.xlsx', sheet_name=0)
B = []
def bul(k, h, det=''):
    B.append((k, h)); print(f'[{k}] {h}')
    if det: print(det)

print('=' * 92); print('12. i1-i47 KODLAMALARI — kaynak kolonla birebir mi'); print('=' * 92)
esles = {'i1':'Histolojik_Yeni','i2':'ER_i','i3':'PR_i','i4':'HER2_i','i5':'Molekuler_i',
 'i6':'Ki_67_i','i7':'Tubul_i','i8':'Nukleer_i','i9':'Mitotik_i','i10':'HistolojikG_i',
 'i11':'Ecaderin_i','i12':'TIL_i','i13':'Metasdaz_Durumu_i','i14':'Metastaz_Yeri_i',
 'i15':'NAC_Once_Evre_i','i46':'Rejim_i','i47':'Kur_Yogunluk_i','i16':'Hangi_Meme',
 'i17':'VKI_Sinifi_17','i18':'Yas_Grubu','i19':'Kan_Grup_i','i20':'Menapoz_i',
 'i45':'Gunesten_Yararlanma','i21':'HT','i22':'DM','i23':'KOAH','i24':'Sigara',
 'i25':'Ailede_MemeCA','i26':'Tiroid','i27':'Retinopati','i28':'Nöropati',
 'i29':'Osteoporoz','i30':'Depresyon','i31':'ALP_Durumu','i32':'ALT_Durumu',
 'i33':'AST_Durumu','i34':'BUN_Durumu','i35':'CA15-3_Durumu','i36':'CEA_Durumu',
 'i37':'CRP_Durumu','i38':'GGT_Durumu','i39':'Glukoz_Durumu','i40':'HbA1c_Durumu',
 'i41':'Kreatinin_Durumu','i42':'LDH_Durumu','i43':'TSH_Durumu','i44':'e-GFR_Durumu'}
kirik = []
for kod, kay in esles.items():
    if kod not in d.columns or kay not in d.columns: continue
    t = pd.crosstab(d[kay].astype(str), d[kod])
    # her kaynak kategori tek bir koda gitmeli ve her kod tek kaynaga
    ok = (t > 0).sum(axis=1).max() == 1 and (t > 0).sum(axis=0).max() == 1
    if not ok:
        kirik.append(kod)
        bul('E1', f'{kod} <-> {kay} birebir DEGIL', t.to_string())
if not kirik: print('i1-i47 kodlamalarinin tamami kaynak kolonla birebir')

print('\n' + '=' * 92); print('13. HANGI MEME vs HASTA_ID SONEKI'); print('=' * 92)
sonek = d.Hasta_ID.astype(str).str.strip().str[-1]
print(pd.crosstab(d.Hangi_Meme, sonek).to_string())
uy = d[((sonek == 'R') & (d.Hangi_Meme.astype(str).str.contains('Sol', na=False))) |
       ((sonek == 'L') & (d.Hangi_Meme.astype(str).str.contains('Sağ', na=False)))]
if len(uy): bul('H1', f'Hasta_ID soneki ile Hangi_Meme celisiyor: {len(uy)} olgu '
                      f'{uy.Hasta_ID.tolist()[:15]}')

print('\n' + '=' * 92); print('14. ALLRED SKORU vs ER/PR ETIKETI'); print('=' * 92)
for rec, allred, etik in [('ER', 'ER_Allred_imputed', 'ER_i'), ('PR', 'PR_Allred_imputed', 'PR_i')]:
    if allred not in d.columns: continue
    a = pd.to_numeric(d[allred], errors='coerce')
    print(f'--- {rec} (Allred 0-8; >=3 pozitif kabul edilir)')
    print(d.groupby(d[etik].astype(str))[allred].apply(
        lambda s: f'n={s.notna().sum()} min={pd.to_numeric(s,errors="coerce").min():.0f} '
                  f'max={pd.to_numeric(s,errors="coerce").max():.0f}').to_string())
    kot = d[(d[etik].astype(str) == 'Negatif') & (a >= 3)]
    if len(kot): bul('A1', f'{rec} negatif etiketli ama Allred >=3: {len(kot)} olgu '
                           f'{kot.Hasta_ID.tolist()[:10]}')
    kot = d[(d[etik].astype(str).str.contains('Pozitif', na=False)) & (a < 3)]
    if len(kot): bul('A2', f'{rec} pozitif etiketli ama Allred <3: {len(kot)} olgu '
                           f'{kot.Hasta_ID.tolist()[:10]}')
    kot = d[a > 8]
    if len(kot): bul('A3', f'{rec} Allred > 8 (olanaksiz): {kot.Hasta_ID.tolist()}')

print('\n' + '=' * 92); print('15. Ki-67 SAYISAL vs KATEGORI ve LUMINAL A/B AYRIMI'); print('=' * 92)
k = pd.to_numeric(d.Ki_67, errors='coerce')
print(d.groupby(d.Ki_67_i.astype(str)).apply(
    lambda g: f'n={len(g)} ham min={pd.to_numeric(g.Ki_67,errors="coerce").min()} '
              f'max={pd.to_numeric(g.Ki_67,errors="coerce").max()}').to_string())
la = d.Molekuler.astype(str) == 'Luminal A'
lb = d.Molekuler.astype(str) == 'Luminal B (HER2 Negatif)'
print(f'Luminal A  Ki-67: n={k[la].notna().sum()} medyan={k[la].median()} max={k[la].max()}')
print(f'Luminal B- Ki-67: n={k[lb].notna().sum()} medyan={k[lb].median()} min={k[lb].min()}')
kot = d[la & (k > 30)]
if len(kot): bul('K2', f'Luminal A ama Ki-67 > %30: {len(kot)} olgu '
                       f'{kot.Hasta_ID.tolist()[:10]} (Ki-67={k[kot.index].tolist()[:10]})')
kot = d[lb & (k < 10)]
if len(kot): bul('K3', f'Luminal B (HER2-) ama Ki-67 < %10: {len(kot)} olgu '
                       f'{kot.Hasta_ID.tolist()[:10]}')

print('\n' + '=' * 92); print('16. E-KADERIN vs HISTOLOJIK TIP'); print('=' * 92)
print(pd.crosstab(d.Histolojik_Yeni, d.Ecaderin_i, dropna=False).to_string())
lob = d.Histolojik_Yeni.astype(str).str.contains('Lobüler|Lobuler', na=False)
kot = d[lob & (d.Ecaderin_i.astype(str).str.contains('Pozitif', na=False))]
if len(kot): bul('EC1', f'Lobuler karsinom ama E-kaderin POZITIF: {len(kot)} olgu '
                        f'{kot.Hasta_ID.tolist()[:12]}')
duk = d.Histolojik_Yeni.astype(str).str.contains('Duktal|Invaziv Duktal', na=False)
kot = d[duk & (d.Ecaderin_i.astype(str).str.contains('Negatif', na=False))]
if len(kot): bul('EC2', f'Duktal karsinom ama E-kaderin NEGATIF: {len(kot)} olgu '
                        f'{kot.Hasta_ID.tolist()[:12]}')

print('\n' + '=' * 92); print('17. KOMORBIDITE BLOGU'); print('=' * 92)
kom = ['HT','DM','KOAH','Sigara','Ailede_MemeCA','Tiroid','Retinopati','Nöropati',
       'Osteoporoz','Depresyon']
for c in kom:
    vc = d[c].value_counts(dropna=False).to_dict()
    print(f'  {c:16s} {vc}')
def var(c): return d[c].astype(str).str.contains('Var|Evet|1', na=False)
kot = d[var('Retinopati') & ~var('DM')]
if len(kot): bul('KM1', f'DM yok ama retinopati var: {len(kot)} olgu {kot.Hasta_ID.tolist()[:10]}')
kot = d[var('Osteoporoz') & (d.Menapoz_i.astype(str) == 'Yok')]
if len(kot): bul('KM2', f'premenopozal ama osteoporoz: {len(kot)} olgu {kot.Hasta_ID.tolist()[:10]}')
tsh_an = ~d.TSH_Durumu.astype(str).str.contains('Normal', na=False)
kot = d[var('Tiroid') & ~tsh_an]
print(f'  tiroid hastaligi var + TSH normal: {int((var("Tiroid") & ~tsh_an).sum())} '
      f'(tedaviyle normalize olabilir, bilgi amacli)')
kot = d[~var('Tiroid') & tsh_an]
if len(kot): bul('KM3', f'tiroid hastaligi kaydi yok ama TSH anormal: {len(kot)} olgu')

print('\n' + '=' * 92); print('18. RCB BILESENLERI (%CA, %CIS, N, M, d1xd2)'); print('=' * 92)
for c in ['d1xd2', '%CA', '%CIS', 'N', 'M', 'RCB_Skor']:
    v = pd.to_numeric(d[c], errors='coerce')
    print(f'  {c:10s} dolu={v.notna().sum():3d} min={v.min()} max={v.max()} '
          f'sifir={int((v==0).sum())}')
if pd.to_numeric(d.d1xd2, errors='coerce').notna().sum() == 0:
    bul('RB1', 'd1xd2 (tumor boyutu) kolonu tamamen bos — RCB formulunun bir bileseni yok')
n = pd.to_numeric(d.N, errors='coerce'); m = pd.to_numeric(d.M, errors='coerce')
kot = d[(n.fillna(0) == 0) & (m.fillna(0) > 0)]
if len(kot): bul('RB2', f'N=0 ama M>0 (tutulan nod yok ama metastaz capi var): {len(kot)} olgu '
                        f'{kot.Hasta_ID.tolist()[:10]}')
kot = d[(n.fillna(0) > 0) & (m.fillna(0) == 0)]
if len(kot): bul('RB3', f'N>0 ama M=0: {len(kot)} olgu (nod tutulmus ama metastaz capi sifir)')
print('  RCB_ML kolonu:', d.RCB_ML.value_counts(dropna=False).to_dict())

print('\n' + '=' * 92); print('19. REJIM ve KUR YOGUNLUGU (modelde yok, kayit icin)'); print('=' * 92)
print(d.Rejim.value_counts(dropna=False).to_string())
print(); print(d.Kur_Yogunluk.value_counts(dropna=False).to_string())
print(); print(pd.crosstab(d.Kur_Yogunluk_i, d.RCB_Kategorize).to_string())

print('\n' + '=' * 92); print('20. DEMOGRAFI KALANI'); print('=' * 92)
for c in ['Kan_Grup', 'Gunesten_Yararlanma', 'Yas_Grubu', 'VKI_Sinifi_17', 'Hangi_Meme']:
    print(f'  {c:22s} {d[c].value_counts(dropna=False).to_dict()}')
b_raw = pd.to_numeric(d.Boy, errors='coerce'); k_raw = pd.to_numeric(d.Kilo, errors='coerce')
print(f'  Boy ham eksik={int(b_raw.isna().sum())} | Kilo ham eksik={int(k_raw.isna().sum())} '
      f'(imputasyonla dolduruldu)')

print('\n' + '=' * 92); print('21. HISTOLOJIK TIP ve TIL'); print('=' * 92)
print(d.Histolojik_Yeni.value_counts(dropna=False).to_string())
print(); print(pd.crosstab(d.TIL_i, d.Molekuler_i, dropna=False).to_string())

print('\n' + '=' * 92); print(f'2. TUR TOPLAM BULGU: {len(B)}')
for k, h in B: print(f'  [{k}] {h}')
