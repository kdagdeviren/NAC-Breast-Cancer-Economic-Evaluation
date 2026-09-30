"""KOHORT DUZELTME HATTI — kararlar sirayla uygulanir, her adim gunluge yazilir.
Girdi : Kagan_TEZ.xlsx (ham)
Cikti : Kohort_v2.xlsx + Duzeltme_Gunlugu.md
"""
import pandas as pd, numpy as np, warnings
warnings.filterwarnings('ignore')

SRC = '../../data/Kagan_TEZ.xlsx'
OUT = '../../outputs/'
SURUM = 'v17'   # her kabul edilen karar setinden sonra artirilir
d = pd.read_excel(SRC, sheet_name=0)
LOG = []
def log(madde, ne, kac, karar):
    LOG.append((madde, ne, kac, karar))
    print(f'[{madde}] {ne} — {kac} kayit')

GRUP = {'NAC_Once_Evre': ['NAC_Once_Evre', 'NAC_Once_Evre_15', 'NAC_Once_Evre_i'],
        'Metastaz_Yeri': ['Metastaz_Yeri', 'Metastaz_Yeri_14', 'Metastaz_Yeri_i'],
        'Metasdaz_Durumu': ['Metasdaz_Durumu', 'Metasdaz_Durumu_13', 'Metasdaz_Durumu_i']}

def duzelt(anahtar, hastalar, etiket, eski_id):
    """Ham / kodlu / imputasyonlu tum surumleri birlikte gunceller."""
    ref = eski_id.isin(hastalar)
    ornek = d[anahtar + '_i'] == etiket
    for c in [c for c in GRUP[anahtar] if c in d.columns]:
        if d[c].dtype == object:
            d.loc[ref, c] = etiket
        else:
            kod = d.loc[ornek, c].mode()
            if len(kod): d.loc[ref, c] = kod.iloc[0]

eski = d['Hasta_ID'].astype(str).str.strip()

# ---------------------------------------------------------------- ONCEKI TURLAR
duzelt('NAC_Once_Evre', ['111R', '143R', '294R'], 'Evre IV', eski)
log('0a', 'uzak metastaz var, evre yanlis kodlanmis -> Evre IV', 3,
    '111R, 143R, 294R')
duzelt('Metastaz_Yeri', ['309L'], 'Lenf Nodu', eski)
duzelt('Metasdaz_Durumu', ['309L'], 'Var', eski)
log('0b', 'metastaz yeri lenf nodu olarak duzeltildi; evre, NAK oncesi raporlarin '
         'yeniden incelenmesiyle Evre IV olarak teyit edildi (degisiklik yok)',
    1, '309L')
duzelt('Metasdaz_Durumu', ['186R', '213R', '218L', '210L'], 'Yok', eski)
duzelt('Metastaz_Yeri', ['210L'], 'Yok', eski)
log('0c', 'metastaz yok olarak duzeltildi (Durumu/Yeri celiskisi)', 4,
    '186R, 213R, 218L, 210L')


# ---------------------------------------------------------------- MADDE 2a
# Molekuler_i, ham Molekuler ile uyumsuz 3 kayitta duzeltildi (i5 zaten ham degeri tasiyor)
M2A = {'036R': 'Luminal B (HER2 Pozitif)', '077R': 'Luminal B (HER2 Negatif)',
       '090R': 'Luminal B (HER2 Negatif)'}
for hid, etiket in M2A.items():
    d.loc[eski == hid, 'Molekuler_i'] = etiket
kalan = int(((~(d.Molekuler.isna() | d.Molekuler.astype(str).str.contains('Bilinmiyor', na=False)))
             & (d.Molekuler.astype(str) != d.Molekuler_i.astype(str))).sum())
log('2a', f'Molekuler_i ham degere gore duzeltildi (kalan uyumsuzluk: {kalan})', 3,
    '036R, 077R, 090R')


# ---------------------------------------------------------------- MADDE 2b (A ve B)
# Ilke: her kayitta ham (dolu) kolon esas alinir, imputasyon ciktisi ona uydurulur.
# A) Alt tip ham veri, HER2 bos idi -> HER2 alt tipten turetildi
M2B_A = ['018R', '084R', '163L', '262R', '278L']
for c in ['HER2_i']:
    d.loc[eski.isin(M2B_A), c] = 'Pozitif'
if 'i4' in d.columns:
    kod = d.loc[d.HER2_i == 'Pozitif', 'i4'].mode()
    if len(kod): d.loc[eski.isin(M2B_A), 'i4'] = kod.iloc[0]
log('2b-A', 'alt tip "Luminal B (HER2 Pozitif)" gercek, HER2 bostu -> HER2_i = Pozitif',
    5, ', '.join(M2B_A))

# B) HER2 ham veri (Negatif), alt tip bos idi -> alt tip HER2'ye gore duzeltildi
M2B_B = ['109R', '122L']
d.loc[eski.isin(M2B_B), 'Molekuler_i'] = 'Luminal B (HER2 Negatif)'
if 'i5' in d.columns:
    kod = d.loc[d.Molekuler_i == 'Luminal B (HER2 Negatif)', 'i5'].mode()
    if len(kod): d.loc[eski.isin(M2B_B), 'i5'] = kod.iloc[0]
log('2b-B', 'HER2 "Negatif" gercek, alt tip bostu -> Luminal B (HER2 Negatif)',
    2, ', '.join(M2B_B))


# ---------------------------------------------------------------- MADDE 2b (C)
# HER2 IHK 2+ ekvokal: ISH dogrulanmadigi icin pozitif/negatif atanmadi.
# Ayri alt tip kategorisi olarak tutulur ve karar kuralinin uygunluk kapsamina alinmaz.
ekv = d.HER2_i.astype(str) == 'Ekvokal'
d['Molekuler_onceki'] = d.Molekuler_i
d.loc[ekv, 'Molekuler_i'] = 'HER2 ekvokal (ISH doğrulanmamış)'
if 'i5' in d.columns:
    d.loc[ekv, 'i5'] = 7
log('2b-C', 'HER2 ekvokal olgular ayri alt tip kategorisine alindi (i5 = 7), '
            'uygunluk kapsami disinda', int(ekv.sum()),
    'ISH doğrulanmadan pozitif/negatif atanmadi')


# ---------------------------------------------------------------- MADDE 2c (Grup 2)
# 2c-A: HER2 ham verisi "HER2-düşük" (gercek), alt tip bostu -> ER-/PR-/HER2-low = triple negatif
M2C_A = ['205L', '207L', '233R']
d.loc[eski.isin(M2C_A), 'Molekuler_i'] = 'Triple Negatif'
if 'i5' in d.columns:
    kod = d.loc[d.Molekuler_i == 'Triple Negatif', 'i5'].mode()
    if len(kod): d.loc[eski.isin(M2C_A), 'i5'] = kod.iloc[0]
log('2c-A', 'HER2 ham "HER2-düşük" gercek, alt tip bostu -> Triple Negatif (HER2-low)',
    3, ', '.join(M2C_A))

# 2c-B: ER- / PR "zayif pozitif" / HER2+ olgular PR Allred skoruna gore ayristirildi (esik >=3)
M2C_B = ['316L']                       # Allred 5 -> gercek PR pozitif
d.loc[eski.isin(M2C_B), 'Molekuler_i'] = 'Luminal B (HER2 Pozitif)'
if 'i5' in d.columns:
    kod = d.loc[d.Molekuler_i == 'Luminal B (HER2 Pozitif)', 'i5'].mode()
    if len(kod): d.loc[eski.isin(M2C_B), 'i5'] = kod.iloc[0]
log('2c-B', 'kaynak kayittaki PR Allred skoruna gore (114R: 2, 257R: 0 -> etiket "zayif '
            'pozitif" olsa da negatif; 316L: 5 -> gercek pozitif) alt tip belirlendi',
    1, '316L -> Luminal B (HER2 Pozitif)')


# ---------------------------------------------------------------- MADDE 2d (Grup 3)
# Ham veride Luminal A/B ayrimi istisnasiz Ki-67'ye dayaniyor (Dusuk -> A, Orta/Yuksek -> B;
# ham Molekuler dolu 145 kayitta tek istisna yok, PR belirleyici degil).
# Bu uc kayitta alt tip imputasyonla atanmis ve kurali ihlal etmis.
M2D = {'237R': 'Luminal B (HER2 Negatif)',   # Ki-67 Yuksek
       '132R': 'Luminal A',                   # Ki-67 Dusuk
       '294R': 'Luminal A'}                   # Ki-67 Dusuk
for hid, etiket in M2D.items():
    d.loc[eski == hid, 'Molekuler_i'] = etiket
    if 'i5' in d.columns:
        kod = d.loc[d.Molekuler_i == etiket, 'i5'].mode()
        if len(kod): d.loc[eski == hid, 'i5'] = kod.iloc[0]
log('2d', 'alt tip imputasyonu Ki-67 kuralini ihlal etmisti; ham Ki-67 esas alindi',
    3, '237R -> Luminal B (HER2-), 132R ve 294R -> Luminal A')


# ---------------------------------------------------------------- MADDE 2e
# Alt tipi hala ML imputasyonuyla atanmis kayitlarda alt tip, mevcut IHK panelinden
# kohortun kendi kuraliyla turetilir. Kural once ham degeri mevcut kayitlarda dogrulanir.
def alt_tip_kurali(er, pr, her, ki):
    hr = ('Pozitif' in str(er)) or ('Pozitif' in str(pr))
    if her == 'Ekvokal':    return 'HER2 ekvokal (ISH doğrulanmamış)'
    if her == 'Pozitif':    return 'Luminal B (HER2 Pozitif)' if hr else 'HER2-Zengin'
    if her == 'HER2-düşük': return 'HER2-Düşük' if hr else 'Triple Negatif'
    if not hr:              return 'Triple Negatif'
    return 'Luminal A' if ki == 'Düşük' else 'Luminal B (HER2 Negatif)'

kural = [alt_tip_kurali(*r) for r in
         d[['ER_i', 'PR_i', 'HER2_i', 'Ki_67_i']].astype(str).itertuples(index=False)]
ham_var = ~(d.Molekuler.isna() | d.Molekuler.astype(str).str.contains('Bilinmiyor', na=False))
ekv = d.Molekuler_i.astype(str).str.contains('ekvokal')
gecerli = ham_var & ~ekv
uyum = (pd.Series(kural)[gecerli.values].values == d.loc[gecerli, 'Molekuler_i'].values).mean()
print(f'   kural dogrulamasi: ham alt tipi mevcut {int(gecerli.sum())} kayitta '
      f'uyum %{100*uyum:.1f}')

hedef = (~ham_var) & (~ekv) & (~eski.isin(
    ['109R', '122L', '205L', '207L', '233R', '237R', '132R', '294R']))
degisen = int((pd.Series(kural)[hedef.values].values != d.loc[hedef, 'Molekuler_i'].values).sum())
d.loc[hedef, 'Molekuler_i'] = pd.Series(kural)[hedef.values].values
for et in d.loc[hedef, 'Molekuler_i'].unique():
    kod = d.loc[(d.Molekuler_i == et) & ~hedef, 'i5'].mode()
    if len(kod): d.loc[hedef & (d.Molekuler_i == et), 'i5'] = kod.iloc[0]

# koken kolonu — makalede alt tipin nereden geldigi izlenebilsin
d['Molekuler_kaynak'] = np.where(ham_var, 'ham patoloji kaydı',
                         np.where(hedef, 'İHK panelinden kural ile türetildi',
                                  'düzeltilmiş / ekvokal'))
log('2e', f'alt tipi eksik kayitlarda ML imputasyonu yerine IHK kurali kullanildi '
          f'(kural ham veriyle %{100*uyum:.1f} uyumlu; etiket degisen: {degisen})',
    int(hedef.sum()), 'Molekuler_kaynak kolonu eklendi')


# ---------------------------------------------------------------- MADDE 3
# HER2-dusuk bir tedavi kategorisidir, icsel molekuler alt tip degildir (St. Gallen / WHO).
# Alt tip Ki-67'ye gore luminal kategorilere dagitilir; HER2-dusuk durumu ayri ikili degisken.
low = d.HER2_i.astype(str) == 'HER2-düşük'
d['HER2_Dusuk'] = np.where(low, 'Var', 'Yok')
d['i65'] = np.where(low, 1, 0)
# yalnizca hormon reseptoru POZITIF olanlar luminal gruba dagitilir;
# ER-/PR- olan HER2-dusuk olgular triple negatif (HER2-low) olarak kalir
hr = (d.ER_i.astype(str).str.contains('Pozitif') | d.PR_i.astype(str).str.contains('Pozitif'))
low = low & hr
yeni = np.where(d.Ki_67_i.astype(str) == 'Düşük', 'Luminal A', 'Luminal B (HER2 Negatif)')
d.loc[low, 'Molekuler_i'] = yeni[low.values]
for et in ['Luminal A', 'Luminal B (HER2 Negatif)']:
    kod = d.loc[(d.Molekuler_i == et) & ~low, 'i5'].mode()
    if len(kod): d.loc[low & (d.Molekuler_i == et), 'i5'] = kod.iloc[0]
log('3', 'HER2-Düşük alt tip kategorisi kaldirildi; Ki-67 ile luminal gruplara dagitildi, '
         'HER2-dusuk durumu ayri ikili degisken oldu (HER2_Dusuk / i65)',
    int(low.sum()), 'St. Gallen / WHO sınıflamasında ayrı alt tip yok')

# ---------------------------------------------------------------- MADDE 1
# 1a. 262 cifti: ilk kayit sag, ikinci kayit sol
idx = d.index[eski.isin(['262L', '262R'])].tolist()
d.loc[idx[0], 'Hangi_Meme'] = 'Sağ'
d.loc[idx[1], 'Hangi_Meme'] = 'Sol'
if 'i16' in d.columns:
    sag_kod = d.loc[d.Hangi_Meme == 'Sağ', 'i16'].mode()
    sol_kod = d.loc[d.Hangi_Meme == 'Sol', 'i16'].mode()
    if len(sag_kod): d.loc[idx[0], 'i16'] = sag_kod.iloc[0]
    if len(sol_kod): d.loc[idx[1], 'i16'] = sol_kod.iloc[0]
log('1a', '262 ciftinde iki kayit da "Sol" idi; ilki "Sağ" olarak duzeltildi', 1,
    'ilk kayit -> Sağ')

# 1b. Iki katmanli kimlik: Hasta_ID (kisi) + Kayit_ID (satir)
kok = eski.str[:-1]
sira = {k: f'P{i+1:03d}' for i, k in enumerate(dict.fromkeys(kok))}
d.insert(0, 'Kayit_ID', [f'K{i+1:03d}' for i in range(len(d))])
d.insert(1, 'Hasta_ID_yeni', kok.map(sira).values)
d.insert(2, 'Eski_ID', eski.values)
d = d.drop(columns=['Hasta_ID']).rename(columns={'Hasta_ID_yeni': 'Hasta_ID'})
n_hasta = d.Hasta_ID.nunique(); n_bil = int((d.Hasta_ID.value_counts() > 1).sum())
log('1b', f'iki katmanli kimlik olusturuldu: {len(d)} kayit / {n_hasta} hasta '
          f'({n_bil} hastada bilateral)', len(d), 'Kayit_ID + Hasta_ID + Eski_ID')
log('1c', 'yas ile tani yili arasindaki 2 yillik fark birakildi (yas kullaniliyor)',
    32, 'degisiklik yok')


# ---------------------------------------------------------------- MADDE 6
# Kurativ niyetli karar problemi: uzak metastatik olgular hedef popülasyonda degil.
# Kayitlar silinmez; analiz kolonu ile isaretlenir (taban senaryo = cM0, duyarlilik = tumu).
cm1 = (d.Metastaz_Yeri_i == 'Uzak Metastaz') | (d.NAC_Once_Evre_i == 'Evre IV')
d['Analiz_Kohortu'] = np.where(cm1, 'duyarlılık (cM1)', 'taban senaryo (cM0)')
log('6', 'uzak metastatik olgular taban senaryodan cikarildi, duyarlilik analizinde '
         'raporlanacak (Analiz_Kohortu kolonu)', int(cm1.sum()),
    'taban senaryo n=%d' % int((~cm1).sum()))


# ---------------------------------------------------------------- MADDE 7
# Patolojik kanitli kanserde BI-RADS 1/2 olamaz. Kural: taniya en yakin, biyopsi oncesi
# son goruntulemenin BI-RADS'i gecerlidir. Degerler radyoloji kayitlarindan teyit edildi.
M7 = {'107R': ('4A', 'Solid kitle'), '028L': ('4B', 'Kalsifikasyon'),
      '035R': ('4C', 'Solid kitle'), '141R': ('5', 'Architectural Distortion'),
      '191L': ('4A', 'Asimetri'),
      '258L': ('1', 'Kalsifikasyon')}   # kayittan geldi; ic tutarsizlik icin bkz. asagi
for hid, (br, lez) in M7.items():
    d.loc[eski == hid, 'BI-RADS_i'] = br
    d.loc[eski == hid, 'Lezyon_Turu_i'] = lez
    for kol, deg in [('i48_e', br), ('i51_e', lez)]:
        if kol in d.columns:
            kaynak = 'BI-RADS_i' if kol == 'i48_e' else 'Lezyon_Turu_i'
            kod = d.loc[d[kaynak].astype(str) == deg, kol].mode()
            if len(kod): d.loc[eski == hid, kol] = kod.iloc[0]
# 258L: BI-RADS 1 ("bulgu yok") ile kalsifikasyon tarifi ve kanitli kanser bir arada
# olamaz; deger kayittan geldigi icin korunuyor ama cozulmemis olarak isaretleniyor.
d['Veri_Notu'] = ''
d.loc[eski == '258L', 'Veri_Notu'] = ('BI-RADS 1 ile kalsifikasyon ve kanitli kanser '
                                      'celisiyor; goruntulemenin tarihi teyit edilmeli')
log('7', 'BI-RADS ve lezyon turu radyoloji kayitlarindan teyit edilerek duzeltildi '
         '(imputasyon uc kayitta lezyon turunu yanlis atamisti); 258L cozulmemis',
    6, ', '.join(M7))


# ---------------------------------------------------------------- MADDE 8
# RCB-0 tanim geregi tumor yatáginda canli invaziv tumor bulunmamasi demektir; %CA != 0
# olamaz. RCB sinifi esas alindi. (Kohortta tek ornek; diger 89 RCB-0 olgusunda %CA = 0.)
d.loc[eski == '143R', '%CA'] = 0
log('8', 'RCB-0 olgusunda %CA = 20 idi; RCB sinifi esas alinarak %CA = 0 yapildi',
    1, '143R')


# ---------------------------------------------------------------- MADDE 10
# Allred skorlari ne bu makalede ne ML/ILR hattinda kullaniliyor; ayrica imputasyonlu
# surumleri bozuk ("Güçlü Pozitif" etiketli olgularda minimum 0). Kolonlar cikarildi.
# Not: 2c-B kararinin dayanagi olan uc olgunun skorlari duzeltme gunlugune yazildi.
ALLRED = [c for c in d.columns if 'Allred' in c]
d = d.drop(columns=ALLRED)
log('10', 'Allred kolonlari cikarildi (kullanilmiyor; imputasyonlu surumleri bozuk)',
    len(ALLRED), ', '.join(ALLRED))


# ---------------------------------------------------------------- MADDE 12
# Ham metin kolonlarindaki bicim varyantlari. Kodlanmis (_i / iX) surumler bunlari zaten
# tekillestirmis, yani model etkilenmemisti; duzeltme paylasilan veri seti icin.
duzeltilen = []
for c in d.columns:
    v = [x for x in d[c].dropna().unique() if isinstance(x, str)]
    if len(v) < 2: continue
    grup = {}
    for x in v: grup.setdefault(x.strip().lower(), []).append(x)
    for esas, varyantlar in grup.items():
        if len(varyantlar) > 1 or any(x != x.strip() for x in varyantlar):
            # en sik gecen bicimi standart kabul et
            std = d[c][d[c].isin(varyantlar)].value_counts().idxmax().strip()
            d[c] = d[c].replace({x: std for x in varyantlar})
            duzeltilen.append(f'{c}: {varyantlar} -> {std}')
# d1xd2 icinde carpma isareti iki farkli karakterle yazilmis (x ve ×) -> tekillestir
if 'd1xd2' in d.columns:
    d['d1xd2'] = (d['d1xd2'].astype(str).str.strip()
                  .str.replace('×', 'x', regex=False).str.replace(' X ', ' x ', regex=False))
# Tani Yili: "2022(2015)" -> 2022; kolon sayisala cevrildi
if 'Tanı Yılı' in d.columns:
    d['Tanı Yılı'] = pd.to_numeric(
        d['Tanı Yılı'].astype(str).str.extract(r'(\d{4})')[0], errors='coerce').astype('Int64')
log('12', 'ham metin kolonlarindaki bicim varyantlari tekillestirildi; Tani Yili sayisala '
          'cevrildi ("2022(2015)" -> 2022); d1xd2 carpma isareti standartlastirildi',
    len(duzeltilen) + 1, ' | '.join(duzeltilen[:6]))


# ---------------------------------------------------------------- MADDE 14
# Rejim bilgisi "Bilinmiyor" olan 31 olgu, ek hasta dosyalarindan antrasiklin + taksan
# olarak teyit edildi. (Imputasyon 28'ini dogru, 3'unu "Sadece Antrasiklin" atamisti.)
bil = d.Rejim.astype(str).str.contains('Bilinmiyor', na=False)
SADECE_ANT = ['011R', '020R', '303L']          # ek dosyalardan teyit
d.loc[bil, 'Rejim'] = 'Antrasiklin + Taksan'
d.loc[bil & eski.isin(SADECE_ANT), 'Rejim'] = 'Sadece Antrasiklin'
for kolon, kod_kaynagi in [('Rejim_46', 'Rejim_46'), ('i46', 'i46')]:
    if kolon in d.columns:
        for et in ['Antrasiklin + Taksan', 'Sadece Antrasiklin']:
            kod = d.loc[(d.Rejim_i == et) & ~bil, kod_kaynagi].mode()
            if len(kod): d.loc[bil & (d.Rejim == et), kolon] = kod.iloc[0]
log('14', 'rejim bilgisi eksik olgular ek hasta dosyalarindan teyit edildi; imputasyon '
          '31/31 dogru bilmis (28 antrasiklin+taksan, 3 sadece antrasiklin)', int(bil.sum()),
    'ham Rejim kolonu dolduruldu; imputasyonlu surumler zaten dogruydu')



# ---------------------------------------------------------------- MADDE 17
# ISH/FISH sonuclari patoloji arsivinden temin edildi. Ekvokal (IHK 2+) olgularin
# dokuzunda HER2 negatif dogrulandi: yedisinde raporda acikca "HER2 negatif",
# ikisinde (173R, 287L) "non-amplifiye" ifadesi. 129R'nin sonucuna ulasilamadi.
# Bu olgular artik HR+/HER2- biyolojide olup karar kuralinin uygunluk kapsamina girer.
ISH_NEG = ['039L', '053L', '112L', '173R', '189R', '203L', '239R', '271R', '287L']
for hid in ISH_NEG:
    m = eski == hid
    if not m.any(): continue
    d.loc[m, 'HER2_i'] = 'Negatif'
    ki = d.loc[m, 'Ki_67_i'].iloc[0]
    yeni = 'Luminal A' if str(ki) == 'Düşük' else 'Luminal B (HER2 Negatif)'
    d.loc[m, 'Molekuler_i'] = yeni
    d.loc[m, 'Veri_Notu'] = 'HER2 durumu ISH/FISH ile doğrulandı (negatif)'
log('17', 'ekvokal olgularda ISH/FISH sonuclari temin edildi; dokuz olgu HER2 negatif '
          'olarak dogrulandi ve HR+/HER2- kapsamina alindi (129R belirsiz kaldi)',
    len(ISH_NEG), ', '.join(ISH_NEG))

# ---------------------------------------------------------------- MADDE 16
# Butunluk onarimi: yapilan tum duzeltmelerden sonra kodlanmis kolonlarin (iX ve _NN)
# etiket kolonlariyla birebir olmasi saglanir. Madde 0 duzeltmeleri iX kolonlarina
# yansimamisti (i13, i14, i15); ayrica bazi _NN ara kodlari eski degerde kalmisti.
IX = {'i1':'Histolojik_Yeni','i2':'ER_i','i3':'PR_i','i4':'HER2_i','i5':'Molekuler_i',
 'i6':'Ki_67_i','i7':'Tubul_i','i8':'Nukleer_i','i9':'Mitotik_i','i10':'HistolojikG_i',
 'i11':'Ecaderin_i','i12':'TIL_i','i13':'Metasdaz_Durumu_i','i14':'Metastaz_Yeri_i',
 'i15':'NAC_Once_Evre_i','i16':'Hangi_Meme','i17':'VKI_Sinifi_17','i18':'Yas_Grubu',
 'i19':'Kan_Grup_i','i20':'Menapoz_i','i21':'HT','i22':'DM','i23':'KOAH','i24':'Sigara',
 'i25':'Ailede_MemeCA','i26':'Tiroid','i27':'Retinopati','i28':'Nöropati','i29':'Osteoporoz',
 'i30':'Depresyon','i31':'ALP_Durumu','i32':'ALT_Durumu','i33':'AST_Durumu','i34':'BUN_Durumu',
 'i35':'CA15-3_Durumu','i36':'CEA_Durumu','i37':'CRP_Durumu','i38':'GGT_Durumu',
 'i39':'Glukoz_Durumu','i40':'HbA1c_Durumu','i41':'Kreatinin_Durumu','i42':'LDH_Durumu',
 'i43':'TSH_Durumu','i44':'e-GFR_Durumu','i45':'Gunesten_Yararlanma','i46':'Rejim_i',
 'i47':'Kur_Yogunluk_i'}
onarim = []
for kod, etiket in IX.items():
    if kod not in d.columns or etiket not in d.columns: continue
    # her etiket icin en sik gecen kodu standart kabul et, sapan kayitlari duzelt
    harita = d.groupby(d[etiket].astype(str))[kod].agg(lambda s: s.mode().iloc[0])
    yeni = d[etiket].astype(str).map(harita)
    sapan = int((d[kod] != yeni).sum())
    if sapan:
        onarim.append(f'{kod}({sapan})')
        d[kod] = yeni
# _NN ara kodlari: 0 "ham kayitta yoktu" demektir ve korunur; sifir olmayan sapmalar duzeltilir
ARA = [('Molekuler_5','Molekuler_i'), ('HER2_4','HER2_i'), ('BI-RADS_48','BI-RADS_i'),
       ('Lezyon_Turu_51','Lezyon_Turu_i'), ('Metastaz_Yeri_14','Metastaz_Yeri_i'),
       ('NAC_Once_Evre_15','NAC_Once_Evre_i'), ('Metasdaz_Durumu_13','Metasdaz_Durumu_i'),
       ('Rejim_46','Rejim_i')]
for ara, etiket in ARA:
    if ara not in d.columns or etiket not in d.columns: continue
    gecerli = d[ara] != 0
    if not gecerli.any(): continue
    harita = (d[gecerli].groupby(d.loc[gecerli, etiket].astype(str))[ara]
              .agg(lambda s: s.mode().iloc[0]))
    yeni = d[etiket].astype(str).map(harita)
    sapan = gecerli & yeni.notna() & (d[ara] != yeni)
    if sapan.any():
        onarim.append(f'{ara}({int(sapan.sum())})')
        d.loc[sapan, ara] = yeni[sapan]
log('16', 'kodlanmis kolonlar etiket kolonlariyla yeniden esitlendi (Madde 0 duzeltmeleri '
          'iX kolonlarina yansimamisti)', len(onarim), ', '.join(onarim) or 'sapma yok')

# ---------------------------------------------------------------- CIKTI
kimlik = ['HASTA ADI', 'HASTA NO', 'Doğum Yeri', 'Enlem', 'Boylam']
d.drop(columns=[c for c in kimlik if c in d.columns]).to_excel(
    OUT + f'Kohort_{SURUM}.xlsx', index=False)

with open(OUT + f'Duzeltme_Gunlugu_{SURUM}.md', 'w') as f:
    f.write(f'# KOHORT DÜZELTME GÜNLÜĞÜ — {SURUM}\n\n')
    f.write(f'Kaynak: `Kagan_TEZ.xlsx` · Çıktı: `Kohort_{SURUM}.xlsx` · '
            f'{len(d)} kayıt / {n_hasta} hasta\n\n')
    f.write('| Madde | Düzeltme | Kayıt | Karar |\n|---|---|---|---|\n')
    for m, ne, kac, kr in LOG:
        f.write(f'| {m} | {ne} | {kac} | {kr} |\n')
    f.write('\n## Kimlik yapısı\n\n')
    f.write('- `Kayit_ID` — satır anahtarı (K001–K%03d), her meme bir kayıt\n' % len(d))
    f.write('- `Hasta_ID` — kişi anahtarı (P001–P%03d), bilateral olgularda iki '
            'kayıtta aynı\n' % n_hasta)
    f.write('- `Eski_ID` — özgün numaralandırma (izlenebilirlik için korundu)\n')
    f.write('- `Hangi_Meme` — taraf bilgisinin tek kaynağı\n\n')
    f.write('Kimlik bilgisi içeren kolonlar çıkarıldı: ' +
            ', '.join(f'`{c}`' for c in kimlik) + '\n')

print(f'\nKohort_{SURUM}.xlsx yazildi: {d.shape}')
print(f'bilateral hasta: {n_bil} | tekil hasta: {n_hasta} | kayit: {len(d)}')
print(d[['Kayit_ID', 'Hasta_ID', 'Eski_ID', 'Hangi_Meme']].head(6).to_string(index=False))
print('...')
print(d.loc[d.Hasta_ID.isin(d.Hasta_ID[d.Hasta_ID.duplicated()].head(3)),
            ['Kayit_ID', 'Hasta_ID', 'Eski_ID', 'Hangi_Meme']].to_string(index=False))
