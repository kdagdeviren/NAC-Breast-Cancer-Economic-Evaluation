"""KOHORT — DERIN KLINIK TUTARLILIK DENETIMI (Kohort v16)
Her kural: klinik/literatur gerekcesi + ihlal eden kayitlar.
"""
import pandas as pd, numpy as np, re, warnings
warnings.filterwarnings('ignore')
pd.set_option('display.width', 220)
d = pd.read_excel('../../outputs/Kohort_v16.xlsx')
S = lambda c: d[c].astype(str)
BULGU = []
def kural(kod, baslik, mask, gerekce, goster=('Kayit_ID','Eski_ID')):
    n = int(mask.sum())
    dur = '✗' if n else '✓'
    print(f'\n[{kod}] {dur} {baslik}  → {n} kayıt')
    print(f'      gerekçe: {gerekce}')
    if n:
        BULGU.append((kod, baslik, n))
        kol = [c for c in goster if c in d.columns]
        print(d.loc[mask, kol].head(12).to_string(index=False))
    return mask

print('=' * 100); print('1. RADYOLOJİ — LEZYON TÜRÜ İLE TANIMLAYICILARIN UYUMU'); print('=' * 100)
kitle = S('Lezyon_Turu_i').str.contains('kitle|Kitle', na=False)
kitle_var = ~S('Kitle_Sekli_i').str.contains('Yok|Değerlend', na=False)
kural('R1', 'Lezyon türü kitle değil ama kitle şekli tanımlı',
      (~kitle) & kitle_var, 'Kitle tanımlayıcıları yalnızca kitle varlığında kaydedilir',
      ('Kayit_ID','Eski_ID','Lezyon_Turu_i','Kitle_Sekli_i','Kitle_Konturu_i'))
kural('R2', 'Lezyon türü kitle ama kitle şekli tanımsız',
      kitle & (~kitle_var), 'Kitle bildirilmişse şekli de tanımlanmalı',
      ('Kayit_ID','Eski_ID','Lezyon_Turu_i','Kitle_Sekli_i'))
kalsif_m = ~S('Kalsifikasyon_Morfolojisi_i').str.contains('Yok|Değerlend', na=False)
kalsif_d = ~S('Kalsifikasyon_Dagilimi_i').str.contains('Yok|Değerlend', na=False)
kural('R3', 'Kalsifikasyon morfolojisi yok ama dağılımı tanımlı', (~kalsif_m) & kalsif_d,
      'Dağılım ancak kalsifikasyon varsa tanımlanır',
      ('Kayit_ID','Eski_ID','Kalsifikasyon_Morfolojisi_i','Kalsifikasyon_Dagilimi_i'))
kural('R4', 'Lezyon türü kalsifikasyon ama morfoloji tanımsız',
      S('Lezyon_Turu_i').eq('Kalsifikasyon') & (~kalsif_m),
      'Kalsifikasyon lezyonunda morfoloji zorunlu tanımlayıcıdır',
      ('Kayit_ID','Eski_ID','Lezyon_Turu_i','Kalsifikasyon_Morfolojisi_i'))
kural('R5', 'Lezyon türü asimetri ama asimetri tipi tanımsız',
      S('Lezyon_Turu_i').eq('Asimetri') & S('Asimetri_i').str.contains('Yok|Değerlend', na=False),
      'Asimetri lezyonunda tipi tanımlanmalı',
      ('Kayit_ID','Eski_ID','Lezyon_Turu_i','Asimetri_i'))
kural('R6', 'BI-RADS 4C/5 ama iki yıldır stabil',
      S('BI-RADS_i').isin(['4C','5']) & S('2_Yildir_Stabil_i').str.contains('Var|Evet|Stabil', na=False),
      'İki yıl stabil lezyon yüksek malignite şüphesi kategorisine girmez',
      ('Kayit_ID','Eski_ID','BI-RADS_i','2_Yildir_Stabil_i'))
kural('R7', 'Kozmetik implant var ama meme dansitesi yağlı',
      S('Kozmetik_Implant').str.contains('Evet|Var', na=False) & S('Meme_Dansite_i').str.contains('Yağlı', na=False),
      'Bilgi amaçlı; implantlı memede dansite değerlendirmesi güvenilirliği düşer',
      ('Kayit_ID','Eski_ID','Kozmetik_Implant','Meme_Dansite_i'))

print('\n' + '=' * 100); print('2. PATOLOJİ — PROLİFERASYON VE GRADE'); print('=' * 100)
kural('P1', 'Grade 1 ama Ki-67 yüksek',
      S('HistolojikG_i').str.contains('1', na=False) & S('Ki_67_i').eq('Yüksek'),
      'Grade 1 düşük proliferasyon demektir; yüksek Ki-67 ile birlikte beklenmez',
      ('Kayit_ID','Eski_ID','HistolojikG_i','Ki_67_i','Mitotik_i'))
kural('P2', 'Grade 3 ama Ki-67 düşük',
      S('HistolojikG_i').str.contains('3', na=False) & S('Ki_67_i').eq('Düşük'),
      'Grade 3 yüksek proliferasyon demektir',
      ('Kayit_ID','Eski_ID','HistolojikG_i','Ki_67_i','Mitotik_i'))
kural('P3', 'Mitotik derece 3 ama Ki-67 düşük',
      S('Mitotik_i').str.contains('3', na=False) & S('Ki_67_i').eq('Düşük'),
      'Mitotik sayı ve Ki-67 aynı biyolojik süreci ölçer',
      ('Kayit_ID','Eski_ID','Mitotik_i','Ki_67_i'))
kural('P4', 'Lobüler karsinom ama grade 3',
      S('Histolojik_Yeni').str.contains('Lobüler|Lobuler', na=False) & S('HistolojikG_i').str.contains('3', na=False),
      'Klasik invaziv lobüler karsinom tipik olarak grade 2\'dir; grade 3 seyrektir',
      ('Kayit_ID','Eski_ID','Histolojik_Yeni','HistolojikG_i','Ecaderin_i'))
kural('P5', 'Triple negatif ama grade 1',
      S('Molekuler_i').eq('Triple Negatif') & S('HistolojikG_i').str.contains('1', na=False),
      'Üçlü negatif tümörler baskın olarak grade 3\'tür',
      ('Kayit_ID','Eski_ID','Molekuler_i','HistolojikG_i','Ki_67_i'))
kural('P6', 'Luminal A ama grade 3',
      S('Molekuler_i').eq('Luminal A') & S('HistolojikG_i').str.contains('3', na=False),
      'Luminal A düşük proliferasyonlu, tipik olarak grade 1-2 tümördür',
      ('Kayit_ID','Eski_ID','Molekuler_i','HistolojikG_i','Ki_67_i'))

print('\n' + '=' * 100); print('3. İMMÜNHİSTOKİMYA VE ALT TİP'); print('=' * 100)
erp = S('ER_i').str.contains('Pozitif'); prp = S('PR_i').str.contains('Pozitif')
kural('I1', 'ER negatif ama PR pozitif', (~erp) & prp,
      'ER−/PR+ fenotipi literatürde büyük ölçüde teknik artefakt kabul edilir',
      ('Kayit_ID','Eski_ID','ER_i','PR_i','HER2_i','Molekuler_i'))
kural('I2', 'HER2-düşük ikili değişkeni ile HER2 immünhistokimyası uyumsuz',
      (S('HER2_Dusuk').eq('Var') & (~S('HER2_i').eq('HER2-düşük'))) |
      (S('HER2_Dusuk').eq('Yok') & S('HER2_i').eq('HER2-düşük')),
      'İkili değişken doğrudan İHK sonucundan türetilmiştir',
      ('Kayit_ID','Eski_ID','HER2_i','HER2_Dusuk'))
kural('I3', 'Triple negatif ama TIL düşük ve grade 1',
      S('Molekuler_i').eq('Triple Negatif') & S('TIL_i').str.contains('<%10') &
      S('HistolojikG_i').str.contains('1', na=False),
      'Bilgi amaçlı: üçlü negatif tümörlerde TIL yoğunluğu tipik olarak yüksektir',
      ('Kayit_ID','Eski_ID','Molekuler_i','TIL_i','HistolojikG_i'))

print('\n' + '=' * 100); print('4. EVRELEME VE CERRAHİ PATOLOJİ'); print('=' * 100)
evre = S('NAC_Once_Evre_i')
kural('E1', 'Evre I ama cerrahide çok sayıda tutulan nod (N>3)',
      evre.str.startswith('Evre I') & ~evre.str.startswith('Evre II') &
      ~evre.str.startswith('Evre III') & ~evre.str.startswith('Evre IV') &
      (pd.to_numeric(d.N, errors='coerce').fillna(0) > 3),
      'Klinik evre I nodal tutulum yok demektir; NAK sonrası 3\'ten fazla tutulu nod uyumsuz',
      ('Kayit_ID','Eski_ID','NAC_Once_Evre_i','Metastaz_Yeri_i','N','M'))
kural('E2', 'Metastaz yeri lenf nodu ama cerrahide hiç tutulan nod yok ve RCB-0 değil',
      S('Metastaz_Yeri_i').eq('Lenf Nodu') & (pd.to_numeric(d.N, errors='coerce').fillna(0) == 0) &
      (d.RCB_Kategorize > 0),
      'Bilgi amaçlı: nodal yanıt olabilir, ancak yüksek rezidüel yükle birlikte beklenmez',
      ('Kayit_ID','Eski_ID','Metastaz_Yeri_i','N','RCB_Kategorize'))
kural('E3', 'Metastaz yok ama cerrahide tutulan nod var',
      S('Metastaz_Yeri_i').eq('Yok') & (pd.to_numeric(d.N, errors='coerce').fillna(0) > 0),
      'Klinik nodal tutulum yokken patolojik tutulum olabilir (okült); yaygınlığı kontrol edilmeli',
      ('Kayit_ID','Eski_ID','Metastaz_Yeri_i','NAC_Once_Evre_i','N'))
kural('E4', 'Cilt çekintisi veya meme başı retraksiyonu var ama evre I-II',
      (S('Cilt_Cekintisi_i').str.contains('Var|Evet', na=False) |
       S('Meme_Basi_Retraksiyonu_i').str.contains('Var|Evet', na=False)) &
      evre.str.match(r'Evre I[AB]?$|Evre II[AB]$'),
      'Cilt/meme başı tutulumu T4b kriteridir ve evre IIIB gerektirir',
      ('Kayit_ID','Eski_ID','Cilt_Cekintisi_i','Meme_Basi_Retraksiyonu_i','NAC_Once_Evre_i'))
kural('E5', 'Dermal lenfatik tutulum ama evre IIIB değil',
      S('Metastaz_Yeri_i').eq('Dermal Lenfatik') & (~evre.isin(['Evre IIIB','Evre IIIC','Evre IV'])),
      'Dermal lenfatik invazyon inflamatuvar meme kanseri (T4d) göstergesidir',
      ('Kayit_ID','Eski_ID','Metastaz_Yeri_i','NAC_Once_Evre_i'))

print('\n' + '=' * 100); print('5. RCB BİLEŞENLERİ VE TEDAVİ YANITI'); print('=' * 100)
def alan(x):
    m = re.findall(r'(\d+(?:[.,]\d+)?)', str(x))
    return float(m[0].replace(',','.')) * float(m[1].replace(',','.')) if len(m) >= 2 else np.nan
A = d.d1xd2.map(alan)
kural('B1', 'RCB-0 ama tümör yatağı alanı sıfır değil',
      (d.RCB_Kategorize == 0) & (A > 0) & (pd.to_numeric(d['%CA'], errors='coerce').fillna(0) > 0),
      'RCB-0 tanım gereği canlı invaziv tümör yokluğudur',
      ('Kayit_ID','Eski_ID','d1xd2','%CA','RCB_Kategorize'))
kural('B2', 'RCB-III ama tümör yatağı çok küçük (<100 mm²) ve nod yok',
      (d.RCB_Kategorize == 3) & (A < 100) & (pd.to_numeric(d.N, errors='coerce').fillna(0) == 0),
      'Yüksek rezidüel yük büyük tümör yatağı veya nodal tutulum gerektirir',
      ('Kayit_ID','Eski_ID','d1xd2','%CA','N','RCB_Skor','RCB_Kategorize'))
kural('B3', '%CIS > %CA (in situ bileşen invazivden büyük) ve RCB-0 değil',
      (pd.to_numeric(d['%CIS'], errors='coerce').fillna(0) >
       pd.to_numeric(d['%CA'], errors='coerce').fillna(0)) & (d.RCB_Kategorize > 0),
      'Bilgi amaçlı: mümkündür ancak seyrektir, kayıt kontrolü değer taşır',
      ('Kayit_ID','Eski_ID','%CA','%CIS','RCB_Kategorize'))
kural('B4', 'Eksik kür ama RCB-0 (tam yanıt)',
      S('Kur_Yogunluk_i').str.contains('Eksik', na=False) & (d.RCB_Kategorize == 0),
      'Bilgi amaçlı: erken tam yanıt nedeniyle tedavi kesilmiş olabilir',
      ('Kayit_ID','Eski_ID','Kur_Yogunluk_i','RCB_Kategorize','Rejim_i'))

print('\n' + '=' * 100); print('6. DEMOGRAFİ VE ANTROPOMETRİ'); print('=' * 100)
yas = pd.to_numeric(d.Tani_Yasi, errors='coerce')
kural('D1', 'Yaş grubu etiketi ile sayısal yaş uyumsuz',
      (S('Yas_Grubu').str.contains('Genç', na=False) & (yas >= 40)) |
      (S('Yas_Grubu').str.contains('Yaşlı', na=False) & (yas < 65)),
      'Etiket ile sayısal değer aynı kaydın iki gösterimi',
      ('Kayit_ID','Eski_ID','Yas_Grubu','Tani_Yasi'))
vki = pd.to_numeric(d.VKI, errors='coerce')
kural('D2', 'VKİ sınıfı ile sayısal VKİ uyumsuz',
      (S('VKI_Sinifi_17').str.contains('Zayıf', na=False) & (vki >= 18.5)) |
      (S('VKI_Sinifi_17').str.contains('Obez', na=False) & (vki < 30)) |
      (S('VKI_Sinifi_17').str.contains('Normal', na=False) & ((vki < 18.5) | (vki >= 25))),
      'Etiket ile sayısal değer aynı kaydın iki gösterimi',
      ('Kayit_ID','Eski_ID','VKI_Sinifi_17','VKI'))
kural('D3', 'Postmenopozal ama 40 yaş altı', S('Menapoz_i').str.contains('Var|Post', na=False) & (yas < 40),
      'Doğal menopoz 40 yaş altında seyrektir',
      ('Kayit_ID','Eski_ID','Menapoz_i','Tani_Yasi'))
kural('D4', 'Premenopozal ama 60 yaş üstü', S('Menapoz_i').str.contains('Yok|Pre', na=False) & (yas > 60),
      '60 yaş üstünde premenopozal durum seyrektir',
      ('Kayit_ID','Eski_ID','Menapoz_i','Tani_Yasi'))

print('\n' + '=' * 100); print('7. LABORATUVAR İÇ TUTARLILIĞI'); print('=' * 100)
kre = pd.to_numeric(d.Kreatinin_Değer, errors='coerce')
egfr = pd.to_numeric(d['e-GFR (CKD-EPI)_Değer'], errors='coerce') if 'e-GFR (CKD-EPI)_Değer' in d.columns else pd.Series(np.nan, index=d.index)
kural('L1', 'Kreatinin yüksek ama e-GFR normal', (kre > 1.2) & (egfr > 90),
      'Yüksek kreatinin ile korunmuş e-GFR birlikte beklenmez',
      ('Kayit_ID','Eski_ID','Kreatinin_Değer','e-GFR (CKD-EPI)_Değer'))
ast = pd.to_numeric(d.AST_Değer, errors='coerce'); alt_ = pd.to_numeric(d.ALT_Değer, errors='coerce')
kural('L2', 'AST veya ALT > 200 U/L', (ast > 200) | (alt_ > 200),
      'Bu düzeyde transaminaz yüksekliği kemoterapi başlangıcında beklenmez',
      ('Kayit_ID','Eski_ID','AST_Değer','ALT_Değer'))
hb = pd.to_numeric(d.HbA1c_Değer, errors='coerce'); glu = pd.to_numeric(d.Glukoz_Değer, errors='coerce')
kural('L3', 'HbA1c > 10 ama glukoz normal', (hb > 10) & (glu < 110),
      'Belirgin yüksek HbA1c ile normal açlık glukozu birlikte beklenmez',
      ('Kayit_ID','Eski_ID','HbA1c_Değer','Glukoz_Değer'))
ca = pd.to_numeric(d['CA15-3_Değer'], errors='coerce')
kural('L4', 'CA15-3 > 100 U/mL ama uzak metastaz yok', (ca > 100) & (~S('Metastaz_Yeri_i').eq('Uzak Metastaz')),
      'Belirgin yüksek CA15-3 metastatik hastalık göstergesi olabilir',
      ('Kayit_ID','Eski_ID','CA15-3_Değer','Metastaz_Yeri_i','NAC_Once_Evre_i'))

print('\n' + '=' * 100); print(f'ÖZET — ihlal saptanan kural sayısı: {len(BULGU)}')
for k, b, n in BULGU: print(f'  [{k}] {b} — {n} kayıt')
