"""KALEM KALEM MALIYET MODELI — NAK yolu vs cerrahi-once yolu
Tum birim fiyatlar: SUT EK-2 (29.06.2026) + TITCK/EK-4A kamu kurum iskontosu uygulanmis.
Standart hasta: HR+/HER2-, BSA 1,7 m2, kur atlanmadan tamamlanan tedavi.
"""
import pandas as pd

BSA = 1.7

# ---------------------------------------------------------------- BIRIM FIYATLAR
SUT = {  # EK-2/B ve EK-2/C, dogrudan TL
 'kemoterapi_uygulama':1482.02, 'gunduz_yatak':74.96, 'port':1321.34,
 'hemogram':33.27, 'onkoloji_muayene':356.00, 'genel_cerrahi_muayene':356.00,
 'radyasyon_onk_muayene':288.00, 'ekokardiyografi':245.55,
 'alt':12.15, 'ast':11.08, 'kan_kulturu':111.00, 'periferik_yayma':11.08,
 'standart_yatak':374.97, 'steril_oda':1674.76,
 'mamografi':188.76, 'meme_mr_kontrastli':1932.32, 'igne_biyopsi':503.13,
 'ihk':312.05, 'fish':1561.30, 'kemik_sintigrafi':1158.85,
 # cerrahi (EK-2/C)
 'mkc':17064.645, 'mkc_slnb':22226.451, 'mkc_alnd':22226.451,
 'mastektomi_basit':21955.219, 'mastektomi_mrm':71932.77,
 # radyoterapi (konformal, 16 fraksiyon + tasarim/planlama/doz)
 'rt_tasarim':3084.96, 'rt_planlama':9042.06, 'rt_doz':2768.94, 'rt_seans':1076.49,
}
ILAC = {  # SGK odemesi (kamu kurum iskontosu uygulanmis)
 'doksorubisin50':371.75, 'doksorubisin10':201.19, 'siklofosfamid500':486.53,
 'dosetaksel80':4690.08, 'dosetaksel40':2523.47, 'dosetaksel20':2713.07,
 'paklitaksel100':1774.75, 'paklitaksel30':160.26,
 'palonosetron':576.62, 'fosaprepitant':1791.46, 'aprepitant_3gun':711.48,
 'deksametazon_amp':43.84, 'deksametazon_8mg_20tb':411.67,
 'filgrastim_5enj':2589.83, 'pegfilgrastim':8001.89,
 'pipTazo_45g':258.89,
 'kapesitabin_500_120tb':2448.68,
 'anastrozol_28tb':486.68, 'letrozol_30tb':414.93, 'tamoksifen_30tb':264.77,
}

def yaz(baslik): print('\n' + '=' * 78); print(baslik); print('=' * 78)
def sat(ad, adet, birim, not_=''):
    t = adet * birim
    print(f'  {ad:44s} {adet:6.1f} × {birim:9.2f} = {t:11.2f}  {not_}')
    return t

# ---------------------------------------------------------------- A. NAK YOLU
yaz('A. NEOADJUVAN KEMOTERAPİ YOLU — sitotoksik ilaç')
ac_kur = 2*ILAC['doksorubisin50'] + ILAC['doksorubisin10'] + 3*ILAC['siklofosfamid500']
tax_kur = ILAC['dosetaksel80'] + ILAC['dosetaksel40'] + ILAC['dosetaksel20']
top = 0
top += sat('AC kürü (doks 110 mg + siklo 1500 mg)', 4, ac_kur)
top += sat('Dosetaksel kürü (140 mg)', 4, tax_kur, '(en ucuz eşdeğer)')
ilac_toplam = top

yaz('A. NEOADJUVAN KEMOTERAPİ YOLU — destek tedavi (kür başına)')
d = 0
d += sat('Palonosetron 0,25 mg IV', 8, ILAC['palonosetron'])
d += sat('Aprepitant 125/80/80 mg (3 gün)', 8, ILAC['aprepitant_3gun'])
d += sat('Deksametazon 8 mg amp (gün 1)', 8, ILAC['deksametazon_amp'])
d += sat('Deksametazon tb (taksan premedikasyonu)', 4, ILAC['deksametazon_8mg_20tb']/2,
         '(yarım kutu/kür)')
destek = d

yaz('A. NEOADJUVAN KEMOTERAPİ YOLU — uygulama ve izlem')
u = 0
u += sat('Günübirlik kemoterapi uygulaması', 8, SUT['kemoterapi_uygulama'])
u += sat('Port kateter yerleştirilmesi', 1, SUT['port'])
u += sat('Hemogram (her kür öncesi)', 8, SUT['hemogram'])
u += sat('ALT + AST', 8, SUT['alt']+SUT['ast'])
u += sat('Tıbbi onkoloji muayenesi', 8, SUT['onkoloji_muayene'])
u += sat('Ekokardiyografi (antrasiklin izlemi)', 2, SUT['ekokardiyografi'])
uygulama = u

yaz('A. NEOADJUVAN KEMOTERAPİ YOLU — febril nötropeni epizodu (olay başına)')
f = 0
f += sat('Yatış (standart yatak, 7 gün)', 7, SUT['standart_yatak'])
f += sat('Piperasilin-tazobaktam 4,5 g ×3/gün', 21, ILAC['pipTazo_45g'])
f += sat('Filgrastim 30 MIU (5 enjeksiyon)', 1, ILAC['filgrastim_5enj'])
f += sat('Kan kültürü (2 set)', 2, SUT['kan_kulturu'])
f += sat('Hemogram + periferik yayma (günlük)', 7, SUT['hemogram']+SUT['periferik_yayma'])
fn_epizod = f
print(f'\n  FN epizod maliyeti: {fn_epizod:.2f} TL')
for p in [0.10, 0.15, 0.20]:
    print(f'    beklenen maliyet (%{p*100:.0f} insidans): {p*fn_epizod:9.2f} TL')

yaz('NAK YOLU — TOPLAM (FN %15 insidans varsayımıyla)')
nak_toplam = ilac_toplam + destek + uygulama + 0.15*fn_epizod
for ad, v in [('Sitotoksik ilaç', ilac_toplam), ('Destek tedavi', destek),
              ('Uygulama ve izlem', uygulama), ('FN beklenen maliyeti', 0.15*fn_epizod)]:
    print(f'  {ad:44s} {v:11.2f}  (%{100*v/nak_toplam:4.1f})')
print(f'  {"TOPLAM":44s} {nak_toplam:11.2f}')

# ---------------------------------------------------------------- B. CERRAHI
yaz('B. CERRAHİ — iki kolda dağılım farkı')
print('  Birim fiyatlar:')
for k in ['mkc','mkc_slnb','mkc_alnd','mastektomi_basit','mastektomi_mrm']:
    print(f'    {k:22s} {SUT[k]:11.2f}')
print('\n  Literatür: HR+/HER2− NAK sonrası meme koruyucu cerrahi oranı %34,5 (Z1071)')
print('             NCDB: NAK sonrası %42,3 · neoadjuvan endokrin sonrası %64,0')
print('  ⬛ Cerrahi-önce kolunda MKC oranı kohorttan/literatürden belirlenecek')
for mkc_nak, mkc_cer in [(0.42, 0.30), (0.42, 0.35), (0.345, 0.25)]:
    c_nak = mkc_nak*SUT['mkc_slnb'] + (1-mkc_nak)*SUT['mastektomi_mrm']
    c_cer = mkc_cer*SUT['mkc_alnd'] + (1-mkc_cer)*SUT['mastektomi_mrm']
    print(f'    MKC oranı NAK %{100*mkc_nak:.1f} / cerrahi-önce %{100*mkc_cer:.0f}  ->  '
          f'NAK {c_nak:10.2f} · cerrahi-önce {c_cer:10.2f} · fark {c_cer-c_nak:+10.2f}')

# ---------------------------------------------------------------- C. RADYOTERAPI
yaz('C. RADYOTERAPİ (meme koruyucu cerrahi sonrası, konformal 16 fraksiyon)')
r = SUT['rt_tasarim'] + SUT['rt_planlama'] + SUT['rt_doz'] + 16*SUT['rt_seans'] \
    + SUT['radyasyon_onk_muayene']
print(f'  tasarım + planlama + doz + 16 seans + muayene = {r:11.2f} TL')

# ---------------------------------------------------------------- D. ADJUVAN
yaz('D. ADJUVAN ENDOKRİN TEDAVİ (5 yıl)')
for ad, kutu in [('Anastrozol (jenerik, 28 tb)', ILAC['anastrozol_28tb']),
                 ('Letrozol (jenerik, 30 tb)', ILAC['letrozol_30tb']),
                 ('Tamoksifen (30 tb)', ILAC['tamoksifen_30tb'])]:
    yillik = kutu * (365/28 if '28' in ad else 365/30)
    print(f'  {ad:36s} yıllık {yillik:10.2f} · 5 yıl {5*yillik:11.2f}')
print(f'\n  Kapesitabin (CREATE-X, 8 kür): 500 mg 120 tb kutu {ILAC["kapesitabin_500_120tb"]:.2f} TL')
print('  ⬛ HR+/HER2−de CREATE-X faydası sınırlı; taban senaryoda uygulanmayacak')
