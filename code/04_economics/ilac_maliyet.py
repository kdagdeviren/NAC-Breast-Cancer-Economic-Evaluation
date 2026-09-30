"""ILAC MALIYET HESAPLAYICISI — TITCK fiyati + EK-4/A iskontosu -> SGK odemesi
Kullanim: python ilac_maliyet.py
Cikti: ../../outputs/ilac_birim_maliyet.csv
"""
import pandas as pd, numpy as np, math, warnings
warnings.filterwarnings('ignore')

TITCK = '../../data/TİTCK_Detaylı_İlaç_Fiyat_Listesi.xlsx'
EK4A  = '../../data/EK-4A_Bedeli_Ödenecek_İlaçlar_Listesi__29_06_20266_.xlsx'

t = pd.read_excel(TITCK)
t['bar'] = t.BARKODU.astype(str).str.strip()
t['ad']  = t['ILAC ADI'].astype(str).str.strip()

e = pd.read_excel(EK4A, skiprows=1)
e.columns = ['kamu_no','barkod','ad_e','eski_barkod','esdeger','terapotik','giris','aktif',
             'pasif','durum','i151','i100','i52','i52a','ozel','eczaci']
e['bar'] = e.barkod.astype(str).str.strip()
for c in ['i151','i100','i52','i52a']:
    e[c] = pd.to_numeric(e[c], errors='coerce')

d = t.merge(e[['bar','durum','esdeger','i151','i100','i52','i52a']], on='bar', how='left')

def band_orani(fiyat, r):
    """EK-4/A: iskonto orani depocuya satis fiyati bandina gore secilir."""
    if pd.isna(fiyat): return np.nan
    if fiyat >= 151.25:                  return r.i151
    if 100.38 <= fiyat <= 151.24:        return r.i100
    if 52.44  <= fiyat <= 100.37:        return r.i52
    return r.i52a

d['isk_ek4a'] = [band_orani(f, r) for f, r in zip(d.FIYATI, d.itertuples())]
d['isk_titck'] = pd.to_numeric(d.KURUM, errors='coerce') / 100
# TITCK KURUM sutunu urune ozgu fiilen uygulanan orandir; EK-4/A bandi capraz kontrol
d['isk'] = d.isk_titck.fillna(d.isk_ek4a)
d['uyum'] = (d.isk_titck.round(3) == d.isk_ek4a.round(3))

# SGK odemesi: firma satis fiyati uzerinden kamu kurum iskontosu
d['sgk_tl'] = (d.FIYATI * (1 - d.isk)).round(2)
# duyarlilik: dagitim kari + KDV eklenmis hali (perakende esasli ust sinir)
KAR, KDV = 0.219, 0.10
d['sgk_ust_tl'] = (d.FIYATI * (1 + KAR) * (1 + KDV) * (1 - d.isk)).round(2)

d[['bar','ad','FIYATI','isk','isk_titck','isk_ek4a','uyum','durum','sgk_tl','sgk_ust_tl']] \
    .to_csv('../../outputs/ilac_birim_maliyet.csv', index=False)

# ---------------------------------------------------------------- rejim hesabi
def bul(desen, n=1):
    m = d[d.ad.str.contains(desen, case=False, regex=True, na=False)]
    return m.nsmallest(n, 'sgk_tl') if len(m) else None

def flakon_coz(hedef_mg, secenekler):
    """Hedef dozu en ucuz flakon kombinasyonuyla karsila (artik dahil, eksik yok)."""
    en_iyi = None
    for a in range(0, math.ceil(hedef_mg / secenekler[0][0]) + 2):
        for b in range(0, math.ceil(hedef_mg / secenekler[-1][0]) + 2):
            mg = a * secenekler[0][0] + b * secenekler[-1][0]
            if mg < hedef_mg: continue
            tl = a * secenekler[0][1] + b * secenekler[-1][1]
            if en_iyi is None or tl < en_iyi[0]: en_iyi = (tl, a, b, mg)
    return en_iyi

BSA = 1.7
print('=' * 78); print('SGK ODEMESI — kamu kurum iskontosu uygulanmis (BSA %.1f m2)' % BSA)
print('=' * 78)

ILAC = {
 'Doksorubisin 50 mg': 'ADRIMISIN 50|DOXORUBICIN 50',
 'Doksorubisin 10 mg': 'ADRIMISIN 10|DOXORUBICIN 10',
 'Siklofosfamid 500 mg': 'ENDOXAN 500',
 'Dosetaksel 80 mg': 'TADOCEL 80',
 'Dosetaksel 20 mg': 'TADOCEL 20',
 'Karboplatin 450 mg': 'CARBOPLATIN KOCAK 450',
 'Filgrastim 30 MIU 5 enj': 'ZARZIO 30|LEUCOSTIM 30',
 'Pegfilgrastim 6 mg': 'ZIEXTENZO 6|NEULASTIM 6',
 'Aprepitant 125/80/80': 'EMEND 125',
 'Granisetron 3 mg amp': 'KYTRIL 3 MG 3 ML 1',
 'Ondansetron 8 mg amp': 'ZOFER 8 MG 4 ML',
 'Deksametazon 8 mg amp': 'DEKORT 8 MG 2 ML|DEKSAMETAZON-PF 8',
 'Tamoksifen 20 mg 30 tb': 'TAMOXIFEN 20',
 'Anastrozol 1 mg 28 tb': 'ARIMIDEX 1 MG',
 'Letrozol 2,5 mg 30 tb': 'FEMARA 2.5',
 'Eksemestan 25 mg 30 drj': 'AROMASIN 25',
 'Trastuzumab 150 mg': 'HERCEPTIN 150|HERZUMA 150',
}
print(f'{"Ilac":26s} {"liste TL":>10s} {"isk":>6s} {"SGK TL":>10s}  durum')
bulunan = {}
for ad, des in ILAC.items():
    r = bul(des)
    if r is None or not len(r):
        print(f'{ad:26s} {"— BULUNAMADI":>10s}'); continue
    r = r.iloc[0]; bulunan[ad] = r.sgk_tl
    print(f'{ad:26s} {r.FIYATI:10.2f} {r.isk:6.1%} {r.sgk_tl:10.2f}  {r.durum}')

print('\n' + '=' * 78); print('AC -> DOSETAKSEL REJIMI (4+4 kur)'); print('=' * 78)
if all(k in bulunan for k in ['Doksorubisin 50 mg','Doksorubisin 10 mg','Siklofosfamid 500 mg',
                              'Dosetaksel 80 mg','Dosetaksel 20 mg']):
    dox = flakon_coz(60*BSA, [(50, bulunan['Doksorubisin 50 mg']), (10, bulunan['Doksorubisin 10 mg'])])
    sik = flakon_coz(600*BSA, [(500, bulunan['Siklofosfamid 500 mg']), (500, bulunan['Siklofosfamid 500 mg'])])
    tax = flakon_coz(75*BSA, [(80, bulunan['Dosetaksel 80 mg']), (20, bulunan['Dosetaksel 20 mg'])])
    print(f'  Doksorubisin {60*BSA:.0f} mg -> {dox[1]}×50 + {dox[2]}×10 mg = {dox[3]} mg : {dox[0]:10.2f} TL')
    print(f'  Siklofosfamid {600*BSA:.0f} mg -> {math.ceil(600*BSA/500)}×500 mg          : '
          f'{math.ceil(600*BSA/500)*bulunan["Siklofosfamid 500 mg"]:10.2f} TL')
    print(f'  Dosetaksel {75*BSA:.1f} mg -> {tax[1]}×80 + {tax[2]}×20 mg = {tax[3]} mg : {tax[0]:10.2f} TL')
    ac = dox[0] + math.ceil(600*BSA/500) * bulunan['Siklofosfamid 500 mg']
    print(f'\n  AC kur basi        {ac:10.2f} TL  × 4 = {4*ac:11.2f} TL')
    print(f'  Dosetaksel kur basi{tax[0]:10.2f} TL  × 4 = {4*tax[0]:11.2f} TL')
    print(f'  ILAC TOPLAMI                        = {4*ac+4*tax[0]:11.2f} TL')
    UYG, PORT, HEM, MUA = 1482.02, 1321.34, 33.27, 356.00
    ek = 8*UYG + PORT + 8*HEM + 8*MUA
    print(f'  + uygulama/port/izlem (SUT)         = {ek:11.2f} TL')
    print(f'  TOPLAM (destek ilaclari haric)      = {4*ac+4*tax[0]+ek:11.2f} TL')

print('\nNOT: SGK odemesi firma satis fiyati uzerinden kamu kurum iskontosu ile hesaplandi;')
print('     dagitim kari ve KDV eklenmemistir (muhafazakar). Ust sinir icin sgk_ust_tl sutunu.')
print('TITCK KURUM ile EK-4/A band orani uyumu: %.1f%%' % (100*d.uyum.mean()))
