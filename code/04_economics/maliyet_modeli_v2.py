"""MALIYET MODELI v2 — Kohort v17
Yenilikler:
  1. Adjuvan kolda dusuk riskli hastalarda TC rejimi (4 kur dosetaksel+siklofosfamid)
  2. Neoadjuvan kolda reziduel hastalikta kapesitabin (CREATE-X)
  3. Cerrahi fiyatlarina yatis/anestezi/ameliyathane DAHIL oldugu dogrulandi (SUT 2.2.2.B-1)
  4. Ulusal olceklemede NAK oranı %10-%35 araligi
"""
import numpy as np, pandas as pd
BSA = 1.7

# ---------------------------------------------------------------- BIRIM FIYATLAR (SGK odemesi)
I = dict(dox50=371.75, dox10=201.19, siklo500=486.53,
         doce80=4690.08, doce40=2523.47, doce20=2713.07,
         palono=576.62, aprepitant=711.48, deksa_amp=43.84, deksa_tb=411.67,
         filgrastim=2589.83, pipTazo=258.89, kapesitabin_kutu=2448.68)
S = dict(uygulama=1482.02, port=1321.34, hemogram=33.27, muayene=356.00,
         alt_ast=23.23, eko=245.55, yatak=374.97, kan_kulturu=111.00, yayma=11.08,
         mkc_slnb=22226.451, mkc_alnd=22226.451, mrm=71932.77,
         rt_toplam=32407.80)

def yaz(b): print('\n'+'='*80); print(b); print('='*80)

# ---------------------------------------------------------------- REJIMLER
def ac_taksan():
    """4×AC + 4×dosetaksel — neoadjuvan ve yuksek riskli adjuvan standardi"""
    ac = 2*I['dox50'] + I['dox10'] + 3*I['siklo500']
    tax = 2*I['doce80']                       # 160 mg >= 127,5 mg, en ucuz kombinasyon
    return 4*ac + 4*tax, 8
def tc():
    """4×TC (dosetaksel + siklofosfamid) — dusuk riskli adjuvan secenegi"""
    kur = 2*I['doce80'] + 3*I['siklo500']
    return 4*kur, 4
def destek(n_kur):
    return n_kur*(I['palono'] + I['aprepitant'] + I['deksa_amp']) + (n_kur/2)*I['deksa_tb']/2
def uygulama(n_kur, port=True):
    return n_kur*(S['uygulama'] + S['hemogram'] + S['muayene'] + S['alt_ast']) \
         + (S['port'] if port else 0) + 2*S['eko']
def fn_epizod():
    return 7*S['yatak'] + 21*I['pipTazo'] + I['filgrastim'] + 2*S['kan_kulturu'] \
         + 7*(S['hemogram'] + S['yayma'])

FN_INS = 0.15
yaz('1. REJİM MALİYETLERİ (hasta başına, SGK ödemesi)')
for ad, (ilac, n) in [('AC→dosetaksel (8 kür)', ac_taksan()), ('TC (4 kür)', tc())]:
    d_ = destek(n); u_ = uygulama(n); f_ = FN_INS*fn_epizod()
    print(f'\n  {ad}')
    print(f'    sitotoksik ilaç       {ilac:11,.2f}')
    print(f'    destek tedavi         {d_:11,.2f}')
    print(f'    uygulama ve izlem     {u_:11,.2f}')
    print(f'    FN (beklenen, %15)    {f_:11,.2f}')
    print(f'    TOPLAM                {ilac+d_+u_+f_:11,.2f}')

AC_TOP = sum(ac_taksan()[:1]) + destek(8) + uygulama(8) + FN_INS*fn_epizod()
TC_TOP = tc()[0] + destek(4) + uygulama(4) + FN_INS*fn_epizod()

yaz('2. KAPESİTABİN (CREATE-X) — yalnızca neoadjuvan kolda doğar')
gun_mg = 2*1250*BSA
kur_tb = np.ceil(gun_mg*14/500)
kutu = kur_tb/120
print(f'  1250 mg/m² × 2 = {gun_mg:.0f} mg/gün · 14 gün = {gun_mg*14:.0f} mg/kür')
print(f'  500 mg tablet: {kur_tb:.0f} tablet/kür = {kutu:.2f} kutu')
print(f'  kutu bedeli {I["kapesitabin_kutu"]:,.2f} TL -> kür başı {kutu*I["kapesitabin_kutu"]:,.2f} TL')
KAPE_8 = 8*np.ceil(kutu)*I['kapesitabin_kutu']
print(f'  8 kür toplam: {KAPE_8:,.2f} TL')
print('  ⚠️ HR+/HER2−de CREATE-X faydası sınırlı; taban senaryoda uygulanan hasta oranı düşük')

yaz('3. İKİ KOLUN KARŞILAŞTIRMASI')
ADJ_KT = 0.886          # cerrahi-once kolunda adjuvan KT alan oran
TC_PAY = 0.30           # adjuvan KT alanlarin dusuk riskli olup TC alan orani ⬛
KAPE_PAY = 0.10         # NAK kolunda kapesitabin verilen reziduel hastalik orani ⬛
MKC_NAK, MKC_ONCE = 0.345, 0.270

nak_sistemik = AC_TOP + KAPE_PAY*KAPE_8
nak_cerrahi  = MKC_NAK*S['mkc_slnb'] + (1-MKC_NAK)*S['mrm']
nak_rt       = MKC_NAK*S['rt_toplam']
once_sistemik = ADJ_KT*((1-TC_PAY)*(AC_TOP-S['port']) + TC_PAY*(TC_TOP-S['port'])) + S['port']*ADJ_KT
once_cerrahi = MKC_ONCE*S['mkc_alnd'] + (1-MKC_ONCE)*S['mrm']
once_rt      = MKC_ONCE*S['rt_toplam']

print(f'{"":26s} {"NAK yolu":>14s} {"Cerrahi-önce":>14s} {"fark":>14s}')
for ad, a, b in [('sistemik tedavi', nak_sistemik, once_sistemik),
                 ('cerrahi (yatış dahil)', nak_cerrahi, once_cerrahi),
                 ('radyoterapi', nak_rt, once_rt)]:
    print(f'  {ad:24s} {a:14,.0f} {b:14,.0f} {b-a:+14,.0f}')
print(f'  {"TOPLAM":24s} {nak_sistemik+nak_cerrahi+nak_rt:14,.0f} '
      f'{once_sistemik+once_cerrahi+once_rt:14,.0f} '
      f'{(once_sistemik+once_cerrahi+once_rt)-(nak_sistemik+nak_cerrahi+nak_rt):+14,.0f}')

yaz('4. TC ORANI ve KAPESİTABİN ORANINA DUYARLILIK (12 aylık fark, TL)')
print(f'{"TC payı":>9s}', ''.join(f'{f"kape %{int(100*k)}":>14s}' for k in [0,0.05,0.10,0.20]))
for tcp in [0, 0.15, 0.30, 0.50]:
    sat = f'{tcp:9.0%}'
    for kp in [0, 0.05, 0.10, 0.20]:
        ns = AC_TOP + kp*KAPE_8
        os_ = ADJ_KT*((1-tcp)*(AC_TOP-S['port']) + tcp*(TC_TOP-S['port'])) + S['port']*ADJ_KT
        sat += f'{(os_+once_cerrahi+once_rt)-(ns+nak_cerrahi+nak_rt):14,.0f}'
    print(sat)

yaz('5. ULUSAL ÖLÇEKLEME — NAK oranı aralığı %10–%35')
YENI=25080; UYGUN=188/320; ISARET=108/188; TASARRUF=2526
print(f'{"NAK oranı":>10s} {"NAK adayı":>11s} {"uygun":>9s} {"işaretlenen":>12s} {"yıllık tasarruf (tam benimseme)":>34s}')
for nak in [0.10, 0.20, 0.35]:
    n=YENI*nak; u=n*UYGUN; i=u*ISARET
    print(f'{nak:10.0%} {n:11,.0f} {u:9,.0f} {i:12,.0f} {u*TASARRUF:34,.0f} TL')
print('\n  ⚠️ NAK oranı Türkiye için doğrulanmış bir kaynağa dayanmamaktadır;')
print('     ulusal kayıt %9,3 ile uzman görüşü %35 arasında geniş bir aralık bildirilmiştir.')
print('  ⚠️ Uygunluk oranı (%58,8) NAK ALMIŞ hastalardan oluşan kohortumuza aittir;')
print('     NAK ADAYI popülasyonda alt tip dağılımı farklı olabilir.')
