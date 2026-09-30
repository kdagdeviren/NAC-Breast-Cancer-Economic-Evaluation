"""KARAR AGACI + ICER — S2 (karar kurali) vs S1 (alt tip kurali) vs S0 (standart bakim)
Kohort v16, cM0, esik 0,80. Tum birim fiyatlar SUT EK-2 (29.06.2026) ve
TITCK/EK-4A kamu kurum iskontosu uygulanmis degerlerdir.

UYARI: Bu ilk surumdur. Isaretli ⬛ parametreler henuz kohorttan/kayittan alinamamistir
ve varsayim olarak girilmistir; duyarlilik analizinde genis araliklarla taranir.
"""
import numpy as np, pandas as pd

# ---------------------------------------------------------------- 1. MODEL PARAMETRELERI
P = dict(
    # --- karar kuralinin cikitlari (Kohort v16, esik 0,80, hasta-gruplu CV)
    kohort=320, uygun=179, isaretlenen=105, dogru=95, yanlis=10,
    s1_isaretlenen=59, s1_dogru=53, s1_yanlis=6,
    # --- maliyet (TL, SGK perspektifi)
    nak_yolu=79448.62,              # ilac + destek + uygulama + FN beklenen
    adjuvan_kt=79448.62,            # ayni rejim adjuvan verilirse (port hari� ~ayni)
    cerrahi_nak=54784.09,           # MKC %34,5 varsayimi
    cerrahi_once=59506.19,          # MKC %25 varsayimi
    karar_destegi=155.00,           # veri girisi is gucu, hasta basi ⬛
    # --- klinik parametreler
    adjuvan_kt_orani=0.886,         # 105 isaretlenenin 93'u yine KT alir (12 kacinir)
    # --- fayda (utility) — Rautenberg 2021
    u_hastaliksiz=0.88, u_kemoterapi=0.86, u_ilk_yil=0.84,
    kt_suresi_yil=0.46,             # 8 kur x 21 gun ~ 168 gun
    # --- disutility
    du_fn=0.150, fn_sure_yil=12/365, fn_insidans=0.15,
    du_mastektomi=0.04,             # MKC 0,79 vs mastektomi 0,75
    # --- ufuk ve indirim
    ufuk_yil=1.0, indirim=0.03,
)

def yaz(b): print('\n'+'='*76); print(b); print('='*76)

# ---------------------------------------------------------------- 2. KOL MALIYETLERI
def hasta_maliyeti(kol, dogru_mu):
    """kol: 'nak' veya 'cerrahi_once'. dogru_mu: model ongorusu dogru cikti mi."""
    if kol == 'nak':
        return P['nak_yolu'] + P['cerrahi_nak']
    m = P['cerrahi_once']
    m += P['adjuvan_kt_orani'] * P['adjuvan_kt']
    return m

def hasta_qaly(kol):
    if kol == 'nak':
        q = P['u_kemoterapi']*P['kt_suresi_yil'] + P['u_ilk_yil']*(1-P['kt_suresi_yil'])
        q -= P['du_fn']*P['fn_sure_yil']*P['fn_insidans']
        return q
    kt = P['adjuvan_kt_orani']
    q = kt*(P['u_kemoterapi']*P['kt_suresi_yil'] + P['u_ilk_yil']*(1-P['kt_suresi_yil'])) \
        + (1-kt)*P['u_ilk_yil']
    q -= P['du_fn']*P['fn_sure_yil']*P['fn_insidans']*kt
    q -= P['du_mastektomi']*(0.345-0.25)      # cerrahi kapsam farki
    return q

yaz('1. HASTA BAŞINA MALİYET VE QALY (12 aylık ufuk)')
for kol in ['nak','cerrahi_once']:
    print(f'  {kol:14s} maliyet {hasta_maliyeti(kol,True):10.2f} TL   QALY {hasta_qaly(kol):.4f}')
fark_m = hasta_maliyeti('cerrahi_once',True) - hasta_maliyeti('nak',True)
fark_q = hasta_qaly('cerrahi_once') - hasta_qaly('nak')
print(f'\n  Yönlendirilen hasta başına: maliyet {fark_m:+10.2f} TL · QALY {fark_q:+.4f}')

# ---------------------------------------------------------------- 3. STRATEJI DUZEYI
yaz('2. STRATEJİ DÜZEYİNDE (320 kayıtlık kohort)')
def strateji(n_isaret, ad):
    """Isaretlenen hastalar cerrahi-once koluna gider, digerleri NAK alir."""
    n = P['kohort']
    mal = n_isaret*hasta_maliyeti('cerrahi_once',True) + (n-n_isaret)*hasta_maliyeti('nak',True)
    if ad != 'S0': mal += n*P['karar_destegi'] if ad == 'S2' else 0
    qal = n_isaret*hasta_qaly('cerrahi_once') + (n-n_isaret)*hasta_qaly('nak')
    return mal, qal

S = {}
S['S0 standart bakım']       = strateji(0, 'S0')
S['S1 alt tip kuralı']       = strateji(P['s1_isaretlenen'], 'S1')
S['S2 karar kuralı']         = strateji(P['isaretlenen'], 'S2')
print(f'{"Strateji":26s} {"Toplam maliyet":>16s} {"Toplam QALY":>13s} {"Hasta başı":>12s}')
for ad,(m,q) in S.items():
    print(f'{ad:26s} {m:16,.0f} {q:13.2f} {m/P["kohort"]:12,.0f}')

# ---------------------------------------------------------------- 4. ICER
yaz('3. ARTIMSAL MALİYET-ETKİNLİK ORANI (ICER)')
sira = sorted(S.items(), key=lambda kv: kv[1][0])
onceki = None
for ad,(m,q) in sira:
    if onceki is None:
        print(f'  {ad:26s} referans (en düşük maliyet)')
    else:
        dm, dq = m-onceki[1][0], q-onceki[1][1]
        if dq > 0:
            print(f'  {ad:26s} Δmaliyet {dm:+12,.0f} · ΔQALY {dq:+7.3f} · ICER {dm/dq:12,.0f} TL/QALY')
        elif dq < 0 and dm < 0:
            print(f'  {ad:26s} Δmaliyet {dm:+12,.0f} · ΔQALY {dq:+7.3f} · daha ucuz ama daha az QALY')
        else:
            print(f'  {ad:26s} Δmaliyet {dm:+12,.0f} · ΔQALY {dq:+7.3f} · DOMİNE EDİLİYOR')
    onceki = (ad,(m,q))

yaz('4. S2 vs S1 — DOĞRUDAN KARŞILAŞTIRMA (asıl soru)')
m2,q2 = S['S2 karar kuralı']; m1,q1 = S['S1 alt tip kuralı']
dm, dq = m2-m1, q2-q1
print(f'  Δ maliyet (kohort)        {dm:+14,.0f} TL')
print(f'  Δ maliyet (hasta başı)    {dm/P["kohort"]:+14,.0f} TL')
print(f'  Δ QALY (kohort)           {dq:+14.3f}')
print(f'  Δ QALY (hasta başı)       {dq/P["kohort"]:+14.5f}')
if dq != 0:
    icer = dm/dq
    print(f'\n  ICER = {icer:,.0f} TL/QALY')
    if dm < 0 and dq > 0: print('  -> DOMİNANT: daha ucuz ve daha fazla QALY')
    elif dm < 0 and dq < 0: print('  -> daha ucuz ama QALY kaybı; eşik altında kalırsa kabul edilebilir')
    elif dm > 0 and dq > 0: print('  -> ek maliyetle ek fayda; eşikle karşılaştırılmalı')
    else: print('  -> DOMİNE EDİLİYOR: hem pahalı hem daha az QALY')

# ---------------------------------------------------------------- 5. DUYARLILIK
yaz('5. TEK YÖNLÜ DUYARLILIK — adjuvan kemoterapi oranı')
print(f'{"adjuvan KT oranı":>18s} {"Δmaliyet (kohort)":>20s} {"hasta başı":>14s}')
saklı = P['adjuvan_kt_orani']
for o in [0.5,0.6,0.7,0.8,0.886,0.95,1.0]:
    P['adjuvan_kt_orani']=o
    m2b,_ = strateji(P['isaretlenen'],'S2'); m1b,_ = strateji(P['s1_isaretlenen'],'S1')
    m2b += P['kohort']*P['karar_destegi']
    print(f'{o:18.3f} {m2b-m1b:+20,.0f} {(m2b-m1b)/P["kohort"]:+14,.0f}')
P['adjuvan_kt_orani']=saklı

yaz('6. TEK YÖNLÜ DUYARLILIK — cerrahi kapsam farkı')
print(f'{"cerrahi fark (TL)":>18s} {"Δmaliyet (kohort)":>20s} {"hasta başı":>14s}')
saklı2 = P['cerrahi_once']
for fark in [0, 3479, 4722, 5965, 10000]:
    P['cerrahi_once']=P['cerrahi_nak']+fark
    m2b,_ = strateji(P['isaretlenen'],'S2'); m1b,_ = strateji(P['s1_isaretlenen'],'S1')
    m2b += P['kohort']*P['karar_destegi']
    print(f'{fark:18,d} {m2b-m1b:+20,.0f} {(m2b-m1b)/P["kohort"]:+14,.0f}')
P['cerrahi_once']=saklı2

yaz('7. ULUSAL ÖLÇEKLEME (yıllık)')
YENI=25080; NAK_ORANI=0.20; UYGUN=179/320; ISARET=105/179; BENIMSEME=0.30
n = YENI*NAK_ORANI
print(f'  Yıllık yeni meme kanseri tanısı            {YENI:>10,}')
print(f'  NAK adayı (%20)                            {n:>10,.0f}')
print(f'  Uygun popülasyon (HR+/HER2−, %{100*UYGUN:.1f})       {n*UYGUN:>10,.0f}')
print(f'  Karar kuralınca işaretlenen (%{100*ISARET:.1f})       {n*UYGUN*ISARET:>10,.0f}')
print(f'  Benimseme %30 sonrası fiilen yönlendirilen {n*UYGUN*ISARET*BENIMSEME:>10,.0f}')
hasta_basi = dm/P['kohort']
print(f'\n  Hasta başına etki {hasta_basi:+,.0f} TL')
print(f'  YILLIK BÜTÇE ETKİSİ {n*UYGUN*ISARET*BENIMSEME*hasta_basi/(105/320):+,.0f} TL')
print('  (işaretlenen hasta başına etkiye çevrildi)')
