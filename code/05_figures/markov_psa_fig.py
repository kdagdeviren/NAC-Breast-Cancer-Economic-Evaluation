"""MARKOV MODELI + OLASILIKSAL DUYARLILIK ANALIZI
S2 (karar kurali) vs S1 (alt tip kurali) vs S0 (standart bakim)

MODEL MANTIGI — kritik:
  Karar kurali RCB'yi DEGISTIRMEZ, ongorur. Dolayisiyla hastanin temel prognozu
  (uzak nuks riski) iki kolda AYNIDIR. EBCTCG 2018: neoadjuvan ile adjuvan kemoterapi
  arasinda uzak nuks, meme kanseri mortalitesi ve genel sagkalim farki YOKTUR.
  Iki kol arasindaki gercek farklar:
    1. Lokorejyonel nuks: neoadjuvan kolda 15 yilda mutlak +%5,5 (EBCTCG)
    2. Cerrahi kapsami: cerrahi-once kolunda daha az meme koruyucu cerrahi,
       daha cok aksiller diseksiyon -> daha cok lenfodem
    3. Kemoterapiden kacinma: isaretlenen hastalarin %11,4'u (kohorttan)
    4. Tedavi sirasi ve maliyet zamanlamasi
"""
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'ps.fonttype':42,'pdf.fonttype':42,'savefig.dpi':600,'figure.dpi':600,'font.size':9,'font.family':'DejaVu Sans',
                     'axes.spines.top':False,'axes.spines.right':False})
rng = np.random.default_rng(2026)
O = '../../figures/'
C = {'ana':'#3d5a80','iyi':'#2a9d8f','kotu':'#e76f51','notr':'#adb5bd','vurgu':'#f4a261'}

# ================================================================= PARAMETRELER
BAZ = dict(
    # --- karar kurali (Kohort v16, esik 0,80, hasta-gruplu CV, bootstrap GA'lari)
    kohort=320, uygun=188,
    s2_n=108, s2_dogru=98, s2_n_ga=(92,125),
    s1_n=61,  s1_dogru=55,
    s3_n=118, s3_dogru=107,   # ek analiz: S1 ∪ S2 (post-hoc)
    # --- kisa donem maliyet (TL, SGK, SUT EK-2 29.06.2026 + TITCK/EK-4A)
    c_nak_sistemik=77262.65,        # AC→taksan + kapesitabin (%10 hasta, CREATE-X)
    c_adjuvan_sistemik=72134.29,    # %70 AC→taksan + %30 TC (dusuk riskli adjuvan secenegi)
    c_cerrahi_nak=54784.09,        # MKC %34,5
    c_cerrahi_once=58512.31,        # MKC %27,0
    c_karar_destegi=155.00,
    c_rt=32407.80,
    # --- uzun donem yillik maliyet
    c_izlem=1200.00,               # yillik muayene + mamografi (SUT)
    c_lrr=180000.00,               # lokorejyonel nuks epizodu (cerrahi+RT+sistemik) ⬛
    c_met_yillik=250000.00,        # metastatik hastalik yillik ⬛ genis PSA araligi
    c_lenfodem_yillik=8000.00,     # konservatif tedavi, kompresyon ⬛
    # --- klinik
    adjuvan_kt_orani=0.886,
    p_lenfodem_slnb=0.06, p_lenfodem_alnd=0.25,
    # meme cerrahisi: HR+/HER2- ozgu NAK sonrasi MKC %34,5 (Z1071);
    # cerrahi-once EBCTCG BCS oran orani (RR 1,28) ile turetildi: 34,5/1,28 = %27,0
    mkc_nak=0.345, mkc_once=0.270,
    # aksiller cerrahi: NCDB (n=792.581, ER+/HER2-) NAK kolunda ALND olasiligi
    # DAHA YUKSEK (OR 1,18). Onceki varsayim tersineydi ve duzeltildi.
    alnd_once=0.45, alnd_nak=0.491,
    lrr_ek_neoadjuvan=0.055,       # EBCTCG 15 yil mutlak fark
    efs10_rcb01=0.835, efs10_rcb23=0.635,   # Yau 2022, HR+/HER2-
    pay_uzak=0.70,                 # olaylarin uzak nuks payi
    # --- fayda
    u_dfs=0.88, u_kt=0.86, u_ilk_yil=0.84, u_lrr=0.73, u_met=0.73,
    du_lenfodem=0.04, du_mastektomi=0.04, du_fn=0.150,
    kt_suresi=0.46, fn_insidans=0.15, fn_sure=12/365,
    # --- model
    ufuk=40, indirim=0.03, yas=50,
)

def yillik_oran(s10):
    """10 yillik sagkalimdan yillik olay oranina"""
    return -np.log(s10)/10

def markov(P, kol, rcb23_orani):
    """Yasam boyu maliyet ve QALY.

    KRITIK MODEL KURALI (EBCTCG 2018):
      Uzak nuks ve mortalite IKI KOLDA AYNIDIR. Karar kurali prognozu degistirmez.
      Iki kol arasindaki tek uzun donem farki:
        - lokorejyonel nuks sikligi (neoadjuvan kolda daha yuksek)
        - lenfodem (cerrahi-once kolunda daha yuksek, aksiller diseksiyon nedeniyle)
      Lokorejyonel nuks GECICI bir durumdur: bir yil surer, maliyet ve fayda kaybi
      dogurur, ancak uzak nukse ilerlemeyi ve mortaliteyi ARTIRMAZ. Aksi halde
      iki kolun sagkalimi ayrisir ve bu EBCTCG bulgusuyla celisirdi.
    """
    n = P['ufuk']; disk = 1/(1+P['indirim'])**np.arange(n)
    r = rcb23_orani*yillik_oran(P['efs10_rcb23']) + (1-rcb23_orani)*yillik_oran(P['efs10_rcb01'])
    p_uzak = 1-np.exp(-r*P['pay_uzak'])            # IKI KOLDA AYNI
    p_lrr_baz = 1-np.exp(-r*(1-P['pay_uzak']))
    # EBCTCG: 15 yilda mutlak +%5,5 -> yillik goreli fazla
    ek = P['lrr_ek_neoadjuvan']/15
    p_lrr = p_lrr_baz + (ek if kol=='nak' else 0)
    p_met_olum = 0.28; p_fon = 0.008
    # lenfodem — aksiller cerrahi tipine bagli, kalici
    alnd = P['alnd_nak'] if kol=='nak' else P['alnd_once']
    p_lenf = alnd*P['p_lenfodem_alnd'] + (1-alnd)*P['p_lenfodem_slnb']
    # durumlar: 0 DFS, 1 LRR (gecici, 1 yil), 2 MET, 3 OLUM
    d = np.zeros((n,4)); d[0] = [1,0,0,0]
    mal = np.zeros(n); qal = np.zeros(n)
    for t in range(n):
        dfs,lrr,met,olu = d[t]
        mal[t] = dfs*(P['c_izlem'] + p_lenf*P['c_lenfodem_yillik']) \
               + lrr*(P['c_lrr'] + P['c_izlem']) + met*P['c_met_yillik']
        qal[t] = dfs*(P['u_dfs'] - p_lenf*P['du_lenfodem']) \
               + lrr*(P['u_lrr'] - p_lenf*P['du_lenfodem']) + met*P['u_met']
        if t == n-1: break
        y = np.zeros(4)
        y[1] = dfs*p_lrr                                   # DFS -> LRR
        y[2] = (dfs+lrr)*p_uzak                            # uzak nuks: IKI KOLDA AYNI oran
        y[3] = olu + (dfs+lrr+met)*p_fon + met*p_met_olum
        y[0] = dfs*(1-p_lrr-p_uzak-p_fon) + lrr*(1-p_uzak-p_fon)   # LRR bir yil sonra DFS'e doner
        y[2] += met*(1-p_met_olum-p_fon)
        d[t+1] = np.clip(y,0,None); d[t+1] /= d[t+1].sum()
    return (mal*disk).sum(), (qal*disk).sum()

def kisa_donem(P, kol):
    """Ilk 12 ay: sistemik tedavi + cerrahi + RT"""
    if kol == 'nak':
        m = P['c_nak_sistemik'] + P['c_cerrahi_nak'] + P['mkc_nak']*P['c_rt']
        q = P['u_kt']*P['kt_suresi'] + P['u_ilk_yil']*(1-P['kt_suresi'])
        q -= P['du_fn']*P['fn_sure']*P['fn_insidans']
        q -= P['du_mastektomi']*(1-P['mkc_nak'])
    else:
        kt = P['adjuvan_kt_orani']
        m = P['c_cerrahi_once'] + kt*P['c_adjuvan_sistemik'] + P['mkc_once']*P['c_rt']
        q = kt*(P['u_kt']*P['kt_suresi'] + P['u_ilk_yil']*(1-P['kt_suresi'])) + (1-kt)*P['u_ilk_yil']
        q -= P['du_fn']*P['fn_sure']*P['fn_insidans']*kt
        q -= P['du_mastektomi']*(1-P['mkc_once'])
    return m, q

def strateji(P, n_isaret, karar_maliyeti):
    """Kohort duzeyinde toplam maliyet ve QALY"""
    n = P['kohort']; rcb23 = 0.622
    m_nak, q_nak = kisa_donem(P,'nak');   M_nak, Q_nak = markov(P,'nak',rcb23)
    m_onc, q_onc = kisa_donem(P,'once');  M_onc, Q_onc = markov(P,'once',rcb23)
    mal = n_isaret*(m_onc+M_onc) + (n-n_isaret)*(m_nak+M_nak) + karar_maliyeti*n
    qal = n_isaret*(q_onc+Q_onc) + (n-n_isaret)*(q_nak+Q_nak)
    return mal, qal

def calistir(P):
    return {'S0': strateji(P,0,0),
            'S1': strateji(P,P['s1_n'],0),
            'S2': strateji(P,P['s2_n'],P['c_karar_destegi']),
            'S3': strateji(P,P['s3_n'],P['c_karar_destegi'])}

# ================================================================= 1. TABAN
print('='*76); print('TABAN SENARYO — yaşam boyu ufuk, %3 indirim'); print('='*76)
baz = calistir(BAZ)
print(f'{"Strateji":8s} {"Maliyet (kohort)":>18s} {"QALY":>10s} {"TRY per patient":>15s} {"QALY per patient":>16s}')
for k,(m,q) in baz.items():
    print(f'{k:8s} {m:18,.0f} {q:10.1f} {m/BAZ["kohort"]:15,.0f} {q/BAZ["kohort"]:16.3f}')
dm = baz['S2'][0]-baz['S1'][0]; dq = baz['S2'][1]-baz['S1'][1]
print(f'\nS2 − S1: Δmaliyet {dm:+,.0f} TL · ΔQALY {dq:+.3f}')
print(f'         hasta başı: {dm/BAZ["kohort"]:+,.0f} TL · {dq/BAZ["kohort"]:+.5f} QALY')
if dq!=0: print(f'         ICER {dm/dq:+,.0f} TL/QALY')

# ================================================================= 2. PSA
print('\n'+'='*76); print('OLASILIKSAL DUYARLILIK ANALİZİ (10.000 yineleme)'); print('='*76)
N = 10000
def gamma_par(ort, sd):
    k = (ort/sd)**2; th = sd**2/ort; return k, th
def beta_par(ort, sd):
    v = sd**2; a = ort*(ort*(1-ort)/v - 1); b = a*(1-ort)/ort; return max(a,0.1), max(b,0.1)

res = {'S0':[], 'S1':[], 'S2':[], 'S3':[]}
for i in range(N):
    P = dict(BAZ)
    for k,sd_orani in [('c_nak_sistemik',.15),('c_adjuvan_sistemik',.15),
                       ('c_cerrahi_nak',.20),('c_cerrahi_once',.20),('c_rt',.15),
                       ('c_lrr',.35),('c_met_yillik',.40),('c_lenfodem_yillik',.40),
                       ('c_izlem',.30),('c_karar_destegi',.50)]:
        a,b = gamma_par(BAZ[k], BAZ[k]*sd_orani); P[k] = rng.gamma(a,b)
    for k,sd in [('u_dfs',.03),('u_kt',.04),('u_ilk_yil',.04),('u_lrr',.05),('u_met',.05)]:
        a,b = beta_par(BAZ[k], sd); P[k] = rng.beta(a,b)
    for k,sd in [('du_lenfodem',.015),('du_mastektomi',.015),('du_fn',.03)]:
        P[k] = max(rng.normal(BAZ[k], sd), 0)
    for k,sd in [('adjuvan_kt_orani',.06),('p_lenfodem_alnd',.05),('p_lenfodem_slnb',.02),
                 ('mkc_nak',.05),('mkc_once',.05),('fn_insidans',.04),
                 ('alnd_nak',.06),('alnd_once',.06),
                 ('efs10_rcb01',.03),('efs10_rcb23',.04)]:
        a,b = beta_par(BAZ[k], sd); P[k] = rng.beta(a,b)
    P['lrr_ek_neoadjuvan'] = max(rng.normal(0.055, 0.016), 0)   # EBCTCG %95 GA 2,4–8,6
    P['s2_n'] = int(np.clip(rng.normal(108, 8.4), 60, 188))
    out = calistir(P)
    for k in res: res[k].append(out[k])

df = {k: np.array(v) for k,v in res.items()}
for k in ['S0','S1','S2','S3']:
    m,q = df[k][:,0], df[k][:,1]
    print(f'  {k}: maliyet {m.mean():12,.0f} ({np.percentile(m,2.5):,.0f}–{np.percentile(m,97.5):,.0f})'
          f'  QALY {q.mean():7.1f} ({np.percentile(q,2.5):.1f}–{np.percentile(q,97.5):.1f})')

dmv = df['S2'][:,0]-df['S1'][:,0]; dqv = df['S2'][:,1]-df['S1'][:,1]
d3m = df['S3'][:,0]-df['S2'][:,0]; d3q = df['S3'][:,1]-df['S2'][:,1]
print(f'\n  S3−S2 Δmaliyet {d3m.mean():+,.0f} ({np.percentile(d3m,2.5):+,.0f} – {np.percentile(d3m,97.5):+,.0f})')
print(f'  S3−S2 ΔQALY    {d3q.mean():+.3f} ({np.percentile(d3q,2.5):+.3f} – {np.percentile(d3q,97.5):+.3f})')
print(f'  S3 daha ucuz olma olasılığı (S2 ye göre): %{100*(d3m<0).mean():.1f}')
print(f'\n  S2−S1 Δmaliyet {dmv.mean():+,.0f} ({np.percentile(dmv,2.5):+,.0f} – {np.percentile(dmv,97.5):+,.0f})')
print(f'  S2−S1 ΔQALY    {dqv.mean():+.3f} ({np.percentile(dqv,2.5):+.3f} – {np.percentile(dqv,97.5):+.3f})')
print(f'  S2 daha ucuz olma olasılığı: %{100*(dmv<0).mean():.1f}')
print(f'  S2 daha fazla QALY olasılığı: %{100*(dqv>0).mean():.1f}')

# ================================================================= 3. CEAC
esikler = np.linspace(0, 3_000_000, 121)
ceac = {k: [] for k in ['S0','S1','S2','S3']}
ceac3 = {k: [] for k in ['S0','S1','S2']}   # birincil: yalnizca onceden tanimli stratejiler
for lam in esikler:
    nmb = np.column_stack([lam*df[k][:,1]-df[k][:,0] for k in ['S0','S1','S2','S3']])
    kazanan = nmb.argmax(axis=1)
    for i,k in enumerate(['S0','S1','S2','S3']): ceac[k].append((kazanan==i).mean())
    nmb3 = nmb[:,:3]; kaz3 = nmb3.argmax(axis=1)
    for i,k in enumerate(['S0','S1','S2']): ceac3[k].append((kaz3==i).mean())
print('\n  Kabul edilebilirlik (S2 en iyi olma olasılığı):')
for lam in [0, 300_000, 600_000, 900_000, 1_800_000]:
    i = np.argmin(abs(esikler-lam))
    print(f'    λ={lam:>9,} -> [4 strateji] S3 %{100*ceac["S3"][i]:5.1f} · S2 %{100*ceac["S2"][i]:5.1f} · S1 %{100*ceac["S1"][i]:5.1f} · S0 %{100*ceac["S0"][i]:5.1f}')
    print(f'    {"":>13}   [3 strateji] S2 %{100*ceac3["S2"][i]:5.1f} · S1 %{100*ceac3["S1"][i]:5.1f} · S0 %{100*ceac3["S0"][i]:5.1f}')

# ================================================================= 4. SEKILLER
fig,ax = plt.subplots(figsize=(4.50,3.38))
ax.scatter(dqv/BAZ['kohort'], dmv/BAZ['kohort'], s=3, alpha=.18, color=C['ana'], edgecolors='none')
ax.axhline(0,color='#495057',lw=.9); ax.axvline(0,color='#495057',lw=.9)
x = np.linspace(dqv.min()/BAZ['kohort'], dqv.max()/BAZ['kohort'], 10)
ax.plot(x, 600_000*x, ls='--', color=C['kotu'], lw=1.2, label='λ = 600.000 TL/QALY')
ax.scatter([dq/BAZ['kohort']],[dm/BAZ['kohort']], s=60, color=C['vurgu'],
           edgecolors='#22333b', zorder=5, label='taban senaryo')
ax.set_xlabel('Incremental QALY (per patient)'); ax.set_ylabel('Incremental cost (per patient, TRY)')
ax.legend(fontsize=7.5, frameon=False, loc='upper left')
plt.tight_layout(); plt.savefig(O+'Fig10.png',bbox_inches='tight'); plt.savefig(O+'Fig10.eps',bbox_inches='tight',format='eps'); plt.savefig(O+'Fig10.pdf',bbox_inches='tight'); plt.close()

fig,ax = plt.subplots(figsize=(4.50,2.89))
for k,c,ls,ad in [('S2',C['vurgu'],'-','S2 prediction rule'),
                  ('S1',C['iyi'],'--','S1 subtype rule'),('S0',C['notr'],'-.','S0 standard care')]:
    ax.plot(esikler/1000, ceac3[k], lw=1.8, color=c, ls=ls, label=ad)
ax.axvline(600, ls=':', lw=1.0, color='#6c757d')
ax.set_xlabel('Threshold (thousand TRY per QALY)'); ax.set_ylabel('Probability of being the optimal strategy')
ax.set_ylim(0,1); ax.legend(fontsize=7.5, frameon=False)
plt.tight_layout(); plt.savefig(O+'Fig8.png',bbox_inches='tight'); plt.savefig(O+'Fig8.eps',bbox_inches='tight',format='eps'); plt.close()

# ================================================================= 5. TORNADO
print('\n'+'='*76); print('TEK YÖNLÜ DUYARLILIK (tornado)'); print('='*76)
TEST = [('adjuvan_kt_orani',0.70,1.00,'Adjuvant chemotherapy rate'),
        ('c_nak_sistemik',60000,100000,'Neoadjuvant systemic therapy cost'),
        ('mkc_once',0.20,0.40,'Breast conservation rate, surgery-first arm'),
        ('mkc_nak',0.28,0.45,'Breast conservation rate, neoadjuvant arm'),
        ('alnd_nak',0.35,0.60,'Axillary dissection rate, neoadjuvant arm'),
        ('alnd_once',0.35,0.60,'Axillary dissection rate, surgery-first arm'),
        ('c_lrr',110000,250000,'Locoregional recurrence cost'),
        ('p_lenfodem_alnd',0.15,0.35,'Lymphoedema after axillary dissection'),
        ('du_mastektomi',0.01,0.08,'Mastectomy disutility'),
        ('c_met_yillik',150000,400000,'Metastatic disease, annual cost'),
        ('lrr_ek_neoadjuvan',0.024,0.086,'Excess locoregional recurrence'),
        ('c_karar_destegi',50,500,'Decision support cost')]
tor = []
for k,lo,hi,ad in TEST:
    v = []
    for x in [lo,hi]:
        P = dict(BAZ); P[k] = x; o = calistir(P)
        v.append((o['S2'][0]-o['S1'][0])/BAZ['kohort'])
    tor.append((ad, v[0], v[1], abs(v[1]-v[0])))
tor.sort(key=lambda z: z[3])
print(f'{"Parametre":34s} {"alt":>12s} {"üst":>12s} {"aralık":>10s}')
for ad,a,b,g in reversed(tor): print(f'{ad:34s} {a:+12,.0f} {b:+12,.0f} {g:10,.0f}')

fig,ax = plt.subplots(figsize=(6.50,3.71))
y = np.arange(len(tor)); taban = dm/BAZ['kohort']
for i,(ad,a,b,g) in enumerate(tor):
    ax.barh(i, b-taban, left=taban, height=.6, color=C['kotu'], alpha=.75)
    ax.barh(i, a-taban, left=taban, height=.6, color=C['iyi'], alpha=.75)
ax.axvline(taban, color='#22333b', lw=1.2)
ax.set_yticks(y); ax.set_yticklabels([t[0] for t in tor], fontsize=8)
ax.set_xlabel('S2 − S1 incremental cost (per patient, TRY)')
plt.tight_layout(); plt.savefig(O+'Fig9.png',bbox_inches='tight'); plt.savefig(O+'Fig9.eps',bbox_inches='tight',format='eps'); plt.close()

pd.DataFrame({'esik':esikler,'S0':ceac['S0'],'S1':ceac['S1'],'S2':ceac['S2'],'S3':ceac['S3']}) \
  .to_csv('../../outputs/kabul_edilebilirlik_egrisi_S3dahil.csv', index=False)
pd.DataFrame({'esik':esikler,'S0':ceac3['S0'],'S1':ceac3['S1'],'S2':ceac3['S2']}) \
  .to_csv('../../outputs/kabul_edilebilirlik_egrisi.csv', index=False)
# ESM sekli: S3 dahil dort strateji
fig,ax=plt.subplots(figsize=(4.5,3.4))
for k,c,ad in [('S3','#7b2cbf','S3 combined (post hoc)'),('S2',C['vurgu'],'S2 prediction rule'),
               ('S1',C['iyi'],'S1 subtype rule'),('S0',C['notr'],'S0 standard care')]:
    ax.plot(esikler/1000, ceac[k], lw=2.2, color=c, label=ad)
ax.axvline(600, ls=':', lw=1.2, color='#6c757d')
ax.set_xlabel('Threshold (thousand TRY per QALY)'); ax.set_ylabel('Probability of being the optimal strategy')
ax.set_ylim(0,1); ax.legend(fontsize=8, frameon=False, loc='center right')
plt.tight_layout(); plt.savefig('../../figures/ESM_Fig_S3_CEAC.png',bbox_inches='tight'); plt.savefig('../../figures/ESM_Fig_S3_CEAC.pdf',bbox_inches='tight',facecolor='white'); plt.savefig('../../figures/ESM_Fig_S3_CEAC.eps',bbox_inches='tight',format='eps'); plt.close()
pd.DataFrame(tor, columns=['parametre','alt','ust','aralik']) \
  .to_csv('../../outputs/tornado.csv', index=False)
print('\n3 şekil ve 2 CSV üretildi.')
