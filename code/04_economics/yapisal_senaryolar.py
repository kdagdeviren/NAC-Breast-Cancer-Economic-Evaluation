# -*- coding: utf-8 -*-
"""Hakem talebi: ekonomik modelin yapısal varsayımlarının stres testi."""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
s=open('../../outputs/markov_psa.py').read()
exec(s[:s.index('# ================================================================= 1. TABAN')])

# uzak metastaz esitligi varsayimini GEVSETEN markov surumu
import types
_markov_orj = markov
def markov_or(P, kol, rcb23, met_or_once=1.0):
    """cerrahi-önce kolunda uzak metastaz olasılığını met_or_once ile çarpar"""
    Q=dict(P)
    if kol=='once' and met_or_once!=1.0:
        Q['efs10_rcb01']=1-(1-P['efs10_rcb01'])*met_or_once
        Q['efs10_rcb23']=1-(1-P['efs10_rcb23'])*met_or_once
    return _markov_orj(Q,kol,rcb23)

def calistir(P, met_or=1.0):
    """S0/S1/S2 maliyet ve QALY — hasta başına"""
    n=P['kohort']; sonuc={}
    for ad,ni in [('S0',0),('S1',P['s1_n']),('S2',P['s2_n'])]:
        m_nak,q_nak=kisa_donem(P,'nak'); m_onc,q_onc=kisa_donem(P,'once')
        M_nak,Q_nak=_markov_orj(P,'nak',0.622); M_onc,Q_onc=markov_or(P,'once',0.622,met_or)
        M=( (n-ni)*(m_nak+M_nak) + ni*(m_onc+M_onc) )/n + (P['c_karar_destegi'] if ad=='S2' else 0)
        Q=( (n-ni)*(q_nak+Q_nak) + ni*(q_onc+Q_onc) )/n
        sonuc[ad]=(M,Q)
    return sonuc

TABAN=dict(BAZ)
b=calistir(TABAN)
print(f'Taban doğrulama: S0−S2 = {b["S0"][0]-b["S2"][0]:,.0f} TL · S1−S2 = {b["S1"][0]-b["S2"][0]:,.0f} TL\n')

SEN=[]
def ekle(ad, P=None, met_or=1.0, aciklama=''):
    Q=dict(TABAN); 
    if P: Q.update(P)
    r=calistir(Q,met_or)
    SEN.append(dict(senaryo=ad, aciklama=aciklama,
        S0_S2=r['S0'][0]-r['S2'][0], S1_S2=r['S1'][0]-r['S2'][0],
        dQALY_S1=r['S2'][1]-r['S1'][1]))

ekle('Taban senaryo', aciklama='tüm varsayımlar yürürlükte')
ekle('Kapesitabin çıkarıldı', {'c_nak_sistemik':TABAN['c_nak_sistemik']-0.10*19589},
     aciklama='CREATE-X yararı esas olarak üçlü negatifte; HR+/HER2−’de kapesitabin verilmediği varsayımı')
ekle('Kapesitabin iki katına', {'c_nak_sistemik':TABAN['c_nak_sistemik']+0.10*19589},
     aciklama='kullanım oranı %10 → %20')
ekle('ALND oranları eşitlendi', {'alnd_nak':TABAN['alnd_once']},
     aciklama='neoadjuvanın aksiller dezavantajı kaldırıldı')
ekle('LRR farkı sıfır', {'lrr_ek_neoadjuvan':0.0},
     aciklama='EBCTCG lokorejyonel nüks fazlası uygulanmadı')
ekle('LRR farkı yarıya', {'lrr_ek_neoadjuvan':0.0275},
     aciklama='modern cerrahi/radyoterapiye uyarlanmış')
ekle('Aynı sistemik rejim', {'c_adjuvan_sistemik':TABAN['c_nak_sistemik']-0.10*19589},
     aciklama='cerrahi-önce kolunda daha ucuz TC seçeneği yok')
ekle('Adjuvan KT oranı %100', {'adjuvan_kt_orani':1.0},
     aciklama='cerrahi-önce kolunda hiç kimse kemoterapiden kaçınmıyor')
ekle('Metastaz OR 1,10 (cerrahi-önce aleyhine)', met_or=1.10,
     aciklama='EBCTCG eşitlik varsayımı gevşetildi')
ekle('Metastaz OR 1,20 (cerrahi-önce aleyhine)', met_or=1.20,
     aciklama='EBCTCG eşitlik varsayımı güçlü biçimde gevşetildi')
ekle('Nüks epizodu 110.000 TL', {'c_lrr':110000}, aciklama='varsayım alt sınırı')
ekle('Nüks epizodu 250.000 TL', {'c_lrr':250000}, aciklama='varsayım üst sınırı')
ekle('İlaç: ortalama eşdeğer (+%40)', {'c_nak_sistemik':TABAN['c_nak_sistemik']*1.4,
     'c_adjuvan_sistemik':TABAN['c_adjuvan_sistemik']*1.4}, aciklama='en ucuz eşdeğer yerine ortalama fiyat')
ekle('Karar desteği 500 TL', {'c_karar_destegi':500}, aciklama='işletme maliyeti üst sınırı')
ekle('İç içe CV işaretleme (105)', {'s2_n':105}, aciklama='iyimserlik düzeltmesi sonrası')

D=pd.DataFrame(SEN)
print(f'{"Senaryo":42s} {"S0−S2 (TL)":>12s} {"S1−S2 (TL)":>12s} {"ΔQALY":>9s}')
for _,r in D.iterrows():
    print(f'{r.senaryo:42s} {r.S0_S2:12,.0f} {r.S1_S2:12,.0f} {r.dQALY_S1:9.5f}')
D.to_csv('../../outputs/yapisal_senaryolar.csv',index=False)
print(f'\nS1−S2 aralığı: {D.S1_S2.min():,.0f} — {D.S1_S2.max():,.0f} TL')
print(f'Tasarrufun kaybolduğu senaryo sayısı: {(D.S1_S2<0).sum()} / {len(D)}')
