"""ONLINE RESOURCE 4 — PSA parametre tablosunu markov_psa.py'den otomatik uretir."""
import re, math
s=open('../../outputs/markov_psa.py').read()
exec(s[:s.index('# ================================================================= 1. TABAN')])
P=BAZ
KUR, PPP = 48.3727, 15.351

# betikteki PSA bloklarindan dagilim atamalarini oku
psa=s[s.index('for i in range(N):'):s.index("df = {k:")]
gama=dict(re.findall(r"\('(\w+)',\.?(\d+)\)", psa[psa.index("gamma_par"):psa.index("for k,sd in")] if 'gamma_par' in psa else ''))
def blok(bas,son):
    try: return psa[psa.index(bas):psa.index(son)]
    except ValueError: return ''
b_gama = blok("for k,sd_orani in", "for k,sd in [('u_")
b_beta_u = blok("for k,sd in [('u_", "for k,sd in [('du_")
b_normal = blok("for k,sd in [('du_", "for k,sd in [('adjuvan")
b_beta_p = blok("for k,sd in [('adjuvan", "P['lrr_ek_neoadjuvan']")

G=dict((k,float('0.'+v)) for k,v in re.findall(r"\('([\w_]+)',\.(\d+)\)", b_gama))
BU=dict((k,float('0.'+v)) for k,v in re.findall(r"\('([\w_]+)',\.(\d+)\)", b_beta_u))
NO=dict((k,float('0.'+v)) for k,v in re.findall(r"\('([\w_]+)',\.(\d+)\)", b_normal))
BP=dict((k,float('0.'+v)) for k,v in re.findall(r"\('([\w_]+)',\.(\d+)\)", b_beta_p))

AD={ # degisken: (aciklama, birim, kaynak)
'c_nak_sistemik':('Neoadjuvan sistemik tedavi (kapesitabin dahil)','TL','SUT EK-2/B; TİTCK; EK-4/A'),
'c_adjuvan_sistemik':('Adjuvan sistemik tedavi (%70 AC→taksan, %30 TC)','TL','SUT EK-2/B; TİTCK; EK-4/A'),
'c_cerrahi_nak':('Cerrahi — neoadjuvan kolu (MKC %34,5)','TL','SUT EK-2/C'),
'c_cerrahi_once':('Cerrahi — cerrahi-önce kolu (MKC %27,0)','TL','SUT EK-2/C'),
'c_rt':('Radyoterapi (tam kür)','TL','SUT EK-2/B'),
'c_lrr':('Lokorejyonel nüks epizodu','TL','Varsayım ⬛'),
'c_met_yillik':('Metastatik hastalık — yıllık','TL','Varsayım ⬛'),
'c_lenfodem_yillik':('Lenfödem — yıllık tedavi','TL','Varsayım ⬛'),
'c_izlem':('Yıllık izlem (muayene, görüntüleme)','TL','SUT EK-2/A, EK-2/B'),
'c_karar_destegi':('Karar desteği — hasta başı işletme maliyeti','TL','Varsayım ⬛'),
'u_dfs':('Hastalıksız sağkalım','utility','Literatür'),
'u_kt':('Kemoterapi sırasında','utility','Literatür'),
'u_ilk_yil':('Tedavi sonrası ilk yıl','utility','Literatür'),
'u_lrr':('Lokorejyonel nüks','utility','Literatür'),
'u_met':('Metastatik hastalık','utility','Literatür'),
'du_lenfodem':('Lenfödem','disutility','Literatür'),
'du_mastektomi':('Mastektomi (MKC’ye göre)','disutility','Literatür'),
'du_fn':('Febril nötropeni epizodu','disutility','Literatür'),
'adjuvan_kt_orani':('Cerrahi-önce kolunda adjuvan kemoterapi alan oran','oran','Literatür (NCDB)'),
'p_lenfodem_alnd':('Lenfödem — aksiller diseksiyon sonrası','olasılık','Literatür'),
'p_lenfodem_slnb':('Lenfödem — sentinel biyopsi sonrası','olasılık','Literatür'),
'mkc_nak':('Meme koruyucu cerrahi — neoadjuvan kolu','oran','Z1071 (Boughey et al. 2014)'),
'mkc_once':('Meme koruyucu cerrahi — cerrahi-önce kolu','oran','EBCTCG oran oranı ile türetildi'),
'alnd_nak':('Aksiller diseksiyon — neoadjuvan kolu','oran','NCDB (Douglas et al. 2026)'),
'alnd_once':('Aksiller diseksiyon — cerrahi-önce kolu','oran','NCDB (Douglas et al. 2026)'),
'fn_insidans':('Febril nötropeni insidansı','olasılık','Literatür'),
'efs10_rcb01':('10 yıllık olaysız sağkalım — RCB-0/I','oran','Yau et al. 2022'),
'efs10_rcb23':('10 yıllık olaysız sağkalım — RCB-II/III','oran','Yau et al. 2022'),
'lrr_ek_neoadjuvan':('Neoadjuvan kola bağlı lokorejyonel nüks fazlası (15 yıl)','mutlak fark','EBCTCG 2018'),
's2_n':('Karar kuralının işaretlediği olgu sayısı','sayı','Bu çalışma — bootstrap'),
}

def sat(k, dag, par):
    ad,birim,kay = AD.get(k,(k,'',''))
    v=P[k]
    if birim=='TL':
        deger=f'{v:,.0f}'.replace(',','.')
        ek=f' ({v/PPP:,.0f})'.replace(',','.')
    elif birim in ('oran','olasılık','mutlak fark'):
        deger=f'{v:.3f}'.replace('.',',')
        ek=''
    elif birim=='sayı':
        deger=f'{v:.0f}'; ek=''
    else:
        deger=f'{v:.3f}'.replace('.',','); ek=''
    return f'| {ad} | {deger}{ek} | {birim} | {dag} | {par} | {kay} |'

BAS='| Parametre | Merkezi değer | Birim | Dağılım | Dağılım parametresi | Kaynak |\n|---|---|---|---|---|---|'
out=[]
out.append('# ONLINE RESOURCE 4 — OLASILIKSAL DUYARLILIK ANALİZİ PARAMETRELERİ\n')
out.append('**CHEERS 2022, madde 22.** Aşağıdaki tablo, 10.000 yinelemeli olasılıksal duyarlılık')
out.append('analizinde kullanılan tüm parametrelerin merkezi değerini, dağılım biçimini, dağılım')
out.append('parametresini ve kaynağını göstermektedir. Tablo analiz betiğinden otomatik üretilmiştir')
out.append('(`markov_psa.py`).\n')
out.append('> Maliyetler 2026 fiyatlarıyla Türk lirası cinsindendir; parantez içindeki değerler')
out.append('> OECD satın alma gücü paritesiyle (2025: 1 uluslararası dolar = 15,351 TL) dönüştürülmüştür.')
out.append('> ⬛ işaretli kalemler resmî fiyat listelerinden kalem kalem kurulamamış, varsayım olarak')
out.append('> girilmiş ve duyarlılık analizinde geniş aralıklarla taranmıştır.\n')

out.append('## Maliyet parametreleri — gama dağılımı\n')
out.append('Gama dağılımı, maliyetlerin negatif olamaması ve sağa çarpık olması nedeniyle seçilmiştir.\n')
out.append(BAS)
for k,sd in G.items(): out.append(sat(k,'Gama',f'SS = merkezi değerin %{sd*100:.0f}’i'))

out.append('\n## Fayda değerleri — beta dağılımı\n')
out.append('Beta dağılımı, utility değerlerinin 0–1 aralığında sınırlı olması nedeniyle seçilmiştir.\n')
out.append(BAS)
for k,sd in BU.items(): out.append(sat(k,'Beta',f'SS = {sd:.3f}'.replace('.',',')))

out.append('\n## Disutility değerleri — normal dağılım\n')
out.append('Sıfırda kesilmiş normal dağılım kullanılmıştır.\n')
out.append(BAS)
for k,sd in NO.items(): out.append(sat(k,'Normal (≥0)',f'SS = {sd:.3f}'.replace('.',',')))

out.append('\n## Oran ve olasılıklar — beta dağılımı\n')
out.append(BAS)
for k,sd in BP.items(): out.append(sat(k,'Beta',f'SS = {sd:.3f}'.replace('.',',')))

out.append('\n## Özel dağılımlar\n')
out.append(BAS)
out.append(sat('lrr_ek_neoadjuvan','Normal (≥0)','SS = 0,016 · %95 GA 0,024–0,086 (EBCTCG)'))
out.append(sat('s2_n','Normal (kırpılmış)','SS = 8,4 · sınırlar 60–188 · bootstrap dağılımından'))

out.append('\n---\n')
out.append('## Sabit tutulan parametreler\n')
out.append('Aşağıdaki parametreler duyarlılık analizinde sabit tutulmuştur:\n')
out.append('| Parametre | Değer | Gerekçe |\n|---|---|---|')
out.append(f'| Zaman ufku | {P["ufuk"]} yıl | Yaşam boyu; hasta yaşı ve fon mortalitesiyle uyumlu |')
out.append(f'| İndirim oranı | %{P["indirim"]*100:.0f} | Türkiye farmakoekonomik kılavuz önerisi |')
out.append(f'| Başlangıç yaşı | {P["yas"]} | Kohort medyanına yakın |')
out.append(f'| Kemoterapi süresi (yıl payı) | {P["kt_suresi"]:.2f} | Sekiz kür × 3 hafta |'.replace('.',','))
out.append(f'| Karar kuralının uygun popülasyonu | {P["uygun"]} olgu | Bu çalışma |')
out.append(f'| Alt tip kuralının işaretlediği olgu | {P["s1_n"]} olgu | Bu çalışma |')

out.append('\n---\n')
out.append('## Yeniden üretilebilirlik\n')
out.append('Tüm parametreler `markov_psa.py` betiğindeki `BAZ` sözlüğünde tanımlıdır. Olasılıksal')
out.append('duyarlılık analizi `numpy.random.default_rng(2026)` tohumuyla çalıştırılmıştır; betik aynı')
out.append('tohumla yeniden çalıştırıldığında sonuçlar birebir tekrarlanabilir.')

open('../../outputs/ESM_4_PSA_parametreleri.md','w').write('\n'.join(out))
print('yazildi ·', len(G)+len(BU)+len(NO)+len(BP)+2, 'parametre')
