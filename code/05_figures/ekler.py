# -*- coding: utf-8 -*-
"""ESM_2 (veri kalitesi) ve ESM_3 (kod ve veri) — mevcut belgelerden turetilir."""
import re
D='../../outputs/'
gun=open(D+'Duzeltme_Gunlugu_v17.md').read()

# gunluk tablosunu al ve anonimlestir/duzelt
satir=[l for l in gun.split('\n') if l.startswith('| ') and '---' not in l and not l.startswith('| Madde')]
def tr(t):
    R={'duzeltildi':'düzeltildi','yanlis':'yanlış','kodlanmis':'kodlanmış','celiskisi':'çelişkisi',
    'degere':'değere','gore':'göre','kalan':'kalan','uyumsuzluk':'uyumsuzluk','gercek':'gerçek',
    'bostu':'boştu','ayri':'ayrı','kategorisine':'kategorisine','alindi':'alındı','kapsami':'kapsamı',
    'disinda':'dışında','dogrulanmadan':'doğrulanmadan','atanmadi':'atanmadı','dusuk':'düşük',
    'kayittaki':'kayıttaki','skoruna':'skoruna','etiket':'etiket','zayif':'zayıf','belirlendi':'belirlendi',
    'imputasyonu':'imputasyonu','kuralini':'kuralını','ihlal':'ihlal','etmisti':'etmişti','esas':'esas',
    'alindi':'alındı','eksik':'eksik','kayitlarda':'kayıtlarda','yerine':'yerine','kullanildi':'kullanıldı',
    'kural':'kural','uyumlu':'uyumlu','degisen':'değişen','kolonu':'kolonu','eklendi':'eklendi',
    'kategorisi':'kategorisi','kaldirildi':'kaldırıldı','luminal':'luminal','gruplara':'gruplara',
    'dagitildi':'dağıtıldı','durumu':'durumu','degisken':'değişken','oldu':'oldu','siniflamasinda':'sınıflamasında',
    'ciftinde':'çiftinde','iki':'iki','kayit':'kayıt','ilki':'ilki','katmanli':'katmanlı','kimlik':'kimlik',
    'olusturuldu':'oluşturuldu','hastada':'hastada','yas':'yaş','tani':'tanı','yili':'yılı',
    'arasindaki':'arasındaki','yillik':'yıllık','fark':'fark','birakildi':'bırakıldı','kullaniliyor':'kullanılıyor',
    'degisiklik':'değişiklik','metastatik':'metastatik','olgular':'olgular','taban':'taban','senaryodan':'senaryodan',
    'cikarildi':'çıkarıldı','duyarlilik':'duyarlılık','analizinde':'analizinde','raporlanacak':'raporlanacak',
    'turu':'türü','radyoloji':'radyoloji','kayitlarindan':'kayıtlarından','teyit':'teyit','edilerek':'edilerek',
    'uc':'üç','lezyon':'lezyon','atamisti':'atamıştı','cozulmemis':'çözülmemiş','olgusunda':'olgusunda',
    'sinifi':'sınıfı','yapildi':'yapıldı','kolonlari':'kolonları','kullanilmiyor':'kullanılmıyor',
    'imputasyonlu':'imputasyonlu','surumleri':'sürümleri','bozuk':'bozuk','metin':'metin','bicim':'biçim',
    'varyantlari':'varyantları','tekillestirildi':'tekilleştirildi','sayisala':'sayısala','cevrildi':'çevrildi',
    'carpma':'çarpma','isareti':'işareti','standartlastirildi':'standartlaştırıldı','rejim':'rejim',
    'bilgisi':'bilgisi','dosyalarindan':'dosyalarından','dogru':'doğru','bilmis':'bilmiş','antrasiklin':'antrasiklin',
    'sadece':'sadece','ham':'ham','dolduruldu':'dolduruldu','zaten':'zaten','dogruydu':'doğruydu',
    'ekvokal':'ekvokal','olgularda':'olgularda','sonuclari':'sonuçları','temin':'temin','edildi':'edildi',
    'dokuz':'dokuz','olgu':'olgu','negatif':'negatif','olarak':'olarak','dogrulandi':'doğrulandı',
    'kapsamina':'kapsamına','alinmis':'alınmış','belirsiz':'belirsiz','kaldi':'kaldı',
    'kodlanmis':'kodlanmış','kolonlar':'kolonlar','kolonlariyla':'kolonlarıyla','yeniden':'yeniden',
    'esitlendi':'eşitlendi','duzeltmeleri':'düzeltmeleri','kolonlarina':'kolonlarına','yansimamisti':'yansımamıştı'}
    for a,b in R.items(): t=re.sub(r'\b'+a+r'\b',b,t)
    return t

out=['# ONLINE RESOURCE 2 — VERİ KALİTESİ DENETİMİ VE DÜZELTME GÜNLÜĞÜ\n',
'Model kurulmadan önce ham veri seti sistematik bir kalite denetiminden geçirilmiştir.',
'Denetim sırasında saptanan her tutarsızlık, düzeltme kararı ve etkilenen kayıt sayısı',
'aşağıda belgelenmiştir. Düzeltmeler sürüm numaralı bir betikle (`kohort_duzeltme.py`)',
'uygulanmış, her sürüm yeniden üretilebilir biçimde kaydedilmiştir.\n',
'**Kaynak:** ham veri seti · **Çıktı:** analiz kohortu (328 kayıt / 312 hasta)\n',
'---\n','## Düzeltme kayıtları\n',
'| Madde | Düzeltme | Etkilenen kayıt | Karar ve gerekçe |','|---|---|---|---|']
for l in satir: out.append(tr(l))

out+= ['\n---\n','## Kimlik yapısı\n',
'İki katmanlı kimlik kullanılmıştır; bilateral olgularda aynı hastanın iki kaydı bulunduğu için',
'tüm çapraz doğrulama işlemleri hasta düzeyinde gruplandırılmıştır.\n',
'| Alan | İşlev |','|---|---|',
'| Kayıt anahtarı | Satır anahtarı; her meme bir kayıt |',
'| Hasta anahtarı | Kişi anahtarı; bilateral olgularda iki kayıtta aynı |',
'| Özgün numara | İzlenebilirlik için korundu |',
'| Taraf | Sağ/sol bilgisinin tek kaynağı |\n',
'Kimlik belirleyici tüm alanlar (ad, hasta numarası, doğum yeri, coğrafi koordinatlar)',
'veri setinden çıkarılmıştır.\n',
'---\n','## En kritik bulgu: sonuç sızıntısı\n',
'Radyolojik değişkenlerin iki farklı kodlama sürümü karşılaştırıldığında, bir sürümde bazı',
'hücrelere radyolojik bulgu yerine **hastanın sonucu** yazıldığı saptanmıştır (390 hücre; ham',
'değerlerin %11,0’i). Bulgu şu testle doğrulanmıştır: etkilenen satırlarda radyoloji bloğunun',
'tek başına ayırt gücü 0,92, etkilenmemiş satırlarda 0,54’tür. Kontamine kodlama kullanıldığında',
'tam modelin ayırt gücü 0,904, temiz kodlamayla 0,806 çıkmaktadır; aradaki fark biyolojik',
'sinyalden değil sızıntıdan kaynaklanmaktadır. Tüm analizler temiz kodlama üzerinde',
'yürütülmüştür.\n',
'---\n','## Hedef değişkenin doğrulanması\n',
'Rezidüel kanser yükü indeksi, kayıtlı sınıflandırmadan bağımsız olarak özgün formülle yeniden',
'hesaplanmış ve 328 olgunun 328’inde kayıtlı sınıfla uyumlu bulunmuştur.\n',
'---\n','## Çözülmemiş konular\n',
'Şeffaflık gereği, denetimde saptanan ancak çözülemeyen konular aşağıda listelenmiştir:\n',
'| Konu | Durum |','|---|---|',
'| Tanı anındaki evre ile nodal tutulum arasında kayıtlarda tutarsızlık | Kolon kayıtlı haliyle kullanıldı; metinde “kayıtlarda bildirilen evre” olarak tanımlandı |',
'| Bir olguda lezyon türü çelişkisi | Çözülemedi; imputasyonlu değer korundu |',
'| Bir ekvokal olguda ISH sonucuna ulaşılamadı | Uygunluk kapsamı dışında bırakıldı |',
'| Kreatinin kolonunda 13 kayıtta birim karışıklığı | İlgili blok bu çalışmada kullanılmadığı için etkisiz |']
open(D+'ESM_2_Veri_kalitesi.md','w').write('\n'.join(out))

# ── ESM_3
ESM3 = """# ONLINE RESOURCE 3 — ANALİZ KODU VE VERİ SETİ

## Veri erişilebilirliği

Bu çalışmanın bulgularını destekleyen anonimleştirilmiş veri seti, makul talep üzerine sorumlu
yazardan temin edilebilir. Bireysel hasta kayıtları kimlik belirleyici olabilecek klinik bilgi
içerdiğinden kamuya açık bir depoda paylaşılamamaktadır. Kimlik belirleyici tüm alanlar analiz
öncesinde veri setinden çıkarılmıştır.

## Kod erişilebilirliği

Analiz kodu, yayım sonrasında kalıcı tanımlayıcılı bir depoda erişime açılacaktır. Depo
aşağıdaki bileşenleri içerecektir:

| Bileşen | İşlev |
|---|---|
| Kohort düzeltme betiği | Ham veriden analiz kohortunu üretir; her düzeltme maddesini sürüm numaralı olarak kaydeder |
| Taban senaryo betiği | Karar kuralını hasta düzeyinde gruplandırılmış çapraz doğrulamayla çalıştırır |
| Katman analizi betiği | Yorumlanabilir katman, örüntü katmanı, istifleme ve kaskad yapılarını karşılaştırır |
| Önyükleme ve karar eğrisi betiği | Güven aralıklarını ve net fayda eğrilerini hesaplar |
| Eşik eğrisi betiği | Karar eşiğinin 0,60–0,90 aralığında taranmasını sağlar |
| İlaç maliyeti betiği | Fiyat listesi ve iskonto oranlarından ödeyici maliyetini türetir |
| Maliyet modeli betiği | Kalem kalem tedavi maliyetlerini hesaplar |
| Markov ve olasılıksal duyarlılık betiği | Yaşam boyu modeli, 10.000 yinelemeli duyarlılık analizini, tornado ve kabul edilebilirlik eğrisini üretir |
| Şekil betikleri | Makaledeki tüm şekilleri yeniden üretir |

## Yeniden üretilebilirlik

Tüm rastgelelik içeren işlemler sabit tohumlarla çalıştırılmıştır:

| İşlem | Tohum |
|---|---|
| Çapraz doğrulama (5 tekrar) | 300–304 |
| Küme önyüklemesi | 2026 |
| Olasılıksal duyarlılık analizi | 2026 |

Betikler aynı tohumlarla yeniden çalıştırıldığında makalede bildirilen tüm sayısal sonuçlar
birebir tekrarlanmaktadır.

## Yazılım ortamı

Analizler Python 3 ile yürütülmüştür. Kullanılan kütüphaneler: `scikit-learn` (model kurulumu ve
çapraz doğrulama), `pandas` ve `numpy` (veri işleme), `scipy` (istatistiksel testler),
`statsmodels` (tanısal analizler), `matplotlib` (şekiller), `openpyxl` (veri okuma).

## Yapay zekâ araçlarının kullanımı

Çalışmanın analitik içeriği, model kurulumu, parametre seçimi ve sonuçların yorumlanması
yazarlara aittir. Büyük dil modeli tabanlı araçlar yalnızca metin düzenleme ve dil kontrolü
amacıyla kullanılmış; sonuç üretimi, veri analizi veya yorumlama aşamalarında kullanılmamıştır.
"""
open(D+'ESM_3_Kod_ve_veri.md','w').write(ESM3)
print('iki ek uretildi')
