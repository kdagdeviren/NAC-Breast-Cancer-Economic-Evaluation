import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'font.size':9,'font.family':'DejaVu Sans','figure.dpi':200})
C={'ana':'#3d5a80','iyi':'#2a9d8f','kotu':'#e76f51','notr':'#adb5bd','vurgu':'#f4a261','mor':'#7b2cbf'}
O='../../figures/'
def kutu(ax,x,y,w,h,renk,a='22',lw=1.2,ls='-'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.08",lw=lw,ec=renk,fc=renk+a,linestyle=ls))
def yaz(ax,x,y,t,fs=9,w='normal',c='#22333b',ha='center'):
    ax.text(x,y,t,ha=ha,va='center',fontsize=fs,fontweight=w,color=c)
def ok(ax,x1,y1,x2,y2,c='#6c757d',lw=1.1):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=11,lw=lw,color=c,shrinkA=0,shrinkB=0))

# ---------- S1: KASKAD MIMARISI
fig,ax=plt.subplots(figsize=(8.6,5.4)); ax.set_xlim(0,13); ax.set_ylim(0,8.4); ax.axis('off')
kutu(ax,0.3,3.4,2.5,1.6,C['ana'])
yaz(ax,1.55,4.55,'HASTA',10,'bold'); yaz(ax,1.55,4.15,'23 değişken',8.4,c='#4a5568')
yaz(ax,1.55,3.85,'patoloji · onkoloji',7.6,c='#6c757d'); yaz(ax,1.55,3.62,'· radyoloji',7.6,c='#6c757d')
ok(ax,2.8,4.2,3.5,6.0); ok(ax,2.8,4.2,3.5,2.4)
kutu(ax,3.5,5.3,4.3,1.9,C['iyi'])
yaz(ax,5.65,6.85,'KATMAN 1 — Yorumlanabilir',9.2,'bold',C['iyi'])
yaz(ax,5.65,6.45,'Katlama içi tek değişkenli tarama (p<0,10)',7.8,c='#4a5568')
yaz(ax,5.65,6.13,'L2 düzenlenmiş lojistik regresyon (C=0,05)',7.8,c='#4a5568')
yaz(ax,5.65,5.72,'AUC 0,798 · Brier 0,179',8.2,'bold')
kutu(ax,3.5,1.5,4.3,1.9,C['vurgu'])
yaz(ax,5.65,3.05,'KATMAN 2 — Örüntü',9.2,'bold','#c77d20')
yaz(ax,5.65,2.65,'Histogram tabanlı gradyan artırma',7.8,c='#4a5568')
yaz(ax,5.65,2.33,'300 yineleme · 15 yaprak · lr 0,05',7.8,c='#4a5568')
yaz(ax,5.65,1.92,'AUC 0,808 · Brier 0,181',8.2,'bold')
ok(ax,7.8,6.25,8.7,4.9); ok(ax,7.8,2.45,8.7,3.8)
kutu(ax,8.7,3.1,2.6,1.6,C['mor'],'20',1.5)
yaz(ax,10.0,4.25,'ORTALAMA',9.4,'bold',C['mor'])
yaz(ax,10.0,3.85,'(p₁ + p₂) / 2',8.6,c='#4a5568')
yaz(ax,10.0,3.45,'AUC 0,820 · Brier 0,166',8.2,'bold')
ok(ax,11.3,3.9,12.0,3.9)
yaz(ax,12.5,4.15,'p',12,'bold',C['mor']); yaz(ax,12.5,3.6,'[0–1]',7.6,c='#6c757d')
yaz(ax,5.65,0.75,'Katmanlar arası korelasyon 0,814 — aynı sinyali yakalayacak kadar benzer,',8.2,c='#4a5568')
yaz(ax,5.65,0.4,'hataları birbirini götürecek kadar farklı',8.2,c='#4a5568')
plt.tight_layout(); plt.savefig(O+'Sunum_A_kaskad_mimarisi.png',bbox_inches='tight'); plt.close()

# ---------- S2: DOGRULAMA SEMASI
fig,ax=plt.subplots(figsize=(8.6,4.6)); ax.set_xlim(0,13); ax.set_ylim(0,7.2); ax.axis('off')
yaz(ax,6.5,6.85,'5 tekrar × 5 katlama · hasta düzeyinde gruplandırılmış',9.6,'bold')
renk=[C['notr']]*5
for t in range(5):
    y=5.5-t*1.05
    yaz(ax,0.55,y+0.28,f'Tur {t+1}',8.6,'bold',ha='center')
    for g in range(5):
        x=1.6+g*2.1
        sinav = (g==t)
        kutu(ax,x,y,1.95,0.62,C['kotu'] if sinav else C['iyi'],'30' if sinav else '18',1.0)
        yaz(ax,x+0.98,y+0.31,'SINAV' if sinav else 'öğrenme',7.8,'bold' if sinav else 'normal',
            C['kotu'] if sinav else '#4a5568')
yaz(ax,6.5,0.55,'Her turda model dört gruptan öğrenir, beşincide sınanır.',8.4,c='#4a5568')
yaz(ax,6.5,0.18,'Değişken seçimi de her turun İÇİNDE, yalnızca öğrenme verisinde yapılır.',8.4,'bold',c='#22333b')
plt.tight_layout(); plt.savefig(O+'Sunum_B_dogrulama_semasi.png',bbox_inches='tight'); plt.close()

# ---------- S3: SIZINTI ONCESI-SONRASI
fig,ax=plt.subplots(figsize=(6.4,3.6))
gr=['Radyoloji bloğu\ntek başına','Tam model']
kont=[0.764,0.904]; temiz=[0.571,0.806]
import numpy as np
x=np.arange(2); w=0.34
ax.bar(x-w/2,kont,w,label='Kontamine kodlama',color=C['kotu'],alpha=.85)
ax.bar(x+w/2,temiz,w,label='Temiz kodlama',color=C['iyi'],alpha=.85)
for i,(a,b) in enumerate(zip(kont,temiz)):
    ax.text(i-w/2,a+.012,f'{a:.3f}',ha='center',fontsize=8.4,fontweight='bold')
    ax.text(i+w/2,b+.012,f'{b:.3f}',ha='center',fontsize=8.4,fontweight='bold')
    ax.annotate('',xy=(i+w/2,b),xytext=(i-w/2,a),arrowprops=dict(arrowstyle='->',color='#495057',lw=1.2))
    ax.text(i,(a+b)/2+.03,f'−{a-b:.3f}',ha='center',fontsize=8.6,fontweight='bold',color=C['kotu'])
ax.set_xticks(x); ax.set_xticklabels(gr,fontsize=9); ax.set_ylabel('Ayırt gücü (AUC)')
ax.set_ylim(0.5,1.0); ax.legend(fontsize=8,frameon=False,loc='upper left')
ax.spines[['top','right']].set_visible(False)
plt.tight_layout(); plt.savefig(O+'Sunum_C_sizinti.png',bbox_inches='tight'); plt.close()
print('3 sunum sekli uretildi')
