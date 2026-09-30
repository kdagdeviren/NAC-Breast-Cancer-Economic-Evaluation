import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'ps.fonttype':42,'pdf.fonttype':42,'savefig.dpi':600,'figure.dpi':600,'font.size':9,'font.family':'DejaVu Sans'})
C={'ana':'#3d5a80','iyi':'#2a9d8f','kotu':'#e76f51','notr':'#adb5bd','vurgu':'#f4a261'}
fig,ax=plt.subplots(figsize=(6.75,5.54)); ax.set_xlim(0,12); ax.set_ylim(0,10.6); ax.axis('off')

def kutu(x,y,w,h,renk,alpha='22',lw=1.1,ls='-'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.09",lw=lw,ec=renk,
                 fc=renk+alpha,linestyle=ls))
def yaz(x,y,t,fs=9,w='normal',c='#22333b',ha='center'):
    ax.text(x,y,t,ha=ha,va='center',fontsize=fs,fontweight=w,color=c)
def ok(x1,y1,x2,y2,c='#6c757d'):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=10,lw=1.0,
                 color=c,shrinkA=0,shrinkB=0))

# ust: kohort
kutu(4.1,9.5,3.8,0.85,C['ana'])
yaz(6,10.13,'320 records / 304 patients',9.6,'bold')
yaz(6,9.77,'Neoadjuvant candidates, cM0 cohort',8.2,c='#4a5568')
for x in [2.0,6.0,10.0]:
    ok(6,9.5,6,9.15); ok(6,9.15,x,9.15) if x!=6 else None
    ok(x,9.15,x,8.75)
ax.plot([2.0,10.0],[9.15,9.15],color='#6c757d',lw=1.0,zorder=0)

STRAT=[
 (2.0,'S0','Standard care',C['notr'],
  ['All cases','receive neoadjuvant','chemotherapy as per practice'],
  [('Cases redirected','0'),('Additional cost','none')],
  None),
 (6.0,'S1','Molecular subtype rule',C['iyi'],
  ['Surgery-first considered','in Luminal A','cases'],
  [('Flagged','61'),('Concordant / discordant','55 / 6'),('Predictive value','0.902')],
  'No additional cost — panel already performed'),
 (10.0,'S2','Prediction rule',C['vurgu'],
  ['In HR+/HER2− biology,','cases with probability','≥ 0.80 are flagged'],
  [('Flagged','108'),('Concordant / discordant','98 / 10'),('Predictive value','0.907')],
  'Additional cost: data entry and system'),
]
for x,kod,ad,renk,kural,sonuc,not_ in STRAT:
    kutu(x-1.875,8.0,3.75,0.75,renk,'33',1.4)
    yaz(x-1.66,8.55,kod,10.2,'bold',renk,ha='left'); yaz(x+0.28,8.53,ad,7.6,'bold')
    kutu(x-1.875,6.35,3.75,1.45,renk,'11',0.9,'--')
    for i,l in enumerate(kural): yaz(x,7.55-i*0.36,l,7.9,c='#4a5568')
    kutu(x-1.875,4.4,3.75,1.75,renk,'18')
    for i,(k,v) in enumerate(sonuc):
        yaz(x-1.70,5.85-i*0.45,k,7.6,c='#4a5568',ha='left')
        yaz(x+1.70,5.85-i*0.45,v,8.2,'bold',ha='right')
    if not_: yaz(x,4.08,not_,7.0,c='#6c757d')

# karsilastirma kutusu
kutu(1.0,0.5,10.0,3.0,C['ana'],'0e',1.0)
yaz(6,3.15,'S2 VERSUS S1',9.6,'bold')
ok(4.9,4.4,4.9,3.5,C['iyi']); ok(11.1,4.4,11.1,3.7,C['vurgu'])
ax.plot([11.1,11.1],[3.7,3.5],color=C['vurgu'],lw=1.0); ok(11.1,3.62,11.05,3.5,C['vurgu'])
sat=[('Additional concordant','+43','+28 – +58',C['iyi']),
     ('Additional discordant','+4','−1 – +9',C['kotu'])]
yaz(3.2,2.62,'',8); yaz(6.6,2.62,'Observed',8.4,'bold',c='#4a5568')
yaz(9.2,2.62,'95% confidence interval',8.4,'bold',c='#4a5568')
for i,(k,v,ga,c) in enumerate(sat):
    y=2.15-i*0.55
    yaz(1.4,y,k,8.8,'bold',c,ha='left'); yaz(6.6,y,v,9.4,'bold',c); yaz(9.2,y,ga,8.8,c='#4a5568')
yaz(6,0.85,'The gain excludes zero; the additional discordant classifications do not',
    8.3,c='#4a5568')
plt.tight_layout(); plt.savefig('../../figures/Fig2.png',bbox_inches='tight'); plt.savefig('../../figures/Fig2.pdf',bbox_inches='tight',facecolor='white'); plt.savefig('../../figures/Fig2.eps',bbox_inches='tight',format='eps'); plt.savefig('../../figures/Fig2.pdf',bbox_inches='tight')
print('ok')
