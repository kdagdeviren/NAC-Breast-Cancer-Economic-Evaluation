# -*- coding: utf-8 -*-
"""Fig 10 — Structure of the economic model. Journal spec: 600 dpi, EPS, ≤174 mm."""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon
plt.rcParams.update({'ps.fonttype':42,'pdf.fonttype':42,'savefig.dpi':600,'figure.dpi':600,
    'font.size':8,'font.family':'DejaVu Sans'})
C={'ana':'#3d5a80','iyi':'#2a9d8f','kotu':'#e76f51','notr':'#8d99a6','vurgu':'#e08a3c','mor':'#7b2cbf'}
def kutu(ax,x,y,w,h,r,a='16',lw=1.1,ls='-'):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.05",lw=lw,ec=r,fc=r+a,linestyle=ls))
def yaz(ax,x,y,t,fs=8,w='normal',c='#22333b',ha='center'):
    ax.text(x,y,t,ha=ha,va='center',fontsize=fs,fontweight=w,color=c)
def ok(ax,x1,y1,x2,y2,c='#6c757d',lw=1.0):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=8,lw=lw,color=c,
                 shrinkA=0,shrinkB=0))

fig,ax=plt.subplots(figsize=(6.85,5.00))
ax.set_xlim(0,20); ax.set_ylim(0,14.6); ax.axis('off')

# ════════ PANEL A
yaz(ax,0.2,14.2,'a  Short-term decision tree (12 months)',9,'bold',C['ana'],ha='left')
ax.add_patch(Polygon([[0.9,11.9],[1.5,12.3],[2.1,11.9],[1.5,11.5]],closed=True,fc='white',ec=C['ana'],lw=1.2))
yaz(ax,1.5,10.9,'Eligible\npatient',7.5)
for yy,ad,renk,alt in [(13.0,'Neoadjuvant chemotherapy first',C['vurgu'],'8 cycles → surgery → RT / endocrine'),
                       (10.5,'Surgery first',C['iyi'],'surgery → pathology → chemotherapy if indicated')]:
    ok(ax,2.2,11.9,3.3,yy+0.45)
    kutu(ax,3.4,yy,7.2,0.95,renk,'20',1.2)
    yaz(ax,7.0,yy+0.60,ad,8,'bold')
    yaz(ax,7.0,yy+0.24,alt,6.6,c='#4a5568')
    ok(ax,10.7,yy+0.47,11.52,yy+0.47)
kutu(ax,11.60,10.3,4.05,3.4,C['ana'],'10',1.0)
yaz(ax,13.62,13.40,'12-month costs',7.5,'bold',C['ana'])
for k,t in enumerate(['Systemic therapy','Surgery and inpatient stay','Radiotherapy',
                      'Supportive care','Febrile neutropenia','Follow-up visits']):
    yaz(ax,11.84,12.95-k*0.42,'· '+t,6.2,ha='left')
ok(ax,15.75,12.0,16.55,12.0)
kutu(ax,16.65,11.05,3.25,1.9,C['mor'],'14',1.2)
yaz(ax,18.28,12.58,'Enter Markov',7.2,'bold',C['mor'])
yaz(ax,18.28,12.24,'model in the',7.2,'bold',C['mor'])
yaz(ax,18.28,11.80,'disease-free',6.9)
yaz(ax,18.28,11.46,'state',6.9)

ax.plot([0.2,19.8],[9.8,9.8],color='#ccd6dd',lw=0.8)

# ════════ PANEL B
yaz(ax,0.2,9.15,'b  Lifetime Markov model (40 years, 1-year cycles, 3% discounting)',9,'bold',C['mor'],ha='left')
# durum kutulari — genis aralikli
DF=(0.6,6.10,4.0,1.75); LR=(8.0,6.10,4.2,1.75); DM=(8.0,2.30,4.2,1.75); DE=(15.4,2.30,3.4,1.75)
for (x,y,w,h),ad,u,renk in [(DF,'DISEASE-FREE','u = 0.88',C['iyi']),
                            (LR,'LOCOREGIONAL\nRECURRENCE','u = 0.78',C['vurgu']),
                            (DM,'DISTANT\nMETASTASIS','u = 0.73',C['kotu']),
                            (DE,'DEATH','u = 0',C['notr'])]:
    kutu(ax,x,y,w,h,renk,'20',1.3)
    yaz(ax,x+w/2,y+h-0.55 if '\n' in ad else y+h-0.62,ad,7.8,'bold')
    yaz(ax,x+w/2,y+0.42,u,7.4,c='#4a5568')

# oklar — etiketler bosluklarda
ok(ax,4.70,7.42,7.90,7.42,C['vurgu'],1.1)
yaz(ax,6.30,7.72,'recurrence',6.6,c='#6c757d')
ok(ax,7.90,6.55,4.70,6.55,'#adb5bd',0.9)
yaz(ax,6.30,6.28,'returns after 1 year',6.2,c='#9aa5b1')
ok(ax,10.10,6.00,10.10,4.15,C['kotu'],1.1)
yaz(ax,10.45,5.10,'progression',6.6,c='#6c757d',ha='left')
ok(ax,2.60,6.00,7.90,4.30,C['kotu'],1.1)
yaz(ax,4.40,4.72,'progression',6.6,c='#6c757d')
ok(ax,12.30,3.18,15.30,3.18,'#6c757d',1.0)
yaz(ax,13.80,3.62,'all-cause and',6.3,c='#6c757d')
yaz(ax,13.80,3.30,'cancer mortality',6.3,c='#6c757d')

# varsayim kutusu — genis, tasmasiz
kutu(ax,13.10,5.55,6.70,3.35,C['kotu'],'07',1.1,'--')
yaz(ax,16.45,8.55,'Load-bearing assumption',7.6,'bold',C['kotu'])
for k,t in enumerate(['The rule forecasts residual disease burden;',
                      'it does not change it. Distant metastasis',
                      'risk is therefore set equal in both arms',
                      '(EBCTCG 2018). Relaxing this assumption',
                      'reverses the economic result (Table 15).']):
    yaz(ax,13.42,8.10-k*0.42,t,6.4,ha='left',c='#4a5568')

yaz(ax,0.6,1.35,'Perspective: payer  ·  Outcome: quality-adjusted life-years  ·  Costs in 2026 Turkish lira',
    7,c='#4a5568',ha='left')
yaz(ax,0.6,0.85,'u = health state utility. Arm-specific inputs: surgery mix, adjuvant regimen,',6.3,c='#8a94a0',ha='left')
yaz(ax,0.6,0.48,'axillary dissection rate, and excess locoregional recurrence.',6.3,c='#8a94a0',ha='left')
plt.tight_layout()
plt.savefig('../../figures/Fig3.png',bbox_inches='tight',facecolor='white'); plt.savefig('../../figures/Fig3.pdf',bbox_inches='tight',facecolor='white')
plt.savefig('../../figures/Fig3.eps',bbox_inches='tight',facecolor='white',format='eps')
plt.close(); print('ok')
