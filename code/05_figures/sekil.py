import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
plt.rcParams.update({'ps.fonttype':42,'pdf.fonttype':42,'savefig.dpi':600,'figure.dpi':600,'font.size':9,'font.family':'DejaVu Sans','axes.spines.top':False,
                     'axes.spines.right':False,'figure.dpi':200})
C={'ana':'#3d5a80','iyi':'#2a9d8f','kotu':'#e76f51','notr':'#adb5bd','vurgu':'#f4a261'}
O='../../figures/'
import os; os.makedirs(O,exist_ok=True)

# ---------- Sekil 1: hasta akis semasi
fig,ax=plt.subplots(figsize=(6.75,7.12)); ax.set_xlim(0,10); ax.set_ylim(0,11.6); ax.axis('off')
def kutu(x,y,w,h,bas,alt,renk,fs=9):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.08",lw=1.1,
                 ec=renk,fc=renk+'22'))
    ax.text(x+w/2,y+h*0.63,bas,ha='center',va='center',fontsize=fs+0.5,fontweight='bold',color='#22333b')
    ax.text(x+w/2,y+h*0.27,alt,ha='center',va='center',fontsize=fs-1.2,color='#4a5568',linespacing=1.25)
def ok(x1,y1,x2,y2):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=11,
                 lw=1.0,color='#6c757d',shrinkA=0,shrinkB=0))
kutu(2.55,10.2,4.9,1.0,'328 records / 312 patients','Breast cancer cohort receiving\nneoadjuvant chemotherapy',C['ana'])
ok(5,10.2,5,9.55); ok(5,9.85,2.45,9.85)
kutu(0.05,9.45,2.35,0.8,'8 records','Distant metastasis (cM1)\nsensitivity analysis',C['notr'],8.5)
kutu(2.55,8.55,4.9,1.0,'320 records / 304 patients','Base-case cohort (cM0)',C['ana'])
ok(5,8.55,5,7.9); ok(5,8.2,2.45,8.2)
kutu(0.05,7.8,2.35,0.8,'132 records','Ineligible biology\nrule does not operate',C['notr'],8.5)
kutu(2.55,6.9,4.9,1.0,'188 records','HR+/HER2− — rule operates',C['ana'])
ok(5,6.9,5,6.25); ok(5,6.55,2.45,6.55)
kutu(0.05,6.15,2.35,0.8,'80 records','Probability < 0.80\nstandard care',C['notr'],8.5)
kutu(2.55,5.25,4.9,1.0,'108 records flagged','High residual burden predicted (p ≥ 0.80)',C['vurgu'])
ax.add_patch(FancyArrowPatch((5,5.25),(3.3,4.6),arrowstyle='-|>',mutation_scale=11,lw=1.0,color='#6c757d'))
ax.add_patch(FancyArrowPatch((5,5.25),(6.7,4.6),arrowstyle='-|>',mutation_scale=11,lw=1.0,color='#6c757d'))
kutu(1.3,3.6,3.0,0.95,'98 records','Concordant (RCB-II/III)',C['iyi'])
kutu(5.7,3.6,3.0,0.95,'10 records','Discordant (RCB-0/I)',C['kotu'])
ax.text(5,2.95,'Positive predictive value = 0.907 (95% CI 0.849–0.958)',ha='center',fontsize=9,
        fontweight='bold',color='#22333b')
ax.add_patch(FancyBboxPatch((0.6,0.55),8.8,2.05,boxstyle="round,pad=0.12",lw=1.0,
             ec='#adb5bd',fc='#f8f9fa'))
ax.text(5,2.33,'Comparator — molecular subtype rule (S1)',ha='center',fontsize=9,
        fontweight='bold',color='#22333b')
ax.text(5,1.9,'61 flagged · 55 concordant · 6 discordant · PPV 0.902',ha='center',fontsize=8.5,color='#4a5568')
ax.text(5,1.35,'Contribution of the prediction rule over the subtype rule',ha='center',fontsize=8.8,
        color='#22333b',fontweight='bold')
ax.text(5,0.9,'+43 additional concordant (95% CI +28 to +58)   ·   +4 additional discordant (95% CI −1 to +9)',
        ha='center',fontsize=8.3,color='#4a5568')
plt.tight_layout(); plt.savefig(O+'Fig1.png',bbox_inches='tight'); plt.savefig(O+'Fig1.eps',bbox_inches='tight',format='eps'); plt.savefig(O+'Fig1.pdf',bbox_inches='tight'); plt.close()

# ---------- Sekil 2: alt tip x pCR
alt=['Luminal A','Luminal B\n(HER2−)','HER2 ekvokal','Luminal B\n(HER2+)','Triple\nnegatif','HER2-Zengin']
n=[53,129,12,55,41,38]; pcr=[1.9,18.6,33.3,41.8,43.9,52.6]
fig,ax=plt.subplots(figsize=(6.75,3.19))
renk=[C['iyi'] if a.startswith(('Luminal A','Luminal B\n(HER2−)')) else C['notr'] for a in alt]
b=ax.bar(range(6),pcr,color=renk,width=0.62,edgecolor='white')
for i,(v,nn) in enumerate(zip(pcr,n)):
    ax.text(i,v+1.4,f'%{v}',ha='center',fontsize=8.5,fontweight='bold',color='#22333b')
    ax.text(i,-4.2,f'n={nn}',ha='center',fontsize=7.8,color='#6c757d')
ax.set_xticks(range(6)); ax.set_xticklabels(alt,fontsize=8)
ax.set_ylabel('Pathological complete response rate'); ax.set_ylim(0,60)
ax.axhspan(0,0,color='none')
ax.text(0.5,55,'biology in which the rule operates',ha='center',fontsize=8,
        color=C['iyi'],style='italic')
ax.plot([-0.42,1.42],[52,52],color=C['iyi'],lw=1.4)
plt.tight_layout(); plt.savefig(O+'Fig4.png',bbox_inches='tight'); plt.savefig(O+'Fig4.eps',bbox_inches='tight',format='eps'); plt.savefig(O+'Fig4.pdf',bbox_inches='tight'); plt.close()

# ---------- Sekil 3: esik egrisi
e=pd.read_csv('../../outputs/esik_egrisi.csv')
fig,axs=plt.subplots(1,2,figsize=(6.75,2.84))
a=axs[0]
a.plot(e.esik,e.PPV,'o-',color=C['ana'],lw=1.6,ms=4,label='Prediction rule')
a.fill_between(e.esik,e.PPV_alt,e.PPV_ust,color=C['ana'],alpha=0.13)
a.axhline(0.898,ls='--',lw=1.2,color=C['kotu'],label='Subtype rule (S1)')
a.axvline(0.80,ls=':',lw=1.0,color='#6c757d')
a.set_xlabel('Decision threshold'); a.set_ylabel('Positive predictive value'); a.set_ylim(0.70,1.0)
a.legend(fontsize=7.5,frameon=False,loc='lower right')
a=axs[1]
a.plot(e.esik,e.net_dogru,'o-',color=C['iyi'],lw=1.6,ms=4,label='additional concordant')
a.fill_between(e.esik,e.net_dogru_alt,e.net_dogru_ust,color=C['iyi'],alpha=0.15)
a.plot(e.esik,e.net_yanlis,'s-',color=C['kotu'],lw=1.6,ms=4,label='additional discordant')
a.fill_between(e.esik,e.net_yanlis_alt,e.net_yanlis_ust,color=C['kotu'],alpha=0.15)
a.axhline(0,color='#495057',lw=0.9); a.axvline(0.80,ls=':',lw=1.0,color='#6c757d')
a.set_xlabel('Decision threshold'); a.set_ylabel('Additional cases over S1')
a.legend(fontsize=7.5,frameon=False)
plt.tight_layout(); plt.savefig(O+'Fig6.png',bbox_inches='tight'); plt.savefig(O+'Fig6.eps',bbox_inches='tight',format='eps'); plt.savefig(O+'Fig6.pdf',bbox_inches='tight'); plt.close()

# ---------- Sekil 4: karar egrisi
k=pd.read_csv('../../outputs/karar_egrisi.csv')
fig,ax=plt.subplots(figsize=(4.50,2.83))
ax.plot(k.pt,k.karar_kurali,'o-',color=C['ana'],lw=1.8,ms=4,label='Prediction rule')
ax.plot(k.pt,k.S1,'s-',color=C['iyi'],lw=1.5,ms=4,label='Subtype rule (S1)')
ax.plot(k.pt,k.hepsi,'^-',color=C['kotu'],lw=1.3,ms=4,label='Flag all')
ax.axhline(0,color='#495057',lw=1.0,ls='--',label='Flag none')
ax.axvline(0.80,ls=':',lw=1.0,color='#6c757d')
ax.set_xlabel('Threshold probability'); ax.set_ylabel('Net benefit'); ax.set_ylim(-1.2,0.42)
ax.legend(fontsize=7.5,frameon=False,loc='lower left')
plt.tight_layout(); plt.savefig(O+'Fig7.png',bbox_inches='tight'); plt.savefig(O+'Fig7.eps',bbox_inches='tight',format='eps'); plt.savefig(O+'Fig7.pdf',bbox_inches='tight'); plt.close()

# ---------- Sekil 5: kalibrasyon
p=np.load('../../outputs/oof_final.npy')
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
y=(d.RCB_Kategorize>=2).astype(int).values
q=pd.qcut(p,8,duplicates='drop')
g=pd.DataFrame({'p':p,'y':y,'q':q}).groupby('q',observed=True).agg(ort=('p','mean'),goz=('y','mean'),n=('y','size'))
fig,ax=plt.subplots(figsize=(3.30,2.70))
ax.plot([0,1],[0,1],ls='--',color='#adb5bd',lw=1.0)
ax.plot(g.ort,g.goz,'o-',color=C['ana'],lw=1.6,ms=5)
ax.set_xlabel('Predicted probability'); ax.set_ylabel('Observed proportion')
ax.set_xlim(0,1); ax.set_ylim(0,1)
ax.text(0.05,0.92,'Brier = 0,168',fontsize=8.5,color='#22333b')
plt.tight_layout(); plt.savefig(O+'Fig5.png',bbox_inches='tight'); plt.savefig(O+'Fig5.eps',bbox_inches='tight',format='eps'); plt.savefig(O+'Fig5.pdf',bbox_inches='tight'); plt.close()
print('5 sekil uretildi')
