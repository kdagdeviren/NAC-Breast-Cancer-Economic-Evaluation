# -*- coding: utf-8 -*-
"""Table 1 — Baseline characteristics, in English, for the manuscript."""
import pandas as pd, numpy as np, warnings; warnings.filterwarnings('ignore')
d0=pd.read_excel('../../data/Kohort_v17.xlsx')
d=d0[d0.Analiz_Kohortu=='taban senaryo (cM0)'].reset_index(drop=True)
elig=np.isin(d.Molekuler_i.values,['Luminal A','Luminal B (HER2 Negatif)'])
y=(d.RCB_Kategorize>=2).astype(int).values
G={'All cM0 records':np.ones(len(d),bool),'Eligible (HR+/HER2−)':elig,
   'Not eligible':~elig}
def satir(ad, fn):
    return [ad]+[fn(d[m], y[m]) for m in G.values()]
R=[]
def say(kol, deger, etiket):
    def f(dd, yy):
        n=(dd[kol].astype(str)==deger).sum()
        return f'{n} ({100*n/len(dd):.1f})'
    R.append(satir(etiket,f))
def surekli(kol, etiket):
    def f(dd, yy):
        v=pd.to_numeric(dd[kol],errors='coerce').dropna()
        return f'{v.median():.0f} ({v.quantile(.25):.0f}–{v.quantile(.75):.0f})'
    R.append(satir(etiket,f))
R.append(satir('Records, n', lambda dd,yy: f'{len(dd)}'))
R.append(satir('Patients, n', lambda dd,yy: f'{dd.Hasta_ID.nunique()}'))
R.append(['— Age —','','',''])
surekli('Tani_Yasi','Age at diagnosis, years, median (IQR)')
for g,e in [('Genç Erişkin','< 35'),('Erken Orta Yaş','35–44'),('Orta Yaş','45–54'),
            ('Geç Orta Yaş','55–64'),('Yaşlı','65–74'),('İleri Yaşlı','≥ 75')]:
    say('Yas_Grubu',g,f'  {e} years')
R.append(['— Tumour biology —','','',''])
for g in ['Luminal A','Luminal B (HER2 Negatif)','Luminal B (HER2 Pozitif)','Triple Negatif','HER2-Zengin']:
    EN={'Luminal A':'Luminal A','Luminal B (HER2 Negatif)':'Luminal B (HER2-negative)',
        'Luminal B (HER2 Pozitif)':'Luminal B (HER2-positive)','Triple Negatif':'Triple-negative',
        'HER2-Zengin':'HER2-enriched'}
    say('Molekuler_i',g,f'  {EN[g]}')
for g,e in [('Düşük','Low'),('Orta','Intermediate'),('Yüksek','High')]:
    say('Ki_67_i',g,f'  Ki-67 {e}')
R.append(['— Histological type —','','',''])
for g,e in [('İnvaziv Duktal Karsinom','Invasive ductal carcinoma'),
            ('İnvaziv Lobüler Karsinom','Invasive lobular carcinoma'),
            ('Mikst (Duktal + Lobüler)','Mixed ductal and lobular'),
            ('Diğer Nadir Tipler','Other rare types'),('DCIS','Ductal carcinoma in situ')]:
    say('Histolojik_Yeni',g,f'  {e}')
R.append(['— Nodal / regional involvement —','','',''])
for g,e in [('Yok','None'),('Lenf Nodu','Regional lymph node'),('Dermal Lenfatik','Dermal lymphatic')]:
    say('Metastaz_Yeri_i',g,f'  {e}')
R.append(['— Tumour-infiltrating lymphocytes —','','',''])
for g,e in [('<%10','< 10%'),('%10-%50','10–50%'),('>%50','> 50%')]:
    say('TIL_i',g,f'  {e}')
R.append(['— Histology and grade —','','',''])
for g in sorted(d.HistolojikG_i.astype(str).unique()):
    if g not in ('nan',): say('HistolojikG_i',g,f'  {g}')
R.append(['— Stage at diagnosis —','','',''])
for g in sorted(d.NAC_Once_Evre_i.astype(str).unique()):
    if g!='nan': say('NAC_Once_Evre_i',g,'  '+g.replace('Evre','Stage'))
R.append(['— Outcome —','','',''])
R.append(satir('RCB-II/III, n (%)', lambda dd,yy: f'{int(yy.sum())} ({100*yy.mean():.1f})'))
R.append(satir('RCB-0/I, n (%)', lambda dd,yy: f'{int((1-yy).sum())} ({100*(1-yy).mean():.1f})'))

T=pd.DataFrame(R, columns=['Characteristic']+list(G.keys()))
out='### Table 1 Baseline characteristics of the cohort ^a^\n\n'
out+='| Characteristic | '+' | '.join(f'{k}<br>(n = {int(m.sum())})' for k,m in G.items())+' |\n'
out+='|---|---|---|---|\n'
for _,r in T.iterrows():
    if r.Characteristic.startswith('—'):
        out+=f'| **{r.Characteristic.strip("— ")}** | | | |\n'
    else:
        out+=f'| {r.Characteristic} | {r.iloc[1]} | {r.iloc[2]} | {r.iloc[3]} |\n'
out+="""
^a^ Values are n (%) unless otherwise stated. The base-case cohort comprises all records with clinical M0 disease. The eligible population is defined by hormone receptor-positive / HER2-negative biology and is the population in which the prediction rule is intended to operate. Stage is as recorded in the clinical registry. RCB = residual cancer burden; IQR = interquartile range. Because 16 patients had bilateral disease, the number of records exceeds the number of patients.
"""
open('../../outputs/Table1_baseline_EN.md','w').write(out)
print(out)
