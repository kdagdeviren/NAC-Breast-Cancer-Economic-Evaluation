# -*- coding: utf-8 -*-
"""Regularised logistic regression alone versus the two-layer model, both under the
fully nested protocol (Section 3.8, Online Resource 4 Table S4.16).

Requires tam_icice.py in the same directory (its header block defines the data,
the modal imputer and the log-likelihood-ratio screen).
"""
import numpy as np, pandas as pd, warnings; warnings.filterwarnings('ignore')
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, brier_score_loss

src = open('tam_icice.py', encoding='utf-8').read()
exec(src[:src.index('S=list(pd.read_csv')])          # data, imp_modal, llr_p, BLOK, VAR, y, grp, elig, lumA

def fit_lr(Xi, vars_, tr, te, yy):
    du = {v: pd.get_dummies(Xi[v], drop_first=True).values.astype(float) for v in vars_}
    sel = [v for v in vars_ if llr_p(du[v][tr], yy[tr]) < 0.10] or list(vars_)
    m = Pipeline([('oh', OneHotEncoder(handle_unknown='ignore', drop='first', sparse_output=False)),
                  ('m', LogisticRegression(penalty='l2', C=0.05, max_iter=5000))]).fit(Xi[sel].iloc[tr], yy[tr])
    return m.predict_proba(Xi[sel].iloc[te])[:, 1]

KOMB = [('P',), ('P', 'O'), ('P', 'R'), ('P', 'O', 'R')]
oof = np.full(len(y), np.nan)
for r in range(3):                                    # seeds 700-702, as in the main protocol
    for tr, te in StratifiedGroupKFold(5, shuffle=True, random_state=700+r).split(X0, y, groups=grp):
        Xi = imp_modal(X0, tr)
        best, best_auc = None, -1
        for k in KOMB:                                # block selection inside the outer fold
            vs = [v for v in sum([BLOK[b] for b in k], []) if v in VAR]
            ip = np.full(len(tr), np.nan)
            Xtr = X0.iloc[tr].reset_index(drop=True); ytr = y[tr]; gtr = grp[tr]
            for itr, ite in StratifiedGroupKFold(3, shuffle=True, random_state=11).split(Xtr, ytr, groups=gtr):
                ip[ite] = fit_lr(imp_modal(Xtr, itr), vs, itr, ite, ytr)
            a = roc_auc_score(ytr, ip)
            if a > best_auc: best_auc, best = a, k
        vs = [v for v in sum([BLOK[b] for b in best], []) if v in VAR]
        oof[te] = fit_lr(Xi, vs, tr, te, y)
np.save('oof_logistic_only.npy', oof)

def slope(p, yy, g, B=800, seed=2026):
    lp = np.log(np.clip(p, 1e-6, 1-1e-6) / (1 - np.clip(p, 1e-6, 1-1e-6))).reshape(-1, 1)
    m = LogisticRegression(penalty=None, max_iter=1000).fit(lp, yy)
    rng = np.random.default_rng(seed); pts = np.unique(g); E = []
    for _ in range(B):
        idx = np.concatenate([np.where(g == h)[0] for h in rng.choice(pts, len(pts), replace=True)])
        if len(np.unique(yy[idx])) < 2: continue
        E.append(LogisticRegression(penalty=None, max_iter=1000).fit(lp[idx], yy[idx]).coef_[0][0])
    return float(m.coef_[0][0]), np.percentile(E, 2.5), np.percentile(E, 97.5)

def summarise(K, label, thr=0.80, matched=None):
    if matched:                                        # compare at a matched flagging rate
        thr = np.sort(K[elig])[::-1][matched-1]
    m = elig & (K >= thr); tp = int(y[m].sum()); n = int(m.sum())
    eg, lo, hi = slope(K[elig], y[elig], grp[elig])
    pr = y[elig].mean(); bs = brier_score_loss(y[elig], K[elig])
    return dict(model=label, threshold=round(thr, 3), AUC=round(roc_auc_score(y[elig], K[elig]), 3),
                Brier=round(bs, 3), scaled_Brier=round(1 - bs/(pr*(1-pr)), 3),
                slope=round(eg, 2), slope_CI=f'{lo:.2f}-{hi:.2f}', flagged=n, concordant=tp,
                discordant=n-tp, PPV=round(tp/n, 3), delta_vs_S1=tp - int(y[lumA].sum()))

nested = np.load('oof_tam_icice.npy')
out = pd.DataFrame([summarise(nested, 'Combined (two-layer)'),
                    summarise(oof, 'Logistic only'),
                    summarise(oof, 'Logistic only, matched flagging rate', matched=113)])
print(out.to_string(index=False))
out.to_csv('../outputs/simple_model_rival.csv', index=False)
