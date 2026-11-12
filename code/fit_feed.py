"""Fit the fermenter feed-temperature relation to measured raw influent temperatures (weighted least squares, each site weighted
equally), choose the air-temperature clamp by leave-one-site-out cross-validation, and bootstrap the coefficients by site."""
import json, numpy as np, pandas as pd
from paths import VALID, OUT
EXCL = ['Linkoping WWTP', 'Jiangsu AAO WWTP', 'Gabal El-Asfar Cairo', 'Singapore WWTPs']  # not raw influent (model input, unspecified point, whole-system water, process water)
v_all = pd.read_csv(OUT / 'validation_table.csv'); v = v_all[~v_all.site.isin(EXCL)].reset_index(drop=True); clim = {k: np.array(x) for k, x in json.load(open(VALID / 'site_clim.json')).items()}
def months(spec): return list(range(1, 13)) if spec == '1-12' else [int(x) for x in str(spec).split(';')]
def feats(c_clamp):
    X = []
    for _, r in v.iterrows():
        c = np.maximum(clim[r.site], c_clamp); ms = [m-1 for m in months(r.months)]
        X.append([1.0, c.mean(), c[ms].mean() - c.mean()])
    return np.array(X)
y = v.T_obs.values; w = v.w.values
def wls(X, y, w): W = np.sqrt(w)[:, None]; return np.linalg.lstsq(X*W, y*np.sqrt(w), rcond=None)[0]
def wrmse(pred): return np.sqrt(np.sum(w*(pred - y)**2)/np.sum(w))
# current model
old = np.array([8 + 0.67*np.mean([clim[r.site][m-1] for m in months(r.months)]) for _, r in v.iterrows()])
print('current model: weighted RMSE', round(wrmse(old), 2), 'bias', round(np.sum(w*(old - y))/np.sum(w), 2))
best = None
for cc in [-50, 0, 2, 4, 6, 8]:
    X = feats(cc); b = wls(X, y, w); pred = X@b
    # leave-one-site-out
    cv = np.zeros_like(y)
    for s in v.site.unique():
        m = (v.site == s).values; bb = wls(X[~m], y[~m], w[~m]); cv[m] = X[m]@bb
    print(f'clamp {cc:>4}: alpha {b[0]:.2f} beta {b[1]:.3f} gamma {b[2]:.3f}  RMSE {wrmse(pred):.2f}  LOSO RMSE {wrmse(cv):.2f}')
    if best is None or wrmse(cv) < best[0]: best = (wrmse(cv), cc, b, pred, cv)
# simple linear a + b*Tair on monthly (no clamp) for comparison
X2 = np.array([[1, np.mean([clim[r.site][m-1] for m in months(r.months)])] for _, r in v.iterrows()])
b2 = wls(X2, y, w); print('simple linear fit a,b', b2.round(3), 'RMSE', round(wrmse(X2@b2), 2))
_, cc, b, pred, cv = best
# bootstrap by site for parameter uncertainty
rng = np.random.default_rng(3); S = v.site.unique(); X = feats(cc); bs = []
for _ in range(2000):
    pick = rng.choice(S, len(S), replace=True); idx = np.concatenate([np.where(v.site == s)[0] for s in pick])
    bs.append(wls(X[idx], y[idx], w[idx]))
bs = np.array(bs)
out = dict(clamp=cc, alpha=b[0], beta=b[1], gamma=b[2], boot=bs.tolist(), rmse=wrmse(pred), loso=wrmse(cv), old_rmse=wrmse(old))
# full 20-site fit for sensitivity
Xa = np.array([[1.0, np.maximum(clim[r.site], cc).mean(), np.maximum(clim[r.site], cc)[[m-1 for m in months(r.months)]].mean() - np.maximum(clim[r.site], cc).mean()] for _, r in v_all.iterrows()])
wa = (1/v_all.groupby('site').site.transform('count')).values
out['full_fit'] = wls(Xa, v_all.T_obs.values, wa).tolist()
json.dump(out, open(OUT / 'feed_fit.json', 'w'))
print('best clamp', cc, 'boot 5-95%', np.percentile(bs, [5, 95], axis=0).round(3))
v['pred_old'] = old; v['pred_new'] = pred; v['pred_loso'] = cv
pa = Xa@np.array(out['full_fit'])
v_all = v_all.drop(columns=[c for c in ['pred_old','pred_new','pred_loso','in_fit','pred_strict'] if c in v_all.columns])
v_all['in_fit'] = ~v_all.site.isin(EXCL)
v_all['pred_old'] = [8 + 0.67*np.mean([clim[r.site][m-1] for m in months(r.months)]) for _, r in v_all.iterrows()]
v_all['pred_strict'] = (np.array([[1.0, np.maximum(clim[r.site], cc).mean(), np.maximum(clim[r.site], cc)[[m-1 for m in months(r.months)]].mean() - np.maximum(clim[r.site], cc).mean()] for _, r in v_all.iterrows()])@b)
v_all = v_all.merge(v[['site','months','pred_loso']], on=['site','months'], how='left')
v_all.to_csv(OUT / 'validation_table.csv', index=False)
print(v[['site', 'months', 'T_obs', 'pred_old', 'pred_new', 'pred_loso']].round(1).to_string())
