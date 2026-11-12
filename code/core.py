"""Core analysis: unheated dark-fermentation (DF) operating windows at HydroWASTE plants, today and at global warming levels."""
import json, numpy as np, pandas as pd
from paths import WORK, OUT
FF = json.load(open(OUT / 'feed_fit.json')); BOOT = np.array(FF['boot'])
p = pd.read_csv(WORK / 'plants.csv'); p.loc[p.CNTRY_ISO == 'TWN', 'income'] = 'HIC'; p['income'] = p.income.fillna('UMC')
dl = pd.read_csv(WORK / 'gwl_deltas.csv', dtype={'gwl': str})
n = len(p); pop = p.POP_SERVED.values.astype(float)
Tair = p[[f'Tair{i}' for i in range(1, 13)]].values
DAYS = np.array([31,28,31,30,31,30,31,31,30,31,30,31])
fw_kg = p.income.map({'HIC': 81., 'UMC': 88., 'LMC': 86., 'LIC': 86.}).values   # UNEP FWI 2024 household, kg/cap/yr
hic = (p.income == 'HIC').values
L_PER_KG = 11126.  # L H2 per kg at 0 C and 1 atm
GWLS = ['1.5', '2', '3', '4']
regs = p.region.values; ureg = np.unique(regs); ridx = np.searchsorted(ureg, regs)
models = {g: sorted(dl[dl.gwl.astype(str) == g].model.unique()) for g in GWLS}
def delta_matrix(g, model=None):
    d = dl[dl.gwl.astype(str) == g]
    d = d[d.model == model] if model else d.groupby(['region', 'mon']).delta.median().reset_index()
    t = d.pivot_table(index='region', columns='mon', values='delta').reindex(ureg)
    t = t.fillna(t.median())  # regions missing in a model: fill with model median
    return t.values[ridx]                                            # (n,12)
CENTRAL = dict(fa=FF['alpha'], fb=FF['beta'], fg=FF['gamma'], clamp=FF['clamp'], a=None, b=None, T0=12.5, T1=31.5, theta=1.07, hrt=2.0, eta37=0.9, cap=0.25, fw_ts=0.25, fw_vsts=0.90, y_fw=55., ps_hic=35., ps_lmic=25., ps_vsts=0.7, y_ps=3.)
def draw(rng):
    """One Monte Carlo parameter set. Draws are taken in a fixed order so that the random stream (seed 7) reproduces the published runs."""
    tri = rng.triangular
    bb = BOOT[rng.integers(len(BOOT))]                 # bootstrap set of feed-temperature coefficients
    T0 = rng.uniform(10, 15); T1 = rng.uniform(28, 35)
    theta = rng.uniform(1.05, 1.10)
    rng.uniform(1.0, 3.0)   # retention-time draw kept for stream compatibility; HRT is fixed at the 2-day design because
                            # k37 is set from eta37 at the design HRT, so a sampled HRT has no independent effect
    eta37 = rng.uniform(0.85, 0.95)
    cap = tri(0.10, 0.25, 0.50); fw_ts = rng.uniform(0.2, 0.3); fw_vsts = rng.uniform(0.85, 0.95); y_fw = tri(30, 55, 100)
    ps_hic = tri(25, 35, 50); ps_lmic = tri(15, 25, 35); ps_vsts = rng.uniform(0.6, 0.8); y_ps = tri(0.5, 3, 10)
    return dict(fa=bb[0], fb=bb[1], fg=bb[2], clamp=FF['clamp'], a=None, b=None, T0=T0, T1=T1, theta=theta, hrt=2.0, eta37=eta37,
                cap=cap, fw_ts=fw_ts, fw_vsts=fw_vsts, y_fw=y_fw, ps_hic=ps_hic, ps_lmic=ps_lmic, ps_vsts=ps_vsts, y_ps=y_ps)
def evaluate(q, dT):
    if q.get('a') is not None:                       # linear sensitivity variant
        Tf = np.maximum(q['a'] + q['b']*(Tair + dT), 4.0)
    else:                                            # validated relation: annual mean + damped seasonal deviation
        Ta = np.maximum(Tair + dT, q['clamp']); Tm = Ta.mean(1, keepdims=True)
        Tf = q['fa'] + q['fb']*Tm + q['fg']*(Ta - Tm)
    Tf = Tf + q.get('offset', 0.0)                   # fermenter minus influent temperature (sensitivity)
    f = np.clip((Tf - q['T0'])/(q['T1'] - q['T0']), 0, 1)
    # first-order CSTR conversion at fixed design HRT; k37 set so that conversion at 37 C equals eta37
    k37 = q['eta37']/(1 - q['eta37'])/q.get('hrt_ref', q['hrt'])   # rate fixed by design HRT
    k = k37*q['theta']**(Tf - 37.0)
    eta = k*q['hrt']/(1 + k*q['hrt'])
    mode = q.get('mode', 'both')                     # temperature penalty: yield factor, kinetics, or both
    if mode == 'yield': eta = np.full_like(eta, q['eta37'])
    elif mode == 'kin': f = np.ones_like(f)
    f = f*eta
    vs = fw_kg*1000/365*q['cap']*q['fw_ts']*q['fw_vsts']*q['y_fw'] + np.where(hic, q['ps_hic'], q['ps_lmic'])*q['ps_vsts']*q['y_ps']
    h2m = vs[:, None]/1000*f*pop[:, None]*DAYS                       # L H2 per month per plant
    h2 = h2m.sum(1)                                                  # L H2/yr per plant
    h2max = vs/1000*pop*365                                          # batch potential at 37 C, full conversion
    evaluate.h2m = h2m; evaluate.q = (Tf >= 25)
    return Tf, f, h2, h2max
def summarise(Tf, h2, h2max):
    o = {}
    for th in (20, 25, 30):
        mo = (Tf >= th).sum(1)
        o[f'm{th}'] = mo
    o['h2'] = h2; o['h2max'] = h2max
    return o
TH = (20, 25, 30)
def agg(o):
    r = {'H2_kt': o['h2'].sum()/L_PER_KG/1e6, 'H2_frac': o['h2'].sum()/o['h2max'].sum(),
         'H2q_kt': (evaluate.h2m*evaluate.q).sum()/L_PER_KG/1e6}
    for th in TH:
        m = o[f'm{th}']
        r[f'pop_any{th}'] = pop[m >= 1].sum()/1e6; r[f'pop_6m{th}'] = pop[m >= 6].sum()/1e6; r[f'pop_12m{th}'] = pop[m == 12].sum()/1e6
        r[f'plants_any{th}'] = int((m >= 1).sum()); r[f'pop_months{th}'] = (m*pop).sum()/pop.sum()
    return r
if __name__ == '__main__':
    res = {'central': {}, 'model_spread': {}, 'mc': {}}
    # 1. central params, ensemble-median deltas -> per-plant fields
    zero = np.zeros_like(Tair)
    Tf, f, h2, h2max = evaluate(CENTRAL, zero); base = summarise(Tf, h2, h2max)
    res['central']['present'] = agg(base)
    h2q = (evaluate.h2m*evaluate.q).sum(1)/L_PER_KG; qdays = (evaluate.q*DAYS).sum(1)
    per = pd.DataFrame({'WASTE_ID': p.WASTE_ID, 'm25_now': base['m25'], 'm30_now': base['m30'], 'm20_now': base['m20'],
                        'h2_now_kg': h2/L_PER_KG, 'h2max_kg': h2max/L_PER_KG, 'Tf_mean_now': (Tf*DAYS).sum(1)/365,
                        'h2q_now_kg': h2q, 'kgd_q_now': np.where(qdays > 0, h2q/np.maximum(qdays, 1), 0.0),
                        'kgd_annual_now': h2/L_PER_KG/365})
    for g in GWLS:
        Tf, f, h2, _ = evaluate(CENTRAL, delta_matrix(g)); o = summarise(Tf, h2, h2max)
        res['central'][g] = agg(o)
        per[f'm25_{g}'] = o['m25']; per[f'm30_{g}'] = o['m30']; per[f'h2_{g}_kg'] = h2/L_PER_KG
        h2q_g = (evaluate.h2m*evaluate.q).sum(1)/L_PER_KG; qd_g = (evaluate.q*DAYS).sum(1)
        per[f'kgd_q_{g}'] = np.where(qd_g > 0, h2q_g/np.maximum(qd_g, 1), 0.0)
        # model agreement: share of models in which plant gains >=1 month at 25 C
        gain = np.zeros(n); spread = []
        for m in models[g]:
            Tf_m, _, h2_m, _ = evaluate(CENTRAL, delta_matrix(g, m)); om = summarise(Tf_m, h2_m, h2max)
            gain += (om['m25'] > base['m25']); spread.append(agg(om))
        per[f'agree25_{g}'] = gain/len(models[g])
        sp = pd.DataFrame(spread)
        res['model_spread'][g] = {k: [float(sp[k].quantile(.05)), float(sp[k].median()), float(sp[k].quantile(.95))] for k in sp.columns}
    per.to_csv(OUT / 'plants_windows.csv', index=False)
    # 2. full Monte Carlo: process params x random climate model
    rng = np.random.default_rng(7); N = 400; rows = []
    for i in range(N):
        q = draw(rng)
        Tf, f, h2, h2max = evaluate(q, zero); r0 = agg(summarise(Tf, h2, h2max)); r0['gwl'] = 'present'; r0['i'] = i; rows.append(r0)
        for g in GWLS:
            m = models[g][rng.integers(len(models[g]))]
            Tf, f, h2, _ = evaluate(q, delta_matrix(g, m)); r = agg(summarise(Tf, h2, h2max)); r['gwl'] = g; r['i'] = i; rows.append(r)
    mc = pd.DataFrame(rows); mc.to_csv(OUT / 'mc_runs.csv', index=False)
    # paired change vs present
    for g in GWLS:
        a = mc[mc.gwl == g].set_index('i'); b = mc[mc.gwl == 'present'].set_index('i')
        res['mc'][g] = {k: [float(np.percentile(a[k], q)) for q in (5, 50, 95)] for k in a.columns if k != 'gwl'}
        res['mc'][g]['dH2_pct'] = [float(np.percentile((a.H2_kt/b.H2_kt - 1)*100, q)) for q in (5, 50, 95)]
        res['mc'][g]['dpop_any25'] = [float(np.percentile(a.pop_any25 - b.pop_any25, q)) for q in (5, 50, 95)]
    b = mc[mc.gwl == 'present']
    res['mc']['present'] = {k: [float(np.percentile(b[k], q)) for q in (5, 50, 95)] for k in b.columns if k not in ('gwl',)}
    json.dump(res, open(OUT / 'core_results.json', 'w'), indent=1)
    c = res['central']
    for k in ['H2_kt', 'H2q_kt', 'H2_frac', 'pop_any25', 'pop_6m25', 'pop_12m25', 'pop_any30', 'pop_12m30', 'pop_months25', 'plants_any25']:
        print(k, round(c['present'][k], 3), [round(c[g][k], 3) for g in GWLS])
    print('MC dH2%', {g: [round(x, 1) for x in res['mc'][g]['dH2_pct']] for g in GWLS})
    print('MC H2 present', [round(x, 1) for x in res['mc']['present']['H2_kt']])
