"""Screening greenhouse-gas balance, robustness to the feed-temperature relation, and country table with SDG 7.1 electricity access."""
import json, numpy as np, pandas as pd
import core as C
from paths import RAW, OUT
rng = np.random.default_rng(11)
p = C.p; pop = C.pop; per = pd.read_csv(OUT / 'plants_windows.csv')
# --- screening GHG: hydrogen-attributable credit vs diversion credit (reported separately)
rows = []
mc = pd.read_csv(OUT / 'mc_runs.csv')
for i in range(400):
    q = C.draw(rng)
    fw_t = (C.fw_kg*q['cap']*pop).sum()/1e3                        # t wet FW/yr diverted
    h2 = mc[(mc.gwl == 'present') & (mc.i == i)].H2_kt.values[0]*1e3  # t H2/yr (paired draw index)
    ef_h2 = rng.uniform(9.9, 11.9); ef_lf = rng.triangular(0.3, 0.6, 1.0)
    rows.append(dict(fw_Mt=fw_t/1e6, h2_credit_Mt=h2*ef_h2/1e6, landfill_credit_Mt=fw_t*ef_lf/1e6))
g = pd.DataFrame(rows); g['ratio'] = g.landfill_credit_Mt/g.h2_credit_Mt
ghg = {k: [float(g[k].quantile(x)) for x in (.05, .5, .95)] for k in g.columns}
ghg['fw_Mt_central'] = float((C.fw_kg*C.CENTRAL['cap']*pop).sum()/1e9)
ghg['h2_credit_Mt_central'] = float(pd.read_csv(OUT / 'plants_windows.csv').h2_now_kg.sum()/1e3*10.9/1e6)
# --- robustness (central params): air-to-feed variants and thresholds, present and +2 C
rob = []
import numpy as _np
_p5, _p95 = _np.percentile(C.BOOT, [5, 95], axis=0)
_lo = C.BOOT[_np.argmin(_np.abs(C.BOOT[:, 0] + 15*C.BOOT[:, 1] - _np.percentile(C.BOOT[:, 0] + 15*C.BOOT[:, 1], 5)))]
_hi = C.BOOT[_np.argmin(_np.abs(C.BOOT[:, 0] + 15*C.BOOT[:, 1] - _np.percentile(C.BOOT[:, 0] + 15*C.BOOT[:, 1], 95)))]
for lab, a, b, fx in [('validated (central fit)', None, None, None), ('validated, cool fit (5th pct)', None, None, _lo),
                      ('validated, warm fit (95th pct)', None, None, _hi), ('validated, all 20 sites', None, None, C.FF['full_fit']), ('linear 8 + 0.67 T_air', 8, .67, None), ('feed = air', 0, 1, None)]:
    q = dict(C.CENTRAL, a=a, b=b)
    if fx is not None: q.update(fa=fx[0], fb=fx[1], fg=fx[2])
    for gl, dT in [('present', np.zeros_like(C.Tair)), ('2', C.delta_matrix('2'))]:
        Tf, f, h2, h2m = C.evaluate(q, dT); r = C.agg(C.summarise(Tf, h2, h2m)); r.update(variant=lab, gwl=gl); rob.append(r)
rob = pd.DataFrame(rob); rob.to_csv(OUT / 'robustness.csv', index=False)
# --- country overlay with SDG 7.1 (electricity access, World Bank)
el = {d['countryiso3code']: d['value'] for d in json.load(open(RAW / 'wb_elec.json'))[1] if d['value'] is not None}
x = p[['CNTRY_ISO', 'COUNTRY', 'income', 'POP_SERVED']].join(per.drop(columns='WASTE_ID'))
cty = x.groupby(['CNTRY_ISO', 'COUNTRY', 'income']).apply(lambda d: pd.Series({
    'plants': len(d), 'pop_M': d.POP_SERVED.sum()/1e6,
    'months25_now': (d.m25_now*d.POP_SERVED).sum()/d.POP_SERVED.sum(), 'months25_2C': (d['m25_2']*d.POP_SERVED).sum()/d.POP_SERVED.sum(),
    'months25_4C': (d['m25_4']*d.POP_SERVED).sum()/d.POP_SERVED.sum(),
    'h2_now_t': d.h2_now_kg.sum()/1e3, 'h2_2C_t': d['h2_2_kg'].sum()/1e3, 'h2max_t': d.h2max_kg.sum()/1e3}), include_groups=False).reset_index()
cty['elec_access'] = cty.CNTRY_ISO.map(el)
cty.to_csv(OUT / 'country_table.csv', index=False)
json.dump(dict(ghg=ghg), open(OUT / 'ghg_results.json', 'w'), indent=1)
print(json.dumps(ghg, indent=0))
print(rob[['variant', 'gwl', 'H2_kt', 'pop_any25', 'pop_12m25', 'pop_months25']].round(1).to_string())
big = cty[cty.pop_M > 5].sort_values('months25_2C', ascending=False)
print(big[['COUNTRY', 'income', 'plants', 'pop_M', 'months25_now', 'months25_2C', 'months25_4C', 'h2_now_t', 'elec_access']].round(1).head(25).to_string())
print('low access (<95%) countries with plants:', cty[cty.elec_access < 95][['COUNTRY', 'plants', 'pop_M', 'months25_now', 'months25_2C', 'elec_access']].round(1).sort_values('pop_M', ascending=False).head(20).to_string())
