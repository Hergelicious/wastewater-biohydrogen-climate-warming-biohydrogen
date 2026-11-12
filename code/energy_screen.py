"""Screening energy balance for unheated dark fermentation, by plant output class (central case, present climate)."""
import numpy as np, pandas as pd
import core as C
from paths import OUT
per = pd.read_csv(OUT / 'plants_windows.csv')
assert len(per) == len(C.pop)
q = C.CENTRAL
fw_t = C.fw_kg*q['cap']*C.pop/1000                     # t wet food waste per year
h2 = per.h2_now_kg.values                               # kg H2 per year (unheated, central)
LHV = 33.3                                              # kWh per kg H2 (lower heating value)
NM3 = 11.1                                              # Nm3 per kg H2
TS_FEED = 0.10                                          # assumed slurry solids after dilution with sludge or effluent
slurry_m3 = fw_t*q['fw_ts']/TS_FEED                     # m3 slurry per year (density ~1 t/m3)
vol = slurry_m3/365*q['hrt']                            # reactor working volume, m3
mix_kwh = vol*5.0*8760/1000                             # 5 W/m3 continuous mixing (US EPA lower bound)
upg_lo, upg_hi = 0.296*NM3, 0.400*NM3                   # kWh per kg H2, membrane upgrading (two-stage)
PARASITIC = 59.1                                        # kWh/t, rest-of-plant electricity at a full-scale food-waste AD plant
gross = h2*LHV
lean = mix_kwh + h2*upg_lo
full = fw_t*PARASITIC + h2*upg_hi
cls = pd.cut(per.kgd_q_now, [1, 10, 100, np.inf], right=False, labels=['1 to 10', '10 to 100', '≥100'])
rows = []
for c in cls.cat.categories:
    m = (cls == c).values
    rows.append({'Output class (kg H2/d, warm months)': c, 'Plants': int(m.sum()),
                 'Median food waste (t/yr)': np.median(fw_t[m]), 'Median H2 (t/yr)': np.median(h2[m])/1000,
                 'H2 energy (kWh/t food waste)': gross[m].sum()/fw_t[m].sum(),
                 'Mixing + upgrading (kWh/t)': lean[m].sum()/fw_t[m].sum(),
                 'Full-plant auxiliaries (kWh/t)': full[m].sum()/fw_t[m].sum()})
t = pd.DataFrame(rows)
t['Ratio, lean'] = t['H2 energy (kWh/t food waste)']/t['Mixing + upgrading (kWh/t)']
t['Ratio, full plant'] = t['H2 energy (kWh/t food waste)']/t['Full-plant auxiliaries (kWh/t)']
t.to_csv(OUT / 'energy_screen.csv', index=False)
pd.set_option('display.width', 250); print(t.round(2).to_string(index=False))
print('all plants: H2 kWh/t', gross.sum()/fw_t.sum(), 'mix kWh/t', mix_kwh.sum()/fw_t.sum())
