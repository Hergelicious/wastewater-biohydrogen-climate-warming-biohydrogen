"""Sensitivity: fermenter offset from influent, temperature-penalty structure, HRT range (central parameters)."""
import numpy as np, pandas as pd
import core as C
from paths import OUT
Z = np.zeros_like(C.Tair); G2 = C.delta_matrix('2')
def run(lab, grp, **kw):
    out = []
    for gl, dT in [('present', Z), ('2', G2)]:
        q = dict(C.CENTRAL, **kw)
        Tf, f, h2, h2m = C.evaluate(q, dT); r = C.agg(C.summarise(Tf, h2, h2m)); r.update(group=grp, case=lab, gwl=gl); out.append(r)
    return out
rows = []
for o in (-3, -2, 0, 2, 3): rows += run(f'{o:+d} °C' if o else 'as influent (central)', 'Fermenter vs influent', offset=float(o))
for m, lab in [('both', 'yield factor and kinetics (central)'), ('yield', 'yield factor only'), ('kin', 'kinetics only')]: rows += run(lab, 'Temperature penalty', mode=m)
for h in (1, 2, 3, 4, 6): rows += run(f'{h} d' + (' (central)' if h == 2 else ''), 'Retention time (rate fixed at 2 d design)', hrt=float(h), hrt_ref=2.0)
r = pd.DataFrame(rows); r.to_csv(OUT / 'sens_kinetics.csv', index=False)
w = r.pivot_table(index=['group', 'case'], columns='gwl', values=['pop_any25', 'pop_12m25', 'H2_kt', 'H2q_kt'], sort=False)
w[('dH2', '%')] = (w[('H2_kt', '2')]/w[('H2_kt', 'present')] - 1)*100
print(w.round(1).to_string())
