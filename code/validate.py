"""Match each measured influent temperature (data/validation/obs.csv, sources in ESI Table S5) to NASA POWER air temperature at the site.
The site climatology used in the paper is shipped in data/validation/site_clim.json; set REFETCH = True to download it again."""
import json, urllib.request, time, numpy as np, pandas as pd
from paths import VALID, OUT
REFETCH = False
o = pd.read_csv(VALID / 'obs.csv')
sites = o[['site', 'lat', 'lon']].drop_duplicates('site')
M = ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']
clim = {}
cache = VALID / 'site_clim.json'
if cache.exists() and not REFETCH:
    clim = {k: np.array(v) for k, v in json.load(open(cache)).items()}
for _, s in sites.iterrows():
    if s.site in clim: continue
    url = f"https://power.larc.nasa.gov/api/temporal/climatology/point?parameters=T2M&community=RE&longitude={s.lon}&latitude={s.lat}&format=JSON"
    for k in range(4):
        try:
            d = json.load(urllib.request.urlopen(url, timeout=120)); break
        except Exception as e: time.sleep(4)
    t = d['properties']['parameter']['T2M']; clim[s.site] = np.array([t[m] for m in M])
json.dump({k: list(map(float, v)) for k, v in clim.items()}, open(cache, 'w'))
def months(spec): return list(range(1, 13)) if spec == '1-12' else [int(x) for x in str(spec).split(';')]
rows = []
for _, r in o.iterrows():
    c = clim[r.site]; ms = months(r.months); Tann = c.mean()
    rows.append(dict(site=r.site, country=r.country, type=r.type, T_obs=r.T_obs, Tair=c[[m-1 for m in ms]].mean(), Tann=Tann, w=None, months=r.months))
v = pd.DataFrame(rows); v['w'] = 1/v.groupby('site').site.transform('count')
v['dev'] = v.Tair - v.Tann
v.to_csv(OUT / 'validation_table.csv', index=False)
print(v.round(1).to_string())
