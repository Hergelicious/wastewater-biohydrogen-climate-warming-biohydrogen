"""Monthly regional land warming (delta method) at 1.5, 2, 3 and 4 C global warming for each CMIP6 model, relative to 2001 to 2025,
and assignment of each plant to an IPCC AR6 reference region."""
import glob, os, numpy as np, pandas as pd, geopandas as gpd
from paths import ATLAS, WORK
A = str(ATLAS)
gwl = pd.read_csv(f'{A}/warming-levels/CMIP6_Atlas_WarmingLevels.csv')
def load(fn):
    d = pd.read_csv(fn, comment='#'); d['year'] = d.date.str[:4].astype(int); d['mon'] = d.date.str[5:7].astype(int); return d
src = f'{A}/datasets-aggregated-regionally/data/CMIP6/CMIP6_tas_land'
out = []
for _, r in gwl.iterrows():
    m = r.model_run
    fh = f'{src}/CMIP6_{m.split("_")[0]}_historical_{m.split("_")[1]}.csv'
    if not os.path.exists(fh): print('missing hist', m); continue
    h = load(fh)
    sc = {}
    for s in ['ssp245', 'ssp585', 'ssp370']:
        f = fh.replace('historical', s)
        if os.path.exists(f): sc[s] = load(f)
    if 'ssp245' not in sc: print('no ssp245', m); continue
    regs = [c for c in h.columns if c not in ('date', 'year', 'mon')]
    base = pd.concat([h[h.year >= 2001], sc['ssp245'][sc['ssp245'].year <= 2025]])
    bm = base.groupby('mon')[regs].mean()
    for lev in ['1.5', '2', '3', '4']:
        y, s = None, None
        for cand in (['ssp245', 'ssp585'] if lev != '4' else ['ssp585', 'ssp370']):
            v = r.get(f'{lev}_{cand}')
            if pd.notna(v) and v != 'NA' and cand in sc: y, s = int(v), cand; break
        if y is None: continue
        full = pd.concat([h, sc[s]]); w = full[(full.year >= y - 9) & (full.year <= y + 10)]
        if w.year.nunique() < 20: continue
        dm = w.groupby('mon')[regs].mean() - bm
        for mo in range(1, 13):
            for rg in regs:
                out.append((m, lev, s, y, mo, rg, dm.loc[mo, rg]))
dl = pd.DataFrame(out, columns=['model', 'gwl', 'ssp', 'year', 'mon', 'region', 'delta'])
dl.to_csv(WORK / 'gwl_deltas.csv', index=False)
print(dl.groupby('gwl').model.nunique())
print(dl[dl.region.isin(['MED', 'SAS', 'NEU', 'WAF', 'ARP'])].groupby(['gwl', 'region']).delta.median().unstack().round(2))
# plant -> region
reg = gpd.read_file(f'{A}/reference-regions/IPCC-WGI-reference-regions-v4.geojson')
print(reg.columns.tolist()[:6])
p = pd.read_csv(WORK / 'plants.csv')
pts = gpd.GeoDataFrame(p[['WASTE_ID']], geometry=gpd.points_from_xy(p.LON_WWTP, p.LAT_WWTP), crs=4326)
reg = reg[reg.Type.str.contains('Land')] if 'Type' in reg else reg
j = gpd.sjoin(pts, reg[['Acronym', 'geometry']], how='left', predicate='within').drop_duplicates('WASTE_ID')
miss = j.Acronym.isna()
if miss.any():
    jn = gpd.sjoin_nearest(pts[miss.values].to_crs(3857), reg[['Acronym', 'geometry']].to_crs(3857), how='left').drop_duplicates('WASTE_ID')
    j.loc[miss, 'Acronym'] = jn.set_index('WASTE_ID').loc[j.loc[miss, 'WASTE_ID'], 'Acronym'].values
print('nearest-assigned', miss.sum())
p['region'] = j.Acronym.values
p.to_csv(WORK / 'plants.csv', index=False)
print(p.groupby('region').POP_SERVED.sum().sort_values(ascending=False).head(12))
