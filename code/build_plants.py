"""Assemble the plant inventory: HydroWASTE plants in use, monthly air temperature at the nearest POWER grid cell, and World Bank income group."""
import json, glob, numpy as np, pandas as pd
from scipy.spatial import cKDTree
from paths import RAW, WORK
hw = pd.read_csv(RAW / 'HydroWASTE_v10.csv', encoding='latin-1')
keep = ['Not Reported', 'Operational', 'Construction Completed']
hw = hw[hw.STATUS.isin(keep) & (hw.POP_SERVED > 0)].copy()
# climate grid
rows = []
for fn in sorted(glob.glob(str(RAW / 'power' / '*.json'))):
    for f in json.load(open(fn))['features']:
        lon, lat = f['geometry']['coordinates'][:2]
        t = f['properties']['parameter']['T2M']
        m = [t[k] for k in ['JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC']]
        if min(m) < -900: continue
        rows.append([lat, lon] + m)
g = pd.DataFrame(rows, columns=['lat','lon'] + [f'T{i}' for i in range(1,13)]).drop_duplicates(['lat','lon'])
def xyz(lat, lon):
    la, lo = np.radians(lat), np.radians(lon)
    return np.c_[np.cos(la)*np.cos(lo), np.cos(la)*np.sin(lo), np.sin(la)]
tree = cKDTree(xyz(g.lat.values, g.lon.values))
dist, idx = tree.query(xyz(hw.LAT_WWTP.values, hw.LON_WWTP.values))
hw['clim_dist_km'] = dist*6371
for i in range(1,13): hw[f'Tair{i}'] = g[f'T{i}'].values[idx]
# income group
wb = json.load(open(RAW / 'wb_income.json'))[1]
inc = {c['id']: c['incomeLevel']['id'] for c in wb}
hw['income'] = hw.CNTRY_ISO.map(inc)
print('unmatched income:', hw[hw.income.isna() | (hw.income=='INX')].CNTRY_ISO.value_counts().head(15).to_dict())
hw.to_csv(WORK / 'plants.csv', index=False)
print(len(hw), 'plants;', round(hw.POP_SERVED.sum()/1e9,3), 'bn people;', 'max clim dist km', round(hw.clim_dist_km.max(),1))
print(hw.groupby('income').POP_SERVED.agg(['count','sum']))
