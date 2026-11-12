"""Download NASA POWER monthly 2-m air temperature climatology (T2M, POWER climatology API) for every
10 x 10 degree tile that contains a HydroWASTE plant."""
import os, time, urllib.request, concurrent.futures as cf
import numpy as np, pandas as pd
from paths import RAW
out = RAW / 'power'; out.mkdir(exist_ok=True)
hw = pd.read_csv(RAW / 'HydroWASTE_v10.csv', encoding='latin-1')
tiles = pd.DataFrame({'tlat': (np.floor(hw.LAT_WWTP/10)*10).astype(int), 'tlon': (np.floor(hw.LON_WWTP/10)*10).astype(int)}).drop_duplicates()
def get(r):
    la, lo = int(r.tlat), int(r.tlon)
    fn = out / f't_{la}_{lo}.json'
    if fn.exists() and fn.stat().st_size > 1000: return fn, 'cached'
    la2 = min(la + 10, 90); lo2 = min(lo + 10, 180); la1 = max(la, -90); lo1 = max(lo, -180)
    if la2 - la1 < 2: la1 = la2 - 2
    if lo2 - lo1 < 2: lo1 = lo2 - 2
    url = (f'https://power.larc.nasa.gov/api/temporal/climatology/regional?latitude-min={la1}&latitude-max={la2}'
           f'&longitude-min={lo1}&longitude-max={lo2}&parameters=T2M&community=RE&format=JSON')
    err = ''
    for k in range(5):
        try:
            with urllib.request.urlopen(url, timeout=180) as f: b = f.read()
            fn.write_bytes(b); return fn, 'ok'
        except Exception as e:
            err = str(e); time.sleep(5*(k + 1))
    return fn, 'FAIL ' + err
with cf.ThreadPoolExecutor(4) as ex:
    for fn, s in ex.map(get, tiles.itertuples()):
        if s != 'cached': print(s, fn.name, flush=True)
print(len(tiles), 'tiles')
