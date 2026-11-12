"""Download the public inputs: HydroWASTE v1.0, World Bank income groups and electricity access,
and the IPCC AR6 WGI Atlas regional CMIP6 series (sparse git checkout at a fixed commit)."""
import subprocess, urllib.request
from paths import RAW, ATLAS
FILES = {
    # HydroWASTE v1.0 (Ehalt Macedo et al. 2022, CC BY 4.0); identical copy of the figshare release
    'HydroWASTE_v10.csv': 'https://raw.githubusercontent.com/rfordatascience/tidytuesday/master/data/2022/2022-09-20/HydroWASTE_v10.csv',
    'wb_income.json': 'https://api.worldbank.org/v2/country?format=json&per_page=400',
    'wb_elec.json': 'https://api.worldbank.org/v2/country/all/indicator/EG.ELC.ACCS.ZS?format=json&mrnev=1&per_page=20000',
}
for name, url in FILES.items():
    fn = RAW / name
    if not fn.exists():
        print('downloading', name); urllib.request.urlretrieve(url, fn)
COMMIT = 'ed29da33d38ae767d19caafd4cde06042cbe00f4'
PATHS = ['datasets-aggregated-regionally/data/CMIP6/CMIP6_tas_land', 'reference-regions', 'warming-levels']
if not ATLAS.exists():
    subprocess.run(['git', 'clone', '--filter=blob:none', '--sparse', 'https://github.com/IPCC-WG1/atlas', str(ATLAS)], check=True)
    subprocess.run(['git', '-C', str(ATLAS), 'sparse-checkout', 'set'] + PATHS, check=True)
    subprocess.run(['git', '-C', str(ATLAS), 'checkout', COMMIT], check=True)
print('inputs ready in', RAW)
