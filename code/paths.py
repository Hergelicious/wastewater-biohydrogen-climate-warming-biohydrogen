"""Folder layout used by all scripts (relative to the package root)."""
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / 'data' / 'raw'              # downloaded public inputs (get_data.py)
ATLAS = RAW / 'atlas'                    # IPCC AR6 WGI Atlas repository (sparse checkout)
VALID = ROOT / 'data' / 'validation'     # compiled measured influent temperatures and site climatology
WORK = ROOT / 'work'                     # intermediate files
OUT = ROOT / 'results'                   # results reported in the paper
for _d in (RAW, WORK, OUT): _d.mkdir(parents=True, exist_ok=True)
