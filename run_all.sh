#!/usr/bin/env bash
# Reproduce all results from public inputs. Run from the package root: bash run_all.sh
set -euo pipefail
cd "$(dirname "$0")/code"
python get_data.py        # HydroWASTE, World Bank, IPCC AR6 Atlas (git)
python fetch_power.py     # NASA POWER T2M climatology tiles
python build_plants.py    # plant inventory with monthly air temperature and income group
python deltas.py          # CMIP6 warming-level deltas and AR6 region of each plant
python validate.py        # measured influent temperatures matched to site air temperature
python fit_feed.py        # feed-temperature relation, cross-validation, bootstrap
python core.py            # windows, hydrogen, model spread, Monte Carlo
python ghg_sdg.py         # greenhouse-gas screening, robustness, country table with SDG 7.1
python sens_kinetics.py   # fermenter offset, temperature-penalty structure, retention time
python energy_screen.py   # screening energy balance by plant output class
echo "Done. Results are in results/."
