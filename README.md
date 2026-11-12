# Unheated dark fermentation of household food waste at wastewater treatment plants

Analysis code for the manuscript *Climate warming expands operating windows for unheated biohydrogen production at wastewater treatment plants*.

The code works out, for each of 56,506 wastewater treatment plants in HydroWASTE, the months in which an unheated dark fermentation reactor fed with household food waste and primary sludge would run warm enough for hydrogen production. It does this for the present climate and for global warming levels of 1.5, 2, 3 and 4 °C. It then estimates hydrogen output and its uncertainty, runs a screening greenhouse-gas balance and a screening energy balance, and overlays country electricity access (SDG 7.1).

The package holds the analysis code only. Code for figures, the graphical abstract and document formatting is not included. Every number in the manuscript, the ESI tables and the figure data comes from the files the scripts write to `results/`.

## Contents

```
run_all.sh                 runs the full pipeline from download to results
requirements.txt           Python package versions used
code/
  paths.py                 folder locations, relative to the package root
  get_data.py              downloads HydroWASTE, World Bank data and the IPCC AR6 Atlas
  fetch_power.py           downloads NASA POWER air temperature climatology
  build_plants.py          plant inventory with monthly air temperature and income group
  deltas.py                CMIP6 monthly warming at each warming level; AR6 region of each plant
  validate.py              matches measured influent temperatures to site air temperature
  fit_feed.py              fits the feed temperature relation; cross-validation; bootstrap
  core.py                  operating windows, hydrogen output, model spread, Monte Carlo
  ghg_sdg.py               greenhouse-gas screen, robustness runs, country table with SDG 7.1
  sens_kinetics.py         fermenter offset, temperature penalty structure, retention time
  energy_screen.py         screening energy balance by plant output class
data/validation/
  obs.csv                  measured raw influent temperatures (38 records, 20 sites; sources in the ESI)
  site_clim.json           NASA POWER monthly climatology at the 20 validation sites
reference_results/         the outputs used in the manuscript, for comparison with a re-run
```

Running the pipeline creates `data/raw/` (downloaded inputs), `work/` (intermediate files) and `results/` (outputs).

## How to run

Requirements: Python 3.13 with the packages in `requirements.txt`, `git`, and internet access for the first run.

```
pip install -r requirements.txt
bash run_all.sh
```

The scripts can also be run one at a time from inside `code/`, in the order listed in `run_all.sh`. Each step reads only what the earlier steps wrote. Downloads are skipped if the files already exist in `data/raw/`. The first run downloads about 220 NASA POWER tiles and a sparse copy of the IPCC Atlas repository; after that the analysis takes a few minutes on a laptop.

## Input data

| Input | Version or access | Licence | Used for |
|---|---|---|---|
| HydroWASTE (Ehalt Macedo et al., 2022, Earth Syst. Sci. Data 14, 559) | v1.0, `HydroWASTE_v10.csv`, md5 `585b1d51401f09592963abd283f933db`, taken from the TidyTuesday mirror of the figshare release | CC BY 4.0 | plant locations, population served, status |
| NASA POWER climatology API | `temporal/climatology`, parameter T2M, community RE, 10° × 10° regional requests | public, no restriction | monthly 2 m air temperature |
| World Bank API | country metadata (income group); indicator EG.ELC.ACCS.ZS, most recent non-empty value | CC BY 4.0 | income group; SDG 7.1 electricity access |
| IPCC AR6 WGI Atlas (Iturbide et al.) | github.com/IPCC-WG1/atlas, commit `ed29da33d38ae767d19caafd4cde06042cbe00f4` | CC BY 4.0 | CMIP6 regional land temperature, warming level years, reference regions |
| Measured influent temperatures | `data/validation/obs.csv`, compiled from the literature (full sources in the ESI) | values cited from published sources | fitting and testing the feed temperature relation |

Plants are kept if their status is operational, construction completed or not reported and they serve more than zero people (56,506 plants serving 2.27 billion people). Each plant takes the monthly air temperature of the nearest POWER grid cell. Taiwan is assigned high income, and plants without a World Bank income group are assigned upper middle income.

## Method in brief

**Climate change.** For each CMIP6 model in the Atlas, the monthly land warming of each AR6 region is the mean over the 20-year window centred on the year the model reaches the warming level, minus the 2001 to 2025 mean (historical plus SSP2-4.5). Levels of 1.5, 2 and 3 °C use SSP2-4.5 where available, otherwise SSP5-8.5; 4 °C uses SSP5-8.5 or SSP3-7.0. This gives 34, 34, 33 and 20 models at 1.5, 2, 3 and 4 °C. Plants outside every land region are assigned to the nearest region. The central case uses the median warming across models.

**Feed temperature.** Monthly feed temperature follows

T_feed = α + β · mean(T\*) + γ · (T\*_m − mean(T\*)), with T\* = max(T_air + ΔT, c)

where T_air is monthly air temperature, ΔT the warming and c a lower clamp chosen by leave-one-site-out cross-validation (c = 0 °C). The fit uses 16 raw influent sites (four records that are not raw influent are kept out of the fit and used in a sensitivity run), with each site given equal weight. Result: α = 8.88, β = 0.823, γ = 0.362; weighted RMSE 2.23 °C, leave-one-site-out RMSE 2.53 °C. Parameter uncertainty comes from 2,000 bootstrap fits resampled by site (seed 3).

**Hydrogen.** A month counts as warm when feed temperature is at least 25 °C (20 and 30 °C are also reported). The yield factor rises linearly from 0 at T0 to 1 at T1. Conversion follows a first-order continuously stirred reactor, η = k·HRT / (1 + k·HRT), with k = k37 · θ^(T − 37). k37 is set so that conversion at 37 °C equals η37 at the 2-day design retention time. Hydrogen comes from captured household food waste (UNEP Food Waste Index 2024 household rates by income group) and primary sludge, and is converted at 11,126 L per kg.

**Central parameters** (`CENTRAL` in `core.py`) and **Monte Carlo ranges** (`draw()`; U uniform, Tri triangular with min, mode, max):

| Parameter | Central | Range |
|---|---|---|
| Yield factor zero point T0 (°C) | 12.5 | U(10, 15) |
| Yield factor full point T1 (°C) | 31.5 | U(28, 35) |
| Temperature coefficient θ | 1.07 | U(1.05, 1.10) |
| Conversion at 37 °C, η37 | 0.90 | U(0.85, 0.95) |
| Design retention time (d) | 2 | fixed (see note) |
| Food waste capture share | 0.25 | Tri(0.10, 0.25, 0.50) |
| Food waste total solids | 0.25 | U(0.20, 0.30) |
| Food waste VS/TS | 0.90 | U(0.85, 0.95) |
| Food waste yield (mL H2/g VS) | 55 | Tri(30, 55, 100) |
| Primary sludge, high income (g TS/person/d) | 35 | Tri(25, 35, 50) |
| Primary sludge, other income groups (g TS/person/d) | 25 | Tri(15, 25, 35) |
| Primary sludge VS/TS | 0.70 | U(0.60, 0.80) |
| Primary sludge yield (mL H2/g VS) | 3 | Tri(0.5, 3, 10) |
| Feed temperature coefficients | fitted | one bootstrap set per draw |

The Monte Carlo has 400 draws (seed 7). Each draw is paired with one randomly chosen climate model at each warming level, and changes are computed against the present climate within the same draw. Model agreement is the share of models in which a plant gains at least one warm month.

*Note on retention time.* Because k37 is set from η37 at the design retention time, a sampled retention time would cancel out. `draw()` therefore keeps the retention time fixed at 2 days. It still takes one uniform draw at that point so that the random stream, and so the published Monte Carlo runs, are reproduced exactly. The effect of a different retention time with the rate fixed at the 2-day design is tested in `sens_kinetics.py` (1 to 6 days).

**Greenhouse-gas screen** (`ghg_sdg.py`, seed 11, 400 draws). Hydrogen credit uses 9.9 to 11.9 kg CO2e per kg H2 for displaced grey hydrogen; the credit for keeping food waste out of landfill uses Tri(0.3, 0.6, 1.0) t CO2e per t. The two credits are reported separately. Each draw takes its hydrogen output from the Monte Carlo run with the same index.

**Robustness** (`ghg_sdg.py`). Central parameters with other feed temperature relations: the central fit, cool and warm bootstrap fits (5th and 95th percentiles of feed temperature at 15 °C air), a fit to all 20 sites, a linear relation 8 + 0.67 T_air, and feed equal to air.

**Process sensitivity** (`sens_kinetics.py`). Fermenter temperature −3, −2, 0, +2 and +3 °C relative to influent; the temperature penalty applied through the yield factor only, the kinetics only, or both; retention time 1, 2, 3, 4 and 6 days with the rate fixed at the 2-day design. Each case is run for the present climate and 2 °C.

**Energy screen** (`energy_screen.py`). Central case, present climate, plants grouped by mean daily hydrogen output in warm months (1 to 10, 10 to 100 and at least 100 kg/d). Hydrogen energy at 33.3 kWh/kg (lower heating value). Lean auxiliaries: mixing at 5 W/m³ of reactor volume (feed diluted to 10% total solids, 2-day retention) plus membrane upgrading at 0.296 kWh/Nm³. Full-plant auxiliaries: 59.1 kWh per tonne of food waste plus upgrading at 0.400 kWh/Nm³. These are screening bounds, not a plant design.

## Outputs

| File in `results/` | Content | Used in |
|---|---|---|
| `feed_fit.json` | feed temperature coefficients, clamp, RMSE, cross-validated RMSE, 2,000 bootstrap sets, 20-site fit | Section 3.1, ESI Section S5 and Fig. S2 |
| `validation_table.csv` | each measured record with matched air temperature and predictions | ESI Table S5 and Fig. S2 |
| `plants_windows.csv` | per plant: warm months at 20, 25 and 30 °C, hydrogen (kg/yr), daily output in warm months, mean feed temperature, model agreement, for the present and each warming level | Fig. 1, Fig. 3, input to the energy screen and country table |
| `core_results.json` | global totals for the central case, 5th/50th/95th percentiles across climate models, Monte Carlo percentiles and paired changes | Abstract, Sections 3.2 to 3.4, Fig. 2, ESI Table S2 |
| `mc_runs.csv` | all 2,000 Monte Carlo rows (400 draws × present and four warming levels) | uncertainty ranges in Sections 3.3 and 3.4, Fig. 2 |
| `robustness.csv` | totals under each feed temperature relation, present and 2 °C | Section 3.9, ESI Table S3 and Fig. S1 |
| `ghg_results.json` | greenhouse-gas screen, 5th/50th/95th percentiles and central values | Section 3.7 |
| `country_table.csv` | per country: plants, population served, population-weighted warm months, hydrogen, electricity access | Sections 3.6 and 3.8, Fig. 4, ESI Table S6 |
| `sens_kinetics.csv` | process sensitivity cases, present and 2 °C | Section 3.9, ESI Table S7 |
| `energy_screen.csv` | screening energy balance by output class | Section 3.10, Table 1 |

Headline central values (present climate, then 1.5, 2, 3 and 4 °C), threshold 25 °C:

| | Present | 1.5 °C | 2 °C | 3 °C | 4 °C |
|---|---|---|---|---|---|
| Hydrogen (kt/yr) | 22.7 | 23.7 | 24.9 | 27.3 | 29.5 |
| People served by plants with at least one warm month (million) | 1,132 | 1,231 | 1,335 | 1,501 | 1,630 |
| People served by plants warm all year (million) | 376 | 401 | 447 | 570 | 612 |
| Warm months, population-weighted mean | 3.55 | 3.85 | 4.23 | 5.01 | 5.62 |

Monte Carlo: present hydrogen 27.4 kt/yr (5th to 95th percentile 12.8 to 51.9); paired change at 2 °C +10.1% (4.6 to 17.0), at 4 °C +30.0% (21.8 to 40.9).

## Verification

The package was re-run from a clean copy with the inputs listed above, and every file in `results/` was compared with `reference_results/`. All ten files matched to a relative tolerance of 1e-9 (in practice identical), including the 56,506 plant rows and the 2,000 Monte Carlo rows. Small numerical differences may appear with other versions of NumPy or pandas, or if NASA POWER revises its climatology; the HydroWASTE md5 and the Atlas commit fix the other inputs.

`validate.py` uses the shipped site climatology in `data/validation/site_clim.json` so the fit is reproduced exactly. Set `REFETCH = True` in that script to download it again from NASA POWER.

## Software

Python 3.13.18, NumPy 2.5.3, pandas 3.0.5, SciPy 1.18.1, GeoPandas 1.2.0, Shapely 2.2.0, git.
