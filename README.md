# Shifting coal frontiers: compound air pollution and heat risk around sub-Saharan Africa's coal power plants

Code and derived data for *Shifting coal frontiers: compound air pollution and heat risk around sub-Saharan Africa's coal power plants and the 2023 Hwange expansion in Zimbabwe* (Nyoni, Moyo, Chinwadzimba, Sibanda, Tigere; manuscript under review).

The study screens 25 coal power plants in sub-Saharan Africa with open satellite and population data. For the 17 utility-scale plants with a detectable NO2 signal it combines hazard (TROPOMI NO2 and SO2 excess, GHAP PM2.5, ERA5-Land hot days and MODIS land surface heat), exposure (WorldPop population within 10, 25 and 50 km) and vulnerability (children, older adults and travel time to health care) into a compound risk index, and evaluates the 2023 commissioning of Hwange Units 7 and 8 in Zimbabwe with an interrupted time series and a placebo-in-space difference in differences. No survey, clinical or individual-level data are used.

## Contents

| Folder | What it holds |
|---|---|
| `analysis/` | Python scripts, run in numerical order (see below) |
| `results/` | Derived datasets: plant screen, annual and monthly NO2 and SO2 near-zone and background columns, heat, population, PM2.5 and access indicators, compound risk index and statistical results (CSV and JSON) |
| `data/` | World Bank crude death rates used in the analysis and geoBoundaries provenance with SHA-256 hashes |

Raster extracts (about 200 MB) and geoBoundaries files are not stored here; `analysis/08` and `analysis/09` reacquire them.

## Reproducing the analysis

Requirements: Python 3.13 with `earthengine-api`, `pandas`, `numpy`, `statsmodels`, `geopandas`, `rasterio`, `pyproj`, `shapely`, `matplotlib` and `cartopy`. Earth Engine steps need your own registered Cloud project; no credentials are distributed.

1. `python analysis/02_africa_coal_screen.py --project YOUR_PROJECT`: plant screen (annual NO2, SO2, heat, population)
2. `python analysis/03_monthly_plant_series.py --project YOUR_PROJECT`: monthly near-zone and background columns
3. `python analysis/04_risk_index.py --project YOUR_PROJECT`: PM2.5, health care access, compound risk index
4. `python analysis/05_statistics.py`: interrupted time series, difference in differences, plant change
5. `python analysis/06_figures.py --project YOUR_PROJECT`: statistical plots of the results
6. `python analysis/07_union_population.py --project YOUR_PROJECT`: population counted once across buffers
7. `python analysis/08_fetch_boundaries.py`: geoBoundaries at commit 9469f09
8. `python analysis/09_plant_rasters.py --project YOUR_PROJECT`: per-plant raster windows
9. `python analysis/10_gis_maps.py`: maps of the plants, analysis zones and per-plant surfaces

`01_feasibility.py` is the initial Hwange feasibility check and is kept for transparency.

## Data sources and licences

Inputs keep their own licences; please cite the original sources: Sentinel-5P TROPOMI NO2 and SO2 (Copernicus), ERA5-Land (Copernicus Climate Change Service), MODIS MOD11A2 and FIRMS (NASA), WorldPop, the Malaria Atlas Project travel time to health care, GHAP PM2.5 (Wei et al. 2023), the Global Power Plant Database v1.3.0 (WRI), geoBoundaries gbOpen (Runfola et al. 2020), World Bank World Development Indicators and Natural Earth.

## Licence

Code is released under the MIT licence (`LICENSE`). Derived datasets are released under CC BY 4.0 (`LICENSE-DATA`).

## Contact

Handsome Bongani Nyoni (nyonih@staff.msu.ac.zw) and Perez Livias Moyo (perezlivias1997@gmail.com).
