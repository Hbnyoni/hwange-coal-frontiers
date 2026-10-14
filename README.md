# Shifting coal frontiers: compound air pollution and heat risk around sub-Saharan Africa's coal power plants

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23367004.svg)](https://doi.org/10.5281/zenodo.23367004)

Code and derived data for *Shifting coal frontiers: compound air pollution and heat risk around coal power plants in Zimbabwe, Botswana and South Africa, and the 2023 Hwange expansion*

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
9. `python analysis/10_gis_maps.py`: population maps and the per-plant atlas
10. `python analysis/12_thematic_maps.py`: NO2, heat and cell-level compound risk maps
11. `python analysis/13_regional_maps.py`: regional frontier map and ranked compound risk map
12. `python analysis/14_charts.py`: plant change and Hwange time series charts

`gis_common.py` holds the shared cartography used by steps 9 to 12.

`01_feasibility.py` is the initial Hwange feasibility check and is kept for transparency.

## Where to get the data

Derived datasets produced by this study are in `results/`; each GitHub release is archived on Zenodo with a DOI
(see "Citation" below). All inputs are public and are reacquired by the scripts:

| Input | Access | Licence / citation |
|---|---|---|
| Sentinel-5P TROPOMI tropospheric NO2, offline L3 | Earth Engine `COPERNICUS/S5P/OFFL/L3_NO2`, https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S5P_OFFL_L3_NO2 | Copernicus open licence |
| Sentinel-5P TROPOMI SO2, offline L3 | Earth Engine `COPERNICUS/S5P/OFFL/L3_SO2`, https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S5P_OFFL_L3_SO2 | Copernicus open licence |
| ERA5-Land daily and monthly aggregates | Earth Engine `ECMWF/ERA5_LAND/DAILY_AGGR`, `ECMWF/ERA5_LAND/MONTHLY_AGGR`, https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_DAILY_AGGR | Copernicus C3S licence; Muñoz-Sabater et al. 2021 |
| MODIS land surface temperature, MOD11A2 v6.1 | Earth Engine `MODIS/061/MOD11A2`, https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD11A2 | NASA open data; Wan 2014 |
| FIRMS active fires | Earth Engine `FIRMS`, https://developers.google.com/earth-engine/datasets/catalog/FIRMS | NASA open data |
| WorldPop 2020 population and age-sex structure | Earth Engine `WorldPop/GP/100m/pop`, `WorldPop/GP/100m/pop_age_sex`, https://www.worldpop.org | CC BY 4.0; Tatem 2017 |
| Travel time to health care, 2019 | Earth Engine `projects/malariaatlasproject/assets/accessibility/accessibility_to_healthcare/2019`, https://malariaatlas.org | CC BY 4.0; Weiss et al. 2020 |
| GHAP PM2.5 (1 km, 2017 to 2022) | Earth Engine Community Catalog, https://gee-community-catalog.org/projects/ghap/ | Wei et al. 2023 |
| Global Power Plant Database v1.3.0 | Earth Engine `WRI/GPPD/power_plants`, https://datasets.wri.org/dataset/globalpowerplantdatabase | CC BY 4.0 |
| geoBoundaries gbOpen ADM0 and ADM1, commit 9469f09 | https://github.com/wmgeolab/geoBoundaries/tree/9469f09 (fetched by `analysis/08_fetch_boundaries.py`) | Runfola et al. 2020 |
| World Bank crude death rates | https://data.worldbank.org/indicator/SP.DYN.CDRT.IN (copy in `data/`) | CC BY 4.0 |
| Natural Earth physical layers | https://www.naturalearthdata.com (downloaded by cartopy) | Public domain |

Earth Engine access needs a free registered account and Cloud project, passed to the scripts with `--project`.

## Licence

Code is released under the MIT licence (`LICENSE`). Derived datasets are released under CC BY 4.0 (`LICENSE-DATA`).

## Citation

Please cite the archived release and the associated article:

Nyoni, H. B., Moyo, P. L., Chinwadzimba, Z., Sibanda, R., & Tigere, D. (2026). *Shifting coal frontiers: Compound air
pollution and heat risk around sub-Saharan African coal power plants: Code and derived data* (Version 1.0.0)
[Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23367004

## Contact

Handsome Bongani Nyoni (nyonih@staff.msu.ac.zw)
