"""Pollution, heat and compound risk maps for every detectable plant, in the population-map layout (14 October 2026).

map_no2.png       mean tropospheric NO2 column, July 2018 to September 2026 excluding the March to June 2023 transition
                  (month-weighted mean of the pre and post rasters), one log colour scale for all plants
map_heat.png      night land surface temperature anomaly, 2019 to 2024: each 1 km cell minus the mean of the plant's
                  40 to 80 km background ring (MODIS MOD11A2)
map_compound.png  cell-level compound index on the 1 km grid of populated cells (at least 1 person per km2):
                  NO2, heat anomaly, population density and travel time to health care are converted to percentiles
                  across all populated cells within 50 km of the 17 plants; hazard is the mean of the NO2 and heat
                  percentiles, exposure the density percentile and vulnerability the travel time percentile; the
                  index is their geometric mean, as in the plant-level index
results/compound_cells_by_plant.csv  people within 50 km living in cells in the top decile of the pooled index
Usage: python analysis/12_thematic_maps.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm, Normalize
from gis_common import RES, order, read, window_distance, small_multiples

PRE_MONTHS, POST_MONTHS = 56, 39  # July 2018 to February 2023; July 2023 to September 2026
dist = window_distance()
ring = (dist >= 40000) & (dist <= 80000)
within50 = dist <= 50000
chg = pd.read_csv(RES / 'plant_change.csv')
pct_no2 = chg[chg.outcome == 'no2_excess'].set_index('id')['pct_change']

def no2_mean(pid):
    pre, ext = read(pid, 'no2_pre'); post, _ = read(pid, 'no2_post')
    return (pre * PRE_MONTHS + post * POST_MONTHS) / (PRE_MONTHS + POST_MONTHS), ext

def heat_anomaly(pid):
    lst, ext = read(pid, 'lst_night')
    return lst - np.nanmean(lst[ring]), ext

def density_1km(pid):
    pop, _ = read(pid, 'pop_500m')
    pop = np.nan_to_num(pop).reshape(180, 2, 180, 2).sum(axis=(1, 3))  # people per 1 km cell = people per km2
    return pop

# pooled percentiles over populated cells within 50 km of every plant
layers = {}
for pid in order.id:
    no2, ext = no2_mean(pid); heat, _ = heat_anomaly(pid); acc, _ = read(pid, 'access')
    layers[pid] = {'no2': no2, 'heat': heat, 'dens': density_1km(pid), 'access': acc, 'ext': ext}
def pooled(key):
    v = np.concatenate([L[key][within50 & (L['dens'] >= 1)] for L in layers.values()])
    return np.sort(v[np.isfinite(v)])
ref = {k: pooled(k) for k in ('no2', 'heat', 'dens', 'access')}
def pct(key, a):
    return np.searchsorted(ref[key], a, side='right') / len(ref[key])

for pid, L in layers.items():
    hazard = (pct('no2', L['no2']) + pct('heat', L['heat'])) / 2
    idx = (hazard * pct('dens', L['dens']) * pct('access', L['access'])).clip(min=1e-6) ** (1 / 3)
    idx[(L['dens'] < 1) | ~np.isfinite(L['no2']) | ~np.isfinite(L['heat']) | ~np.isfinite(L['access'])] = np.nan
    L['compound'] = idx
cut = np.nanpercentile(np.concatenate([L['compound'][within50] for L in layers.values()]), 90)
rows = []
for pid, L in layers.items():
    top = within50 & (L['compound'] >= cut)
    rows.append({'id': pid, 'people_top_decile_cells_50km': float(L['dens'][top].sum()),
                 'people_50km': float(L['dens'][within50].sum()), 'top_decile_threshold': float(cut)})
cells = pd.DataFrame(rows).set_index('id')
cells['share_top_decile'] = cells.people_top_decile_cells_50km / cells.people_50km
cells.round(4).to_csv(RES / 'compound_cells_by_plant.csv')
print(cells.join(order.set_index('id')['name']).sort_values('people_top_decile_cells_50km', ascending=False).round(3).to_string())

# NO2
NO2_CMAP = plt.get_cmap('YlOrRd').copy()
small_multiples('map_no2.png', no2_mean, NO2_CMAP, LogNorm(10, 250),
                'Tropospheric NO$_2$ column, 2018 to 2026 (µmol m$^{-2}$)',
                lambda r: f"Excess {r.no2_excess:.0f} µmol m$^{{-2}}$ | {pct_no2[r.id]:+.0f}% after 2023", extend='both',
                cbar_ticks=[10, 20, 50, 100, 200], cbar_ticklabels=['10', '20', '50', '100', '200'])
# heat
small_multiples('map_heat.png', heat_anomaly, plt.get_cmap('RdBu_r'), TwoSlopeNorm(0, -4, 4),
                'Night land surface temperature anomaly, 2019 to 2024 (°C)',
                lambda r: f"{r.hot_days_popw_25km:.0f} days > 35 °C a year | night {r.lst_night_excess:+.1f} °C".replace('-0.0', '0.0').replace('+0.0', '0.0'), extend='both')
# compound
CMP = plt.get_cmap('YlOrRd').copy()
small_multiples('map_compound.png', lambda pid: (layers[pid]['compound'], layers[pid]['ext']), CMP, Normalize(0.2, 0.9),
                'Cell-level compound index (hazard × exposure × vulnerability)',
                lambda r: f"Index {r.compound_index:.2f} | {cells.at[r.id, 'people_top_decile_cells_50km'] / 1000:,.1f}k in top-decile cells",
                extend='both')
print('thematic maps written')
