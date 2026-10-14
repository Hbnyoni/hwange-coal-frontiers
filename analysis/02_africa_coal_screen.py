"""Screen sub-Saharan African coal plants for compound air pollution and heat risk (13 October 2026).

Plants: WRI Global Power Plant Database (WRI/GPPD/power_plants, fuel1 = Coal), plus Medupi,
which the database omits. For each plant:
  hazard     NO2 and SO2 excess (10 km disc minus a 40 to 80 km ring with 15 km around other
             plants masked), annual 2019 to 2025; MODIS land surface temperature excess (5 km vs ring);
             ERA5-Land days with Tmax > 35 C (2015 to 2024) and the 1981 to 2024 Tmax trend
  exposure   WorldPop 2020 population within 10, 25 and 50 km
  vulnerab.  children under 5 and adults 65+ within 25 km; mean travel time to health care (MAP 2019)
Usage: python analysis/02_africa_coal_screen.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse
import ee
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)

SSA = ['AGO', 'BWA', 'ZAF', 'ZWE', 'ZMB', 'MOZ', 'MWI', 'TZA', 'KEN', 'NGA', 'SEN', 'GHA', 'MDG', 'NAM',
       'SWZ', 'LSO', 'COD', 'ETH', 'CIV', 'MUS', 'NER', 'GIN', 'SDN', 'TCD', 'CMR', 'BEN', 'TGO', 'UGA']
gppd = ee.FeatureCollection('WRI/GPPD/power_plants').filter(ee.Filter.eq('fuel1', 'Coal')).filter(ee.Filter.inList('country', SSA))
plants = pd.DataFrame([f['properties'] for f in gppd.getInfo()['features']])[['country', 'name', 'capacitymw', 'latitude', 'longitude']]
plants = pd.concat([plants, pd.DataFrame([{'country': 'ZAF', 'name': 'Medupi', 'capacitymw': 4764.0,
                                           'latitude': -23.7050, 'longitude': 27.5550}])], ignore_index=True)
plants['id'] = [f'P{i:02d}' for i in range(len(plants))]

pts = {r.id: ee.Geometry.Point([r.longitude, r.latitude]) for r in plants.itertuples()}
BOX = ee.Geometry.Rectangle([5, -35, 50, 20])
all_pts = ee.FeatureCollection([ee.Feature(g) for g in pts.values()])
near_plant = ee.Image.constant(1).clip(all_pts.map(lambda f: f.buffer(15000))).mask().Not()  # 0 within 15 km of any plant
def fc(kind):
    feats = []
    for k, g in pts.items():
        geom = {'d5': g.buffer(5000), 'd10': g.buffer(10000), 'd25': g.buffer(25000), 'd50': g.buffer(50000),
                'ring': g.buffer(80000, 500).difference(g.buffer(40000, 500), 500)}[kind]
        feats.append(ee.Feature(geom, {'id': k}))
    return ee.FeatureCollection(feats)

def table(img, kind, scale, reducer=ee.Reducer.mean()):
    if kind == 'ring':
        img = img.updateMask(near_plant)  # background excludes 15 km around every plant
    print('  reducing', kind, img.bandNames().getInfo(), flush=True)
    got = img.reduceRegions(fc(kind), reducer, scale).getInfo()['features']
    return pd.DataFrame([f['properties'] for f in got]).set_index('id')

# Hazard 1: NO2 and SO2 excess by year.
NO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').filterBounds(BOX).select('tropospheric_NO2_column_number_density')
SO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2').filterBounds(BOX).select('SO2_column_number_density')
rows = []
for y in ([] if (OUT / 'africa_coal_air_annual.csv').exists() else range(2019, 2026)):
    img = NO2.filterDate(f'{y}-01-01', f'{y + 1}-01-01').mean().rename('no2').addBands(
          SO2.filterDate(f'{y}-01-01', f'{y + 1}-01-01').mean().rename('so2')).multiply(1e6)  # umol m-2
    near, ring = table(img, 'd10', 1113.2), table(img, 'ring', 1113.2)
    for k in near.index:
        rows.append({'id': k, 'year': y, 'no2_near': near.at[k, 'no2'], 'no2_ring': ring.at[k, 'no2'],
                     'so2_near': near.at[k, 'so2'], 'so2_ring': ring.at[k, 'so2']})
air = pd.DataFrame(rows) if rows else pd.read_csv(OUT / 'africa_coal_air_annual.csv').drop(columns=['no2_excess', 'so2_excess'])
air['no2_excess'] = air.no2_near - air.no2_ring; air['so2_excess'] = air.so2_near - air.so2_ring
air.to_csv(OUT / 'africa_coal_air_annual.csv', index=False)

# Hazard 2: heat.
lst = ee.ImageCollection('MODIS/061/MOD11A2').filterDate('2019-01-01', '2025-01-01')
lst_img = lst.select('LST_Day_1km').mean().multiply(0.02).subtract(273.15).rename('lst_day').addBands(
          lst.select('LST_Night_1km').mean().multiply(0.02).subtract(273.15).rename('lst_night'))
lst_near, lst_ring = table(lst_img, 'd5', 1000), table(lst_img, 'ring', 1000)
daily = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').select('temperature_2m_max')
hot_days = ee.ImageCollection([daily.filterDate(f'{y}-01-01', f'{y + 1}-01-01').map(lambda i: i.gt(308.15)).sum()
                               for y in range(2015, 2025)]).mean().rename('hot_days_35c')
monthly = ee.ImageCollection('ECMWF/ERA5_LAND/MONTHLY_AGGR').select('temperature_2m_max')
annual_t = ee.ImageCollection([monthly.filterDate(f'{y}-01-01', f'{y + 1}-01-01').mean().subtract(273.15)
                               .rename('t').addBands(ee.Image.constant(y).float().rename('year')) for y in range(1981, 2025)])
trend = annual_t.select(['year', 't']).reduce(ee.Reducer.linearFit()).select('scale').multiply(10).rename('tmax_trend_per_decade')
heat = table(hot_days.addBands(trend), 'd10', 5000)

# Exposure and vulnerability.
pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('year', 2020)).mosaic()
age = ee.ImageCollection('WorldPop/GP/100m/pop_age_sex').filter(ee.Filter.eq('year', 2020)).mosaic()
u5 = age.select(['M_0', 'F_0', 'M_1', 'F_1']).reduce(ee.Reducer.sum()).rename('under5')
o65 = age.select(['M_65', 'F_65', 'M_70', 'F_70', 'M_75', 'F_75', 'M_80', 'F_80']).reduce(ee.Reducer.sum()).rename('over65')
popn = {f'pop_{k}km': table(pop, f'd{k}', 100, ee.Reducer.sum())['sum'] for k in (10, 25, 50)}
vul = table(u5.addBands(o65), 'd25', 100, ee.Reducer.sum())
try:
    access = table(ee.Image('Oxford/MAP/accessibility_to_healthcare_2019').select('accessibility').rename('health_travel_min'), 'd25', 1000)
except ee.EEException as e:  # the catalogue asset was reported corrupted on 13 October 2026
    print('health care access skipped:', e)
    access = pd.DataFrame({'health_travel_min': np.nan}, index=plants.id)

air_mean = air.groupby('id')[['no2_near', 'no2_ring', 'no2_excess', 'so2_excess']].mean()
air_sd = air.groupby('id')['no2_excess'].std()
early = air[air.year <= 2022].groupby('id')[['no2_excess', 'so2_excess']].mean()
late = air[air.year >= 2023].groupby('id')[['no2_excess', 'so2_excess']].mean()
s = plants.set_index('id').join(air_mean)
s['no2_excess_t'] = air_mean.no2_excess / (air_sd / np.sqrt(7))
s['no2_change_2023_25_vs_2019_22_pct'] = 100 * (late.no2_excess / early.no2_excess - 1)
s['so2_change_2023_25_vs_2019_22_pct'] = 100 * (late.so2_excess / early.so2_excess - 1)
s['lst_day_excess'] = lst_near.lst_day - lst_ring.lst_day
s['lst_night_excess'] = lst_near.lst_night - lst_ring.lst_night
s = s.join(heat[['hot_days_35c', 'tmax_trend_per_decade']]).join(pd.DataFrame(popn)).join(vul[['under5', 'over65']]).join(access[['health_travel_min']])
s.round(3).to_csv(OUT / 'africa_coal_screen.csv')
print(s[['country', 'name', 'capacitymw', 'no2_excess', 'no2_excess_t', 'so2_excess', 'no2_change_2023_25_vs_2019_22_pct',
         'lst_day_excess', 'lst_night_excess', 'hot_days_35c', 'tmax_trend_per_decade', 'pop_25km', 'under5', 'over65',
         'health_travel_min']].round(2).to_string())
