"""Compound pollution and heat risk index for the screened coal plants (14 October 2026).

Adds to results/africa_coal_screen.csv:
  GHAP V1 PM2.5, 2017 to 2022 mean: 10 km disc, ring and population weighted within 25 km
  conditional PM2.5 mortality contrast for the near-plant excess (log-linear, RR 1.08 per 10 ug m-3,
  95% CI 1.06 to 1.09; Chen and Hoek 2020), using World Bank crude death rates (mean 2017 to 2019)
  population weighted ERA5-Land days with Tmax > 35 C within 25 km (2015 to 2024)
  travel time to health care within 25 km (Malaria Atlas Project 2019, Weiss et al. 2020)
Then builds hazard, exposure and vulnerability scores (min-max across detectable plants) and a
geometric-mean compound index, with arithmetic-mean and leave-one-indicator-out sensitivity ranks.
Usage: python analysis/04_risk_index.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json
import ee
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)

s = pd.read_csv(OUT / 'africa_coal_screen.csv').set_index('id')
pts = {k: ee.Geometry.Point([r.longitude, r.latitude]) for k, r in s.iterrows()}
near_plant = ee.Image.constant(1).clip(ee.FeatureCollection([ee.Feature(g.buffer(15000)) for g in pts.values()])).mask().Not()

def fc(kind):
    geoms = {'d10': lambda g: g.buffer(10000), 'd25': lambda g: g.buffer(25000),
             'ring': lambda g: g.buffer(80000, 500).difference(g.buffer(40000, 500), 500)}[kind]
    return ee.FeatureCollection([ee.Feature(geoms(g), {'id': k}) for k, g in pts.items()])

def table(img, kind, scale, reducer=ee.Reducer.mean()):
    if kind == 'ring':
        img = img.updateMask(near_plant)
    got = img.reduceRegions(fc(kind), reducer, scale).getInfo()['features']
    return pd.DataFrame([f['properties'] for f in got]).set_index('id')

pm = ee.ImageCollection('projects/sat-io/open-datasets/GHAP/GHAP_Y1K_PM25').mean().rename('pm')  # 2017 to 2022
pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('year', 2020)).mosaic().rename('pop')
daily = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').select('temperature_2m_max')
hot = ee.ImageCollection([daily.filterDate(f'{y}-01-01', f'{y + 1}-01-01').map(lambda i: i.gt(308.15)).sum()
                          for y in range(2015, 2025)]).mean().rename('hot')
access = ee.Image('projects/malariaatlasproject/assets/accessibility/accessibility_to_healthcare/2019').rename('access')

s['pm25_near'] = table(pm, 'd10', 1000)['mean']
s['pm25_ring'] = table(pm, 'ring', 1000)['mean']
s['pm25_excess'] = s.pm25_near - s.pm25_ring
weighted = table(pop.addBands(pm.multiply(pop).rename('pm_pop')).addBands(hot.multiply(pop).rename('hot_pop'))
                 .addBands(access.multiply(pop).rename('acc_pop')), 'd25', 100, ee.Reducer.sum())
s['pm25_popw_25km'] = weighted.pm_pop / weighted['pop']
s['hot_days_popw_25km'] = weighted.hot_pop / weighted['pop']
s['health_travel_min_popw_25km'] = weighted.acc_pop / weighted['pop']

# Conditional mortality contrast for the near-plant PM2.5 excess, applied to people within 10 km.
s['pop_10km'] = s['pop_10km']
cdr = json.load(open(ROOT / 'data/worldbank_crude_death_rate.json'))[1]
cdr = pd.DataFrame([{'country': r['countryiso3code'], 'year': int(r['date']), 'cdr': r['value']} for r in cdr])
cdr = cdr[cdr.year <= 2019].groupby('country').cdr.mean()
s['crude_death_rate'] = s.country.map(cdr)
deaths = s.pop_10km * s.crude_death_rate / 1000
for tag, rr in (('', 1.08), ('_lo', 1.06), ('_hi', 1.09)):
    af = 1 - np.exp(-np.log(rr) / 10 * s.pm25_excess.clip(lower=0))
    s['pm25_excess_deaths_10km' + tag] = deaths * af

# Detectable utility-scale plants: capacity >= 300 MW, mean annual NO2 excess > 5 umol m-2 and its t value > 3.
s['detectable'] = (s.capacitymw >= 300) & (s.no2_excess > 5) & (s.no2_excess_t > 3)
d = s[s.detectable].copy()
d['vulnerable_share'] = (d.under5 + d.over65) / d.pop_25km
ind = {'hazard': ['no2_excess', 'so2_excess', 'pm25_excess', 'hot_days_popw_25km', 'lst_night_excess'],
       'exposure': ['pop_25km'],
       'vulnerability': ['vulnerable_share', 'health_travel_min_popw_25km']}
def mm(x): return (x - x.min()) / (x.max() - x.min())
norm = {c: mm(np.log1p(d[c]) if c in ('pop_25km', 'so2_excess', 'no2_excess') else d[c]) for v in ind.values() for c in v}
def index(drop=None, geometric=True):
    dom = {k: pd.concat([norm[c] for c in v if c != drop], axis=1).mean(axis=1) for k, v in ind.items()}
    dom = pd.DataFrame(dom)
    return (dom.clip(lower=0.01).prod(axis=1) ** (1 / 3) if geometric else dom.mean(axis=1)), dom
d['compound_index'], dom = index()
d[['hazard_score', 'exposure_score', 'vulnerability_score']] = dom.values
d['rank'] = d.compound_index.rank(ascending=False).astype(int)
d['rank_arithmetic'] = index(geometric=False)[0].rank(ascending=False).astype(int)
loo = pd.DataFrame({c: index(drop=c)[0].rank(ascending=False) for v in ind.values() for c in v if len(v) > 1})
d['rank_loo_min'] = loo.min(axis=1).astype(int); d['rank_loo_max'] = loo.max(axis=1).astype(int)

s.round(4).to_csv(OUT / 'africa_coal_screen_full.csv')
d.sort_values('rank').round(4).to_csv(OUT / 'africa_coal_risk_index.csv')
print(d.sort_values('rank')[['name', 'country', 'compound_index', 'hazard_score', 'exposure_score', 'vulnerability_score',
                             'rank_arithmetic', 'rank_loo_min', 'rank_loo_max', 'pm25_excess', 'pm25_popw_25km',
                             'hot_days_popw_25km', 'health_travel_min_popw_25km', 'pm25_excess_deaths_10km']].round(2).to_string())
