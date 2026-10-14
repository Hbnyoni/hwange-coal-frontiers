"""Feasibility check for the Hwange expansion study (13 October 2026).

Monthly TROPOMI NO2 and SO2 over the Hwange plant against a regional background ring,
VIIRS night lights for Hwange, Bulawayo and Harare, and FIRMS fire counts near the plant,
July 2018 to September 2026. Units 7 and 8 were synchronised from 20 March 2023.
Usage: python analysis/01_feasibility.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json
import ee
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'; OUT.mkdir(exist_ok=True)
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)

NO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').select('tropospheric_NO2_column_number_density')
SO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2').select('SO2_column_number_density')
VIIRS = ee.ImageCollection('NOAA/VIIRS/DNB/MONTHLY_V1/VCMSLCFG').select('avg_rad')
FIRE = ee.ImageCollection('FIRMS').select('T21')

# 1. Locate the plant NO2 hotspot: maximum of the 2019 to 2025 mean within 25 km of Hwange town.
search = ee.Geometry.Point([26.49, -18.37]).buffer(25000)
mean_no2 = NO2.filterDate('2019-01-01', '2026-01-01').mean()
ll = mean_no2.rename('no2').addBands(ee.Image.pixelLonLat())
peak = ll.reduceRegion(ee.Reducer.max(3), search, 1113.2, maxPixels=1e8).getInfo()
lon, lat = peak['max1'], peak['max2']  # max(3) maximises the first band (NO2) and returns the paired lon and lat
hot = ee.Geometry.Point([lon, lat])
print('hotspot', round(lat, 4), round(lon, 4), 'mean NO2', peak['max'])

regions = {
    'plant_10km': hot.buffer(10000),
    'ring_40_80km': hot.buffer(80000).difference(hot.buffer(40000)),
}
lights = {
    'hwange_town': ee.Geometry.Point([26.50, -18.36]).buffer(8000),
    'bulawayo': ee.Geometry.Point([28.58, -20.15]).buffer(15000),
    'harare': ee.Geometry.Point([31.05, -17.83]).buffer(20000),
}
months = pd.date_range('2018-07-01', '2026-09-01', freq='MS')

def monthly(ic, geoms, reducer, scale, name):
    fc = ee.FeatureCollection([ee.Feature(g, {'region': k}) for k, g in geoms.items()])
    rows = []
    for y in sorted(set(months.year)):
        ms = [m for m in months if m.year == y]
        def one(d):
            d = ee.Date(d); img = ic.filterDate(d, d.advance(1, 'month'))
            img = ee.Image(ee.Algorithms.If(img.size().gt(0), img.reduce(reducer), ee.Image.constant(0).updateMask(0)))
            return img.rename('v').reduceRegions(fc, ee.Reducer.mean(), scale).map(lambda f: f.set('month', d.format('YYYY-MM')))
        got = ee.FeatureCollection(ee.List([m.strftime('%Y-%m-%d') for m in ms]).map(one)).flatten().getInfo()
        rows += [{'month': f['properties']['month'], 'region': f['properties']['region'], name: f['properties'].get('mean')} for f in got['features']]
    return pd.DataFrame(rows)

no2 = monthly(NO2, regions, ee.Reducer.mean(), 1113.2, 'no2')
so2 = monthly(SO2, regions, ee.Reducer.mean(), 1113.2, 'so2')
ntl = monthly(VIIRS, lights, ee.Reducer.mean(), 463.8, 'avg_rad')
fire_area = hot.buffer(80000)
def fire_month(d):
    d = ee.Date(d); img = FIRE.filterDate(d, d.advance(1, 'month')).map(lambda i: i.gt(0).unmask(0)).sum()
    v = img.reduceRegion(ee.Reducer.sum(), fire_area, 1000, maxPixels=1e9).get('T21')
    return ee.Feature(None, {'month': d.format('YYYY-MM'), 'fire_pixel_days_80km': v})
got = ee.FeatureCollection(ee.List([m.strftime('%Y-%m-%d') for m in months]).map(fire_month)).getInfo()
fire = pd.DataFrame([f['properties'] for f in got['features']])

air = no2.merge(so2, on=['month', 'region'])
air.to_csv(OUT / 'feasibility_air_monthly.csv', index=False)
ntl.to_csv(OUT / 'feasibility_nightlights_monthly.csv', index=False)
fire.to_csv(OUT / 'feasibility_fires_monthly.csv', index=False)
json.dump({'hotspot_lat': lat, 'hotspot_lon': lon, 'mean_no2_2019_2025_mol_m2': peak['max']},
          open(OUT / 'feasibility_hotspot.json', 'w'), indent=2)
print('saved', len(air), len(ntl), len(fire))
