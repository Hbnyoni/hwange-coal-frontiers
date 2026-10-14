"""Monthly NO2 and SO2 near every screened coal plant and in its background ring (14 October 2026).

Same plants and geometry as 02_africa_coal_screen.py: 10 km disc against a 40 to 80 km ring with
15 km around every plant masked. July 2018 (NO2) and December 2018 (SO2) to September 2026.
Usage: python analysis/03_monthly_plant_series.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse
import ee
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)

plants = pd.read_csv(OUT / 'africa_coal_screen.csv')[['id', 'country', 'name', 'latitude', 'longitude']]
pts = {r.id: ee.Geometry.Point([r.longitude, r.latitude]) for r in plants.itertuples()}
BOX = ee.Geometry.Rectangle([5, -35, 60, 20])
near_plant = ee.Image.constant(1).clip(ee.FeatureCollection([ee.Feature(g.buffer(15000)) for g in pts.values()])).mask().Not()
fc = ee.FeatureCollection([ee.Feature(g.buffer(10000), {'id': k, 'kind': 'near'}) for k, g in pts.items()] +
                          [ee.Feature(g.buffer(80000, 500).difference(g.buffer(40000, 500), 500), {'id': k, 'kind': 'ring'})
                           for k, g in pts.items()])
NO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').filterBounds(BOX).select('tropospheric_NO2_column_number_density')
SO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_SO2').filterBounds(BOX).select('SO2_column_number_density')

def month_image(d):
    d = ee.Date(d); e = d.advance(1, 'month')
    def mean_or_masked(ic, name):
        sub = ic.filterDate(d, e)
        return ee.Image(ee.Algorithms.If(sub.size().gt(0), sub.mean(), ee.Image.constant(0).updateMask(0))).rename(name)
    img = mean_or_masked(NO2, 'no2').addBands(mean_or_masked(SO2, 'so2')).multiply(1e6)
    # ring features read the masked bands, near features the unmasked ones
    img = img.addBands(img.updateMask(near_plant).rename(['no2_bg', 'so2_bg']))
    return img.reduceRegions(fc, ee.Reducer.mean(), 1113.2).map(lambda f: f.set('month', d.format('YYYY-MM')))

out = OUT / 'africa_coal_monthly.csv'
done = pd.read_csv(out) if out.exists() else pd.DataFrame()
months = pd.date_range('2018-07-01', '2026-09-01', freq='MS')
for q in range(0, len(months), 6):
    chunk = [m for m in months[q:q + 6] if done.empty or m.strftime('%Y-%m') not in set(done.month)]
    if not chunk:
        continue
    got = ee.FeatureCollection(ee.List([m.strftime('%Y-%m-%d') for m in chunk]).map(month_image)).flatten().getInfo()
    rows = []
    for f in got['features']:
        pr = f['properties']; ring = pr['kind'] == 'ring'
        rows.append({'month': pr['month'], 'id': pr['id'], 'kind': pr['kind'],
                     'no2': pr.get('no2_bg' if ring else 'no2'), 'so2': pr.get('so2_bg' if ring else 'so2')})
    done = pd.concat([done, pd.DataFrame(rows)], ignore_index=True)
    done.to_csv(out, index=False)
    print('done through', chunk[-1].strftime('%Y-%m'), flush=True)
print('EXIT complete', len(done))
