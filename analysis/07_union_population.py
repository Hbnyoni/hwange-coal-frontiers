"""Population within 10, 25 and 50 km of any detectable utility-scale plant, counting each person once (14 October 2026).
Usage: python analysis/07_union_population.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json
import ee
import pandas as pd
ROOT = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)
r = pd.read_csv(ROOT / 'results/africa_coal_risk_index.csv')
pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('year', 2020)).mosaic()
age = ee.ImageCollection('WorldPop/GP/100m/pop_age_sex').filter(ee.Filter.eq('year', 2020)).mosaic()
img = pop.rename('pop').addBands(age.select(['M_0', 'F_0', 'M_1', 'F_1']).reduce(ee.Reducer.sum()).rename('under5')).addBands(
      age.select(['M_65', 'F_65', 'M_70', 'F_70', 'M_75', 'F_75', 'M_80', 'F_80']).reduce(ee.Reducer.sum()).rename('over65'))
out = {}
for km in (10, 25, 50):
    u = ee.FeatureCollection([ee.Feature(ee.Geometry.Point([x.longitude, x.latitude]).buffer(km * 1000)) for x in r.itertuples()]).union(100).geometry()
    out[f'{km}km'] = img.reduceRegion(ee.Reducer.sum(), u, 100, maxPixels=1e10).getInfo()
json.dump(out, open(ROOT / 'results/union_population.json', 'w'), indent=2)
print(json.dumps(out, indent=1))
