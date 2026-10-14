"""Download per-plant raster windows for the GIS atlas (14 October 2026).

For each detectable plant, a 180 x 180 km window centred on the plant in its local UTM zone:
  pop_100m.tif    WorldPop 2020 unconstrained population on a 100 m UTM grid (people per cell)
  pop_500m.tif    the same summed to 500 m cells
  no2_pre.tif     TROPOMI tropospheric NO2, mean July 2018 to February 2023, 1 km (umol m-2)
  no2_post.tif    TROPOMI tropospheric NO2, mean July 2023 to September 2026, 1 km (umol m-2)
  lst_night.tif   MODIS MOD11A2 night land surface temperature, mean 2019 to 2024, 1 km (deg C)
  access.tif      Malaria Atlas Project travel time to the nearest health facility, 2019, 1 km (minutes)
Usage: python analysis/09_plant_rasters.py --project YOUR_EARTH_ENGINE_PROJECT
"""
from pathlib import Path
import argparse, json
from concurrent.futures import ThreadPoolExecutor
import ee
import numpy as np
import rasterio
import pandas as pd
from pyproj import Transformer

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/plant_rasters'; OUT.mkdir(parents=True, exist_ok=True)
p = argparse.ArgumentParser(); p.add_argument('--project', required=True); args = p.parse_args()
ee.Initialize(project=args.project)
HALF = 90000

risk = pd.read_csv(ROOT / 'results/africa_coal_risk_index.csv')
pop = ee.ImageCollection('WorldPop/GP/100m/pop').filter(ee.Filter.eq('year', 2020)).mosaic()
no2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2').select('tropospheric_NO2_column_number_density')
lst = ee.ImageCollection('MODIS/061/MOD11A2').filterDate('2019-01-01', '2025-01-01').select('LST_Night_1km')
layers = {
    'no2_pre': (no2.filterDate('2018-07-01', '2023-03-01').mean().multiply(1e6), 1000),
    'no2_post': (no2.filterDate('2023-07-01', '2026-10-01').mean().multiply(1e6), 1000),
    'lst_night': (lst.mean().multiply(0.02).subtract(273.15), 1000),
    'access': (ee.Image('projects/malariaatlasproject/assets/accessibility/accessibility_to_healthcare/2019'), 1000),
}
def fetch(r):
    zone = int((r.longitude + 180) // 6) + 1; epsg = 32700 + zone if r.latitude < 0 else 32600 + zone
    x, y = Transformer.from_crs(4326, epsg, always_xy=True).transform(r.longitude, r.latitude)
    x, y = round(x / 1000) * 1000, round(y / 1000) * 1000  # snap so 500 m and 1 km grids align
    d = OUT / r.id; d.mkdir(exist_ok=True)
    todo = {'pop_100m': (pop, 100)}  # summed to 500 m locally below
    todo.update(layers)
    for name, (img, scale) in todo.items():
        f = d / f'{name}.tif'
        if f.exists():
            continue
        n = 2 * HALF // scale
        data = ee.data.computePixels({
            'expression': img.rename('v').toFloat().unmask(-9999),
            'fileFormat': 'GEO_TIFF',
            'grid': {'dimensions': {'width': n, 'height': n}, 'crsCode': f'EPSG:{epsg}',
                     'affineTransform': {'scaleX': scale, 'shearX': 0, 'translateX': x - HALF,
                                         'shearY': 0, 'scaleY': -scale, 'translateY': y + HALF}}})
        f.write_bytes(data)
    if not (d / 'pop_500m.tif').exists():
        with rasterio.open(d / 'pop_100m.tif') as src:
            a = src.read(1); a = np.where(a < -9000, 0.0, a).clip(min=0)
            agg = a.reshape(a.shape[0] // 5, 5, a.shape[1] // 5, 5).sum(axis=(1, 3))
            prof = src.profile; prof.update(width=agg.shape[1], height=agg.shape[0], transform=src.transform * src.transform.scale(5, 5))
            with rasterio.open(d / 'pop_500m.tif', 'w', **prof) as dst: dst.write(agg.astype('float32'), 1)
    print(r.id, r.name, epsg, flush=True)
    return {'id': r.id, 'name': r.name, 'epsg': epsg, 'x': x, 'y': y, 'half_width_m': HALF}

with ThreadPoolExecutor(2) as pool:  # two at a time stays within the noncommercial concurrency limit
    meta = list(pool.map(fetch, risk.itertuples()))
json.dump(meta, open(OUT / 'windows.json', 'w'), indent=1)
