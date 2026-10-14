"""Download geoBoundaries gbOpen ADM0 and ADM1 for the mapped countries at the pinned commit (14 October 2026).

Same source and commit as the Midlands paper (wmgeolab/geoBoundaries 9469f09). Files and SHA-256 hashes go to
data/boundaries and data/boundaries/provenance.json. Simplified geometries are used for ADM1 to keep files small.
"""
from pathlib import Path
import hashlib, json, urllib.request

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'data/boundaries'; OUT.mkdir(parents=True, exist_ok=True)
COMMIT = '9469f09'
ISO = ['ZAF', 'ZWE', 'ZMB', 'BWA', 'NAM', 'MOZ', 'LSO', 'SWZ', 'MWI', 'AGO', 'TZA', 'COD', 'MDG', 'MUS', 'NER']
prov = {'source': 'geoBoundaries gbOpen (Runfola et al. 2020), CC BY 4.0 / CC BY 3.0 IGO by country', 'commit': COMMIT, 'files': []}
for iso in ISO:
    for level, suffix in (('ADM0', '_simplified'), ('ADM1', '_simplified')):
        name = f'geoBoundaries-{iso}-{level}{suffix}.geojson'
        url = f'https://raw.githubusercontent.com/wmgeolab/geoBoundaries/{COMMIT}/releaseData/gbOpen/{iso}/{level}/{name}'
        dest = OUT / name
        if not dest.exists():
            data = urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent': 'boundary-fetch/1.0'}), timeout=120).read()
            if data[:200].lstrip().startswith(b'version https://git-lfs'):
                lfs = url.replace('raw.githubusercontent.com', 'media.githubusercontent.com/media')
                data = urllib.request.urlopen(urllib.request.Request(lfs, headers={'User-Agent': 'boundary-fetch/1.0'}), timeout=300).read()
                url = lfs
            dest.write_bytes(data)
        prov['files'].append({'iso': iso, 'level': level, 'file': name, 'source_url': url,
                              'sha256': hashlib.sha256(dest.read_bytes()).hexdigest(), 'bytes': dest.stat().st_size})
        print(iso, level, dest.stat().st_size)
json.dump(prov, open(OUT / 'provenance.json', 'w'), indent=1)
