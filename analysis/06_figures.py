"""Manuscript figures for the African coal plant compound risk paper (14 October 2026).

Figure 1  map of the screened plants coloured by compound risk index
Figure 2  Hwange monthly NO2 and SO2 excess with commissioning dates, and monthly hot days
Figure 3  change in NO2 and SO2 excess after mid 2023 at every detectable plant
Figure 4  compound risk: hazard against vulnerability, sized by population, and domain scores
Usage: python analysis/06_figures.py --project YOUR_EARTH_ENGINE_PROJECT   (the project is used once to cache Hwange hot days)
"""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = Path(__file__).resolve().parents[1]
RES, FIG = ROOT / 'results', ROOT / 'figures'
FIG.mkdir(exist_ok=True)
p = argparse.ArgumentParser(); p.add_argument('--project'); args = p.parse_args()
plt.rcParams.update({'font.family': 'Arial', 'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False})
HW = 'P12'; HWC = '#b2182b'; GREY = '0.55'
U7, U8 = pd.Timestamp('2023-03-20'), pd.Timestamp('2023-06-09')

risk = pd.read_csv(RES / 'africa_coal_risk_index.csv')
full = pd.read_csv(RES / 'africa_coal_screen_full.csv')
mon = pd.read_csv(RES / 'africa_coal_monthly_excess.csv', parse_dates=['date'])
chg = pd.read_csv(RES / 'plant_change.csv')
short = lambda n: n.replace(' Coal Power Plant Zimbabwe', '').replace(' Power Station Unit 1', '').replace('Morpule', 'Morupule')

# Hwange monthly hot days (ERA5-Land daily Tmax > 35 C), cached.
hot_csv = RES / 'hwange_monthly_hot_days.csv'
if not hot_csv.exists():
    import ee
    ee.Initialize(project=args.project)
    r = full.set_index('id').loc[HW]
    g = ee.Geometry.Point([r.longitude, r.latitude]).buffer(10000)
    daily = ee.ImageCollection('ECMWF/ERA5_LAND/DAILY_AGGR').select('temperature_2m_max')
    def one(d):
        d = ee.Date(d); img = daily.filterDate(d, d.advance(1, 'month')).map(lambda i: i.gt(308.15)).sum()
        return ee.Feature(None, {'month': d.format('YYYY-MM'),
                                 'hot_days': img.reduceRegion(ee.Reducer.mean(), g, 5000).get('temperature_2m_max')})
    ms = [m.strftime('%Y-%m-%d') for m in pd.date_range('2018-07-01', '2026-09-01', freq='MS')]
    got = ee.FeatureCollection(ee.List(ms).map(one)).getInfo()
    pd.DataFrame([f['properties'] for f in got['features']]).to_csv(hot_csv, index=False)
hot = pd.read_csv(hot_csv); hot['date'] = pd.to_datetime(hot.month)

# Figure 1: map.
fig = plt.figure(figsize=(7, 6))
ax = fig.add_subplot(projection=ccrs.PlateCarree()); ax.set_extent([14, 36, -30, -14])
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '50m'), facecolor='#f2f0eb', edgecolor='none')
ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '50m'), facecolor='none', edgecolor='0.4', lw=0.6)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '50m'), facecolor='#cfe3f0', edgecolor='none')
other = full[~full.id.isin(risk.id)]
ax.scatter(other.longitude, other.latitude, s=28, marker='x', color='0.25', lw=1.2, zorder=4, transform=ccrs.PlateCarree(), label='Not detectable or < 300 MW')
sc = ax.scatter(risk.longitude, risk.latitude, s=np.sqrt(risk.capacitymw) * 3, c=risk.compound_index, cmap='YlOrRd',
                vmin=0, vmax=risk.compound_index.max(), edgecolor='k', lw=0.5, transform=ccrs.PlateCarree(), zorder=3)
for _, r in risk.iterrows():
    if r.id in (HW, 'P00', 'P16', 'P24', 'P13') or r['rank'] <= 3:
        off = {'P00': (-48, -4), 'P07': (8, 6)}.get(r.id, (6, 4))
        ax.annotate(short(r['name']), (r.longitude, r.latitude), xytext=off, textcoords='offset points', fontsize=8,
                    transform=ccrs.PlateCarree(), fontweight='bold' if r.id == HW else 'normal')
mb = full[full.id == 'P13'].iloc[0]
ax.annotate('Maamba', (mb.longitude, mb.latitude), xytext=(6, 4), textcoords='offset points', fontsize=8)
for c, (x, y) in {'ZIMBABWE': (29.5, -19.2), 'ZAMBIA': (26.5, -15.0), 'BOTSWANA': (23.6, -22.3), 'SOUTH AFRICA': (22.5, -29.2),
                  'NAMIBIA': (16.5, -21.0), 'MOZAMBIQUE': (33.6, -21.0)}.items():
    ax.text(x, y, c, fontsize=7, color='0.45', ha='center')
cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02); cb.set_label('Compound risk index (0 to 1)')
for mw in (300, 1500, 4000):
    ax.scatter([], [], s=np.sqrt(mw) * 3, facecolor='white', edgecolor='k', lw=0.5, label=f'{mw:,} MW')
ax.legend(loc='upper left', fontsize=7, frameon=True, framealpha=0.9)
gl = ax.gridlines(draw_labels=True, lw=0.3, color='0.8'); gl.top_labels = gl.right_labels = False
plt.close(fig)  # superseded: Figure 1 is drawn by 10_gis_maps.py with geoBoundaries

# Figure 2: Hwange time series.
h = mon[mon.id == HW].sort_values('date')
fig, axs = plt.subplots(3, 1, figsize=(7, 7), sharex=True, gridspec_kw={'height_ratios': [1, 1, 0.6]})
for ax, y, lab in [(axs[0], 'no2_excess', 'NO$_2$ excess\n(µmol m$^{-2}$)'), (axs[1], 'so2_excess', 'SO$_2$ excess\n(µmol m$^{-2}$)')]:
    ax.plot(h.date, h[y], color=GREY, lw=0.8)
    ax.plot(h.date, h[y].rolling(12, center=True, min_periods=9).mean(), color=HWC, lw=2, label='12 month running mean')
    for per, col in (('pre', '#2166ac'), ('post', HWC)):
        v = h.loc[h.period == per, y].mean(); d = h.loc[h.period == per, 'date']
        ax.hlines(v, d.min(), d.max(), color=col, ls=':', lw=1.2)
    ax.set_ylabel(lab)
axs[0].legend(frameon=False, fontsize=8, loc='upper left')
axs[2].bar(hot.date, hot.hot_days, width=25, color='#ef8a62')
axs[2].set_ylabel('Days with\nTmax > 35 °C')
for ax in axs:
    ax.axvspan(U7, pd.Timestamp('2023-07-01'), color='0.9', zorder=0)
    ax.axvline(U7, color='k', ls='--', lw=0.8); ax.axvline(U8, color='k', ls='--', lw=0.8)
axs[0].text(U7, axs[0].get_ylim()[1] * 0.97, 'Unit 7 ', fontsize=8, va='top', ha='right'); axs[0].text(U8, axs[0].get_ylim()[1] * 0.97, ' Unit 8', fontsize=8, va='top')
fig.align_ylabels(axs); fig.tight_layout()
(FIG / 'draft_hwange_timeseries.png').unlink(missing_ok=True); fig.savefig(FIG / 'draft_hwange_timeseries.png', dpi=300); plt.close(fig)

# Figure 3: change by plant.
fig, axs = plt.subplots(1, 2, figsize=(7.5, 5), sharey=True)
order = chg[chg.outcome == 'no2_excess'].sort_values('pct_change').name.tolist()
for ax, y, lab in [(axs[0], 'no2_excess', 'NO$_2$'), (axs[1], 'so2_excess', 'SO$_2$')]:
    c = chg[chg.outcome == y].set_index('name').reindex(order)
    pos = np.arange(len(c))
    pct_lo, pct_hi = 100 * c.ci_low / c.pre_mean, 100 * c.ci_high / c.pre_mean
    col = [HWC if 'Hwange' in n else '0.3' for n in c.index]
    ax.hlines(pos, pct_lo, pct_hi, color=col, lw=1.2); ax.scatter(c['pct_change'], pos, color=col, s=22, zorder=3)
    ax.axvline(0, color='0.6', lw=0.8); ax.set_xlabel(f'Change in {lab} excess after mid 2023 (%)')
axs[0].set_yticks(np.arange(len(order))); axs[0].set_yticklabels([short(n) for n in order])
fig.tight_layout(); (FIG / 'draft_plant_change.png').unlink(missing_ok=True); fig.savefig(FIG / 'draft_plant_change.png', dpi=300); plt.close(fig)

# Figure 4: compound risk.
fig, axs = plt.subplots(1, 2, figsize=(8, 4.2), gridspec_kw={'width_ratios': [1.1, 1]})
ax = axs[0]
ax.scatter(risk.hazard_score, risk.vulnerability_score, s=risk.pop_25km / 4000 + 10, c=risk.compound_index, cmap='YlOrRd',
           edgecolor='k', lw=0.5, vmin=0, vmax=risk.compound_index.max())
for _, r in risk.iterrows():
    off = {'P24': (-10, -11), 'P04': (-38, 4), 'P07': (4, -10)}.get(r.id, (4, 3))
    ax.annotate(short(r['name']), (r.hazard_score, r.vulnerability_score), xytext=off, textcoords='offset points',
                fontsize=7, color=HWC if r.id == HW else 'k', fontweight='bold' if r.id == HW else 'normal')
ax.set_xlabel('Hazard score'); ax.set_ylabel('Vulnerability score'); ax.set_title('a', loc='left', fontweight='bold')
ax = axs[1]
r = risk.sort_values('compound_index')
pos = np.arange(len(r))
for i, (col, c) in enumerate([('hazard_score', '#d6604d'), ('exposure_score', '#4393c3'), ('vulnerability_score', '#878787')]):
    ax.barh(pos + (i - 1) * 0.27, r[col], height=0.27, color=c, label=col.split('_')[0].title())
ax.scatter(r.compound_index, pos, color='k', marker='D', s=14, zorder=3, label='Compound index')
ax.set_yticks(pos); ax.set_yticklabels([short(n) for n in r['name']], fontsize=7)
ax.set_xlabel('Score (0 to 1)'); ax.legend(fontsize=7, frameon=False, loc='lower right'); ax.set_title('b', loc='left', fontweight='bold')
fig.tight_layout(); (FIG / 'draft_compound_scores.png').unlink(missing_ok=True); fig.savefig(FIG / 'draft_compound_scores.png', dpi=300); plt.close(fig)
print('figures written')
