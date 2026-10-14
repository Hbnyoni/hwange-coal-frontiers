"""GIS maps for the coal frontier paper (14 October 2026).

Figure 1  regional map: geoBoundaries countries and provinces, plants coloured by compound risk, Africa locator,
          north arrow and scale bar (Lambert azimuthal equal-area projection centred on southern Africa)
Figure 2  analysis zones for every detectable plant in one layout: population density (WorldPop 2020, 500 m),
          10, 25 and 50 km buffers, 40 to 80 km background ring, boundaries; identical 180 x 180 km extent, local UTM
S1 Fig    atlas, one page per plant in a fixed layout: population, NO2 change after mid 2023, night land surface
          temperature, with the same zones and a summary of the plant's indicators
Usage: python analysis/10_gis_maps.py
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm, Normalize
from matplotlib.patches import Circle, Polygon as MPolygon, Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from shapely.geometry import box, Point
import cartopy.crs as ccrs
import cartopy.feature as cfeature

ROOT = Path(__file__).resolve().parents[1]
RES, FIG, RAS, BND = ROOT / 'results', ROOT / 'figures', ROOT / 'data/plant_rasters', ROOT / 'data/boundaries'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6})
NAVY = '#1f2d3d'; HW = 'P12'
risk = pd.read_csv(RES / 'africa_coal_risk_index.csv')
full = pd.read_csv(RES / 'africa_coal_screen_full.csv')
from pyproj import Transformer
def _window(r):  # same grid as 09_plant_rasters.py
    zone = int((r.longitude + 180) // 6) + 1; epsg = 32700 + zone if r.latitude < 0 else 32600 + zone
    x, y = Transformer.from_crs(4326, epsg, always_xy=True).transform(r.longitude, r.latitude)
    return {'epsg': epsg, 'x': round(x / 1000) * 1000, 'y': round(y / 1000) * 1000, 'half_width_m': 90000}
windows = {r.id: _window(r) for r in risk.itertuples()}
short = lambda n: n.replace(' Coal Power Plant Zimbabwe', '').replace(' Power Station Unit 1', '').replace('Morpule', 'Morupule')
CNAME = {'ZAF': 'South Africa', 'ZWE': 'Zimbabwe', 'BWA': 'Botswana', 'ZMB': 'Zambia', 'NAM': 'Namibia'}

adm0 = pd.concat([gpd.read_file(f) for f in sorted(BND.glob('*-ADM0_simplified.geojson'))], ignore_index=True)
adm1 = pd.concat([gpd.read_file(f) for f in sorted(BND.glob('*-ADM1_simplified.geojson'))], ignore_index=True)
adm0 = gpd.GeoDataFrame(adm0, crs=4326); adm1 = gpd.GeoDataFrame(adm1, crs=4326)

def north_arrow(ax, x, y, size=0.06):
    """Compass-style arrow in axes fraction coordinates, drawn natively."""
    t = ax.transAxes
    ax.add_patch(MPolygon([[x, y + size], [x - size * 0.4, y - size * 0.6], [x, y - size * 0.3]], closed=True, transform=t,
                          facecolor=NAVY, edgecolor=NAVY, lw=0.6, zorder=20))
    ax.add_patch(MPolygon([[x, y + size], [x + size * 0.4, y - size * 0.6], [x, y - size * 0.3]], closed=True, transform=t,
                          facecolor='white', edgecolor=NAVY, lw=0.6, zorder=20))
    ax.text(x, y + size * 1.15, 'N', transform=t, ha='center', va='bottom', fontsize=8, fontweight='bold', color=NAVY, zorder=20)

def scale_bar(ax, x0, y0, length_m, label, height=None, segments=4):
    """Alternating scale bar in data (metre) coordinates."""
    height = height or length_m * 0.04
    step = length_m / segments
    for k in range(segments):
        ax.add_patch(Rectangle((x0 + k * step, y0), step, height, facecolor=NAVY if k % 2 == 0 else 'white',
                               edgecolor=NAVY, lw=0.5, zorder=20))
    for v, s in ((0, '0'), (length_m, label)):
        ax.text(x0 + v, y0 - height * 0.8, s, ha='center', va='top', fontsize=7, color=NAVY, zorder=20)

# ---------------- Figure 1: regional map ----------------
LAEA = ccrs.LambertAzimuthalEqualArea(central_longitude=25, central_latitude=-22)
fig = plt.figure(figsize=(7.1, 5.6))
ax = fig.add_axes([0.07, 0.07, 0.70, 0.90], projection=LAEA)
ax.set_extent([13.5, 36.5, -31, -12.5], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m'), facecolor='#dcebf5', edgecolor='none', zorder=0)
a0 = adm0.to_crs(LAEA.proj4_init); a1 = adm1.to_crs(LAEA.proj4_init)
a0.plot(ax=ax, facecolor='#f4f2ee', edgecolor='none', zorder=1, transform=LAEA)
a1.boundary.plot(ax=ax, color='#b9b4ab', lw=0.3, zorder=2, transform=LAEA)
a0.boundary.plot(ax=ax, color='#5b5b5b', lw=0.7, zorder=3, transform=LAEA)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '50m'), facecolor='#b9d7ea', edgecolor='none', zorder=2)
other = full[~full.id.isin(risk.id)]
xs = ax.scatter(other.longitude, other.latitude, s=26, marker='x', color='0.25', lw=1.1, zorder=5, transform=ccrs.PlateCarree(),
           label='Not detectable or < 300 MW')
sc = ax.scatter(risk.longitude, risk.latitude, s=np.sqrt(risk.capacitymw) * 3, c=risk.compound_index, cmap='YlOrRd',
                vmin=0, vmax=risk.compound_index.max(), edgecolor='k', lw=0.5, zorder=6, transform=ccrs.PlateCarree())
labels = {'P12': (8, 4), 'P00': (-46, -10), 'P24': (9, 0), 'P16': (8, 2)}
for _, r in risk[risk.id.isin(labels)].iterrows():
    ax.annotate(short(r['name']), (r.longitude, r.latitude), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax),
                xytext=labels[r.id], textcoords='offset points', fontsize=8, zorder=7,
                fontweight='bold' if r.id == HW else 'normal')
ax.text(30.6, -25.3, 'Highveld\ncluster', fontsize=7.5, transform=ccrs.PlateCarree(), zorder=7)
mb = full[full.id == 'P13'].iloc[0]
ax.annotate('Maamba', (mb.longitude, mb.latitude), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax), xytext=(6, 4),
            textcoords='offset points', fontsize=8, zorder=7)
for c, (x, y) in {'ZIMBABWE': (30.2, -19.6), 'ZAMBIA': (27.5, -14.2), 'BOTSWANA': (23.8, -22.0), 'SOUTH AFRICA': (22.8, -30.0),
                  'NAMIBIA': (17.0, -21.5), 'MOZAMBIQUE': (34.6, -18.6), 'ANGOLA': (17.5, -13.5)}.items():
    ax.text(x, y, c, fontsize=7, color='0.4', ha='center', transform=ccrs.PlateCarree(), zorder=4)
gl = ax.gridlines(draw_labels=True, lw=0.3, color='0.75', x_inline=False, y_inline=False)
gl.top_labels = gl.right_labels = False; gl.xlabel_style = gl.ylabel_style = {'size': 7}
x0, x1, y0, y1 = ax.get_extent(LAEA)
scale_bar(ax, x0 + 0.05 * (x1 - x0), y0 + 0.05 * (y1 - y0), 400000, '400 km', height=(y1 - y0) * 0.008)
na = fig.add_axes([0.80, 0.20, 0.10, 0.10]); na.axis('off'); north_arrow(na, 0.5, 0.35, size=0.45)
from matplotlib.lines import Line2D
leg = [Line2D([], [], marker='o', ls='', mfc='white', mec='k', mew=0.5, ms=np.sqrt(np.sqrt(mw) * 3), label=f'{mw:,} MW')
       for mw in (300, 1500, 4000)]
ax.legend(handles=[xs] + leg, loc='lower right', fontsize=7, frameon=True, framealpha=0.95, labelspacing=1.1, borderpad=0.8)
cax = fig.add_axes([0.81, 0.36, 0.025, 0.30]); cb = fig.colorbar(sc, cax=cax); cb.set_label('Compound risk index', fontsize=8)
cb.ax.tick_params(labelsize=7)
ins = fig.add_axes([0.79, 0.72, 0.19, 0.25], projection=ccrs.LambertAzimuthalEqualArea(central_longitude=20, central_latitude=0))
ins.set_extent([-19, 52, -36, 38], crs=ccrs.PlateCarree())
ins.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m'), facecolor='#d9d5cc', edgecolor='none')
ins.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '110m'), facecolor='none', edgecolor='white', lw=0.3)
ins.plot([13.5, 36.5, 36.5, 13.5, 13.5], [-31, -31, -12.5, -12.5, -31], color='#b2182b', lw=1, transform=ccrs.PlateCarree())
for sp in ins.spines.values(): sp.set_edgecolor('0.5'); sp.set_linewidth(0.5)
fig.text(0.07, 0.008, 'Boundaries: geoBoundaries gbOpen ADM0 and ADM1 (commit 9469f09); ocean, lakes and locator: Natural Earth. '
         '\nProjection: Lambert azimuthal equal area (25°E, 22°S).', fontsize=6, color='0.35')
fig.savefig(FIG / 'Figure1_map.png', dpi=300, bbox_inches='tight', pad_inches=0.05); plt.close(fig)

# ---------------- per-plant helpers ----------------
def window_layers(pid):
    w = windows[pid]; epsg = w['epsg']; x, y, H = w['x'], w['y'], w['half_width_m']
    win = box(x - H, y - H, x + H, y + H)
    b0 = adm0.to_crs(epsg).clip(win); b1 = adm1.to_crs(epsg).clip(win)
    return w, b0, b1

def read(pid, name):
    with rasterio.open(RAS / pid / f'{name}.tif') as src:
        a = src.read(1).astype(float); b = src.bounds
    a[a < -9000] = np.nan  # computePixels fill value
    return a, [b.left, b.right, b.bottom, b.top]

def draw_zones(ax, w, b0, b1, ring=True):
    x, y = w['x'], w['y']
    b1.boundary.plot(ax=ax, color='0.55', lw=0.35, zorder=4)
    b0.boundary.plot(ax=ax, color='0.15', lw=0.8, zorder=5)
    if ring:
        ax.add_patch(Circle((x, y), 80000, facecolor='none', edgecolor='#2166ac', lw=0.7, ls=(0, (3, 2)), zorder=6))
        ax.add_patch(Circle((x, y), 40000, facecolor='none', edgecolor='#2166ac', lw=0.7, ls=(0, (3, 2)), zorder=6))
    for rad, lw in ((10000, 1.1), (25000, 0.9), (50000, 0.8)):
        ax.add_patch(Circle((x, y), rad, facecolor='none', edgecolor='white', lw=lw + 0.8, zorder=7))
        ax.add_patch(Circle((x, y), rad, facecolor='none', edgecolor='k', lw=lw, zorder=8))
    ax.plot(x, y, marker='*', ms=9, mfc='#ffd92f', mec='k', mew=0.6, zorder=9)
    H = w['half_width_m']; ax.set_xlim(x - H, x + H); ax.set_ylim(y - H, y + H)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal')

POP_NORM = LogNorm(1, 5000); POP_CMAP = plt.get_cmap('magma_r').copy(); POP_CMAP.set_under('#f7f7f7')
def pop_density(pid):
    a, ext = read(pid, 'pop_500m'); dens = a / 0.25  # people per km2 in 500 m cells
    return np.where(dens > 0, dens, np.nan), ext

order = risk.sort_values('rank')
order = order[[(RAS / i / 'lst_night.tif').exists() and (RAS / i / 'pop_500m.tif').exists() for i in order.id]]  # only downloaded windows
# ---------------- Figure 2: zones and population for every plant ----------------
ncol, nrow = 4, 5
fig, axs = plt.subplots(nrow, ncol, figsize=(7.1, 9.4))
plt.subplots_adjust(left=0.01, right=0.99, top=0.975, bottom=0.01, wspace=0.04, hspace=0.16)
for ax, (_, r) in zip(axs.flat, order.iterrows()):
    w, b0, b1 = window_layers(r.id)
    dens, ext = pop_density(r.id)
    im = ax.imshow(dens, extent=ext, origin='upper', cmap=POP_CMAP, norm=POP_NORM, interpolation='nearest', zorder=1)
    ax.set_facecolor('#f7f7f7')
    draw_zones(ax, w, b0, b1)
    ax.set_title(f"{int(r['rank'])}. {short(r['name'])} ({r.country}, {r.capacitymw:,.0f} MW)", fontsize=7.5, pad=2,
                 fontweight='bold' if r.id == HW else 'normal', color='#b2182b' if r.id == HW else 'k')
    ax.text(0.03, 0.03, f"{r.pop_25km / 1000:,.0f}k within 25 km", transform=ax.transAxes, fontsize=6.5,
            bbox=dict(facecolor='white', edgecolor='none', alpha=0.85, pad=1.2), zorder=10)
for ax in axs.flat[len(order):]: ax.axis('off')
leg_ax, cb_ax, sb_ax = axs.flat[len(order)], axs.flat[len(order) + 1], axs.flat[len(order) + 2]
# legend panel
from matplotlib.lines import Line2D
handles = [Line2D([], [], marker='*', ls='', mfc='#ffd92f', mec='k', ms=9, label='Power plant'),
           Line2D([], [], color='k', lw=1, label='10, 25 and 50 km buffers'),
           Line2D([], [], color='#2166ac', lw=0.8, ls=(0, (3, 2)), label='40 to 80 km background ring'),
           Line2D([], [], color='0.15', lw=0.8, label='National boundary'),
           Line2D([], [], color='0.55', lw=0.4, label='Province boundary')]
leg_ax.legend(handles=handles, loc='center', fontsize=6.5, frameon=False)
cax = cb_ax.inset_axes([0.1, 0.45, 0.8, 0.12]); cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend='both')
cb.set_label('Population density, 2020 (people per km$^2$)', fontsize=6.5); cb.ax.tick_params(labelsize=6)
w0 = windows[order.iloc[0].id]
sb_ax.set_xlim(0, 180000); sb_ax.set_ylim(0, 180000); sb_ax.set_aspect('equal')
scale_bar(sb_ax, 20000, 70000, 100000, '100 km', height=6000)
north_arrow(sb_ax, 0.85, 0.55, size=0.14)
sb_ax.text(0.5, 0.15, 'All panels: 180 × 180 km,\nlocal UTM projection', transform=sb_ax.transAxes, ha='center', fontsize=6.5)
fig.savefig(FIG / 'Figure2_plant_zones.png', dpi=300); plt.close(fig)

# ---------------- S1 Fig: one page per plant ----------------
with PdfPages(FIG / 'S1_Fig_plant_atlas.pdf') as pdf:
    for _, r in order.iterrows():
        w, b0, b1 = window_layers(r.id)
        fig = plt.figure(figsize=(8.27, 5.4))
        fig.text(0.02, 0.955, f"{short(r['name'])}, {CNAME.get(r.country, r.country)} ({r.capacitymw:,.0f} MW): compound risk rank "
                 f"{int(r['rank'])} of {len(risk)}", fontsize=11, fontweight='bold')
        dens, ext = pop_density(r.id)
        pre, ext2 = read(r.id, 'no2_pre'); post, _ = read(r.id, 'no2_post'); lst, ext3 = read(r.id, 'lst_night')
        dno2 = post - pre
        panels = [('a  Population density, 2020', dens, ext, POP_CMAP, POP_NORM, 'people per km$^2$'),
                  ('b  NO$_2$ change, after minus before mid 2023', dno2, ext2, 'RdBu_r',
                   TwoSlopeNorm(0, -max(5, np.nanpercentile(abs(dno2), 99)), max(5, np.nanpercentile(abs(dno2), 99))), 'µmol m$^{-2}$'),
                  ('c  Night land surface temperature, 2019 to 2024', lst, ext3, 'inferno',
                   Normalize(np.nanpercentile(lst, 2), np.nanpercentile(lst, 98)), '°C')]
        for i, (title, arr, e, cmap, norm, unit) in enumerate(panels):
            ax = fig.add_axes([0.02 + i * 0.325, 0.37, 0.30, 0.52])
            im = ax.imshow(arr, extent=e, origin='upper', cmap=cmap, norm=norm, interpolation='nearest', zorder=1)
            ax.set_facecolor('#f7f7f7'); draw_zones(ax, w, b0, b1, ring=(i > 0))
            ax.set_title(title, fontsize=8, loc='left', pad=3)
            cax = fig.add_axes([0.04 + i * 0.325, 0.325, 0.26, 0.016]); cb = fig.colorbar(im, cax=cax, orientation='horizontal')
            cb.set_label(unit, fontsize=7); cb.ax.tick_params(labelsize=6.5)
            if i == 0:
                scale_bar(ax, w['x'] - 82000, w['y'] - 84000, 50000, '50 km', height=2500)
                north_arrow(ax, 0.92, 0.84, size=0.07)
        rows = [('NO$_2$ excess, 2019 to 2025', f"{r.no2_excess:.1f} µmol m$^{{-2}}$"), ('SO$_2$ excess, 2019 to 2025', f"{r.so2_excess:.0f} µmol m$^{{-2}}$"),
                ('PM$_{2.5}$, population weighted, 2017 to 2022', f"{r.pm25_popw_25km:.1f} µg m$^{{-3}}$"),
                ('Days above 35 °C a year, 2015 to 2024', f"{r.hot_days_popw_25km:.0f}"),
                ('Night land surface heat excess', f"{r.lst_night_excess:+.1f} °C"),
                ('Population within 10 / 25 / 50 km', f"{r.pop_10km:,.0f} / {r.pop_25km:,.0f} / {r.pop_50km:,.0f}"),
                ('Children under 5 / adults 65+ within 25 km', f"{r.under5:,.0f} / {r.over65:,.0f}"),
                ('Travel time to health care', f"{r.health_travel_min_popw_25km:.0f} min"),
                ('Hazard / exposure / vulnerability scores', f"{r.hazard_score:.2f} / {r.exposure_score:.2f} / {r.vulnerability_score:.2f}")]
        for j, (k, v) in enumerate(rows):
            col, row = divmod(j, 5)
            fig.text(0.03 + col * 0.49, 0.215 - row * 0.036, k, fontsize=7.5, color='0.3')
            fig.text(0.31 + col * 0.49, 0.215 - row * 0.036, v, fontsize=7.5, fontweight='bold')
        fig.text(0.02, 0.012, 'Zones: buffers at 10, 25 and 50 km (black); background ring 40 to 80 km (blue dashed, panels b and c). '
                 'Data: WorldPop, TROPOMI, MODIS, ERA5-Land, geoBoundaries. Local UTM projection.', fontsize=6, color='0.35')
        pdf.savefig(fig)
        (FIG / 'atlas').mkdir(exist_ok=True)
        fig.savefig(FIG / 'atlas' / f"{int(r['rank']):02d}_{short(r['name']).replace(' ', '_')}.png", dpi=220)
        plt.close(fig)
print('maps written')
