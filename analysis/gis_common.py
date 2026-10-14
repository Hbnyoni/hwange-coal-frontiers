"""Shared cartography for the coal frontier maps (14 October 2026).

Plant windows (same grid as 09_plant_rasters.py), geoBoundaries layers, raster reading, analysis zones, north arrow,
scale bar and the common multi-panel layout used for every per-plant figure (one panel per detectable plant, ordered
by compound risk rank, identical 180 x 180 km extent in the local UTM projection).
"""
from pathlib import Path
import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Polygon as MPolygon, Rectangle
from pyproj import Transformer
from shapely.geometry import box

ROOT = Path(__file__).resolve().parents[1]
RES, FIG, RAS, BND = ROOT / 'results', ROOT / 'figures', ROOT / 'data/plant_rasters', ROOT / 'data/boundaries'
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6})
NAVY = '#1f2d3d'; HW = 'P12'; HALF = 90000
risk = pd.read_csv(RES / 'africa_coal_risk_index.csv')
full = pd.read_csv(RES / 'africa_coal_screen_full.csv')
order = risk.sort_values('rank')
short = lambda n: n.replace(' Coal Power Plant Zimbabwe', '').replace(' Power Station Unit 1', '').replace('Morpule', 'Morupule')
CNAME = {'ZAF': 'South Africa', 'ZWE': 'Zimbabwe', 'BWA': 'Botswana', 'ZMB': 'Zambia', 'NAM': 'Namibia'}

def _window(r):
    zone = int((r.longitude + 180) // 6) + 1; epsg = 32700 + zone if r.latitude < 0 else 32600 + zone
    x, y = Transformer.from_crs(4326, epsg, always_xy=True).transform(r.longitude, r.latitude)
    return {'epsg': epsg, 'x': round(x / 1000) * 1000, 'y': round(y / 1000) * 1000, 'half_width_m': HALF}
windows = {r.id: _window(r) for r in risk.itertuples()}

adm0 = gpd.GeoDataFrame(pd.concat([gpd.read_file(f) for f in sorted(BND.glob('*-ADM0_simplified.geojson'))], ignore_index=True), crs=4326)
adm1 = gpd.GeoDataFrame(pd.concat([gpd.read_file(f) for f in sorted(BND.glob('*-ADM1_simplified.geojson'))], ignore_index=True), crs=4326)

def save(fig, path, **kw):
    """Save a figure, deleting any old file first (OneDrive can refuse an in-place overwrite)."""
    Path(path).unlink(missing_ok=True)
    fig.savefig(path, **kw)

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

def window_layers(pid):
    w = windows[pid]; x, y, H = w['x'], w['y'], w['half_width_m']
    win = box(x - H, y - H, x + H, y + H)
    return w, adm0.to_crs(w['epsg']).clip(win), adm1.to_crs(w['epsg']).clip(win)

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
        for rad in (80000, 40000):
            ax.add_patch(Circle((x, y), rad, facecolor='none', edgecolor='#2166ac', lw=0.7, ls=(0, (3, 2)), zorder=6))
    for rad, lw in ((10000, 1.1), (25000, 0.9), (50000, 0.8)):
        ax.add_patch(Circle((x, y), rad, facecolor='none', edgecolor='white', lw=lw + 0.8, zorder=7))
        ax.add_patch(Circle((x, y), rad, facecolor='none', edgecolor='k', lw=lw, zorder=8))
    ax.plot(x, y, marker='*', ms=9, mfc='#ffd92f', mec='k', mew=0.6, zorder=9)
    H = w['half_width_m']; ax.set_xlim(x - H, x + H); ax.set_ylim(y - H, y + H)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect('equal')

def window_distance():
    """Distance (m) from the plant to every 1 km cell centre of a window (identical for all plants)."""
    n = 2 * HALF // 1000
    c = (np.arange(n) + 0.5) * 1000 - HALF
    xx, yy = np.meshgrid(c, -c)
    return np.hypot(xx, yy)

def small_multiples(fname, layer, cmap, norm, cbar_label, label, extend='neither', cbar_ticks=None, cbar_ticklabels=None,
                    background='#f0f0f0'):
    """One panel per detectable plant in compound-risk order, with shared legend, colour bar, scale and north arrow.

    layer(pid) -> (array, extent); label(row) -> text for the lower-left box of each panel.
    """
    fig, axs = plt.subplots(5, 4, figsize=(7.1, 9.4))
    plt.subplots_adjust(left=0.01, right=0.99, top=0.975, bottom=0.01, wspace=0.04, hspace=0.16)
    for ax, (_, r) in zip(axs.flat, order.iterrows()):
        w, b0, b1 = window_layers(r.id)
        arr, ext = layer(r.id)
        im = ax.imshow(arr, extent=ext, origin='upper', cmap=cmap, norm=norm, interpolation='nearest', zorder=1)
        ax.set_facecolor(background)
        draw_zones(ax, w, b0, b1)
        ax.set_title(f"{int(r['rank'])}. {short(r['name'])} ({r.country}, {r.capacitymw:,.0f} MW)", fontsize=7.5, pad=2,
                     fontweight='bold' if r.id == HW else 'normal', color='#b2182b' if r.id == HW else 'k')
        ax.text(0.03, 0.03, label(r), transform=ax.transAxes, fontsize=6.3, va='bottom',
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.88, pad=1.2), zorder=10)
    for ax in axs.flat[len(order):]: ax.axis('off')
    leg_ax, cb_ax, sb_ax = axs.flat[len(order)], axs.flat[len(order) + 1], axs.flat[len(order) + 2]
    handles = [Line2D([], [], marker='*', ls='', mfc='#ffd92f', mec='k', ms=9, label='Power plant'),
               Line2D([], [], color='k', lw=1, label='10, 25 and 50 km buffers'),
               Line2D([], [], color='#2166ac', lw=0.8, ls=(0, (3, 2)), label='40 to 80 km background ring'),
               Line2D([], [], color='0.15', lw=0.8, label='National boundary'),
               Line2D([], [], color='0.55', lw=0.4, label='Province boundary')]
    leg_ax.legend(handles=handles, loc='center', fontsize=6.5, frameon=False)
    cax = cb_ax.inset_axes([0.0, 0.45, 1.0, 0.12] if cbar_ticklabels is not None else [0.1, 0.45, 0.8, 0.12])
    cb = fig.colorbar(im, cax=cax, orientation='horizontal', extend=extend)
    if cbar_ticks is not None:
        cb.set_ticks(cbar_ticks)
        if cbar_ticklabels is not None:
            cb.set_ticklabels(cbar_ticklabels); cb.ax.tick_params(length=0)
    cb.set_label(cbar_label, fontsize=6.5); cb.ax.tick_params(labelsize=5.8 if cbar_ticklabels is not None else 6)
    sb_ax.set_xlim(0, 180000); sb_ax.set_ylim(0, 180000); sb_ax.set_aspect('equal')
    scale_bar(sb_ax, 20000, 70000, 100000, '100 km', height=6000)
    north_arrow(sb_ax, 0.85, 0.55, size=0.14)
    sb_ax.text(0.5, 0.15, 'All panels: 180 × 180 km,\nlocal UTM projection', transform=sb_ax.transAxes, ha='center', fontsize=6.5)
    save(fig, FIG / fname, dpi=600)
    save(fig, FIG / fname.replace('.png', '.pdf'), dpi=300); plt.close(fig)
