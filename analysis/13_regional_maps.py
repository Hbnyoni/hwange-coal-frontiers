"""Regional maps (14 October 2026).

map_region            the shifting coal frontier: every screened plant coloured by the change in NO2 excess after mid
                      2023 and sized by capacity, with Africa locator, north arrow and scale bar
map_compound_ranked   compound risk of all 17 detectable stations on one map, numbered by rank and sized by the
                      population within 25 km, with a Highveld zoom and ranked bars of the index and its hazard,
                      exposure and vulnerability scores
Both use geoBoundaries ADM0 and ADM1 in a Lambert azimuthal equal-area projection (25 E, 22 S); 600 dpi PNG and PDF.
Usage: python analysis/13_regional_maps.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from gis_common import RES, FIG, HW, risk, full, short, adm0, adm1, north_arrow, scale_bar, save

LAEA = ccrs.LambertAzimuthalEqualArea(central_longitude=25, central_latitude=-22)
PC = ccrs.PlateCarree()
EXTENT = [13.5, 36.5, -31, -12.5]
HV = [27.55, 30.45, -27.45, -25.55]  # Highveld zoom window
a0 = adm0.to_crs(LAEA.proj4_init); a1 = adm1.to_crs(LAEA.proj4_init)
SOURCE = ('Boundaries: geoBoundaries gbOpen ADM0 and ADM1 (commit 9469f09); ocean, lakes and locator: Natural Earth. '
          'Lambert azimuthal equal-area projection (25°E, 22°S).')
COUNTRIES = {'ZIMBABWE': (30.2, -19.6), 'ZAMBIA': (27.5, -14.2), 'BOTSWANA': (23.8, -22.0), 'SOUTH AFRICA': (22.8, -30.0),
             'NAMIBIA': (17.0, -21.5), 'MOZAMBIQUE': (34.6, -18.6), 'ANGOLA': (17.5, -13.5)}

def base(ax, extent, labels=True, grid=True):
    ax.set_extent(extent, crs=PC)
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'ocean', '50m'), facecolor='#dcebf5', edgecolor='none', zorder=0)
    a0.plot(ax=ax, facecolor='#f4f2ee', edgecolor='none', zorder=1, transform=LAEA)
    a1.boundary.plot(ax=ax, color='#c4bfb5', lw=0.3, zorder=2, transform=LAEA)
    a0.boundary.plot(ax=ax, color='#5b5b5b', lw=0.7, zorder=3, transform=LAEA)
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', '50m'), facecolor='#b9d7ea', edgecolor='none', zorder=2)
    if labels:
        for c, (x, y) in COUNTRIES.items():
            ax.text(x, y, c, fontsize=7, color='0.45', ha='center', transform=PC, zorder=4)
    if grid:
        gl = ax.gridlines(draw_labels=True, lw=0.3, color='0.78', x_inline=False, y_inline=False)
        gl.top_labels = gl.right_labels = False; gl.xlabel_style = gl.ylabel_style = {'size': 7}

def locator(fig, rect):
    ins = fig.add_axes(rect, projection=ccrs.LambertAzimuthalEqualArea(central_longitude=20, central_latitude=0))
    ins.set_extent([-19, 52, -36, 38], crs=PC)
    ins.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '110m'), facecolor='#d9d5cc', edgecolor='none')
    ins.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', '110m'), facecolor='none',
                    edgecolor='white', lw=0.3)
    x0, x1, y0, y1 = EXTENT
    ins.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], color='#b2182b', lw=1, transform=PC)
    for sp in ins.spines.values(): sp.set_edgecolor('0.5'); sp.set_linewidth(0.5)

def map_scale(ax, km=400):
    x0, x1, y0, y1 = ax.get_extent(LAEA)
    scale_bar(ax, x0 + 0.05 * (x1 - x0), y0 + 0.05 * (y1 - y0), km * 1000, f'{km} km', height=(y1 - y0) * 0.008)

other = full[~full.id.isin(risk.id)]
rk = risk.set_index('id')

# ---------------- the shifting coal frontier ----------------
chg = pd.read_csv(RES / 'plant_change.csv').query("outcome == 'no2_excess'").set_index('id')['pct_change'].rename('chg')
det = rk.join(chg)
fig = plt.figure(figsize=(7.1, 5.6))
ax = fig.add_axes([0.07, 0.08, 0.70, 0.89], projection=LAEA); base(ax, EXTENT)
xs = ax.scatter(other.longitude, other.latitude, s=26, marker='x', color='0.3', lw=1.1, zorder=5, transform=PC,
                label='Not detectable or < 300 MW')
sc = ax.scatter(det.longitude, det.latitude, s=np.sqrt(det.capacitymw) * 3, c=det.chg, cmap='RdBu_r',
                norm=TwoSlopeNorm(0, -60, 80), edgecolor='k', lw=0.5, zorder=6, transform=PC)
for pid, off in {'P12': (9, 3), 'P00': (-62, -16), 'P24': (9, 0), 'P16': (9, 2), 'P14': (-58, 9), 'P23': (-66, 2)}.items():
    r = det.loc[pid]
    ax.annotate(f"{short(r['name'])} ({r.chg:+.0f}%)", (r.longitude, r.latitude), xycoords=PC._as_mpl_transform(ax),
                xytext=off, textcoords='offset points', fontsize=7.5, zorder=7, fontweight='bold' if pid == HW else 'normal')
mb = full[full.id == 'P13'].iloc[0]
ax.annotate('Maamba (Phase II from 2026)', (mb.longitude, mb.latitude), xycoords=PC._as_mpl_transform(ax), xytext=(6, 4),
            textcoords='offset points', fontsize=7, zorder=7, color='0.25')
map_scale(ax)
na = fig.add_axes([0.80, 0.17, 0.10, 0.10]); na.axis('off'); north_arrow(na, 0.5, 0.35, size=0.45)
leg = [Line2D([], [], marker='o', ls='', mfc='white', mec='k', mew=0.5, ms=np.sqrt(np.sqrt(mw) * 3), label=f'{mw:,} MW')
       for mw in (300, 1500, 4000)]
ax.legend(handles=[xs] + leg, loc='lower right', fontsize=7, frameon=True, framealpha=0.95, labelspacing=1.1, borderpad=0.8)
cax = fig.add_axes([0.81, 0.34, 0.025, 0.32]); cb = fig.colorbar(sc, cax=cax, extend='both')
cb.set_label('Change in NO$_2$ excess after mid 2023 (%)', fontsize=8); cb.ax.tick_params(labelsize=7)
locator(fig, [0.79, 0.72, 0.19, 0.25])
fig.text(0.07, 0.01, SOURCE, fontsize=6, color='0.35')
save(fig, FIG / 'map_region.png', dpi=600, bbox_inches='tight', pad_inches=0.05)
save(fig, FIG / 'map_region.pdf', bbox_inches='tight', pad_inches=0.05); plt.close(fig)

# ---------------- compound risk of every station, ranked ----------------
CMAP = plt.get_cmap('YlOrRd'); NORM = Normalize(0.1, 0.45)
size = lambda pop: 18 + np.sqrt(pop) * 0.55
def numbered(ax, rows, scale=1.0):
    ax.scatter(rows.longitude, rows.latitude, s=size(rows.pop_25km) * scale, c=rows.compound_index, cmap=CMAP, norm=NORM,
               edgecolor='k', lw=0.6, zorder=6, transform=PC)

fig = plt.figure(figsize=(7.1, 5.0))
ax = fig.add_axes([0.035, 0.15, 0.56, 0.83], projection=LAEA); base(ax, EXTENT, grid=False)
ax.scatter(other.longitude, other.latitude, s=20, marker='x', color='0.35', lw=1, zorder=5, transform=PC)
numbered(ax, rk)
hv = rk.longitude.between(HV[0], HV[1]) & rk.latitude.between(HV[2], HV[3])
for pid, r in rk[~hv].iterrows():
    off = {'P24': (-11, 9), 'P01': (9, -9)}.get(pid)
    if off is None:
        ax.text(r.longitude, r.latitude, f"{int(r['rank'])}", transform=PC, ha='center', va='center', fontsize=7,
                fontweight='bold', color='white' if r.compound_index > 0.3 else 'k', zorder=8)
    else:
        ax.annotate(f"{int(r['rank'])}", (r.longitude, r.latitude), xycoords=PC._as_mpl_transform(ax), xytext=off,
                    textcoords='offset points', fontsize=7, fontweight='bold', ha='center', va='center', zorder=8)
ax.annotate('Hwange', (rk.loc[HW].longitude, rk.loc[HW].latitude), xycoords=PC._as_mpl_transform(ax), xytext=(12, 6),
            textcoords='offset points', fontsize=8, fontweight='bold', color='#b2182b', zorder=8)
ax.add_patch(Rectangle((HV[0], HV[2]), HV[1] - HV[0], HV[3] - HV[2], transform=PC, facecolor='none', edgecolor='k', lw=0.8, zorder=7))
ax.text(HV[1] + 0.25, HV[3], 'a', transform=PC, fontsize=8, fontweight='bold', zorder=7)
map_scale(ax)
north_arrow(ax, 0.07, 0.16, size=0.05)
z = ax.inset_axes([0.015, 0.52, 0.44, 0.40], projection=LAEA); base(z, HV, labels=False, grid=False)
z.scatter(rk[hv].longitude, rk[hv].latitude, s=70, c=rk[hv].compound_index, cmap=CMAP, norm=NORM, edgecolor='k', lw=0.6,
          zorder=6, transform=PC)  # equal sizes in the zoom; ranks set beside the symbols with leader lines
ZOOM_OFFSETS = {'P07': (-15, -9), 'P08': (17, -3), 'P03': (2, 14), 'P04': (-6, -16), 'P05': (-15, -4), 'P14': (-8, 12), 'P02': (13, 9)}
for pid, r in rk[hv].iterrows():
    z.annotate(f"{int(r['rank'])}", (r.longitude, r.latitude), xycoords=PC._as_mpl_transform(z), xytext=ZOOM_OFFSETS.get(pid, (9, 6)),
               textcoords='offset points', ha='center', va='center', fontsize=6.5, fontweight='bold', zorder=8,
               arrowprops=dict(arrowstyle='-', lw=0.4, color='0.3', shrinkA=0, shrinkB=3))
z.text(0.03, 0.91, 'a  Highveld', transform=z.transAxes, fontsize=7, fontweight='bold', zorder=9)
for sp in z.spines.values(): sp.set_edgecolor('k'); sp.set_linewidth(0.8)
handles = [plt.scatter([], [], s=size(p), facecolor='white', edgecolor='k', lw=0.6, label=f'{p / 1000:,.0f}k')
           for p in (50000, 200000, 700000)]
ax.legend(handles=handles, loc='lower right', title='People within 25 km', title_fontsize=6.5, fontsize=6.5, frameon=True,
          framealpha=0.95, labelspacing=1.3, borderpad=0.9)
b = fig.add_axes([0.715, 0.15, 0.27, 0.80])
o = risk.sort_values('rank'); y = np.arange(len(o))[::-1]
b.barh(y, o.compound_index, color=CMAP(NORM(o.compound_index)), edgecolor='k', lw=0.4, height=0.72)
for col, mk, lab in (('hazard_score', '^', 'Hazard'), ('exposure_score', 'o', 'Exposure'), ('vulnerability_score', 's', 'Vulnerability')):
    b.scatter(o[col], y, marker=mk, s=13, facecolor='white' if mk == 'o' else '0.15', edgecolor='0.15', lw=0.6, zorder=3, label=lab)
b.set_yticks(y); b.set_yticklabels([f"{int(r['rank'])}. {short(r['name'])}" for _, r in o.iterrows()], fontsize=7)
for t, (_, r) in zip(b.get_yticklabels(), o.iterrows()):
    if r.id == HW: t.set_color('#b2182b'); t.set_fontweight('bold')
b.set_xlim(0, 1.0); b.set_xlabel('Score (0 to 1)', fontsize=7.5); b.tick_params(axis='x', labelsize=7)
b.grid(axis='x', color='0.88', lw=0.5); b.set_axisbelow(True)
for sp in ('top', 'right'): b.spines[sp].set_visible(False)
b.legend(loc='lower right', fontsize=6.3, frameon=True, framealpha=0.95, handletextpad=0.3, borderpad=0.5)
b.set_title('b  Index (bars) and domain scores', fontsize=8, loc='left', pad=4, fontweight='bold')
cax = fig.add_axes([0.07, 0.10, 0.22, 0.016])
cb = fig.colorbar(plt.cm.ScalarMappable(norm=NORM, cmap=CMAP), cax=cax, orientation='horizontal', extend='both')
cb.set_label('Compound risk index', fontsize=7); cb.ax.xaxis.set_label_position('top'); cb.ax.tick_params(labelsize=6.5)
fig.text(0.035, 0.008, 'Numbers: compound risk rank. Symbol size: people within 25 km (main map; equal sizes in the zoom).\n'
         + SOURCE, fontsize=5.5, color='0.35')
save(fig, FIG / 'map_compound_ranked.png', dpi=600)
save(fig, FIG / 'map_compound_ranked.pdf'); plt.close(fig)
print('regional maps written')
