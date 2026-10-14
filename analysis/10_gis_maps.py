"""GIS maps for the coal frontier paper (14 October 2026).

Figure 1  regional map: geoBoundaries countries and provinces, plants coloured by compound risk, Africa locator,
          north arrow and scale bar (Lambert azimuthal equal-area projection centred on southern Africa)
Figure 2  analysis zones for every detectable plant in one layout: population density (WorldPop 2020, 500 m),
          10, 25 and 50 km buffers, 40 to 80 km background ring, boundaries; identical 180 x 180 km extent, local UTM
S1 Fig    atlas, one page per plant in a fixed layout: population, NO2 change after mid 2023, night land surface
          temperature, with the same zones and a summary of the plant's indicators
Usage: python analysis/10_gis_maps.py
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm, Normalize
from matplotlib.backends.backend_pdf import PdfPages
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from gis_common import (FIG, RAS, NAVY, HW, risk, full, order, short, CNAME, adm0, adm1, north_arrow, scale_bar, save,
                        window_layers, read, draw_zones, small_multiples)

# Regional maps (map_region, map_compound_ranked) are drawn by 13_regional_maps.py.

POP_NORM = LogNorm(1, 5000); POP_CMAP = plt.get_cmap('magma_r').copy(); POP_CMAP.set_under('#f7f7f7')
def pop_density(pid):
    a, ext = read(pid, 'pop_500m'); dens = a / 0.25  # people per km2 in 500 m cells
    return np.where(dens > 0, dens, np.nan), ext

# ---------------- Figure: zones and population for every plant ----------------
small_multiples('map_population.png', pop_density, POP_CMAP, POP_NORM,
                'Population density, 2020 (people per km$^2$)', lambda r: f"{r.pop_25km / 1000:,.0f}k within 25 km", extend='both')

# ---------------- S1 Fig: one page per plant ----------------
(FIG / 'S1_Fig_plant_atlas.pdf').unlink(missing_ok=True)
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
        save(fig, FIG / 'atlas' / f"{int(r['rank']):02d}_{short(r['name']).replace(' ', '_')}.png", dpi=220)
        plt.close(fig)
print('maps written')
