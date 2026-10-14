"""Journal-style charts (14 October 2026), matching the map style; 600 dpi PNG and vector PDF.

plot_plant_change       change in NO2 and SO2 excess after mid 2023 at every detectable plant (dot and whisker)
plot_hwange_timeseries  Hwange monthly NO2 and SO2 excess with the interrupted time series estimates, and hot days
Usage: python analysis/14_charts.py
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from gis_common import RES, FIG, HW, risk, short, save

plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.linewidth': 0.6, 'xtick.major.width': 0.6,
                     'ytick.major.width': 0.6, 'xtick.major.size': 3, 'ytick.major.size': 3,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False})
INK = '#1f2d3d'; RED = '#b2182b'; BLUE = '#2166ac'; GREY = '#8c8c8c'
COUNTRY = {'ZAF': ('#2c5d8a', 'South Africa'), 'ZWE': (RED, 'Zimbabwe'), 'BWA': ('#1b9e77', 'Botswana')}
U7, U8, POST = pd.Timestamp('2023-03-20'), pd.Timestamp('2023-06-09'), pd.Timestamp('2023-07-01')

def letter(ax, s, x=-0.02, y=1.02):
    ax.text(x, y, s, transform=ax.transAxes, fontsize=10, fontweight='bold', ha='right', va='bottom')

# ---------------- change after mid 2023, every plant ----------------
chg = pd.read_csv(RES / 'plant_change.csv').merge(risk[['id', 'country']], on='id', suffixes=('', '_r'))
no2 = chg[chg.outcome == 'no2_excess'].sort_values('pct_change')
order = list(no2.id)
fig, axs = plt.subplots(1, 2, figsize=(7.1, 4.6), sharey=True, gridspec_kw={'wspace': 0.08})
for ax, (y, lab, lim) in zip(axs, (('no2_excess', 'NO$_2$', (-100, 120)), ('so2_excess', 'SO$_2$', (-140, 160)))):
    c = chg[chg.outcome == y].set_index('id').loc[order]
    pos = np.arange(len(c))
    for k in pos[::2]: ax.axhspan(k - 0.5, k + 0.5, color='#f3f3f3', zorder=0, lw=0)
    lo, hi = 100 * c.ci_low / c.pre_mean, 100 * c.ci_high / c.pre_mean
    for k, (pid, r) in enumerate(c.iterrows()):
        col = COUNTRY[r.country][0]
        ax.plot([lo.iloc[k], hi.iloc[k]], [k, k], color=col, lw=1.4 if pid == HW else 1.0, solid_capstyle='round', zorder=2)
        ax.scatter(r['pct_change'], k, s=42 if pid == HW else 24, color=col, edgecolor='white', lw=0.6, zorder=3)
        if pid == HW:
            ax.annotate(f"{r['pct_change']:+.0f}%", (hi.iloc[k], k), xytext=(4, 0), textcoords='offset points',
                        va='center', fontsize=7.5, fontweight='bold', color=RED)
    ax.axvline(0, color=INK, lw=0.7, zorder=1)
    ax.set_xlim(*lim); ax.set_ylim(-0.6, len(c) - 0.4)
    ax.set_xlabel(f'Change in {lab} excess after mid 2023 (%)')
    ax.grid(axis='x', color='#e3e3e3', lw=0.5, zorder=0); ax.set_axisbelow(True)
    ax.spines['left'].set_visible(False); ax.tick_params(axis='y', length=0)
    letter(ax, 'a' if y == 'no2_excess' else 'b')
axs[0].set_yticks(np.arange(len(order)))
axs[0].set_yticklabels([short(n) for n in no2.set_index('id').loc[order, 'name']])
for t, pid in zip(axs[0].get_yticklabels(), order):
    if pid == HW: t.set_color(RED); t.set_fontweight('bold')
handles = [Line2D([], [], marker='o', ls='-', color=col, mec='white', ms=6, lw=1, label=name) for col, name in COUNTRY.values()]
axs[1].legend(handles=handles, loc='lower right', fontsize=7, handlelength=1.6)
fig.text(0.5, -0.005, 'Points: mean change between the same calendar months before (Jul 2018 to Feb 2023) and after '
         '(Jul 2023 to Sep 2026); lines: 95% CI.', ha='center', fontsize=6.5, color='0.35')
save(fig, FIG / 'plot_plant_change.png', dpi=600, bbox_inches='tight', pad_inches=0.04)
save(fig, FIG / 'plot_plant_change.pdf', bbox_inches='tight', pad_inches=0.04); plt.close(fig)

# ---------------- Hwange time series ----------------
mon = pd.read_csv(RES / 'africa_coal_monthly_excess.csv', parse_dates=['date'])
h = mon[mon.id == HW].sort_values('date')
its = pd.read_csv(RES / 'hwange_its.csv').query("model == 'primary'").set_index('outcome')
hot = pd.read_csv(RES / 'hwange_monthly_hot_days.csv'); hot['date'] = pd.to_datetime(hot.month)
fig, axs = plt.subplots(3, 1, figsize=(7.1, 6.6), sharex=True, gridspec_kw={'height_ratios': [1, 1, 0.55], 'hspace': 0.12})
start, end = h.date.min() - pd.Timedelta(days=20), h.date.max() + pd.Timedelta(days=20)
for ax, (y, lab, unit) in zip(axs[:2], (('no2_excess', 'NO$_2$ excess', 'µmol m$^{-2}$'), ('so2_excess', 'SO$_2$ excess', 'µmol m$^{-2}$'))):
    r = its.loc[y]
    pre, post = h[h.date < U7], h[h.date >= POST]
    ax.scatter(h.date, h[y], s=9, color=GREY, alpha=0.75, lw=0, zorder=2, label='Monthly excess')
    ax.plot(h.date, h[y].rolling(12, center=True, min_periods=9).mean(), color=INK, lw=1.3, zorder=3, label='12 month running mean')
    ax.hlines(r.pre_mean, pre.date.min(), pre.date.max(), color=BLUE, lw=1.8, zorder=4, label='Pre-expansion mean')
    lvl, lo, hi = r.pre_mean + r.post_effect, r.pre_mean + r.ci_low, r.pre_mean + r.ci_high
    ax.fill_between([post.date.min(), post.date.max()], lo, hi, color=RED, alpha=0.15, lw=0, zorder=1)
    ax.hlines(lvl, post.date.min(), post.date.max(), color=RED, lw=1.8, zorder=4, label='Post-expansion level (95% CI)')
    ax.text(post.date.min() + pd.Timedelta(days=30), hi, f"+{r.post_effect:.1f} {unit} "
            f"({r.ci_low:.1f} to {r.ci_high:.1f}), +{r['pct_change']:.0f}%", fontsize=7.5, color=RED, va='bottom', fontweight='bold')
    ax.set_ylabel(f'{lab}\n({unit})')
    ax.grid(axis='y', color='#e8e8e8', lw=0.5); ax.set_axisbelow(True)
axs[0].legend(loc='upper left', fontsize=6.8, ncol=2, columnspacing=1.2, handlelength=1.6)
axs[2].bar(hot.date, hot.hot_days, width=24, color='#e8743b', lw=0)
axs[2].set_ylabel('Days with\nT$_{max}$ > 35 °C')
axs[2].grid(axis='y', color='#e8e8e8', lw=0.5); axs[2].set_axisbelow(True)
for ax, s in zip(axs, 'abc'):
    ax.axvspan(U7, POST, color='#ededed', zorder=0, lw=0)
    for d in (U7, U8): ax.axvline(d, color=INK, ls=(0, (3, 2)), lw=0.7, zorder=1)
    ax.set_xlim(start, end); letter(ax, s, x=-0.075)
top = axs[0].get_ylim()[1]
axs[0].text(U7, top, 'Unit 7 ', ha='right', va='top', fontsize=7); axs[0].text(U8, top, ' Unit 8', ha='left', va='top', fontsize=7)
fig.align_ylabels(axs)
save(fig, FIG / 'plot_hwange_timeseries.png', dpi=600, bbox_inches='tight', pad_inches=0.04)
save(fig, FIG / 'plot_hwange_timeseries.pdf', bbox_inches='tight', pad_inches=0.04); plt.close(fig)
print('charts written')
