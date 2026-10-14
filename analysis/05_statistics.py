"""Before and after statistics for Hwange and the other screened plants (14 October 2026).

Monthly excess = 10 km disc minus background ring (results/africa_coal_monthly.csv).
Periods: pre = July 2018 to February 2023; transition = March to June 2023 (Unit 7 synchronised
20 March 2023, Unit 8 about 9 to 13 June 2023); post = July 2023 to September 2026.

1. Hwange interrupted time series: excess ~ post + calendar month, Newey-West (6 lags);
   sensitivity with a linear trend, with the transition included as post, and placebo break dates.
2. Difference in differences with placebo in space: log(near/ring) ~ plant + month + treated x post,
   using the 17 detectable utility-scale plants; every control plant is re-assigned as the treated
   plant in turn and the Hwange estimate is ranked against that distribution.
Outputs results/hwange_its.csv, results/did_placebo.csv and results/plant_change.csv.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results'
HWANGE = 'P12'
m = pd.read_csv(OUT / 'africa_coal_monthly.csv')
w = m.pivot_table(index=['id', 'month'], columns='kind', values=['no2', 'so2']).reset_index()
w.columns = ['id', 'month', 'no2_near', 'no2_ring', 'so2_near', 'so2_ring']
w['no2_excess'] = w.no2_near - w.no2_ring
w['so2_excess'] = w.so2_near - w.so2_ring
w['no2_logratio'] = np.log(w.no2_near / w.no2_ring)
w['date'] = pd.to_datetime(w.month)
w['cal'] = w.date.dt.month.astype(str)
w['t'] = (w.date - pd.Timestamp('2018-07-01')).dt.days / 365.25
def period(d, start='2023-03-01', end='2023-07-01'):
    return np.where(d < pd.Timestamp(start), 'pre', np.where(d < pd.Timestamp(end), 'transition', 'post'))
w['period'] = period(w.date)
w.to_csv(OUT / 'africa_coal_monthly_excess.csv', index=False)

risk = pd.read_csv(OUT / 'africa_coal_risk_index.csv')
plants = pd.read_csv(OUT / 'africa_coal_screen_full.csv').set_index('id')
detect = list(risk.id)

# 1. Hwange interrupted time series.
h = w[w.id == HWANGE].dropna(subset=['no2_excess']).copy()
rows = []
def its(df, y, label, formula_extra='', post_col='post'):
    fit = smf.ols(f'{y} ~ {post_col} + C(cal){formula_extra}', df).fit(cov_type='HAC', cov_kwds={'maxlags': 6})
    pre_mean = df.loc[df[post_col] == 0, y].mean()
    b, (lo, hi) = fit.params[post_col], fit.conf_int().loc[post_col]
    rows.append({'model': label, 'outcome': y, 'n_months': int(fit.nobs), 'pre_mean': pre_mean, 'post_effect': b,
                 'ci_low': lo, 'ci_high': hi, 'p': fit.pvalues[post_col], 'pct_change': 100 * b / pre_mean})
for y in ('no2_excess', 'so2_excess'):
    d = h[h.period != 'transition'].dropna(subset=[y]).assign(post=lambda x: (x.period == 'post').astype(int))
    its(d, y, 'primary')
    its(d, y, 'with linear trend', ' + t')
    its(h.dropna(subset=[y]).assign(post=lambda x: (x.period != 'pre').astype(int)), y, 'transition counted as post')
    for brk in ('2020-07-01', '2021-07-01', '2022-01-01'):  # placebo breaks inside the pre period
        d0 = h[(h.date < pd.Timestamp('2023-03-01'))].dropna(subset=[y]).copy()
        d0['post'] = (d0.date >= pd.Timestamp(brk)).astype(int)
        its(d0, y, f'placebo break {brk[:7]}')
pd.DataFrame(rows).round(4).to_csv(OUT / 'hwange_its.csv', index=False)
print(pd.DataFrame(rows).round(2).to_string())

# 2. Difference in differences with placebo in space (NO2 log ratio and SO2 excess).
panel = w[w.id.isin(detect) & (w.period != 'transition')].copy()
panel['post'] = (panel.period == 'post').astype(int)
res = []
for y in ('no2_logratio', 'so2_excess'):
    dd = panel.dropna(subset=[y])
    for k in detect:
        dd = dd.assign(treat_post=((dd.id == k) & (dd.post == 1)).astype(int))
        fit = smf.ols(f'{y} ~ treat_post + C(id) + C(month)', dd).fit()
        res.append({'outcome': y, 'id': k, 'name': plants.at[k, 'name'], 'estimate': fit.params['treat_post']})
did = pd.DataFrame(res)
did['rank_desc'] = did.groupby('outcome').estimate.rank(ascending=False)
did['n_plants'] = did.groupby('outcome').id.transform('count')
did['permutation_p'] = did.rank_desc / did.n_plants
did.round(4).to_csv(OUT / 'did_placebo.csv', index=False)
print(did[did.id == HWANGE].round(3).to_string())

# 3. Change in mean excess by plant, post vs pre (same calendar months).
chg = []
for k, g in w[w.id.isin(detect)].groupby('id'):
    for y in ('no2_excess', 'so2_excess'):
        c = g[g.period != 'transition'].groupby(['cal', 'period'])[y].mean().unstack()
        if {'pre', 'post'} <= set(c.columns):
            diff = (c.post - c.pre).dropna()
            se = diff.std() / np.sqrt(len(diff))
            chg.append({'id': k, 'name': plants.at[k, 'name'], 'country': plants.at[k, 'country'], 'outcome': y,
                        'pre_mean': c.pre.mean(), 'post_mean': c.post.mean(), 'change': diff.mean(),
                        'ci_low': diff.mean() - 1.96 * se, 'ci_high': diff.mean() + 1.96 * se,
                        'pct_change': 100 * diff.mean() / c.pre.mean()})
pd.DataFrame(chg).round(3).to_csv(OUT / 'plant_change.csv', index=False)
print(pd.DataFrame(chg).query("outcome == 'no2_excess'").sort_values('pct_change', ascending=False).round(1).to_string())
