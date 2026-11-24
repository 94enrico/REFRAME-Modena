"""Baseline characteristics and period prevalence in the same people under two follow-up rules."""

from pathlib import Path
from argparse import ArgumentParser
import sys
import numpy as np
import pandas as pd
from scipy import stats
from threadpoolctl import threadpool_limits

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'E1_prevalenza/Scripts'))
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from standardizzazione import (
    fit_standardisation,
    prevalence_table,
    contrast_table,
)
from coorti import STATES, living, all_three, baseline, random_period

FIELDS = [
    ('Age (years)', 'access_age', 1),
    ('Women', 'female', None),
    ('Postmenopausal at the visit (women)', 'postmenopausal', None),
    ('BMI (kg/m²)', 'BMI', 1),
    ('Waist (cm)', 'waist_circumference__cm', 1),
    ('Waist-to-height ratio', 'WHtR', 2),
    ('VAT (cm²)', 'VAT_dexa', 1),
    ('HIV duration (years)', 'hiv_duration_yrs', 1),
    ('CD4 count (cells/µL)', 'CD4_absolute__cells/uL', 0),
    ('CD4 nadir (cells/µL)', 'Nadir_CD4__cells/uL', 0),
    ('CD4/CD8 ratio', 'CD4/CD8__ratio', 2),
    ('Suppressed HIV viral load', 'vl_suppressed', None),
    ('INSTI', 'INSTI', None),
    ('TAF', 'TAF', None),
    ('Protease inhibitor', 'IP', None),
    ('Physical activity', 'phys_act_bin', None),
    ('Diabetes', 'comorbidity_DM', None),
    ('Hypertension', 'comorbidity_HTN', None),
    ('Dyslipidaemia', 'comorbidity_Dislipidemia', None),
]


def baseline_pvalue(samples, binary):
    if not binary:
        return stats.kruskal(*samples).pvalue
    yes = np.array([int(group.sum()) for group in samples])
    totals = np.array([len(group) for group in samples])
    cells = np.column_stack([yes, totals - yes])
    test = stats.chi2_contingency(cells, correction=False)
    if test.expected_freq.min() >= 5:
        return test.pvalue
    # Conditional Monte Carlo test with fixed group sizes and number of positives.
    rng = np.random.default_rng(52019)
    draws = rng.multivariate_hypergeometric(totals, int(yes.sum()), size=20000)
    simulated = np.stack([draws, totals - draws], axis=-1)
    scores = ((simulated - test.expected_freq) ** 2 / test.expected_freq).sum(axis=(1, 2))
    return (1 + np.count_nonzero(scores >= test.statistic - 1e-12)) / (len(draws) + 1)


def baseline_table(d):
    entry = baseline(d).copy()
    entry['female'] = entry.sex.eq('F')
    groups = {'Overall': entry, **{name: entry[entry.state.eq(s)] for s, name in STATES.items()}}
    rows = [
        dict(Panel='A. Baseline', Measure='Participants', **{k: len(v) for k, v in groups.items()})
    ]
    for label, field, decimals in FIELDS:
        row = dict(Panel='A. Baseline', Measure=label)
        for name, group in groups.items():
            obs = group[field].dropna()
            if decimals is None:
                row[name] = (
                    f'{int(obs.sum())}/{len(obs)} ({100 * obs.mean():.1f}%)' if len(obs) else '0/0'
                )
            else:
                q1, med, q3 = obs.quantile([0.25, 0.5, 0.75])
                row[name] = f'{med:.{decimals}f} ({q1:.{decimals}f}, {q3:.{decimals}f})'
        observed = [entry.loc[entry.state.eq(st), field].dropna() for st in STATES]
        p = baseline_pvalue(observed, binary=decimals is None)
        row['P value'] = '<0.001' if p < 0.001 else f'{p:.3f}'
        rows.append(row)
    return rows


def output(d, model, target):
    rows = baseline_table(d)
    p = prevalence_table(model)
    delta = contrast_table(model)
    for method in ['Crude', 'Standardised']:
        for period, block in p.groupby('Period', sort=False):
            row = dict(Panel='B. ' + method, Measure=period, Overall=int(block.N.iloc[0]))
            for _, r in block.iterrows():
                count = f'{r.n}/{r.N}, ' if method == 'Crude' else ''
                row[r.State] = (
                    count
                    + f"{r[method + ' %']:.1f}% "
                    + f"({r[method + ' lower 95%']:.1f}, {r[method + ' upper 95%']:.1f})"
                )
            rows.append(row)
        row = dict(Panel='B. ' + method, Measure='Last minus first period (pp)')
        for r in delta[delta.Method.eq(method)].itertuples():
            row[r.State] = f'{r.Difference_pp:+.1f} ({r.Lower95_pp:+.1f}, {r.Upper95_pp:+.1f})'
        rows.append(row)
    target.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(target, index=False)


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--rifitta', action='store_true')
    parser.add_argument('--coorte', choices=['entrambe', 'continua', 'rientri'], default='entrambe')
    args = parser.parse_args()
    folder = Path(__file__).resolve().parents[1]
    for tag, returns, seed in [('continua', False, 20261140), ('rientri', True, 20261144)]:
        if args.coorte not in ['entrambe', tag]:
            continue
        d = all_three(living(returns))
        n = d.patient_id.nunique()
        model = folder / f'Modelli/{tag}.npz'
        if args.rifitta or not model.exists():
            with threadpool_limits(1):
                fit_standardisation(random_period(d), model, seed)
        output(d, model, folder / f'Tables/{tag}_{n}.csv')
