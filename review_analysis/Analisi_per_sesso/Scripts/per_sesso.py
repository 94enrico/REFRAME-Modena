"""One manuscript table for the reviewer requests on sex-specific descriptions and transitions."""

from pathlib import Path
import sys
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import PERIODS, STATES, baseline, last_period, random_period, retained

panel = retained()
living = panel[~panel.death]
entry = baseline(living)
periods = random_period(living)
last = last_period(living)
sexes = {'Women': 'F', 'Men': 'M'}
rows = []


def anthropometry(frame, panel_name, group):
    for label, field in [
        ('Waist circumference (cm)', 'waist_circumference__cm'),
        ('Waist-to-hip ratio', 'hip_ratio'),
    ]:
        row = dict(Panel=panel_name, Group=group, Measure=label)
        digits = 1 if field == 'waist_circumference__cm' else 2
        for name, sex in sexes.items():
            observed = frame.loc[frame.sex.eq(sex), field].dropna()
            q1, median, q3 = observed.quantile([0.25, 0.5, 0.75])
            row[name] = f'{median:.{digits}f} ({q1:.{digits}f}, {q3:.{digits}f}), n={len(observed)}'
        rows.append(row)


anthropometry(entry, 'A. Baseline', 'Overall')
for state, name in STATES.items():
    anthropometry(entry[entry.state.eq(state)], 'A. Baseline', name)
for period in PERIODS:
    frame = periods[periods.period.eq(period)]
    anthropometry(last[last.period.eq(period)], 'B. Anthropometry at last visit', period)
    for state, name in STATES.items():
        row = dict(Panel='C. Prevalence at random visit', Group=period, Measure=name)
        for label, sex in sexes.items():
            subset = frame[frame.sex.eq(sex)]
            n, total = subset.state.eq(state).sum(), len(subset)
            interval = binomtest(n, total).proportion_ci(method='wilson')
            row[label] = (
                f'{n}/{total}, {100 * n / total:.1f}% ({100 * interval.low:.1f}, {100 * interval.high:.1f})'
            )
        rows.append(row)

# Count pairs of consecutive observations. Death is observed at its exact date.
following = panel.groupby('patient_id').state.shift(-1)
for origin, first_label in STATES.items():
    origin_rows = panel.state.eq(origin) & following.notna()
    for dest, last_label in {**STATES, 4: 'Death'}.items():
        row = dict(Panel='D. Consecutive observations', Group=first_label, Measure=last_label)
        for label, sex in sexes.items():
            endpoints = following[origin_rows & panel.sex.eq(sex)]
            n, total = endpoints.eq(dest).sum(), len(endpoints)
            row[label] = f'{n}/{total} ({100 * n / total:.1f}%)'
        rows.append(row)
(ROOT / "Analisi_per_sesso/Tables").mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(ROOT / 'Analisi_per_sesso/Tables/per_sesso.csv', index=False)
