"""Summarise the status recorded at each selected visit among women."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import living, baseline, last_period, STATES

d = living()
women = d[d.sex.eq('F')]
entry = baseline(women)
groups = [('Baseline', 'Overall', entry)]
groups += [('Baseline', name, entry[entry.state.eq(s)]) for s, name in STATES.items()]
groups += [
    ('Last visit in period', period, sub) for period, sub in last_period(women).groupby('period')
]
rows = []
for visit, group, sub in groups:
    status = sub.postmenopausal
    rows.append(
        dict(
            Visit=visit,
            Group=group,
            Women=len(sub),
            Yes=int(status.sum()),
            Known=status.count(),
            Missing=status.isna().sum(),
            Percent=round(100 * status.mean(), 1),
        )
    )
(ROOT / "R3_10_menopausa/Tables").mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(ROOT / 'R3_10_menopausa/Tables/menopausa.csv', index=False)
