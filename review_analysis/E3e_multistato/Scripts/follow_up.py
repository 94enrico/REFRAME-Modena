"""Observation intervals and covariate imputation in the two multistate cohorts."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import retained

rows = []
for cohort, methods in [('Historical', False), ('Methods', True)]:
    panel = retained(methods)
    multi = panel[panel.groupby('patient_id').state.transform('size').ge(2)]
    day_panel = retained(methods, calendar_days=True)
    added = day_panel.loc[day_panel.index.difference(panel.index)]
    duration = multi.groupby('patient_id').hiv_duration_yrs.transform('first')
    gaps = multi.groupby('patient_id').age_days.diff().dropna()
    q1, median, q3 = gaps.quantile([0.25, 0.50, 0.75])
    rows.append(
        dict(
            Cohort=cohort,
            People=multi.patient_id.nunique(),
            Observations=len(multi),
            Deaths=multi.death.sum(),
            Intervals=len(gaps),
            Median_days=median,
            Q1_days=q1,
            Q3_days=q3,
            Extra_living_visits_using_calendar_days=int((~added.death).sum()),
            Physical_activity_imputed=multi.phys_act_bin.isna().sum(),
            HIV_duration_imputed=duration.isna().sum(),
        )
    )
print(pd.DataFrame(rows).to_string(index=False))
