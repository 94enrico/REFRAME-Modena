"""Count distinct periods with an eligible, classifiable living visit."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import living

rows = []
for rule, returns in [
    ('Stop at first gap >3 years', False),
    ('Allow returns after >3 years', True),
]:
    periods = living(returns).groupby('patient_id').period.nunique()
    for n in (1, 2, 3):
        count = int(periods.eq(n).sum())
        rows.append(
            dict(
                Follow_up=rule,
                Periods=n,
                Participants=count,
                Total=len(periods),
                Percent=round(100 * count / len(periods), 1),
            )
        )
(ROOT / "R2_1_R3_4_tre_periodi/Tables").mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(Path(__file__).resolve().parents[1] / 'Tables/presenze.csv', index=False)
