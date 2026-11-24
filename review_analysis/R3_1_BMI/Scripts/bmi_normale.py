"""Classification and supporting adiposity criteria in people with BMI <25."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import living, STATES
from agreement_draw import select_one

visits = living().dropna(subset=['bmi_cat'])
sample = select_one(visits)
normal = sample.loc[sample.bmi_cat.eq(1)].copy()
vat = normal.VAT_dexa.ge(150)
anthro = normal.anthro_pair
criteria = {
    'Both criteria met': vat & anthro,
    'VAT met, anthropometry below threshold': vat & normal.anthro_known & ~anthro,
    'VAT met, anthropometry incomplete': vat & ~normal.anthro_known,
    'Anthropometry met, VAT below threshold': anthro & normal.VAT_dexa.lt(150),
    'Anthropometry met, VAT unavailable': anthro & normal.VAT_dexa.isna(),
}
rows = []
for s, label in STATES.items():
    group = normal.state.eq(s)
    row = dict(
        State=label,
        People=int(group.sum()),
        Normal_BMI=len(normal),
        Percent=round(100 * group.mean(), 1),
    )
    for name, flag in criteria.items():
        row[name] = int((group & flag).sum())
    if s != 1 and sum(row[k] for k in criteria) != row['People']:
        raise ValueError('Incomplete attribution of the adiposity criteria')
    rows.append(row)
(ROOT / "R3_1_BMI/Tables").mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(ROOT / 'R3_1_BMI/Tables/BMI_normale.csv', index=False)
