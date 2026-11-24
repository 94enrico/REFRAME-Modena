"""Current use of first- and second-generation INSTIs at the baseline visit."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import living, baseline, STATES

entry = baseline(living())
if not entry.g1.isna().equals(entry.g2.isna()):
    raise ValueError('Partial generation coding requires separate denominators')
rows = []
for label, group in [('Overall', entry)] + [
    (name, entry[entry.state.eq(s)]) for s, name in STATES.items()
]:
    known = group.dropna(subset=['g1', 'g2'])
    rows.append(
        dict(
            Group=label,
            Participants=len(group),
            Known=len(known),
            Missing=len(group) - len(known),
            First_generation=int(known.g1.sum()),
            First_percent=round(100 * known.g1.mean(), 1),
            Second_generation=int(known.g2.sum()),
            Second_percent=round(100 * known.g2.mean(), 1),
        )
    )
(ROOT / "R3_14_INSTI/Tables").mkdir(parents=True, exist_ok=True)
pd.DataFrame(rows).to_csv(ROOT / 'R3_14_INSTI/Tables/INSTI.csv', index=False)
