"""Compare clinical obesity prevalence for each pair of calendar periods."""

from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'E1_prevalenza/Scripts'))
from standardizzazione import contrast_table

rows = []
for label, name in [('Random visit', 'casuale'), ('Last visit', 'ultima')]:
    model = ROOT / f'E1_prevalenza/Modelli/{name}.npz'
    for earlier, later in [(0, 1), (1, 2), (0, 2)]:
        estimates = contrast_table(model, earlier, later)
        clinical = estimates.loc[estimates.State.eq('Clinical obesity')].drop(columns='State')
        clinical.insert(0, 'Selection', label)
        rows.append(clinical)

table = pd.concat(rows, ignore_index=True)
(ROOT / "R2_2_eta/Tables").mkdir(parents=True, exist_ok=True)
table.to_csv(ROOT / 'R2_2_eta/Tables/differenza_clinica.csv', index=False, float_format='%.2f')
