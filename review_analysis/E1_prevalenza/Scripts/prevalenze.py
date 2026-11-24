"""Crude, conditional and marginal prevalence for the random and last visit in each period."""

from argparse import ArgumentParser
from pathlib import Path
import sys
import pandas as pd
from threadpoolctl import threadpool_limits

FOLDER = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(FOLDER.parent / "Paper/Scripts"))
from standardizzazione import fit_standardisation, prevalence_table, conditional_table
from coorti import living, random_period, last_period

MODELS = {'Random visit': 'casuale', 'Last visit': 'ultima'}


def export():
    (FOLDER / "Tables").mkdir(parents=True, exist_ok=True)
    crude, conditional, standardised = [], [], []
    for label, filename in MODELS.items():
        model = FOLDER / f'Modelli/{filename}.npz'
        values = prevalence_table(model)
        for method, result in [('Crude', crude), ('Standardised', standardised)]:
            identifiers = (
                ['Period', 'State', 'n', 'N'] if method == 'Crude' else ['Period', 'State']
            )
            names = [method + ' %', method + ' lower 95%', method + ' upper 95%']
            sub = values[identifiers + names].copy()
            sub = sub.rename(
                columns=dict(zip(names, ['Prevalence_pct', 'Lower95_pct', 'Upper95_pct']))
            )
            sub.insert(0, 'Selection', label)
            result.append(sub.round(2))
        fixed = conditional_table(model)
        fixed.insert(0, 'Selection', label)
        conditional.append(fixed)
    for rows, name in [
        (crude, 'grezze'),
        (conditional, 'condizionate'),
        (standardised, 'standardizzate'),
    ]:
        pd.concat(rows, ignore_index=True).to_csv(FOLDER / f'Tables/{name}.csv', index=False)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__)
    p.add_argument('--rifitta', action='store_true')
    args = p.parse_args()
    for name, select in [('casuale', random_period), ('ultima', last_period)]:
        model = FOLDER / f'Modelli/{name}.npz'
        if args.rifitta or not model.exists():
            with threadpool_limits(1):
                fit_standardisation(select(living()), model, seed=20261171)
    export()
