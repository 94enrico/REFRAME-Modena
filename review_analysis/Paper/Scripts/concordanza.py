"""Supplementary Table 2 and the BMI agreement figures. Bootstrap unit is the participant."""

import numpy as np
import pandas as pd
from coorti import ROOT, PERIODS, STATES, living, last_period
from agreement_draw import select_one

BMI = ['<25', '25 to <30', '≥30']
METRICS = [
    'Exact agreement',
    'Quadratic weighted kappa',
    'REFRAME category higher than BMI',
    'REFRAME category lower than BMI',
]
WEIGHTS = (np.arange(3)[:, None] - np.arange(3)[None, :]) ** 2


def cells(frame):
    return (frame.state.to_numpy(dtype=int) - 1) * 3 + frame.bmi_cat.to_numpy(dtype=int) - 1


def matrix(frame):
    return np.bincount(cells(frame), minlength=9).reshape(3, 3)


def measures(table):
    joint = table / table.sum()
    expected = np.outer(joint.sum(axis=1), joint.sum(axis=0))
    kappa = 1 - (WEIGHTS * joint).sum() / (WEIGHTS * expected).sum()
    return np.array([np.trace(joint), kappa, np.tril(joint, -1).sum(), np.triu(joint, 1).sum()])


def bootstrap(frame, latest, seed):
    ids = pd.Index(frame.patient_id.drop_duplicates())
    n = len(ids)
    sizes = frame.groupby('patient_id', sort=True).size().to_numpy()
    offsets = np.r_[0, np.cumsum(sizes)[:-1]]
    codes = cells(frame)
    # Each person's four tables contain their last visits by period and all visits pooled.
    clustered = np.zeros((n, 4, 9), dtype=int)
    np.add.at(clustered, (ids.get_indexer(frame.patient_id), np.full(len(frame), 3), codes), 1)
    period = pd.Categorical(latest.period, categories=PERIODS).codes
    np.add.at(clustered, (ids.get_indexer(latest.patient_id), period, cells(latest)), 1)
    rng = np.random.default_rng(seed)
    replicates = np.empty((2000, 5, len(METRICS)))
    for b in range(len(replicates)):
        sample = rng.integers(0, n, size=n)
        # The random-visit analysis redraws one visit per resampled participant.
        selected = offsets[sample] + rng.integers(sizes[sample])
        replicates[b, 0] = measures(np.bincount(codes[selected], minlength=9).reshape(3, 3))
        frequencies = np.bincount(sample, minlength=n)
        totals = (frequencies @ clustered.reshape(n, -1)).reshape(4, 3, 3)
        for j, table in enumerate(totals, start=1):
            replicates[b, j] = measures(table)
    return np.quantile(replicates, [0.025, 0.975], axis=0)


def analyse():
    frame = living().dropna(subset=['bmi_cat'])
    random_visit = select_one(frame)
    latest = last_period(frame)
    rows, figure_data = [], []
    for j, (group, sex) in enumerate([('Overall', None), ('Men', 'M'), ('Women', 'F')]):

        def subgroup(data):
            return data if sex is None else data[data.sex.eq(sex)]

        all_visits, last = subgroup(frame), subgroup(latest)
        selections = {'One random visit': subgroup(random_visit)}
        selections.update({p: last[last.period.eq(p)] for p in PERIODS})
        selections['All BMI-complete visits'] = all_visits
        lower, upper = bootstrap(all_visits, last, 52201 + j)
        for k, (selection, data) in enumerate(selections.items()):
            table = matrix(data)
            estimates = measures(table)
            for m, name in enumerate(METRICS):
                rows.append(
                    {
                        'Selection': selection,
                        'Group': group,
                        'Participants': data.patient_id.nunique(),
                        'Observations': len(data),
                        'Measure': name,
                        'Estimate': estimates[m],
                        'Lower 95%': lower[k, m],
                        'Upper 95%': upper[k, m],
                    }
                )
            for state, bmi in np.ndindex(3, 3):
                figure_data.append(
                    dict(
                        Selection=selection,
                        Group=group,
                        State=STATES[state + 1],
                        BMI=BMI[bmi],
                        n=table[state, bmi],
                    )
                )
    return pd.DataFrame(rows), pd.DataFrame(figure_data)


if __name__ == '__main__':
    from figure import agreement_figures

    results, figure_data = analyse()
    (ROOT / "Paper/Tables").mkdir(parents=True, exist_ok=True)
    results.to_csv(ROOT / 'Paper/Tables/Supplementary_table_2.csv', index=False)
    agreement_figures(figure_data, results)
