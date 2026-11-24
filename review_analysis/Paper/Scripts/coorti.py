"""Build the historical and Methods cohorts from the shared visit-level extract."""

from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PERIODS = ['2005–2014', '2015–2019', '2020–2025']
STATES = {1: 'No excess adiposity', 2: 'Preclinical obesity', 3: 'Clinical obesity'}
MODEL_FIELDS = [
    'patient_id',
    'access_age',
    'state',
    'sex',
    'INSTI',
    'TAF',
    'IP',
    'hiv_duration_yrs',
    'phys_act_bin',
    'comorbidity_Dislipidemia',
    'liver_disease',
    'insulin_resistance',
    'sppb_low',
]


def events(methods=True):
    d = pd.read_csv(
        ROOT / 'Dati/visite.csv.gz',
        parse_dates=['date', 'exam_DEXA_date'],
        float_precision='round_trip',
    )
    d['access_age'] = d.age_days / 365.25
    d['access_date'] = d.date
    elasto = d.cap_med.notna() & d.e_med.notna()
    ct = d.FEGMIL.notna() & d.FIB4.notna()
    # Supplementary Table 1 requires both measurements in either diagnostic pair.
    hepatic = (elasto & d.cap_med.ge(248) & d.e_med.ge(8)) | (
        ct & d.FEGMIL.lt(1.1) & d.FIB4.gt(1.45)
    )
    d['liver_disease'] = hepatic.astype(float).where(elasto | ct)
    domains = [
        'comorbidity_COPD',
        'comorbidity_DM',
        'comorbidity_HTN',
        'comorbidity_Kidn',
        'comorbidity_CVD',
        'liver_disease',
    ]
    if d.excess_adiposity.isna().any() or d[domains].isna().all(axis=1).any():
        raise ValueError('Unclassifiable visit in input')
    d['organ_damage'] = d[domains].eq(1).any(axis=1)
    d['state'] = (1 + d.excess_adiposity.astype(int) * (1 + d.organ_damage.astype(int))).where(
        ~d.death, 4
    )
    if methods:
        d = d.loc[d.adult & d.date.dt.year.between(2005, 2025)].copy()
        d = d.loc[d.patient_id.isin(d.loc[~d.death, 'patient_id'])].copy()
    d['study_entry_age'] = d.groupby('patient_id').access_age.transform('min')
    d['period'] = pd.cut(d.date.dt.year, [2004, 2014, 2019, 2025], labels=PERIODS).astype('object')
    d['anthro_known'] = d.WHtR.notna() & d.waist_circumference__cm.notna()
    d['anthro_pair'] = d.WHtR.ge(0.50) & d.waist_circumference__cm.ge(
        np.where(d.sex.eq('M'), 102, 88)
    )
    d['postmenopausal'] = d.menopausa.map({'SI': 1, 'NO': 0}).where(d.sex.eq('F'))
    d['smoking'] = d.smoking.map({'Intense (>10day)': 1, 'Moderate (<=10day)': 1, 'No': 0})
    for drug in ('NRTI', 'NNRTI'):
        d[drug] = (
            d[drug]
            .astype('string')
            .str.upper()
            .map({'TRUE': 1, 'FALSE': 0, 'SI': 1, 'NO': 0, '1': 1, '0': 0})
        )
    viral = d['HIV_viral_load__copies/mL']
    d['vl_suppressed'] = (viral.isna() | viral.le(50)).astype(int)
    return d.sort_values(['patient_id', 'access_age']).reset_index(drop=True)


def retained(methods=True, returns=False, calendar_days=False):
    d = events(methods)
    # The delivered fits use differences in attained age, as in the original program.
    # Integer days are used in the editor's boundary check.
    clock = d.age_days.to_numpy() if calendar_days else d.access_age.to_numpy()
    min_gap = 30 if calendar_days else 30 / 365.25
    max_gap = 3 * 365.25 if calendar_days else 3
    rows = []
    for ix in d.groupby('patient_id', sort=False).indices.values():
        selected = [ix[0]]
        for i in ix[1:]:
            gap = clock[i] - clock[selected[-1]]
            if gap < min_gap:
                continue
            if gap > max_gap and not returns:
                break
            selected.append(i)
        # Single initial observations enter cross-sectional analyses. Repeated
        # trajectories must keep at least two observations after spacing.
        if len(ix) == 1 or len(selected) >= 2:
            rows.extend(selected)
    return d.iloc[rows].copy()


def living(returns=False):
    d = retained(returns=returns)
    return d.loc[~d.death].copy()


def baseline(d):
    return d.drop_duplicates('patient_id', keep='first')


def last_period(d):
    return d.drop_duplicates(['patient_id', 'period'], keep='last')


def random_period(d, seed=52001):
    """Uniform draw per person and period, reproducible in overlapping cohorts."""
    ordered = d.sort_values(['patient_id', 'period', 'date']).reset_index(drop=True)
    rows = []
    for (person, period), positions in ordered.groupby(
        ['patient_id', 'period'], sort=True
    ).indices.items():
        rng = np.random.default_rng([seed, int(person), PERIODS.index(period)])
        rows.append(positions[rng.integers(len(positions))])
    return ordered.iloc[rows].copy()


def all_three(d):
    return d.loc[d.groupby('patient_id').period.transform('nunique').eq(3)].copy()


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser(description='Export the pre-imputation data for one multistate cohort.')
    parser.add_argument('cohort', choices=['storica', 'metodi'])
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    panel = retained(methods=args.cohort == 'metodi')
    multi = panel[panel.groupby('patient_id').state.transform('size').ge(2)]
    multi[MODEL_FIELDS].to_csv(args.output, index=False)
