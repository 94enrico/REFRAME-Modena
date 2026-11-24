"""Compare VAT and total fat on paired adult DXA examinations and at the baseline visit."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import living, baseline, STATES

OUT = ROOT / 'R3_2_DEXA/Tables'
OUT.mkdir(parents=True, exist_ok=True)
KEY = ['patient_id', 'DEXA_date']
FAT = 'dexa_total_fat_mass__g'
WEIGHT = 'dexa_total_weight__kg'
raw = pd.read_csv(
    ROOT / 'Dati/dexa.csv.gz', parse_dates=['DEXA_date'], float_precision='round_trip'
)
components = raw[[FAT, 'dexa_total_lean_mass__g', 'dexa_total_bone_mass__g']].sum(
    axis=1, min_count=3
)
stored = raw[WEIGHT].dropna()
grams = components.loc[stored.index]
if stored.empty or not (stored.gt(0) & grams.gt(0)).all():
    raise ValueError('Incomplete DXA mass components. The unit cannot be verified.')
scales = [s for s in (1, 1000) if np.allclose(s * stored, grams, atol=0.5, rtol=0)]
if len(scales) != 1:
    raise ValueError('DXA total weight has ambiguous or inconsistent units.')
scale = scales[0]

visits = living()
first = baseline(visits)
records = raw[raw.DEXA_date.dt.year.between(2005, 2025)].dropna(subset=['VAT_dexa', FAT, WEIGHT])
records = records.merge(
    first[['patient_id', 'sex', 'date', 'age_days']], on='patient_id', validate='m:1'
)
birth = records.date - pd.to_timedelta(records.age_days, unit='D')
records = records[records.DEXA_date.ge(birth + pd.DateOffset(years=18))]
records = records.drop_duplicates(KEY + ['VAT_dexa', FAT, WEIGHT])
valid = (
    records[WEIGHT].gt(0)
    & records[FAT].gt(0)
    & records[FAT].le(scale * records[WEIGHT])
    & records.VAT_dexa.ge(0)
)
exam = records[valid & ~records.duplicated(KEY, keep=False)].copy()
exam['high_vat'] = exam.VAT_dexa.ge(150)
exam['fat_percent'] = 100 * exam[FAT] / (scale * exam[WEIGHT])
exam['high_fat'] = exam.fat_percent.ge(np.where(exam.sex.eq('M'), 25, 35))

# Sample one paired examination per person for agreement. A different sample is
# used for the figure, which describes repeated observations.
exam = exam.sort_values(KEY).reset_index(drop=True)
sample = exam.groupby('patient_id', sort=True).sample(n=1, random_state=52017)
sample = sample.sort_values('patient_id')


def agreement(values):
    """Absolute-agreement single-measure ICC(A,1) and kappa for the two thresholds."""
    n, k = values.shape
    grand = values.mean()
    subjects = values.mean(axis=1)
    methods = values.mean(axis=0)
    ms_subjects = k * np.square(subjects - grand).sum() / (n - 1)
    ms_methods = n * np.square(methods - grand).sum() / (k - 1)
    residuals = values - subjects[:, None] - methods[None, :] + grand
    ms_error = np.square(residuals).sum() / ((n - 1) * (k - 1))
    icc = (ms_subjects - ms_error) / (
        ms_subjects + (k - 1) * ms_error + k * (ms_methods - ms_error) / n
    )
    above = values >= 1
    observed = (above[:, 0] == above[:, 1]).mean()
    rates = above.mean(axis=0)
    expected = np.prod(rates) + np.prod(1 - rates)
    return np.array([icc, (observed - expected) / (1 - expected)])


def median_iqr(values):
    q1, median, q3 = values.quantile([0.25, 0.5, 0.75])
    return f'{median:.1f} ({q1:.1f}, {q3:.1f})'


def proportion(n, total):
    return f'{n}/{total} ({100 * n / total:.1f}%)'


table = {}
for label, group in [
    ('Overall', sample),
    ('Women', sample[sample.sex.eq('F')]),
    ('Men', sample[sample.sex.eq('M')]),
]:
    thresholds = np.where(group.sex.eq('M'), 25, 35)
    ratios = np.column_stack([group.VAT_dexa / 150, group.fat_percent / thresholds])
    count = pd.crosstab(group.high_vat, group.high_fat).reindex(
        index=[False, True], columns=[False, True], fill_value=0
    )
    point = agreement(ratios)
    rng = np.random.default_rng(52018)
    estimates = [
        agreement(ratios[rng.integers(0, len(group), size=len(group))]) for _ in range(2000)
    ]
    lo, hi = np.quantile(estimates, [0.025, 0.975], axis=0)
    n = len(group)
    a, b = 'A. Paired measurements', 'B. Agreement'
    table[f'{label} (N = {n:,})'] = {
        (a, 'VAT (cm²), median (Q1, Q3)'): median_iqr(group.VAT_dexa),
        (a, 'Total body fat (%), median (Q1, Q3)'): median_iqr(group.fat_percent),
        (a, 'High VAT, n/N (%)'): proportion(group.high_vat.sum(), n),
        (a, 'High total body fat, n/N (%)'): proportion(group.high_fat.sum(), n),
        (b, 'Neither high, n'): count.iat[0, 0],
        (b, 'High VAT only, n'): count.iat[1, 0],
        (b, 'High total body fat only, n'): count.iat[0, 1],
        (b, 'Both high, n'): count.iat[1, 1],
        (b, 'Binary agreement, n/N (%)'): proportion(count.iat[0, 0] + count.iat[1, 1], n),
        (
            b,
            'ICC(A,1) for threshold-scaled values (95% CI)',
        ): f'{point[0]:.3f} ({lo[0]:.3f}, {hi[0]:.3f})',
        (
            b,
            'Cohen kappa for threshold classifications (95% CI)',
        ): f'{point[1]:.3f} ({lo[1]:.3f}, {hi[1]:.3f})',
    }
pd.DataFrame(table).rename_axis(['Panel', 'Measure']).to_csv(OUT / 'confronto_DEXA.csv')

from bubble_plot import plot_repeated

plot_repeated(exam, ROOT / 'R3_2_DEXA/Figures/VAT_massa_grassa.pdf')

# Use the exact examination already linked to the first clinical visit by the paper.
cols = KEY + ['VAT_dexa', FAT, WEIGHT, 'high_fat']
paired = first.merge(
    exam[cols].rename(columns={'DEXA_date': 'exam_DEXA_date', 'VAT_dexa': 'paired_vat'}),
    on=['patient_id', 'exam_DEXA_date'],
    validate='1:1',
)
same = (
    np.isclose(paired.exam_VAT_dexa, paired.paired_vat)
    & np.isclose(paired['exam_' + FAT], paired[FAT])
    & np.isclose(paired['exam_' + WEIGHT], paired[WEIGHT])
)
paired = paired[same].copy()
excess = paired.high_fat | paired.anthro_pair
paired['fat_state'] = np.select(
    [excess, paired.anthro_known], [2 + paired.organ_damage.astype(int), 1], default=0
)
counts = pd.crosstab(paired.state, paired.fat_state).reindex(
    index=[1, 2, 3], columns=[1, 2, 3, 0], fill_value=0
)
counts = counts.rename(index=STATES, columns={**STATES, 0: 'Unclassifiable'})
counts.index.name = 'State with VAT'
counts.to_csv(OUT / 'riclassificazione_basale.csv')
print(f'DXA mass unit check: {len(stored)} rows, {scale} gram(s) per stored unit.')
print(
    f'Paired participants: {len(sample)}. Baseline visits with a matched exam: {len(paired)}/{len(first)}.'
)
