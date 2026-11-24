"""Produce Table 1 at the first visit and Table 2 at the last visit in each period."""

from pathlib import Path
from coorti import living, baseline, last_period
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
PERIODS = ['2005–2014', '2015–2019', '2020–2025']
STATES = ['No excess adiposity', 'Preclinical obesity', 'Clinical obesity']
# Display precision follows the published tables. A blank precision denotes a binary variable.
FIELDS = {
    'age': ('Age, years', 'access_age', 1),
    'entry': ('Age at study entry, years', 'study_entry_age', 1),
    'end': ('Age at the study end date, years', 'study_end_age', 1),
    'male': ('Male sex', 'male', None),
    'female': ('Female sex', 'female', None),
    'smoke': ('Current smoking', 'smoking', None),
    'meno': ('Postmenopausal at the visit (women), n/N (%)', 'postmenopausal', None),
    'pa': ('Physical activity', 'phys_act_bin', None),
    'bmi': ('BMI, kg/m²', 'BMI', 1),
    'waist': ('Waist circumference, cm', 'waist_circumference__cm', 1),
    'whtr': ('Waist-to-height ratio', 'WHtR', 2),
    'cd4': ('Current CD4 count, cells/µL', 'CD4_absolute__cells/uL', 0),
    'nadir': ('Nadir CD4 count, cells/µL', 'Nadir_CD4__cells/uL', 0),
    'ratio': ('CD4/CD8 ratio', 'CD4/CD8__ratio', 2),
    'viral': ('Suppressed HIV viral load', 'vl_suppressed', None),
    'duration': ('Time since HIV diagnosis, years', 'hiv_duration_yrs', 1),
    'nrti': ('Current NRTI use', 'NRTI', None),
    'nnrti': ('Current NNRTI use', 'NNRTI', None),
    'pi': ('Current PI use', 'IP', None),
    'insti': ('Current INSTI use', 'INSTI', None),
    'taf': ('Current TAF use', 'TAF', None),
    'vat': ('Visceral adipose tissue, cm²', 'VAT_dexa', 1),
    'dm': ('Diabetes', 'comorbidity_DM', None),
    'htn': ('Hypertension', 'comorbidity_HTN', None),
    'dys': ('Dyslipidemia', 'comorbidity_Dislipidemia', None),
    'liver': ('Hepatic organ involvement', 'liver_disease', None),
    'ckd': ('Chronic kidney disease', 'comorbidity_Kidn', None),
    'copd': ('Chronic obstructive pulmonary disease', 'comorbidity_COPD', None),
    'low_sppb': ('SPPB score <11', 'sppb_under11', None),
    'cvd': ('Cardiovascular disease', 'comorbidity_CVD', None),
    'ir': ('Insulin resistance', 'insulin_resistance', None),
    'homa': ('HOMA-IR', 'HOMA', 1),
    'sppb': ('SPPB score', 'SPPB', 0),
}
FIELDS['g1'] = ('Current first-generation INSTI use, n/N (%)', 'g1', None)
FIELDS['g2'] = ('Current second-generation INSTI use, n/N (%)', 'g2', None)
for i, t in enumerate(
    ['Normal BMI (<25 kg/m²)', 'Overweight (25 to <30 kg/m²)', 'Obesity (≥30 kg/m²)'], 1
):
    FIELDS[f'bmi{i}'] = (t, f'bmi{i}', None)
for i, t in enumerate(STATES, 1):
    FIELDS[f'state{i}'] = (t, f'state{i}', None)
T1 = [
    (
        'Demographic, anthropometric and lifestyle characteristics',
        'age male meno smoke pa waist whtr',
    ),
    ('HIV-related characteristics', 'cd4 nadir ratio viral duration nrti nnrti pi insti g1 g2 taf'),
    ('Cardiometabolic characteristics and physical performance', 'dys ir homa vat sppb'),
]
T2 = [
    (
        'Demographic, anthropometric and lifestyle characteristics',
        'entry end female meno smoke pa bmi waist whtr',
    ),
    ('HIV-related characteristics', 'duration nadir cd4 ratio viral nrti nnrti pi insti taf'),
    ('Cardiometabolic characteristics', 'vat dm htn dys liver ckd cvd copd low_sppb'),
    ('BMI categories', 'bmi1 bmi2 bmi3'),
    ('REFRAME states', 'state1 state2 state3'),
]


def cell(values, decimals, group_size):
    if decimals is not None:
        if values.empty:
            return 'NA [n=0]'
        q1, median, q3 = values.quantile([0.25, 0.50, 0.75])
        result = f'{median:.{decimals}f} ({q1:.{decimals}f}, {q3:.{decimals}f})'
        return result if len(values) == group_size else f'{result} [n={len(values):,}]'
    if not values.isin([0, 1]).all():
        raise ValueError(f'Unexpected binary value in {values.name}')
    if values.empty:
        return '0/0'
    count = int(values.sum())
    return f'{count:,}/{len(values):,} ({100 * count / len(values):.1f}%)'


def p_value(groups, binary):
    if binary:
        counts = [[x.sum(), len(x) - x.sum()] for x in groups]
        test = stats.chi2_contingency(counts, correction=False)
        if test.expected_freq.min() < 5:
            raise ValueError('The chi-squared approximation is unsuitable for this row.')
    else:
        test = stats.kruskal(*groups)
    return '<0.001' if test.pvalue < 0.001 else f'{test.pvalue:.3f}'


def make(groups, sections, compare=False):
    output = []
    headers = {name: f'{name} (N = {len(group):,})' for name, group in groups.items()}
    if compare:
        overall = len(groups['Overall'])
        for name in STATES:
            n = len(groups[name])
            headers[name] = f'{name} (N = {n:,}, {100 * n / overall:.1f}%)'
    for section, names in sections:
        for key in names.split():
            label, field, decimals = FIELDS[key]
            row = {'Section': section, 'Characteristic': label}
            samples = []
            for group_name, group in groups.items():
                eligible = group[group.sex.eq('F')] if key == 'meno' else group
                available = eligible[field].dropna()
                row[headers[group_name]] = cell(available, decimals, len(eligible))
                if group_name != 'Overall':
                    samples.append(available)
            if compare:
                row['P value'] = p_value(samples, binary=decimals is None)
            output.append(row)
    return pd.DataFrame(output)


def write_tables(destination):
    destination.mkdir(parents=True, exist_ok=True)
    d = living()
    d['male'] = d.sex.eq('M').astype(int)
    d['female'] = 1 - d.male
    d['sppb_under11'] = d.SPPB.lt(11).astype(float).where(d.SPPB.notna())
    for n in (1, 2, 3):
        d[f'bmi{n}'] = d.bmi_cat.eq(n).astype(float).where(d.bmi_cat.notna())
        d[f'state{n}'] = d.state.eq(n).astype(int)
    entry = baseline(d)
    groups = {'Overall': entry}
    groups.update({name: entry[entry.state.eq(n)] for n, name in enumerate(STATES, 1)})
    last = last_period(d)
    periods = {p: last[last.period.eq(p)] for p in PERIODS}
    make(groups, T1, compare=True).to_csv(destination / 'Table_1.csv', index=False)
    make(periods, T2).to_csv(destination / 'Table_2.csv', index=False)


if __name__ == '__main__':
    write_tables(ROOT / 'Paper/Tables')
