"""Age/sex standardisation and participant-level bootstrap for the three REFRAME states."""

from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy.special import softmax, logsumexp
from scipy.linalg import solve

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'Paper/Scripts'))
from coorti import PERIODS, STATES


def design(d, knots):
    """Period indicators, sex, restricted cubic spline in age, and sex–age interactions."""
    age = (d.access_age.to_numpy() - 50) / 10
    knots = (knots - 50) / 10
    male = d.sex.eq('M').to_numpy(dtype=float)

    def positive_cube(knot):
        return np.maximum(age - knot, 0) ** 3

    basis = [age]
    for knot in knots[:-2]:
        last, previous = knots[-1], knots[-2]
        curve = (
            positive_cube(knot)
            - (last - knot) / (last - previous) * positive_cube(previous)
            + (previous - knot) / (last - previous) * positive_cube(last)
        )
        basis.append(curve / (last - knots[0]) ** 2)
    basis = np.column_stack(basis)
    return np.column_stack(
        [
            np.ones(len(age)),
            d.period.eq(PERIODS[1]),
            d.period.eq(PERIODS[2]),
            male,
            basis,
            basis * male[:, None],
        ]
    )


def fit_model(X, Y, weight, start):
    """Weighted multinomial likelihood, with no excess adiposity as the reference state."""
    theta = start.copy()
    reference = np.zeros(len(Y))
    for _ in range(45):
        logits = np.column_stack([reference, X @ theta])
        probability = softmax(logits, axis=1)
        score = X.T @ ((Y - probability[:, 1:]) * weight[:, None])
        if np.max(np.abs(score)) / weight.sum() < 2e-9:
            return theta
        blocks = []
        for i in range(2):
            row = []
            for j in range(2):
                covariance = probability[:, i + 1] * (float(i == j) - probability[:, j + 1])
                row.append(X.T @ ((weight * covariance)[:, None] * X))
            blocks.append(row)
        step = solve(np.block(blocks), score.T.ravel(), assume_a='pos').reshape(2, X.shape[1]).T
        current = (weight * ((Y * logits[:, 1:]).sum(axis=1) - logsumexp(logits, axis=1))).sum()
        for scale in (1, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625):
            trial = theta + scale * step
            linear = X @ trial
            normaliser = logsumexp(np.column_stack([reference, linear]), axis=1)
            likelihood = weight @ ((Y * linear).sum(axis=1) - normaliser)
            if likelihood >= current - 1e-8:
                theta = trial
                break
        else:
            raise ArithmeticError('The multinomial likelihood did not increase.')
    raise ArithmeticError('The multinomial score did not converge.')


def estimate(X, Y, periods, weights, start):
    theta = fit_model(X, Y, weights, start)
    period_size = np.bincount(periods, weights=weights, minlength=3)
    reference = weights / (3 * period_size[periods])
    observed = np.column_stack([1 - Y.sum(axis=1), Y])
    crude = []
    standardised = []
    # Predict all reference-population members as though seen in each period.
    baseline = X @ theta - X[:, 1:3] @ theta[1:3]
    for p in range(3):
        in_period = periods == p
        crude.append(weights[in_period] @ observed[in_period] / period_size[p])
        linear = baseline if p == 0 else baseline + theta[p]
        predicted = softmax(np.column_stack([np.zeros(len(periods)), linear]), axis=1)
        standardised.append(reference @ predicted)
    estimates = np.stack([crude, standardised])
    if not np.allclose(estimates.sum(axis=2), 1):
        raise ValueError('Prevalences must sum to one in each period.')
    return theta, estimates


def fit_standardisation(d, output, seed, B=2000):
    d = d.sort_values(['patient_id', 'period']).reset_index(drop=True)
    if d.duplicated(['patient_id', 'period']).any():
        raise ValueError('Each person can contribute only one visit in a period.')
    knots = np.quantile(d.access_age, [0.05, 0.35, 0.65, 0.95])
    X = design(d, knots)
    Y = np.column_stack([d.state.eq(2), d.state.eq(3)]).astype(float)
    periods = pd.Categorical(d.period, categories=PERIODS).codes
    person, ids = pd.factorize(d.patient_id)
    n = len(ids)
    theta, point = estimate(X, Y, periods, np.ones(len(d)), np.zeros((X.shape[1], 2)))
    rng = np.random.default_rng(seed)
    reps = np.empty((B, 2, 3, 3))
    bootstrap_theta = np.empty((B, *theta.shape))
    for b in range(B):
        frequencies = np.bincount(rng.integers(0, n, size=n), minlength=n)
        bootstrap_theta[b], reps[b] = estimate(X, Y, periods, frequencies[person], theta)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        coefficients=theta,
        knots=knots,
        point=point,
        bootstrap=reps,
        bootstrap_coefficients=bootstrap_theta,
        seed=seed,
        participants=n,
        counts=counts(d),
        periods=PERIODS,
        states=list(STATES.values()),
        terms=[
            'Intercept',
            'Period 2',
            'Period 3',
            'Male',
            'Age',
            'Spline 1',
            'Spline 2',
            'Male × age',
            'Male × spline 1',
            'Male × spline 2',
        ],
    )


def counts(d):
    return (
        pd.crosstab(d.period, d.state)
        .reindex(index=PERIODS, columns=[1, 2, 3], fill_value=0)
        .to_numpy()
    )


def prevalence_table(model):
    with np.load(model) as f:
        sample = f['counts']
        point = 100 * f['point']
        low, high = 100 * np.quantile(f['bootstrap'], [0.025, 0.975], axis=0)
    rows = []
    for p, period in enumerate(PERIODS):
        for state, label in STATES.items():
            s = state - 1
            row = dict(Period=period, State=label, n=sample[p, s], N=sample[p].sum())
            for m, name in enumerate(['Crude', 'Standardised']):
                row[name + ' %'] = point[m, p, s]
                row[name + ' lower 95%'] = low[m, p, s]
                row[name + ' upper 95%'] = high[m, p, s]
            rows.append(row)
    return pd.DataFrame(rows)


def contrast_table(model, earlier=0, later=2):
    """Compare two periods using matched participant-bootstrap replicates."""
    with np.load(model) as f:
        prevalence = f['point']
        replicates = f['bootstrap']
        periods = f['periods']
        difference = 100 * (prevalence[:, later] - prevalence[:, earlier])
        changes = replicates[:, :, later] - replicates[:, :, earlier]
        low, high = 100 * np.quantile(changes, [0.025, 0.975], axis=0)
    rows = []
    for m, method in enumerate(['Crude', 'Standardised']):
        for s, name in enumerate(STATES.values()):
            rows.append(
                dict(
                    Method=method,
                    State=name,
                    Earlier_period=periods[earlier],
                    Earlier_prevalence_pct=100 * prevalence[m, earlier, s],
                    Later_period=periods[later],
                    Later_prevalence_pct=100 * prevalence[m, later, s],
                    Difference_pp=difference[m, s],
                    Lower95_pp=low[m, s],
                    Upper95_pp=high[m, s],
                )
            )
    return pd.DataFrame(rows)


def conditional_table(model, ages=(40, 50, 60)):
    """Predict state probabilities at fixed ages and sex, without averaging across profiles."""
    profiles = pd.DataFrame(
        [(p, sex, age) for p in PERIODS for sex in ('F', 'M') for age in ages],
        columns=['Period', 'Sex', 'Age_years'],
    )
    frame = profiles.rename(columns={'Period': 'period', 'Sex': 'sex', 'Age_years': 'access_age'})
    with np.load(model) as fit:
        X = design(frame, fit['knots'])

        def predict(coefficients):
            linear = X @ coefficients
            zeros = np.zeros((*linear.shape[:-1], 1))
            return softmax(np.concatenate([zeros, linear], axis=-1), axis=-1)

        point = 100 * predict(fit['coefficients'])
        lower, upper = np.quantile(
            100 * predict(fit['bootstrap_coefficients']), [0.025, 0.975], axis=0
        )
    result = profiles.copy()
    for j, state in enumerate(STATES.values()):
        result[state + ' % (95% CI)'] = [
            f'{p:.1f} ({lo:.1f}, {hi:.1f})'
            for p, lo, hi in zip(point[:, j], lower[:, j], upper[:, j])
        ]
    return result
