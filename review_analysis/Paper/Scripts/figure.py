"""Draw the numerical figures used in the article and its supplement from the recalculated tables."""

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'Paper/Figures'
TAB = ROOT / 'Paper/Tables'
STATES = ['No excess adiposity', 'Preclinical obesity', 'Clinical obesity']
BMI = ['<25', '25 to <30', '≥30']
PERIOD = ['One random visit', '2005–2014', '2015–2019', '2020–2025']
GROUPS = ['Overall', 'Men', 'Women']
COLORS = ['#9AC5AB', '#EDB372', '#923C53']
COLOR_BMI = ['#AECCBE', '#D8AE70', '#8C4E63']
plt.rcParams.update(
    {
        'font.family': 'DejaVu Sans',
        'pdf.fonttype': 42,
        'font.size': 9,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.titlesize': 10,
    }
)


def save(fig, name):
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        OUT / (name + '.pdf'), facecolor='white', metadata={'Title': name.replace('_', ' ')}
    )
    plt.close(fig)


def count_table(d, s, g):
    return (
        d[(d.Selection == s) & (d.Group == g)]
        .pivot(index='State', columns='BMI', values='n')
        .loc[STATES, BMI]
        .to_numpy()
    )


def stacked_percentages(ax, t, colors, labels, show_small_top=False):
    v = 100 * t / t.sum(0)
    bottom = np.zeros(v.shape[1])
    for i, (color, rank) in enumerate(zip(colors, v)):
        ax.bar(range(len(bottom)), rank, bottom=bottom, width=0.65, color=color)
        for j, p in enumerate(rank):
            if p >= 4.4:
                ax.text(
                    j,
                    bottom[j] + p / 2,
                    f'{p:.1f}',
                    ha='center',
                    va='center',
                    fontsize=8,
                    color='white' if i == 2 else '#182326',
                )
            elif show_small_top and i == 2:
                ax.text(j, 101.5, f'{p:.1f}', ha='center', va='bottom', fontsize=8, color=colors[2])
        bottom += rank
    ax.set_xticks(range(len(labels)), labels)
    ax.set_ylim(0, 107 if show_small_top else 100)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.tick_params(axis='both', length=0)
    ax.spines[['left', 'bottom']].set_color('#CCCFCF')


def legend(fig, items, colors, y):
    fig.legend(
        handles=[Patch(facecolor=c, label=s) for s, c in zip(items, colors)],
        loc='lower center',
        bbox_to_anchor=(0.5, y),
        ncol=len(items),
        frameon=False,
    )


def agreement_figures(source, stats):
    fig, axes = plt.subplots(1, 4, figsize=(12, 3.85), sharey=True)
    for i, (ax, p) in enumerate(zip(axes, PERIOD)):
        t = count_table(source, p, 'Overall')
        stacked_percentages(
            ax, t, COLORS, [f'{s}\n(n={n:,})' for s, n in zip(['<25', '25–<30', '≥30'], t.sum(0))]
        )
        ax.set_title(
            ('A   ' if i == 0 else 'B   ' if i == 1 else '') + p, loc='left', fontweight='bold'
        )
    axes[0].set_ylabel('Participants within BMI category (%)')
    fig.suptitle(
        'REFRAME states within BMI categories', x=0.062, ha='left', fontsize=16, fontweight='bold'
    )
    fig.text(0.52, 0.15, 'BMI (kg/m²)', ha='center')
    legend(fig, STATES, COLORS, 0.015)
    fig.subplots_adjust(left=0.065, right=0.99, bottom=0.3, top=0.79, wspace=0.18)
    save(fig, 'Figure_1_BMI')
    cmap = LinearSegmentedColormap.from_list('joint', ['#F4F7F4', '#A7CFC4', '#3F8B81', '#194F50'])
    fig, axes = plt.subplots(3, 4, figsize=(12.1, 8.8))
    for i, g in enumerate(GROUPS):
        for j, p in enumerate(PERIOD):
            ax = axes[i, j]
            t = count_table(source, p, g)
            fr = 100 * t / t.sum()
            ax.imshow(fr, cmap=cmap, vmin=0, vmax=65, aspect='auto')
            for h, k in np.ndindex(3, 3):
                ax.text(
                    k,
                    h,
                    f'{fr[h, k]:.1f}',
                    ha='center',
                    va='center',
                    color='white' if fr[h, k] > 35 else '#172B30',
                    fontsize=9,
                )
            vals = stats[(stats.Selection == p) & (stats.Group == g)].set_index('Measure').Estimate
            ax.set_title(f'{p}  ·  n={t.sum():,}', fontsize=10)
            ax.set_xticks(range(3), ['<25', '25–<30', '≥30'], fontsize=8)
            ax.set_yticks(
                range(3),
                ['No excess', 'Preclinical', 'Clinical'] if not j else ['', '', ''],
                fontsize=8,
            )
            ax.set_xlabel(
                f"Agreement {100 * vals['Exact agreement']:.1f}%    "
                f"κ {vals['Quadratic weighted kappa']:.3f}\n"
                f"REFRAME higher / lower   {100 * vals['REFRAME category higher than BMI']:.1f}% / "
                f"{100 * vals['REFRAME category lower than BMI']:.1f}%",
                fontsize=8,
            )
            ax.tick_params(length=0)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if not j:
                ax.set_ylabel(g, fontsize=11, fontweight='bold')
    fig.suptitle(
        'Agreement between REFRAME and BMI categories',
        x=0.095,
        ha='left',
        fontsize=17,
        fontweight='bold',
    )
    fig.text(
        0.095,
        0.025,
        'Cells show the percentage of each panel in both categories. '
        'Darker shading indicates a larger percentage. Columns show BMI in kg/m².',
        fontsize=9,
    )
    fig.subplots_adjust(left=0.1, right=0.99, top=0.9, bottom=0.12, wspace=0.14, hspace=0.8)
    save(fig, 'Figure_2_agreement')
    for number, state in zip([2, 3, 4], [2, 1, 0]):
        fig, ax = plt.subplots(figsize=(7.1, 4.25))
        t = np.column_stack([count_table(source, p, 'Overall')[state] for p in PERIOD])
        stacked_percentages(
            ax,
            t,
            COLOR_BMI,
            [f'{p}\n(n={n:,})' for p, n in zip(PERIOD, t.sum(0))],
            show_small_top=state == 0,
        )
        ax.set_ylabel('Participants within REFRAME state (%)')
        ax.set_title(
            'BMI among participants with ' + STATES[state].lower(),
            loc='left',
            fontweight='bold',
            fontsize=13,
            pad=12,
        )
        legend(fig, ['BMI <25', 'BMI 25 to <30', 'BMI ≥30'], COLOR_BMI, 0.012)
        fig.subplots_adjust(left=0.12, right=0.98, bottom=0.23, top=0.86)
        save(fig, f'Supplementary_Figure_{number}_BMI')


def multistate_figure():
    h = pd.read_csv(TAB / 'HR.csv').query('Model == 3')
    cols = ['1-2', '1-3', '1-4', '2-3', '2-4', '3-4', '2-1', '3-1', '3-2']
    vars_ = [
        'access_age',
        'hiv_tertileMedium',
        'hiv_tertileLong',
        'phys_act_bin',
        'INSTI',
        'TAF',
        'IP',
        'comorbidity_Dislipidemia',
        'sex',
    ]
    labels = [
        'Age (per year)',
        'HIV duration, middle vs lowest',
        'HIV duration, highest vs lowest',
        'Physical activity',
        'INSTI use',
        'TAF use',
        'PI use',
        'Dyslipidemia',
        'Male vs female sex',
    ]
    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.axis('off')
    ax.set_xlim(-2.95, 8.55)
    ax.set_ylim(9.75, -1.5)
    for j, e in enumerate(cols):
        ax.text(
            j,
            -0.62,
            e.replace('-', ' → ').replace('4', 'D'),
            va='center',
            ha='center',
            fontweight='bold',
        )
        for i, v in enumerate(vars_):
            part = h[(h.Transition == e) & (h.Variable == v)]
            if part.empty:
                s = '—'
                color = '#F4F5F5'
            else:
                r = part.iloc[0]
                if pd.isna(r.HR):
                    s = 'NE'
                    color = '#F4F5F5'
                else:
                    s = f'{r.HR:.2f}\n({r.Lower95:.2f}–{r.Upper95:.2f})'
                    color = (
                        ('#D8EADF' if e in ['2-1', '3-1', '3-2'] else '#F1DCD8')
                        if r.Lower95 > 1 or r.Upper95 < 1
                        else '#FAFAF8'
                    )
            ax.add_patch(
                plt.Rectangle((j - 0.485, i - 0.46), 0.97, 0.92, fc=color, ec='white', lw=0.6)
            )
            ax.text(j, i, s, ha='center', va='center', fontsize=8.2, color='#172B30')
    for i, l in enumerate(labels):
        ax.text(-0.58, i, l, ha='right', va='center', fontsize=9.5)
    fig.suptitle(
        'Adjusted hazard ratios for transitions between REFRAME states',
        x=0.04,
        ha='left',
        y=0.98,
        fontsize=16,
        fontweight='bold',
    )
    fig.text(
        0.04,
        0.095,
        '1   No excess adiposity       2   Preclinical obesity       3   Clinical obesity       D   Death',
        fontsize=9.2,
    )
    note = (
        'Each cell shows the HR and 95% CI. Shading marks intervals excluding 1 using '
        'unrounded estimates. A dash marks a term outside that transition model.'
    )
    if h.HR.isna().any():
        note += '  NE indicates an unestimable HR.'
    fig.text(0.04, 0.045, note, fontsize=8.1)
    fig.subplots_adjust(left=0.01, right=0.99, top=0.89, bottom=0.095)
    save(fig, 'Figure_3_multistate')


if __name__ == '__main__':
    multistate_figure()
