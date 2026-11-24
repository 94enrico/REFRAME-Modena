"""Paired raw-scale measurements from people with at least two eligible DXA dates."""

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np


def plot_repeated(exam, destination):
    data = exam[exam.groupby('patient_id').DEXA_date.transform('nunique').ge(2)].copy()
    # Put each clinical cutoff on a bin edge.
    data['fat_bin'] = 2 * np.floor((data.fat_percent - 1) / 2) + 2
    data['vat_bin'] = 20 * np.floor((data.VAT_dexa - 10) / 20) + 20
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 5.5), sharey=True, sharex=True)
    colors = ['#007773', '#63539a']
    for ax, (sex, title, cutoff), color in zip(
        axes, [('F', 'Women', 35), ('M', 'Men', 25)], colors
    ):
        d = data[data.sex.eq(sex)]
        cells = d.groupby(['fat_bin', 'vat_bin']).size().rename('exams').reset_index()
        ax.grid(axis='y', color='#e6e9e9', linewidth=0.65)
        ax.set_axisbelow(True)
        ax.axhline(150, color='#59666a', linestyle='--', linewidth=0.9)
        ax.axvline(cutoff, color='#59666a', linestyle='--', linewidth=0.9)
        ax.scatter(
            cells.fat_bin,
            cells.vat_bin,
            s=1.65 * cells.exams,
            c=color,
            alpha=0.57,
            linewidths=0,
            zorder=3,
        )
        ax.set_title(
            f'{title}   ·   {d.patient_id.nunique():,} people, {len(d):,} exams',
            loc='left',
            fontsize=10.5,
            fontweight='bold',
            pad=11,
        )
        ax.set_xlabel('Total body fat (%)')
        ax.spines[['right', 'top']].set_visible(False)
        ax.spines[['left', 'bottom']].set_color('#a1a8a9')
        ax.set_yscale('symlog', linthresh=150, linscale=1.2, base=2)
        ax.minorticks_off()
    top = max(300, 20 * (np.floor(data.VAT_dexa.max() / 20) + 1))
    axes[0].set_ylim(0, top * 1.1)
    ticks = [y for y in (0, 50, 100, 150, 300, 600, 1200, 2400, 4800) if y <= top * 1.1]
    axes[0].set_yticks(ticks, [str(y) for y in ticks])
    axes[0].set_xlim(0, 10 * np.ceil((data.fat_bin.max() + 1) / 10))
    axes[0].set_ylabel('Visceral adipose tissue area (cm²)')
    handles = [
        axes[1].scatter(
            [],
            [],
            s=1.65 * n,
            color='#49676a',
            alpha=0.62,
            linewidths=0,
            label=f'{n} exam' if n == 1 else f'{n} exams',
        )
        for n in (1, 20, 100)
    ]
    handles.append(
        Line2D([0], [0], color='#59666a', linestyle='--', linewidth=1, label='Adiposity thresholds')
    )
    fig.legend(
        handles=handles,
        loc='lower center',
        bbox_to_anchor=(0.5, 0.073),
        ncol=4,
        frameon=False,
        fontsize=9.5,
    )
    fig.suptitle(
        'Visceral and total adiposity on repeated DXA examinations',
        x=0.085,
        y=0.98,
        ha='left',
        fontsize=17,
        fontweight='bold',
    )
    fig.text(
        0.085,
        0.035,
        'Bubble area is proportional to the number of exams in each 2-percentage-point × 20-cm² bin.',
        fontsize=9.2,
        color='#344649',
    )
    fig.text(
        0.085,
        0.005,
        'VAT is shown on a linear scale to 150 cm² and on a logarithmic scale above it.',
        fontsize=9.2,
        color='#344649',
    )
    fig.subplots_adjust(left=0.085, right=0.98, bottom=0.245, top=0.84, wspace=0.12)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        destination, metadata={'Title': 'Visceral and total adiposity on repeated DXA examinations'}
    )
    plt.close(fig)
