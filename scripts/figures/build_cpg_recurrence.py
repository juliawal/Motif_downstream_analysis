#!/usr/bin/env python3
"""Render TP63/CELSR3 from the two manuscript data tables; no pipeline inputs."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import argparse
import os
import tempfile

os.environ.setdefault('MPLCONFIGDIR', str(Path(tempfile.gettempdir()) / 'motif_f4_figure_mpl'))
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from figure_typography import readable_text

ROOT = Path(__file__).resolve().parents[2]


def draw_figure(tp, cel, output_dir):
    BLUE = '#0072B2'
    ORANGE = '#D55E00'
    GREY = '#8F9BA4'
    DARK = '#24313C'
    MID = '#60717D'
    GRID = '#E5EAF0'
    mpl.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 7.5,
        'axes.labelsize': 7.4, 'axes.titlesize': 9,
        'xtick.labelsize': 7, 'ytick.labelsize': 7.2,
        'legend.fontsize': 6.8, 'pdf.fonttype': 42, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.linewidth': .6, 'axes.edgecolor': '#A4AFB7',
        'text.color': DARK, 'axes.labelcolor': DARK,
        'xtick.color': MID, 'ytick.color': MID,
    })

    fig = plt.figure(figsize=(156 / 25.4, 104 / 25.4))
    fig.text(.025, .97, 'A  TP63 in CESC and ESCA', fontsize=9, weight='bold', va='top')
    fig.text(.63, .97, 'B  CELSR3–cg22702772', fontsize=9, weight='bold', va='top')
    ax = fig.add_axes([.13, .22, .43, .65])
    # Mix drawing order without moving points, so no support class is always on top.
    draw = tp.iloc[np.random.default_rng(42).permutation(len(tp))]
    colours = draw.support_count.map({0: GREY, 1: ORANGE, 2: BLUE})
    def darker(colour):
        return tuple(channel * .68 for channel in mpl.colors.to_rgb(colour))

    ax.scatter(draw.a_CESC, draw.a_ESCA, s=14, c=colours, alpha=1,
               edgecolors=[darker(colour) for colour in colours], linewidths=.35, zorder=3)
    ax.axhline(0, color=MID, lw=.6, ls='--', zorder=0)
    ax.axvline(0, color=MID, lw=.6, ls='--', zorder=0)
    ax.set_xlim(-1.03, 1.03)
    ax.set_ylim(-1.03, 1.03)
    ax.set_xticks([-1, -.5, 0, .5, 1])
    ax.set_yticks([-1, -.5, 0, .5, 1])
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('CESC coefficient $a$', linespacing=1.45)
    ax.set_ylabel('ESCA coefficient $a$', linespacing=1.45)
    handles = [Line2D([], [], marker='o', linestyle='none', color=col,
                      markeredgecolor=darker(col), markeredgewidth=.35, markersize=4,
                      label=f'{label} (n={count})')
               for col, label, count in [(BLUE, 'Both cohorts', int(tp.support_count.eq(2).sum())),
                                         (ORANGE, 'One cohort only', int(tp.support_count.eq(1).sum())),
                                         (GREY, 'Neither cohort', int(tp.support_count.eq(0).sum()))]]
    ax.legend(handles=handles, loc='upper left', frameon=False,
              title='F4 support in:', title_fontsize=6.9,
              borderaxespad=.3, handletextpad=.45, handlelength=.9, labelspacing=.35)
    bx = fig.add_axes([.77, .22, .20, .65])
    cohort_labels = ['Bladder\n(BLCA)', 'Breast\n(BRCA)', 'Head and neck\n(HNSC)',
                     'Lung adeno-\ncarcinoma (LUAD)', 'Prostate\n(PRAD)']
    for i, row in cel.iterrows():
        y = 4 - i
        bx.plot([0, row.a_min], [y, y], color='#CED6DD', lw=1)
        bx.scatter([row.a_min], [y], s=28, color=BLUE, zorder=3)
        bx.text(row.a_min + .027, y, f'+{row.a_min:.3f}', va='center', fontsize=7.1)
    bx.set_yticks(range(5), cohort_labels[::-1])
    bx.set_ylim(-.65, 4.65)
    bx.set_xlim(0, .99)
    bx.set_xticks([0, .5, 1])
    bx.grid(axis='x', color=GRID, lw=.5)
    bx.spines['left'].set_visible(False)
    bx.tick_params(axis='y', length=0, pad=6)
    bx.set_xlabel('Coefficient $a$', linespacing=1.45)

    readable_text(fig, minimum=8.5)

    fig.savefig(output_dir / 'cpg_recurrence.pdf', dpi=400, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--data-dir', type=Path, default=ROOT / 'analysis_results')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'figures/main')
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    tp = pd.read_csv(args.data_dir / 'tp63_cpg_recurrence.tsv', sep='\t')
    cel = pd.read_csv(args.data_dir / 'celsr3_cpg_recurrence.tsv', sep='\t')
    draw_figure(tp, cel, output_dir)
    print('Finished CpG recurrence figure.')


if __name__ == '__main__':
    main()
