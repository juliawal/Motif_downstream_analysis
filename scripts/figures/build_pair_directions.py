#!/usr/bin/env python3
"""Draw descriptive MF--CpG directions from saved analysis tables.

No F4-support filter, significance encoding, background-gene comparison or
model fitting is used. See analysis/prepare_direction_results.py for preparation.
"""
import argparse
import json
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'utils'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
NEG, POS, MIX = '#0072B2', '#D55E00', '#8E5BA4'
INK, GREY, GRID = '#24333E', '#E4E8EC', '#EFF2F4'


def plot_pair_directions(input_dir, output_dir):
    pairs = pd.read_csv(input_dir / 'plotted_pair_directions.tsv', sep='\t')
    counts = pd.read_csv(input_dir / 'pair_direction_counts.tsv', sep='\t').set_index('direction')
    cohort_counts = pd.read_csv(input_dir / 'pair_cohort_counts.tsv', sep='\t')
    summary = json.loads((input_dir / 'pair_direction_summary.json').read_text())
    total = summary['repeated_tested_pairs']
    if int(counts.n_pairs.sum()) != total or int(cohort_counts.n_pairs.sum()) != total:
        raise ValueError('Direction totals disagree with the saved summary.')

    fig = plt.figure(figsize=(156 / 25.4, 132 / 25.4))
    fig.text(.025, .98, 'A  MF–CpG association directions across tested cohorts', weight='bold', fontsize=10, va='top')
    fig.text(.025, .924,
             f'{total} pairs tested in ≥ 2 cohorts; '
             f'{summary["repeated_percent_of_tested"]:.2f}% of all '
             f'{summary["tested_distinct_pairs"]} tested pairs', fontsize=8.5)

    ax = fig.add_axes([.025, .835, .95, .055])
    categories = [('negative', NEG, 'negative throughout'),
                  ('positive', POS, 'positive throughout'), ('mixed', MIX, 'mixed signs')]
    if counts.loc['zero', 'n_pairs']:
        categories.append(('zero', GREY, 'zero without both signs'))
    left = 0
    handles = []
    for key, colour, label in categories:
        value = int(counts.loc[key, 'n_pairs'])
        ax.barh(0, value, left=left, color=colour, height=1)
        handles.append(Patch(facecolor=colour, label=f'{value} {label}'))
        left += value
    ax.set(xlim=(0, total), ylim=(-.5, .5))
    ax.axis('off')
    fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(.018, .818),
               ncol=3 if len(categories) == 3 else 2, frameon=False, fontsize=8.1,
               columnspacing=1.35, handlelength=1, handletextpad=.5, borderaxespad=0)
    fig.text(.025, .746,
             f'Same nonzero sign in every tested cohort: {summary["same_sign_pairs"]} pairs '
             f'({summary["observed_same_sign_percent"]:.2f}%)', fontsize=8.5)

    tx = fig.add_axes([.025, .638, .95, .075])
    tx.axis('off')
    table = tx.table(cellText=[
        ['Cohorts tested per pair'] + cohort_counts.n_tested_cohorts.astype(str).tolist(),
        ['Number of pairs'] + [f'{n}' for n in cohort_counts.n_pairs],
    ], cellLoc='center', colWidths=[.32] + [.68 / len(cohort_counts)] * len(cohort_counts),
        bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_linewidth(1)
        cell.set_facecolor(GRID if row == 0 else 'white')
        cell.get_text().set_color(INK)
        if col == 0:
            cell.get_text().set_ha('left')
            cell.PAD = .03

    # Every finite tested coefficient uses the same text style. Grey cells have
    # no tested estimate; they are not a significance or support category.
    fig.text(.025, .567, 'B  MF–CpG association for selected pairs', weight='bold', fontsize=10)
    keys = pairs[['MF', 'CpG', 'pair_order']].drop_duplicates().sort_values('pair_order')
    tested = pairs.loc[pairs.tested].copy()
    preferred = ['CESC', 'LUAD', 'COAD', 'STAD', 'LGG', 'PAAD']
    present = set(tested.cohort)
    cohorts = [c for c in preferred if c in present] + sorted(present - set(preferred))
    row_lookup = {order: i for i, order in enumerate(keys.pair_order)}
    matrix = np.full((len(keys), len(cohorts)), np.nan)
    for row in tested.itertuples():
        matrix[row_lookup[row.pair_order], cohorts.index(row.cohort)] = row.a
    ax = fig.add_axes([.34, .272, .635, .235])
    cmap = LinearSegmentedColormap.from_list('mf_signed', [NEG, '#FBFAF7', POS])
    cmap.set_bad(GREY)
    im = ax.imshow(matrix, cmap=cmap, vmin=-1, vmax=1, aspect='auto')
    ax.set_yticks(range(len(keys)), [f'{m}  {c}' for m, c in zip(keys.MF, keys.CpG)], fontsize=8.6)
    ax.set_xticks(range(len(cohorts)), cohorts, fontsize=8.5)
    ax.tick_params(length=0, pad=5)
    ax.set_xticks(np.arange(-.5, len(cohorts), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(keys), 1), minor=True)
    ax.grid(which='minor', color='white', lw=1)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in tested.itertuples():
        ax.text(cohorts.index(row.cohort), row_lookup[row.pair_order],
                f'{row.a:+.2f}'.replace('-', '−'), ha='center', va='center',
                fontsize=9, weight='normal', color='white' if abs(row.a) > .58 else INK)
    cx = fig.add_axes([.395, .148, .525, .019])
    cb = fig.colorbar(im, cax=cx, orientation='horizontal', ticks=[-1, -.5, 0, .5, 1])
    cb.outline.set_visible(False)
    cb.set_label('Standardised coefficient $a$', fontsize=8.5, labelpad=3)
    fig.legend(handles=[Patch(facecolor=GREY, label='No tested estimate in that cohort')],
               loc='lower left', bbox_to_anchor=(.018, .015), frameon=False,
               fontsize=8.5, handlelength=1, borderaxespad=0)
    for extension in ['pdf', 'png']:
        fig.savefig(output_dir / f'pair_directions.{extension}', dpi=240, facecolor='white')
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / 'analysis_results')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'figures/main')
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        'font.family': 'DejaVu Sans', 'font.size': 8.5, 'text.color': INK,
        'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
        'pdf.fonttype': 42, 'axes.linewidth': .6,
    })
    plot_pair_directions(args.input_dir.resolve(), output_dir)
    print(f'Wrote pair_directions.pdf and pair_directions.png to {output_dir}')


if __name__ == '__main__':
    main()
