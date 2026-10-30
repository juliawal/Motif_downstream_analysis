#!/usr/bin/env python3
"""Draw the pair-direction figure from saved analysis tables.

Inputs: analysis_results/paper_data/.
Outputs: PDF in figures/main/ and plotted data in
analysis_results/paper_data/. No models or significance tests are rerun.
"""
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import project_layout as artifacts
from project_layout import figure_data_dir

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
NEG = "#0072B2"
POS = "#D55E00"
MIX = "#8E5BA4"
INK = "#24333E"
GREY = "#E4E8EC"
GRID = "#EFF2F4"


def save_figure(fig, output_dir, name):
    artifacts.save_figure(fig, output_dir / f"{name}.pdf", dpi=320, facecolor="white")
    plt.close(fig)


def plot_pair_directions(input_dir, output_dir, source_dir):
    """Show cross-cohort signs, their F4-supported subsets and selected pairs."""
    pairs = pd.read_csv(input_dir / 'plotted_pair_directions.tsv', sep='\t')
    repeated = pd.read_csv(input_dir / 'all_repeated_pairs.tsv', sep='\t')
    strict = repeated[repeated.strict_recurrent_supported]
    summary = {'repeated_tested_pairs': len(repeated),
               'observed_same_sign_percent': 100 * repeated.same_sign_all_tested.mean(),
               'strict_negative': int(strict.n_positive.eq(0).sum()),
               'strict_positive': int(strict.n_negative.eq(0).sum())}
    summary.update({'repeated_' + sign: int(repeated.direction_all_tested.eq(sign).sum())
                    for sign in ['negative', 'positive', 'mixed']})
    supported_counts = (repeated.loc[repeated.n_supported_cohorts.ge(2)]
        .groupby('direction_all_tested').size())
    subsets = pd.DataFrame([{
        'direction': key,
        'n_tested_pairs': summary['repeated_' + key],
        'n_pairs_with_F4_support_in_at_least_two_cohorts': int(supported_counts.get(key, 0)),
    } for key in ['negative', 'positive', 'mixed']])
    artifacts.write_dataframe(subsets, source_dir / 'pair_support_subsets.tsv', sep='\t', index=False)
    cohort_counts = (repeated.groupby('n_tested_cohorts').size()
        .reindex(range(2, int(repeated.n_tested_cohorts.max()) + 1), fill_value=0)
        .rename('n_pairs').reset_index())
    artifacts.write_dataframe(cohort_counts, source_dir / 'pair_cohort_counts.tsv', sep='\t', index=False)
    keys = pairs[['MF', 'CpG', 'pair_order']].drop_duplicates().sort_values('pair_order')
    preferred = ['CESC', 'COAD', 'ESCA', 'STAD', 'LUAD', 'LUSC', 'SARC',
                 'SKCM_METASTATIC', 'SKCM_PRIMARY', 'LGG', 'PAAD']
    cohorts = [c for c in preferred if c in set(pairs.loc[pairs.tested, 'cohort'])]
    matrix = np.full((len(keys), len(cohorts)), np.nan)
    fig = plt.figure(figsize=(156 / 25.4, 156 / 25.4))

    # A: all repeated pairs and the supported subsets, on a shared count scale.
    # Separate the panel letter so the subtitle shares the title's exact left edge.
    fig.text(.025, .98, 'A', weight='bold', fontsize=10, va='top')
    fig.text(.055, .98, 'Expression–methylation associations\nof identical gene–CpG pairs across cohorts',
             weight='bold', fontsize=10, va='top', linespacing=1.2)
    fig.text(.055, .9, f'{summary["repeated_tested_pairs"]} pairs tested in at least two cohorts', fontsize=9)
    ax = fig.add_axes([.025, .786, .95, .07])
    left = 0
    for key, col in [('negative', NEG), ('positive', POS), ('mixed', MIX)]:
        value = summary['repeated_' + key]
        ax.barh(0, value, left=left, color=col, height=.8)
        if key != 'mixed':
            ax.text(left + value / 2, 0, f'{value} consistently {key}', ha='center',
                    va='center', color='white', fontsize=8.5, weight='bold')
        # Negative subsets align left; positive and mixed subsets align right.
        supported = int(supported_counts[key])
        subset_left = left if key == 'negative' else left + value - supported
        ax.barh(1, supported, left=subset_left, color=col, height=.8)
        if key != 'mixed':  # The purple subset is too narrow for its count label.
            ax.text(subset_left + supported / 2, 1, str(supported), ha='center',
                    va='center', color='white', fontsize=8.5, weight='bold')
        for edge in [subset_left, subset_left + supported]:
            ax.plot([edge, edge], [.4, .6], color=col, lw=.65, alpha=.65)
        left += value
    gap_left = summary['strict_negative']
    gap_right = summary['repeated_negative'] + summary['repeated_positive'] - summary['strict_positive']
    ax.text((gap_left + gap_right) / 2, 1, 'Subset with F4 support in ≥ 2 cohorts',
            ha='center', va='center', fontsize=8.5)
    ax.set(xlim=(0, summary['repeated_tested_pairs']), ylim=(-.5, 1.6))
    ax.axis('off')
    same_sign = summary['repeated_positive'] + summary['repeated_negative']
    fig.text(.025, .768, f'Same sign: {same_sign} pairs ({summary["observed_same_sign_percent"]:.2f}%)',
             fontsize=8.5)
    fig.text(.53, .768,
             f'Opposite signs: {summary["repeated_mixed"]} pairs ({100 - summary["observed_same_sign_percent"]:.2f}%)',
             fontsize=8.5, color=MIX)
    tx = fig.add_axes([.025, .689, .95, .064])
    tx.axis('off')
    counts_table = tx.table(cellText=[
        ['Cohorts tested per pair'] + cohort_counts.n_tested_cohorts.astype(str).tolist(),
        ['Number of pairs'] + cohort_counts.n_pairs.astype(str).tolist()],
        cellLoc='center', colWidths=[.32] + [.68 / len(cohort_counts)] * len(cohort_counts),
        bbox=[0, 0, 1, 1])
    counts_table.auto_set_font_size(False)
    counts_table.set_fontsize(8)
    for (row, col), cell in counts_table.get_celld().items():
        cell.set_edgecolor('white')
        cell.set_linewidth(1)
        cell.set_facecolor(GRID if row == 0 else 'white')
        cell.get_text().set_color(INK)
        if col == 0:
            cell.get_text().set_ha('left')
            cell.PAD = .03

    # B: coefficients and support states of the selected gene--CpG pairs.
    fig.text(.025, .66, 'B  Examples of specific gene–CpG pairs', weight='bold', fontsize=9.3)
    ax = fig.add_axes([.335, .267, .64, .36])
    for row in pairs[pairs.tested].itertuples():
        matrix[row.pair_order, cohorts.index(row.cohort)] = row.a
    cmap = LinearSegmentedColormap.from_list('mf_signed', [NEG, '#FBFAF7', POS])
    cmap.set_bad(GREY)
    im = ax.imshow(matrix, cmap=cmap, vmin=-1, vmax=1, aspect='auto')
    ax.set_yticks(range(len(keys)), [f'{m}  {c}' for m, c in zip(keys.MF, keys.CpG)], fontsize=8.2)
    cohort_labels = [c.replace('SKCM_METASTATIC', 'SKCM-MET').replace('SKCM_PRIMARY', 'SKCM-PRI')
                     for c in cohorts]
    ax.set_xticks(range(len(cohorts)), cohort_labels, rotation=55, ha='right', fontsize=8.1)
    ax.tick_params(length=0, pad=4)
    ax.set_xticks(np.arange(-.5, len(cohorts), 1), minor=True)
    ax.set_yticks(np.arange(-.5, len(keys), 1), minor=True)
    ax.grid(which='minor', color='white', lw=1)
    ax.tick_params(which='minor', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for row in pairs[pairs.tested].itertuples():
        x = cohorts.index(row.cohort)
        y = row.pair_order
        label = f'{row.a:+.2f}'.replace('+0.', '+.').replace('-0.', '−.')
        ax.text(x, y, label, ha='center', va='center', fontsize=8.0,
                weight='bold' if row.supported else 'normal',
                color='white' if abs(row.a) > .58 else INK)
        if row.supported_all_four:
            ax.add_patch(Rectangle((x - .46, y - .45), .92, .9, fill=False, edgecolor=INK, lw=.8))
    cx = fig.add_axes([.40, .142, .50, .018])
    cb = fig.colorbar(im, cax=cx, orientation='horizontal', ticks=[-1, -.5, 0, .5, 1])
    cb.outline.set_visible(False)
    cb.set_label('Expression–methylation association (standardised F4 $a$)', fontsize=8, labelpad=4)
    fig.text(.025, .063, 'Bold number: F4 support in the union analysis.', fontsize=8)
    fig.text(.025, .044, 'Bold number + border: also supported in each setting (S1, S2 and S3).', fontsize=8)
    fig.text(.025, .025, 'Regular number: tested in the union analysis, without F4 support.', fontsize=8)
    fig.text(.025, .006, 'Grey cell: no tested estimate in that cohort.', fontsize=8)
    save_figure(fig, output_dir, 'pair_directions')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir", type=Path,
        default=ROOT / "analysis_results/paper_data",
        help="Directory containing the prepared direction tables",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=ROOT / "figures/main",
        help="Figure directory; source tables follow the shared project layout",
    )
    args = parser.parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    source_dir = figure_data_dir(output_dir)
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5, "text.color": INK,
        "axes.labelcolor": INK, "xtick.color": INK, "ytick.color": INK,
        "pdf.fonttype": 42, "axes.linewidth": .6,
    })
    plot_pair_directions(input_dir, output_dir, source_dir)
    print(f"Wrote pair_directions.pdf to {output_dir}")
    print(f"Plotted source data: {source_dir}")


if __name__ == "__main__":
    main()
