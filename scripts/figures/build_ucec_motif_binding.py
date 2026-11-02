#!/usr/bin/env python3
"""Build the UCEC FIMO manuscript figure from saved F1 results.

From the project root, in the Motif environment, run:
    python scripts/figures/build_ucec_motif_binding.py

From another working directory, use the absolute path to this script.
Default input and output paths are resolved relative to the script.
Defaults read analysis_results/ and write ucec_motif_binding.pdf
to figures/main/, with plotted source rows in analysis_results/.
Use --evidence-dir and --output-dir to select other locations. No models,
hypothesis tests or multiple-testing corrections are recomputed.
Panels A and B are stacked at 156 x 179 mm for a portrait manuscript page.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from project_layout import figure_data_dir
import argparse
import csv
import os

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from matplotlib.transforms import ScaledTranslation

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    with path.open() as f:
        return [dict(r, source_line=i) for i, r in enumerate(csv.DictReader(f, delimiter='\t'), 2)]


def write_table(path, rows):
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)


def draw_portrait(union, lookup, settings_order, groups):
    """Stack both panels at final print size (156 x 179 mm)."""
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9,
                         'axes.labelsize': 9, 'xtick.labelsize': 8.5,
                         'ytick.labelsize': 9, 'pdf.fonttype': 42,
                         'text.color': '#24313C', 'axes.labelcolor': '#24313C',
                         'axes.edgecolor': '#A8B0B6', 'xtick.color': '#4B5963'})
    blue, grey, dark = '#0072B2', '#8B949C', '#24313C'
    cmap = LinearSegmentedColormap.from_list(
        'F1_difference', ['#C45D19', '#FBFBFA', '#0072B2'])
    norm = Normalize(-55, 55)
    fig = plt.figure(figsize=(156 / 25.4, 179 / 25.4))

    # Millimetres keep the two row grids and their text at final print size.
    def axes(left, bottom, width, height):
        return fig.add_axes([left / 156, bottom / 179, width / 156, height / 179])

    def heading(x, y, text, **kwargs):
        return fig.text(x / 156, y / 179, text, **kwargs)

    labels_a = axes(1, 106, 19, 49)
    rates = axes(22, 106, 66, 49)
    targets = axes(94, 106, 26, 49)
    backgrounds = axes(123, 106, 32, 49)
    labels_b = axes(1, 14, 19, 65)
    settings = axes(22, 14, 133, 65)
    for ax in [labels_a, rates, targets, backgrounds, labels_b, settings]:
        ax.set_ylim(len(union) - .5, -.5)
        ax.set_yticks([])
        ax.spines[:].set_visible(False)
        ax.set_xticks([])
    for ax in [labels_a, labels_b, targets, backgrounds]:
        ax.set_xlim(0, 1)
    rates.set_xlim(-1, 63)
    rates.set_xticks([0, 20, 40, 60])
    rates.spines['bottom'].set_visible(True)
    rates.grid(axis='x', color='#E4E8EB', linewidth=.5)
    rates.set_axisbelow(True)
    rates.set_xlabel('Percentage of analysed CpG windows with at least one motif hit', labelpad=3)
    settings.set_xlim(-.5, 2.5)

    heading(1, 175, 'A  Union target and background hit rates in UCEC',
            fontsize=10.5, weight='bold')
    heading(1, 169.5, f'UCEC · FIMO · all {len(union)} evaluated factors', fontsize=9)
    fig.legend(handles=[
        Line2D([], [], color=blue, marker='o', linestyle='none', markersize=4.5,
               label='Target'),
        Line2D([], [], markeredgecolor=grey, markerfacecolor='white', marker='s',
               linestyle='none', markersize=4.5, label='Background')],
        loc='center left', bbox_to_anchor=(22 / 156, 162 / 179),
        borderaxespad=0, frameon=False, ncol=2, fontsize=8.5,
        handletextpad=.5, columnspacing=1.3)
    for x, text in [(107, 'Target'), (139, 'Background')]:
        heading(x, 158.5, text, ha='center', fontsize=9)
    heading(124, 163, 'hits / analysed CpG windows', ha='center', fontsize=8.5)

    heading(1, 89.5, 'B  UCEC hit rate differences across settings', fontsize=10.5, weight='bold')
    for x, (analysis, label) in enumerate(zip(settings_order, ['S1', 'S2', 'S3'])):
        centre = 22 + (x + .5) * 133 / 3
        heading(centre, 84.5, label, ha='center', fontsize=9.5, weight='bold')
        heading(centre, 80.5, f"p = {float(groups[analysis]['p_value_sign_test']):.3g}",
                ha='center', fontsize=9)

    for yy, row in enumerate(union):
        for ax in [labels_a, labels_b]:
            ax.text(0, yy, row['MF'], va='center', weight='normal', fontsize=9)
        target, background = 100 * float(row['target_hit_rate']), 100 * float(row['bg_hit_rate'])
        rates.plot([background, target], [yy, yy], color='#B1BBC2', lw=1.25, zorder=2)
        rates.scatter(background, yy, s=24, marker='s', facecolor='white',
                      edgecolor=grey, linewidth=1, zorder=4)
        rates.scatter(target, yy, s=22, color=blue, linewidth=.4, zorder=5)
        for ax, hit_key, count_key in [
                (targets, 'n_target_hits', 'n_target_cpgs'),
                (backgrounds, 'n_background_hits', 'n_background_cpgs')]:
            ax.text(.5, yy, f'{int(row[hit_key])}/{int(row[count_key])}',
                    ha='center', va='center', fontsize=9)
        for x, analysis in enumerate(settings_order):
            effect = 100 * float(lookup[(analysis, row['MF_ensembl'])]['difference'])
            settings.add_patch(Rectangle((x - .49, yy - .47), .98, .94,
                                         facecolor=cmap(norm(effect)), edgecolor='white', lw=.35))
            label = '0.0' if effect == 0 else f'{effect:+.1f}'
            # Centre the visible numeral outlines, including their signs.
            # A typographic line box reserves descender space below digits,
            # which otherwise makes va='center' look too high in the cell.
            font = FontProperties(family='DejaVu Sans', weight='normal', size=9.5)
            bounds = TextPath((0, 0), label, prop=font).get_extents()
            offset = ScaledTranslation(-(bounds.x0 + bounds.width / 2) / 72,
                                       -(bounds.y0 + bounds.height / 2) / 72,
                                       fig.dpi_scale_trans)
            settings.text(x, yy, label, transform=settings.transData + offset,
                          ha='left', va='baseline', fontproperties=font,
                          color='white' if abs(effect) > 35 else dark)
    cax = axes(47, 9, 83, 1.6)
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax,
                      orientation='horizontal', ticks=[-50, 0, 50])
    cb.ax.tick_params(labelsize=8.5, length=2, pad=1)
    cb.set_label('Target hit rates (in %) − background hit rates (in %)',
                 fontsize=9, labelpad=1)
    cb.outline.set_visible(False)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--evidence-dir', type=Path, default=ROOT / 'analysis_results',
                        help='Directory containing the saved binding measurements and original group tests')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'figures/main',
                        help='Directory for PDF; source tables use analysis_results/ for project outputs, otherwise data/')
    args = parser.parse_args()
    OUT = args.output_dir.resolve()
    DATA = figure_data_dir(OUT)
    DATA.mkdir(parents=True, exist_ok=True)
    BASE = args.evidence_dir.resolve()
    MEASUREMENTS = BASE / 'binding_measurements.tsv'
    GROUPS = BASE / 'binding_group_tests.tsv'
    all_rows = read(MEASUREMENTS)
    rows = [r for r in all_rows if r['cohort']=='UCEC' and r['track']=='FIMO']
    lookup = {(r['analysis'], r['MF_ensembl']):r for r in rows}
    union = sorted([r for r in rows if r['analysis']=='union' and r['status']=='tested'],
                   key=lambda r:-float(r['difference']))
    settings_order = ['rank_1','rank_2','rank_3']
    group_rows = [r for r in read(GROUPS) if r['cohort']=='UCEC' and r['track']=='FIMO'
                  and r['analysis'] in settings_order]
    groups = {r['analysis']:r for r in group_rows}
    source_rows = [lookup[(a,r['MF_ensembl'])] for r in union for a in ['union',*settings_order]]
    # Only saved p-values are used; no sign test or multiplicity correction is run.
    write_table(DATA / 'ucec_motif_binding_measurements.tsv',source_rows)
    write_table(DATA / 'ucec_motif_binding_tests.tsv',[groups[a] for a in settings_order])

    fig = draw_portrait(union, lookup, settings_order, groups)
    fig.savefig(OUT / 'ucec_motif_binding.pdf', dpi=300, facecolor='white')
    plt.close(fig)
    print('Finished figure and data export.')


if __name__ == '__main__':
    main()
