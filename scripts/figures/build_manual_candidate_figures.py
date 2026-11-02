"""Render binding/network and cohort F4 figures for the current 13-gene panel.

F1 uses the saved section 2.9 evaluation; F2 retains the original cohort points.
No group evaluation is performed.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import argparse
import json
import os
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from figure_typography import readable_text
from cohort_f4_figure import draw_cohort_f4

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir', type=Path, default=ROOT / 'analysis_results')
parser.add_argument('--binding-results', type=Path, help='Override the prepared section 2.9 binding measurements.')
parser.add_argument('--output-dir', type=Path, default=ROOT / 'figures/supplement')
parser.add_argument('--only-overview', action='store_true', help='Render only the binding/network panels.')
args = parser.parse_args()
SOURCE = args.data_dir.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, GREY = '#0072B2', '#C66D00', '#78848D'
PALE, POS, NEG = '#F3F6F8', '#B45C23', '#0072B2'
mpl.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8,
    'axes.labelsize':8, 'axes.titlesize':9, 'xtick.labelsize':7.2,
    'ytick.labelsize':7.6, 'legend.fontsize':7.2, 'axes.linewidth':.6,
    'pdf.fonttype':42, 'savefig.dpi':240})

def read(name):
    return pd.read_csv(SOURCE / name, sep='\t')


config = json.loads((ROOT / 'scripts/paper_cohorts_and_genes.json').read_text())
expected_genes = set(config['manual_prioritised'])
binding_results = args.binding_results or SOURCE / 'manual_binding_measurements.tsv'
f1 = pd.read_csv(binding_results, sep='\t').sort_values('MF').reset_index(drop=True)
if not f1.MF.is_unique or set(f1.MF) != expected_genes:
    raise ValueError('S3 requires saved F1 results for exactly the current manual_prioritised genes.')
genes = f1.MF.tolist()
contexts = read('combined_evidence/union_evidence_contexts.tsv')
contexts = contexts.loc[contexts.MF.isin(genes)].copy()
if contexts.duplicated(['MF', 'cohort']).any() or set(contexts.MF) != expected_genes:
    raise ValueError('S3 requires unique original cohort contexts covering every current gene.')
for measure in ['harmonic', 'pagerank']:
    values = contexts[measure + '_degree_bin_excess']
    if not np.isfinite(values).all():
        raise ValueError(f'Missing or non-finite saved {measure} differences in S3 cohort contexts.')
    contexts[measure + '_difference_sign'] = np.select(
        [values.gt(0), values.lt(0)], ['positive', 'negative'], default='zero')
f2_sign_counts = {measure: {sign: int(contexts[measure + '_difference_sign'].eq(sign).sum())
    for sign in ['positive', 'negative', 'zero']} for measure in ['harmonic', 'pagerank']}


def clean(ax):
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=5)
    ax.tick_params(axis='x',length=3,width=.6)
    ax.set_axisbelow(True)

def rows(ax,labels=True):
    for i in range(len(genes)):
        if i%2==0:ax.axhspan(i-.5,i+.5,color=PALE,zorder=0)
    ax.set(ylim=(len(genes)-.4,-.6),yticks=range(len(genes)),yticklabels=genes if labels else [])
    clean(ax)


def f1_panel(ax):
    rows(ax)
    for i,r in f1.iterrows():
        if r.status=='tested':
            color=BLUE if r.track=='FIMO' else ORANGE
            ax.plot([100*r.bg_hit_rate,100*r.target_hit_rate],[i,i],color='#BDC7CE',lw=1.4)
            ax.scatter(100*r.bg_hit_rate,i,marker='s',s=14,color=GREY,zorder=3)
            ax.scatter(100*r.target_hit_rate,i,s=24,color=color,edgecolor='white',lw=.4,zorder=4)
        else:
            ax.text(5,i,'Not measured',va='center',color=GREY,fontsize=7.2)
    ax.set(xlim=(-3,80),xticks=[0,20,40,60,80],xlabel='CpG windows with a hit (%)')
    ax.grid(axis='x',color='#E1E7EB',lw=.5)

f1handles=[Line2D([],[],marker='o',ls='',color=BLUE,label='FIMO target'),Line2D([],[],marker='o',ls='',color=ORANGE,label='ChIP target'),Line2D([],[],marker='s',ls='',color=GREY,label='Background')]
f2handles=[Line2D([],[],marker='o',ls='',color=POS,label='Positive'),Line2D([],[],marker='o',ls='',color=NEG,label='Negative')]
if any(counts['zero'] for counts in f2_sign_counts.values()):
    f2handles.append(Line2D([],[],marker='o',ls='',color=GREY,label='Zero'))

def f2_panel(axes):
    for m,ax in zip(['harmonic','pagerank'],axes):
        rows(ax,labels=(m=='harmonic'))
        for i,g in enumerate(genes):
            d=contexts.loc[contexts.MF.eq(g)].sort_values('cohort')
            offset=np.linspace(-.22,.22,len(d)) if len(d)>1 else np.array([0.])
            values=d[m+'_degree_bin_excess'].to_numpy()*(1e4 if m=='pagerank' else 1)
            colours=np.where(values>0,POS,np.where(values<0,NEG,GREY))
            ax.scatter(values,i+offset,s=14,color=colours,edgecolor='white',lw=.25,zorder=3)
        ax.axvline(0,color=GREY,ls='--',lw=.8)
        ax.set_xlabel('Proximity difference'+('\n(× 10⁻⁴)' if m=='pagerank' else ''))
        if m=='harmonic':ax.set(xlim=(-.28,1.62),xticks=[0,.5,1,1.5])
        else:ax.set(xlim=(-.34,4.8),xticks=[0,2,4])


def render_overview():
    # Supplementary two-panel figure; quantitative results are unchanged.
    with mpl.rc_context({'font.size': 9.5, 'axes.labelsize': 9,
                         'xtick.labelsize': 9, 'ytick.labelsize': 9.5,
                         'legend.fontsize': 9}):
        width, height = 156, 110
        fig = plt.figure(figsize=(width / 25.4, height / 25.4))

        def axes(left, bottom, axis_width, axis_height):
            return fig.add_axes([left / width, bottom / height,
                                 axis_width / width, axis_height / height])

        a = axes(20, 20, 54, 75)
        h = axes(101, 20, 20, 75)
        pr = axes(126, 20, 29, 75)
        f1_panel(a)
        f2_panel([h, pr])
        fig.text(1 / width, 106 / height, 'A  F1 Binding Site Evaluation', weight='bold', fontsize=10)
        fig.text(82 / width, 106 / height, 'B  F2 PPI Proximity', weight='bold', fontsize=10)
        a.set_xlabel('CpG windows with a hit (%)', labelpad=4)
        h.set_xticks([0, 1.5])
        h.set_xlabel('Harmonic\ndifference', labelpad=4)
        pr.set_xticks([0, 2, 4])
        pr.set_xlabel('PageRank\ndifference (× 10⁻⁴)', labelpad=4)
        legend_style = dict(frameon=False, fontsize=8.5, handlelength=.8,
                            handletextpad=.35, columnspacing=.8, borderaxespad=0)
        fig.legend(handles=f1handles, loc='center left',
                   bbox_to_anchor=(1 / width, 101 / height), ncol=3, **legend_style)
        fig.legend(handles=f2handles, loc='center left',
                   bbox_to_anchor=(82 / width, 101 / height),
                   ncol=len(f2handles), **legend_style)
        readable_text(fig, minimum=8.5)
        fig.savefig(OUT / 'manual_candidate_evidence.pdf', facecolor='white')
        plt.close(fig)


render_overview()

# Detailed F4 counts/directions at final manuscript size.
# Keep the fixed canvas: tight bounding boxes would change the embedded font size.
def render_f4_cohort_details():
    fig = draw_cohort_f4(read('manual_prioritised_contexts.tsv'))
    fig.savefig(OUT / 'manual_candidate_associations.pdf', facecolor='white')
    plt.close(fig)

if not args.only_overview:
    with mpl.rc_context():
        render_f4_cohort_details()

print('Saved supplementary manual-candidate figure(s).')
