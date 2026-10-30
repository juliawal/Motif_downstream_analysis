"""Render both manual-panel figures solely from analysis_results/paper_data/."""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import project_layout as artifacts
from project_layout import figure_output_dir
import argparse
import os
os.environ.setdefault('MPLCONFIGDIR', '/tmp/motif_26_figure_mpl')
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from figure_typography import readable_text

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data-dir', type=Path, default=ROOT / 'analysis_results/paper_data')
args = parser.parse_args()
SOURCE = args.data_dir.resolve()
OUT = ROOT / 'figures/main'
DATA = ROOT / 'analysis_results/paper_data'
OUT.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)
BLUE, ORANGE, PURPLE, GREY = '#0072B2', '#C66D00', '#79549D', '#78848D'
PALE, POS, NEG = '#F3F6F8', '#B45C23', '#0072B2'
mpl.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8,
    'axes.labelsize':8, 'axes.titlesize':9, 'xtick.labelsize':7.2,
    'ytick.labelsize':7.6, 'legend.fontsize':7.2, 'axes.linewidth':.6,
    'pdf.fonttype':42, 'savefig.dpi':240})

def read(name):
    return pd.read_csv(SOURCE / name, sep='\t')


f1 = read('manual_binding_measurements.tsv').sort_values('MF').reset_index(drop=True)
genes = f1.MF.tolist()
f2g = read('manual_network_group_tests.tsv').set_index('measure')
terms = read('manual_enrichment_terms.tsv')
matrix = read('manual_enrichment_memberships.tsv').set_index('MF').loc[genes, terms.plot_id].to_numpy()
contexts = read('manual_network_contexts.tsv').merge(
    read('manual_association_contexts.tsv'), on=['MF', 'MF_ensembl', 'cohort'], validate='one_to_one')
f2_sign_counts = {measure: {sign: int(contexts[measure + '_difference_sign'].eq(sign).sum())
    for sign in ['positive', 'negative', 'zero']} for measure in ['harmonic', 'pagerank']}


def save(fig,name):
    artifacts.save_figure(fig, figure_output_dir(OUT, name)/f'{name}.pdf',bbox_inches='tight',pad_inches=.09,facecolor='white')
    plt.close(fig)

def clean(ax):
    ax.spines[['top','right','left']].set_visible(False)
    ax.tick_params(axis='y',length=0,pad=5)
    ax.tick_params(axis='x',length=3,width=.6)
    ax.set_axisbelow(True)

def rows(ax,labels=True):
    for i in range(17):
        if i%2==0:ax.axhspan(i-.5,i+.5,color=PALE,zorder=0)
    ax.set(ylim=(16.6,-.6),yticks=range(17),yticklabels=genes if labels else [])
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


def f3_panel(ax):
    for i in range(17):
        for j in range(6):
            ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor=PURPLE if matrix[i,j] else '#F1EDF5',edgecolor='white',lw=1.1))
    ax.set(xlim=(-.55,5.55),ylim=(16.6,-.6),xticks=range(6),xticklabels=terms.plot_id,yticks=range(17),yticklabels=genes)
    ax.tick_params(length=0)
    for s in ax.spines.values():s.set_visible(False)
    for j,n in enumerate(matrix.sum(axis=0)):ax.text(j,18.35,str(n),ha='center',fontsize=7.1)
    ax.set_xlabel('Six selected GO terms (purple = membership)',labelpad=22,fontsize=7.5)

def f4_panel(ax):
    rows(ax)
    for i,g in enumerate(genes):
        selected=contexts.loc[contexts.MF.eq(g)].sort_values('cohort')
        tested=selected.loc[selected.n_tested_pairs.gt(0)]
        offset=np.linspace(-.22,.22,len(tested)) if len(tested)>1 else np.array([0.])
        ax.scatter(tested.f4_support_percent,i+offset,s=23,color=BLUE,edgecolor='white',lw=.4,zorder=3)
        ax.text(1.14,i,f'{len(tested)}/{len(selected)}',transform=ax.get_yaxis_transform(),ha='center',va='center',fontsize=7.3)
    ax.set(xlim=(-3,103),xticks=[0,25,50,75,100],xlabel='Supported / tested MF–CpG pairs (%)')
    ax.grid(axis='x',color='#E1E7EB',lw=.5)


def render_overview():
    # Final manuscript width. Each evaluation repeats the same gene order;
    # the two network measures share labels within panel B.
    with mpl.rc_context({'font.size':9, 'axes.labelsize':8.5,
                         'xtick.labelsize':8.5, 'ytick.labelsize':9,
                         'legend.fontsize':8.5}):
        fig = plt.figure(figsize=(156 / 25.4, 154 / 25.4))

        def axes(left, bottom, width, height):
            return fig.add_axes([left / 156, bottom / 154, width / 156, height / 154])

        a = axes(19, 92, 55, 50)
        h = axes(99, 92, 22, 50)
        pr = axes(126, 92, 29, 50)
        c = axes(19, 14, 55, 50)
        d = axes(99, 14, 43, 50)
        f1_panel(a)
        f2_panel([h, pr])
        f3_panel(c)
        f4_panel(d)
        c.tick_params(axis='y', pad=5)

        for x, y, label in [(1, 150, 'A  Binding'), (81, 150, 'B  Enzyme proximity'),
                            (1, 74, 'C  Function'), (81, 74, 'D  Associations')]:
            fig.text(x / 156, y / 154, label, weight='bold', fontsize=10)

        a.set_xticks([0, 20, 40, 60, 80])
        a.set_xlabel('CpG windows with a hit (%)', labelpad=3)
        h.set_xticks([0, 1.5])
        h.set_xlabel('Harmonic\ndifference', labelpad=3)
        pr.set_xticks([0, 2, 4])
        pr.set_xlabel('PageRank difference\n(× 10⁻⁴)', labelpad=3)
        d.set_xticks([0, 50, 100])
        d.set_xlabel('Supported / tested\npairs (%)', labelpad=3)
        d.text(1.14, 1.01, 'Testable /\nselected', transform=d.transAxes,
               fontsize=8.5, ha='center', va='bottom')

        # Put GO identifiers and counts together so they cannot be mistaken
        # for additional data rows; preserve the original six-term order.
        for artist in list(c.texts):
            artist.remove()
        c.set_xticklabels([f'{term}\n{count}' for term, count in
                           zip(terms.plot_id, matrix.sum(axis=0))])
        c.set_xlabel('GO term (genes out of 17)', labelpad=3, fontsize=8.5)

        legend_style = dict(frameon=False, fontsize=8.5, handlelength=.8,
                            handletextpad=.35, columnspacing=.8, borderaxespad=0)
        fig.legend(handles=f1handles, loc='center left', bbox_to_anchor=(1 / 156, 145 / 154),
                   ncol=3, **legend_style)
        fig.legend(handles=f2handles, loc='center left', bbox_to_anchor=(81 / 156, 145 / 154),
                   ncol=len(f2handles), **legend_style)
        fig.legend(handles=[Patch(facecolor=PURPLE, label='Annotated'),
                            Patch(facecolor='#F1EDF5', label='Not annotated')],
                   loc='center left', bbox_to_anchor=(19 / 156, 69 / 154),
                   ncol=2, **legend_style)
        readable_text(fig, minimum=8.5)
        # Fixed canvas keeps the font sizes when included at manuscript width.
        artifacts.save_figure(fig, OUT / 'manual_candidate_evidence.pdf', facecolor='white')
        plt.close(fig)


render_overview()

# Detailed F4 counts/directions at final manuscript size (approved 29 October 2026).
# Keep the fixed canvas: tight bounding boxes would change the embedded font size.
def render_f4_cohort_details():
    d = contexts.sort_values(['MF', 'cohort']).reset_index(drop=True)
    BLUE, BROWN, GREY, PALE, INK = '#0072B2','#B45C23','#D9DFE4','#F3F6F8','#25333D'
    mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'pdf.fonttype':42,'axes.linewidth':.6,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
    genes=sorted(d.MF.unique())
    def sequence(group):
     seq=[]
     for g in group:
      seq.append(('header',g,None))
      for _,r in d.loc[d.MF.eq(g)].iterrows():seq.append(('row',g,r))
     return seq
    sequences=[sequence(genes[:8]),sequence(genes[8:])]
    height=max(map(len,sequences)); short={'SKCM_METASTATIC':'SKCM-M','SKCM_PRIMARY':'SKCM-P'}
    fig=plt.figure(figsize=(156/25.4,166/25.4))
    text_artists=[]
    for left,seq in zip([.105,.605],sequences):
     ax=fig.add_axes([left,.15,.235,.72])
     labels=[];ticks=[]
     for y,(kind,g,r) in enumerate(seq):
      if kind=='header':
       text_artists.append(ax.text(-.03,y,g,transform=ax.get_yaxis_transform(),ha='right',va='center',weight='bold',fontsize=9.0))
       ax.axhspan(y-.43,y+.43,color=PALE,zorder=0,xmin=-.03,xmax=1.62,clip_on=False)
      else:
       labels.append(short.get(r.cohort,r.cohort));ticks.append(y)
       if r.n_tested_pairs>0:
        neg=100*r.n_a_supported_negative/r.n_tested_pairs;pos=100*r.n_a_supported_positive/r.n_tested_pairs
        unsupported=100*(r.n_tested_pairs-r.n_supported_pairs)/r.n_tested_pairs
        ax.barh(y,neg,height=.64,color=BLUE)
        ax.barh(y,pos,left=neg,height=.64,color=BROWN)
        ax.barh(y,unsupported,left=neg+pos,height=.64,color=GREY)
        text_artists.append(ax.text(104,y,f'{int(r.n_supported_pairs):,}/{int(r.n_tested_pairs):,}',va='center',fontsize=8.5))
       else:
        text_artists.append(ax.text(2,y,'Not testable',va='center',fontsize=8.5,color='#62727D'))
        text_artists.append(ax.text(104,y,'—',va='center',fontsize=8.5,color='#62727D'))
     ax.set(yticks=ticks,yticklabels=labels,ylim=(height-.2,-.8),xlim=(0,100),xticks=[0,50,100])
     ax.spines[['top','right','left']].set_visible(False)
     ax.tick_params(axis='y',length=0,pad=4)
     ax.tick_params(axis='x',length=3,width=.6,pad=3)
     text_artists.append(ax.text(104,-.9,'Supported /\ntested',va='bottom',fontsize=8.5,linespacing=1.0))
    fig.text(.04,.985,'F4 · Cohort-specific support and association direction\nfor all 17 candidates · Union cohort analyses',fontsize=10,weight='bold',va='top',linespacing=1.3)
    fig.text(.5,.092,'Share of tested MF–CpG pairs (%)',ha='center',fontsize=8.5)
    fig.legend(handles=[Patch(facecolor=BLUE,label='Supported: a < 0'),Patch(facecolor=BROWN,label='Supported: a > 0'),Patch(facecolor=GREY,label='Tested, not supported')],loc='lower center',bbox_to_anchor=(.5,.040),ncol=3,frameon=False,columnspacing=1.2,handlelength=1.1,handletextpad=.4)
    fig.text(.04,.010,f'{len(d)} selected contexts; {int(d.n_tested_pairs.eq(0).sum())} not testable. Original F4 results.',fontsize=8.5,va='bottom')
    artifacts.save_figure(fig, figure_output_dir(OUT, 'manual_candidate_associations') / 'manual_candidate_associations.pdf', facecolor='white')
    plt.close(fig)

with mpl.rc_context():
    render_f4_cohort_details()

print('Saved both manual-panel figures.')
