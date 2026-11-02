#!/usr/bin/env python3
"""Render the F3 manuscript figure from the prepared enrichment tables.

Run in the Motif environment:
    python scripts/figures/build_functional_enrichment.py

Reads the term table, cohort ranking and category summary in analysis_results/.
Writes functional_enrichment.pdf and the plotted gene memberships.

Use --evidence-dir and --output-dir for alternative locations.
Saved tests and memberships are retained; no new inference is run.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from project_layout import figure_data_dir
import argparse
import os

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Patch

from figure_typography import readable_text
from enrichment_data import CATEGORIES, ROOT, short_cohort

INK = '#243746'
GREY = '#7A8791'
PALE = '#F1F4F6'


def figure(terms, rank, summary, output_dir):
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':7, 'pdf.fonttype':42,
        'axes.linewidth':.5, 'savefig.facecolor':'white',
        'text.color':INK, 'axes.labelcolor':INK, 'xtick.color':INK, 'ytick.color':INK})
    fig = plt.figure(figsize=(156/25.4, 177/25.4))
    def title(y, letter, text):
        fig.text(.025,y,letter,fontsize=9,weight='bold',va='top')
        fig.text(.059,y,text,fontsize=8,weight='bold',va='top')
    title(.976,'A',f'Functional enrichment across {len(rank)} cohorts')
    fig.text(.025,.938,'Union · all cohorts with enriched terms',fontsize=6.8,color=GREY)
    fig.text(.025,.902,'Cohort',fontsize=6.4,color=GREY)
    fig.text(.920,.921,'Lowest\nadjusted p-value\nover all terms',fontsize=8,color=GREY,ha='center',va='top')
    cohorts=rank.cohort.tolist()
    for k,(category,(label,header,color)) in enumerate(CATEGORIES.items()):
        left=.205+k*.168
        ax=fig.add_axes([left,.49,.145,.37])
        ax.set_xscale('log');ax.set_xlim(1,160);ax.set_ylim(len(rank)-.4,-.6)
        ax.set_xticks([1,10,100],labels=['1','10','100']);ax.minorticks_off()
        ax.tick_params(axis='x',labelsize=6,length=2,width=.5,pad=2)
        ax.set_yticks([])
        for spine in ['top','left','right']:ax.spines[spine].set_visible(False)
        ax.grid(axis='x',color='#DDE4E8',lw=.5,zorder=0)
        display_header = ['DNA binding /\ntranscription', 'Development\n/ cell identity',
                          'Immune\nsignalling /\ncell state', 'Other\nfunctions'][k]
        fig.text(left+.0725,.920,display_header,color=color,fontsize=8,ha='center',va='top',weight='bold')
        for i,cohort in enumerate(cohorts):
            if i%2==0:ax.axhspan(i-.5,i+.5,color=PALE,zorder=-1)
            row=summary[(summary.cohort==cohort)&(summary.display_category==category)]
            if row.empty:
                ax.text(10,i,'–',ha='center',va='center',color='#B3BEC6',fontsize=7)
            else:
                r=row.iloc[0]
                ax.plot([r.minimum_fold_enrichment,r.maximum_fold_enrichment],[i,i],color=color,lw=1.2,zorder=2)
                ax.scatter(r.median_fold_enrichment,i,s=17,color=color,zorder=3)
                ax.annotate(str(int(r.n_terms)),(r.median_fold_enrichment,i),xytext=(0,5),textcoords='offset points',ha='center',fontsize=5.8,color=color)
            if k==0:
                yy=.49+.37*(len(rank)-.4-i)/(len(rank)+.2)
                fig.text(.025,yy,short_cohort(cohort),fontsize=6.5,va='center',weight='normal')
                value=rank.loc[rank.cohort.eq(cohort),'minimum_adjusted_p_value'].iloc[0]
                fig.text(.920,yy,f'{value:.2g}',fontsize=6.3,ha='center',va='center',color=GREY)
    fig.text(.53,.454,'Fold enrichment (log scale)',ha='center',fontsize=7)
    fig.text(.025,.431,'Point: median  ·  Line: minimum–maximum  ·  Number of enriched terms',fontsize=6.0,color=GREY)
    membership=[]
    def matrix(cohort, specific_id, letter, top, bottom):
        broad=terms[(terms.cohort==cohort)&(terms.term_id=='GO:0000977')].iloc[0]
        specific=terms[(terms.cohort==cohort)&(terms.term_id==specific_id)].iloc[0]
        a=set(broad.mf_symbols_in_term.split(';'));b=set(specific.mf_symbols_in_term.split(';'))
        order=sorted(b)+sorted(a-b)
        title(top,letter,f'{cohort}: shared genes between terms')
        cell_width=.76/21
        fig.text(.205+cell_width*len(order),top-.025,f'{len(b)} shared genes',fontsize=6.8,ha='right',color=GREY)
        ax=fig.add_axes([.205,bottom,cell_width*len(order),.04])
        ax.set_xlim(-.5,len(order)-.5);ax.set_ylim(1.5,-.5)
        ax.set_yticks([0,1],labels=[f'DNA binding ({len(a)}/{int(broad.n_mf_symbols)})',f'{"MET" if cohort=="SARC" else "Beta-cell dev."} ({len(b)}/{int(specific.n_mf_symbols)})'])
        ax.tick_params(axis='y',length=0,labelsize=6.4,pad=4)
        ax.set_xticks(range(len(order)),labels=order);ax.tick_params(axis='x',length=0,pad=2,labelsize=6.4)
        plt.setp(ax.get_xticklabels(),rotation=60,ha='right',rotation_mode='anchor')
        for spine in ax.spines.values():spine.set_visible(False)
        for j,gene in enumerate(order):
            for k,(r,genes) in enumerate([(broad,a),(specific,b)]):
                yes=gene in genes
                color=['#0072B2','#D55E00'][k] if yes else '#E8EDF0'
                ax.add_patch(Rectangle((j-.46,k-.43),.92,.86,facecolor=color,edgecolor='white',lw=.4))
                membership.append(dict(cohort=cohort,term_id=r.term_id,gene=gene,annotated_in_term=yes,n_candidates=int(r.n_mf_symbols)))
        ax.axvline(len(b)-.5,color='white',lw=2)
    matrix('SARC','GO:0060231','B',.390,.320)
    matrix('UCEC','R-HSA-186712','C',.230,.160)
    handles=[Patch(color=v[2],label=label) for v,label in zip(CATEGORIES.values(),
        ['DNA binding / transcription','Development / cell identity','Immune signalling / cell state','Other functions'])]
    handles.append(Patch(color='#E8EDF0',label='B–C: not annotated in this term'))
    fig.legend(handles=handles,loc='lower left',bbox_to_anchor=(.025,.003),ncol=2,frameon=False,fontsize=6.1,handlelength=.9,columnspacing=1.3)
    readable_text(fig)
    fig.savefig(output_dir / 'functional_enrichment.pdf', dpi=300)
    plt.close(fig)
    membership = pd.DataFrame(membership)
    return membership


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--evidence-dir', type=Path, default=ROOT / 'analysis_results',
                        help='Directory containing the prepared enrichment term, ranking and summary tables')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'figures/main',
                        help='Directory for PDF; accompanying tables use analysis_results/ for project outputs, otherwise data/')
    args = parser.parse_args()
    source = args.evidence_dir.resolve()
    terms, rank, summary = [
        pd.read_csv(source / f'{name}.tsv', sep='\t', float_precision='round_trip')
        for name in ['union_enrichment_terms', 'enrichment_cohort_ranking', 'enrichment_category_summary']
    ]
    output_dir = args.output_dir.resolve()
    data_dir = figure_data_dir(output_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    membership = figure(terms, rank, summary, output_dir)
    membership.to_csv(data_dir / 'enrichment_gene_memberships.tsv', sep='\t', index=False)
    print('Finished figure and data export.')


if __name__ == '__main__':
    main()
