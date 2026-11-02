#!/usr/bin/env python3
"""Render the approved PPI figure for Results section 2.5 from saved results.

Run from the project root:
    python scripts/figures/build_network_proximity.py

Defaults read analysis_results/ and write
network_proximity.pdf to figures/. Plotted tables
are written to analysis_results/. Use --evidence-dir and --output-dir to override.
The original nominal p-values are retained. No models, network scores,
hypothesis tests or multiple-testing corrections are recomputed.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from project_layout import figure_data_dir
import argparse
import os

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from figure_typography import readable_text

ROOT = Path(__file__).resolve().parents[2]

BLUE='#0072B2'; ORANGE='#D55E00'; TEAL='#009E73'; INK='#243746'; GREY='#7A8791'; PALE='#F1F4F6'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':7,'axes.titlesize':8,
 'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6.5,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.5,
 'xtick.major.width':.5,'ytick.major.width':.5,'pdf.fonttype':42,'savefig.facecolor':'white','text.color':INK,
 'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})


def save(fig,stem,output_dir):
    # Explicit physical dimensions match the final manuscript text width.
    readable_text(fig)
    fig.savefig(output_dir / f'{stem}.pdf', dpi=300)
    plt.close(fig)

def title(fig,x,y,letter,text,fs=8,va='top'):
    fig.text(x,y,letter,weight='bold',fontsize=9,va=va)
    fig.text(x+.034,y,text,weight='bold',fontsize=fs,va=va)


def main():
    parser=argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--evidence-dir',type=Path,default=ROOT/'analysis_results',
                        help='Directory with the three saved PPI evidence tables')
    parser.add_argument('--output-dir',type=Path,default=ROOT/'figures/main',
                        help='Figure directory; plotted tables use analysis_results/ for project outputs, otherwise data/')
    args=parser.parse_args()
    SOURCE=args.evidence_dir.resolve();OUT=args.output_dir.resolve();DATA=figure_data_dir(OUT)
    DATA.mkdir(parents=True,exist_ok=True)
    source_names=['union_network_effects.tsv','network_positive_contributions.tsv','network_contribution_concentration.tsv']

    effects=pd.read_csv(SOURCE/source_names[0],sep='\t');effects=effects[effects.analysis.eq('union')].copy()
    shares=pd.read_csv(SOURCE/source_names[1],sep='\t')
    concentration=pd.read_csv(SOURCE/source_names[2],sep='\t')
    cohorts=sorted(effects.cohort.unique())
    prad=shares[shares.cohort.eq('PRAD')].copy()

    fig=plt.figure(figsize=(156/25.4,173/25.4))
    positions=[(.19,.235),(.48,.205),(.735,.245)]
    for (left,_),letter,heading in zip(positions,'ABC',['Harmonic','PageRank','Top three candidates']):
        title(fig,left,.982,letter,heading,fs=8)
    for measure,(x,_) in zip(['harmonic','pagerank'],positions[:2]):
        n=int(effects.loc[effects.measure.eq(measure),'p_value'].le(.05).sum())
        fig.text(x,.951,f'{n}/{len(cohorts)} cohorts significant',fontsize=6.2,color=GREY,ha='left')
    fig.legend(handles=[Line2D([],[],marker='o',lw=0,color=BLUE,markersize=3,label='Harmonic'),Line2D([],[],marker='s',lw=0,color=ORANGE,markersize=3,label='PageRank')],loc='lower left',bbox_to_anchor=(.728,.936),ncol=2,frameon=False,fontsize=6.2,handlelength=.6,columnspacing=.7,handletextpad=.4)
    axs=[]
    labels=[c.replace('SKCM_METASTATIC','SKCM-MET').replace('SKCM_PRIMARY','SKCM-PRI') for c in cohorts]
    for j,(left,width) in enumerate(positions):
        ax=fig.add_axes([left,.33,width,.595]);axs.append(ax)
        ax.set_ylim(-.7,len(cohorts)-.3);ax.invert_yaxis();ax.set_axisbelow(True)
        for i,co in enumerate(cohorts):
            if i%2==0:ax.axhspan(i-.5,i+.5,color=PALE,zorder=0)
        ax.grid(axis='x',color='#D9E1E6',lw=.45)
        ax.spines['left'].set_visible(False)
        ax.set_yticks(range(len(cohorts)),labels=labels if j==0 else ['']*len(cohorts));ax.tick_params(axis='y',length=0,pad=3,labelsize=6.6)
        if j<2:
            measure=['harmonic','pagerank'][j];color=[BLUE,ORANGE][j]
            dd=effects[effects.measure.eq(measure)].set_index('cohort').loc[cohorts]
            ax.set_xscale('log')
            ax.axvline(.05,lw=.7,color=GREY,ls='--')
            for i,r in enumerate(dd.itertuples()):
                ax.plot([.05,r.p_value],[i,i],color=color,alpha=.25,lw=.8)
                ax.scatter(r.p_value,i,s=12,facecolors=color if r.p_value<=.05 else 'white',edgecolors=color,lw=.8,zorder=3)
            ax.set_xlim(.0001,1.15)
            ax.set_xticks([.0001,.005,.05,1],labels=['0.0001','0.005','0.05','1'])
            ax.minorticks_off()
            ax.set_xlabel('p-value (log scale)',labelpad=2,fontsize=6.5)
        else:
            cc=concentration.pivot(index='cohort',columns='measure',values='top3_positive_mass_fraction').loc[cohorts]*100
            for measure,color in [('harmonic',BLUE),('pagerank',ORANGE)]:
                ax.axvline(cc[measure].median(),color=color,ls='--',lw=.8,zorder=2.5)
            for i,r in enumerate(cc.itertuples()):
                ax.plot([r.harmonic,r.pagerank],[i,i],color='#B5C0C8',lw=.8)
            ax.scatter(cc.harmonic,np.arange(len(cohorts)),c=BLUE,s=10,zorder=3)
            ax.scatter(cc.pagerank,np.arange(len(cohorts)),c=ORANGE,s=10,zorder=3,marker='s')
            ax.set_xlim(0,103);ax.set_xticks([0,50,100]);ax.set_xlabel('Share of the sum of positive\ndifferences contributed by\nthe top three candidates (%)',labelpad=2,fontsize=6.5)
            fig.text(.735,.225,'Medians: ',fontsize=6.4,color=GREY)
            fig.text(.823,.225,f'{cc.harmonic.median():.1f}%',fontsize=6.4,color=BLUE)
            fig.text(.900,.225,f'{cc.pagerank.median():.1f}%',fontsize=6.4,color=ORANGE)
    fig.legend(handles=[Line2D([],[],marker='o',lw=0,color=INK,markersize=3,label='p ≤ 0.05'),Line2D([],[],marker='o',lw=0,color=INK,markerfacecolor='white',markersize=3,label='p > 0.05')],
     loc='lower left',bbox_to_anchor=(.025,.239),frameon=False,ncol=2,fontsize=6,handletextpad=.4,columnspacing=.8)

    title(fig,.025,.210,'D','PRAD: contributions by gene',fs=8,va='baseline')
    ax=fig.add_axes([.19,.115,.66,.070])
    genes=['WT1','ERG','SYCP2L','PAX8'];colors=['#7563A8',TEAL,'#56B4E9','#E69F00','#D4DCE2']
    bar_rows=[]
    for i,measure in enumerate(['harmonic','pagerank']):
        dd=prad[prad.measure.eq(measure)].set_index('MF')
        vals=[float(dd.loc[g,'positive_mass_fraction']*100) for g in genes]
        vals.append(100-sum(vals))
        start=0
        for gene,value,color in zip(genes+['Others'],vals,colors):
            ax.barh(i,value,left=start,height=.56,color=color,edgecolor='white',lw=.5)
            if value>=6:ax.text(start+value/2,i,f'{value:.1f}',ha='center',va='center',fontsize=6.1,color='white' if gene in ['WT1','ERG'] else INK)
            bar_rows.append(dict(measure=measure,gene=gene,share_percent=value))
            start+=value
        cc=concentration[(concentration.cohort=='PRAD')&(concentration.measure==measure)].iloc[0]
        ax.text(103,i,f'Top three\n{cc.top3_positive_mass_fraction*100:.1f}%',va='center',fontsize=6.3)
    ax.set_xlim(0,100);ax.set_ylim(-.5,1.5);ax.invert_yaxis()
    ax.set_yticks([0,1],labels=['Harmonic','PageRank']);ax.tick_params(axis='y',length=0,labelsize=6.3)
    ax.set_xticks([0,25,50,75,100]);ax.set_xlabel('Gene contribution (%)',labelpad=2,fontsize=6.4)
    ax.spines['left'].set_visible(False)
    fig.legend(handles=[Patch(color=c,label=g) for g,c in zip(genes+['Others'],colors)],loc='lower left',bbox_to_anchor=(.175,.012),ncol=5,frameon=False,fontsize=6.2,handlelength=.9,columnspacing=1.05)
    save(fig,'network_proximity',OUT)
    pd.DataFrame(bar_rows).to_csv(DATA / 'prad_network_contribution_shares.tsv', sep='\t', index=False)

    print('Finished figure and data export.')


if __name__=='__main__':
    main()
