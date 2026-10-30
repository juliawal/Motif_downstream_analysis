#!/usr/bin/env python3
"""Draw the literature comparison from the saved literature comparison counts."""
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
PAPERS = [('LevyJurgenson2019', 'Levy-Jurgenson', '2019'),
          ('Detilleux2022', 'Detilleux', '2022'),
          ('Silva2019_ELMER', 'Silva / ELMER', '2019')]
CRITERIA = [('F1', 'F1\nBinding'), ('F2_harmonic', 'F2\nHarmonic'),
            ('F2_PageRank', 'F2\nPageRank'), ('F3', 'F3\nGene-set\nenrichment'),
            ('F4', 'F4\nMediation\nsupport')]
BLUE, GRAY, INK, EDGE = '#0072B2', '#D4DADF', '#243746', '#8C99A3'


def main():
    counts = pd.read_csv(ROOT / 'analysis_results/literature_comparison.tsv',
                         sep='\t', keep_default_na=False)
    out = ROOT / 'figures/main'
    out.mkdir(parents=True, exist_ok=True)
    stem = 'literature_comparison'
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':8,
                         'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none',
                         'hatch.linewidth':.55, 'text.color':INK, 'savefig.facecolor':'white'})
    fig = plt.figure(figsize=(156/25.4, 78/25.4))
    ax = fig.add_axes([0,0,1,1])
    ax.set(xlim=(0,156), ylim=(0,78))
    ax.axis('off')
    centers = [47,71,95,119,143]
    bar_width, bar_height = 20, 3.5
    ax.text(2,75,'Shared candidate\nMFs',ha='left',va='top',fontsize=8.5,weight='bold',linespacing=1.2)
    for x,(_,label) in zip(centers,CRITERIA):
        ax.text(x,75,label,ha='center',va='top',fontsize=8.5,weight='bold',linespacing=1.2)
    ax.plot([2,154],[63,63],color='#CBD3D9',lw=.6)
    for y,(paper,name,year) in zip([53,38,23],PAPERS):
        dd = counts[counts.paper.eq(paper)].set_index('criterion')
        n_shared = int(dd.shared_genes.iloc[0])
        ax.text(2,y+3,name,ha='left',va='center',weight='bold',fontsize=8.2)
        ax.text(2,y-2,f'{year} · {n_shared} genes',ha='left',va='center',fontsize=8,color='#61717E')
        for x,(criterion,_) in zip(centers,CRITERIA):
            row = dd.loc[criterion]
            ax.text(x,y+3,row.label,ha='center',va='center',fontsize=9,weight='bold')
            left, bottom = x-bar_width/2, y-3.8
            for fraction,color,hatch in [(row.bar_met_fraction,BLUE,None),
                                         (row.bar_not_met_fraction,GRAY,None),
                                         (row.bar_unavailable_fraction,'white','////')]:
                width = bar_width*fraction
                if width:
                    ax.add_patch(Rectangle((left,bottom),width,bar_height,
                                           facecolor=color,edgecolor=EDGE if hatch else color,
                                           hatch=hatch,linewidth=.35 if hatch else 0))
                left += width
            ax.add_patch(Rectangle((x-bar_width/2,bottom),bar_width,bar_height,
                                   fill=False,edgecolor=EDGE,linewidth=.35))
    fig.legend(handles=[Patch(facecolor=BLUE,label='Criterion met'),
                        Patch(facecolor=GRAY,label='Assessed, not met'),
                        Patch(facecolor='white',edgecolor=EDGE,hatch='////',label='Not assessable')],
               loc='lower center',bbox_to_anchor=(.5,.05),ncol=3,frameon=False,
               fontsize=8,handlelength=1.4,columnspacing=1.5)
    ax.text(78,2.8,'Labels: criterion met / assessable. Bars: all shared candidate MFs.',
            ha='center',va='center',fontsize=8,color='#61717E')
    for extension in ['pdf', 'svg', 'png']:
        fig.savefig(out / f'{stem}.{extension}', dpi=300)
    plt.close(fig)
    print(f'Literature comparison figure: {out / stem}')


if __name__ == '__main__':
    main()
