"""Shared two-column F4 layout for supplementary Figures S4 and S5."""
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def draw_cohort_f4(data, *, height_mm=176):
    d = data.sort_values(['MF', 'cohort']).reset_index(drop=True)
    assert not d.duplicated(['MF', 'cohort']).any(), 'Duplicate current-panel contexts'
    BLUE, BROWN, GREY, PALE, INK = '#0072B2','#B45C23','#D9DFE4','#F3F6F8','#25333D'
    mpl.rcParams.update({'font.family':'DejaVu Sans','font.size':8.5,'axes.labelsize':8.5,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'legend.fontsize':8.5,'pdf.fonttype':42,'axes.linewidth':.6,'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})
    genes=sorted(d.MF.unique())
    def sequence(group):
     seq=[]
     for g in group:
      seq.append(('header',g,None))
      for _,r in d.loc[d.MF.eq(g)].iterrows():seq.append(('row',g,r))
     return seq
    split = min(range(1, len(genes)), key=lambda i: abs(len(sequence(genes[:i]))-len(sequence(genes[i:]))))
    sequences=[sequence(genes[:split]),sequence(genes[split:])]
    height=max(map(len,sequences)); short={'SKCM_METASTATIC':'SKCM-MET','SKCM_PRIMARY':'SKCM-PRI'}
    label_padding = .015 if d.cohort.isin(short).any() else 0
    fig=plt.figure(figsize=(156/25.4,height_mm/25.4))
    text_artists=[]
    for left,seq in zip([.105,.605],sequences):
     ax=fig.add_axes([left+label_padding,.15,.235-label_padding,.72])
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
    fig.text(.04,.985,f'F4 · Cohort-specific support and association direction\nfor all {len(genes)} candidates · Union cohort analyses',fontsize=10,weight='bold',va='top',linespacing=1.3)
    fig.text(.5,.092,'Share of tested MF–CpG pairs (%)',ha='center',fontsize=8.5)
    fig.legend(handles=[Patch(facecolor=BLUE,label='Supported: a < 0'),Patch(facecolor=BROWN,label='Supported: a > 0'),Patch(facecolor=GREY,label='Tested, not supported')],loc='lower center',bbox_to_anchor=(.5,.040),ncol=3,frameon=False,columnspacing=1.2,handlelength=1.1,handletextpad=.4)
    fig.text(.04,.010,f'{len(d)} MF–cohort combinations; {int(d.n_tested_pairs.eq(0).sum())} not testable. Original F4 results.',fontsize=8.5,va='bottom')
    return fig

