#!/usr/bin/env python3
"""Candidate coverage and descriptive sign baselines from the prepared pair results.

No new patient-level tests. Screening cutoffs are exploratory coverage rules,
not significance thresholds or validated candidate scores.
"""
import argparse
from pathlib import Path
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--analysis-dir', type=Path, default=ROOT/'analysis_results')
args = parser.parse_args()
DIRECTIONS=args.analysis_dir/'pair_directions'
CANDIDATES=args.analysis_dir/'additional_candidates'
read=lambda p:pd.read_csv(p,sep='\t',low_memory=False)
write=lambda d,p:d.to_csv(p,sep='\t',index=False,na_rep='NA')
d=read(DIRECTIONS/'all_pair_cohort_evidence.tsv.gz')
old=read(DIRECTIONS/'all_exact_pair_recurrence.tsv.gz')
allraw=d[['cohort','MF_ensembl','MF','CpG','tested','supported','a']].copy()
allraw[['tested','supported']]=allraw[['tested','supported']].astype(int)
allraw['sign']=np.where(allraw.tested,np.sign(allraw.a),0).astype(int)
t=allraw[allraw.tested.eq(1)].copy();t['positive']=t.sign.gt(0)
repeat=old[old.n_tested_cohorts.ge(2)]
x=t.merge(repeat[['MF_ensembl','CpG']],on=['MF_ensembl','CpG'])
# Descriptive reference used in the manuscript: signs among repeated pairs.
marginal=x.groupby(['MF_ensembl','cohort']).agg(p_positive=('positive','mean')).reset_index()
z=x.merge(marginal,on=['MF_ensembl','cohort'])
z['p_negative']=1-z.p_positive
e=z.groupby(['MF_ensembl','CpG']).agg(MF=('MF','first'),p_all_positive=('p_positive','prod'),p_all_negative=('p_negative','prod')).reset_index()
e['expected_same_sign']=e.p_all_positive+e.p_all_negative
e=e.merge(repeat[['MF_ensembl','CpG','same_sign_all_tested']],on=['MF_ensembl','CpG'])
write(e,DIRECTIONS/'sign_frequency_reference.tsv')

u=read(CANDIDATES/'union_context_evidence.tsv');g=read(CANDIDATES/'additional_candidate_evidence.tsv')
u=u[~u.manual17].copy()
u['n_favourable_layers']=u[['f1_positive','f2_positive','f3_member_significant_min3_term','f4_supported']].sum(axis=1)
u['screen_all4']=u.all_four_positive
u['screen_three_layers_20']=u.n_favourable_layers.ge(3)&u.n_supported_pairs.ge(20)
u['screen_fraction_stable']=u.n_tested_pairs.ge(50)&u.supported_pair_fraction.ge(.5)&u.n_supported_all_four_pairs.ge(10)
u['screen_large1000']=u.n_supported_pairs.ge(1000)
flags=['screen_all4','screen_three_layers_20','screen_fraction_stable','screen_large1000']
g=g.merge(u.groupby('MF_ensembl')[flags].any().reset_index(),on='MF_ensembl')
g['screen_recurrent5']=g.n_strict_recurrent_pairs.ge(5)
flags+=['screen_recurrent5']
g['screen_any']=g[flags].any(axis=1)
g['screen_core']=g.screen_all4|g.screen_three_layers_20|g.screen_recurrent5
g['flag_reasons']=g.apply(lambda z:';'.join(c for c in flags if z[c]),axis=1)
write(g,CANDIDATES/'candidate_coverage_flags.tsv')
counts=[]
for k in [1,2,5,10,25]:counts.append({'criterion':'strict_recurrent_pairs','threshold':k,'n_genes':int(g.n_strict_recurrent_pairs.ge(k).sum())})
for k in [5,10,20,50]:counts.append({'criterion':'supported_pairs_with_3_favourable_layers','threshold':k,'n_genes':int(u.loc[u.n_favourable_layers.ge(3)&u.n_supported_pairs.ge(k),'MF'].nunique())})
write(pd.DataFrame(counts),CANDIDATES/'screen_threshold_sensitivity.tsv')
# The upstream MF--CpG question does not require a supported downstream G' arm.
# Flag large standardized a values even where the original joint F4 rule fails.
# This is an effect-size diagnostic, not a new significance test or causal claim.
t['large_a']=t.a.abs().ge(.5)
t['large_a_without_F4_support']=t.large_a&t.supported.eq(0)
upstream=t.groupby(['MF_ensembl','MF','cohort']).agg(
    n_tested=('CpG','size'),n_large_a=('large_a','sum'),
    n_large_a_without_F4=('large_a_without_F4_support','sum'),
    median_abs_a=('a',lambda v:v.abs().median())).reset_index()
upstream=upstream[upstream.MF_ensembl.isin(g.MF_ensembl)].copy()
upstream['fraction_large_a_without_F4']=upstream.n_large_a_without_F4/upstream.n_tested
upstream['a_only_review_flag']=upstream.n_large_a_without_F4.ge(20)&upstream.fraction_large_a_without_F4.ge(.5)
write(upstream[upstream.a_only_review_flag].sort_values('n_large_a_without_F4',ascending=False),CANDIDATES/'expression_methylation_review_queue.tsv')

print(f'{len(repeat)} repeated pairs; {g.screen_core.sum()} core candidates; {g.screen_any.sum()} broadly screened candidates.')
