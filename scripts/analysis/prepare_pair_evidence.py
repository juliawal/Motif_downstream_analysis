#!/usr/bin/env python3
"""Descriptive, reproducible prioritisation from saved F0/F1--F4 outputs.

No models or original hypothesis tests are refitted. Union counts each pair once
per cohort; the three settings are sensitivity analyses on the same patients.
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from prepare_direction_results import classify_pair_directions

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--pipeline-dir', type=Path, default=ROOT / 'pipeline_results')
parser.add_argument('--analysis-dir', type=Path, default=ROOT / 'analysis_results')
args = parser.parse_args()
PIPELINE, ANALYSIS = args.pipeline_dir.resolve(), args.analysis_dir.resolve()
DIRECTIONS = ANALYSIS / 'pair_directions'
CANDIDATES = ANALYSIS / 'candidate_evidence'
for directory in [DIRECTIONS, CANDIDATES]:
    directory.mkdir(parents=True, exist_ok=True)
K = ['MF_ensembl', 'CpG']


def read(p, **kwargs):
    return pd.read_csv(p, sep='\t', **kwargs)


def write(d, path):
    d.to_csv(path, sep='\t', index=False, na_rep='NA',
             compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None)


u = read(ANALYSIS / 'combined_evidence/union_evidence_contexts.tsv')
pairdir = PIPELINE / 'f4_methylation_factor_cpg_pairs'
tables = []
pair_paths = sorted(pairdir.glob('union_*.tsv.gz'))
for n, p in enumerate(pair_paths, 1):
    d = read(p)
    cohort = d.cohort.iloc[0]
    d['tested'] = d.n_tested_triplets.gt(0)
    d['supported'] = d.n_supported_triplets.gt(0)
    d['a'] = d.a_min
    d['negative'] = d.tested & d.a.lt(0)
    d['positive'] = d.tested & d.a.gt(0)
    d['supported_negative'] = d.supported & d.a.lt(0)
    d['supported_positive'] = d.supported & d.a.gt(0)
    for arm in ['rank_1','rank_2','rank_3']:
        x = read(pairdir/f'{arm}_{cohort}.tsv.gz', usecols=K+['n_tested_triplets','n_supported_triplets','a_min'])
        x['selected'] = True
        x['tested'] = x.n_tested_triplets.gt(0)
        x['supported'] = x.n_supported_triplets.gt(0)
        x = x[K+['selected','tested','supported','a_min']].rename(columns={c:f'{arm}_{c}' for c in ['selected','tested','supported','a_min']})
        d = d.merge(x, on=K, how='left')
        for c in ['selected','tested','supported']:
            d[f'{arm}_{c}'] = d[f'{arm}_{c}'].eq(True)
        w = read(PIPELINE / 'original_pipeline_outputs/selected_settings' / arm / cohort / 'evaluation_input/mf_cpg_pairs.tsv.gz', usecols=K+['weight'])
        d = d.merge(w.rename(columns={'weight':f'{arm}_weight'}),on=K,how='left')
        d = d.drop(columns=f'{arm}_a_min')
    d['selected_all_three'] = d[[f'{a}_selected' for a in ['rank_1','rank_2','rank_3']]].all(axis=1)
    d['supported_all_four'] = d.supported & d[[f'{a}_supported' for a in ['rank_1','rank_2','rank_3']]].all(axis=1)
    d['n_settings_selected'] = d[[f'{a}_selected' for a in ['rank_1','rank_2','rank_3']]].sum(axis=1)
    d['n_settings_supported'] = d[[f'{a}_supported' for a in ['rank_1','rank_2','rank_3']]].sum(axis=1)
    tables.append(d)
    print(f'{n}/{len(pair_paths)} {cohort}: {len(d):,} pairs',flush=True)
pairs = pd.concat(tables,ignore_index=True)

# Classify each exact pair once using every tested cohort, irrespective of F4 support.
rec = pairs.groupby(K,sort=False).agg(MF=('MF','first'),n_input_cohorts=('cohort','size'),
    n_tested_cohorts=('tested','sum'),n_supported_cohorts=('supported','sum'),
    n_negative=('negative','sum'),n_positive=('positive','sum'),
    n_supported_negative=('supported_negative','sum'),n_supported_positive=('supported_positive','sum'),
    n_selected_all_three=('selected_all_three','sum'),n_supported_all_four=('supported_all_four','sum'),
    a_min=('a','min'),a_max=('a','max'),a_median=('a','median')).reset_index()
rec['mixed_tested_signs'] = rec.n_negative.gt(0)&rec.n_positive.gt(0)
rec['direction_all_tested'] = classify_pair_directions(rec)
rec['same_sign_all_tested'] = rec.direction_all_tested.isin(['negative', 'positive'])
rec['negative_all_tested'] = rec.direction_all_tested.eq('negative')
rec['positive_all_tested'] = rec.direction_all_tested.eq('positive')
rec['min_abs_a_all_tested'] = np.where(rec.mixed_tested_signs,0,np.minimum(rec.a_min.abs(),rec.a_max.abs()))
write(rec,DIRECTIONS/'all_exact_pair_recurrence.tsv.gz')
write(rec[rec.same_sign_all_tested].sort_values(['n_tested_cohorts','min_abs_a_all_tested','MF_ensembl','CpG'],ascending=[False,False,True,True]),DIRECTIONS/'consistent_direction_pairs.tsv')
write(rec[rec.mixed_tested_signs],DIRECTIONS/'pairs_with_mixed_directions.tsv')

# Fully auditable context matrix for all selected genes, including missing evidence.
mc = pairs.groupby(['MF_ensembl','cohort']).agg(n_selected_all_three_pairs=('selected_all_three','sum'),n_supported_all_four_pairs=('supported_all_four','sum')).reset_index()
u = u.merge(mc,on=['MF_ensembl','cohort'])
write(u,CANDIDATES/'union_context_evidence.tsv')
mr = rec.groupby(['MF_ensembl','MF']).agg(n_unique_input_pairs=('CpG','size'),
    n_repeated_tested_pairs=('n_tested_cohorts',lambda x:int(x.ge(2).sum())),
    n_consistent_direction_pairs=('same_sign_all_tested','sum'),
    n_consistent_negative_pairs=('negative_all_tested','sum'),
    n_consistent_positive_pairs=('positive_all_tested','sum'),
    n_mixed_direction_pairs=('mixed_tested_signs','sum')).reset_index()
write(mr.sort_values(['MF','MF_ensembl']),DIRECTIONS/'pair_direction_counts_by_gene.tsv')
genes = u.groupby(['MF_ensembl','MF']).agg(n_selected_cohorts=('cohort','size'),
    n_tested_cohorts=('f4_testable','sum'),n_supported_cohorts=('f4_supported','sum'),
    n_f1_available=('f1_available','sum'),n_f1_positive=('f1_positive','sum'),
    n_f2_available=('f2_available','sum'),n_f2_positive=('f2_positive','sum'),
    n_f3_positive=('f3_member_significant_min3_term','sum'),n_all_four_positive=('all_four_positive','sum'),
    n_all_four_positive_all_settings=('all_four_positive_in_union_and_all_settings','sum'),
    n_selected_all_three_contexts=('selected_all_three_settings','sum'),
    n_tested_pairs=('n_tested_pairs','sum'),n_supported_pairs=('n_supported_pairs','sum'),
    n_supported_all_four_pairs=('n_supported_all_four_pairs','sum'),
    median_cohort_support_fraction=('supported_pair_fraction','median'),
    n_supported_negative=('n_a_supported_negative','sum'),n_supported_positive=('n_a_supported_positive','sum'),
    median_harmonic_excess=('harmonic_degree_bin_excess','median')).reset_index().merge(mr,on=['MF_ensembl','MF'],validate='one_to_one')
genes['pooled_support_fraction'] = genes.n_supported_pairs/genes.n_tested_pairs.replace(0,np.nan)
genes = genes.sort_values(['n_all_four_positive','n_supported_cohorts','n_consistent_direction_pairs','n_supported_pairs'],ascending=False)
write(genes,CANDIDATES/'all_candidate_evidence.tsv')

# Case-level F1 annotation is exported only for HNF1B/UCEC, ERG/PRAD, RUNX1/THCA.
binding = read(PIPELINE / 'cpg_binding_annotations/target_cpg_hits.tsv.gz')
binding = binding[binding.analysis.eq('union')]
binding = binding[['cohort']+K+['track','cpg_passed_filters','binding_hit']]
pairs = pairs.merge(binding,on=['cohort']+K,how='left')
pairs = pairs.merge(u[['MF_ensembl','cohort','n_samples','f1_available','f1_positive','f2_positive','f3_member_significant_min3_term','all_four_positive']],on=['MF_ensembl','cohort'])
write(pairs,DIRECTIONS/'all_pair_cohort_evidence.tsv.gz')
print(f'Prepared {len(rec):,} distinct pairs and evidence for {len(genes)} MFs.',flush=True)
