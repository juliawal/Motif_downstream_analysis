#!/usr/bin/env python3
"""Join the current-paper evidence using cohort, analysis and canonical Ensembl ID."""
import argparse
from pathlib import Path
import sys
sys.dont_write_bytecode = True

import pandas as pd

KEY = ['analysis', 'cohort', 'MF_ensembl']
ARMS = ['rank_1', 'rank_2', 'rank_3', 'union']
CASES = {('ERG', 'PRAD'), ('TCL1A', 'TGCT'), ('RUNX1', 'THCA'), ('HNF1B', 'UCEC')}


def join_contexts(left, right):
    columns = KEY + [col for col in right if col not in left]
    return left.merge(right[columns], on=KEY, how='left')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[2] / 'analysis_results')
    args = parser.parse_args()
    root = args.root
    e, f, out = root / 'binding_network_enrichment', root / 'associations', root / 'combined_evidence'
    selection = root / 'candidate_selection'
    out.mkdir(parents=True, exist_ok=True)
    read = lambda path: pd.read_csv(path, sep='\t')
    members = read(selection / 'selected_gene_cohort_memberships.tsv')
    binding = read(e / 'binding_measurements.tsv')
    ppi = read(e / 'network_individual_scores.tsv')
    f4 = read(f / 'association_gene_summary.tsv')
    f3 = read(e / 'enrichment_gene_contexts.tsv')
    f3counts = read(e / 'enrichment_test_counts.tsv')
    groups = read(e / 'network_group_tests.tsv')
    cross = join_contexts(members, binding)
    cross = join_contexts(cross, ppi)
    cross = join_contexts(cross, f4)
    cross = join_contexts(cross, f3)
    cross['n_significant_min3_terms'] = cross.n_significant_min3_terms.fillna(0).astype(int)
    cross['f3_themes'] = cross.f3_themes.fillna('')
    cross = cross.merge(f3counts[['analysis', 'cohort', 'n_tested_terms']].rename(columns={'n_tested_terms': 'n_f3_tests'}),
                        on=['analysis', 'cohort'], how='left')
    gp = groups.pivot(index=['analysis', 'cohort'], columns='measure', values='p_value').add_prefix('p_value_').reset_index()
    cross = cross.merge(gp, on=['analysis', 'cohort'], how='left')
    cross['f1_available'] = cross.status.eq('tested')
    cross['f1_positive'] = cross.f1_available & cross.difference.gt(0)
    cross['f2_available'] = cross.mapping_status.eq('in_network') & cross.harmonic_degree_bin_excess.notna()
    cross['f2_positive'] = cross.f2_available & cross.harmonic_degree_bin_excess.gt(0)
    cross['f2_harmonic_excess_positive'] = cross.f2_positive
    cross['f2_pagerank_excess_positive'] = cross.f2_available & cross.pagerank_degree_bin_excess.gt(0)
    cross['f3_available'] = cross.n_f3_tests.gt(0)
    cross['f3_member_significant_min3_term'] = cross.f3_available & cross.n_significant_min3_terms.gt(0)
    cross['f4_testable'] = cross.n_tested_pairs.gt(0)
    cross['f4_supported'] = cross.f4_testable & cross.n_supported_pairs.gt(0)
    cross['all_four_available'] = cross[['f1_available', 'f2_available', 'f3_available', 'f4_testable']].all(axis=1)
    other_positive = cross[['f2_positive', 'f3_member_significant_min3_term', 'f4_supported']].all(axis=1)
    cross['all_four_positive'] = cross.f1_positive & other_positive
    cross['descriptive_all_four_positive'] = cross.all_four_positive
    cross['all_four_positive_and_pagerank_excess'] = cross.all_four_positive & cross.f2_pagerank_excess_positive
    cross['f2_harmonic_group_nominal'] = cross.p_value_harmonic.le(.05)

    union = cross[cross.analysis.eq('union')].copy()
    settings_flags = ['f1_positive', 'f2_positive', 'f3_member_significant_min3_term', 'f4_supported', 'all_four_positive']
    for arm in ARMS[:3]:
        x = cross[cross.analysis.eq(arm)][['cohort', 'MF_ensembl'] + settings_flags + ['n_tested_pairs', 'n_supported_pairs']].copy()
        x['selected'] = True
        x = x.rename(columns={col: f'{arm}_{col}' for col in x if col not in ['cohort', 'MF_ensembl']})
        union = union.merge(x, on=['cohort', 'MF_ensembl'], how='left')
        for flag in ['selected'] + settings_flags:
            union[f'{arm}_{flag}'] = union[f'{arm}_{flag}'].eq(True)
    for outflag, suffix in [('selected_all_three_settings', 'selected'), ('f1_positive_all_three_settings', 'f1_positive'),
                           ('f2_positive_all_three_settings', 'f2_positive'), ('f3_membership_all_three_settings', 'f3_member_significant_min3_term'),
                           ('f4_support_all_three_settings', 'f4_supported')]:
        union[outflag] = union[[f'{a}_{suffix}' for a in ARMS[:3]]].all(axis=1)
    union['all_four_positive_in_union_and_all_settings'] = union.all_four_positive & union[[f'{a}_all_four_positive' for a in ARMS[:3]]].all(axis=1)

    positive = union[union.all_four_positive].sort_values(['cohort', 'MF'])
    for data, name in [(cross, 'combined_evidence_contexts.tsv'), (union, 'union_evidence_contexts.tsv'),
                       (positive, 'union_positive_contexts.tsv')]:
        data.to_csv(out / name, sep='\t', index=False, na_rep='NA')
    cases = cross[[pair in CASES for pair in zip(cross.MF, cross.cohort)]]
    cases.to_csv(out / 'case_evaluation_criteria.tsv', sep='\t', index=False, na_rep='NA')
    print(f'Combined evidence: {len(cross)} contexts, {len(positive)} positive union contexts.')


if __name__ == '__main__':
    main()
