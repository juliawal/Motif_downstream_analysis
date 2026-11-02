#!/usr/bin/env python3
"""Summarize the saved F4 export for the paper.

One branch is read at a time, retaining at most four cohort pair tables and the
small union columns required for cross-cohort exact-CpG recurrence. Original
triplet p values and support flags are never refitted or rethresholded.
"""
from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path
import sys
sys.dont_write_bytecode = True

import numpy as np
import pandas as pd

ARMS = ['rank_1', 'rank_2', 'rank_3', 'union']
CASES = [('ERG', 'PRAD'), ('TCL1A', 'TGCT'), ('RUNX1', 'THCA'), ('HNF1B', 'UCEC')]
KEY = ['MF_ensembl', 'CpG']
COUNT_COLS = ['n_input_pairs', 'n_candidate_pairs', 'n_tested_pairs', 'n_supported_pairs',
              'n_candidate_triplets', 'n_tested_triplets', 'n_supported_triplets',
              'n_nominal_joint_triplets', 'n_no_candidate_pairs', 'n_only_untestable_pairs',
              'n_tested_negative_pairs', 'n_untestable_triplets', 'n_negative_triplets']


def ratio(n, d):
    return float(n / d) if d else np.nan


def median(values):
    """Preserve missing estimates without NumPy's empty-slice warning."""
    present = values.dropna()
    return present.median() if len(present) else np.nan


def add_fractions(frame):
    for name, numerator, denominator in [
        ('pair_support_fraction', 'n_supported_pairs', 'n_tested_pairs'),
        ('triplet_support_fraction', 'n_supported_triplets', 'n_tested_triplets'),
        ('supported_positive_a_fraction', 'n_a_supported_positive', 'n_supported_pairs'),
        ('supported_mean_b', 'supported_b_sum', 'n_supported_triplets'),
        ('supported_mean_abs_b', 'supported_b_abs_sum', 'n_supported_triplets'),
        ('supported_mean_abs_distance', 'supported_abs_distance_sum', 'n_supported_triplets'),
    ]:
        if numerator in frame:
            frame[name] = frame[numerator] / frame[denominator].replace(0, np.nan)
    # Explicit aliases match the preceding paper's descriptive tables.
    frame['supported_pair_fraction'] = frame.pair_support_fraction
    frame['supported_triplet_fraction'] = frame.triplet_support_fraction
    return frame


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path(__file__).resolve().parents[2] / 'pipeline_results')
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parents[2] / 'analysis_results/associations')
    args = parser.parse_args()
    root, out = args.input.resolve(), args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)

    def read(relative, **kwargs):
        path = root / relative
        return pd.read_csv(path, sep='\t', **kwargs)

    def read_json(relative):
        path = root / relative
        return json.loads(path.read_text())

    def save(frame, name):
        frame.to_csv(out / name, sep='\t', index=False, na_rep='NA',
                     compression={'method': 'gzip', 'mtime': 0} if name.endswith('.gz') else None)

    def branch(arm, cohort):
        base = 'original_pipeline_outputs/union' if arm == 'union' else f'original_pipeline_outputs/selected_settings/{arm}'
        return f'{base}/{cohort}'

    metas = [read_json(p.relative_to(root)) for p in sorted((root / 'f4_methylation_factor_cpg_pairs').glob('*.json'))]
    cohorts = sorted({m['cohort'] for m in metas})
    metadata = {(m['analysis'], m['cohort']): m for m in metas}
    ann = read('reference_annotations/gene_annotation_file/gene_annotation.v22.tsv.gz')
    symbols = ann.set_index('ensembl_id').selected_gene_symbol
    case_ids = {}
    for symbol in [c[0] for c in CASES] + ['TDRD1']:
        ids = ann.loc[ann.selected_gene_symbol.eq(symbol), 'ensembl_id']
        case_ids[symbol] = ids.iloc[0]
    cohort_rows, contexts = [], []
    pairwise_rows, pattern_rows, union_frames = [], [], []
    common_tested, common_support = Counter(), Counter()

    for cohort in cohorts:
        tables = {}
        for arm in ARMS:
            meta = metadata[arm, cohort]
            relative = f'f4_methylation_factor_cpg_pairs/{arm}_{cohort}.tsv.gz'
            z = read(relative)
            n_c, n_t, n_s = (z[c] for c in ['n_candidate_triplets', 'n_tested_triplets', 'n_supported_triplets'])
            source_base = branch(arm, cohort)
            # The JSON preserves the complete original summary as strings.
            original = pd.Series(meta['original_summary'])
            original['testable'] = original.testable == 'True'
            selected = read(f'{source_base}/selected_mfs/mfs.tsv.gz', usecols=['MF_ensembl'])
            selected['MF'] = selected.MF_ensembl.map(symbols)
            counts = dict(n_input_pairs=len(z), n_candidate_pairs=int(n_c.gt(0).sum()),
                          n_tested_pairs=int(n_t.gt(0).sum()), n_supported_pairs=int(n_s.gt(0).sum()),
                          n_candidate_triplets=int(n_c.sum()), n_tested_triplets=int(n_t.sum()),
                          n_supported_triplets=int(n_s.sum()), n_nominal_joint_triplets=int(z.n_nominal_joint_triplets.sum()),
                          n_no_candidate_pairs=int(n_c.eq(0).sum()), n_only_untestable_pairs=int((n_c.gt(0) & n_t.eq(0)).sum()),
                          n_tested_negative_pairs=int((n_t.gt(0) & n_s.eq(0)).sum()),
                          n_untestable_triplets=int((n_c - n_t).sum()), n_negative_triplets=int((n_t - n_s).sum()))
            row = dict(analysis=arm, cohort=cohort, n_samples=int(original.n_samples), testable=bool(original.testable),
                       n_selected_mfs=len(selected), n_input_mfs=z.MF_ensembl.nunique(),
                       n_candidate_mfs=z.loc[n_c.gt(0), 'MF_ensembl'].nunique(),
                       n_tested_mfs=z.loc[n_t.gt(0), 'MF_ensembl'].nunique(),
                       n_supported_mfs=z.loc[n_s.gt(0), 'MF_ensembl'].nunique(),
                       fdr_alpha=float(original.fdr_alpha), max_cpg_tss_distance_bp=int(original.max_cpg_tss_distance_bp),
                       age_and_sex_columns=original.age_and_sex_columns, source_summary=f'f4_methylation_factor_cpg_pairs/{arm}_{cohort}.json#original_summary',
                       source_pairs=relative, source_triplets=meta['source'], **counts)
            cohort_rows.append(row)
            z['candidate_pair'] = n_c.gt(0).astype(int)
            z['tested_pair'] = n_t.gt(0).astype(int)
            z['supported_pair'] = n_s.gt(0).astype(int)
            z['a_supported'] = z.a_min.where(n_s.gt(0))
            z['a_tested'] = z.a_min.where(n_t.gt(0))
            for population, mask in [('supported', n_s.gt(0)), ('tested', n_t.gt(0))]:
                for sign, comparison in [('positive', z.a_min.gt(0)), ('negative', z.a_min.lt(0)), ('zero', z.a_min.eq(0))]:
                    z[f'a_{population}_{sign}'] = (mask & comparison).astype(int)
            specifications = {
                'n_input_pairs': ('CpG', 'size'), 'n_candidate_pairs': ('candidate_pair', 'sum'),
                'n_tested_pairs': ('tested_pair', 'sum'), 'n_supported_pairs': ('supported_pair', 'sum'),
                'a_supported_median': ('a_supported', 'median'), 'a_tested_median': ('a_tested', 'median'),
                'a_supported_q25': ('a_supported', lambda v: v.quantile(.25)),
                'a_supported_q75': ('a_supported', lambda v: v.quantile(.75)),
                'triplet_q_min': ('joint_q_min', 'min')}
            for column in ['n_candidate_triplets', 'n_tested_triplets', 'n_supported_triplets', 'n_nominal_joint_triplets',
                           'n_supported_b_positive', 'n_supported_b_negative', 'supported_b_sum',
                           'supported_b_abs_sum', 'supported_abs_distance_sum']:
                specifications[column] = (column, 'sum')
            for population in ['tested', 'supported']:
                for sign in ['positive', 'negative', 'zero']:
                    specifications[f'n_a_{population}_{sign}'] = (f'a_{population}_{sign}', 'sum')
            ag = z.groupby('MF_ensembl').agg(**specifications).reset_index()
            ag = selected.merge(ag, on='MF_ensembl', how='left')
            for col in [c for c in ag if c.startswith('n_')]:
                ag[col] = ag[col].fillna(0).astype('int64')
            for col in ['supported_b_sum', 'supported_b_abs_sum', 'supported_abs_distance_sum']:
                ag[col] = ag[col].fillna(0.)
            ag['analysis'], ag['cohort'], ag['n_samples'], ag['cohort_testable'] = arm, cohort, int(original.n_samples), bool(original.testable)
            ag['n_no_candidate_pairs'] = ag.n_input_pairs - ag.n_candidate_pairs
            ag['n_only_untestable_pairs'] = ag.n_candidate_pairs - ag.n_tested_pairs
            ag['n_tested_negative_pairs'] = ag.n_tested_pairs - ag.n_supported_pairs
            ag['n_untestable_triplets'] = ag.n_candidate_triplets - ag.n_tested_triplets
            ag['n_negative_triplets'] = ag.n_tested_triplets - ag.n_supported_triplets
            ag['context_status'] = np.select([ag.n_supported_pairs.gt(0), ag.n_tested_pairs.gt(0), ag.n_candidate_pairs.gt(0), ag.n_input_pairs.gt(0)],
                                             ['tested_with_support', 'tested_without_support', 'candidate_not_testable',
                                              'no_candidate_triplets' if original.testable else 'cohort_below_min_samples'], default='no_input_pairs')
            contexts.append(add_fractions(ag))
            tables[arm] = z
            if arm == 'union':
                union_frames.append(z[['cohort', 'MF', *KEY, 'n_candidate_triplets', 'n_tested_triplets', 'n_supported_triplets', 'a_min']].copy())

        # Exact pair identities: bit order is S1, S2, S3, union. Missing input
        # cannot accidentally be interpreted as a tested negative result.
        wide = None
        for arm, z in tables.items():
            small = z.set_index(KEY)[['candidate_pair', 'tested_pair', 'supported_pair', 'a_min']].copy()
            small['input_pair'] = 1
            small.columns = [f'{arm}_{c}' for c in small]
            wide = small if wide is None else wide.join(small, how='outer')
        for population in ['input', 'candidate', 'tested', 'supported']:
            cols = [f'{arm}_{population}_pair' for arm in ARMS]
            wide[cols] = wide[cols].fillna(0).astype('int8')
            pattern = wide[cols].astype(str).agg(''.join, axis=1)
            wide[f'{population}_pattern'] = pattern
        all_test = wide.tested_pattern.eq('1111')
        all_support = wide.supported_pattern.eq('1111')
        common_tested.update(wide.index.get_level_values('MF_ensembl')[all_test])
        common_support.update(wide.index.get_level_values('MF_ensembl')[all_support])
        pats = wide.reset_index().groupby(['MF_ensembl', 'input_pattern', 'tested_pattern', 'supported_pattern']).size().reset_index(name='n_pairs')
        pats.insert(0, 'cohort', cohort); pats['MF'] = pats.MF_ensembl.map(symbols)
        pattern_rows.append(pats)
        for symbol, cc in CASES:
            if cohort == cc:
                casewide = wide.reset_index()
                casewide = casewide[casewide.MF_ensembl.eq(case_ids[symbol])].copy()
                casewide.insert(0, 'cohort', cohort); casewide.insert(1, 'MF', symbol)
                casewide = casewide.rename(columns={'supported_pattern': 'support_pattern'})
                if symbol == 'TCL1A':
                    save(casewide, 'tcl1a_pair_patterns.tsv')
        for a, b in itertools.combinations(ARMS, 2):
            ai, bi = (wide[f'{arm}_input_pair'].eq(1) for arm in [a, b])
            at, bt = (wide[f'{arm}_tested_pair'].eq(1) for arm in [a, b])
            ap, bp = (wide[f'{arm}_supported_pair'].eq(1) for arm in [a, b])
            shared = at & bt
            diff = (wide[f'{a}_a_min'] - wide[f'{b}_a_min']).abs()
            pairwise_rows.append(dict(cohort=cohort, analysis_a=a, analysis_b=b,
                n_input_a=int(ai.sum()), n_input_b=int(bi.sum()), n_input_shared=int((ai & bi).sum()),
                n_tested_a=int(at.sum()), n_tested_b=int(bt.sum()), n_tested_shared=int(shared.sum()),
                n_supported_a=int(ap.sum()), n_supported_b=int(bp.sum()), n_supported_both=int((ap & bp).sum()),
                n_supported_a_only_shared_tested=int((ap & ~bp & shared).sum()),
                n_supported_b_only_shared_tested=int((bp & ~ap & shared).sum()),
                n_supported_a_absent_input_b=int((ap & ~bi).sum()), n_supported_b_absent_input_a=int((bp & ~ai).sum()),
                n_supported_a_input_b_untested=int((ap & bi & ~bt).sum()), n_supported_b_input_a_untested=int((bp & ai & ~at).sum()),
                max_abs_a_difference_shared_tested=diff[shared].max(),
                supported_jaccard=ratio((ap & bp).sum(), (ap | bp).sum())))
        print(f'Summarized {cohort}', flush=True)

    cd = add_fractions(pd.DataFrame(cohort_rows))
    md = pd.concat(contexts, ignore_index=True)
    save(cd, 'association_cohort_summary.tsv'); save(md, 'association_gene_summary.tsv')
    save(pd.DataFrame(pairwise_rows), 'pairwise_setting_overlap.tsv')
    patterns = pd.concat(pattern_rows, ignore_index=True)
    save(patterns, 'all_four_pair_patterns_by_context.tsv')
    save(patterns.groupby(['input_pattern', 'tested_pattern', 'supported_pattern'], as_index=False).n_pairs.sum(), 'all_four_pair_patterns_global.tsv')
    robustness = []
    for ens, g in patterns.groupby('MF_ensembl'):
        common = g.tested_pattern.eq('1111')
        all3 = g.supported_pattern.str[:3].eq('111')
        us = g.supported_pattern.str[3].eq('1')
        row = dict(MF_ensembl=ens, MF=symbols[ens],
                   n_common_input_pair_cohort_units=int(g.loc[g.input_pattern.eq('1111'), 'n_pairs'].sum()),
                   n_common_tested_pair_cohort_units=int(g.loc[common, 'n_pairs'].sum()),
                   n_supported_all_three_settings_common_tested=int(g.loc[common & all3, 'n_pairs'].sum()),
                   n_supported_all_four_common_tested=int(g.loc[common & all3 & us, 'n_pairs'].sum()),
                   n_union_supported_without_all_three_common_tested=int(g.loc[common & ~all3 & us, 'n_pairs'].sum()),
                   n_all_three_supported_without_union_common_tested=int(g.loc[common & all3 & ~us, 'n_pairs'].sum()))
        for i, arm in enumerate(ARMS):
            row[f'n_supported_{arm}_common_tested'] = int(g.loc[common & g.supported_pattern.str[i].eq('1'), 'n_pairs'].sum())
        robustness.append(row)
    robust = pd.DataFrame(robustness)
    casecontexts = md[[ (mf, cohort) in CASES for mf, cohort in zip(md.MF, md.cohort) ]]
    save(casecontexts, 'case_association_summary.tsv')
    global_rows = []
    for arm, c in cd.groupby('analysis', sort=False):
        m = md[md.analysis.eq(arm)]
        row = dict(analysis=arm, n_cohorts=len(c), n_testable_cohorts=int(c.testable.sum()),
                   n_selected_contexts=len(m), n_selected_mfs_distinct=m.MF_ensembl.nunique(),
                   n_tested_contexts=int(m.n_tested_pairs.gt(0).sum()), n_supported_contexts=int(m.n_supported_pairs.gt(0).sum()),
                   n_tested_mfs_distinct=m.loc[m.n_tested_pairs.gt(0), 'MF_ensembl'].nunique(),
                   n_supported_mfs_distinct=m.loc[m.n_supported_pairs.gt(0), 'MF_ensembl'].nunique())
        row.update({col:int(c[col].sum()) for col in COUNT_COLS})
        global_rows.append(row)
    global_frame = add_fractions(pd.DataFrame(global_rows))
    save(global_frame, 'association_global_summary.tsv')
    recurrence = []
    for (arm, ens), g in md.groupby(['analysis', 'MF_ensembl']):
        tested, supported = g.n_tested_pairs.gt(0), g.n_supported_pairs.gt(0)
        row = dict(analysis=arm, MF_ensembl=ens, MF=symbols[ens], n_selected_cohorts=len(g),
                   n_testable_cohorts=int(tested.sum()), n_supported_cohorts=int(supported.sum()),
                   n_tested_no_support_cohorts=int((tested & ~supported).sum()), n_untestable_cohorts=int((~tested).sum()),
                   n_negative_median_a_supported_cohorts=int(g.a_supported_median.lt(0).sum()),
                   n_positive_median_a_supported_cohorts=int(g.a_supported_median.gt(0).sum()),
                   n_zero_median_a_supported_cohorts=int(g.a_supported_median.eq(0).sum()),
                   selected_cohorts='|'.join(sorted(g.cohort)), tested_cohorts='|'.join(sorted(g.loc[tested, 'cohort'])),
                   supported_cohorts='|'.join(sorted(g.loc[supported, 'cohort'])),
                   pair_support_fraction_min=g.pair_support_fraction.min(), pair_support_fraction_median=median(g.pair_support_fraction),
                   pair_support_fraction_max=g.pair_support_fraction.max(),
                   a_supported_median_across_contexts=median(g.a_supported_median))
        row.update({col:int(g[col].sum()) for col in COUNT_COLS})
        recurrence.append(row)
    rec = add_fractions(pd.DataFrame(recurrence))
    save(rec, 'gene_recurrence.tsv')
    u = pd.concat(union_frames, ignore_index=True)
    u['tested'] = u.n_tested_triplets.gt(0)
    u['supported'] = u.n_supported_triplets.gt(0)
    u['supported_positive_a'] = u.supported & u.a_min.gt(0)
    u['supported_negative_a'] = u.supported & u.a_min.lt(0)
    cp = u.groupby(KEY).agg(n_input_cohorts=('cohort', 'size'), n_tested_cohorts=('tested', 'sum'),
                            n_supported_cohorts=('supported', 'sum'),
                            n_supported_positive_a_cohorts=('supported_positive_a', 'sum'),
                            n_supported_negative_a_cohorts=('supported_negative_a', 'sum')).reset_index()
    cp['MF'] = cp.MF_ensembl.map(symbols)
    cp['mixed_supported_a_signs'] = cp.n_supported_positive_a_cohorts.gt(0) & cp.n_supported_negative_a_cohorts.gt(0)
    cp['supported_once_tested_multiple_cohorts'] = cp.n_supported_cohorts.eq(1) & cp.n_tested_cohorts.ge(2)
    cp['supported_once_tested_one_cohort'] = cp.n_supported_cohorts.eq(1) & cp.n_tested_cohorts.eq(1)
    save(cp, 'union_same_cpg_recurrence.tsv.gz')
    gene_cp = cp.groupby('MF_ensembl').agg(
        n_unique_input_cpgs=('CpG', 'size'), n_unique_tested_cpgs=('n_tested_cohorts', lambda v: int(v.gt(0).sum())),
        n_unique_supported_cpgs=('n_supported_cohorts', lambda v: int(v.gt(0).sum())),
        n_cpgs_tested_in_at_least_two_cohorts=('n_tested_cohorts', lambda v: int(v.ge(2).sum())),
        n_cpgs_supported_in_at_least_two_cohorts=('n_supported_cohorts', lambda v: int(v.ge(2).sum())),
        max_supported_cohort_recurrence=('n_supported_cohorts', 'max'),
        n_cpgs_supported_once_tested_multiple_cohorts=('supported_once_tested_multiple_cohorts', 'sum'),
        n_cpgs_supported_once_tested_one_cohort=('supported_once_tested_one_cohort', 'sum'),
        n_recurrent_supported_cpgs_mixed_a_signs=('mixed_supported_a_signs', 'sum')).reset_index()
    ur = rec[rec.analysis.eq('union')].merge(gene_cp, on='MF_ensembl', how='left')
    for c in gene_cp.columns.drop('MF_ensembl'):
        ur[c] = ur[c].fillna(0).astype(int)
    ur['n_tested_same_pair_cohort_units_union_and_all_three_settings'] = ur.MF_ensembl.map(common_tested).fillna(0).astype(int)
    ur['n_supported_same_pair_cohort_units_union_and_all_three_settings'] = ur.MF_ensembl.map(common_support).fillna(0).astype(int)
    ur = ur.merge(robust.drop(columns='MF'), on='MF_ensembl', how='left')
    save(ur, 'union_gene_recurrence.tsv')
    recurrence_pairs = []
    for ens, g in u.groupby('MF_ensembl', sort=False):
        by_cohort = {c:v.set_index('CpG') for c,v in g.groupby('cohort')}
        for a,b in itertools.combinations(sorted(by_cohort), 2):
            ga, gb = by_cohort[a], by_cohort[b]
            ta, tb = set(ga.index[ga.tested]), set(gb.index[gb.tested])
            sa, sb = set(ga.index[ga.supported]), set(gb.index[gb.supported])
            shared = sorted(sa & sb)
            signs_a, signs_b = np.sign(ga.loc[shared, 'a_min']), np.sign(gb.loc[shared, 'a_min'])
            recurrence_pairs.append(dict(MF_ensembl=ens, MF=symbols[ens], cohort_a=a, cohort_b=b,
                n_tested_a=len(ta), n_tested_b=len(tb), n_tested_both=len(ta & tb), n_supported_a=len(sa),
                n_supported_b=len(sb), n_supported_both=len(shared), supported_jaccard=ratio(len(shared), len(sa | sb)),
                n_same_a_sign=int(signs_a.eq(signs_b).sum()), n_opposite_a_sign=int((signs_a * signs_b).lt(0).sum())))
    save(pd.DataFrame(recurrence_pairs), 'union_pairwise_cohort_recurrence.tsv')

    # Retain binding annotations for the separate ERG--TDRD1 triplet table.
    binding = read('cpg_binding_annotations/target_cpg_hits.tsv.gz')
    bkeys = ['analysis', 'cohort', *KEY]
    tss_annotation = read('reference_annotations/gene_tss_annotations_file/gene_annot.tsv.gz')
    tdrd1_strands = tss_annotation.loc[tss_annotation.gene_id.eq(case_ids['TDRD1']), 'gene_strand'].dropna().unique()
    tdrd1_strand = tdrd1_strands[0]
    erg_frames = []
    for arm in ARMS:
        e = read(f'prad_erg_tdrd1_triplet_case_study/{arm}_PRAD_ERG_TDRD1_all_triplets.tsv.gz')
        e.insert(0, 'cohort', 'PRAD'); e.insert(0, 'analysis', arm)
        e = e.merge(binding[binding.track.eq('ChIP_seq')].drop(columns='MF'), on=bkeys, how='left')
        e['binding_record_status'] = e.f1_status.fillna('not_exported_for_this_case')
        # Saved F4 distance is genomic TSS minus CpG, without strand orientation.
        # Convert it explicitly before applying the promoter window convention.
        e['G_prime_strand'] = tdrd1_strand
        e['G_prime_oriented_cpg_distance_bp'] = e.G_prime_distance_bp * (-1 if tdrd1_strand == '+' else 1)
        e['promoter_proximal'] = e.G_prime_oriented_cpg_distance_bp.between(-1500, 199)
        erg_frames.append(e)
    save(pd.concat(erg_frames, ignore_index=True), 'erg_tdrd1_binding.tsv')
    tp, cel = cpg_figure_data(root / "f4_methylation_factor_cpg_pairs")
    save(tp, "tp63_cpg_recurrence.tsv")
    save(cel, "celsr3_cpg_recurrence.tsv")
    print(global_frame.to_string(index=False), flush=True)


def cpg_figure_data(pair_dir):

    def pair_source(cohort):
        path = pair_dir / f'union_{cohort}.tsv.gz'
        data = pd.read_csv(path, sep='\t', low_memory=False)
        return data

    tp = []
    for cohort in ['CESC', 'ESCA']:
        data = pair_source(cohort)
        data = data[data.MF_ensembl.eq('ENSG00000073282') & data.n_tested_triplets.gt(0)].copy()
        data = data[['MF_ensembl', 'CpG', 'a_min', 'n_tested_triplets', 'n_supported_triplets']]
        tp.append(data.rename(columns={x: x + '_' + cohort for x in data if x not in ['MF_ensembl', 'CpG']}))
    tp = tp[0].merge(tp[1], on=['MF_ensembl', 'CpG'])
    tp['support_count'] = tp.n_supported_triplets_CESC.gt(0).astype(int) + tp.n_supported_triplets_ESCA.gt(0).astype(int)
    tp['a_CESC'] = tp.a_min_CESC
    tp['a_ESCA'] = tp.a_min_ESCA

    cel = []
    for cohort in ['BLCA', 'BRCA', 'HNSC', 'LUAD', 'PRAD']:
        data = pair_source(cohort)
        data = data[data.MF_ensembl.eq('ENSG00000008300') & data.CpG.eq('cg22702772')].copy()
        cel.append(data)
    cel = pd.concat(cel, ignore_index=True)

    return tp, cel


if __name__ == '__main__':
    main()
