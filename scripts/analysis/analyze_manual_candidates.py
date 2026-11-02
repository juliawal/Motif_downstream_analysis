#!/usr/bin/env python3
"""Summarise the saved F1–F3 evaluations of the two current manuscript groups.

All p/q-values come from the server. Cohort profiles and F4 support come from
the original union results; no evaluation or statistical test is rerun here.
"""
import argparse
import json
import shutil
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def read(path):
    return pd.read_csv(path if path.is_file() else path.with_suffix(path.suffix + '.gz'), sep='\t')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'pipeline_results')
    parser.add_argument('--analysis', type=Path, default=ROOT / 'analysis_results')
    args = parser.parse_args()
    config = json.loads((ROOT / 'scripts/paper_cohorts_and_genes.json').read_text())
    union = read(args.analysis / 'combined_evidence/union_evidence_contexts.tsv')
    summaries = []
    for prefix, genes_key, section in [('manual', 'manual_prioritised', '2.9'),
                                       ('additional', 'additional_candidates', '2.10')]:
        genes = config[genes_key]
        raw = args.input / 'panel_evaluations' / config[prefix + '_evaluation'] / 'results/manual_mfs'
        selected = read(raw / 'selected_mfs/mfs.tsv')
        binding = read(raw / 'evaluation/binding_sites_per_mf.tsv').sort_values('MF')
        for frame in (selected, binding):
            if not frame.MF.is_unique or set(frame.MF) != set(genes):
                raise ValueError(f'Section {section}: saved panel membership differs from {genes_key}.')
        f1 = read(raw / 'evaluation/binding_sites.tsv')
        f2 = read(raw / 'evaluation/ppi_proximity.tsv')
        scores = read(raw / 'evaluation/ppi_proximity_per_mf.tsv').merge(
            selected, on='MF_ensembl', validate='many_to_one')
        enrichment = read(raw / 'evaluation/gene_set_enrichment.tsv')
        terms = enrichment.loc[enrichment.adj_p_value.le(.05) &
                               enrichment.term_has_min_mf_overlap].sort_values('adj_p_value', kind='stable')
        binding.loc[binding.status.ne('tested'), ['target_hit_rate', 'bg_hit_rate']] = np.nan
        binding['difference_percentage_points'] = 100 * (binding.target_hit_rate - binding.bg_hit_rate)
        memberships = terms[['library_name', 'term_name', 'adj_p_value', 'mf_symbols_in_term']].copy()
        memberships['MF'] = memberships.pop('mf_symbols_in_term').str.split(';')
        memberships = memberships.explode('MF').sort_values(['MF', 'adj_p_value'])
        binding['n_reportable_F3_terms'] = binding.MF.map(memberships.MF.value_counts()).fillna(0).astype(int)

        folder = args.analysis / (prefix + '_candidates')
        folder.mkdir(parents=True, exist_ok=True)
        # Preserve original statistical output bytes, including decimal precision.
        for source, name in [('binding_sites.tsv', 'binding_group_tests'),
                             ('ppi_proximity.tsv', 'network_group_tests'),
                             ('gene_set_enrichment.tsv', 'enrichment_all_tests')]:
            shutil.copyfile(raw / 'evaluation' / source, folder / f'{prefix}_{name}.tsv')
        frames = {'binding_measurements': binding, 'pooled_network_scores': scores,
                  'enrichment_terms': terms, 'enrichment_memberships': memberships}
        for name, frame in frames.items():
            frame.to_csv(folder / f'{prefix}_{name}.tsv', sep='\t', index=False)

        summary = {'section': section, 'n_genes': len(genes), 'genes': ';'.join(genes),
                   'saved_evaluation': config[prefix + '_evaluation']}
        for track, label in [('FIMO', 'motif'), ('ChIP_seq', 'chip')]:
            row = f1.set_index('track').loc[track]
            for field in ['n_mf', 'n_positive_differences', 'n_nonzero_differences', 'p_value_sign_test']:
                summary[f'f1_{label}_{field}'] = row[field]
        for measure in ['harmonic', 'pagerank']:
            row = f2.set_index('measure').loc[measure]
            summary[f'f2_{measure}_n_mapped'] = row.n_mfs_in_network
            summary[f'f2_{measure}_p_value'] = row.p_value
        summary['f3_n_reportable_terms'] = len(terms)
        summary['f3_n_go_terms'] = int(terms.library_name.str.startswith('GO_').sum())
        summary['f3_n_reactome_terms'] = int(terms.library_name.str.startswith('Reactome').sum())
        summaries.append(summary)

        if prefix == 'manual':
            panel = pd.DataFrame({'MF': genes})
            contexts = panel.merge(union, on='MF', validate='one_to_many')
            if set(contexts.MF) != set(genes):
                raise ValueError('Missing original cohort results for the current biological panel.')
            contexts = contexts.sort_values(['MF', 'cohort']).reset_index(drop=True)
            contexts.to_csv(folder / 'manual_prioritised_contexts.tsv', sep='\t', index=False)
            panel.to_csv(folder / 'manual_panel_membership.tsv', sep='\t', index=False)
            for measure in ['harmonic', 'pagerank']:
                excess = contexts[measure + '_degree_bin_excess']
                contexts[measure + '_difference_sign'] = np.select(
                    [excess.gt(0), excess.lt(0)], ['positive', 'negative'], default='zero')
            f2cols = ['MF', 'MF_ensembl', 'cohort', 'entrez_id'] + [
                measure + suffix for measure in ['pagerank', 'harmonic']
                for suffix in ['_score', '_degree_bin_null_mean', '_degree_bin_excess', '_difference_sign']]
            contexts[f2cols].to_csv(folder / 'manual_network_contexts.tsv', sep='\t', index=False)
        print(f'Section {section}: {len(genes)} genes; {len(terms)} saved reportable F3 terms.')
    pd.DataFrame(summaries).to_csv(args.analysis / 'section_evaluation_summary.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
