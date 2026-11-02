#!/usr/bin/env python3
"""Prepare the paper's additional-candidate tables from saved evidence.

Membership comes from paper_cohorts_and_genes.json. No screening thresholds,
priority tags, model fitting or statistical tests are applied here.
"""
import argparse
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--analysis-dir', type=Path, default=ROOT / 'analysis_results')
    args = parser.parse_args()
    analysis = args.analysis_dir.resolve()
    config = json.loads((ROOT / 'scripts/paper_cohorts_and_genes.json').read_text())
    candidates = pd.DataFrame({'MF': config['additional_candidates']})
    if candidates.MF.duplicated().any():
        raise ValueError('Duplicate additional candidates in the paper configuration.')

    evidence = analysis / 'candidate_evidence'
    genes = pd.read_csv(evidence / 'all_candidate_evidence.tsv', sep='\t')
    shortlist = candidates.merge(genes, on='MF', validate='one_to_one').sort_values('MF')
    if len(shortlist) != len(candidates):
        raise ValueError('Missing saved evidence for an additional candidate.')
    contexts = pd.read_csv(evidence / 'union_context_evidence.tsv', sep='\t')
    contexts = contexts.merge(candidates, on='MF', validate='many_to_one')
    contexts = contexts.sort_values(['MF', 'cohort'])
    if set(contexts.MF) != set(candidates.MF):
        raise ValueError('Missing cohort contexts for an additional candidate.')

    output = analysis / 'additional_candidates'
    output.mkdir(parents=True, exist_ok=True)
    for frame, name in [(shortlist, 'candidate_shortlist'), (contexts, 'candidate_contexts')]:
        frame.to_csv(output / f'{name}.tsv', sep='\t', index=False, na_rep='NA')
    print(f'Prepared {len(shortlist)} additional candidates and {len(contexts)} cohort contexts.')


if __name__ == '__main__':
    main()
