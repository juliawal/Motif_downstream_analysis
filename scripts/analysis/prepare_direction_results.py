#!/usr/bin/env python3
"""Summarise saved union MF--CpG coefficient signs, irrespective of F4 support.

This descriptive analysis uses the existing regression results. It does not
refit models, test significance, or compare candidates with background genes.
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
EXAMPLE_PAIRS = [
    ('POU2AF2', 'cg03078767'),
    ('ZNF256', 'cg00541956'),
    ('ZNF256', 'cg03430717'),
    ('LHX5', 'cg26843324'),
]
DIRECTIONS = ['negative', 'positive', 'mixed', 'zero']


def classify_pair_directions(pairs):
    """Classify signs across all tested cohorts; F4 support is not an input."""
    tested, negative, positive = pairs.n_tested_cohorts, pairs.n_negative, pairs.n_positive
    if (tested.lt(0) | negative.lt(0) | positive.lt(0) | (negative + positive).gt(tested)).any():
        raise ValueError('Invalid counts of tested cohorts or coefficient signs.')
    return pd.Series(np.select([
        tested.lt(2),
        negative.eq(tested),
        positive.eq(tested),
        negative.gt(0) & positive.gt(0),
    ], ['not_repeated', 'negative', 'positive', 'mixed'], default='zero'), index=pairs.index)


def prepare_direction_results(input_dir, output_dir):
    """Count every repeated exact pair once; retain all its tested cohorts."""
    input_dir, output_dir = Path(input_dir), Path(output_dir)
    tested_distinct = 0
    repeated_chunks = []
    for chunk in pd.read_csv(input_dir / 'all_exact_pair_recurrence.tsv.gz',
                             sep='\t', chunksize=50000):
        tested_distinct += int(chunk.n_tested_cohorts.gt(0).sum())
        repeated_chunks.append(chunk.loc[chunk.n_tested_cohorts.ge(2)].copy())
    repeated = pd.concat(repeated_chunks, ignore_index=True)
    if repeated.empty or repeated.duplicated(['MF_ensembl', 'CpG']).any():
        raise ValueError('Expected unique repeated MF--CpG pairs.')
    signed = repeated.n_negative + repeated.n_positive
    if (signed.gt(repeated.n_tested_cohorts) | repeated.n_negative.lt(0)
            | repeated.n_positive.lt(0)).any():
        raise ValueError('Invalid counts of coefficient signs.')
    # Use the common classification from prepare_pair_evidence.py for both
    # the overview figure and the candidate/reference tables.
    if not repeated.direction_all_tested.isin(DIRECTIONS).all():
        raise ValueError('Invalid direction classification for a repeated pair.')
    if not repeated.same_sign_all_tested.eq(
            repeated.direction_all_tested.isin(DIRECTIONS[:2])).all():
        raise ValueError('Direction labels and consistency flags disagree.')
    counts = (repeated.groupby('direction_all_tested').size()
              .reindex(DIRECTIONS, fill_value=0).rename('n_pairs')
              .rename_axis('direction').reset_index())
    counts['percent'] = 100 * counts.n_pairs / len(repeated)
    cohort_counts = (pd.crosstab(repeated.n_tested_cohorts, repeated.direction_all_tested)
                     .reindex(index=range(2, int(repeated.n_tested_cohorts.max()) + 1),
                              columns=DIRECTIONS, fill_value=0))
    cohort_counts['n_pairs'] = cohort_counts.sum(axis=1)
    cohort_counts = cohort_counts.reset_index()
    count_map = dict(zip(counts.direction, counts.n_pairs))
    same_sign = count_map['negative'] + count_map['positive']
    summary = {
        'scope': 'Selected union MF--CpG pairs with saved tested coefficients a',
        'f4_support_filter': False,
        'comparison_genes': False,
        'tested_distinct_pairs': tested_distinct,
        'repeated_tested_pairs': len(repeated),
        'repeated_percent_of_tested': 100 * len(repeated) / tested_distinct,
        **{'repeated_' + key: int(value) for key, value in count_map.items()},
        'same_sign_pairs': int(same_sign),
        'observed_same_sign_percent': 100 * same_sign / len(repeated),
        'pairs_with_any_zero': int(signed.lt(repeated.n_tested_cohorts).sum()),
        'zero_category': 'Zero coefficients without both positive and negative signs',
    }

    # Fixed manuscript illustrations; selection does not consult support fields.
    # Read only the small set of relevant columns, in chunks, to limit memory.
    key = pd.DataFrame(EXAMPLE_PAIRS, columns=['MF', 'CpG'])
    key['pair_order'] = np.arange(len(key))
    columns = ['MF_ensembl', 'MF', 'CpG', 'cohort', 'tested', 'a', 'n_samples']
    example_chunks = []
    for chunk in pd.read_csv(input_dir / 'all_pair_cohort_evidence.tsv.gz',
                             sep='\t', usecols=columns, chunksize=50000):
        example_chunks.append(chunk.merge(key, on=['MF', 'CpG'], validate='many_to_one'))
    examples = pd.concat(example_chunks, ignore_index=True).sort_values(['pair_order', 'cohort'])
    if examples.duplicated(['MF_ensembl', 'CpG', 'cohort']).any():
        raise ValueError('Duplicate example pair--cohort estimates.')
    tested = examples.loc[examples.tested]
    if set(tested.pair_order) != set(key.pair_order) or not np.isfinite(tested.a).all():
        raise ValueError('A manuscript example has no finite tested coefficient.')

    output_dir.mkdir(parents=True, exist_ok=True)
    for frame, name in [
        (repeated, 'all_repeated_pairs.tsv'),
        (examples, 'plotted_pair_directions.tsv'),
        (counts, 'pair_direction_counts.tsv'),
        (cohort_counts, 'pair_cohort_counts.tsv'),
    ]:
        frame.to_csv(output_dir / name, sep='\t', index=False, na_rep='NA')
    (output_dir / 'pair_direction_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(f'Direction-only results: {len(repeated):,} repeated pairs; '
          f'{same_sign:,} with the same nonzero sign ({summary["observed_same_sign_percent"]:.2f}%).')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / 'analysis_results/pair_directions')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'analysis_results')
    args = parser.parse_args()
    prepare_direction_results(args.input_dir, args.output_dir)


if __name__ == '__main__':
    main()
