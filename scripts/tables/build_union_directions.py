#!/usr/bin/env python3
"""Generate the data behind the manuscript's Union F4 direction summary.

Reads the finished F4 summaries. Counts have already been calculated during
analysis; the original Union pair tables are not read again here.
Default tables: tables/supplement/.
No figures are generated.
Paths are resolved relative to this script, not the working directory.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import argparse
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--association-dir', type=Path, default=ROOT / 'analysis_results/associations')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'tables/supplement')
    args = parser.parse_args()
    F4, OUT = args.association_dir.resolve(), args.output_dir.resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    contexts = pd.read_csv(F4 / 'association_gene_summary.tsv', sep='\t').query("analysis == 'union'").copy()
    cohorts = pd.read_csv(F4 / 'association_cohort_summary.tsv', sep='\t').query("analysis == 'union'").copy()
    contexts['negative_a_percent'] = 100 * contexts.n_a_tested_negative / contexts.n_tested_pairs.replace(0, np.nan)
    contexts['positive_a_percent'] = 100 * contexts.n_a_tested_positive / contexts.n_tested_pairs.replace(0, np.nan)
    contexts['direction_summary'] = np.select(
        [contexts.n_tested_pairs.eq(0), contexts.negative_a_percent.gt(50),
         contexts.positive_a_percent.gt(50), contexts.n_a_tested_negative.eq(contexts.n_a_tested_positive)],
        ['not_tested', 'majority_negative', 'majority_positive', 'equal_negative_positive'],
        default='no_direction_majority')
    contexts = contexts.sort_values(['cohort', 'MF', 'MF_ensembl'])
    tested_contexts = contexts[contexts.n_tested_pairs.gt(0)]

    cols = ['analysis', 'cohort', 'MF', 'MF_ensembl', 'n_samples', 'n_input_pairs', 'n_tested_pairs',
            'n_a_tested_negative', 'n_a_tested_positive', 'n_a_tested_zero',
            'negative_a_percent', 'positive_a_percent', 'a_tested_median',
            'n_supported_pairs', 'direction_summary']
    contexts[cols].to_csv(OUT / 'union_directions_by_gene.tsv', sep='\t', index=False, na_rep='NA')

    counts = tested_contexts.groupby('cohort').agg(
        n_tested_MFs=('MF_ensembl', 'size'), n_negative_a_pairs=('n_a_tested_negative', 'sum'),
        n_positive_a_pairs=('n_a_tested_positive', 'sum'),
        n_MFs_majority_negative=('negative_a_percent', lambda x: int(x.gt(50).sum())),
        n_MFs_majority_positive=('positive_a_percent', lambda x: int(x.gt(50).sum())),
        n_MFs_equal_directions=('direction_summary', lambda x: int(x.eq('equal_negative_positive').sum())))
    coverage = cohorts[['analysis', 'cohort', 'n_samples', 'testable', 'n_selected_mfs', 'n_tested_pairs']].merge(
        counts, left_on='cohort', right_index=True, how='left').fillna(0)
    for col in coverage.columns:
        if col.startswith('n_'):
            coverage[col] = coverage[col].astype(int)
    coverage.to_csv(OUT / 'union_directions_by_cohort.tsv', sep='\t', index=False)

    print(f'Union F4 direction tables: {OUT}')

if __name__ == "__main__":
    main()
