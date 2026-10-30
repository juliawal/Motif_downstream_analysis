#!/usr/bin/env python3
"""Export the saved pair directions and additional candidate tables.

The manuscript table rows are rendered separately by scripts/build_tables.py.
"""
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parents[2]
CANDIDATES = ROOT / 'analysis_results/additional_candidates'
DIRECTIONS = ROOT / 'analysis_results/pair_directions'
OUTPUT = ROOT / 'analysis_results'
OUTPUT.mkdir(parents=True, exist_ok=True)

# Export complete supporting tables and every plotted pair context.
for directory, names in [
    (DIRECTIONS, ['all_repeated_pairs.tsv','plotted_pair_directions.tsv',
                  'strict_recurrent_supported_pairs.tsv','pairs_with_supported_sign_reversals.tsv',
                  'sign_frequency_reference.tsv']),
    (CANDIDATES, ['candidate_shortlist.tsv','candidate_contexts.tsv',
                  'candidate_coverage_flags.tsv','expression_methylation_review_queue.tsv',
                  'screen_threshold_sensitivity.tsv']),
]:
    for name in names:
        shutil.copy2(directory/name,OUTPUT/name)
print(f'Pair and candidate tables: {OUTPUT}')
