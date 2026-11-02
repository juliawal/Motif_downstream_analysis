#!/usr/bin/env python3
"""Render one F4 figure for the additional candidates using the S4 layout.

Read saved cohort results only, including non-testable contexts.
"""
import argparse
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'utils'))
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from cohort_f4_figure import draw_cohort_f4

ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir', type=Path, default=ROOT / 'analysis_results',
                        help='Directory containing the prepared candidate tables')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'figures/supplement')
    args = parser.parse_args()
    data = pd.read_csv(args.input_dir / 'candidate_contexts.tsv', sep='\t')
    shortlist = pd.read_csv(args.input_dir / 'candidate_shortlist.tsv', sep='\t')
    if not shortlist.MF.is_unique or set(data.MF) != set(shortlist.MF):
        raise ValueError('Cohort contexts must cover exactly the additional-candidate shortlist.')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    # The same two-column layout as S4, with a taller canvas for all 66 contexts.
    with plt.rc_context():
        fig = draw_cohort_f4(data, height_mm=196)
        fig.savefig(output / 'candidate_evidence.pdf', facecolor='white')
        plt.close(fig)
    print(f'Wrote one F4 figure for {data.MF.nunique()} candidates in {len(data)} cohort contexts.')


if __name__ == '__main__':
    main()
