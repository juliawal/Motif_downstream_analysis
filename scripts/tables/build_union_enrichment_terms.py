#!/usr/bin/env python3
"""Prepare the shared data for the functional-enrichment table and figure.

From the project root:
    python scripts/tables/build_union_enrichment_terms.py

Reads analysis_results/binding_network_enrichment/enriched_terms.tsv.
Writes the full term records and category mapping in analysis_results/.
The files supply tab:f3-enriched-terms and functional_enrichment.pdf. The manuscript
rows are rendered by scripts/build_tables.py.

Use --evidence-dir and --output-dir for alternative folders. Data preparation
uses enrichment_data.py once; the figure reads these prepared tables. No working drafts,
existing output tables or graphics are needed. No tests are recomputed.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import argparse

from enrichment_data import ROOT, prepare


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--evidence-dir', type=Path, default=ROOT / 'analysis_results/binding_network_enrichment',
                        help='Directory containing enriched_terms.tsv')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'analysis_results',
                        help='Directory for the prepared enrichment term, category and cohort tables')
    args = parser.parse_args()
    source = args.evidence_dir.resolve() / 'enriched_terms.tsv'
    terms, rank, summary, mapping = prepare(source)
    output_dir = args.output_dir.resolve()
    data_dir = output_dir
    data_dir.mkdir(parents=True, exist_ok=True)
    terms.to_csv(data_dir / 'union_enrichment_terms.tsv', sep='\t', index=False)
    mapping.to_csv(data_dir / 'enrichment_term_categories.tsv', sep='\t', index=False)
    rank.to_csv(data_dir / 'enrichment_cohort_ranking.tsv', sep='\t', index=False)
    summary.to_csv(data_dir / 'enrichment_category_summary.tsv', sep='\t', index=False)
    print(f'F3 term tables: {data_dir}')


if __name__ == '__main__':
    main()
