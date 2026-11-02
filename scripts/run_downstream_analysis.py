#!/usr/bin/env python3
"""Analyse the pipeline results and prepare the paper datasets in analysis_results/."""
import argparse
from pathlib import Path
import shutil
import subprocess
import sys
sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parents[1]

ANALYSES = [
    'analyze_selection_and_evidence.py',
    'analyze_associations.py',
    'integrate_evidence.py',
    'analyze_manual_candidates.py',
    'prepare_pair_evidence.py',
    'prepare_direction_results.py',
    'prepare_candidate_results.py',
    'compare_literature_genes.py',
]

# Retain the raw analysis tables and expose the paper inputs in one place.
ANALYSIS_TABLES = {
    'candidate_selection': [
        'selected_gene_cohort_memberships', 'selection_overlap', 'selection_recurrence',
        'stability_by_cohort', 'stability_setting_ranking', 'stability_selected_settings',
        'stability_selected_by_cohort',
    ],
    'binding_network_enrichment': [
        'binding_measurements', 'binding_group_tests', 'network_individual_scores',
        'network_positive_contributions', 'network_contribution_concentration',
        'network_setting_agreement', 'enriched_terms', 'enrichment_test_counts',
        'enrichment_setting_agreement', 'evaluation_coverage_and_tests', 'union_network_effects',
    ],
    'associations': [
        'association_gene_summary', 'association_global_summary',
        'tp63_cpg_recurrence', 'celsr3_cpg_recurrence',
    ],
    'combined_evidence': ['combined_evidence_contexts', 'union_evidence_contexts'],
    'manual_candidates': [
        'manual_binding_measurements', 'manual_binding_group_tests', 'manual_network_contexts',
        'manual_network_group_tests', 'manual_pooled_network_scores', 'manual_enrichment_all_tests',
        'manual_enrichment_terms',
        'manual_enrichment_memberships',
        'manual_panel_membership', 'manual_prioritised_contexts',
    ],
    'additional_candidates': [
        'additional_binding_measurements', 'additional_binding_group_tests',
        'additional_network_group_tests', 'additional_pooled_network_scores',
        'additional_enrichment_all_tests', 'additional_enrichment_terms',
        'additional_enrichment_memberships',
        'candidate_shortlist', 'candidate_contexts',
    ],
    'pair_directions': ['consistent_direction_pairs', 'pairs_with_mixed_directions',
                        'pair_direction_counts_by_gene'],
    'literature_comparison': [
        'literature_gene_comparison', 'literature_paper_contexts', 'literature_gene_summary',
        'literature_evaluation_summary', 'literature_shortlist_comparison', 'literature_comparison',
    ],
}
DATA_MODULES = [
    'tables/build_manual_reference_means.py',
    'tables/build_union_enrichment_terms.py',
    'tables/build_manuscript_tables.py',
]


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    for number, script in enumerate(ANALYSES, 1):
        print(f"\nAnalyse {number}/{len(ANALYSES)}: {script}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts/analysis" / script)],
                       cwd=ROOT, check=True)

    output = ROOT / "analysis_results"
    output.mkdir(parents=True, exist_ok=True)
    for folder, names in ANALYSIS_TABLES.items():
        for name in names:
            shutil.copyfile(ROOT / "analysis_results" / folder / (name + ".tsv"),
                            output / (name + ".tsv"))
    for number, module in enumerate(DATA_MODULES, 1):
        print(f"\ntables {number}/{len(DATA_MODULES)}: {module}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts" / module)],
                       cwd=ROOT, check=True)
    print("\nFinished: analysis_results/", flush=True)


if __name__ == "__main__":
    main()
