#!/usr/bin/env python3
"""Die Pipeline-Ergebnisse auswerten und die 89 Paper-Datentabellen aufbereiten."""
import argparse
import os
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
    'analyze_candidate_coverage_and_signs.py',
    'prepare_candidate_results.py',
    'compare_literature_genes.py',
]

# Fachmodule teilen dieselben Daten- und Darstellungsregeln. Der Ausgabemodus
# trennt Tabellen und PDFs, ohne wissenschaftliche Berechnungen zu duplizieren.
# Complete retained analysis tables are copied once, before preparing views.
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
        'association_gene_summary', 'association_global_summary', 'binding_case_summary',
        'binding_context_strata', 'tp63_cpg_recurrence', 'celsr3_cpg_recurrence',
    ],
    'combined_evidence': ['combined_evidence_contexts', 'union_evidence_contexts'],
    'manual_candidates': [
        'manual_binding_measurements', 'manual_binding_group_tests', 'manual_network_contexts',
        'manual_network_group_tests', 'manual_enrichment_term_selection', 'manual_enrichment_terms',
        'manual_enrichment_memberships', 'manual_association_contexts',
    ],
    'additional_candidates': ['candidate_evidence'],
    'literature_comparison': [
        'literature_gene_comparison', 'literature_paper_contexts', 'literature_gene_summary',
        'literature_evaluation_summary', 'literature_shortlist_comparison', 'literature_comparison',
    ],
}
DATA_MODULES = [
    'tables/build_union_enrichment_terms.py',
    'tables/build_union_directions.py',
    'tables/build_manuscript_tables.py',
    'tables/build_pair_and_candidate_tables.py',
    'figures/build_candidate_selection.py',
    'figures/build_ucec_motif_binding.py',
    'figures/build_functional_enrichment.py',
    'figures/build_network_proximity.py',
    'figures/build_association_figures.py',
    'figures/build_binding_and_direction.py',
    'figures/build_stability_figures.py',
    'figures/build_pair_directions.py',
]



def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    for number, script in enumerate(ANALYSES, 1):
        print(f"\nAnalyse {number}/{len(ANALYSES)}: {script}", flush=True)
        subprocess.run([sys.executable, "-B", str(ROOT / "scripts/analysis" / script)],
                       cwd=ROOT, check=True)

    output = ROOT / "tables/supplement"
    output.mkdir(parents=True, exist_ok=True)
    for folder, names in ANALYSIS_TABLES.items():
        for name in names:
            shutil.copyfile(ROOT / "analysis_results" / folder / (name + ".tsv"),
                            output / (name + ".tsv"))
    env = dict(os.environ, MOTIF_OUTPUT_KIND="tables", PYTHONDONTWRITEBYTECODE="1")
    for number, module in enumerate(DATA_MODULES, 1):
        print(f"\ntables {number}/{len(DATA_MODULES)}: {module}", flush=True)
        subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / module)],
                       cwd=ROOT, env=env, check=True)
    print("\nFertig: analysis_results/ und tables/supplement/", flush=True)


if __name__ == "__main__":
    main()
