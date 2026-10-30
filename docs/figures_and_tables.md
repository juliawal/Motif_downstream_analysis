# Figures and tables

`run_downstream_analysis.py` saves results in `analysis_results/`, including 89 datasets in `analysis_results/paper_data/`. The table and figure scripts read these prepared datasets without requiring the original Motif run, `pipeline_results/` or other analysis folders.

Run from this repository:

```bash
python scripts/build_tables.py
python scripts/build_figures.py
```

These commands leave the datasets and `Masterarbeit/` unchanged. For the full workflow, see the [README](../README.md).

## Figures

`build_figures.py` generates 17 PDFs. Output paths below are relative to this repository; modules are in `scripts/figures/` and data files are in `analysis_results/paper_data/`.

| Output | Module | Data files |
| --- | --- | --- |
| `figures/main/pipeline_overview.pdf` | `build_pipeline_overview.py`, called by `build_association_figures.py` | Diagram defined in the script |
| `figures/supplement/setting_overlap_and_stability.pdf` | `build_stability_figures.py` | `selected_gene_cohort_memberships.tsv`, `selection_overlap.tsv`, `stability_selected_by_cohort.tsv` |
| `figures/main/candidate_selection_and_recurrence.pdf` | `build_candidate_selection.py` | `selected_gene_cohort_memberships.tsv`, `selection_overlap.tsv`, `selection_recurrence.tsv`, `stability_selected_by_cohort.tsv` |
| `figures/main/ucec_motif_binding.pdf` | `build_ucec_motif_binding.py` | `binding_measurements.tsv`, `binding_group_tests.tsv` |
| `figures/main/network_proximity.pdf` | `build_network_proximity.py` | `union_network_effects.tsv`, `network_positive_contributions.tsv`, `network_contribution_concentration.tsv` |
| `figures/main/functional_enrichment.pdf` | `build_functional_enrichment.py` | `enriched_terms.tsv` |
| `figures/main/association_coverage.pdf` | `build_association_figures.py` | `association_coverage.tsv` |
| `figures/main/combined_evidence.pdf` | `build_association_figures.py` | `union_evidence_contexts.tsv` |
| `figures/main/pair_directions.pdf` | `build_pair_directions.py` | `all_repeated_pairs.tsv`, `plotted_pair_directions.tsv` |
| `figures/main/manual_candidate_evidence.pdf` | `build_manual_candidate_figures.py` | `manual_binding_measurements.tsv`, `manual_network_group_tests.tsv`, `manual_enrichment_terms.tsv`, `manual_enrichment_memberships.tsv`, `manual_network_contexts.tsv`, `manual_association_contexts.tsv` |
| `figures/main/erg_and_tcl1a_associations.pdf` | `build_association_figures.py` | `erg_tdrd1_triplets.tsv`, `case_association_table.tsv` |
| `figures/main/binding_and_direction.pdf` | `build_binding_and_direction.py` | `binding_case_summary.tsv` |
| `figures/main/candidate_evidence.pdf` | `build_candidate_evidence.py` | `candidate_evidence.tsv` |
| `figures/main/literature_comparison.pdf` | `build_literature_comparison.py` | `literature_comparison.tsv` |
| `figures/supplement/stability_across_settings.pdf` | `build_stability_figures.py` | `stability_by_cohort.tsv`, `stability_setting_ranking.tsv`, `stability_selected_settings.tsv` |
| `figures/main/cpg_recurrence.pdf` | `build_cpg_recurrence.py` | `tp63_cpg_recurrence.tsv`, `celsr3_cpg_recurrence.tsv` |
| `figures/supplement/manual_candidate_associations.pdf` | `build_manual_candidate_figures.py` | `manual_association_contexts.tsv`, `manual_network_contexts.tsv` |

The literature figure also has PNG and SVG versions. Figure placement and captions are controlled by `Masterarbeit/main.tex`.

## Supplementary tables

`build_tables.py` writes all 13 LaTeX table bodies to `tables/supplements/manuscript_rows.tex`. The current manuscript has no main-text tables; `tables/main/` is reserved for them. Data files below are in `analysis_results/paper_data/`.

| Table | Content | Data file |
| --- | --- | --- |
| S1 | Union candidate selection and evaluation coverage | `union_evaluation_coverage.tsv` |
| S2 | Binding identity exclusions and ChIP-alias sensitivity | `binding_identity_sensitivity.tsv` |
| S3 | Supported union enrichment terms | `union_enrichment_terms.tsv` |
| S4 | Association coverage across cohorts and settings | `association_coverage.tsv` |
| S5 | Union contexts meeting all four evidence criteria | `combined_positive_contexts.tsv` |
| S6 | Index of association overlap datasets | `overlap_index.tsv` |
| S7 | ERG–CpG–TDRD1 triplets in PRAD | `erg_tdrd1_triplets.tsv` |
| S8 | HNF1B binding and expression–methylation directions in UCEC | `hnf1b_binding_directions.tsv` |
| S9 | Genomic-context sensitivity of the HNF1B direction contrast | `hnf1b_genomic_context_sensitivity.tsv` |
| S10 | Association summaries for the four paper cases | `case_association_table.tsv` |
| S11 | Separate evaluation criteria for the four paper cases | `case_criteria_table.tsv` |
| S12 | Additional candidates: evidence and limitations | `candidate_shortlist.tsv` |
| S13 | Aggregate pair evidence for the additional candidates | `candidate_shortlist.tsv` |

S6 lists nine entries referring to eight datasets in `analysis_results/paper_data/`; entries D and G use the same file. Table captions, layouts and numbering remain in `Masterarbeit/main.tex`. The generated bodies do not automatically replace its inline tables.
