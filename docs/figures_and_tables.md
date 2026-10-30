# Figures and tables

All figure and table filenames use the descriptive names in the manuscript folder `Masterarbeit/`. Names are lowercase English words separated by underscores and do not encode figure or table order. The current manuscript has 17 figures and 13 supplementary tables, numbered S1–S13 by LaTeX. F1–F4 refer to evaluation methods; S1–S3 also identify the three analysis settings.

Paths below are relative to this repository. The workflow is `scripts/prepare_pipeline_results.py` → `scripts/run_downstream_analyses.py` → `scripts/build_tables.py` → `scripts/build_figures.py`. `scripts/reproduce.py` runs these four steps in order.

## Figures

Figure modules below are relative to `scripts/figures/`. All numeric figure inputs are among the 89 retained files in `tables/supplement/`. The figure command works without `analysis_results/` or `pipeline_results/` and does not rewrite the TSV files. The downstream analysis command prepares display summaries through the same modules with figure output disabled.

| Figure file | Figure module | Input under `tables/supplement/` |
| --- | --- | --- |
| `figures/main/pipeline_overview.pdf` | `build_pipeline_overview.py`, called by `build_association_figures.py` | Pipeline methods and layout defined in the script |
| `figures/supplement/setting_overlap_and_stability.pdf` | `build_stability_figures.py` | Selection and stability tables |
| `figures/main/candidate_selection_and_recurrence.pdf` | `build_candidate_selection.py` | Selection and stability tables |
| `figures/main/ucec_motif_binding.pdf` | `build_ucec_motif_binding.py` | `binding_measurements.tsv`, `binding_group_tests.tsv` |
| `figures/main/network_proximity.pdf` | `build_network_proximity.py` | `union_network_effects.tsv`, `network_positive_contributions.tsv`, `network_contribution_concentration.tsv` |
| `figures/main/functional_enrichment.pdf` | `build_functional_enrichment.py` | `enriched_terms.tsv`; shared preparation in `scripts/utils/enrichment_data.py` |
| `figures/main/association_coverage.pdf` | `build_association_figures.py` | `association_coverage.tsv` |
| `figures/main/combined_evidence.pdf` | `build_association_figures.py` | `union_evidence_contexts.tsv` |
| `figures/main/pair_directions.pdf` | `build_pair_directions.py` | `all_repeated_pairs.tsv`, `plotted_pair_directions.tsv` |
| `figures/main/manual_candidate_evidence.pdf` | `build_manual_candidate_figures.py` | Manual binding, network, enrichment and association tables |
| `figures/main/erg_and_tcl1a_associations.pdf` | `build_association_figures.py` | `erg_tdrd1_triplets.tsv`, `case_association_table.tsv` |
| `figures/main/binding_and_direction.pdf` | `build_binding_and_direction.py` | `binding_case_summary.tsv` |
| `figures/main/candidate_evidence.pdf` | `build_candidate_evidence.py` | `candidate_evidence.tsv` |
| `figures/main/literature_comparison.pdf` | `build_literature_comparison.py` | `literature_comparison.tsv` |
| `figures/supplement/stability_across_settings.pdf` | `build_stability_figures.py` | `stability_by_cohort.tsv`, `stability_setting_ranking.tsv`, `stability_selected_settings.tsv` |
| `figures/main/cpg_recurrence.pdf` | `build_cpg_recurrence.py` | `tp63_cpg_recurrence.tsv`, `celsr3_cpg_recurrence.tsv` |
| `figures/supplement/manual_candidate_associations.pdf` | `build_manual_candidate_figures.py` | `manual_association_contexts.tsv` and `manual_network_contexts.tsv` |

The literature figure also has PNG and SVG versions named `literature_comparison.png` and `literature_comparison.svg`.

## Manuscript tables

All paths in the last column are relative to `tables/supplement/`. This is a flat directory, matching the manuscript's table directory. `scripts/build_tables.py` renders the 13 LaTeX bodies into `tables/supplement/manuscript_rows.tex`. Captions, layouts, table numbering and their placement remain controlled by the manuscript. Use the LaTeX labels when matching tables after changes to their order.

| Table | Content | LaTeX label | Data file |
| --- | --- | --- | --- |
| S1 | Union candidate selection and evaluation coverage | `tab:rerun-coverage` | `union_evaluation_coverage.tsv` |
| S2 | Binding identity exclusions and ChIP-alias sensitivity | `tab:rerun-identity` | `binding_identity_sensitivity.tsv` |
| S3 | Supported union enrichment terms | `tab:f3-supported-terms` | `union_enrichment_terms.tsv` |
| S4 | Association coverage across cohorts and settings | `tab:rerun-f4-coverage` | `association_coverage.tsv` |
| S5 | Union contexts meeting all four evidence criteria | `tab:rerun-cross-layer` | `combined_positive_contexts.tsv` |
| S6 | Index of complete association overlap records | `tab:rerun-f4-overlap` | `overlap_index.tsv` |
| S7 | ERG–CpG–TDRD1 triplets in PRAD | `tab:rerun-erg-triplets` | `erg_tdrd1_triplets.tsv` |
| S8 | HNF1B binding and expression–methylation directions in UCEC | `tab:rerun-hnf1b-direction` | `hnf1b_binding_directions.tsv` |
| S9 | Genomic-context sensitivity of the HNF1B direction contrast | `tab:rerun-hnf1b-context` | `hnf1b_genomic_context_sensitivity.tsv` |
| S10 | Association summaries for the four paper cases | `tab:rerun-case-f4` | `case_association_table.tsv` |
| S11 | Separate evaluation criteria for the four paper cases | `tab:rerun-case-criteria` | `case_criteria_table.tsv` |
| S12 | Additional candidates: evidence and qualifications | `tab:candidate-priorities` | `candidate_shortlist.tsv` |
| S13 | Aggregate pair evidence for the additional candidates | `tab:candidate-pair-counts` | `candidate_shortlist.tsv` |

The technical reference inventory is no longer a manuscript table.

| Table module, relative to `scripts/tables/` | Outputs in `tables/supplement/` |
| --- | --- |
| `build_manuscript_tables.py` | Coverage, identity, combined-evidence, case and HNF1B tables; `overlap_index.tsv` and its eight data files (nine entries; D/G share a file) |
| `build_union_enrichment_terms.py` | `union_enrichment_terms.tsv`, `enrichment_term_categories.tsv` |
| `build_union_directions.py` | `union_directions_by_gene.tsv`, `union_directions_by_cohort.tsv` |
| `build_pair_and_candidate_tables.py` | Candidate shortlist, contexts and review tables; repeated pairs and direction data |
| `../build_tables.py` | `manuscript_rows.tex`, using the files above |

The `file` entries in `overlap_index.tsv` are relative to the directory containing that index. Table rendering reads `candidate_shortlist.tsv` from the same manuscript data directory. Run `python scripts/build_tables.py` to regenerate the 13 table bodies without any analysis inputs. S6 includes row counts, file checksums, analysis source paths and the twelve-column selection for entry D. HNF1B strata are retained in `binding_context_strata.tsv` without a separate HNF1B copy.

## Preparation and intermediate files

The complete workflow is `prepare_pipeline_results.py` → `run_downstream_analyses.py` → `build_tables.py` → `build_figures.py`. Only the first two steps need the original run or the prepared pipeline inputs. `run_downstream_analyses.py` runs the eight analyses, copies 41 retained analysis tables, and prepares the remaining manuscript views. It ends with the 89 TSV/TSV.GZ files listed in `Masterarbeit/data_inventory.tsv`. The next two commands independently render the 13 LaTeX table bodies and the 17 PDF figures from those data; the downstream analysis command does not invoke them.

| Folder under `analysis_results/` | Module under `scripts/analysis/` |
| --- | --- |
| `candidate_selection/`, `binding_network_enrichment/` | `analyze_selection_and_evidence.py` |
| `associations/` (including prepared TP63/CELSR3 figure data) | `analyze_associations.py` |
| `combined_evidence/` | `analyze_selection_and_evidence.py`, `integrate_evidence.py` |
| `manual_candidates/` (eight final manuscript tables) | `analyze_manual_candidates.py` |
| `pair_directions/`, `additional_candidates/` | `prepare_pair_evidence.py`, `analyze_candidate_coverage_and_signs.py`, `prepare_candidate_results.py` |
| `literature_comparison/` | `compare_literature_genes.py` |

Some intermediate files are still needed by later analyses or table preparation, for example `network_group_tests.tsv`, `enrichment_gene_contexts.tsv`, `additional_candidate_evidence.tsv` and `binding_context_sensitivity.tsv`. They stay under `analysis_results/` and are not copied to the manuscript. Complete test-list exports, unused diagnostics and duplicate final views are no longer generated. Original complete tests remain in the prepared pipeline inputs. Historical reproduction exports are stored separately in `Reproduktionsarchiv_Masterarbeit/`.

For a fresh reproduction, start with new output directories; scripts do not delete unrelated or historical files from an existing directory. To refresh the manuscript, copy the 89 final data files and the figures to the corresponding manuscript paths. The final S6 index already contains the consolidated D/G mapping and current checksums. Scientific values in the inline manuscript tables are updated separately from generating `manuscript_rows.tex`.
