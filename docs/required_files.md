# Required files

`prepare_pipeline_results.py` reads the completed sibling `Motif/` run and writes `pipeline_results/`. Keep the original Motif code, environment, configuration and references available. The export reads saved evaluation tables; it does not import Motif's binding or network helpers.

The run must contain the three selected settings, their union and the manual evaluation for the paper's [34 cohorts and 17 genes](../scripts/paper_cohorts_and_genes.json). See the [README](../README.md) for running the workflow.

## Saved results

Paths below are relative to `Motif/`. `<cohort>` means each configured cohort.

| Path | Contents |
| --- | --- |
| `config.yaml` | Run settings and reference paths |
| `results/stability_comparison/selection.tsv`, `results/stability_comparison/ranking.tsv` | Selected settings and ranking |
| `results/<cohort>/prepared/regulator_names.tsv` | Tested genes |
| `results/<cohort>/stability/grid_stability.tsv` | Stability results |

Each `results/selected_settings/rank_<1–3>/<cohort>/` and `results/union/<cohort>/` directory must contain:

| Path within each directory | Contents |
| --- | --- |
| `selected_mfs/mfs.tsv` | Selected candidates |
| `evaluation_input/mf_cpg_pairs.tsv` | Candidate–CpG pairs |
| `evaluation/binding_sites.tsv`, `evaluation/binding_sites_per_mf.tsv` | Binding results |
| `evaluation/binding_evidence/target_cpg_hits.tsv` | Saved target-CpG eligibility and binding hits; needed for PRAD, THCA and UCEC |
| `evaluation/ppi_proximity.tsv` | Group network results |
| `evaluation/ppi_proximity_per_mf.tsv` | Per-MF scores, background means and contributions for PageRank and harmonic centrality |
| `evaluation/gene_set_enrichment.tsv` | Enrichment results |
| `evaluation/mediation.tsv`, `evaluation/mediation_triplets.tsv` | Association summaries and full triplets |

`ppi_proximity_per_mf.tsv` must contain `MF_ensembl`, `measure`, `contribution`, `score` and `background_mean`. It is copied during preparation. Downstream reads the scores and background means and computes their differences; the saved union contributions are used for the top-one and top-three summaries. MFs without a saved result retain missing scores and background means.

Group p-values are read from `ppi_proximity.tsv`, and observed group means are calculated from the saved scores once per distinct Entrez node. The unsaved random-group background mean and relative group excess are omitted. Entrez IDs come from the gene annotation, and evaluated MFs are identified from the saved per-MF rows. No network is loaded and no degree bins, comparison pools, proximity scores or random sets are recalculated. Graph degree, comparison-pool size and network-background counts are omitted.

F1 writes `target_cpg_hits.tsv` in its existing `binding_evidence/` output directory. Preparation reads the saved rows for ERG/PRAD, RUNX1/THCA and HNF1B/UCEC; sequences, GC content and peak overlaps are not recalculated. Completed runs made before this export was added need their F1 stage rerun to produce it.

F1 availability follows the saved `status == tested`; positive evidence requires a target hit rate above the background hit rate. No additional motif/peak identity filter is applied downstream.

Under `results/manual_mfs/`, keep `selected_mfs/mfs.tsv` and, in `evaluation/`, these four files: `binding_sites.tsv`, `binding_sites_per_mf.tsv`, `ppi_proximity.tsv` and `gene_set_enrichment.tsv`.

## Reference files

Keep the files named by these keys in `config.yaml`. Relative paths are resolved from `Motif/`.

| Configuration keys | Contents |
| --- | --- |
| `gene_annotation_file`, `gene_tss_annotations_file`, `coding_gene_list_file` | Gene annotations |
| `cpg_coordinates_file`, `cpg_gene_annotations_file` | CpG annotations |

The preparation step requires `pipeline_results/` to be absent. Once prepared, `run_downstream_analysis.py` uses this folder. Rebuilding tables and figures needs only `analysis_results/`; see [figures and tables](figures_and_tables.md).
