# Required files

`prepare_pipeline_results.py` reads the completed sibling `Motif/` run and writes `pipeline_results/`. Keep the original Motif code, environment, configuration and references available. The export uses Motif's binding and network helpers.

The run must contain the three selected settings, their union and the manual evaluation for the paper's [34 cohorts and 17 genes](../scripts/paper_cohorts_and_genes.json). See the [README](../README.md) for running the workflow.

## Saved results

Paths below are relative to `Motif/`. `<cohort>` means each configured cohort.

| Path | Contents |
| --- | --- |
| `config.yaml` | Run settings and reference paths |
| `results/stability_comparison/selection.tsv`, `results/stability_comparison/ranking.tsv` | Selected settings and ranking |
| `results/<cohort>/prepared/regulator_names.tsv` | Tested genes |
| `results/<cohort>/stability/grid_stability.tsv` | Stability results |
| `results/<cohort>/prepared/target_names.tsv` | Prepared CpGs; needed for PRAD, THCA and UCEC |

Each `results/selected_settings/rank_<1–3>/<cohort>/` and `results/union/<cohort>/` directory must contain:

| Path within each directory | Contents |
| --- | --- |
| `selected_mfs/mfs.tsv` | Selected candidates |
| `evaluation_input/mf_cpg_pairs.tsv` | Candidate–CpG pairs |
| `evaluation/binding_sites.tsv`, `evaluation/binding_sites_per_mf.tsv` | Binding results |
| `evaluation/ppi_proximity.tsv` | Network results |
| `evaluation/gene_set_enrichment.tsv` | Enrichment results |
| `evaluation/mediation.tsv`, `evaluation/mediation_triplets.tsv` | Association summaries and full triplets |

Keep `evaluation/binding_evidence/fimo_motifs.meme` where produced; it is copied when present. `fimo_hits_with_cpg_ids.tsv` in that directory is needed when the PRAD/ERG, THCA/RUNX1 or UCEC/HNF1B case uses FIMO.

Under `results/manual_mfs/`, keep `selected_mfs/mfs.tsv` and, in `evaluation/`, these four files: `binding_sites.tsv`, `binding_sites_per_mf.tsv`, `ppi_proximity.tsv` and `gene_set_enrichment.tsv`.

## Reference files

Keep the files named by these keys in `config.yaml`. Relative paths are resolved from `Motif/`.

| Configuration keys | Contents |
| --- | --- |
| `gene_annotation_file`, `gene_tss_annotations_file`, `coding_gene_list_file` | Gene annotations |
| `cpg_coordinates_file`, `cpg_gene_annotations_file`, `cpg_annotations_file` | CpG annotations |
| `genome_fasta`, `genome_fasta_index` | Reference genome and index |
| `hippie_network_file` | Interaction network |
| `chipseq_dir` | `peak_files.json` and its referenced peak files; older runs may use `manifest.json` |

The preparation step requires `pipeline_results/` to be absent. Once prepared, `run_downstream_analysis.py` uses this folder. Rebuilding tables and figures needs only `analysis_results/`; see [figures and tables](figures_and_tables.md).
