# Required files

All paths below are relative to the completed sibling `Motif/` run. Configuration resource paths are resolved relative to that run, not the caller's working directory.

For entries such as `genome_fasta`, look up the file path in `Motif/config.yaml`. For example, `genome_fasta: "resources/hg38.fa"` means that the required file is `Motif/resources/hg38.fa`.

| Input | Produced by current Motif | Used by Motif_downstream_analysis |
| --- | --- | --- |
| `config.yaml` | Run configuration | Cohorts, settings, seeds and reference paths |
| `results/stability_comparison/{selection,ranking}.tsv` | Stage D | Selected settings and stability ranking |
| `results/<cohort>/prepared/regulator_names.tsv` | Stage A | Tested gene universe and F2 background |
| `results/<cohort>/prepared/target_names.tsv` | Stage A | Binding eligibility in PRAD, THCA and UCEC |
| `results/<cohort>/stability/grid_stability.tsv` | Stage D | Stability figures and summaries |
| `results/selected_settings/rank_<1–3>/<cohort>/selected_mfs/mfs.tsv` | Stage E | Candidate membership |
| `.../evaluation_input/mf_cpg_pairs.tsv` | Stage F0 | Original pair sets, including untestable pairs |
| `.../evaluation/binding_sites.tsv`, `binding_sites_per_mf.tsv` | Stage F1 | Original group tests and individual binding measurements |
| `.../evaluation/binding_evidence/fimo_motifs.meme` | Stage F1 when motifs are assigned | Motif identity audit |
| `.../evaluation/binding_evidence/fimo_hits_with_cpg_ids.tsv` | Stage F1 | CpG-specific motif matches for the UCEC case |
| `.../evaluation/ppi_proximity.tsv` | Stage F2 | Original PPI group results |
| `.../evaluation/gene_set_enrichment.tsv` | Stage F3 | Functional enrichment and membership |
| `.../evaluation/mediation.tsv`, `mediation_triplets.tsv` | Stage F4 | Original summary, every tested/untestable triplet and support decision |
| `results/union/<cohort>/...` with the same files | `motif union_evaluation` | Union analyses alongside S1–S3 |
| `results/manual_mfs/selected_mfs/mfs.tsv` and `evaluation/{binding_sites,binding_sites_per_mf,ppi_proximity,gene_set_enrichment}.tsv` | `motif manually_evaluation --mfs ...` | Pooled 17-gene F1–F3 analysis; cohort F4 is taken from ordinary results |
| `<chipseq_dir>/peak_files.json`, `<chipseq_dir>/peaks/*.bed` | `prepare_chipseq_resources` | Selected ChIP resource identity and interval overlaps |
| `gene_annotation_file`, `gene_tss_annotations_file`, `coding_gene_list_file` | Gene references configured in `config.yaml` | Gene identities, TSS positions and coding-gene lists |
| `cpg_coordinates_file`, `cpg_gene_annotations_file`, `cpg_annotations_file` | CpG references configured in `config.yaml` | CpG coordinates, binding eligibility and genomic context |
| `genome_fasta`, `genome_fasta_index` | Reference genome and FASTA index configured in `config.yaml` | Sequence eligibility for binding annotations |
| `hippie_network_file` | PPI network configured in `config.yaml` | PPI effects |
