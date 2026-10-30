# Motif downstream analysis

Generate the paper's results, tables and figures from a completed [Motif](https://github.com/juliawal/Motif) run.

## Requirements

Use the existing `motif` environment. Keep both repositories on the same server, next to each other:

```text
project/
├── Motif/
└── Motif_downstream_analysis/
```

The Motif run must include the three selected settings, union evaluation and manual-gene evaluation for the paper's 34 cohorts and 17 manual genes (scripts/paper_cohorts_and_genes.json).

## Run

Replace the path, account and partition with your server settings:

```bash
cd /path/to/project/Motif_downstream_analysis
conda activate motif
python scripts/reproduce.py --submit --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

This submits all four steps below to Slurm. To run directly in an allocated server session, use `python scripts/reproduce.py`.

## Data flow

All scripts are in `scripts/`; all outputs stay inside this repository.

```text
Completed Motif run
    │ prepare_pipeline_results.py
    ▼
pipeline_results/
    │ run_downstream_analysis.py
    ▼
analysis_results/
    │
    ├── build_tables.py
    │   → tables/supplements/
    │
    └── build_figures.py
        → figures/main/ and figures/supplement/
```