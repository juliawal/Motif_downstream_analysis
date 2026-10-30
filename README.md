# Motif_downstream_analysis

Generate the paper's results, tables and figures from a completed [Motif](https://github.com/juliawal/Motif) run.

## Requirements

Use the existing `motif` environment. Keep both repositories on the same server, next to each other:

```text
project/
├── Motif/
└── Motif_downstream_analysis/
```

The Motif run must include the three selected settings, union evaluation and manual-gene evaluation for the paper's [34 cohorts and 17 manual genes](scripts/paper_cohorts_and_genes.json). Keep its configuration, results and reference files in place; see [required files](docs/required_files.md).

## Run

Replace the path, account and partition with your server settings:

```bash
cd /path/to/project/Motif_downstream_analysis
conda activate motif
python scripts/reproduce.py --submit --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

This submits all four steps below to Slurm using 4 CPUs, 24 GB and 24 hours. Logs are saved in `logs/`. To run directly in an allocated server session, use `python scripts/reproduce.py`.

## Data flow

All scripts are in `scripts/`; all outputs stay inside this repository.

```text
Completed Motif run
    │ prepare_pipeline_results.py
    ▼
pipeline_results/ — prepared inputs
    │ run_downstream_analysis.py
    ▼
analysis_results/ — 89 datasets, with intermediate data in subfolders
    │
    ├── build_tables.py
    │   → tables/supplements/ — 13 supplementary tables
    │
    └── build_figures.py
        → figures/main/ and figures/supplement/ — 17 figures
```

The analysis produces data; table and figure generation are separate steps. All current tables (S1–S13) are supplementary. Their LaTeX bodies are saved together in `tables/supplements/manuscript_rows.tex`. `tables/main/` is reserved for main-text tables and is currently empty.

## Rebuild tables and figures

Once the 89 datasets are available in `analysis_results/`:

```bash
python scripts/build_tables.py
python scripts/build_figures.py
```

Both commands read only these prepared datasets. They leave `Masterarbeit` unchanged. See the [figure and table guide](docs/figures_and_tables.md) for individual outputs and data sources.
