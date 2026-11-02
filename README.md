# Motif downstream analysis

Generate the paper's analyses, tables and figures from saved Motif results.
Use the existing `motif` Python environment.

## Reproduce the paper outputs

Place this project next to the completed `Motif/` run on the server. From this project, in an allocated server session:

```bash
conda activate motif
python scripts/reproduce.py
```

This runs four steps:

1. `prepare_pipeline_results.py` copies the required saved results and reference files into `pipeline_results/`.
2. `run_downstream_analysis.py` prepares the analyses and manuscript table inputs.
3. `build_tables.py` renders the LaTeX tables.
4. `build_figures.py` renders each figure once and saves its plotted data.

The first run requires `pipeline_results/` to be absent and leaves the original Motif run unchanged. The run must contain the paper's 34 cohorts, three selected settings, Union and a manual evaluation matching the additional evaluated gene sets in `scripts/paper_cohorts_and_genes.json`.

If `pipeline_results/` is already prepared, use `python scripts/reproduce.py --from-export` to repeat the downstream steps.

## Download the results for a local run (optional)

To run the downstream analysis on your local computer, prepare the results from the server project in an allocated server session. Skip this preparation if `pipeline_results/` is already available:

```bash
conda activate motif
python scripts/prepare_pipeline_results.py --run-root ../Motif --cpus 4
```

Download the prepared folder from the server. Run this from your local copy of this project and replace `USER`, `SERVER` and the path with your server details:

```bash
scp -r USER@SERVER:/path/to/Motif_downstream_analysis/pipeline_results ./
```

This copies `pipeline_results/` with all its subfolders and saved evaluations. From the local project directory, run:

```bash
conda activate motif
python scripts/reproduce.py --from-export
```

The run uses this project and the installed Python environment. It does not need the original Motif folder, its genome files or a network connection.

## Folders

| Folder | Purpose |
| --- | --- |
| `pipeline_results/` | Complete server export; retained unchanged |
| `analysis_results/` | Raw analysis tables, manuscript inputs and plotted values |
| `tables/` | Generated LaTeX table bodies and reference table |
| `figures/` | Main-text and supplementary figures |
| `scripts/` | Export, analysis, table and figure code |
| `logs/` | Logs from the latest run |