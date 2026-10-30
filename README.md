# Motif_downstream_analysis

Generate the paper's figures and tables from a completed [Motif](https://github.com/juliawal/Motif) run on the same server.

## Data flow

All scripts shown below are in `scripts/`.

```text
Completed Motif run
        │
        │ prepare_pipeline_results.py
        ▼
pipeline_results/
Prepared inputs from the original run
        │
        │ run_downstream_analyses.py
        │   1. Analysis → analysis_results/ (intermediate results)
        │   2. Preparation of manuscript data
        ▼
tables/supplement/
The 89 required data tables
        │
        ├── build_tables.py
        │       → 13 LaTeX table bodies
        │
        └── build_figures.py
                → 17 figures
```

`run_downstream_analyses.py` performs the analyses and prepares the manuscript data. `build_tables.py` and `build_figures.py` are separate subsequent steps. Once the 89 data tables are available, these two rendering steps can run independently of the analysis and pipeline directories.

## Prepare Motif

Place both repositories in the same parent directory, using these exact names:

```text
project/
├── Motif/
└── Motif_downstream_analysis/
```

Follow `Motif/README.md` to install the Motif environment and prepare its data, references and `config.yaml`. This analysis requires the paper's 34 cohorts, three selected settings, their union and the 17 manual genes listed in [paper_cohorts_and_genes.json](scripts/paper_cohorts_and_genes.json).

Run these commands in order, letting each finish successfully. Skip them if all three evaluations are already complete:

```bash
cd /path/to/project/Motif
conda activate motif
motif all
motif union_evaluation
motif manually_evaluation --mfs BMP2 CELSR3 ERG FOXA2 HNF1B HOXA9 HOXB13 HOXB9 IRX2 MNX1 PRDM16 PTPRN2 RUNX1 SOX2 TCL1A TP63 WT1
```

Keep the pipeline code, configuration, results and references in place. See [required files](docs/required_files.md). Reproducing the paper's numbers requires the same input data, references, settings and seeds.

## Run the analysis

Replace `/path/to/project` with your parent directory. Use the existing Motif environment; no additional packages are needed:

```bash
cd /path/to/project/Motif_downstream_analysis
conda activate motif
python scripts/reproduce.py --submit --account YOUR_ACCOUNT --partition YOUR_PARTITION
```

Replace the account and partition with your Slurm settings. The default request is **4 CPUs, 24 GB and 24 hours**; adjust with `--cpus`, `--mem` and `--time`. Wait for the job to finish; its log is `logs/slurm-JOBID.log`.

Alternatively, run `python scripts/reproduce.py` directly in an allocated server session. The script automatically finds the neighbouring `Motif` directory and runs all four steps from the beginning.

## Outputs

All outputs are written inside `Motif_downstream_analysis`:

| Folder | Contents |
| --- | --- |
| `pipeline_results/` | Prepared inputs from the completed pipeline run |
| `analysis_results/` | Calculated results and intermediate tables needed by later steps |
| `tables/supplement/` | Exactly 89 manuscript data tables and the 13 LaTeX table bodies |
| `figures/main/`, `figures/supplement/` | All 17 PDF figures, including the generated pipeline overview |

The [figure and table guide](docs/figures_and_tables.md) lists the scripts and data behind each figure and maps supplementary Tables S1–S13 to their LaTeX labels and data files. Generated files use the same descriptive names as the manuscript, such as `pipeline_overview.pdf`, `functional_enrichment.pdf` and `union_enrichment_terms.tsv`. All exported tables are stored directly in `tables/supplement/`.

The workflow scripts are `scripts/prepare_pipeline_results.py`, `scripts/run_downstream_analyses.py`, `scripts/build_tables.py` and `scripts/build_figures.py`. `scripts/reproduce.py` runs them in that order. Figure and table filenames are independent of their position in the paper.

To run the same four steps individually:

```bash
python scripts/prepare_pipeline_results.py
python scripts/run_downstream_analyses.py
python scripts/build_tables.py
python scripts/build_figures.py
```

`run_downstream_analyses.py` produces the data only. `reproduce.py` explicitly runs table rendering and figure generation as the next two steps.

## Redraw from manuscript data only

Once the 89 data tables are in `tables/supplement/`, regenerate the table bodies and figures with:

```bash
python scripts/build_tables.py
python scripts/build_figures.py
```

These two commands read only the manuscript data tables. They do not require `pipeline_results/`, `analysis_results/`, or the Motif repository. Figure generation leaves the data tables unchanged. The table bodies are written to `tables/supplement/manuscript_rows.tex`; they do not automatically replace the inline tables in `Masterarbeit/main.tex`.

The full workflow prepares the data first. It no longer writes the removed duplicate views or complete test-list exports. Original tests remain in the pipeline inputs; earlier archived exports remain in the separate reproduction archive. Table S6 has nine entries referring to eight files: D and G share `union_gene_recurrence.tsv`. The HNF1B strata remain in `binding_context_strata.tsv`.
