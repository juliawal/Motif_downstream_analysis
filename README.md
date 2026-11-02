# Motif downstream analysis

Generate the paper's analyses, tables and figures from saved Motif results.
Use the existing `motif` Python environment. These scripts summarise completed
server evaluations; they do not rerun F1–F4 or calculate new significance tests.

## Reproduce the paper outputs

Place this project next to the completed `Motif/` run. The run must contain the
paper's 34 cohorts, three selected settings and Union. Also provide the completed
F1–F3 evaluations for the 13 genes in section 2.9 and the 17 genes in section 2.10:

```text
Motif/
Motif_downstream_analysis/
panel_evaluations/
  section29_13genes_F1-F3_2026-11-02/results/manual_mfs/...
  section210_17genes_F1-F3_2026-11-02/results/manual_mfs/...
```

The two group-result archives can be extracted into `panel_evaluations/` with
`tar -xzf /path/to/ARCHIVE_results.tar.gz -C ../panel_evaluations` after creating
that directory. The gene lists and directory names are defined in
`scripts/paper_cohorts_and_genes.json`.

From this project, in an allocated server session:

```bash
conda activate motif
python scripts/reproduce.py
```

For results stored elsewhere, pass `--run-root /path/to/Motif` and
`--panel-results /path/to/panel_evaluations`.

The command runs four steps:

1. `prepare_pipeline_results.py` copies saved cohort results, both group evaluations and required reference files into `pipeline_results/`.
2. `run_downstream_analysis.py` prepares descriptive analyses and manuscript table inputs.
3. `build_tables.py` renders the four supplementary table bodies and the main-text reference table.
4. `build_figures.py` renders the 16 figures used by the manuscript and saves their data.

The first run requires `pipeline_results/` to be absent and leaves the source
runs unchanged. Each group is checked against its configured gene list.
The section summaries and original group p-values are saved in
`analysis_results/section_evaluation_summary.tsv`. F3 retains the original BH
correction and reports terms with adjusted p-value ≤ 0.05 and at least three
candidate genes. F4 and the cohort-specific F2 points use the original Union
results. `manual17` records the historical gene set used for the original
974-gene screen; it is not the current section 2.10 group.

## Run from a downloaded export

If needed, prepare the export in an allocated server session:

```bash
conda activate motif
python scripts/prepare_pipeline_results.py --run-root ../Motif --panel-results ../panel_evaluations --cpus 4
```

Download the entire prepared directory, including `panel_evaluations/`. Run the
following from your local copy of this project, replacing the server details:

```bash
scp -r USER@SERVER:/path/to/Motif_downstream_analysis/pipeline_results ./
conda activate motif
python scripts/reproduce.py --from-export
```

This mode uses only the export and the installed Python environment. It does not
need the original Motif folder, genome files or a network connection. The same
command regenerates the outputs when `pipeline_results/` is already available.

## Output folders

| Folder | Contents |
| --- | --- |
| `pipeline_results/` | Saved server results and reference files; retained unchanged |
| `analysis_results/` | Descriptive analyses, manuscript inputs and plotted values |
| `tables/main/` | Main-text reference table |
| `tables/supplement/` | Four supplementary table bodies in `manuscript_rows.tex` |
| `figures/main/`, `figures/supplement/` | Figures referenced by the manuscript |
| `logs/` | Logs from the latest run |

The supplementary captions and column layouts, and the manually written criteria
table, remain in the manuscript. Generated data, figures, logs and caches are
excluded from Git. The repository contains only the scripts, their configuration,
this README and `.gitignore`.
