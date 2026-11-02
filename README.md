# Motif downstream analysis

Generate the paper's analyses, tables and figures from saved Motif results.
Use the existing `motif` Python environment.

## Folders

| Folder | Purpose |
| --- | --- |
| `pipeline_results/` | Complete server export; retained unchanged |
| `analysis_results/` | Raw analysis tables, manuscript inputs and plotted values |
| `tables/` | Generated LaTeX table bodies and reference table |
| `figures/` | Main-text and supplementary figures |
| `scripts/` | Export, analysis, table and figure code |
| `logs/` | Logs from the latest run |

## Additional group evaluations

Place this project next to the completed `Motif/` run on the server.

After completing the main pipeline and union evaluation, run these commands
on the server from `Motif_downstream_analysis/`. Each group is evaluated and
saved separately. Existing group results are not overwritten.

```bash
(
set -e
conda activate motif
cd ../Motif

group1="../panel_evaluations/section29_13genes_F1-F3_2026-11-02"
group2="../panel_evaluations/section210_17genes_F1-F3_2026-11-02"

test ! -e results/manual_mfs
test ! -e "$group1"
test ! -e "$group2"
mkdir -p "$group1/results" "$group2/results"

motif manually_evaluation --mfs \
    BMP2 ERG FOXA2 GATA3 HNF1A HOXB4 PRDM16 RUNX1 RUNX2 SOX2 TCL1A TP63 WT1
mv results/manual_mfs "$group1/results/"

motif manually_evaluation --mfs \
    BEX2 CELSR3 DLX6 HNF1B HOXA9 HOXB13 HOXB9 IRX2 KNDC1 LHX5 LMX1B \
    MNX1 PITX2 PLEKHG4B PTPRN2 RIMBP2 SIM2
mv results/manual_mfs "$group2/results/"
)
```

After both evaluations finish successfully, continue with the downstream
reproduction commands.

## Reproduce the paper outputs

From this project, in an allocated server session:

```bash
conda activate motif
python scripts/reproduce.py
```

This runs four steps:

1. `prepare_pipeline_results.py` copies the required saved results and reference files into `pipeline_results/`.
2. `run_downstream_analysis.py` prepares the analyses and manuscript table inputs.
3. `build_tables.py` renders the LaTeX tables.
4. `build_figures.py` renders each figure once and saves its plotted data.

The first run requires `pipeline_results/` to be absent and leaves the original
Motif run unchanged. The Motif run must contain the paper's 34 cohorts, three
selected settings and Union. The two completed group evaluations defined in
`scripts/paper_cohorts_and_genes.json` must be available in `../panel_evaluations/`.

If `pipeline_results/` is already prepared, use
`python scripts/reproduce.py --from-export` to repeat the downstream steps.

## Download the results for a local run (optional)

To run the downstream analysis on your local computer, prepare the results
from the server project in an allocated server session. Skip this preparation
if `pipeline_results/` is already available:

```bash
conda activate motif
python scripts/prepare_pipeline_results.py --run-root ../Motif --cpus 4
```

Download the prepared folder from the server. Run this from your local copy
of this project and replace `USER`, `SERVER` and the path with your server details:

```bash
scp -r USER@SERVER:/path/to/Motif_downstream_analysis/pipeline_results ./
```

This copies `pipeline_results/` with all its subfolders and saved evaluations.
From the local project directory, run:

```bash
conda activate motif
python scripts/reproduce.py --from-export
```

The run uses this project and the installed Python environment. It does not
need the original Motif folder, its genome files or a network connection.