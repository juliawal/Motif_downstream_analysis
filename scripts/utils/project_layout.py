"""Shared output locations; alternative output folders remain self-contained."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

def figure_data_dir(output):
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    try:
        output.relative_to(ROOT / 'figures')
    except ValueError:
        target = output / 'data'
    else:
        target = ROOT / 'tables/supplement'
    target.mkdir(parents=True, exist_ok=True)
    return target

def figure_output_dir(output, name):
    output = Path(output).resolve()
    if output == ROOT / 'figures/main':
        if name == 'manual_candidate_associations':
            output = ROOT / 'figures/supplement'
    output.mkdir(parents=True, exist_ok=True)
    return output


# run_downstream_analyses.py and build_figures.py select their own outputs. Shared calculations are
# identical in both modes. No pandas or matplotlib functions are monkeypatched.
import os
OUTPUT_KIND = os.environ.get("MOTIF_OUTPUT_KIND", "figures")
TABLES = OUTPUT_KIND != "figures"
FIGURES = OUTPUT_KIND != "tables"


def write_dataframe(frame, path, **kwargs):
    if TABLES:
        frame.to_csv(path, **kwargs)


def write_bytes(path, content):
    if TABLES:
        path.write_bytes(content)


def save_figure(figure, path, **kwargs):
    if FIGURES:
        figure.savefig(path, **kwargs)


