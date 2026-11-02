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
        target = ROOT / 'analysis_results'
    target.mkdir(parents=True, exist_ok=True)
    return target
