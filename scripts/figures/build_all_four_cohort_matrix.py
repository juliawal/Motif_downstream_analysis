#!/usr/bin/env python3
"""Plot candidates meeting all four criteria within the same union cohort.

Uses saved union evidence flags; no models or statistical tests are rerun.
Only genes and cohorts with at least one complete four-criterion result are
shown. White cells also include unselected or incompletely assessed contexts.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BLUE, INK, GRID = "#0072B2", "#202A33", "#CBD3D9"
COHORT_LABELS = {"SKCM_METASTATIC": "SKCM-\nMET", "SKCM_PRIMARY": "SKCM-\nPRI"}
CRITERIA = ["f1_positive", "f2_positive", "f3_member_significant_min3_term", "f4_supported"]


def load_matrix(path):
    evidence = pd.read_csv(path, sep="\t", low_memory=False)
    union = evidence.loc[evidence.analysis.eq("union")].copy()
    if union.empty:
        raise ValueError("No union contexts in the input table.")
    if union.duplicated(["MF", "cohort"]).any():
        raise ValueError("Gene symbols do not uniquely identify union contexts.")
    for column in CRITERIA + ["all_four_positive"]:
        values = union[column].astype(str).str.lower()
        if not values.isin(["true", "false"]).all():
            raise ValueError(f"Missing or invalid evidence flag: {column}")
        union[column] = values.eq("true")
    if not union[CRITERIA].all(axis=1).equals(union.all_four_positive):
        raise ValueError("Saved combined evidence disagrees with the four criteria.")

    positive = union.loc[union.all_four_positive]
    if positive.empty:
        raise ValueError("No union context meets all four criteria.")
    genes = sorted(positive.MF.unique())
    cohorts = sorted(positive.cohort.unique())
    matrix = (positive.assign(value=1)
              .pivot(index="cohort", columns="MF", values="value")
              .reindex(index=cohorts, columns=genes).fillna(0).astype(int))
    matrix.columns.name = None
    return matrix


def draw_matrix(matrix, output):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "savefig.facecolor": "white",
    })
    n_cohorts, n_genes = matrix.shape
    fig, ax = plt.subplots(figsize=(156 / 25.4, 146 / 25.4), layout="constrained")
    ax.imshow(matrix.to_numpy().T, cmap=ListedColormap(["white", BLUE]),
              vmin=0, vmax=1, origin="lower", interpolation="none", aspect="auto")
    ax.set_xticks(np.arange(n_cohorts), [COHORT_LABELS.get(c, c) for c in matrix.index],
                  rotation=0, fontstyle="normal")
    ax.set_yticks(np.arange(n_genes), matrix.columns, fontstyle="normal")
    ax.set_xticks(np.arange(-.5, n_cohorts, 1), minor=True)
    ax.set_yticks(np.arange(-.5, n_genes, 1), minor=True)
    ax.grid(which="minor", color=GRID, linewidth=.45)
    ax.tick_params(which="both", length=0, pad=4)
    ax.set_xlabel("TCGA cohort", labelpad=6)
    ax.set_ylabel("Candidate MF", labelpad=7)
    for spine in ax.spines.values():
        spine.set_color(GRID)
        spine.set_linewidth(.45)
    ax.legend(handles=[Patch(facecolor=BLUE, label="All four criteria met (Union)"),
                       Patch(facecolor="white", edgecolor=GRID, label="Not met / unavailable")],
              loc="lower left", bbox_to_anchor=(0, 1.015), borderaxespad=0,
              ncol=2, frameon=False, fontsize=8, handlelength=1.15,
              columnspacing=1.3)
    output.mkdir(parents=True, exist_ok=True)
    stem = output / "all_four_cohort_matrix"
    for extension in ("pdf", "svg", "png"):
        fig.savefig(stem.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "analysis_results")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/main")
    parser.add_argument("--matrix-output", type=Path,
                        help="TSV export path (default: data directory); 1 means all four criteria met in Union.")
    args = parser.parse_args()
    matrix = load_matrix(args.data_dir / "union_evidence_contexts.tsv")
    draw_matrix(matrix, args.output_dir)
    matrix_output = args.matrix_output or args.data_dir / "all_four_cohort_matrix.tsv"
    matrix_output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(matrix_output, sep="\t")
    print(f"All-four matrix: {matrix.shape[1]} genes, {matrix.shape[0]} cohorts, "
          f"{matrix.to_numpy().sum()} filled cells; output: {args.output_dir}")


if __name__ == "__main__":
    main()
