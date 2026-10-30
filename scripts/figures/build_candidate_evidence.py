#!/usr/bin/env python3
"""Render the additional-candidate figure from prepared evaluation tables.

Internal module for scripts/build_figures.py. Defaults are
relative to this script; --input-dir and --output-dir override the paths.
The full figure caption belongs to the manuscript.
"""
import argparse
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import project_layout as artifacts

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
NAME = "candidate_evidence"
NEGATIVE = "#0072B2"
POSITIVE = "#D55E00"
TEAL = "#007F73"
INK = "#24333E"
MID = "#62727D"
GREY = "#E4E8EC"
GRID = "#EFF2F4"


def draw_figure(data):
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 8.5,
        "text.color": INK, "axes.labelcolor": INK,
        "xtick.color": INK, "ytick.color": INK,
        "pdf.fonttype": 42, "axes.linewidth": .6,
    })
    n = len(data)
    y = np.arange(n)
    fig = plt.figure(figsize=(156 / 25.4, 142 / 25.4))
    fig.text(.025, .98, "Evidence for eight additional candidate MFs",
             fontsize=10.2, weight="bold", va="top", linespacing=1.2)
    fig.text(.025, .892, "A  Expression–methylation associations", fontsize=9.5, weight="bold")
    fig.text(.64, .892, "B  Evaluation results", fontsize=9.5, weight="bold")
    ax = fig.add_axes([.20, .245, .425, .57])
    bx = fig.add_axes([.65, .245, .325, .57])

    tested = data.n_tested_pairs.to_numpy()
    negative = 100 * data.n_a_supported_negative.to_numpy() / tested
    positive = 100 * data.n_a_supported_positive.to_numpy() / tested
    ax.barh(y, negative, color=NEGATIVE, height=.68)
    ax.barh(y, positive, left=negative, color=POSITIVE, height=.68)
    ax.barh(y, 100 - negative - positive, left=negative + positive, color=GREY, height=.68)
    labels = []
    for row in data.itertuples():
        cohort = row.cohort.replace("SKCM_METASTATIC", "SKCM-MET").replace("SKCM_PRIMARY", "SKCM-PRI")
        labels.append(f"{row.MF} / {cohort}")
    ax.set_yticks(y, labels, fontsize=8)
    ax.set(xlim=(0, 159), ylim=(n - .5, -.5), xticks=[0, 50, 100], xlabel="F4-tested pairs (%)")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_bounds(0, 100)
    ax.tick_params(axis="y", length=0, pad=5)
    ax.grid(axis="x", color=GRID, lw=.5)
    ax.set_axisbelow(True)
    ax.text(104, -.93, "F4 pairs\nSupported / tested", fontsize=7.6,
            ha="left", va="bottom", linespacing=1.2)
    for i, row in enumerate(data.itertuples()):
        ax.text(104, i, f"{row.n_supported_pairs:,}/{row.n_tested_pairs:,}", va="center", fontsize=8)

    bx.set(xlim=(-.45, 6.1), ylim=(n - .5, -.5))
    bx.axis("off")
    positions = [0, .8, 1.6, 3.0, 5.15]
    headers = ["F1", "F2", "F3", "Pairs in\nall settings", "Recurrent\npairs"]
    for x, label in zip(positions, headers):
        bx.text(x, -.93, label, ha="center", va="bottom", fontsize=7.3, linespacing=1.2)

    previous_mf = None
    for i, row in enumerate(data.itertuples()):
        available = [row.f1_available, row.f2_available, row.f3_available]
        met = [row.f1_positive, row.f2_positive, row.f3_member_significant_min3_term]
        for x, evaluated, criterion_met in zip(positions[:3], available, met):
            if evaluated:
                bx.scatter(x, i, s=20, facecolor=TEAL if criterion_met else "white",
                           edgecolor=TEAL if criterion_met else MID, lw=.8)
            else:
                bx.text(x, i, "×", color=MID, ha="center", va="center", fontsize=10)
        bx.text(positions[3], i, str(row.n_supported_all_four_pairs),
                ha="center", va="center", fontsize=8.5)
        if row.MF != previous_mf:
            # This count belongs to the MF across cohorts, so centre it on its group.
            group_rows = np.flatnonzero(data.MF.eq(row.MF).to_numpy())
            centre = (group_rows[0] + group_rows[-1]) / 2
            bx.text(positions[4], centre, str(int(row.n_strict_recurrent_pairs)),
                    ha="center", va="center", fontsize=8.5, weight="bold")
            if previous_mf is not None:
                for axis in [ax, bx]:
                    axis.axhline(i - .5, color="#BBC5CC", lw=.55)
        previous_mf = row.MF

    handles = [
        Patch(facecolor=NEGATIVE, label="$a<0$, F4-supported"),
        Patch(facecolor=POSITIVE, label="$a>0$, F4-supported"),
        Patch(facecolor=GREY, label="Tested, no F4 support"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .117),
               ncol=3, frameon=False, fontsize=8, handlelength=1.1, columnspacing=1.0)
    handles = [
        Line2D([], [], marker="o", ls="none", ms=4, mfc=TEAL, mec=TEAL, label="Criterion met"),
        Line2D([], [], marker="o", ls="none", ms=4, mfc="white", mec=MID, label="Criterion not met"),
        Line2D([], [], marker="x", ls="none", ms=4, color=MID, label="Not tested"),
    ]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(.5, .078),
               ncol=3, frameon=False, fontsize=8, handlelength=1.1, columnspacing=1.0)
    fig.text(.025, .055,
             "All settings: identical pairs supported in the union and each of S1–S3 within a cohort.",
             fontsize=7.9)
    fig.text(.025, .031, "Recurrent pairs: support in ≥2 cohorts and the same sign in every tested cohort.", fontsize=7.9)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path,
                        default=ROOT / "analysis_results",
                        help="Directory containing the prepared candidate tables")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/main",
                        help="Figure output directory")
    args = parser.parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.input_dir / "candidate_evidence.tsv", sep="\t")

    fig = draw_figure(data)
    artifacts.save_figure(fig, output_dir / f"{NAME}.pdf", dpi=320, facecolor="white")
    plt.close(fig)
    print(f"Wrote candidate figure for {data.MF.nunique()} MFs in {len(data)} cohort contexts to {output_dir}")


if __name__ == "__main__":
    main()
