#!/usr/bin/env python3
"""Render the overview, F4 coverage and combined-evidence figures.

No models, hypothesis tests, or multiplicity corrections are recomputed here.
The saved figure-data tables specify every display.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from project_layout import figure_data_dir

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextPath
from figure_typography import readable_text
from build_pipeline_overview import render_overview
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
WORK = PROJECT / "analysis_results"
ARMS = ["rank_1", "rank_2", "rank_3", "union"]
BLUE, INK, GRAY, LIGHT = "#0072B2", "#202A33", "#6E7881", "#E7ECEF"
OUT = PROJECT / "figures/main"


def read(path, **kwargs):
    kwargs.setdefault("low_memory", False)
    return pd.read_csv(path, sep="\t", **kwargs)


def data(frame, name):
    frame.to_csv(figure_data_dir(OUT) / f'{name}.tsv', sep='\t', index=False)


def save(fig, name):
    if name in {"association_coverage", "combined_evidence"}:
        readable_text(fig)
    fig.savefig(OUT / f'{name}.pdf', dpi=300, facecolor='white')
    plt.close(fig)


def title(ax, label):
    ax.set_title(label, loc="left", fontsize=8.2, fontweight="bold", pad=10)


def cohort_label(value):
    return value.replace("SKCM_METASTATIC", "SKCM-MET").replace("SKCM_PRIMARY", "SKCM-PRI")


def yes(value):
    if pd.isna(value):
        return None
    if isinstance(value, str):
        return {"true": True, "false": False, "1": True, "0": False}.get(value.lower())
    return bool(value)


def coverage_label(ax, x, y, label, **kwargs):
    # Centre the visible glyphs, excluding unused font ascent/descent space.
    font = FontProperties(family="DejaVu Sans", size=8)
    bounds = TextPath((0, 0), label, prop=font).get_extents()
    offset = (-(bounds.x0 + bounds.x1) / 2, -(bounds.y0 + bounds.y1) / 2)
    return ax.annotate(label, (x, y), xytext=offset, textcoords="offset points",
                       ha="left", va="baseline", fontproperties=font, **kwargs)


def main_coverage(f4dir, integration):
    frame = read(f4dir / "association_coverage.tsv")
    cross = read(integration / "union_evidence_contexts.tsv").query("analysis == 'union'")
    cohorts = sorted(frame.cohort.unique())
    tested = frame.pivot(index="cohort", columns="analysis", values="n_tested_pairs").reindex(index=cohorts, columns=ARMS)
    rates = frame.pivot(index="cohort", columns="analysis", values="pair_support_fraction").reindex(index=cohorts, columns=ARMS) * 100
    rates = rates.where(tested.gt(0))
    # Extra row height leaves padding around the 8 pt cell labels.
    fig = plt.figure(figsize=(156/25.4, 175/25.4), layout="constrained")
    upper = fig.add_gridspec(1, 2, width_ratios=[.85, .15], wspace=.025)
    ax, text = fig.add_subplot(upper[0]), fig.add_subplot(upper[1])
    cmap = LinearSegmentedColormap.from_list("support", ["#F1F6FA", "#B8D9EC", BLUE]); cmap.set_bad("#ECECEC")
    im = ax.imshow(np.ma.masked_invalid(rates.to_numpy()), cmap=cmap, vmin=0, vmax=100, aspect="auto")
    title(ax, "F4 association support across cohorts")
    ax.set_xticks(range(4), ["S1", "S2", "S3", "Union"]); ax.xaxis.tick_top()
    ax.set_yticks(range(len(cohorts)), [cohort_label(x) for x in cohorts]); ax.tick_params(length=0, labelsize=7)
    text.set(xlim=(0, 1), ylim=(len(cohorts)-.5, -.5)); text.axis("off"); text.set_title("Union\nsupport", fontsize=7.5)
    for i, cohort in enumerate(cohorts):
        for j, arm in enumerate(ARMS):
            value = rates.loc[cohort, arm]
            coverage_label(ax, j, i, str(int(tested.loc[cohort, arm])) if pd.notna(value) else "NT", color="white" if pd.notna(value) and value > 70 else INK)
        value = rates.loc[cohort, "union"]
        coverage_label(text, .5, i, f"{value:.2f}%" if pd.notna(value) else "NT")
    ax.set_xticks(np.arange(-.5, 4), minor=True); ax.set_yticks(np.arange(-.5, len(cohorts)), minor=True)
    ax.grid(which="minor", color="white", lw=1); ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    cb = fig.colorbar(im, ax=[ax, text], orientation="horizontal", fraction=.025, pad=.025, aspect=50)
    cb.set_label("Supported / tested pairs (%); NT = no tested pair", fontsize=7); cb.outline.set_visible(False)
    save(fig, "association_coverage")

    fig, bottom = plt.subplots(figsize=(156/25.4, 51/25.4), layout="constrained")
    bottom.axis("off"); bottom.set(xlim=(0, 1), ylim=(-.3, 1.5))
    title(bottom, "Evaluation overlap by MF–cohort combination")
    available = int(cross.all_four_available.map(yes).sum()); positive = int(cross.all_four_positive.map(yes).sum())
    for y, count, total, label in [(1, available, len(cross), "Assessed by all four evaluations"), (.2, positive, available, "Meeting all four criteria")]:
        bottom.text(0, y+.27, label, fontsize=8.5, va="center")
        bottom.add_patch(Rectangle((0, y-.07), .76, .19, facecolor=LIGHT))
        bottom.add_patch(Rectangle((0, y-.07), .76*count/total if total else 0, .19, facecolor=BLUE))
        bottom.text(.79, y+.025, f"{count}/{total}", fontsize=9, va="center", fontweight="bold")
    save(fig, "combined_evidence")


def main():
    global WORK, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=WORK,
                        help="Directory containing the manuscript data tables")
    parser.add_argument("--output-dir", type=Path,
                        help="Override the figure output directory")
    args = parser.parse_args()
    WORK = args.data_dir.resolve()
    OUT = args.output_dir.resolve() if args.output_dir else PROJECT / "figures/main"
    OUT.mkdir(parents=True, exist_ok=True);figure_data_dir(OUT)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":7.2,"axes.labelsize":7,"axes.titlesize":8.2,"xtick.labelsize":6.5,"ytick.labelsize":7,"text.color":INK,"axes.labelcolor":INK,"axes.edgecolor":GRAY,"axes.linewidth":.5,"pdf.fonttype":42,"figure.constrained_layout.w_pad":.04,"figure.constrained_layout.h_pad":.04})
    f4dir = integration = WORK
    render_overview(OUT / "pipeline_overview.pdf")
    main_coverage(f4dir, integration)
    print(f"Wrote the overview, F4 coverage and combined evidence figures to {OUT}")


if __name__ == "__main__":
    main()
