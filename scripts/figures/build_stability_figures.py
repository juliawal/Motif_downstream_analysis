#!/usr/bin/env python3
"""Render setting overlap and seed-split stability from saved results.

No models, hypothesis tests or multiple-testing corrections are rerun.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from project_layout import figure_data_dir

import build_association_figures as shared
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
# Match the blue/orange palette in the current F1, F2 and F4 main figures.
BLUE, ORANGE, PALE = "#0072B2", "#D55E00", "#FBFBFA"
STABILITY_CMAP = LinearSegmentedColormap.from_list("paper_stability", [PALE, BLUE])
OVERLAP_CMAP = LinearSegmentedColormap.from_list("paper_overlap", [PALE, ORANGE])


def cell_text_color(rgba):
    """Choose black or white by contrast against the actual cell colour."""
    rgb = np.asarray(rgba[:3])
    linear = np.where(rgb <= .04045, rgb / 12.92, ((rgb + .055) / 1.055) ** 2.4)
    luminance = float(linear @ np.array([.2126, .7152, .0722]))
    return "white" if 1.05 / (luminance + .05) > (luminance + .05) / .05 else "black"




def stability_grid(evidence):
    grid = shared.read(evidence / "stability_by_cohort.tsv")
    ranking = shared.read(evidence / "stability_setting_ranking.tsv")
    selected = shared.read(evidence / "stability_selected_settings.tsv")
    n_cohorts = grid.cohort.nunique()
    summary = grid.groupby(["setting_id", "configuration_name", "k_edges"]).agg(
        n_cohorts=("cohort", "nunique"), n_observed=("mean_kuncheva_stability", "count"),
        median_kuncheva=("mean_kuncheva_stability", "median")).reset_index()
    compare = summary.merge(ranking[["setting_id", "kuncheva_median", "selection_rank"]],
                            on="setting_id")
    score_order = {"max": 0, "mean": 1, "sum": 2}
    names = sorted(grid.configuration_name.unique(),
                   key=lambda name: (score_order.get(name.rsplit("_", 1)[-1], 3), name))
    cuts = sorted(grid.k_edges.unique())
    values = summary.pivot(index="configuration_name", columns="k_edges", values="median_kuncheva").loc[names, cuts]
    fig, ax = plt.subplots(figsize=(150 / 25.4, 183 / 25.4), layout="constrained")
    norm = Normalize(vmin=0, vmax=.6)
    cmap = STABILITY_CMAP
    im = ax.imshow(values.to_numpy(), cmap=cmap, norm=norm, aspect="auto")
    ax.set(xticks=range(len(cuts)), xticklabels=[f"{n:,}" for n in cuts], yticks=range(len(names)),
           yticklabels=names, xlabel="Retained edges per run", ylabel="Aggregation and MF score")
    shared.title(ax, f"Seed stability across all {len(ranking)} settings")
    ax.tick_params(length=0, labelsize=7)
    ax.set_xticks(np.arange(-.5, len(cuts)), minor=True)
    ax.set_yticks(np.arange(-.5, len(names)), minor=True)
    ax.grid(which="minor", color="white", lw=.5)
    ax.tick_params(which="minor", length=0)
    setting_labels = {(r.configuration_name, r.k_edges): r.paper_label for r in selected.itertuples()}
    for i, name in enumerate(names):
        for j, cut in enumerate(cuts):
            value = values.loc[name, cut]
            setting_label = setting_labels.get((name, cut))
            ax.text(j, i, f"{value:.3f}" + (f"  {setting_label}" if setting_label else ""),
                    va="center", ha="center", fontsize=6.8, weight="bold" if setting_label else "normal",
                    color=cell_text_color(cmap(norm(value))))
            if setting_label:
                ax.add_patch(Rectangle((j-.47, i-.43), .94, .86, fill=False, ec="#18252B", lw=1.2))
    fig.colorbar(im, ax=ax, fraction=.04, pad=.025, shrink=.7).set_label(
        "Median of cohort mean Kuncheva indices", fontsize=7)
    shared.data(compare, "stability_medians")
    shared.save(fig, "stability_across_settings")
    return {"settings": len(compare), "cohorts": n_cohorts, "row_order": names,
            "top_three_medians": compare.sort_values("selection_rank").median_kuncheva.head(3).tolist()}


def setting_overlap(evidence):
    members = shared.read(evidence / "selected_gene_cohort_memberships.tsv")
    overlap = shared.read(evidence / "selection_overlap.tsv").set_index("cohort")
    seed = shared.read(evidence / "stability_selected_by_cohort.tsv")
    cohorts = sorted(overlap.index)
    sets = {(c, a): set(d.MF_ensembl) for (c, a), d in members.groupby(["cohort", "analysis"])}
    rows = []
    for cohort in cohorts:
        candidates = [sets[cohort, f"rank_{i}"] for i in (1, 2, 3)]
        for i, j in [(1, 2), (1, 3), (2, 3)]:
            shared_n = len(candidates[i-1] & candidates[j-1])
            union_n = len(candidates[i-1] | candidates[j-1])
            value = shared_n / union_n
            rows.append({"cohort": cohort, "comparison": f"S{i}–S{j}", "shared_mfs": shared_n,
                         "union_mfs": union_n, "jaccard_percent": 100 * value})
    pairs = pd.DataFrame(rows)
    left = pairs.pivot(index="cohort", columns="comparison", values="jaccard_percent").loc[cohorts]
    right = seed.pivot(index="cohort", columns="selection_rank", values="mean_kuncheva_stability").loc[cohorts, [1, 2, 3]]
    fig, axes = plt.subplots(1, 2, figsize=(170 / 25.4, 215 / 25.4), sharey=True,
                             layout="constrained", gridspec_kw={"wspace": .12})
    panels = [(left, list(left.columns), OVERLAP_CMAP, 0, 100, "A  Pairwise setting overlap", "Jaccard index (%)", ".1f"),
              (right, ["S1", "S2", "S3"], STABILITY_CMAP, -.01, 1, "B  Stability across seed splits", "Mean Kuncheva index", ".3f")]
    for ax, (values, labels, cmap, low, high, title, unit, fmt) in zip(axes, panels):
        norm = Normalize(vmin=low, vmax=high)
        ax.imshow(values.to_numpy(), cmap=cmap, norm=norm, aspect="auto")
        ax.set_title(title, loc="left", pad=28, fontsize=8.5, weight="bold")
        ax.text(.5, 1.027, unit, transform=ax.transAxes, ha="center", fontsize=7.4)
        ax.set(xticks=range(len(labels)), xticklabels=labels, yticks=range(len(cohorts)),
               yticklabels=[shared.cohort_label(c) for c in cohorts], ylim=(len(cohorts)-.5, -.5))
        ax.xaxis.tick_top()
        ax.tick_params(length=0, labelsize=7.4)
        ax.set_xticks(np.arange(-.5, len(labels)), minor=True)
        ax.set_yticks(np.arange(-.5, len(cohorts)), minor=True)
        ax.grid(which="minor", color="white", lw=.5)
        ax.tick_params(which="minor", length=0)
        ax.spines[:].set_visible(False)
        for i in range(len(cohorts)):
            for j in range(len(labels)):
                value = values.iloc[i, j]
                ax.text(j, i, format(value, fmt), ha="center", va="center", fontsize=7.5,
                        color=cell_text_color(cmap(norm(value))))
    fig.supxlabel("A: selected lists from all ten runs. B: mean over 20 splits of the same ten runs.", fontsize=7.4)
    shared.data(pairs, "setting_overlap")
    shared.save(fig, "setting_overlap_and_stability")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=ROOT / "analysis_results")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/supplement")
    args = parser.parse_args()
    source = args.data_dir.resolve()
    shared.OUT = args.output_dir.resolve()
    (figure_data_dir(shared.OUT)).mkdir(parents=True, exist_ok=True)
    evidence = source
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.2, "axes.labelsize": 7,
        "axes.titlesize": 8.2, "xtick.labelsize": 6.5, "ytick.labelsize": 7, "text.color": shared.INK,
        "axes.labelcolor": shared.INK, "axes.edgecolor": shared.GRAY, "axes.linewidth": .5,
        "pdf.fonttype": 42, "figure.constrained_layout.w_pad": .04, "figure.constrained_layout.h_pad": .04})
    stability_grid(evidence)
    setting_overlap(evidence)
    print("Rendered two retained stability supplements and accompanying data.")


if __name__ == "__main__":
    main()
