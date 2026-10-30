#!/usr/bin/env python3
"""Render the overview, separate F4/combined-evidence figures and ERG/TCL1A.

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
import project_layout as artifacts
from project_layout import figure_data_dir, figure_output_dir
import tempfile

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "motif_v22_figure_mpl"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from figure_typography import readable_text
from build_pipeline_overview import render_overview
import numpy as np
import pandas as pd

PROJECT = Path(__file__).resolve().parents[2]
WORK = PROJECT / "tables/supplement"
ARMS = ["rank_1", "rank_2", "rank_3", "union"]
LABELS = ["S1", "S2", "S3", "Union"]
BLUE, TEAL, ORANGE, INK, GRAY, LIGHT = "#0072B2", "#009E73", "#D55E00", "#202A33", "#6E7881", "#E7ECEF"
OUT = PROJECT / "figures/main"


def read(path, **kwargs):
    kwargs.setdefault("low_memory", False)
    return pd.read_csv(path, sep="\t", **kwargs)


def data(frame, name):
    artifacts.write_dataframe(frame, figure_data_dir(OUT) / f"{name}.tsv", sep="\t", index=False)


def save(fig, name):
    if name in {"association_coverage", "combined_evidence", "erg_and_tcl1a_associations"}:
        readable_text(fig)
    artifacts.save_figure(fig, figure_output_dir(OUT, name) / f"{name}.pdf", dpi=300, facecolor="white")
    plt.close(fig)


def figure(height):
    return plt.figure(figsize=(150 / 25.4, height / 25.4), layout="constrained")


def title(ax, label):
    ax.set_title(label, loc="left", fontsize=8.2, fontweight="bold", pad=10)


def clean(ax):
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color=LIGHT, lw=.5, zorder=0)


def limits(values, *, zero=True, pad=.15):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if zero:
        x = np.append(x, 0.)
    if not len(x):
        return (-1, 1)
    low, high = float(x.min()), float(x.max())
    span = high - low if high > low else max(abs(high), 1.)
    return low - span * pad, high + span * pad


def cohort_label(value):
    return value.replace("SKCM_METASTATIC", "SKCM met.").replace("SKCM_PRIMARY", "SKCM primary")


def yes(value):
    if pd.isna(value):
        return None
    if isinstance(value, str):
        return {"true": True, "false": False, "1": True, "0": False}.get(value.lower())
    return bool(value)


def pct(numerator, denominator):
    return f"{100 * numerator / denominator:.2f}%" if denominator else "NA"


def method_overview():
    render_overview(OUT / "pipeline_overview.pdf")


def main_coverage(f4dir, integration):
    frame = read(f4dir / "association_coverage.tsv")
    cross = read(integration / "union_evidence_contexts.tsv").query("analysis == 'union'")
    cohorts = sorted(frame.cohort.unique())
    tested = frame.pivot(index="cohort", columns="analysis", values="n_tested_pairs").reindex(index=cohorts, columns=ARMS)
    rates = frame.pivot(index="cohort", columns="analysis", values="pair_support_fraction").reindex(index=cohorts, columns=ARMS) * 100
    rates = rates.where(tested.gt(0))
    fig = plt.figure(figsize=(156/25.4, 150/25.4), layout="constrained")
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
            ax.text(j, i, f"{int(tested.loc[cohort, arm]):,}" if pd.notna(value) else "NT", fontsize=6.7, ha="center", va="center", color="white" if pd.notna(value) and value > 70 else INK)
        value = rates.loc[cohort, "union"]
        text.text(.5, i, f"{value:.2f}%" if pd.notna(value) else "NT", ha="center", va="center", fontsize=7)
    ax.set_xticks(np.arange(-.5, 4), minor=True); ax.set_yticks(np.arange(-.5, len(cohorts)), minor=True)
    ax.grid(which="minor", color="white", lw=1); ax.tick_params(which="minor", length=0)
    for spine in ax.spines.values(): spine.set_visible(False)
    cb = fig.colorbar(im, ax=[ax, text], orientation="horizontal", fraction=.025, pad=.025, aspect=50)
    cb.set_label("Supported / tested pairs (%); NT = no tested pair", fontsize=7); cb.outline.set_visible(False)
    save(fig, "association_coverage")

    fig, bottom = plt.subplots(figsize=(156/25.4, 51/25.4), layout="constrained")
    bottom.axis("off"); bottom.set(xlim=(0, 1), ylim=(-.3, 1.5))
    title(bottom, "Evidence overlap by MF–cohort combination")
    available = int(cross.all_four_available.map(yes).sum()); positive = int(cross.all_four_positive.map(yes).sum())
    for y, count, total, label in [(1, available, len(cross), "Assessed by all four evaluations"), (.2, positive, available, "Meeting all four criteria")]:
        bottom.text(0, y+.27, label, fontsize=8.5, va="center")
        bottom.add_patch(Rectangle((0, y-.07), .76, .19, facecolor=LIGHT))
        bottom.add_patch(Rectangle((0, y-.07), .76*count/total if total else 0, .19, facecolor=BLUE))
        bottom.text(.79, y+.025, f"{count:,}/{total:,}", fontsize=9, va="center", fontweight="bold")
    save(fig, "combined_evidence")


def main_cases(f4dir):
    all_erg = read(f4dir / "erg_tdrd1_triplets.tsv")
    erg = all_erg.query("analysis == 'union'").sort_values(["G_prime_distance_bp", "CpG"]).reset_index(drop=True)
    cases = read(f4dir / "case_association_table.tsv")
    tcl = cases.query("MF == 'TCL1A' and cohort == 'TGCT'").set_index("analysis").reindex(ARMS)
    data(tcl.reset_index(), "tcl1a_associations")
    fig = figure(145); outer = fig.add_gridspec(2, 1, height_ratios=[.68, .32], hspace=.2)
    top = outer[0].subgridspec(1, 5, width_ratios=[.23, .27, .27, .11, .12], wspace=.06)
    axes = [fig.add_subplot(top[i]) for i in range(5)]
    count = max(len(erg), 1)
    for ax in axes:
        ax.set_ylim(count-.35, -1.1)
        for i in range(0, count, 2): ax.axhspan(i-.5, i+.5, color="#F5F7F9", lw=0)
    for ax in [axes[0], axes[3], axes[4]]: ax.set_xlim(0,1); ax.axis("off")
    patients = sorted(set(erg.n_samples.dropna().astype(int)))
    n_label = str(patients[0]) if len(patients)==1 else "variable per triplet" if patients else "NA"
    axes[0].set_title(f"A  ERG–CpG–TDRD1 in PRAD (n = {n_label})", loc="left", fontsize=8.1, fontweight="bold", pad=29)
    axes[0].text(0, -.8, "CpG", fontsize=7, color=GRAY)
    axes[3].set_title("F4\nsupport", fontsize=7, pad=10); axes[4].set_title("Binding\noverlap", fontsize=7, pad=10)
    if not len(erg): axes[0].text(0, 0, "No current case triplets", fontsize=7)
    for i, row in erg.iterrows():
        axes[0].text(0, i, row.CpG, fontsize=7, va="center")
        tested = row.triplet_status in {"significant", "not_significant"}
        supported, hit = yes(row.mediation_supported), yes(row.binding_hit)
        axes[3].text(.5, i, "Yes" if supported else "No" if tested else "NT", ha="center", va="center", fontsize=7, color=BLUE if supported else GRAY)
        axes[4].text(.5, i, "NA" if hit is None else "Yes" if hit else "No", ha="center", va="center", fontsize=7, color=BLUE if hit else GRAY)
    coef_limit = limits(erg[["a_estimate", "b_estimate"]].to_numpy().ravel(), pad=.24) if len(erg) else (-1,1)
    for ax, field, label in [(axes[1], "a_estimate", "a: ERG–methylation"), (axes[2], "b_estimate", "b: methylation–TDRD1")]:
        ax.axvline(0, color=GRAY, lw=.65, ls="--")
        for i, row in erg.iterrows():
            value = row[field]
            if pd.isna(value): ax.text(0, i, "NA", fontsize=7); continue
            supported = yes(row.mediation_supported)
            ax.plot([0,value],[i,i], color="#AEBFC9", lw=.8)
            ax.scatter(value, i, s=23, facecolor=BLUE if supported else "white", edgecolor=BLUE if supported else GRAY, zorder=3)
            ax.annotate(f"{value:.3f}", (value,i), xytext=(0,7), textcoords="offset points", ha="center", fontsize=6.1, color=GRAY)
        ax.set(xlim=coef_limit, yticks=[], xlabel="Coefficient"); ax.set_title(label.replace("–", "–\n", 1), fontsize=8, pad=10); clean(ax)
    b = fig.add_subplot(outer[1]); b.axis("off")
    title(b, "B  TCL1A–CpG pairs in TGCT")
    rows = [[label, f"{int(row.n_tested_pairs):,}", f"{int(row.n_supported_pairs):,}",
             pct(row.n_supported_pairs, row.n_tested_pairs)]
            for label, row in zip(LABELS, tcl.itertuples())]
    table = b.table(cellText=rows, colLabels=["Analysis", "Tested", "Supported", "Supported / tested"],
                    cellLoc="center", colWidths=[.17,.23,.23,.37], bbox=[0,.06,1,.86])
    table.auto_set_font_size(False); table.set_fontsize(8)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_facecolor(LIGHT if row == 0 else "#F5F7F9" if row % 2 else "white")
        if row == 0: cell.set_text_props(weight="bold")
    save(fig,"erg_and_tcl1a_associations")


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
    method_overview()
    main_coverage(f4dir, integration)
    main_cases(f4dir)
    print(f"Wrote the overview, F4 coverage, combined evidence and ERG/TCL1A figures to {OUT}")


if __name__ == "__main__":
    main()
