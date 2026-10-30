"""Reproduce the adopted Union binding-direction figure from saved F4 counts."""
import argparse
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
import project_layout as artifacts
from project_layout import figure_data_dir
import os

os.environ.setdefault("MPLCONFIGDIR", "/tmp/motif_2_7_figure_mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from figure_typography import readable_text

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source", type=Path, default=ROOT / "analysis_results/paper_data/binding_case_summary.tsv")
parser.add_argument("--output-dir", type=Path, default=ROOT / "figures/main")
args = parser.parse_args()
SOURCE = args.source.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
DATA = figure_data_dir(OUT)
DATA.mkdir(exist_ok=True)
data = pd.read_csv(SOURCE, sep="\t")
data = data[(data.analysis == "union") & (data.population == "tested")]
data = data.set_index("MF").loc[["HNF1B", "ERG", "RUNX1"]].reset_index()

BLUE = "#0072B2"
INK = "#24333E"
GRAY = "#62727D"
GRID = "#E4E9EE"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "text.color": INK, "axes.labelcolor": INK,
                     "xtick.color": INK, "ytick.color": INK,
                     "pdf.fonttype": 42})


def save(fig, name):
    readable_text(fig)
    artifacts.save_figure(fig, OUT / f"{name}.pdf", dpi=240, facecolor="white")
    plt.close(fig)


def setup_axis(ax, ylim, xlabel):
    ax.set(xlim=(0, 100), ylim=ylim, xticks=[0, 25, 50, 75, 100],
           yticks=[], xlabel=xlabel)
    ax.set_axisbelow(True)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.spines[["left", "top", "right"]].set_visible(False)
    ax.spines["bottom"].set_color(GRAY)
    ax.tick_params(axis="x", labelsize=8, length=3, color=GRAY)
    ax.xaxis.labelpad = 10


def bars(ax, row, positions, labels):
    for group, y, label in zip(("hit", "no_hit"), positions, labels):
        negative = int(row[f"n_negative_a_{group}"])
        total = int(row[f"n_a_available_{group}"])
        percentage = 100 * negative / total
        ax.barh(y, percentage, height=0.66, color=BLUE, zorder=3)
        ax.text(percentage - 2, y, f"{percentage:.1f} %", ha="right", va="center",
                fontsize=9, color="white", fontweight="bold")
        ax.text(-0.035, y, label, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8.5)
        ax.text(-0.035, y - 0.40, f"{total:,} CpG sites examined",
                transform=ax.get_yaxis_transform(), ha="right", va="center",
                fontsize=8, color=GRAY)


def all_cases():
    fig = plt.figure(figsize=(156 / 25.4, 112 / 25.4))
    title = "Binding annotations and expression–methylation direction"
    fig.text(0.035, 0.955, title, fontsize=10, fontweight="bold", va="top", linespacing=1.35)
    ax = fig.add_axes([0.34, 0.17, 0.62, 0.63])
    xlabel = ("CpGs with negative expression–methylation associations (%)\n"
              "(F4 coefficient $a < 0$)")
    setup_axis(ax, (-0.9, 8.45), xlabel)
    names = {"HNF1B": "A  HNF1B · UCEC", "ERG": "B  ERG · PRAD", "RUNX1": "C  RUNX1 · THCA"}
    for k, row in data.iterrows():
        top = 7 - 3 * k
        motif = row.MF == "HNF1B"
        group_labels = (["With motif match", "Without motif match"] if motif
                        else ["With peak overlap", "Without peak overlap"])
        bars(ax, row, [top, top - 1], group_labels)
        y_fig = fig.transFigure.inverted().transform(ax.transData.transform((0, top + 0.85)))[1]
        fig.text(0.035, y_fig, names[row.MF], fontsize=10, fontweight="bold", va="center")
    save(fig, "binding_and_direction")


all_cases()
artifacts.write_dataframe(data, DATA / "binding_and_direction.tsv", sep="\t", index=False)
print('Finished binding-direction figure and data export.')
