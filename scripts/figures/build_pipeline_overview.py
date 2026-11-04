#!/usr/bin/env python3
"""Draw the pipeline overview from editable labels, boxes and arrows."""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PROJECT = Path(__file__).resolve().parents[2]
OUTPUT = PROJECT / "figures/main/pipeline_overview.pdf"

# Dimensions in PDF points preserve the original figure's size in the manuscript.
WIDTH, HEIGHT = 430.726, 265.049
INK, ARROW = "#25313A", "#5C6770"
GREEN, PALE_GREEN = "#007E5C", "#EBF7F4"
STAGES = [
    ("Input: Matched TCGA data", "#4A525A", "#F2F3F4"),
    ("A: Prepare data", "#005B8E", "#EBF4F9"),
    ("B: Infer gene–CpG edges", "#667221", "#F5F6EE"),
    ("C1: Filter sex-dependent associations", "#B87F00", "#FDF7EB"),
    ("C2: Filter negative associations near the gene’s TSS", "#B87F00", "#FDF7EB"),
    ("D: Compare selection stability across random seeds", "#A36186", "#FBF4F8"),
    ("E: Select MFs with the three most stable settings", "#AA4B00", "#FCF2EB"),
    ("F0: Prepare inputs for each setting and their union", GREEN, PALE_GREEN),
]
EVALUATIONS = [
    "F1: Binding site",
    "F2: PPI proximity",
    "F3: Gene-set\nenrichment",
    "F4: Supported\nassociations",
]


def render_overview(output: Path = OUTPUT, *, visible_stages=None) -> None:
    """Write a vector PDF without reading an existing figure or requiring LaTeX."""
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "font.weight": "bold",
        "pdf.fonttype": 42,
        "text.color": INK,
        "text.usetex": False,
        "savefig.bbox": None,
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(WIDTH / 72, HEIGHT / 72))
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set(xlim=(0, WIDTH), ylim=(HEIGHT, 0))
        ax.axis("off")

        def box(x, y, width, height, label, border, fill):
            ax.add_patch(FancyBboxPatch(
                (x, y), width, height,
                boxstyle="round,pad=0,rounding_size=1.4",
                linewidth=0.65, edgecolor=border, facecolor=fill,
                zorder=2,
            ))
            if label:
                ax.text(x + width / 2, y + height / 2, label,
                        ha="center", va="center", linespacing=1.25, zorder=4)

        def arrow(x, start, end):
            ax.add_patch(FancyArrowPatch(
                (x, start), (x, end), arrowstyle="-|>",
                mutation_scale=4.6, linewidth=0.6, color=ARROW,
                shrinkA=0, shrinkB=0, zorder=3,
            ))

        left, width = 7.014, 416.698
        top, row_height, row_step = 2.316, 15.591, 20.923
        center = WIDTH / 2
        for index, (label, border, fill) in enumerate(STAGES):
            if visible_stages is not None and index >= visible_stages:
                continue
            y = top + index * row_step
            box(left, y, width, row_height, label, border, fill)
            if index < len(STAGES) - 1 and (visible_stages is None or index+1 < visible_stages):
                arrow(center, y + row_height + 0.5, y + row_step - 0.8)

        if visible_stages is not None and visible_stages <= len(STAGES):
            fig.savefig(output, facecolor='white', metadata={'CreationDate': None, 'ModDate': None})
            plt.close(fig)
            return

        box(left, 170.437, width, 67.805, "", "#00674B", "#FAFDFC")
        arrow(center, top + 7 * row_step + row_height + 0.5, 169.5)
        ax.text(center, 183.5, "Evaluate the union and each of the three settings",
                color="#00573F", ha="center", va="center", zorder=4)

        padding, gap = 7.5, 8.7
        evaluation_left = left + padding
        evaluation_width = (width - 2 * padding - gap * (len(EVALUATIONS) - 1)) / len(EVALUATIONS)
        evaluation_top, evaluation_height, branch_y = 206.6, 25.8, 197.7
        centers = [evaluation_left + i * (evaluation_width + gap) + evaluation_width / 2
                   for i in range(len(EVALUATIONS))]
        ax.plot([center, center], [191, branch_y], color=ARROW, linewidth=0.6, zorder=3)
        ax.plot([centers[0], centers[-1]], [branch_y, branch_y],
                color=ARROW, linewidth=0.6, zorder=3)
        for i, (x, label) in enumerate(zip(centers, EVALUATIONS)):
            if visible_stages is not None and 8+i >= visible_stages:
                continue
            box(x - evaluation_width / 2, evaluation_top, evaluation_width,
                evaluation_height, label, GREEN, PALE_GREEN)
            arrow(x, branch_y, evaluation_top - 0.8)
            if visible_stages is None or visible_stages >= 13:
                arrow(x, evaluation_top + evaluation_height + 0.5, 247.1)

        if visible_stages is None or visible_stages >= 13:
            box(left, 247.966, width, row_height,
                "G: Summarise evaluation results", "#74669A", "#F6F4F9")
        try:
            fig.savefig(output, facecolor='white', metadata={'Title': 'Motif pipeline overview', 'CreationDate': None, 'ModDate': None})
        finally:
            plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--reveal-dir", type=Path, help="Export 13 fixed-position build states for the presentation")
    args = parser.parse_args()
    if args.reveal_dir:
        args.reveal_dir.mkdir(parents=True, exist_ok=True)
        for i in range(1,14):
            render_overview(args.reveal_dir / f'pipeline_reveal_{i:02d}.pdf', visible_stages=i)
        print(f"Wrote 13 pipeline reveal states to {args.reveal_dir}")
    else:
        render_overview(args.output)
        print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
