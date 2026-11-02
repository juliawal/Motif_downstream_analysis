"""Render manuscript selection/recurrence from saved memberships; no new inference.

Run: python scripts/figures/build_candidate_selection.py
Outputs: figures/main/candidate_selection_and_recurrence.pdf and plotted TSVs.
"""
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from project_layout import figure_data_dir
import os

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parents[2] / '.cache/matplotlib'))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
from matplotlib.textpath import TextPath
from matplotlib.transforms import ScaledTranslation
import numpy as np
import pandas as pd
from figure_typography import readable_text

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "figures/main"
DATA = figure_data_dir(OUT)
DATA.mkdir(parents=True, exist_ok=True)
BASE = ROOT / "analysis_results"
SOURCES = ["selected_gene_cohort_memberships.tsv",
           "selection_overlap.tsv", "selection_recurrence.tsv",
           "stability_selected_by_cohort.tsv"]
data = {p: pd.read_csv(BASE / p, sep="\t") for p in SOURCES}
m = data[SOURCES[0]]
reference = data[SOURCES[2]].set_index("MF_ensembl")
seed = data[SOURCES[3]].pivot(index="cohort", columns="selection_rank",
                            values="mean_kuncheva_stability")
u = m.loc[m.analysis.eq("union")].copy()
cohorts = sorted(u.cohort.unique())
selection = {(c, a): set(d.MF_ensembl) for (c, a), d in m.groupby(["cohort", "analysis"])}
gene_cohorts = u.groupby("MF_ensembl").cohort.apply(lambda x: set(x))
recur = gene_cohorts.map(len)
summary_rows, detail_rows, rec_rows, pattern_rows = [], [], [], []
max_other = int(recur.max() - 1)

for c in cohorts:
    sets = [selection[c, f"rank_{i}"] for i in [1, 2, 3]]
    union = selection[c, "union"]
    patterns = {p: 0 for p in ["S1", "S2", "S3", "S1;S2", "S1;S3", "S2;S3", "S1;S2;S3"]}
    for g in sorted(union):
        settings = [f"S{i}" for i, s in enumerate(sets, 1) if g in s]
        pattern = ";".join(settings)
        patterns[pattern] += 1
        others = sorted(gene_cohorts[g] - {c})
        eligible = set(reference.loc[g, "eligible_cohorts"].split(";"))
        detail_rows.append({
            "cohort": c, "MF_ensembl": g, "MF": reference.loc[g, "MF"],
            "selected_settings": pattern, "n_selected_settings": len(settings),
            "n_other_selected_cohorts": len(others), "other_selected_cohorts": ";".join(others),
            "n_other_eligible_cohorts": len(eligible - {c}),
            "other_eligible_cohorts": ";".join(sorted(eligible - {c})),
        })
    counts = {n: sum(v for p, v in patterns.items() if len(p.split(";")) == n) for n in [1, 2, 3]}
    row = {"cohort": c, "union_n": len(union), "in_all_three": counts[3],
           "in_exactly_two": counts[2], "in_exactly_one": counts[1],
           "all_three_fraction_union": counts[3] / len(union),
           "all_three_fraction_each_setting": counts[3] / 30}
    for i, j in [(1, 2), (1, 3), (2, 3)]:
        shared = len(sets[i-1] & sets[j-1])
        jac = shared / len(sets[i-1] | sets[j-1])
        row[f"S{i}_S{j}_shared_n"] = shared
        row[f"S{i}_S{j}_jaccard"] = jac
    for i in [1, 2, 3]:
        row[f"S{i}_mean_seed_Kuncheva"] = float(seed.loc[c, i])
    summary_rows.append(row)
    pattern_rows.append({"cohort": c, **patterns})
    hist = {k: sum(recur[g] - 1 == k for g in union) for k in range(max_other + 1)}
    rec_rows.append({"cohort": c, "union_n": len(union),
                     **{f"exactly_{k}_other_cohorts": hist[k] for k in hist},
                     "exactly_9_to_33_other_cohorts": sum(hist[k] for k in hist if k >= 9),
                     "in_at_least_one_other": len(union) - hist[0],
                     "recurrent_fraction_union": (len(union) - hist[0]) / len(union)})

summary = pd.DataFrame(summary_rows).set_index("cohort")
recurrence = pd.DataFrame(rec_rows).set_index("cohort")
details = pd.DataFrame(detail_rows).sort_values(["cohort", "MF", "MF_ensembl"])
patterns = pd.DataFrame(pattern_rows).set_index("cohort")

for name, df in [("cohort_setting_agreement", summary), ("cohort_union_recurrence", recurrence),
                 ("exact_setting_memberships", patterns)]:
    df.to_csv(DATA / f'{name}.tsv', sep='\t')
details.to_csv(DATA / 'union_candidate_details.tsv', sep='\t', index=False)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8,
    "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7.5,
    "ytick.labelsize": 8, "legend.fontsize": 7.2,
    "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": .6, "savefig.dpi": 300,
})
labels = [c.replace("SKCM_METASTATIC", "SKCM-MET").replace("SKCM_PRIMARY", "SKCM-PRI") for c in cohorts]
y = np.arange(len(cohorts))
colors = ["#165A7A", "#75B6D0", "#D8E2E8"]
teal = LinearSegmentedColormap.from_list("count_teal", ["#FFFFFF", "#DAEBE8", "#8CC9BE", "#258979", "#005545"])

fig, axs = plt.subplots(1, 3, figsize=(156/25.4, 165/25.4), sharey=True, layout="constrained",
                        gridspec_kw={"width_ratios": [2.6, .37, 3.45], "wspace": .035})
a, n_ax, b = axs
fig.suptitle("MF selection across settings and cohorts", fontsize=12, fontweight="bold")
a.set_title("A  Agreement across settings", loc="left", pad=31)
left = np.zeros(len(cohorts))
centered_labels = []
for col, color in zip(["in_all_three", "in_exactly_two", "in_exactly_one"], colors):
    values = summary[col].to_numpy()
    a.barh(y, values, left=left, height=.94, color=color, edgecolor="white", linewidth=.35)
    for i, (val, offset) in enumerate(zip(values, left)):
        if val >= 3:
            centered_labels.append(a.text(
                offset + val / 2, i, str(val), ha="center", va="center", fontsize=7,
                color="white" if col == "in_all_three" else "#18313D"))
    left += values
a.set_xlim(0, 56)
a.set_xticks([0, 15, 30, 45])
a.set_xlabel("Number of candidate MFs")
a.set_yticks(y, labels)
a.tick_params(axis="y", length=0)
a.grid(axis="x", color="#D9DFE3", linewidth=.5)
a.set_axisbelow(True)
a.legend([Patch(facecolor=c) for c in colors], ["3", "2", "1"],
         title="Number of settings", title_fontsize=7.5, ncol=3,
         handlelength=1, columnspacing=.8, handletextpad=.35, frameon=False,
         loc="lower center", bbox_to_anchor=(.5, 1.002), borderaxespad=0)
n_ax.set_xlim(0, 1)
n_ax.axis("off")
n_ax.text(.5, 1.014, "Union\nsize", ha="center", va="bottom", fontsize=7.5, transform=n_ax.transAxes)
for i, n in enumerate(summary.union_n):
    n_ax.text(.5, i, str(n), ha="center", va="center", fontweight="bold", fontsize=7.8)

matrix = recurrence[[f"exactly_{k}_other_cohorts" for k in range(max_other + 1)]].to_numpy()
im = b.imshow(matrix, cmap=teal, vmin=0, vmax=int(matrix.max()), aspect="auto", interpolation="nearest")
b.set_title("B  Union-MF recurrence in other cohorts", loc="left", pad=31)
b.text(.5, 1.033, "Number of other cohorts", ha="center", va="bottom", fontsize=7.5, transform=b.transAxes)
b.set_xticks(range(max_other + 1), [str(i) for i in range(max_other + 1)])
b.xaxis.tick_top()
b.tick_params(axis="both", length=0)
b.set_xticks(np.arange(-.5, max_other + 1, 1), minor=True)
b.set_yticks(np.arange(-.5, len(cohorts), 1), minor=True)
b.grid(which="minor", color="#DCE4E5", linewidth=.45)
b.tick_params(which="minor", bottom=False, left=False, top=False)
b.spines[:].set_visible(False)
for i in range(len(cohorts)):
    for k in range(max_other + 1):
        val = int(matrix[i, k])
        centered_labels.append(b.text(
            k, i, str(val) if val else "·", ha="center", va="center", fontsize=7.2,
            color="white" if val >= 29 else ("#203A36" if val else "#A4B4B0")))
b.set_xlabel("Number of union MFs (Cell values are the exact number)", loc="right")
for ax in [a, b]:
    ax.set_ylim(len(cohorts) - .5, -.5)
cb = fig.colorbar(im, ax=b, orientation="horizontal", fraction=.025, pad=.025)
cb.outline.set_visible(False)
fig.align_xlabels([a, b])
readable_text(fig)
# Centre visible glyphs at each segment/cell midpoint after setting final font sizes.
# Font line boxes include unused descender space and otherwise lift the digits.
for label in centered_labels:
    bounds = TextPath((0, 0), label.get_text(), prop=label.get_fontproperties()).get_extents()
    offset = ScaledTranslation(-(bounds.x0 + bounds.x1) / 144,
                               -(bounds.y0 + bounds.y1) / 144, fig.dpi_scale_trans)
    label.set(ha="left", va="baseline")
    label.set_transform(label.get_transform() + offset)
fig.savefig(OUT / 'candidate_selection_and_recurrence.pdf')
plt.close(fig)

