#!/usr/bin/env python3
"""Prepare manuscript data and index the accompanying overlap files.

Reads saved analyses, preserving their values and full columns. The manuscript table rows are rendered separately by scripts/build_tables.py.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path
import pandas as pd

ARMS = ["rank_1", "rank_2", "rank_3", "union"]
CASE_ORDER = {"ERG": 0, "TCL1A": 1, "RUNX1": 2, "HNF1B": 3}


def ordered(frame, keys):
    frame = frame.copy()
    frame["_arm"] = frame.analysis.map({a: i for i, a in enumerate(ARMS)}) if "analysis" in frame else 0
    frame["_case"] = frame.MF.map(CASE_ORDER).fillna(99) if "MF" in frame else 0
    return frame.sort_values(keys, kind="stable").drop(columns=["_arm", "_case"])


def build(analysis, output):
    analysis, output = Path(analysis).resolve(), Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    table_data = output
    table_data.mkdir(exist_ok=True)

    def read(name):
        path = analysis / name
        return pd.read_csv(path, sep="\t")

    def save(frame, name):
        path = table_data / name
        frame.to_csv(path, sep="\t", index=False)
        return frame

    coverage = read("combined_evidence/union_coverage.tsv").sort_values("cohort")
    f4 = ordered(read("associations/association_cohort_summary.tsv"), ["cohort", "_arm"])
    cross = read("combined_evidence/union_positive_contexts.tsv").sort_values(["MF", "cohort"])
    erg = ordered(read("associations/erg_tdrd1_binding.tsv"), ["_arm", "G_prime_distance_bp", "CpG"])
    cases = ordered(read("associations/case_association_summary.tsv"), ["_case", "cohort", "_arm"])
    criteria = ordered(read("combined_evidence/case_evaluation_criteria.tsv"), ["_case", "cohort", "_arm"])
    for df, filename in [(coverage, "union_evaluation_coverage.tsv"), (f4, "association_coverage.tsv"),
                         (cross, "combined_positive_contexts.tsv"),
                         (erg, "erg_tdrd1_triplets.tsv"), (cases, "case_association_table.tsv"),
                         (criteria, "case_criteria_table.tsv")]:
        save(df, filename)

    # Include the saved TCL1A values reported in the candidate results.
    tcl = pd.read_csv(table_data / "case_association_table.tsv", sep="\t", low_memory=False)
    tcl = tcl.query("MF == 'TCL1A' and cohort == 'TGCT'").set_index("analysis").reindex(ARMS)
    save(tcl.reset_index(), "tcl1a_associations.tsv")

    # Record filenames relative to this output directory.
    # Copy original bytes so leading-zero state patterns and compressed records
    # survive unchanged; this step performs no overlap analysis or new testing.
    overlap_specs = [
        ("A", "pairwise_setting_overlap.tsv", "Setting comparisons within each cohort"),
        ("B", "all_four_pair_patterns_by_context.tsv", "All four-analysis patterns by MF and cohort"),
        ("C", "all_four_pair_patterns_global.tsv", "All four-analysis patterns pooled across contexts"),
        ("D", "union_gene_recurrence.tsv", "Shared tested and supported pairs by MF"),
        ("E", "union_same_cpg_recurrence.tsv.gz", "Exact pair recurrence across cohorts in the union"),
        ("F", "union_pairwise_cohort_recurrence.tsv", "Pairwise cohort overlaps and coefficient directions"),
        ("G", "union_gene_recurrence.tsv", "Union gene recurrence and exact-pair summaries"),
        ("H", "gene_recurrence.tsv", "Gene recurrence for each setting and the union"),
        ("I", "tcl1a_pair_patterns.tsv", "Every TCL1A input pair and its four-analysis states"),
    ]
    overlap_dir = output
    overlap_dir.mkdir(exist_ok=True)
    overlap_records = []
    copied = {}
    robustness_columns = [
        "MF_ensembl", "MF", "n_common_input_pair_cohort_units",
        "n_common_tested_pair_cohort_units", "n_supported_all_three_settings_common_tested",
        "n_supported_all_four_common_tested", "n_union_supported_without_all_three_common_tested",
        "n_all_three_supported_without_union_common_tested",
        *[f"n_supported_{arm}_common_tested" for arm in ARMS],
    ]
    for key, name, description in overlap_specs:
        relative = "associations/" + name
        source = analysis / relative
        target = overlap_dir / name
        if name not in copied:
            shutil.copyfile(source, target)
            copied[name] = (len(pd.read_csv(target, sep="\t", usecols=[0])),
                            hashlib.sha256(target.read_bytes()).hexdigest())
        rows, checksum = copied[name]
        record = {"key": key, "file": name,
                  "location": "analysis_results",
                  "description": description,
                  "rows": rows, "sha256": checksum, "analysis_source": relative,
                  "columns": ";".join(robustness_columns) if key == "D" else "*"}
        overlap_records.append(record)
    pd.DataFrame(overlap_records).to_csv(output / "overlap_index.tsv", sep="\t", index=False)
    print(f"Wrote six manuscript data tables and an index of {len(copied)} overlap files to {output}")


if __name__ == "__main__":
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, default=root / "analysis_results")
    parser.add_argument("--output", type=Path, default=root / "analysis_results")
    args = parser.parse_args()
    build(args.analysis, args.output)
