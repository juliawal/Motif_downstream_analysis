"""Read saved F1 target-CpG results for the three paper cases."""
from __future__ import annotations

import csv
import gzip
from pathlib import Path

from .inventory import ARMS, branch_path


CONTEXTS = {"PRAD": "ERG", "THCA": "RUNX1", "UCEC": "HNF1B"}


def _rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def export_binding(run_root: Path, out_dir: Path, config: dict, cohorts: list[str]) -> None:

    run_root, out_dir = Path(run_root).resolve(), Path(out_dir)
    destination = out_dir / "cpg_binding_annotations"
    destination.mkdir(parents=True, exist_ok=True)
    target_path = destination / "target_cpg_hits.tsv.gz"
    fields = ["analysis", "cohort", "MF", "MF_ensembl", "CpG", "track", "f1_status", "cpg_passed_filters", "binding_hit"]
    n_targets = 0
    with gzip.open(target_path, "wt", newline="") as target_handle:
        writer = csv.DictWriter(target_handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for cohort in cohorts:
            if cohort not in CONTEXTS:
                continue
            mf = CONTEXTS[cohort]
            for arm in ARMS:
                base = branch_path(run_root, arm, cohort)
                saved_path = base / "evaluation/binding_evidence/target_cpg_hits.tsv"
                for row in _rows(saved_path):
                    if row["MF"] == mf:
                        writer.writerow(dict(analysis=arm, cohort=cohort, **row))
                        n_targets += 1
    print(f"Binding export: {n_targets} target CpG rows", flush=True)
