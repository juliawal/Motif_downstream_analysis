#!/usr/bin/env python3
"""Render all main-text and supplementary figures from saved analysis data."""
import argparse
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
# Each figure module runs once and saves its plotted data alongside the analyses.
MODULES = [
    'figures/build_candidate_selection.py',
    'figures/build_ucec_motif_binding.py',
    'figures/build_functional_enrichment.py',
    'figures/build_network_proximity.py',
    'figures/build_association_figures.py',
    'figures/build_all_four_cohort_matrix.py',
    'figures/build_cpg_recurrence.py',
    'figures/build_manual_candidate_figures.py',
    'figures/build_stability_figures.py',
    'figures/build_pair_directions.py',
    'figures/build_candidate_evidence.py',
    'figures/build_literature_comparison.py',
]


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    for number, module in enumerate(MODULES, 1):
        print(f"\nfigures {number}/{len(MODULES)}: {module}", flush=True)
        subprocess.run([sys.executable, str(ROOT / "scripts" / module)],
                       cwd=ROOT, check=True)
    print("\nFinished: figures/", flush=True)


if __name__ == "__main__":
    main()
