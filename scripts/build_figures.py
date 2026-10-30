#!/usr/bin/env python3
"""Alle Abbildungen des Papers und Supplements erzeugen."""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
# Fachmodule teilen dieselben Daten- und Darstellungsregeln. Der Ausgabemodus
# trennt Tabellen und PDFs, ohne wissenschaftliche Berechnungen zu duplizieren.
MODULES = [
    'figures/build_candidate_selection.py',
    'figures/build_ucec_motif_binding.py',
    'figures/build_functional_enrichment.py',
    'figures/build_network_proximity.py',
    'figures/build_association_figures.py',
    'figures/build_cpg_recurrence.py',
    'figures/build_manual_candidate_figures.py',
    'figures/build_binding_and_direction.py',
    'figures/build_stability_figures.py',
    'figures/build_pair_directions.py',
    'figures/build_candidate_evidence.py',
    'figures/build_literature_comparison.py',
]


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    env = dict(os.environ, MOTIF_OUTPUT_KIND="figures", PYTHONDONTWRITEBYTECODE="1")
    for number, module in enumerate(MODULES, 1):
        print(f"\nfigures {number}/{len(MODULES)}: {module}", flush=True)
        subprocess.run([sys.executable, "-B", str(ROOT / "scripts" / module)],
                       cwd=ROOT, env=env, check=True)
    print("\nFertig: figures/", flush=True)


if __name__ == "__main__":
    main()
