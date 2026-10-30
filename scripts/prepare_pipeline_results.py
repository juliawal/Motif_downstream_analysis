#!/usr/bin/env python3
"""Prepare saved results from the sibling Motif pipeline, on the same server."""
import argparse
import os
from pathlib import Path
import sys
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, default=ROOT.parent / 'Motif')
    parser.add_argument('--cpus', type=int, default=4)
    args = parser.parse_args()
    root = args.run_root.expanduser().resolve(strict=True)
    output = ROOT / 'pipeline_results'
    if ROOT.is_relative_to(root) or root.is_relative_to(ROOT):
        parser.error('Motif and Motif_downstream_analysis must be separate sibling directories')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = str(args.cpus)
    sys.path.insert(0, str(ROOT / 'scripts'))
    if output.exists():
        parser.error(f'{output} already exists; move it aside explicitly before preparing a new run')
    from server.export import run_export
    run_export(root, output)


if __name__ == '__main__':
    main()
