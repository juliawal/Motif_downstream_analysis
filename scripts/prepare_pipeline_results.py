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
    parser.add_argument('--panel-results', type=Path, default=ROOT.parent / 'panel_evaluations',
                        help='Directory containing the completed section 2.9 and 2.10 server runs')
    parser.add_argument('--cpus', type=int, default=4)
    parser.add_argument('--check-only', action='store_true',
                        help='Check current pipeline outputs without writing an export')
    args = parser.parse_args()
    root = args.run_root.expanduser().resolve(strict=True)
    output = ROOT / 'pipeline_results'
    if ROOT.is_relative_to(root) or root.is_relative_to(ROOT):
        parser.error('Motif and Motif_downstream_analysis must be separate sibling directories')
    for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
        os.environ[name] = str(args.cpus)
    sys.path.insert(0, str(ROOT / 'scripts'))
    if args.check_only:
        from server.inventory import load_config
        from server.validation import validate_run, validate_panels
        try:
            paths = validate_run(root, load_config(root)) + validate_panels(args.panel_results)
        except ValueError as error:
            parser.error(str(error))
        print(f'Current Motif output contract satisfied: {len(paths)} files.')
        return
    if output.exists():
        parser.error(f'{output} already exists; move it aside explicitly before preparing a new run')
    from server.export import run_export
    run_export(root, output, args.panel_results.expanduser().resolve())


if __name__ == '__main__':
    main()
