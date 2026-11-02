#!/usr/bin/env python3
"""Generate the paper outputs from a Motif run or an existing local export."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
STEPS = ['prepare_pipeline_results.py', 'run_downstream_analysis.py',
         'build_tables.py', 'build_figures.py']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group()
    source.add_argument('--run-root', type=Path, default=ROOT.parent / 'Motif',
                        help='Completed Motif run to export (default: sibling Motif directory)')
    source.add_argument('--from-export', action='store_true',
                        help='Use local pipeline_results/ without accessing or exporting Motif')
    parser.add_argument('--panel-results', type=Path, default=ROOT.parent / 'panel_evaluations',
                        help='Directory containing the completed section 2.9 and 2.10 server runs')
    parser.add_argument('--submit', action='store_true', help='Submit the selected workflow to Slurm')
    parser.add_argument('--account')
    parser.add_argument('--partition')
    parser.add_argument('--cpus', type=int, default=4)
    parser.add_argument('--mem', default='24G')
    parser.add_argument('--time', default='24:00:00')
    args = parser.parse_args()
    steps = STEPS[1:] if args.from_export else STEPS
    if args.from_export:
        required = [ROOT / 'pipeline_results/run_configuration/config.yaml',
                    ROOT / 'pipeline_results/source_file_index.tsv']
        missing = [str(path.relative_to(ROOT)) for path in required if not path.is_file()]
        if missing:
            parser.error('--from-export requires a completed current export in pipeline_results/. '
                         'Missing: ' + ', '.join(missing) + '. Copy the complete export first; '
                         'this mode does not access Motif or create missing exports.')
        from server.validation import validate_panels
        try:
            validate_panels(ROOT / 'pipeline_results/panel_evaluations')
        except ValueError as error:
            parser.error(str(error))
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg',
               MPLCONFIGDIR=str(ROOT / '.cache/matplotlib'), PYTHONHASHSEED='0',
               OMP_NUM_THREADS=str(args.cpus), OPENBLAS_NUM_THREADS=str(args.cpus),
               MKL_NUM_THREADS=str(args.cpus))
    command = [sys.executable, '-u']
    if args.submit:
        (ROOT / 'logs').mkdir(exist_ok=True)
        invocation = command + [str(Path(__file__).resolve()), '--cpus', str(args.cpus)]
        if args.from_export:
            invocation.append('--from-export')
        else:
            invocation += ['--run-root', str(args.run_root.expanduser().resolve()),
                           '--panel-results', str(args.panel_results.expanduser().resolve())]
        sbatch = ['sbatch', '--job-name=motif_downstream_analysis', '--nodes=1', '--ntasks=1',
                  f'--cpus-per-task={args.cpus}', f'--mem={args.mem}', f'--time={args.time}',
                  '--no-requeue', '--export=ALL', '--chdir=' + str(ROOT),
                  '--output=' + str(ROOT / 'logs/slurm-%j.log')]
        for flag, value in [('--account', args.account), ('--partition', args.partition)]:
            if value:
                sbatch.append(flag + '=' + value)
        subprocess.run(sbatch, input='#!/bin/bash\nset -euo pipefail\nexec ' + shlex.join(invocation) + '\n',
                       text=True, env=env, check=True)
        return
    (ROOT / 'logs').mkdir(exist_ok=True)
    for number, script in enumerate(steps, 1):
        step_command = command + [str(ROOT / 'scripts' / script)]
        if script == 'prepare_pipeline_results.py':
            step_command += ['--run-root', str(args.run_root.expanduser().resolve()), '--cpus', str(args.cpus),
                             '--panel-results', str(args.panel_results.expanduser().resolve())]
        print(f'\nStep {number}/{len(steps)}: {script}', flush=True)
        with (ROOT / 'logs' / (script.removesuffix('.py') + '.log')).open('w') as log:
            result = subprocess.run(step_command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.name + ':\n' + Path(log.name).read_text()[-10000:], file=sys.stderr)
            raise SystemExit(result.returncode)
    print('Finished. Figures: figures/; tables: tables/.', flush=True)


if __name__ == '__main__':
    main()
