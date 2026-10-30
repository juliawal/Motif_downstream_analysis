#!/usr/bin/env python3
"""Generate analyses, tables and figures from the completed sibling Motif run."""
import argparse
import os
from pathlib import Path
import shlex
import subprocess
import sys

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
STEPS = ['prepare_pipeline_results.py', 'run_downstream_analyses.py',
         'build_tables.py', 'build_figures.py']


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root', type=Path, default=ROOT.parent / 'Motif')
    parser.add_argument('--submit', action='store_true', help='Submit the complete workflow to Slurm')
    parser.add_argument('--account')
    parser.add_argument('--partition')
    parser.add_argument('--cpus', type=int, default=4)
    parser.add_argument('--mem', default='24G')
    parser.add_argument('--time', default='24:00:00')
    args = parser.parse_args()
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1', MPLBACKEND='Agg',
               MPLCONFIGDIR=str(ROOT / '.cache/matplotlib'), PYTHONHASHSEED='0',
               OMP_NUM_THREADS=str(args.cpus), OPENBLAS_NUM_THREADS=str(args.cpus),
               MKL_NUM_THREADS=str(args.cpus))
    command = [sys.executable, '-B', '-u']
    if args.submit:
        (ROOT / 'logs').mkdir(exist_ok=True)
        invocation = command + [str(Path(__file__).resolve()), '--run-root',
                               str(args.run_root.expanduser().resolve()),
                               '--cpus', str(args.cpus)]
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
    for number, script in enumerate(STEPS, 1):
        step_command = command + [str(ROOT / 'scripts' / script)]
        if number == 1:
            step_command += ['--run-root', str(args.run_root.expanduser().resolve()), '--cpus', str(args.cpus)]
        print(f'\nStep {number}/{len(STEPS)}: {script}', flush=True)
        with (ROOT / 'logs' / (script.removesuffix('.py') + '.log')).open('w') as log:
            result = subprocess.run(step_command, cwd=ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            print(log.name + ':\n' + Path(log.name).read_text()[-10000:], file=sys.stderr)
            raise SystemExit(result.returncode)
    print('Finished. Figures: figures/; tables: tables/.', flush=True)


if __name__ == '__main__':
    main()
