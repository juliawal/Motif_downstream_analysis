"""Check the current Motif output contract before starting the paper export."""
import csv
import hashlib
import json
from pathlib import Path

from .inventory import ARMS, EVALUATIONS, RESOURCE_COPY_KEYS, STABILITY, PANEL_FILES, branch_path, configured_path

PAPER = Path(__file__).resolve().parents[1] / 'paper_cohorts_and_genes.json'
CASE_COHORTS = {'PRAD', 'THCA', 'UCEC'}
DETAIL_COLUMNS = {
    'ppi_proximity_per_mf.tsv': {'MF_ensembl', 'measure', 'contribution', 'score', 'background_mean'},
    'target_cpg_hits.tsv': {'MF', 'MF_ensembl', 'CpG', 'track', 'f1_status', 'cpg_passed_filters', 'binding_hit'},
    'mediation_triplets.tsv': {'MF', 'MF_ensembl', 'CpG', 'G_prime_ensembl', 'G_prime_distance_bp',
                              'a_estimate', 'a_pvalue', 'b_estimate', 'b_pvalue',
                              'joint_adj_p_value', 'triplet_status', 'mediation_supported'},
}


def required_paths(root, config):
    paths = {root / 'config.yaml'}
    paths.update(configured_path(root, config, key) for key in RESOURCE_COPY_KEYS)
    paths.update(root / 'results/stability_comparison' / name for name in STABILITY)
    for cohort in config['cohorts']:
        paths.update(root / 'results' / cohort / name for name in
                     ['prepared/regulator_names.tsv', 'stability/grid_stability.tsv'])
        for arm in ARMS:
            branch = branch_path(root, arm, cohort)
            names = ['selected_mfs/mfs.tsv', 'evaluation_input/mf_cpg_pairs.tsv',
                     *['evaluation/' + name for name in EVALUATIONS],
                     'evaluation/ppi_proximity.tsv', 'evaluation/mediation.tsv',
                     'evaluation/mediation_triplets.tsv']
            if cohort in CASE_COHORTS:
                names.append('evaluation/binding_evidence/target_cpg_hits.tsv')
            paths.update(branch / name for name in names)
    return sorted({path.resolve() for path in paths})


def validate_run(root, config, *, paper=None):
    root = Path(root).resolve()
    paper = json.loads(PAPER.read_text()) if paper is None else paper
    errors = []
    cohorts = config.get('cohorts', [])
    if len(cohorts) != len(set(cohorts)) or set(cohorts) != set(paper['cohorts']):
        errors.append('config.yaml must specify the paper cohorts in scripts/paper_cohorts_and_genes.json.')
    if config.get('n_selected_stable_settings') != 3:
        errors.append('The paper requires n_selected_stable_settings: 3.')
    if config.get('should_aggregate_cpgs', False):
        errors.append('The paper requires individual CpGs: should_aggregate_cpgs must be false.')
    missing_keys = set(RESOURCE_COPY_KEYS) - config.keys()
    if missing_keys:
        errors.append('Missing reference settings: ' + ', '.join(sorted(missing_keys)))
    if errors:
        raise ValueError('\n'.join(errors))
    paths = required_paths(root, config)
    for path in paths:
        if not path.is_file():
            errors.append(f'Missing file: {path}')
            continue
        if path.stat().st_size == 0:
            errors.append(f'Empty file: {path}')
            continue
        expected = DETAIL_COLUMNS.get(path.name)
        if expected:
            with path.open(newline='') as stream:
                columns = next(csv.reader(stream, delimiter='\t'), [])
            missing = expected - set(columns)
            if missing:
                errors.append(f'Missing columns in {path}: {", ".join(sorted(missing))}')
    if errors:
        shown = '\n'.join(errors[:30])
        if len(errors) > 30:
            shown += f'\n... and {len(errors) - 30} more input errors.'
        raise ValueError('Incomplete current Motif run. Complete motif all, union_evaluation and '
                         'the saved cohort evaluations before exporting.\n' + shown)
    return paths


def write_source_index(paths, output):
    """Record actual source bytes, including full saved pair/triplet inputs."""
    with (Path(output) / 'source_file_index.tsv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=['source', 'bytes', 'sha256'], delimiter='\t')
        writer.writeheader()
        for path in paths:
            digest = hashlib.sha256()
            with path.open('rb') as source:
                for block in iter(lambda: source.read(1024 * 1024), b''):
                    digest.update(block)
            writer.writerow(dict(source=str(path), bytes=path.stat().st_size, sha256=digest.hexdigest()))


def validate_panels(panel_root, *, paper=None):
    """Require the completed current groups before exporting or analysing data."""
    paper = json.loads(PAPER.read_text()) if paper is None else paper
    paths = []
    for prefix, genes_key in [('manual', 'manual_prioritised'),
                              ('additional', 'additional_candidates')]:
        base = Path(panel_root) / paper[prefix + '_evaluation'] / 'results/manual_mfs'
        for name in PANEL_FILES:
            path = base / name
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f'Missing saved group result: {path}. Supply the completed '
                                 'section 2.9 and 2.10 server evaluations; no evaluation is run here.')
            paths.append(path.resolve())
        for name in ('selected_mfs/mfs.tsv', 'evaluation/binding_sites_per_mf.tsv'):
            with (base / name).open(newline='') as stream:
                rows = list(csv.DictReader(stream, delimiter='\t'))
            genes = [row['MF'] for row in rows]
            if len(genes) != len(set(genes)) or set(genes) != set(paper[genes_key]):
                raise ValueError(f'{base / name} does not match {genes_key}.')
        with (base / 'evaluation/ppi_proximity.tsv').open(newline='') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        if (len(rows) != 2 or {row['measure'] for row in rows} != {'harmonic', 'pagerank'}
                or any(int(row['n_mfs_nominated']) != len(genes) for row in rows)):
            raise ValueError(f'The saved network group results do not match {genes_key}.')
    return paths
