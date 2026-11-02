#!/usr/bin/env python3
"""Summarise saved union evidence for the manual-candidate reference table.

Uses the availability flags supplied by integrate_evidence.py. Cohort-level
means weight each evaluable MF--cohort combination equally. Direction-consistency
counts weight each eligible MF equally; tested-cohort counts weight each repeated
MF--CpG pair equally. No evaluation is run here.
"""
import argparse
import csv
import gzip
from collections import defaultdict
import math
from pathlib import Path
from statistics import mean

ROOT = Path(__file__).resolve().parents[2]
FIELDS = ['unit', 'stage', 'measure', 'n_units', 'mean', 'display']
REQUIRED = {
    'analysis', 'cohort', 'MF_ensembl', 'f1_available', 'track', 'difference',
    'target_hit_rate', 'bg_hit_rate',
    'f2_available', 'harmonic_degree_bin_excess', 'harmonic_score',
    'harmonic_degree_bin_null_mean', 'pagerank_score', 'pagerank_degree_bin_null_mean',
    'n_tested_pairs', 'n_supported_pairs', 'n_a_supported_negative', 'n_a_supported_positive',
}


def iter_tsv(path):
    path = Path(path)
    with (gzip.open(path, 'rt') if path.suffix == '.gz' else path.open()) as stream:
        yield from csv.DictReader(stream, delimiter='\t')


def read_tsv(path):
    return list(iter_tsv(path))


def number(row, key):
    value = row[key]
    if value is None or str(value).strip().lower() in {'', 'na', 'nan', 'none'}:
        return None
    value = float(value)
    return value if math.isfinite(value) else None


def ratio(row, numerator, denominator):
    a, b = number(row, numerator), number(row, denominator)
    return a / b if a is not None and b is not None and b > 0 else None


def summarize_reference_means(rows, recurrence_rows=None):
    if not rows:
        raise ValueError('The union evidence table is empty.')
    missing = REQUIRED - rows[0].keys()
    if missing:
        raise ValueError(f'Missing union evidence columns: {sorted(missing)}')
    seen = set()
    for row in rows:
        key = (row['cohort'], row['MF_ensembl'])
        if row['analysis'] != 'union' or not all(key) or key in seen:
            raise ValueError(f'Expected one row per union gene--cohort combination: {key}')
        seen.add(key)
        for flag in ['f1_available', 'f2_available']:
            if row[flag] not in {'True', 'False'}:
                raise ValueError(f'Invalid saved {flag} for {key}: {row[flag]!r}')

    references = []

    def add(stage, label, values, suffix='', decimals=2, scale=1, offset=0, unit='MF--cohort'):
        values = [scale * (v + offset) for v in values if v is not None]
        average = mean(values) if values else None
        references.append(dict(unit=unit, stage=stage, measure=label,
                               n_units=len(values),
                               mean=average if average is not None else 'NA',
                               display=f'{average:.{decimals}f}{suffix}' if average is not None else '--'))

    for track, label in [('FIMO', 'Motif hit rate minus background rate'),
                         ('ChIP_seq', 'ChIP overlap rate minus background rate')]:
        selected = [r for r in rows if r['f1_available'] == 'True' and r['track'] == track]
        add('F1', label, [number(r, 'difference') for r in selected], ' pp', scale=100)
    motif = [r for r in rows if r['f1_available'] == 'True' and r['track'] == 'FIMO']
    for column, label in [('target_hit_rate', 'Motif hit rate at selected CpGs'),
                          ('bg_hit_rate', 'Motif hit rate at background CpGs')]:
        add('F1', label, [number(r, column) for r in motif], r'\%', scale=100)
    selected = [r for r in rows if r['f2_available'] == 'True']
    add('F2', 'Harmonic score: percentage above matched mean',
        [ratio(r, 'harmonic_score', 'harmonic_degree_bin_null_mean') for r in selected],
        r'\%', scale=100, offset=-1)
    add('F2', 'PageRank: percentage above matched mean',
        [ratio(r, 'pagerank_score', 'pagerank_degree_bin_null_mean') for r in selected],
        r'\%', scale=100, offset=-1)
    testable = [r for r in rows if number(r, 'n_tested_pairs') is not None
                and number(r, 'n_tested_pairs') > 0]
    for column, label in [('n_tested_pairs', 'Tested MF--CpG pairs'),
                          ('n_supported_pairs', 'Supported MF--CpG pairs')]:
        add('F4', label, [number(r, column) for r in testable], decimals=1)
    add('F4', 'Supported pairs among tested pairs',
        [ratio(r, 'n_supported_pairs', 'n_tested_pairs') for r in rows], r'\%', 1, scale=100)
    for sign in ['negative', 'positive']:
        add('F4', sign.capitalize() + ' directions among supported pairs',
            [ratio(r, 'n_a_supported_' + sign, 'n_supported_pairs') for r in rows], r'\%', 1, scale=100)
    if recurrence_rows is not None:
        add_direction_references(recurrence_rows, rows, add)
    return references


def add_direction_references(recurrence_rows, contexts, add):
    # Every eligible MF counts once, including zero consistent-direction pairs.
    # Eligibility requires at least one exact pair tested in >=2 cohorts.
    genes = {r['MF_ensembl'] for r in contexts}
    counts = defaultdict(lambda: dict(total=0, negative=0, positive=0, mixed=0))
    tested_cohorts = []
    seen = set()
    for row in recurrence_rows:
        if int(row['n_tested_cohorts']) < 2:
            continue
        key = (row['MF_ensembl'], row['CpG'])
        if key[0] not in genes or key in seen:
            raise ValueError(f'Unexpected or duplicate recurrent-pair record: {key}')
        seen.add(key)
        gene = counts[key[0]]
        consistent = row['same_sign_all_tested'] == 'True'
        gene['total'] += consistent
        gene['negative'] += row['direction_all_tested'] == 'negative'
        gene['positive'] += row['direction_all_tested'] == 'positive'
        gene['mixed'] += row['direction_all_tested'] == 'mixed'
        tested_cohorts.append(int(row['n_tested_cohorts']))
    for key, label in [('total', 'Pairs with consistent association directions'),
                       ('negative', 'Pairs with negative directions throughout'),
                       ('positive', 'Pairs with positive directions throughout'),
                       ('mixed', 'Pairs with opposite directions across tested cohorts')]:
        add('F4', label, [r[key] for r in counts.values()], unit='MF')
    add('F4', 'Cohorts in which the same MF--CpG pair was tested', tested_cohorts,
        unit='MF--CpG')


def build_reference_means(source, output, recurrence_source=None):
    recurrence_rows = iter_tsv(recurrence_source) if recurrence_source is not None else None
    references = summarize_reference_means(read_tsv(source), recurrence_rows)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=FIELDS, delimiter='\t')
        writer.writeheader()
        writer.writerows(references)
    return references


def render_reference_table(source, output):
    references = read_tsv(source)
    body = []
    previous_unit = None
    for row in references:
        if row['unit'] != previous_unit:
            if previous_unit is not None:
                body.append(r'\addlinespace')
            body.append(r'\multicolumn{4}{@{}l}{\textit{Per ' + row['unit']
                        + (' combination' if row['unit'] == 'MF--cohort' else ' pair' if row['unit'] == 'MF--CpG' else '')
                        + r'}}\\')
            previous_unit = row['unit']
        body.append(f"{row['stage']} & {row['measure']} & {row['n_units']} & {row['display']} " + r'\\')
    table = r'''% Generated by scripts/build_tables.py from saved reference means.
\begin{table}[!htbp]
\centering
\small
\caption{Reference means for the measures reported in \Cref{sec:results-manual17,sec:results-candidate-priorities}.}\label{tab:manual17-reference}
\setlength{\tabcolsep}{4pt}
\renewcommand{\arraystretch}{1.15}
\begin{tabularx}{\linewidth}{@{}lXrr@{}}
\toprule
Stage & Measure & $n$ & Mean \\
\midrule
''' + '\n'.join(body) + r'''
\bottomrule
\end{tabularx}
\par\smallskip
\begin{minipage}{\linewidth}
\small The reference means summarise the available union results across all candidate MFs in the 34 cohorts, including the candidates discussed in \Cref{sec:results-manual17,sec:results-candidate-priorities}. They provide a comparison for the individual candidate results. Means are calculated across MF--cohort combinations, individual MFs or MF--CpG pairs, as indicated. $n$ is the number of evaluable units for each measure; each unit is weighted equally.
\par\smallskip
For MF--cohort means, the same MF is counted separately in each cohort. Percentages are calculated per combination before averaging. Missing measurements are excluded. F4 pair counts and support percentages require at least one tested pair; direction percentages require at least one supported pair.
\par\smallskip
The per-MF direction counts include MFs with at least one pair tested in two or more cohorts, including MFs with no consistent pairs. A pair has a consistent direction if its coefficient $a$ has the same nonzero sign in every tested cohort, regardless of F4 support. The final row reports the mean number of tested cohorts across all pairs tested in two or more cohorts. pp: percentage points; matched mean: mean score of degree-matched background proteins.
\end{minipage}
\end{table}
'''
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(table)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=ROOT / 'analysis_results/union_evidence_contexts.tsv')
    parser.add_argument('--output', type=Path, default=ROOT / 'analysis_results/manual_candidate_reference_means.tsv')
    parser.add_argument('--recurrence', type=Path, default=ROOT / 'analysis_results/pair_directions/all_exact_pair_recurrence.tsv.gz')
    args = parser.parse_args()
    rows = build_reference_means(args.input, args.output, args.recurrence)
    print(f'Wrote {len(rows)} reference means from saved union evidence: {args.output}')


if __name__ == '__main__':
    main()
