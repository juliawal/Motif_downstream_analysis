#!/usr/bin/env python3
"""Prepare descriptive manuscript results from the audited MF--CpG exports.

No model fitting, new p-values, or inferred experimental effects. Priority groups
are qualitative editorial choices, with limitations retained explicitly.
"""
import argparse
import json
from pathlib import Path
import pandas as pd
import numpy as np
from prepare_direction_results import prepare_direction_results

ROOT = Path(__file__).resolve().parents[2]
DIRECTIONS = ROOT / 'analysis_results/pair_directions'
CANDIDATES = ROOT / 'analysis_results/additional_candidates'
read = lambda path: pd.read_csv(path, sep='\t', low_memory=False)
write = lambda frame, path: frame.to_csv(path, sep='\t', index=False, na_rep='NA')

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--candidates-only', action='store_true',
                    help='Update candidate tables without recalculating the separate direction-only analysis')
args = parser.parse_args()

# An editorial shortlist, not an estimated probability of being a regulator.
# Strong individual findings and cross-cohort evidence are both represented.
curation = [
 ('SIM2','Main','Evidence across tumour types','Evidenz über Tumorarten',
  'F2--F4 in four cohorts; 32 recurrent pairs, including 12 supported in at least three cohorts.',
  'No usable binding measurement; LUAD absent in S3; F3 not met in every setting.'),


 ('PITX2','Main','Integrated evidence with recurrence','Kombinierte Evidenz mit Paarrekurrenz',
  '18 recurrent positive pairs; all four criteria in every setting in SARC and both SKCM cohorts.',
  'Only 11/532 supported in SKCM-PRI; PageRank below background; 12 recurrent pairs link SARC and SKCM-MET.'),
 ('LMX1B','Main','Recurrent positive associations','Wiederkehrende positive Paare',
  '36 recurrent positive pairs; 266/507 supported in LUAD and 1236/1754 in LUSC; F2--F4 in both.',
  'Small LUAD binding excess; LUSC binding difference negative; LUSC absent in S3.'),
 ('KNDC1','Main','Strong recurrent positive associations','Starke rekurrente positive Assoziationen',
  '65 recurrent pairs; 2758/3187 supported in BRCA and 127/152 in LIHC; both network differences positive.',
  'No usable binding measurement or favourable F3; LIHC absent in S3.'),
 ('BEX2','Main','Broad recurrent negative associations','Breite rekurrente negative Assoziationen',
  '733 recurrent pairs between ESCA and STAD; 154 supported in both cohorts in all settings.',
  '733/2740 repeatedly tested pairs recurrent; no usable F1 or favourable F3; STAD PageRank below background.'),

 ('LHX5','Main','Context-dependent integrated evidence','Kontextabhängige kombinierte Evidenz',
  'All four criteria in every setting in LGG; 77/101 supported pairs; 35 in all settings.',
  'Only 1/71 supported in PAAD; one exact pair supported with opposite signs in LGG and PAAD.'),
 ('DLX6','Main','Evidence across several cohorts','Evidenz über mehrere Kohorten',
  'F2--F4 in HNSC, LUAD and SARC; additional F4 support in ESCA; both network differences positive.',
  'No motif hits; only 3/51 repeatedly tested pairs recurrent; ESCA selected only in S3.'),
 ('RIMBP2','Main','Pair recurrence across four cohorts','Paarrekurrenz über vier Kohorten',
  '34 recurrent pairs; cg06403845 supported in four cohorts and every setting.',
  'Binding and network differences negative; no favourable F3.'),
 ('PLEKHG4B','Main','High recurrence among repeatedly tested pairs','Hohe Rekurrenz wiederholt getesteter Paare',
  'All 26 repeatedly tested pairs recurrent; 32/41 supported in HNSC and 28/29 in LUSC.',
  'No usable F1; network differences negative and no favourable F3; limited pair numbers.'),
 ('CELSR3', 'Main', 'Recurrent associations across five cohorts', 'Rekurrente Assoziationen über fünf Kohorten', 'F2 and F4 in five cohorts; 55/86 repeatedly tested pairs recurrent; cg22702772 supported positively in all five.', 'No usable F1; PageRank below background in all five cohorts; no F3 membership.'),
 ('HNF1B', 'Main', 'Integrated evidence in UCEC', 'Kombinierte Evidenz in UCEC', 'All four criteria in the union and S1--S3; motif excess 24.6 percentage points; 490/929 supported pairs.', 'Selected in one cohort; same patients across settings; no demonstrated DNMT/TET effect in the cited study.'),
 ('HOXA9', 'Main', 'Integrated evidence with recurrence', 'Kombinierte Evidenz mit Rekurrenz', 'All four criteria in the union and S1--S3 in UCEC; 20/38 repeatedly tested pairs recurrent.', 'Only 1/29 supported pairs in PCPG; HNSC binding difference negative; cited mechanism concerns histone methylation.'),
 ('HOXB13', 'Main', 'Integrated evidence across cohorts', 'Kombinierte Evidenz über Kohorten', 'All four criteria in five union cohorts and all settings in four; 86/201 repeatedly tested pairs recurrent.', 'Only 10/1618 supported in SKCM-P; the cited study demonstrates methylation-sensitive DNA binding.'),
 ('HOXB9', 'Main', 'Integrated evidence in UCEC', 'Kombinierte Evidenz in UCEC', 'All four criteria in the union and S1--S3; 86/126 supported pairs; both network scores above background.', 'Selected in one cohort; the cited methylation association concerns its own locus.'),
 ('IRX2', 'Main', 'Evidence across several cohorts', 'Evidenz über mehrere Kohorten', 'F2--F4 in four cohorts, F2/F4 in two more; 10/42 repeatedly tested pairs recurrent.', 'No usable F1; harmonic excess only 0.9--1.2%; five pairs have opposite directions in supported cohorts.'),
 ('MNX1', 'Main', 'Multi-stage evidence with recurrence', 'Mehrstufige Evidenz mit Rekurrenz', 'F2--F4 in five cohorts and F2/F4 in KIRP; 54/95 repeatedly tested pairs recurrent.', 'No motif hits; BRCA absent in S3; cited experiments concern histones rather than a DNMT/TET mechanism.'),
 ('PTPRN2', 'Main', 'Recurrent associations across settings', 'Rekurrente Assoziationen über Settings', '157/349 repeatedly tested pairs recurrent; 132 supported in at least two cohorts in the union and every setting.', 'Only F4 met; no usable F1, both network measures below background, and no F3 membership.'),
]
curated = pd.DataFrame(curation, columns=['MF','display_group','evidence_profile','evidence_profile_de',
    'data_rationale','main_limitation'])
expected = json.loads((ROOT / 'scripts/paper_cohorts_and_genes.json').read_text())['additional_candidates']
if not curated.MF.is_unique or set(curated.MF) != set(expected):
    raise ValueError('The curated shortlist must match additional_candidates in the paper configuration.')
curated = curated.sort_values('MF').reset_index(drop=True)
curated['display_order'] = np.arange(len(curated))
genes = read(CANDIDATES/'all_candidate_evidence.tsv')
short = curated.merge(genes, on='MF', validate='one_to_one')
assert len(short) == len(curated), 'Missing curated candidate evidence'
write(short, CANDIDATES/'candidate_shortlist.tsv')
contexts = read(CANDIDATES/'union_context_evidence.tsv').merge(curated[['MF','display_order','display_group']], on='MF')
contexts = contexts.sort_values(['display_order','cohort'])
write(contexts, CANDIDATES/'candidate_contexts.tsv')
main_contexts = contexts[contexts.display_group.eq('Main') & contexts.n_tested_pairs.gt(0)].copy()
write(main_contexts.merge(short[['MF', 'n_strict_recurrent_pairs']], on='MF', how='left'),
      CANDIDATES/'candidate_evidence.tsv')


if not args.candidates_only:
    prepare_direction_results(DIRECTIONS, ROOT / 'analysis_results')
print(f'Prepared {len(short)} candidates and {len(main_contexts)} plotted contexts.')
