#!/usr/bin/env python3
"""Prepare descriptive manuscript results from the audited MF--CpG exports.

No model fitting, new p-values, or inferred experimental effects. Priority groups
are qualitative editorial choices, with limitations retained explicitly.
"""
from pathlib import Path
import pandas as pd
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
DIRECTIONS = ROOT / 'analysis_results/pair_directions'
CANDIDATES = ROOT / 'analysis_results/additional_candidates'
read = lambda path: pd.read_csv(path, sep='\t', low_memory=False)
write = lambda frame, path: frame.to_csv(path, sep='\t', index=False, na_rep='NA')

# An editorial shortlist, not an estimated probability of being a regulator.
# The three descriptive groups do not imply a statistically calibrated ordering.
curation = [
 ('LMX1B','Main','Recurrent positive associations','Wiederkehrende positive Paare',
  '36 recurrent pairs; four favourable layers in LUAD.',
  'Small LUAD binding excess; LUSC binding difference negative; no LUSC S3 support.'),
 ('POU2AF2','Main','Recurrent negative associations','Wiederkehrende negative Paare',
  '7 recurrent pairs, including a pair supported in every setting in CESC and LUAD.',
  'No usable binding measurement.'),
 ('ZNF256','Main','Both CpG directions','Beide CpG-Richtungen',
  '3 recurrent pairs; positive binding and network differences in STAD.',
  'Mixed directions; only 5/36 supported pairs in COAD.'),
 ('HOXB6','Main','Strong single-cohort evidence','Starker Einzelkohortenbefund',
  '105/125 supported in COAD; 82 pairs supported in union and S1--S3.',
  'No recurrent supported pair.'),
 ('ZNF518B','Main','Recurrent predominantly negative associations','Wiederkehrende überwiegend negative Paare',
  '6 recurrent pairs; mostly negative in COAD and ESCA.',
  'No usable binding measurement; no COAD pair supported in all settings.'),
 ('LHX5','Main','Context-dependent integrated evidence','Kontextabhängige kombinierte Evidenz',
  'Four favourable layers in every setting in LGG; 35 pairs supported in all settings.',
  'One supported exact pair reverses direction in PAAD.'),
 ('FOXD1','Main','Context-dependent integrated evidence','Kontextabhängige kombinierte Evidenz',
  'Four favourable layers in every setting in LGG; 281/350 supported pairs.',
  'No recurrent supported pair; small binding and network differences.'),
 ('PITX2','Main','Integrated evidence with recurrence','Kombinierte Evidenz mit Paarrekurrenz',
  '18 recurrent positive pairs; all four layers in every setting in three contexts.',
  'Primary/metastatic SKCM may not be independent; 12 recurrent pairs also link SARC to SKCM-MET.'),
 ('RIMBP2','Alternative','Strong pair recurrence','Starke Paarrekurrenz',
  '34 recurrent pairs; cg06403845 supported in four cohorts and every setting.',
  'No positive network or functional layer.'),
 ('KNDC1','Alternative','Strong pair recurrence','Starke Paarrekurrenz',
  '65 recurrent pairs; predominantly positive; network excess in BRCA and LIHC.',
  'No usable binding evidence or favourable functional layer; LIHC absent in S3.'),
 ('BEX2','Alternative','Broad recurrent negative signal','Breites wiederkehrendes negatives Signal',
  '733 recurrent pairs between ESCA and STAD.',
  'Broad, almost uniform sign; no usable binding or favourable F3.'),
 ('HOXC9','Reserve','Strong LUAD context','Starker LUAD-Kontext',
  '125/158 supported pairs and four favourable layers in LUAD.',
  'No recurrent supported pair and none supported in all settings.'),
 ('TBX1','Reserve','Site-dependent directions','Zielstellenabhängige Richtungen',
  '6 recurrent pairs; distinct exact pairs show positive and negative directions.',
  'Some favourable contexts contain few supported pairs; no global MF sign.'),
 ('PBX1','Reserve','Strong melanoma context','Starker Melanom-Kontext',
  'Four favourable layers in every setting in SKCM-MET; 232/480 supported pairs.',
  'No recurrent supported pair; small binding difference.'),
 ('TFAP2B','Reserve','Integrated BRCA context','Kombinierte Evidenz in BRCA',
  'Four favourable layers in BRCA; 45/100 supported pairs.',
  'No recurrent supported pair; only 6 pairs supported in all settings; mixed signs.'),
 ('INSM1','Reserve','Predominantly positive PRAD context','Überwiegend positiver PRAD-Kontext',
  '223/311 supported pairs; positive binding and network differences.',
  'No recurrent supported pair.'),
]
curated = pd.DataFrame(curation, columns=['MF','display_group','evidence_profile','evidence_profile_de',
    'data_rationale','main_limitation'])
curated['display_order'] = np.arange(len(curated))
genes = read(CANDIDATES/'additional_candidate_evidence.tsv')
short = curated.merge(genes, on='MF')
r = read(DIRECTIONS/'all_exact_pair_recurrence.tsv.gz')
write(short, CANDIDATES/'candidate_shortlist.tsv')
contexts = read(CANDIDATES/'union_context_evidence.tsv').merge(curated[['MF','display_order','display_group']], on='MF')
contexts = contexts.sort_values(['display_order','cohort'])
write(contexts, CANDIDATES/'candidate_contexts.tsv')
main_contexts = contexts[contexts.display_group.eq('Main') & contexts.n_tested_pairs.gt(0)].copy()
write(main_contexts.merge(short[['MF', 'n_strict_recurrent_pairs']], on='MF', how='left'),
      CANDIDATES/'candidate_evidence.tsv')

repeated = r[r.n_tested_cohorts.ge(2)].copy()
repeated['direction_all_tested'] = np.select(
    [repeated.n_negative.eq(0), repeated.n_positive.eq(0)], ['positive','negative'], default='mixed')
write(repeated, DIRECTIONS/'all_repeated_pairs.tsv')

d = read(DIRECTIONS/'all_pair_cohort_evidence.tsv.gz')
# Choose one PITX2 example requiring supported SARC and SKCM-MET, not solely
# primary/metastatic melanoma. Prioritise support in all settings in both.
pit = d[d.MF.eq('PITX2') & d.supported & d.cohort.isin(['SARC','SKCM_METASTATIC'])]
pc = pit.groupby('CpG').agg(n=('cohort','nunique'),all_settings=('supported_all_four','sum'),min_abs_a=('a',lambda x:x.abs().min()))
pc = pc[pc.n.eq(2)].sort_values(['all_settings','min_abs_a'],ascending=False)
pit_cpg = pc.index[0]
plot_pairs = [('LMX1B','cg14325930'),('POU2AF2','cg03078767'),('ZNF256','cg00541956'),
    ('ZNF256','cg03430717'),('HOXB6','cg22660299'),('ZNF518B','cg07850527'),
    ('PITX2',pit_cpg),('FOXD1','cg04968473'),('LHX5','cg26843324')]
key = pd.DataFrame(plot_pairs,columns=['MF','CpG']);key['pair_order']=np.arange(len(key))
plot = d.merge(key,on=['MF','CpG']).merge(
    r[['MF_ensembl','CpG','n_tested_cohorts','n_supported_cohorts','strict_recurrent_supported']],on=['MF_ensembl','CpG'])
write(plot.sort_values(['pair_order','cohort']), DIRECTIONS/'plotted_pair_directions.tsv')
print(f'Prepared {len(short)} candidates, {len(main_contexts)} plotted contexts and {len(plot)} pair rows.')
