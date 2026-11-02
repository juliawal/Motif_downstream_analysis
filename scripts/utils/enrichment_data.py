"""Prepare the saved F3 results for the manuscript figure and enrichment table.

Retains the original statistics and the approved display categories and ordering.
No files are written and no hypothesis tests are run by this module.
"""
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
FOLD = 'observed_to_expected_term_overlap_ratio'
CATEGORIES = {
    'DNA_binding_transcription': ('DNA binding / transcription', 'DNA binding /\ntranscription', '#0072B2'),
    'development_cell_identity': ('Development / cell identity', 'Development /\ncell identity', '#D55E00'),
    'immune_cell_state': ('Immune signalling / cell state', 'Immune signalling /\ncell state', '#009E73'),
    'other': ('Other functions', 'Other\nfunctions', '#8064A2'),
}
# Explicit display-only refinement of the saved name-keyword classification.
# Existing categories remain in descriptive_theme; no source tables are edited.
IMMUNE_IDS = {
    'GO:0038093', 'GO:0038094', 'R-HSA-388841', 'R-HSA-2424491',
    'R-HSA-202433', 'R-HSA-389948', 'R-HSA-202427', 'R-HSA-202403',
    'R-HSA-202430',
}


def short_cohort(cohort):
    return cohort.replace('SKCM_METASTATIC', 'SKCM-MET').replace('SKCM_PRIMARY', 'SKCM-PRI')


def prepare(source):
    all_terms = pd.read_csv(source, sep='\t')
    union = all_terms[all_terms.analysis.eq('union') & all_terms.adj_p_value.le(.05)].copy()
    terms = union[union.term_has_min_mf_overlap].copy()
    terms['term_id'] = terms.term_name.str.extract(r'(GO:\d+|R-HSA-\d+)', expand=False)
    terms['display_category'] = terms.descriptive_theme
    terms.loc[terms.term_id.isin(IMMUNE_IDS), 'display_category'] = 'immune_cell_state'
    terms['category_label'] = terms.display_category.map(lambda x: CATEGORIES[x][0])
    terms['category_order'] = terms.display_category.map({k:i for i,k in enumerate(CATEGORIES)})
    rank = terms.groupby('cohort', as_index=False).agg(
        minimum_adjusted_p_value=('adj_p_value', 'min'), n_enriched_terms=('term_id', 'size'))
    rank = rank.sort_values(['minimum_adjusted_p_value', 'cohort']).reset_index(drop=True)
    rank['cohort_rank'] = np.arange(1, len(rank)+1)
    terms = terms.merge(rank[['cohort', 'cohort_rank']], on='cohort')
    terms = terms.sort_values(['cohort_rank', 'category_order', 'adj_p_value', 'library_name', 'term_id']).reset_index(drop=True)
    summary = terms.groupby(['cohort', 'display_category'], as_index=False).agg(
        n_terms=('term_id', 'size'), median_fold_enrichment=(FOLD, 'median'),
        minimum_fold_enrichment=(FOLD, 'min'), maximum_fold_enrichment=(FOLD, 'max'),
        minimum_term_adjusted_p_value=('adj_p_value', 'min'))
    summary = summary.merge(rank[['cohort', 'cohort_rank']], on='cohort')
    summary['category_order'] = summary.display_category.map({k:i for i,k in enumerate(CATEGORIES)})
    summary = summary.sort_values(['cohort_rank', 'category_order']).reset_index(drop=True)
    mapping = terms[['library_name','term_id','term_name','descriptive_theme','display_category','category_label']].drop_duplicates()
    mapping['assignment_basis'] = np.where(mapping.descriptive_theme.eq(mapping.display_category),
        'Saved descriptive category', 'Explicit immune-receptor/signalling term-ID assignment for display')
    mapping['reference_url'] = mapping.term_id.map(lambda x: 'https://reactome.org/content/detail/'+x if x.startswith('R-HSA') else 'https://amigo.geneontology.org/amigo/term/'+x)
    return terms, rank, summary, mapping

