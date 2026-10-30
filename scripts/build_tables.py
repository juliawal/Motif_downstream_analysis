"""Render the 13 supplementary table bodies from analysis_results/.

Captions, column layouts and explanations belong to the manuscript. Numbers,
selection order and the curated candidate descriptions come from the analysis and table builders.
All current manuscript tables belong to the supplement; there are no main-text tables.
Only the 89 prepared datasets are needed, not the other analysis intermediates.
"""
import csv
import math
import re
from pathlib import Path

TABLES = [
    ('tab:rerun-coverage', 'Coverage'),
    ('tab:rerun-identity', 'Identity'),
    ('tab:f3-supported-terms', 'FthreeTerms'),
    ('tab:rerun-f4-coverage', 'FfourCoverage'),
    ('tab:rerun-cross-layer', 'CrossLayer'),
    ('tab:rerun-f4-overlap', 'Overlap'),
    ('tab:rerun-erg-triplets', 'ErgTriplets'),
    ('tab:rerun-hnf1b-direction', 'HnfDirections'),
    ('tab:rerun-hnf1b-context', 'HnfContext'),
    ('tab:rerun-case-f4', 'CaseFfour'),
    ('tab:rerun-case-criteria', 'CaseCriteria'),
    ('tab:candidate-priorities', 'CandidatePriorities'),
    ('tab:candidate-pair-counts', 'CandidateCounts'),
]
ARMS = {"rank_1": "S1", "rank_2": "S2", "rank_3": "S3", "union": "Union"}


def integer(value):
    return f"{int(float(value)):,}".replace(",", r"\,")


def yes(value):
    if isinstance(value, str):
        value = value == "True"
    return "Yes" if value else "No"


def fixed(value, digits=4, scale=1):
    if value == "" or not math.isfinite(float(value)):
        return "--"
    return f"{float(value) * scale:.{digits}f}"


def pvalue(value, significant=False):
    value = float(value)
    if value == 1:
        return "1"
    if 0 < value < .001:
        mantissa, exponent = f"{value:.2e}".split('e')
        return rf"${mantissa}\times10^{{{int(exponent)}}}$"
    return f"{value:.3g}" if significant else f"{value:.4f}"


def cohort(value):
    return value.replace("SKCM_METASTATIC", "SKCM-MET").replace("SKCM_PRIMARY", "SKCM-PRI")


def escape(value):
    replacements = {'&': r'\&', '%': r'\%', '_': r'\_', '#': r'\#', '$': r'\$'}
    return ''.join(replacements.get(c, c) for c in str(value))


def render(root):
    root = Path(root)
    data = root / 'analysis_results'
    def read(path):
        with path.open() as handle:
            return list(csv.DictReader(handle, delimiter='\t'))
    def source(name):
        return read(data / (name + '.tsv'))
    def ratio(row, a, b):
        return integer(row[a]) + '/' + integer(row[b])
    def lines(rows):
        return ''.join(' & '.join(row) + r' \\' + '\n' for row in rows)
    bodies = {}
    bodies['Coverage'] = lines([
        [cohort(r['cohort']), integer(r['n_union_selected']), integer(r['n_in_all_three_settings']),
         integer(r['n_formal_F1']), ratio(r,'n_usable_FIMO','n_usable_ChIP'), integer(r['n_PPI_mapped']),
         integer(r['n_F4_samples']), yes(r['F4_testable'])]
        for r in source('union_evaluation_coverage')])
    bodies['FfourCoverage'] = lines([
        [cohort(r['cohort']), ARMS[r['analysis']], integer(r['n_input_pairs']), integer(r['n_candidate_pairs']),
         ratio(r,'n_supported_pairs','n_tested_pairs') if r['testable']=='True' else 'NT',
         ratio(r,'n_supported_triplets','n_tested_triplets') if r['testable']=='True' else 'NT',
         integer(r['n_untestable_triplets']), fixed(r['pair_support_fraction'],2,100), fixed(r['triplet_support_fraction'],2,100)]
        for r in source('association_coverage')])
    bodies['CrossLayer'] = lines([
        [r['MF'],cohort(r['cohort']),r['track'].replace('ChIP_seq','ChIP'), ratio(r,'n_target_hits','n_target_cpgs'),
         ratio(r,'n_background_hits','n_background_cpgs'), fixed(r['harmonic_degree_bin_excess']), integer(r['n_significant_min3_terms']),
         ratio(r,'n_supported_pairs','n_tested_pairs'),
         yes(all(r[k]=='True' for k in ['selected_all_three_settings','f1_positive_all_three_settings','f3_membership_all_three_settings','f4_support_all_three_settings'])),
         yes(r['all_four_positive_in_union_and_all_settings'])]
        for r in source('combined_positive_contexts')])
    bodies['Identity'] = lines([
        [r['MF'],cohort(r['cohort']),ARMS[r['analysis']]] + [yes(r[k]) for k in [
          'f1_confirmed_motif_mismatch','f1_unresolved_chip_identity','f1_available',
          'f1_available_including_unresolved_chip','all_four_available_including_unresolved_chip','all_four_positive_including_unresolved_chip']]
        for r in source('binding_identity_sensitivity')])
    bodies['ErgTriplets'] = lines([
        [ARMS[r['analysis']],r['CpG'],integer(r['G_prime_distance_bp']),fixed(r['a_estimate'],3),pvalue(r['a_pvalue']),
         fixed(r['b_estimate'],3),pvalue(r['b_pvalue']),pvalue(r['joint_adj_p_value']),yes(r['mediation_supported']),yes(r['binding_hit'])]
        for r in source('erg_tdrd1_triplets')])
    bodies['CaseFfour'] = lines([
        [r['MF'],cohort(r['cohort']),ARMS[r['analysis']],integer(r['n_input_pairs']),ratio(r,'n_supported_pairs','n_tested_pairs'),
         ratio(r,'n_supported_triplets','n_tested_triplets'),integer(r['n_untestable_triplets']),fixed(r['a_supported_median'],3),
         ratio(r,'n_a_supported_negative','n_supported_pairs'),fixed(r['pair_support_fraction'],2,100)]
        for r in source('case_association_table')])
    bodies['CaseCriteria'] = lines([
        [r['MF'],cohort(r['cohort']),ARMS[r['analysis']],fixed(r['difference']) if r['f1_available']=='True' else 'NT',fixed(r['harmonic_degree_bin_excess']),
         integer(r['n_significant_min3_terms']),yes(r['f4_supported']),yes(r['all_four_available']),yes(r['all_four_positive'])]
        for r in source('case_criteria_table')])
    bodies['HnfDirections'] = lines([
        [ARMS[r['analysis']], {'tested':'All tested','supported':'Supported'}[r['population']],
         ratio(r,'n_negative_a_hit','n_a_available_hit'),ratio(r,'n_negative_a_no_hit','n_a_available_no_hit'),
         fixed(r['negative_a_fraction_hit'],1,100),fixed(r['negative_a_fraction_no_hit'],1,100),
         fixed(float(r['negative_a_fraction_hit'])-float(r['negative_a_fraction_no_hit']),1,100)]
        for r in source('hnf1b_binding_directions')])
    bodies['HnfContext'] = lines([
        [{'cgi':'CpG-island context','region':'Gene region','cgi_x_region':'Context + region'}[r['stratification']],
         integer(r['n_common_cells']),ratio(r,'n_retained_binding','n_total_binding'),ratio(r,'n_retained_no_binding','n_total_no_binding'),
         fixed(r['standardized_difference'],1,100)] for r in source('hnf1b_genomic_context_sensitivity')])
    bodies['Overlap'] = lines([[r['key'],escape(r['description']),integer(r['rows'])]
         for r in source('overlap_index')])
    term_rows=[]; last_cohort=None
    for r in source('union_enrichment_terms'):
        if last_cohort is not None and last_cohort != r['cohort']:
            term_rows.append('\\midrule\n')
        last_cohort=r['cohort']
        name=re.sub(r'\s*\(?'+re.escape(r['term_id'])+r'\)?\s*$','',r['term_name'])
        library={'GO_Molecular_Function_2023':'GO MF','GO_Biological_Process_2023':'GO BP','GO_Cellular_Component_2023':'GO CC'}.get(r['library_name'],'Reactome')
        term_rows.append(lines([[r['cohort'].replace('SKCM_METASTATIC','SKCM met.').replace('SKCM_PRIMARY','SKCM primary'),
            escape(r['category_label']),escape(name)+r'\newline '+r['term_id'],library,
            str(int(r['n_mf_symbols_in_term']))+'/'+str(int(r['n_mf_symbols'])),
            f"{int(r['n_term_genes']):,}/{int(r['n_background_genes']):,}",
            fixed(r['observed_to_expected_term_overlap_ratio'],2),pvalue(r['adj_p_value'],significant=True),
            ', '.join(sorted(set(r['mf_symbols_in_term'].split(';'))))]]))
    bodies['FthreeTerms']=''.join(term_rows)
    candidates=source('candidate_shortlist')
    bodies['CandidatePriorities']=lines([[r'\textbf{'+r['MF']+r'}\newline '+r['display_group'],
        escape(r['data_rationale']),escape(r['main_limitation'])] for r in candidates])
    bodies['CandidateCounts']=lines([[r['MF']]+[integer(r[k]) for k in ['n_tested_pairs','n_supported_pairs',
        'n_supported_negative','n_supported_positive','n_supported_all_four_pairs','n_strict_recurrent_pairs']] for r in candidates])
    return bodies


def build(root):
    bodies=render(root)
    output=Path(root)/'tables/supplements/manuscript_rows.tex'
    output.parent.mkdir(parents=True,exist_ok=True)
    content='% Generated by scripts/build_tables.py; edit source analyses, not these rows.\n'
    for _,name in TABLES:
        content += '\\newcommand{\\PaperTable'+name+'}{%\n'+bodies[name]+'}\n'
    output.write_text(content)
    print(f'Wrote all {len(TABLES)} manuscript table bodies: {output}')


if __name__ == '__main__':
    build(Path(__file__).resolve().parents[1])
