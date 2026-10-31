#!/usr/bin/env python3
"""Summarize the manual MF panel as eight manuscript data tables.

Read saved F1–F3 pipeline results and the main analysis's cohort/F4 summaries.
Keep original p/q-values; the F4 summaries describe the existing cohort runs,
not a separate F4 test of the manual panel. No archive or report files are needed.
"""
import argparse
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "utils"))
from pipeline_file_layout import PipelineFiles
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=ROOT / "pipeline_results", help="Pipeline export directory")
    parser.add_argument("--analysis", type=Path, default=ROOT / "analysis_results", help="Analysis result root containing the F4 and joint evaluation folders")
    parser.add_argument("--output", type=Path, default=ROOT / "analysis_results/manual_candidates", help="Directory for the eight manual-panel data tables")
    args = parser.parse_args()
    raw = args.input.resolve() / "original_pipeline_outputs/manual_mfs"
    files = PipelineFiles(args.input)
    read = lambda path: pd.read_csv(files.path(path) if Path(path).is_absolute() else path, sep="\t")
    evidence = read(raw / "evaluation/binding_sites_per_mf.tsv.gz").sort_values("MF")
    terms = read(raw / "evaluation/gene_set_enrichment.tsv.gz")
    reportable = terms[terms.adj_p_value.lt(.05) & terms.term_has_min_mf_overlap]
    reportable = reportable.sort_values("adj_p_value", kind="stable")
    # A missing binding measurement must not appear as a measured zero.
    evidence.loc[evidence.status.ne("tested"), ["target_hit_rate", "bg_hit_rate"]] = np.nan
    evidence["difference_percentage_points"] = 100 * (evidence.target_hit_rate - evidence.bg_hit_rate)
    memberships = reportable.mf_symbols_in_term.str.split(";").explode().value_counts()
    evidence["n_reportable_F3_terms"] = evidence.MF.map(memberships).fillna(0).astype(int)
    f2 = read(raw / "evaluation/ppi_proximity.tsv.gz")
    f2g = f2.set_index('measure')
    f1, f3 = evidence, reportable
    genes = f1.MF.tolist()
    all_context = read(args.analysis / 'combined_evidence/union_evidence_contexts.tsv')
    contexts = all_context.loc[all_context.MF.isin(genes)].sort_values(['MF','cohort']).reset_index(drop=True)
    for measure in ['harmonic','pagerank']:
        excess = contexts[measure+'_degree_bin_excess']
        contexts[measure+'_difference_sign'] = np.select([excess.gt(0), excess.lt(0)], ['positive', 'negative'], default='zero')
    contexts['f4_support_percent']=np.where(contexts.n_tested_pairs>0,100*contexts.n_supported_pairs/contexts.n_tested_pairs,np.nan)
    contexts['f4_positive_a_percent']=np.where(contexts.n_tested_pairs>0,100*contexts.n_a_supported_positive/contexts.n_tested_pairs,np.nan)
    contexts['f4_negative_a_percent']=np.where(contexts.n_tested_pairs>0,100*contexts.n_a_supported_negative/contexts.n_tested_pairs,np.nan)
    contexts['f4_not_supported_percent']=np.where(contexts.n_tested_pairs>0,100*(contexts.n_tested_pairs-contexts.n_supported_pairs)/contexts.n_tested_pairs,np.nan)
    spec=[('GO:0000987','Cis-regulatory region sequence-specific DNA binding'),
          ('GO:0006355','Regulation of DNA-templated transcription'),
          ('GO:1990837','Sequence-specific double-stranded DNA binding'),
          ('GO:0001067','Transcription regulatory region nucleic acid binding'),
          ('GO:0009952','Anterior/posterior pattern specification'),
          ('GO:0045596','Negative regulation of cell differentiation')]
    terms=[]
    for go,label in spec:
        row=f3.loc[f3.term_name.str.contains(go,regex=False)].iloc[0].copy()
        row['plot_label']=label;row['go_id']=go;terms.append(row)
    terms=pd.DataFrame(terms).reset_index(drop=True)
    terms['plot_id']=[f'T{i+1}' for i in range(6)]
    term_selection = f3.copy()
    term_selection['go_id'] = term_selection.term_name.str.extract(r'(GO:\d+)')[0]
    term_selection['plot_id'] = term_selection.go_id.map(dict(zip(terms.go_id, terms.plot_id))).fillna('')
    term_selection['shown_in_figure'] = term_selection.plot_id.ne('')
    term_selection['overlap_rank_among_18'] = term_selection.n_mf_symbols_in_term.rank(method='min', ascending=False).astype(int)
    term_selection = term_selection.sort_values(['n_mf_symbols_in_term','adj_p_value'], ascending=[False,True])
    matrix=np.array([[g in r.mf_symbols_in_term.split(';') for _,r in terms.iterrows()] for g in genes],int)

    f2cols=['MF','MF_ensembl','cohort','entrez_id','pagerank_score','pagerank_degree_bin_null_mean','pagerank_degree_bin_excess','pagerank_difference_sign','harmonic_score','harmonic_degree_bin_null_mean','harmonic_degree_bin_excess','harmonic_difference_sign']
    f4cols=['MF','MF_ensembl','cohort','n_samples','context_status','n_tested_pairs','n_supported_pairs','n_a_supported_negative','n_a_supported_positive','f4_support_percent','f4_negative_a_percent','f4_positive_a_percent','f4_not_supported_percent']

    tables = args.output
    tables.mkdir(parents=True, exist_ok=True)
    for frame, name, index in [
        (f1, 'manual_binding_measurements.tsv', False),
        (read(raw / 'evaluation/binding_sites.tsv.gz'), 'manual_binding_group_tests.tsv', False),
        (contexts[f2cols], 'manual_network_contexts.tsv', False),
        (f2g, 'manual_network_group_tests.tsv', True),
        (terms, 'manual_enrichment_terms.tsv', False),
        (term_selection, 'manual_enrichment_term_selection.tsv', False),
        (pd.DataFrame(matrix, index=genes, columns=terms.plot_id).rename_axis('MF'), 'manual_enrichment_memberships.tsv', True),
        (contexts[f4cols], 'manual_association_contexts.tsv', False),
    ]:
        frame.to_csv(tables / name, sep="\t", index=index)
    print(f"Wrote eight manuscript data tables for {len(genes)} manual MFs to {tables}")


if __name__ == "__main__":
    main()
