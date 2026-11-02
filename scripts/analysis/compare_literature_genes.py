#!/usr/bin/env python3
"""Compare fixed literature gene lists with the saved Motif union evaluations."""
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "analysis_results"
OUT = RESULTS / "literature_comparison"

# Published gene labels and tumour contexts; no PDF extraction is needed.
# Levy-Jurgenson: Figure 3; historical labels are kept as published.
# Detilleux: expression-supported main-text candidates and Figure S2,
# not an exhaustive list from the accompanying web resource.
# Silva/ELMER: Supplementary Table S3. Motif-only labels are excluded.
PAPER_GENES = {
    "LevyJurgenson2019": {
        "BRCA;PRAD": """
            ALG8 AMZ2 B3GAT2 BRCC3 CAMK2A COL13A1 CRB1 ELAVL1 ERI2 FOXL1 FOXO3B GLUD1 GORASP1
            GRPEL1 HAND1 ICA1 IMPA1 KRTAP21-2 LMLN LOC100130987 LRRC42 MRPS18B NODAL NPR1 NTF3
            ODF1 OSBP2 OTOP2 RAB13 RABGAP1 RABIF RASL11B RNASE12 RPS6KB1 SLC25A3 SNW1 SRPK2 SSTR4
            ST7 TRMU USP2 VPS72 XRCC5 ZBTB5 ZKSCAN5 ZNF621
        """,
    },
    "Detilleux2022": {
        "BLCA": """
            BATF EGR1 EGR2 EGR3 GRHL2 KLF10 MYCL NR2F6 PATZ1 PRDM1 SNAI2 SOX2 SOX4 SREBF1 TBX2
            TBX3 TCF3 TP53 USF1 USF2 ZBTB4
        """,
        "BRCA": """
            ARNT2 BATF BHLHE40 CREB3 CREB3L1 CREB3L4 EGR1 EGR2 ESR1 FOXA1 GATA3 HOXA5 IRX5 KLF10
            KLF11 KLF13 KLF15 KLF2 KLF4 KLF5 RARA RFX5 SREBF1 TCF3 TFAP2A TFAP2C THRB USF1 USF2
            ZBTB4
        """,
        "CHOL": """
            ATF6 BATF CEBPB CREB3 CTCF ELF3 FOSL1 FOSL2 FOXA2 JUNB KLF6 KLF9 MAFK NFE2L2 NFIB
            NFKB1 THRA
        """,
        "COAD": """
            ATF6B BATF EGR1 ELF3 EPAS1 FOSL1 FOXQ1 HOXA9 ISX JUN JUNB KLF10 KLF2 KLF3 KLF4 KLF6
            MEF2D SPIB SRF ZBTB7B
        """,
        "HNSC": """
            ARNTL2 ATF6B BATF CEBPB CREB3 EGR1 ETS1 FOSL1 HES2 IRF1 IRF3 IRX4 MSC RUNX3 SNAI2
            SREBF1 USF1 YY1 ZBTB7B
        """,
        "KIRC": """
            ARID5B BACH1 BATF BATF3 BHLHE40 BHLHE41 CEBPB CEBPD EGR1 EMX1 EPAS1 IRF1 IRF7 JUN
            MAFB MXI1 MYC NFIL3 NFIX NFKB2 NR3C1 SPI1 SREBF1 STAT1 USF2
        """,
        "KIRP": """
            ATF4 BATF BHLHE41 CEBPB CEBPD CREB5 IRF5 IRX3 IRX5 KLF6 MLXIPL MYC NFIB NR1H4 SOX4
            SP1 SREBF1 STAT1 TGIF1 TP53 USF2
        """,
        "LIHC": """
            ATF6B CEBPA CEBPG CREB3 EGR1 ESRRA FOXA1 FOXA2 HNF4A HSF1 IRF6 MLXIPL MSC NFE2L1 RARA
            RFX5 SOX9 SREBF2 STAT1 TBX3 USF1 USF2 ZBTB18 ZNF787
        """,
        "LUAD": """
            BATF CREB3L1 CREB3L4 EGR1 EGR2 EPAS1 ERG ETS1 ETV1 FLI1 FOSL1 GATA2 HOXA5 HSF1 KLF10
            SOX18 SPI1 SREBF1 TP53 USF1 ZBTB4
        """,
        "LUSC": """
            ARNTL2 ATF4 EGR1 EGR2 ELF1 ELK3 ERG ETV1 FLI1 FOSL2 GATA2 GATA6 GRHL1 GRHL2 HOXA5
            HSF1 KLF9 MYC MYCL NFE2L2 RARG SNAI2 SPI1 SREBF1 SREBF2 TCF3 TFAP2A TFAP2C TGIF2 TP53
            USF2 ZBTB4
        """,
        "PRAD": """
            CREB3L1 EGR1 EMX2 FOXA1 HOXB8 JUND KLF10 KLF5 MYC NR2F6 SPDEF SREBF1 STAT3 TWIST1
            USF2
        """,
        "UCEC": """
            ASCL2 BHLHE40 CEBPB CREB3L1 CREB3L4 DLX5 DLX6 EGR1 EGR2 ESRRA FOXA1 FOXA2 FOXJ1 FOXM1
            HMG20B HOXB5 KLF10 KLF11 KLF15 LEF1 MSX2 MYBL2 PAX8 PITX1 RXRG SOX12 SOX17 SOX4 SOX9
            SREBF1 TCF3 TP53 USF1 ZNF263 ZNF787
        """,
    },
    "Silva2019_ELMER": {
        "BRCA": """
            ADNP ALX4 AR ARID2 ATOH7 ATOH8 BARX1 BATF BATF3 BCL11A BCL6 BHLHE41 BNC1 CEBPB CEBPG
            DBX2 DLX4 DLX6 DMBX1 E2F3 E2F5 E4F1 EGR1 EGR2 EGR3 ELF5 ELK1 EMX1 EN1 EPAS1 ERG ESR1
            ETS2 ETV3 ETV3L ETV5 ETV6 ETV7 FOSL1 FOXA1 FOXC1 FOXD2 FOXH1 FOXJ1 FOXK2 FOXM1 FOXN1
            FOXO6 FOXP1 FOXP3 GATA2 GATA3 GBX2 GLI2 GLI3 HAND2 HIF3A HIVEP3 HKR1 HNF1A HOMEZ
            HOXB1 HOXB13 HOXB2 HOXB3 HOXB5 HOXB6 HOXC10 HOXC11 HOXC13 IKZF5 IRF7 IRX5 ISL2 ISX
            KDM5B KLF11 KLF14 KLF15 KLF2 KLF4 KLF5 KLF6 KLF9 LBX2 LEF1 LHX1 LHX2 LMX1B MAZ MEIS3
            MNX1 MSX2 MYB MYBL1 MYBL2 MYC MZF1 NEUROG2 NFATC4 NFE2L3 NFIB NFIL3 NKX1-2 NKX2-2
            NKX2-5 NPAS3 NPAS4 NR2E3 NR2F1 NR2F6 NR3C1 NR3C2 OSR1 OTP OTX1 OVOL1 OVOL2 OVOL3
            PATZ1 PAX7 PAX9 PBX1 PGR PITX1 PLAGL1 POU2F1 POU4F1 POU5F1 PPARA PPARG PRDM16 PRDM6
            PROX1 PURA RARA RARB RELB RFX6 RHOXF1 RREB1 RUNX3 RXRA SALL2 SALL4 SIM2 SIX2 SIX4
            SNAI1 SOX11 SOX12 SOX13 SOX6 SOX8 SOX9 SPDEF SPIB STAT5B TAL1 TAL2 TBX15 TCF7L1
            TCF7L2 TEAD4 TFAP2B THRA THRB TLX1 TLX3 TRERF1 TRPS1 TSHZ1 TWIST2 VEZF1 WT1 ZBTB2
            ZBTB4 ZBTB42 ZBTB45 ZBTB47 ZBTB7A ZBTB7B ZBTB7C ZFHX3 ZFP14 ZFP2 ZFP30 ZFP64 ZFP82
            ZFX ZIC1 ZIC4 ZKSCAN1 ZKSCAN2 ZNF124 ZNF138 ZNF157 ZNF175 ZNF20 ZNF213 ZNF232 ZNF24
            ZNF25 ZNF263 ZNF28 ZNF280A ZNF280B ZNF280D ZNF281 ZNF286A ZNF286B ZNF296 ZNF30 ZNF302
            ZNF329 ZNF333 ZNF350 ZNF366 ZNF37A ZNF397 ZNF423 ZNF44 ZNF442 ZNF446 ZNF467 ZNF468
            ZNF500 ZNF502 ZNF512 ZNF516 ZNF521 ZNF525 ZNF546 ZNF552 ZNF559 ZNF563 ZNF586 ZNF587
            ZNF587B ZNF595 ZNF599 ZNF613 ZNF615 ZNF619 ZNF620 ZNF639 ZNF652 ZNF665 ZNF667 ZNF671
            ZNF672 ZNF680 ZNF687 ZNF689 ZNF701 ZNF710 ZNF721 ZNF75A ZNF763 ZNF764 ZNF768 ZNF776
            ZNF785 ZNF799 ZNF814 ZNF816 ZNF823 ZNF844 ZSCAN16 ZSCAN29 ZSCAN32
        """,
    },
}


CRITERIA = {
    'F1': ('f1_available', 'f1_positive'),
    'F2_harmonic': ('f2_available', 'f2_positive'),
    'F2_PageRank': ('f2_available', 'f2_pagerank_excess_positive'),
    'F3': ('f3_available', 'f3_member_significant_min3_term'),
    'F4': ('f4_testable', 'f4_supported'),
    'F1_F4_joint': ('all_four_available', 'all_four_positive'),
}


def joined(values):
    return ';'.join(sorted(set(values)))


def counts(data):
    result = {'n_genes': data.MF.nunique(), 'n_contexts': len(data)}
    for label, (available, positive) in CRITERIA.items():
        for state, flag in [('available', available), ('positive', positive)]:
            result[f'{label}_{state}_contexts'] = int(data[flag].sum())
            result[f'{label}_{state}_genes'] = data.loc[data[flag], 'MF'].nunique()
    stable = data[data.all_four_positive_in_union_and_all_settings]
    result['all_four_all_settings_contexts'] = len(stable)
    result['all_four_all_settings_genes'] = stable.MF.nunique()
    return result


def main():
    union = pd.read_csv(RESULTS / 'combined_evidence/union_evidence_contexts.tsv', sep='\t')
    shortlist = pd.read_csv(RESULTS / 'additional_candidates/candidate_shortlist.tsv', sep='\t')
    rows = [(paper, gene, cohort) for paper, cohorts in PAPER_GENES.items()
            for cohort, genes in cohorts.items() for gene in genes.split()]
    papers = pd.DataFrame(rows, columns=['paper', 'gene', 'paper_cohorts'])
    papers = papers.groupby(['paper', 'gene'], as_index=False).agg(paper_cohorts=('paper_cohorts', joined))
    selected = union.groupby('MF').cohort.agg(joined)
    papers['motif_selected'] = papers.gene.isin(selected.index)
    papers['motif_cohorts'] = papers.gene.map(selected).fillna('')
    papers['same_cohorts'] = [joined(set(p.split(';')) & set(m.split(';')))
                              for p, m in zip(papers.paper_cohorts, papers.motif_cohorts)]

    shared = union[union.MF.isin(papers.gene)].copy().sort_values(['MF', 'cohort'])
    for metric in ['harmonic', 'pagerank']:
        shared[f'{metric}_relative_excess_pct'] = (
            100 * shared[f'{metric}_degree_bin_excess'] / shared[f'{metric}_degree_bin_null_mean']
        ).where(shared.f2_available)
    shared['F1_delta_percentage_points'] = (100 * shared.difference).where(shared.f1_available)
    shared['F4_supported_pairs_pct'] = (100 * shared.n_supported_pairs / shared.n_tested_pairs).where(shared.f4_testable)
    for label, (available, positive) in CRITERIA.items():
        shared[f'{label}_state'] = 'not_available'
        shared.loc[shared[available], f'{label}_state'] = 'not_positive'
        shared.loc[shared[positive], f'{label}_state'] = 'positive'
    shared['papers'] = shared.MF.map(papers.groupby('gene').paper.agg(joined))
    # Group p-values are not gene-level p-values.
    shared = shared.drop(columns=['p_value_harmonic', 'p_value_pagerank', 'f2_harmonic_group_nominal'])
    contexts = papers[papers.motif_selected].merge(shared, left_on='gene', right_on='MF')
    contexts['same_paper_cohort'] = [c in cs.split(';') for c, cs in zip(contexts.cohort, contexts.paper_cohorts)]

    summaries, genes, cells = [], [], []
    for paper in PAPER_GENES:
        data = contexts[contexts.paper.eq(paper)]
        n_paper = int(papers.paper.eq(paper).sum())
        summary = dict(paper=paper, n_literature_genes=n_paper,
                       n_not_selected=n_paper-data.MF.nunique(), **counts(data))
        summaries.append(summary)
        for label in ['F1', 'F2_harmonic', 'F2_PageRank', 'F3', 'F4']:
            total = summary['n_genes']
            available = summary[f'{label}_available_genes']
            positive = summary[f'{label}_positive_genes']
            cells.append(dict(paper=paper, criterion=label, shared_genes=total,
                assessable_genes=available, criterion_met_genes=positive,
                assessable_not_met_genes=available-positive, not_assessable_genes=total-available,
                label=f'{positive}/{available}' if available else 'n/a',
                bar_met_fraction=positive/total, bar_not_met_fraction=(available-positive)/total,
                bar_unavailable_fraction=(total-available)/total))
        for gene, group in data.groupby('MF'):
            row = dict(paper=paper, gene=gene, motif_cohorts=joined(group.cohort), **counts(group))
            row['same_paper_cohorts'] = joined(group.loc[group.same_paper_cohort, 'cohort'])
            for label, (available, positive) in CRITERIA.items():
                row[f'{label}_positive_cohorts'] = joined(group.loc[group[positive], 'cohort'])
                row[f'{label}_not_positive_cohorts'] = joined(group.loc[group[available] & ~group[positive], 'cohort'])
                row[f'{label}_unavailable_cohorts'] = joined(group.loc[~group[available], 'cohort'])
            for metric in ['harmonic', 'pagerank']:
                row[f'{metric}_excess_pct_min'] = group[f'{metric}_relative_excess_pct'].min()
                row[f'{metric}_excess_pct_max'] = group[f'{metric}_relative_excess_pct'].max()
            genes.append(row)

    short = shortlist[['MF', 'n_consistent_direction_pairs',
                       'n_supported_negative', 'n_supported_positive']].copy()
    short['motif_cohorts'] = short.MF.map(selected)
    for paper in PAPER_GENES:
        lookup = papers[papers.paper.eq(paper)].set_index('gene').paper_cohorts
        short[paper] = short.MF.map(lookup).fillna('not in extracted candidate set')
    outputs = {
        'literature_gene_comparison': papers,
        'literature_shortlist_comparison': short,
        'literature_evaluation_summary': pd.DataFrame(summaries),
        'literature_gene_summary': pd.DataFrame(genes),
        'literature_paper_contexts': contexts,
        'literature_comparison': pd.DataFrame(cells),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    for name, frame in outputs.items():
        frame.to_csv(OUT / f'{name}.tsv', sep='\t', index=False, na_rep='NA')
    print(f'Literature comparison: {shared.MF.nunique()} shared genes; results in {OUT}')


if __name__ == '__main__':
    main()
