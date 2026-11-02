#!/usr/bin/env python3
"""Summarize selection, stability and F1–F3 from saved pipeline results.
Original p/q-values are retained. Individual effects and contributions are descriptive.
The output root contains separate selection/stability and F1–F3 result folders;
the cross-stage Union coverage is saved alongside the joint F1–F4 evaluation.
"""
from __future__ import annotations

import argparse
import gzip
import itertools
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import pandas as pd
import yaml

ARMS = ["rank_1", "rank_2", "rank_3", "union"]
KEY = ["analysis", "cohort"]
PPI_ORIGINAL_COLUMNS = [
    "measure", "n_random_sets", "p_value", "n_mfs_nominated", "n_mfs_in_network",
    "network_release", "network_file", "mapping_file", "background_gene_ids_file",
    "n_graph_nodes", "n_graph_edges", "degree_bin_min_size", "random_seed",
    "pagerank_alpha", "dnmt_tet_gene_symbols", "dnmt_tet_entrez_ids",
]
SELECTION_TABLES = {
    "stability_selected_settings.tsv", "stability_setting_ranking.tsv",
    "stability_by_cohort.tsv", "stability_selected_by_cohort.tsv",
    "selected_gene_cohort_memberships.tsv", "selection_overlap.tsv",
    "selection_recurrence.tsv",
}


def joined(values):
    return ";".join(sorted({str(value) for value in values if pd.notna(value) and str(value)}))


def tokens(value, separator=";"):
    return set(str(value).split(separator)) - {"", "nan"} if pd.notna(value) else set()


def theme(name):
    """The paper's original descriptive vocabulary; classes have no own test."""
    s = name.lower()
    if any(v in s for v in ["methylation", "demethylation", "methyltransferase"]):
        return "DNA_methylation_named"
    if any(v in s for v in ["immune", "lymphocyte", "t cell", "t-cell", "b cell", "leukocyte", "antigen", "thym", "hematop", "haemopo", "cytokine", "interleukin"]):
        return "immune_cell_state"
    if any(v in s for v in ["dna binding", "dna-binding", "transcription", "cis-regulatory"]):
        return "DNA_binding_transcription"
    if any(v in s for v in ["development", "differentiation", "morphogenesis", "pattern", "neuron", "neural", "neurogen", "embry", "axon", "limb", "brain", "epithelial", "kidney", "mesenchym", "skeletal", "organ"]):
        return "development_cell_identity"
    return "other"


class Analysis:
    def __init__(self, source, output):
        self.source = Path(source).resolve()
        self.output = Path(output).resolve()
        self.selection_output = self.output / "candidate_selection"
        self.evaluation_output = self.output / "binding_network_enrichment"
        self.joint_output = self.output / "combined_evidence"
        for directory in (self.selection_output, self.evaluation_output, self.joint_output):
            directory.mkdir(parents=True, exist_ok=True)
        self.config = yaml.safe_load(self.text("run_configuration/config.yaml"))
        self.cohorts = sorted(self.config["cohorts"])
        self.annotation = self.read(next((self.source / "reference_annotations/gene_annotation_file").glob("*.tsv.gz")))
        self.symbols = self.annotation.set_index("ensembl_id").selected_gene_symbol.to_dict()
        self.regulators = {c: set(self.text(f"original_pipeline_outputs/{c}/prepared/regulator_names.tsv.gz").splitlines()) for c in self.cohorts}
        # The export retains every original F2 field in this combined table.
        self.ppi_effects = self.read("ppi_network_proximity/group_effects.tsv")

    def text(self, path):
        path = self.source / path
        return gzip.decompress(path.read_bytes()).decode() if path.suffix == ".gz" else path.read_text()

    def read(self, path, **kwargs):
        return pd.read_csv(self.source / path, sep="\t", **kwargs)

    def save(self, frame, name):
        if name in SELECTION_TABLES:
            directory = self.selection_output
        elif name == "union_coverage.tsv":
            directory = self.joint_output
        else:
            directory = self.evaluation_output
        frame.to_csv(directory / name, sep="\t", index=False,
                     compression={"method": "gzip", "mtime": 0} if name.endswith(".gz") else None)

    def base(self, arm, cohort):
        return Path("original_pipeline_outputs") / ("union" if arm == "union" else f"selected_settings/{arm}") / cohort

    def stability(self):
        prefix = "original_pipeline_outputs/stability_comparison/"
        selected = self.read(prefix + "selection.tsv.gz")
        ranking = self.read(prefix + "ranking.tsv.gz")
        grid_rows = []
        for cohort in self.cohorts:
            grid = self.read(f"original_pipeline_outputs/{cohort}/stability/grid_stability.tsv.gz")
            grid["setting_id"] = grid.k_edges.astype(str) + "__" + grid.configuration_name
            grid_rows.append(grid)
        grid = pd.concat(grid_rows, ignore_index=True)
        eligible = grid[grid.enough_mfs.eq(True)]
        selected_cohort = grid[["cohort", "k_edges", "configuration_name", "n_mfs",
                               "mean_kuncheva_stability", "setting_id"]].merge(
            selected[["setting_id", "selection_rank"]], on="setting_id")
        selected_cohort = selected_cohort.sort_values(["selection_rank", "cohort"], kind="stable")
        selected = selected.copy()
        selected["analysis"] = "rank_" + selected.selection_rank.astype(str)
        selected["paper_label"] = "S" + selected.selection_rank.astype(str)
        selected["descriptive_label"] = [f"S{r.selection_rank}: {r.configuration_name}; k={r.k_edges:,}" for r in selected.itertuples()]
        self.save(selected, "stability_selected_settings.tsv")
        self.save(ranking, "stability_setting_ranking.tsv")
        self.save(grid, "stability_by_cohort.tsv")
        self.save(selected_cohort, "stability_selected_by_cohort.tsv")
        return dict(n_settings=len(ranking), n_cohorts=eligible.cohort.nunique(), n_eligible_rows=len(eligible), selected_settings=selected.fillna("").to_dict("records"), status="PASS")

    def run(self):
        self.stability()
        selections, scopes, overlaps, f1_frames, f1_groups, f3_sig, f3_counts, ppi_original = [], [], [], [], [], [], [], []
        for cohort in self.cohorts:
            pairsets, mfsets = {}, {}
            bg_symbols = {self.symbols[g] for g in self.regulators[cohort] if g in self.symbols and pd.notna(self.symbols[g])}
            for arm in ARMS:
                base = self.base(arm, cohort)
                ids = self.read(base / "selected_mfs/mfs.tsv.gz").MF_ensembl
                pair_path = (Path("f4_methylation_factor_cpg_pairs") / f"union_{cohort}.tsv.gz"
                             if arm == "union" else base / "evaluation_input/mf_cpg_pairs.tsv.gz")
                pairs = self.read(pair_path, usecols=["MF", "MF_ensembl", "CpG"])
                mfsets[arm] = set(ids)
                pairsets[arm] = set(zip(pairs.MF_ensembl, pairs.CpG))
                counts = pairs.groupby("MF_ensembl").size()
                for gene in ids:
                    selections.append(dict(analysis=arm, cohort=cohort, MF_ensembl=gene, MF=self.symbols[gene], n_input_pairs=int(counts.get(gene, 0)), n_unique_cpg=int(counts.get(gene, 0)), in_saved_regulator_universe=gene in self.regulators[cohort], source=str(base / "selected_mfs/mfs.tsv.gz")))
                scopes.append(dict(analysis=arm, cohort=cohort, n_selected=len(ids), n_with_pairs=pairs.MF_ensembl.nunique(), n_unique_pairs=len(pairs), n_rows_pairs=len(pairs), n_unique_cpg=pairs.CpG.nunique(), n_predictor_universe=len(self.regulators[cohort]), selected_source=str(base / "selected_mfs/mfs.tsv.gz"), pairs_source=str(pair_path)))
                symbol_ids = pairs[["MF", "MF_ensembl"]].drop_duplicates()
                f1path = base / "evaluation/binding_sites_per_mf.tsv.gz"
                binding = self.read(f1path).merge(symbol_ids, on="MF")
                binding["analysis"], binding["cohort"], binding["source"] = arm, cohort, str(f1path)
                binding["difference"] = binding.n_target_hits / binding.n_target_cpgs - binding.n_background_hits / binding.n_background_cpgs
                tested = binding[binding.status.eq("tested")]
                f1_frames.append(binding)
                group_path = base / "evaluation/binding_sites.tsv.gz"
                group = self.read(group_path)
                group["analysis"], group["cohort"], group["source"], group["n_selected_mfs"] = arm, cohort, str(group_path), len(ids)
                for index, g in group.iterrows():
                    bt = tested[tested.track.eq(g.track)]
                    # Reconstruct differences from saved integer numerators
                    # to avoid decimal-serialization noise at exact ties.
                    delta = bt.n_target_hits / bt.n_target_cpgs - bt.n_background_hits / bt.n_background_cpgs
                    pos, neg, zero = int(delta.gt(0).sum()), int(delta.lt(0).sum()), int(delta.eq(0).sum())
                    for name, value in dict(n_positive_reconstructed=pos, n_negative_reconstructed=neg, n_zero_reconstructed=zero,
                                            median_target_rate=bt.target_hit_rate.median(), median_background_rate=bt.bg_hit_rate.median(),
                                            n_input_pairs_measured=bt.n_target_input_cpgs.sum(), n_retained_target_cpg_memberships=bt.n_target_cpgs.sum(),
                                            n_target_hits_memberships=bt.n_target_hits.sum(), n_background_hits_memberships=bt.n_background_hits.sum(),
                                            n_background_cpg_memberships=bt.n_background_cpgs.sum()).items():
                        group.loc[index, name] = value
                f1_groups.append(group)
                ppi = self.ppi_effects.loc[self.ppi_effects.analysis.eq(arm) &
                                           self.ppi_effects.cohort.eq(cohort), PPI_ORIGINAL_COLUMNS].copy()
                ppi["analysis"], ppi["cohort"], ppi["source"] = arm, cohort, f"ppi_network_proximity/group_effects.tsv#analysis={arm}&cohort={cohort}"
                ppi_original.append(ppi)

                term_path = base / "evaluation/gene_set_enrichment.tsv.gz"
                terms = self.read(term_path)
                terms["analysis"], terms["cohort"], terms["source"] = arm, cohort, str(term_path)
                selected_symbols = set(pairs.MF)
                terms["descriptive_theme"] = terms.term_name.map(theme)
                sig = terms[terms.adj_p_value.le(.05)].copy()
                f3_sig.append(sig)
                f3_counts.append(dict(analysis=arm, cohort=cohort, n_tested_terms=len(terms), n_background_genes=len(bg_symbols), n_mf_symbols=len(selected_symbols), n_nominal_terms=int(terms.p_value.lt(.05).sum()), n_bh_terms=len(sig), n_bh_min3_terms=int(sig.term_has_min_mf_overlap.sum()), min_p=terms.p_value.min(), min_q=terms.adj_p_value.min(), source=str(term_path)))
            overlap = dict(cohort=cohort, rank1=len(mfsets["rank_1"]), rank2=len(mfsets["rank_2"]), rank3=len(mfsets["rank_3"]), intersection_all3=len(set.intersection(*(mfsets[a] for a in ARMS[:3]))), union=len(mfsets["union"]), union_exact_original_union=True, union_added_to_rank1=len(mfsets["union"] - mfsets["rank_1"]), pairs_union=len(pairsets["union"]), pairs_union_exact_original_union=True)
            for a, b in itertools.combinations(ARMS[:3], 2):
                overlap[f"{a}_{b}_jaccard"] = len(mfsets[a] & mfsets[b]) / len(mfsets[a] | mfsets[b])
            overlaps.append(overlap)
            print(f"EVIDENCE_COHORT_PASS {cohort}", flush=True)
        members, scopes, overlap = pd.DataFrame(selections), pd.DataFrame(scopes), pd.DataFrame(overlaps)
        binding, groups, terms, term_counts = pd.concat(f1_frames, ignore_index=True), pd.concat(f1_groups, ignore_index=True), pd.concat(f3_sig, ignore_index=True), pd.DataFrame(f3_counts)
        original_ppi = pd.concat(ppi_original, ignore_index=True)
        for frame, name in [(members, "selected_gene_cohort_memberships.tsv"), (overlap, "selection_overlap.tsv"), (binding, "binding_measurements.tsv"), (groups, "binding_group_tests.tsv"), (terms, "enriched_terms.tsv"), (term_counts, "enrichment_test_counts.tsv"), (original_ppi, "network_group_tests.tsv")]:
            self.save(frame, name)
        self.summaries(members, scopes, overlap, binding, groups, terms, term_counts)

    def summaries(self, members, scopes, overlap, binding, groups, terms, term_counts):
        recurrence = []
        for gene, rows in members.groupby("MF_ensembl"):
            sets = {a: set(rows[rows.analysis.eq(a)].cohort) for a in ARMS}
            all3 = set.intersection(*(sets[a] for a in ARMS[:3]))
            eligible = {c for c, universe in self.regulators.items() if gene in universe}
            rec = dict(MF_ensembl=gene, MF=self.symbols[gene], n_eligible_cohorts=len(eligible), eligible_cohorts=joined(eligible), n_union_cohorts=len(sets["union"]), union_fraction_eligible=len(sets["union"]) / len(eligible), union_cohorts=joined(sets["union"]), n_all3_same_cohorts=len(all3), all3_cohorts=joined(all3), all3_fraction_union=len(all3) / len(sets["union"]), all3_fraction_eligible=len(all3) / len(eligible), n_union_input_pairs=int(rows[rows.analysis.eq("union")].n_input_pairs.sum()))
            for arm in ARMS[:3]:
                rec[f"n_{arm}_cohorts"], rec[f"{arm}_cohorts"] = len(sets[arm]), joined(sets[arm])
            recurrence.append(rec)
        recurrence = pd.DataFrame(recurrence).sort_values(["n_union_cohorts", "n_all3_same_cohorts", "MF"], ascending=[False, False, True])
        self.save(recurrence, "selection_recurrence.tsv")
        term_members = []
        lookup = members.groupby(KEY + ["MF"])["MF_ensembl"].agg(list).to_dict()
        for row in terms.itertuples():
            for symbol in sorted(tokens(row.mf_symbols_in_term)):
                for gene in lookup[row.analysis, row.cohort, symbol]:
                    term_members.append(dict(analysis=row.analysis, cohort=row.cohort, MF_ensembl=gene, MF=symbol, library=row.library_name, term=row.term_name, q_term=row.adj_p_value, term_has_min_mf_overlap=row.term_has_min_mf_overlap, theme=row.descriptive_theme, source=row.source))
        term_members = pd.DataFrame(term_members)
        term_settings = []
        for (cohort, library, term), tt in terms.groupby(["cohort", "library_name", "term_name"]):
            arms = set(tt.analysis)
            min3_arms = set(tt[tt.term_has_min_mf_overlap].analysis)
            term_settings.append(dict(cohort=cohort, library_name=library, term_name=term, significant_arms=joined(arms), significant_min3_arms=joined(min3_arms), significant_all_four=set(ARMS) <= arms, significant_min3_all_four=set(ARMS) <= min3_arms, significant_in_union="union" in arms, significant_min3_in_union="union" in min3_arms))
        self.save(pd.DataFrame(term_settings), "enrichment_setting_agreement.tsv")
        f3_gene_context = term_members[term_members.term_has_min_mf_overlap].groupby(KEY + ["MF_ensembl", "MF"]).agg(n_significant_min3_terms=("term", "size"), f3_themes=("theme", joined)).reset_index()
        self.save(f3_gene_context, "enrichment_gene_contexts.tsv")
        effects = self.ppi_effects
        per_mf = self.read("ppi_network_proximity/per_mf.tsv")
        self.save(effects[effects.analysis.eq("union")], "union_network_effects.tsv")
        self.save(per_mf, "network_individual_scores.tsv")
        self.ppi_contributions(per_mf)
        concordance = []
        for measure, ee in effects.groupby("measure"):
            pivot = ee.pivot(index="cohort", columns="analysis", values="p_value").reindex(columns=ARMS)
            for cohort, row in pivot.iterrows():
                concordance.append(dict(cohort=cohort, measure=measure, all_four_nominal=bool(row.le(.05).all()), union_nominal=bool(row["union"] <= .05), **{f"p_{arm}": row[arm] for arm in ARMS}))
        self.save(pd.DataFrame(concordance), "network_setting_agreement.tsv")
        coverage = scopes.copy()
        for track in ["FIMO", "ChIP_seq"]:
            gg = groups[groups.track.eq(track)].set_index(KEY)
            for name in ["n_mf", "n_selected_mfs", "n_positive_reconstructed", "n_negative_reconstructed", "n_zero_reconstructed", "median_target_rate", "median_background_rate", "median_difference", "p_value_sign_test"]:
                coverage = coverage.merge(gg[[name]].rename(columns={name: f"F1_{track}_{name}"}), on=KEY)
        for measure in ["harmonic", "pagerank"]:
            ee = effects[effects.measure.eq(measure)].set_index(KEY)
            names = ["n_mfs_nominated", "n_mfs_in_network", "p_value", "observed_mean_score"]
            coverage = coverage.merge(ee[names].rename(columns={n: f"F2_{measure}_{n}" for n in names}), on=KEY)
        coverage = coverage.merge(term_counts.rename(columns={"n_tested_terms": "F3_n_tested_terms", "n_bh_terms": "F3_n_adjusted_p_le_05", "n_bh_min3_terms": "F3_n_adjusted_p_le_05_min3_mfs"}), on=KEY)
        self.save(coverage, "evaluation_coverage_and_tests.tsv")
        st1 = []
        for c in self.cohorts:
            ff = binding[binding.analysis.eq("union") & binding.cohort.eq(c)]
            mm = per_mf[per_mf.analysis.eq("union") & per_mf.cohort.eq(c)]
            mediation = pd.Series(json.loads(self.text(f"f4_methylation_factor_cpg_pairs/union_{c}.json")))
            st1.append(dict(cohort=c, n_union_selected=len(members[members.analysis.eq("union") & members.cohort.eq(c)]), n_in_all_three_settings=int(overlap.set_index("cohort").loc[c, "intersection_all3"]), n_formal_F1=int(ff.status.eq("tested").sum()), n_usable_FIMO=int((ff.status.eq("tested") & ff.track.eq("FIMO")).sum()), n_usable_ChIP=int((ff.status.eq("tested") & ff.track.eq("ChIP_seq")).sum()), n_PPI_mapped=int(mm.mapping_status.eq("in_network").sum()), n_F4_samples=int(mediation.n_samples), F4_testable=bool(mediation.testable)))
        self.save(pd.DataFrame(st1), "union_coverage.tsv")
        print(f"Selection and F1–F3 summaries: {len(self.cohorts)} cohorts", flush=True)


    def ppi_contributions(self, per_mf):
        """Summarize saved pipeline contributions per distinct union network node.

        These shares do not decompose the group null or assign gene significance.
        Ensembl aliases for one Entrez node receive one combined row.
        """
        contributions, concentration = [], []
        union = per_mf[per_mf.analysis.eq("union") & per_mf.eligible_for_group_test]
        for cohort in self.cohorts:
            saved_contributions = self.read(
                self.base("union", cohort) / "evaluation/ppi_proximity_per_mf.tsv.gz",
                float_precision="round_trip")
            for measure in ("harmonic", "pagerank"):
                rows = union[union.cohort.eq(cohort)]
                shares = saved_contributions[saved_contributions.measure.eq(measure)]
                rows = rows.merge(shares[["MF_ensembl", "contribution"]], on="MF_ensembl", how="left")
                nodes = rows.groupby("entrez_id", sort=False).agg(
                    MF=("MF", joined), MF_ensembl=("MF_ensembl", joined),
                    score=(f"{measure}_score", "first"),
                    degree_pool_mean=(f"{measure}_degree_bin_null_mean", "first"),
                    excess=(f"{measure}_degree_bin_excess", "first"),
                    contribution=("contribution", "first")).reset_index()
                nodes["positive_excess"] = nodes.excess.clip(lower=0)
                total = nodes.positive_excess.sum()
                nodes["positive_mass_fraction"] = nodes.pop("contribution")
                nodes["cohort"], nodes["measure"] = cohort, measure
                nodes = nodes.sort_values(["positive_mass_fraction", "MF"], ascending=[False, True])
                contributions.append(nodes)
                positive = nodes[nodes.positive_mass_fraction.gt(0)]
                concentration.append(dict(cohort=cohort, measure=measure, n_in_network=len(nodes),
                    n_positive_excess=len(positive), top1_gene=joined(positive.head(1).MF),
                    top3_genes=";".join(positive.head(3).MF),
                    top1_positive_mass_fraction=positive.head(1).positive_mass_fraction.sum(min_count=1),
                    top3_positive_mass_fraction=positive.head(3).positive_mass_fraction.sum(min_count=1),
                    positive_excess_mass=total))
        self.save(pd.concat(contributions, ignore_index=True), "network_positive_contributions.tsv")
        self.save(pd.DataFrame(concentration), "network_contribution_concentration.tsv")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path(__file__).resolve().parents[2] / "pipeline_results")
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[2] / "analysis_results",
                        help="Analysis result root, e.g. analysis_results; creates separate result subfolders")
    args = parser.parse_args()
    Analysis(args.input, args.output).run()
