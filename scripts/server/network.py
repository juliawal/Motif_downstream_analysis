"""Export saved PPI results without loading or recalculating the network."""

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .inventory import ARMS, branch_path, configured_path


def _ids(path):
    frame = pd.read_csv(path, sep="\t", dtype={"MF_ensembl": str})
    return frame, frame.MF_ensembl.drop_duplicates().tolist()


def _summary(saved_per_mf, measure, mapping):
    scores = []
    seen_nodes = set()
    for (gene, saved_measure), row in saved_per_mf.iterrows():
        node = mapping.get(gene)
        if saved_measure == measure and node not in seen_nodes:
            scores.append(row["score"])
            seen_nodes.add(node)
    return {"observed_mean_score": float(np.mean(scores))}


def _per_mf_rows(arm, cohort, selected, pair_ids, pair_symbols, annotation_symbols,
                 memberships, mapping, nodes, saved_per_mf):
    multiplicity = Counter(mapping.get(gene) for gene in pair_ids if mapping.get(gene) in nodes)
    pair_set = set(pair_ids)
    rows = []
    for gene in selected:
        node = mapping.get(gene)
        present = node in nodes
        in_input = gene in pair_set
        if present:
            status = "in_network"
        elif not in_input:
            status = "not_in_evaluation_input"
        elif node is None:
            status = "missing_entrez_mapping"
        else:
            status = "mapped_but_absent_network"
        row = {
            "analysis": arm, "cohort": cohort,
            "MF_ensembl": gene, "MF": pair_symbols.get(gene, annotation_symbols.get(gene, "")),
            "selection_ranks": ";".join(memberships.get(gene, [])),
            "entrez_id": node, "mapping_status": status,
            "in_evaluation_input": in_input,
            "eligible_for_group_test": bool(in_input and present),
            "n_group_entrez_in_network": len(nodes),
            "entrez_ensembl_multiplicity": multiplicity.get(node, 0) if in_input else 0,
        }
        for measure in ("pagerank", "harmonic"):
            observed = np.nan
            mean = np.nan
            if (gene, measure) in saved_per_mf.index:
                observed = saved_per_mf.loc[(gene, measure), "score"]
                mean = saved_per_mf.loc[(gene, measure), "background_mean"]
            row[f"{measure}_score"] = observed
            excess = observed - mean
            row.update({
                f"{measure}_degree_bin_null_mean": mean,
                f"{measure}_degree_bin_excess": excess,
            })
        rows.append(row)
    return rows


def export_ppi(run_root: Path, out_dir: Path, config: dict, cohorts: list[str]) -> None:
    """Export group effects and individual contributions for all four arms.

    The caller should invoke sequentially, use a new output directory, and only
    publish the whole export after all modules succeed.
    """
    run_root = Path(run_root).resolve()
    destination = Path(out_dir) / "ppi_network_proximity"
    annotation = configured_path(run_root, config, "gene_annotation_file").resolve()
    ann = pd.read_csv(annotation, sep="\t", dtype={"ensembl_id": str})
    symbols = dict(zip(ann.ensembl_id, ann.selected_gene_symbol.fillna("")))
    mapped = ann.dropna(subset=["entrez_id"])
    mapping = dict(zip(mapped.ensembl_id, mapped.entrez_id.astype(int)))
    groups, individual = [], []
    for cohort in cohorts:
        memberships = {}
        arm_tables = {}
        for arm in ARMS:
            branch = branch_path(run_root, arm, cohort)
            selection_path = branch / "selected_mfs/mfs.tsv"
            pairs_path = branch / "evaluation_input/mf_cpg_pairs.tsv"
            original_path = branch / "evaluation/ppi_proximity.tsv"
            selected_frame, selected = _ids(selection_path)
            pairs_frame, pair_ids = _ids(pairs_path)
            if arm != "union":
                for gene in selected:
                    memberships.setdefault(gene, []).append(arm.removeprefix("rank_"))
            arm_tables[arm] = (selected, pairs_frame, pair_ids, original_path)
        for arm, (selected, pairs_frame, pair_ids, original_path) in arm_tables.items():
            original = pd.read_csv(original_path, sep="\t", float_precision="round_trip",
                                   dtype={"network_release": str})
            saved_per_mf = pd.read_csv(original_path.with_name("ppi_proximity_per_mf.tsv"),
                                       sep="\t", float_precision="round_trip")
            saved_per_mf = saved_per_mf.set_index(["MF_ensembl", "measure"])
            pair_symbols = dict(zip(pairs_frame.MF_ensembl, pairs_frame.MF.fillna("")))
            saved_genes = set(saved_per_mf.index.get_level_values("MF_ensembl"))
            nodes = []
            for gene in pair_ids:
                node = mapping.get(gene)
                if gene in saved_genes and node not in nodes:
                    nodes.append(node)
            common = {"analysis": arm, "cohort": cohort, "n_selected_mfs": len(selected),
                      "original_source": str(original_path), "n_mapped_ensembl_ids": len(saved_genes),
                      "mf_order": "first appearance in this arm's saved F0 pair table; distinct Entrez IDs preserve that order"}
            original_by_measure = original.set_index("measure")
            for measure in ("pagerank", "harmonic"):
                if not nodes:
                    group = dict.fromkeys(original.columns, np.nan)
                    group.update(measure=measure, n_mfs_nominated=len(pair_ids), n_mfs_in_network=0,
                                 status="no_mfs_in_network", observed_mean_score=np.nan)
                else:
                    group = original_by_measure.loc[measure].to_dict()
                    group.update(measure=measure, status="tested", **_summary(saved_per_mf, measure, mapping))
                groups.append({**group, **common})
            individual.extend(_per_mf_rows(
                arm, cohort, selected, pair_ids, pair_symbols, symbols, memberships,
                mapping, nodes, saved_per_mf))
            print(f"PPI {arm}/{cohort}: {len(pair_ids)} nominated IDs, {len(nodes)} distinct network nodes", flush=True)
    destination.mkdir(parents=True)
    products = {
        "group_effects.tsv": groups, "per_mf.tsv": individual,
    }
    for filename, rows in products.items():
        pd.DataFrame(rows).to_csv(destination / filename, sep="\t", index=False)
    print(f"PPI export: {len(groups)} group rows and {len(individual)} individual rows", flush=True)
