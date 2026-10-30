"""Export PPI effects using the completed run's own F2 implementation.

Original group p-values are read from every saved setting and union branch. New null draws reproduce the original RNG sequence; they are
not a new significance correction. Individual values remain descriptive.
"""

from collections import Counter
from contextlib import contextmanager
import importlib
from pathlib import Path
import random
import sys

import numpy as np
import pandas as pd

from .inventory import ARMS, branch_path, configured_path


@contextmanager
def _run_implementation(run_root):
    """Import only the requested run, restoring callers' module cache afterward.

    Python caches imports by package name, so merely prepending sys.path would
    silently reuse a previously loaded implementation from another run. This context is
    deliberately sequential and is not safe to call concurrently in threads.
    """
    prefixes = ("pipeline", "evaluation", "data_processing")

    def related(name):
        return any(name == p or name.startswith(p + ".") for p in prefixes)

    saved = {name: module for name, module in sys.modules.items() if related(name)}
    old_path = sys.path[:]
    for name in saved:
        del sys.modules[name]
    sys.path.insert(0, str(run_root))
    importlib.invalidate_caches()
    try:
        runner = importlib.import_module("pipeline.evaluate_ppi_network_proximity")
        yield runner
    finally:
        for name in list(sys.modules):
            if related(name):
                del sys.modules[name]
        sys.modules.update(saved)
        sys.path[:] = old_path


def _ids(path):
    frame = pd.read_csv(path, sep="\t", dtype={"MF_ensembl": str})
    return frame, frame.MF_ensembl.drop_duplicates().tolist()


def _summary(scores, nodes, samples):
    observed = float(np.mean([scores.get(node, 0.0) for node in nodes]))
    null = np.asarray([np.mean([scores.get(node, 0.0) for node in sample]) for sample in samples])
    mean = float(null.mean())
    return {
        "observed_mean_score": observed,
        "random_set_null_mean": mean,
        "relative_excess_percent": 100 * (observed - mean) / mean if mean else float("nan"),
    }


def _per_mf_rows(arm, cohort, selected, pair_ids, pair_symbols, annotation_symbols,
                 memberships, graph, mapping, nodes, matched, scores, n_background):
    multiplicity = Counter(mapping.get(gene) for gene in pair_ids if mapping.get(gene) in graph)
    pair_set = set(pair_ids)
    rows = []
    for gene in selected:
        node = mapping.get(gene)
        present = node in graph
        in_input = gene in pair_set
        pool = matched.get(node, []) if in_input else []
        status = "in_network" if present else "missing_entrez_mapping" if node is None else "mapped_but_absent_network"
        row = {
            "analysis": arm, "cohort": cohort,
            "MF_ensembl": gene, "MF": pair_symbols.get(gene, annotation_symbols.get(gene, "")),
            "selection_ranks": ";".join(memberships.get(gene, [])),
            "entrez_id": node, "mapping_status": status,
            "in_evaluation_input": in_input,
            "eligible_for_group_test": bool(in_input and present),
            "comparison_status": "tested" if pool else "not_in_evaluation_input" if not in_input else status,
            "graph_degree": graph.degree(node) if present else np.nan,
            "n_group_entrez_in_network": len(nodes),
            "n_background_entrez_in_network": n_background,
            "entrez_ensembl_multiplicity": multiplicity.get(node, 0) if in_input else 0,
            "n_degree_bin_comparators": len(pool),
        }
        for measure, score in scores.items():
            observed = score.get(node, 0.0) if present else np.nan
            row[f"{measure}_score"] = observed
            values = np.asarray([score.get(g, 0.0) for g in pool])
            mean = float(values.mean()) if pool else np.nan
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
    network = configured_path(run_root, config, "hippie_network_file").resolve()
    annotation = configured_path(run_root, config, "gene_annotation_file").resolve()
    n_random = int(config["n_random_sets"])
    bin_size = int(config["degree_bin_min_size"])
    seed = int(config["proximity_random_seed"])
    alpha = float(config["proximity_pagerank_alpha"])
    groups, individual = [], []
    with _run_implementation(run_root) as runner:
        graph = runner.load_hippie_network(network)
        mapping = runner.load_ensembl_to_entrez_mapping(annotation)
        enzymes = dict(runner.dnmt_tet_genes)
        ann = pd.read_csv(annotation, sep="\t", dtype={"ensembl_id": str})
        symbols = dict(zip(ann.ensembl_id, ann.selected_gene_symbol.fillna("")))
        # These graph-wide quantities are identical in all cohorts and arms.
        scores = None
        for cohort in cohorts:
            background_path = run_root / "results" / cohort / "prepared/regulator_names.tsv"
            gene_ids = background_path.read_text().splitlines()
            background = runner.get_background_genes(graph, gene_ids, mapping)
            degrees = runner.get_genes_degrees_in_ppi(graph, background)
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
            current_cache = {}
            for arm, (selected, pairs_frame, pair_ids, original_path) in arm_tables.items():
                original = pd.read_csv(original_path, sep="\t", float_precision="round_trip",
                                       dtype={"network_release": str})
                pair_symbols = dict(zip(pairs_frame.MF_ensembl, pairs_frame.MF.fillna("")))
                nodes = runner.get_mf_entrez_ids_in_ppi(graph, pair_ids, mapping)
                parameters = {
                    "n_random_sets": n_random, "n_mfs_nominated": len(pair_ids),
                    "n_mfs_in_network": len(nodes), "network_release": str(config["hippie_release"]),
                    # Pipeline provenance records the configured argument's
                    # basename, which can differ from a symlink target's name.
                    "network_file": Path(config["hippie_network_file"]).name,
                    "mapping_file": Path(config["gene_annotation_file"]).name,
                    "background_gene_ids_file": background_path.name,
                    "n_graph_nodes": graph.number_of_nodes(), "n_graph_edges": graph.number_of_edges(),
                    "degree_bin_min_size": bin_size, "random_seed": seed, "pagerank_alpha": alpha,
                    "dnmt_tet_gene_symbols": ",".join(enzymes),
                    "dnmt_tet_entrez_ids": ",".join(str(e) for e in enzymes.values()),
                }
                common = {"analysis": arm, "cohort": cohort, **parameters,
                          "n_selected_mfs": len(selected), "n_background_genes_in_network": len(background),
                          "original_source": str(original_path), "n_mapped_ensembl_ids": sum(mapping.get(g) in graph for g in pair_ids),
                          "mf_order": "first appearance in this arm's saved F0 pair table; distinct Entrez IDs preserve that order"}
                if scores is None and (nodes or any(mapping.get(gene) in graph for gene in selected)):
                    scores = {
                        "pagerank": runner.personalized_pagerank(graph, enzymes, alpha),
                        "harmonic": runner.harmonic_centrality_from_enzymes(graph, enzymes),
                    }
                if not nodes:
                    matched = {}
                    for measure in ("pagerank", "harmonic"):
                        empty_stats = {key: np.nan for key in (
                            "observed_mean_score", "random_set_null_mean", "relative_excess_percent", "p_value")}
                        empty = {**common, "measure": measure,
                                 "status": "no_mfs_in_network", **empty_stats}
                        groups.append(empty)
                    # A completely unmapped cohort needs neither PageRank nor BFS.
                    local_scores = scores or {"pagerank": {}, "harmonic": {}}
                else:
                    local_scores = scores
                    key = tuple(nodes)
                    if key not in current_cache:
                        bins = runner.build_degree_bins(degrees, bin_size)
                        matched = runner.get_degree_matched_genes(nodes, degrees, bins, enzymes)
                        rng = random.Random(seed)
                        samples = [runner.sample_degree_matched_genes(matched, nodes, rng) for _ in range(n_random)]
                        stats = {measure: _summary(score, nodes, samples) for measure, score in scores.items()}
                        current_cache[key] = (matched, stats)
                    matched, stats = current_cache[key]
                    original_by_measure = original.set_index("measure")
                    for measure, values in stats.items():
                        stored_p = float(original_by_measure.loc[measure, "p_value"])
                        groups.append({**common, "measure": measure, "status": "tested",
                                       **values,
                                       "p_value": stored_p})
                individual.extend(_per_mf_rows(
                    arm, cohort, selected, pair_ids, pair_symbols, symbols, memberships, graph, mapping,
                    nodes, matched, local_scores, len(background)))
                print(f"PPI {arm}/{cohort}: {len(pair_ids)} nominated IDs, {len(nodes)} distinct network nodes", flush=True)
    destination.mkdir(parents=True)
    products = {
        "group_effects.tsv": groups, "per_mf.tsv": individual,
    }
    for filename, rows in products.items():
        pd.DataFrame(rows).to_csv(destination / filename, sep="\t", index=False)
    print(f"PPI export: {len(groups)} group rows and {len(individual)} individual rows", flush=True)
