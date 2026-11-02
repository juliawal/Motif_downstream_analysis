"""Saved pipeline results and reference sources consumed by the paper scripts."""
from pathlib import Path
import yaml
from .common import copy_file

ARMS = ("rank_1", "rank_2", "rank_3", "union")
EVALUATIONS = ("binding_sites.tsv", "binding_sites_per_mf.tsv",
               "ppi_proximity_per_mf.tsv", "gene_set_enrichment.tsv")
PANEL_FILES = ("selected_mfs/mfs.tsv", "evaluation/binding_sites.tsv",
               "evaluation/binding_sites_per_mf.tsv", "evaluation/ppi_proximity.tsv",
               "evaluation/ppi_proximity_per_mf.tsv", "evaluation/gene_set_enrichment.tsv")
STABILITY = ("selection.tsv", "ranking.tsv")
RESOURCE_COPY_KEYS = ("gene_annotation_file", "gene_tss_annotations_file", "cpg_coordinates_file",
                      "cpg_gene_annotations_file", "coding_gene_list_file")


def configured_path(root, config, key):
    path = Path(config[key])
    return path if path.is_absolute() else root / path


def branch_path(root, arm, cohort):
    return root / "results" / ("union" if arm == "union" else f"selected_settings/{arm}") / cohort


def load_config(root):
    return yaml.safe_load((root / "config.yaml").read_text())


def copy_inputs(root, payload, config):
    def copy(path, relative, *, compress=True):
        copy_file(path, payload / (str(relative) + (".gz" if compress else "")),
                  compress=compress)

    copy(root / "config.yaml", "run_configuration/config.yaml", compress=False)
    for key in RESOURCE_COPY_KEYS:
        path = configured_path(root, config, key)
        copy(path, Path("reference_annotations") / key / path.name, compress=path.suffix != ".gz")
    for name in STABILITY:
        relative = Path("results/stability_comparison") / name
        copy(root / relative, Path("original_pipeline_outputs") / relative.relative_to("results"))
    for cohort in config["cohorts"]:
        for name in ("prepared/regulator_names.tsv", "stability/grid_stability.tsv"):
            relative = Path("results") / cohort / name
            copy(root / relative, Path("original_pipeline_outputs") / relative.relative_to("results"))
        for arm in ARMS:
            base = branch_path(root, arm, cohort)
            pair_inputs = ("evaluation_input/mf_cpg_pairs.tsv",) if arm != "union" else ()
            for name in ("selected_mfs/mfs.tsv", *pair_inputs,
                         *("evaluation/" + name for name in EVALUATIONS)):
                copy(base / name, Path("original_pipeline_outputs") / base.relative_to(root / "results") / name)


def copy_panels(panel_root, payload, paper):
    """Keep the two completed group evaluations byte-for-byte in the export."""
    for key in ('manual_evaluation', 'additional_evaluation'):
        relative = Path(paper[key]) / 'results/manual_mfs'
        for name in PANEL_FILES:
            copy_file(panel_root / relative / name,
                      payload / 'panel_evaluations' / relative / name, compress=False)
