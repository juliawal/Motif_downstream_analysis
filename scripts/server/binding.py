"""Recover case F1 target identities with the completed run's own helpers/resources.

No FIMO, background resampling or inferential testing is performed. Eligibility
uses the same prepared-CpG/coordinate/sequence/GC rules as F1. Matching-category
construction only excludes CpGs with missing GC; its quantile bins do not change
which targets are eligible. Only case target windows are read from the genome.
"""
from __future__ import annotations

import csv
import gzip
from pathlib import Path
import types

from .inventory import ARMS, branch_path, configured_path


CONTEXTS = {"PRAD": "ERG", "THCA": "RUNX1", "UCEC": "HNF1B"}


def _rows(path: Path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def _module(path: Path, name: str):
    # Source execution avoids creating __pycache__ inside the read-only run.
    module = types.ModuleType(name)
    module.__file__ = str(path)
    exec(compile(path.read_bytes(), str(path), "exec"), module.__dict__)
    return module


def _load_helpers(run_root: Path):
    paths = [run_root / "evaluation/binding_site/find_chipseq_overlaps.py",
             run_root / "evaluation/binding_site/prepare_cpgs_for_binding_site_evaluation.py"]
    return tuple(_module(path, f"_paper_export_binding_helper_{i}") for i, path in enumerate(paths))


def resolve_current_chipseq_sources(run_root, config, mf_symbols):
    from .chipseq import resolve_peaks
    resolved, _ = resolve_peaks(run_root, config, sorted(set(mf_symbols)))
    return resolved



def export_binding(run_root: Path, out_dir: Path, config: dict, cohorts: list[str]) -> None:

    run_root, out_dir = Path(run_root).resolve(), Path(out_dir)
    peaks, sequences = _load_helpers(run_root)
    coords_path = configured_path(run_root, config, "cpg_annotations_file")
    genome = configured_path(run_root, config, "genome_fasta")
    fai = configured_path(run_root, config, "genome_fasta_index")
    coordinates = sequences.load_cpg_coordinates(coords_path)
    window = int(config["window_around_cpg"])
    destination = out_dir / "cpg_binding_annotations"
    destination.mkdir(parents=True, exist_ok=True)
    target_path = destination / "target_cpg_hits.tsv.gz"
    fields = ["analysis", "cohort", "MF", "MF_ensembl", "CpG", "track", "f1_status", "eligible", "exclusion_reason", "binding_hit", "chrom", "window_start", "window_end", "peak_file", "jaspar_ids"]
    n_targets = 0
    resolved_peak_cache = {}
    with gzip.open(target_path, "wt", newline="") as target_handle:
        writer = csv.DictWriter(target_handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for cohort in cohorts:
            if cohort not in CONTEXTS:
                continue
            mf = CONTEXTS[cohort]
            names_path = run_root / "results" / cohort / "prepared/target_names.tsv"
            prepared = set(names_path.read_text().splitlines())
            for arm in ARMS:
                base = branch_path(run_root, arm, cohort)
                pair_path = base / "evaluation_input/mf_cpg_pairs.tsv"
                summary_path = base / "evaluation/binding_sites_per_mf.tsv"
                selected_pairs = [row for row in _rows(pair_path) if row["MF"] == mf]
                ids = {}
                for row in selected_pairs:
                    ids[row["CpG"]] = row["MF_ensembl"]
                original = [row for row in _rows(summary_path) if row["MF"] == mf]
                if not ids:
                    continue
                original = original[0]
                available = set(ids) & prepared & set(coordinates.index)
                target_coordinates = coordinates.loc[sorted(available)]
                seqs = sequences.extract_cpg_sequences(target_coordinates, genome, fai, window_around_cpg=window)
                gc = sequences.compute_gc_fraction_by_cpg(seqs)
                eligible = set(gc[gc.notna()].index)
                track, status = original["track"], original["status"]
                hit_cpgs, selected_peak = set(), None
                if track == "FIMO":
                    hit_path = base / "evaluation/binding_evidence/fimo_hits_with_cpg_ids.tsv"
                    pwm_ids = set(original["jaspar_ids"].split(",")) - {""}
                    hit_cpgs = {row["cpg_id"] for row in _rows(hit_path) if row["motif_id"] in pwm_ids}
                elif track == "ChIP_seq":
                    if mf not in resolved_peak_cache:
                        mapping = resolve_current_chipseq_sources(run_root, config, [mf])
                        resolved_peak_cache.update(mapping)
                    selected_peak = resolved_peak_cache[mf]
                    if selected_peak is not None:
                        hit_cpgs = peaks.cpgs_with_peak_overlap(coordinates.loc[sorted(eligible)], selected_peak, window_around_cpg=window) if eligible else set()
                for cpg in sorted(ids):
                    if cpg not in prepared:
                        reason = "not_in_prepared_targets"
                    elif cpg not in coordinates.index:
                        reason = "missing_coordinates"
                    elif cpg not in seqs:
                        reason = "sequence_more_than_10_percent_N"
                    elif cpg not in eligible:
                        reason = "missing_gc_fraction"
                    else:
                        reason = ""
                    coordinate = coordinates.loc[cpg] if cpg in coordinates.index else None
                    row = dict(analysis=arm, cohort=cohort, MF=mf, MF_ensembl=ids[cpg], CpG=cpg, track=track, f1_status=status,
                               eligible=not bool(reason), exclusion_reason=reason,
                               binding_hit=(cpg in hit_cpgs) if not reason and status == "tested" else "",
                               chrom=coordinate["CpG_chrm"] if coordinate is not None else "",
                               window_start=max(0, int(coordinate["CpG_beg"])-window) if coordinate is not None else "",
                               window_end=int(coordinate["CpG_end"])+window if coordinate is not None else "",
                               peak_file=str(selected_peak.absolute()) if selected_peak else "", jaspar_ids=original.get("jaspar_ids", ""))
                    writer.writerow(row)
                    n_targets += 1
    print(f"Binding export: {n_targets} target CpG rows", flush=True)
