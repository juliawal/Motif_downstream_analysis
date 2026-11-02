"""Reduce saved F4 triplets in chunks; retain all F0 pairs and ERG–TDRD1 rows.

Support and q-values come from the pipeline. No models or tests are fitted here.
"""
from collections import Counter
from contextlib import ExitStack
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
from .common import json_write

KEY = ["MF_ensembl", "CpG"]
PAIR_COLUMNS = [
    "analysis", "cohort", "MF", "MF_ensembl", "CpG",
    "n_candidate_triplets", "n_tested_triplets", "n_supported_triplets",
    "n_nominal_joint_triplets", "a_min", "a_max", "a_pvalue_min", "joint_q_min",
    "n_supported_b_positive", "n_supported_b_negative", "supported_b_sum",
    "supported_b_abs_sum", "supported_abs_distance_sum", "pair_status",
]
TESTED = {"significant", "not_significant"}
COUNT_COLUMNS = [
    "n_candidate_triplets", "n_tested_triplets", "n_supported_triplets",
    "n_nominal_joint_triplets", "n_supported_b_positive", "n_supported_b_negative",
]
SUM_COLUMNS = ["supported_b_sum", "supported_b_abs_sum", "supported_abs_distance_sum"]
MIN_COLUMNS = ["a_min", "a_pvalue_min", "joint_q_min"]
ERG_ID = "ENSG00000157554"
TDRD1_ID = "ENSG00000095627"




def export_f4_branch(run_root, branch_dir, out_dir, arm, cohort, config, *, chunksize=100_000):
    run_root, branch_dir, out_dir = map(Path, (run_root, branch_dir, out_dir))
    input_path = branch_dir / "evaluation_input/mf_cpg_pairs.tsv"
    summary_path = branch_dir / "evaluation/mediation.tsv"
    triplet_path = branch_dir / "evaluation/mediation_triplets.tsv"
    inputs = pd.read_csv(input_path, sep="\t", usecols=["MF", *KEY], dtype=str, keep_default_na=False)
    original = pd.read_csv(summary_path, sep="\t", dtype=str, keep_default_na=False).iloc[0].to_dict()
    expected = {name: int(original[name]) for name in
                ["n_samples", "n_input_pairs", "n_candidate_triplets", "n_tested_triplets", "n_supported_triplets"]}
    alpha = float(original["fdr_alpha"])
    max_distance = int(original["max_cpg_tss_distance_bp"])
    min_samples = int(config["mediation_min_samples"])
    n_pairs = len(inputs)
    pair_index = pd.MultiIndex.from_frame(inputs[KEY])
    if not pair_index.is_unique:
        raise ValueError(f'Duplicate MF-CpG input pairs: {input_path}')
    values = {name: np.zeros(n_pairs, dtype=np.int64) for name in COUNT_COLUMNS}
    values.update({name: np.zeros(n_pairs, dtype=float) for name in SUM_COLUMNS})
    values.update({name: np.full(n_pairs, np.nan) for name in MIN_COLUMNS + ["a_max"]})
    statuses = Counter()
    n_rows, n_cases = 0, 0
    case_path = out_dir / "prad_erg_tdrd1_triplet_case_study" / f"{arm}_{cohort}_ERG_TDRD1_all_triplets.tsv.gz"
    with ExitStack() as stack:
        case_handle = None
        if cohort == "PRAD":
            case_path.parent.mkdir(parents=True, exist_ok=True)
            case_handle = stack.enter_context(gzip.open(case_path, "wt"))
            pd.read_csv(triplet_path, sep="\t", nrows=0).to_csv(case_handle, sep="\t", index=False)
        reader = stack.enter_context(pd.read_csv(triplet_path, sep="\t", dtype=str,
                                                keep_default_na=False, chunksize=chunksize))
        for chunk in reader:
            if chunk.empty:
                continue
            indexes = pair_index.get_indexer(pd.MultiIndex.from_frame(chunk[KEY]))
            if (indexes < 0).any():
                raise ValueError(f'A saved triplet has no corresponding F0 input pair: {triplet_path}')
            tested = chunk.triplet_status.isin(TESTED)
            supported = chunk.mediation_supported.eq("True")
            numeric = {name: pd.to_numeric(chunk[name].replace("", np.nan)) for name in
                       ["G_prime_distance_bp", "a_estimate", "a_pvalue", "b_estimate", "b_pvalue", "joint_adj_p_value"]}
            distances = numeric["G_prime_distance_bp"]
            a, b = numeric["a_estimate"], numeric["b_estimate"]
            work = pd.DataFrame({
                "pair": indexes,
                "n_candidate_triplets": np.ones(len(chunk), dtype=np.int64),
                "n_tested_triplets": tested.to_numpy(dtype=np.int64),
                "n_supported_triplets": supported.to_numpy(dtype=np.int64),
                "n_nominal_joint_triplets": (tested & numeric["a_pvalue"].lt(.05)
                                              & numeric["b_pvalue"].lt(.05)).to_numpy(dtype=np.int64),
                "a_min": a.to_numpy(), "a_max": a.to_numpy(),
                "a_pvalue_min": numeric["a_pvalue"].to_numpy(),
                "joint_q_min": numeric["joint_adj_p_value"].to_numpy(),
                "n_supported_b_positive": (supported & b.gt(0)).to_numpy(dtype=np.int64),
                "n_supported_b_negative": (supported & b.lt(0)).to_numpy(dtype=np.int64),
                "supported_b_sum": b.where(supported, 0).to_numpy(),
                "supported_b_abs_sum": b.abs().where(supported, 0).to_numpy(),
                "supported_abs_distance_sum": distances.abs().where(supported, 0).to_numpy(),
            })
            operations = {name: "sum" for name in COUNT_COLUMNS + SUM_COLUMNS}
            operations.update({name: "min" for name in MIN_COLUMNS})
            operations["a_max"] = "max"
            reduced = work.groupby("pair", sort=False).agg(operations)
            positions = reduced.index.to_numpy()
            for name in COUNT_COLUMNS + SUM_COLUMNS:
                values[name][positions] += reduced[name].to_numpy()
            for name in MIN_COLUMNS:
                values[name][positions] = np.fmin(values[name][positions], reduced[name].to_numpy())
            values["a_max"][positions] = np.fmax(values["a_max"][positions], reduced.a_max.to_numpy())
            statuses.update({str(k): int(v) for k, v in chunk.triplet_status.value_counts().items()})
            if case_handle is not None:
                cases = chunk[chunk.MF_ensembl.eq(ERG_ID) & chunk.G_prime_ensembl.eq(TDRD1_ID)]
                cases.to_csv(case_handle, sep="\t", index=False, header=False)
                n_cases += len(cases)
            n_rows += len(chunk)


    observed = {name: int(values[name].sum()) for name in
                ["n_candidate_triplets", "n_tested_triplets", "n_supported_triplets"]}
    actual_counts = dict(n_input_pairs=n_pairs, **observed)
    if any(actual_counts[key] != expected[key] for key in actual_counts):
        raise ValueError(f'Saved F4 summary and detail rows disagree: {summary_path}')
    testable = original["testable"] == "True"
    result = inputs.copy()
    result.insert(0, "cohort", cohort)
    result.insert(0, "analysis", arm)
    for name, vector in values.items():
        result[name] = vector
    candidate = result.n_candidate_triplets.gt(0)
    tested_pair = result.n_tested_triplets.gt(0)
    supported_pair = result.n_supported_triplets.gt(0)
    result["pair_status"] = np.select(
        [~candidate, ~tested_pair, ~supported_pair],
        ["no_candidate_triplets", "candidate_not_testable", "tested_without_support"],
        default="tested_with_support")
    result = result[PAIR_COLUMNS]
    pair_path = out_dir / "f4_methylation_factor_cpg_pairs" / f"{arm}_{cohort}.tsv.gz"
    metadata_path = pair_path.with_suffix("").with_suffix(".json")
    pair_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(pair_path, sep="\t", index=False, compression={"method": "gzip", "mtime": 0})
    summary = {
        "analysis": arm, "cohort": cohort, "source": str(triplet_path),
        "n_rows": n_rows, "n_input_pairs": n_pairs, "n_exported_pairs": n_pairs,
        "n_candidate_pairs": int(candidate.sum()), "n_tested_pairs": int(tested_pair.sum()),
        "n_supported_pairs": int(supported_pair.sum()),
        "n_no_candidate_pairs": int((~candidate).sum()),
        "n_candidate_mfs": int(result.loc[candidate, "MF_ensembl"].nunique()),
        "n_tested_mfs": int(result.loc[tested_pair, "MF_ensembl"].nunique()),
        "n_supported_mfs": int(result.loc[supported_pair, "MF_ensembl"].nunique()),
        "n_input_mfs": int(inputs.MF_ensembl.nunique()),
        "status_counts": dict(sorted(statuses.items())),
        "pair_status_counts": {str(k): int(v) for k, v in result.pair_status.value_counts().items()},
        "n_samples": expected["n_samples"], "testable": testable,
        "min_samples": min_samples, "fdr_alpha": alpha,
        "max_cpg_tss_distance_bp": max_distance,
        "age_and_sex_columns": original["age_and_sex_columns"],
        "original_summary": original, **observed,
        "pair_scope": "all_F0_input_pairs",
        "n_case_triplets": n_cases,
    }
    json_write(metadata_path, summary)
    return summary
