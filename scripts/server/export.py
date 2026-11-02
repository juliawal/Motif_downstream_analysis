#!/usr/bin/env python3
"""Prepare a completed pipeline run for the paper analyses on the same server."""
import sys
import json
sys.dont_write_bytecode = True

from .inventory import ARMS, branch_path, copy_inputs, copy_panels, load_config
from .associations import export_f4_branch
from .network import export_ppi
from .binding import export_binding
from .validation import PAPER, validate_run, validate_panels, write_source_index


def run_export(root, output, panel_root):
    root, output = root.resolve(strict=True), output.resolve()
    if output.is_relative_to(root) or root.is_relative_to(output):
        raise ValueError("Use an output directory outside the completed pipeline run")
    config = load_config(root)
    paper = json.loads(PAPER.read_text())
    sources = validate_run(root, config) + validate_panels(panel_root, paper=paper)
    output.mkdir(parents=True, exist_ok=False)
    payload = output
    copy_inputs(root, payload, config)
    copy_panels(panel_root, payload, paper)
    for cohort in config["cohorts"]:
        for arm in ARMS:
            print(f"F4 {arm}/{cohort}", flush=True)
            export_f4_branch(root, branch_path(root, arm, cohort), payload, arm, cohort, config)
    export_ppi(root, payload, config, config["cohorts"])
    export_binding(root, payload, config, config["cohorts"])
    write_source_index(sources, payload)
    print(f"Export complete: {payload}", flush=True)
