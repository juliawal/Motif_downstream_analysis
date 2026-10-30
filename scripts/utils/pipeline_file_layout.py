"""Read pipeline files, including previously compacted inputs with a TSV index."""
import csv
import os
from pathlib import Path

INDEX = Path("original_pipeline_outputs/shared_file_index.tsv")


class PipelineFiles:
    """Resolve an original context path to its stored file."""

    def __init__(self, root):
        self.root = Path(root).resolve()
        self.rows = {}
        index = self.root / INDEX
        if not index.exists():
            return  # Uncompacted server exports remain readable.
        with index.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.rows = {row["context_path"]: row for row in reader}

    def path(self, path):
        candidate = Path(path)
        if not candidate.is_absolute():
            candidate = self.root / candidate
        candidate = Path(os.path.abspath(candidate))
        try:
            context = candidate.relative_to(self.root).as_posix()
        except ValueError:
            return candidate  # Analysis outputs may be outside pipeline_results.
        row = self.rows.get(context)
        return self.root / row["data_path"] if row else candidate

