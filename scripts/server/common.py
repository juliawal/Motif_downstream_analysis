"""Copy inputs and write the data needed by the downstream analysis."""
import gzip
import json
from pathlib import Path
import shutil


def json_write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n")


def copy_file(source, destination, *, compress=False):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with Path(source).open("rb") as src, destination.open("xb") as dst:
        if compress:
            with gzip.GzipFile(filename="", mode="wb", fileobj=dst, mtime=0) as compressed:
                shutil.copyfileobj(src, compressed)
        else:
            shutil.copyfileobj(src, dst)
