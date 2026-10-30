"""Read current peak_files.json and the archived manifest layout without downloads."""
import json
from urllib.parse import unquote
from .inventory import configured_path


def read_peak_mapping(root, config):
    directory = configured_path(root, config, 'chipseq_dir')
    current = directory / 'peak_files.json'
    if current.is_file():
        mapping = json.loads(current.read_text())
        entries = {}
        for mf, value in mapping.items():
            if value is None:
                entries[mf] = dict(status='absent', path=None, symbol=None)
                continue
            path = directory / value
            symbol = unquote(path.name.removesuffix('.gz').removesuffix('.bed'))
            entries[mf] = dict(status='available', path=value, symbol=symbol)
        return {'entries': entries, 'source_format': 'peak_files.json',
                'source_mapping': str(current)}, [current]
    manifest = directory / 'manifest.json'
    result = json.loads(manifest.read_text())
    result['source_format'] = 'legacy_manifest'
    sources = [manifest]
    return result, sources


def resolve_peaks(root, config, symbols):
    directory = configured_path(root, config, 'chipseq_dir')
    manifest, sources = read_peak_mapping(root, config)
    result = {}
    for symbol in symbols:
        entry = manifest['entries'][symbol]
        if entry['status'] == 'absent':
            result[symbol] = None
        elif entry['status'] == 'available':
            path = directory / entry['path']
            result[symbol] = path
            sources.append(path)
    return result, sources
