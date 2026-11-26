#!/usr/bin/env python3
"""Extract identifiers and graph edges from a pinned SO OBO file."""
import hashlib
import json
import sys
from pathlib import Path


def build(path):
    raw = Path(path).read_bytes()
    terms = {}
    version = None
    term = None
    for line in raw.decode('utf-8').splitlines():
        if line.startswith('data-version: '):
            version = line.split(': ', 1)[1]
        if line.startswith('['):
            term = {} if line == '[Term]' else None
        elif term is not None and ': ' in line:
            key, value = line.split(': ', 1)
            if key == 'id':
                terms[value] = term
            elif key in ('name', 'is_obsolete'):
                term[key] = value
            elif key in ('is_a', 'alt_id', 'replaced_by', 'consider'):
                term.setdefault(key, []).append(value.split()[0])
            elif key == 'relationship' and value.startswith('part_of '):
                term.setdefault('part_of', []).append(value.split()[1])
    if not version or 'SO:0000110' not in terms:
        raise ValueError('Missing SO version or sequence_feature root')
    return {'version': version, 'source_sha256': hashlib.sha256(raw).hexdigest(), 'terms': terms}


if __name__ == '__main__':
    data = build(sys.argv[1])
    Path(sys.argv[2]).write_text(json.dumps(data, sort_keys=True, separators=(',', ':'))+'\n')
