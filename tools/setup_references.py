#!/usr/bin/env python3
"""Build local reference indexes from user-acquired files; no downloads."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import zipfile
from build_so_snapshot import build

ROOT = Path(__file__).resolve().parents[1]
SO_HASH = '22a8f3ec2b49125dbb6cee8456f0bca86dd8c98f433165ffa4b554da4f155204'
JAR_HASH = '722b23c530447623597b58f5d3a929bcffe953062480ff3d2fddd6fc864c66cc'

def digest(raw):
    return hashlib.sha256(raw).hexdigest()

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--obo', required=True, type=Path)
    p.add_argument('--jar', required=True, type=Path)
    p.add_argument('--dbxref-tags', required=True, type=Path, help='User-supplied registry tags, one per line')
    a = p.parse_args()
    so = build(a.obo)
    if so['version'] != '2026-08-07' or so['source_sha256'] != SO_HASH:
        p.error('OBO does not match pinned source version and checksum')
    with zipfile.ZipFile(a.jar) as z:
        if 'feature-mapping.tsv' in z.namelist():
            raw = a.jar.read_bytes()
        else:
            entries = [n for n in z.namelist() if n.startswith('BOOT-INF/lib/gff3tools-') and n.endswith('.jar')]
            if len(entries) != 1:
                p.error('Expected exactly one embedded gff3tools JAR')
            raw = z.read(entries[0])
    if digest(raw) != JAR_HASH:
        p.error('gff3tools JAR does not match pinned checksum')
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        rows = csv.DictReader(io.StringIO(z.read('feature-mapping.tsv').decode('utf-8')), delimiter='\t')
        mapping = {r['SOID'].strip(): r['SO term'].strip() for r in rows if r['SOID'].strip() and r['SO term'].strip()}
    tags_raw = a.dbxref_tags.read_bytes()
    tags = sorted({s.strip() for s in tags_raw.decode('utf-8').splitlines() if s.strip() and not s.startswith('#')})
    if not tags or any(':' in t or '\t' in t for t in tags):
        p.error('Expected nonempty database tags, one per line, without identifiers')
    vocab = {'mapping': mapping, 'dbxref_tags': tags, 'jar_sha256': digest(raw), 'dbxref_source_sha256': digest(tags_raw)}
    for name, data in [('so_snapshot.json', so), ('local_vocabularies.json', vocab)]:
        path = ROOT / name
        tmp = path.with_suffix('.tmp')
        tmp.write_text(json.dumps(data, sort_keys=True, separators=(',', ':'))+'\n')
        tmp.replace(path)
    print('Local reference indexes generated; do not commit them.')

if __name__ == '__main__':
    main()
