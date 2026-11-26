"""Extract local Dbxref tags from a user-supplied registry HTML file."""
import argparse
import hashlib
from html.parser import HTMLParser
from pathlib import Path
import re
import sys

EXAMPLE = re.compile(r'db_xref\s*=\s*[“”"\s]*([^:\r\n]+):')


class Tables(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.tables = []
        self.table = self.row = self.cell = None

    def handle_starttag(self, tag, attrs):
        if tag == 'table':
            if self.table is not None:
                raise ValueError('Nested table is unsupported')
            self.table = []
        elif self.table is not None:
            if tag == 'tr':
                self.row = []
            elif tag in ('td', 'th') and self.row is not None:
                self.cell = []
            elif tag in ('br', 'p') and self.cell is not None:
                self.cell.append(' ')

    def handle_data(self, data):
        if self.cell is not None:
            self.cell.append(data)

    def handle_endtag(self, tag):
        if tag in ('td', 'th') and self.cell is not None:
            self.row.append(''.join(self.cell).strip())
            self.cell = None
        elif tag == 'tr' and self.row is not None:
            self.table.append(self.row)
            self.row = None
        elif tag == 'table' and self.table is not None:
            self.tables.append(self.table)
            self.table = None


def extract(html):
    parser = Tables()
    parser.feed(html)
    parser.close()
    candidates = [t for t in parser.tables if any(
        len(r) >= 3 and EXAMPLE.search(r[2]) for r in t)]
    if len(candidates) != 1 or parser.table is not None:
        raise ValueError('Expected one complete registry table')
    tags, issues = set(), []
    for row in candidates[0]:
        if len(row) not in (4, 5) or any(row[3:]) or not row[0] or not row[1]:
            raise ValueError('Unexpected registry row structure')
        examples = {x.strip() for x in EXAMPLE.findall(row[2])}
        if not examples or any(not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_./() -]*', x) for x in examples):
            raise ValueError(f'Missing or malformed example for {row[0]!r}')
        if examples != {row[0]}:
            issues.append((row[0], sorted(examples)))
        if tags & examples:
            raise ValueError('Duplicate example tag across rows')
        tags.update(examples)
    if not tags:
        raise ValueError('Empty registry')
    return tags, issues


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('html', type=Path)
    ap.add_argument('output', type=Path)
    ap.add_argument('--accept-example-tags', action='store_true',
                    help='After review, retain all example tags including discrepancies')
    args = ap.parse_args()
    try:
        raw = args.html.read_bytes()
        tags, issues = extract(raw.decode('utf-8'))
        print('HTML SHA-256:', hashlib.sha256(raw).hexdigest())
        for label, examples in issues:
            print(f'Review: {label!r} -> {examples!r}', file=sys.stderr)
        if issues and not args.accept_example_tags:
            raise ValueError('Review discrepancies before using --accept-example-tags')
        data = ('\n'.join(sorted(tags)) + '\n').encode('utf-8')
        with args.output.open('xb') as out:
            out.write(data)
        print(f'{len(tags)} tags; SHA-256: {hashlib.sha256(data).hexdigest()}')
    except (ValueError, OSError) as exc:
        ap.exit(1, f'Error: {exc}\n')


if __name__ == '__main__':
    main()
