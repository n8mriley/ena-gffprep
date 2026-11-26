"""Pinned SO graph; ontology evidence is distinct from INSDC acceptance."""
import json
from functools import lru_cache
from pathlib import Path

path = Path(__file__).with_name('so_snapshot.json')
if not path.is_file():
    raise ValueError('Missing local ontology index; run tools/setup_references.py (see README)')
DATA = json.loads(path.read_text())
TERMS = DATA['terms']
NAMES = {v['name']: k for k,v in TERMS.items()}
ALIASES = {a:k for k,v in TERMS.items() for a in v.get('alt_id', [])}


def resolve(value):
    return value if value in TERMS else ALIASES.get(value, NAMES.get(value))


@lru_cache(maxsize=None)
def ancestors(identifier, relations=('is_a',)):
    seen = set()
    pending = [identifier]
    while pending:
        current = pending.pop()
        if current in seen or current not in TERMS:
            continue
        seen.add(current)
        if TERMS[current].get('is_obsolete') == 'true':
            continue
        for relation in relations:
            pending.extend(TERMS[current].get(relation, []))
    return frozenset(seen)


def is_a(value, root):
    return resolve(root) in ancestors(resolve(value))


def relationship(child, parent):
    """Positive part_of evidence only; absence is not a contradiction."""
    c, p = resolve(child), resolve(parent)
    if not c or not p:
        return 'unresolved_term'
    if any(target in ancestors(p) for a in ancestors(c)
           for target in TERMS[a].get('part_of', [])):
        return 'ontology_supported'
    return 'not_established'


def feature_status(value):
    identifier = resolve(value)
    if not identifier:
        return 'unknown', None
    term = TERMS[identifier]
    if term.get('is_obsolete') == 'true':
        return 'obsolete', identifier
    if 'SO:0000110' not in ancestors(identifier, ('is_a','part_of')):
        return 'not_sequence_feature', identifier
    return 'sequence_feature', identifier
