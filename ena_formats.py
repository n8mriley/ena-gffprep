"""Input syntax helpers and prediction provenance detection.

Detection records observed sources, never guesses versions from filenames or
model names. Structural normalization is implemented in ena_normalize and never
selects rules by predictor/source name.
"""
import re
from collections import Counter

FORMATS = {'auto', 'generic', 'annevo', 'tiberius', 'carbon-a', 'oriongeno', 'helixer', 'braker', 'gtf', 'gff3'}
SOURCES = {'tiberius': ('tiberius', 'Tiberius'), 'carbon-a': ('carbon-a', 'Carbon-A'),
           'oriongeno': ('oriongeno', 'OrionGeno'), 'helixer': ('helixer', 'Helixer')}


def detect(path, opener):
    """Scan the entire input for dialect signatures and source composition.

    BRAKER component signatures identify a compatible dialect, not proof of the
    pipeline name/version. Explicit predictor metadata takes priority outside
    the strict legacy ANNEVO path.
    """
    counts = Counter(); versions = set(); annevo = False; syntax = set(); models = set()
    with opener(path) as handle:
        for line in handle:
            if line.startswith('#'):
                if re.search(r'\b(?:generated\s+(?:with|by)\s+ANNEVO|ANNEVO\s+is\s+an?\s+ab\s+initio)', line, re.I):
                    annevo = True
                m = re.search(r'\bANNEVO\s+(?:version\s*|v\s*)?([0-9]+(?:\.[0-9]+)+)', line, re.I)
                if m: versions.add(m[1]); annevo = True
                if '.h5' in line and 'Helixer' in line:
                    models.add(line.split()[-1].split('/')[-1])
                continue
            if not line.strip(): continue
            c = line.rstrip('\r\n').split('\t')
            if len(c) != 9: continue  # parser supplies precise diagnostics
            counts[c[1]] += 1
            syntax.add('gtf' if re.match(r'\s*\w+\s+"', c[8]) else 'gff3')
    if len(versions) > 1: raise ValueError('Conflicting ANNEVO versions in input header')
    lower = {token.casefold() for source in counts for token in source.split(',')}
    fmt, predictor = 'generic', None
    if annevo: fmt, predictor = 'annevo', 'ANNEVO'
    elif lower & {'augustus', 'genemark.hmm3', 'gmst', 'stringtie2utr'}: fmt = 'braker'
    else:
        matches = [v for k,v in SOURCES.items() if k in lower]
        if len(matches) == 1: fmt, predictor = matches[0]
        elif syntax == {'gtf'}: fmt = 'gtf'
    assisted = bool(lower & {'transdecoder', 'gmst', 'genemark.hmm3', 'stringtie2utr'}) or fmt == 'braker'
    return {'format':fmt, 'predictor':predictor, 'version':next(iter(versions),None),
            'sources':dict(sorted(counts.items())), 'syntax':sorted(syntax),
            'evidence_assisted_or_mixed':assisted, 'model_hints':sorted(models)}


def gtf_attributes(text):
    """Parse quoted GTF fields without splitting semicolons inside quoted values.

    Repeated attributes are retained as lists; malformed or unquoted fields fail.
    Escaped quotes and backslashes are decoded, without interpreting arbitrary
    Python escape sequences.
    """
    result = {}; pos = 0
    pattern = re.compile(r'\s*([A-Za-z][A-Za-z0-9_]*)\s+"((?:[^"\\]|\\.)*)"\s*(?:;|$)')
    while pos < len(text):
        if not text[pos:].strip(): break
        m = pattern.match(text, pos)
        if not m: raise ValueError('Malformed quoted GTF attribute near '+repr(text[pos:pos+60]))
        value = re.sub(r'\\(["\\])', r'\1', m[2])
        result.setdefault(m[1], []).append(value); pos = m.end()
    return result


def gtf_identity(a, feature, line):
    """Construct explicit GFF3 relationships from GTF gene/transcript identifiers.

    Missing containers are not fabricated: normal graph checks require them.
    Line-based leaf IDs are deterministic for the exact input; the report records
    the source checksum so that input-dependent identifiers can be reproduced.
    """
    def one(key):
        """Require one nonempty GTF identifier for this relationship."""
        vals = a.get(key, [])
        if len(vals) != 1 or not vals[0]: raise ValueError(f'Line {line}: GTF requires one {key}')
        return vals[0]
    if feature == 'gene': a['ID'] = [one('gene_id')]
    elif feature in {'transcript','mRNA'}:
        a['ID'] = [one('transcript_id')]; a['Parent'] = [one('gene_id')]
    else:
        parent = one('transcript_id'); a['Parent'] = [parent]
        a['ID'] = [f'{parent}.ena_{feature}.{line}']
    return a


def adapt(rows, scalar, require):
    """Compatibility entry point for the shared structural normalizer.

    The main CLI calls ena_normalize.normalize directly for every source.
    This helper retains the earlier in-place API for programmatic callers.
    """
    from ena_normalize import normalize
    normalized, counts = normalize(rows, scalar, require)
    rows[:] = normalized
    return counts


def source_inference(detection, predictor=None, version=None):
    """Return an explicit policy for directly identified ab initio sources.

    These five direct sources are classified as ab initio by the submitter.
    Combined sources and other pipeline components receive no new inference.
    An unknown version is omitted, never replaced by a pipeline/model version.
    """
    known = {'annevo':'ANNEVO', 'helixer':'Helixer', 'tiberius':'Tiberius',
             'carbon-a':'Carbon-A', 'oriongeno':'OrionGeno'}
    mapping = {}
    for source in detection['sources']:
        tool = known.get(source.casefold())
        value = 'ab initio prediction:'+tool if tool else None
        if tool and predictor and predictor.casefold() == tool.casefold() and version:
            value += ':'+str(version)
        mapping[source] = value
    return mapping
