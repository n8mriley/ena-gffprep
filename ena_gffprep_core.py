#!/usr/bin/env python3
"""GFF3 preparation engine for ENA annotation workflows.

Rows are represented as dictionaries containing columns (c), decoded attribute
lists (a), integer start/end positions and original line numbers. Repeated IDs
represent intervals of one logical feature, including discontinuous CDSs.
Unsupported biological representations are rejected rather than guessed."""
import argparse
import csv
import gzip
import hashlib
import json
import os
import tempfile
from contextlib import contextmanager
import re
import sys
from collections import defaultdict
from pathlib import Path
from urllib.parse import quote, unquote

VERSION = '0.6.1'
GENES = {'gene', 'pseudogene'}
GENIC = {'mRNA', 'transcript', 'CDS', 'exon', 'intron', 'five_prime_UTR',
         'three_prime_UTR', 'tRNA', 'rRNA', 'ncRNA', 'lnc_RNA', 'miRNA',
         'snRNA', 'snoRNA', 'misc_RNA', 'primary_transcript', 'prim_transcript',
         'tmRNA', 'precursor_RNA'}
# NCBI genetic-code registry checked 2026-10-15; IDs are not a continuous range.
TRANSLATION_TABLES = {1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23,24,25,26,27,28,29,30,31,32,33}
SUPPORTED_TYPES = {'gene','exon','intron','CDS','mRNA','rRNA','tRNA','ncRNA_gene',
    'protein_coding_gene','pseudogene','repeat_region','gap','five_prime_UTR',
    'three_prime_UTR','recoded_codon','stop_codon','anticodon'}
EXCEPTIONS = {'ribosomal_slippage','trans_splicing','low-quality sequence region',
    'heterogeneous population sequenced','RNA editing','reasons given in citation',
    'rearrangement required for product','annotated by transcript or proteomic data','circular_RNA'}
GAP_TYPES = {'between scaffolds','centromere','contamination','heterochromatin',
    'repeat within scaffold','short arm','telomere','within scaffold'}

class Invalid(ValueError):
    """An input or preparation constraint failed before successful publication."""
    pass

def require(ok, message):
    """Raise Invalid with the supplied message when a preparation check fails."""
    if not ok:
        raise Invalid(message)

def opened(path, mode='rt'):
    """Open ASCII text, choosing gzip from the .gz filename suffix."""
    return gzip.open(path, mode, encoding='ascii') if str(path).endswith('.gz') else open(path, mode, encoding='ascii')

def sha256(path):
    """Return the SHA-256 digest of raw file bytes, including gzip compression."""
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


@contextmanager
def staged_outputs(paths):
    """Stage named output files and publish them without overwriting destinations.
    
    Yield a mapping of output names to adjacent temporary paths. Preserve .gz
    suffixes so writers select compression correctly. On failure, remove temporary
    files and outputs published by this call. Publication uses exclusive hard links;
    multiple files do not constitute a filesystem-wide atomic transaction."""
    staged, published = {}, []
    try:
        for name in paths:
            destination = Path(name)
            fd, temp = tempfile.mkstemp(prefix='.ena_gffprep-', suffix='.gz' if str(name).endswith('.gz') else '.tmp', dir=destination.parent)
            os.close(fd)
            staged[name] = Path(temp)
        yield staged
        for name in paths:
            os.link(staged[name], name)  # exclusive publication: never overwrite
            published.append(name)
    except BaseException:
        for name in published:
            if Path(name).exists() and os.path.samefile(name, staged[name]):
                Path(name).unlink()
        raise
    finally:
        for temp in staged.values():
            temp.unlink(missing_ok=True)


def decode(s):
    """Decode a GFF3 escaped value after checking percent escapes and ASCII controls."""
    require(not re.search(r'%(?![0-9A-Fa-f]{2})', s), f'Malformed percent escape: {s!r}')
    value = unquote(s, encoding='ascii', errors='strict')
    require(all(32 <= ord(c) <= 126 for c in value), 'Encoded control/non-ASCII character')
    return value

def encode(s):
    """Percent-escape a GFF3 value while retaining the permitted safe characters."""
    return quote(s, safe=':._-+|/@()[] ')

def attributes(s):
    """Parse column 9 into attribute names mapped to lists of decoded values.
    
    Split comma-separated values before decoding so an escaped comma stays literal.
    Reject duplicate keys and empty values; a dot represents no attributes."""
    out = {}
    if s == '.':
        return out
    for pair in s.rstrip(';').split(';'):
        require('=' in pair, f'Expected GFF3 key=value, got {pair!r}; GTF/GFF2 requires a separate adapter')
        key, value = pair.split('=', 1)
        require(bool(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', key)), f'Invalid attribute key {key!r}')
        require(key not in out and value != '', f'Duplicate/empty attribute {key!r}')
        values = [decode(x) for x in value.split(',')]
        require(all(values), f'Empty list member in {key}')
        out[key] = values
    return out

def render(a):
    """Render attributes with ID, Parent and locus_tag first, escaping each value."""
    keys = [k for k in ('ID', 'Parent', 'locus_tag') if k in a]
    keys += sorted(set(a) - set(keys))
    return ';'.join(k + '=' + ','.join(encode(v) for v in a[k]) for k in keys) or '.'

def scalar(a, key, default=None):
    """Return a single-valued attribute or default; reject multiple values."""
    if key not in a:
        return default
    require(len(a[key]) == 1, f'{key} must have one value')
    return a[key][0]

def parse(path):
    """Parse supported GFF3 into (rows, sequence_regions, comments).
    
    Each row retains columns, decoded attribute lists, integer coordinates and its
    physical input line number. Reject malformed records and embedded FASTA.
    This parser alone does not establish biological validity or sequence identity."""
    rows, regions, comments = [], {}, []
    with opened(path) as f:
        for line_no, line in enumerate(f, 1):
            line = line.rstrip('\r\n')
            if not line:
                continue
            require(all((32 <= ord(c) < 127) or c == '\t' for c in line), f'Line {line_no}: control character')
            require(not line.startswith('##FASTA'), 'Embedded FASTA detected. Supply annotation-only GFF3; nucleotide FASTA is forbidden and peptide FASTA is not handled by this version.')
            directive = line.split()[0] if line.startswith('#') else ''
            if directive == '##sequence-region':
                fields = line.split()
                require(len(fields) == 4, f'Line {line_no}: invalid sequence-region')
                seq, lo, hi = decode(fields[1]), int(fields[2]), int(fields[3])
                require(lo == 1 and hi >= lo, 'This version requires sequence regions beginning at 1')
                require(seq not in regions, f'Duplicate sequence-region {seq}')
                regions[seq] = hi
            elif directive == '##gff-version':
                require(bool(re.fullmatch(r'##gff-version\s+3(?:\.[0-9]+){0,2}', line)), f'Line {line_no}: only GFF3 is supported')
            elif directive == '##species':
                require(len(line.split()) == 2, f'Line {line_no}: malformed species directive')
                comments.append(line)  # retained for metadata consistency checks, not copied verbatim
            elif directive == '#!insdc-gff-version':
                require(bool(re.fullmatch(r'#!insdc-gff-version\s+[0-9]+(?:\.[0-9]+){0,2}', line)), f'Line {line_no}: invalid INSDC version')
            elif line == '###':
                continue
            elif line.startswith('#'):
                comments.append(line)
            else:
                c = line.split('\t')
                require(len(c) == 9, f'Line {line_no}: expected nine TAB-separated columns')
                c[0] = decode(c[0])
                require(c[0] and not re.search(r'\s', c[0]), f'Line {line_no}: invalid seqid')
                require(bool(re.fullmatch(r'[0-9]+', c[3])) and bool(re.fullmatch(r'[0-9]+', c[4])), f'Line {line_no}: coordinates must be integers')
                require(c[1] and c[2] and not re.search(r'\s', c[1] + c[2]), f'Line {line_no}: invalid source/type')
                start, end = int(c[3]), int(c[4])
                require(1 <= start <= end, f'Line {line_no}: invalid coordinates')
                require(c[6] in {'+', '-', '.', '?'}, f'Line {line_no}: invalid strand')
                require(c[7] in {'0', '1', '2'} if c[2] == 'CDS' else c[7] in {'.','0','1','2'} if c[2] in {'recoded_codon','stop_codon'} else c[7] == '.', f'Line {line_no}: invalid phase')
                if c[5] != '.':
                    require(__import__('math').isfinite(float(c[5])), f'Line {line_no}: invalid score')
                a = attributes(c[8])
                for key in ('ID', 'Parent', 'Derives_from', 'transcript_id', 'protein_id'):
                    if key in a:
                        require(all(not re.search(r'\s', v) for v in a[key]), f'Line {line_no}: whitespace in {key}')
                for key in ('ID','transcript_id','protein_id'):
                    scalar(a, key)
                rows.append({'c': c, 'a': a, 'start': start, 'end': end, 'line': line_no})
    require(rows, 'No features found')
    return rows, regions, comments

def annevo_lengths(comments):
    """Extract sequence lengths from recognised ANNEVO comments; reject conflicts."""
    sizes = {}
    pattern = re.compile(r'^#\s*-+\s*prediction\s+on\s+sequence\s+number\s+[0-9]+\s*\(\s*length\s*=\s*([0-9]+)\s*,\s*name\s*=\s*([^)]*?\S)\s*\)\s*-+\s*$', re.I)
    for comment in comments:
        match = pattern.fullmatch(comment)
        if match:
            size, seq = int(match[1]), match[2]
            require(size > 0 and (seq not in sizes or sizes[seq] == size), f'Conflicting ANNEVO lengths: {seq}')
            sizes[seq] = size
    return sizes


def adapt_annevo(rows):
    """Normalise numbered ANNEVO CDS segments in place into logical CDS groups.
    
    Segments for one transcript share a CDS ID; original row identifiers are retained
    in original_id. Reject incompatible metadata, collisions and references that
    would make this one-CDS-per-transcript interpretation ambiguous."""
    ids = {scalar(r['a'], 'ID') for r in rows if 'ID' in r['a']}
    refs = {p for r in rows for p in r['a'].get('Parent', [])}
    referenced_other = {v for r in rows for k in ('Derives_from',) for v in r['a'].get(k, [])}
    transcripts = {scalar(r['a'], 'ID') for r in rows if r['c'][2] == 'mRNA'}
    groups = defaultdict(list)
    for r in rows:
        if r['c'][2] == 'CDS':
            parent = scalar(r['a'], 'Parent')
            require(parent in transcripts, f'ANNEVO CDS parent must be one mRNA: line {r["line"]}')
            groups[parent].append(r)
    changed = 0
    for parent, parts in groups.items():
        old_ids = {scalar(r['a'], 'ID') for r in parts}
        # Already-normalized groups are retained, including a single existing interval.
        if len(old_ids) == 1 and old_ids == {parent + '.CDS'}:
            continue
        require(all(old is not None and re.fullmatch(re.escape(parent) + r'\.CDS\.[0-9]+', old) for old in old_ids),
                f'Unexpected ANNEVO CDS IDs for {parent}; refusing to infer product grouping')
        require(len(old_ids) == len(parts), f'Duplicate ANNEVO segment IDs for {parent}')
        require(not old_ids & (refs | referenced_other), f'ANNEVO segment IDs are referenced by other features: {parent}')
        new_id = parent + '.CDS'
        require(new_id not in ids, f'ANNEVO normalized ID collision: {new_id}')
        metadata = [{k: v for k, v in r['a'].items() if k not in {'ID', 'original_id'}} for r in parts]
        require(all(a == metadata[0] for a in metadata), f'Conflicting CDS attributes across ANNEVO segments: {parent}')
        for r in parts:
            require('original_id' not in r['a'], f'original_id already present: {parent}')
            r['a']['original_id'] = r['a']['ID']
            r['a']['ID'] = [new_id]
            changed += 1
    return {'cds_transcripts': len(groups), 'cds_rows_renamed': changed,
            'assumption': 'Explicit ANNEVO adapter: one CDS product per mRNA; original segment IDs retained in original_id.'}


def audit(path):
    """Return descriptive feature counts and CDS arithmetic without preparing output.
    
    CDS intervals are grouped by parent for this audit. Short-length thresholds
    are advisory observations, not evidence that a prediction should be deleted."""
    rows, regions, comments = parse(path)
    counts, per_seq = defaultdict(int), defaultdict(int)
    groups, exons = defaultdict(list), defaultdict(list)
    for r in rows:
        counts[r['c'][2]] += 1
        per_seq[r['c'][0]] += 1
        if r['c'][2] == 'CDS':
            groups[scalar(r['a'], 'Parent', scalar(r['a'], 'ID', 'unparented'))].append(r)
        if r['c'][2] == 'exon':
            exons[scalar(r['a'], 'Parent', 'unparented')].append(r)
    phase_errors, short, nontriplet, nonzero = [], [], [], []
    exon_differences = []
    for transcript, parts in groups.items():
        parts.sort(key=lambda r: r['start'], reverse=parts[0]['c'][6] == '-')
        initial = int(parts[0]['c'][7])
        if initial:
            nonzero.append(transcript)
        for a, b in zip(parts, parts[1:]):
            expected = (3 - ((a['end'] - a['start'] + 1 - int(a['c'][7])) % 3)) % 3
            if expected != int(b['c'][7]):
                phase_errors.append({'transcript': transcript, 'line': b['line'], 'expected': expected, 'observed': int(b['c'][7])})
        total = sum(r['end'] - r['start'] + 1 for r in parts)
        if total < 90:
            short.append({'transcript': transcript, 'seqid': parts[0]['c'][0], 'strand': parts[0]['c'][6], 'cds_length_nt': total, 'segments': len(parts)})
        if (total - initial) % 3:
            nontriplet.append(transcript)
        if sorted((r['start'],r['end']) for r in parts) != sorted((r['start'],r['end']) for r in exons[transcript]):
            exon_differences.append(transcript)
    headers = annevo_lengths(comments)
    return {'tool_version': VERSION, 'mode': 'audit_only_no_submission_metadata',
            'input_name': Path(path).name, 'input_sha256': sha256(path),
            'feature_rows': len(rows), 'feature_counts': dict(counts), 'annotated_sequences': len(per_seq),
            'feature_rows_by_sequence': dict(per_seq), 'sequence_region_directives': len(regions),
            'annevo_header_sequences': len(headers), 'annevo_header_lengths': headers,
            'gene_rows_with_locus_tag': sum(r['c'][2] == 'gene' and 'locus_tag' in r['a'] for r in rows),
            'cds_rows_with_product': sum(r['c'][2] == 'CDS' and 'product' in r['a'] for r in rows),
            'cds_rows_with_transl_table': sum(r['c'][2] == 'CDS' and 'transl_table' in r['a'] for r in rows),
            'transcripts_with_cds': len(groups), 'phase_continuity_errors': phase_errors,
            'nonzero_initial_phase_transcripts': nonzero, 'nontriplet_cds_transcripts': nontriplet,
            'cds_exon_interval_difference_transcripts': exon_differences,
            'shortest_cds_segment_nt': min((r['end']-r['start']+1 for p in groups.values() for r in p), default=None),
            'cds_under_90_nt_count': len(short), 'cds_under_90_nt': short,
            'cautions': ['The 90 nt threshold is a descriptive review threshold, not an ENA validation rule.',
                         'Comment lengths are unverified against the assembly.',
                         'Phase arithmetic does not establish correct start/stop codons or protein translation.',
                         'Interval comparison does not infer UTRs or validate whether stop codons are included.',
                         'Audit mode is descriptive; run preparation for hierarchy/bounds/tag checks and official Webin for acceptance.']}


def validate_attributes(r):
    """Check the locally supported attributes on one parsed row.
    
    Reject unsupported specialist representations rather than silently rewriting
    their biological meaning. This is a local subset check, not official validation."""
    a, c = r['a'], r['c']
    for key in ('locus_tag','transl_table','pseudo','pseudogene','product','amino_acid','codon_redefined'):
        scalar(a,key)
    if 'pseudo' in a:
        require(a['pseudo'] == ['true'], f'Line {r["line"]}: pseudo must be true')
    if 'exception' in a:
        require(set(a['exception']) <= EXCEPTIONS, f'Line {r["line"]}: unrecognized exception value; use trans_splicing with an underscore')
    require(not ({'part','location'} & set(a)) and not ({'trans_splicing','circular_RNA'} & set(a.get('exception', []))),
            f'Line {r["line"]}: explicit interval ordering/remote/circular annotation requires specialist handling; refusing to reorder')
    if c[2] == 'gap':
        estimated = scalar(a, 'estimated_length')
        require(estimated is not None and (estimated == 'unknown' or bool(re.fullmatch(r'[0-9]+',estimated))), f'Line {r["line"]}: gap requires estimated_length=integer or unknown')
        if 'gap_type' in a:
            require(scalar(a,'gap_type') in GAP_TYPES, f'Line {r["line"]}: invalid gap_type')
    if c[2] == 'recoded_codon':
        require('codon_redefined' in a and 'Parent' in a, f'Line {r["line"]}: recoded_codon needs codon_redefined and Parent')
    if c[2] == 'anticodon':
        require('amino_acid' in a and 'Parent' in a, f'Line {r["line"]}: anticodon needs amino_acid and Parent')
    if 'transcript_id' in a:
        require(c[2] in {'mRNA','transcript','primary_transcript','rRNA','tRNA','ncRNA','lnc_RNA','miRNA','snRNA','snoRNA','misc_RNA'}, f'Line {r["line"]}: transcript_id outside an RNA feature requires source-specific review')
    if 'protein_id' in a:
        require(c[2] == 'CDS', f'Line {r["line"]}: protein_id outside CDS requires source-specific review')


# A repeated CDS ID denotes one discontinuous feature, not duplicate predictions.
def review_cds(nodes, pseudo_nodes, children):
    """Return advisory CDS checks using coordinates, phases and graph metadata.
    
    Internal segment phases describe the reading frame and are not lengths to trim
    from every interval. No FASTA is translated here; specialist exception groups
    are handled conservatively rather than interpreted as ordinary collinear CDSs."""
    result = {'logical_cds_count': 0, 'phase_continuity_errors': [],
              'nontriplet_without_partial_end': [], 'cds_under_90_nt': [],
              'phase_check_skipped_exception_or_pseudogene': [], 'multiple_CDS_ids_per_parent': []}
    per_parent = defaultdict(list)
    for fid, parts in nodes.items():
        if parts[0]['c'][2] != 'CDS':
            continue
        result['logical_cds_count'] += 1
        parent = scalar(parts[0]['a'],'Parent')
        per_parent[parent].append(fid)
        exceptions = {x for r in parts for x in r['a'].get('exception',[])}
        if pseudo_nodes[fid] or exceptions:
            result['phase_check_skipped_exception_or_pseudogene'].append(fid)
            continue
        ordered = sorted(parts,key=lambda r:r['start'],reverse=parts[0]['c'][6] == '-')
        for first, second in zip(ordered,ordered[1:]):
            expected = (3 - ((first['end']-first['start']+1-int(first['c'][7])) % 3)) % 3
            if expected != int(second['c'][7]):
                result['phase_continuity_errors'].append({'CDS':fid,'line':second['line'],'expected':expected,'observed':int(second['c'][7])})
        coding_length = sum(r['end']-r['start']+1 for r in ordered)-int(ordered[0]['c'][7])
        partials = {x for r in parts for x in r['a'].get('partial',[])}
        if coding_length % 3 and 'end' not in partials:
            result['nontriplet_without_partial_end'].append(fid)
        if coding_length < 90:
            result['cds_under_90_nt'].append({'CDS':fid,'locus_tag':scalar(parts[0]['a'],'locus_tag'),
                 'coding_length_nt':coding_length,'short_review_threshold_nt':90,
                 'has_evidence_or_inference': any('inference' in r['a'] or 'experiment' in r['a'] for r in parts)})
    for parent, ids in per_parent.items():
        if parent is not None and len(ids)>1:
            result['multiple_CDS_ids_per_parent'].append({'parent':parent,'CDS_IDs':ids})
    # Validate placement of the supported specialized child types.
    for fid, parts in nodes.items():
        kind = parts[0]['c'][2]
        parent = scalar(parts[0]['a'],'Parent')
        if kind in {'recoded_codon','anticodon'}:
            expected = 'CDS' if kind == 'recoded_codon' else 'tRNA'
            require(parent in nodes and nodes[parent][0]['c'][2] == expected, f'{kind} {fid} requires a {expected} parent')
            require(sum(r['end']-r['start']+1 for r in parts) == 3, f'{kind} {fid} must span three bases in total')
            if kind == 'anticodon':
                exons = [r for child in children[parent] if nodes[child][0]['c'][2] == 'exon' for r in nodes[child]]
                require(exons and all(any(e['start']<=r['start'] and r['end']<=e['end'] for e in exons) for r in parts), f'Anticodon {fid} must be inside tRNA exon intervals')
    return result


def lengths(args, regions):
    """Combine sequence aliases and lengths from regions, FAI and assembly report.
    
    Return alias mappings, sizes and source metadata for coordinate checks.
    Matching names and lengths do not prove sequence identity with an ENA assembly."""
    aliases, sizes, fai_seen = {}, {}, set()
    report_accessions = []
    if args.assembly_report:
        with opened(args.assembly_report) as f:
            for line in f:
                if line.startswith('#'):
                    match = re.match(r'# (?:GenBank|RefSeq) assembly accession:\s+(GC[AF]_[0-9]+\.[0-9]+)', line)
                    if match:
                        report_accessions.append(match[1])
                    continue
                if not line.strip():
                    continue
                c = line.rstrip('\r\n').split('\t')
                require(len(c) >= 9, 'Assembly report needs standard NCBI columns')
                target = c[4] if args.accession_column == 'genbank' else c[6]
                require(target != 'na', f'No {args.accession_column} accession for {c[0]}')
                size = int(c[8])
                require(size > 0 and target not in sizes, f'Duplicate target or invalid length: {target}')
                sizes[target] = size
                for alias in [c[0], c[4], c[6]] + ([c[9]] if len(c) > 9 else []):
                    if alias != 'na':
                        require(alias not in aliases or aliases[alias] == target, f'Ambiguous alias {alias}')
                        aliases[alias] = target
        if args.assembly_accession and report_accessions:
            require(args.assembly_accession in report_accessions, 'Assembly accession disagrees with assembly-report header')
    if args.fai:
        with opened(args.fai) as f:
            for line in f:
                c = line.rstrip('\r\n').split('\t')
                require(len(c) >= 2 and c[0] and not re.search(r'\s',c[0]), 'Malformed FAI row')
                require(bool(re.fullmatch(r'[0-9]+',c[1])), 'FAI length must be an integer')
                if args.assembly_report:
                    require(c[0] in aliases, f'FAI ID absent from assembly report: {c[0]}')
                seq, size = aliases.get(c[0], c[0]), int(c[1])
                require(size > 0 and seq not in fai_seen, f'Duplicate FAI ID or invalid length: {seq}')
                require(seq not in sizes or sizes[seq] == size, f'Length mismatch for {seq}')
                sizes[seq] = size
                fai_seen.add(seq)
    external = bool(args.assembly_report or args.fai)
    for old, size in regions.items():
        seq = aliases.get(old, old)
        if external:
            require(seq in sizes, f'Sequence metadata absent from external length source: {seq}')
            if args.fai:
                require(seq in fai_seen, f'Sequence metadata absent from FAI: {seq}')
        require(seq not in sizes or sizes[seq] == size, f'sequence-region length disagrees for {seq}')
        if not external:
            sizes[seq] = size
    require(sizes, 'No sequence lengths available')
    return aliases, sizes, {'fai_ids': fai_seen, 'assembly_report_accessions': report_accessions}


def prepare(args):
    """Prepare GFF3, locus tags, manifest and reports from an argument namespace.
    
    Optionally remove NCBI-classified organelle predictions before adaptation.
    Validate the retained feature graph, apply supported representation changes, propagate
    submission metadata and stage outputs. Preserve prediction coordinates and
    retain biological QC findings as warnings unless strict warning mode is chosen.
    Return the preparation report; do not authenticate or submit to ENA."""
    require(bool(re.fullmatch('[A-Z][A-Z0-9]{2,11}', args.locus_tag_prefix)), 'Prefix must be 3-12 uppercase alphanumeric characters, beginning with a letter')
    require(args.width > 0 and args.start > 0 and args.step > 0, 'width/start/step must be positive')
    requested_taxon = getattr(args, 'taxon_id', None)
    require(requested_taxon is None or (isinstance(requested_taxon, int) and not isinstance(requested_taxon, bool) and requested_taxon > 0), 'Taxon ID must be a positive integer when supplied')
    inference = getattr(args, 'inference', None)
    if inference is not None:
        require(bool(inference.strip()) and all(32 <= ord(c) <= 126 for c in inference), '--inference must be nonempty ASCII without controls')
    require(args.translation_table in TRANSLATION_TABLES, 'Translation table is not in the supported NCBI registry')
    if args.assembly_accession:
        require(bool(re.fullmatch(r'GC[AF]_[0-9]+\.[0-9]+', args.assembly_accession)), 'Assembly accession must be GCA_/GCF_ accession.version')
    require(isinstance(getattr(args, 'remove_organelles', False),bool), 'remove_organelles must be a boolean')
    remove_organelles = bool(getattr(args, 'remove_organelles', False))
    organelle_report = (str(getattr(args, 'organelle_report', None) or Path(args.report).with_suffix('.organelles.tsv'))
                        if remove_organelles else None)
    outputs = [args.output, args.write_locus_tag_map, args.report] + ([args.manifest] if args.manifest else [])
    if organelle_report:
        outputs.append(organelle_report)
    inputs = [args.input, args.locus_tag_map, args.assembly_report, args.fai]
    if remove_organelles and getattr(args, 'organelle_metadata', None):
        inputs.append(args.organelle_metadata)
    resolved = [Path(p).resolve() for p in outputs]
    require(len(set(resolved)) == len(resolved), 'Output paths must be distinct')
    require(not set(resolved) & {Path(p).resolve() for p in inputs if p}, 'Outputs must not overwrite inputs; use a new mapping version')
    require(not any(p.exists() for p in resolved), 'An output already exists; choose new paths')
    require(all(p.parent.is_dir() for p in resolved), 'Every output parent directory must exist')
    if args.manifest:
        require(bool(args.assembly_accession), '--assembly-accession required with --manifest')
        for value in [args.name, args.description, args.assembly_accession, str(Path(args.output).resolve())]:
            require(not any(c in value for c in '\r\n\t'), 'Manifest values must be single fields')
    rows, regions, comments = parse(args.input)
    require(not any(c.startswith('#!transl_table') for c in comments), 'File-level transl_table pragma requires explicit per-CDS translation tables before preparation; this version does not interpret organelle scope')
    input_species = []
    input_taxa = set()
    for comment in comments:
        if comment.startswith('##species'):
            match = re.search(r'(?:/taxonomy/|/tax-id/|[?&]id=)([0-9]+)(?:/|$|&)', comment)
            require(match is not None and int(match[1]) > 0, 'Input species URI is unrecognized or has an invalid taxon ID')
            input_taxa.add(int(match[1]))
            input_species.append(comment)
    require(len(input_taxa) <= 1, 'Conflicting input species directives')
    require(requested_taxon is None or not input_taxa or input_taxa == {requested_taxon},
            'Input species URI conflicts with --taxon-id')
    effective_taxon = requested_taxon if requested_taxon is not None else next(iter(input_taxa), None)
    species_directive = (f'##species https://www.ncbi.nlm.nih.gov/datasets/taxonomy/{requested_taxon}/'
                         if requested_taxon is not None else input_species[0] if input_species else None)
    taxonomy_source = 'explicit' if requested_taxon is not None else 'input_header' if input_species else 'not_provided'
    comments = [c for c in comments if not c.startswith('##species')]
    source_format = getattr(args, 'source_format', 'generic')
    comment_sizes = annevo_lengths(comments) if getattr(args, 'use_annevo_lengths', False) else {}
    if getattr(args, 'use_annevo_lengths', False):
        require(source_format == 'annevo' and comment_sizes, '--use-annevo-lengths requires --source-format annevo and recognized headers')
        for seq, size in comment_sizes.items():
            require(seq not in regions or regions[seq] == size, f'ANNEVO header/directive length mismatch: {seq}')
            regions[seq] = size
    organelle_filter, excluded_assignments = None, []
    if remove_organelles:
        from ena_organelles import filter_rows
        # Filter before adapting CDSs or allocating new locus tags.
        rows, regions, comments, organelle_filter, excluded_assignments = filter_rows(rows, regions, comments, args)
    adapter = adapt_annevo(rows) if source_format == 'annevo' else {}
    aliases, sizes, length_checks = lengths(args, regions)
    if organelle_filter:
        for seq, expected in organelle_filter['matched_sequence_lengths'].items():
            require(sizes.get(aliases.get(seq,seq)) == expected, 'External/GFF length disagrees with NCBI metadata: '+seq)
    warnings, notes = [], []
    if organelle_filter and organelle_filter['unclassified_retained_sequences']:
        warnings.append('Organelle status is unspecified for retained primary-assembly scaffolds: '+', '.join(organelle_filter['unclassified_retained_sequences']))
    removed_comments = []
    if comment_sizes and not args.assembly_report and not args.fai:
        warnings.append('Sequence lengths came only from ANNEVO comments; confirm against the exact submitted assembly before submission.')
    if source_format == 'annevo':
        # Reordering invalidates position-dependent prediction comments. Retain only provenance.
        removed_comments = [c for c in comments if not c.startswith(('# This output', '# ANNEVO', '# Citation:'))]
        comments = [c for c in comments if c.startswith(('# This output', '# ANNEVO', '# Citation:'))]
    for comment in comments:
        if comment.startswith(('##','#!')) and not comment.startswith(('##attribute-ontology','##source-ontology')):
            notes.append(f'Uninterpreted input directive retained: {comment}')
    missing_products = 0
    for r in rows:
        c, a = r['c'], r['a']
        validate_attributes(r)
        if c[2] == 'CDS' and inference is not None:
            require('inference' not in a or a['inference'] == [inference], f'Line {r["line"]}: existing inference differs; refusing to overwrite')
            a['inference'] = [inference]
        if args.assembly_report:
            require(c[0] in aliases, f'No assembly report mapping for {c[0]}')
        c[0] = aliases.get(c[0], c[0])
        if args.fai:
            require(c[0] in length_checks['fai_ids'], f'Annotated sequence missing from FAI: {c[0]}')
        require(c[0] in sizes, f'Missing true sequence length for {c[0]}; supply .fai, assembly report, or sequence-region')
        require(r['end'] <= sizes[c[0]], f'Line {r["line"]}: coordinate beyond sequence length; circular virtual coordinates require specialist handling')
        if c[2] in {'ncRNA_gene', 'protein_coding_gene'}:
            c[2] = 'gene'
        if c[2] == 'pseudogene':
            c[2] = 'gene'
            a.setdefault('pseudo', ['true'])
        if c[2] not in SUPPORTED_TYPES and c[2] not in {'region','chromosome','contig','scaffold'}:
            warnings.append(f'Line {r["line"]}: feature type {c[2]} needs official vocabulary validation')
        if 'partial' in a:
            require(set(a['partial']) <= {'start', 'end'} and len(a['partial']) == len(set(a['partial'])), f'Line {r["line"]}: ambiguous partial flag; provide biological start/end explicitly')
    # Build graph by logical feature, retaining legal repeated IDs for split intervals.
    nodes = defaultdict(list)
    for i, r in enumerate(rows):
        fid = scalar(r['a'], 'ID', f'__anonymous_{i}')
        require(not fid.startswith('__anonymous_') or 'ID' not in r['a'], 'Reserved ID prefix')
        r['key'] = fid
        nodes[fid].append(r)
    referenced = {p for r in rows for p in r['a'].get('Parent', [])}
    for fid, parts in nodes.items():
        signatures = {(r['c'][0], r['c'][2], r['c'][6], tuple(r['a'].get('Parent', []))) for r in parts}
        require(len(signatures) == 1, f'Inconsistent repeated ID {fid}; remote/trans-spliced features require specialist handling')
        for attr in ('protein_id','transcript_id','transl_table','product','pseudo','pseudogene','exception'):
            values = {tuple(r['a'][attr]) for r in parts if attr in r['a']}
            require(len(values) <= 1, f'Conflicting {attr} across repeated ID {fid}')
            if values:
                value = list(next(iter(values)))
                for part in parts:
                    part['a'].setdefault(attr, list(value))
        intervals = sorted((r['start'], r['end']) for r in parts)
        require(all(a[1] < b[0] for a, b in zip(intervals, intervals[1:])), f'Overlapping/duplicate intervals for repeated ID {fid}')
    derives_refs = {p for r in rows for p in r['a'].get('Derives_from', [])}
    # Shared leaves can be copied safely; shared internal subtrees need a dedicated adapter.
    newnodes = {}
    splits = 0
    for fid, parts in nodes.items():
        parents = parts[0]['a'].get('Parent', [])
        if len(parents) <= 1:
            newnodes[fid] = parts
            continue
        require(len(set(parents)) == len(parents), f'Duplicate Parent on {fid}')
        require(fid not in referenced and fid not in derives_refs, f'Multiple-parent internal or referenced feature {fid}: cannot safely duplicate it')
        for n, parent in enumerate(parents, 1):
            new_id = f'{fid}.ena_parent{n}'
            require(new_id not in nodes and new_id not in newnodes, f'Generated ID collision: {new_id}')
            copies = []
            for r in parts:
                copy = dict(r, c=list(r['c']), a={k:list(v) for k,v in r['a'].items()}, key=new_id)
                copy['a']['Parent'] = [parent]
                if 'ID' in copy['a']:
                    copy['a']['ID'] = [new_id]
                copies.append(copy)
            newnodes[new_id] = copies
        splits += 1
    nodes = newnodes
    children = defaultdict(list)
    indegree = {fid: 0 for fid in nodes}
    for fid, parts in nodes.items():
        parent = scalar(parts[0]['a'], 'Parent')
        if parent:
            require(parent in nodes, f'Missing Parent {parent} for {fid}')
            p = nodes[parent][0]
            require(p['c'][0] == parts[0]['c'][0] and (p['c'][6] == parts[0]['c'][6] or p['c'][6] in '.?'), f'Parent/child seqid or strand mismatch: {fid}')
            lo, hi = min(r['start'] for r in nodes[parent]), max(r['end'] for r in nodes[parent])
            require(all(lo <= r['start'] <= r['end'] <= hi for r in parts), f'Feature {fid} exceeds parent bounds')
            if p['c'][2] == 'CDS':
                require(all(any(x['start'] <= r['start'] and r['end'] <= x['end'] for x in nodes[parent]) for r in parts), f'Feature {fid} is in an intron of its CDS parent')
            children[parent].append(fid)
            indegree[fid] = 1
    def key(fid):
        """Return a deterministic genomic ordering key for a logical feature ID."""
        parts = nodes[fid]
        return (parts[0]['c'][0], min(r['start'] for r in parts), max(r['end'] for r in parts), fid)
    roots = sorted((fid for fid,d in indegree.items() if d == 0), key=key)
    ordered = []
    stack = list(reversed(roots))
    while stack:
        fid = stack.pop()
        ordered.append(fid)
        stack.extend(reversed(sorted(children[fid], key=key)))
    require(len(ordered) == len(nodes), 'Parent cycle detected')
    for fid, parts in nodes.items():
        for ref in parts[0]['a'].get('Derives_from', []):
            require(ref in nodes, f'Missing Derives_from target {ref} for {fid}')
    mapping, used, used_numbers = {}, {}, {}
    pattern = re.compile(re.escape(args.locus_tag_prefix) + r'_([A-Za-z0-9]+)')
    def register(gene, tag):
        """Register a gene-to-locus-tag assignment and reserve its suffix.
        
        Reject collisions so existing tags remain stable and cannot be reallocated."""
        require(gene and not re.search(r'\s', gene), 'Map gene_id must be nonempty without whitespace')
        require(bool(pattern.fullmatch(tag)), f'Tag {tag} must use the registered prefix, an underscore and an alphanumeric suffix (no decimal versions)')
        require(gene not in mapping or mapping[gene] == tag, f'Conflicting locus tags for gene {gene}')
        require(tag not in used or used[tag] == gene, f'Locus tag {tag} reused by different genes')
        suffix = pattern.fullmatch(tag).group(1)
        # Preserve alphanumeric legacy suffixes verbatim; generate only numeric ones.
        if suffix.isdigit():
            number = int(suffix)
            require(number > 0 and (number not in used_numbers or used_numbers[number] == tag), f'Ambiguous numeric tag suffix: {tag}')
            used_numbers[number] = tag
        mapping[gene], used[tag] = tag, gene
    if args.locus_tag_map:
        with opened(args.locus_tag_map) as f:
            reader = csv.DictReader(f, delimiter='\t')
            require({'gene_id', 'locus_tag'} <= set(reader.fieldnames or []), 'Map requires gene_id and locus_tag columns')
            for r in reader:
                register(r['gene_id'], r['locus_tag'])
    # Retain tags already carried by excluded genes even without an input map.
    for gene, tag in excluded_assignments:
        register(gene, tag)
    owners = {}
    for fid in ordered:
        r = nodes[fid][0]
        parent = scalar(r['a'], 'Parent')
        owner = fid if r['c'][2] == 'gene' else owners.get(parent)
        owners[fid] = owner
        require(r['c'][2] != 'gene' or 'ID' in r['a'], 'Every gene needs an explicit stable ID')
        require(r['c'][2] not in GENIC or owner is not None, f'Orphan gene component {fid}: gene parent missing; use a source-specific adapter')
        for part in nodes[fid]:
            existing = scalar(part['a'], 'locus_tag')
            if existing:
                require(part['c'][2] != 'repeat_region', f'repeat_region {fid} must not carry locus_tag; remove that attribute explicitly')
                require(owner is not None, f'Locus tag outside gene hierarchy: {fid}')
                register(owner, existing)
    # Reserve every supplied mapping entry, including genes absent from this annotation.
    # Allocate above the numeric high-water mark; never recycle retired identifiers.
    next_number = max([args.start - args.step] + list(used_numbers)) + args.step
    active = [fid for fid in ordered if nodes[fid][0]['c'][2] == 'gene']
    for fid in active:
        if fid not in mapping:
            register(fid, f'{args.locus_tag_prefix}_{next_number:0{args.width}d}')
            next_number += args.step
    transcript_identifiers, protein_identifiers = {}, {}
    for fid in ordered:
        for r in nodes[fid]:
            for attr, registry in [('transcript_id',transcript_identifiers), ('protein_id',protein_identifiers)]:
                value = scalar(r['a'], attr)
                if value:
                    require(value not in registry or registry[value] == fid, f'{attr} {value} belongs to multiple logical features')
                    registry[value] = fid
    pseudo_nodes = {}
    for fid in ordered:
        a = nodes[fid][0]['a']
        pseudo_nodes[fid] = ('pseudogene' in a or a.get('pseudo') == ['true'] or pseudo_nodes.get(scalar(a, 'Parent'), False))
    for fid in ordered:
        owner = owners[fid]
        for r in nodes[fid]:
            a, c = r['a'], r['c']
            if pseudo_nodes[fid] and 'pseudogene' not in a:
                a['pseudo'] = ['true']
            # DDBJ excludes repeat_region even when it is nested under a gene.
            if owner and c[2] != 'repeat_region':
                a['locus_tag'] = [mapping[owner]]
            if c[2] == 'CDS':
                pseudo = pseudo_nodes[fid]
                if pseudo:
                    a.pop('transl_table', None)
                elif 'transl_table' not in a:
                    a['transl_table'] = [str(args.translation_table)]
                else:
                    require(str(scalar(a, 'transl_table')).isdigit() and int(scalar(a,'transl_table')) in TRANSLATION_TABLES, f'Invalid transl_table for {fid}')
                if 'product' not in a:
                    if args.default_product:
                        a['product'] = [args.default_product]
                    else:
                        missing_products += 1
    if missing_products:
        notes.append(f'{missing_products} CDS rows lack product. No product fallback was requested; acceptance by the intended validator remains unconfirmed. No products were invented.')
    require(active, 'No gene features found')
    active_set = set(active)
    biological_review = review_cds(nodes, pseudo_nodes, children)
    if biological_review['phase_continuity_errors']:
        warnings.append(f"{len(biological_review['phase_continuity_errors'])} CDS phase-continuity errors; no phases changed")
    if biological_review['cds_under_90_nt']:
        warnings.append(f"{len(biological_review['cds_under_90_nt'])} CDS features have less than 90 coding nucleotides (descriptive review threshold only, not a confirmed decoupled-pilot rejection rule); no predictions removed")
    if biological_review['multiple_CDS_ids_per_parent']:
        warnings.append(f"{len(biological_review['multiple_CDS_ids_per_parent'])} parents have multiple CDS IDs; check source-specific grouping")
    if biological_review['nontriplet_without_partial_end']:
        warnings.append(f"{len(biological_review['nontriplet_without_partial_end'])} CDS lengths are not multiples of three without an indicated partial end")
    warnings = sorted(set(warnings))
    report = {'tool_version': VERSION, 'input_sha256': sha256(args.input),
              'assembly_accession': args.assembly_accession, 'source_format': source_format, 'adapter': adapter,
              'sequence_lengths_source': 'assembly_report' if args.assembly_report else 'fai' if args.fai else 'annevo_comments' if comment_sizes else 'sequence_region',
              'genes': len(active),
              'feature_rows': sum(len(v) for v in nodes.values()), 'shared_features_split': splits,
              'reserved_retired_gene_ids': sorted(set(mapping) - active_set), 'warnings': warnings, 'notes': notes,
              'submission_metadata': {'taxon_id': effective_taxon, 'requested_taxon_id': requested_taxon,
                  'taxonomy_source': taxonomy_source, 'species_directive': species_directive, 'assembly_accession': args.assembly_accession,
                  'locus_tag_policy': {'generated_suffix': 'zero-padded positive integer', 'minimum_width': args.width,
                  'existing_suffix': 'ASCII alphanumeric, preserved verbatim', 'numbering_scope': 'whole annotation across sequences and gene types',
                  'retired_tags': 'reserved from supplied mapping; never automatically reused', 'repeat_region': 'no locus_tag'},
                  'locus_tag_prefix': args.locus_tag_prefix, 'translation_table_fallback': args.translation_table,
                  'predictor_model': getattr(args, 'predictor_model', None), 'CDS_inference': inference},
              'input_files_sha256': {str(p): sha256(p) for p in inputs if p},
              'sequence_lengths': {seq: sizes[seq] for seq in sorted({r['c'][0] for r in rows})},
              'feature_types_in_supported_subset': sorted({r['c'][2] for r in rows} & SUPPORTED_TYPES),
              'feature_types_needing_external_vocabulary_check': sorted({r['c'][2] for r in rows} - SUPPORTED_TYPES),
              'biological_review': biological_review,
              'removed_comments': removed_comments,
              'product_policy': 'explicit_fallback' if args.default_product else 'preserve_existing_only',
              'checks_not_performed': ['FASTA-based CDS translation, start/stop codons, internal stops and biological partiality',
                  'Complete SO/INSDC feature and Dbxref vocabulary validation',
                  'All INSDC minimum annotation specifications and overlap rules',
                  'Prefix registration, study ownership, BioSample/taxon/code agreement',
                  'Sequence/assembly association and development-environment availability',
                  'Official gff3tools or Webin validation; no submission performed'],
              'validation_scope': 'Structural preparation only; official vocabulary, sequence translations, minimum biological specifications and Webin acceptance are not validated.'}
    if organelle_filter is not None:
        report['organelle_filter'] = organelle_filter
    if args.strict_warnings:
        require(not warnings, '\n'.join(warnings))
    with staged_outputs(outputs) as staged:
        if organelle_filter is not None:
            from ena_organelles import write_exclusions
            write_exclusions(staged[organelle_report], organelle_filter)
        with opened(staged[args.output], 'wt') as f:
            f.write('##gff-version 3\n#!insdc-gff-version 0.5.1\n')
            if species_directive is not None:
                f.write(species_directive + '\n')
            for seq in sorted({r['c'][0] for parts in nodes.values() for r in parts}):
                f.write(f'##sequence-region {encode(seq)} 1 {sizes[seq]}\n')
            for comment in comments:
                # A single explicit translation table is written to each translating CDS.
                if not comment.startswith('#!transl_table'):
                    f.write(comment + '\n')
            for fid in ordered:
                for r in sorted(nodes[fid], key=lambda r: (r['start'], r['end'])):
                    c = list(r['c'])
                    c[0], c[8] = encode(c[0]), render(r['a'])
                    f.write('\t'.join(c) + '\n')
        with opened(staged[args.write_locus_tag_map], 'wt') as f:
            writer = csv.writer(f, delimiter='\t', lineterminator='\n')
            writer.writerow(['gene_id', 'locus_tag', 'status'])
            for gene in sorted(mapping, key=lambda g: (0, int(pattern.fullmatch(mapping[g]).group(1)))
                               if pattern.fullmatch(mapping[g]).group(1).isdigit() else (1, mapping[g])):
                writer.writerow([gene, mapping[gene], 'active' if gene in active_set else 'reserved'])
        Path(staged[args.report]).write_text(json.dumps(report, indent=2) + '\n', encoding='ascii')
        if args.manifest:
            require(bool(args.assembly_accession), '--assembly-accession required with --manifest')
            with opened(staged[args.manifest], 'wt') as f:
                for k,v in [('NAME',args.name), ('GFF3',os.path.relpath(Path(args.output).resolve(), Path(args.manifest).resolve().parent)), ('ANALYSIS_TYPE','DECOUPLED_ANNOTATION'), ('PRIMARY_ID',args.assembly_accession), ('DESCRIPTION',args.description)]:
                    require(not any(x in v for x in '\r\n\t'), f'Manifest {k} must be a single field')
                    f.write(k + '\t' + v + '\n')
    return report

