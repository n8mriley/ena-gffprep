"""Retrieve assembly sequence metadata and exclude explicitly classified organelles.

This module supplies shared metadata for organelle filtering and length fallback.
It uses the ENA sequence report HTTPS endpoint, not nucleotide downloads.
Unmatched sequences and partial responses are errors. Organelle filtering also
requires supported compartment classifications; length-only lookup does not.
"""
import csv
import hashlib
import os
import re
import tempfile
import time
from collections import Counter, defaultdict
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

# Explicit cellular-compartment classifications, not sequence-name heuristics.
ORGANELLES = {'mitochondrion', 'chloroplast', 'plastid', 'chromoplast', 'apicoplast',
              'cyanelle', 'leucoplast', 'proplastid', 'kinetoplast', 'hydrogenosome',
              'chromatophore', 'nucleomorph', 'mitochondrion:kinetoplast',
              'plastid:chloroplast', 'plastid:apicoplast', 'plastid:chromoplast',
              'plastid:cyanelle', 'plastid:leucoplast', 'plastid:proplastid'}
NON_ORGANELLES = {'chromosome', 'nuclear', 'nucleus', 'plasmid', 'linkagegroup',
                 'linkage group', 'scaffold', 'contig', 'unplaced', 'unlocalized',
                 'unplaced-scaffold', 'unlocalized-scaffold', 'unlocalized scaffold',
                 'unplaced scaffold', 'na', 'un', 'unknown'}
from ena_vocabularies import GENE_FEATURES
GENE_TYPES = GENE_FEATURES
API_BASE = 'https://ftp.ebi.ac.uk/pub/databases/ena/assembly/'


def check(condition, message):
    """Raise a preparation-compatible error if a metadata constraint fails."""
    if not condition:
        raise ValueError(message)


FIELDS = ('accession','sequence-name','sequence-length','sequence-role','replicon-name','replicon-type','assembly-unit')


def report_url(accession):
    """Build the ENA sequence-report URL for an exact versioned GCA accession."""
    check(bool(re.fullmatch(r'GCA_[0-9]{9}\.[0-9]+', accession or '')), 'ENA metadata requires a versioned GCA accession')
    return f'{API_BASE}{accession[:7]}/{accession[:10]}/{accession}_sequence_report.txt'


def validate_records(records, accession, require_classification=True):
    """Validate report rows and derive aliases, lengths and compartment classifications."""
    check(bool(records), 'ENA sequence report is empty')
    result, seen = [], set()
    for row in records:
        check(all(k in row and isinstance(row[k], str) for k in FIELDS), 'Malformed ENA sequence report row')
        canonical = row['accession']
        check(bool(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*\.[0-9]+', canonical)), 'Invalid versioned sequence accession: '+canonical)
        check(canonical not in seen, 'Duplicate ENA sequence record: '+canonical)
        seen.add(canonical)
        check(row['sequence-length'].isdigit() and int(row['sequence-length']) > 0, 'Invalid ENA sequence length: '+canonical)
        name = row['sequence-name']
        aliases = [canonical] + ([name] if name.casefold() not in {'','na'} and name != canonical else [])
        check(all(not re.search(r'\s', alias) for alias in aliases), 'Invalid sequence alias: '+canonical)
        location = row['replicon-type'].strip()
        key = location.casefold()
        unit = row['assembly-unit'].strip()
        role = row['sequence-role'].casefold().replace('unlocalised','unlocalized')
        unclassified = key in {'', 'unknown', 'na', 'un'}
        if require_classification:
            if unclassified:
                check(unit.casefold() == 'primary assembly' and role in {'unplaced-scaffold','unlocalized-scaffold'}, 'Missing/uncertain organelle classification for '+canonical)
            else:
                check(key in ORGANELLES | NON_ORGANELLES, 'Unrecognised ENA compartment: '+location)
            check(unit.casefold() != 'non-nuclear' or key in ORGANELLES | {'plasmid'}, 'Non-nuclear sequence lacks a recognised compartment: '+canonical)
            check(key not in ORGANELLES or unit.casefold() == 'non-nuclear', 'Contradictory organelle classification: '+canonical)
        result.append(dict(canonical=canonical, aliases=aliases, length=int(row['sequence-length']), location=location or 'unclassified', organelle=key in ORGANELLES, unclassified=unclassified, assembly_unit=unit))
    return result


def parse_report(raw, accession, require_classification=True):
    """Parse the ENA TSV, rejecting malformed headers and incomplete rows."""
    check(raw.endswith(b'\n'), 'ENA report lacks final newline; possible incomplete transfer')
    reader = csv.DictReader(raw.decode('utf-8-sig').splitlines(), delimiter='\t')
    check(reader.fieldnames and len(reader.fieldnames)==len(set(reader.fieldnames)) and set(FIELDS)<=set(reader.fieldnames), 'Expected ENA sequence report TSV header')
    rows = list(reader)
    check(all(None not in row for row in rows), 'Malformed ENA sequence report row')
    return validate_records(rows, accession, require_classification)


def request_report(url):
    """Fetch report bytes, retry transient failures and check transfer length."""
    request = Request(url, headers={'Accept':'text/plain','User-Agent':'ena-gffprep/0.8.0'})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                raw = response.read()
                size = response.headers.get('Content-Length')
                check(size is None or len(raw)==int(size), 'Incomplete ENA report transfer')
                return raw
        except HTTPError as error:
            if error.code not in {429,500,502,503,504} or attempt==2:
                raise ValueError(f'ENA metadata request failed (HTTP {error.code})') from error
        except (URLError, TimeoutError, OSError) as error:
            if attempt==2:
                raise ValueError('ENA metadata request failed; supply --sequence-metadata for offline preparation') from error
        time.sleep(2**attempt)


def load_metadata(args, require_classification=True):
    """Resolve offline/cache/remote metadata, validate it and publish fresh cache atomically."""
    accession = args.assembly_accession
    url = report_url(accession)
    offline = getattr(args,'sequence_metadata',None) or getattr(args,'organelle_metadata',None) or getattr(args,'assembly_report',None)
    refresh = getattr(args,'refresh_sequence_cache',False) or getattr(args,'refresh_organelle_cache',False)
    check(not offline or not refresh, 'Offline sequence metadata cannot be combined with cache refresh')
    if offline:
        path = Path(offline)
        check(path.name == accession+'_sequence_report.txt', 'Offline report filename must match assembly accession/version: '+accession)
        raw = path.read_bytes(); source = 'offline_file'
    else:
        default = Path(os.environ.get('XDG_CACHE_HOME',str(Path.home()/'.cache')))/'ena-gffprep'/'ena-sequences'
        directory = Path(getattr(args,'sequence_cache_dir',None) or getattr(args,'organelle_cache_dir',None) or default)
        path = directory/(accession+'_sequence_report.txt')
        if path.exists() and not refresh:
            raw = path.read_bytes(); source = 'cache'
        else:
            raw = request_report(url); source = 'ena_ftp'
    records = parse_report(raw,accession,require_classification)
    if source == 'ena_ftp':
        directory.mkdir(parents=True,exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile('wb',dir=directory,delete=False) as handle:
                temporary = Path(handle.name); handle.write(raw)
            os.replace(temporary,path)
        finally:
            if temporary and temporary.exists(): temporary.unlink()
    return records, dict(source=source, metadata_file=str(path.resolve()), metadata_sha256=hashlib.sha256(raw).hexdigest(), source_urls=[url] if not offline else [], assembly_accession=accession, metadata_sequence_count=len(records))


def manual_sequence_ids(value):
    """Validate optional exact sequence IDs; booleans retain automatic/off modes."""
    check(isinstance(value, (bool, str)), 'remove_organelles must be a boolean or comma-separated sequence IDs')
    if isinstance(value, bool):
        return None
    names = [part.strip() for part in value.split(',')]
    check(all(names) and all(all(33 <= ord(c) <= 126 for c in name) for name in names),
          'Manual organelle list contains an empty or invalid sequence ID')
    check(len(names) == len(set(names)), 'Manual organelle list contains duplicate sequence IDs')
    return set(names)


def filter_rows(rows, regions, comments, args):
    """Return retained rows/regions/comments plus summary and excluded rows.

    Match all annotated sequences before filtering. Preserve exact versioned IDs;
    reject cross-boundary feature IDs and references. Keep unclassified primary
    scaffolds with an explicit report entry, rather than declaring them nuclear.
    """
    manual = manual_sequence_ids(getattr(args, 'remove_organelles', True))
    if manual is not None:
        check(not any(getattr(args, key, None) for key in
                      ['organelle_metadata', 'organelle_cache_dir', 'refresh_organelle_cache']),
              'Manual organelle IDs cannot be combined with organelle metadata/cache options')
        observed = {r['c'][0] for r in rows} | set(regions)
        check(not manual - observed, 'Manual organelle sequence IDs absent from input: '+', '.join(sorted(manual-observed)))
        # Local exact IDs only: do not invent ENA classification, aliases or lengths.
        records = [dict(aliases=[seq], canonical=seq, length=regions.get(seq),
                        organelle=seq in manual, location='manual exclusion' if seq in manual else 'unknown',
                        unclassified=seq not in manual) for seq in sorted(observed)]
        provenance = dict(source='manual', requested_sequence_ids=sorted(manual),
                          classification_scope='Only explicitly listed input IDs are excluded; retained sequences are not classified')
    else:
        records, provenance = load_metadata(args)
        args._ena_sequence_metadata = (records, provenance)  # shared within this preparation only
    aliases = defaultdict(list)
    for record in records:
        for alias in record['aliases']:aliases[alias].append(record)
    matched = {}
    for seq in sorted(({r['c'][0] for r in rows} | manual) if manual is not None else {r['c'][0] for r in rows}):
        candidates = aliases.get(seq,[])
        candidates = [r for r in candidates if r['canonical'] == seq] or candidates
        check(len(candidates)==1, ('Ambiguous' if candidates else 'Unmatched')+' GFF sequence in ENA metadata: '+seq)
        matched[seq] = candidates[0]
        check(seq not in regions or regions[seq] == candidates[0]['length'],
              'GFF/ENA sequence length mismatch: '+seq)
    removed = [r for r in rows if matched[r['c'][0]]['organelle']]
    kept = [r for r in rows if not matched[r['c'][0]]['organelle']]
    check(kept, 'Organelle filtering would remove all annotation features')
    removed_ids = {v for r in removed for v in r['a'].get('ID',[])}
    kept_ids = {v for r in kept for v in r['a'].get('ID',[])}
    check(not removed_ids & kept_ids, 'A feature ID is shared by retained and excluded sequences')
    for row in kept:
        check(not any(ref in removed_ids for key in ['Parent','Derives_from'] for ref in row['a'].get(key,[])),
              f"Retained feature on line {row['line']} references an excluded organelle feature")
    # For removed records carrying old tags, resolve their original gene ownership.
    by_id = {v:r for r in removed for v in r['a'].get('ID',[])}
    assignments = []
    for row in removed:
        tag = row['a'].get('locus_tag')
        if not tag:continue
        check(len(tag)==1, 'Excluded feature has multiple locus tags')
        current, seen = row, set()
        while current['c'][2] not in GENE_TYPES:
            parents = current['a'].get('Parent',[])
            check(len(parents)==1 and parents[0] in by_id and parents[0] not in seen,
                  'Cannot reserve excluded locus tag: ambiguous/missing gene ancestry')
            seen.add(parents[0]); current=by_id[parents[0]]
        ids = current['a'].get('ID',[])
        check(len(ids)==1, 'Excluded tagged gene needs one stable ID')
        assignments.append((ids[0],tag[0]))
    counts = Counter(r['c'][0] for r in removed)
    genes = Counter(seq for seq, fid in {(r['c'][0], tuple(r['a'].get('ID',[str(r['line'])])))
                                         for r in removed if r['c'][2] in GENE_TYPES})
    excluded = [dict(seqid=seq, ena_sequence=record['canonical'] if manual is None else '', location=record['location'],
                     feature_rows=counts[seq], genes=genes[seq]) for seq,record in sorted(matched.items()) if record['organelle']]
    # Prediction comments can be positional. Retain only global directives/provenance.
    retained_comments = [c for c in comments if c.startswith(('##species','##attribute-ontology','##source-ontology',
                                                              '# This output','# ANNEVO','# Citation:'))]
    kept_regions = {seq:length for seq,length in regions.items() if seq not in aliases or len(aliases[seq])==1 and not aliases[seq][0]['organelle']}
    summary = dict(provenance, input_feature_rows=len(rows), retained_feature_rows=len(kept),
                   removed_feature_rows=len(removed), removed_genes=sum(genes.values()),
                   removed_by_feature_type=dict(Counter(r['c'][2] for r in removed)),
                   excluded_sequences=excluded,
                   organelles_in_assembly=[dict(sequence=r['canonical'],location=r['location'],aliases=r['aliases']) for r in records if r['organelle']],
                   matched_sequences={seq:r['canonical'] for seq,r in sorted(matched.items())},
                   matched_sequence_lengths={seq:r['length'] for seq,r in sorted(matched.items()) if not r['organelle'] and r['length'] is not None},
                   unclassified_retained_sequences=[seq for seq,r in sorted(matched.items()) if r['unclassified']],
                   discarded_comment_count=len(comments)-len(retained_comments),
                   excluded_gene_ids=sorted({v for r in removed if r['c'][2] in GENE_TYPES for v in r['a'].get('ID',[])}))
    return kept, kept_regions, retained_comments, summary, assignments


def write_exclusions(path, summary):
    """Write a header-bearing TSV, including when no annotated organelles occur."""
    with open(path,'w',encoding='ascii',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['seqid','ena_sequence','location','feature_rows','genes'],delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(summary['excluded_sequences'])


def sequence_lengths(args, seqids, regions):
    """Match exact assembly aliases and return fallback aliases, sizes, provenance.

    Reuse metadata already loaded for organelle filtering in this preparation.
    Length-only lookups do not require cellular-compartment classifications.
    Never strip accession versions or infer lengths from feature coordinates.
    """
    shared = getattr(args, '_ena_sequence_metadata', None)
    records, provenance = shared or load_metadata(args, require_classification=False)
    args._ena_sequence_metadata = (records, provenance)
    index = defaultdict(list)
    for record in records:
        for alias in record['aliases']:
            index[alias].append(record)
    aliases, sizes, targets = {}, {}, {}
    for seq in sorted(set(seqids) | set(regions)):
        candidates = index.get(seq, [])
        candidates = [r for r in candidates if r['canonical'] == seq] or candidates
        check(len(candidates) == 1, ('Ambiguous' if candidates else 'Unmatched')+' GFF sequence in ENA metadata: '+seq)
        record = candidates[0]
        target = record['canonical']
        check(target not in sizes or sizes[target] == record['length'], 'Conflicting ENA lengths for '+target)
        if seq in seqids:
            check(target not in targets or targets[target] == seq,
                  'Multiple annotated sequence names map to one ENA sequence: '+str(target))
            targets[target] = seq
        check(seq not in regions or regions[seq] == record['length'], 'GFF/ENA sequence length mismatch: '+seq)
        aliases[seq], sizes[target] = target, record['length']
    return aliases, sizes, dict(provenance, reused_organelle_metadata=bool(shared),
                               matched_sequences=aliases,
                               matched_sequence_lengths={seq:sizes[target] for seq,target in aliases.items()})
