"""Retrieve assembly sequence metadata and exclude explicitly classified organelles.

Only this opt-in module makes network requests. It uses the NCBI Datasets v2
sequence report API, not nucleotide sequence downloads. Unknown classifications,
unmatched sequences and partial responses are errors rather than nuclear guesses.
"""
import csv
import hashlib
import json
import os
import re
import tempfile
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import urlencode, quote
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
GENE_TYPES = {'gene', 'pseudogene', 'protein_coding_gene', 'ncRNA_gene'}
API_BASE = 'https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'


def check(condition, message):
    """Raise a preparation-compatible error if a metadata constraint fails."""
    if not condition:
        raise ValueError(message)


def field(record, snake, camel):
    """Read equivalent API snake_case and Datasets JSONL camelCase fields."""
    a, b = record.get(snake), record.get(camel)
    check(a is None or b is None or a == b, 'Conflicting NCBI field spellings: '+snake)
    return a if a is not None else b


def validate_records(records, accession):
    """Normalise and check one complete assembly report; reject wrong versions.

    Missing location on a primary-assembly unplaced/unlocalized scaffold is
    retained as explicitly unclassified, not asserted to be nuclear. Other
    missing or unknown location values require review.
    """
    check(isinstance(records, list) and records, 'NCBI sequence report is empty or invalid')
    result = []
    identities = set()
    for record in records:
        check(isinstance(record, dict), 'Invalid NCBI sequence record')
        get = lambda snake, camel: field(record, snake, camel)
        check(get('assembly_accession', 'assemblyAccession') == accession,
              'NCBI metadata assembly accession/version does not match '+accession)
        aliases = [get('genbank_accession','genbankAccession'),
                   get('refseq_accession','refseqAccession'),
                   get('sequence_name','sequenceName'), get('ucsc_style_name','ucscStyleName')]
        aliases = [a for a in aliases if a and a != 'na']
        check(aliases and all(isinstance(a,str) and not re.search(r'\s',a) for a in aliases),
              'Invalid or missing sequence aliases in NCBI metadata')
        aliases = sorted(set(aliases))
        gb, rs = get('genbank_accession','genbankAccession'), get('refseq_accession','refseqAccession')
        canonical = gb if gb and gb != 'na' else rs if rs and rs != 'na' else aliases[0]
        check(canonical not in identities, 'Duplicate NCBI sequence record: '+canonical)
        identities.add(canonical)
        length = record.get('length')
        check(isinstance(length,int) and not isinstance(length,bool) and length > 0,
              'Invalid NCBI sequence length: '+canonical)
        location = get('assigned_molecule_location_type','assignedMoleculeLocationType') or ''
        check(isinstance(location,str), 'Invalid NCBI molecule location')
        unit = get('assembly_unit','assemblyUnit') or ''
        role = record.get('role','')
        key = location.strip().casefold()
        unclassified = key in {'', 'unknown', 'na', 'un'}
        if unclassified:
            check(str(unit).casefold() == 'primary assembly' and role in {'unplaced-scaffold','unlocalized-scaffold'},
                  'Missing/uncertain organelle classification for '+canonical)
        else:
            check(key in ORGANELLES | NON_ORGANELLES,
                  'Unrecognised NCBI molecule location '+repr(location)+' for '+canonical)
        check(not (str(unit).casefold() == 'non-nuclear' and key not in ORGANELLES | {'plasmid'}),
              'Non-nuclear sequence lacks a recognised compartment: '+canonical)
        result.append(dict(canonical=canonical, aliases=aliases, length=length,
                           location=location or 'unclassified', organelle=key in ORGANELLES,
                           unclassified=unclassified, assembly_unit=unit))
    return result


def request_json(url):
    """Fetch one HTTPS page with a 30-second timeout and three bounded attempts."""
    request = Request(url, headers={'Accept':'application/json', 'User-Agent':'ena-gffprep/0.6.0'})
    for attempt in range(3):
        try:
            with urlopen(request, timeout=30) as response:
                return json.load(response)
        except HTTPError as error:
            if error.code not in {429,500,502,503,504} or attempt == 2:
                raise ValueError(f'NCBI metadata request failed (HTTP {error.code}); use --organelle-metadata for offline preparation') from error
        except (URLError, TimeoutError, OSError) as error:
            if attempt == 2:
                raise ValueError('NCBI metadata request failed; check HTTPS access or supply --organelle-metadata') from error
        except (ValueError, UnicodeError) as error:
            raise ValueError('NCBI returned invalid JSON metadata') from error
        time.sleep(2 ** attempt)


def fetch_report(accession):
    """Fetch every page, checking totals/tokens before accepting the report."""
    reports, urls, tokens = [], [], set()
    token = None
    total = None
    while True:
        query = {'page_size':1000}
        if token:
            query['page_token'] = token
        url = API_BASE+quote(accession,safe='')+'/sequence_reports?'+urlencode(query)
        page = request_json(url)
        check(isinstance(page,dict) and isinstance(page.get('reports'),list), 'Invalid NCBI sequence-report response')
        reports.extend(page['reports']); urls.append(url)
        current_total = page.get('total_count',page.get('totalCount'))
        if current_total is not None:
            check(str(current_total).isdigit(), 'Invalid NCBI sequence count')
            check(total is None or total == int(current_total), 'NCBI sequence count changed during pagination')
            total = int(current_total)
        token = page.get('next_page_token',page.get('nextPageToken'))
        if not token:
            break
        check(isinstance(token,str) and token not in tokens, 'Repeated/invalid NCBI pagination token')
        tokens.add(token)
    check(total is None or total == len(reports), 'Incomplete NCBI sequence report; filtering refused')
    validate_records(reports,accession)
    return dict(schema_version=1, assembly_accession=accession, reports=reports,
                retrieved_at=datetime.now(timezone.utc).isoformat(), source_urls=urls)


def read_report(path):
    """Read API JSON, the tool's cache envelope or Datasets sequence JSONL."""
    text = Path(path).read_text(encoding='utf-8-sig')
    try:
        data = json.loads(text)
    except json.JSONDecodeError:
        data = [json.loads(line) for line in text.splitlines() if line.strip()]
    if isinstance(data,list):
        return {'reports':data}
    check(isinstance(data,dict), 'Metadata must be NCBI sequence JSON or JSONL')
    if 'reports' not in data and ('assemblyAccession' in data or 'assembly_accession' in data):
        return {'reports':[data]}
    check(isinstance(data.get('reports'),list), 'Metadata has no sequence report records')
    check(not data.get('next_page_token',data.get('nextPageToken')), 'Offline metadata contains only a page; supply a complete report')
    total = data.get('total_count',data.get('totalCount'))
    check(total is None or str(total).isdigit() and int(total) == len(data['reports']), 'Offline sequence report count is incomplete')
    return data


def load_metadata(args):
    """Load offline metadata, cached data or NCBI pages, with source provenance.

    Caches use versioned assembly names and atomic replacement. Cache validation
    failures stop the run; explicit refresh is required to replace corrupt data.
    """
    accession = args.assembly_accession
    check(bool(re.fullmatch(r'GC[AF]_[0-9]+\.[0-9]+',accession or '')), 'Organelle removal needs an exact versioned assembly accession')
    offline = getattr(args,'organelle_metadata',None)
    refresh = getattr(args,'refresh_organelle_cache',False)
    check(not offline or not refresh, 'Offline organelle metadata cannot be combined with cache refresh')
    if offline:
        path = Path(offline)
        data = read_report(path)
        source = 'offline_file'
    else:
        default = Path(os.environ.get('XDG_CACHE_HOME',str(Path.home()/'.cache')))/'ena-gffprep'/'ncbi-sequences'
        directory = Path(getattr(args,'organelle_cache_dir',None) or default)
        path = directory/(accession+'.json')
        if path.exists() and not refresh:
            data = read_report(path); source = 'cache'
        else:
            data = fetch_report(accession)
            directory.mkdir(parents=True,exist_ok=True)
            # Publish a complete validated report only; no partial network response is cached.
            temporary = None
            try:
                with tempfile.NamedTemporaryFile('w',encoding='utf-8',dir=directory,delete=False) as handle:
                    temporary = Path(handle.name)
                    json.dump(data,handle,indent=2); handle.write('\n')
                os.replace(temporary,path)
            finally:
                if temporary and temporary.exists():temporary.unlink()
            source = 'ncbi_api'
    check(data.get('assembly_accession',accession) == accession, 'Metadata cache assembly/version mismatch')
    records = validate_records(data['reports'],accession)
    return records, dict(source=source, metadata_file=str(path.resolve()),
                         metadata_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                         retrieved_at=data.get('retrieved_at'), source_urls=data.get('source_urls',[]),
                         assembly_accession=accession, metadata_sequence_count=len(records))


def filter_rows(rows, regions, comments, args):
    """Return retained rows/regions/comments plus summary and excluded rows.

    Match all annotated sequences before filtering. Preserve exact versioned IDs;
    reject cross-boundary feature IDs and references. Keep unclassified primary
    scaffolds with an explicit report entry, rather than declaring them nuclear.
    """
    records, provenance = load_metadata(args)
    aliases = defaultdict(list)
    for record in records:
        for alias in record['aliases']:aliases[alias].append(record)
    matched = {}
    for seq in sorted({r['c'][0] for r in rows}):
        candidates = aliases.get(seq,[])
        check(len(candidates)==1, ('Ambiguous' if candidates else 'Unmatched')+' GFF sequence in NCBI metadata: '+seq)
        matched[seq] = candidates[0]
        check(seq not in regions or regions[seq] == candidates[0]['length'],
              'GFF/NCBI sequence length mismatch: '+seq)
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
    excluded = [dict(seqid=seq, ncbi_sequence=record['canonical'], location=record['location'],
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
                   matched_sequence_lengths={seq:r['length'] for seq,r in sorted(matched.items()) if not r['organelle']},
                   unclassified_retained_sequences=[seq for seq,r in sorted(matched.items()) if r['unclassified']],
                   discarded_comment_count=len(comments)-len(retained_comments),
                   excluded_gene_ids=sorted({v for r in removed if r['c'][2] in GENE_TYPES for v in r['a'].get('ID',[])}))
    return kept, kept_regions, retained_comments, summary, assignments


def write_exclusions(path, summary):
    """Write a header-bearing TSV, including when no annotated organelles occur."""
    with open(path,'w',encoding='ascii',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['seqid','ncbi_sequence','location','feature_rows','genes'],delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(summary['excluded_sequences'])
