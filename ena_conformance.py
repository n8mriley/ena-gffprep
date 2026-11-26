"""Additional submission-profile checks; no biological repairs."""

from ena_vocabularies import load_vocabularies, PSEUDOGENE_TYPES


def feature_name(value, require, line):
    from ena_ontology import feature_status, TERMS
    FEATURE_TERMS, SO_ACCESSIONS, _ = load_vocabularies()
    status, identifier = feature_status(value)
    # Preserve explicitly supported legacy names and validator mapping aliases.
    name = SO_ACCESSIONS.get(value, value)
    if status == 'obsolete':
        replacements = TERMS[identifier].get('replaced_by', []) + TERMS[identifier].get('consider', [])
        require(False, f"Line {line}: obsolete SO term {value}; review replacements {replacements}")
    require(status != 'not_sequence_feature', f"Line {line}: {value} is not an SO sequence feature")
    if identifier and name == value and value not in FEATURE_TERMS:
        name = TERMS[identifier]['name']
    require(name in FEATURE_TERMS, f"Line {line}: feature type {value} not in local INSDC vocabulary ({status}); specialist review required")
    return name


AMINO_ACIDS = frozenset('alanine arginine asparagine aspartic_acid cysteine glutamine glutamic_acid glycine histidine isoleucine leucine lysine methionine phenylalanine proline serine threonine tryptophan tyrosine valine selenocysteine pyrrolysine'.split())


def validate_row(row, require):
    """Check explicit amino acid attributes and cross-reference syntax."""
    a = row['a']
    for key in ('amino_acid', 'codon_redefined'):
        if key in a:
            require(all(v in AMINO_ACIDS for v in a[key]), f'Line {row["line"]}: invalid {key} amino acid name')
    _, _, DBXREF_TAGS = load_vocabularies()
    for value in a.get('Dbxref', []):
        tag, sep, identifier = value.partition(':')
        require(bool(tag and sep and identifier and identifier.strip()), f'Line {row["line"]}: invalid Dbxref; expected DBTAG:ID')

        require(tag in DBXREF_TAGS, f'Line {row["line"]}: unrecognized Dbxref database {tag}')
    if "pseudogene" in a:
        require(all(v in PSEUDOGENE_TYPES for v in a["pseudogene"]), f'Line {row["line"]}: invalid pseudogene subtype')


def validate_partials(fid, parts, require):
    """Check explicit partial flags on multipart features in biological order."""
    if len(parts) < 2 or not any(r['a'].get('partial') for r in parts):
        return
    strand = parts[0]['c'][6]
    require(strand in ('+', '-'), f'Partial multipart feature {fid} needs known strand')
    ordered = sorted(parts, key=lambda r: (r['start'], r['end']), reverse=strand == '-')
    for i, row in enumerate(ordered):
        for boundary in row['a'].get('partial', []):
            require((boundary == 'start' and i == 0) or (boundary == 'end' and i == len(ordered)-1), f'Internal partial boundary on {fid}, line {row["line"]}')


def validate_relationship(fid, parts, parent_parts, require):
    """Enforce unambiguous relationships in conventional gene models."""
    kind, parent_kind = parts[0]['c'][2], parent_parts[0]['c'][2]
    if kind == 'stop_codon' and parent_kind == 'CDS':
        strand = parent_parts[0]['c'][6]
        require(strand in ('+', '-'), f'Exceptional stop codon {fid} needs known strand')
        edge = 'end' if strand == '+' else 'start'
        extreme = max if strand == '+' else min
        require(extreme(r[edge] for r in parts) == extreme(r[edge] for r in parent_parts), f'Exceptional stop codon {fid} must be at biological CDS end')
        require(1 <= sum(r['end']-r['start']+1 for r in parts) <= 3, f'Exceptional stop codon {fid} must span 1–3 bases')
    if kind == 'mRNA':
        require(parent_kind == 'gene', f'Invalid parent type {parent_kind} for mRNA {fid}')
    if kind == 'gene':
        require(parent_kind not in {'gene', 'mRNA', 'CDS', 'exon', 'intron'}, f'Invalid parent type {parent_kind} for gene {fid}')
    if kind == 'CDS':
        require(parent_kind not in {'CDS', 'exon', 'intron', 'tRNA', 'rRNA'}, f'Invalid parent type {parent_kind} for CDS {fid}')
    if kind in {'mRNA', 'CDS', 'exon', 'intron'}:
        require(parent_kind not in {'repeat_region', 'gap', 'recoded_codon', 'anticodon', 'stop_codon'}, f'Invalid parent type {parent_kind} for {kind} {fid}')
