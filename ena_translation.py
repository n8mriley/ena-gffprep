"""Resolve explicit genetic-code pragmas without guessing organelle codes."""
import re
from collections import defaultdict


def parse_pragma(line, tables, require):
    match = re.fullmatch(r'#!transl_table\s+([A-Za-z_][A-Za-z0-9_]*:[0-9]+(?:,[A-Za-z_][A-Za-z0-9_]*:[0-9]+)*)', line)
    require(match is not None, 'Invalid transl_table pragma')
    result = {}
    for item in match[1].split(','):
        key, value = item.split(':')
        require(key not in result, 'Duplicate transl_table compartment: '+key)
        require(int(value) in tables, 'Invalid transl_table for compartment '+key)
        result[key] = int(value)
    return result


def resolve_compartments(args, seqids, sizes, require):
    from ena_organelles import load_metadata
    shared = getattr(args, '_ena_sequence_metadata', None)
    require(shared or getattr(args, 'assembly_accession', None), 'compartment-specific transl_table requires ENA assembly metadata')
    require(shared or getattr(args, 'sequence_metadata', None) or getattr(args, 'assembly_report', None) or getattr(args, 'organelle_metadata', None) or getattr(args, 'ena_lengths', False), 'compartment-specific transl_table requires --sequence-metadata or enabled ENA lookup')
    records, provenance = shared or load_metadata(args, require_classification=False)
    args._ena_sequence_metadata = (records, provenance)
    index = defaultdict(list)
    for record in records:
        for alias in record['aliases']:
            index[alias].append(record)
    result = {}
    for seq in sorted(seqids):
        candidates = index[seq]
        candidates = [r for r in candidates if r['canonical'] == seq] or candidates
        require(len(candidates) == 1, 'Unmatched or ambiguous translation compartment: '+seq)
        r = candidates[0]
        require(r['length'] == sizes[seq], 'Translation compartment metadata length mismatch: '+seq)
        location = r['location'].casefold()
        unit = r['assembly_unit'].casefold()
        if unit == 'primary assembly' and location in {'chromosome','nuclear','nucleus','linkagegroup','linkage group'}:
            key = 'primary'
        elif unit == 'non-nuclear' and r['organelle']:
            key = location.replace(':', '_')
        else:
            key = None
        result[seq] = key
    return result, provenance
