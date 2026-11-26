"""Local reference data; original preparation policy is separate."""
import json
from pathlib import Path
from functools import lru_cache

PSEUDOGENE_TYPES = frozenset(('allelic', 'processed', 'unitary', 'unknown', 'unprocessed'))

@lru_cache(maxsize=1)
def reference_data():
    path = Path(__file__).with_name('local_vocabularies.json')
    if not path.is_file():
        raise ValueError('Missing local vocabulary index; run tools/setup_references.py (see README)')
    return json.loads(path.read_text())

def effective_dbxref_tags(tags):
    """Apply explicit project compatibility policy without rewriting input."""
    tags = set(tags)
    for pair in (("BOLD", "Bold"), ("InterPro", "Interpro")):
        if tags.intersection(pair):
            tags.update(pair)
    tags.discard("Phytozme")
    return frozenset(tags)

@lru_cache(maxsize=1)
def load_vocabularies():
    data = reference_data()
    mapping = dict(data['mapping'])
    mapping.update({'SO:0000145': 'recoded_codon', 'SO:0000319': 'stop_codon'})
    names = frozenset(mapping.values()) | frozenset('anticodon chromosome contig lnc_RNA misc_RNA prim_transcript protein_coding_gene recoded_codon region scaffold start_codon stop_codon'.split())
    return names, mapping, effective_dbxref_tags(data['dbxref_tags'])

PSEUDOGENE_FEATURES = {"processed_pseudogene": "processed",
                       "non_processed_pseudogene": "unprocessed",
                       "unitary_pseudogene": "unitary",
                       "allelic_pseudogene": "allelic"}
GENE_FEATURES = frozenset({"gene", "pseudogene", "protein_coding_gene", "ncRNA_gene"} | set(PSEUDOGENE_FEATURES))
