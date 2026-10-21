# GFF preparation examples

All accessions, prefixes and annotations are placeholders or synthetic fixtures.
Do not submit them.

- `annotations.tsv`: edit paths, assembly accessions, optional taxa and registered prefixes. Paths resolve relative to this TSV; blank cells inherit defaults.
- `nuclear_defaults.json`: shared nuclear-code fallback and output root. Paths resolve relative to the JSON file.
- `pilot.gff3`: synthetic structural GFF3 fixture. Use explicit `--no-inference` because it has no real prediction provenance, and supply your preparation metadata.

The optional `.fai` input is a sequence-name/length index used for GFF bounds and headers, not a sequence-preparation workflow.
