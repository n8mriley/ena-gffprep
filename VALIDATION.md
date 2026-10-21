# Release verification — 0.5.5

## Current scope and checks

This release contains GFF preparation and GFF-based review tools. The local
bundled-validator helper accepts a GFF file only.

All 77 retained automated test methods pass:

```bash
python3 -m unittest discover -s tests -v
```

Coverage includes GFF parsing and escaping, ANNEVO adaptation, feature graph
checks, phase arithmetic, stable locus-tag allocation and reservation,
alphanumeric suffixes, repeat-region exclusion, RNA component propagation,
optional taxonomy, inference provenance, package publication, config and batch
handling, and mapping validator diagnostics to the exact GFF line numbers.
The embedded-FASTA rejection test remains because it checks a GFF input boundary.

FAI name/length index support remains. It is useful for GFF coordinate checking
and sequence-region declarations without reading or preparing nucleotide FASTA.

## Historical integration evidence

Earlier releases were tested using decoupled annotation instructions and
Webin-CLI 10.0.0-beta (gff3tools 5.3.2).
A full ilConVacc3 annotation preparation contained 17,451 genes and 282,576 rows.
Local bundled validation returned exit 0 with no errors and 10,725 warnings:
10,693 EXON_LENGTH and 32 CDS_INTRON_LENGTH. The optional-taxonomy regression
also passed with those counts when the species directive was omitted.

These full-data results are inherited evidence, not fresh validation of 0.5.5.
This release was verified with the retained automated suite; no authenticated
Webin validation or submission was performed. The current locus-tag policy has
local regression coverage but no new full-data Webin result.

## External checks

Registered prefix ownership, exact archived sequence association and development
availability, study/sample relationships and archive processing remain external
checks. Unsupported circular or exceptional models need separate review.
Preparation does not establish biological correctness or submission acceptance.
