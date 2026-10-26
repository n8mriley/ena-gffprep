# Excluding organelle annotations

Use `--remove-organelles` when predictions made with a nuclear annotation model
should be excluded from organelle sequences. It removes all annotated feature
rows on those sequences; it does not change the assembly, reannotate organelles
or substitute another genetic code.

Add the flag to an ordinary preparation command:

```bash
python3 ena_gffprep.py \
  --input annotation.gff3 \
  --assembly-accession GCA_123456789.1 \
  --locus-tag-prefix YOURPREFIX \
  --remove-organelles
```

The normal format and provenance requirements still apply. Replace the example
accession/prefix and add predictor or inference options when required. The flag
does not add a filtering suffix to the output name: choose a new output directory
or label if an earlier package exists.

## Metadata and matching

The exact assembly's NCBI Datasets sequence report is loaded from the shared
cache, fetched over HTTPS or read from a complete offline report. It is also
reused if sequence-length fallback needs it. See the README's
[sequence metadata section](README.md#sequence-lengths-and-ncbi-metadata) for
cache locations, refresh, offline formats and general options.

The general `--sequence-metadata`, `--sequence-cache-dir` and
`--refresh-sequence-cache` options work with filtering. The filtering-specific
`--organelle-metadata`, `--organelle-cache-dir` and `--refresh-organelle-cache`
names remain supported and require the removal flag. General and specific paths
must agree when both are supplied. Offline metadata and refresh are incompatible.

All annotated sequence IDs must match exactly one metadata record using its
GenBank accession, RefSeq accession, original name or UCSC-style alias. Accession
versions are not stripped. Wrong assembly versions, unmatched/ambiguous aliases,
conflicting declared lengths and retained external-length disagreements fail.
A paired GCA/GCF accession is not silently substituted for the requested one.

## Classification policy

Classification uses `assigned_molecule_location_type` (or camelCase equivalent),
not sequence names, length, gene content or `non-nuclear` alone.

Supported organelle locations include mitochondrion, chloroplast, plastid,
chromoplast, apicoplast, cyanelle, leucoplast, proplastid, kinetoplast,
hydrogenosome, chromatophore and nucleomorph. Supported compound forms are
listed in `ena_organelles.ORGANELLES`.

An unplaced/unlocalized scaffold in Primary Assembly with an unspecified
location is retained and reported as unclassified, with a warning. Other
missing or unsupported compartment classifications fail. NCBI misclassification
or absent labels can therefore leave organelles undetected; the filter does not
certify every retained sequence as nuclear.

An organelle present in the assembly but absent from the annotation is reported,
but contributes no removed records. A valid report with no annotated organelles
produces a successful no-op.

## Feature relationships and locus tags

Parsing and relationship construction run first. Filtering then precedes shared
structural normalization and new tag allocation. Fundamentally
malformed input can still fail before exclusion.

The filter removes corresponding output sequence declarations and positional
comments. Global provenance and selected global ontology/species directives are
retained. It refuses feature IDs shared across the retained/excluded boundary,
retained Parent/Derives_from links into excluded features and removal of every
annotation feature.

Excluded genes receive no new tags. On updates, supply the complete previous
locus-tag map: all historical entries remain reserved, including removed genes.
Tags already present on excluded features are also reserved when gene ancestry
can be resolved safely. See [LOCUS_TAGS.md](LOCUS_TAGS.md).

## Reports

Normal package mode adds `organelle_exclusions.tsv`:

| Column | Meaning |
| --- | --- |
| `seqid` | Original GFF sequence ID. |
| `ncbi_sequence` | Matched metadata record's canonical accession, preferring GenBank. |
| `location` | Organelle classification. |
| `feature_rows` | Number of removed rows on that sequence. |
| `genes` | Number of removed logical genes. |

The file always has a header. A header-only file means no annotated organelle
features were removed, not necessarily that the assembly lacks organelles.

The `organelle_filter` section of `preparation.json` records metadata provenance,
matched aliases/lengths, removed counts, excluded gene IDs, all classified
organelles in the assembly and retained unclassified scaffolds. The text summary
includes removal totals. Explicit-output mode writes
`REPORT_STEM.organelles.tsv` beside the JSON report instead.

The manifest continues to identify the full assembly. Filtering does not change
the default NAME/DESCRIPTION. Cache files can remain after a failed preparation,
but incomplete output packages are not published. No Webin credentials are used.

## Config, batch and audit

Use the same option names with underscores in JSON or TSV, for example
`remove_organelles`, `sequence_metadata` and `sequence_cache_dir`. Boolean values
are `true`/`false`; paths follow the normal config/row rules. Each offline report
must match that job's assembly. Audit mode does not filter or fetch metadata.

Validate each prepared annotation using [VALIDATION.md](VALIDATION.md).

References:

- [NCBI genome sequence reports](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/)
- [NCBI Datasets REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)
