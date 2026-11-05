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

## Manual sequence exclusions

Supply one or more exact input sequence IDs immediately after the flag:

```bash
python3 ena_gffprep.py \
  --input annotation.gff3 \
  --assembly-accession GCA_123456789.1 \
  --locus-tag-prefix YOURPREFIX \
  --no-inference \
  --remove-organelles OV123456.1,OV123457.1 \
  --no-ncbi-lengths
```

All accessions above are examples. Replace them with the assembly and sequence
IDs actually used by your annotation. `--no-inference` preserves existing evidence
without adding inference; alternatively supply accurate predictor provenance.

Matching uses exact column-1 IDs or IDs declared in `##sequence-region`, including
accession versions. A header-only sequence can therefore be excluded and reported
with zero feature rows. Empty entries, duplicate IDs and names absent from both
features and sequence declarations are errors. Quote a list if it contains spaces
around commas; whitespace within an ID is invalid.

Manual mode performs no organelle metadata request and does not infer aliases,
classifications or whether retained sequences are nuclear. Supply the complete
list of sequences you intend to exclude; unlisted sequences are retained. The
shared relationship and tag safeguards below apply to both modes.

Sequence lengths are still required for retained annotations. Missing lengths can
trigger an independent NCBI lookup even in manual mode. `--no-ncbi-lengths`
prevents that lookup: use valid input lengths or supply `--fai` or
`--assembly-report` from the actual reference assembly. Complete offline
`--sequence-metadata` is also supported for length fallback when fallback is
not disabled.

Manual IDs cannot be combined with `--organelle-metadata`,
`--organelle-cache-dir` or `--refresh-organelle-cache`. General sequence metadata
and cache options remain available for length lookup. Changing a suppressed
assembly accession to a newer version is not a substitute for confirming that
the annotation matches that reference.

## Automatic metadata and matching

With bare `--remove-organelles`, the exact assembly's NCBI Datasets sequence report is loaded from the shared
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

## Automatic classification policy

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

## Network failures

Automatic filtering requires a complete usable report from cache, an offline
file or NCBI. A network/server error or incomplete response stops preparation;
it is not treated as evidence that the assembly contains no organelles. Retry
later, supply a complete offline sequence report, or use manually confirmed
sequence IDs. Cached reports remain useful during outages; refresh requires a
working service. Manual filtering alone does not disable length lookup.

## Reports

Normal package mode adds `organelle_exclusions.tsv`:

| Column | Meaning |
| --- | --- |
| `seqid` | Original GFF sequence ID. |
| `ncbi_sequence` | Canonical metadata accession in automatic mode; empty in manual mode. |
| `location` | Metadata classification in automatic mode; `manual exclusion` in manual mode. |
| `feature_rows` | Number of removed rows on that sequence. |
| `genes` | Number of removed logical genes. |

The file always has a header. Automatic mode can produce a header-only file
when no annotated organelles are matched. Manual mode can include entries with
zero counts for excluded header-only sequences. Consult removal counts rather
than the presence of rows to determine whether features were removed.

In automatic mode, the `organelle_filter` section of `preparation.json` records
metadata provenance,
matched aliases/lengths, removed counts, excluded gene IDs, all classified
organelles in the assembly and retained unclassified scaffolds. The text summary
includes removal totals. Explicit-output mode writes
`REPORT_STEM.organelles.tsv` beside the JSON report instead.

Manual reports set `organelle_filter.source` to `manual` and record
`requested_sequence_ids` and `classification_scope`. Removed counts and gene IDs
use the same fields as automatic filtering. Retained sequences are reported as
unclassified; the manual list supplies no evidence of their cellular compartment.
The `organelles_in_assembly` field in manual mode lists only the supplied IDs,
not an independent inventory of the assembly's organelles.

The manifest continues to identify the full assembly. Filtering does not change
the default NAME/DESCRIPTION. Cache files can remain after a failed preparation,
but incomplete output packages are not published. No Webin credentials are used.

## Config, batch and audit

Use the same option names with underscores in JSON or TSV, for example
`remove_organelles`, `sequence_metadata` and `sequence_cache_dir`. Set
`remove_organelles` to JSON `true` for automatic filtering, `false` for no
filtering, or a string such as `"OV123456.1,OV123457.1"` for manual filtering.
In TSV, use `true`, `false` or the comma-separated IDs in the corresponding
cell. Paths follow the normal config/row rules. Each offline report
must match that job's assembly. Audit mode does not filter or fetch metadata.

Validate each prepared annotation using [VALIDATION.md](VALIDATION.md).

References:

- [NCBI genome sequence reports](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/)
- [NCBI Datasets REST API](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/)
