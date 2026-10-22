# Remove organelle predictions

`--remove-organelles` excludes all annotation features on sequences explicitly
classified by NCBI as organelles for the exact assembly accession/version.
It is intended for predictions made by nuclear annotation models on assemblies
that also contain mitochondria or plastids. It does not change the assembly,
reannotate organelles or substitute another translation table.

## Online use

```bash
python3 ena_gffprep.py \
  --input /path/to/Abraxas_sylvata.GCA_964187985.1.annevo_2.2.1.Insecta.gff \
  --assembly-accession GCA_964187985.1 \
  --locus-tag-prefix YOURPREFIX \
  --predictor ANNEVO \
  --predictor-version 2.2.1 \
  --predictor-model Insecta \
  --remove-organelles \
  --output-dir /path/to/new-filtered-package
```

Use the real registered prefix. The existing automatic label/filename convention
is unchanged; the flag does not append a filtering suffix. Choose a distinct
`--label` or `--output-dir` if an earlier package already exists. The model is
optional. Normal source-format, provenance and length requirements still apply.

The tool queries:

```text
https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/ACCESSION/sequence_reports
```

All pages are retrieved and counts checked when the API provides a total. A
response for a different accession/version, including a paired GCF/GCA assembly,
is rejected. No nucleotide sequences are downloaded. The request uses Python's
standard library with a 30-second timeout per attempt and up to three attempts
for transient network/server failures (including rate limiting). An initial
online run needs outbound HTTPS access. There is no new Python-package, NCBI
CLI, Java or credential dependency. Ordinary preparation without the flag makes
no organelle metadata request.

## Cache and refresh

Downloaded metadata is cached as `ACCESSION.json` under:

```text
$XDG_CACHE_HOME/ena-gffprep/ncbi-sequences/
```

If XDG_CACHE_HOME is unset, `~/.cache` is used. Use an absolute XDG cache root or
set an explicit directory:

```bash
--remove-organelles --organelle-cache-dir /path/to/ncbi-cache
```

Only complete validated responses are published into the cache. Matching cached
metadata is reused without a network call and has no automatic expiry. The report
records source, retrieval time, file hash and assembly version. Force refresh:

```bash
--remove-organelles --refresh-organelle-cache
```

A refresh failure stops preparation; it does not silently fall back to old data.
A corrupt cache also stops preparation; refresh it explicitly or supply a valid
local metadata file. Metadata can change without a new sequence version, so
review the cache provenance when deciding to refresh.

## Offline metadata

```bash
--remove-organelles --organelle-metadata /path/to/sequence_report.jsonl
```

Accepted formats are a complete NCBI sequence-report API JSON response, a
complete Datasets `sequence_report.jsonl`, or this tool's cached JSON file.
Snake_case API fields and camelCase Datasets fields are accepted. Use the
sequence report, not the assembly-level summary. Each sequence record must carry
the exact assembly accession/version and a positive integer length.

For example, if the optional NCBI Datasets CLI is available on another machine:

```bash
datasets summary genome accession GCA_964187985.1 \
  --report sequence --as-json-lines > sequence_report.jsonl
```

Transfer that complete file to the preparation machine. Partial API pages with
continuation tokens or inconsistent counts are rejected. JSONL does not contain
a completeness marker: supply the complete report, not a hand-filtered list of
organelle sequences. Every annotated GFF sequence must be matched before any rows
are removed. Offline metadata is not supplemented by a network request.

`--organelle-metadata` and `--refresh-organelle-cache` are mutually exclusive.
Metadata/cache options require `--remove-organelles`. Existing `--assembly-report`
still handles aliases and length checks, but is not an offline classification
input for this release; use the JSON/JSONL option for that purpose.

## Classification and matching policy

- Use NCBI's `assigned_molecule_location_type` / `assignedMoleculeLocationType`.
  Explicit supported classes include mitochondrion, chloroplast, plastid,
  chromoplast, apicoplast, cyanelle, leucoplast, proplastid, kinetoplast,
  hydrogenosome, chromatophore and nucleomorph, plus the explicitly listed
  compound classes in `ena_organelles.ORGANELLES`.
- Match exact GenBank accession, RefSeq accession, original sequence name or
  UCSC-style alias to decoded GFF column 1. Accession versions are not stripped.
  Existing output sequence-name mapping behaviour is unchanged.
- Do not infer organelle identity from `MT`, `chrM`, size, gene names or
  `non-nuclear` alone. NCBI records with ambiguous aliases or unknown compartment
  classes fail clearly. Unmatched annotated sequences also fail.
- An unplaced/unlocalized scaffold explicitly assigned to Primary Assembly may
  lack a molecule location. It is retained with a warning and listed as
  unclassified; the tool does not certify that it is nuclear. Other missing
  classifications fail. Misclassified or unlabelled organelle sequence cannot
  be discovered by this metadata-only filter.
- GFF sequence declarations and retained external lengths must agree with NCBI
  lengths. Matching names/lengths still does not prove nucleotide identity.
- Organelle records present in the assembly but absent from the input annotation
  are reported separately and do not count as removed features. A valid report
  containing no annotated organelles produces a successful no-op.

## Annotation and locus-tag handling

All feature rows on matched organelle sequences are excluded, including genes,
transcripts, CDSs, exons and RNA features. Filtering precedes the source adapter,
feature-attribute conversion and new locus-tag allocation. Normal GFF parsing
still runs first, so fundamentally malformed input remains an error even on a
sequence that would subsequently be excluded.

Corresponding sequence-region declarations are omitted from the output.
Position-dependent input comments are discarded; recognised global provenance
and selected global ontology/species directives are retained. Dropped comment
counts are recorded. Retained Parent or Derives_from links into excluded
features, and feature IDs shared across the retained/excluded boundary, fail
rather than produce a broken annotation. If filtering removes all features,
preparation fails and no package is published.

On updates, supply the complete previous `--locus-tag-map`. Its historical entries
remain reserved, including excluded genes. Existing tags found on excluded
features are also reserved by resolving their gene ancestry, even without a map.
An ambiguous tagged ancestry or conflicting assignment fails. New untagged
excluded genes do not receive tags. Do not omit a historical map merely because
you are now removing organelle predictions.

The manifest PRIMARY_ID continues to identify the full assembly. DESCRIPTION
uses the current Annotator / Version / Model wording; filtering details are in
the preparation and exclusion reports rather than automatically changing NAME.

## Outputs

Normal package mode adds `organelle_exclusions.tsv` with these columns:

```text
seqid    ncbi_sequence    location    feature_rows    genes
```

The actual separators are tabs. The file always has a header, including when no
annotated organelles were removed. `seqid` is the original matched GFF ID, while
`ncbi_sequence` is the report's preferred GenBank accession (or available alias).

`preparation.json` adds `organelle_filter`, containing:

- Metadata source (`ncbi_api`, `cache` or `offline_file`), path, SHA-256, retrieval
  timestamp when available and source URLs when available.
- Original, retained and removed feature counts, removed logical gene counts
  and removed counts by feature type.
- Matched sequence aliases and retained metadata lengths.
- Removed annotated organelles, all organelles identified in the assembly,
  excluded gene IDs and unclassified retained scaffold IDs.

`preparation.txt` includes removed feature/gene counts. Explicit-output mode
writes the TSV beside the requested JSON report as `REPORT_STEM.organelles.tsv`.
The sidecar participates in output staging and overwrite checks. Neither mode
creates a separate excluded GFF file. Cached downloads can remain after a failed
preparation for reuse, but incomplete output packages are not published.

## Config, batch and audit

JSON keys are `remove_organelles` (boolean), `organelle_metadata`,
`organelle_cache_dir` and `refresh_organelle_cache` (boolean). File paths resolve
relative to the configuration file. The same names are accepted batch columns;
TSV booleans use `true` / `false`, and row paths resolve relative to the TSV.
The normal config < CLI < nonblank-row precedence applies. A shared cache can
serve many assemblies; each row must supply metadata for its own assembly if
using offline files. Audit-only mode does not accept organelle filtering.

## Verification scope

Automated tests use synthetic metadata and mocked network responses, so the test
suite remains offline. Real-data preparation additionally exercises an NCBI
response for the supplied Abraxas annotation. No request uses Webin credentials,
and metadata retrieval never submits anything to ENA.

NCBI references:

- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-reports/genome-sequence/
- https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/
