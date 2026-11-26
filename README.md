# ENA GFF preparation tools

[![Python tests](https://github.com/n8mriley/ena-gffprep/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/n8mriley/ena-gffprep/actions/workflows/tests.yml)

[![Release](https://img.shields.io/github/v/release/n8mriley/ena-gffprep)](https://github.com/n8mriley/ena-gffprep/releases) [![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE) [![Python 3.9+](https://img.shields.io/badge/Python-3.9%2B-blue.svg)](#requirements)

Prepare existing genome annotations for ENA decoupled GFF3 submission. Version
**0.8.0** supports several annotation formats, creates a manifest and locus-tag
map, and records the preparation settings. It does not run annotation software,
authenticate with ENA or submit data.

Retained prediction coordinates, strands, scores and CDS phases are preserved.
Representation changes can include feature types, identifiers, attributes and
sequence-name mappings. Optional organelle removal excludes whole sequences'
annotation records. Short predictions are not removed to silence validation.

## Scope and maturity

The primary application is formatting existing deep-learning and ab initio
nuclear genome annotations for ENA submission (e.g. ANNEVO, Helixer, Tiberius
ab initio output, OrionGeno and Carbon-A). Format recognition and shared
normalization are not guarantees of end-to-end submission compatibility.

**Experimental compatibility:** BRAKER and other evidence-based or mixed
annotations have not been fully tested. Their existing evidence and product
relationships need careful review; use the appropriate explicit provenance
policy and validate each output. More comprehensive evidence-based annotation
support and testing are planned. Do not treat successful preparation as proof
that evidence semantics or biological models are correct.

## Quick start

### Install

Clone the initial public release and keep the Python modules together. No
package installation, GPU or Java is needed for preparation.

```bash
git clone https://github.com/n8mriley/ena-gffprep.git
cd ena-gffprep
python3 ena_gffprep.py --version
```

### Try a local example

Run the included synthetic fixture without network access or credentials:

```bash
python3 ena_gffprep.py \
  --input examples/pilot.gff3 \
  --assembly-accession GCA_123456789.1 \
  --bioproject PRJEB123456 \
  --locus-tag-prefix EXAMPLE \
  --no-inference --no-ena-lengths \
  --output-dir example-preparation
```

Inspect `example-preparation/preparation.txt` and `manifest.txt`. The accessions
and prefix are placeholders: **do not submit this package**. Choose a new output
directory for each run.

### Prepare a real annotation

For an automatically recognized annotation with a supported provenance policy:

```bash
python3 ena_gffprep.py \
  --input annotation.gff3 \
  --assembly-accession GCA_123456789.1 \
  --bioproject PRJEB123456 \
  --locus-tag-prefix YOURPREFIX
```

Replace the example accessions and prefix with the exact GCA assembly version,
the study your submission account can reference, and your registered prefix. Input may be plain text or gzip-compressed (`.gz`).
The `.gff` or `.gff3` extension does not determine its syntax.

Add metadata when it is known, for example:

```bash
--predictor Helixer --predictor-version 0.3.3 --predictor-model invertebrate
```

These are additions to the preparation command. The predictor and version must
describe the software actually used. The model is optional and requires both;
it records a model label, not the organism or a model to run. Versions are never
guessed from filenames, model names or the GFF version directive.

For an experimental BRAKER-compatible preparation, identify the pipeline explicitly and choose an
inference policy, for example preserving existing evidence without additions
(replace ACTUAL_VERSION with the version used):

```bash
python3 ena_gffprep.py \
  --input braker.gff3 \
  --assembly-accession GCA_123456789.1 \
  --bioproject PRJEB123456 \
  --locus-tag-prefix YOURPREFIX \
  --predictor BRAKER --predictor-version ACTUAL_VERSION \
  --no-inference
```

For ANNEVO, a recognized versioned header can supply predictor/version. Otherwise
supply them explicitly, or choose an accurate `--inference` or `--no-inference`.
ANNEVO header and inference handling is retained, but structural normalization is now the same for every source.

| Information | Requirement |
| --- | --- |
| Annotation, versioned GCA assembly accession, BioProject, locus-tag prefix | Required for preparation; supply via CLI, configuration or batch row. BioProject syntax is checked locally; study existence and account permissions are checked by Webin. |
| Prediction provenance | Automatic for supported direct sources, with the ANNEVO exception above; other inputs need an explicit policy. See [provenance](#prediction-provenance). |
| Sequence lengths | Required, but an index is optional. Local lengths or exact-assembly ENA fallback can supply them. |
| Previous locus-tag map | Supply the complete cumulative map when updating an annotation. |
| Taxon ID, model, output directory, product fallback | Optional. |

After preparation, review the reports and [validate the generated manifest](VALIDATION.md).

## Documentation

The toolkit ends at package preparation; Webin validation and submission are
separate user-controlled steps.

| Guide | Purpose |
| --- | --- |
| [Examples](examples/README.md) | Synthetic fixture and batch templates |
| [Input syntax and normalization](#input-syntax-and-structural-normalization) | Feature graph and CDS grouping policies |
| [Prediction provenance](#prediction-provenance) | Inference and evidence handling |
| [Organelle filtering](ORGANELLES.md) | Automatic and manual exclusions |
| [Locus tags](LOCUS_TAGS.md) | Stable allocation and annotation updates |
| [Configuration and batches](#configuration-and-batches) | Reproducible multi-annotation runs |
| [Validation](VALIDATION.md) | Webin, local checks and troubleshooting |
| [Compatibility](COMPATIBILITY.md) | Evidence levels and real-data reporting |
| [Changelog](CHANGELOG.md) | Development history |

Standards and archive resources: [ENA submission documentation](https://ena-docs.readthedocs.io/),
[Sequence Ontology GFF3 specification](https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md)
and [INSDC](https://www.insdc.org/). The general GFF3 specification alone does
not establish compliance with ENA's decoupled annotation submission rules.

## Requirements

- Python **3.9+**, using only the standard library; no pip installation or GPU.
  Keep the distributed Python scripts together.
- Supported GFF3 or quoted GTF with explicit gene/transcript containers. Arbitrary
  GFF2, embedded sequences and some specialist locations are not supported.
  Text must meet the scripts' ASCII checks.
- Outbound HTTPS for an uncached ENA lookup, or a complete offline sequence
  report. No external CLI or credentials are needed.
- Java **17+** and a separately obtained compatible Webin JAR for validation.
  Neither is required for Python preparation; the JAR is not included.

## Input syntax and structural normalization

Every source now uses one graph-based normalization engine. It examines feature
relationships, IDs and metadata; column 2 and predictor names never choose its
transformation rules. Unknown annotators receive the same structural treatment
as recognized ones. Source detection remains separate for provenance, model hints
and header metadata.

| Option | Values / behavior |
| --- | --- |
| `--input-format FORMAT` | `auto` (default), `gff3` or `gtf`. Chooses attribute syntax, not annotator. Auto examines the actual rows; mixed GTF/GFF3 attribute syntax is rejected. |
| `--normalization MODE` | `standard` (default) or `strict`. Standard applies shared transformations. Strict skips type/phase/CDS-group rewrites, but still expands shared leaves for single-parent output and runs structural checks. |
| `--cds-grouping POLICY` | `auto` (default), `single-product` or `preserve`. Explicit overrides require standard normalization; see below. |
| `--source-format HINT` | Compatibility option for earlier commands. Named annotators no longer select adapters. Explicit `gtf`/`gff3` supply a syntax default; explicit `generic` supplies a strict-normalization default. Explicit new options take precedence. |

Accepted compatibility hints remain `auto`, `annevo`, `tiberius`, `carbon-a`,
`oriongeno`, `helixer`, `braker`, `gtf`, `gff3` and `generic`. Existing working
commands can be retained; new commands should use the syntax/policy options.
An automatically detected unknown source labeled generic does **not** select
strict mode. Only an explicit generic hint supplies that compatibility default.

GFF3 and quoted GTF parsers produce the same internal graph. GTF needs explicit
gene and transcript container rows; missing containers are not fabricated.
Relationships are constructed from gene/transcript identifiers, with deterministic
leaf IDs for the exact input. Arbitrary GFF2 and embedded sequence sections remain
unsupported.

In standard mode the engine:

- Converts coding `transcript` features to `mRNA`.
- Gives compatible segments one logical CDS ID, retaining renamed IDs in
  `original_id` without discarding earlier provenance.
- Normalizes numeric non-CDS phases to `.`, retaining `original_phase`.
- Checks redundant gene/transcript identifiers against their ancestors and moves
  child transcript identifiers to `original_transcript_id`.
- Duplicates unreferenced shared leaves into single-parent features, retaining
  their original identifiers where present. Shared internal features are refused.

Coordinates, strands, scores, CDS phases, isoforms and record sources are retained.
Exon, intron, UTR and split-codon intervals are preserved. No exons or missing
container features are invented, no exon is extended to cover UTRs, and no
stop-codon bases are appended. Labels such as incomplete/no_stop do not become
biological partial flags. Filtering occurs before these transformations.

### CDS grouping policies

| Policy | Behavior |
| --- | --- |
| `auto` | Keep existing logical groups. Merge compatible, nonoverlapping uniquely identified segment rows under one mRNA; all-anonymous segment sets can receive a logical ID. Preserve groups carrying distinct explicit protein IDs. Refuse ambiguous mixtures of existing multipart groups. |
| `single-product` | Explicitly treat compatible CDS intervals under a transcript as one product, including multiple pre-grouped sets. Conflicting protein/product/evidence metadata, overlaps and referenced segment IDs still fail. |
| `preserve` | Keep supplied CDS IDs and grouping; explicit IDs are required. Other standard normalization still applies. |

Automatic grouping assumes ordinary CDS segment rows with compatible metadata
represent one product per coding transcript. Parent membership alone is not proof:
sequence/strand, intervals, product/evidence attributes, references and conflicting
custom metadata are checked. Distinct products identified by unique protein IDs
are preserved. Use preserve when separate IDs intentionally represent separate
products, or curate ambiguous input. Single-product is not a way to override
biological conflicts. Per-segment `cds_type` metadata can be retained without
preventing grouping.

New grouped IDs use `TRANSCRIPT.CDS`. Existing valid logical IDs, including directly gene-parented CDSs, are retained;
regrouping requires a coding transcript parent. Numbered `TRANSCRIPT.CDS.N` segments, including a singleton, normalize to the
established `TRANSCRIPT.CDS` form irrespective of source. Valid ANNEVO output
conventions are preserved without an ANNEVO transformation branch. Conversely,
already-grouped IDs are no longer rejected just because they do not match an
annotator-specific naming pattern.

`preparation.json` records the shared engine, syntax/policies and rule counts
under `normalization`. The old adapter field is retained as a compatibility
summary of those counts, not as evidence of a selected source adapter. Original
attributes and input hashes support tracing changes and reproducible reruns.

## Prediction provenance

In `--prediction-method auto` mode, direct source labels **ANNEVO, Helixer,
Tiberius, Carbon-A and OrionGeno** are treated as ab initio. In supported
non-ANNEVO provenance policies, CDSs from those direct sources receive:

```text
ab initio prediction:TOOL[:VERSION]
```

An unknown version is omitted. A supplied predictor version is used for a
component only when the predictor name matches that component. Existing
inference is preserved. Combined labels such as `Tiberius,transdecoder` and
other components receive no new inference automatically. Mixed files can
therefore contain both automatically annotated and unchanged evidence fields.

ANNEVO's compatibility provenance policy uses its detected/supplied predictor and version to
construct one inference for CDSs, and rejects conflicting existing values.
Without a usable version, choose explicit inference or no-inference. Detected
ANNEVO provenance conflicts are rejected.

Choose an override when automatic treatment is inappropriate:

| Option | Effect |
| --- | --- |
| `--no-inference` | Add nothing; preserve existing evidence. |
| `--inference "TEXT"` | Apply one accurate inference to every CDS; reject existing conflicting values. |
| `--inference-map evidence.json` | Apply explicit values by exact column-2 source name; preserve evidence for unmapped sources. |
| `--prediction-method ab-initio` | Disable source-specific auto policy and allow a global inference from predictor/version. |
| `--prediction-method evidence-assisted` or `unspecified` | Require an explicit inference, inference map or no-inference choice. |

The three inference override options are mutually exclusive. An inference map
is a nonempty JSON object; each value is an accurate inference string or `null`
to add nothing. Unknown source keys and conflicting existing values are errors:

```json
{
  "Tiberius": "ab initio prediction:Tiberius:ACTUAL_VERSION",
  "transdecoder": null
}
```

Use only the version actually run. The map is hashed and recorded in the report.
Inference describes the prediction method; it does not prove expression or
biological correctness. Choosing no-inference may leave short-CDS evidence
requirements unresolved in the external validator.

## Sequence lengths and ENA metadata

ENA is the only remote assembly metadata provider. Exact versioned GCA sequence
reports are downloaded over HTTPS from the ENA assembly directory. No NCBI
API, JSON/JSONL, RefSeq accession mapping or NCBI assembly report is supported.
Local FAI, sequence-region directives and recognized ANNEVO comments remain
supported; automatic ENA length lookup is needed when local lengths are missing.
Explicit offline reports also enable accession mapping. Missing or ambiguous
aliases and conflicting lengths stop preparation. ENA chromosome names may differ
from original submitted scaffold names; use versioned INSDC accessions directly.

| Option | Default / effect |
| --- | --- |
| `--fai PATH` | Optional local sequence length index. |
| `--assembly-report PATH` | Local ENA TSV report; compatibility name. |
| `--sequence-metadata PATH` | Complete local ENA TSV named `GCA_..._sequence_report.txt` for the exact version. |
| `--no-ena-lengths` | Disable automatic ENA length lookup, not automatic organelle filtering. |
| `--use-annevo-lengths` | Recognized ANNEVO length comments. |
| `--sequence-cache-dir PATH` | Override shared ENA metadata cache. |
| `--refresh-sequence-cache` | Refresh when lookup is required. |

The default cache is `$XDG_CACHE_HOME/ena-gffprep/ena-sequences/` (or
`~/.cache/ena-gffprep/ena-sequences/`). It stores the original TSV with exact
versioned filename. Old NCBI caches are not read. Offline filenames are checked,
but the TSV does not independently encode the assembly accession or total count;
use the original complete ENA download, not a renamed or subset report.
Reports are checked for schema, unique versioned accessions, positive lengths,
and all required annotated IDs. A truncated report ending on a complete row
may not be detectable without an independent checksum or record count.
Automatic filtering and length fallback share one report. Provenance records
path, SHA-256, source and matches. Local metadata cannot establish nucleotide identity.

## Optional annotation metadata and filtering

| Option | Default / effect |
| --- | --- |
| `--taxon-id INTEGER` | Optional positive taxon ID; adds a canonical `##species` URI. Otherwise a recognized input URI is preserved, or the directive is omitted. Taxonomy conflicts fail. |
| `--translation-table INTEGER` | `1`; fallback for translating CDSs without a table. Existing supported per-CDS tables are preserved; pseudogenic CDSs omit them. |
| `--default-product "TEXT"` | No default; fills missing CDS products only when explicitly requested. Existing products are retained. |
| `--remove-organelles [SEQID[,SEQID...]]` | Off when omitted. With no value, identify organelles from ENA metadata. With a value, exclude the exact input sequence IDs listed, without an organelle metadata lookup. Both modes run before normalization and new tag allocation. |
| `--organelle-metadata PATH` | Filtering-specific name for an offline sequence report; requires automatic mode (bare `--remove-organelles`). |
| `--organelle-cache-dir PATH` | Filtering-specific cache override; requires automatic mode (bare `--remove-organelles`). |
| `--refresh-organelle-cache` | Filtering-specific refresh flag; requires automatic mode (bare `--remove-organelles`). |
| `--strict-warnings` | Off; warnings in the preparation report become fatal. Notes are separate, and external validators can still report their own warnings. |

Table 1 is a nuclear assumption, not a universal genetic code. No code is inferred
from taxonomy or sequence names. Supported table IDs are 1, 2, 3, 4, 5, 6, 9,
10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 and 33.
Mixed-code inputs need correct per-CDS values
or another deliberate strategy. Taxonomy lookup is not inferred from a species
name: Helixer's name-only directive is omitted with a warning if no taxon is
provided, while an explicit taxon produces a canonical URI.

For automatic/manual exclusion examples, classification rules, offline operation
and tag reservation, see [ORGANELLES.md](ORGANELLES.md). The assembly and
manifest PRIMARY_ID remain unchanged.

## Locus tags

| Option | Default / effect |
| --- | --- |
| `--locus-tag-prefix PREFIX` | Required; `[A-Z][A-Z0-9]{2,11}`: 3–12 uppercase alphanumeric characters, beginning with a letter. Registration is not checked. |
| `--locus-tag-map PATH` | Existing cumulative TSV map for preserving assignments on updates. |
| `--width INTEGER` | `6`; minimum numeric suffix width. |
| `--start INTEGER` | `1`; lower starting point for new allocation, raised above existing numeric tags. |
| `--step INTEGER` | `1`; increment between newly generated tags. |

Width, start and step must be positive integers.

A gene and its supported components share one `PREFIX_suffix` tag across all
isoforms; `repeat_region` is excluded. Allocation spans the whole annotation.
Keep the cumulative map permanently, including reserved tags. Stable matching
uses gene IDs, not coordinates: renamed, split or merged genes need curated
reconciliation. See [LOCUS_TAGS.md](LOCUS_TAGS.md) for accepted existing suffixes,
registration and update examples.

## Output package and naming

The default directory is `prepared_annotations/LABEL/`. It is staged before
publication; an existing destination is refused. Choose a new label/directory
for a rerun. The input annotation is not overwritten.

| File | Contents |
| --- | --- |
| `LABEL.ena.gff3.gz` | Prepared annotation. |
| `manifest.txt` | Tab-separated NAME, STUDY, GFF3, ANALYSIS_TYPE, PRIMARY_ID and DESCRIPTION fields. The GFF3 path is relative to the manifest. |
| `locus_tags.tsv` | Cumulative `gene_id`, `locus_tag`, `status` map, including reserved entries. |
| `preparation.json` | Input hashes, settings, source detection, transformations, lengths, provenance, warnings and advisory checks. |
| `preparation.txt` | Human-readable preparation summary. |
| `VALIDATE.txt` | Command guidance for the separate validator. |
| `organelle_exclusions.tsv` | Present when filtering was requested; removed sequence and feature/gene counts. |

The default label is `ACCESSION.PREDICTOR[.VERSION][.MODEL]`; unavailable optional
components and their separating dots are omitted. Without a predictor, it is
`ACCESSION.INPUT_STEM`. Predictor spelling is preserved and unsafe filename
characters become underscores. A model requires predictor/version even if no
inference is being added.

NAME defaults to `LABEL decoupled annotation`. DESCRIPTION starts with either:

- `Ab initio genome annotation prepared for ENA decoupled submission.` for the
  resolved ab-initio method; or
- `Genome annotation prepared for ENA decoupled submission.` for other methods.

Available metadata is appended as `Annotator: NAME Version: VERSION Model: MODEL`.
| Option | Default / effect |
| --- | --- |
| `--predictor NAME` | Detected where supported; explicit pipeline metadata can override component names, except conflicting ANNEVO headers. |
| `--predictor-version VERSION` | Optional for direct-source auto policy outside ANNEVO; needed for model metadata and global predictor/version inference. |
| `--predictor-model MODEL` | Optional label recorded in default filename, NAME, DESCRIPTION and report; does not change inference or run a model. |
| `--label LABEL` | Override filename stem/default directory and default NAME. Must start with a letter/digit and contain only letters, digits, dots, underscores or hyphens. |
| `--out-root PATH` | `prepared_annotations`; parent of automatically named packages. |
| `--output-dir PATH` | Exact new package directory; overrides out-root/label directory placement, not the filename. |
| `--name "TEXT"` | Override manifest NAME. |
| `--description "TEXT"` | Override DESCRIPTION completely; metadata is not automatically appended to custom text. |
| `--bioproject PRJEB123456` | Required BioProject accession, written as STUDY. Accepted prefixes: PRJEB, PRJNA, PRJDB. |
| `--assembly-accession ACCESSION` | Required versioned assembly accession; written as PRIMARY_ID. |

Compressed GFF bytes are reproducible for identical input/settings because gzip
headers omit timestamps and embedded filenames. Reports contain actual paths
and can differ between directories. The header includes `##gff-version 3`,
`#!insdc-gff-version 0.5.1`, optional species metadata and sequence declarations.
The profile version is separate from the converter version.

## Configuration and batches

JSON defaults use CLI option names with underscores. Unknown keys are rejected.
Use JSON booleans for boolean settings. `remove_organelles` also accepts a
comma-separated string of exact sequence IDs for manual filtering (see
[ORGANELLES.md](ORGANELLES.md#config-batch-and-audit)). For example:

```json
{
  "translation_table": 1,
  "out_root": "../prepared_annotations",
  "sequence_cache_dir": "../ena-cache"
}
```

```bash
python3 ena_gffprep.py --config defaults.json --input annotation.gff3 \
  --assembly-accession GCA_123456789.1 --bioproject PRJEB123456 --locus-tag-prefix YOURPREFIX
```

A batch TSV needs an `input` column. Required job values can be inherited from
configuration/CLI defaults or supplied in columns. See [examples](examples/README.md)
for a table containing different annotators.

```bash
python3 ena_gffprep.py \
  --config examples/nuclear_defaults.json \
  --batch examples/annotations.tsv \
  --batch-summary batch_summary.json
```

| Rule | Behavior |
| --- | --- |
| Precedence | Config < explicit CLI defaults < nonblank batch cells. Blank cells inherit; they cannot clear a default. |
| Paths | CLI paths are relative to the working directory; config paths to the JSON file; row paths to the TSV. Batch-summary is relative to the working directory. |
| Booleans | JSON uses `true`/`false`; TSV uses the strings `true`/`false`. CLI switches enable flags; omit true defaults from a config when a single run needs them off. |
| Shared metadata | Only use shared predictor/version/model defaults when accurate for every inheriting job. |
| Safety | Repeated output directories or prefixes are refused. Prepare an assembly together; independent per-chromosome tag allocations can collide. |
| Failures | Jobs run sequentially; individual preparation failures are recorded and later jobs continue. Choose a new summary path for reruns. |

`--config PATH`, `--batch PATH` and `--batch-summary PATH` are mode options, not
keys in a preparation config. Batch-summary defaults to `batch_summary.json`.
Explicit-output options are not valid batch columns. Avoid a shared output-dir
for multiple jobs; use out-root and individual labels instead.

## Audit and explicit outputs

Audit requires only an input and a new report path (no BioProject):

```bash
python3 ena_gffprep.py --input annotation.gff3 --audit-only --report audit.json
```

It is offline and descriptive: no prepared GFF, manifest, locus tags or ENA
lookup. Filtering and ENA metadata options cannot be used. `--input-format` can restrict the parsed syntax; audit does not run the normalizer. CDS arithmetic is
grouped by parent and can conflate distinct products. The 90-nt review threshold
is not an asserted ENA rejection rule. The report's parent directory must exist.

Advanced callers can choose individual output paths instead of a package:

| Option | Requirement |
| --- | --- |
| `--output PATH` | Prepared GFF; required with report and map output. `.gz` selects compression. |
| `--report PATH` | Preparation JSON; required with output and map output. Also used alone with audit. |
| `--write-locus-tag-map PATH` | New cumulative mapping; required with output and report. Different from input locus-tag-map. |
| `--manifest PATH` | Optional manifest in explicit-output mode. |

Outside audit, any of these four options selects explicit-output mode; output,
report and map output must then all be supplied. All normal preparation
requirements still apply. Paths must be new, distinct, have existing parent
directories and not overwrite inputs. Package directory options do not relocate
explicit files. This mode cannot be combined with batch, and does not create
preparation.txt or VALIDATE.txt. Filtering adds `REPORT_STEM.organelles.tsv`.

`--help`/`-h` lists arguments; `--version` prints the converter version. Exit codes
are **0** for success, **1** for failed batch jobs and **2** for caught top-level
preparation/usage errors.

## Validation and troubleshooting

Follow [VALIDATION.md](VALIDATION.md) for Webin commands, local validation,
report-review helpers and test coverage. Preparation checks structure and
metadata; external validation and account/archive checks remain separate.

| Symptom | Action |
| --- | --- |
| Cannot determine predictor/version | Supply accurate metadata or choose explicit inference/no-inference; inspect detected format/provenance in the report. |
| Model requires predictor/version | Provide both; model names do not establish software versions. |
| Missing/ambiguous sequence or length conflict | Check the exact assembly, aliases and supplied length files; use complete offline ENA metadata if needed. |
| ENA request/cache failure | Check HTTPS/cache permissions, explicitly refresh a corrupt cache, or supply offline metadata. |
| Unexpected output name/description | Check label, name, description and explicit-output overrides. |
| Existing output | Choose a new directory, label or file path. |
| Invalid/conflicting locus tag | Check prefix syntax, stable gene IDs and the complete cumulative map. |
| Species-name warning | Supply the correct taxon ID if a canonical species directive is wanted. |
| Short-CDS or other biological validation finding | Review evidence and the prediction; do not invent inference, partial flags or exceptions to silence the rule. |

Origin crossings, remote locations, trans-splicing, circular RNA, explicit
interval ordering and some slippage representations require specialist handling.
The tool does not establish nucleotide identity, translation correctness,
biological partiality, all vocabularies, prefix ownership or archive acceptance.
Product names and protein IDs are not invented by default. See also
[LOCUS_TAGS.md](LOCUS_TAGS.md), [ORGANELLES.md](ORGANELLES.md) and [CHANGELOG.md](CHANGELOG.md).

## Code layout

| Component | Responsibility |
| --- | --- |
| `ena_gffprep.py` | CLI, configuration and batch resolution, provenance and package naming. |
| `ena_gffprep_core.py` | Parse feature graphs, validate metadata and relationships, allocate tags, write GFF3 and reports. |
| `ena_formats.py` | Detect input syntax and predictor provenance; construct GTF relationships. |
| `ena_normalize.py` | Source-independent structural transformations. |
| `ena_organelles.py` | ENA sequence reports, caching, alias matching and sequence exclusion. |
| `ena_common.py` | Shared I/O, validation and staged publication helpers. |
| `run_bundled_validator.py` | Run local rules from a separately supplied Webin JAR. |
| `ena_validation_review.py` | Associate validator diagnostics with features. |
| `ena_add_inference.py` | Optional explicit inference addition to a separate GFF3 copy. |

See [VALIDATION.md](VALIDATION.md#automated-tests-and-scope) for the regression suite.

## Authors and citation

Developed by Nathan Riley as part of the Tree of Life programme at the
Wellcome Sanger Institute. To cite this software, use the citation metadata
in [CITATION.cff](CITATION.cff) or GitHub's "Cite this repository" feature.

## License

Copyright (c) 2026 Genome Research Ltd.

This software is distributed under the [MIT License](LICENSE).

## Submission conformance

See [CONFORMANCE.md](CONFORMANCE.md) for the targeted INSDC Submission profile, local checks, unsupported representations and external validation requirements.

The development conformance checks use locally generated feature and Dbxref vocabulary indexes; see [CONFORMANCE.md](CONFORMANCE.md) for provenance and limitations. Valid specialist terms outside the snapshot require review.

### Compartment-specific genetic codes (development)

For mixed nuclear and organellar annotations, use explicit CDS `transl_table` attributes or a file-level pragma such as `#!transl_table primary:1,mitochondrion:2` with verified ENA sequence metadata (`--sequence-metadata`). The example codes are not universal. See [CONFORMANCE.md](CONFORMANCE.md) for supported compartment keys and unresolved-classification safeguards.

Specialist feature support is conservative: explicitly declared one-base ribosomal-slippage CDS overlaps and supported pseudogenic component types are handled; circular virtual coordinates require authoritative sequence topology and remain unsupported. See [CONFORMANCE.md](CONFORMANCE.md).

The development conformance checks use a pinned local Sequence Ontology graph
alongside the INSDC-oriented feature mapping. See [CONFORMANCE.md](CONFORMANCE.md)
for provenance, legacy aliases and relationship assessment limitations.

### Local reference setup (development branch)

Before preparing annotations or running the full regression suite, obtain the
SO OBO release 2026-08-07, the compatible Webin CLI JAR containing gff3tools
5.5.1, and a plain-text list of current INSDC Dbxref database tags (one per line,
case-sensitive, no identifiers). See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

```bash
python3 tools/setup_references.py --obo /path/to/so.obo \
  --jar /path/to/webin-cli.jar --dbxref-tags /path/to/dbxref-tags.txt
```

The exact OBO bytes and embedded JAR are checksum-pinned. The tag list is
user-maintained; use the same list for reproducible results. The regression
suite expects the Stage 5 reference version and tags including GeneID, GO and
ATCC(in host). Generated indexes are ignored by Git and must not be published.
The setup command performs no network access. Original software remains MIT;
external reference materials retain their upstream terms.
