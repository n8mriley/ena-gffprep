# ENA GFF preparation tools

[![Python tests](https://github.com/n8mriley/ena-gffprep/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/n8mriley/ena-gffprep/actions/workflows/tests.yml)

Prepare existing ab initio annotations for ENA decoupled GFF3 submission. This
README describes **tool version 0.6.0**. The scripts prepare local files; they do
not authenticate to ENA or submit data.

By default, the converter preserves prediction coordinates and intervals while normalising
supported GFF representations and adding submission metadata. It does not remove
short predictions to silence biological warnings. Supported feature-type
normalisation, shared-leaf duplication and ANNEVO CDS grouping can change the
representation and identifiers. Optional `--remove-organelles` excludes predictions on NCBI-classified organelle sequences before conversion. Retained prediction coordinates are preserved.

## Requirements

- Python **3.9 or newer**, using only the standard library. No pip installation,
  GPU or Python environment is required.
- Keep the scripts together: `ena_gffprep.py` imports `ena_gffprep_core.py` and
  `ena_common.py`. Organelle filtering also uses `ena_organelles.py`.
- `--remove-organelles` needs outbound HTTPS for an uncached lookup, or a complete offline NCBI sequence report. No extra Python packages or NCBI CLI are required.
- Input must be supported GFF3, plain text or gzip-compressed with a `.gz` suffix.
  A `.gff` extension is accepted if the contents are GFF3. GTF/GFF2 is not supported.
- Input text and metadata must be compatible with the scripts' ASCII checks.
- Java **17 or newer** and a separately obtained compatible Webin JAR are needed
  for validation, not for Python preparation. The JAR is not distributed here.

Check the checkout:

```bash
python3 ena_gffprep.py --version
python3 -m unittest discover -s tests -q
```

Version 0.6.0 has 120 automated tests. See [VALIDATION.md](VALIDATION.md) for their
scope and the distinction between historical integration evidence and fresh
Webin validation.

## What is required for a normal preparation?

Requirements apply to the resolved settings: values may come from the command
line, a JSON configuration file or a batch row.

| Information | Required? | How to provide it |
| --- | --- | --- |
| Input annotation | **Yes** | `--input` |
| Assembly accession with version | **Yes** | `--assembly-accession`, e.g. `GCA_948150665.1`; the CLI also accepts GCF syntax, but verify the correct ENA primary accession for your submission |
| Locus-tag prefix | **Yes** | `--locus-tag-prefix`; use the prefix registered for the annotation |
| True sequence lengths | **Yes**, but a separate file is not always required | `--fai`, `--assembly-report`, existing `##sequence-region` declarations, or recognised ANNEVO length comments |
| Prediction provenance policy | **Yes**, but explicit flags are not always required | Automatic ANNEVO header detection, both `--predictor` and `--predictor-version`, a full `--inference`, or explicit `--no-inference` |
| Predictor model | **No** | `--predictor-model`; records the actual model in default naming, manifest description and reports; requires a predictor and version |
| Organelle exclusion | **No**; off by default | `--remove-organelles`; fetch/cache exact-assembly sequence metadata or provide `--organelle-metadata` |
| Taxonomy ID | **No** | Optional `--taxon-id`; otherwise preserve a valid existing `##species`, or omit it when absent |
| Output path | **No** | Automatic package layout unless overridden |
| Translation table | **No explicit flag required** | Existing per-CDS tables are preserved; missing translating-CDS tables use fallback 1 unless overridden |
| Previous locus-tag map | **Needed to protect historical assignments on updates** | `--locus-tag-map`; not needed for the first allocation |
| Product fallback | **No** | Existing products are preserved; no new product name is invented by default |
| Webin credentials, study or sample flags | **Not required by this preparation CLI** | Account and archive metadata checks belong to the separate Webin workflow |

A successful preparation does not verify prefix registration, archived sequence
identity, development-environment availability or acceptance of the annotation.

## Quick start: ANNEVO

### Explicit predictor, version and model

This example works with an unversioned ANNEVO header when recognised
sequence-length comments are present:

```bash
python3 ena_gffprep.py \
  --input /path/to/Abraxas_sylvata.GCA_964187985.1.annevo_2.2.1.Insecta.gff \
  --assembly-accession GCA_964187985.1 \
  --locus-tag-prefix YOURPREFIX \
  --predictor ANNEVO \
  --predictor-version 2.2.1 \
  --predictor-model Insecta
```

Replace the paths, accession and prefix for your annotation. `YOURPREFIX` is a
placeholder, not a registered prefix. `Insecta` is the model used to generate
this example annotation; supply the actual model for each input, or omit the
option if you are not recording one. It is not the organism name or taxon ID.

`--predictor ANNEVO` selects ANNEVO processing in automatic source-format mode.
You do not also need `--source-format annevo`. If the input lacks usable length
comments, add an independent length index:

```bash
--fai /path/to/matching-assembly.fa.fai
```

A FAI is recommended for checking sequence names and lengths against the exact
assembly, even when comment lengths are available. The accession flag does not
download sequence metadata.

### Minimum run with versioned ANNEVO headers

If recognised input comments already identify ANNEVO, its version and sequence
lengths, only these three options need to be supplied:

```bash
python3 ena_gffprep.py \
  --input /path/to/annotation.gff \
  --assembly-accession GCA_948150665.1 \
  --locus-tag-prefix YOURPREFIX
```

Optionally add `--predictor-model Invertebrate` if that was the model actually
used. Model metadata is never inferred from a filename or guessed from taxonomy.

Versioned headers supply predictor/version automatically. Unversioned headers
identify ANNEVO but cannot provide its version: normally supply
`--predictor-version` explicitly. An accurate full inference or explicit
`--no-inference` is an alternative provenance policy, but recording a model still
requires predictor and version. Explicit metadata conflicting with a detected
version is rejected.

When neither FAI nor assembly report is supplied, ANNEVO mode attempts to read
recognised length comments automatically, including supported spacing/case
variants. Missing or unrecognised comments require an external source. Lengths
are never estimated from the last annotated feature. Comment-only lengths
produce a warning; compare them with the exact submitted assembly.

## Exclude predictions on organelle sequences

For annotations generated by a nuclear model on an assembly containing organelles,
add `--remove-organelles` to the normal preparation command:

```bash
python3 ena_gffprep.py \
  --input /path/to/Abraxas_sylvata.GCA_964187985.1.annevo_2.2.1.Insecta.gff \
  --assembly-accession GCA_964187985.1 \
  --locus-tag-prefix YOURPREFIX \
  --predictor ANNEVO --predictor-version 2.2.1 \
  --predictor-model Insecta \
  --remove-organelles \
  --output-dir /path/to/new-filtered-package
```

The tool matches exact NCBI sequence aliases for this assembly/version and removes
all feature types on explicitly identified organelles, plus their output sequence
declarations. NCBI metadata is cached; the input GFF and assembly are not edited.
The flag does not change automatic output labels, so use a new output directory
if a previous package exists. The manifest still refers to the full assembly.

The package adds `organelle_exclusions.tsv` and an `organelle_filter` section in
`preparation.json`. Inspect them to see which sequences and genes were excluded.
Unmatched sequences, ambiguous aliases and broken cross-boundary references fail.
Unclassified primary-assembly unplaced/unlocalized scaffolds are retained with a
warning, not certified nuclear. Use the latest cumulative locus-tag map on updates
so removed genes remain reserved. See [ORGANELLES.md](ORGANELLES.md) for full policy,
cache/offline handling, classification limits and reporting.

## Other annotation sources

For GFF3 with a supported gene hierarchy and correctly grouped CDS intervals:

```bash
python3 ena_gffprep.py \
  --input /path/to/genes.gff3.gz \
  --source-format generic \
  --predictor ACTUAL_TOOL \
  --predictor-version ACTUAL_VERSION \
  --predictor-model ACTUAL_MODEL \
  --fai /path/to/assembly.fa.fai \
  --assembly-accession GCA_123456789.1 \
  --locus-tag-prefix YOURPREFIX
```

Use actual provenance and omit `--predictor-model` when it is not applicable.
A model may be a named trained checkpoint, model family or other model label
used by the predictor; it does not select or run that model. Alternatively supply a complete, accurate inference value:

```bash
--inference "ab initio prediction:ACTUAL_TOOL:ACTUAL_VERSION"
```

To preserve existing evidence without adding inference, explicitly use:

```bash
--no-inference
```

This does not delete existing inference. Existing differing inference is not
silently overwritten. Adding inference records the prediction method; it does
not demonstrate expression or biological correctness.

`generic` does not repair every predictor dialect. Intervals of a single CDS
must already share an appropriate logical ID and consistent Parent metadata.
Orphan components, shared internal subtrees and ambiguous CDS products require
an appropriate adapter or prior correction.

## Naming, model metadata and output package

The default GFF filename follows:

```text
ACCESSION.PREDICTOR.PREDICTOR-VERSION.PREDICTOR-MODEL.ena.gff3.gz
```

If no model is supplied, its component and separating dot are omitted:

```text
ACCESSION.PREDICTOR.PREDICTOR-VERSION.ena.gff3.gz
```

For the explicit ANNEVO example above:

| Item | With `--predictor-model Insecta` | Without a model |
| --- | --- | --- |
| Label / package directory name | `GCA_964187985.1.ANNEVO.2.2.1.Insecta` | `GCA_964187985.1.ANNEVO.2.2.1` |
| GFF filename | `GCA_964187985.1.ANNEVO.2.2.1.Insecta.ena.gff3.gz` | `GCA_964187985.1.ANNEVO.2.2.1.ena.gff3.gz` |
| Manifest NAME | `GCA_964187985.1.ANNEVO.2.2.1.Insecta decoupled annotation` | `GCA_964187985.1.ANNEVO.2.2.1 decoupled annotation` |
| Manifest DESCRIPTION model text | Includes `Model: Insecta` in the default description | No Model field |

With the model specified, the default output is:

```text
prepared_annotations/GCA_964187985.1.ANNEVO.2.2.1.Insecta/
    GCA_964187985.1.ANNEVO.2.2.1.Insecta.ena.gff3.gz
    manifest.txt
    locus_tags.tsv
    preparation.json
    preparation.txt
    VALIDATE.txt
```

The default manifest description is:

```text
Ab initio genome annotation prepared for ENA decoupled submission. Annotator: ANNEVO Version: 2.2.1 Model: Insecta
```

Model metadata is also recorded in `preparation.json` under
`submission_metadata.predictor_model` and, for package runs,
`workflow.source_detection.model`. The CDS inference remains
`ab initio prediction:ANNEVO:2.2.1`; this option does not append the model to
inference or add a new per-feature model attribute.

### Naming overrides

| Option | Effect |
| --- | --- |
| `--label LABEL` | Replaces the automatic filename stem and default directory name. The default manifest NAME then follows that label. |
| `--output-dir PATH` | Changes the exact package directory, leaving the label-based GFF filename unchanged. |
| `--out-root PATH` | Changes the parent of automatically named package directories. |
| `--name TEXT` | Replaces the manifest NAME, including its automatic model-bearing form. |
| `--description TEXT` | Replaces the entire generated description; no model text is appended. |
| `--output PATH` | In advanced explicit-output mode, uses the exact requested GFF path instead of automatic naming. |

Avoid filename/name overrides if you want the standard pattern. Predictor case
is preserved, so `ANNEVO` and `annevo` produce different label spelling. Unsafe
filename characters, such as spaces or slashes in a model, become underscores;
the original model text remains in description/report metadata. Use concise
filename-safe model names for predictable labels. The model must be nonempty
printable ASCII when provided; a blank setting is treated as omitted.

Without predictor metadata (for example, generic `--no-inference` preparation),
the fallback label is accession + input filename stem. That fallback cannot be
combined with a model unless predictor and version are also supplied.

### Output files

| File | Purpose |
| --- | --- |
| `*.ena.gff3.gz` | Prepared annotation |
| `manifest.txt` | Annotation manifest with a relative path to the GFF3 |
| `locus_tags.tsv` | Cumulative gene-to-tag mapping, including reserved historical entries; keep permanently |
| `preparation.json` | Detailed checks, metadata, hashes and warnings |
| `preparation.txt` | Short preparation summary |
| `VALIDATE.txt` | Example validation command for the generated package |
| `organelle_exclusions.tsv` | Added only with organelle removal; original sequence ID, matched NCBI sequence, compartment, removed rows and genes |

Existing output directories are refused. Package files are staged privately and
published after successful preparation. Changing a model label does not assign
new biological gene identities: use the appropriate cumulative locus-tag map
when preparing an annotation update.

## Complete option reference

The tables below describe **normal package preparation** unless another mode is
specified. Defaults are resolved by the workflow, so `--help` may not display all
of them. Paths supplied directly on the CLI resolve from your working directory.

### Input, sequence lengths and source format

| Option | Required? / default | Behaviour |
| --- | --- | --- |
| `--input PATH` | **Required** for preparation and audit; each batch job needs an input | Plain GFF3 or `.gz` input. No extension conversion is needed for valid GFF3 stored as `.gff`. |
| `--source-format FORMAT` | Optional; `auto` | `auto` uses the ANNEVO adapter when the supplied predictor is ANNEVO (case-insensitive) or recognised headers identify it; otherwise it uses `generic`. Explicit format choices override this selection. `annevo` applies the numbered-CDS adapter; `generic` expects suitable grouping already. Setting `annevo` alone does not supply a predictor version. |
| `--fai PATH` | Optional file, but true lengths are required | Reads sequence IDs and lengths from an existing index; does not read nucleotide sequence. Checks names, duplicate IDs and length consistency. |
| `--assembly-report PATH` | Optional alternative or supplement to FAI | Reads a standard NCBI assembly report for lengths and alias-to-accession mapping. Use the exact assembly version. Header accession conflicts are rejected when that metadata is present. |
| `--accession-column COLUMN` | Optional; `genbank` | Assembly-report target accessions: use `genbank` or `refseq`. It has no useful effect without `--assembly-report`. |
| `--use-annevo-lengths` | Optional; automatically enabled for ANNEVO without FAI/report | Reads recognised ANNEVO sequence-length comments. Requires the ANNEVO adapter and matching comments. Can also be requested with an external source for consistency checks. |

A generic input may use existing `##sequence-region` declarations without an
external length file. The ANNEVO automatic-comment path still requires recognised
length comments if no FAI/report is supplied. FAI and assembly-report data must
agree when used together. Sequence declarations are checked against external
metadata rather than used to hide missing external entries. An assembly report
can remap recognised sequence aliases; a FAI alone does not infer alias mappings.
No combination of names and lengths proves archived nucleotide identity.

### Assembly and taxonomy metadata

| Option | Required? / default | Behaviour |
| --- | --- | --- |
| `--assembly-accession ACCESSION` | **Required** for normal and explicit-output preparation | Versioned `GCA_` or `GCF_` syntax. Written as manifest `PRIMARY_ID` when a manifest is produced; does not itself translate sequence names in column 1. |
| `--taxon-id INTEGER` | Optional; absent | Positive taxonomy ID. Adds a canonical `##species` URI and rejects disagreement with existing input taxonomy. |
| `--name TEXT` | Optional; label + ` decoupled annotation` | Manifest `NAME`. Quote values containing spaces. Does not rename genes or sequences. |
| `--description TEXT` | Optional; generated Annotator / Version / Model description | Manifest `DESCRIPTION`. Quote values containing spaces. Does not add feature products. |

Without `--taxon-id`, a recognised existing `##species` directive is preserved;
if there is none, no species directive is added. Conflicting input taxonomies are
rejected. To omit a taxon supplied by a config, remove it from that config; a
blank batch cell inherits defaults rather than clearing them. Omitting the flag
does not strip an existing species directive.

### Inference, products and genetic code

| Option | Required? / default | Behaviour |
| --- | --- | --- |
| `--predictor NAME` | Conditional | ANNEVO selects the ANNEVO adapter in auto mode, case-insensitively. With `--predictor-version`, constructs `ab initio prediction:NAME:VERSION`. Optional when detected from ANNEVO or when using full inference / explicit no-inference. |
| `--predictor-model MODEL` | Optional; absent | Records the specific prediction model. Included after predictor/version in the default GFF filename, package label and manifest NAME, and recorded as Model in the default DESCRIPTION. An explicit --description is preserved unchanged. Does not change inference. Requires predictor and version, explicit or detected; no model is guessed. Supported in JSON as `predictor_model` and as a batch column. |
| `--predictor-version VERSION` | Conditional | Version of the actual predictor; do not use the converter version. Both predictor fields are needed to construct inference for an unrecognised source. |
| `--inference TEXT` | Optional alternative to predictor pair | Complete nonempty printable ASCII inference value, quoted as one shell argument. Applied to CDS rows; conflicting existing values are rejected. Detected ANNEVO provenance must agree. |
| `--no-inference` | Optional; off | Explicitly disables adding new inference while retaining existing evidence. Cannot be combined with `--inference`. |
| `--translation-table INTEGER` | Optional; `1` | Fallback for translating CDSs missing a table. Existing valid per-CDS tables remain unchanged. Pseudogenic CDSs omit translation tables. |
| `--default-product TEXT` | Optional; absent | Explicit fallback for CDSs missing a product; existing products remain unchanged. No product is invented when omitted; missing products are noted. |

Table 1 is a **nuclear assumption**, not a universal choice. The tool does not
infer organelle codes from names or taxonomy. Mixed-code annotations need correct
existing per-CDS values or another deliberate preparation strategy. An explicit
fallback does not override existing tables.

### Locus tags

| Option | Required? / default | Behaviour |
| --- | --- | --- |
| `--locus-tag-prefix PREFIX` | **Required** | 3–12 uppercase alphanumeric characters, starting with a letter: `[A-Z][A-Z0-9]{2,11}`. The tool checks syntax, not registration or ownership. |
| `--locus-tag-map PATH` | Optional for first allocation; important for updates | Existing TSV containing `gene_id` and `locus_tag`. Preserves exact assignments and reserves entries absent from the current annotation. |
| `--width INTEGER` | Optional; `6` | Minimum numeric suffix width. `--width 5` gives tags such as `WSITEST_00001`. Numbers expand beyond this width rather than wrap. Must be positive. |
| `--start INTEGER` | Optional; `1` | Lower starting point for new numeric allocation. Supplied tags can raise the allocation above their numeric high-water mark. Must be positive. |
| `--step INTEGER` | Optional; `1` | Positive increment between newly allocated numeric tags. Does not renumber existing assignments. |

New tags use `PREFIX_number`, with the underscore separator. Existing ASCII
alphanumeric suffixes are accepted and preserved. Decimal-version suffixes are
rejected. Differently padded versions of the same numeric suffix remain rejected
as ambiguous by this tool.

Allocation is shared across the annotation's sequences and gene types. A gene's
supported components share its tag; repeat_region is excluded. Historical numeric
tags contribute to the high-water mark; alphanumeric suffixes are reserved but
do not contribute a numeric maximum. Unused gaps are not automatically filled.

For updates, supply the **complete latest cumulative map**, including reserved
entries. Stable matching is by source gene ID, not coordinates or similarity;
split, merged or renamed genes need curated reconciliation. See
[LOCUS_TAGS.md](LOCUS_TAGS.md) for examples and registration guidance.

### Package naming and warning policy

| Option | Required? / default | Behaviour |
| --- | --- | --- |
| `--label LABEL` | Optional; derived | Output GFF filename stem and default package directory name. Start with a letter/digit; use only letters, digits, dots, underscores and hyphens. |
| `--out-root PATH` | Optional; `prepared_annotations` | Parent directory for automatically named packages. |
| `--output-dir PATH` | Optional | Exact new package directory; overrides the directory derived from output root + label. The label still controls the GFF filename. |
| `--remove-organelles` | Optional; off | Remove all feature rows on NCBI-classified organelle sequences before conversion and tag allocation. Requires exact-assembly metadata, fetched automatically or supplied offline. Does not change translation tables to repair organelle predictions. |
| `--organelle-metadata PATH` | Optional with removal | Complete NCBI sequence JSON/JSONL or tool cache file for offline use. Requires `--remove-organelles`; cannot combine with refresh. |
| `--organelle-cache-dir PATH` | Optional with removal | Override the default `$XDG_CACHE_HOME/ena-gffprep/ncbi-sequences` cache (`~/.cache` if unset). Requires `--remove-organelles`. |
| `--refresh-organelle-cache` | Optional; off | Fetch fresh metadata, replacing the validated cached report. Requires `--remove-organelles`; fails rather than silently reusing old data on a failed refresh. |
| `--strict-warnings` | Optional; off | Treat findings in the preparation report's warnings list as fatal. Notes are distinct from warnings; this does not guarantee a warning-free Webin run. |

There is no overwrite flag. Choose a new directory or label for each preparation.

### Configuration, batch, audit and information

| Option | Required? / default | Behaviour |
| --- | --- | --- |
| `--config PATH` | Optional | JSON defaults using option names with underscores; unknown keys are rejected. Preparation paths in it resolve relative to the JSON file. |
| `--batch PATH` | Optional | Run jobs from a tab-separated table; cannot combine with audit or explicit-output options. Details below. |
| `--batch-summary PATH` | Optional; `batch_summary.json` | JSON summary for batch mode only. Its path resolves from the working directory, not the batch file. Existing summaries are refused. |
| `--audit-only` | Optional | Descriptive inspection requiring only input and a new report path; does not prepare a submission package. |
| `--help`, `-h` | Optional | Display the CLI argument list and exit. |
| `--version` | Optional | Display the converter version and exit. |

### Explicit-output mode (advanced)

These options are for scripts that need individual output paths. For most users,
`--output-dir` is simpler.

| Option | Requirement | Behaviour |
| --- | --- | --- |
| `--output PATH` | Required together with report and map output | Prepared GFF3; `.gz` suffix selects gzip compression. |
| `--report PATH` | Required together with output and map output; also used by audit | Preparation JSON, or descriptive audit JSON in audit mode. |
| `--write-locus-tag-map PATH` | Required together with output and report | New cumulative map; distinct from the input `--locus-tag-map`. |
| `--manifest PATH` | Optional in explicit-output mode | Write the annotation manifest. Supplying this alone also selects explicit-output mode and therefore requires the other three output paths. |

Outside audit mode, any of these four settings selects explicit-output mode.
All normal preparation requirements still apply, including assembly accession.
Output paths must be distinct, new, and must not overwrite inputs. Their parent
directories must already exist. Package naming options do not relocate these
explicit paths. This mode does not create `preparation.txt` or `VALIDATE.txt`.
It cannot be combined with batch mode. With organelle removal it also writes `REPORT_STEM.organelles.tsv` beside the JSON report.

```bash
mkdir -p /path/to/new-output
python3 ena_gffprep.py \
  --input /path/to/annotation.gff \
  --fai /path/to/assembly.fa.fai \
  --assembly-accession GCA_948150665.1 \
  --locus-tag-prefix YOURPREFIX \
  --predictor ACTUAL_TOOL --predictor-version ACTUAL_VERSION \
  --predictor-model ACTUAL_MODEL \
  --output /path/to/new-output/annotation.ena.gff3.gz \
  --report /path/to/new-output/preparation.json \
  --write-locus-tag-map /path/to/new-output/locus_tags.tsv \
  --manifest /path/to/new-output/manifest.txt
```

## Shared configuration and batch runs

Example JSON defaults:

```json
{
  "translation_table": 1,
  "out_root": "../prepared_annotations"
}
```

Use JSON booleans (`true`/`false`), not quoted boolean strings. Do not put mode
switches such as `batch` or `audit_only` in this preparation-settings object.
CLI flags can enable boolean settings; use a separate config when you need to
unset a true config boolean for a single run.

A batch file must have an `input` column. Each job also needs all required values,
which may be inherited from config or CLI defaults. A model can be set in JSON:

```json
{
  "predictor": "ANNEVO",
  "predictor_version": "2.2.1",
  "predictor_model": "Insecta",
  "out_root": "../prepared_annotations"
}
```

Use shared predictor/version/model defaults only when they accurately describe
all jobs inheriting them. For mixed annotations, use per-row columns. This
Python example creates a table with real tab separators; edit the placeholder
paths and prefixes before running preparation:

```bash
python3 - <<'PYTHON'
import csv
with open("annotations.tsv", "w", newline="") as handle:
    writer = csv.writer(handle, delimiter="\t")
    writer.writerow(["input", "assembly_accession", "locus_tag_prefix",
                     "predictor", "predictor_version", "predictor_model"])
    writer.writerow(["/path/to/abraxas.gff", "GCA_964187985.1", "YOURPREFIX1",
                     "ANNEVO", "2.2.1", "Insecta"])
    writer.writerow(["/path/to/other.gff", "GCA_123456789.1", "YOURPREFIX2",
                     "ANNEVO", "2.2.3", ""])
PYTHON
```

The second row has no model only if no model is supplied by config/CLI defaults.
A blank cell inherits rather than clears a default. For mixed model/no-model
runs, omit a shared `predictor_model` and populate only the applicable row cells.

Paths in nonblank row cells resolve relative to the TSV. Optional columns include
`fai`, `assembly_report`, `taxon_id`, `locus_tag_map` and `label`. Column names use
underscores, not CLI hyphens. Explicit-output options are not valid batch columns.

```bash
python3 ena_gffprep.py \
  --config examples/nuclear_defaults.json \
  --batch examples/annotations.tsv \
  --batch-summary my_batch_summary.json
```

Precedence is **config < explicit CLI defaults < nonblank batch cells**. Blank
cells inherit; they do not erase a supplied value. Boolean TSV cells use `true`
or `false`. Avoid a shared `output_dir` for multiple jobs; use an output root and
individual labels or distinct per-row directories.

Jobs run sequentially. Individual preparation failures are recorded and later
jobs continue. Repeated output directories or locus-tag prefixes are refused:
independent jobs must not allocate colliding tags. Prepare all chromosomes of one
assembly together, or coordinate shared-prefix updates deliberately.

Exit statuses: **0** for success, **1** for one or more failed batch jobs,
**2** for caught top-level preparation/usage errors. Inspect the batch summary
for each job's outcome and use a fresh summary path when rerunning.

For batch filtering, add `remove_organelles` (`true`/`false`),
`organelle_metadata`, `organelle_cache_dir` and/or `refresh_organelle_cache`
columns as needed. Config uses real JSON booleans. Offline metadata must match
each row's assembly. Supplied paths follow the normal config/TSV rules.

## Audit without submission metadata

```bash
python3 ena_gffprep.py \
  --input /path/to/annotation.gff3.gz \
  --audit-only \
  --report source.audit.json
```

Only input and report are required. Organelle filtering is unavailable in audit mode. No accession, prefix, taxon, predictor or
external length file is required for this descriptive mode. The report path must
be new and its parent directory must exist. It does not produce a prepared GFF,
manifest or locus-tag map. CDS counts/arithmetic grouped by parent are descriptive
and can conflate multiple products under one parent. The 90-nt threshold is a
review threshold, not an asserted ENA rejection rule.

## Validate the prepared package

Change into the generated package directory. Check `preparation.txt`, the JSON
warnings and `VALIDATE.txt` before running the separately obtained Webin build.

With credentials already set in your environment:

```bash
java -Xmx4g -XX:ActiveProcessorCount=2 \
  -jar /absolute/path/to/webin-cli-10.0.0-beta.jar \
  -context annotation \
  -manifest manifest.txt \
  -username "$WEBIN_USERNAME" \
  -passwordEnv WEBIN_PASSWORD \
  -validate
```

`-passwordEnv WEBIN_PASSWORD` names the environment variable; it does not take
`"$WEBIN_PASSWORD"` as its value. The password variable must be exported so Java
can read it. Keep credentials out of the repository and manifests.

Follow your ENA contact's instructions for the intended environment and JAR;
the supplied pilot command omits `-test`. Successful preparation does not check
archive routing. A successful validation can still contain warnings; inspect the
detailed report. No Python command in this repository submits data.

## Additional GFF helpers

Run from the toolkit directory and choose new output paths:

```bash
# Associate supported diagnostic lines with the exact validated GFF:
python3 ena_validation_review.py \
  --gff /path/to/prepared.ena.gff3.gz \
  --webin-report /path/to/exact-file.report \
  --outdir webin_review

# Local bundled gff3tools check without account/registry checks:
python3 run_bundled_validator.py \
  --jar /path/to/webin-cli-10.0.0-beta.jar \
  --gff /path/to/prepared.ena.gff3.gz \
  --log local-validation.log
```

The report mapper needs the exact GFF that produced the log, including line
order. It recognises supported diagnostic syntax rather than every possible log
format. The bundled check is not a full authenticated Webin validation.

`ena_add_inference.py` is a separate helper for explicit provenance additions to
already-prepared files. Normal preparation already handles inference; use the
helper's `--help` before choosing a separate workflow.

## Common mistakes

| Symptom | What to check |
| --- | --- |
| NCBI metadata lookup failed | Check HTTPS access, use a validated cache, or supply a complete offline sequence report. Filtering does not silently proceed without metadata. |
| Unmatched/ambiguous GFF sequence | Check exact assembly version and sequence aliases; accession versions are not stripped. |
| Cannot determine predictor/version | Supply both predictor flags, a full inference, or explicit no-inference; confirm ANNEVO provenance comments are present. |
| Model requires predictor/version | Supply both fields explicitly or use a recognised versioned header. A model does not supply the predictor version. |
| Model missing from a custom filename or NAME | `--label`, `--name` and explicit `--output` override automatic naming; remove them or include the desired components yourself. |
| No lengths / missing true sequence length | Supply the correct FAI/report, valid declarations for generic input, or recognised ANNEVO comments. Do not substitute maximum feature coordinates for sequence lengths. |
| Missing parent / orphan component | Input needs a supported gene hierarchy; generic mode cannot invent one safely. |
| Existing output directory/file | Choose a new label, directory or explicit file path. |
| Invalid/conflicting locus tag | Check registered prefix syntax and use the latest compatible cumulative map. |
| Legacy mode requires three outputs | Use `--output-dir` for a package; otherwise supply output, report and map output together. Audit mode is the exception for standalone `--report`. |
| Webin cannot read the GFF file | Check the manifest path, actual filename and working directory. Generated package manifests point to a sibling GFF file. |
| Webin cannot read password environment variable | Export the variable and pass its name to `-passwordEnv`, not the expanded password. |

## Scope and limitations

Supported preparation focuses on ordinary explicit gene hierarchies, known
sequence bounds, safe multipart IDs and single-parent output. Shared leaf features
may be duplicated per parent; shared internal subtrees are refused. Repeated CDS
IDs identify intervals of one discontinuous feature, not necessarily an error.

Origin crossings, remote locations, trans-splicing, explicit interval ordering,
circular RNA and some ribosomal-slippage representations remain unsupported or
refused. Embedded sequence sections are rejected. The tool does not establish
translation correctness, biological partiality, all controlled vocabularies or
archive acceptance. No product names or protein IDs are invented by default.

The GFF header includes `##gff-version 3`, `#!insdc-gff-version 0.5.1`, optional
species metadata and sequence-region declarations. The profile version in the
header is separate from the converter's software version.

The preparation report uses `feature_types_in_supported_subset` for its local
feature vocabulary summary. The corresponding internal set is `SUPPORTED_TYPES`.

See [ORGANELLES.md](ORGANELLES.md), [LOCUS_TAGS.md](LOCUS_TAGS.md), [VALIDATION.md](VALIDATION.md) and
[CHANGELOG.md](CHANGELOG.md) for further details. Start with `ena_gffprep.py` to
follow the CLI workflow, `ena_gffprep_core.py` for transformations and
`ena_common.py` for shared I/O and package publication.
