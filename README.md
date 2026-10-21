# ENA GFF preparation tools 0.5.5
\n[![Python tests](https://github.com/n8mriley/ena-gffprep/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/n8mriley/ena-gffprep/actions/workflows/tests.yml)


Prepare existing ab initio annotations for ENA’s decoupled GFF3 pilot.
The tools preserve predictions and do not submit anything or authenticate to ENA.

`ena_gffprep.py` is the main CLI; keep it alongside `ena_gffprep_core.py`
and `ena_common.py`. Python 3.9+ and the standard library are sufficient.
Java 17+ is needed only for validation with the supplied Webin build.

This distribution focuses on GFF preparation and GFF-based report review.
It contains no FASTA preparation, header conversion or sequence QC tools.
FAI index support remains available to check sequence identifiers and lengths
and generate `##sequence-region` declarations; nucleotide sequence is not read.

## Quick start: an ANNEVO annotation

```bash
python3 ena_gffprep.py \
  --input /path/to/ilConVacc3.annevo_2.2.3.Invertebrate.gff \
  --fai /path/to/ilConVacc3.fa.fai \
  --assembly-accession GCA_948150665.1 \
  --locus-tag-prefix YOURPREFIX
```

Replace YOURPREFIX with the prefix registered for your annotation. This example
uses the supplied ilConVacc3 metadata; other assemblies need their own metadata.
TEST is only a local test namespace, not a real registered prefix.

The tool detects ANNEVO v2.2.3 from the **file header**, selects its adapter and
adds `inference=ab initio prediction:ANNEVO:2.2.3` to every CDS interval. It does
not infer versions from filenames. A conflicting explicit predictor/version or
ANNEVO inference is rejected.

Default output:

```
prepared_annotations/GCA_948150665.1.annevo_2.2.3/
    GCA_948150665.1.annevo_2.2.3.ena.gff3.gz
    manifest.txt
    locus_tags.tsv
    preparation.json
    preparation.txt
    VALIDATE.txt
```

`--taxon-id` is optional. When supplied, the tool writes a matching
`##species` directive. When omitted, no directive is added unless the input already
contains a valid one, which is preserved. Explicit/header conflicts, conflicting
input directives and nonpositive supplied IDs are still rejected. The report
records whether taxonomy was explicit, preserved from the header or not provided.
If a config supplies `taxon_id`, remove that config entry to omit it; blank batch
cells inherit defaults. A header is not silently stripped just because a CLI
argument is absent.

The manifest points to the GFF3 in the same directory, ready for validation.
`preparation.txt` is the short summary; JSON contains detailed checks and provenance.
Existing output directories are refused. Failed preparations are not published.
Keep the latest `locus_tags.tsv` permanently.

### Defaults and overrides

| Setting | Behaviour |
| --- | --- |
| Output root | `prepared_annotations/` under the current working directory |
| Label | Assembly accession + detected predictor and version; generic inputs use accession + input stem |
| `--label` | Choose a safe filename/directory label for a distinct annotation build |
| `--output-dir` | Set the exact new package directory instead of the default |
| `--out-root` | Choose a common parent for automatically named packages |
| Genetic code | Fallback **1**, a nuclear assumption, explicitly recorded in the report |
| `--translation-table` | Override the fallback; existing per-CDS tables remain unchanged |
| Length source | Supplied FAI, assembly report or sequence-region declarations |
| ANNEVO without external lengths | Automatically use recognised length comments, with a warning |
| Biological QC | Advisory; no short gene/exon/intron is deleted or altered |
| `--strict-warnings` | Optional: make warning findings fatal |

Do not apply code 1 indiscriminately to organelles. The tool cannot infer the
right genetic code from sequence names. For mixed inputs, organelle CDSs need
their correct existing code or separate preparation. Pseudogenic CDSs omit tables.

A FAI is recommended as an independent length/membership check, but is not always
required. ANNEVO comment lengths are convenient but not an independent assembly
check. An NCBI assembly report can replace `--fai` via `--assembly-report` and map
standard aliases to accessions. Use the report for the exact assembly version.
No length source proves nucleotide identity with the archived reference.

### Other GFF3 predictors

For a standard GFF3 with correct multipart CDS grouping:

```bash
python3 ena_gffprep.py \
  --input /path/to/genes.gff3.gz --source-format generic \
  --predictor ACTUAL_TOOL --predictor-version ACTUAL_VERSION \
  --fai /path/to/genome.fa.fai \
  --assembly-accession GCA_123456789.1 \
  --locus-tag-prefix YOURPREFIX
```

Use the real values, or supply an accurate full `--inference` string. To explicitly
preserve existing evidence without adding an inference, use `--no-inference`.
It does not remove existing attributes. Differing existing inference is never
silently overwritten; it needs a deliberate merge decision.

`auto` detects recognised ANNEVO headers and otherwise selects the generic adapter.
Generic is **not** a universal dialect converter: one protein's CDS intervals must
already share one ID and have consistent Parent relationships. Missing gene
parents, ambiguous numbered CDS products, GTF/GFF2 and specialised biology require
another adapter. No tool/version is fabricated for unrecognised sources.

## Many annotations: shared config and a TSV

A config is a JSON object of CLI setting names with underscores, e.g.
`translation_table`, `out_root`, `predictor_version`. Only documented preparation
settings are accepted. Paths in a config resolve relative to that config file.

`examples/nuclear_defaults.json` is a small starting point:

```json
{
  "translation_table": 1,
  "out_root": "../prepared_annotations"
}
```

Edit `examples/annotations.tsv` with actual input paths, assembly accessions,
optional taxa and registered prefixes, then run:

```bash
python3 ena_gffprep.py \
  --config examples/nuclear_defaults.json \
  --batch examples/annotations.tsv \
  --batch-summary my_batch_summary.json
```

TSV column names are CLI option names with underscores. Only `input` must be a
column; all required values must ultimately come from config, CLI or a row.
Precedence is **config < explicit CLI defaults < nonblank row values**. Row paths
resolve relative to the TSV; relative CLI paths resolve from your current directory.
Boolean TSV values must be `true` or `false`. Empty cells inherit defaults.

Each row gets its own package. A failed row is reported and the others continue.
The process exits 1 if any row failed, 0 if all succeeded, or 2 for top-level usage
errors. Read the batch JSON for individual results. Use a new summary path on repeat
runs. All rows are processed sequentially, with no scheduler or network dependency.

The initial batch implementation refuses repeated locus-tag prefixes in one batch:
independent runs under the same prefix could allocate duplicate identifiers. Use
distinct registered prefixes or manage shared-prefix allocation as coordinated
updates separately. Do not generate independent maps with the same prefix and
assume they are globally unique.

## Updating an annotation

Add the latest cumulative map to your next preparation:

    --locus-tag-map /path/to/previous-package/locus_tags.tsv

Choose a new `--label` or `--output-dir`. Existing source gene IDs retain tags;
new genes get unused numbers and retired tags remain reserved. Do not reuse a TEST
map for a different registered namespace. Stable mapping follows **gene IDs**;
renumbered, split or merged models need identity reconciliation first.

The older explicit output-path options (`--output`, `--write-locus-tag-map`,
`--report`, optional `--manifest`) remain available for scripting. Supply the first
three together. They bypass the automatic package layout; use normal metadata and
source options too. These legacy path options are not accepted in batch rows.

## Validation and submission

Enter a generated package directory. `VALIDATE.txt` contains the command for its
context. Set credentials in your environment using your normal account procedure.

```bash
java -Xmx4g -XX:ActiveProcessorCount=2 \
  -jar /absolute/path/to/webin-cli-10.0.0-beta.jar \
  -context annotation -manifest manifest.txt \
  -username "$WEBIN_USERNAME" -passwordEnv WEBIN_PASSWORD \
  -validate
```

These commands follow
the supplied pilot instructions by omitting `-test`. This particular custom JAR
uses development submission URLs on both branches but different paths/account
routes; follow your ENA contact's routing instructions rather than assuming
ordinary release defaults.

Success ends with `Submission(s) validated successfully.` Read the detailed
report as well: warning-only validation can succeed. Check development availability
of the exact sequence accessions and use registered metadata before submitting.
If you choose to submit, replace `-validate` with `-submit` in the appropriate
command. Under the supplied pilot instructions, share the returned analysis ID
with ENA so they can manually start processing. No toolkit command performs that
submission or sends messages.

## Optional advisory QC

Run these from the toolkit directory; choose new output report paths/directories.

```bash
# Source inspection without submission metadata:
python3 ena_gffprep.py --input /path/to/genes.gff \
  --audit-only --report source.audit.json

# Map an exact Webin report to genes/locus tags:
python3 ena_validation_review.py --gff /path/to/prepared.ena.gff3.gz \
  --webin-report /path/to/exact-file.report --outdir webin_review

# Bundled validator, without account/registry checks:
python3 run_bundled_validator.py --jar /path/to/webin-cli-10.0.0-beta.jar \
  --gff /path/to/prepared.ena.gff3.gz --log local-validation.log
```

The report mapper needs the exact validated GFF3, with matching line order.
The 90-nt review threshold is descriptive, not an ENA rejection threshold.
Audit-only can conflate multiple CDS products under one Parent; use it for the
supported one-CDS-per-transcript cases. This distribution does not perform
sequence-based translation or splice-motif QC.

`ena_add_inference.py` remains a separate helper for already-prepared files;
normal runs now add provenance automatically/explicitly in the main preparer.

## Limits and testing

GFF preparation retains the tested scope: ordinary explicit gene hierarchies,
known bounds, single-parent representation and safe multipart IDs. It does not
provide a complete biological validator or all SO/qualifier vocabularies. Some
source gene feature types/pseudogene attributes are normalised as documented in
the engine; genes, coordinates and intervals are not biologically curated.

Origin crossings, remote locations, trans-splicing, explicit interval ordering,
circular RNA, shared internal subtrees and some ribosomal-slippage representations
remain unsupported or refused. Optional embedded peptide FASTA is also refused.
No prefix registration, sequence availability, study ownership or archived base
identity is verified. No product names or protein IDs are invented by default.

```bash
python3 -m unittest discover -s tests -v
```

Keep the entire directory together; scripts import the shared modules. See
VALIDATION.md for release checks, the exact supplied binary and the distinction
between local test results and full submission acceptance.

## Reading the code

Start with `ena_gffprep.py` for annotation CLI workflows. `ena_gffprep_core.py` implements the GFF3 transformations; `ena_common.py` holds shared I/O and publication helpers. Module docstrings explain scope, and function docstrings describe inputs, results and limitations. QC helpers produce advisory reports and do not submit data.

See [LOCUS_TAGS.md](LOCUS_TAGS.md) for registration, syntax, suffix conventions, feature propagation and safe annotation updates.

## Report compatibility in 0.5.5

The feature vocabulary summary is named `feature_types_in_supported_subset` in
`preparation.json`. Update downstream report readers to use this key. The internal
feature set is now `SUPPORTED_TYPES`. Neither change alters feature selection or
GFF output. Product and exception diagnostic wording has also been clarified.
