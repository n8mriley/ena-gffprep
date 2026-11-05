# Validation

Preparation and ENA acceptance are separate steps. Review `preparation.txt`,
`preparation.json` and any `organelle_exclusions.tsv`, then validate the generated
manifest with the intended Webin build. A successful preparation is not a
submission receipt or a biological quality assessment.

## Webin validation

Webin CLI 11.0.0 requires `STUDY` in the annotation manifest. Pass `--bioproject PRJEB123456` during preparation to populate `STUDY` (use your real accession). Older packages without `STUDY` must be regenerated with this option before submission. This option changes only manifest metadata, not GFF3 coordinates or predictions.


Change into the prepared package directory. `VALIDATE.txt` contains the command
pattern. With credentials already exported in the environment:

```bash
java -Xmx4g -XX:ActiveProcessorCount=2 \
  -jar /absolute/path/to/webin-cli-dev-11.0.0.jar \
  -context annotation -manifest manifest.txt \
  -username "$WEBIN_USERNAME" -passwordEnv WEBIN_PASSWORD \
  -validate
```

`-passwordEnv WEBIN_PASSWORD` takes the variable's **name**, not its expanded
value. The variable must be exported for Java to read it. Keep credentials out
of manifests and source files. The generated GFF3 path is relative to its manifest;
check that the referenced file exists if Webin cannot read it.

Use the JAR and environment instructions supplied by ENA. This example follows
the supplied development workflow and omits `-test`; do not assume that every
build routes requests identically. No Python command here submits data.

Read the detailed reports even after a success message. Warnings may remain.
Short CDSs and other biological findings can require accurate evidence or review;
do not invent products, inference, partial flags or exceptions to suppress them.

If submission later fails, inspect the complete `webin-cli.report` and server
receipt. Account acknowledgements or archive errors can appear under a generic
authentication/password summary. Local GFF validation does not check every
account or submission-service requirement.

## Local validator and report review

The optional bundled check extracts gff3tools from your separately supplied JAR.
It runs local rules without authentication or archive/registry checks:

```bash
python3 run_bundled_validator.py \
  --jar /path/to/webin-cli-dev-11.0.0.jar \
  --gff /path/to/prepared.ena.gff3.gz \
  --log local-validation.log
```

To associate supported diagnostic lines with features:

```bash
python3 ena_validation_review.py \
  --gff /path/to/prepared.ena.gff3.gz \
  --webin-report /path/to/exact-file.report \
  --outdir webin_review
```

Use the exact GFF that generated the log, including line order. The helper
recognizes supported diagnostic syntax, not every possible Java/Webin log format.
It writes finding and coordinate-based review tables plus a summary. Choose new
log/output paths. A report-review failure can leave a partial review directory;
the prepared GFF is unchanged.

`ena_add_inference.py` is an optional helper for adding explicitly supplied
provenance to a separate copy of an already prepared file. Normal preparation
already handles inference. Consult `--help`; the helper refuses existing CDS
inference rather than merging it automatically:

```bash
python3 ena_add_inference.py \
  --input prepared.ena.gff3.gz --output prepared.with-inference.gff3.gz \
  --inference "ab initio prediction:ACTUAL_TOOL:ACTUAL_VERSION"
```

Replace the provenance with the true method/version; this is an optional separate
workflow, not an additional required preparation step.

## Automated tests and scope

```bash
python3 -m unittest discover -s tests -v
```

Tests use temporary files, synthetic annotations and mocked network responses.
They need no NCBI connection, Webin credentials, Java or extra Python packages.
`release_verification.json` records the current test result.

| Area | Checks |
| --- | --- |
| Input formats | Parsing, escaping, GTF relationships, source-independent normalization, syntax overrides, explicit grouping policies, distinct products, isoforms, phase normalization and output compatibility. |
| Provenance/metadata | Direct-source and explicit inference policies, conflicts, model naming, descriptions, optional taxonomy and genetic-code preservation. |
| Structure | Parent graphs, single-parent output, multipart IDs, shared leaves, coordinate bounds and advisory CDS arithmetic. |
| Locus tags | Stable allocation, collisions, suffix syntax, propagation and historical reservations. |
| NCBI metadata | Exact assemblies, aliases, local-length precedence, bounds, offline data, pagination, caching/refresh, failure handling and shared filtering lookups. |
| Filtering/workflows | Automatic and manual organelle exclusions, optional CLI values, exact-ID validation, header-only exclusions, network bypass, retained references, tag reservation, config/batch behavior, output protection and reproducible compression. |
| Diagnostics | Supported validator-log associations and review tables. |

Local checks do not establish nucleotide identity with the archived assembly,
translation correctness, expression, biological partiality, all controlled
vocabularies, prefix ownership or archive acceptance. Metadata-only filtering
can miss unlabelled or misclassified organelles. Manual filtering excludes only
the supplied IDs and does not establish that retained sequences are nuclear. Synthetic fixtures and local
validator checks, including checks using synthetic bounds, do not certify a real
submission. Use true assembly metadata and the intended Webin workflow for each
prepared annotation.
