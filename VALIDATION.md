# Validation and test coverage

## Automated checks

Run the bundled automated tests with Python's standard-library
`unittest` runner:

```bash
python3 -m unittest discover -s tests -v
```

Tests use temporary files, synthetic annotations and mocked network responses.
They require no NCBI connection, Webin credentials, Java or external packages.
The compact result for this package is recorded in `release_verification.json`.

| Area | Coverage |
| --- | --- |
| GFF structure | Parsing/escaping, supported attributes, feature relationships, multipart IDs, shared leaves and coordinate bounds |
| Predictor metadata | ANNEVO compatibility, source/dialect detection, direct-source and mixed-source inference policy, GTF relationships, CDS grouping, non-CDS phase normalization, model naming and manifest metadata |
| Taxonomy and genetic code | Optional taxonomy, conflict checks, supported table identifiers and preservation of existing per-CDS tables |
| Locus tags | Allocation, stable assignments, collisions, alphanumeric suffixes, propagation and reserved historical tags |
| Organelle filtering | Exact-assembly metadata, aliases, excluded records, references, retained unclassified scaffolds, tag reservation and reports |
| Length fallback | Source precedence, exact-assembly matching, aliases, local/NCBI disagreements, bounds, offline data, disabled lookup and sharing with organelle removal |
| Metadata transport | JSON/JSONL input, pagination, response completeness, cache reuse/refresh and network failure handling |
| Workflows | Configuration, batch handling, audit mode, output protection and package publication |
| Diagnostics | Supported validator-log mapping and coordinate-based advisory checks |

## Validate a prepared annotation

Preparation performs local structural and metadata checks. Review
`preparation.txt`, `preparation.json` and, if filtering was enabled,
`organelle_exclusions.tsv`. Then follow `VALIDATE.txt` or the README to run the
intended Webin build against the generated manifest.

The optional `run_bundled_validator.py` helper runs the bundled local GFF
validator. It does not perform the complete authenticated Webin workflow.
Warnings may remain even when Webin reports successful validation. Review the
detailed GFF report rather than treating a success message as biological proof.

If submission fails, inspect the full `webin-cli.report` and server receipt.
Authentication, account acknowledgements and archive processing are separate
from local preparation and annotation validation.

## Limits of the checks

The test suite does not establish nucleotide identity with the archived assembly,
translation correctness, expression, biological partiality, all controlled
vocabularies, prefix registration, account permissions or final archive acceptance.
The organelle filter relies on NCBI metadata; unlabelled or misclassified sequences
can remain. See [ORGANELLES.md](ORGANELLES.md) and the README for supported scope.

A passing local suite is not a submission receipt. This package does
not claim authenticated Webin validation or submission.

## Cross-format verification

The supplied annotation examples were exercised independently; they are not
biological comparisons between genomes. ANNEVO output was compared against
0.6.1: decompressed GFF, locus-tag map and manifest matched exactly. Gzip header
metadata is now deterministic, so compressed bytes can differ from older builds.

Local validator checks use the supplied development JAR. Where true sequence
lengths were unavailable, temporary synthetic bounds were used only to exercise
conversion and local rules. Those outputs are not submission packages. No example
is evidence of archive sequence association or authenticated Webin acceptance.
Short-CDS evidence requirements can remain after a structurally correct conversion,
especially when choosing `--no-inference`. Do not invent evidence or partial flags
to suppress them.
