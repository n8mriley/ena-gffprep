# Changelog

## 0.6.0 — Optional organelle exclusion

- Add --remove-organelles with direct NCBI sequence-metadata lookup and exact-version matching.
- Add validated caching, explicit refresh and offline JSON/JSONL metadata inputs; no new dependencies.
- Exclude organelle features before conversion and tag allocation, with reference checks and cumulative tag reservation.
- Add exclusion TSVs and preparation-report provenance/counts; retain unclassified primary scaffolds with warnings.
- Support CLI, config, batch and advanced explicit-output workflows; audit-only remains descriptive.
- Preserve the current Annotator / Version / Model description and custom-description handling.
- All 120 offline tests pass. Full Abraxas filtering removes 50 rows / 7 genes and passes the bundled local validator; no authenticated submission performed.


## 0.5.6 — Predictor selection and model provenance

- Infer ANNEVO source format from an explicit ANNEVO predictor in auto mode.
- Recognise unversioned ANNEVO headers without inventing a predictor version.
- Accept spacing/case variants in ANNEVO length comments.
- Add --predictor-model to CLI, config, batch and reports.
- Use accession.predictor.version[.model] default labels and manifest names; append model to descriptions.
- Preserve explicit source-format and naming overrides.
- Add 11 regression methods; all 88 tests pass.


This history covers the GFF preparation tool.

## 0.5.5 — Repository documentation cleanup

- Refresh documentation, comments and diagnostic wording for repository publication.
- Rename the report field to `feature_types_in_supported_subset` and the internal feature set to `SUPPORTED_TYPES`.
- Preserve GFF transformation, validation checks and locus-tag allocation behaviour.
- All 77 automated tests pass; no new authenticated Webin validation or submission performed.

## 0.5.4 — Focused GFF distribution

- Retain GFF preparation, inference and validator-report helpers, and FAI name/length checks.
- Simplify the distribution, examples and validation helper around GFF workflows.
- All 77 retained tests pass; GFF preparation and locus-tag policy are unchanged from 0.5.3.

## 0.5.3 — Locus-tag handling

- Accept existing ASCII alphanumeric suffixes while generating numeric suffixes.
- Exclude repeat_region from propagation and reject explicitly tagged repeats.
- Extend gene-component checks to tmRNA, precursor_RNA and prim_transcript.
- Preserve six-digit generation, historical reservations and global numbering.
- Add locus-tag policy to reports, dedicated guidance and regression tests.

## 0.5.2 — Code documentation

- Add module descriptions and function/class docstrings, including nested helpers.
- Explain side effects, error handling and supported scope.
- Document inference provenance, optional taxonomy and discontinuous CDS IDs.

## 0.5.1 — Optional annotation taxonomy

- Make --taxon-id optional in single, batch and explicit-output preparation.
- Preserve valid input taxonomy when no explicit taxon is supplied.
- Reject conflicting directives, explicit/header mismatches and nonpositive IDs.
- Record requested/effective taxonomy and its source in reports.

## 0.5.0 — Streamlined GFF workflow

- Separate the conversion engine into ena_gffprep_core.py.
- Add predictable package layout, relative manifest paths, summaries and validation instructions.
- Detect ANNEVO provenance from headers and reject conflicting metadata.
- Add JSON defaults and sequential TSV batch preparation with per-row results.
- Record genetic-code fallback and automatic ANNEVO comment lengths.
- Preserve prediction data and retain diagnostic and local-validator helpers.
